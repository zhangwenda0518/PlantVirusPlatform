#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
_server_collect.py — 服务器侧：把一个 MMPV-RNA discovery 输出目录打成交接包
================================================================================

这是 adopt_mmpv.py 的**服务器侧半边**。平台侧 `main.py submit-adopt --from-server`
会把它 scp 到服务器 /tmp 下执行，再把产出的 tar.gz 拉回本地解包，然后走本地
的适配器流程（注册运行 → 建表 → 导出）。

为什么需要它
------------
discovery 跑在 Linux 服务器上（实测数据根 /home/zhangwenda/data-test/out<N>），
平台在本地 Windows。两者的衔接此前只有人工拷文件。真正需要跨网络搬的只是几
个小文件（FASTA ~百 KB、分类表 ~KB、CheckV 表 ~KB、provenance.json），
而 discovery 侧**知道**这些文件叫什么、放在哪 —— 所以采集逻辑放这一侧。

只依赖标准库（服务器实测 Python 3.10.0，不保证有 pandas）。

用法
----
  python3 _server_collect.py <数据集目录> [-o /tmp/handoff.tar.gz] [--max-mb 200]

  例: python3 _server_collect.py /home/zhangwenda/data-test/out10

产出 tar.gz 内含（相对包根，括号内为真实数据里出现过的等价位置）：
  manifest.json                    本次采集的清单与计数（供下游核对）
  <数据集名>/08_Rescue/HQ_plant_viruses.fasta
      （或 09b_Analysis_Verify/class_KEEP.fasta、all_plant_viruses.fasta）
  <数据集名>/08_Rescue/suvtk.taxonomy_output/taxonomy.tsv（或 suvtk_taxonomy/）
  <数据集名>/08_Rescue/checkv/<组>/completeness.tsv（或 07_Checkv/<组>/，跳过 tmp/）
  <数据集名>/09_Virome_Analysis/<集合>_info.tsv（同上，09a_/09_ 两个名字都收；
      或 All_plant.viruses_info.tsv、ref_info.tsv —— out5 的富表就是后者，21 列）
  <数据集名>/provenance.json / run_config.json / pipeline_config.yaml
      （pipeline_config.yaml 也接受 MMPV 仓库根处的那份）

实测包体积：out10 = 38 KB，out5 = 147 KB。若超出 --max-mb 会报告并拒绝打包。
"""

import argparse
import datetime
import hashlib
import io
import json
import os
import subprocess
import sys
import tarfile
from pathlib import Path

# ── 与 adopt_mmpv.py 保持一致的产物位置约定 ──
# 09 分析阶段目录名：新名 09a_（与 00a/00b、03a/03b、09b_ 同规）优先，旧名 09_ 回退
ANALYSIS_DIR_NAMES = ('09a_Virome_Analysis', '09_Virome_Analysis')

# 候选 FASTA，顺序即优先级：KEEP 子集 > HQ 全集 > 最宽的 all_plant_viruses
# （实测 out10 有 HQ_plant_viruses.fasta，out5 只有 all_plant_viruses.fasta）
FASTA_CANDIDATES = (
    ('class_KEEP', '09b_Analysis_Verify/class_KEEP.fasta'),
    ('class_KEEP', '09b_ACVirus_Analysis/class_KEEP.fasta'),
    ('HQ_plant_viruses', '08_Rescue/HQ_plant_viruses.fasta'),
    ('HQ_plant_viruses', 'HQ_plant_viruses.fasta'),
    ('all_plant_viruses', '08_Rescue/all_plant_viruses.fasta'),
) + tuple(
    (label, f'{nm}/{fn}')
    for nm in ANALYSIS_DIR_NAMES
    for label, fn in (('all_plant_viruses', 'all_plant_viruses.fasta'),
                      ('all_plant_viruses', 'all_plant_viruses_filtered.fasta'))
)

REF_INFO_REL = tuple(f'{nm}/ref_info.tsv' for nm in ANALYSIS_DIR_NAMES)
PUBLIC_META_REL = ('metadata_output/Global_Unified_Metadata_Full.tsv',
                   'metadata_output/Global_Unified_Metadata_Core14.tsv')
SUVTK_TAX_REL = tuple(
    f'{d}/{sub}/taxonomy.tsv'
    for d, sub in [('08_Rescue', 'suvtk.taxonomy_output'),
                   ('08_Rescue', 'suvtk_taxonomy')]
    + [(nm, 'suvtk.taxonomy_output') for nm in ANALYSIS_DIR_NAMES]
    + [(nm, 'suvtk_taxonomy') for nm in ANALYSIS_DIR_NAMES]
)
ROOT_JSONS = ('provenance.json', 'run_config.json', 'pipeline_config.yaml')
# CheckV 输出位置：08_Rescue/checkv/ 或 07_Checkv/（实测 out5 用的是后者）
CHECKV_DIRS = ('08_Rescue/checkv', '07_Checkv')

# 富表 info.tsv：09 分析目录下与待提交 FASTA 同名者，外加固定兜底名。
# 注意 out5 的富表叫 ref_info.tsv（21 列），已由 REF_INFO_REL 覆盖。
INFO_FALLBACK = ('All_plant.viruses_info.tsv',)
# CheckV：只收 completeness.tsv，跳过 tmp/（实测一个组目录里 tmp/ 有 diamond/aai 中间文件）
CHECKV_SKIP_DIRS = ('tmp',)


def _norm_group(name):
    """CheckV 组目录名归一化：out5 上 07_Checkv/Unknown 与 08_Rescue/checkv/unknown
    是同一批数据的两套命名，归一到 'unknown' 以便去重。"""
    s = str(name or '').strip().lower()
    return s[len('skipped_'):] if s.startswith('skipped_') else s


def _first(dataset, rels):
    for rel in rels:
        p = dataset / rel
        if p.is_file() and p.stat().st_size > 20:
            return p
    return None


def collect(dataset):
    """扫描数据集目录，返回 (文件列表, 说明 dict, 待提交 FASTA 路径)。"""
    picked = []
    notes = {}

    fasta = None
    kind = None
    for label, rel in FASTA_CANDIDATES:
        p = dataset / rel
        if p.is_file() and p.stat().st_size > 100:
            fasta, kind = p, label
            break
    notes['fasta_kind'] = kind if fasta else None

    # 富表 info.tsv（真实数据的主分类来源，Species 常为完整双名）
    # 09 分析目录新旧名都看；顺序上通用兜底先、与 FASTA 同名者后（后者更贴切）
    info = []
    for nm_dir in ANALYSIS_DIR_NAMES:
        for nm in INFO_FALLBACK:
            p = dataset / nm_dir / nm
            if p.is_file() and p not in info:
                info.append(p)
    if fasta is not None:
        for nm_dir in ANALYSIS_DIR_NAMES:
            cand = dataset / nm_dir / f'{fasta.stem}_info.tsv'
            if cand.is_file() and cand not in info:
                info.append(cand)

    for p in (fasta, *info):
        if p is not None:
            picked.append((p, p.relative_to(dataset).as_posix()))

    suvtk = _first(dataset, SUVTK_TAX_REL)
    if suvtk:
        picked.append((suvtk, suvtk.relative_to(dataset).as_posix()))

    ref_info = _first(dataset, REF_INFO_REL)
    if ref_info:
        picked.append((ref_info, ref_info.relative_to(dataset).as_posix()))

    pub = _first(dataset, PUBLIC_META_REL)
    if pub:
        picked.append((pub, pub.relative_to(dataset).as_posix()))

    # CheckV：<数据集>/<CHECKV_DIRS>/<组>/completeness.tsv，跳过 tmp/
    # 同组在两种命名下重复出现（实测 out5 有 07_Checkv/Unknown 与
    # 08_Rescue/checkv/unknown）→ 按归一化组名去重，减小包体积
    n_ck = 0
    seen_groups = set()
    for rel in CHECKV_DIRS:
        ck_root = dataset / rel
        if not ck_root.is_dir():
            continue
        for grp in sorted(ck_root.iterdir()):
            if not grp.is_dir() or grp.name in CHECKV_SKIP_DIRS:
                continue
            g = _norm_group(grp.name)
            if g in seen_groups:
                continue
            p = grp / 'completeness.tsv'
            if p.is_file():
                seen_groups.add(g)
                picked.append((p, p.relative_to(dataset).as_posix()))
                n_ck += 1
    notes['checkv_groups'] = n_ck

    for nm in ROOT_JSONS:
        p = dataset / nm
        if p.is_file():
            picked.append((p, nm))

    # pipeline_config.yaml 常不在数据集目录，而在 MMPV 仓库根（实测服务器上如此）。
    # 它给出各 profile 的组装器名，配合 run_config.json 的 profile 才能定
    # cmt-Assembly_Method；没有它只能回落平台默认值（版本号则来自 provenance.json）。
    cfg_in_ds = dataset / 'pipeline_config.yaml'
    if not cfg_in_ds.is_file():
        for cand in (Path.home() / 'MMPV-RNA' / 'pipeline_config.yaml',
                     Path('/home/zhangwenda/MMPV-RNA/pipeline_config.yaml')):
            if cand.is_file():
                picked.append((cand, 'pipeline_config.yaml'))
                break

    # 去重（同一文件可能命中多条规则）
    seen, uniq = set(), []
    for p, rel in picked:
        if p in seen:
            continue
        seen.add(p)
        uniq.append((p, rel))

    notes['files'] = len(uniq)
    notes['total_bytes'] = sum(p.stat().st_size for p, _ in uniq)
    notes['fasta_sequences'] = _count_fasta(fasta) if fasta else 0
    return uniq, notes, fasta


def _count_fasta(path):
    n = 0
    with open(path, 'r', encoding='utf-8', errors='replace') as f:
        for line in f:
            if line.startswith('>'):
                n += 1
    return n


# manifest 哈希：超过该体积的文件不读盘（实测交接包都是小文件，此阈值只为兜底）
_HASH_MAX_BYTES = 500 * 1024 * 1024


def _sha256(path, max_bytes=_HASH_MAX_BYTES):
    """文件 sha256 → (hash, note)。> max_bytes 跳过并标注，读失败也标注。"""
    try:
        size = path.stat().st_size
        if size > max_bytes:
            return '', f'skipped_size>{max_bytes // 1024 ** 2}MB'
        h = hashlib.sha256()
        with open(path, 'rb') as f:
            for chunk in iter(lambda: f.read(1024 * 1024), b''):
                h.update(chunk)
        return h.hexdigest(), ''
    except OSError as e:
        return '', f'error_{type(e).__name__}'


def _git_rev(cands=()):
    """尽力取 git 版本（rev-parse --short HEAD）；非仓库/无 git/失败一律返回 ''。"""
    for d in cands:
        try:
            r = subprocess.run(
                ['git', '-C', str(d), 'rev-parse', '--short', 'HEAD'],
                capture_output=True, text=True, timeout=10)
            if r.returncode == 0 and (r.stdout or '').strip():
                return r.stdout.strip()
        except (OSError, subprocess.SubprocessError):
            continue
    return ''


def main():
    ap = argparse.ArgumentParser(description='打包 MMPV discovery 产物供平台收养')
    ap.add_argument('dataset', help='discovery 数据集目录，如 /home/zhangwenda/data-test/out10')
    ap.add_argument('-o', '--output', default=None,
                    help='输出 tar.gz（默认 /tmp/mmpv_handoff_<数据集名>.tar.gz）')
    ap.add_argument('--max-mb', type=float, default=200.0,
                    help='包体积上限 MB（默认 200，超出则拒绝打包）')
    args = ap.parse_args()

    dataset = Path(args.dataset).expanduser().resolve()
    if not dataset.is_dir():
        raise SystemExit(f'[x] 数据集目录不存在: {dataset}')

    picked, notes, fasta = collect(dataset)
    if fasta is None:
        raise SystemExit(
            f'[x] {dataset} 下找不到 class_KEEP.fasta 或 HQ_plant_viruses.fasta；'
            f'请确认这是 discovery 输出根（含 08_Rescue/ 或 09b_Analysis_Verify/）')
    if not picked:
        raise SystemExit('[x] 没有可打包的文件')

    mb = notes['total_bytes'] / 1024 / 1024
    if mb > args.max_mb:
        raise SystemExit(
            f'[x] 采集内容 {mb:.1f} MB 超过上限 {args.max_mb} MB，已拒绝打包。\n'
            f'    文件清单:\n' +
            '\n'.join(f'      {rel}  {p.stat().st_size/1024:.1f} KB' for p, rel in picked))

    out = Path(args.output) if args.output else \
        Path('/tmp') / f'mmpv_handoff_{dataset.name}.tar.gz'
    ds_name = dataset.name
    with tarfile.open(out, 'w:gz') as tf:
        for p, rel in picked:
            # 解引用符号链接：MMPV 的 ensure_taxonomy 用 os.symlink 把 suvtk 输出
            # 链到 08_Rescue/suvtk.taxonomy_output/，实测 out5 就是这样。若按 symlink
            # 打包，接收端安全解包会（正确地）拒绝 —— 这里存实体内容，入口名不变。
            src = p.resolve()
            if not src.is_file():
                raise SystemExit(f'[x] {p} 解析后不是普通文件: {src}')
            tf.add(str(src), arcname=f'{ds_name}/{rel}')
        man = {
            'dataset': ds_name,
            'source_dir': str(dataset),
            'fasta_kind': notes['fasta_kind'],
            'fasta_sequences': notes['fasta_sequences'],
            'checkv_groups': notes['checkv_groups'],
            'files': [rel for _, rel in picked],
            'total_bytes': notes['total_bytes'],
            'hostname': os.uname().nodename if hasattr(os, 'uname') else '',
            # ── 溯源信息：抓取时间(UTC) / 代码版本 / 每文件 sha256 ──
            'created_utc': datetime.datetime.now(
                datetime.timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ'),
            'git_rev': _git_rev((dataset, Path.home() / 'MMPV-RNA',
                                 Path('/home/zhangwenda/MMPV-RNA'))),
            'file_hashes': [],
        }
        for p, rel in picked:
            _h, _note = _sha256(p)
            entry = {'path': rel, 'bytes': p.stat().st_size, 'sha256': _h}
            if _note:
                entry['note'] = _note           # 大文件跳过/读失败时如实标注
            man['file_hashes'].append(entry)
        raw = json.dumps(man, ensure_ascii=False, indent=1).encode('utf-8')
        ti = tarfile.TarInfo(f'{ds_name}/manifest.json')
        ti.size = len(raw)
        tf.addfile(ti, io.BytesIO(raw))

    print(f'[OK] {out}  ({out.stat().st_size/1024:.1f} KB)')
    print(f'    数据集 {ds_name}: {notes["fasta_sequences"]} 条序列, '
          f'{notes["fasta_kind"]}, CheckV {notes["checkv_groups"]} 组, '
          f'{notes["files"]} 个文件')
    for _, rel in picked:
        print(f'      {rel}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
