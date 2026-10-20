# -*- coding: utf-8 -*-
"""矫正鉴定库 ref_info 的 Segment 值（原位规范化，带备份）。

⚠️ 仅数据库维护者使用（维护者一次性工具，非平台运行链路）：
   默认 --data-dir 指向维护者机器上的 plant_virus_db_pipeline/docs/data
   （用于按 accession 取 family/molecule 上下文）。普通用户运行平台不需要它。

对 databases/virusref_db/{pv_full,pv_complete}/reference.ref_info.tsv：
  - Segment 列原位规范化（307 种原始写法 → 受控词表）；
  - 家族/分子类型上下文从 docs/data 全库 Info 按 accession 精确取；
  - 原文件备份为 .bak_segmentstd。

用法: python scripts/normalize_ref_info.py [--data-dir <docs/data 目录>]
"""

import argparse
import os
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
from Virus_Platform_Core.segment_norm import norm_segment_ctx  # noqa: E402
from Virus_Platform_Core.config import DIRS  # noqa: E402

# 默认 = 维护者机器上的路径；换机/换目录请用 --data-dir 覆盖
# 维护者默认路径可用环境变量 VP_DBPIPE_DATA_DIR 覆盖
DEFAULT_DATA_DIR = os.environ.get('VP_DBPIPE_DATA_DIR') or os.path.join('D:', os.sep, '桌面', 'C-host_classify',
                                'plant_virus_db_pipeline', 'docs', 'data')
LIBS = [
    os.path.join(DIRS['virus_src'], 'pv_full', 'reference.ref_info.tsv'),
    os.path.join(DIRS['virus_src'], 'pv_complete', 'reference.ref_info.tsv'),
]


def build_acc_context(data_dir):
    """accession(去版本号) → (family, molecule)。"""
    # 2026-10-20：polars → pandas（与 kv_identify/kv_filter 同步，去掉
    # polars 依赖）。dtype=str 等价于 polars 的 infer_schema_length=0
    # （全部按字符串读，避免 accession 被推断成数字丢前导零）。
    import pandas as pd
    ctx = {}
    src = os.path.join(data_dir, 'Plant_Virus_Full.Info.tsv')
    if os.path.isfile(src):
        df = pd.read_csv(src, sep='\t', dtype=str, engine='python',
                         on_bad_lines='skip')
        for acc, fam, mol in zip(df['Accession'].tolist(),
                                 df['Family'].tolist(),
                                 df['Molecule_type'].tolist()):
            key = (acc or '').split('.')[0]
            if key and key not in ctx:
                ctx[key] = (fam or '', mol or '')
    return ctx


def main(data_dir: str) -> int:
    # 入口先校验输入目录：默认路径只存在于维护者机器上，普通用户误运行时
    # 应给出明确指引，而不是静默输出 0 条上下文后返回 0。
    if not os.path.isdir(data_dir):
        print(f'✘ 数据目录不存在: {data_dir}')
        print('  本脚本仅供数据库维护者使用；请用 --data-dir 指向'
              ' plant_virus_db_pipeline 的 docs/data 目录')
        return 1
    acc_ctx = build_acc_context(data_dir)
    print(f'accession 上下文条目: {len(acc_ctx)}')
    for path in LIBS:
        if not os.path.isfile(path):
            print('跳过（不存在）:', path)
            continue
        bak = path + '.bak_segmentstd'
        if not os.path.isfile(bak):
            shutil.copy2(path, bak)
        else:
            # 原位改写型脚本的二跑陷阱（2026-10-04 审查 P3）：备份存在但
            # 与当前文件不同 = 当前已是上一轮改写产物，不得覆盖首跑原稿
            import filecmp
            if not filecmp.cmp(path, bak, shallow=False):
                print('  备份保留首跑原稿（当前文件视为已改写产物）:', bak)
        with open(path, encoding='utf-8', errors='replace') as f:
            lines = f.read().split('\n')
        while lines and lines[-1] == '':
            lines.pop()
        header = lines[0].split('\t')
        low = [h.strip().lower() for h in header]
        i_acc = low.index('accession')
        i_seg = low.index('segment')
        out = [lines[0]]
        n_fixed = 0
        for line in lines[1:]:
            cols = line.split('\t')
            if len(cols) <= i_seg:
                out.append(line)
                continue
            acc = cols[i_acc].strip()
            fam, mol = acc_ctx.get(acc.split('.')[0], ('', ''))
            fixed = norm_segment_ctx(cols[i_seg].strip(), fam, mol)
            if fixed != cols[i_seg].strip():
                n_fixed += 1
            cols[i_seg] = fixed
            out.append('\t'.join(cols))
        with open(path, 'w', encoding='utf-8', newline='') as f:
            f.write('\n'.join(out) + '\n')
        name = os.path.basename(os.path.dirname(path))
        print(f'OK {name}: 矫正 {n_fixed} 个 Segment 值')
    return 0


if __name__ == '__main__':
    _p = argparse.ArgumentParser(
        description='原位规范化 ref_info 的 Segment 列（仅数据库维护者使用）')
    _p.add_argument('--data-dir', default=DEFAULT_DATA_DIR,
                    help='维护者机器上的 docs/data 目录'
                         f'（默认: {DEFAULT_DATA_DIR}）')
    _args = _p.parse_args()
    sys.exit(main(_args.data_dir))
