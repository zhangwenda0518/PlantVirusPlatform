# -*- coding: utf-8 -*-
"""
⑥b 层2：HMM / 结构域注释（序列同源层之后的功能兜底）。

- pyhmmer 流式扫描 **Pfam 病毒库**（内存恒定 ~50MB，与库大小无关）：
    Pfam  databases/annot/hmm/pfam/Pfam-A-Viruses.hmm   1,074 profiles
  路径由 Virus_Platform_Core.config.db_path() 解析（见下方 refresh_paths），不硬编码 ——
  数据库重组或自定义数据库根下都能对上。库可缺：available_libs()
  只返回可用的库，缺则整个 HMM 层跳过（不报错）。
  **库可用 = 主文件存在，或 press 产物齐全（`.h3m`+`.h3i`）** —— 后者允许
  删掉 107MB 的纯文本 `.hmm` 只留派生件省空间（见 _hmm_lib_usable）。
  注意入口名必须仍是 `.hmm`（Easel 拿它挂同名 `.h3m`，直接传 `.h3m` 会报错）。
- 命中过滤采用 rosekantor/viral_fams 口径双门槛：
    域级 i-Evalue ≤ 1e-3  且  HMM 模型覆盖率 ≥ 0.5
  hmmscan 调用侧的序列级 E 只作**漏斗**（SCAN_PREFILTER_E=10），不参与
  判定——序列级 E 设成 1e-3 会先滤掉整条序列，远严于域级口径。
  一个命中内的全部合规域都保留（多域蛋白不再只留 best_domain），
  product/category 归属仍取最优域。
- 功能/分类归属：Pfam-A-Viruses 为**自描述 HMM**（NAME/ACC/DESC 内嵌），
  描述直接取命中自带 DESC，无需外部元数据文件。

> 历史：曾支持 VOGDB(vog_all.hmm) / RVDB(rvdb_v32_prot.hmm) 两路，三者都
> 不随平台分发。2026-09-10 起**只保留 Pfam 病毒库**这一路（不再有缺库提示），
> 相关库位、元数据加载与合并分支已一并移除。
> 曾并行支持 CDD 结构域层（mmseqs2 搜索 NCBI Cdd）：因 CDD prefilter 索引
> 一次性分配 ~5GB、小内存机器反复 OOM，2026-09-29 移除；同日改接 **BFVD**
> 病毒域层（bfvd_search.search_orfs，与 SUVTK 提交链同库同命名），
> 层2 合并入口为 merge_layer2(hmm + bfvd)。

对外接口：
  hmm_available()            pyhmmer 是否可用
  available_libs()           实际存在的 HMM 库 [(名, 路径)]
  annotate_orfs_hmm(...)     逐库扫描 → {orf_id: {'hits': [...], 'best': hit}}
  merge_layer2(...)          层2（HMM + BFVD）与层1结果合并 → product/category/family 补全
"""
import os

from .config import get_config, db_path

# i-Evalue / 覆盖率门槛（viral_fams 口径）
DOM_IEVAL_MAX = 1e-3
HMM_COV_MIN = 0.5
#: hmmscan 的**序列级**预过滤 E 值（漏斗作用，不是判定门槛）。
#: 判定门槛始终是上面的域级 DOM_IEVAL_MAX：hmmscan 的 E 是整条序列的
#: E-value，若也设 1e-3，长序列/多域蛋白的合规域会被序列级门槛先滤掉
#: （序列级 E 通常比域级 i-E 大得多），远严于文档口径。
#: 故预过滤放宽到 10，合规判定仍走 :_scan_library 内双门槛（2026-10-03 审计 P3）。
SCAN_PREFILTER_E = 10.0

# ── 库路径（全部走 db_path 注册表；DB_LAYOUT 见 Virus_Platform_Core/config.py）─────────
HMM_DIR = db_path('annot', 'hmm')


# HMM 库「可用」判定：主文件存在，或 press 产物齐全（.h3m+.h3i）即视为可用。
#   —— 为什么要支持「只有 press 产物」：press 后 .hmm 纯文本可删以省空间，
#      但 **HMMFile 的入口名必须仍是 .hmm 主文件名**（Easel 拿它去挂同名
#      .h3m/.h3i，见 _scan_library）；直接传 .h3m 反而 is_pressed=False 报错。
#      故库位保存的是「逻辑入口名」（磁盘上可能只有派生件）。
_PRESS_EXTS = ('.h3m', '.h3i')


def _hmm_lib_usable(entry_path):
    """该 HMM 库位是否可用（满足其一）：
      a) entry_path 本体存在（纯文本 .hmm，未 press 或已 press 的主文件）
      b) entry_path + '.h3m' 且 entry_path + '.h3i' 都存在且非空（仅留 press 产物）
    预压进行中（.h3m 有但 .h3i 空）不算可用——由 _press_in_progress 再拦一层。
    """
    if os.path.isfile(entry_path):
        return True
    return all(os.path.isfile(entry_path + e) and os.path.getsize(entry_path + e) > 0
               for e in _PRESS_EXTS)


# 唯一的 HMM 库：annot/hmm/pfam/Pfam-A-Viruses.hmm（Pfam 病毒谱）
#   —— virsorter2 的病毒 Pfam 子集，1,074 profiles，纯文本 HMMER3/f，
#      **自带 NAME/ACC/DESC**（不需要外部元数据，见 merge_layer2）；
#      实测 2.2s 扫 1,074 profiles、零噪音。
# 目录 2026-09-10 由 hmm/vfam/ 更名为 hmm/pfam/（名字与内容一致）。
# VOGDB(vog_all.hmm) / RVDB(rvdb_v32_prot.hmm) / vFam(vFam-B_2014.hmm) 三路
# 均已移除：库不随平台分发，留着只会产生「缺库」提示。**HMM 层只保留 Pfam**。
#
# 2026-09-10: 支持「只保留 press 产物」——配上 .h3m/.h3i/.h3f/.h3p 后，
#   107MB 的 .hmm 纯文本可删（省空间），库仍判定为可用（见 _hmm_lib_usable）。
PFAM_PATH = os.path.join(HMM_DIR, 'pfam', 'Pfam-A-Viruses.hmm')
# 该库是「自描述」HMM（DESC 内嵌在 HMM 文件里）——Pfam-A-Viruses 属此类
PFAM_SELF_DESC = os.path.basename(PFAM_PATH).startswith('Pfam-A-')

HMM_LIBS = [
    ('pfam', PFAM_PATH),
]
PFAM_ANNOT_DIR = os.path.join(HMM_DIR, 'pfam', 'annot')

_PFAM_META_CACHE = None


def refresh_paths():
    """重算模块级库路径常量（切换「设置 → 数据库目录」后由 config 回调）。

    本模块在 import 时快照了 HMM_DIR 及其派生常量。
    若不重算，`db-migrate` 之后长驻的 Web 进程仍指向旧库位置：
    旧路径已不存在 -> available_libs() 返回空 -> HMM 功能注释
    **静默失效**（不报错，只是结果变空）。故注册进
    Virus_Platform_Core/config.py::_refresh_module_paths()。
    """
    global HMM_DIR, HMM_LIBS
    global PFAM_ANNOT_DIR
    global PFAM_PATH, PFAM_SELF_DESC
    global _PFAM_META_CACHE

    HMM_DIR = db_path('annot', 'hmm')
    PFAM_PATH = os.path.join(HMM_DIR, 'pfam', 'Pfam-A-Viruses.hmm')
    PFAM_SELF_DESC = os.path.basename(PFAM_PATH).startswith('Pfam-A-')
    HMM_LIBS = [
        ('pfam', PFAM_PATH),
    ]
    PFAM_ANNOT_DIR = os.path.join(HMM_DIR, 'pfam', 'annot')
    # 缓存按旧库内容构建，一并失效
    _PFAM_META_CACHE = None


def hmm_available():
    """pyhmmer 可导入即视为 HMM 层可用（库文件缺失时逐库跳过）。"""
    try:
        import pyhmmer  # noqa: F401
        return True
    except ImportError:
        return False


def available_libs():
    """返回 (库名, 路径) 列表（仅可用的库）。

    可用 = 主文件存在，**或** press 产物齐全（.h3m+.h3i）。
    后者支持「删掉纯文本 .hmm、只留 press 产物」以省空间（见 _hmm_lib_usable）。
    """
    return [(n, p) for n, p in HMM_LIBS if _hmm_lib_usable(p)]


def lib_trace():
    """HMM 库可追溯信息 → (db_files, db_versions)。

    清单含 press 派生件（.h3m/.h3i/.h3f/.h3p，可能只有它们存在于磁盘）
    与主文件，文件名+size+mtime；供 ⑥b summary.json 记录本次扫描所用库
    （2026-10-03 审计 P2）。
    """
    from .utils import db_trace
    paths = []
    for _name, p in HMM_LIBS:
        paths.append(p)
        paths += [p + e for e in _PRESS_EXTS + ('.h3f', '.h3p')]
    return db_trace(paths, n_max=40)


# ------------------------------------------------------------------
# Pfam 元数据
# ------------------------------------------------------------------
def _strip(s):
    return (s or '').strip()


def load_pfam_meta(fam_ids):
    """按需读取旧 vFam 的 annot/<fam>_annotations.txt →
    {fam: {'family': 主科, 'genus': 主属, 'desc': 代表产物名}}。
    格式：FAMILIES/GENERA 计数字典 + 成员 FASTA 标题。"""
    import ast
    meta = {}
    for fam in fam_ids:
        path = os.path.join(PFAM_ANNOT_DIR, f'{fam}_annotations.txt')
        d = {}
        try:
            with open(path, encoding='utf-8', errors='replace') as f:
                titles = []
                in_titles = False
                for line in f:
                    line = line.rstrip('\n')
                    if line.startswith('FAMILIES\t'):
                        try:
                            fams = ast.literal_eval(
                                line.split('\t', 1)[1].strip())
                            d['family'] = (max(fams.items(),
                                               key=lambda x: x[1])[0]
                                           if fams else '')
                        except (ValueError, SyntaxError):
                            pass
                    elif line.startswith('GENERA\t'):
                        try:
                            gens = ast.literal_eval(
                                line.split('\t', 1)[1].strip())
                            d['genus'] = (max(gens.items(),
                                              key=lambda x: x[1])[0]
                                          if gens else '')
                        except (ValueError, SyntaxError):
                            pass
                    elif line.startswith('FASTA SEQUENCE TITLES'):
                        in_titles = True
                    elif in_titles and '|' in line:
                        titles.append(line)
                for t in titles:
                    # 标题形如 gi|..|ref|YP_xxx.1|vFam_1000| 产物 [物种]（仅旧 vFam 库需要）
                    body = t.split('|', 5)[-1] if t.count('|') >= 5 else t
                    prod = body.split(' [')[0].strip()
                    low = prod.lower()
                    if prod and 'hypothetical' not in low \
                            and 'unknown' not in low:
                        d['desc'] = prod
                        break
                if 'desc' not in d and titles:
                    first = titles[0].split('|', 5)[-1]
                    d['desc'] = first.split(' [')[0].strip()
        except OSError:
            pass
        meta[fam] = d
    return meta


# ------------------------------------------------------------------
# pyhmmer 扫描
# ------------------------------------------------------------------
_ASCII_LINK = None


def _ascii_path(path):
    """Easel(ANSI fopen) 打不开含中文的路径：经 %TEMP%\\vp_ascii_root 目录
    junction（无需管理员，建一次）把平台内路径转成纯 ASCII 别名。"""
    try:
        path.encode('ascii')
        return path
    except (UnicodeEncodeError, AttributeError):
        pass
    global _ASCII_LINK
    if _ASCII_LINK is None:
        import tempfile
        import subprocess
        from .config import PLATFORM_ROOT
        link = os.path.join(tempfile.gettempdir(), 'vp_ascii_root')
        if not os.path.isdir(link):
            subprocess.run(['cmd', '/c', 'mklink', '/J', link, PLATFORM_ROOT],
                           capture_output=True)
        _ASCII_LINK = link if os.path.isdir(link) else ''
    if not _ASCII_LINK:
        return path
    from .config import PLATFORM_ROOT
    try:
        rel = os.path.relpath(path, PLATFORM_ROOT)
    except ValueError:
        return path
    cand = os.path.join(_ASCII_LINK, rel)
    return cand if str(cand).isascii() else path


def _press_in_progress(lib_path):
    """hmmpress 预压进行中（h3m 已建但索引未写完）。预压期间同文件的其他
    读取会报格式错误，故跳过该库待下次运行自动启用。"""
    h3i = lib_path + '.h3i'
    return (os.path.isfile(lib_path + '.h3m')
            and (not os.path.isfile(h3i) or os.path.getsize(h3i) == 0))


def _scan_library(faa, lib_name, lib_path, threads, logger=None):
    """单库扫描（流式，内存恒定）。返回 {orf_id: [hit, ...]}，
    hit: {lib,target,acc,desc,evalue,score,cov,hmm_from,hmm_to}。"""
    import pyhmmer
    alpha = pyhmmer.easel.Alphabet.amino()
    lib_path = _ascii_path(lib_path)
    faa = _ascii_path(faa)
    if _press_in_progress(lib_path):
        if logger:
            logger.log(f"  [{lib_name}] hmmpress 预压进行中，本次跳过"
                       f"（完成后自动启用，更快）", "WARN")
        return {}
    pressed = os.path.isfile(lib_path + '.h3m') and \
        os.path.getsize(lib_path + '.h3m') > 1024
    by_q = {}
    with pyhmmer.easel.SequenceFile(faa, digital=True, alphabet=alpha) as sf:
        seqs = list(sf)
    if not seqs:
        return by_q
    with pyhmmer.plan7.HMMFile(lib_path) as hf:
        targets = hf.optimized_profiles() if pressed else hf
        for hits in pyhmmer.hmmer.hmmscan(seqs, targets, cpus=threads,
                                          E=SCAN_PREFILTER_E):
            qname = hits.query.name
            qname = qname.decode() if isinstance(qname, bytes) else qname
            found = []
            for hit in hits:
                tname = hit.name
                tname = tname.decode() if isinstance(tname, bytes) else tname
                # 遍历该命中的全部 domain（不再只取 best_domain——多域蛋白
                # 的其余合规域此前被整体丢弃，2026-10-03 审计 P3）；同坐标
                # 重复域去重（保留 i-Evalue 更小者）。
                doms = list(hit.domains or [])
                if hit.best_domain is not None:
                    doms.append(hit.best_domain)
                seen_dom = {}
                for dom in doms:
                    if dom is None:
                        continue
                    ie = dom.i_evalue
                    aln = dom.alignment
                    hmm_len = aln.hmm_length or 0
                    cov = ((aln.hmm_to - aln.hmm_from + 1) / hmm_len
                           if hmm_len else 0.0)
                    # viral_fams 双门槛：域级 i-Evalue + 模型覆盖率
                    if ie > DOM_IEVAL_MAX or cov < HMM_COV_MIN:
                        continue
                    key = (aln.hmm_from, aln.hmm_to)
                    prev = seen_dom.get(key)
                    if prev is not None and prev['evalue'] <= ie:
                        continue
                    seen_dom[key] = {
                        'lib': lib_name, 'target': tname,
                        'acc': (getattr(hit, 'accession', '') or ''),
                        'desc': (getattr(hit, 'description', '') or ''),
                        'evalue': ie, 'score': dom.score,
                        'cov': round(cov, 2),
                        'hmm_from': aln.hmm_from, 'hmm_to': aln.hmm_to,
                    }
                found.extend(seen_dom.values())
            if found:
                found.sort(key=lambda h: h['evalue'])
                by_q[qname] = found[:5]
    if logger:
        n_orf = len(by_q)
        logger.log(f"  [{lib_name}] {n_orf} 个 ORF 获得合规 HMM 命中")
    return by_q


def annotate_orfs_hmm(faa, threads=None, logger=None, progress=None,
                      frac_from=0.0, frac_to=1.0):
    """逐库顺序扫描（当前仅 pfam）。返回
    ({orf_id: {'hits': [...], 'best': hit}}, [库名])。"""
    libs = available_libs()
    if not libs:
        return {}, []
    merged = {}
    done_libs = []
    for i, (lib_name, path) in enumerate(libs):
        if progress:
            progress(frac_from + (frac_to - frac_from) * i / len(libs),
                     f'HMM 扫描 {lib_name}')
        try:
            by_q = _scan_library(faa, lib_name, path,
                                 threads or get_config().threads, logger)
        except Exception as e:
            if logger:
                logger.log(f"  [{lib_name}] 扫描失败: {e}", "WARN")
            continue
        done_libs.append(lib_name)
        for q, hits in by_q.items():
            ent = merged.setdefault(q, {'hits': [], 'best': None})
            ent['hits'].extend(hits)
            cur = ent['best']
            if cur is None or hits[0]['evalue'] < cur['evalue']:
                ent['best'] = hits[0]
        if progress:
            progress(frac_from + (frac_to - frac_from) * (i + 1) / len(libs),
                     f'HMM 扫描 {lib_name} 完成')
    return merged, done_libs


# ------------------------------------------------------------------
# 与层1合并
# ------------------------------------------------------------------
def format_hit_cell(hits, max_targets=5, max_doms_per_target=3):
    """命中列表 → 紧凑单元格文本。
    HMM：PF00680(2e-40,cov0.93)|PF00946(...)；同一模型的多个合规域
    （多域蛋白/重复域）在同一 target 内以 ';' 连接（2026-10-03 审计 P3，
    此前只保留 best_domain）；
    BFVD：target 是 AlphaFold 模型代号（可读性差），用 product 名替代
    （与 suvtk featuretable 的 product 同源）。"""
    order, per_target = [], {}
    for h in (hits or []):
        name = h['target']
        if h.get('lib') == 'bfvd' and h.get('desc'):
            name = h['desc']
        if name not in per_target:
            if len(order) >= max_targets:
                continue
            per_target[name] = []
            order.append(name)
        if len(per_target[name]) < max_doms_per_target:
            per_target[name].append(
                f"{name}({h['evalue']:.0e},cov{h['cov']})")
    return '|'.join(';'.join(per_target[n]) for n in order)


def merge_layer2(rows, hmm_by_q, bfvd_by_q, classify_product):
    """把层2结果并回层1行。返回 (n_hmm_only, n_bfvd_only)。

    优先级：层1 序列命中（informative）> HMM（Pfam 自描述 DESC）> BFVD 域名。
    rows 为 orf_annotation 行 dict（就地更新 product/organism/family/category/
    informative/evidence，并新增 hmm_hits/bfvd_hits 列数据）。"""
    n_hmm_only = n_bfvd_only = 0
    for r in rows:
        q = r['orf_id']
        hh = (hmm_by_q or {}).get(q)
        bb = (bfvd_by_q or {}).get(q)
        r['hmm_hits'] = format_hit_cell(hh['hits']) if hh else ''
        r['bfvd_hits'] = format_hit_cell(bb) if bb else ''
        if r.get('informative') == 'Y':
            r['evidence'] = 'seq'
            continue
        if hh:
            best = hh['best']
            target = best['target']
            # HMM 层当前仅 pfam（Pfam-A-Viruses 自描述，DESC 直接来自 HMM 文件）
            hit_desc = (best.get('desc') or '').strip()
            if hit_desc:
                # 自描述库（Pfam-A-Viruses）：DESC 直接来自 HMM 文件，
                # 干净可读（如 "RNA-dependent RNA polymerase"），
                # 无需 annot/ 外部元数据。
                desc = hit_desc
                meta = {}
            else:
                meta = load_pfam_meta([target]).get(target, {})
                # 兜底标签用命中代号本身（Pfam 库自描述，正常走不到这里）
                desc = meta.get('desc') or target
            family = meta.get('family', '')
            cat = classify_product(desc)
            r['product'] = desc or target
            r['category'] = cat
            if not r.get('family'):
                r['family'] = family
            r['informative'] = 'Y'
            r['evidence'] = 'hmm'
            n_hmm_only += 1
            continue
        if bb:
            best = bb[0]
            name = best.get('desc') or best['target']
            r['product'] = name
            r['organism'] = r['organism'] or best.get('organism', '')
            r['category'] = classify_product(name)
            r['informative'] = 'Y'
            r['evidence'] = 'bfvd'
            n_bfvd_only += 1
            continue
        r['evidence'] = 'none'
    return n_hmm_only, n_bfvd_only
