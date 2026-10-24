#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
adopt_mmpv.py — MMPV-RNA discovery 产物 → 平台 NCBI 提交项目（一键收养）
=========================================================================

背景
----
MMPV-RNA 的发现流水线（virome_discovery_pipeline）与平台提交模块
（Virus_Platform_Core.ncbi_submit）此前完全解耦：没有共享代码、没有编排脚本，
唯一的衔接是人工把 class_KEEP.fasta 拷过来。本模块补上这一层。

衔接契约（上游产物 → 平台输入）
------------------------------
  09b_Analysis_Verify/class_KEEP.fasta
      → viral_contigs.fasta + sequence_name
        （KEEP 高置信子集；缺失时回退 08_Rescue/HQ_plant_viruses.fasta 全集）
  09_Virome_Analysis/ref_info.tsv
      → taxonomy.tsv + organism
        （12 列：Accession/Length/Species/Genus/Family/Realm/suvtk_Species/
          suvtk_Genus/suvtk_Family/CDS_Count/tRNA_Count/Primary_Tool）
  08_Rescue/checkv/<宿主>/completeness.tsv（可选）→ near_complete
  public_metadata_pipeline 的 Global_Unified_Metadata_Full.tsv（可选）
      → collection_date / geo_loc_name / tissue / cultivar / dev_stage /
        collected_by / sra（按 contig 名内的 Run 号匹配；本地测序序列无 Run 号）
  pipeline_config.yaml（可选）→ cmt-Assembly_Method 提示

产出
----
1. run/tool_runs/contigs_mmpv_<数据集>/     ← 注册为平台"分类运行"
     virus_classification.tsv  contig/species/genus/family/taxon/host/near_complete
     viral_contigs.fasta       class_KEEP.fasta 的副本
   作用：出现在 Web「从运行导入」下拉（store.list_contig_runs 只认 contigs_ 前缀），
         且 store.infer_source_fasta 能自动找到序列 FASTA。
2. run/submissions/<项目名>/
     unified_metadata.csv      平台全部列（UNIFIED_COLUMNS，28 列），能填的都填
     taxonomy.tsv              contig<TAB>organism（回填用，兼作来源留痕）
     discovery_source.json     上游路径/计数留痕
     ＋ store.export_files 生成的全套产物
       （source.src / biosample_template.tsv / miuvig.tsv / assembly.tsv /
         report.html / validation_report.txt / authorset.sbt）

有机体名（organism）走平台的 store._infer_organism，规则与平台其它入口一致：
  species（双名）→ Genus sp. → Family sp.，未培养病毒落到属级符合 NCBI UViG 惯例。

诚实性约定
----------
上游确实没有的字段（作者、你自己的 BioProject/BioSample 等）一律留平台占位符
（YYYY-MM-DD / PRJNAXXXXXX / SAMNXXXXXXXX / Country:Region / XX.XX N XXX.XX E），
绝不编造。`submit-validate` 会精确列出还差什么。

用法
----
  # 单样本（一个数据集共享一套采集信息）
  python main.py submit-adopt --out-dir /path/RNA-Lycium_chinense_out \\
      --name chinense2024 --dataset chinense \\
      --geo-loc China:Ningxia --lat-lon "37.48 N 105.68 E" \\
      --collection-date 2024 --host "Lycium chinense" \\
      --authors "Zhang, Wenda" --bioproject PRJNA123456

  # 多样本（显式样本表按 contig 前缀匹配）
  python main.py submit-adopt --out-dir /path/out --name multi \\
      --samples samples.tsv --bioproject PRJNA123456

  # 只看会生成什么，不落盘
  python main.py submit-adopt --out-dir /path/out --name dry --dry-run
"""

import argparse
import fnmatch
import json
import math
import os
import re
import shutil
import sys
from pathlib import Path

if __package__ in (None, ''):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    from Virus_Platform_Core.ncbi_submit import store
    from Virus_Platform_Core.ncbi_submit.unified_metadata import geo_to_latlon
    from Virus_Platform_Core.config import DIRS
    from Virus_Platform_Core.utils import check_path, safe_open
else:
    from . import store
    from .unified_metadata import geo_to_latlon
    from ..config import DIRS
    from ..utils import check_path, safe_open


# ══════════════════════════════════════════════════════════════
# 上游产物定位与解析
# ══════════════════════════════════════════════════════════════

# 09 分析阶段目录名：新名 09a_（与 00a/00b、03a/03b、09b_ 同规）优先，旧名 09_ 回退。
# 上游把 09_ 改成 09a_ 后两个名字都用同一套解析，新旧输出树都能读。
ANALYSIS_DIR_NAMES = ('09a_Virome_Analysis', '09_Virome_Analysis')

# KEEP 子集（首选）→ HQ 全集 → 最宽的 all_plant_viruses（实测 out5 只有后者，
# 且其中混有参考序列，靠 REF_ACC_RE 剔除）。与 discovery 侧
# virome_pipeline.py:2468 / submission_pipeline.py:162 的探测口径一致。
FASTA_CANDIDATES = [
    ('class_KEEP', '09b_Analysis_Verify/class_KEEP.fasta'),
    ('class_KEEP(旧目录名)', '09b_ACVirus_Analysis/class_KEEP.fasta'),
    ('HQ_plant_viruses', '08_Rescue/HQ_plant_viruses.fasta'),
    ('HQ_plant_viruses', 'HQ_plant_viruses.fasta'),
    ('all_plant_viruses', '08_Rescue/all_plant_viruses.fasta'),
] + [(f'all_plant_viruses({nm})', f'{nm}/all_plant_viruses.fasta')
     for nm in ANALYSIS_DIR_NAMES] \
  + [(f'all_plant_viruses(过滤版,{nm})', f'{nm}/all_plant_viruses_filtered.fasta')
     for nm in ANALYSIS_DIR_NAMES]

REF_INFO_REL = [Path(nm) / 'ref_info.tsv' for nm in ANALYSIS_DIR_NAMES]
PUBLIC_META_REL = [
    Path('metadata_output') / 'Global_Unified_Metadata_Full.tsv',
    Path('metadata_output') / 'Global_Unified_Metadata_Core14.tsv',
]
CONFIG_NAME = 'pipeline_config.yaml'

# CheckV 输出位置：08_Rescue/checkv/<组>/ 与 07_Checkv/<组>/ 都出现过
CHECKV_DIRS = ('08_Rescue/checkv', '07_Checkv')

# NCBI 参考序列 accession（AA_123456 / AA123456.1）。上游 fix_sqn.py:39 用它
# 「剔除参考序列 (NC_/PV_/HM_ 等 accession, 不可提交)」——实测 out5 的
# all_plant_viruses.fasta 里就混着 >OR489165.1，必须同样剔除，否则会把参考
# 序列当成自己的 contig 提交上去。
REF_ACC_RE = re.compile(r'^[A-Z]{2}_?\d+(\.\d+)?$')

TAG_OUR = 'our|'          # cluster_pipeline.py:90 给自家 contig 加的前缀
RUN_ID_RE = re.compile(r'^((?:SRR|ERR|DRR|CRR|SRP|ERP|DRP|PRJ[A-Z]{2})\d+)')

# 平台占位符（store.PLACEHOLDER_RE 认得的那几个）——确未知时写这些，
# 让 submit-validate 精确报出来，绝不编造真值。
# 注意 authors 用空串而不是 'Last, First'：后者不在 store.PLACEHOLDER_RE 内，
# 校验会误判为已填。
PLACEHOLDER = {
    'collection_date': 'YYYY-MM-DD',
    'bioproject': 'PRJNAXXXXXX',
    'biosample': 'SAMNXXXXXXXX',
    'geo_loc': 'Country:Region',
    'lat_lon': 'XX.XX N XXX.XX E',
}

# 组装器名（pipeline_config.yaml 的 assembly.assembler 值）→ source.src 用字符串。
# 只映射认识的名字，未知值原样带过并提示人工确认，不猜版本号。
ASSEMBLER_HINT = {
    'megahit': 'MEGAHIT',
    'rnaviralspades': 'SPAdes rnaviral',
    'spades': 'SPAdes',
    'rnaspades': 'rnaSPAdes',
    'trinity': 'Trinity',
    'idba': 'IDBA',
}

# provenance.json 里 tools 键 → cmt-Assembly_Method 用的名字。
# 实测 out10：{"megahit": "MEGAHIT v1.2.9",
#              "rnaviralspades": "SPAdes genome assembler v4.2.0 [rnaSPAdes mode]"}
PROVENANCE_TOOLS = {
    'megahit': 'megahit', 'rnaviralspades': 'rnaviralspades',
    'spades': 'spades', 'rnaspades': 'rnaspades',
    'trinity': 'trinity', 'idba': 'idba',
}

# suvtk taxonomy 的目录命名变体：点号=submission_pipeline 写入，下划线=virome_analysis
# 写入；09 分析目录新旧名都收（见 ANALYSIS_DIR_NAMES）。
SUVTK_TAX_REL = tuple(
    Path(d) / sub / 'taxonomy.tsv'
    for d, sub in [('08_Rescue', 'suvtk.taxonomy_output'),
                   ('08_Rescue', 'suvtk_taxonomy')]
    + [(nm, 'suvtk.taxonomy_output') for nm in ANALYSIS_DIR_NAMES]
    + [(nm, 'suvtk_taxonomy') for nm in ANALYSIS_DIR_NAMES]
)
PROVENANCE_NAME = 'provenance.json'
RUN_CONFIG_NAME = 'run_config.json'

# CheckV 组目录的优先级：contig 可能同时出现在多个组里，用靠前的组定稿
CHECKV_GROUP_ORDER = ('Plant', 'all', 'no_rescue', 'known', 'rescued')

# ref_info.tsv 的判别列：真 discovery ref_info 至少含其一。
# 用于排除 data-test/ref_info.tsv 那类 4 列参考库表（Accession/Length/Species/Segment）。
REF_INFO_MARKER_COLS = ('Genus', 'Family', 'suvtk_Species', 'suvtk_Genus', 'suvtk_Family')


def _strip_our(seq_id):
    """去掉 cluster_pipeline 的 'our|' 前缀，还原 contig 原名。"""
    s = str(seq_id or '').strip()
    return s[len(TAG_OUR):] if s.startswith(TAG_OUR) else s


def _open_text(path):
    """读文本，容忍 BOM / 非 UTF-8 尾巴（上游 TSV 多为 utf-8-sig）。"""
    return open(path, 'r', encoding='utf-8-sig', errors='replace')


def _read_tsv(path):
    """TSV → list[dict]（表头做 strip；缺列返回空串由调用方兜底）。"""
    try:
        from ..io_utils import read_dicts
    except ImportError:  # 直跑脚本时的兜底
        from Virus_Platform_Core.io_utils import read_dicts
    return read_dicts(path, sep='\t', default=[])


def _read_fasta_ids(path):
    """流式读 FASTA 的表头与长度，返回 [(id, length)]，按文件顺序。

    不依赖 Biopython：平台侧只关心 ID 与长度。
    """
    out = []
    cur_id = None
    cur_len = 0
    with _open_text(path) as f:
        for line in f:
            if line.startswith('>'):
                if cur_id is not None:
                    out.append((cur_id, cur_len))
                head = line[1:].strip()
                cur_id = head.split()[0] if head else ''
                cur_len = 0
            elif cur_id is not None:
                cur_len += len(line.strip())
    if cur_id is not None:
        out.append((cur_id, cur_len))
    return out


def _find_first(root, rels, min_size=100):
    """在 root 及其上两级里依次探测相对路径。"""
    bases = [root, root.parent, root.parent.parent]
    for base in bases:
        for rel in rels:
            p = base / rel
            if p.is_file() and p.stat().st_size > min_size:
                return p, base
    return None, None


def locate_discovery(out_dir):
    """定位 discovery 产物，返回 dict（找不到的项目为 None）。

    out_dir 可以是 $OUT 本身，也可以是 $OUT/08_Rescue（与上游
    submission_pipeline.py 的 --work-dir 习惯一致）。
    """
    root = Path(out_dir).expanduser().resolve()
    if not root.is_dir():
        raise FileNotFoundError(f'discovery 输出目录不存在: {root}')

    fasta = None
    kind = None
    bases = [root, root.parent, root.parent.parent]
    for label, rel in FASTA_CANDIDATES:
        for b in bases:
            p = b / rel
            if p.is_file() and p.stat().st_size > 100:
                fasta, kind = p, label
                base = b
                break
        if fasta is not None:
            break
    if fasta is None:
        raise FileNotFoundError(
            f'在 {root} 及上级找不到待提交 FASTA；请确认这是 discovery 输出根'
            f'（含 09b_Analysis_Verify/ 或 08_Rescue/）')

    ref_info, _ = _find_first(root, REF_INFO_REL)
    pub_meta, _ = _find_first(root, PUBLIC_META_REL)

    # CheckV：<root>/08_Rescue/checkv/<组>/ 或 <root>/07_Checkv/<组>/
    # 组目录可能很多，按 CHECKV_GROUP_ORDER 排序让最贴合的组最后读入、覆盖其它组。
    # 同名组在两种命名下会重复出现（实测 out5 同时有 07_Checkv/Unknown 与
    # 08_Rescue/checkv/unknown、Animal 与 skipped_Animal）→ 按归一化组名去重，
    # 只保留排名靠前的那份，避免重复解析几千行的表。
    checkv = []
    seen_paths = set()
    seen_groups = set()
    for b in [base, root, root.parent]:
        if b is None:
            continue
        for rel in CHECKV_DIRS:
            ck = b / rel
            if not ck.is_dir():
                continue
            found = [p for p in ck.rglob('completeness.tsv') if p.is_file()]

            def _rank(p):
                grp = p.parent.name
                return (CHECKV_GROUP_ORDER.index(grp)
                        if grp in CHECKV_GROUP_ORDER else len(CHECKV_GROUP_ORDER),
                        str(p))

            for p in sorted(found, key=_rank):
                g = _norm_group(p.parent.name)
                if p in seen_paths or g in seen_groups:
                    continue
                seen_paths.add(p)
                seen_groups.add(g)
                checkv.append(p)

    # 富表 info.tsv：09(a)_Virome_Analysis/ 下与 FASTA 同名的 <stem>_info.tsv
    # （实测 out10：HQ_plant_viruses.fasta → HQ_plant_viruses_info.tsv，17 列，
    #  含 Species 双名 / Genus / Family / Segment / Molecule_type，是真实数据里
    #  质量最高的分类来源），以及 All_plant.viruses_info.tsv 作为兜底。
    info_tsv = []
    # 顺序即优先级（后读覆盖先读）：通用兜底先，与 FASTA 同名的那份最后 —— 后者
    # 对应本次要提交的序列集合，比全量集合表更贴切。
    info_names = ['All_plant.viruses_info.tsv', f'{Path(fasta).stem}_info.tsv']
    for b in [base, root, root.parent]:
        if b is None:
            continue
        for nm_dir in ANALYSIS_DIR_NAMES:
            fa_set = b / nm_dir
            for nm in info_names:
                p = fa_set / nm
                if p.is_file() and p not in info_tsv:
                    info_tsv.append(p)

    # suvtk taxonomy.tsv（点号/下划线两种命名都在真实数据里出现过）
    suvtk_tax = None
    for b in [base, root, root.parent]:
        if b is None:
            continue
        for rel in SUVTK_TAX_REL:
            p = b / rel
            if p.is_file() and p.stat().st_size > 20:
                suvtk_tax = p
                break
        if suvtk_tax:
            break

    # pipeline_config.yaml / provenance.json / run_config.json：从输出根向上找
    def _up(name):
        for b in [base, root, root.parent, root.parent.parent]:
            if b is None:
                continue
            c = b / name
            if c.is_file():
                return c
        return None

    return {
        'root': base or root,
        'out_dir': root,
        'fasta': fasta,
        'fasta_kind': kind,
        'ref_info': ref_info,
        'info_tsv': info_tsv,
        'suvtk_tax': suvtk_tax,
        'public_meta': pub_meta,
        'checkv': checkv,
        'config': _up(CONFIG_NAME),
        'provenance': _up(PROVENANCE_NAME),
        'run_config': _up(RUN_CONFIG_NAME),
    }


def read_ref_info(path, log=None):
    """ref_info.tsv → {accession: {...}}；列名容错（大小写/同义词）。

    只接受**真 discovery ref_info**：表头须含 Genus/Family/suvtk_* 之一。
    同名但只有 Accession/Length/Species/Segment 的文件是参考病毒库
    （实测 data-test/ref_info.tsv 15202 行 NC_ 登录号），误读会把参考序列当成
    我们的 contig，这里直接拒绝。
    """
    if not path:
        return {}
    rows = _read_tsv(path)
    if not rows:
        return {}
    header = set(rows[0].keys())
    if not any(c in header for c in REF_INFO_MARKER_COLS):
        if log:
            log(f'[!] {path} 表头 {sorted(header)} 不像 discovery ref_info'
                f'（缺 {"/".join(REF_INFO_MARKER_COLS)} 任一），已忽略')
        return {}
    out = {}
    for r in rows:
        acc = ''
        for k in ('Accession', 'accession', 'contig_id', 'contig', 'seq_name'):
            if r.get(k):
                acc = _strip_our(r[k])
                break
        if not acc:
            continue
        out[acc] = r
    return out


def read_info_tsv(paths, log=None):
    """09_Virome_Analysis/<集合>_info.tsv → {contig: {...}}（真实数据的主分类来源）。

    实测 17 列：Accession Length Species Genus Family Realm Kingdom Class Order
    suvtk_Species suvtk_Genus suvtk_Family CDS_Count Primary_Tool Confidence
    Molecule_type Segment

    比 build_ref_info.py 的 12 列多了 Kingdom/Class/Order、Confidence、
    Molecule_type、Segment —— Species 常是完整双名（Betacytorhabdovirus lycii），
    而 suvtk taxonomy.tsv 往往只到科（Rhabdoviridae sp.）。Segment 更是
    source.src 里原先只能留空的字段。
    """
    out = {}
    for p in paths or []:
        try:
            rows = _read_tsv(p)
        except OSError:
            continue
        n = 0
        for r in rows:
            acc = ''
            for k in ('Accession', 'accession', 'contig_id', 'contig', 'seq_name'):
                if r.get(k):
                    acc = _strip_our(r[k])
                    break
            if not acc:
                continue
            out[acc] = r
            n += 1
        if log and n:
            log(f'[*] 富表分类: {Path(p).name} → {n} 条')
    return out


def read_suvtk_taxonomy(path, log=None):
    """suvtk taxonomy.tsv（contig<TAB>taxonomy）→ {contig: taxonomy}。

    真实数据里这里常是科/属级（'Rhabdoviridae sp.' / 'Potexvirus sp.'），
    作为富表 info.tsv 缺失时的兜底。表头名大小写与别名都容错。
    """
    if not path:
        return {}
    out = {}
    try:
        rows = _read_tsv(path)
    except OSError:
        return {}
    for r in rows:
        cid = t = ''
        for k, v in r.items():
            kl = str(k).lower()
            if kl in ('contig', 'seq_id', 'sequence_id', 'accession') and not cid:
                cid = _strip_our(v)
            elif kl in ('taxonomy', 'tax') and not t:
                t = _last_segment(v)
        if not cid and len(r) >= 2:                      # 无表头：按位置取
            vals = list(r.values())
            cid, t = _strip_our(vals[0]), _last_segment(vals[1])
        if cid:
            out[cid] = t
    if log and out:
        log(f'[*] suvtk taxonomy: {Path(path).name} → {len(out)} 条')
    return out


def read_provenance(path, log=None):
    """provenance.json → 真实工具版本（MMPV 每次运行都会写）。

    实测内容含 {"pipeline": "MMPV-RNA v2.3",
                "tools": {"megahit": "MEGAHIT v1.2.9",
                          "rnaviralspades": "SPAdes genome assembler v4.2.0 [rnaSPAdes mode]",
                          "checkv": "CheckV, version 1.1.1", ...}}
    用它给出 cmt-Assembly_Method / cmt-Annotation_Pipeline 的真实版本串，
    胜过按组装器名猜版本。
    """
    if not path:
        return {}
    try:
        import json
        with _open_text(path) as f:
            data = json.load(f)
    except (OSError, ValueError):
        return {}
    tools = data.get('tools') or {}
    if log and tools:
        log(f'[*] provenance: pipeline={data.get("pipeline") or "?"}, '
            f'{len(tools)} 个工具版本')
    return {
        'pipeline': str(data.get('pipeline') or ''),
        'timestamp': str(data.get('timestamp') or ''),
        'hostname': str(data.get('hostname') or ''),
        'tools': {str(k).lower(): str(v) for k, v in tools.items()},
    }


def assembler_string(provenance, assembler_name, log=None):
    """真实组装器串 → (值, 来源)。

    来源：'provenance'（真实版本，如 'MEGAHIT v1.2.9'）/ 'config'（只有名字）
    / 'provenance(推断)' / ''（都没有）。调用方据此如实报告来源，不误标。
    """
    tools = (provenance or {}).get('tools') or {}
    for key, pkey in PROVENANCE_TOOLS.items():
        if key in tools and (not assembler_name or
                             key in str(assembler_name).lower() or
                             pkey in str(assembler_name).lower()):
            return tools[key], 'provenance'
    if assembler_name:
        return assembler_name, 'config'
    if tools:
        # 名字未知时按常见顺序取第一个命中的组装器，并提示确认
        for key in ('megahit', 'rnaviralspades', 'spades', 'rnaspades'):
            if key in tools:
                if log:
                    log(f'[!] 组装器名未确定，按 provenance 取 {key}：{tools[key]}，请人工确认')
                return tools[key], 'provenance(推断)'
    return '', ''


def _last_segment(val):
    """取分类串的末段（防上游带 'Viruses; Riboviria; ...' 谱系前缀）。"""
    s = str(val or '').strip()
    if ';' in s:
        parts = [p.strip() for p in s.split(';') if p.strip()]
        s = parts[-1] if parts else ''
    return re.sub(r'\s+', ' ', s)


def taxon_triple(tax_row):
    """分类行（富表 info.tsv 或 ref_info.tsv）→ (species, genus, family)。

    优先 Species/Genus/Family（富表里 Species 常是完整双名），回退 suvtk_*。
    与上游 build_ref_info.py:142-150 的「共识 > suvtk」口径一致。
    """
    if not tax_row:
        return '', '', ''
    species = _last_segment(tax_row.get('Species') or tax_row.get('suvtk_Species'))
    genus = _last_segment(tax_row.get('Genus') or tax_row.get('suvtk_Genus'))
    family = _last_segment(tax_row.get('Family') or tax_row.get('suvtk_Family'))
    # suvtk 的最小可定级别是 '<Taxon> sp.'，属级名混在 Species 列时下沉；
    # 但分类等级本身混进 Species 的脏值（实测 suvtk 会输出 'Riboviria sp.'，那是
    # realm）必须先判掉——原实现把"清高等级"放在下沉之后，导致该判不可达，
    # realm 噪音最终漏进 genus。这里合并为一次判断：命中 _HIGH_RANKS 一律清空。
    if species.endswith(' sp.'):
        base = species[:-4].strip()
        if _is_high_rank(base):
            species = ''                     # 高等级噪音：整段清空（不漏到 genus）
        elif not genus:
            genus, species = base, ''        # 属级名下沉策略保持
    # genus/species 列本身命中高等级名（多来源合并时也会出现）→ 一律清空该字段
    if _is_high_rank(genus):
        genus = ''
    if _is_high_rank(species):
        species = ''
    return species, genus, family


# 高于科的分类单元名（这些落到 '<X> sp.' 时是 suvtk 的噪音，实测见过 'Riboviria sp.'）
_HIGH_RANKS = {'riboviria', 'orthornavirae', 'monjiviricetes', 'mononegavirales',
               'viruses', 'duplornaviricota', 'pisuviricota',
               'negarnaviricota', 'kitrinoviricota', 'lenarviricota'}


def _norm_group(name):
    """CheckV 组目录名归一化：out5 上 07_Checkv/Unknown 与 08_Rescue/checkv/unknown
    是同一批数据的两套命名，归一到 'unknown' 以便去重。"""
    s = str(name or '').strip().lower()
    return s[len('skipped_'):] if s.startswith('skipped_') else s


def _is_high_rank(name):
    return str(name or '').strip().lower() in _HIGH_RANKS


def read_checkv(checkv_paths, log=None):
    """CheckV completeness → {contig_id: {'completeness': float|None, 'near': str}}。

    真实 completeness.tsv 的完整度列叫 **aai_completeness**（实测 CheckV 1.1.1，
    15 列：contig_id contig_length viral_length aai_expected_length aai_completeness
    aai_confidence aai_error aai_num_hits aai_top_hit aai_id aai_af
    hmm_completeness_lower hmm_completeness_upper hmm_num_hits kmer_freq），
    且**没有** near_complete 列 —— 早期版本按 'completeness' 找会静默取不到值。
    """
    out = {}
    for p in checkv_paths or []:
        try:
            rows = _read_tsv(p)
        except OSError:
            continue
        n = 0
        for r in rows:
            cid = _strip_our(r.get('contig_id') or r.get('contig') or
                             r.get('seq_name') or '')
            if not cid:
                continue
            raw = (r.get('aai_completeness') or r.get('completeness')
                   or r.get('completeness_score') or '')
            comp_f = None
            try:
                comp_f = float(raw)
                if not math.isfinite(comp_f):
                    comp_f = None
            except (TypeError, ValueError):
                comp_f = None
            nc = (r.get('near_complete') or r.get('is_near_complete') or '').strip()
            if not nc and comp_f is not None:
                nc = 'Yes' if comp_f >= 90 else ('Partial' if comp_f >= 50 else 'No')
            out[cid] = {'completeness': comp_f, 'near': nc}
            n += 1
        if log and n:
            log(f'[*] CheckV: {Path(p).parent.name}/{Path(p).name} → {n} 条')
    return out


def read_run_config(path, log=None):
    """run_config.json → {'profile': 'default'|'plant'|...}。

    MMPV 每次运行都写它（实测 out10: profile=default）。profile 决定该用哪个
    组装器：default→megahit，plant→rnaviralspades。
    """
    if not path:
        return {}
    try:
        import json
        with _open_text(path) as f:
            return {'profile': str(json.load(f).get('profile') or '')}
    except (OSError, ValueError):
        return {}


def read_assembler(config_path, profile='', log=None):
    """从 pipeline_config.yaml 读某 profile 的组装器名。

    优先 PyYAML（平台环境有），按 run_config 的 profile 取 assembly.assembler；
    profile 缺该键时回落 default profile（downstream/analysis 这类只跑下游的
    profile 不设 assembler）。无 PyYAML 时退回正则取第一个 assembler 值。
    """
    if not config_path:
        return ''
    profs = None
    try:
        import yaml
        with _open_text(config_path) as f:
            profs = (yaml.safe_load(f) or {}).get('profiles') or {}
    except Exception:                                        # noqa: BLE001
        profs = None

    def _asm(name):
        if not profs or name not in profs:
            return ''
        return str(((profs[name] or {}).get('assembly') or {})
                   .get('assembler') or '').strip()

    raw = _asm(profile) or _asm('default')
    if not raw:
        try:
            with _open_text(config_path) as f:
                txt = f.read()
        except OSError:
            return ''
        m = re.search(r'^\s*assembler\s*:\s*([A-Za-z0-9_.\-]+)\s*$', txt, re.M)
        raw = m.group(1).strip() if m else ''
    if not raw:
        return ''
    hint = ASSEMBLER_HINT.get(raw.lower())
    if not hint and log:
        log(f'[warn] pipeline_config 的 assembler="{raw}" 未在已知表内，'
            f'请人工确认 cmt-Assembly_Method')
    return hint or raw


def read_software_versions(root, log=None, provenance=None):
    """拼注释流程串（cmt-Annotation_Pipeline / MIUVIG vir_ident_software）。

    优先级：provenance.json 的 pipeline + 真实工具版本 > SOFTWARE_VERSIONS.txt
    的 suvtk 版本 > 中性串。不编造版本号。
    """
    parts = []
    if provenance:
        if provenance.get('pipeline'):
            parts.append(provenance['pipeline'])
        tools = provenance.get('tools') or {}
        for key in ('checkv', 'genomad', 'diamond', 'mmseqs', 'blastn'):
            if tools.get(key):
                parts.append(f"{key}:{tools[key]}")
                break
    if not parts:
        for base in (root, root.parent, root.parent.parent):
            if base is None:
                continue
            p = Path(base) / 'SOFTWARE_VERSIONS.txt'
            if not p.is_file():
                continue
            try:
                with _open_text(p) as f:
                    txt = f.read()
            except OSError:
                continue
            suvtk = re.search(r'suvtk[^\n]*?v?(\d+\.\d+(?:\.\d+)?)', txt, re.I)
            parts.append('MMPV-RNA')
            if suvtk:
                parts.append(f'suvtk v{suvtk.group(1)}')
            break
    if not parts:
        parts.append('MMPV-RNA virome_discovery')
    return ' + '.join(parts)


# ══════════════════════════════════════════════════════════════
# 样本元数据：显式样本表 / 公共元数据 / CLI 兜底
# ══════════════════════════════════════════════════════════════

_SAMPLE_HEADER_ALIASES = {
    'sample_name': 'sample', 'name': 'sample', 'dataset': 'sample',
    'match': 'match', 'contig_prefix': 'match', 'prefix': 'match',
    'pattern': 'match', 'glob': 'match',
    'geo_loc_name': 'geo_loc_name', 'src-geo_loc_name': 'geo_loc_name',
    'location': 'geo_loc_name', 'geo': 'geo_loc_name',
    'lat_lon': 'lat_lon', 'src-lat_lon': 'lat_lon',
    'collection_date': 'collection_date', 'date': 'collection_date',
    'host': 'host', 'src-host': 'host',
    'bioproject': 'bioproject', 'biosample': 'biosample', 'sra': 'sra',
    'tissue': 'tissue', 'src-tissue_type': 'tissue',
    'cultivar': 'cultivar', 'src-cultivar': 'cultivar',
    'dev_stage': 'dev_stage', 'src-dev_stage': 'dev_stage',
    'collected_by': 'collected_by', 'src-collected_by': 'collected_by',
    'isolation_source': 'isolation_source',
    'isolate_prefix': 'isolate_prefix',
    'authors': 'authors', 'title': 'title',
    'assembler': 'assembler', 'sequencer': 'sequencer', 'coverage': 'coverage',
}


def load_sample_map(path):
    """样本表（TSV/CSV）→ list[dict]，表头名容错（平台列名或短名都认）。

    必需列：match（contig 前缀/glob，决定哪些 contig 归这个样本）。
    """
    if not path:
        return []
    rows = []
    sep = '\t' if str(path).lower().endswith(('.tsv', '.txt')) else ','
    with _open_text(path) as f:
        header = None
        for line in f:
            line = line.rstrip('\n').rstrip('\r')
            if not line.strip():
                continue
            cols = [c.strip() for c in line.split(sep)]
            if header is None:
                header = [_SAMPLE_HEADER_ALIASES.get(c.strip().lower(),
                                                     c.strip().lower())
                          for c in cols]
                continue
            d = {k: (cols[i] if i < len(cols) else '') for i, k in enumerate(header)}
            if d.get('match') or d.get('sample'):
                rows.append(d)
    if rows and not any(r.get('match') for r in rows):
        raise ValueError(f'样本表 {path} 缺 match 列（contig 前缀或 glob）')
    return rows


def match_sample(contig, samples):
    """按 'match' 前缀/glob 匹配样本，最长匹配优先；无匹配返回 None。"""
    bare = _strip_our(contig)
    best, best_len = None, -1
    for s in samples:
        pat = str(s.get('match') or '').strip()
        if not pat:
            continue
        if any(ch in pat for ch in '*?['):
            hit = fnmatch.fnmatch(bare, pat)
        else:
            hit = bare.startswith(pat)
        if hit and len(pat) > best_len:
            best, best_len = s, len(pat)
    return best


def load_public_metadata(path):
    """Global_Unified_Metadata_Full.tsv → {run_id: {...}}。

    与平台 unified_metadata.load_metadata 的字段口径保持一致
    （CollectionDate/Location/Tissue/Source/Age_GrowthStage/CenterName/
      Platform/BioProject/BioSample/ScientificName/TaxID）。
    """
    if not path:
        return {}
    out = {}
    try:
        rows = _read_tsv(path)
    except OSError:
        return {}
    for r in rows:
        run = (r.get('Run') or r.get('query_id') or '').strip()
        if not run:
            continue
        out[run] = {
            'collection_date': r.get('CollectionDate', '').strip(),
            'geo_loc_name': r.get('Location', '').strip(),
            'tissue': r.get('Tissue', '').strip(),
            'cultivar': r.get('Source', '').strip(),
            'dev_stage': r.get('Age_GrowthStage', '').strip(),
            'collected_by': r.get('CenterName', '').strip(),
            'platform': r.get('Platform', '').strip(),
            'bioproject': r.get('BioProject', '').strip(),
            'biosample': r.get('BioSample', '').strip(),
            'host': r.get('ScientificName', '').strip(),
            'taxid': r.get('TaxID', '').strip(),
        }
    return out


def extract_run_id(contig):
    """contig 名开头的 Run 号（公共数据才有）。

    比平台 unified_metadata.extract_sra 宽：后者是 re.match(r'([SC]RR\\d+)')，
    只认 SRR/CRR；这里补上 ERR/DRR/SRP/ERP/DRP/PRJ*。
    """
    m = RUN_ID_RE.match(_strip_our(contig))
    return m.group(1) if m else ''


def _first_nonempty(*vals):
    """取第一个非空且不是 'nan'/'not_provided' 的值。"""
    for v in vals:
        s = str(v or '').strip()
        if s and s.lower() not in ('nan', 'not_provided', 'n/a', 'na', 'none'):
            return s
    return ''


# ══════════════════════════════════════════════════════════════
# 注册为平台「分类运行」
# ══════════════════════════════════════════════════════════════

def run_name_for(dataset, run_name=None):
    """运行目录名。必须以 contigs_ 开头才出现在 store.list_contig_runs()。"""
    if run_name:
        base = run_name
    elif dataset:
        base = f'mmpv_{dataset}'
    else:
        base = 'mmpv'
    safe = re.sub(r'[^A-Za-z0-9_\-]', '_', base)
    if not safe.startswith('contigs_'):
        safe = f'contigs_{safe}'
    return safe


def register_run(run_name, contigs, taxonomy, checkv, fasta_path,
                 host_default='', dry_run=False, log=print):
    """把 discovery 产物注册为平台 tool_runs/<run_name>/。

    写 virus_classification.tsv（平台 import_from_run 的输入口径：
    contig + species/genus/family/taxon + host + near_complete + segment）并复制
    FASTA 为 viral_contigs.fasta（store.infer_source_fasta 首选名）。

    taxonomy: {contig: 分类行} —— 由 merge_taxonomy() 合并富表/ref_info/suvtk 后给出。
    checkv:   {contig: {'completeness': float|None, 'near': str}}
    """
    root = DIRS.get('tool_runs')
    run_dir = Path(root) / run_name
    tsv = run_dir / 'virus_classification.tsv'
    fa = run_dir / 'viral_contigs.fasta'
    if dry_run:
        log(f'  [dry-run] 将注册运行 {run_name}：{tsv}')
        log(f'  [dry-run] 将复制 FASTA → {fa}')
        return run_dir, tsv

    run_dir.mkdir(parents=True, exist_ok=True)
    cols = ['contig', 'species', 'genus', 'family', 'taxon', 'host',
            'near_complete', 'completeness', 'length', 'cds', 'trna', 'segment']
    with safe_open(str(tsv), 'wt') as f:
        f.write('\t'.join(cols) + '\n')
        for cid, length in contigs:
            row = taxonomy.get(cid, {})
            species, genus, family = taxon_triple(row)
            ck = checkv.get(cid) or {}
            comp_f = ck.get('completeness')
            f.write('\t'.join(str(v) for v in (
                _strip_our(cid), species, genus, family, '', host_default,
                ck.get('near', ''),
                '' if comp_f is None else f'{comp_f:g}',
                length, row.get('CDS_Count', ''), row.get('tRNA_Count', ''),
                row.get('Segment', ''))) + '\n')

    # FASTA 副本：平台内部的读路径受 check_path 约束，副本留在平台内最稳
    shutil.copyfile(str(fasta_path), str(fa))
    return run_dir, tsv


def merge_taxonomy(info_tsv_map, ref_info_map, suvtk_map, log=print):
    """按质量合并三个分类来源 → {contig: 分类行}，并统计各级覆盖。

    优先级（实测依据）：
      1. 富表 info.tsv      — Species 常为完整双名，另有 Segment/Molecule_type
      2. ref_info.tsv       — R 共识 + suvtk 双栏（12 列）
      3. suvtk taxonomy.tsv — 常只到科/属（'Rhabdoviridae sp.'）
    """
    out = {}
    for cid, t in (suvtk_map or {}).items():
        out[cid] = {'suvtk_Species': t}
    for cid, r in (ref_info_map or {}).items():
        out.setdefault(cid, {}).update(r)
    for cid, r in (info_tsv_map or {}).items():
        out.setdefault(cid, {}).update(r)

    # 富表没给 Species 时，用 suvtk 的科/属级补上（哪怕只到科也比空着强）
    filled = 0
    for cid, row in out.items():
        if not row.get('Species') and row.get('suvtk_Species'):
            row['Species'] = row['suvtk_Species']
            filled += 1
    if log:
        log(f'[*] 合并分类来源: 共 {len(out)} 条'
            f'（富表 {len(info_tsv_map or {})} / ref_info {len(ref_info_map or {})}'
            f' / suvtk {len(suvtk_map or {})}），其中 {filled} 条靠 suvtk 补位')
    return out


# ══════════════════════════════════════════════════════════════
# 平台全列填充（UNIFIED_COLUMNS，28 列）
# ══════════════════════════════════════════════════════════════

def norm_id(s, max_len=45):
    """MMPV 惯用的本地 ID 规范化：非字母数字→下划线；超长截断并加 md5 短尾保唯一。

    与 virome_submission_pipeline/sync_sqn_from_csv.py:52 的 norm_id 同规则，
    避免同一批数据在两侧得到不同的 Isolate 名。
    """
    s = re.sub(r'[^A-Za-z0-9]+', '_', str(s or '')).strip('_')
    if len(s) > max_len:
        import hashlib
        s = s[:max_len - 9] + '_' + hashlib.md5(s.encode()).hexdigest()[:8]
    return s


def build_rows(details, taxonomy, checkv, samples, pub_meta, args, log=print):
    """按平台 UNIFIED_COLUMNS 口径填充每一行的全部列（28）。

    details 来自 store.import_from_run 的返回（已含 organism/sequence_name），
    这里只补齐平台导入不负责的提交元数据列。

    优先级：显式样本表 > 公共元数据（按 Run 号）> CLI > 平台占位符。
    绝不编造：确实未知的一律写平台占位符，由 submit-validate 报出来。
    """
    rows = []
    n_seg = 0
    for d in details:
        cid = str(d.get('contig') or '').strip()
        if not cid:
            continue
        organism = str(d.get('organism') or '').strip() or 'unclassified plant virus'
        near_complete = str(d.get('near_complete') or '').strip()
        ck = checkv.get(cid) or {}
        comp_f = ck.get('completeness')
        # 片段号：真实富表 info.tsv 的 Segment 列（实测 'Unsegmented' / 'RNA' / '3,'）。
        # 非分段病毒留空，分段病毒原样带入 —— 这是原先只能留空的字段。
        seg_raw = str((taxonomy.get(cid) or {}).get('Segment') or '').strip()
        segment = '' if seg_raw.lower() in ('unsegmented', 'na', 'nan', 'none') else seg_raw
        if segment:
            n_seg += 1

        smp = match_sample(cid, samples)
        run_id = extract_run_id(cid)
        pub = pub_meta.get(run_id, {}) if run_id else {}

        def pick(key, cli_val):
            """样本表 > 公共元数据 > CLI。"""
            return _first_nonempty((smp or {}).get(key), pub.get(key), cli_val)

        collection_date = pick('collection_date', args.collection_date)
        geo_loc = pick('geo_loc_name', args.geo_loc)
        host = pick('host', args.host)
        tissue = pick('tissue', args.tissue)
        cultivar = pick('cultivar', args.cultivar)
        dev_stage = pick('dev_stage', args.dev_stage)
        collected_by = pick('collected_by', args.collected_by)
        sra = pick('sra', '') or run_id
        authors = _first_nonempty((smp or {}).get('authors'), args.authors)

        lat_lon = pick('lat_lon', args.lat_lon) or geo_to_latlon(geo_loc) or ''

        # BioProject/BioSample：你自己的提交号优先，公共号只留痕
        ref_bp = pick('bioproject', '')
        bioproject = args.bioproject or ''
        bs_prefix = args.biosample_prefix or ''
        ref_bs = pick('biosample', '')
        if bs_prefix:
            biosample = (f'{bs_prefix}_{sra}' if sra else bs_prefix)[:50]
        else:
            biosample = ref_bs or ''

        # Isolate：同一病毒的所有片段必须同号 → 用 样本+organism 派生，不按行递增
        # （NCBI 要求分段病毒各片段共享同一 Isolate；逐行递增会违反这一点）
        prefix = _first_nonempty((smp or {}).get('isolate_prefix'),
                                 args.isolate_prefix, args.dataset, 'MMPV')
        isolate = norm_id(f'{prefix}_{organism}', max_len=45)

        # Note：留上游来源与公共号，便于 NCBI 审核与自查
        note_parts = []
        if not smp and run_id:
            note_parts.append(f'source_run={run_id}')
        if ref_bp and ref_bp != bioproject:
            note_parts.append(f'original_BioProject={ref_bp}')
        if ref_bs and ref_bs != biosample:
            note_parts.append(f'original_BioSample={ref_bs}')
        if near_complete:
            note_parts.append(f'checkv={near_complete}')
        if comp_f is not None:
            note_parts.append(f'completeness={comp_f:g}%')
        note_parts.append(f'origin=MMPV-RNA discovery/{args.dataset or "-"}')
        note = '; '.join(note_parts)

        se = (smp or {}).get('sequencer') or args.sequencer
        asm = (smp or {}).get('assembler') or args.assembler
        cov = (smp or {}).get('coverage') or args.coverage

        rows.append({
            'organism': organism,
            'sequence_name': cid,
            'authors': authors,     # 未知留空（会占位符告警），不写 'Last, First' 假值
            'collection_date': collection_date or PLACEHOLDER['collection_date'],
            'bioproject': bioproject or PLACEHOLDER['bioproject'],

            'src-Isolate': isolate,
            'src-geo_loc_name': geo_loc or PLACEHOLDER['geo_loc'],
            'src-Lat_Lon': lat_lon or PLACEHOLDER['lat_lon'],
            'src-Host': host,
            'src-Segment': segment,
            'src-Isolation-source': (smp or {}).get('isolation_source') or
                                    args.isolation_source or '',
            'src-Note': note,
            'src-Tissue_type': tissue,
            'src-Collected_by': collected_by,
            'src-Cultivar': cultivar,
            'src-Dev_stage': dev_stage,

            'gb-sample_name': isolate[:50],
            'gb-title': (smp or {}).get('title') or args.title or
                        f'{organism} genome sequencing',
            'sra': sra,
            'biosample': biosample or PLACEHOLDER['biosample'],

            'cmt-Assembly_Method': asm,
            'cmt-Sequencing_Technology': se,
            'cmt-Genome_Coverage': cov,
            'cmt-Annotation_Pipeline': args.annotation_pipeline,

            'bs-isolate': isolate,
            'bs-geo_loc_name': geo_loc or PLACEHOLDER['geo_loc'],
            'bs-host': host,
            'bs-isolation_source': (smp or {}).get('isolation_source') or
                                   args.isolation_source or '',
        })
    if n_seg:
        log(f'[*] src-Segment: {n_seg}/{len(rows)} 条从上游富表带入片段号')
    return rows


# ══════════════════════════════════════════════════════════════
# 服务器拉取（discovery 在 Linux 服务器，平台在本地）
# ══════════════════════════════════════════════════════════════

# 与 MMPV-RNA submission_gui/sync_client.py 的约定一致；可用环境变量覆盖
DEFAULT_SSH_HOST = os.environ.get('MMPV_SSH_HOST', 'zhangwenda@202.119.189.246')
DEFAULT_SERVER_DATA_ROOT = os.environ.get('MMPV_SERVER_DATA_ROOT',
                                          '/home/zhangwenda/data-test')
SERVER_COLLECT_NAME = '_server_collect.py'
# 数据集名/相对路径白名单：避免拼进远端 shell 造成注入
_REMOTE_SAFE_RE = re.compile(r'^[A-Za-z0-9_][A-Za-z0-9_./\-]*$')


def _server_collect_path():
    """服务器侧采集脚本（与本模块同目录，随平台一起分发）。"""
    return Path(__file__).resolve().parent / SERVER_COLLECT_NAME


def _run(cmd, log, timeout=600):
    """跑一条外部命令；返回 (rc, stdout, stderr)。不抛异常，由调用方判断。"""
    import subprocess
    try:
        r = subprocess.run(cmd, capture_output=True, text=True,
                           timeout=timeout, encoding='utf-8', errors='replace')
    except (OSError, subprocess.SubprocessError) as e:
        return 1, '', f'{type(e).__name__}: {e}'
    for line in (r.stdout or '').splitlines():
        if line.strip():
            log(f'      {line.strip()}')
    if r.returncode != 0 and (r.stderr or '').strip():
        log(f'      [stderr] {r.stderr.strip()[-800:]}')
    return r.returncode, r.stdout or '', r.stderr or ''


def _safe_extract(tar_path, dest, max_bytes=2 * 1024 ** 3, max_members=200000):
    """解包 tar.gz，拒绝绝对路径与 .. 穿越（包来自网络，必须校验）。

    同时做 tar 炸弹防护：成员数上限 max_members、解压总量上限 max_bytes
    （默认 2GB）——超限直接拒绝，绝不 extractall 到磁盘。
    """
    import tarfile
    dest = Path(dest).resolve()
    total = 0
    with tarfile.open(tar_path, 'r:gz') as tf:
        members = tf.getmembers()
        if len(members) > max_members:
            raise ValueError(f'包内成员数 {len(members)} 超上限 {max_members}，拒绝解包')
        for m in members:
            if not (m.isfile() or m.isdir()):
                raise ValueError(
                    f'包内含非普通文件/目录成员: {m.name}'
                    f'（符号链接？服务器侧 _server_collect.py 应在打包时解引用）')
            target = (dest / m.name).resolve()
            if target != dest and dest not in target.parents:
                raise ValueError(f'包内路径越界: {m.name}')
            if m.isfile():
                total += max(0, int(m.size or 0))
                if total > max_bytes:
                    raise ValueError(
                        f'包内文件解压总量超上限 {max_bytes / 1024 ** 3:.1f}GB，'
                        f'拒绝解包（疑似 tar 炸弹或交接包过大）')
        tf.extractall(dest)
    return dest


def pull_from_server(args, log=print):
    """ssh/scp 把服务器上的 discovery 数据集拉成本地 staging 目录。

    流程（沿 MMPV-RNA sync_client.py 的既有约定）：
      1. scp _server_collect.py → <host>:/tmp/
      2. ssh 执行它，在服务器 /tmp 打出交接 tar.gz（只含小文件）
      3. scp tar.gz 回本地 <平台>/run/uploads/mmpv_handoff/
      4. 安全解包，返回解包后的数据集目录（供后续走本地适配流程）

    返回解包后的数据集目录；失败抛 RuntimeError。
    """
    collect = _server_collect_path()
    if not collect.is_file():
        raise RuntimeError(f'缺服务器侧脚本 {collect}（应随平台一起分发）')
    if not _REMOTE_SAFE_RE.match(args.from_server or ''):
        raise RuntimeError(f'--from-server 只接受相对数据集名/路径，'
                           f'收到 {args.from_server!r}')

    host = args.server_host or DEFAULT_SSH_HOST
    root = (args.server_data_root or DEFAULT_SERVER_DATA_ROOT).rstrip('/')
    remote_ds = f'{root}/{args.from_server.strip("/")}'
    ds_name = Path(remote_ds).name
    remote_tar = f'/tmp/mmpv_handoff_{ds_name}.tar.gz'
    ssh_base = ['ssh', '-o', 'BatchMode=yes', '-o', 'ConnectTimeout=15', host]
    scp_base = ['scp', '-o', 'BatchMode=yes', '-o', 'ConnectTimeout=15']

    log(f'[*] 服务器: {host}')
    log(f'[*] 远端数据集: {remote_ds}')

    log('[*] [1/4] 上传采集脚本 → /tmp/')
    if _run(['scp', '-o', 'BatchMode=yes', '-o', 'ConnectTimeout=15',
             str(collect), f'{host}:/tmp/{SERVER_COLLECT_NAME}'], log)[0] != 0:
        raise RuntimeError(f'上传 {SERVER_COLLECT_NAME} 失败（检查 ssh 免密是否可用）')

    log('[*] [2/4] 远端打包 discovery 产物')
    rc, out, err = _run(ssh_base + [
        f"python3 /tmp/{SERVER_COLLECT_NAME} '{remote_ds}' -o '{remote_tar}'"], log)
    if rc != 0:
        raise RuntimeError(
            f'远端采集失败（exit={rc}）。上方 [stderr] 是服务器侧报错原文；'
            f'若确属“找不到产物”，检查 {remote_ds} 是否为 discovery 输出根'
            f'（含 08_Rescue/ 或 09b_Analysis_Verify/）')

    from ..config import PLATFORM_ROOT
    staging = Path(DIRS.get('uploads')
                   or os.path.join(PLATFORM_ROOT, 'run', 'uploads')) / 'mmpv_handoff'
    staging.mkdir(parents=True, exist_ok=True)
    local_tar = staging / f'{ds_name}.tar.gz'
    log(f'[*] [3/4] 拉回交接包 → {local_tar}')
    if _run(scp_base + [f'{host}:{remote_tar}', str(local_tar)], log)[0] != 0:
        raise RuntimeError('拉回交接包失败')

    log('[*] [4/4] 解包')
    # 包内结构是 <数据集名>/...，直接解到 staging 下，得到 staging/<数据集名>/...
    ds_dir = staging / ds_name
    if ds_dir.exists():
        shutil.rmtree(ds_dir)
    _safe_extract(local_tar, staging)
    if not ds_dir.is_dir():
        ds_dir = staging          # 包内未套数据集目录时退回解包根
    man = ds_dir / 'manifest.json'
    if man.is_file():
        try:
            info = json.loads(man.read_text(encoding='utf-8'))
            log(f'[*] 清单: {info.get("fasta_sequences")} 条序列, '
                f'{info.get("fasta_kind")}, CheckV {info.get("checkv_groups")} 组, '
                f'{info.get("files") and len(info["files"])} 个文件')
        except (OSError, ValueError):
            pass
    return ds_dir


# ══════════════════════════════════════════════════════════════
# 主流程
# ══════════════════════════════════════════════════════════════

def adopt(args, log=print):
    """执行收养，返回汇总 dict。"""
    if getattr(args, 'from_server', None):
        log('[*] 从服务器拉取 discovery 产物 …')
        args.out_dir = str(pull_from_server(args, log))
    found = locate_discovery(args.out_dir)
    log(f'[*] discovery 输出根: {found["root"]}')
    log(f'[*] KEEP/HQ FASTA: {found["fasta"]}  [{found["fasta_kind"]}]')
    log(f'[*] ref_info.tsv : {found["ref_info"] or "（无）"}')
    log(f'[*] 富表 info.tsv: {len(found["info_tsv"])} 个'
        + (f' — {", ".join(Path(p).name for p in found["info_tsv"])}'
           if found['info_tsv'] else '（无）'))
    log(f'[*] suvtk taxonomy: {found["suvtk_tax"] or "（无）"}')
    log(f'[*] 公共元数据   : {found["public_meta"] or "（无）"}')
    log(f'[*] CheckV       : {len(found["checkv"])} 个 completeness.tsv')
    log(f'[*] 管道配置     : {found["config"] or "（无）"}')
    log(f'[*] provenance   : {found["provenance"] or "（无）"}')

    contigs = _read_fasta_ids(found['fasta'])
    if not contigs:
        raise ValueError(f'{found["fasta"]} 里没有序列')

    # 剔除 NCBI 参考序列：all_plant_viruses.fasta 这类集合会混入参考 accession
    # （实测 out5 含 >OR489165.1）。上游 fix_sqn.py 同样剔除（"不可提交"）。
    refs = [c for c, _ in contigs if REF_ACC_RE.match(_strip_our(c))]
    if refs:
        contigs = [(c, n) for c, n in contigs if not REF_ACC_RE.match(_strip_our(c))]
        log(f'[!] 剔除 {len(refs)} 条参考序列 accession（不可提交）: '
            f'{", ".join(refs[:5])}{" …" if len(refs) > 5 else ""}')

    if args.min_length:
        before = len(contigs)
        contigs = [(c, n) for c, n in contigs if n >= args.min_length]
        log(f'[*] --min-length {args.min_length}: {before} → {len(contigs)} 条')
    if not contigs:
        raise ValueError('剔除参考序列并过滤长度后没有剩余序列可提交')
    log(f'[*] 待提交序列: {len(contigs)} 条')

    # 分类三来源合并：富表 info.tsv > ref_info.tsv > suvtk taxonomy.tsv
    ref_info = read_ref_info(found['ref_info'], log)
    info_map = read_info_tsv(found['info_tsv'], log)
    suvtk_map = read_suvtk_taxonomy(found['suvtk_tax'], log)
    taxonomy = merge_taxonomy(info_map, ref_info, suvtk_map, log)
    provenance = read_provenance(found['provenance'], log)

    checkv = read_checkv(found['checkv'], log)
    samples = load_sample_map(args.samples)
    if samples:
        log(f'[*] 样本表: {len(samples)} 个样本规则')
    pub_meta = load_public_metadata(args.public_metadata or found['public_meta'])
    if pub_meta:
        log(f'[*] 公共元数据索引: {len(pub_meta)} 个 Run')

    # 交叉校验：分类覆盖情况（上游分类表覆盖 HQ 全集，KEEP 是其子集）
    covered = [c for c, _ in contigs if _strip_our(c) in taxonomy]
    if taxonomy and len(covered) < len(contigs):
        miss = [c for c, _ in contigs if _strip_our(c) not in taxonomy]
        log(f'[!] {len(miss)}/{len(contigs)} 条序列无任何分类记录'
            f'（organism 将回落 unclassified）: '
            f'{", ".join(_strip_our(m) for m in miss[:3])}')
    if taxonomy:
        len_mismatch = [c for c, n in contigs
                        if str((taxonomy.get(_strip_our(c)) or {})
                               .get('Length', '')).strip().isdigit()
                        and int(taxonomy[_strip_our(c)]['Length']) != n]
        if len_mismatch:
            log(f'[!] {len(len_mismatch)} 条长度与分类表不一致（分类表可能过期）: '
                f'{", ".join(_strip_our(m) for m in len_mismatch[:3])}')

    # ── 组装默认值 ──
    dataset = args.dataset or Path(found['root']).name
    if not args.dataset:
        log(f'[*] 未给 --dataset，从输出目录名推断: {dataset}')
    args.dataset = dataset          # 供 build_rows 的 Isolate 前缀兜底与 Note 使用
    if not args.assembler:
        run_cfg = read_run_config(found['run_config'], log)
        if run_cfg.get('profile'):
            log(f'[*] 运行 profile: {run_cfg["profile"]}')
        hint = read_assembler(found['config'], run_cfg.get('profile', ''), log)
        if hint:
            log(f'[*] 从 pipeline_config 读到组装器: {hint}')
        # provenance.json 有真实版本串（实测 'MEGAHIT v1.2.9'），优先于裸名字
        real, asm_src = assembler_string(provenance, hint, log)
        if asm_src == 'provenance':
            args.assembler = real
            log(f'[*] cmt-Assembly_Method ← provenance 真实版本: {real}')
        elif asm_src.startswith('provenance'):
            args.assembler = real
            log(f'[*] cmt-Assembly_Method ← provenance: {real}')
        elif hint:
            args.assembler = hint
            log(f'[!] provenance 无版本号，cmt-Assembly_Method 用组装器名 {hint}'
                f'（来自 pipeline_config），版本请人工补（或用 --assembler 指定）')
        else:
            args.assembler = 'SPAdes;4.3.0;rnaviral'   # 平台 submit-init 同款默认
            log('[!] 未能从上游确定组装方法，cmt-Assembly_Method 暂用平台默认值，'
                '请人工确认（或用 --assembler 指定）')
    if not args.sequencer:
        args.sequencer = 'Illumina NovaSeq 6000'
    if not args.annotation_pipeline:
        args.annotation_pipeline = read_software_versions(found['root'], log, provenance)
    # 未显式给 --host 时，尝试样本表/输出目录名里带出的宿主
    if not args.host and samples and len(samples) == 1:
        args.host = samples[0].get('host') or ''

    if not args.host and not samples:
        log('[!] 未给 --host，src-Host 将留空（NCBI 接受，但建议补；'
            '或用 --auto-host 从病毒名推断）')

    # ── 注册平台运行 ──
    # ── 建表（先占位，避免后面注册了运行却在重名处失败留下残留） ──
    if args.dry_run:
        log(f'  [dry-run] 将新建提交项目: {args.name}')
        log('  [dry-run] 未落盘，结束')
        return {'dry_run': True, 'sequences': len(contigs),
                'run': run_name_for(dataset, args.run_name)}

    try:
        store.create_table(args.name)
    except FileExistsError:
        if not args.force:
            raise SystemExit(
                f'提交项目 {args.name} 已存在；换 --name，或加 --force 覆盖（会清空重建）')
        d = store.table_dir(args.name)
        log(f'[!] --force：清空并重建 {d}')
        shutil.rmtree(d)
        store.create_table(args.name)
    log(f'[*] 提交项目: {store.table_dir(args.name)}')

    # ── 注册平台分类运行（建表成功后才落盘） ──
    run_name = run_name_for(dataset, args.run_name)
    log(f'[*] 注册平台分类运行: {run_name}')
    register_run(run_name, contigs, taxonomy, checkv, found['fasta'],
                 host_default=args.host or '', dry_run=False, log=log)

    # ── 导入（organism 交给平台 _infer_organism，口径统一） ──
    imp = store.import_from_run(args.name, run_name)
    log(f'[*] 导入完成: 新增 {imp["added"]} 行，跳过 {imp["skipped"]} 行'
        f'（organism 分级列: {", ".join(imp["lineage_cols"]) or "无"}）')

    rows = build_rows(imp['rows'], taxonomy, checkv, samples, pub_meta, args, log)
    store.save_table(args.name, rows)
    log(f'[*] 全列已写入: {len(rows)} 行 × {len(store.UNIFIED_COLUMNS)} 列')

    # ── taxonomy.tsv 留痕（contig<TAB>organism，可直接被 submit-init 复用） ──
    tax_path = os.path.join(store.table_dir(args.name), 'taxonomy.tsv')
    with safe_open(tax_path, 'wt') as f:
        f.write('contig\ttaxonomy\n')
        for r in rows:
            f.write(f'{r["sequence_name"]}\t{r["organism"]}\n')
    log(f'[*] taxonomy.tsv → {tax_path}')

    # ── 来源留痕 ──
    src_json = os.path.join(store.table_dir(args.name), 'discovery_source.json')
    with safe_open(src_json, 'wt') as f:
        json.dump({
            'discovery_root': str(found['root']),
            'fasta': str(found['fasta']), 'fasta_kind': found['fasta_kind'],
            'ref_info': str(found['ref_info']) if found['ref_info'] else None,
            'info_tsv': [str(p) for p in found['info_tsv']],
            'suvtk_taxonomy': str(found['suvtk_tax']) if found['suvtk_tax'] else None,
            'public_meta': str(found['public_meta']) if found['public_meta'] else None,
            'checkv_files': [str(p) for p in found['checkv']],
            'config': str(found['config']) if found['config'] else None,
            'provenance': str(found['provenance']) if found['provenance'] else None,
            'provenance_pipeline': provenance.get('pipeline') or None,
            'provenance_timestamp': provenance.get('timestamp') or None,
            'provenance_host': provenance.get('hostname') or None,
            'registered_run': run_name,
            'sequences': len(contigs),
            'taxonomy_covered': len(covered),
            'dataset': dataset,
        }, f, ensure_ascii=False, indent=1)

    # ── 平台产物 ──
    if not args.no_export:
        out = store.export_files(args.name, assembler=args.assembler,
                                 sequencer=args.sequencer,
                                 enrichment=args.enrichment, log=log)
        log('[*] 平台产物:')
        for k, v in out.items():
            log(f'      {k:20s} {v}')

    # ── 校验报告 ──
    issues = store.validate_table(args.name)
    if issues:
        log('')
        log(f'[!] 还需人工补全 {len(issues)} 处（submit-validate 同口径）:')
        for it in issues:
            kind = it.get('kind', '')
            n = it.get('count')
            if kind == 'missing_col' or n == -1:
                log(f'      [缺少列] {it["column"]} — {it.get("desc", "")}')
                continue
            ex = ', '.join(repr(e) for e in (it.get('examples') or [])[:2])
            log(f'      [{kind or "问题"} {n} 条] {it["column"]} '
                f'— {it.get("desc", "")}  例: {ex}')
    else:
        log('')
        log('[✓] 所有必填字段已填写，可以提交')

    return {'name': args.name, 'sequences': len(rows), 'run': run_name,
            'issues': len(issues), 'table_dir': store.table_dir(args.name)}


def _build_parser():
    p = argparse.ArgumentParser(
        prog='adopt_mmpv',
        description='MMPV-RNA discovery 产物 → 平台 NCBI 提交项目（class_KEEP.fasta '
                    '+ ref_info.tsv → taxonomy.tsv + 全部列 unified_metadata.csv）',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='例: python main.py submit-adopt --out-dir out/ --name chinense2024 '
               '--dataset chinense --geo-loc China:Ningxia --host "Lycium chinense"')
    p.add_argument('--out-dir',
                   help='discovery 输出根（含 09b_Analysis_Verify/ 或 08_Rescue/），'
                        '也可直接给 08_Rescue/。给了 --from-server 时可省')
    p.add_argument('--from-server', metavar='DATASET',
                   help='从服务器拉取 discovery 产物（相对 <server-data-root> 的'
                        '数据集名，如 out10）；沿 MMPV sync_client.py 的 ssh/scp 约定')
    p.add_argument('--server-host', help=f'ssh 目标（默认 {DEFAULT_SSH_HOST}）')
    p.add_argument('--server-data-root',
                   help=f'服务器数据根（默认 {DEFAULT_SERVER_DATA_ROOT}）')
    p.add_argument('--name', required=True, help='提交项目名（字母/数字/_/-）')
    p.add_argument('--dataset', help='数据集名（默认取输出目录名），用作 Isolate 前缀')
    p.add_argument('--samples', help='样本元数据表（TSV/CSV，需 match 列=contig 前缀/glob）')
    p.add_argument('--run-name', help='平台分类运行名（默认 contigs_mmpv_<dataset>）')
    p.add_argument('--force', action='store_true', help='项目已存在时清表重建')
    p.add_argument('--dry-run', action='store_true', help='只报会生成什么，不落盘')
    p.add_argument('--no-export', action='store_true',
                   help='不跑 store.export_files（只出 unified_metadata.csv）')
    p.add_argument('--min-length', type=int, default=0,
                   help='过滤短于该长度的序列（0=不过滤；class_KEEP 已按上游阈值筛过）')
    p.add_argument('--isolate-prefix', help='Isolate 前缀（默认用 --dataset）')
    p.add_argument('--public-metadata',
                   help='公共元数据表；默认自动探测 <root>/metadata_output/...')

    g = p.add_argument_group('提交元数据（单样本或作为多样本的兜底）')
    g.add_argument('--geo-loc', help='采集地点 "China:Ningxia"')
    g.add_argument('--lat-lon', help='经纬度 "37.48 N 105.68 E"（缺省按 geo-loc 推断）')
    g.add_argument('--collection-date', help='采集日期 YYYY / YYYY-MM / YYYY-MM-DD')
    g.add_argument('--host', help='宿主物种名')
    g.add_argument('--tissue', help='组织类型')
    g.add_argument('--cultivar', help='栽培品种')
    g.add_argument('--dev-stage', help='发育阶段')
    g.add_argument('--collected-by', help='采集人/机构')
    g.add_argument('--isolation-source', help='分离来源')
    g.add_argument('--authors', help='作者 "Last, First; ..."')
    g.add_argument('--title', help='提交标题')
    g.add_argument('--bioproject', help='你自己的 BioProject（PRJNA...）')
    g.add_argument('--biosample-prefix',
                   help='你自己的 BioSample 前缀（SAMN...），按 <前缀>_<Run> 编号')
    g.add_argument('--assembler', help='组装方法（缺省从 pipeline_config 读）')
    g.add_argument('--sequencer', help='测序平台（默认 Illumina NovaSeq 6000）')
    g.add_argument('--coverage', help='基因组覆盖度（如 42.5x）')
    g.add_argument('--annotation-pipeline', help='注释流程串（默认从 SOFTWARE_VERSIONS.txt 拼）')
    g.add_argument('--enrichment', default='rRNA depletion',
                   help='miuvig 富集方式（默认 rRNA depletion）')
    return p


def main(argv=None):
    args = _build_parser().parse_args(argv)
    if not args.out_dir and not args.from_server:
        raise SystemExit('[x] 需要 --out-dir（本地 discovery 输出根）'
                         '或 --from-server（从服务器拉取）')
    if args.samples:
        args.samples = check_path(args.samples, must_exist=True)
    if args.public_metadata:
        args.public_metadata = check_path(args.public_metadata, must_exist=True)
    try:
        res = adopt(args)
    except (FileNotFoundError, ValueError) as e:
        raise SystemExit(f'[x] {e}')
    if res.get('dry_run'):
        print(f'\n[dry-run] 预计 {res["sequences"]} 条序列，运行名 {res["run"]}')
        return 0
    print(f'\n[✓] 提交项目 {res["name"]}：{res["sequences"]} 条序列'
          f'（源运行 {res["run"]}）')
    print(f'    目录: {res["table_dir"]}')
    if res['issues']:
        print(f'    下一步: python main.py submit-validate --name {res["name"]}  '
              f'（还有 {res["issues"]} 处待补）')
    return 0


if __name__ == '__main__':
    sys.exit(main())
