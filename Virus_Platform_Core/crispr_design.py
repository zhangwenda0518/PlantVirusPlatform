# -*- coding: utf-8 -*-
"""CRISPR/Cas12a 现场检测试剂设计（crRNA + RPA/PCR 引物协同）。

链路：输入多分离物 FASTA（可多分体）→ 按片段分组 → MAFFT 比对 →
逐列保守性 + 多数共识 → TTTV PAM 扫描 → crRNA 漏斗（保守/GC/同聚物/
poly-T/RNAfold MFE）→ 位点侧翼引物设计（primer3 主路 + 枚举兜底）→
双库特异性筛查（宿主基因组 + 近缘病毒；k-mer 索引路 / 大基因组 BLAST 路）
→ crRNA×引物对 复合评分（CASPER 式 21 特征效用加权）→ 群体覆盖率回填 →
输出候选表 / 排序组合 / 订购表 / L1–L5 验证清单 / 参数快照。

规则来源（净室实现，规则与参数不受版权保护、代码未复制，见
docs/CRISPR检测设计模块规划_20261006.md §1）：
- PAM/容错 PAM/保守阈值/RNAfold MFE≥-5/GC 30–70：AutoCORDSv2（GPL-2.0，
  仅参数参考）
- TTTA/TTTC/TTTG 优先、同聚物≥5 剔除、近缘 ≥3 错配口径、RNA 病毒放松档
  （黄病毒先例）：PathoGD, Commun Biol 2025（GPL-3.0，仅论文参数）
- 复合评分特征体系与初始权重：CASPER（iGEM 2025 Lambert GA，无许可证，
  仅公开配置参考；训练数据非植物病毒体系，报告注明）
- RPA primer3 参数（30–35nt / Tm 60–70 / SELF_ANY 8 等）：
  RPA-Primer-Design（MIT）/ CRISPRPrimerDesigner（无许可证）参数口径
- T7 适配：T7 启动子 + DR-spacer 体外转录模板（PathoGD 论文适配方式）

证据边界（ACS Sensors 2026 综述口径，写入每份输出）：
计算设计仅支持候选优先级排序，不构成 crRNA 活性、错配区分度、扩增兼容
性或基质性能的证据；实验结论须按 L1–L5 验证清单逐级建立。

工程约定（对齐 dsrna_pipeline）：
- 全程离线，仅复用已纳管工具：MAFFT（phylo._run_mafft）、primer3
  （primer3_runtime 惰性加载）、ViennaRNA（import RNA，3rd/python 内置）、
  BLAST（大宿主基因组兜底路）；不引入新外部依赖。
- 进程可取消（cancel 置位即抛 RuntimeError）；中文输出统一 utf-8 落盘。
"""
import json
import math
import os
import re
import shutil
import tempfile
import time

import numpy as np

from .config import DIRS, current_host_genome, get_config
from .utils import check_path, iter_fasta, run_cmd, safe_open

ENGINE_VERSION = '1.0.0'

# ---------------------------------------------------------------
# 常量：PAM / 重复序列 / T7
# ---------------------------------------------------------------
# 严格档：TTTV（V=A/C/G，排除 TTTT——LbCas12a 对 TTTT 活性偏低）
_PAM_STRICT_RE = re.compile(r'^TTT[ACG]$')
# 容错档：TTTV 前 3 位容 1 错配（AutoCORDSv2 四组正则的并集）
_PAM_TOLERANT_RE = re.compile(
    r'^(?:[ACG]TT[ACG]|T[ACG]T[ACG]|TT[ACG][ACG]|TTT[ACG])$')

# LbCas12a crRNA 直接重复（handle）：RNA / 对应 DNA 编码
DR_RNA = 'UAAUUUCUACUAAGUGUAGAU'
DR_DNA = 'TAAATTCTACTAAGTGTAGAT'
T7_PROMOTER = 'TAATACGACTCACTATAGGG'      # 与 dsrna_scoring.T7_PROMOTER 一致

# seed 区长度（PAM 近端 spacer 5' 端；Cas12a 错配敏感区）
SEED_LEN = 8


class NucleaseProfile:
    """效应酶参数档案（借鉴 crisprDesign CrisprNuclease 抽象的净室轻量版，
    2026-10-06 二期第一项：引擎逻辑全部从档案取参，新增酶 = 新增档案条目）。

    - requires_pam=False 时（Cas13a，PFS 宽松）：不做 PAM 扫描，保守窗口
      即候选位点，pam 置空、pam_priority 恒 1（排序键中性）；
    - t7_polyt_hard=True（Cas13a 体外转录路线）：spacer 含 TTTT 直接否决
      （T7 转录终止），而非 Cas12a 的仅标记；
    - spacer 长度界由档案决定（Cas12a 20–24 / Cas13a 27–30）。
    """

    def __init__(self, key, label, requires_pam, pam_strict=None,
                 pam_tolerant=None, pfs_note='', spacer_min=20, spacer_max=24,
                 spacer_default=20, dr_rna='', dr_dna='', seed_len=8,
                 t7_polyt_hard=False, target='DNA', note=''):
        self.key = key
        self.label = label
        self.requires_pam = requires_pam
        self.pam_strict = pam_strict          # re 模式字符串（match 全串）
        self.pam_tolerant = pam_tolerant
        self.pfs_note = pfs_note
        self.spacer_min, self.spacer_max = spacer_min, spacer_max
        self.spacer_default = spacer_default
        self.dr_rna, self.dr_dna = dr_rna, dr_dna
        self.seed_len = seed_len
        self.t7_polyt_hard = t7_polyt_hard
        self.target = target
        self.note = note

    def spacer_clamp(self, v):
        try:
            n = int(v)
        except (TypeError, ValueError):
            n = self.spacer_default
        if n <= 0:                       # 0/None/空 = 用档案默认
            n = self.spacer_default
        return max(self.spacer_min, min(self.spacer_max, n))


# LwaCas13a crRNA 直接重复（36 nt，Kellner et al. 2017 Nat Protoc
# SHERLOCK 口径；Lsh/LwaCas13a 共用此 DR 序列）
_CAS13A_DR_DNA = 'GATTTAGACTACCCCAAAAACGAAGGGGACTAAAAC'
_CAS13A_DR_RNA = _CAS13A_DR_DNA.replace('T', 'U')

NUCLEASES = {
    'cas12a': NucleaseProfile(
        key='cas12a', label='LbCas12a', requires_pam=True,
        pam_strict=r'TTT[ACG]',
        pam_tolerant=(r'(?:[ACG]TT[ACG]|T[ACG]T[ACG]|TT[ACG][ACG]|TTT[ACG])'),
        pfs_note='TTTV PAM（严格档）；容错档前 3 位容 1 错配，活性需实验验证',
        spacer_min=20, spacer_max=24, spacer_default=20,
        dr_rna=DR_RNA, dr_dna=DR_DNA, seed_len=8, t7_polyt_hard=False,
        target='DNA', note='trans-cleavage 诊断主力；DR-spacer 化学合成或 T7 转录皆可'),
    'cas13a': NucleaseProfile(
        key='cas13a', label='LwaCas13a', requires_pam=False,
        pfs_note='无 PAM；PFS 宽松（LwaCas13a 报告无强 3′ 偏好）',
        spacer_min=27, spacer_max=30, spacer_default=30,
        dr_rna=_CAS13A_DR_RNA, dr_dna=_CAS13A_DR_DNA, seed_len=8,
        t7_polyt_hard=True, target='RNA',
        note='直接识别 RNA（免逆转录路线）；crRNA 一律 T7 体外转录，'
             'spacer 过 poly-T≥4 硬检查（转录终止）'),
}

# ---------------------------------------------------------------
# 参数档位（严格 / 放松；放松档来自 PathoGD 黄病毒实测先例）
# ---------------------------------------------------------------
PROFILES = {
    'strict': {
        'min_cons_crrna': 0.95,      # crRNA spacer 平均保守分
        'min_cons_primer': 0.95,     # 引物区平均保守分
        'seed_cons_min': 1.0,        # 种子区逐列保守下限（PAM 近端
                                     # 8nt 错配容忍极低，硬约束）
        'product_min': 100,
        'product_max': 200,
        'pam_tolerant': False,       # 仅 TTTV
        'max_sites': 30,             # 每片段进入引物设计的 crRNA 位点数
    },
    'relaxed': {
        'min_cons_crrna': 0.85,
        'min_cons_primer': 0.80,
        'seed_cons_min': 0.95,
        'product_min': 90,
        'product_max': 300,
        'pam_tolerant': True,        # 前 3 位容 1 错配（需实验验证活性）
        'max_sites': 40,
    },
}

# RPA primer3 参数（净室口径：MangoBank / CRISPRPrimerDesigner 参数并集）
RPA_P3_GLOBAL = {
    'PRIMER_OPT_SIZE': 31, 'PRIMER_MIN_SIZE': 30, 'PRIMER_MAX_SIZE': 35,
    'PRIMER_OPT_TM': 65.0, 'PRIMER_MIN_TM': 60.0, 'PRIMER_MAX_TM': 70.0,
    'PRIMER_MIN_GC': 30.0, 'PRIMER_MAX_GC': 70.0,
    'PRIMER_GC_CLAMP': 1,
    'PRIMER_MAX_POLY_X': 4,
    'PRIMER_MAX_SELF_ANY': 8, 'PRIMER_MAX_SELF_END': 4,
    'PRIMER_PAIR_MAX_COMPL_ANY': 8, 'PRIMER_PAIR_MAX_COMPL_END': 4,
    'PRIMER_MAX_HAIRPIN_TH': 20.0,
    'PRIMER_MAX_NS_ACCEPTED': 0,
    'PRIMER_SALT_MONOVALENT': 50.0, 'PRIMER_DNA_CONC': 50.0,
    'PRIMER_NUM_RETURN': 5,
}
# PCR 模式沿用平台引物设计口径（primer_design.CORE_DEFAULTS 附近值，
# 保持与 /primer 页可比），仅放开产物上限到 RPA 同区间由档位控制
PCR_P3_GLOBAL = {
    'PRIMER_OPT_SIZE': 20, 'PRIMER_MIN_SIZE': 18, 'PRIMER_MAX_SIZE': 25,
    'PRIMER_OPT_TM': 60.0, 'PRIMER_MIN_TM': 57.0, 'PRIMER_MAX_TM': 63.0,
    'PRIMER_MIN_GC': 40.0, 'PRIMER_MAX_GC': 60.0,
    'PRIMER_MAX_POLY_X': 4,
    'PRIMER_MAX_SELF_ANY': 8, 'PRIMER_MAX_SELF_END': 4,
    'PRIMER_PAIR_MAX_COMPL_ANY': 8, 'PRIMER_PAIR_MAX_COMPL_END': 4,
    'PRIMER_MAX_NS_ACCEPTED': 0,
    'PRIMER_SALT_MONOVALENT': 50.0, 'PRIMER_DNA_CONC': 50.0,
    'PRIMER_NUM_RETURN': 5,
}

# 宿主 k-mer 索引路的 FASTA 体量上限（MB）。
# 2026-10-07 真实番茄宿主（810MB）实测：np.unique 峰值约 bp×40B（双链
# 编码 + 排序去重），810MB → ~40GB，64GB 机在可用内存不足时必然分配
# 失败。250MB → ~10GB 峰值留足余量；更大的宿主一律走 BLAST 路（磁盘
# 库、内存友好），k-mer 构建失败也自动降级 BLAST，绝不静默放弃筛查
HOST_KMER_CAP_MB = 250
# 特异性筛查的 k-mer 长度（20 为 Cas12a spacer 主流长度）
SCREEN_K = 20


class CrisprDesignError(RuntimeError):
    """输入/配置错误（对用户可读，任务侧直接呈现）。"""


# ---------------------------------------------------------------
# ViennaRNA（惰性加载；失败时结构过滤自动降级为跳过，不阻塞设计）
# ---------------------------------------------------------------
_RNA_MOD = None
_RNA_TRIED = False


def _rna():
    global _RNA_MOD, _RNA_TRIED
    if not _RNA_TRIED:
        _RNA_TRIED = True
        try:
            import RNA                      # noqa: F401  ViennaRNA 2.x
            _RNA_MOD = RNA
        except Exception:
            _RNA_MOD = None
    return _RNA_MOD


def rna_mfe(seq_rna):
    """RNA 序列 MFE（kcal/mol）。不可用时返回 None（调用方按未过滤处理）。"""
    rna = _rna()
    if rna is None or not seq_rna:
        return None
    try:
        md = rna.md()
        fc = rna.fold_compound(seq_rna, md)
        struct, mfe = fc.mfe()
        return float(mfe), str(struct)
    except Exception:
        return None


def dna_mfe(seq):
    """DNA 参数集的 MFE（引物兜底路用）。不可用时返回 None。"""
    rna = _rna()
    if rna is None or not seq:
        return None
    try:
        md = rna.md()
        md.dna = True
        fc = rna.fold_compound(seq, md)
        struct, mfe = fc.mfe()
        return float(mfe)
    except Exception:
        return None


# ---------------------------------------------------------------
# 输入解析与片段分组
# ---------------------------------------------------------------
# 常见多分体命名：RNA1/RNA-2/RNA_3/segment 4/S1/M2/L1（含 DNA/A/B 分段）。
# 锚点用字母负向后顾而不是 \b：'A_RNA1' 里 '_' 与 'R' 都是单词字符，
# \b 在 R 前不成立（实测漏分组）；(?<![A-Za-z]) 允许 '_RNA1'。
# 分段标签两种命名（2026-10-07 真实数据教训）：
#   数字式：RNA1 / RNA-2 / segment RNA3 / S1 / M2 / L1 / DNA-A
#   字母式：RNA L / RNA M / RNA S（负义 RNA 病毒 TSWV/INSV 等分段命名）
# 分段标签（两种命名，2026-10-07 真实数据两轮教训）：
#   数字式 RNA1 / RNA-2 / segment RNA3；字母式 RNA L / RNA M /
#   RNA S（TSWV 等负义分段）与 DNA-A/B（双生病毒）。
#   字母/数字必须紧跟 RNA|DNA|SEG(MENT) 前缀——裸单字母分支会把
#   isolate 编号（A-25/Far-25）误判为分段（ToMV 真实数据实测误伤）
# 分段标签（两种命名，2026-10-07 真实数据两轮教训）：
#   数字式 RNA1 / RNA-2 / segment RNA3；字母式 RNA L / RNA M / RNA S
#   （TSWV 等负义分段）与 DNA-A/B（双生病毒）。
#   字母/数字必须紧跟 RNA|DNA|SEG(MENT) 前缀——裸单字母分支会把
#   isolate 编号（A-25/Far-25）误判为分段（ToMV 真实数据实测误伤）
_SEG_RE = re.compile(
    r'(?<![A-Za-z])(?:(?:RNA|DNA|SEG(?:MENT)?)[\s\-_?]'
    r'?([0-9]{1,2}|[IVX]{1,4}|[LSMAB]))(?![0-9A-Za-z])',
    re.IGNORECASE)


def read_targets(path):
    """读目标 FASTA → [(name, seq)]；同名记录自动加后缀防撞。

    name 保留**完整 header**（截 80 字符）：分段标签（segment RNA1 /
    RNA L 等）常在真实 NCBI 标题中段而非首词 accession，只取首词会让
    多分体分组失效（2026-10-07 真实 CMV 冒烟教训：24 条全进 main 组）。
    MSA 写盘仍按 first word 取唯一名（_write_tmp_fasta / _read_aln 兜底
    按位置对齐），报告里的分离株名反而更可读。
    """
    recs = []
    for h, s in iter_fasta(check_path(path, must_exist=True)):
        name = (h or 'seq').strip()[:80] or 'seq'
        s = re.sub(r'[^A-Za-z]', '', (s or '').upper())
        if s:
            recs.append((name, s))
    if not recs:
        raise CrisprDesignError('目标 FASTA 中没有可用序列')
    seen = {}
    out = []
    for name, s in recs:
        if name in seen:
            seen[name] += 1
            out.append((f'{name}_{seen[name]}', s))
        else:
            seen[name] = 0
            out.append((name, s))
    return out


def group_segments(records, group_by_segment=True):
    """按 header 中的片段标签分组；无标签/关闭分组 → 全部进 'main'。

    组键 = 匹配到的完整标签（大写、去分隔符），如 RNA1 / RNA-1 / RNA_1
    → 'segRNA1'；A1 分段命名 → 'segA1'（与 RNA1 天然不混组）。
    返回 {seg_key: [(name, seq), ...]}，保持输入顺序。
    """
    groups = {}
    for name, seq in records:
        key = 'main'
        if group_by_segment:
            m = _SEG_RE.search(name)
            if m:
                label = re.sub(r'[\s\-_]', '', m.group(0)).upper()
                label = label.replace('SEGMENT', 'SEG')
                key = 'seg' + label
        groups.setdefault(key, []).append((name, seq))
    return groups


# ---------------------------------------------------------------
# MSA → 逐列保守性 + 多数共识
# ---------------------------------------------------------------
def _write_tmp_fasta(records, path):
    with safe_open(path, 'wt') as f:
        for name, seq in records:
            f.write('>%s\n%s\n' % (name, seq))


def _read_aln(path, names):
    aln = []
    for h, s in iter_fasta(path):
        first = h.split()[0] if h.split() else ''
        aln.append((first, re.sub(r'[^A-Za-z\-.]', '', (s or '').upper())
                    .replace('.', '-')))
    # MAFFT --inputorder 保证顺序；按输入名兜底对齐（防头被截断）
    if [n for n, _ in aln] != names:
        if len(aln) == len(names):
            aln = [(names[i], s) for (__, s), i in zip(aln, range(len(names)))]
        else:                          # 罕见：解析异常，按名字映射
            by = {n: s for n, s in aln}
            aln = [(n, by.get(n, '')) for n in names]
    return [(n, s) for n, s in aln if s]


def consensus_and_conservation(records, work_dir, threads=None, logger=None,
                               cancel=None):
    """返回 (consensus, cons_list, n_used, msa_used)。

    单条记录：保守性恒 1.0；多条：MAFFT 比对后逐列统计。
    cons_list[k] = 共识第 k 位（0 基）对应的保守分 ∈ [0,1]。
    """
    if len(records) == 1:
        seq = records[0][1]
        return seq, [1.0] * len(seq), 1, False

    from .phylo import _run_mafft
    in_fa = os.path.join(work_dir, 'in.fa')
    aln_fa = os.path.join(work_dir, 'aln.fa')
    _write_tmp_fasta(records, in_fa)
    if cancel is not None and cancel.is_set():
        raise RuntimeError('用户取消')
    _run_mafft(in_fa, aln_fa, threads=threads, logger=logger)
    if cancel is not None and cancel.is_set():
        raise RuntimeError('用户取消')
    aln = _read_aln(aln_fa, [n for n, _ in records])
    if not aln:
        raise CrisprDesignError('MAFFT 比对结果为空')

    n = len(aln)
    cons, cols = [], []
    n_ambig = 0
    for col in zip(*(s for _, s in aln)):
        bases = [c for c in col if c in 'ACGT']
        if not bases:
            cons.append('N')
            cols.append(0.0)
            continue
        counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
        for b in bases:
            counts[b] += 1
        best, cnt = max(counts.items(), key=lambda kv: (kv[1], kv[0]))
        if cnt * 2 > len(bases):
            cons.append(best)          # 严格多数才定碱基
        else:
            # 平票不猜：2 碱基平票记 IUPAC 兼并码（W/S/M/K/R/Y，保住引物
            # 模板区），≥3 碱基分歧记 N。crRNA/引物序列的纯 ACGT 约束由
            # 候选扫描与引物后置守卫另行把关
            tied = sorted(b for b, c in counts.items() if c == cnt)
            cons.append(_IUPAC_PAIRS.get(''.join(tied), 'N')
                        if len(tied) == 2 else 'N')
            n_ambig += 1
        cols.append(cnt / float(len(bases)))
    if logger and n_ambig:
        # 评审 C3：让「候选窗骤缩」可归因——区分真不保守与平票技术性 N
        logger.log('共识含 %d 个非多数列（平票记 IUPAC / ≥3 分歧记 N）：'
                   '落在这些列上的 PAM/引物扫描会被跳过' % n_ambig, 'WARN')
    return ''.join(cons), cols, n, True


# ---------------------------------------------------------------
# 序列特征小工具
# ---------------------------------------------------------------
def max_homo_run(seq):
    """最长同碱基连续 run（如 GGGG → 4）。"""
    best = run = 0
    prev = ''
    for c in seq:
        run = run + 1 if c == prev else 1
        prev = c
        best = max(best, run)
    return best


def gc_frac(seq):
    if not seq:
        return 0.0
    gc = sum(1 for c in seq if c in 'GC')
    acgt = sum(1 for c in seq if c in 'ACGT')
    return gc / acgt if acgt else 0.0


# IUPAC 兼并码（平票共识）+ 完整互补表（负链窗口含兼并码也要正确翻链）
_IUPAC_PAIRS = {'AT': 'W', 'CG': 'S', 'AC': 'M', 'GT': 'K',
                'AG': 'R', 'CT': 'Y'}
_RC = str.maketrans('ACGTURYSWKMBDHVN', 'TGCAAYRSWMKVHDBN')
_PURE_ACGT = set('ACGT')


def revcomp(seq):
    return seq.translate(_RC)[::-1]


def complementarity_pct(a, b):
    """a 与 revcomp(b) 全部错位重叠里最高 Watson-Crick 配对比例（%）。

    二聚体倾向的粗口径（0–100）：比例越高越容易二聚。纯 Python 滑动，
    引物 30–35nt 规模下每次调用 <1ms。
    """
    rb = revcomp(b)
    best = 0.0
    for off in range(-(len(rb) - 1), len(a)):
        s0 = max(0, -off)
        e0 = min(len(a), len(rb) - off)
        if e0 - s0 < 8:                       # <8bp 重叠不构成二聚体核
            continue
        matched = 0
        for i in range(s0, e0):
            if a[i] + rb[i + off] in ('AT', 'TA', 'GC', 'CG'):
                matched += 1
        best = max(best, matched / float(e0 - s0))
    return round(best * 100.0, 1)


def three_prime_gc(seq):
    """3′ 端最后 3 碱基全 G/C（等温扩增规范要求规避，防非特异延伸）。"""
    return len(seq) >= 3 and set(seq[-3:]) <= {'G', 'C'}


def three_prime_self_run(seq, tail=12):
    """3' 末端同碱基连续 run（CASPER 口径：自最末碱基起数）。"""
    s = seq[-tail:]
    if not s:
        return 0
    last = s[-1]
    run = 0
    for ch in reversed(s):
        if ch == last:
            run += 1
        else:
            break
    return run


def cross_3p_run(fp, rp, tail=8):
    """fp 与 rp 的 3'/3' 交叉配对 run（CASPER 口径）。

    max k 使 fp[-k:] == revcomp(rp[-tail:])[:k] —— 即两条引物 3' 端
    可互相退火延伸的连续配对长度。
    """
    fp_tail = (fp or '')[-tail:].upper()
    rp_tail_rc = revcomp((rp or '')[-tail:])
    best = 0
    for k in range(1, tail + 1):
        if fp_tail[-k:] == rp_tail_rc[:k]:
            best = k
    return best


# ---------------------------------------------------------------
# crRNA 候选生成与漏斗
# ---------------------------------------------------------------
def find_pam_sites(consensus, cons, spacer_len, tolerant, profile):
    """按酶档案扫描候选位点（Cas12a：PAM+spacer；Cas13a：保守窗口即位点）。

    site_start/site_end = 共识**正链** 0 基闭区间（PAM+spacer 或 spacer）；
    负链位点 spacer 取自负链 5'→3'，保守分按负链坐标取值后映射。
    """
    L = len(consensus)
    out = []
    if profile.requires_pam:
        pam_re = _pam_regex(profile, tolerant)
        # 正链：PAM [i,i+4)，spacer [i+4, i+4+cl)
        for i in range(0, L - 4 - spacer_len + 1):
            pam = consensus[i:i + 4]
            if not pam_re.match(pam):
                continue
            spacer = consensus[i + 4:i + 4 + spacer_len]
            if set(spacer) - _PURE_ACGT:
                continue
            out.append(_make_site(spacer, '+', pam, i, cons, spacer_len,
                                  profile))
        # 负链：在 revcomp 上扫（cons 翻转），坐标映射回正链
        rseq = revcomp(consensus)
        rcons = cons[::-1]
        for j in range(0, L - 4 - spacer_len + 1):
            pam = rseq[j:j + 4]
            if not pam_re.match(pam):
                continue
            spacer = rseq[j + 4:j + 4 + spacer_len]
            if set(spacer) - _PURE_ACGT:
                continue
            # rseq[j : j+4+cl] ↔ 正链区间 [L-(j+4+cl), L-j)
            out.append(_make_site(spacer, '-', pam, L - (j + 4 + spacer_len),
                                  rcons, spacer_len, profile,
                                  spacer_start=j + 4))
    else:
        # Cas13a：无 PAM，保守窗口即候选（每位置一个窗口，漏斗过滤）
        for i in range(0, L - spacer_len + 1):
            spacer = consensus[i:i + spacer_len]
            if set(spacer) - _PURE_ACGT:
                continue
            out.append(_make_site(spacer, '+', '', i, cons, spacer_len,
                                  profile))
        rseq = revcomp(consensus)
        rcons = cons[::-1]
        for j in range(0, L - spacer_len + 1):
            spacer = rseq[j:j + spacer_len]
            if set(spacer) - _PURE_ACGT:
                continue
            out.append(_make_site(spacer, '-', '', L - (j + spacer_len),
                                  rcons, spacer_len, profile,
                                  spacer_start=j))
    # 去重（同 spacer 同位点）
    seen = set()
    dedup = []
    for c in sorted(out, key=lambda d: -d['cons_score']):
        key = (c['spacer'], c['site_start'], c['strand'])
        if key not in seen:
            seen.add(key)
            dedup.append(c)
    return dedup


def _pam_regex(profile, tolerant):
    import re as _re
    pat = profile.pam_tolerant if tolerant else profile.pam_strict
    return _re.compile('^(?:%s)$' % pat)


def _make_site(spacer, strand, pam, site_start, strand_cons, spacer_len,
               profile, spacer_start=None):
    """组装候选。spacer_start = spacer 在**该链**坐标里的起点（0 基）；
    正链 = site_start+PAM 长，负链由调用方按 revcomp 坐标显式给出。"""
    if spacer_start is None:
        spacer_start = site_start + (4 if profile.requires_pam else 0)
    cs = strand_cons[spacer_start:spacer_start + spacer_len]
    cons_score = (sum(cs) / len(cs)) if cs else 0.0
    if not profile.requires_pam:
        pri = 1                                  # 无 PAM：排序键中性
    else:
        strict = bool(_PAM_STRICT_RE.match(pam))
        if strict and pam != 'TTTT':
            pri = 2                              # TTTA/TTTC/TTTG 最优
        elif strict:
            pri = 1                              # TTTT 次选（活性偏低）
        else:
            pri = 0                              # 容错 PAM：需实验验证
    sl = strand_cons[spacer_start:spacer_start + profile.seed_len]
    return {
        'spacer': spacer, 'strand': strand, 'pam': pam,
        'pam_priority': pri, 'cons_score': round(cons_score, 4),
        'seed': spacer[:profile.seed_len],
        'seed_cons': round(min(sl), 4) if sl else 0.0,
        'site_start': site_start,
        'site_end': site_start
        + (4 if profile.requires_pam else 0) + spacer_len - 1,
    }


_DR_MFE_CACHE = {}


def _dr_mfe(nuclease):
    """DR 骨架单独折叠 MFE（按酶档案缓存一次）。"""
    if nuclease.key not in _DR_MFE_CACHE:
        f = rna_mfe(nuclease.dr_rna)
        _DR_MFE_CACHE[nuclease.key] = f[0] if f else None
    return _DR_MFE_CACHE[nuclease.key]


def crrna_funnel(candidates, profile, params, nuclease=None, logger=None,
                 cancel=None):
    """crRNA 漏斗：保守 / GC / 同聚物 / poly-T / RNAfold MFE。

    profile = 参数档（PROFILES 条目）；nuclease = NucleaseProfile。
    Cas12a：poly-T≥4 仅标记；Cas13a（t7_polyt_hard）：直接否决。
    就地给候选补 mfe / t7_unsafe 标记；返回通过列表与漏斗统计。
    """
    nuc = nuclease or NUCLEASES['cas12a']
    stats = {'total': len(candidates), 'by_cons': 0, 'by_gc': 0,
             'by_homopolymer': 0, 'by_seed': 0, 'by_mfe': 0,
             'by_polyt': 0, 'by_dr_pair': 0}
    kept = []
    min_cons = float(profile['min_cons_crrna'])
    for c in candidates:
        if cancel is not None and cancel.is_set():
            raise RuntimeError('用户取消')
        sp = c['spacer']
        if c['cons_score'] < min_cons:
            stats['by_cons'] += 1
            continue
        # 种子区（PAM 近端 seed_len 位）逐列保守硬约束：错配容忍极低，
        # 不允许群体 SNV 落在种子区
        seed_thr = float(profile.get('seed_cons_min', 0.99))
        if c.get('seed_cons', 1.0) < seed_thr:
            stats['by_seed'] += 1
            continue
        gc = gc_frac(sp)
        if not (0.30 <= gc <= 0.70):
            stats['by_gc'] += 1
            continue
        if max_homo_run(sp) >= 5:
            stats['by_homopolymer'] += 1
            continue
        c['gc'] = round(gc * 100.0, 1)
        # poly-T ≥4：T7 体外转录终止风险。Cas12a 标记不剔除（订购表警示）；
        # Cas13a 一律体外转录 → 硬否决
        c['t7_unsafe'] = 'TTTT' in sp
        if nuc.t7_polyt_hard and c['t7_unsafe']:
            stats['by_polyt'] += 1
            continue
        # RNAfold：spacer 作为 RNA 折叠，MFE < -5 kcal/mol 剔除
        fold = rna_mfe(sp.replace('T', 'U'))
        if fold is None:
            c['mfe'] = None            # ViennaRNA 不可用：跳过结构过滤
        else:
            c['mfe'] = round(fold[0], 2)
            c['seed_unpaired_frac'] = round(
                sum(1 for ch in fold[1][:nuc.seed_len] if ch == '.')
                / float(nuc.seed_len), 3)
            if fold[0] < -5.0:
                stats['by_mfe'] += 1
                continue
            # DR+spacer 联合折叠（ΔMFE 闯门）：DR 自身带设计性茎环，全长
            # 绝对阈值会误杀；spacer 与 DR 骨架强配对/自折叠时全长 MFE
            # 显著低于 DR 单独折叠 → 判 DR 相互作用，剔除
            dr_mfe = _dr_mfe(nuc)
            full = rna_mfe(nuc.dr_rna + sp.replace('T', 'U'))
            if full is not None and dr_mfe is not None:
                c['crna_mfe'] = round(full[0], 2)
                if full[0] < dr_mfe - 3.0:
                    stats['by_dr_pair'] += 1
                    continue
        c.setdefault('seed_unpaired_frac', None)
        kept.append(c)
    stats['kept'] = len(kept)
    if logger:
        logger.log('crRNA 漏斗: %d 个候选位点 → 保守-%d 种子-%d GC-%d '
                   '同聚物-%d polyT-%d MFE-%d DR互作-%d → 通过 %d'
                   % (stats['total'], stats['by_cons'], stats['by_seed'],
                      stats['by_gc'], stats['by_homopolymer'],
                      stats['by_polyt'], stats['by_mfe'],
                      stats['by_dr_pair'], stats['kept']))
    # 排序：PAM 优先级 → 保守性 → MFE（越不脆越好）
    kept.sort(key=lambda d: (-d['pam_priority'], -d['cons_score'],
                             99.0 if d['mfe'] is None else d['mfe']))
    return kept, stats


# ---------------------------------------------------------------
# 引物设计（primer3 主路 + 枚举兜底）
# ---------------------------------------------------------------
def _p3():
    from .primer3_runtime import get_module
    try:
        return get_module()
    except RuntimeError as e:
        raise RuntimeError('primer3 不可用：CRISPR 引物设计需要 primer3（%s）'
                           % e)


def _amp_global(amp_mode, profile):
    base = dict(RPA_P3_GLOBAL if amp_mode == 'RPA' else PCR_P3_GLOBAL)
    base['PRIMER_PRODUCT_SIZE_RANGE'] = [
        [int(profile['product_min']), int(profile['product_max'])]]
    return base


def design_primers_for_site(consensus, cons, site, amp_mode, profile,
                            params, work_dir, logger=None, cancel=None):
    """对一个 crRNA 位点设计引物对（扩增子完整含 PAM+spacer，两侧留边距）。

    返回 [{'f_seq','f_start','f_len','f_tm','f_gc','r_seq','r_start','r_len',
           'r_tm','r_gc','product','route'}]（坐标 = 共识 1 基）。
    """
    spacer_len = len(site['spacer'])
    flank = int(params.get('flank') or 250)
    margin = int(params.get('site_margin') or 25)
    L = len(consensus)
    w0 = max(0, site['site_start'] - flank)
    w1 = min(L, site['site_end'] + 1 + flank)
    window = consensus[w0:w1]
    if len(window) < int(profile['product_min']) + 40:
        return []
    # 位点在窗口内的坐标（0 基 → primer3 亦 0 基起）
    s_in = site['site_start'] - w0
    e_in = site['site_end'] - w0

    pairs = []
    # ---- 路线 A：primer3 ----
    try:
        opts = _amp_global(amp_mode, profile)
        opts['PRIMER_PRODUCT_SIZE_RANGE'] = [[
            max(int(profile['product_min']), 80),
            min(int(profile['product_max']), len(window) - 1)]]
        # 生产修复（2026-10-07 真实 PVY 冒烟抓到）：primer3 模板只认
        # ACGTN，IUPAC 兼并码（M/K/R/Y…平票共识）会整路报错降级。把
        # 兼并位转 N：模板区长度保住，PRIMER_MAX_NS_ACCEPTED=0 自动
        # 阻止引物落在该位（兼并码的本意），坐标不变
        window_p3 = re.sub(r'[^ACGTN]', 'N', window)
        res = _p3().bindings.design_primers(
            {'SEQUENCE_ID': 'win', 'SEQUENCE_TEMPLATE': window_p3}, opts)
        for i in range(res.get('PRIMER_PAIR_NUM_RETURNED', 0)):
            f_pos = res.get('PRIMER_LEFT_%d' % i, (0, 0))
            r_pos = res.get('PRIMER_RIGHT_%d' % i, (0, 0))
            # 模板可含 IUPAC 兼并码（平票共识），但引物交付序列必须纯 ACGT
            fseq_a = res.get('PRIMER_LEFT_%d_SEQUENCE' % i, '')
            rseq_a = res.get('PRIMER_RIGHT_%d_SEQUENCE' % i, '')
            if (set(fseq_a) | set(rseq_a)) - _PURE_ACGT:
                continue
            if three_prime_gc(fseq_a) or three_prime_gc(rseq_a):
                continue
            amp_a = f_pos[0]
            amp_z = r_pos[0]                    # primer3: 正链最右碱基（0 基）
            if amp_a + margin > s_in or e_in + margin > amp_z:
                continue                        # 扩增子须完整含位点且留边距
            pairs.append({
                'f_seq': res.get('PRIMER_LEFT_%d_SEQUENCE' % i, ''),
                'f_start': f_pos[0] + 1 + w0, 'f_len': f_pos[1],
                'f_tm': round(res.get('PRIMER_LEFT_%d_TM' % i, 0.0), 1),
                'f_gc': round(res.get('PRIMER_LEFT_%d_GC_PERCENT' % i, 0.0), 1),
                'r_seq': res.get('PRIMER_RIGHT_%d_SEQUENCE' % i, ''),
                'r_start': r_pos[0] + 1 + w0, 'r_len': r_pos[1],
                'r_tm': round(res.get('PRIMER_RIGHT_%d_TM' % i, 0.0), 1),
                'r_gc': round(res.get('PRIMER_RIGHT_%d_GC_PERCENT' % i, 0.0), 1),
                'product': res.get('PRIMER_PAIR_%d_PRODUCT_SIZE' % i, 0),
                'route': 'primer3',
            })
    except (OSError, RuntimeError) as e:
        if logger:
            logger.log('primer3 主路失败（%s），转枚举兜底' % e, 'WARN')

    # ---- 路线 B：枚举兜底（AutoCORDSv2 语义的规则过滤）----
    if len(pairs) < int(params.get('min_pairs_per_site') or 2):
        pairs = _enumerate_primers(window, cons[w0:w1], s_in, e_in,
                                   amp_mode, profile, margin, pairs, w0=w0)
    return pairs


def _enumerate_primers(window, wcons, s_in, e_in, amp_mode, profile, margin,
                       existing, w0=0):
    """滑动枚举 30–35nt（PCR 18–25nt）候选 + 规则过滤 + 简单配对。

    输出坐标与路线 A 同口径：共识 1 基全局坐标（w0 = 窗口在共识上的偏移）。
    """
    lo, hi = ((30, 35) if amp_mode == 'RPA' else (18, 25))
    pmin, pmax = int(profile['product_min']), int(profile['product_max'])
    min_cons = float(profile['min_cons_primer'])
    f_cands, r_cands = [], []

    def _ok_primer(seq, c0):
        # 模板（window/consensus）可含 IUPAC 兼并码（平票共识，见 build_consensus），
        # 但**引物交付序列必须纯 ACGT** —— 与路线 A（primer3 分支）和探针设计
        # 同口径。2026-10-07 审查 P0-2：原实现只挡字面 'N'，W/S/M/K/R/Y 全部放行，
        # 路线 B 会把含兼并码的引物一路送到 ranked_sets / order_sheet（无法合成）。
        if set(seq) - _PURE_ACGT:
            return False
        gc = gc_frac(seq)
        if not (0.30 <= gc <= 0.70):
            return False
        if max_homo_run(seq) >= 5 or three_prime_gc(seq):
            return False
        csl = wcons[c0:c0 + len(seq)]
        if not csl or sum(csl) / len(csl) < min_cons:
            return False
        mfe = dna_mfe(seq)
        if mfe is not None and mfe < -5.0:
            return False
        return True

    # 正向：3' 端止于 s_in - margin（保守取 f1 ≤ f_end_zone）
    f_end_zone = s_in - margin
    for f1 in range(max(0, f_end_zone - 60), f_end_zone + 1):
        for flen in range(lo, hi + 1):
            f0 = f1 - flen + 1 - 1 + 1 - 1      # f1 为 3' 末碱基下标
            f0 = f1 - flen + 1
            if f0 < 0:
                continue
            seq = window[f0:f1 + 1]
            if _ok_primer(seq, f0):
                f_cands.append((f0, f1, seq))
    # 反向：5' 端起于 e_in + 1 + margin
    r_start_zone = e_in + 1 + margin
    for r0 in range(r_start_zone, min(len(window), r_start_zone + 60)):
        for rlen in range(lo, hi + 1):
            r1 = r0 + rlen - 1
            if r1 >= len(window):
                continue
            seq = revcomp(window[r0:r1 + 1])
            if _ok_primer(seq, r0):
                r_cands.append((r0, r1, seq))
    if not f_cands or not r_cands:
        return existing

    p3 = None
    try:
        p3 = _p3()
    except RuntimeError:
        pass

    def _tm(seq):
        if p3 is not None:
            try:
                return round(float(p3.bindings.calc_tm(seq)), 1)
            except Exception:
                pass
        # Wallace 兜底
        return round(2.0 * (seq.count('A') + seq.count('T'))
                     + 4.0 * (seq.count('G') + seq.count('C')), 1)

    out = list(existing)
    seen = {(p['f_start'], p['product']) for p in existing}
    scored = []
    for f0, f1, fseq in f_cands:
        for r0, r1, rseq in r_cands:
            product = r1 - f0 + 1
            if not (pmin <= product <= pmax):
                continue
            ftm, rtm = _tm(fseq), _tm(rseq)
            if abs(ftm - rtm) > 10:
                continue
            key = (f0 + 1, product)
            if key in seen:
                continue
            scored.append((abs(ftm - rtm), f0, f1, fseq, ftm, r0, r1, rseq,
                           rtm, product))
    scored.sort(key=lambda t: t[0])
    for t in scored[:5]:
        _, f0, f1, fseq, ftm, r0, r1, rseq, rtm, product = t
        out.append({
            'f_seq': fseq, 'f_start': f0 + 1 + w0, 'f_len': f1 - f0 + 1,
            'f_tm': ftm, 'f_gc': round(gc_frac(fseq) * 100.0, 1),
            'r_seq': rseq, 'r_start': r0 + 1 + w0, 'r_len': r1 - r0 + 1,
            'r_tm': rtm, 'r_gc': round(gc_frac(rseq) * 100.0, 1),
            'product': product, 'route': 'enumerate',
        })
    return out


# ---------------------------------------------------------------
# 特异性双库筛查
# ---------------------------------------------------------------
def _fasta_mb(path):
    try:
        return os.path.getsize(path) / 1048576.0
    except OSError:
        return 0.0


def build_exclusion_index(name, fasta, k, cache_dir, logger=None, cancel=None):
    """建 k-mer 排序数组索引（复用 dsrna_offtarget，双链入集 + 磁盘缓存）。

    返回 {'name','arr','bp','n_records','route':'kmer'}；失败/超限返回
    {'name','skipped':原因}。
    """
    from .dsrna_offtarget import build_species_index
    try:
        st = build_species_index(fasta, k=k, cache_dir=cache_dir,
                                 log=logger, cancel=cancel)
    except Exception as e:
        return {'name': name, 'skipped': f'索引构建失败: {e}'}
    return {'name': name, 'arr': st['arr'], 'bp': st['bp'],
            'n_records': st['n_records'], 'route': 'kmer',
            'cached': bool(st.get('cached'))}


def screen_windows(seq, index, k=SCREEN_K, max_mm=2):
    """序列全部 k-mer 窗口对单个索引的筛查。

    返回 {'min_tier': 0/1/2/None, 'n_mm0','n_mm1','n_mm2',
          'three_prime_mm1': bool}（3' 端 15nt 内是否存在 mm1 窗口）。
    """
    from .dsrna_offtarget import _best_tier, encode_kmers
    arr = index['arr']
    if arr is None or arr.size == 0 or len(seq) < k:
        return {'min_tier': None, 'n_mm0': 0, 'n_mm1': 0, 'n_mm2': 0,
                'three_prime_mm1': False}
    kmers = encode_kmers(seq, k)
    if kmers.size == 0:
        return {'min_tier': None, 'n_mm0': 0, 'n_mm1': 0, 'n_mm2': 0,
                'three_prime_mm1': False}
    counts = {0: 0, 1: 0, 2: 0}
    min_tier = None
    tp_mm1 = False
    # 窗口起点 ≥ len(seq)-k-15+1 → 与 3' 端 15nt 重叠
    tp_zone = max(0, len(seq) - k - 14)
    for i in range(int(kmers.size)):
        t = _best_tier(int(kmers[i]), arr, k, max_mm)
        if t is None:
            continue
        counts[t] += 1
        if min_tier is None or t < min_tier:
            min_tier = t
        if t <= 1 and i >= tp_zone:
            tp_mm1 = True
    return {'min_tier': min_tier, 'n_mm0': counts[0], 'n_mm1': counts[1],
            'n_mm2': counts[2], 'three_prime_mm1': tp_mm1}


def screen_crrna_candidate(cand, indexes, k=SCREEN_K):
    """crRNA 对全部排除库：任一库 mm0 → veto；mm1 → flag。"""
    res = {}
    veto = False
    worst = None
    for name, idx in indexes.items():
        r = screen_windows(cand['spacer'], idx, k=k)
        res[name] = r
        if r['min_tier'] == 0:
            veto = True
        if r['min_tier'] is not None and (worst is None or r['min_tier'] < worst):
            worst = r['min_tier']
    return {'veto': veto, 'worst_tier': worst, 'per_lib': res}


def screen_primer_pair(pair, indexes, k=SCREEN_K):
    """引物对筛查：任一引物 mm0 → veto；3' 端 mm1 → 重罚 flag。"""
    res = {}
    veto = False
    worst = None
    for tag in ('f', 'r'):
        res[tag] = {}
        for name, idx in indexes.items():
            r = screen_windows(pair[tag + '_seq'], idx, k=k)
            res[tag][name] = r
            if r['min_tier'] == 0:
                veto = True
            if r['min_tier'] is not None and (worst is None
                                              or r['min_tier'] < worst):
                worst = r['min_tier']
    return {'veto': veto, 'worst_tier': worst, 'per_primer': res}


def blast_screen_large(seqs_by_id, db_prefix, work_dir, logger=None):
    """BLAST 兜底路（超大宿主基因组）：blastn-short，命中≈mm≤1 判近缘。

    返回 {query_id: {'min_tier': 0|1|None, 'hits': [(sseqid, mismatch)]}}。
    精度口径：查询覆盖 ≥ len(query)-1 且错配 ≤1 → tier0/1。同一
    (query, subject) 的多条 HSP 先按查询坐标合并再判（评审 C2：单个
    错配可能把短 query 拆成两条 HSP，各自覆盖不足会被整行过滤 → 近缘
    命中漏报，方向危险）。重叠 HSP 的错配有重复计入的可能，偏保守。
    """
    cfg = get_config()
    blastn = cfg.tool('blastn')
    makeblastdb = cfg.tool('makeblastdb')
    import glob as _g
    if not _g.glob(db_prefix + '.n??'):
        if logger:
            logger.log('构建宿主 BLAST 库（仅首次，大基因组需较长时间）…',
                       'WARN')
        host_fa = os.path.splitext(db_prefix)[0] + '.fa'
        run_cmd([makeblastdb, '-in', host_fa, '-dbtype', 'nucl',
                 '-out', db_prefix], logger=logger)
    q_fa = os.path.join(work_dir, 'screen_query.fa')
    with safe_open(q_fa, 'wt') as f:
        for qid, seq in seqs_by_id.items():
            f.write('>%s\n%s\n' % (qid, seq))
    out_tsv = os.path.join(work_dir, 'screen_hits.tsv')
    run_cmd([blastn, '-task', 'blastn-short', '-query', q_fa,
             '-db', db_prefix, '-word_size', '7', '-dust', 'no',
             '-evalue', '100', '-max_target_seqs', '5',
             '-outfmt', '6 qseqid sseqid length mismatch qstart qend',
             '-out', out_tsv], logger=logger)
    res = {qid: {'min_tier': None, 'hits': []}
           for qid in seqs_by_id}
    qlen = {qid: len(s) for qid, s in seqs_by_id.items()}
    # 按 (query, subject) 聚合 HSP → 查询坐标合并后统一判层
    agg = {}
    with open(out_tsv, encoding='utf-8', errors='replace') as f:
        for line in f:
            p = line.rstrip('\n').split('\t')
            if len(p) < 6:
                continue
            try:
                qid, sseqid = p[0], p[1]
                mism = int(p[3])
                qstart, qend = int(p[4]), int(p[5])
            except ValueError:
                continue
            if qstart > qend:               # 负链命中：坐标翻转
                qstart, qend = qend, qstart
            agg.setdefault((qid, sseqid), []).append((qstart, qend, mism))
    for (qid, sseqid), hsps in agg.items():
        merged = []
        for qs, qe, m in sorted(hsps):
            if merged and qs <= merged[-1][1] + 1:
                s0, e0, m0 = merged[-1]
                merged[-1] = (s0, max(e0, qe), m0 + m)
            else:
                merged.append((qs, qe, m))
        cur = res[qid]
        need = min(qlen.get(qid, 20), SCREEN_K) - 1
        for s0, e0, m0 in merged:
            if (e0 - s0 + 1) < need:
                continue
            tier = 0 if m0 == 0 else (1 if m0 <= 1 else None)
            if tier is None:
                continue
            if cur['min_tier'] is None or tier < cur['min_tier']:
                cur['min_tier'] = tier
            cur['hits'].append((sseqid, m0))
    return res


# ---------------------------------------------------------------
# 复合评分（CASPER 式效用加权；初始权重经归一）
# ---------------------------------------------------------------
FEATURE_WEIGHTS = {
    'crrna_offtarget_imp':  0.08,
    'fp_offtarget_imp':     0.07,
    'rp_offtarget_imp':     0.07,
    'fp_3p_self_run':       0.045,
    'rp_3p_self_run':       0.045,
    'fp_rp_3p_cross_run':   0.050,
    'fp_self_dimer_pct':    0.10,
    'rp_self_dimer_pct':    0.10,
    'fp_rp_cross_dimer_pct': 0.10,
    'fp_gc_pct':            0.025,
    'rp_gc_pct':            0.025,
    'amplicon_gc_pct':      0.025,
    'seed_max_run':         0.035,
    'seed_gc_pct':          0.025,
    'guide_mfe_kcal':       0.015,
    'guide_seed_unpaired_frac': 0.03,
    'delta_tm_C':           0.08,
    'fp_tm_C':              0.03,
    'rp_tm_C':              0.03,
}
_W_SUM = float(sum(FEATURE_WEIGHTS.values()))

# 特征 → 评估类别（CASPER CATEGORY_FEATURES 口径；报告/前端按类展示）。
# 原「交叉干扰」类（overlap_protospacer/pam 恒 0 死特征）已移除（评审 C4）
CATEGORY_GROUPS = {
    '脱靶': ('crrna_offtarget_imp', 'fp_offtarget_imp', 'rp_offtarget_imp'),
    '二聚体': ('fp_self_dimer_pct', 'rp_self_dimer_pct', 'fp_rp_cross_dimer_pct'),
    '3′端稳定': ('fp_3p_self_run', 'rp_3p_self_run', 'fp_rp_3p_cross_run'),
    'Tm 平衡': ('delta_tm_C', 'fp_tm_C', 'rp_tm_C'),
    'GC 含量': ('fp_gc_pct', 'rp_gc_pct', 'amplicon_gc_pct'),
    '种子与结构': ('seed_max_run', 'seed_gc_pct', 'guide_mfe_kcal',
                    'guide_seed_unpaired_frac'),
}


def category_scores(utilities):
    """按类别聚合效用（类内加权平均，0–100），供雷达/堆叠条形可视化。"""
    out = {}
    for cat, feats in CATEGORY_GROUPS.items():
        wsum = usum = 0.0
        for f in feats:
            u = utilities.get(f)
            if u is None:
                continue
            w = FEATURE_WEIGHTS[f]
            wsum += w
            usum += w * u
        if wsum > 0:
            out[cat] = round(usum / wsum * 100.0, 1)
    return out


def downsample_profile(values, cap=2000):
    """保守性剖面降采样（分箱均值，保端点），控制 report.json 体积。"""
    n = len(values)
    if n <= cap:
        return [round(float(v), 3) for v in values]
    out = []
    for i in range(cap):
        a = int(i * n / cap)
        b = max(a + 1, int((i + 1) * n / cap))
        bin_ = values[a:b]
        out.append(round(sum(bin_) / len(bin_), 3))
    return out


def _u_smaller_lin(x, lo, hi):
    if x is None:
        return 0.5
    return max(0.0, min(1.0, 1.0 - (x - lo) / max(hi - lo, 1e-9)))


def _u_bigger_lin(x, lo, hi):
    if x is None:
        return 0.5
    return max(0.0, min(1.0, (x - lo) / max(hi - lo, 1e-9)))


def _u_smaller_exp(x, half_life=1.0):
    if x is None:
        return 1.0
    tau = max(half_life, 1e-9) / math.log(2.0)
    return max(0.0, min(1.0, math.exp(-(x) / tau)))


def _u_laplace(x, t, s):
    if x is None:
        return 0.5
    return math.exp(-abs(x - t) / max(s, 1e-9))


def _u_plateau(x, a, b, falloff):
    if x is None:
        return 0.5
    if a <= x <= b:
        return 1.0
    d = (a - x) if x < a else (x - b)
    return max(0.0, min(1.0, 1.0 - d / max(falloff, 1e-9)))


def compute_features(crrna, pair, consensus, off_imp, amplicon_gc):
    """对一个 (crRNA, 引物对) 组合算全部特征值（供评分与明细表）。

    off_imp = {'crrna': float, 'fp': float, 'rp': float}——三方各自的
    错配加权脱靶重要度（mm1=2/库，mm2=1/库；由筛查结果预先归并）。
    """
    seed = crrna['seed']
    feats = {
        'crrna_offtarget_imp': off_imp.get('crrna', 0.0),
        'fp_offtarget_imp': off_imp.get('fp', 0.0),
        'rp_offtarget_imp': off_imp.get('rp', 0.0),
        'fp_3p_self_run': three_prime_self_run(pair['f_seq']),
        'rp_3p_self_run': three_prime_self_run(pair['r_seq']),
        'fp_rp_3p_cross_run': cross_3p_run(pair['f_seq'], pair['r_seq']),
        'fp_self_dimer_pct': complementarity_pct(pair['f_seq'],
                                                 pair['f_seq']),
        'rp_self_dimer_pct': complementarity_pct(pair['r_seq'],
                                                 pair['r_seq']),
        'fp_rp_cross_dimer_pct': complementarity_pct(pair['f_seq'],
                                                     pair['r_seq']),
        'fp_gc_pct': pair.get('f_gc'),
        'rp_gc_pct': pair.get('r_gc'),
        'amplicon_gc_pct': amplicon_gc,
        'seed_max_run': max_homo_run(seed),
        'seed_gc_pct': round(gc_frac(seed) * 100.0, 1),
        'guide_mfe_kcal': crrna.get('mfe'),
        'guide_seed_unpaired_frac': crrna.get('seed_unpaired_frac'),
        'delta_tm_C': round(abs((pair.get('f_tm') or 0)
                                - (pair.get('r_tm') or 0)), 1),
        'fp_tm_C': pair.get('f_tm'),
        'rp_tm_C': pair.get('r_tm'),
    }
    return feats


def score_features(feats):
    """特征 → 加权效用总分（0–100）+ 分类效用明细。"""
    u = {
        'crrna_offtarget_imp': _u_smaller_exp(feats['crrna_offtarget_imp']),
        'fp_offtarget_imp': _u_smaller_exp(feats['fp_offtarget_imp']),
        'rp_offtarget_imp': _u_smaller_exp(feats['rp_offtarget_imp']),
        'fp_3p_self_run': _u_smaller_lin(feats['fp_3p_self_run'], 0, 5),
        'rp_3p_self_run': _u_smaller_lin(feats['rp_3p_self_run'], 0, 5),
        'fp_rp_3p_cross_run': _u_smaller_lin(feats['fp_rp_3p_cross_run'], 0, 8),
        'fp_self_dimer_pct': _u_smaller_lin(feats['fp_self_dimer_pct'], 0, 100),
        'rp_self_dimer_pct': _u_smaller_lin(feats['rp_self_dimer_pct'], 0, 100),
        'fp_rp_cross_dimer_pct': _u_smaller_lin(
            feats['fp_rp_cross_dimer_pct'], 0, 100),
        'fp_gc_pct': _u_plateau(feats['fp_gc_pct'], 40, 60, 10),
        'rp_gc_pct': _u_plateau(feats['rp_gc_pct'], 40, 60, 10),
        'amplicon_gc_pct': _u_plateau(feats['amplicon_gc_pct'], 35, 65, 10),
        'seed_max_run': _u_smaller_lin(feats['seed_max_run'], 0, 5),
        'seed_gc_pct': _u_plateau(feats['seed_gc_pct'], 35, 60, 10),
        'guide_mfe_kcal': _u_bigger_lin(feats['guide_mfe_kcal'], -60, 0),
        'guide_seed_unpaired_frac': _u_bigger_lin(
            feats['guide_seed_unpaired_frac'], 0, 1),
        'delta_tm_C': _u_laplace(feats['delta_tm_C'], 0.0, 1.5),
        'fp_tm_C': _u_plateau(feats['fp_tm_C'], 55, 65, 7),
        'rp_tm_C': _u_plateau(feats['rp_tm_C'], 55, 65, 7),
    }
    total = sum(FEATURE_WEIGHTS[k] * u[k] for k in FEATURE_WEIGHTS) / _W_SUM
    return round(total * 100.0, 1), {k: round(v, 3) for k, v in u.items()}


# ---------------------------------------------------------------
# 群体覆盖率回填
# ---------------------------------------------------------------
def isolate_index(seq, k=SCREEN_K):
    """单条分离株的正向 k-mer 集合（小序列，Python set 足够）。"""
    from .dsrna_offtarget import encode_kmers
    return set(int(x) for x in encode_kmers(seq, k))


def check_coverage(spacer, isolate_kmer_sets, k=SCREEN_K):
    """spacer 在各分离株中的存在性（mm0 精确 / mm1 单错配）。

    口径注（评审 C5）：k-mer 编码跳过含非 ACGT（N/兼并码）的窗口，
    含 N 分离株的覆盖率系统性偏保守（安全侧），报告解读时注意。
    """
    from .dsrna_offtarget import encode_kmers
    kmers = [int(x) for x in encode_kmers(spacer, k)]
    cov0 = cov1 = 0
    per_isolate = []
    for name, kset in isolate_kmer_sets:
        if all(x in kset for x in kmers):
            cov0 += 1
            per_isolate.append((name, 0))
        elif any(x in kset for x in kmers):
            cov1 += 1
            per_isolate.append((name, 1))
        else:
            per_isolate.append((name, None))
    n = len(isolate_kmer_sets)
    return {'n_isolates': n, 'mm0': cov0, 'mm0_or_1': cov0 + cov1,
            'pct_mm0': round(100.0 * cov0 / n, 1) if n else 0.0,
            'pct_mm0_or_1': round(100.0 * (cov0 + cov1) / n, 1) if n else 0.0,
            'uncovered': [nm for nm, t in per_isolate if t is None]}


# ---------------------------------------------------------------
# 输出
# ---------------------------------------------------------------
EVIDENCE_NOTE = ('证据边界（ACS Sensors 2026 综述口径）：本结果为序列级候选'
                 '优先级排序，不构成 crRNA 活性、错配区分度、扩增兼容性或'
                 '基质性能的证据；实验结论须按 validation_checklist.md 的'
                 ' L1–L5 逐级建立。')

_CHECKLIST_TMPL = """# CRISPR 检测试剂验证清单（L1–L5）

目标：{target_label}
设计时间：{ts}
组合数：{n_sets}（详见 ranked_sets.tsv）

> {evidence}

## L1 纯化靶标
- [ ] 合成含 protospacer 的 dsDNA（或体外转录 RNA），梯度稀释测 aLoD
- [ ] time-to-positive：记录不同拷贝数下荧光/侧流出峰时间
- [ ] crRNA 活性初筛：与 1–2 个备选 crRNA 并排比较
- [ ] 容错 PAM（如有）：确认非 TTTV PAM 位点的实际活性

## L2 加标基质
- [ ] 健康植株粗提液（纸片提取 / 碱裂解两种前处理各做一版）加标梯度
- [ ] 测回收率与抑制情况（与纯化靶标 aLoD 对比）
- [ ] **近缘交叉反应**：近缘病毒核酸加标，验证不检出
- [ ] 宿主核酸背景：健康植株提取液不检出（对应设计中的宿主筛查）

## L3 回顾性田间样本
- [ ] 用平台历史阳性样品（02_kvsuite 定量结果召回清单）盲测
- [ ] 以 RT-qPCR 为对照，算一致性（敏感度 / 特异度 / kappa）
- [ ] 低滴度样本单独统计（Ct>30 段）

## L4 前瞻性田间试验
- [ ] 多地点 / 多品种 / 不同生育期
- [ ] 操作者为非开发人员；记录单样 <20 min 出结果的比例
- [ ] 现场温度（非恒温实验室条件）下复测

## L5 实际使用
- [ ] 田间站 / 检测点实际人员操作
- [ ] 储存稳定性：37℃ 加速老化 + 常温运输模拟
- [ ] 批次一致性（≥3 个独立合成批次）
- [ ] 单样有效结果成本核算
- [ ] 结果的处置决策记录（诊断结论如何指导防控）

---
引用：Wang H. et al. ACS Sensors 2026（框架）；Low S.J. et al. Commun Biol
2025 PathoGD（群体覆盖口径）；CASPER（初始评分权重）；AutoCORDSv2 /
RPA-Primer-Design / CRISPRPrimerDesigner（规则与参数参考，净室实现）。
"""


def _fmt_order_formats(cand, nuclease):
    """crRNA 订购格式：spacer / DR-spacer / T7 模板（DR 取自酶档案）。"""
    sp = cand['spacer']
    dr_sp = nuclease.dr_dna + sp
    t7 = T7_PROMOTER + nuclease.dr_dna + sp
    t7_risk = cand.get('t7_unsafe') or 'TTTT' in dr_sp
    if nuclease.t7_polyt_hard:
        note = ('Cas13a 路线：T7 体外转录（spacer 已过 poly-T 硬检查；'
                '转录本 5′ 带启动子来源的 GGG 前导）')
    else:
        note = ('⚠ spacer 或 DR-spacer 含 ≥4 连续 T，体外转录可能提前'
                '终止，建议化学合成 crRNA' if t7_risk else
                'T7 体外转录可用（转录本 5′ 带启动子来源的 GGG 前导）')
    return {
        'crRNA_spacer_DNA': sp,
        'crRNA_DR_spacer_DNA': dr_sp,
        'T7_template_DNA': t7,
        'T7_note': note,
    }


def write_outputs(out_dir, result, params, logger=None):
    """落盘全部产物；返回文件清单 dict。"""
    os.makedirs(out_dir, exist_ok=True)
    files = {}

    # 1. crrna_candidates.tsv
    p = os.path.join(out_dir, 'crrna_candidates.tsv')
    with safe_open(p, 'wt') as f:
        f.write('\t'.join([
            'segment', 'rank', 'spacer', 'PAM', 'strand', 'cons_start',
            'cons_end', 'cons_score', 'GC%', 'MFE', 'PAM_priority',
            't7_unsafe', 'host_worst_tier', 'near_worst_tier',
            'veto', 'note']) + '\n')
        for seg, seg_data in result['segments'].items():
            for i, c in enumerate(seg_data['candidates_all'], 1):
                pl = c.get('screening', {}).get('per_lib', {})

                def _worst(prefix):
                    ts = [pl[k].get('min_tier') for k in pl
                          if k.startswith(prefix)
                          and pl[k].get('min_tier') is not None]
                    return min(ts) if ts else ''
                hs = _worst('__host')
                ns = _worst('__near')
                f.write('\t'.join(str(x) for x in [
                    seg, i, c['spacer'], c['pam'], c['strand'],
                    c['site_start'] + 1, c['site_end'] + 1,
                    c['cons_score'], c.get('gc', ''), c.get('mfe', ''),
                    c['pam_priority'],
                    '是' if c.get('t7_unsafe') else '否',
                    '' if hs is None else hs,
                    '' if ns is None else ns,
                    '是' if c.get('vetoed') else '否',
                    c.get('note', '')]) + '\n')
    files['candidates'] = p

    # 2. ranked_sets.tsv
    p = os.path.join(out_dir, 'ranked_sets.tsv')
    with safe_open(p, 'wt') as f:
        f.write('\t'.join([
            'rank', 'score', 'segment', 'crRNA_spacer', 'PAM', 'strand',
            'site', 'cons_score', 'coverage_mm0%', 'forward_primer',
            'f_start', 'f_len', 'f_Tm', 'f_GC%', 'reverse_primer', 'r_start',
            'r_len', 'r_Tm', 'r_GC%', 'product_bp', 'route', 'probe',
            'crrna_offtarget', 'fp_offtarget', 'rp_offtarget',
            'fp_self_dimer%', 'rp_self_dimer%', 'cross_dimer%',
            'guide_MFE', 'flags']) + '\n')
        for i, s in enumerate(result['sets'], 1):
            f.write('\t'.join(str(x) for x in [
                i, s['score'], s['segment'], s['crrna']['spacer'],
                s['crrna']['pam'], s['crrna']['strand'],
                f"{s['crrna']['site_start'] + 1}-{s['crrna']['site_end'] + 1}",
                s['crrna']['cons_score'], s.get('coverage', {}).get('pct_mm0', ''),
                s['pair']['f_seq'], s['pair']['f_start'], s['pair']['f_len'],
                s['pair']['f_tm'], s['pair']['f_gc'],
                s['pair']['r_seq'], s['pair']['r_start'], s['pair']['r_len'],
                s['pair']['r_tm'], s['pair']['r_gc'], s['pair']['product'],
                s['pair']['route'],
                (s.get('probe') or {}).get('marked')
                or (s.get('probe') or {}).get('skipped') or '-',
                s['feats']['crrna_offtarget_imp'],
                s['feats']['fp_offtarget_imp'],
                s['feats']['rp_offtarget_imp'],
                s['feats']['fp_self_dimer_pct'],
                s['feats']['rp_self_dimer_pct'],
                s['feats']['fp_rp_cross_dimer_pct'],
                s['crrna'].get('mfe', ''),
                '; '.join(s['flags']) or '-']) + '\n')
    files['ranked'] = p

    # 2b. panel.tsv（每片段最优组合，多联检视图）
    p = os.path.join(out_dir, 'panel.tsv')
    with safe_open(p, 'wt') as f:
        f.write(chr(9).join([
            'segment', 'score', 'crRNA_spacer', 'PAM', 'site',
            'coverage_mm0%', 'forward_primer', 'reverse_primer',
            'product_bp', 'probe', 'IMM']) + chr(10))
        for s in result.get('panel', []):
            pr = s.get('probe') or {}
            imm = s['crrna'].get('imm')
            f.write(chr(9).join(str(x) for x in [
                s['segment'], s['score'], s['crrna']['spacer'],
                s['crrna']['pam'] or '-',
                str(s['crrna']['site_start'] + 1) + '-'
                + str(s['crrna']['site_end'] + 1),
                s.get('coverage', {}).get('pct_mm0', ''),
                s['pair']['f_seq'], s['pair']['r_seq'],
                s['pair']['product'],
                pr.get('marked') or '-',
                imm['text'] if imm else '-']) + chr(10))
    files['panel'] = p

    # 3. order_sheet.xlsx（openpyxl 可用）/ .tsv 兜底
    try:
        import openpyxl
        p = os.path.join(out_dir, 'order_sheet.xlsx')
        wb = openpyxl.Workbook()
        ws = wb.active
        ws.title = 'crRNA'
        ws.append(['rank', 'segment', '用途', '序列（5′→3′，DNA）',
                   '长度', '备注'])
        for i, s in enumerate(result['sets'], 1):
            fmts = s['order_formats']
            ws.append([i, s['segment'], 'crRNA spacer（合成参照）',
                       fmts['crRNA_spacer_DNA'],
                       len(fmts['crRNA_spacer_DNA']), ''])
            ws.append([i, s['segment'], 'crRNA DR+spacer（合成参照）',
                       fmts['crRNA_DR_spacer_DNA'],
                       len(fmts['crRNA_DR_spacer_DNA']),
                       'LbCas12a handle: ' + DR_RNA])
            ws.append([i, s['segment'], 'T7 体外转录模板',
                       fmts['T7_template_DNA'],
                       len(fmts['T7_template_DNA']), fmts['T7_note']])
            pr = s.get('probe') or {}
            if pr.get('marked'):
                ws.append([i, s['segment'], 'RPA 探针',
                           pr['marked'], pr['len'], pr.get('note', '')])
        ws2 = wb.create_sheet('primers')
        ws2.append(['rank', 'segment', '名称', '序列（5′→3′）', '长度',
                    'Tm', 'GC%', '产物bp'])
        for i, s in enumerate(result['sets'], 1):
            ws2.append([i, s['segment'], 'forward',
                        s['pair']['f_seq'], s['pair']['f_len'],
                        s['pair']['f_tm'], s['pair']['f_gc'],
                        s['pair']['product']])
            ws2.append([i, s['segment'], 'reverse',
                        s['pair']['r_seq'], s['pair']['r_len'],
                        s['pair']['r_tm'], s['pair']['r_gc'], ''])
            if s['pair'].get('f_t7'):
                ws2.append([i, s['segment'],
                            'forward+T7（Cas13a DNA 体系：扩增后体外转录'
                            '生成靶标 RNA）',
                            s['pair']['f_t7'], len(s['pair']['f_t7']), '', ''])
        ws3 = wb.create_sheet('panel')
        ws3.append(['segment', 'score', 'crRNA(DR+spacer)', 'forward',
                    'reverse', 'product', 'probe'])
        for s in result.get('panel', []):
            fmts = s['order_formats']
            pr = s.get('probe') or {}
            ws3.append([s['segment'], s['score'],
                        fmts['crRNA_DR_spacer_DNA'], s['pair']['f_seq'],
                        s['pair']['r_seq'], s['pair']['product'],
                        pr.get('marked', '-')])
        wb.save(p)
        files['order'] = p
    except ImportError:
        p = os.path.join(out_dir, 'order_sheet.tsv')
        with safe_open(p, 'wt') as f:
            f.write('rank\tsegment\t用途\t序列\t长度\t备注\n')
            for i, s in enumerate(result['sets'], 1):
                fmts = s['order_formats']
                f.write(f"{i}\t{s['segment']}\tcrRNA DR+spacer\t"
                        f"{fmts['crRNA_DR_spacer_DNA']}\t"
                        f"{len(fmts['crRNA_DR_spacer_DNA'])}\t\n")
                f.write(f"{i}\t{s['segment']}\tT7模板\t"
                        f"{fmts['T7_template_DNA']}\t"
                        f"{len(fmts['T7_template_DNA'])}\t"
                        f"{fmts['T7_note']}\n")
                f.write(f"{i}\t{s['segment']}\tforward\t{s['pair']['f_seq']}"
                        f"\t{s['pair']['f_len']}\t\n")
                f.write(f"{i}\t{s['segment']}\treverse\t{s['pair']['r_seq']}"
                        f"\t{s['pair']['r_len']}\t\n")
        files['order'] = p

    # 4. validation_checklist.md
    p = os.path.join(out_dir, 'validation_checklist.md')
    with safe_open(p, 'wt') as f:
        f.write(_CHECKLIST_TMPL.format(
            target_label=result['target_label'],
            ts=time.strftime('%Y-%m-%d %H:%M:%S'),
            n_sets=len(result['sets']), evidence=EVIDENCE_NOTE))
    files['checklist'] = p

    # 5. report.json（前端渲染） + crispr_design.json（参数快照）
    # 位点轨道可视化数据：漏斗存活位点 + top 组合标记（rank 0 = 未入选）
    rank_of = {}
    for i, s in enumerate(result['sets'], 1):
        key = (s['segment'], s['crrna']['site_start'], s['crrna']['strand'])
        rank_of[key] = min(rank_of.get(key, 999), i)
    seg_viz = {}
    for seg, d in result['segments'].items():
        sites = []
        for c in d.get('candidates_all', []):
            sites.append({
                's': c['site_start'] + 1, 'e': c['site_end'] + 1,
                'strand': c['strand'], 'pam': c['pam'],
                'cons': c['cons_score'], 'spacer': c['spacer'],
                'veto': bool(c.get('vetoed')),
                'top': rank_of.get((seg, c['site_start'], c['strand']), 0),
            })
        seg_viz[seg] = {
            'n_isolates': d['n_isolates'],
            'consensus_len': d['consensus_len'],
            'n_pam_sites': d['funnel']['total'],
            'n_after_funnel': d['funnel']['kept'],
            'n_sites_designed': d['n_designed'],
            'conservation_profile': d.get('cons_profile'),
            'sites': sites,
        }
    report = {
        'summary': result['summary'],
        'segments': seg_viz,
        'sets': [{
            'rank': i, 'score': s['score'], 'segment': s['segment'],
            'spacer': s['crrna']['spacer'], 'pam': s['crrna']['pam'],
            'strand': s['crrna']['strand'],
            'site': f"{s['crrna']['site_start'] + 1}"
                    f"-{s['crrna']['site_end'] + 1}",
            'cons_score': s['crrna']['cons_score'],
            'mfe': s['crrna'].get('mfe'),
            'coverage_pct': s.get('coverage', {}).get('pct_mm0'),
            'coverage_pct_1mm': s.get('coverage', {}).get('pct_mm0_or_1'),
            'f_seq': s['pair']['f_seq'], 'f_start': s['pair']['f_start'],
            'f_tm': s['pair']['f_tm'], 'f_gc': s['pair']['f_gc'],
            'r_seq': s['pair']['r_seq'], 'r_start': s['pair']['r_start'],
            'r_tm': s['pair']['r_tm'], 'r_gc': s['pair']['r_gc'],
            'product': s['pair']['product'], 'route': s['pair']['route'],
            'probe': s.get('probe'),
            'imm': s['crrna'].get('imm'),
            'f_t7': s['pair'].get('f_t7'),
            'order_formats': s['order_formats'],
            'feats': s['feats'],
            'utilities': s['utilities'],
            'category_scores': category_scores(s['utilities']),
            'flags': s['flags'],
        } for i, s in enumerate(result['sets'], 1)],
        'evidence_note': EVIDENCE_NOTE,
        'exclusions': result['exclusions_meta'],
        'funnel_total': result['funnel_total'],
        'panel': [{'segment': s['segment'], 'score': s['score'],
                   'spacer': s['crrna']['spacer'], 'pam': s['crrna']['pam'],
                   'coverage_pct': s.get('coverage', {}).get('pct_mm0'),
                   'f_seq': s['pair']['f_seq'], 'r_seq': s['pair']['r_seq'],
                   'product': s['pair']['product'],
                   'probe': (s.get('probe') or {}).get('marked') or '-'}
                  for s in result.get('panel', [])],
        'category_weights': {
            cat: round(sum(FEATURE_WEIGHTS[f] for f in feats) / _W_SUM, 4)
            for cat, feats in CATEGORY_GROUPS.items()},
    }
    p = os.path.join(out_dir, 'report.json')
    with safe_open(p, 'wt') as f:
        json.dump(report, f, ensure_ascii=False, indent=1)
    files['report'] = p

    p = os.path.join(out_dir, 'crispr_design.json')
    with safe_open(p, 'wt') as f:
        json.dump({'engine_version': ENGINE_VERSION,
                   'ts': time.strftime('%Y-%m-%d %H:%M:%S'),
                   'params': params,
                   'summary': result['summary'],
                   'files': {k: os.path.basename(v) for k, v in files.items()},
                   'evidence_note': EVIDENCE_NOTE,
                   'citations': [
                       'Wang H. et al. ACS Sensors 2026. '
                       'DOI:10.1021/acssensors.6c02636',
                       'Low S.J. et al. Commun Biol 2025 (PathoGD). '
                       'DOI:10.1038/s42003-025-07591-1',
                       'CASPER (iGEM 2025 Lambert GA) — 初始评分权重',
                       'AutoCORDSv2 / RPA-Primer-Design / '
                       'CRISPRPrimerDesigner — 规则与参数参考（净室实现）',
                   ]}, f, ensure_ascii=False, indent=1)
    files['snapshot'] = p

    # 6. report.html（plotly 交互三图 + top 表 + 证据边界）
    p = os.path.join(out_dir, 'report.html')
    data_json = (json.dumps(report, ensure_ascii=False)
                 .replace('</', '<\\/'))      # 防 </script> 提前闭合
    with safe_open(p, 'wt') as f:
        f.write(_REPORT_HTML_TMPL.replace('__DATA__', data_json))
    files['html'] = p

    if logger:
        logger.log('产物已落盘: ' + ', '.join(
            os.path.basename(v) for v in files.values()))
    return files


_REPORT_HTML_TMPL = """<!DOCTYPE html>
<html lang="zh"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>CRISPR 检测试剂设计报告</title>
<script src="/static/plotly.min.js"
        onerror="window.__NO_PLOTLY=1"></script>
<style>
body{font-family:'Microsoft YaHei',sans-serif;margin:20px;color:#222;
max-width:1180px;background:#fff}
h1{font-size:20px} h2{font-size:15px;margin:22px 0 6px}
.chips span{display:inline-block;background:#f0f4f8;border-radius:12px;
padding:2px 12px;margin-right:8px;font-size:13px}
.note{background:#fff8e6;border-left:4px solid #e6a23c;padding:8px 12px;
font-size:13px}
table{border-collapse:collapse;font-size:12px;margin-top:8px;width:100%}
th,td{border:1px solid #bbb;padding:3px 7px;text-align:left}
th{background:#f0f4f8;position:sticky;top:0}
code{font-family:Consolas,monospace}
.plot{width:100%}
.noplot{background:#fafafa;color:#888;padding:18px;font-size:13px}
.tw{max-height:420px;overflow:auto}
</style></head><body>
<h1>🧬 CRISPR/Cas12a 检测试剂设计报告</h1>
<div class="chips" id="chips"></div>
<div class="note" id="evd"></div>
<h2>① 各片段保守性剖面与 crRNA 位点</h2>
<p class="hint" style="font-size:12px;color:#666">曲线 = 逐列保守性（降采样）；
标记 = 漏斗存活 crRNA 位点（▲/▽ = 正/负链；红 = 入选 top 组合，灰 = 未入选，
空心 = 被排除库否决）。悬停查看 spacer/PAM/保守分/名次。</p>
<div id="segplots"></div>
<h2>② Top 组合评分构成（加权效用堆叠）</h2>
<p class="hint" style="font-size:12px;color:#666">各段长度 = 类别权重 × 类别
效用；条形总长 = 覆盖率折减前的基础分。悬停查看各类明细。</p>
<div class="plot" id="stack"></div>
<h2>③ 评分 vs 群体覆盖率</h2>
<div class="plot" id="scatter"></div>
<h2>④ Top 组合明细</h2>
<div class="tw"><table id="tbl"></table></div>
<script>
const REP = __DATA__;
const CAT_COLORS = {'脱靶':'#d62728','二聚体':'#ff7f0e','3′端稳定':'#8c564b',
  'Tm 平衡':'#e377c2','GC 含量':'#17becf','种子与结构':'#2ca02c',
  '交叉干扰':'#7f7f7f'};
const fmt = v => (v===null||v===undefined) ? '—' : String(v);
document.getElementById('chips').innerHTML =
  ['目标：' + (REP.summary.target||''), '片段 ' + (REP.summary.n_segments||0),
   '组合 ' + (REP.summary.n_sets_output||0) + '/' + (REP.summary.n_sets_total||0),
   '最优评分 ' + (REP.summary.best_score||0),
   'PAM 位点 ' + (REP.funnel_total || REP.summary.funnel_total || 0),
   (REP.summary.amp_mode||'') + ' · ' + (REP.summary.profile||'') + ' 档',
   '耗时 ' + (REP.summary.elapsed_sec||0) + 's']
  .map(s=>'<span>'+s+'</span>').join('');
document.getElementById('evd').textContent = '⚠ ' + (REP.evidence_note||'');
// ---- 表格 ----
const cols = [['rank','#'],['score','评分'],['segment','片段'],
  ['spacer','crRNA spacer'],['pam','PAM'],['site','位点'],
  ['cons_score','保守'],['coverage_pct','覆盖%'],['f_seq','Forward'],
  ['r_seq','Reverse'],['product','产物bp'],['flags','标记']];
document.getElementById('tbl').innerHTML = '<tr>'+cols.map(c=>
  '<th>'+c[1]+'</th>').join('')+'</tr>' + REP.sets.slice(0,20).map(s=>
  '<tr>'+cols.map(c=>'<td>'+fmt(s[c[0]])+'</td>').join('')+'</tr>').join('');
if (window.__NO_PLOTLY || typeof Plotly === 'undefined') {
  ['segplots','stack','scatter'].forEach(id => {
    document.getElementById(id).innerHTML =
      '<div class="noplot">📊 交互图需要通过平台服务打开本报告'
      + '（/tool_runs/&lt;run&gt;/crispr/report.html）以加载 plotly；'
      + '直接双击打开的本地副本仅显示表格。表格数据不受影响。</div>';
  });
} else {
  const cfg = {responsive:true, displayModeBar:false};
  // ---- ① 每片段：保守性曲线 + 位点 ----
  const segHost = document.getElementById('segplots');
  Object.entries(REP.segments||{}).forEach(([seg, d]) => {
    const box = document.createElement('div');
    box.innerHTML = '<h2 style="font-size:13px;color:#555">' + seg +
      ' · ' + d.n_isolates + ' 株 · 共识 ' + d.consensus_len + 'bp · 漏斗后 ' +
      d.n_after_funnel + ' 位点</h2><div class="plot" id="seg_' + seg +
      '" style="height:220px"></div>';
    segHost.appendChild(box);
    const prof = d.conservation_profile || [];
    const L = d.consensus_len || prof.length;
    const stride = L / Math.max(1, prof.length);
    const x = prof.map((_, i) => Math.round((i + 0.5) * stride));
    const traces = [{x:x, y:prof, mode:'lines', name:'保守性',
      line:{color:'#9ecae1', width:1}, hoverinfo:'skip'}];
    const sites = d.sites || [];
    [['top', s=>s.top>0], ['pass', s=>!s.top && !s.veto],
     ['veto', s=>s.veto]].forEach(([nm, pred]) => {
      const ss = sites.filter(pred);
      if (!ss.length) return;
      traces.push({x:ss.map(s=>(s.s+s.e)/2), y:ss.map(()=>1.08),
        mode:'markers', name:nm==='top'?'入选组合':(nm==='veto'?'已否决':'未入选'),
        marker:{symbol:ss.map(s=>s.strand==='+'?'triangle-up':'triangle-down'),
          size:ss.map(s=>s.top===1?11:7),
          color:ss.map(s=>nm==='top'?'#d62728':(nm==='veto'?'#c7c7c7':'#7f7f7f')),
          line:ss.map(s=>s.veto?{width:1,color:'#999'}:{width:0})},
        text:ss.map(s=>'#'+(s.top||'—')+' '+s.spacer+' '+s.pam+
          (s.strand==='+'?' 正链':' 负链')+' 保守'+s.cons+
          (s.veto?' [已否决]':'')),
        hovertemplate:'%{text}<extra></extra>'});
    });
    Plotly.newPlot('seg_'+seg, traces, {margin:{l:56,r:16,t:6,b:46},
      showlegend:true, legend:{orientation:'h', y:1.28, font:{size:10}},
      xaxis:{title:{text:'共识位置 (bp)', font:{size:11}}},
      yaxis:{title:{text:'保守性', font:{size:11}},
        range:[Math.min(0.7, Math.min(...prof, 1)), 1.15],
        tickvals:[0.8,0.9,1.0]}}, cfg);
  });
  // ---- ② 评分构成堆叠 ----
  const cats = Object.keys(REP.category_weights||{});
  const cw = REP.category_weights||{};
  const y = REP.sets.map(s=>'#'+s.rank+' '+(s.spacer||'').slice(0,14)+
    (s.coverage_pct!==null&&s.coverage_pct!==undefined&&s.coverage_pct<100
      ?'（覆盖'+s.coverage_pct+'%）':''));
  Plotly.newPlot('stack', cats.map(c => ({y:y.slice().reverse(),
    x:REP.sets.slice().reverse().map(s => {
      const cs = (s.category_scores||{})[c];
      return cs===undefined?0:(cw[c]*cs/100.0);
    }),
    name:c, type:'bar', orientation:'h', textposition:'none',
    marker:{color:CAT_COLORS[c]||'#999'},
    text:REP.sets.slice().reverse().map(s =>
      c+' '+fmt((s.category_scores||{})[c])+' (w' +
      Math.round((cw[c]||0)*100)+')'),
    hovertemplate:'%{y}<br>%{text}<extra></extra>'})), {
    barmode:'stack', margin:{l:150,r:34,t:8,b:52},
    xaxis:{title:{text:'加权效用（条形总长 = 基础分/100）',
      font:{size:11}}},
    yaxis:{automargin:true, tickfont:{size:10}},
    legend:{orientation:'h', y:1.12, font:{size:10}}}, cfg);
  // ---- ③ 评分 vs 覆盖率 ----
  const segsAll = [...new Set(REP.sets.map(s=>s.segment))];
  const palette = ['#1f77b4','#d62728','#2ca02c','#9467bd','#8c564b',
    '#e377c2','#17becf','#bcbd22'];
  Plotly.newPlot('scatter', segsAll.map((sg,i) => {
    const ss = REP.sets.filter(s=>s.segment===sg &&
      s.coverage_pct!==null && s.coverage_pct!==undefined);
    return {x:ss.map(s=>s.coverage_pct), y:ss.map(s=>s.score),
      mode:'markers', name:sg, text:ss.map(s=>'#'+s.rank+' '+s.spacer+
        '<br>'+s.pam+' '+s.f_seq.slice(0,16)+'…<br>产物 '+s.product+'bp'),
      hovertemplate:'%{text}<br>覆盖 %{x}% · 评分 %{y}<extra></extra>',
      marker:{size:ss.map(s=>Math.max(7, s.product/8)), opacity:0.8,
        color:palette[i%palette.length]}};
  }), (() => {
    const xs = REP.sets.filter(s => s.coverage_pct !== null &&
      s.coverage_pct !== undefined).map(s => s.coverage_pct);
    let x0 = 0, x1 = 100;
    if (xs.length) {
      const lo = Math.min(...xs), hi = Math.max(...xs);
      const span = Math.max(8, (hi - lo) * 0.15);
      x0 = Math.max(-4, lo - span); x1 = Math.min(104, hi + span);
    }
    return {margin:{l:56,r:20,t:8,b:52},
      xaxis:{title:{text:'群体覆盖率 mm0 (%)', font:{size:11}},
        range:[x0, x1]},
      yaxis:{title:{text:'综合评分（含覆盖率因子）', font:{size:11}}},
      legend:{orientation:'h', y:1.15, font:{size:10}}};
  })(), cfg);
}
</script>
</body></html>
"""


# ---------------------------------------------------------------
# RPA 探针（EXO / NFO）与 intentional mismatch 建议（二期）
# ---------------------------------------------------------------
def design_probe(consensus, cons, pair, profile, probe_type,
                 probe_len=(30, 35), logger=None):
    """在扩增子内部（引物结合区之间）设计 EXO/NFO 探针（TwistAmp 口径）。

    - EXO：5′-FAM — THF(中位) — BHQ1-3′（ExoIII 切 THF 释放荧光）
    - NFO：5′-FAM — THF(中位) — C3-3′（EndoIV 切割，侧流读出，
      需配套一端 5′ 生物素标记的引物）
    约束：纯 ACGT、GC 30–70、同聚物 <4、保守分 ≥ 引物档阈值、
    Tm ≥ 60℃（calc_tm），尽量居中。返回 dict 或 {'skipped': 原因}。
    """
    amp0 = pair['f_start'] - 1
    amp1 = amp0 + pair['product'] - 1
    gap = 5
    inner0 = amp0 + pair['f_len'] + gap
    inner1 = amp1 - pair['r_len'] - gap
    if inner1 - inner0 + 1 < probe_len[0]:
        return {'skipped': f'扩增子 {pair["product"]}bp 内部空间不足'
                           f'（引物区之间仅 {inner1 - inner0 + 1}bp）'}
    min_cons = float(profile['min_cons_primer'])
    cands = []
    center = (amp0 + amp1) / 2.0
    for plen in range(probe_len[0], probe_len[1] + 1):
        for s in range(inner0, inner1 - plen + 2):
            seq = consensus[s:s + plen]
            if set(seq) - _PURE_ACGT:
                continue
            gc = gc_frac(seq)
            if not (0.30 <= gc <= 0.70):
                continue
            if max_homo_run(seq) >= 4 or three_prime_gc(seq):
                continue
            csl = cons[s:s + plen]
            if not csl or sum(csl) / len(csl) < min_cons:
                continue
            thf = plen // 2                       # 1 基中位
            cands.append((abs(s + plen / 2.0 - center), s, plen, seq,
                          round(gc * 100.0, 1), thf))
    if not cands:
        return {'skipped': '内部窗口无纯 ACGT 且达保守/GC 阈值的候选'}
    cands.sort(key=lambda t: t[0])
    _, s0, plen, seq, gc, thf = cands[0]
    tm = None
    try:
        tm = round(float(_p3().bindings.calc_tm(seq)), 1)
    except Exception:
        pass
    h1, h2 = seq[:thf - 1], seq[thf:]
    if probe_type == 'exo':
        marked = ('5′-[FAM]-' + h1 + '-[THF]-' + h2 + '-[BHQ1]-3′')
        note = 'EXO：ExoIII 切 THF 释放荧光（实时荧光读出）'
    else:
        marked = ('5′-[FAM]-' + h1 + '-[THF]-' + h2 + '-[C3]-3′')
        note = ('NFO：EndoIV 切割侧流读出；需配套一端引物 5′ 生物素标记'
                '（订购表 primers 页任选其一加 Biotin）')
    return {'probe_seq': seq, 'start': s0 + 1, 'len': plen, 'tm': tm,
            'gc': gc, 'thf_pos': thf, 'marked': marked, 'note': note}


def load_fasta_records(path, cap_records=200):
    """读近缘库记录供 IMM 定位扫描（库通常 MB 级，直接内存扫）。"""
    out = []
    for h, s in iter_fasta(check_path(path, must_exist=True)):
        out.append(((h.split()[0] if h.split() else 'seq'),
                    re.sub(r'[^A-Za-z]', '', (s or '').upper())))
        if len(out) >= cap_records:
            break
    return out


def locate_near_mismatch(spacer, near_recs, max_h=2):
    """定位 spacer 在近缘库中的最佳比对窗口与天然错配位置。

    返回 {'h', 'name', 'mm': [0基错配位…], 'win'} 或 None（无 ≤max_h 窗口）。
    """
    best = None
    L = len(spacer)
    for name, seq in near_recs:
        if not seq:
            continue
        for strand in (seq, revcomp(seq)):
            for i in range(0, len(strand) - L + 1):
                w = strand[i:i + L]
                mm = [j for j in range(L) if spacer[j] != w[j]]
                if len(mm) <= max_h and (best is None
                                          or len(mm) < best['h']):
                    best = {'h': len(mm), 'name': name, 'mm': mm, 'win': w}
                    if best['h'] == 0:
                        return best
    return best


def suggest_intentional_mismatch(cand, near_recs, nuclease):
    """近缘 ≤1–2 错配候选的 intentional mismatch 建议（综述 §guide eng.）。

    策略：天然差异位 q（近缘已有）之外的 PAM 远端位 p 引入人工错配 →
    靶标 1 个人工错配 vs 近缘 ≥2 错配，差分去稳定。位置与序列依赖，
    **必须实验验证区分度增益与活性损失**（报告与此建议一并注明）。
    """
    loc = locate_near_mismatch(cand['spacer'], near_recs)
    if not loc or loc['h'] == 0 or loc['h'] > 2:
        return None
    L = len(cand['spacer'])
    lo = max(nuclease.seed_len + 2, L // 2)       # 避开 PAM 近端 seed
    distal = [p for p in range(lo, L - 1) if p not in loc['mm']]
    if not distal:
        return None
    p = distal[len(distal) // 2]
    cur = cand['spacer'][p]
    near_b = loc['win'][p]
    alts = [b for b in 'ACGT' if b != cur and b != near_b]
    alt = alts[0] if alts else ('A' if cur != 'A' else 'C')
    return {
        'natural_pos': [q + 1 for q in loc['mm']],
        'near': loc['name'], 'h': loc['h'],
        'pos': p + 1, 'from': cur, 'to': alt,
        'text': (f"近缘 {loc['name']} 与 spacer 仅 {loc['h']} 错配"
                 f"（天然差异位 {','.join(str(q + 1) for q in loc['mm'])}）；"
                 f"建议在 PAM 远端第 {p + 1} 位引入人工错配 "
                 f"{cur}->{alt}，使近缘合计 ≥{loc['h'] + 1} 错配、"
                 f"靶标仅 1 个人工错配（差分去稳定）。序列与位置依赖，"
                 f"需实验验证区分度增益与活性损失"),
    }


# ---------------------------------------------------------------
# 主链路
# ---------------------------------------------------------------
class _LogAdapter:
    """统一日志适配器：TaskLogger 型（.log 方法）与任务回调（可调用）皆可。

    phylo._run_mafft / utils.run_cmd 以 logger.log(msg[, level]) 调用；
    dsrna_offtarget._log 兼容可调用与 .log 两种。TaskLogger 本身**不可
    callable**（CLI 直跑实测 TypeError），必须走 .log 分支。
    """

    def __init__(self, fn):
        self.fn = fn

    def __call__(self, msg, level='INFO'):
        fn = self.fn
        if fn is None:
            return
        if hasattr(fn, 'log'):
            try:
                fn.log(msg, level)
            except TypeError:
                try:
                    fn.log(msg)
                except TypeError:
                    pass
            return
        try:
            fn(msg)
        except TypeError:
            try:
                fn(msg, level)
            except TypeError:
                pass

    def log(self, msg, level='INFO'):
        self(msg, level)


def run_chain(targets, params, out_dir, log=None, prog=None, cancel=None):
    """CRISPR 设计全链路入口。

    targets : 目标 FASTA 路径（多分离物，可多分体，header 可含 RNA1 等标签）
    params  : 见 PROFILES + {'amp_mode','profile','spacer_len','flank',
              'site_margin','top_n','use_host','near_fasta',
              'group_by_segment','threads'}
    out_dir : 产物根目录（写 out_dir/crispr/）
    返回 (summary, report_dict)。
    """
    _log = _LogAdapter(log)

    def _prog(stage, frac, msg=''):
        if prog is not None:
            prog(stage, frac, msg)

    t0 = time.time()
    profile_name = params.get('profile') or 'strict'
    profile = dict(PROFILES.get(profile_name) or PROFILES['strict'])
    nuc_key = str(params.get('nuclease') or 'cas12a').strip().lower()
    nuclease = NUCLEASES.get(nuc_key) or NUCLEASES['cas12a']
    # 显式参数覆盖档位
    for key in ('product_min', 'product_max', 'max_sites'):
        v = params.get(key)
        if v not in (None, '', 0, '0'):
            profile[key] = int(v)
    if (nuclease.requires_pam
            and str(params.get('pam_tolerant', '')).lower()
            in ('1', 'true', 'yes')):
        profile['pam_tolerant'] = True
    amp_mode = 'PCR' if str(params.get('amp_mode') or 'RPA').upper() == 'PCR' \
        else 'RPA'
    probe_type = str(params.get('probe_type') or 'none').strip().lower()
    if probe_type not in ('exo', 'nfo'):
        probe_type = 'none'
    if amp_mode != 'RPA' and probe_type != 'none':
        probe_type = 'none'               # TaqMan 探针不在本期口径内
    # spacer 长度按酶档案钳制（Cas12a 20–24 / Cas13a 27–30）：筛查与
    # 覆盖检查的 k-mer 口径为 k=20，<20 无法参与（空窗口假象）
    spacer_len = nuclease.spacer_clamp(params.get('spacer_len'))
    top_n = max(1, min(200, int(params.get('top_n') or 20)))

    out_root = os.path.join(out_dir, 'crispr')
    work_dir = os.path.join(out_root, '_work')
    os.makedirs(work_dir, exist_ok=True)
    # 数据库根一律走 DIRS['databases']：分发版由 config._detect_database_root()
    # 重定向到同级 VirusPlatform-Database，与 PLATFORM_ROOT/databases 不等价。
    # DIRS['databases'] 由 config 保证为非空绝对路径（apply_database_root 的
    # 两个分支都非空），故不需要 `or` 回退。2026-10-07 审查 P1-2。
    cache_dir = os.path.join(DIRS['databases'], 'crispr', 'index_cache')

    # ---------- ① 输入与分组 ----------
    _prog('读取', 0.03, '读取目标 FASTA')
    records = read_targets(targets)
    groups = group_segments(
        records, group_by_segment=str(params.get('group_by_segment', '1'))
        .lower() not in ('0', 'false', 'no'))
    _log('输入 %d 条记录 → %d 个片段组: %s' % (
        len(records), len(groups), ', '.join(
            f'{k}({len(v)})' for k, v in groups.items())))

    # ---------- ② 排除库索引 ----------
    _prog('建索引', 0.08, '构建特异性筛查索引')
    indexes = {}
    exclusions_meta = []
    near_paths = []
    near_in = (params.get('near_fasta') or '').strip()
    if near_in:
        near_paths.append(('near', check_path(near_in, must_exist=True)))
    # databases/crispr/near_neighbors/*.fa 自动纳入
    # 走 DIRS['databases']（与上方 cache_dir 同口径）—— 2026-10-07 审查 P1-2
    nn_dir = os.path.join(DIRS['databases'], 'crispr', 'near_neighbors')
    if os.path.isdir(nn_dir):
        for fn in sorted(os.listdir(nn_dir)):
            if fn.lower().endswith(('.fa', '.fasta', '.fna')):
                near_paths.append(
                    ('near', os.path.join(nn_dir, fn)))

    host_blast_db = None
    host_fa = None
    if str(params.get('use_host', '1')).lower() not in ('0', 'false', 'no'):
        host_fa = current_host_genome()
        if host_fa and os.path.isfile(host_fa):
            mb = _fasta_mb(host_fa)
            host_kmer_ok = False
            if mb <= HOST_KMER_CAP_MB:
                idx = build_exclusion_index('__host__', host_fa, SCREEN_K,
                                            cache_dir, logger=_log,
                                            cancel=cancel)
                if 'arr' in idx:
                    indexes['__host__'] = idx
                    host_kmer_ok = True
                    exclusions_meta.append(
                        {'name': '宿主基因组', 'key': '__host__',
                         'bp': idx['bp'], 'n_records': idx['n_records'],
                         'route': 'k-mer(k=%d)' % SCREEN_K,
                         'cached': idx.get('cached', False)})
                else:
                    # k-mer 失败（典型：内存不足）→ 降级 BLAST 路，绝不
                    # 静默放弃宿主筛查（2026-10-07 真实数据教训）
                    _log('宿主 k-mer 索引失败（%s）→ 降级 BLAST 筛查路'
                         % idx.get('skipped', ''), 'WARN')
            if not host_kmer_ok:
                # 大基因组 / k-mer 失败 → BLAST 路。库建在 %TEMP% 纯 ASCII
                # 目录（2026-10-07 真实数据实测：BLAST+ 的 LMDB 在含中文
                # 的平台根路径下 mdb_env_open 直接失败，与 MAFFT 同坑同
                # 方案 _ascii_work_base 中转）；TEMP 逐用户稳定，可跨运行
                # 缓存，host.fa 按大小比对决定是否重拷
                from .assembly import _ascii_work_base
                import shutil as _sh
                bdir = _ascii_work_base('vp_crispr_blast')
                os.makedirs(bdir, exist_ok=True)
                staged_fa = os.path.join(bdir, 'host.fa')
                if (not os.path.isfile(staged_fa)
                        or os.path.getsize(staged_fa)
                        != os.path.getsize(host_fa)):
                    _sh.copyfile(host_fa, staged_fa)
                host_blast_db = os.path.join(bdir, 'host')
                exclusions_meta.append(
                    {'name': '宿主基因组', 'key': '__host__',
                     'route': 'BLAST(blastn-short)', 'note':
                     f'FASTA {mb:.0f}MB（k-mer 上限 {HOST_KMER_CAP_MB}MB；'
                     'BLAST 库缓存在本机 TEMP，首次运行构建较慢）'})
        else:
            exclusions_meta.append(
                {'name': '宿主基因组', 'key': '__host__',
                 'skipped': '未配置宿主基因组（host-db）；设计未做宿主筛查'})
    else:
        exclusions_meta.append(
            {'name': '宿主基因组', 'key': '__host__', 'skipped': '用户关闭'})

    near_recs = []
    for kind, path in near_paths:         # IMM 定位扫描用（≤200 条）
        try:
            near_recs.extend(load_fasta_records(path))
        except Exception as e:
            _log(f'近缘记录读取失败（IMM 建议将跳过）: {e}', 'WARN')

    for i, (kind, path) in enumerate(near_paths):
        key = '__near__' if len(near_paths) == 1 else f'__near{i}__'
        idx = build_exclusion_index(key, path, SCREEN_K, cache_dir,
                                    logger=_log, cancel=cancel)
        if 'arr' in idx:
            indexes[key] = idx
            exclusions_meta.append(
                {'name': f'近缘病毒({os.path.basename(path)})', 'key': key,
                 'bp': idx['bp'], 'n_records': idx['n_records'],
                 'route': 'k-mer(k=%d)' % SCREEN_K,
                 'cached': idx.get('cached', False)})
        else:
            exclusions_meta.append(
                {'name': f'近缘病毒({os.path.basename(path)})', 'key': key,
                 'skipped': idx['skipped']})

    # ---------- ③ 各片段组：MSA → crRNA → 引物 ----------
    segments_out = {}
    funnel_total = 0
    all_sets = []
    seg_names = sorted(groups)
    for si, seg in enumerate(seg_names):
        recs = groups[seg]
        base_frac = 0.15 + 0.55 * (si / max(1, len(seg_names)))
        _prog('比对', base_frac, f'片段 {seg}：比对与保守性')
        if cancel is not None and cancel.is_set():
            raise RuntimeError('用户取消')
        lens = [len(s) for _, s in recs]
        if max(lens) / max(1, min(lens)) > 1.35 and len(recs) > 1:
            _log(f'片段 {seg} 序列长度差异大（{min(lens)}-{max(lens)}bp），'
                 f'可能混入不同分子，保守性与共识请人工复核', 'WARN')
        consensus, cons, n_used, msa_used = consensus_and_conservation(
            recs, work_dir if len(recs) > 1 else work_dir,
            threads=int(params.get('threads') or 0) or None,
            logger=_log, cancel=cancel)
        _log(f'片段 {seg}: {len(recs)} 条 → 共识 {len(consensus)}bp'
             f'（MSA={"是" if msa_used else "单序列直通"}；'
             f'效应酶 {nuclease.label}）')

        # crRNA 生成 + 漏斗
        _prog('crRNA筛选', base_frac + 0.05, f'片段 {seg}：crRNA 漏斗')
        cands = find_pam_sites(consensus, cons, spacer_len,
                               profile['pam_tolerant'], nuclease)
        funnel_total += len(cands)
        kept, fstats = crrna_funnel(cands, profile, params, nuclease=nuclease,
                                    logger=_log, cancel=cancel)
        if not kept:
            segments_out[seg] = {
                'n_isolates': len(recs), 'consensus_len': len(consensus),
                'funnel': fstats, 'candidates_all': [], 'n_designed': 0,
                'cons_profile': downsample_profile(cons)}
            _log(f'片段 {seg}: 漏斗后无候选（可试放松档）', 'WARN')
            continue

        # 特异性筛查（crRNA 级）
        _prog('特异性', base_frac + 0.10, f'片段 {seg}：crRNA 特异性筛查')
        if host_blast_db is not None:
            # BLAST 路（超大宿主）：批量查 spacer；近缘 k-mer 结果**合并**
            # 进同一 per_lib，不能整包替换（否则近缘筛查被覆盖丢失）
            qmap = {f'c{i}': c['spacer'] for i, c in enumerate(kept)}
            bl = blast_screen_large(qmap, host_blast_db, work_dir,
                                    logger=_log)
            for i, c in enumerate(kept):
                base = (screen_crrna_candidate(c, indexes, SCREEN_K)
                        if indexes else {'veto': False, 'worst_tier': None,
                                         'per_lib': {}})
                r = bl.get(f'c{i}', {'min_tier': None, 'hits': []})
                base['per_lib']['__host__'] = {
                    'min_tier': r['min_tier'], 'n_mm0': 0, 'n_mm1': 0,
                    'n_mm2': 0, 'three_prime_mm1': False,
                    'blast_hits': r['hits'][:5]}
                base['veto'] = base['veto'] or r['min_tier'] == 0
                wt = base.get('worst_tier')
                ht = r['min_tier']
                if ht is not None and (wt is None or ht < wt):
                    base['worst_tier'] = ht
                c['screening'] = base
        else:
            for c in kept:
                if cancel is not None and cancel.is_set():
                    raise RuntimeError('用户取消')
                c['screening'] = (screen_crrna_candidate(c, indexes, SCREEN_K)
                                  if indexes else {'veto': False,
                                                   'worst_tier': None,
                                                   'per_lib': {}})
        n_veto = 0
        for c in kept:
            if c['screening']['veto']:
                c['vetoed'] = True
                c['note'] = '排除库精确命中（mm0），剔除'
                n_veto += 1
        pool = [c for c in kept if not c.get('vetoed')]
        # 排序已含 PAM 优先级；mm1 近缘命中降权（note 标记）
        for c in pool:
            wt = c['screening'].get('worst_tier')
            if wt == 1:
                c['note'] = '近缘库 ≤1 错配命中：区分度风险，实验需验证'
                if near_recs:
                    imm = suggest_intentional_mismatch(
                        c, near_recs, nuclease)
                    if imm:
                        c['imm'] = imm
        segments_out[seg] = {
            'n_isolates': len(recs), 'consensus_len': len(consensus),
            'funnel': fstats, 'candidates_all': kept, 'n_designed': 0,
            'cons_profile': downsample_profile(cons)}

        # 引物设计（top max_sites 个位点）
        _prog('引物', base_frac + 0.18, f'片段 {seg}：引物设计')
        design_pool = pool[:int(profile['max_sites'])]
        iso_kmer_sets = [(nm, isolate_index(sq, SCREEN_K)) for nm, sq in recs]
        n_designed = 0
        for c in design_pool:
            if cancel is not None and cancel.is_set():
                raise RuntimeError('用户取消')
            pairs = design_primers_for_site(
                consensus, cons, c, amp_mode, profile, params, work_dir,
                logger=_log, cancel=cancel)
            for pair in pairs:
                amp0 = pair['f_start'] - 1
                amp1 = amp0 + pair['product'] - 1
                amp_gc = round(gc_frac(consensus[amp0:amp1 + 1]) * 100.0, 1)
                scr_crrna = c['screening']
                fp_scr = (screen_primer_pair(
                    {'f_seq': pair['f_seq'], 'r_seq': pair['r_seq']},
                    indexes, SCREEN_K) if indexes else None)
                if fp_scr is not None and fp_scr['veto']:
                    continue          # 引物 mm0 脱靶 → 整对剔除
                feats = compute_features(
                    c, pair, consensus,
                    {'crrna': _off_imp(scr_crrna.get('per_lib', {})),
                     'fp': (_off_imp(fp_scr['per_primer']['f'])
                            if fp_scr is not None else 0.0),
                     'rp': (_off_imp(fp_scr['per_primer']['r'])
                            if fp_scr is not None else 0.0)},
                    amp_gc)
                score, utils = score_features(feats)
                cov = check_coverage(c['spacer'], iso_kmer_sets, SCREEN_K)
                # 群体覆盖率乘性因子（PathoGD prevalence 口径）：序列特征
                # 再好的候选，覆盖不了目标群体就不是可用试剂。多分离株输入
                # 时 总分 × (0.5 + 0.5×mm0覆盖率%)；单分离株不折减。
                if cov['n_isolates'] > 1:
                    score = round(score * (0.5 + 0.5 * cov['pct_mm0']
                                           / 100.0), 1)
                flags = []
                if c['screening'].get('worst_tier') == 1:
                    flags.append('近缘≤1mm')
                if c.get('t7_unsafe'):
                    flags.append('T7 风险(≥4T)')
                if fp_scr is not None:
                    for tag in ('f', 'r'):
                        for lk, r in fp_scr['per_primer'][tag].items():
                            if r.get('three_prime_mm1'):
                                flags.append(f"{'F' if tag == 'f' else 'R'}"
                                              f"引物3'端≤1mm@{lk}")
                if (cov['pct_mm0'] or 0) < 90 and cov['n_isolates'] > 1:
                    flags.append(f"覆盖{cov['pct_mm0']}%")
                # Cas13a + DNA 扩增体系：前向引物 5′ 加 T7 启动子，
                # 扩增后体外转录生成靶标 RNA 供 Cas13a 识别
                if nuclease.key == 'cas13a':
                    pair['f_t7'] = T7_PROMOTER + pair['f_seq']
                probe = None
                if probe_type != 'none':
                    if amp_mode != 'RPA':
                        probe = {'skipped': 'PCR 模式不设计 RPA 探针'}
                    else:
                        probe = design_probe(consensus, cons, pair,
                                             profile, probe_type,
                                             logger=_log)
                    if probe and 'skipped' not in probe:
                        flags.append('探针' + probe_type.upper())
                all_sets.append({
                    'segment': seg, 'crrna': c, 'pair': pair, 'feats': feats,
                    'utilities': utils, 'score': score, 'flags': flags,
                    'coverage': cov, 'probe': probe,
                    'order_formats': _fmt_order_formats(c, nuclease),
                })
                n_designed += 1
        segments_out[seg]['n_designed'] = n_designed
        _log(f'片段 {seg}: {len(design_pool)} 位点 → {n_designed} 组合'
             f'（crRNA 筛查剔除 {n_veto}）')

    # ---------- ④ 排序与输出 ----------
    _prog('评分', 0.75, '复合评分排序')
    all_sets.sort(key=lambda s: -s['score'])
    # 每个位点最多留 2 组合，避免单位点霸榜
    per_site = {}
    ranked = []
    for s in all_sets:
        k = (s['segment'], s['crrna']['site_start'], s['crrna']['strand'])
        if per_site.get(k, 0) >= 2:
            continue
        per_site[k] = per_site.get(k, 0) + 1
        ranked.append(s)
    ranked = ranked[:top_n]
    # 跨片段 panel：每片段取评分最高的一套组合（多联检视图）
    panel, seen_seg = [], set()
    for s in ranked:
        if s['segment'] in seen_seg:
            continue
        seen_seg.add(s['segment'])
        panel.append(s)

    _prog('产出', 0.92, '写出产物')
    target_label = os.path.basename(str(targets))
    # 摘要必须先于 write_outputs 组装：report.json 的 chips/前端要用
    # （此前 summary 在 write_outputs 之后才填，report.json 里嵌的是 {}）
    summary = {
        'tool': 'crispr',
        'ok': True,
        'engine_version': ENGINE_VERSION,
        'target': target_label,
        'n_input_records': len(records),
        'n_segments': len(groups),
        'segments': {seg: {
            'n_isolates': d['n_isolates'],
            'consensus_len': d['consensus_len'],
            'n_pam_sites': d['funnel']['total'],
            'n_after_funnel': d['funnel']['kept'],
            'n_designed': d['n_designed'],
        } for seg, d in segments_out.items()},
        'exclusions': exclusions_meta,
        'n_sets_total': len(all_sets),
        'n_sets_output': len(ranked),
        'best_score': ranked[0]['score'] if ranked else 0.0,
        'amp_mode': amp_mode,
        'profile': profile_name,
        'nuclease': nuclease.label,
        'nuclease_key': nuclease.key,
        'pfs_note': nuclease.pfs_note,
        'out_dir': out_root,
        'elapsed_sec': round(time.time() - t0, 1),
        'evidence_note': EVIDENCE_NOTE,
    }
    result = {
        'target_label': target_label,
        'segments': segments_out,
        'sets': ranked,
        'panel': panel,
        'exclusions_meta': exclusions_meta,
        'funnel_total': funnel_total,
        'summary': summary,
    }
    files = write_outputs(out_root, result, {
        'profile': profile_name, 'amp_mode': amp_mode,
        'nuclease': nuclease.key, 'pfs_note': nuclease.pfs_note,
        'probe_type': probe_type,
        'spacer_len': spacer_len, 'flank': params.get('flank') or 250,
        'site_margin': params.get('site_margin') or 25,
        'top_n': top_n, 'use_host': params.get('use_host', '1'),
        'near_fasta': near_in, 'group_by_segment':
            params.get('group_by_segment', '1'),
        'screen_k': SCREEN_K,
        'engine': 'VirusPlatform crispr_design ' + ENGINE_VERSION,
    }, logger=_log)

    summary['files'] = {k: os.path.relpath(v, out_dir)
                        for k, v in files.items()}
    summary['elapsed_sec'] = round(time.time() - t0, 1)
    result['summary'] = summary
    with safe_open(os.path.join(out_root, 'summary.json'), 'wt') as f:
        json.dump(summary, f, ensure_ascii=False, indent=1)
    try:
        shutil.rmtree(work_dir, ignore_errors=True)
    except OSError:
        pass
    _prog('完成', 1.0, '完成')
    _log('CRISPR 设计完成：%d 组合输出（耗时 %.1fs）'
         % (len(ranked), summary['elapsed_sec']))
    return summary, result


def _off_imp(per_lib):
    """错配加权脱靶重要度：mm0=4/库（本应 veto，兜底计），mm1=2，mm2=1。"""
    v = 0.0
    for r in per_lib.values():
        t = r.get('min_tier') if isinstance(r, dict) else None
        if t == 0:
            v += 4.0
        elif t == 1:
            v += 2.0
        elif t == 2:
            v += 1.0
    return v


def _merge_primer_feats(feats, fp_scr, indexes):
    """把引物对筛查结果并进特征（fp/rp 分别的脱靶重要度）。"""
    for tag, fk in (('f', 'fp_offtarget_imp'), ('r', 'rp_offtarget_imp')):
        feats[fk] = _off_imp(fp_scr['per_primer'][tag])
