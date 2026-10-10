"""
PAmiRDB v3 — TAPIR 引擎（erbon7/tapir 的 tapir_fasta 忠实移植版）。

移植自 https://github.com/erbon7/tapir （TAPIR v1.2, Bonnet et al. 2010,
Bioinformatics 26:1566 — 植物miRNA 靶标预测，含 target mimics）。

tapir_fasta 口径（与 src/tapir_fasta 逐条对应）：
  1. 比对：fasta34（miRNA 反向互补链 vs 目标，-E 150），query 链须为
     反向（strand == -1，即反平行配对）——本引擎对 miRanda 比对列直接
     判定，天然反平行；
  2. 逐列分类（q = miRNA 反向互补字符，h = 目标字符）：
       q == h                        → 配对（不扣分）
       (q,h) = (A,G) / (C,T)         → G:U 摆动（gu）
       其余                           → 错配（mismatch）
       任一侧 '-'                     → 缺口（gap）
  3. 种子区 = 比对列里 miRNA 第 2-12 位（跨缺口计数，get_upper_lower）；
  4. score = gap + mismatch + seed_mismatch*2 + seed_gap*2 + seed_gu
             + gu*0.5 ；截止 score ≤ 3.5；
  5. mfe ratio：约束折叠（linker "GCGGGGACGC"，约束 "(((xxxx)))"）
       ratio = MFE(mir+linker+target) / MFE(mir+linker+revcomp(mir))
     截止 ratio ≥ 0.7。
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Optional


# ── 参考实现缺省参数 ────────────────────────────────────────────────
SCORE_CUTOFF = 3.5            # tapir_fasta --score 缺省
MFE_RATIO_CUTOFF = 0.7        # tapir_fasta --mfe_ratio 缺省
LINKER = 'GCGGGGACGC'         # 约束折叠的环区 linkers
LINKER_CS = '(((xxxx)))'      # linker 的结构约束
SEED_POS_START = 2            # 种子区（miRNA 原始 5'→3' 帧）起点
SEED_POS_END = 12             # 种子区终点（get_upper_lower 口径 2-12）


@dataclass
class TAPIRResult:
    """TAPIR target prediction result."""
    mirna_id: str
    target_id: str
    seed_score: float            # 种子区罚分（seed_mismatch*2+seed_gap*2+seed_gu）
    accessibility_score: float    # mfe ratio × 100（≥70 通过）
    energy_score: float           # mir+linker+target 约束折叠 MFE
    tapir_score: float            # TAPIR 原生罚分（score，越低越好）
    target_region_structure: str  # 约束折叠结构（dot-bracket）
    up_energy: float              # self 对照 MFE（mir+linker+revcomp(mir)）
    passes: bool
    consensus_votes: int = 0


def _revcomp_rna(seq: str) -> str:
    comp = {'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G', '-': ''}
    return ''.join(comp.get(b, 'N') for b in reversed(seq))


# ── score() 子例程的忠实移植 ────────────────────────────────────────

def _classify_column(q: str, h: str):
    """q = miRNA 反向互补字符，h = 目标字符（tapir_fasta score() 的分类）。

    返回 'gap' | 'match' | 'gu' | 'mismatch'。"""
    if q == '-' or h == '-':
        return 'gap'
    if q == h:
        return 'match'
    if (q, h) in (('A', 'G'), ('C', 'T'), ('A', 'U'), ('C', 'U')):
        # tapir 以 DNA 字母判 (A,G)/(C,T)；我们把 T/U 归一后再判一次
        if (q, h) in (('A', 'G'), ('C', 'T'), ('A', 'U'), ('C', 'U')):
            return 'gu'
    return 'mismatch'


def _norm_rna(ch: str) -> str:
    return ch.upper().replace('T', 'U')


def score_alignment_columns(mline: str, tline: str):
    """对 miRanda 比对的两行按 tapir_fasta::score 打分。

    **列语义（2026-10-03 修复 P1-2）**：mline/tline 必须是**比对列**
    （miRanda 书写约定：上行 miRNA 3'→5'、下行 target 5'→3'，含缺口），
    不能传「无缺口 miRNA vs 原始目标片段」——反平行位点按同向两条原始
    序列逐位比对时列语义整体错位（错位后的列既不是配对也不是缺口）。
    tapir 的 q = miRNA 在该列的互补字符，q == h 即配对。

    种子区 = miRNA 原始 5' 第 2-12 位（get_upper_lower 口径）。上行是
    3'→5' 帧，故先翻回原始 5'→3' 帧（mline[::-1]）再自左向右数非缺口
    字符；原始帧第 k 个非缺口字符对应的列下标 i 换算回上行：L-1-i 的
    翻转位（本函数直接记录翻转后的列下标，等价于原来自右向左数）。

    两行长度不一致（target 短于 miRNA 等）时按缺口补齐到等长
    ——原实现对短 tline 直接索引会 IndexError（2026-10-03 修复）。

    返回 (score, gap, seed_gap, gu, seed_gu, match, mismatch, seed_mismatch,
          seed_start_col, seed_stop_col, aln_str)。"""
    L = max(len(mline), len(tline))
    if not L:
        return (0.0, 0, 0, 0, 0, 0, 0, 0, -1, -1, '')
    mline = mline.ljust(L, '-')
    tline = tline.ljust(L, '-')

    # 种子列定位：翻回 miRNA 原始 5'→3' 帧后自左向右数第 2..12 个非缺口
    # 字符；记录其在上行中的列下标（自右向左数第 2..12 个的等价写法）。
    seed_cols = []
    k = 0
    for i_fwd, ch in enumerate(mline[::-1]):
        if ch != '-':
            k += 1
            if SEED_POS_START <= k <= SEED_POS_END:
                seed_cols.append(L - 1 - i_fwd)
    seed_start_col = min(seed_cols) if seed_cols else -1
    seed_stop_col = max(seed_cols) if seed_cols else -1

    gap = seed_gap = gu = seed_gu = match = mismatch = seed_mismatch = 0
    aln_chars = []
    for i in range(L):
        q = _norm_rna(_revcomp_char(mline[i]))
        h = _norm_rna(tline[i])
        cls = _classify_column(q, h)
        in_seed = (seed_start_col != -1 and seed_stop_col != -1 and
                   seed_start_col <= i <= seed_stop_col)
        if cls == 'gap':
            if in_seed:
                seed_gap += 1
            else:
                gap += 1
            aln_chars.append('.')
        elif cls == 'match':
            match += 1
            aln_chars.append('|')
        elif cls == 'gu':
            if in_seed:
                seed_gu += 1
            else:
                gu += 1
            aln_chars.append('o')
        else:
            if in_seed:
                seed_mismatch += 1
            else:
                mismatch += 1
            aln_chars.append('.')
    score = (gap + mismatch + seed_mismatch * 2 + seed_gap * 2 +
             seed_gu + gu * 0.5)
    return (score, gap, seed_gap, gu, seed_gu, match, mismatch,
            seed_mismatch, seed_start_col, seed_stop_col, ''.join(aln_chars))


def _revcomp_char(m: str) -> str:
    """miRNA 字符的互补字符（比对照配对方向的反向链字符）。"""
    comp = {'A': 'T', 'U': 'A', 'T': 'A', 'G': 'C', 'C': 'G', '-': '-'}
    return comp.get(m.upper(), 'N')


# ── constraint_fold 的忠实移植（ViennaRNA Python 绑定） ─────────────

def _constraint_fold(mir: str, reftarget: str):
    """(mfe, structure)：mir + linker + reftarget 的约束折叠。

    construct = mir + linker + target（T→U），constraint = '.'*（mir/target）
    + '(((xxxx)))'（linker）。对应 tapir_fasta 的 constraint_fold 子例程。"""
    import RNA
    construct = (mir + LINKER + reftarget).upper().replace('T', 'U')
    constraint = ('.' * len(mir)) + LINKER_CS + ('.' * len(reftarget))
    fc = RNA.fold_compound(construct)
    fc.constraints_add(constraint)
    ss, mfe = fc.mfe()
    return mfe, ss


def check_mfe_ratio(alq: str, alh: str) -> float:
    """tapir_fasta::check_mfe：约束折叠 self/daemon 两个 MFE 的比值。

    alq = miRNA 比对串（含缺口），alh = 目标比对串。比值 = mfe2/mfe1，
    两值皆为负，故比值越高（→1）代表 target 区配对越接近 self 对照。"""
    mir = alq.replace('-', '')
    revmir = _revcomp_rna(mir)
    hit = alh.replace('-', '')

    mfe1, _ = _constraint_fold(mir, revmir)
    mfe2, _ = _constraint_fold(mir, hit)
    if mfe1 == 0:
        return 0.0
    return mfe2 / mfe1


def tapir_predict(mirna_id: str, mirna_seq: str,
                  target_id: str, target_seq: str,
                  site_start: int, site_end: int,
                  ) -> Optional[TAPIRResult]:
    """在给定 target 区间上运行 tapir_fasta 口径（评分 + mfe ratio）。

    列语义（2026-10-03 修复）：site 是 target 5'→3' 子串，与 miRNA 反平行
    ——site 第 i 位对应 miRNA 原始位置 L-i+1，故把 miRNA 反转后同向逐位比
    （q = 该位 miRNA 字符的互补字符，q == h 即配对），与 score_alignment_columns
    的列语义一致。原实现直接拿 miRNA 5'→3' 原向与 site 同向比，等于找
    comp(miRNA) 镜像位点（与 P0-4 同类错误），已修正。
    """
    if site_end <= site_start or site_end > len(target_seq):
        return None
    site = target_seq[site_start - 1:site_end]
    mirna = mirna_seq.upper().replace('T', 'U')
    mir_rev = mirna[::-1]        # 反转帧：与 target 5'→3' 反平行
    L = len(mirna)
    if not site or not mirna:
        return None

    # 评分（列对：miRNA 互补字符 vs 目标字符；列 i ↔ miRNA 原始位置 L-i+1）
    score = gap = seed_gap = gu = seed_gu = match = mismatch = seed_mismatch = 0
    aln_chars = []
    for i, (m, h) in enumerate(zip(mir_rev, site), start=1):
        q = _revcomp_char(m)
        h_rna = h.upper().replace('T', 'U')
        in_seed = SEED_POS_START <= (L - i + 1) <= SEED_POS_END
        if q == '-' or h == '-':
            cls = 'gap'
        elif q == h_rna:
            cls = 'match'
        elif (q, h_rna) in (('A', 'G'), ('C', 'U')):
            cls = 'gu'
        else:
            cls = 'mismatch'
        if cls == 'gap':
            if in_seed:
                seed_gap += 1
            else:
                gap += 1
            aln_chars.append('.')
        elif cls == 'match':
            match += 1
            aln_chars.append('|')
        elif cls == 'gu':
            if in_seed:
                seed_gu += 1
            else:
                gu += 1
            aln_chars.append('o')
        else:
            if in_seed:
                seed_mismatch += 1
            else:
                mismatch += 1
            aln_chars.append('.')
    score = (gap + mismatch + seed_mismatch * 2 + seed_gap * 2 +
             seed_gu + gu * 0.5)

    # mfe ratio（constraint_fold）。ViennaRNA 不可用时 mfe1/mfe2 归零，
    # 避免 try 外引用未绑定变量（原实现在 ImportError 时 NameError）。
    mfe1 = mfe2 = 0.0
    try:
        mfe1, _ = _constraint_fold(mirna, revcomp_of(mirna))
        mfe2, _ = _constraint_fold(mirna, site)
        ratio = (mfe2 / mfe1) if mfe1 else 0.0
    except Exception:
        ratio = 0.0

    passes = (score <= SCORE_CUTOFF) and (ratio >= MFE_RATIO_CUTOFF)
    return TAPIRResult(
        mirna_id=mirna_id, target_id=target_id,
        seed_score=float(seed_mismatch * 2 + seed_gap * 2 + seed_gu),
        accessibility_score=round(ratio * 100, 1),
        energy_score=float(mfe2),
        tapir_score=float(score),
        target_region_structure='',
        up_energy=float(mfe1),
        passes=passes,
        consensus_votes=0,
    )


def revcomp_of(seq: str) -> str:
    comp = {'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G'}
    return ''.join(comp.get(b, 'N') for b in reversed(seq))


def validate_with_tapir(miranda_result, mirna_id: str,
                         target_id: str, target_seq: str,
                         ) -> Optional[TAPIRResult]:
    """Validate a miRanda hit with TAPIR scoring.

    对 miRanda 比对列直接套 tapir_fasta 口径（评分 + mfe ratio）。

    2026-10-03 修复 P1-2①：传入**比对列**（miRanda 三行格式的上下两行，
    含缺口，上行 miRNA 3'→5'、下行 target 5'→3'），而不是「无缺口 miRNA
    vs 原始目标片段」——后者列语义整体错位（反平行位点被当同向序列比）。
    兼容旧格式（上行 5'→3'）：此时把上行反转后再按新帧口径打分。
    """
    mirna_seq = miranda_result.mirna_seq
    t_start = miranda_result.target_start
    t_end = miranda_result.target_end
    if t_end <= t_start or t_start < 1:
        return None
    site = target_seq[t_start - 1:t_end]
    mirna = mirna_seq.upper().replace('T', 'U')
    revmir = _revcomp_rna(mirna)
    if not site:
        return None

    lines = (getattr(miranda_result, 'alignment_str', '') or '').split('\n')
    if len(lines) < 3:
        return None
    top_label = lines[0][:12]
    mline = (lines[0].replace("miRNA  3' ", "").replace("miRNA  5' ", "")
             .replace(" 5'", "").replace(" 3'", "").strip())
    tline = (lines[2].replace("target 5' ", "").replace("target 3' ", "")
             .replace(" 5'", "").replace(" 3'", "").strip())
    if not mline or not tline:
        return None
    if "5'" in top_label:                # 旧格式：上行原为 5'→3'（镜像帧）
        mline = mline[::-1]
    # 长度不足时 score_alignment_columns 内部按缺口补齐（不再 IndexError）

    score, gap, seed_gap, gu, seed_gu, match, mismatch, seed_mismatch, \
        _s_col, _e_col, _aln = score_alignment_columns(mline, tline)

    # mfe ratio（constraint_fold：self 对照 vs 目标区）
    mfe1, _st1 = _constraint_fold(mirna, revmir)
    mfe2, _st2 = _constraint_fold(mirna, site)
    ratio = (mfe2 / mfe1) if mfe1 else 0.0

    passes = (score <= SCORE_CUTOFF) and (ratio >= MFE_RATIO_CUTOFF)
    return TAPIRResult(
        mirna_id=mirna_id, target_id=target_id,
        seed_score=float(seed_mismatch * 2 + seed_gap * 2 + seed_gu),
        accessibility_score=round(ratio * 100, 1),
        energy_score=float(mfe2),
        tapir_score=float(score),
        target_region_structure='',
        up_energy=float(mfe1),
        passes=passes,
        consensus_votes=0,
    )


# ── 自检 ────────────────────────────────────────────────────────────

if __name__ == '__main__':
    mirna = "UUUGGAUUGAAGGGAGCUCUA"
    target = "GUGAAACCUAACUUCCCUCGAGAUUUCAACGAUU"
    r = validate_with_tapir(
        type('MR', (), {'mirna_seq': mirna, 'target_start': 6,
                        'target_end': 26})(),
        'ath-miR159a', 'test_target', target)
    print('passes:', r.passes if r else None,
          '| score:', r.tapir_score if r else None,
          '| ratio:', r.accessibility_score if r else None)
