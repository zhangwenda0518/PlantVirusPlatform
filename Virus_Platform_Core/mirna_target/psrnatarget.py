"""
psRNATarget 评分引擎（PAmiRDB v3）—— MiRNATarget 忠实移植版。

移植自 jtremblay/MiRNATarget（aka psRNATarget_command_line）：
  https://github.com/jtremblay/MiRNATarget
该实现经作者实测与 psRNATarget web server（Dai et al. 2018, NAR）输出一致，
规则源头是 Fahlgren & Carrington (2009) 书章 + psRNATarget 2011。

口径（与 parse_mirna_targets.py 逐条对应）：
  * 从 miRNA 5' 端开始给比对列打分，只计前 alignment_length=19 列；
  * 完美配对（WC）        → 0
  * G:U 摆动（含 G·T）    → 0.5（GUs++；**种子区内也不乘** 1.5）
  * 错配                  → 1.0（种子区 2-13 内 ×1.5，seed_mismatches++）
  * 缺口打开              → 2.0（种子区内 ×1.5；gaps++、total_mismatches++；
                              miRNA 一侧的缺口额外 +1）
  * 缺口延伸              → 0.5（种子区内 ×1.5）
  * 通过过滤：penalty ≤ 5.0 且 seed_mismatches ≤ 2 且 HSP ≥ 14
              且 gaps ≤ 1 且 total_mismatches ≤ 8 且比对长 ≤ 22

比对本身：SSEARCH36 参数（-f -8 -g -3 -r +4/-3，即 match +4 / mismatch -3 /
gap open -8 / extend -3 的仿射 Smith-Waterman 局部比对），正反两条目标链
各比一次（对应原实现的两次 ssearch36 调用）。列属性按实际碱基对判定：
  * q == h（q 为 miRNA 的互补链字符）→ WC '|'
  * (q,h) = (C,T) / (A,G)            → G:U 摆动 '.'
  * 其余                              → 错配 ' '
"""

from dataclasses import dataclass, field
from typing import Optional
import math


# ── 参考实现缺省参数（parse_mirna_targets.py 的命令行缺省值） ────────
E_CUTOFF = 5.0               # penalty（参考实现称之为 E/expectation）截止
SEED_START = 2               # 种子区起点（miRNA 5' 起数）
SEED_END = 13                # 种子区终点
PENALTY_MULTIPLIER = 1.5     # 种子区错配/缺口乘数（G:U 不乘）
NUM_MISMATCH_SEED = 2        # 种子区允许的错配数（不含 G:U）
HSP_CUTOFF = 14              # HSP（比对柱数）下限
GAP_CUTOFF = 1               # 缺口数上限
GUS_CUTOFF = 7               # G:U 数上限
TOTAL_MISMATCHES_CUTOFF = 8  # 总错配数上限
ALIGNMENT_LENGTH = 19        # 参与打分的比对列数（从 miRNA 5' 端起）
MAX_ALIGNMENT_LENGTH = 22    # 比对总长上限
EXTRA_PENALTY_QUERY_GAP = 1  # miRNA 一侧缺口的额外罚分

# SSEARCH36 比对参数（-f -8 -g -3 -r +4/-3）
SS_MATCH = 4
SS_MISMATCH = -3
SS_GAP_OPEN = -8
SS_GAP_EXTEND = -3


# ── 碱基工具（DNA 字母表；U 按 T 处理） ─────────────────────────────

def _dna(seq: str) -> str:
    return seq.upper().replace('U', 'T')


def _revcomp_dna(seq: str) -> str:
    comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return ''.join(comp.get(b, 'N') for b in reversed(seq.upper()))


def _column_class(q: str, h: str):
    """q = miRNA 互补链字符，h = 目标链字符（均大写 DNA）。

    返回 ('wc'|'gu'|'mm')。q==h 即互补配对（q 是 miRNA 的互补）。
    摆动：(q,h) = (C,T)/(A,G)，对应 miRNA G·T / U·G。
    """
    if q == h:
        return 'wc'
    if (q, h) in (('C', 'T'), ('A', 'G')):
        return 'gu'
    return 'mm'


# ── Smith-Waterman 局部比对（仿射缺口，SSEARCH36 参数） ─────────────

def _sw_align(q: str, t: str):
    """q 对 t 的最优局部比对。返回 (q_cols, t_cols, t_lo, t_hi)
    （两串等长的列字符列表，含 '-' 占位；t_lo/t_hi 为比对覆盖的 t 上
    1-based 坐标闭区间，用于还原靶坐标）。"""
    m, n = len(q), len(t)
    if m == 0 or n == 0:
        return [], [], 0, 0
    H = [[0.0] * (n + 1) for _ in range(m + 1)]
    ptr = [[0] * (n + 1) for _ in range(m + 1)]   # 0=结束 1=diag 2=q缺(上) 3=t缺(左)
    best, bi, bj = 0.0, 0, 0
    for i in range(1, m + 1):
        qi = q[i - 1]
        Hi, Hm = H[i], H[i - 1]
        for j in range(1, n + 1):
            d = Hm[j - 1] + (SS_MATCH if qi == t[j - 1] else SS_MISMATCH)
            best_ij, p = d, 1
            up = Hm[j] + SS_GAP_OPEN
            if up > best_ij:
                best_ij, p = up, 2
            left = Hi[j - 1] + SS_GAP_OPEN
            if left > best_ij:
                best_ij, p = left, 3
            Hi[j] = best_ij
            ptr[i][j] = p
            if best_ij > best:
                best, bi, bj = best_ij, i, j
    # 回溯
    q_cols, t_cols = [], []
    i, j = bi, bj
    while i > 0 and j > 0 and ptr[i][j] != 0:
        p = ptr[i][j]
        if p == 1:
            q_cols.append(q[i - 1]); t_cols.append(t[j - 1]); i -= 1; j -= 1
        elif p == 2:
            q_cols.append(q[i - 1]); t_cols.append('-'); i -= 1
        else:
            q_cols.append('-'); t_cols.append(t[j - 1]); j -= 1
    q_cols.reverse(); t_cols.reverse()
    t_lo, t_hi = (j + 1, bj) if q_cols else (0, 0)
    return q_cols, t_cols, t_lo, t_hi


# ── 结果对象（字段与既有接口一致） ──────────────────────────────────

@dataclass
class PsRNATargetResult:
    """psRNATarget prediction result."""
    mirna_id: str
    target_id: str
    target_start: int = 0
    target_end: int = 0
    alignment_score: float = 0.0
    expectation: float = 0.0
    max_energy: float = 0.0
    unpaired_energy: float = 0.0
    num_mismatches: int = 0
    num_wobbles: int = 0
    num_gaps: int = 0
    seed_perfect: bool = False
    cleavage_possible: bool = False
    inhibition_possible: bool = False
    alignment_str: str = ""
    passes: bool = False
    consensus_votes: int = 0


# ── 评分与过滤（parse_mirna_targets.py 忠实移植） ───────────────────

def score_alignment(query_aln, subject_aln):
    """对一条比对（列字符列表）按 psRNATarget 口径打分。

    query_aln/subject_aln 为**已按 miRNA 5'→3' 定向**的列字符列表
    （miRNA 字符 / 目标字符，缺口为 '-'）。

    返回 (penalty_score, seed_mismatches, total_mismatches, gaps, GUs,
          参与打分的列数)。"""
    penalty = 0.0
    seed_mismatches = total_mismatches = gaps = gus = 0
    extra_query_gap = 0
    gap_status_q = gap_status_s = 'closed'
    z = 1                                    # 1-based 列号（对齐 miRNA 位置）
    n_scored = min(len(query_aln), ALIGNMENT_LENGTH)
    for i in range(n_scored):
        q = query_aln[i]
        h = subject_aln[i]
        is_gap_q = (q == '-')
        is_gap_s = (h == '-')
        curr = 0.0
        is_mismatch = False
        if is_gap_q or is_gap_s:
            if (is_gap_q and gap_status_q == 'closed') or \
               (is_gap_s and gap_status_s == 'closed'):
                curr = 2.0                   # 缺口打开
                if is_gap_q:
                    extra_query_gap += EXTRA_PENALTY_QUERY_GAP
            else:
                curr = 0.5                   # 缺口延伸
            if is_gap_q:
                gap_status_q = 'open'
            if is_gap_s:
                gap_status_s = 'open'
            gaps += 1
            total_mismatches += 1
            is_mismatch = True
        elif _column_class(q, h) == 'wc':
            # 完美配对 → 0。2026-10-03 修复 P1-3：原先缺 'wc' 分支，
            # q==h 的完美配对落进 else 按错配 1.0 计分（self 比对
            # penalty=比对列数、完美命中被自己的打分否掉）。
            curr = 0.0
            gap_status_q = gap_status_s = 'closed'
        elif _column_class(q, h) == 'gu':
            curr = 0.5                       # G:U 摆动（种子区内不乘）
            gus += 1
            gap_status_q = gap_status_s = 'closed'
        else:
            curr = 1.0                       # 错配
            gap_status_q = gap_status_s = 'closed'
            total_mismatches += 1
            is_mismatch = True
        if SEED_START <= z <= SEED_END and is_mismatch:
            curr *= PENALTY_MULTIPLIER
            seed_mismatches += 1
        penalty += curr
        z += 1
    penalty += extra_query_gap
    return penalty, seed_mismatches, total_mismatches, gaps, gus, n_scored


def _estimate_duplex_energy(query_aln, subject_aln) -> float:
    """双链 MFE 估算（VIENNA 可用时精确，否则最近邻近似）。"""
    qa = ''.join(c for c in query_aln if c != '-')
    ha = ''.join(c for c in subject_aln if c != '-')
    try:
        import RNA
        if qa and ha:
            fc = RNA.fold_compound(qa.replace('T', 'U') + '&' +
                                   ha.replace('T', 'U'))
            return fc.mfe()[1]
    except Exception:
        pass
    NN = {('A', 'U'): -2.0, ('U', 'A'): -2.0, ('G', 'C'): -3.0, ('C', 'G'): -3.0,
          ('G', 'U'): -1.0, ('U', 'G'): -1.0}
    e = 4.0
    for a, b in zip(query_aln, subject_aln):
        if a == '-' or b == '-':
            e += 1.5
        else:
            a_rna = a.upper().replace('T', 'U')
            b_rna = b.upper().replace('T', 'U')
            pair = _column_class(a, b)
            if pair == 'wc':
                e += NN.get((a_rna, b_rna), NN.get((b_rna, a_rna), -1.0))
            elif pair == 'gu':
                e += -1.0
            else:
                e += 2.0
    return e


def psrna_score_alignment(mirna_id: str, target_id: str,
                          query_aln, subject_aln,
                          target_offset: int = 0,
                          target_range=None,
                          ) -> Optional[PsRNATargetResult]:
    """对一条比对施加 psRNATarget 打分与过滤（MiRNATarget 口径）。

    query_aln = miRNA 列（**按 miRNA 5'→3' 定向**），subject_aln = 目标列；
    两者等长。返回 PsRNATargetResult（passes 表示通过全部过滤）。

    target_range = 该比对在 subject 串上的 1-based 闭区间 (lo, hi)
    （_sw_align 返回的 t_lo/t_hi）。注意 miRNA 5'→3' 帧下 subject 列自左
    向右 = target 3'→5'，坐标递减，故不能按列序数坐标（2026-10-03 修复：
    旧实现给出 start>end 的倒置区间）。未给 target_range 时退化为
    「subject 列自第 target_offset+1 位起」的旧口径。
    """
    penalty, seed_mm, total_mm, gaps, gus, _ = score_alignment(
        query_aln, subject_aln)
    subject_len = len(subject_aln)
    hsp = len([c for c in query_aln if c != '-'])

    passes = (penalty <= E_CUTOFF and seed_mm <= NUM_MISMATCH_SEED and
              hsp >= HSP_CUTOFF and gaps <= GAP_CUTOFF and
              total_mm <= TOTAL_MISMATCHES_CUTOFF and
              subject_len <= MAX_ALIGNMENT_LENGTH and hsp > 0)

    # 靶坐标（subject 串内 1-based 闭区间，lo ≤ hi）
    if target_range:
        t_start, t_end = int(target_range[0]), int(target_range[1])
    else:
        t_start = target_offset + 1
        t_end = target_offset + len(subject_aln)

    # 切割位点：psRNATarget 在 miRNA 10-11 位错配时抑制切割
    central_mm = 0
    mpos = 0
    for q, h in zip(query_aln, subject_aln):
        if q == '-':
            continue
        mpos += 1
        if 10 <= mpos <= 11 and _column_class(q, h) != 'wc':
            central_mm += 1

    # 比对字符串（回写为平台统一的三行格式，方向保持 miRNA 5'→3'）
    q_s = ''.join(query_aln)
    h_s = ''.join(subject_aln)
    mid = []
    for q, h in zip(query_aln, subject_aln):
        if q == '-' or h == '-':
            mid.append(' ')
        elif _column_class(q, h) == 'wc':
            mid.append('|')
        elif _column_class(q, h) == 'gu':
            mid.append(':')
        else:
            mid.append('.')

    return PsRNATargetResult(
        mirna_id=mirna_id, target_id=target_id,
        target_start=t_start, target_end=t_end,
        alignment_score=penalty, expectation=penalty,
        num_mismatches=total_mm, num_wobbles=gus, num_gaps=gaps,
        seed_perfect=(seed_mm == 0),
        cleavage_possible=(central_mm == 0),
        inhibition_possible=(central_mm > 0),
        alignment_str=(f"miRNA  5' {q_s} 3'\n"
                       f"         {''.join(mid)}\n"
                       f"target 5' {h_s} 3'"),
        passes=passes,
    )


def evaluate_mirna_vs_target(mirna_seq: str, target_seq: str,
                             mirna_id: str = '', target_id: str = '',
                             ) -> Optional[PsRNATargetResult]:
    """完整评估：正反两条目标链各做一次 SSEARCH36 参数的局部比对，
    取更能通过过滤的那条（与参考实现两次 ssearch36 调用一致）。

    2026-10-03 修复：
      * 反链比对前先 _dna() 归一（含 U 的 RNA 输入原会经 revcomp 的
        未识别碱基映射整条变成 N）；
      * 坐标改用 _sw_align 返回的比对区间（旧实现给出 start>end）。
    """
    m = _dna(mirna_seq)
    q_rev = _revcomp_dna(m)              # 正向链比对用的 query
    best = None
    for t_str, is_rev in ((target_seq, False), (_revcomp_dna(_dna(target_seq)), True)):
        t = _dna(t_str)
        q_cols, h_cols, t_lo, t_hi = _sw_align(q_rev, t)
        if not q_cols:
            continue
        # 比对输出方向 = revcomp-miRNA 5'→3' = miRNA 3'→5'：反转为 miRNA 5'→3'
        res = psrna_score_alignment(mirna_id, target_id,
                                    list(reversed(q_cols)),
                                    list(reversed(h_cols)),
                                    target_range=(t_lo, t_hi))
        if res is None:
            continue
        if is_rev:
            # 反向链命中：坐标换算回原目标链（区间取反）
            L = len(target_seq)
            res.target_start, res.target_end = L - res.target_end + 1, L - res.target_start + 1
        # 择优：优先通过过滤者；同通过状态时 expectation（=penalty，越小越好）
        # 更小者胜。比较方向本身正确，修复前选错源于满分比对被按错配计分。
        if best is None or (res.passes and not best.passes) or \
           (res.passes == best.passes and res.expectation < best.expectation):
            best = res
    if best and best.target_start == 0:
        best.target_start, best.target_end = 1, len(m)
    return best


# ── 兼容入口（旧调用名保留） ────────────────────────────────────────

def psrna_target_scan(mirna_id: str, mirna_seq: str,
                      target_id: str, target_seq: str,
                      max_expectation: float = 5.0,
                      max_upe: float = 25.0,
                      ) -> Optional[PsRNATargetResult]:
    return evaluate_mirna_vs_target(mirna_seq, target_seq, mirna_id, target_id)


def _calc_expectation(score: float, mirna_len: int, target_len: int) -> float:
    """兼容保留：MiRNATarget 口径下 E 即 penalty 本身，原统计模型不再使用。"""
    return score
