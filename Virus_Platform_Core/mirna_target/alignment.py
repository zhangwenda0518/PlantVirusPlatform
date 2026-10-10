"""
Core miRNA-target alignment engine for PAmiRDB v2.
Enhanced with: position-specific binding, similarity%, matching status.

Satish et al. (2019) Scientific Reports 9:4627
"""

from dataclasses import dataclass, field
from typing import Optional

SCORE = {('A','U'):5,('U','A'):5,('G','C'):5,('C','G'):5,('G','U'):2,('U','G'):2}
GAP_OPEN = -9
GAP_EXTEND = -4

# 种子区（miRNA 5' 起算）。本模块口径 2-8（miRanda 经典种子），与
# tapir 的 2-12、psRNATarget 的 2-13 各自参考实现口径不同，互不混用。
SEED_START = 2
SEED_END = 8


def _pair_score(a: str, b: str) -> int:
    a, b = a.upper().replace('T','U'), b.upper().replace('T','U')
    return SCORE.get((a,b), -3)


def _is_wc(a: str, b: str) -> bool:
    """Watson-Crick pair?"""
    a, b = a.upper().replace('T','U'), b.upper().replace('T','U')
    return (a,b) in {('A','U'),('U','A'),('G','C'),('C','G')}


def _is_wobble(a: str, b: str) -> bool:
    """G:U wobble pair?"""
    a, b = a.upper().replace('T','U'), b.upper().replace('T','U')
    return (a,b) in {('G','U'),('U','G')}


@dataclass
class AlignmentResult:
    mirna_id: str
    mirna_seq: str
    target_id: str
    target_start: int
    target_end: int
    score: float
    energy: float
    num_mismatches: int
    num_gaps: int
    seed_perfect: bool
    query_cover: float
    identity: float
    similarity: float = 0.0           # WC + wobble %
    alignment_str: str = ""
    passes_filters: bool = False
    tier_strict: bool = False
    tier_moderate: bool = False
    tier_loose: bool = False
    
    # New PAmiRDB v2 fields
    query_start: int = 1              # Query_S: miRNA 5' start position (1-indexed)
    query_end: int = 0                # Query_E: miRNA 5' end position
    alignment_length: int = 0
    target_protein: str = ""          # Viral protein targeted
    matching_status: str = ""         # Withinmatch / Overlapping
    position_binding: dict = field(default_factory=dict)  # {pos: 'WC'|'wobble'|'mismatch'|'gap'}


def _get_position_binding(aln_a: list, aln_b: list, pipes: list) -> dict:
    """
    Map miRNA positions (1-indexed from 5') to binding status.
    Returns {pos: 'WC'|'wobble'|'mismatch'|'gap'}

    aln_a 为反转帧（miRNA 3'→5' 书写，左端 = miRNA 3' 端），故位置计数
    必须自右向左，才能得到文档承诺的「5' 起算」位置（2026-10-03 修复 P0-4）。
    """
    pos_map = {}
    mirna_pos = 0
    for a, b, p in reversed(list(zip(aln_a, aln_b, pipes))):
        if a != '-':
            mirna_pos += 1
            if p == '|':
                if _is_wc(a, b):
                    pos_map[mirna_pos] = 'WC'
                elif _is_wobble(a, b):
                    pos_map[mirna_pos] = 'wobble'
                else:
                    pos_map[mirna_pos] = 'mismatch'
            elif b == '-':
                pos_map[mirna_pos] = 'gap'
            else:
                pos_map[mirna_pos] = 'mismatch'
    return pos_map


def _check_position_constraint(pos_binding: dict, required_positions: list) -> bool:
    """
    Check if all required positions have complementary binding (WC or wobble).
    required_positions: list of 1-indexed miRNA positions
    """
    for pos in required_positions:
        status = pos_binding.get(pos)
        if status not in ('WC', 'wobble'):
            return False
    return True


def _check_au_at_positions(pos_binding: dict, aln_a: list, positions: list) -> bool:
    """
    Check if specific miRNA positions are A or U.
    positions: list of 1-indexed miRNA positions (5' 起算)
    """
    # Build mapping of miRNA position to nucleotide（aln_a 为反转帧，
    # 自右向左数才是 miRNA 5'→3' 顺序；2026-10-03 修复 P0-4）
    mirna_pos = 0
    nuc_at_pos = {}
    for a in reversed(aln_a):
        if a != '-':
            mirna_pos += 1
            nuc_at_pos[mirna_pos] = a.upper()
    
    for pos in positions:
        if nuc_at_pos.get(pos, '') not in ('A', 'U'):
            return False
    return True


def align_mirna_to_target(
    mirna_id: str, mirna_seq: str,
    target_id: str, target_seq: str,
    required_positions: list = None,  # positions that MUST be paired
    check_au_positions: list = None,  # positions that MUST be A or U
) -> Optional[AlignmentResult]:
    """Run miRanda-style alignment with optional position constraints.

    方向（2026-10-03 修复 P0-4）：miRNA 与靶标是**反平行**配对，物理位点是
    target 5'→3' 上的一段 revcomp(miRNA)。因此 DP 让 miRNA 的反转帧
    （mir_rev = mirna[::-1]，即 miRNA 3'→5' 书写）沿 target 5'→3' 扫描，
    每列的 _pair_score(mir_rev[i], target[j]) 即真实碱基配对；输出比对串与
    miRanda 官方输出的书写约定一致（上行 miRNA 3'→5'、下行 target 5'→3'，
    已用官方 miranda.exe v3.3a 实测输出核对）。

    修复前 DP 直接用 miRNA 5'→3' 原向做同索引互补匹配，等价于找"同向镜像"
    的 comp(miRNA) 段，把物理位点判反（审计三重验证：含 comp(miRNA) 的 target
    identity=100%，含 revcomp(miRNA) 的 target 仅 57%；ViennaRNA duplexfold
    判 revcomp 为物理位点 MFE -37.3 vs -9.4）。

    坐标约定：反转帧位置 k（1..L）对应 miRNA 原始 5'→3' 位置 L-k+1；
    query_start / query_end 输出原始帧（5' 起算、start<end）；
    target_start / target_end 为 target 5'→3' 上的 1-based 闭区间。

    种子区：miRNA 原始 5' 第 SEED_START-SEED_END 位（本模块口径 2-8），
    在反转帧中即第 L-SEED_END+1 .. L-SEED_START+1 位。

    打分矩阵（WC +5 / 错配 -3 / 缺口 -9/-4）非 miRanda 原值（miRanda 原值
    为 match +5 / mismatch -4 / gap open -9 / extend -4 且带 scaling），
    分数字段仅供排序，不可与 miRanda 官方 score 直接比较。
    """
    mirna = mirna_seq.upper().replace('T','U')
    target = target_seq.upper().replace('T','U')
    mir_rev = mirna[::-1]        # 反转帧：上行 = miRNA 3'→5'，与 target 5'→3' 反平行
    L = len(mirna)
    m, n = len(mir_rev), len(target)

    # DP forward（mir_rev 沿 target 5'→3' 扫描）
    H = [[0.0]*(n+1) for _ in range(m+1)]
    best_score, best_pos = 0.0, (0,0)
    for i in range(1,m+1):
        for j in range(1,n+1):
            diag = H[i-1][j-1]+_pair_score(mir_rev[i-1],target[j-1])
            up = H[i-1][j]+(GAP_EXTEND if i>1 and H[i-1][j]>0 else GAP_OPEN+GAP_EXTEND)
            left = H[i][j-1]+(GAP_EXTEND if j>1 and H[i][j-1]>0 else GAP_OPEN+GAP_EXTEND)
            H[i][j]=max(0.0,diag,up,left)
            if H[i][j]>best_score:
                best_score=H[i][j]; best_pos=(i,j)

    if best_score==0:
        return None

    # Traceback
    i,j=best_pos; aln_a,aln_b,pipes=[],[],[]
    mismatches=gaps=0
    while i>0 and j>0 and H[i][j]>0:
        cur=H[i][j]
        if cur==H[i-1][j-1]+_pair_score(mir_rev[i-1],target[j-1]):
            s=_pair_score(mir_rev[i-1],target[j-1])
            aln_a.append(mir_rev[i-1]); aln_b.append(target[j-1])
            pipes.append('|' if s>0 else '.')
            if s<=0: mismatches+=1
            i-=1;j-=1
        elif cur==H[i-1][j]+(GAP_EXTEND if i>1 and H[i-1][j]>0 else GAP_OPEN+GAP_EXTEND):
            aln_a.append(mir_rev[i-1]); aln_b.append('-'); pipes.append(' '); gaps+=1; i-=1
        else:
            aln_a.append('-'); aln_b.append(target[j-1]); pipes.append(' '); gaps+=1; j-=1

    aln_a.reverse(); aln_b.reverse(); pipes.reverse()
    # 回溯止步处 (i, j) 是最后一个未参与比对的列：target 真实区间为
    # [j+1, best_pos[1]]（1-based 闭区间）；修复前写作 t_start=j，整体左偏 1。
    t_start=j+1; t_end=best_pos[1]

    # Position binding map（内部按原始 5' 帧换算，见函数 docstring）
    pos_binding = _get_position_binding(aln_a, aln_b, pipes)

    # If position constraints specified, check them
    if required_positions:
        if not _check_position_constraint(pos_binding, required_positions):
            return None
    if check_au_positions:
        if not _check_au_at_positions(pos_binding, aln_a, check_au_positions):
            return None

    # Seed check：按 miRNA 原始 5' 第 2-8 位判定（wobble 由 '|' 计入通过，
    # 与原逻辑一致）。反转帧自左向右第 k 个非缺口字符 = 原始位置 L-k+1。
    seed_ok=True
    k=0
    for a,b,p in zip(aln_a,aln_b,pipes):
        if a!='-':
            k+=1
            if SEED_START<=L-k+1<=SEED_END and p!='|':
                seed_ok=False; break

    # Energy（在原始 5'→3' 帧上做 cofold：两链各按 5'→3' 读入，
    # 与比对列的反平行配对语义一致；仅能量项依赖此帧）
    energy=_estimate_energy(aln_a[::-1],aln_b)

    # Query_S and Query_E：比对覆盖的 miRNA 位置（原始 5'→3' 帧，1-based）
    k_first=0; k_last=0; k=0
    for a in aln_a:
        if a!='-':
            k+=1
            if k_first==0: k_first=k
            k_last=k
    query_s=L-k_last+1 if k_last else 1
    query_e=L-k_first+1 if k_first else 0

    # Alignment length
    aln_len=len(aln_a)

    # Metrics
    aln_bases=len([c for c in aln_a if c!='-'])
    query_cover=aln_bases/len(mirna_seq)*100 if mirna_seq else 0

    # Identity: WC pairs / miRNA length
    wc_pairs=sum(1 for a,b in zip(aln_a,aln_b) if _is_wc(a,b))
    identity=wc_pairs/len(mirna_seq)*100 if mirna_seq else 0

    # Similarity: (WC + wobble) / miRNA length
    wobble_pairs=sum(1 for a,b in zip(aln_a,aln_b) if _is_wobble(a,b))
    similarity=(wc_pairs+wobble_pairs)/len(mirna_seq)*100 if mirna_seq else 0

    # Consecutive mismatch check
    max_cons=0; cur_cons=0
    for ch in pipes:
        if ch!='|': cur_cons+=1; max_cons=max(max_cons,cur_cons)
        else: cur_cons=0

    # Score threshold from original miRanda C source (all.c)
    # Default: 50, strict (i-miRNA): 140
    SCORE_THRESHOLD = 50
    
    # Tiered filter: support progressive sensitivity
    tier_strict  = (seed_ok and mismatches<=4 and max_cons<=2 and energy<=-20 and best_score>=SCORE_THRESHOLD)
    tier_moderate = (seed_ok and mismatches<=6 and max_cons<=3 and energy<=-15 and best_score>=SCORE_THRESHOLD*0.7)
    tier_loose   = (mismatches<=10 and max_cons<=4 and energy<=-10 and best_score>=SCORE_THRESHOLD*0.4)

    # 书写约定与 miRanda 官方输出一致：上行 miRNA 3'→5'（= 反转帧），
    # 下行 target 5'→3'（官方 miranda.exe 实测格式，见函数 docstring）
    aln_str=(f"miRNA  3' {''.join(aln_a)} 5'\n"
             f"         {''.join(pipes)}\n"
             f"target 5' {''.join(aln_b)} 3'")

    return AlignmentResult(
        mirna_id=mirna_id, mirna_seq=mirna_seq,
        target_id=target_id, target_start=t_start, target_end=t_end,
        score=best_score, energy=energy,
        num_mismatches=mismatches, num_gaps=gaps,
        seed_perfect=seed_ok, query_cover=query_cover,
        identity=identity, similarity=similarity,
        alignment_str=aln_str, passes_filters=tier_strict,
        tier_strict=tier_strict, tier_moderate=tier_moderate, tier_loose=tier_loose,
        query_start=query_s, query_end=query_e,
        alignment_length=aln_len,
        position_binding=pos_binding,
    )


def _estimate_energy(aln_a: list, aln_b: list) -> float:
    """双链 MFE 估算。

    aln_a 需为 miRNA 原始 5'→3' 帧（调用方传入 aln_a[::-1]）——cofold 语义
    是两条链各按 5'→3' 读入，故能量项必须用原始帧，不可用反转帧。
    最近邻回退分支与帧无关（逐列求和）。
    """
    try:
        import RNA
        seq_a=''.join(c for c in aln_a if c!='-')
        seq_b=''.join(c for c in aln_b if c!='-')
        if seq_a and seq_b:
            fc=RNA.fold_compound(seq_a+'&'+seq_b)
            return fc.mfe()[1]
    except ImportError: pass

    NN={('A','U'):-2.0,('U','A'):-2.0,('G','C'):-3.0,('C','G'):-3.0,
        ('G','U'):-1.0,('U','G'):-1.0}
    e=4.0
    for a,b in zip(aln_a,aln_b):
        if a=='-' or b=='-': e+=1.0
        else: e+=NN.get((a.upper(),b.upper()),2.0)
    return e


# Re-export helper for pipeline
def get_position_binding(aln_a, aln_b, pipes):
    return _get_position_binding(list(aln_a), list(aln_b), list(pipes))


def check_position_constraint(pos_binding, positions):
    return _check_position_constraint(pos_binding, positions)
