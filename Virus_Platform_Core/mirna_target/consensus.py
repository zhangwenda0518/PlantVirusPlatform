"""PAmiRDB v3 Consensus Scoring System.
Integrates miRanda + psRNATarget + RNAhybrid via alignment validation.
"""
from __future__ import annotations
import sys, math, logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

log = logging.getLogger(__name__)

try:  # 包内相对导入（Virus_Platform_Core.mirna_target）
    from .alignment import align_mirna_to_target, AlignmentResult
    from .psrnatarget import (PsRNATargetResult, SEED_START, SEED_END,
                              PENALTY_MULTIPLIER, E_CUTOFF,
                              NUM_MISMATCH_SEED, HSP_CUTOFF, GAP_CUTOFF,
                              TOTAL_MISMATCHES_CUTOFF, ALIGNMENT_LENGTH)
    from .rnahybrid import RNAhybridResult
    from .rna22 import RNA22Result, validate_with_rna22
    from .tapir import TAPIRResult, validate_with_tapir
    from .psrobot import PsRobotResult, validate_with_psrobot
except ImportError:  # 独立脚本方式运行（原 PAmiRDB 布局）
    sys.path.insert(0, str(Path(__file__).parent))
    from alignment import align_mirna_to_target, AlignmentResult
    from psrnatarget import (PsRNATargetResult, SEED_START, SEED_END,
                             PENALTY_MULTIPLIER, E_CUTOFF,
                             NUM_MISMATCH_SEED, HSP_CUTOFF, GAP_CUTOFF,
                             TOTAL_MISMATCHES_CUTOFF, ALIGNMENT_LENGTH)
    from rnahybrid import RNAhybridResult
    from rna22 import RNA22Result, validate_with_rna22
    from tapir import TAPIRResult, validate_with_tapir
    from psrobot import PsRobotResult, validate_with_psrobot


# ── Seed Region Quality Scoring (from miRNA-seed-matching-counter) ───

def compute_seed_quality(alignment_str: str) -> dict:
    """
    Analyze seed region (positions 2-8) complementarity from alignment.

    比对串为 miRanda 书写约定（上行 miRNA 3'→5'、下行 target 5'→3'，
    见 alignment.align_mirna_to_target），故位置计数必须自右向左才是
    miRNA 5' 起算（2026-10-03 修复 P0-4 联动）。逐列按碱基对分类——
    本模块比对串的 '|' 同时覆盖 WC 与 G:U wobble，只看杠线会把 wobble
    误计成 WC。

    Returns:
      seed_wc:    number of Watson-Crick pairs in seed (0-7)
      seed_wobble: number of G:U wobble pairs in seed (0-7)
      seed_mismatch: number of mismatches/gaps in seed (0-7)
      seed_score:    weighted quality (0-100), WC=15pt, wobble=5pt each
      is_seed_perfect: True if all 7 positions are WC pairs
    """
    lines = alignment_str.split('\n')
    if len(lines) < 3:
        return {"seed_wc": 0, "seed_wobble": 0, "seed_mismatch": 7,
                "seed_score": 0, "is_seed_perfect": False}

    # Parse the miRNA and target lines from miRanda format
    pipe_line = lines[1].strip()
    mirna_line = lines[0].replace("miRNA  3' ", "").replace(" 5'", "").strip()
    target_line = lines[2].replace("target 5' ", "").replace(" 3'", "").strip()

    # Count miRNA positions through the alignment (skip gaps in miRNA)
    seed_wc = 0
    seed_wobble = 0
    seed_mismatch = 0
    mirna_pos = 0

    for m_char, t_char, p_char in reversed(list(zip(mirna_line, target_line,
                                                    pipe_line))):
        if m_char == '-':
            continue
        mirna_pos += 1
        if not (2 <= mirna_pos <= 8):
            continue
        if t_char == '-':
            seed_mismatch += 1
        elif _is_wc(m_char, t_char):
            seed_wc += 1
        elif _is_wobble(m_char, t_char):
            seed_wobble += 1
        else:
            seed_mismatch += 1

    # Score: WC = 15pt each (max 105, capped at 100), wobble = 5pt each
    seed_score = min(100, seed_wc * 15 + seed_wobble * 5)

    return {
        "seed_wc": seed_wc,
        "seed_wobble": seed_wobble,
        "seed_mismatch": seed_mismatch + max(0, 7 - seed_wc - seed_wobble - seed_mismatch),
        "seed_score": seed_score,
        "is_seed_perfect": (seed_wc == 7),
    }


def _is_wc(a: str, b: str) -> bool:
    a, b = a.upper().replace('T', 'U'), b.upper().replace('T', 'U')
    return (a, b) in {('A', 'U'), ('U', 'A'), ('G', 'C'), ('C', 'G')}


def _is_wobble(a: str, b: str) -> bool:
    a, b = a.upper().replace('T', 'U'), b.upper().replace('T', 'U')
    return (a, b) in {('G', 'U'), ('U', 'G')}


def _validate_psrna(mr, mirna_id: str, target_id: str) -> Optional[PsRNATargetResult]:
    """Validate miRanda alignment with psRNATarget scoring.

    2026-09-19 改为 MiRNATarget（jtremblay/psRNATarget_command_line）忠实
    移植口径：从 miRNA 5' 端起逐列打分（WC 0 / G:U 0.5 不乘种子系数 /
    错配 1 / 缺口开 2、延 0.5，种子区 2-13 内错配与缺口 ×1.5，miRNA 一侧
    缺口额外 +1），只计前 19 列；过滤：penalty ≤ 5 且种子错配 ≤ 2 且
    HSP ≥ 14 且 gap ≤ 1 且总错配 ≤ 8 且比对长 ≤ 22。切割判定 = miRNA
    10-11 位无错配。

    2026-10-03 修复 P0-4 联动：比对串改为 miRanda 官方书写约定（上行
    miRNA 3'→5'），「自 miRNA 5' 端起」= 自右向左遍历比对列。"""
    lines = mr.alignment_str.split('\n')
    if len(lines) < 3: return None
    mline = lines[0].replace("miRNA  3' ", "").replace(" 5'", "").strip()
    tline = lines[2].replace("target 5' ", "").replace(" 3'", "").strip()

    # ---- MiRNATarget 逐列打分（z = 比对列号，含缺口列，前 19 列） ----
    penalty = 0.0
    seed_mm = 0
    total_mm = 0
    gaps = 0
    gus = 0
    extra_query_gap = 0
    gap_status_q = gap_status_s = 'closed'
    hsp = 0                       # miRNA 非缺口字符数（HSP 判定用）
    ctr = 0                       # 10-11 位错配计数
    col = 0                       # 比对列号（1-based，自 miRNA 5' 端起）
    for m, t in zip(mline[::-1], tline[::-1]):
        col += 1
        if col > ALIGNMENT_LENGTH:
            break
        is_gap_q = (m == '-')
        is_gap_s = (t == '-')
        curr = 0.0
        is_mismatch = False
        if is_gap_q or is_gap_s:
            if (is_gap_q and gap_status_q == 'closed') or \
               (is_gap_s and gap_status_s == 'closed'):
                curr = 2.0
                if is_gap_q:
                    extra_query_gap += 1
            else:
                curr = 0.5
            if is_gap_q:
                gap_status_q = 'open'
            if is_gap_s:
                gap_status_s = 'open'
            gaps += 1
            total_mm += 1
            is_mismatch = True
        else:
            if not is_gap_q:
                hsp += 1
            m_rna = m.upper().replace('T', 'U')
            t_rna = t.upper().replace('T', 'U')
            if (m_rna, t_rna) in {('A', 'U'), ('U', 'A'), ('G', 'C'), ('C', 'G')}:
                pass                              # WC → 0
            elif (m_rna, t_rna) in {('G', 'U'), ('U', 'G')}:
                curr = 0.5
                gus += 1
                gap_status_q = gap_status_s = 'closed'
            else:
                curr = 1.0
                total_mm += 1
                is_mismatch = True
                gap_status_q = gap_status_s = 'closed'
            if 10 <= col <= 11 and curr != 0.0:
                ctr += 1
        if SEED_START <= col <= SEED_END and is_mismatch:
            curr *= PENALTY_MULTIPLIER
            seed_mm += 1
        penalty += curr
    penalty += extra_query_gap

    seed = (seed_mm == 0)
    passes = (penalty <= E_CUTOFF and seed_mm <= NUM_MISMATCH_SEED and
              hsp >= HSP_CUTOFF and gaps <= GAP_CUTOFF and
              total_mm <= TOTAL_MISMATCHES_CUTOFF)

    return PsRNATargetResult(
        mirna_id=mirna_id, target_id=target_id,
        target_start=mr.target_start, target_end=mr.target_end,
        alignment_score=penalty, expectation=penalty, max_energy=mr.energy,
        num_mismatches=total_mm, num_wobbles=gus, num_gaps=gaps,
        seed_perfect=seed,
        cleavage_possible=(ctr == 0), inhibition_possible=(ctr >= 1),
        alignment_str=mr.alignment_str, passes=passes,
    )


def _validate_rnahybrid(mr, mirna_id: str, target_id: str) -> Optional[RNAhybridResult]:
    """Validate miRanda alignment with RNAhybrid MFE criteria.

    2026-10-03 修复 P0-4 联动：比对串为 miRanda 官方书写约定（上行
    miRNA 3'→5'），种子区（miRNA 5' 第 2-8 位）须自右向左数。"""
    lines = mr.alignment_str.split('\n')
    if len(lines) < 3: return None
    mline = lines[0].replace("miRNA  3' ", "").replace(" 5'", "").strip()
    tline = lines[2].replace("target 5' ", "").replace(" 3'", "").strip()
    
    mfe = 4.0; swc = 0; mpos = 0
    for m, t in zip(mline[::-1], tline[::-1]):
        if m != '-' and t != '-':
            if _is_wc(m, t): mfe += -2.0
            elif _is_wobble(m, t): mfe += -1.0
            else: mfe += 1.5
        else: mfe += 1.0
        if m != '-':
            mpos += 1
            if 2 <= mpos <= 8 and _is_wc(m, t): swc += 1
    
    seed_ok = swc >= 6
    passes = mfe <= -20.0 and seed_ok
    
    return RNAhybridResult(
        mirna_id=mirna_id, target_id=target_id, mfe=mfe,
        target_start=mr.target_start, target_end=mr.target_end,
        seed_perfect=seed_ok, alignment_str=mr.alignment_str, passes=passes,
    )


@dataclass
class ConsensusResult:
    mirna_id: str; mirna_seq: str; target_id: str
    target_name: str = ""; target_family: str = ""
    miranda: Optional[AlignmentResult] = None
    psrna: Optional[PsRNATargetResult] = None
    rnahyb: Optional[RNAhybridResult] = None
    rna22: Optional[RNA22Result] = None
    tapir: Optional[TAPIRResult] = None
    psrobot: Optional[PsRobotResult] = None
    num_algorithms: int = 0; consensus_level: str = ""
    best_energy: float = 0.0; best_identity: float = 0.0
    seed_quality_score: float = 0.0; seed_perfect_wc: int = 0; seed_wobble_count: int = 0
    consensus_detail: str = ""; passes_pamirdb_filter: bool = False
    engines_actual: list = field(default_factory=list)   # 实际跑成并投票的引擎键名
    engines_skipped: list = field(default_factory=list)  # 失败/异常引擎 + 原因
    engine_notes: list = field(default_factory=list)     # 引擎口径旁注（诚实化声明）


# 引擎键名 ↔ 投票名（algorithms 参数用键名，consensus_detail 用投票名）
ENGINE_KEYS = ('miranda', 'psrna', 'rnahybrid', 'rna22', 'tapir', 'psrobot')


def run_consensus_prediction(mirna_id, mirna_seq, target_id, target_name, target_family, target_seq,
                              fast_mode=True, algorithms=None):
    """
    Run consensus prediction with configurable algorithm depth.

    fast_mode=True (default): only miRanda+psRNATarget+RNAhybrid (core 3, fastest)
    fast_mode=False: full 6-algorithm chain with RNA22/TAPIR/psRobot
    algorithms: 引擎键名子集（ENGINE_KEYS），给出时优先生效——miRanda 恒为基座
    （无 miRanda 命中直接返回 None），其余只跑勾选的校验器。None = 沿用 fast_mode。
    """
    if algorithms is not None:
        selected = {a for a in algorithms if a in ENGINE_KEYS}
        want = {k: (k == 'miranda' or k in selected) for k in ENGINE_KEYS}
    else:
        want = None  # 原始 fast_mode 语义

    result = ConsensusResult(mirna_id=mirna_id, mirna_seq=mirna_seq,
                              target_id=target_id, target_name=target_name, target_family=target_family)
    algos = []

    r1 = align_mirna_to_target(mirna_id, mirna_seq, target_id, target_seq)
    if r1 and r1.passes_filters:
        result.miranda = r1; algos.append('miRanda')

    r2 = r3 = r4 = r5 = r6 = None
    if 'miRanda' in algos:
        result.engines_actual.append('miranda')
        run_slow = (want is None and not fast_mode) or \
                   (want is not None and (want['rna22'] or want['tapir'] or want['psrobot']))

        if want is None or want['psrna']:
            r2 = _validate_psrna(r1, mirna_id, target_id)
            if r2 and r2.passes:
                result.psrna = r2; r2.consensus_votes = 1; algos.append('psRNATarget')
                result.engines_actual.append('psrna')

        if want is None or want['rnahybrid']:
            r3 = _validate_rnahybrid(r1, mirna_id, target_id)
            if r3 and r3.passes:
                result.rnahyb = r3; algos.append('RNAhybrid')
                result.engines_actual.append('rnahybrid')

        # Slow algorithms — only when requested (or in original full mode)
        if run_slow:
            if want is None or want['rna22']:
                r4 = validate_with_rna22(r1, mirna_id, target_id, target_seq)
                if r4 and r4.passes:
                    result.rna22 = r4; algos.append('RNA22')
                    result.engines_actual.append('rna22')
                # 口径诚实化：本 RNA22 移植只看靶标结构特征（GC 富集种子 +
                # hairpin），完全不使用 miRNA 序列——命中不代表 miRNA 相关。
                result.engine_notes.append(
                    'RNA22：仅靶标结构特征（GC 富集种子 + hairpin 结构），'
                    '未使用 miRNA 序列，命中不能视为 miRNA 特异性证据')

            if want is None or want['tapir']:
                try:
                    r5 = validate_with_tapir(r1, mirna_id, target_id, target_seq)
                    if r5 and r5.passes:
                        result.tapir = r5; algos.append('TAPIR')
                        result.engines_actual.append('tapir')
                except Exception as e:
                    log.warning("TAPIR 校验失败跳过: %r", e)
                    result.engines_skipped.append(
                        {'engine': 'tapir', 'reason': f'{type(e).__name__}: {e}'})

            if want is None or want['psrobot']:
                try:
                    r6 = validate_with_psrobot(r1, mirna_id, target_id, target_seq)
                    if r6 and r6.passes:
                        result.psrobot = r6; algos.append('psRobot')
                        result.engines_actual.append('psrobot')
                except Exception as e:
                    log.warning("psRobot 校验失败跳过: %r", e)
                    result.engines_skipped.append(
                        {'engine': 'psrobot', 'reason': f'{type(e).__name__}: {e}'})


    if not algos: return None

    result.num_algorithms = len(algos)
    result.consensus_detail = '+'.join(algos)
    LEVEL_MAP = {1:'LOW', 2:'LOW+', 3:'MEDIUM', 4:'MEDIUM+', 5:'HIGH',
                 6:'HIGH+', 7:'VERY HIGH', 8:'VERY HIGH'}
    result.consensus_level = LEVEL_MAP.get(len(algos), 'LOW')

    # best_energy：只在「通过过滤的引擎」能量中取 min。被过滤掉的引擎
    # 的 mfe/max_energy 从未参与投票，混进来会把 best_energy 拉低到
    # 非命中值（2026-10-03 修复 P3-15）。
    es = [r1.energy] if (r1 and result.miranda) else []
    if result.psrna: es.append(getattr(result.psrna, 'max_energy', 0))
    if result.rnahyb: es.append(getattr(result.rnahyb, 'mfe', 0))
    result.best_energy = min(es) if es else None   # 通过数=0 时不伪造 0.0
    result.best_identity = r1.identity if r1 else 0.0
    result.passes_pamirdb_filter = result.miranda is not None

    # Seed quality scoring (from miRNA-seed-matching-counter methodology)
    if r1 and r1.alignment_str:
        sq = compute_seed_quality(r1.alignment_str)
        result.seed_quality_score = sq["seed_score"]
        result.seed_perfect_wc = sq["seed_wc"]
        result.seed_wobble_count = sq["seed_wobble"]
    return result


def to_dict(cr: ConsensusResult) -> dict:
    d = {"mirna_id": cr.mirna_id, "mirna_seq": cr.mirna_seq,
         "target_virus_id": cr.target_id, "target_virus_name": cr.target_name,
         "target_virus_family": cr.target_family,
         "num_algorithms": cr.num_algorithms, "consensus_level": cr.consensus_level,
         "consensus_detail": cr.consensus_detail,
         "best_energy": cr.best_energy, "best_identity": cr.best_identity,
         "seed_quality_score": cr.seed_quality_score,
         "seed_perfect_wc": cr.seed_perfect_wc,
         "seed_wobble_count": cr.seed_wobble_count,
         # 实际参与的引擎（engines_actual）与失败跳过的引擎（engines_skipped）
         # —— 后者用于区分「引擎跑失败」与「引擎跑了但未通过过滤」，
         # 避免用请求引擎列表冒充实际结果（2026-10-03 修复 P1-4）。
         "engines_actual": list(cr.engines_actual),
         "engines_skipped": list(cr.engines_skipped),
         "engine_notes": list(cr.engine_notes),
         "passes_pamirdb_filter": cr.passes_pamirdb_filter}
    if cr.miranda:
        r = cr.miranda
        d.update(seed_perfect=r.seed_perfect, query_cover=r.query_cover,
                 query_start=r.query_start, query_end=r.query_end,
                 alignment_length=r.alignment_length,
                 target_start=r.target_start, target_end=r.target_end,
                 miranda_score=r.score, miranda_energy=r.energy,
                 miranda_identity=r.identity, miranda_similarity=r.similarity,
                 miranda_mismatches=r.num_mismatches, miranda_gaps=r.num_gaps,
                 alignment=r.alignment_str)
    if cr.psrna:
        r = cr.psrna
        d.update(psrna_score=r.alignment_score, psrna_expectation=r.expectation,
                 psrna_energy=r.max_energy, psrna_cleavage=r.cleavage_possible,
                 psrna_inhibition=r.inhibition_possible, psrna_alignment=r.alignment_str)
    if cr.rnahyb:
        r = cr.rnahyb
        d.update(rnahybrid_mfe=r.mfe, rnahybrid_seed=r.seed_perfect,
                 rnahybrid_alignment=r.alignment_str)
    if cr.rna22:
        r = cr.rna22
        d.update(rna22_mfe=r.mfe, rna22_pattern_score=r.pattern_score,
                 rna22_seed=r.seed_perfect, rna22_alignment=r.alignment_str,
                 # 口径诚实化：该引擎判定未使用 miRNA 序列
                 rna22_miRNA_independent=getattr(r, 'miRNA_independent', True))
    if cr.tapir:
        r = cr.tapir
        d.update(tapir_score=r.tapir_score, tapir_accessibility=r.accessibility_score,
                 tapir_seed=r.seed_score, tapir_up_energy=r.up_energy)
    if cr.psrobot:  # 平台移植补丁：原版 to_dict 漏导出 psrobot，前端算法表显示 SKIP
        r = cr.psrobot
        d.update(psrobot_score=r.score, psrobot_mismatches=r.num_mismatches,
                 psrobot_wobbles=r.num_wobbles, psrobot_alignment=r.alignment_str)
    return d


# ── Union/Intersection/Individual stratified output (Paper 2) ─────

def stratified_consensus(results: list[dict]) -> dict:
    """
    Generate three-layer consensus output per paper 2 methodology:
      - Individual:  per-algorithm predictions
      - Union:       all unique miRNA-virus pairs (high sensitivity)
      - Intersection: pairs found by >=2, >=3, >=4 tools (high specificity)
    """
    # Individual
    individual = {}
    for r in results:
        detail = r.get('consensus_detail', 'miRanda')
        for tool in detail.split('+'):
            tool = tool.strip()
            individual.setdefault(tool, []).append(r['mirna_id'] + '|' + r.get('target_virus_id', r.get('virus_id', '')))

    # Union (deduplicated miRNA-virus pairs)
    union = list(set(r['mirna_id'] + '|' + r.get('target_virus_id', r.get('virus_id', ''))
                     for r in results))

    # Intersection by algorithm count
    intersection = {}
    for n in [2, 3, 4, 5]:
        pairs = [r for r in results if r.get('num_algorithms', 0) >= n]
        intersection[f'>={n}_tools'] = {
            'count': len(pairs),
            'pairs': [r['mirna_id'] + '|' + r.get('target_virus_id', r.get('virus_id', ''))
                       for r in pairs],
        }

    return {
        'individual': {tool: {'count': len(pairs)} for tool, pairs in individual.items()},
        'union': {'count': len(union)},
        'intersection': intersection,
    }


if __name__ == '__main__':
    cr = run_consensus_prediction(
        "ath-miR408-5p", "CGGGGAACAGGCAGAGCAUGG",
        "ACMV", "African cassava mosaic virus", "geminiviridae",
        "GTCCCCTTCGTCCCTACGAGCT")
    if cr:
        s = {"HIGH":"***","MEDIUM":"** ","LOW":"*  "}.get(cr.consensus_level,"?  ")
        print(f"{s} [{cr.consensus_level}] {cr.mirna_id} -> {cr.target_name}")
        print(f"   Algorithms: {cr.consensus_detail} ({cr.num_algorithms}/3)")
        print(f"   Energy: {cr.best_energy:.1f}")
    else:
        print("No consensus")
