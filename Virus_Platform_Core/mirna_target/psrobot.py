"""
psRobot-style plant miRNA target prediction engine.
Based on: Wu et al. (2012) psRobot: a web-based tool for plant small RNA analysis.

Key features:
  - Simple complementarity scoring (0-5, lower = better, default cutoff=2.5)
  - Designed for plant miRNAs where near-perfect complementarity is expected
  - Target penalty score based on mismatches and G:U wobbles
  - No complex free energy calculation needed
  - Fast: suitable for genome-wide screening

Scoring scheme (from psRobot documentation):
  - Perfect WC match:     0.0
  - G:U wobble:            0.5 penalty
  - Mismatch:              1.0 penalty
  - Cutoff:                ≤2.5 (default)

  ⚠️ 2026-10-03 docstring 修正：原文档声称「Gap: 2.0 penalty」，但本实现是
  等长滑窗逐列比对（zip(miRNA, 片段)），**不产生缺口也从不计 gap 罚分**——
  gap=2.0 从未生效，已删除该声明（对齐真实能力；num_gaps 恒为 0）。
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import Optional


def _is_wc(a: str, b: str) -> bool:
    a, b = a.upper().replace('T', 'U'), b.upper().replace('T', 'U')
    return (a, b) in {('A', 'U'), ('U', 'A'), ('G', 'C'), ('C', 'G')}


def _is_wobble(a: str, b: str) -> bool:
    a, b = a.upper().replace('T', 'U'), b.upper().replace('T', 'U')
    return (a, b) in {('G', 'U'), ('U', 'G')}


@dataclass
class PsRobotResult:
    """psRobot target prediction result."""
    mirna_id: str
    target_id: str
    score: float                # target penalty (0-5, lower=better)
    target_start: int
    target_end: int
    alignment_str: str
    num_mismatches: int
    num_wobbles: int
    num_gaps: int
    passes: bool                 # score <= 2.5


def psrobot_scan(mirna_id: str, mirna_seq: str,
                 target_id: str, target_seq: str,
                 max_score: float = 2.5,
                 ) -> list[PsRobotResult]:
    """
    psRobot-style plant miRNA target prediction.
    
    Slides miRNA across target, computes complementarity penalty score.
    Returns all hits with score <= max_score.
    """
    mirna = mirna_seq.upper().replace('T', 'U')
    target = target_seq.upper().replace('T', 'U')
    results = []
    
    for offset in range(len(target) - len(mirna) + 1):
        t_seg = target[offset:offset + len(mirna)]
        score = 0.0
        mm = wb = gaps = 0
        aln_a, aln_b, pipes = [], [], []
        
        for a, b in zip(mirna, t_seg):
            aln_a.append(a)
            aln_b.append(b)
            if _is_wc(a, b):
                pipes.append('|')
            elif _is_wobble(a, b):
                pipes.append(':')
                score += 0.5
                wb += 1
            else:
                pipes.append('.')
                score += 1.0
                mm += 1
        
        if score <= max_score:
            aln_str = (f"miRNA  5' {''.join(aln_a)} 3'\n"
                       f"         {''.join(pipes)}\n"
                       f"target 3' {''.join(aln_b)} 5'")
            results.append(PsRobotResult(
                mirna_id=mirna_id, target_id=target_id,
                score=score, target_start=offset,
                target_end=offset + len(mirna),
                alignment_str=aln_str,
                num_mismatches=mm, num_wobbles=wb, num_gaps=gaps,
                passes=score <= max_score,
            ))
    
    return results


def validate_with_psrobot(miranda_result, mirna_id: str,
                           target_id: str, target_seq: str,
                           max_score: float = 2.5,
                           ) -> Optional[PsRobotResult]:
    """Validate a miRanda hit with psRobot scoring."""
    t_start = max(0, miranda_result.target_start - 10)
    t_end = min(len(target_seq), miranda_result.target_end + 10)
    region = target_seq[t_start:t_end]
    
    results = psrobot_scan(mirna_id, miranda_result.mirna_seq,
                            target_id, region, max_score)
    if results:
        best = min(results, key=lambda x: x.score)
        best.target_start += t_start
        best.target_end += t_start
        return best
    return None


# ── miRTarBase validation data ────────────────────────────────────

MIRTARBASE_PLANT_MIRNAS = {
    # Experimentally validated plant miRNA-target interactions
    # Source: miRTarBase v8.0 (PMID: 29126174)
    # Filtered for plant antiviral / stress-related
    "ath-miR156": {
        "targets": ["SPL3", "SPL4", "SPL5", "SPL9", "SPL10", "SPL15"],
        "evidence": "Degradome-seq, 5'-RACE",
        "pmid": "19763152",
    },
    "ath-miR159": {
        "targets": ["MYB33", "MYB65", "MYB101"],
        "evidence": "Degradome-seq",
        "pmid": "18676820",
    },
    "ath-miR164": {
        "targets": ["CUC1", "CUC2", "NAC1"],
        "evidence": "5'-RACE, Degradome-seq",
        "pmid": "15269176",
    },
    "ath-miR165/166": {
        "targets": ["PHB", "PHV", "REV", "ATHB8", "ATHB15"],
        "evidence": "Degradome-seq, 5'-RACE",
        "pmid": "15200956",
    },
    "ath-miR168": {
        "targets": ["AGO1"],
        "evidence": "Degradome-seq",
        "pmid": "16669758",
    },
    "ath-miR398": {
        "targets": ["CSD1", "CSD2", "CCS"],
        "evidence": "5'-RACE",
        "pmid": "15937229",
    },
    "ath-miR399": {
        "targets": ["PHO2", "UBC24"],
        "evidence": "5'-RACE, Degradome-seq",
        "pmid": "15937229",
    },
    "ath-miR408": {
        "targets": ["LAC3", "LAC12", "LAC13", "Plantacyanin"],
        "evidence": "Degradome-seq",
        "pmid": "19763152",
    },
    "osa-miR156": {
        "targets": ["SPL3", "SPL12", "SPL14", "IPA1"],
        "evidence": "Degradome-seq, 5'-RACE",
        "pmid": "20547592",
    },
    "osa-miR159": {
        "targets": ["GAMYB", "GAMYBL1"],
        "evidence": "Degradome-seq",
        "pmid": "20547592",
    },
    "osa-miR168": {
        "targets": ["AGO1a", "AGO1b"],
        "evidence": "5'-RACE",
        "pmid": "16669758",
    },
    "osa-miR398": {
        "targets": ["CSD1", "CSD2", "SOD"],
        "evidence": "5'-RACE",
        "pmid": "20547592",
    },
    "osa-miR528": {
        "targets": ["AAO", "LAC", "F-box protein"],
        "evidence": "Degradome-seq",
        "pmid": "20547592",
    },
}

# psRobot scoring table (for reference)
PSROBOT_SCORING = {
    "WC_match": 0.0,
    "GU_wobble": 0.5,
    "mismatch": 1.0,
    "gap": 2.0,          # ⚠️ 未实现：等长滑窗不产生缺口（2026-10-03 如实标注）
    "default_cutoff": 2.5,
}
