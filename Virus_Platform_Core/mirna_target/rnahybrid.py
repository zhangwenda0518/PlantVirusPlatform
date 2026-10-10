"""
RNAhybrid-style MFE engine for PAmiRDB v3.
Fast miRNA-target duplex free energy calculation.
"""
from dataclasses import dataclass
from typing import Optional


def _is_wc(a: str, b: str) -> bool:
    a, b = a.upper().replace('T', 'U'), b.upper().replace('T', 'U')
    return (a, b) in {('A', 'U'), ('U', 'A'), ('G', 'C'), ('C', 'G')}


def _is_wobble(a: str, b: str) -> bool:
    a, b = a.upper().replace('T', 'U'), b.upper().replace('T', 'U')
    return (a, b) in {('G', 'U'), ('U', 'G')}


NN_ENERGY = {
    ('A', 'U'): -2.0, ('U', 'A'): -2.0,
    ('G', 'C'): -3.0, ('C', 'G'): -3.0,
    ('G', 'U'): -1.0, ('U', 'G'): -1.0,
}


@dataclass
class RNAhybridResult:
    mirna_id: str
    target_id: str
    mfe: float
    target_start: int
    target_end: int
    seed_perfect: bool
    alignment_str: str
    passes: bool
    consensus_votes: int = 0


def rnahybrid_scan(mirna_id: str, mirna_seq: str,
                   target_id: str, target_seq: str,
                   max_energy: float = -20.0,
                   min_seed_match: int = 7,
                   ) -> Optional[RNAhybridResult]:
    """Fast MFE-based target validation.

    ⚠️ 2026-10-03 标注 deprecated（死函数）：返回对象的 passes 恒为 False
    （见函数末尾），调用方拿到的永远是「永不通过」的结果；平台实际使用的
    是 consensus._validate_rnahybrid（对 miRanda 比对列做 MFE 复核）。
    保留仅为兼容旧调用名，请勿在新代码中使用。
    """
    mirna = mirna_seq.upper().replace('T', 'U')
    target = target_seq.upper().replace('T', 'U')
    
    if len(target) < len(mirna):
        return None

    best_mfe = float('inf')
    best_off = 0

    # Slide and compute complementarity energy
    for off in range(len(target) - len(mirna) + 1):
        t_seg = target[off:off + len(mirna)]
        e = 4.0  # initiation
        for a, b in zip(mirna, t_seg):
            if _is_wc(a, b):
                e += NN_ENERGY.get((a, b), -2.0)
            elif _is_wobble(a, b):
                e += -1.0
            else:
                e += 1.5
        if e < best_mfe:
            best_mfe = e
            best_off = off

    if best_mfe > max_energy:
        return None

    # Build alignment
    t_seg = target[best_off:best_off + len(mirna)]
    aln_a, aln_b, pipes = [], [], []
    for a, b in zip(mirna, t_seg):
        aln_a.append(a); aln_b.append(b)
        if _is_wc(a, b): pipes.append('|')
        elif _is_wobble(a, b): pipes.append(':')
        else: pipes.append('.')

    # Seed check
    seed_ok = True
    wc_count = 0
    for i in range(min(len(aln_a), 8)):
        if _is_wc(aln_a[i], aln_b[i]):
            wc_count += 1
    seed_ok = wc_count >= min_seed_match

    aln_str = (f"miRNA  5' {''.join(aln_a)} 3'\n"
               f"         {''.join(pipes)}\n"
               f"target 3' {''.join(aln_b)} 5'")

    return RNAhybridResult(
        mirna_id=mirna_id, target_id=target_id,
        mfe=best_mfe, target_start=best_off,
        target_end=best_off + len(mirna),
        seed_perfect=seed_ok, alignment_str=aln_str,
        passes=False,
    )
