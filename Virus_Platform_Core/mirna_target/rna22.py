"""
PAmiRDB v3 — RNA22 引擎（sRNA-target-prediction-windows 参考实现移植版）。

用户指定替换源（2026-09-19）：
  https://github.com/every-algorithm/python/blob/
  db08a602d78bd0fc026de39a4400a41886028e79/bioinformatics/rna22.py

参考实现的口径（逐函数对应，见下方 vendored 区）：
  1. find_seed：在序列上找全部 7-nt 窗口（不含 N）；
  2. score_seed：种子的 GC 含量计数（revcomp 不改变 GC 组成）；
  3. 合格种子 = score_seed > 4（7 nt 中 ≥5 个 GC）；
  4. hairpin 验证：种子后随窗口内存在「两段互补 4-nt 茎 + 4-9 nt 环」。

本引擎将其作为 miRanda 命中的**验证器**：在 miRanda 命中区间附近按参考
口径扫描，存在合格位点且与命中区间重叠 → passes。
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Optional


# ══ 参考实现（vendored from every-algorithm/python bioinformatics/rna22.py） ══

def reverse_complement(seq):
    complement = {'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G'}
    return ''.join(complement.get(base, 'N') for base in reversed(seq))


def find_seed(seq, seed_len=7):
    """Find all subsequences of length seed_len without ambiguous bases."""
    seeds = []
    for i in range(len(seq) - seed_len + 1):
        seg = seq[i:i + seed_len]
        if 'N' not in seg:
            seeds.append((i, seg))
    return seeds


def score_seed(seg):
    """Simple scoring: count G and C bases."""
    return sum(1 for base in seg if base in 'GC')


def find_hairpin(seq, min_loop=4, max_loop=9):
    """Detect simple hairpin structures: two complementary 4-nt stems separated
    by a loop of length between min_loop and max_loop."""
    hairpins = []
    for i in range(len(seq) - 4):
        stem1 = seq[i:i + 4]
        rc = reverse_complement(stem1)
        for j in range(i + min_loop, len(seq) - 3):
            stem2 = seq[j:j + 4]
            if stem2 == rc:
                hairpins.append((i, j, i + 4, j + 4))
    return hairpins


def RNA22(seq):
    """Main function to find potential microRNA binding sites."""
    sites = []
    seeds = find_seed(seq)
    for idx, seg in seeds:
        rc = reverse_complement(seg)
        # Validate with hairpin check（2026-10-03 合并嵌套 if，语义不变）
        if score_seed(rc) > 4 and find_hairpin(seq[idx:idx + len(seg) + 5]):
            sites.append((idx, seg))
    return sites


# ══ 平台验证器接口 ═══════════════════════════════════════════════════

@dataclass
class RNA22Result:
    """RNA22 prediction result.

    miRNA_independent 恒 True（2026-10-03 口径诚实化）：本移植版参考实现
    只分析**靶标序列自身**的结构特征（GC 富集 7-nt 种子 + 下游 hairpin），
    miRNA 序列仅用于展示与种子错配计数，不参与任何筛选/打分——命中不代表
    miRNA 特异性证据，汇总时必须带此声明。
    """
    mirna_id: str
    target_id: str
    mfe: float                    # 参考实现无 MFE 概念，恒 0（保留字段兼容）
    pattern_score: float          # 合格种子的 GC 计数（参考实现的 score_seed）
    target_start: int
    target_end: int
    num_mismatches: int
    num_wobbles: int
    seed_perfect: bool
    alignment_str: str
    passes: bool
    consensus_votes: int = 0
    miRNA_independent: bool = True


def _to_rna(seq: str) -> str:
    return seq.upper().replace('T', 'U')


def validate_with_rna22(miranda_result, mirna_id: str,
                         target_id: str, target_seq: str,
                         ) -> Optional[RNA22Result]:
    """按参考实现口径验证 miRanda 命中：

    在 miRanda 命中区间附近（±8nt）扫描目标序列，存在「GC 富集 7-nt 种子
    （score_seed > 4）+ 下游 hairpin」的合格位点且与命中区间重叠 → 通过。

    2026-10-03 口径诚实化：本引擎**不使用 miRNA 序列**判定（miRNA 只在
    输出里用于展示错配数），故返回对象标注 miRNA_independent=True；
    passes=True 的含义是「靶标自身存在合格结构位点」，不是 miRNA 特异性证据。
    """
    if not target_seq:
        return None
    target_rna = _to_rna(target_seq)
    t_start = max(0, miranda_result.target_start - 1)
    t_end = min(len(target_rna), miranda_result.target_end)
    if t_end <= t_start:
        return None

    # 合格位点全集（参考实现 RNA22(target) 的输出）
    sites = RNA22(target_rna)
    if not sites:
        return None

    hit_lo, hit_hi = t_start, t_end
    best = None                       # 与命中区间重叠的位点里取 GC 计数最高者
    for idx, seg in sites:
        site_lo, site_hi = idx, idx + len(seg)
        if site_hi <= hit_lo or site_lo >= hit_hi:
            continue                  # 不与 miRanda 命中区间重叠
        gc = score_seed(reverse_complement(seg))
        if best is None or gc > best[0]:
            best = (gc, idx, seg)
    if best is None:
        return None

    gc, idx, seg = best
    seed_rna = reverse_complement(seg)
    mism = sum(1 for a, b in zip(_to_rna(miranda_result.mirna_seq)[1:8],
                                 seed_rna) if a != b)
    mid = []
    for a, b in zip(seg, seed_rna):
        mid.append('|' if a == b else ('o' if a in 'GC' and b in 'GC' else '.'))
    aln_str = (f"target 5' {seg} 3'\n"
               f"         {''.join(mid)}\n"
               f"miRNA*  3' {seed_rna} 5'")

    return RNA22Result(
        mirna_id=mirna_id, target_id=target_id,
        mfe=0.0,
        pattern_score=float(gc),
        target_start=idx + 1, target_end=idx + len(seg),
        num_mismatches=mism, num_wobbles=0,
        seed_perfect=(mism == 0),
        alignment_str=aln_str,
        passes=True,
        consensus_votes=0,
        miRNA_independent=True,          # 判定未用 miRNA 序列（口径诚实化）
    )


if __name__ == '__main__':
    # 自检：参考实现的示例行为
    example_seq = "AUGCUAGCUAGCGUAGCUAGCUAGCUAGCUAGCUGAUGC"
    print("Potential sites:", RNA22(example_seq))
