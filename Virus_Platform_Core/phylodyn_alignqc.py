# -*- coding: utf-8 -*-
"""比对质量检查（QC）：长度 / gap% / 简并字符% / identity → 干净集 + 剔除集 + 报告。

来源与定位
----------
口径移植自上游 `MMPV-RNA/virome_phylo_pipeline/utils/align_qc.py`
（MAFFT → 逐序列质量检查（长度/gap%/N%/与参考 identity）→ 去除不完整序列 → QC 报告）。
**为什么平台需要它**：A0（时间与地理推断）的硬前置是「**等长比对**」，
而 explorer 导出的是**未比对**序列 —— 此前只能在 A0 里报"不等长"，
没有任何量化体检，用户不知道是哪几条拖后腿、该剔谁。

与上游的**有意差异**（都在这里写明，不留暗差）
------------------------------------------
  1. **不自己 shell 调 MAFFT**，复用平台 `phylo._run_mafft`
     （它已处理"Windows 原生 mafft 不认中文路径"的 ASCII 中转，以及
     "序列数中等 + 位点数大"时自动从迭代法降级到 FFT-NS-2，避免跑不动）。
  2. 额外给**「是否等长」**与长度分布（A0 的硬前置），上游只给"长度/参考长度"。
  3. identity 默认**对多数票共有序列（consensus）**算，没给参考也能跑；
     给了 `reference` 就按上游口径对该参考算。
  4. 判定阈值默认沿用上游（gap 10% / N 5% / 长度 ≥90% / identity 跌 >2%），
     但**全部可配**并写进报告 —— 阈值是数据相关的，不许藏在代码里。
"""

import io
import os
from collections import Counter, OrderedDict
from typing import Dict, List, Optional, Tuple

# 默认阈值（= 上游 align_qc 的默认值）
DEFAULTS = {'max_gap': 0.10, 'max_n': 0.05, 'min_length': 0.90,
            'max_identity_drop': 0.02}

_CANON = set('ACGTRYSWKMBDHVN')       # DNA（含简并）
_ALPHA = set('ACGTU')                 # 明确的碱基
_GAP_CHARS = '-.?'                    # gap 口径：`-` 与比对软件常用的 `.` / `?`


def read_fasta(path) -> "OrderedDict[str, str]":
    """读 FASTA → 有序 {名字(首词): 序列(大写, 无空白)}。

    ⚠️ 2026-10-03：首词重名**不再静默覆盖** —— 原实现 `seqs[name] = …` 会让
    后写的把先写的顶掉，序列条数凭空少掉、QC 剔除集对不上。第二条起改名
    `名字.2` / `名字.3` … 保留，重名明细见 `read_fasta_detailed`。
    """
    return read_fasta_detailed(path)[0]


def read_fasta_detailed(path):
    """同 `read_fasta`，但额外返回 (重名明细 [{name, kept_as}, …])。

    重名不是小事：explorer 导出的 FASTA 头常是 `>acc 描述`，首词可能撞车；
    丢掉任何一条都会让 QC 的条数与剔除集对不上（静默丢数据）。
    """
    seqs: "OrderedDict[str, str]" = OrderedDict()
    dups: List[Dict] = []
    name = None
    buf: List[str] = []
    with io.open(path, encoding='utf-8', errors='replace') as f:
        for ln in f:
            ln = ln.rstrip('\r\n')
            if ln.startswith('>'):
                if name is not None:
                    _put_seq(seqs, dups, name, ''.join(buf).upper())
                tok = ln[1:].strip().split()
                name = tok[0] if tok else ''
                buf = []
            elif ln.strip():
                buf.append(''.join(ln.split()))
    if name is not None:
        _put_seq(seqs, dups, name, ''.join(buf).upper())
    return seqs, dups


def _put_seq(seqs, dups, name, seq):
    """存一条序列；首词重名时加 `.2`/`.3`… 后缀（不覆盖、并留痕）。"""
    nm = name
    if nm in seqs:
        k = 2
        while '%s.%d' % (name, k) in seqs:
            k += 1
        nm = '%s.%d' % (name, k)
        dups.append({'name': name, 'kept_as': nm})
    seqs[nm] = seq


# ══════════════════════════════════════════════════════════════════════
# 单序列指标
# ══════════════════════════════════════════════════════════════════════

def seq_metrics(name: str, seq: str) -> Dict:
    """一条序列的 QC 指标（不含判定）。

    ⚠️ 口径（2026-10-03 修）：`-` / `.` / `?` **统一算 gap**（原实现把 `.`
    既算 gap 又算简并字符 → `n_frac` 与 `gap_frac` 对同一批 `.` 双计）；
    `n_frac` 只统计 **IUPAC 简并码 + N**（R/Y/S/W/K/M/B/D/H/V/N），
    真正非法的字符（X / `*` / 数字等）仍只进 `degenerate`，不冒充 gap。
    """
    n = len(seq)
    if n == 0:
        return {'name': name, 'length': 0, 'gap': 0, 'gap_frac': 1.0,
                'n_frac': 1.0, 'n_amb': 0, 'degenerate': 0}
    gaps = sum(1 for c in seq if c in _GAP_CHARS)
    amb = sum(1 for c in seq if c in _CANON and c not in _ALPHA)  # 简并/N
    bad = sum(1 for c in seq if c not in _CANON and c not in _GAP_CHARS)
    return {'name': name, 'length': n, 'gap': gaps, 'gap_frac': gaps / n,
            'n_frac': amb / n, 'n_amb': amb, 'degenerate': bad}


def consensus(seqs: "OrderedDict[str, str]") -> str:
    """按列多数票求共有序列（**只统计明确碱基**；并列时取字典序，保证可复现）。"""
    if not seqs:
        return ''
    L = max(len(s) for s in seqs.values())
    out = []
    for i in range(L):
        c: Counter = Counter()
        for s in seqs.values():
            ch = s[i] if i < len(s) else '-'
            if ch in _ALPHA or ch == '-':
                c[ch] += 1
        if not c:
            out.append('N')
            continue
        top = max(c.values())
        out.append(sorted(k for k, v in c.items() if v == top)[0])
    return ''.join(out)


def pairwise_identity(a: str, b: str) -> Optional[float]:
    """两位点均非 gap 时才计入分母（与上游 `pairwise_identity` 同口径）。"""
    if not a or not b:
        return None
    m = min(len(a), len(b))
    same = tot = 0
    for i in range(m):
        x, y = a[i], b[i]
        if x == '-' or y == '-':
            continue
        tot += 1
        if x == y:
            same += 1
    return (same / tot) if tot else None


# ══════════════════════════════════════════════════════════════════════
# 表级 QC
# ══════════════════════════════════════════════════════════════════════

def qc(fasta: str, *, reference: Optional[str] = None,
       max_gap: float = DEFAULTS['max_gap'], max_n: float = DEFAULTS['max_n'],
       min_length: float = DEFAULTS['min_length'],
       max_identity_drop: float = DEFAULTS['max_identity_drop'],
       log=None) -> Dict:
    """对**已有**比对/序列做 QC，返回报告 dict（不写盘、不剔除）。

    判定（每条**只给第一个**命中的理由，便于报告短而准）：
      · `too_short`     长度 < `min_length` × 参考长度（参考 = 众数长度或显式 `reference`）
      · `too_gappy`     gap 比例 > `max_gap`（比对质量差）
      · `too_ambiguous` 简并字符比例 > `max_n`（低质量）
      · `identity_low`  与参考的 identity 跌破（中位数 − `max_identity_drop`）（偏离群体）
      · `ok`
    """
    say = log or (lambda *_: None)
    seqs, dup_names = read_fasta_detailed(fasta)
    if not seqs:
        return {'ok': False, 'error': 'FASTA 是空的', 'n_seq': 0}
    lengths = [len(s) for s in seqs.values()]
    modal = Counter(lengths).most_common(1)[0][0]
    ref_len = len(seqs[reference]) if (reference and reference in seqs) else modal
    aligned = len(set(lengths)) == 1

    ref_seq = seqs[reference] if (reference and reference in seqs) else None
    if ref_seq is None:
        ref_seq = consensus(seqs)
        ref_kind = 'consensus(多数票共有序列)'
    else:
        ref_kind = f'指定参考 {reference}'

    # 先算 identity 的中位数（用于 identity_low 判据）
    # ⚠️ 两个口径都要：
    #   `identity`      与上游一致 —— 只在**双方都有碱基**的列上比；
    #   `core_identity` 只在**多数序列都有碱基**的列上比（≈不用 gap 块）。
    #   为什么必须加后者：当比对里出现大 gap 块（成片插入/缺失）时，
    #   `identity` 会遇到"对角块"假象 —— 每列在场的就是一个小亚群，
    #   多数票自然等于该亚群自己的碱基 → 再离谱的序列也能报 0.99。
    #   实测（用户那批 359bp 类病毒）：`identity` 报 0.94~0.997，
    #   而真实两两比只有 **0.66** —— 差到会让人做出错误判断。
    ids, core_ids, ncol = {}, {}, max(len(s) for s in seqs.values())
    col_n = [0] * ncol
    for s in seqs.values():
        for i, ch in enumerate(s):
            if i < ncol and ch != '-':
                col_n[i] += 1
    core_cols = [i for i in range(ncol) if col_n[i] >= 0.9 * len(seqs)]
    for nm, s in seqs.items():
        if nm == reference:
            continue
        ids[nm] = pairwise_identity(s, ref_seq)
        if core_cols:
            core_ids[nm] = pairwise_identity(
                ''.join(s[i] if i < len(s) else '-' for i in core_cols),
                ''.join(ref_seq[i] if i < len(ref_seq) else '-' for i in core_cols))
    vals_core = sorted(v for v in core_ids.values() if v is not None)
    vals_ids = sorted(v for v in ids.values() if v is not None)
    # ⚠️ 两个中位数**分别算、分别报**（2026-10-03 修）：原实现把同一个 med 同时
    #    写进 `identity_median` 与 `core_identity_median`，报告里两列永远相等，
    #    读者会以为"两种口径一致"，其实只是根本没算第二个。
    #    判 identity_low 的阈值仍用 core 优先、无 core 列时退回 id 口径（行为不变）。
    med_ids = vals_ids[len(vals_ids) // 2] if vals_ids else None
    med_core = vals_core[len(vals_core) // 2] if vals_core else None
    med = med_core if med_core is not None else med_ids
    n_big_gap_cols = sum(1 for c in col_n if c < 0.5 * len(seqs))
    # ⚠️ **gap 判定不能直接用整条序列的 gap 比例**：多条比对的 gap 是"共享"的 ——
    #    只要**少数**序列带插入/移位，MAFFT 就会把比对撑长，而**其余每条**都被追加大段
    #    gap。实测（用户那批 PSTVd 全长的 171 条）：末尾 ~90 列只有极少数序列在场，
    #    于是 171 条**全部**被判 `too_gappy`（每条 ~22%）—— 其实是"1 条离群、170 条无辜"。
    #    所以判定改用 **只统计"≥50% 序列在场"的列** 里的 gap（`gap_frac_core`）：
    #    成片 block 由集合级诊断（`block_hint` + 下面的 outlier 列表）负责，
    #    不该算在每条序列头上。`gap_frac` 仍原样保留在报告里供人工核对。
    present_cols = [i for i in range(ncol) if col_n[i] >= 0.5 * len(seqs)]
    n_present = max(1, len(present_cols))

    rows, keep, drop = [], [], []
    for nm, s in seqs.items():
        m = seq_metrics(nm, s)
        m['identity'] = ids.get(nm)
        m['core_identity'] = core_ids.get(nm)
        m['gap_frac_core'] = (sum(1 for i in present_cols
                                  if i >= len(s) or s[i] == '-') / n_present)
        # 谁在"造 block"：在 majority-gap 列里**有碱基**的序列 —— 它们才是离群者
        m['in_block_cols'] = sum(1 for c in range(ncol)
                                 if col_n[c] < 0.5 * len(seqs)
                                 and c < len(s) and s[c] != '-')
        verdict, reason = 'ok', ''
        if m['length'] < max(1, int(min_length * ref_len)):
            verdict, reason = 'too_short', f"长度 {m['length']} < {min_length:.0%}×参考({ref_len})"
        elif m['gap_frac_core'] > max_gap:
            verdict, reason = 'too_gappy', (
                '核心区 gap %.1f%% > %.0f%%（整条 gap 比例 %.1f%%；'
                '成片 gap 由集合级诊断负责）'
                % (100 * m['gap_frac_core'], 100 * max_gap, 100 * m['gap_frac']))
        elif m['n_frac'] > max_n:
            verdict, reason = 'too_ambiguous', f"简并字符 {m['n_frac']:.1%} > {max_n:.0%}"
        elif (med is not None and m['core_identity'] is not None
              and m['core_identity'] < med - max_identity_drop):
            verdict, reason = 'identity_low', (
                '核心 identity %.3f < 中位 %.3f − %s'
                % (m['core_identity'], med, max_identity_drop))
        m['verdict'], m['reason'] = verdict, reason
        rows.append(m)
        (keep if verdict == 'ok' else drop).append(nm)

    rep = {
        'ok': True, 'fasta': os.path.abspath(fasta), 'n_seq': len(seqs),
        'aligned': aligned,
        'length_min': min(lengths), 'length_max': max(lengths),
        'length_modal': modal,
        'reference_kind': ref_kind, 'reference_len': ref_len,
        'identity_median': med_ids,
        'core_identity_median': med_core,
        # 首词重名（保留为 `名.2`/`名.3`…）：条数与明细必须能对上，
        # 否则用户拿报告核 FASTA 时会发现"条数少了几条"却找不到原因
        'n_dup_names': len(dup_names), 'dup_names': dup_names[:20],
        'n_core_cols': len(core_cols),
        'n_gap_block_cols': n_big_gap_cols,
        'n_present_cols': len(present_cols),
        'block_outliers': sorted(
            ((r['name'], r['in_block_cols']) for r in rows
             if r['in_block_cols'] > 0.3 * max(1, n_big_gap_cols)),
            key=lambda t: -t[1])[:10],
        'n_keep': len(keep), 'n_drop': len(drop),
        'keep': keep, 'drop': drop,
        'drop_reasons': dict(Counter(r['verdict'] for r in rows
                                     if r['verdict'] != 'ok')),
        'per_seq': rows,
        'thresholds': {'max_gap': max_gap, 'max_n': max_n,
                       'min_length': min_length,
                       'max_identity_drop': max_identity_drop},
    }
    say('比对 QC：%d 条，%s；长度 %s~%s（众数 %s）；参考 %s；可留 %d / 剔除 %d'
        % (rep['n_seq'], '已等长' if aligned else '**不等长**',
           rep['length_min'], rep['length_max'], modal, ref_kind,
           rep['n_keep'], rep['n_drop']))
    if dup_names:
        say('比对 QC：首词重名 %d 条（已按 `名.2`/`名.3`… 保留，不会丢序列）：%s'
            % (len(dup_names),
               '、'.join('%s→%s' % (d['name'], d['kept_as'])
                         for d in dup_names[:5])))
        rep['dup_hint'] = (
            'FASTA 里有 %d 条序列的首词重名（如 %s）—— 已按 `名.2` 后缀保留，'
            '**没有丢序列**；但下游按名取序列时仍需注意这是靠后缀区分的。'
            % (len(dup_names),
               '、'.join(d['name'] for d in dup_names[:3])))
    if not aligned:
        rep['advice'] = ('序列**不等长**（A0/建树要求等长比对）→ 先跑一遍序列比对'
                         '（MAFFT），再回来 QC。')
    # ⚠️ 「大 gap 块」是**成片**差异（不是零散 indel）的特征：逐条 gap% 都差不多，
    #    且大量列只有少数序列在场。对类病毒/环状基因组尤其常见（起点不一致），
    #    也常见于"同一份导出混了不同种/株"。这条提示比"gap 太多"有用得多。
    if n_big_gap_cols > 0.05 * max(1, ncol):
        rep['block_hint'] = (
            '比对里**成片 gap**：%d/%d 列有超过一半的序列是 gap —— 这通常意味着'
            '①序列之间是**成片差异**（同科不同种/株，或环状基因组**线性化起点不一致**），'
            '而不是个别序列质量差；②逐条 gap 比例都相近，说明是**全局现象**。'
            '建议：按病毒/株拆分后再比对，或改用能处理环状基因组的流程；'
            '阈值可以放宽（max_gap），但**放宽不会让不同源的序列变成同源**。'
            % (n_big_gap_cols, ncol))
        _ol = rep['block_outliers']
        if _ol:
            rep['block_hint'] += ('　造出这些 block 的主要是 **%d 条**离群序列（如 %s）—— '
                                  '它们各自带插入/移位；其余序列是**被连累**的 '
                                  '（已改用核心区 gap 判定，不再一并判不合格）。'
                                  % (len(_ol), '、'.join(n for n, _c in _ol[:3])))
    if rep['drop_reasons']:
        rep['advice_drop'] = ('被剔除的原因分布：'
                              + '、'.join(f'{k} {v}' for k, v in
                                          rep['drop_reasons'].items()))
    return rep


def write_qc_outputs(fasta: str, report: Dict, out_dir: str,
                     *, write_sets: bool = True) -> Dict:
    """把 QC 结果落盘：`align_qc_report.tsv` +（可选）`clean.fasta` / `removed.fasta`。"""
    os.makedirs(out_dir, exist_ok=True)
    seqs = read_fasta(fasta)
    rp = os.path.join(out_dir, 'align_qc_report.tsv')
    with io.open(rp, 'w', encoding='utf-8', newline='') as f:
        f.write('name\tlength\tgap_frac\tn_frac\tidentity\tverdict\treason\n')
        for m in report['per_seq']:
            ident = '' if m['identity'] is None else f"{m['identity']:.4f}"
            f.write('%s\t%d\t%.4f\t%.4f\t%s\t%s\t%s\n'
                    % (m['name'], m['length'], m['gap_frac'], m['n_frac'],
                       ident, m['verdict'], m['reason']))
    outs = {'report': rp}
    if write_sets:
        clean = os.path.join(out_dir, 'clean.fasta')
        removed = os.path.join(out_dir, 'removed.fasta')
        keep = set(report['keep'])
        with io.open(clean, 'w', encoding='utf-8', newline='') as fc, \
                io.open(removed, 'w', encoding='utf-8', newline='') as fr:
            for nm, s in seqs.items():
                tgt = fc if nm in keep else fr
                tgt.write('>%s\n' % nm)
                for i in range(0, len(s), 60):
                    tgt.write(s[i:i + 60] + '\n')
        outs.update({'clean': clean, 'removed': removed})
    return outs


def align(fasta: str, out_fasta: str, *, threads: int = 4,
          strategy: str = 'auto', log=None) -> str:
    """调平台既有的 MAFFT 封装（含中文路径 ASCII 中转 + 体量降级）。"""
    from Virus_Platform_Core.phylo import _run_mafft
    _run_mafft(fasta, out_fasta, threads=threads, logger=None, strategy=strategy)
    if log:
        log('MAFFT 完成 → %s' % os.path.basename(out_fasta))
    return out_fasta
