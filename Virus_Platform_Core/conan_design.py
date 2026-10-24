# -*- coding: utf-8 -*-
"""CONAN-SWIFT（Cas3 / I 型 I-E）crRNA 现场检测设计（独立工具链）。

为什么独立成模块、而不是往 crispr_design.py 加一个酶档案（与 t-crispr 的区别）：
- 评分模型不同：CONAN 的核心是「负组同源位点双重失配」——负组记录在目标
  共识同源的位点窗口上，要求 PAM 已变 且 种子有效错配 才算背景干净（双双
  保留记 leak）；既有引擎的近缘筛查是 k-mer/BLAST 序列级命中统计，没有
  「同源位点」概念，硬塞进 NucleaseProfile 会把两套语义搅在一起。
- 输出形态不同：CONAN-SWIFT 是 RT-LAMP + Cas3 反式切割 + 侧流层析路线，
  只产出 crRNA spacer（向导），不做 RPA/PCR 引物对复合评分。
- 种子口径不同：I-E Cascade 种子区第 +6 位天然容忍（ignored position），
  既有引擎没有这个概念。
- 输入比对方式不同：本模块把目标组与负组合并进同一次 MAFFT，用比对列
  提取逐记录同源窗口（既有引擎比对只服务共识与保守分）。

链路：多分离物目标 FASTA（可多分体）+ 可选负组 FASTA → 片段分组 →
目标组+负组合并 MAFFT → 目标组逐列多数票共识 + 保守分 → AAG(3nt PAM)
+ 32nt spacer 双链扫描 → 逐候选在每条目标/负组记录的同源比对列上取窗 →
目标规则（PAM 保留 + 种子有效匹配 + 非种子错配 ≤ max）→ 负组规则
（PAM 已变 且 种子有效错配 = 背景干净；双双保留 = leak）→ 过滤
（距末端 ≥ margin / GC 区间 / 种子同碱基 run）→ 评分排序 → 候选表 /
位点核查表 / summary.json / 中文报告。

规则来源（净室实现：仅方法学与参数口径对齐，代码未复制。参考
kimihira-tetsutaro-code/conan-swift-designer 与 conan-crrna-webui，MIT
许可；方法学出处 bioRxiv 2026.08.25.746999 CONAN-SWIFT 论文）：
- PAM = AAG（5′ 前置，可配置其他 I 型 PAM）、spacer 默认 32nt、
  种子 = PAM+1 起 8nt、忽略位 +6（I-E Cascade 容忍位）
- 目标非种子错配上限默认 4；末端排除默认 150bp；GC 默认 35–65%；
  种子同碱基 run 默认 ≥3 剔除
- 评分 = 目标 PAM 保守率 + 目标种子匹配率 − 负组 PAM 残留率
  − 负组种子匹配率（无负组时负项记 0）

证据边界（与平台 CRISPR 检测设计同一口径，写入每份报告）：计算设计仅
支持候选优先级排序，不构成 crRNA 活性、错配区分度、扩增兼容性或基质
性能的证据；实验结论须按平台 L1–L5 验证清单逐级建立。I 型 crRNA 直接
重复（DR 骨架）属实验构造件，不在计算输出范围——订购前按所用系统的
Protocol 补齐 DR-spacer 完整序列。

工程约定（对齐 crispr_design）：全程离线，仅复用已纳管工具 MAFFT
（phylo._run_mafft，内部自带中文路径 ASCII 中转）；进程可取消（cancel
置位即抛 RuntimeError）；中文输出统一 utf-8 落盘。
"""
import json
import os
import re
import time

from .crispr_design import (CrisprDesignError, _read_aln, _write_tmp_fasta,
                            downsample_profile, gc_frac, group_segments,
                            max_homo_run, read_targets, revcomp)
from .phylo import _run_mafft
from .utils import safe_open

ENGINE_VERSION = '1.0.0'

# ---- CONAN-SWIFT 口径常量（均可被 params 覆盖，钳制见 run_chain）----
PAM_DEFAULT = 'AAG'                    # I-E Cascade PAM（5′ 前置）
SPACER_MIN, SPACER_MAX = 24, 36       # CONAN-SWIFT 论文用 32
SPACER_DEFAULT = 32
SEED_LEN_MIN, SEED_LEN_MAX = 4, 12
SEED_LEN_DEFAULT = 8
IGNORED_DEFAULT = '6'                  # 种子忽略位（1 基，逗号分隔）
MAX_NONSEED_MM_DEFAULT = 4
TERMINAL_MARGIN_DEFAULT = 150
GC_RANGE_DEFAULT = (35.0, 65.0)
HOMOPOLYMER_DEFAULT = 3                # 种子同碱基 run 下限

_PURE_ACGT = set('ACGT')
_RC = str.maketrans('ACGTURYSWKMBDHVN', 'TGCAAYRSWMKVHDBN')
_IUPAC_PAIRS = {'AT': 'W', 'CG': 'S', 'AC': 'M', 'GT': 'K',
                'AG': 'R', 'CT': 'Y'}


class _Log:
    """日志适配：log 回调可能为 None；异常吞掉不阻塞主流程。

    同时暴露 .log()（phylo._run_mafft 等平台工具按 TaskLogger 口径调用）。
    """

    def __init__(self, callback):
        self._cb = callback

    def _emit(self, msg, level='INFO'):
        if self._cb is None:
            return
        try:
            self._cb(msg if level == 'INFO' else '[%s] %s' % (level, msg))
        except Exception:
            pass

    def __call__(self, msg, level='INFO'):
        self._emit(msg, level)

    def log(self, msg, level='INFO'):
        self._emit(msg, level)


# ---------------------------------------------------------------
# 小工具
# ---------------------------------------------------------------
def parse_ignored(text):
    """'6' / '6, 9' → {6, 9}（1 基种子位置集合）；空/非法输入 → 空集。"""
    out = set()
    for tok in re.split(r'[,\s;]+', str(text or '')):
        if tok.isdigit() and int(tok) > 0:
            out.add(int(tok))
    return out


def _ungapped(seq):
    return seq.replace('-', '').replace('.', '')


def _pure_acgt(seq):
    return all(c in _PURE_ACGT for c in seq)


def _mismatch_positions(ref, obs, ignored):
    """1 基错配位置列表；ignored = 忽略的 1 基位置集合。"""
    return [i + 1 for i, (a, b) in enumerate(zip(ref, obs))
            if a != b and (i + 1) not in ignored]


# ---------------------------------------------------------------
# 比对与共识（目标组 + 负组合并进同一次 MAFFT，保证同源列对齐）
# ---------------------------------------------------------------
def _align_combined(target_recs, negative_recs, work_dir, threads=None,
                    logger=None, cancel=None):
    """合并 MAFFT。返回 ({'t': [(name, aln)], 'n': [(name, aln)]}, msa_used)。

    单目标且无负组：直通不比对。名字前缀 t0/n0 防两库重名撞车。
    """
    if not negative_recs and len(target_recs) == 1:
        return ({'t': [(target_recs[0][0], target_recs[0][1].upper())],
                 'n': []}, False)
    combined = ([('t%d' % i, s) for i, (_, s) in enumerate(target_recs)]
                + [('n%d' % i, s) for i, (_, s) in enumerate(negative_recs)])
    in_fa = os.path.join(work_dir, 'conan_in.fa')
    aln_fa = os.path.join(work_dir, 'conan_aln.fa')
    _write_tmp_fasta(combined, in_fa)
    if cancel is not None and cancel.is_set():
        raise RuntimeError('用户取消')
    _run_mafft(in_fa, aln_fa, threads=threads, logger=logger)
    if cancel is not None and cancel.is_set():
        raise RuntimeError('用户取消')
    aln = _read_aln(aln_fa, [n for n, _ in combined])
    if not aln:
        raise CrisprDesignError('MAFFT 比对结果为空')
    by_name = dict(aln)
    t_out = [(name, by_name.get('t%d' % i, '')) for i, (name, _) in
             enumerate(target_recs)]
    n_out = [(name, by_name.get('n%d' % i, '')) for i, (name, _) in
             enumerate(negative_recs)]
    return ({'t': [(n, s) for n, s in t_out if s],
             'n': [(n, s) for n, s in n_out if s]}, True)


def _consensus_from_aln(target_aln):
    """目标组比对行 → (共识串, 逐列保守分, 列映射 col_map)。

    与 crispr_design.consensus_and_conservation 同一规则：多数票定碱基，
    两碱基平票记 IUPAC 兼并码（W/S/M/K/R/Y），≥3 分歧记 N；保守分 =
    多数票碱基占该列非 gap 碱基的比例。

    关键差异（负组同框比对的必然要求）：**全目标皆 gap 的列直接跳过**
    （多为负组特有插入），不产生 N 位——否则负组插入列会把目标保守窗
    打断成 N，扫描全军覆没（合成数据实测教训）。共识第 k 位 ↔ 合并
    比对第 col_map[k] 列；目标/负组记录的同源窗口都按 col_map 提取。
    """
    n_cols = max(len(s) for _, s in target_aln)
    rows = [s for _, s in target_aln]
    cons, cols, col_map = [], [], []
    for k in range(n_cols):
        bases = [r[k] for r in rows if k < len(r) and r[k] in 'ACGT']
        if not bases:
            continue                     # 全目标 gap（负组插入列）→ 跳过
        counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
        for b in bases:
            counts[b] += 1
        best, cnt = max(counts.items(), key=lambda kv: (kv[1], kv[0]))
        if cnt * 2 > len(bases):
            cons.append(best)
        else:
            tied = sorted(b for b, c in counts.items() if c == cnt)
            cons.append(_IUPAC_PAIRS.get(''.join(tied), 'N')
                        if len(tied) == 2 else 'N')
        cols.append(cnt / float(len(bases)))
        col_map.append(k)
    return ''.join(cons), cols, col_map


# ---------------------------------------------------------------
# 候选扫描与逐记录同源窗核查
# ---------------------------------------------------------------
def scan_candidates(consensus, pam, spacer_len, seed_len):
    """双链扫描 PAM+spacer 候选（共识串上，0 基半开区间 span）。

    返回 dict 列表：strand / span(共识 0 基 [start, end)) / pam /
    protospacer(搜索取向 PAM 后端) / guide(向导 = protospacer 反义链) /
    seed(种子 = protospacer 前 seed_len 位) / window(搜索取向全窗)。
    去重键 = (strand, window)，与参考实现一致。
    """
    query_len = len(pam) + spacer_len
    L = len(consensus)
    out = []
    seen = set()
    for strand, seq in (('+', consensus), ('-', revcomp(consensus))):
        for j in range(0, max(0, L - query_len + 1)):
            window = seq[j:j + query_len]
            if window[:len(pam)] != pam or not _pure_acgt(window):
                continue
            if (strand, window) in seen:
                continue
            seen.add((strand, window))
            spacer = window[len(pam):]
            start = j if strand == '+' else L - (j + query_len)
            out.append({
                'strand': strand, 'span': (start, start + query_len),
                'pam': pam, 'protospacer': spacer, 'guide': revcomp(spacer),
                'seed': spacer[:seed_len], 'window': window,
            })
    return out


def _homologous_window(aln_seq, span, strand, col_map):
    """按列映射从记录比对行取候选同源窗，翻到搜索取向。

    span 是共识（目标列空间）0 基半开区间；记录在该列无碱基（gap）
    时贡献空串 → 窗口变短 → 上游 clean 判 False。
    """
    start, end = span
    cols = col_map[start:end]
    w = _ungapped(''.join(aln_seq[c] for c in cols if c < len(aln_seq)))
    return revcomp(w) if strand == '-' else w


def evaluate_candidate(cand, target_aln, negative_aln, seed_len, ignored,
                       max_nonseed_mm, col_map):
    """单候选 × 全部目标/负组记录的同源位点核查。

    目标规则：clean（等长纯 ACGT）且 PAM 匹配 且 种子有效匹配 且
    非种子有效错配 ≤ max_nonseed_mm。
    负组规则：clean 且 PAM 已变 且 种子有效错配 → 该记录背景干净；
    PAM 残留且种子匹配 = leak。
    返回 (row, target_rows, negative_rows)。
    """
    seed = cand['seed']
    spacer = cand['protospacer']
    strand, span = cand['strand'], cand['span']

    target_rows, negative_rows = [], []
    t_pam = t_seed = t_ok = 0
    n_pam_res = n_seed_match = n_leak = n_pass = 0
    t_nonseed_total = 0

    for name, aln_seq in target_aln:
        w = _homologous_window(aln_seq, span, strand, col_map)
        obs_pam = w[:len(cand['pam'])]
        obs_spacer = w[len(cand['pam']):]
        obs_seed = obs_spacer[:seed_len]
        clean = len(obs_spacer) == len(spacer) and _pure_acgt(w)
        seed_mm = _mismatch_positions(seed, obs_seed, ignored)
        spacer_mm = _mismatch_positions(spacer, obs_spacer, ignored)
        nonseed_mm = [p for p in spacer_mm if p > seed_len]
        pam_match = clean and obs_pam == cand['pam']
        seed_match = clean and not seed_mm
        ok = (clean and obs_pam == cand['pam'] and not seed_mm
              and len(nonseed_mm) <= max_nonseed_mm)
        t_pam += int(pam_match)
        t_seed += int(seed_match)
        t_ok += int(ok)
        t_nonseed_total += len(nonseed_mm)
        target_rows.append({
            'candidate_id': cand['id'], 'target_record': name,
            'observed_PAM': obs_pam, 'PAM_match': pam_match,
            'observed_seed': obs_seed, 'seed_match': seed_match,
            'effective_seed_mm': ','.join(map(str, seed_mm)),
            'effective_spacer_mm': ','.join(map(str, spacer_mm)),
            'effective_nonseed_mm_count': len(nonseed_mm),
            'target_pass': ok, 'observed_window': w,
        })

    for name, aln_seq in negative_aln:
        w = _homologous_window(aln_seq, span, strand, col_map)
        obs_pam = w[:len(cand['pam'])]
        obs_seed = w[len(cand['pam']):len(cand['pam']) + seed_len]
        clean = (len(obs_pam) == len(cand['pam'])
                 and len(obs_seed) == seed_len and _pure_acgt(obs_pam + obs_seed))
        seed_mm = _mismatch_positions(seed, obs_seed, ignored)
        pam_residual = clean and obs_pam == cand['pam']
        seed_match = clean and not seed_mm
        leak = pam_residual and seed_match
        n_pam_res += int(pam_residual)
        n_seed_match += int(seed_match)
        n_leak += int(leak)
        n_pass += int(clean and not pam_residual and not seed_match)
        negative_rows.append({
            'candidate_id': cand['id'], 'negative_record': name,
            'observed_PAM': obs_pam,
            'PAM_residual': pam_residual,
            'PAM_changed': clean and not pam_residual,
            'observed_seed': obs_seed, 'seed_match': seed_match,
            'effective_seed_mm_count': len(seed_mm),
            'effective_seed_mm': ','.join(map(str, seed_mm)),
            'negative_locus_pass': clean and not pam_residual
            and not seed_match,
            'observed_window': w,
        })

    n_t = max(1, len(target_aln))
    has_neg = bool(negative_aln)
    n_n = max(1, len(negative_aln))
    # 未运行负组评分的字段一律 None（JSON null / TSV 空串），不再用
    # 'not_run' 字符串混进数值列（评审 R8）
    row = {
        'candidate_id': cand['id'], 'strand': strand, 'PAM': cand['pam'],
        'cons_start': span[0] + 1, 'cons_end': span[1],
        'protospacer_5to3_after_PAM': spacer,
        'guide_spacer_5to3_RC': cand['guide'],
        'locus_sequence_5to3': cand['window'],
        'seed': seed,
        'target_PAM_conservation_rate': round(100 * t_pam / n_t, 2),
        'target_seed_match_rate': round(100 * t_seed / n_t, 2),
        'target_coverage_rate': round(100 * t_ok / n_t, 2),
        'negative_PAM_residual_rate':
            round(100 * n_pam_res / n_n, 2) if has_neg else None,
        'negative_seed_match_rate':
            round(100 * n_seed_match / n_n, 2) if has_neg else None,
        'negative_leak_rate':
            round(100 * n_leak / n_n, 2) if has_neg else None,
        'target_pass_count': t_ok, 'target_total_count': len(target_aln),
        'negative_locus_pass_count': n_pass if has_neg else None,
        'negative_total_count': len(negative_aln),
        'all_targets_pass': t_ok == len(target_aln),
        'all_negatives_locus_pass':
            (n_pass == len(negative_aln)) if has_neg else None,
        'target_nonseed_mm_total': t_nonseed_total,
    }
    if has_neg:
        row['score'] = round(
            row['target_PAM_conservation_rate']
            + row['target_seed_match_rate']
            - row['negative_PAM_residual_rate']
            - row['negative_seed_match_rate'], 2)
    else:
        row['score'] = round(
            row['target_PAM_conservation_rate']
            + row['target_seed_match_rate'], 2)
    return row, target_rows, negative_rows


def apply_filters(cand_row, cons_len, params):
    """末端 margin / GC 区间 / 种子同碱基 run 过滤。返回 (通过?, 原因)。"""
    reasons = []
    margin = params['terminal_margin']
    if params['exclude_terminal'] and margin > 0:
        if (cand_row['cons_start'] <= margin
                or cand_row['cons_end'] > cons_len - margin):
            reasons.append('距共识末端<%dbp' % margin)
    gc = round(100 * gc_frac(cand_row['protospacer_5to3_after_PAM']), 2)
    cand_row['spacer_GC_percent'] = gc
    if params['select_gc'] and not (
            params['gc_min'] <= gc <= params['gc_max']):
        reasons.append('GC%% 超出 %g-%g' % (params['gc_min'], params['gc_max']))
    homo = max_homo_run(cand_row['seed'])
    cand_row['seed_max_homo_run'] = homo
    if params['exclude_homopolymer'] and homo >= params['homopolymer']:
        reasons.append('种子含≥%d同碱基run' % params['homopolymer'])
    return (not reasons), '; '.join(reasons)


def _sort_key(row):
    """排序：覆盖↓ / leak↑ / PAM 保守↓ / 种子匹配↓ / score↓ / 非种子错配↑。"""
    leak = row['negative_leak_rate'] if isinstance(
        row['negative_leak_rate'], (int, float)) else 0.0
    neg_pass = row['all_negatives_locus_pass']
    return (-row['target_coverage_rate'], 0 if leak == 0 else 1, leak,
            -row['target_PAM_conservation_rate'],
            -row['target_seed_match_rate'], -row['score'],
            row['target_nonseed_mm_total'],
            0 if neg_pass is True else 1)


# ---------------------------------------------------------------
# 输出落盘
# ---------------------------------------------------------------
_CAND_COLS = [
    'segment', 'rank', 'candidate_id', 'strand', 'PAM', 'cons_start',
    'cons_end', 'locus_sequence_5to3', 'protospacer_5to3_after_PAM',
    'guide_spacer_5to3_RC', 'seed', 'spacer_GC_percent',
    'seed_max_homo_run', 'cons_mean',
    'filter_pass', 'filter_reason', 'target_coverage_rate',
    'target_PAM_conservation_rate', 'target_seed_match_rate',
    'negative_PAM_residual_rate', 'negative_seed_match_rate',
    'negative_leak_rate', 'target_pass_count', 'target_total_count',
    'negative_locus_pass_count', 'negative_total_count',
    'all_targets_pass', 'all_negatives_locus_pass', 'score']


def _write_tsv(path, cols, rows, yes_no=()):
    with safe_open(path, 'wt') as f:
        f.write('\t'.join(cols) + '\n')
        for r in rows:
            vals = []
            for c in cols:
                v = r.get(c, '')
                if c in yes_no and isinstance(v, bool):
                    v = '是' if v else '否'
                vals.append('' if v is None else str(v))
            f.write('\t'.join(vals) + '\n')


def _fmt_rate(v):
    if isinstance(v, (int, float)):
        return '%.2f' % v
    return '' if v is None else str(v)


def write_report(path, summary, params, top_rows):
    """中文 Markdown 报告（含参数快照与证据边界）。"""
    seg_lines = ['| 片段组 | 目标记录 | 负组记录 | 候选(过滤前) | 过滤后 |',
                 '|---|---|---|---|---|']
    for s in summary['segments']:
        if s['n_negatives'] is None:
            neg_disp = '未提供'
        elif s.get('negative_applied') is False:
            neg_disp = '未应用（无同名片段标签，跨片段比对无同源语义）'
        else:
            neg_disp = str(s['n_negatives'])
        seg_lines.append('| %s | %d | %s | %d | %d |' % (
            s['segment'], s['n_targets'], neg_disp, s['n_raw'], s['n_kept']))
    top_lines = []
    if top_rows:
        top_lines = [
            '| # | 候选 | 链 | 位点 | spacer(向导=RC列) | GC%% | 覆盖%% | leak%% | score |',
            '|---|---|---|---|---|---|---|---|---|']
        for i, r in enumerate(top_rows, 1):
            top_lines.append('| %d | %s | %s | %d-%d | `%s` | %s | %s | %s | %s |' % (
                i, r['candidate_id'], r['strand'], r['cons_start'],
                r['cons_end'], r['guide_spacer_5to3_RC'],
                r['spacer_GC_percent'], _fmt_rate(r['target_coverage_rate']),
                _fmt_rate(r['negative_leak_rate']), _fmt_rate(r['score'])))
    with safe_open(path, 'wt') as f:
        f.write('''# CONAN-SWIFT（Cas3 / I 型 I-E）crRNA 设计报告

## 结果概览

%s

## Top 候选（全表见 candidates.tsv）

%s

## 参数快照

- PAM：`%s`（5′ 前置）；spacer：%d nt；种子：PAM+1 起 %d nt（忽略位 %s）
- 目标非种子错配上限：%d；末端排除：%s bp；GC 区间：%s–%s%%；
  种子同碱基 run ≥%d 剔除
- 负组 FASTA：%s
- 引擎：VirusPlatform conan_design %s

## 方法学与规则来源

净室实现（代码未复制，仅方法学与参数口径对齐）：
- CONAN-SWIFT 论文（bioRxiv 2026.08.25.746999）：AAG PAM / 32nt spacer /
  种子 8nt 忽略 +6 / 负组同源位点双重失配规则 / 四项评分公式
- 参考开源 conan-swift-designer / conan-crrna-webui（MIT）：
  过滤参数（末端 150bp / GC 35–65 / 种子 run≥3）

**评分使用边界**：score 为率差线性组合，数值依赖本次目标/负组条数，
仅限单次输入内的候选排序比较，跨数据集不可比。

**负组同源性要求**：负组应为近缘非目标病毒的全基因组序列（与目标整体
同源 ≥70–80%% 为宜）。位点块嵌于无关序列时，比对列配准可能漂移，
负组核查结果不可信。

## 证据边界

计算设计仅支持候选优先级排序，不构成 crRNA 活性、错配区分度、扩增兼容性
或基质性能的证据；实验结论须按平台 L1–L5 验证清单逐级建立。

**订购提示**：本表只输出 spacer 向导序列（guide_spacer_5to3_RC 列）。
I 型 CRISPR crRNA 需 5′ 端接直接重复（DR）骨架构成完整 crRNA，DR 序列
属所用 Cas3 系统的实验构造件——订购前按所用系统的 Protocol 补齐
DR-spacer 完整序列并复核方向。
''' % ('\n'.join(seg_lines), '\n'.join(top_lines) if top_lines
            else '（无候选通过过滤——检查 PAM/参数或负组设置）',
            params['pam'], params['spacer_len'], params['seed_len'],
            params['ignored_text'] or '无', params['max_nonseed_mm'],
            params['terminal_margin'] if params['exclude_terminal'] else '不排除',
            params['gc_min'], params['gc_max'], params['homopolymer'],
            params['negative'] or '未提供（负组评分按 0 计）',
            ENGINE_VERSION))


# ---------------------------------------------------------------
# 结果图（matplotlib 直出；不可用时跳过不阻塞设计）
# ---------------------------------------------------------------
def _write_plots(out_root, segments_out, cons_store, logger):
    """结果可视化：Top 候选柱状图 + 各片段保守性剖面位点总览图。

    字体约定对齐 dfv_plot：临时切 Microsoft YaHei/SimHei，画完还原。
    返回落盘文件 dict（键与 summary['files'] 对齐）。
    """
    def _logw(msg):
        if logger:
            logger(msg, 'WARN')

    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except Exception:
        _logw('matplotlib 不可用，跳过结果图（候选表/报告不受影响）')
        return {}
    old = (plt.rcParams['font.sans-serif'], plt.rcParams['axes.unicode_minus'])
    plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei',
                                       'DejaVu Sans']
    plt.rcParams['axes.unicode_minus'] = False
    files = {}
    try:
        def _num(v):
            return v if isinstance(v, (int, float)) else 0.0

        rows = [r for seg in sorted(segments_out)
                for r in segments_out[seg]['candidates']][:15]
        if rows:
            fig, ax = plt.subplots(
                figsize=(max(8.0, 0.7 * len(rows)), 4.4))
            x = list(range(len(rows)))
            w = 0.38
            ax.bar([i - w / 2 for i in x],
                   [_num(r['target_coverage_rate']) for r in rows], w,
                   label='目标覆盖率%', color='#2e7d32')
            ax.bar([i + w / 2 for i in x],
                   [_num(r['negative_leak_rate']) for r in rows], w,
                   label='负组leak%', color='#c62828')
            ax2 = ax.twinx()
            ax2.plot(x, [_num(r['score']) for r in rows], 'o-',
                     color='#1565c0', label='score')
            ax2.set_ylabel('score', color='#1565c0')
            ax.set_xticks(x)
            ax.set_xticklabels(
                [r['candidate_id'].replace('CN', '') for r in rows],
                rotation=45, ha='right', fontsize=8)
            ax.set_ylabel('率 (%)')
            ax.set_ylim(0, 105)
            h1, l1 = ax.get_legend_handles_labels()
            h2, l2 = ax2.get_legend_handles_labels()
            ax.legend(h1 + h2, l1 + l2, loc='lower left', fontsize=8)
            ax.set_title('Top 候选：目标覆盖率 / 负组 leak / score')
            fig.tight_layout()
            p = os.path.join(out_root, 'top_candidates.png')
            fig.savefig(p, dpi=150)
            for ext in ('svg', 'pdf'):
                try:
                    fig.savefig(p[:-4] + '.' + ext, bbox_inches='tight')
                except Exception:
                    pass      # 附加矢量格式失败不影响 PNG 主产物
            plt.close(fig)
            files['top_plot'] = p

        segs = [seg for seg in sorted(segments_out) if seg in cons_store]
        if segs:
            fig, axes = plt.subplots(
                len(segs), 1, squeeze=False,
                figsize=(10.0, 2.6 * len(segs)))
            for ax, seg in zip(axes[:, 0], segs):
                cons, L = cons_store[seg]
                prof = downsample_profile(cons, 2000)
                xs = [i * L / float(len(prof)) for i in range(len(prof))]
                ax.plot(xs, prof, lw=0.8, color='#455a64')
                ax.set_ylim(0, 1.05)
                ax.set_ylabel('保守性', fontsize=8)
                ax.set_title('%s（共识 %dbp，橙色带 = Top5 候选位点）'
                             % (seg, L), fontsize=9)
                for r in segments_out[seg]['candidates'][:5]:
                    ax.axvspan(r['cons_start'] - 1, r['cons_end'],
                               color='#ff8f00', alpha=0.25)
            axes[-1][0].set_xlabel('共识坐标 (bp)', fontsize=8)
            fig.tight_layout()
            p = os.path.join(out_root, 'locus_overview.png')
            fig.savefig(p, dpi=150)
            for ext in ('svg', 'pdf'):
                try:
                    fig.savefig(p[:-4] + '.' + ext, bbox_inches='tight')
                except Exception:
                    pass      # 附加矢量格式失败不影响 PNG 主产物
            plt.close(fig)
            files['locus_plot'] = p
    except Exception as e:
        _logw('结果图绘制失败（不影响候选表/报告）: %s' % e)
    finally:
        plt.rcParams['font.sans-serif'], plt.rcParams['axes.unicode_minus'] = old
    return files


# ---------------------------------------------------------------
# 全链路入口
# ---------------------------------------------------------------
def run_chain(targets, params, out_dir, log=None, prog=None, cancel=None):
    """CONAN-SWIFT 设计全链路入口。

    targets : 目标 FASTA 路径（多分离物，可多分体，header 可含 RNA1 标签）
    params  : {'negative', 'pam', 'spacer_len', 'seed_len', 'ignored',
               'max_nonseed_mm', 'terminal_margin', 'exclude_terminal',
               'select_gc', 'gc_min', 'gc_max', 'exclude_homopolymer',
               'homopolymer', 'group_by_segment', 'top_n', 'threads'}
    out_dir : 产物根目录（写 out_dir/conan/）
    返回 (summary, report_dict)。
    """
    _log = _Log(log)

    def _prog(stage, frac, msg=''):
        if prog is not None:
            prog(stage, frac, msg)

    def _int(key, default, lo, hi):
        # 仅 None/空串取默认；数值 0 是合法参数（max_nonseed_mm=0 严格
        # 档 / terminal_margin=0 关闭末端排除），不得吞成默认值
        # （评审 R1：实测 0→默认 的静默漂移）
        raw = params.get(key)
        if raw is None or (isinstance(raw, str) and not raw.strip()):
            return default
        try:
            v = int(str(raw).strip())
        except (TypeError, ValueError):
            return default
        return max(lo, min(hi, v))

    def _flag(key, default=True):
        v = params.get(key, '1' if default else '0')
        return str(v).strip().lower() not in ('0', 'false', 'no')

    t0 = time.time()
    pam = re.sub(r'[^ACGT]', '', str(params.get('pam') or PAM_DEFAULT)
                 .upper())
    if not pam:
        raise CrisprDesignError('PAM 必须为 ACGT 序列（CONAN-SWIFT 默认 AAG）')
    spacer_len = _int('spacer_len', SPACER_DEFAULT, SPACER_MIN, SPACER_MAX)
    seed_len = _int('seed_len', SEED_LEN_DEFAULT, SEED_LEN_MIN, SEED_LEN_MAX)
    ignored_text = str(params.get('ignored') or IGNORED_DEFAULT).strip()
    ignored = parse_ignored(ignored_text)
    ignored = {p for p in ignored if p <= spacer_len}
    max_nonseed_mm = _int('max_nonseed_mm', MAX_NONSEED_MM_DEFAULT, 0,
                          spacer_len)
    def _gc(key, default):
        # 同 _int 口径：数值 0（GC 下限 0）不得被 or 吞成默认（评审 R1）
        raw = params.get(key)
        if raw is None or (isinstance(raw, str) and not raw.strip()):
            return default
        try:
            return max(0.0, min(100.0, float(raw)))
        except (TypeError, ValueError):
            return default

    gc_min = _gc('gc_min', GC_RANGE_DEFAULT[0])
    gc_max = _gc('gc_max', GC_RANGE_DEFAULT[1])
    if gc_min > gc_max:
        gc_min, gc_max = gc_max, gc_min
    over_seed = sorted(p for p in parse_ignored(ignored_text) if p > seed_len)
    if over_seed:
        # 忽略位是种子位置（1 基）；超出部分对种子核查无效果，只对
        # spacer 错配计数生效（评审 R9：此前静默截断）
        _log('忽略位 %s 超出种子长度 %d：对种子核查无效，仅对 spacer '
             '错配计数生效' % (over_seed, seed_len), 'WARN')
    cfg = {
        'pam': pam, 'spacer_len': spacer_len, 'seed_len': seed_len,
        'ignored': ignored, 'ignored_text': ignored_text,
        'max_nonseed_mm': max_nonseed_mm,
        'terminal_margin': _int('terminal_margin', TERMINAL_MARGIN_DEFAULT,
                                0, 2000),
        'exclude_terminal': _flag('exclude_terminal', True),
        'select_gc': _flag('select_gc', True),
        'gc_min': gc_min, 'gc_max': gc_max,
        'exclude_homopolymer': _flag('exclude_homopolymer', True),
        'homopolymer': _int('homopolymer', HOMOPOLYMER_DEFAULT, 2, seed_len),
        'group_by_segment': _flag('group_by_segment', True),
        'top_n': _int('top_n', 20, 1, 200),
        'negative': str(params.get('negative') or '').strip() or None,
        'threads': int(params.get('threads') or 0) or None,
    }

    out_root = os.path.join(out_dir, 'conan')
    work_dir = os.path.join(out_root, '_work')
    os.makedirs(work_dir, exist_ok=True)

    # ---------- ① 输入 ----------
    _prog('读取', 0.03, '读取目标 / 负组 FASTA')
    records = read_targets(targets)
    groups = group_segments(records, cfg['group_by_segment'])
    negative_alns_all = None
    if cfg['negative']:
        try:
            neg_recs = read_targets(cfg['negative'])
        except CrisprDesignError as e:
            # read_targets 文案写死「目标」；负组入口改写（评审 R5）
            raise CrisprDesignError(str(e).replace('目标', '负组'))
        neg_groups = group_segments(neg_recs, cfg['group_by_segment'])
        _log('负组 %d 条记录（%d 个片段组）' % (len(neg_recs), len(neg_groups)))
    else:
        neg_recs, neg_groups = [], {}
    _log('输入目标 %d 条记录 → %d 个片段组: %s' % (
        len(records), len(groups), ', '.join(
            '%s(%d)' % (k, len(v)) for k, v in groups.items())))

    # ---------- ② 逐片段组：合并比对 → 扫描 → 核查 → 过滤 → 排序 ----------
    segments_out = {}
    seg_meta = []
    all_t_rows, all_n_rows = [], []
    cons_store = {}                          # 保守性剖面（结果图用）
    seg_names = sorted(groups)
    n_prog = 0.02
    for seg in seg_names:
        t_recs = groups[seg]
        negative_applied = False
        if neg_groups:
            if seg in neg_groups:
                seg_negs = neg_groups[seg]
                negative_applied = True
            elif len(groups) == 1 and len(neg_groups) == 1:
                # 双方都是单片段组：无跨片段错配语义问题，回退全量
                seg_negs = neg_recs
                negative_applied = True
                _log('片段组 %s 在负组 FASTA 中无同名标签（双方均为单一'
                     '片段组）→ 负组全量 %d 条参与该组比对'
                     % (seg, len(seg_negs)), 'WARN')
            else:
                # 多片段组下跨片段「同源位点」实为随机序列比对，leak 率
                # 无意义 → 跳过该组负组评分并在报告标注（评审 R2）
                seg_negs = []
                _log('片段组 %s 在负组 FASTA 中无同名片段标签且任一侧为'
                     '多片段组 → 该组跳过负组评分（跨片段比对无同源语义）'
                     % seg, 'WARN')
        else:
            seg_negs = []
        _prog('比对', n_prog, '片段 %s：目标+负组合并 MAFFT（%d+%d 条）'
              % (seg, len(t_recs), len(seg_negs)))
        aln, msa_used = _align_combined(t_recs, seg_negs, work_dir,
                                        threads=cfg['threads'], logger=_log,
                                        cancel=cancel)
        consensus, cons, col_map = _consensus_from_aln(aln['t'])
        cons_store[seg] = (cons, len(consensus))
        _log('片段 %s: 目标 %d 条 → 共识 %d 列（MSA=%s；负组 %d 条同框比对）'
             % (seg, len(aln['t']), len(consensus),
                '是' if msa_used else '否', len(aln['n'])))

        _prog('扫描', n_prog + 0.06, '片段 %s：%s PAM+%dnt 双链扫描'
              % (seg, cfg['pam'], cfg['spacer_len']))
        cands = scan_candidates(consensus, cfg['pam'], cfg['spacer_len'],
                                cfg['seed_len'])
        cand_rows, t_rows, n_rows = [], [], []
        n_ambig = sum(1 for ch in consensus if ch not in 'ACGT')
        if n_ambig:
            _log('片段 %s: 共识含 %d 个非 ACGT 列（平票 IUPAC / 分歧 N），'
                 '落在其上的候选窗会被跳过' % (seg, n_ambig), 'WARN')
        for i, c in enumerate(cands):
            c['id'] = 'CN%s_%05d' % (seg, i + 1)
            row, tr, nr = evaluate_candidate(
                c, aln['t'], aln['n'], cfg['seed_len'], cfg['ignored'],
                cfg['max_nonseed_mm'], col_map)
            c0, c1 = c['span']
            seg_cons = cons[c0:c1]
            row['cons_mean'] = round(sum(seg_cons) / len(seg_cons), 4) \
                if seg_cons else 0.0
            ok, reason = apply_filters(row, len(consensus), cfg)
            row['segment'] = seg
            row['filter_pass'] = ok
            row['filter_reason'] = reason
            cand_rows.append(row)
            if ok:                     # 位点核查明细只落过滤后候选（控表体量）
                t_rows.extend(tr)
                n_rows.extend(nr)
        all_t_rows.extend(dict(r, segment=seg) for r in t_rows)
        all_n_rows.extend(dict(r, segment=seg) for r in n_rows)
        kept = [r for r in cand_rows if r['filter_pass']]
        kept.sort(key=_sort_key)
        for rank, r in enumerate(kept, 1):
            r['rank'] = rank
        segments_out[seg] = {
            'n_targets': len(aln['t']), 'n_negatives': len(aln['n']),
            'consensus_len': len(consensus), 'n_raw': len(cands),
            'n_kept': len(kept), 'all_targets_pass': sum(
                1 for r in kept if r['all_targets_pass']),
            'background_clear': sum(
                1 for r in kept
                if r['all_negatives_locus_pass'] is True),
            'candidates': kept[:cfg['top_n']],
        }
        seg_meta.append({
            'segment': seg, 'n_targets': len(aln['t']),
            'n_negatives': len(aln['n']) if cfg['negative'] else None,
            'negative_applied': (negative_applied
                                 if cfg['negative'] else None),
            'n_raw': len(cands), 'n_kept': len(kept),
            'all_targets_pass': segments_out[seg]['all_targets_pass'],
            'background_clear': segments_out[seg]['background_clear'],
        })
        n_prog = min(0.88, n_prog + 0.8 / max(1, len(seg_names)))

    # ---------- ③ 落盘 ----------
    _prog('输出', 0.95, '写候选表 / 位点核查 / 结果图 / 报告')
    plot_files = _write_plots(out_root, segments_out, cons_store, _log)
    all_kept = [dict(r, segment=seg) for seg, d in
                sorted(segments_out.items())
                for r in d['candidates']]
    # 位点核查矩阵：逐片段组 候选(Top) × 记录 → 状态格（ok/mmN/pam/seed/
    # leak/gap）。CONAN-SWIFT 图1e 那格是湿实验交叉反应矩阵；这里是同一
    # 布局的计算侧类比，数据源即 target/negative_checks
    top_map = {r['candidate_id']: r for r in all_kept[:cfg['top_n']]}
    matrices = []
    for seg in seg_names:
        t_names = sorted({r['target_record'] for r in all_t_rows
                          if r['segment'] == seg})
        n_names = sorted({r['negative_record'] for r in all_n_rows
                          if r['segment'] == seg})
        rows = []
        for cid, crow in top_map.items():
            if crow['segment'] != seg:
                continue
            cells = {}
            for r in all_t_rows:
                if r['segment'] != seg or r['candidate_id'] != cid:
                    continue
                if r['target_pass']:
                    st = 'ok'
                elif not r['PAM_match']:
                    st = 'pam'
                elif not r['seed_match']:
                    st = 'seed'
                elif r['effective_nonseed_mm_count']:
                    st = 'mm%d' % r['effective_nonseed_mm_count']
                else:
                    st = 'gap'
                cells[r['target_record']] = st
            for r in all_n_rows:
                if r['segment'] != seg or r['candidate_id'] != cid:
                    continue
                if r['PAM_residual'] and r['seed_match']:
                    st = 'leak'
                elif r['PAM_residual']:
                    st = 'pam'
                elif r['seed_match']:
                    st = 'seed'
                elif r['negative_locus_pass']:
                    st = 'ok'
                else:
                    st = 'gap'
                cells[r['negative_record']] = st
            rows.append({'candidate_id': cid, 'rank': crow.get('rank'),
                         'cells': cells})
        if rows:
            matrices.append({
                'segment': seg,
                'cols': ([{'kind': 'T', 'name': n} for n in t_names]
                         + [{'kind': 'N', 'name': n} for n in n_names]),
                'rows': rows,
            })
    _write_tsv(os.path.join(out_root, 'candidates.tsv'), _CAND_COLS, all_kept,
               yes_no=('filter_pass', 'all_targets_pass',
                       'all_negatives_locus_pass'))
    # 位点核查补 rank，便于与候选表直接对位（评审 R4）
    rank_by_id = {r['candidate_id']: r.get('rank', '')
                  for r in all_kept}
    for r in all_t_rows:
        r['rank'] = rank_by_id.get(r['candidate_id'], '')
    for r in all_n_rows:
        r['rank'] = rank_by_id.get(r['candidate_id'], '')
    _write_tsv(os.path.join(out_root, 'target_checks.tsv'),
               ['segment', 'candidate_id', 'rank', 'target_record',
                'observed_PAM',
                'PAM_match', 'observed_seed', 'seed_match',
                'effective_seed_mm', 'effective_spacer_mm',
                'effective_nonseed_mm_count', 'target_pass',
                'observed_window'],
               all_t_rows,
               yes_no=('PAM_match', 'seed_match', 'target_pass'))
    files_extra = {'target_checks': os.path.join(out_root,
                                                 'target_checks.tsv')}
    if cfg['negative']:
        _write_tsv(os.path.join(out_root, 'negative_checks.tsv'),
                   ['segment', 'candidate_id', 'rank', 'negative_record',
                    'observed_PAM', 'PAM_residual', 'PAM_changed',
                    'observed_seed', 'seed_match', 'effective_seed_mm_count',
                    'effective_seed_mm', 'negative_locus_pass',
                    'observed_window'],
                   all_n_rows,
                   yes_no=('PAM_residual', 'PAM_changed', 'seed_match',
                           'negative_locus_pass'))
        files_extra['negative_checks'] = os.path.join(
            out_root, 'negative_checks.tsv')

    summary = {
        'engine': 'VirusPlatform conan_design ' + ENGINE_VERSION,
        'ts': time.strftime('%Y-%m-%d %H:%M:%S'),
        'params': {k: v for k, v in cfg.items() if k != 'ignored'},
        'params_snapshot': {
            'pam': cfg['pam'], 'spacer_len': cfg['spacer_len'],
            'seed_len': cfg['seed_len'], 'ignored': sorted(cfg['ignored']),
            'max_nonseed_mm': cfg['max_nonseed_mm'],
            'terminal_margin': cfg['terminal_margin'],
            'gc': [cfg['gc_min'], cfg['gc_max']],
            'homopolymer': cfg['homopolymer'],
        },
        'segments': seg_meta,
        'n_candidates': len(all_kept),
        # top 候选行（≤top_n）：summary.json 自带，供前端结果面板直读
        'top': [
            {'segment': r['segment'], 'candidate_id': r['candidate_id'],
             'strand': r['strand'], 'cons_start': r['cons_start'],
             'cons_end': r['cons_end'],
             'guide': r['guide_spacer_5to3_RC'],
             'pam_len': len(cfg['pam']),
             'gc': r['spacer_GC_percent'],
             'coverage': r['target_coverage_rate'],
             'leak': r['negative_leak_rate'],
             'score': r['score'],
             'all_targets_pass': r['all_targets_pass'],
             'all_negatives_locus_pass': r['all_negatives_locus_pass']}
            for r in all_kept[:cfg['top_n']]],
        # 位点核查矩阵（候选 × 记录；类比 CONAN-SWIFT 图1e 的可计算版：
        # 论文那格是湿实验交叉反应矩阵，这里是计算侧核查状态）
        'matrices': matrices,
        'elapsed_sec': 0.0,
        'files': {},
    }
    summary_path = os.path.join(out_root, 'summary.json')
    report_path = os.path.join(out_root, 'report.md')
    summary['files'] = {
        'candidates': 'conan/candidates.tsv',
        'report': 'conan/report.md',
        'summary': 'conan/summary.json',
    }
    summary['files'].update({
        k: 'conan/' + os.path.basename(v)
        for k, v in files_extra.items()})
    summary['files'].update({
        k: 'conan/' + os.path.basename(v)
        for k, v in plot_files.items()})
    summary['elapsed_sec'] = round(time.time() - t0, 1)
    with safe_open(summary_path, 'wt') as f:
        json.dump(summary, f, ensure_ascii=False, indent=1)
    write_report(report_path, summary, cfg, all_kept[:cfg['top_n']])
    try:
        import shutil as _sh
        _sh.rmtree(work_dir, ignore_errors=True)
    except OSError:
        pass
    _prog('完成', 1.0, '完成')
    _log('CONAN-SWIFT 设计完成：%d 候选（%d 个片段组，耗时 %.1fs）'
         % (len(all_kept), len(seg_names), summary['elapsed_sec']))
    return summary, {'segments': segments_out}
