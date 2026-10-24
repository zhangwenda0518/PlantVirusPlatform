# -*- coding: utf-8 -*-
"""conan_design 纯函数单测 + UI 契约 + 合成数据端到端。

合成数据带已知答案：
- 目标 3 株埋 engineered 位点 AAG+SP32（应全目标通过、进入 Top）；
- 负组 2 条：一条在同源位点保留 AAG+种子（leak 案）、一条 PAM 已变
  且种子错配（背景干净案）→ 该候选 leak 率 = 50%、背景干净 = 否；
- 端到端跑 run_chain（目标+负组合并 MAFFT），断言产物齐全、
  engineered 候选进表且速率与理论值一致。

MAFFT 缺失时端到端用例自动跳过（其余用例为纯函数，无外部依赖）。
"""
import json
import os
import random
import re
import shutil
import uuid

import pytest

from Virus_Platform_Core import conan_design as cnd

TOOLS = 'webapp/templates/tools.html'
I18N = 'webapp/static/i18n.js'

_SP = 'ACGTTGCACTGGATCAAGTCAAGGTTACGATA'   # 32nt，GC 43.75%，种子无 ≥3 run
assert len(_SP) == 32


def _rand(rng, n):
    return ''.join(rng.choice('ACGT') for _ in range(n))


def _snv(rng, seq, k, exclude):
    """在 exclude（0 基区间列表）之外打 k 个 SNV。"""
    pos = set()
    while len(pos) < k:
        p = rng.randrange(len(seq))
        if any(a <= p <= b for a, b in exclude):
            continue
        pos.add(p)
    out = list(seq)
    for p in pos:
        alt = {'A': 'C', 'C': 'A', 'G': 'T', 'T': 'G'}[out[p]]
        out[p] = alt
    return ''.join(out)


def _write_fasta(path, recs):
    with open(path, 'w', encoding='utf-8', newline='\n') as f:
        for name, seq in recs:
            f.write('>%s\n%s\n' % (name, seq))


def _ident_map(seq):
    """恒等列映射（单测用无 gap 伪比对行）。"""
    return list(range(len(seq)))


def _mk_cand(consensus, strand='+'):
    """取 consensus 上 SP32 位点的候选（无负组扫描即可）。"""
    for c in cnd.scan_candidates(consensus, 'AAG', 32, 8):
        if c['protospacer'] == _SP and c['strand'] == strand:
            c['id'] = 'CNtest_00001'
            return c
    raise AssertionError('合成序列中未找到 engineered 位点')


# ── 纯函数 ──────────────────────────────────────────────────
def test_parse_ignored():
    assert cnd.parse_ignored('6') == {6}
    assert cnd.parse_ignored('6, 9') == {6, 9}
    assert cnd.parse_ignored('') == set()
    assert cnd.parse_ignored('x, 0, -3') == set()


def test_scan_candidates_both_strands():
    rng = random.Random(11)
    cons = _rand(rng, 50) + 'AAG' + _SP + _rand(rng, 50)
    cands = cnd.scan_candidates(cons, 'AAG', 32, 8)
    plus = [c for c in cands if c['protospacer'] == _SP]
    assert plus and plus[0]['strand'] == '+'
    assert plus[0]['guide'] == cnd.revcomp(_SP)
    assert plus[0]['seed'] == _SP[:8]
    cons2 = _rand(rng, 40) + cnd.revcomp('AAG' + _SP) + _rand(rng, 40)
    minus = [c for c in cnd.scan_candidates(cons2, 'AAG', 32, 8)
             if c['protospacer'] == _SP and c['strand'] == '-']
    assert minus, '负链位点未识别'


def _site_consensus():
    """共识 400bp，位点 AAG+SP32 起于 0 基 199（1 基 200–234）。"""
    rng = random.Random(20261020)
    return _rand(rng, 199) + 'AAG' + _SP + _rand(rng, 400 - 199 - 35)


def test_seed_ignored_position6():
    """种子第 +6 位错配不计（I-E 口径）；+5 位错配照常计。"""
    cons = _site_consensus()
    cand = _mk_cand(cons)
    span = cand['span']
    spacer0 = span[0] + 3                      # spacer 0 基起点

    def _eval(mut_cons):
        row, _, _ = cnd.evaluate_candidate(cand, [('t1', mut_cons)], [],
                                           8, {6}, 4, _ident_map(cons))
        return row, _

    # +6 位（0 基 5）错配 → 种子仍算匹配
    mut = list(cons)
    alt = {'A': 'C', 'C': 'A', 'G': 'T', 'T': 'G'}[mut[spacer0 + 5]]
    mut[spacer0 + 5] = alt
    row, _ = _eval(''.join(mut))
    assert row['target_seed_match_rate'] == 100.0
    # +5 位（0 基 4）错配 → 种子失配
    mut2 = list(cons)
    alt2 = {'A': 'C', 'C': 'A', 'G': 'T', 'T': 'G'}[mut2[spacer0 + 4]]
    mut2[spacer0 + 4] = alt2
    row2, _ = _eval(''.join(mut2))
    assert row2['target_seed_match_rate'] == 0.0


def test_negative_dual_fail_scoring():
    """负组规则：PAM+种子双保留 = leak；PAM 已变且种子错配 = 干净。"""
    cons = _site_consensus()
    cand = _mk_cand(cons)
    span = cand['span']
    leak_neg = list(cons)                      # 同源位点原样保留 → leak
    clean_neg = list(cons)
    # 干净案：PAM 首位与种子首位同变 → PAM 断、种子断
    clean_neg[span[0]] = 'A' if clean_neg[span[0]] != 'A' else 'G'
    seed0 = span[0] + 3
    clean_neg[seed0] = 'A' if clean_neg[seed0] != 'A' else 'G'
    row, _, n_rows = cnd.evaluate_candidate(
        cand, [('t1', cons)],
        [('leak', ''.join(leak_neg)), ('clean', ''.join(clean_neg))],
        8, {6}, 4, _ident_map(cons))
    assert row['negative_leak_rate'] == 50.0
    assert row['negative_PAM_residual_rate'] == 50.0
    assert row['negative_seed_match_rate'] == 50.0
    assert row['all_negatives_locus_pass'] is False
    # leak 案标记；干净案通过
    by = {r['negative_record']: r for r in n_rows}
    assert by['leak']['PAM_residual'] and by['leak']['seed_match']
    assert not by['leak']['negative_locus_pass']
    assert by['clean']['PAM_changed'] and not by['clean']['seed_match']
    assert by['clean']['negative_locus_pass']
    # 评分：100 + 100 − 50 − 50 = 100
    assert row['score'] == 100.0


def test_filters_terminal_gc_homopolymer():
    cons = _site_consensus()
    cand = _mk_cand(cons)
    row, _, _ = cnd.evaluate_candidate(cand, [('t1', cons)], [], 8, {6}, 4,
                                       _ident_map(cons))
    # 默认参数（末端 150 / GC 35–65 / 种子 run≥3）→ 全过
    ok, reason = cnd.apply_filters(dict(row), len(cons), {
        'terminal_margin': 150, 'exclude_terminal': True, 'select_gc': True,
        'gc_min': 35, 'gc_max': 65, 'exclude_homopolymer': True,
        'homopolymer': 3, 'seed_len': 8})
    assert ok, reason
    # 末端排除放大到 300 → 位点（200–234）进排除带
    ok2, reason2 = cnd.apply_filters(dict(row), len(cons), {
        'terminal_margin': 300, 'exclude_terminal': True, 'select_gc': False,
        'gc_min': 0, 'gc_max': 100, 'exclude_homopolymer': False,
        'homopolymer': 3, 'seed_len': 8})
    assert not ok2 and '末端' in reason2
    # GC 收窄到 50–60（SP32 GC 43.75%）→ GC 过滤命中
    ok3, reason3 = cnd.apply_filters(dict(row), len(cons), {
        'terminal_margin': 0, 'exclude_terminal': False, 'select_gc': True,
        'gc_min': 50, 'gc_max': 60, 'exclude_homopolymer': False,
        'homopolymer': 3, 'seed_len': 8})
    assert not ok3 and 'GC' in reason3


def test_ui_contract():
    """t-conan 卡 UI 契约：控件 / i18n 中英成对 / JS 接线 / 注册与导航。"""
    src = open(TOOLS, encoding='utf-8').read()
    i = src.index('<section class="card" id="t-conan">')
    j = src.index('</section>', i)
    card = src[i:j]
    for sel in ['id="cn_input"', 'id="cn_negative"', 'id="cn_pam"',
                'id="cn_spacer"', 'id="cn_seed"', 'id="cn_ignored"',
                'id="cn_maxmm"', 'id="cn_margin_on"', 'id="cn_gc_on"',
                'id="cn_homo_on"', 'id="cn_group"', 'id="cn_topn"',
                'id="toolrun-conan"', 'id="cnResult"', 'id="cnTable"',
                'id="cnFlow"', 'id="cnMatrices"']:
        assert sel in card, sel
    for fn in ('runConan', 'loadConanResult', '_conanRate', 'conanStruct'):
        assert ('function ' + fn) in src or (fn + '(') in src, fn
    for key in ('negative:', 'pam:', 'spacer_len:', 'seed_len:', 'ignored:',
                'max_nonseed_mm:', 'terminal_margin:', 'group_by_segment:',
                'top_n:'):
        assert key in src, key
    i18n = open(I18N, encoding='utf-8').read()
    keys = set(re.findall(r'data-i18n="(tk\.cn[A-Za-z]+)"', card))
    keys |= set(re.findall(r'data-i18n-ph="(tk\.cn[A-Za-z]+)"', card))
    assert keys, '卡内未找到 tk.cn* i18n 锚点'
    for k in sorted(keys):
        assert i18n.count("'%s'" % k) >= 2, k + ' 中英不成对'
    # 注册表与导航卡
    from Virus_Platform_Core.web.tools_api import TOOL_REGISTRY
    assert 'conan' in TOOL_REGISTRY
    from Virus_Platform_Core.web.pages import NAV_GROUPS
    ids = [it.get('id') for g in NAV_GROUPS if isinstance(g, dict)
           for it in g.get('items', []) if isinstance(it, dict)]
    assert 't-conan' in ids


# ── 端到端（需 MAFFT）────────────────────────────────────────
@pytest.fixture()
def synth(tmp_path):
    """目标 3 株（位点保留）+ 负组 2 条（1 leak 案 + 1 干净案）。

    负组按真实 CONAN 场景构造：目标骨架的近同源副本（~98% 同一，
    位点区外散布 SNV）。不要给负组全随机侧翼——35bp 位点孤岛嵌在无关
    序列里时 MAFFT 会产生 ±1 列配准漂移（exact-match 块也滑，合成数据
    实测复现）；真实负组（近缘病毒全基因组）高同源无此问题。
    """
    rng = random.Random(20261020)
    base = _site_consensus()
    t1 = base
    t2 = _snv(rng, base, 12, [(199, 233)])
    t3 = _snv(rng, base, 12, [(199, 233)])
    # leak 案：近同源副本、位点原样保留（PAM+种子 双双保留 → leak）
    leak = _snv(rng, base, 6, [(199, 233)])
    # 干净案：近同源副本、位点无歧义销毁（PAM 三位 TTC，零同源可滑；
    # 只改首位 AAG→GAG 会被 MAFFT 用 gap 滑掉——G-run 边界单碱基突变
    # 的比对器固有行为，实测复现）+ 种子首位变异
    clean = list(base)
    clean[199], clean[200], clean[201] = 'T', 'T', 'C'
    clean[202] = 'T'
    clean = _snv(rng, ''.join(clean), 6, [(199, 233)])
    targets = tmp_path / 'targets.fa'
    _write_fasta(targets, [('isoA', t1), ('isoB', t2), ('isoC', t3)])
    negatives = tmp_path / 'negatives.fa'
    _write_fasta(negatives, [('near_leak', leak),
                             ('near_clean', clean)])
    return {'targets': str(targets), 'negatives': str(negatives),
            'base': base}


def test_run_chain_synthetic(synth, tmp_path):
    if not _mafft_ok():
        pytest.skip('MAFFT 不可用')
    out_dir = str(_plat_dir(tmp_path))
    summary, rep = cnd.run_chain(
        synth['targets'],
        {'negative': synth['negatives'], 'threads': 0},
        out_dir)
    # 产物齐全
    root = os.path.join(out_dir, 'conan')
    for fn in ('candidates.tsv', 'target_checks.tsv',
               'negative_checks.tsv', 'summary.json', 'report.md',
               'top_candidates.png', 'locus_overview.png'):
        assert os.path.isfile(os.path.join(root, fn)), fn
    assert os.path.getsize(os.path.join(root, 'candidates.tsv')) > 0
    # 评审修复回归：seed / cons_mean 列、核查 rank 列、summary 时间戳
    with open(os.path.join(root, 'candidates.tsv'), encoding='utf-8') as f:
        header = f.readline().rstrip('\n').split('\t')
    assert 'seed' in header and 'cons_mean' in header
    assert summary.get('ts')
    with open(os.path.join(root, 'target_checks.tsv'), encoding='utf-8') as f:
        assert 'rank' in f.readline().rstrip('\n').split('\t')
    with open(os.path.join(root, 'negative_checks.tsv'), encoding='utf-8') as f:
        assert 'rank' in f.readline().rstrip('\n').split('\t')
    # summary 结构
    assert summary['engine'].startswith('VirusPlatform conan_design')
    assert summary['segments'][0]['n_targets'] == 3
    assert summary['segments'][0]['n_negatives'] == 2
    # engineered 候选进表：全目标通过、leak 50%（以 TSV 全列行为准）
    rows = _read_tsv_candidates(os.path.join(root, 'candidates.tsv'))
    hit = [c for c in rows
           if c.get('protospacer_5to3_after_PAM') == _SP]
    assert hit, 'engineered 位点未进入候选'
    h = hit[0]
    assert float(h['target_coverage_rate']) == 100.0
    assert float(h['negative_leak_rate']) == 50.0
    assert float(h['score']) == 100.0
    assert h['all_targets_pass'] == '是'
    assert h['all_negatives_locus_pass'] == '否'
    assert h['seed'] == _SP[:8]
    assert 0.0 < float(h['cons_mean']) <= 1.0
    # 位点核查矩阵（CONAN-SWIFT 图1e 计算类比）：目标格全 ok、leak 负组
    # 格 = leak、干净负组格 = ok
    assert summary.get('matrices'), 'matrices 缺失'
    m0 = summary['matrices'][0]
    assert m0['rows'] and len(m0['cols']) == 5        # 3 目标 + 2 负组
    mrow = next(r for r in m0['rows']
                if r['candidate_id'] == h['candidate_id'])
    t_vals = [v for k, v in mrow['cells'].items() if k.startswith('iso')]
    assert t_vals == ['ok', 'ok', 'ok']
    assert mrow['cells']['near_leak'] == 'leak'
    assert mrow['cells']['near_clean'] == 'ok'


def test_param_zero_passthrough(tmp_path):
    """评审 R1：数值 0 是合法参数（严格 0 错配档 / 关闭末端排除），
    不得被 or 吞成默认值。"""
    base = _site_consensus()
    d = _plat_dir(tmp_path)
    tf = os.path.join(d, 't.fa')
    _write_fasta(tf, [('isoA', base)])
    summary, _ = cnd.run_chain(tf, {'max_nonseed_mm': 0}, d)
    assert summary['params_snapshot']['max_nonseed_mm'] == 0
    summary2, _ = cnd.run_chain(tf, {'terminal_margin': 0}, d)
    assert summary2['params_snapshot']['terminal_margin'] == 0
    summary3, _ = cnd.run_chain(tf, {'select_gc': '1', 'gc_min': 0,
                                     'gc_max': 100}, d)
    assert summary3['params_snapshot']['gc'] == [0.0, 100.0]
    shutil.rmtree(d, ignore_errors=True)


def test_negative_segment_policy(tmp_path):
    """评审 R2：多片段组目标 × 单片段组负组 → 无同名标签的片段跳过
    负组评分（跨片段比对无同源语义，leak 率无意义）。"""
    if not _mafft_ok():
        pytest.skip('MAFFT 不可用')
    rng = random.Random(7)
    base = _site_consensus()
    leak = _snv(rng, base, 6, [(199, 233)])
    d = _plat_dir(tmp_path)
    tf = os.path.join(d, 't.fa')
    nf = os.path.join(d, 'n.fa')
    _write_fasta(tf, [('isoA_RNA1', base),
                      ('isoA_RNA2', _rand(rng, 300))])
    _write_fasta(nf, [('near_RNA1', leak)])
    summary, _ = cnd.run_chain(tf, {'negative': nf}, d)
    by = {s['segment']: s for s in summary['segments']}
    assert by['segRNA1']['negative_applied'] is True
    assert by['segRNA1']['n_negatives'] == 1
    assert by['segRNA2']['negative_applied'] is False
    assert by['segRNA2']['n_negatives'] == 0
    shutil.rmtree(d, ignore_errors=True)


def _mafft_ok():
    try:
        from Virus_Platform_Core.config import get_config
        return bool(os.path.isfile(str(get_config().tool('mafft'))))
    except Exception:
        return False


def _plat_dir(tmp_path):
    # check_path 强制写路径在平台根/输出根内 → 用平台 run/tests_conan 专目录
    d = os.path.join(os.path.dirname(os.path.dirname(
        os.path.abspath(__file__))), 'run', 'tests_conan',
        uuid.uuid4().hex[:10])
    os.makedirs(d, exist_ok=True)
    return d


def _read_tsv_candidates(path):
    """candidates.tsv → [dict]（补齐 summary['top'] 缺的列）。"""
    import csv
    with open(path, encoding='utf-8') as f:
        return list(csv.DictReader(f, delimiter='\t'))
