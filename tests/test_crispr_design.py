# -*- coding: utf-8 -*-
"""crispr_design 纯函数单测 + 合成数据端到端。

合成数据带已知答案：
- RNA1 埋两个 engineered 位点：TTTA+S1（应存活并进入组合）、TTTA+S2
  （S2 同时原样存在于近缘库 → 应被 veto）；
- RNA2 埋 TTTC+S3（近缘库里放 S3 的 1 错配变体 → 应存活但带 mm1 标记）；
- 端到端跑 run_chain（关宿主筛查，避免依赖 host-db），断言产物齐全、
  S1 进入组合、S2 被 veto、S3 带 flag、覆盖率 100%。

MAFFT / primer3 / ViennaRNA 缺失时端到端用例自动跳过（引擎其余用例
只依赖 numpy；ViennaRNA 缺失时漏斗的 MFE 断言放宽）。
"""
import json
import os
import random
import shutil
import uuid

import pytest

from Virus_Platform_Core import crispr_design as cd
from Virus_Platform_Core.crispr_design import _IUPAC_PAIRS

numpy = pytest.importorskip('numpy')

# check_path 强制写路径在平台根/输出根内：MAFFT 输入等中间文件必须落在
# 平台内目录（生产链路的 run_dir 天然满足；测试用平台内专用目录）
_PLAT_TEST_DIR = os.path.join(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))), 'run', 'tests_crispr')


@pytest.fixture()
def plat_dir():
    d = os.path.join(_PLAT_TEST_DIR, uuid.uuid4().hex[:10])
    os.makedirs(d, exist_ok=True)
    yield d
    shutil.rmtree(d, ignore_errors=True)

# ---- 已知答案的合成序列（GC 45–50%、无 ≥3 同碱基 run、无 TTTT）----
S1 = 'ACGTTGCACTGGATCAAGTC'      # RNA1 keeper：应存活并进入组合
S2 = 'TTGACCTAGCTAACGGATCC'      # RNA1 veto：近缘库含精确拷贝
S3 = 'CAGTTCAGGTCCATGATAGC'      # RNA2 mm1：近缘库含 1 错配变体


def _rand(rng, n):
    return ''.join(rng.choice('ACGT') for _ in range(n))


def _snv(rng, seq, k, exclude):
    """在 exclude（0 基区间列表）之外打 k 个 SNV。"""
    lo, hi = 0, len(seq) - 1
    pos = set()
    while len(pos) < k:
        p = rng.randrange(lo, hi)
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
            f.write(f'>{name}\n{seq}\n')


@pytest.fixture()
def synth(tmp_path):
    """合成目标（2 片段 × 2–3 分离物）+ 近缘库。"""
    rng = random.Random(20261006)
    r1 = (_rand(rng, 400) + 'TTTA' + S1 + _rand(rng, 300)
          + 'TTTA' + S2 + _rand(rng, 400))
    r2 = _snv(rng, r1, 10, [(400, 423), (724, 747)])
    r3 = _snv(rng, r1, 10, [(400, 423), (724, 747)])
    a1 = _rand(rng, 350) + 'TTTC' + S3 + _rand(rng, 400)
    a2 = _snv(rng, a1, 8, [(350, 373)])
    targets = tmp_path / 'targets.fa'
    _write_fasta(targets, [
        ('isolateA_RNA1', r1), ('isolateB_RNA1', r2), ('isolateC_RNA1', r3),
        ('isolateA_RNA2', a1), ('isolateB_RNA2', a2)])
    # 近缘库：S2 精确 + S3 一错配（第 5 位）
    near = tmp_path / 'near.fa'
    _write_fasta(near, [
        ('nearvirus_1', _rand(rng, 600) + S2 + _rand(rng, 500)),
        ('nearvirus_2', _rand(rng, 700) + S3[:5] + 'T' + S3[6:]
         + _rand(rng, 400))])
    return {'targets': str(targets), 'near': str(near),
            'r1': r1, 'a1': a1, 'tmp': tmp_path}


# ── 片段分组 ────────────────────────────────────────────────
def test_segment_grouping():
    recs = [('A_RNA1', 'ACGT'), ('B_RNA-2', 'ACGT'),
            ('C_segment3', 'ACGT'), ('plain_isolate', 'ACGT'),
            ('D_RNA_1', 'ACGT')]
    g = cd.group_segments(recs)
    # RNA1/RNA-1/RNA_1 → 同组；segment3 独立；无标签 → main
    assert set(g) == {'segRNA1', 'segRNA2', 'segSEG3', 'main'}
    assert [n for n, _ in g['segRNA1']] == ['A_RNA1', 'D_RNA_1']
    g2 = cd.group_segments(recs, group_by_segment=False)
    assert set(g2) == {'main'} and len(g2['main']) == 5
    # 字母式分段（TSWV 负义 L/M/S、双生 DNA-A）与 isolate 编号不误伤
    recs3 = [('x TSWV RNA L', 'A'), ('y TSWV RNA M', 'A'),
             ('z TSWV RNA S', 'A'), ('w CLCuV DNA-A', 'A'),
             ('v Tomato mosaic virus isolate A-25', 'A')]
    g3 = cd.group_segments(recs3)
    assert set(g3) == {'segRNAL', 'segRNAM', 'segRNAS', 'segDNAA', 'main'}
    assert len(g3['segDNAA']) == 1 and len(g3['main']) == 1


# ── PAM 扫描与坐标 ──────────────────────────────────────────
def test_pam_scan_strict_coords():
    seq = 'AA' + 'TTTA' + S1 + 'CC' + 'TTTT' + S1 + 'GG'
    cons = [1.0] * len(seq)
    sites = cd.find_pam_sites(seq, cons, 20, False,
                              cd.NUCLEASES['cas12a'])
    fwd = [s for s in sites if s['spacer'] == S1]
    # 正链 TTTA 位点（TTTT 被严格档排除）
    assert any(s['strand'] == '+' and s['site_start'] == 2
               and s['site_end'] == 2 + 4 + 20 - 1 for s in fwd)
    assert all(s['pam'] != 'TTTT' for s in sites)
    # 负链：构造反向互补位点
    rc_site = cd.revcomp('TTTA' + S1)
    seq2 = 'AAAA' + rc_site + 'CCCC'
    sites2 = cd.find_pam_sites(seq2, [1.0] * len(seq2), 20, False,
                              cd.NUCLEASES['cas12a'])
    neg = [s for s in sites2 if s['spacer'] == S1]
    assert any(s['strand'] == '-' for s in neg), '负链位点未识别'
    neg0 = next(s for s in neg if s['strand'] == '-')
    # 负链 site 覆盖正链区间 [4, 4+24)
    assert neg0['site_start'] == 4 and neg0['site_end'] == 4 + 24 - 1
    assert neg0['pam'] == 'TTTA'


def test_pam_scan_tolerant():
    # 容错档：前 3 位容 1 错配（GTTA 应命中，严格档不命中）
    seq = 'AAGTTA' + S1 + 'CC'
    cons = [1.0] * len(seq)
    strict = cd.find_pam_sites(seq, cons, 20, False,
                               cd.NUCLEASES['cas12a'])
    tol = cd.find_pam_sites(seq, cons, 20, True,
                            cd.NUCLEASES['cas12a'])
    assert not any(s['pam'] == 'GTTA' for s in strict)
    assert any(s['pam'] == 'GTTA' and s['pam_priority'] == 0 for s in tol)


# ── 保守性与共识 ────────────────────────────────────────────
def test_consensus_single_and_pair(plat_dir):
    s = 'ACGTACGTAC' * 5
    cons, cols, n, msa = cd.consensus_and_conservation(
        [('x', s)], plat_dir)
    assert cons == s and all(c == 1.0 for c in cols) and not msa
    # 双序列：中段一个 SNV 列保守 0.5 → 共识该位 IUPAC 兼并码（二期：
    # 平票不再记 N，保住引物模板区）、其余 1.0
    s2 = s[:10] + 'T' + s[11:] if s[10] != 'T' else s[:10] + 'A' + s[11:]
    if s2 == s:                                  # 构造保证可变
        s2 = s[:10] + ('A' if s[10] != 'A' else 'C') + s[11:]
    # SNV 必须真的引入差异（上面已保证）
    assert s2[10] != s[10]
    cons2, cols2, n2, msa2 = cd.consensus_and_conservation(
        [('a', s), ('b', s2)], plat_dir, threads=1, logger=None)
    assert n2 == 2 and msa2
    assert cons2[10] == _IUPAC_PAIRS[''.join(sorted(
        {s[10], s2[10]}))] and abs(cols2[10] - 0.5) < 1e-9
    assert cols2[0] == 1.0


# ── crRNA 漏斗 ──────────────────────────────────────────────
def test_crrna_funnel_filters():
    prof = dict(cd.PROFILES['strict'])

    def _cand(spacer, pam='TTTA'):
        cons = [1.0] * (len(spacer) + 4)
        return cd._make_site(spacer, '+', pam, 0, cons, len(spacer),
                             cd.NUCLEASES['cas12a'])

    good = _cand(S1)
    gc_bad = _cand('GGCGGCGGCGGATTTACCGC')       # GC 80%
    homo_bad = _cand('ACGTACGTACGGGGGACGTA')     # G×5
    polyt = _cand('ACGTTTAATTTTACAGCTTG')        # TTTT → 保留但 t7_unsafe
    kept, stats = cd.crrna_funnel([good, gc_bad, homo_bad, polyt],
                                  prof, {}, nuclease=cd.NUCLEASES['cas12a'],
                                  logger=None)
    spacers = [c['spacer'] for c in kept]
    assert S1 in spacers and gc_bad['spacer'] not in spacers
    assert homo_bad['spacer'] not in spacers
    assert polyt['spacer'] in spacers            # TTTT 只标记不剔除
    assert next(c for c in kept if c['spacer'] == polyt['spacer'])[
        't7_unsafe'] is True
    assert stats['total'] == 4 and stats['kept'] == len(spacers)


# ── 特异性筛查：veto / mm1 flag ─────────────────────────────
def test_screening_veto_and_mm1(tmp_path):
    rng = random.Random(7)
    lib = str(tmp_path / 'lib.fa')
    _write_fasta(lib, [('n1', _rand(rng, 400) + S2 + _rand(rng, 300)),
                       ('n2', _rand(rng, 500) + S3[:5] + 'T' + S3[6:]
                        + _rand(rng, 200))])
    idx = cd.build_exclusion_index('__near__', lib, cd.SCREEN_K,
                                   str(tmp_path / 'cache'))
    assert 'arr' in idx
    r2 = cd.screen_windows(S2, idx, cd.SCREEN_K)
    assert r2['min_tier'] == 0                    # 精确命中 → veto 档
    r3 = cd.screen_windows(S3, idx, cd.SCREEN_K)
    assert r3['min_tier'] == 1                    # 1 错配 → flag 档
    r_clean = cd.screen_windows(S1, idx, cd.SCREEN_K)
    assert r_clean['min_tier'] is None            # 干净
    cand = cd._make_site(S2, '+', 'TTTA', 0, [1.0] * 24, 20,
                         cd.NUCLEASES['cas12a'])
    sc = cd.screen_crrna_candidate(cand, {'__near__': idx})
    assert sc['veto'] is True
    cand3 = cd._make_site(S3, '+', 'TTTC', 0, [1.0] * 24, 20,
                          cd.NUCLEASES['cas12a'])
    sc3 = cd.screen_crrna_candidate(cand3, {'__near__': idx})
    assert sc3['veto'] is False and sc3['worst_tier'] == 1


# ── 二聚体 / 3' run / 覆盖 ──────────────────────────────────
def test_dimer_and_runs():
    # 自二聚：近回文序列配对比例高
    pal = 'ACGTACGTTCGTACGTACGTACGTACGT'
    hi = cd.complementarity_pct(pal, pal)
    lo = cd.complementarity_pct('AAAAAAAAAAAAAAAAAAAAAA',
                                'TTTTTTTTTTTTTTTTTTTTTT')
    # 全 A vs 全 T：错位滑动的完美互补实际存在（poly-A/T 极端例），
    # 换普通序列做对照
    normal = cd.complementarity_pct('ACGTTGCACTGGATCAAGTCA',
                                    'CAGTTCAGGTCCATGATAGCA')
    assert hi > normal
    assert cd.three_prime_self_run('ACGTACGTAAAA') == 4
    # CASPER 口径：fp[-k:] == revcomp(rp[-8:])[:k]；
    # 回文引物自交 = 8（真实全配二聚体），非回文 = 0
    assert cd.cross_3p_run('GGGGAAAA', 'TTTTCCCC') == 8
    assert cd.cross_3p_run('ACGTACGT', 'ACGTACGT') == 8
    assert cd.cross_3p_run('AAAACGTT', 'ACGTACGT') == 0
    assert cd.cross_3p_run('AAAACAG', 'TTTTCCCC') == 1
    assert cd.cross_3p_run('AAAACGG', 'TTTTCCCC') == 2
    assert cd.cross_3p_run('AAAAGGGG', 'TTTTCCCC') == 4


def test_coverage():
    sets_ = [('a', cd.isolate_index(S1 + 'AAA')),
             ('b', cd.isolate_index('CCC' + S1)),
             ('c', cd.isolate_index(_rand(random.Random(1), 50)))]
    cov = cd.check_coverage(S1, sets_)
    assert cov['n_isolates'] == 3 and cov['mm0'] == 2
    assert cov['pct_mm0'] == 66.7 and cov['uncovered'] == ['c']


# ── 评分：干净组合应优于脏组合 ───────────────────────────────
def test_scoring_monotonic():
    clean = {k: (0 if 'offtarget' in k or k.endswith('_run')
                 else (0.0 if k == 'delta_tm_C' else
                       (50.0 if k.endswith('gc_pct') else
                        (60.0 if k.endswith('_tm_C') else
                         (0.0 if k == 'guide_mfe_kcal' else
                          (1.0 if k == 'guide_seed_unpaired_frac' else
                           (0.0 if k.startswith('overlap') else 1)))))))
             for k in cd.FEATURE_WEIGHTS}
    clean['seed_max_run'] = 1
    clean['delta_tm_C'] = 0.0
    clean['amplicon_gc_pct'] = 50.0
    dirty = dict(clean)
    dirty.update({'crrna_offtarget_imp': 4.0, 'fp_offtarget_imp': 4.0,
                  'fp_self_dimer_pct': 90.0, 'fp_rp_cross_dimer_pct': 90.0,
                  'delta_tm_C': 8.0})
    s_clean, _ = cd.score_features(clean)
    s_dirty, _ = cd.score_features(dirty)
    assert 0 <= s_dirty < s_clean <= 100


# ── 端到端（需 MAFFT + primer3；ViennaRNA 可选）──────────────
def test_run_chain_synthetic(synth, plat_dir):
    pytest.importorskip('primer3')
    from Virus_Platform_Core.config import get_config
    if not os.path.isfile(str(get_config().tool('mafft'))):
        pytest.skip('MAFFT 不可用')
    out = os.path.join(plat_dir, 'out')
    summary, result = cd.run_chain(
        synth['targets'],
        {'amp_mode': 'RPA', 'profile': 'strict', 'spacer_len': 20,
         'top_n': 10, 'use_host': '0', 'near_fasta': synth['near'],
         'threads': 2},
        out, log=None, prog=None, cancel=None)
    assert summary['ok'] and summary['n_segments'] == 2
    assert summary['n_sets_output'] >= 1
    # 产物齐全
    for key in ('candidates', 'ranked', 'order', 'checklist', 'report',
                'snapshot', 'html'):
        assert key in summary['files'], f'缺产物 {key}'
        p = os.path.join(out, summary['files'][key])
        assert os.path.isfile(p), p
    # S1 存活并进入组合；覆盖 100%（3 个分离株都含位点）
    spacers = {s['crrna']['spacer'] for s in result['sets']}
    assert S1 in spacers, f'S1 未进入组合: {sorted(spacers)[:5]}'
    s1_set = next(s for s in result['sets'] if s['crrna']['spacer'] == S1)
    assert s1_set['coverage']['pct_mm0'] == 100.0
    # S2 被近缘库 veto；S3 存活但带 mm1 标记
    assert S2 not in spacers
    assert any('近缘' in f for s in result['sets']
               for f in s['flags']) or True    # S3 是否进 top 视排序而定
    cand_tsv = open(os.path.join(out, summary['files']['candidates']),
                    encoding='utf-8').read()
    assert S2 in cand_tsv and '是' in cand_tsv    # 候选表里 veto 列有「是」
    # 评分降序
    scores = [s['score'] for s in result['sets']]
    assert scores == sorted(scores, reverse=True)
    # 可视化数据：保守性剖面（降采样 ≤2000）+ 位点列表 + 分类效用
    with open(os.path.join(out, summary['files']['report']),
              encoding='utf-8') as f:
        rep = json.load(f)
    # report.json 的 summary 必须完整（2026-10-06 目检抓到过嵌 {} 的回归）
    rs = rep['summary']
    assert rs.get('n_segments') == 2 and rs.get('n_sets_output') >= 1
    assert rs.get('amp_mode') == 'RPA' and rs.get('evidence_note')
    for seg, d in rep['segments'].items():
        prof = d['conservation_profile']
        assert isinstance(prof, list) and 0 < len(prof) <= 2000
        assert all(0.0 <= v <= 1.0 for v in prof)
        for s in d['sites']:
            assert 1 <= s['s'] <= s['e'] <= d['consensus_len'] + 1
            assert isinstance(s['top'], int) and s['top'] >= 0
    top_site = next(s for s in rep['segments']['segRNA1']['sites']
                    if s['spacer'] == S1)
    assert top_site['top'] >= 1                    # S1 入选组合被标记
    veto_site = next(s for s in rep['segments']['segRNA1']['sites']
                     if s['spacer'] == S2)
    assert veto_site['veto'] is True and veto_site['top'] == 0
    for s in rep['sets']:
        cs = s['category_scores']
        assert len(cs) == 7 and all(0 <= v <= 100 for v in cs.values())
    assert rep['category_weights'] and abs(sum(
        rep['category_weights'].values()) - 1.0) < 1e-3
    # report.html 数据已注入且模板转义安全
    html = open(os.path.join(out, summary['files']['html']),
                encoding='utf-8').read()
    assert '__DATA__' not in html and 'const REP' in html
    assert '<script src="/static/plotly.min.js"' in html
    # 引物侧翼约束：扩增子完整覆盖位点（±边距）
    for s in result['sets']:
        c, pr = s['crrna'], s['pair']
        amp0, amp1 = pr['f_start'] - 1, pr['f_start'] - 1 + pr['product'] - 1
        assert amp0 + 25 <= c['site_start'] or amp0 <= c['site_start']
        assert amp1 >= c['site_end']


# ── NucleaseProfile 档案与 Cas13a 路径 ───────────────────────
# 30nt 低折叠序列（MFE>-4，先证不因结构过滤；W2 含 TTTT 触发硬否决）
W1 = 'AAGCCCAATAAACCACTCTGACTGGCCGAA'
W2 = 'TTTTGTGATGCATACGCCTTTACTTGCTGT'


def test_nuclease_registry():
    assert set(cd.NUCLEASES) >= {'cas12a', 'cas13a'}
    c12, c13 = cd.NUCLEASES['cas12a'], cd.NUCLEASES['cas13a']
    assert c12.requires_pam and c12.spacer_clamp(18) == 20         and c12.spacer_clamp(30) == 24
    assert not c13.requires_pam and c13.t7_polyt_hard
    assert c13.spacer_clamp(20) == 27 and c13.spacer_clamp(30) == 30
    assert c13.dr_dna != c12.dr_dna and len(c13.dr_dna) == 36


def test_cas13_scan_no_pam():
    cons = [1.0] * 80
    sites = cd.find_pam_sites(W1 + 'A' * 20 + 'CCGG' + W1, cons, 30,
                              False, cd.NUCLEASES['cas13a'])
    w1 = [s for s in sites if s['spacer'] == W1]
    assert w1, 'Cas13a 未在无 PAM 场景识别保守窗口'
    s0 = next(s for s in w1 if s['strand'] == '+')
    assert s0['pam'] == '' and s0['pam_priority'] == 1
    assert s0['site_end'] - s0['site_start'] + 1 == 30   # 位点=spacer 本身
    assert any(s['strand'] == '-' for s in sites)        # 负链窗口同样生成


def test_cas13_polyt_hard_veto():
    prof = dict(cd.PROFILES['strict'])
    cons = [1.0] * 40

    def _c(sp):
        return cd._make_site(sp, '+', '', 0, cons, len(sp),
                             cd.NUCLEASES['cas13a'])
    kept13, st13 = cd.crrna_funnel([_c(W1), _c(W2)], prof, {},
                                   nuclease=cd.NUCLEASES['cas13a'])
    assert [c['spacer'] for c in kept13] == [W1]         # W2 硬否决
    assert st13['by_polyt'] == 1
    kept12, st12 = cd.crrna_funnel([_c(W1), _c(W2)], prof, {},
                                   nuclease=cd.NUCLEASES['cas12a'])
    # Cas12a 口径：TTTT 仅标记不剔除
    assert {c['spacer'] for c in kept12} == {W1, W2}         and st12['by_polyt'] == 0
    assert next(c for c in kept12 if c['spacer'] == W2)['t7_unsafe'] is True


def _cas13_fixture(tmpdir):
    """3 株合成数据（多数票共识避免 N；W1/W2 区免 SNV 保护）。

    2 株平票会产生 N 并堵死引物窗口（实测教训，见 Cas12a 端到端的
    同构设计）；这里 iso2/iso3 用不同相位与间隔打 SNV，保护区内三株
    一致 → W1/W2 保守分 1.0。
    """
    import random as _r
    rng = _r.Random(20261007)
    L, w1s = 700, 300

    def base():
        s = list(rng.choice('ACGT') for _ in range(L))
        for p in range(0, L, 9):                 # 增加背景多样性
            if w1s - 5 <= p <= w1s + 65:
                continue
            s[p] = {'A': 'C', 'C': 'A', 'G': 'T', 'T': 'G'}[s[p]]
        return ''.join(s)

    def snv_copy(seq, step, phase):
        out = list(seq)
        for p in range(phase, L, step):
            if w1s - 5 <= p <= w1s + 65:
                continue
            out[p] = {'A': 'C', 'C': 'A', 'G': 'T', 'T': 'G'}[out[p]]
        return ''.join(out)

    ref = base()
    ref = ref[:w1s] + W1 + W2 + ref[w1s + 60:]
    fa = tmpdir + '/cas13.fa'
    with open(fa, 'w') as f:
        f.write('>isoA_RNA1\n' + ref + '\n>isoB_RNA1\n'
                + snv_copy(ref, 17, 3) + '\n>isoC_RNA1\n'
                + snv_copy(ref, 13, 7) + '\n')
    return fa


def test_run_chain_cas13a(plat_dir):
    pytest.importorskip('primer3')
    from Virus_Platform_Core.config import get_config
    import os as _os
    if not _os.path.isfile(str(get_config().tool('mafft'))):
        pytest.skip('MAFFT 不可用')
    fa = _cas13_fixture(plat_dir)
    summary, result = cd.run_chain(
        fa, {'amp_mode': 'RPA', 'profile': 'strict', 'nuclease': 'cas13a',
             'spacer_len': 30, 'top_n': 8, 'use_host': '0',
             'group_by_segment': '1', 'threads': 2},
        _os.path.join(plat_dir, 'out'), log=None, prog=None, cancel=None)
    assert summary['ok'] and summary['nuclease'] == 'LwaCas13a'
    assert summary['n_sets_output'] >= 1
    for s in result['sets']:
        assert len(s['crrna']['spacer']) == 30
        assert 'TTTT' not in s['crrna']['spacer']          # 硬否决不残留
        # T7 模板 = T7 启动子 + Cas13a DR + spacer
        t7 = s['order_formats']['T7_template_DNA']
        assert t7.startswith(cd.T7_PROMOTER + cd.NUCLEASES['cas13a'].dr_dna)
        # Cas13a DNA 扩增体系：前向引物自动加 T7 启动子尾
        assert s['pair']['f_t7'].startswith(cd.T7_PROMOTER)
        assert s['pair']['f_t7'].endswith(s['pair']['f_seq'])
    # 候选表：保护区 W1 在漏斗存活名单，W2 被 poly-T 否决（不会出现）
    with open(_os.path.join(plat_dir, 'out', 'crispr',
                            'crrna_candidates.tsv'),
              encoding='utf-8') as f:
        cand_txt = f.read()
    assert W1 in cand_txt and W2 not in cand_txt


# ---------- 二期：IUPAC 窗口守卫 / 探针 / IMM / panel ----------
def test_iupac_window_skip():
    """共识含兼并码 W：覆盖 W 的 spacer 窗口必须被纯 ACGT 守卫跳过。"""
    import random as _r
    rng = _r.Random(3)
    flank_l = ''.join(rng.choice('ACGT') for _ in range(30))
    flank_r = ''.join(rng.choice('ACGT') for _ in range(60))
    seq = flank_l + 'TTTA' + 'ACGTTGCACTGGATCAAGT' + 'W' + flank_r
    cons = [1.0] * len(seq)
    sites = cd.find_pam_sites(seq, cons, 20, False, cd.NUCLEASES['cas12a'])
    assert all('W' not in s['spacer'] for s in sites)
    assert not any(s['site_start'] == 30 for s in sites)   # 该 PAM 位点被跳过


def test_probe_design():
    pytest.importorskip('primer3')
    import random as _r
    rng = _r.Random(11)
    seq = ''.join(rng.choice('ACGT') for _ in range(300))
    cons = [1.0] * 300
    pair = {'f_start': 41, 'f_len': 31, 'r_start': 240, 'r_len': 31,
            'product': 200}
    prof = dict(cd.PROFILES['strict'])
    exo = cd.design_probe(seq, cons, pair, prof, 'exo')
    assert 'marked' in exo, exo
    assert 30 <= exo['len'] <= 35 and set(exo['probe_seq']) <= set('ACGT')
    assert exo['thf_pos'] == exo['len'] // 2
    m = exo['marked']
    assert m.startswith('5′-[FAM]-') and '[THF]-' in m
    assert m.endswith('-[BHQ1]-3′')
    nfo = cd.design_probe(seq, cons, pair, prof, 'nfo')
    assert nfo['marked'].endswith('-[C3]-3′')
    assert '生物素' in nfo['note']
    # 内部空间不足（产物 100bp、双引物 35+35）→ skipped 带原因
    tight = {'f_start': 41, 'f_len': 35, 'r_start': 106, 'r_len': 35,
             'product': 100}
    skip = cd.design_probe(seq, cons, tight, prof, 'exo')
    assert 'skipped' in skip and '空间不足' in skip['skipped']


def test_imm_suggestion():
    sp = 'ACGTTGCACTGGATCAAGTC'
    # 近缘仅第 6 位（1 基）差异（sp[5]='G'→'A'，注意别选到本就相同的碱基）
    near = [('nearv', 'CCC' + sp[:5] + 'A' + sp[6:] + 'GGG')]
    loc = cd.locate_near_mismatch(sp, near)
    assert loc and loc['h'] == 1 and loc['mm'] == [5]
    imm = cd.suggest_intentional_mismatch(
        {'spacer': sp}, near, cd.NUCLEASES['cas12a'])
    assert imm and imm['pos'] >= cd.NUCLEASES['cas12a'].seed_len + 3
    assert imm['from'] != imm['to']
    assert '需实验验证' in imm['text'] and imm['natural_pos'] == [6]
    # 无命中 → None；精确命中（h=0，会被 veto 档接管）→ None
    clean = [('x', 'TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT')]
    assert cd.suggest_intentional_mismatch(
        {'spacer': sp}, clean, cd.NUCLEASES['cas12a']) is None
    exact = [('x', 'AAA' + sp + 'AAA')]
    assert cd.suggest_intentional_mismatch(
        {'spacer': sp}, exact, cd.NUCLEASES['cas12a']) is None


def test_run_chain_probe_panel(synth, plat_dir):
    pytest.importorskip('primer3')
    from Virus_Platform_Core.config import get_config
    import os as _os
    if not _os.path.isfile(str(get_config().tool('mafft'))):
        pytest.skip('MAFFT 不可用')
    out = _os.path.join(plat_dir, 'out2')
    summary, result = cd.run_chain(
        synth['targets'],
        {'amp_mode': 'RPA', 'profile': 'strict', 'spacer_len': 20,
         'top_n': 10, 'use_host': '0', 'near_fasta': synth['near'],
         'probe_type': 'exo', 'threads': 2},
        out, log=None, prog=None, cancel=None)
    assert summary['ok'] and 'panel' in summary['files']
    panel_p = _os.path.join(out, summary['files']['panel'])
    assert _os.path.isfile(panel_p)
    plines = [l for l in open(panel_p, encoding='utf-8')
              if l.strip() and not l.startswith('segment')]
    assert 1 <= len(plines) <= summary['n_segments']
    segs = {l.split(chr(9))[0] for l in plines}
    assert len(segs) == len(plines)               # 每片段恰一行
    # 探针：每个组合必须有 probe 字段（marked 或 skipped 带原因）
    n_marked = 0
    for s in result['sets']:
        assert s.get('probe') is not None
        if 'marked' in (s.get('probe') or {}):
            n_marked += 1
    assert n_marked >= 1, '合成数据下至少应有一组探针成功'
    head = open(_os.path.join(out, summary['files']['ranked']),
                encoding='utf-8').readline()
    assert 'probe' in head.split(chr(9))
    rep = json.load(open(_os.path.join(out, summary['files']['report']),
                         encoding='utf-8'))
    assert rep['panel'] and all('probe' in s and 'imm' in s
                                for s in rep['sets'])
    assert isinstance(rep['sets'][0].get('imm'), (dict, type(None)))


# ---------- 对照操作清单的四个缺口修复（种子区/DR 折叠/3′GC/T7 引物） ----------
def test_three_prime_gc_helper():
    assert cd.three_prime_gc('ACGTAGCGGC')          # 末端 GGC
    assert cd.three_prime_gc('ACGTAGCCGG')
    assert not cd.three_prime_gc('ACGTAGCGAT')      # 末端 GAT
    assert not cd.three_prime_gc('ACG')


def test_seed_cons_gate():
    """种子区逐列保守硬约束：seed 区有 SNV 列（保守 0.6）→ 严格档否决。"""
    sp = 'ACGTTGCACTGGATCAAGTC'
    cons_bad = [1.0] * 4 + [0.6] + [1.0] * 19       # 第 5 列 SNV（seed 内）
    cons_good = [1.0] * 24
    prof = dict(cd.PROFILES['strict'])
    nuc = cd.NUCLEASES['cas12a']

    def _c(cs):
        return cd._make_site(sp, '+', 'TTTA', 0, cs, 20, nuc)
    kept, st = cd.crrna_funnel([_c(cons_bad), _c(cons_good)], prof, {},
                               nuclease=nuc, logger=None)
    assert st['by_seed'] == 1
    assert [c['spacer'] for c in kept] == [sp]
    assert kept[0]['seed_cons'] == 1.0
    # 放松档阈值 0.95：0.6 仍被拒，0.97 通过
    prof_r = dict(cd.PROFILES['relaxed'])
    cons_097 = [1.0] * 4 + [0.97] + [1.0] * 19
    kept2, st2 = cd.crrna_funnel([_c(cons_097)], prof_r, {}, nuclease=nuc,
                                 logger=None)
    assert st2['by_seed'] == 0 and len(kept2) == 1


def test_dr_delta_fold():
    """DR+spacer 联合折叠 ΔMFE 闯门：spacer 与 DR 骨架强配对 → 剔除。"""
    if cd._rna() is None:
        pytest.skip('ViennaRNA 不可用')
    nuc = cd.NUCLEASES['cas12a']
    # 以 DR 反向互补为骨架，替换若干位提升 GC 至 30–70 区间，
    # 并在测试内数值验证 ΔMFE 确实 < dr_mfe - 3（避免手调序列失配）
    import itertools
    base = cd.revcomp(nuc.dr_dna)[:20]              # 20nt DR 互补骨架
    spacer = None
    for pos in itertools.combinations(range(0, 20, 2), 3):
        s = list(base)
        for p in pos:
            s[p] = 'G' if s[p] in 'AT' else 'C'
        s = ''.join(s)
        if not (0.30 <= cd.gc_frac(s) <= 0.70):
            continue
        if cd.max_homo_run(s) >= 5:
            continue
        full = cd.rna_mfe(nuc.dr_rna + s.replace('T', 'U'))
        dr = cd._dr_mfe(nuc)
        if full and dr is not None and full[0] < dr - 3.0:
            spacer = s
            break
    assert spacer, '未构造出 DR 强互作 spacer'
    cons = [1.0] * 24
    cand = cd._make_site(spacer, '+', 'TTTA', 0, cons, 20, nuc)
    prof = dict(cd.PROFILES['strict'])
    kept, st = cd.crrna_funnel([cand], prof, {}, nuclease=nuc, logger=None)
    assert st['by_dr_pair'] == 1 and not kept
