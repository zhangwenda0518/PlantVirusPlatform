# -*- coding: utf-8 -*-
"""miRNA 靶标比对方向（P0-4）与联动修复的回归测试。

覆盖（不依赖外部引擎的用例为主；ViennaRNA/primer3 缺失时自动跳过）：
  * P0-4 正控/反控：target 含 revcomp(miRNA) 判中、含 comp(miRNA) 同向镜像不判中；
  * 靶标坐标对准插入段、query 坐标回原始 5'→3' 帧；
  * ViennaRNA duplexfold 交叉验证（能量与位点区间一致）；
  * P1-2 tapir 列语义/种子区（第 3 位错配=种子错配、第 16 位=种子保留）
    与「target 短于 miRNA 不 IndexError」；
  * P1-3 psRNATarget score_alignment 完美配对 penalty=0；
  * P3-24 psRNATarget U→T 归一 + 择优 + 坐标区间；
  * P2-9 bin_engines miRanda Q/R 区间解析（含 Score 数字干扰）与 Energy 行窗口；
  * P2-10 primer_thermo calc_end_stability 参数方向（F 3' vs R 原序列）。
"""
import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from Virus_Platform_Core.mirna_target.alignment import (   # noqa: E402
    align_mirna_to_target)
from Virus_Platform_Core.mirna_target import tapir as T     # noqa: E402
from Virus_Platform_Core.mirna_target import psrnatarget as P  # noqa: E402
from Virus_Platform_Core.mirna_target import bin_engines as B  # noqa: E402
from Virus_Platform_Core.mirna_target.consensus import (   # noqa: E402
    compute_seed_quality, run_consensus_prediction, to_dict)

MIRNA = 'UUUGGAUUGAAGGGAGCUCUA'          # ath-miR159a 长度 21
L = len(MIRNA)
_COMP = {'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G'}
FLANK, TAIL = 'GGGGG', 'AAAAA'
# revcomp(MIRNA) 插入段在 target 上的 1-based 闭区间
SITE_LO = len(FLANK) + 1
SITE_HI = len(FLANK) + L


def revcomp(s):
    return ''.join(_COMP[c] for c in reversed(s))


def comp(s):
    return ''.join(_COMP[c] for c in s)


def pos_target():
    """正控 target：flank + revcomp(miRNA)（真实物理位点）+ tail。"""
    return FLANK + revcomp(MIRNA) + TAIL


def neg_target():
    """反控 target：flank + comp(miRNA)（同向互补，镜像位点，非物理位点）。"""
    return FLANK + comp(MIRNA) + TAIL


# ---------------------------------------------------------------- P0-4 方向

def test_positive_control_revcomp_site_hit():
    """正控：target 含 revcomp(miRNA) → identity≈100、seed_ok、坐标对准插入段。"""
    r = align_mirna_to_target('m', MIRNA, 't', pos_target())
    assert r is not None
    assert r.identity > 95.0, f'identity={r.identity}'
    assert r.seed_perfect is True
    assert (r.target_start, r.target_end) == (SITE_LO, SITE_HI)
    # query 坐标为 miRNA 原始 5'→3' 帧且 start < end
    assert (r.query_start, r.query_end) == (1, L)
    # 比对串按 miRanda 书写约定：上行 miRNA 3'→5'
    assert r.alignment_str.startswith("miRNA  3' ")
    assert "target 5' " in r.alignment_str
    assert r.tier_strict is True


def test_negative_control_comp_site_not_hit():
    """反控：target 含 comp(miRNA)（同向互补）→ 不得判中（旧实现此处 identity=100）。"""
    r = align_mirna_to_target('m', MIRNA, 't', neg_target())
    assert r is None or not r.tier_strict, (
        r and f'identity={r.identity} seed={r.seed_perfect}')
    assert r is None or r.identity < 90.0, (r and r.identity)


def test_vienna_duplexfold_cross_check():
    """ViennaRNA 交叉验证：revcomp 位点 MFE 显著低于镜像位点，且与比对能量一致。"""
    RNA = pytest.importorskip('RNA')
    mfe_pos = RNA.duplexfold(MIRNA, revcomp(MIRNA)).energy
    mfe_neg = RNA.duplexfold(MIRNA, comp(MIRNA)).energy
    assert mfe_pos < mfe_neg - 15.0, f'{mfe_pos} vs {mfe_neg}'
    r = align_mirna_to_target('m', MIRNA, 't', pos_target())
    # 位点区间与比对能量一致（cofold 原始 5'→3' 帧）
    assert abs(r.energy - mfe_pos) < 3.0, f'{r.energy} vs {mfe_pos}'
    assert r.energy <= -20.0


def test_seed_quality_new_frame():
    """consensus.compute_seed_quality 在新帧下按 5' 起算，且 wobble 不计成 WC。"""
    r = align_mirna_to_target('m', MIRNA, 't', pos_target())
    sq = compute_seed_quality(r.alignment_str)
    assert sq['seed_wc'] == 7 and sq['seed_wobble'] == 0
    assert sq['is_seed_perfect'] is True

    # 人工构造：种子第 3 位（原始帧）换成 G:U wobble → seed_wobble=1
    mir_rev = MIRNA[::-1]
    row = [comp(c) for c in mir_rev]          # 完美互补的目标行
    i = L - 3                                  # 原始第 3 位 => 新帧列下标
    q = mir_rev[i]
    row[i] = {'G': 'U', 'U': 'G', 'C': 'U', 'A': 'G'}[q]   # wobble 伙伴
    aln = ("miRNA  3' %s 5'\n" % ''.join(mir_rev)
           + '         ' + ''.join(row) + '\n'
           + "target 5' %s 3'" % ''.join(row))
    sq2 = compute_seed_quality(aln)
    assert sq2['seed_wobble'] == 1 and sq2['seed_wc'] == 6
    assert sq2['is_seed_perfect'] is False


# ---------------------------------------------------------------- P1-2 tapir

def _tapir_rows(mismatch_pos=None):
    """新帧比对列（上行 miRNA 3'→5'）；mismatch_pos = 原始 5' 帧位次（None=全配对）。"""
    mir_rev = MIRNA[::-1]
    row = [comp(c) for c in mir_rev]
    if mismatch_pos:
        i = L - mismatch_pos
        q = T._norm_rna(T._revcomp_char(mir_rev[i]))
        for cand in 'ACGU':
            if T._classify_column(q, cand) == 'mismatch':
                row[i] = cand
                break
        else:
            raise AssertionError('no mismatch base')
    return mir_rev, ''.join(row)


def test_tapir_seed_band_positions():
    """第 3 位错配 → 种子错配；第 16 位错配 → 种子保留（种子区 = miRNA 第 2-12 位）。"""
    mir_rev, row = _tapir_rows(3)
    s3 = T.score_alignment_columns(mir_rev, row)
    assert s3[7] == 1 and s3[6] == 0, f'seed_mismatch={s3[7]} nonseed={s3[6]}'
    assert s3[0] >= 2.0

    mir_rev, row = _tapir_rows(16)
    s16 = T.score_alignment_columns(mir_rev, row)
    assert s16[7] == 0 and s16[6] == 1, f'seed_mismatch={s16[7]} nonseed={s16[6]}'

    # 完美比对：满分 0，无错配无缺口
    mir_rev, row = _tapir_rows(None)
    s0 = T.score_alignment_columns(mir_rev, row)
    assert s0[0] == 0.0 and s0[7] == 0 and s0[6] == 0


def test_tapir_short_target_no_indexerror():
    """target 短于 miRNA：按缺口补齐，不再 IndexError。"""
    mir_rev, _ = _tapir_rows(None)
    s = T.score_alignment_columns(mir_rev, 'ACGUACGU')     # 8 < 21
    assert s[0] >= 0.0


def test_tapir_validate_uses_alignment_columns():
    """validate_with_tapir 用比对列（而非原始 miRNA vs 片段）打分：完美位点 score=0。"""
    r = align_mirna_to_target('m', MIRNA, 't', pos_target())
    tr = T.validate_with_tapir(r, 'm', 't', pos_target())
    assert tr is not None
    assert tr.tapir_score == 0.0, tr.tapir_score


# ---------------------------------------------------------------- P1-3 / P3-24

def test_psrna_wc_branch_perfect_pairing():
    """P1-3：完美配对 penalty=0、mismatches=0（原先缺 'wc' 分支按错配计分）。"""
    q = list('TAGAGCTCCCTTCAATCCAAA')          # 19 列内全部 q==h
    penalty, seed_mm, total_mm, gaps, gus, n = P.score_alignment(q, q)
    assert penalty == 0.0
    assert (seed_mm, total_mm, gaps, gus) == (0, 0, 0, 0)
    assert n == 19


def test_psrna_evaluate_coords_and_u_input():
    """P3-24：坐标对准插入段；含 U 的 RNA target 不被 revcomp 打成 N。"""
    tgt = pos_target()
    r = P.evaluate_mirna_vs_target(MIRNA, tgt, 'm', 't')
    assert r is not None and r.passes
    assert (r.target_start, r.target_end) == (SITE_LO, SITE_HI)
    # 全 U 版 target：坐标与通过性一致（U→T 归一后等价）
    r_u = P.evaluate_mirna_vs_target(MIRNA, tgt.replace('T', 'U'), 'm', 't')
    assert r_u is not None and r_u.passes
    assert (r_u.target_start, r_u.target_end) == (SITE_LO, SITE_HI)
    assert P._revcomp_dna(P._dna('ACGU')) == 'ACGT'


# ---------------------------------------------------------------- P2-9 / P2-12

def test_miranda_block_parsing(tmp_path, monkeypatch):
    """miRanda 命中块解析：Q/R 区间不取到 Score 数字，Energy 行在窗口内。"""
    block = (
        "Performing Scan: mir vs tgt\n"
        "   Forward:\tScore: 195.000000  Q:2 to 20  R:6 to 26 "
        "Align Len (18) (100.00%) (100.00%)\n"
        "\n"
        "   Query:    3' auCUCGAGGGAAGUUAGGUUu 5'\n"
        "                  |||||||||||||||||| \n"
        "   Ref:      5' uaGAGCUCCCUUCAAUCCAAa 3'\n"
        "\n"
        "   Energy:  -39.000000 kCal/Mol\n"
    )

    class _Fake:
        stdout = block

    monkeypatch.setattr(B, 'MIRANDA_EXE', __file__)      # 只过 isfile 检查
    monkeypatch.setattr(B, '_run', lambda *a, **k: _Fake())
    res = B.run_miranda_batch([('mir', MIRNA)], [('tgt', pos_target())],
                             str(tmp_path))
    assert res and len(res['hits']) == 1
    h = res['hits'][0]
    assert (h['q_start'], h['q_end']) == (2, 20)         # 不是 (195, 2)
    assert (h['t_start'], h['t_end']) == (6, 26)
    assert h['energy'] == -39.0                          # 第 6 行也在窗口内
    assert h['query_aln'].startswith('Query:')           # 比对三行取到真实内容


def test_bin_engines_u_to_t():
    """P2-12：miRanda/RNAhybrid 输入写盘前 U→T（含小写 u→t），去缺口。"""
    assert B._to_dna('AUGCu-') == 'ATGCt'      # 保留原大小写
    assert B._to_dna('UUU-ggg') == 'TTTggg'


# ---------------------------------------------------------------- P2-10

def test_end_stability_direction():
    """P2-10：calc_end_stability(F, R) 判 3' 互配为结合（负 ΔG），不再传 revcomp(R)。"""
    primer3 = pytest.importorskip('primer3')
    from Virus_Platform_Core import primer_thermo as pt
    if not pt.BIO_OK:
        pytest.skip('Biopython 不可用')
    f, r = 'ACGTA', 'TACGT'                    # F 3'=ACGTA 与 R 3'=TACGT 互配
    out = pt.end_stability_analysis(f, r)
    assert out['fwd_rev_end_dg'] < 0, out
    old = primer3.calc_end_stability(f, str(__import__('Bio.Seq', fromlist=['Seq']).Seq(r).reverse_complement()))
    assert old.dg > 0                          # 旧写法（revcomp(R)）恰为不结合


# ---------------------------------------------------------------- 共识汇总字段

def test_consensus_engines_actual_and_notes():
    """P1-4/P1-5：engines_actual 如实、rna22 口径旁注、best_energy 只取通过引擎。"""
    cr = run_consensus_prediction('m', MIRNA, 't', 't', '', pos_target(),
                                  algorithms=['psrna', 'rna22'])
    assert cr is not None
    d = to_dict(cr)
    assert d['engines_actual'][0] == 'miranda'
    assert 'psrna' in d['engines_actual']
    assert d['engines_skipped'] == []
    assert d['best_energy'] == min(d['miranda_energy'], d['psrna_energy'])
    assert any('未使用 miRNA 序列' in n for n in d['engine_notes'])


def test_rna22_flagged_mirna_independent():
    """P1-5①：RNA22 结果标注 miRNA_independent=True（判定不看 miRNA 序列）。"""
    import dataclasses
    from Virus_Platform_Core.mirna_target.rna22 import (
        validate_with_rna22, RNA22Result)
    names = {f.name for f in dataclasses.fields(RNA22Result)}
    assert 'miRNA_independent' in names

    class _MR:
        mirna_seq = MIRNA
        target_start = 5
        target_end = 16

    # 构造同时含「GC 富集 7-mer 种子 + hairpin」的窗口（该引擎的通过条件）
    res = validate_with_rna22(_MR(), 'm', 't', 'AAAAAAAACCCCGGGGAAAA')
    assert res is not None and res.passes is True
    assert res.miRNA_independent is True


def test_consensus_module_typing_and_docstring():
    """P2-6：Optional 已导入（get_type_hints 不再 NameError）、__doc__ 非 None。"""
    import typing
    from Virus_Platform_Core.mirna_target import consensus as C
    assert C.__doc__, '__doc__ 为 None（docstring 顺序不对）'
    hints = typing.get_type_hints(C.ConsensusResult)
    assert 'miranda' in hints
    assert typing.get_type_hints(C._validate_psrna)


def test_consensus_engine_failure_recorded(monkeypatch):
    """P1-4：校验引擎抛异常 → engines_skipped 如实记录 + WARN，不再静默吞。"""
    from Virus_Platform_Core.mirna_target import consensus as C

    def _boom(*a, **k):
        raise ValueError('模拟引擎失败')

    monkeypatch.setattr(C, 'validate_with_tapir', _boom)
    cr = run_consensus_prediction('m', MIRNA, 't', 't', '', pos_target(),
                                  algorithms=['tapir'])
    assert cr is not None
    d = to_dict(cr)
    assert 'tapir' not in d['engines_actual']
    assert d['engines_skipped'] and d['engines_skipped'][0]['engine'] == 'tapir'
    assert 'ValueError' in d['engines_skipped'][0]['reason']


def test_interactive_primer_specificity_declared(monkeypatch):
    """P2-14④：交互式设计接入可选特异性检查；默认关时 summary 如实声明未检查。"""
    pytest.importorskip('primer3')
    import os
    import random
    from Virus_Platform_Core import primer_design as pd
    from Virus_Platform_Core import primer as P
    from Virus_Platform_Core import config as C
    random.seed(20261003)
    seq = ''.join(random.choice('ACGT') for _ in range(400))

    res0 = pd.run_design(fasta_text='>t1\n' + seq, ptype='PCR')
    assert res0['summary']['specificity_requested'] is False
    assert res0['summary']['specificity_checked'] is False
    assert '未做' in res0['summary']['specificity_note']

    # 打桩：宿主库存在 + blastn 命中，验证扩增子 FASTA 与逐对 host_hit 标注
    seen = {}

    def fake_check(amp_fa, host_genome=None, logger=None):
        txt = open(amp_fa, encoding='utf-8').read()
        seen['amp'] = txt
        names = [l[1:] for l in txt.splitlines() if l.startswith('>')]
        return {names[0]: ['host_ctg1']} if names else {}

    monkeypatch.setattr(C, 'current_host_genome',
                        lambda: os.path.abspath('examples/example_tmv.fasta'))
    monkeypatch.setattr(P, 'check_specificity_blast', fake_check)
    res1 = pd.run_design(fasta_text='>t1\n' + seq, ptype='PCR', specificity=True)
    if res1['summary']['n_pairs'] == 0:
        pytest.skip('本环境未设计出引物，跳过脱靶检查链路')
    assert res1['summary']['specificity_checked'] is True
    hits = [(p.get('host_hit')) for rec in res1['records'] for p in rec['pairs']]
    assert hits.count('YES') == 1 and set(hits) <= {'YES', 'NO'}
    assert seen['amp'].count('>') == res1['summary']['n_pairs']
