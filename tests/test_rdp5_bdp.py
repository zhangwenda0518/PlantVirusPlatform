# -*- coding: utf-8 -*-
"""rdp5_engine 断点分布兜底单测（2026-10-12 修复配套）。

背景：RDP5CL -rbdp 只统计双亲本都已鉴定的事件，事件带 Unknown 亲本时
（真实数据常见）分布 CSV 静默全 0——10 条 CMV RNA3 实测复现，示例比对
删到 3 条序列（迫使次亲本 Unknown）同样复现。这里覆盖退化判定与重建。
"""
import os

from Virus_Platform_Core.rdp5_engine import (
    BDP_CSV,
    breakpoint_distribution_degenerate,
    orf_csv_to_aln_map,
    synthesize_breakpoint_distribution,
)

PLATFORM = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _write_csv(path, lines):
    with open(path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(lines) + '\n')


def test_degenerate_missing_file(tmp_path):
    assert breakpoint_distribution_degenerate(str(tmp_path / 'nope.csv'))


def test_degenerate_all_zero(tmp_path):
    p = str(tmp_path / BDP_CSV)
    _write_csv(p, [
        'Position in alignment, Recombination breakpoint number (200nt win),'
        ' Upper 99% CI, Lower 99% CI',
        '1,0,0,0', '3,0,0,0', '5,0,0,0',
    ])
    assert breakpoint_distribution_degenerate(p)


def test_degenerate_real_data_is_not(tmp_path):
    p = str(tmp_path / BDP_CSV)
    _write_csv(p, ['Position in alignment, Recombination breakpoint number',
                   '1,0.56,0.82,0', '3,0.55,0.82,0'])
    assert not breakpoint_distribution_degenerate(p)


def test_degenerate_platform_example_is_not():
    """examples 里 RDP5 正常产物（有断点密度）必须判为非退化。"""
    p = os.path.join(PLATFORM, 'examples', 'results', 'rdp', 'rdp', BDP_CSV)
    if os.path.isfile(p):
        assert not breakpoint_distribution_degenerate(p)


def test_synthesize_counts_sliding_window(tmp_path):
    p = str(tmp_path / BDP_CSV)
    evs = [{'start': 1126, 'end': 1212},          # 线性事件
           {'bp_start': 2164, 'bp_end': 366}]     # 环状事件（键名兼容）
    n_rows = synthesize_breakpoint_distribution(p, evs, aln_len=2292)
    assert n_rows == (2292 + 1) // 2               # 步长 2 的位置网格
    rows = {}
    with open(p, encoding='utf-8') as f:
        lines = f.read().splitlines()
    assert lines[0].startswith('Position in alignment')
    for ln in lines[1:]:
        pos, v = ln.split(',')
        rows[int(pos)] = int(v)
    # 两断点同窗（[1101,1301) 覆盖 1126 与 1212）→ 计 2
    assert rows[1101] == 2
    # 只覆盖单断点的窗 → 计 1
    assert rows[941] == 1 and rows[1181] == 1
    # 环状事件两个断点按位置直计
    assert rows[2101] == 1 and rows[301] == 1
    # 远离断点处为 0
    assert rows[1601] == 0


def test_synthesize_without_aln_len_uses_last_breakpoint(tmp_path):
    p = str(tmp_path / BDP_CSV)
    n_rows = synthesize_breakpoint_distribution(p, [{'start': 500, 'end': 800}])
    assert n_rows > 0
    last_pos = max(int(ln.split(',')[0])
                   for ln in open(p, encoding='utf-8').read().splitlines()[1:])
    assert last_pos >= 800


def test_synthesize_no_breakpoints_writes_nothing(tmp_path):
    p = str(tmp_path / BDP_CSV)
    assert synthesize_breakpoint_distribution(p, [{}], aln_len=100) == 0
    assert not os.path.isfile(p)


# ── orf_csv_to_aln_map：基因坐标 → 比对坐标 ─────────────────
def _write_aln(path, ref_gaps=None):
    """两条 12nt 序列的迷你比对；ref 行可在指定列插 gap。"""
    ref = 'ACGTACGTACGT'
    alt = 'ACGTACGTACGT'
    if ref_gaps:
        ref = ''.join('-' if i + 1 in ref_gaps else c
                      for i, c in enumerate(ref))
    with open(path, 'w', encoding='utf-8') as f:
        f.write('>RefSeq\n' + ref + '\n>Other\n' + alt + '\n')


def test_orf_map_ref_nt_coords_skip_gaps(tmp_path):
    """ref 第 3、4 位是 gap：ref nt5 → 比对列 5+2=7。"""
    fa = tmp_path / 'a.fasta'
    _write_aln(fa, ref_gaps={3, 4})
    csvp = tmp_path / 'orf.csv'
    csvp.write_text('Gene.symbol,Start,Stop\nMP,1,5\nCP,8,10\n',
                    encoding='utf-8')
    out = orf_csv_to_aln_map(str(csvp), str(fa))
    assert out['ref'].startswith('RefSeq') and out['skipped'] == 0
    assert out['genes'][0] == {'name': 'MP', 'start': 1, 'end': 7}
    assert out['genes'][1] == {'name': 'CP', 'start': 10, 'end': 12}


def test_orf_map_header_case_and_ref_pick(tmp_path):
    """表头大小写不敏感；ref_id 前缀匹配选参考行。"""
    fa = tmp_path / 'a.fasta'
    _write_aln(fa)
    with open(fa, 'a', encoding='utf-8') as f:
        f.write('>NC_001440.1|xyz\nACGTACGTACGT\n')
    csvp = tmp_path / 'orf.csv'
    csvp.write_text('gene.symbol,start,stop\n3A,2,6\n', encoding='utf-8')
    out = orf_csv_to_aln_map(str(csvp), str(fa), ref_id='NC_001440.1')
    assert out['ref'].startswith('NC_001440.1')
    assert out['genes'] == [{'name': '3A', 'start': 2, 'end': 6}]


def test_orf_map_skips_out_of_range_and_junk(tmp_path):
    fa = tmp_path / 'a.fasta'
    _write_aln(fa)
    csvp = tmp_path / 'orf.csv'
    csvp.write_text('Gene.symbol,Start,Stop\n'
                    'OK,2,6\nBAD,99,120\nJUNK,x,y\n', encoding='utf-8')
    out = orf_csv_to_aln_map(str(csvp), str(fa))
    assert [g['name'] for g in out['genes']] == ['OK']
    assert out['skipped'] == 2


def test_orf_map_all_invalid_raises(tmp_path):
    fa = tmp_path / 'a.fasta'
    _write_aln(fa)
    csvp = tmp_path / 'orf.csv'
    csvp.write_text('Gene.symbol,Start,Stop\nBAD,99,120\n', encoding='utf-8')
    import pytest
    with pytest.raises(ValueError):
        orf_csv_to_aln_map(str(csvp), str(fa))
