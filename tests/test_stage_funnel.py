# -*- coding: utf-8 -*-
"""stage_funnel 流程漏斗三件套单测（fixture = tests/fixtures/sample_demo）。

守护语义（MMPV 蓝本退化守卫）：
- 各阶段产物存在 → 逐层提取数字与来源；
- 产物缺失 → 该层 value=None 跳过，绝不报 0；
- 全缺 → build_funnel 返回 None；
- 出图：双图落盘、TSV 行与有效层一致；空目录不产图。
"""
import csv
import json
import os

import pytest

from Virus_Platform_Core import stage_funnel as sf

PLATFORM = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIXTURE = os.path.join(PLATFORM, 'tests', 'fixtures', 'sample_demo')


def test_fixture_funnel_rows():
    res = sf.build_funnel(FIXTURE)
    assert res is not None
    by_key = {r['key']: r for r in res['rows']}
    # fixture 覆盖的核心层必须提取到数字
    for k, expect in (('input_reads', 7200), ('qc_reads', 7150),
                      ('host_free_reads', 7112), ('viral_reads', 1834),
                      ('contigs', 24), ('viral_contigs', 2),
                      ('verified', 2), ('consensus_refs', 2),
                      ('orf_genes', 14)):
        assert by_key[k]['value'] == expect, (k, by_key[k])
        assert by_key[k]['src'], f'{k} 缺数据源标注'
    # 每层 note 不含科学计数残留
    for r in res['rows']:
        assert 'e+' not in str(r['note'])


def test_kpi_cards_from_side_stages():
    kpi = dict((k, v) for k, v, _n in sf.kpi_stats(FIXTURE))
    assert '2 个物种' in kpi.get('已知病毒检出', '')
    assert kpi.get('变异位点') == '27'
    assert '6 对' in kpi.get('引物', '')
    assert '基因组图' in kpi


def test_tsv_written(tmp_path):
    res = sf.build_funnel(FIXTURE, out_dir=str(tmp_path))
    tsv = tmp_path / 'stage_funnel.tsv'
    assert tsv.is_file()
    with open(tsv, encoding='utf-8') as f:
        rows = list(csv.DictReader(f, delimiter='\t'))
    valid = [r for r in rows if r['value']]
    assert len(valid) == 9           # fixture 九层全有
    assert valid[0]['stage'] == '输入 reads'


def test_missing_everything_returns_none(tmp_path):
    empty = tmp_path / 'empty_sample'
    empty.mkdir()
    assert sf.build_funnel(str(empty)) is None


def test_partial_missing_skips_layer(tmp_path):
    d = tmp_path / 'partial'
    (d / '01_host_removal').mkdir(parents=True)
    json.dump({'total_pairs': 100, 'kept_pairs': 90,
               'host_ratio': 0.1, 'single_end': False},
              open(d / '01_host_removal' / 'stats.json', 'w',
                   encoding='utf-8'))
    res = sf.build_funnel(str(d))
    assert res is not None
    by_key = {r['key']: r for r in res['rows']}
    assert by_key['host_free_reads']['value'] == 180   # 对×2
    assert by_key['contigs']['value'] is None          # 缺产物 = None 而非 0
    # 无 fastp 时输入层回退到宿主去除 stats：单位必须与去宿主层一致
    # （双端 对×2 折条；2026-10-13 真实全流程实抓混单位 → 留存率虚高 2 倍）
    assert by_key['input_reads']['value'] == 200
    assert '对×2' in by_key['input_reads']['note']


def test_single_end_units_not_doubled(tmp_path):
    """单端 stats：kept_pairs 本就是单条 reads（host_removal kept1=R1 记录数），
    输入/去宿主两层都不得再 ×2。"""
    d = tmp_path / 'se'
    (d / '01_host_removal').mkdir(parents=True)
    json.dump({'total_pairs': 500, 'kept_pairs': 450,
               'host_ratio': 0.1, 'single_end': True},
              open(d / '01_host_removal' / 'stats.json', 'w',
                   encoding='utf-8'))
    by_key = {r['key']: r for r in sf.build_funnel(str(d))['rows']}
    assert by_key['input_reads']['value'] == 500
    assert by_key['input_reads']['note'] == ''
    assert by_key['host_free_reads']['value'] == 450


def test_plots_render_and_qa(tmp_path):
    res = sf.render(FIXTURE, out_dir=str(tmp_path))
    assert res is not None
    png, pdf = tmp_path / 'stage_funnel.png', tmp_path / 'stage_funnel.pdf'
    rpng = tmp_path / 'stage_funnel_retention.png'
    assert png.is_file() and pdf.is_file() and rpng.is_file()
    # 图形体检（MMPV 蓝本 §1.6）：文字两两无重叠、无越界
    import sys
    dev_tools = os.path.join(PLATFORM, 'dev_tools')
    if dev_tools not in sys.path:
        sys.path.insert(0, dev_tools)
    import check_figures
    with check_figures.fig_quality_guard('test_funnel') as rep:
        sf.plot_funnel(res, str(tmp_path / 'qa_funnel'))
        sf.plot_retention(res, str(tmp_path / 'qa_retention'))
    assert not rep.fails, rep.summary()
