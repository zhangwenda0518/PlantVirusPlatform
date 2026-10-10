# -*- coding: utf-8 -*-
"""t-crispr 卡 UI 契约测试：控件、i18n 键、渲染分支齐备性。

浏览器人工走查之外的可回归保障（2026-10-07 生产审查 ③）：工具卡每次
加参数/列/详情分支，这里同步断言，防止前端与引擎参数面漂移。
"""
import re

import pytest

TOOLS = 'webapp/templates/tools.html'
I18N = 'webapp/static/i18n.js'


@pytest.fixture(scope='module')
def card():
    src = open(TOOLS, encoding='utf-8').read()
    i = src.index('<section class="card" id="t-crispr">')
    j = src.index('</section>', i)
    return src[i:j]


def test_card_inputs(card):
    """八个关键控件齐备（输入/排除库/酶/扩增/探针/档位/宿主/分组）。"""
    for sel in ['id="cr_input"', 'id="cr_near"', 'id="cr_nuclease"',
                'id="cr_amp"', 'id="cr_probe"', 'id="cr_profile"',
                'id="cr_use_host"', 'id="cr_group"', 'id="cr_spacer"',
                'id="cr_topn"', 'id="cr_flank"', 'id="toolrun-crispr"']:
        assert sel in card, sel
    # 酶与探针下拉的选项值
    assert 'value="cas12a"' in card and 'value="cas13a"' in card
    assert 'value="exo"' in card and 'value="nfo"' in card


def test_card_js_wiring(card, tools_src=None):
    src = open(TOOLS, encoding='utf-8').read()
    for fn in ('runCrispr', 'loadCrisprResult', 'renderCrisprTracks',
               'crToggleDetail'):
        assert ('function ' + fn) in src, fn
    # 提交参数与引擎参数面对齐
    for key in ('nuclease:', 'probe_type:', 'use_host:', 'group_by_segment:',
                'amp_mode:', 'profile:', 'spacer_len:', 'top_n:', 'flank:'):
        assert key in src, key
    # 详情面板渲染分支：探针 / IMM / T7 前向引物
    assert 's[i].probe && s[i].probe.marked' in src
    assert 's[i].imm' in src
    assert 's[i].f_t7' in src
    # 结果表列：探针列 + route 缩写
    assert "['probe', '探针']" in src
    assert "v === 'primer3' ? 'P3'" in src


def test_card_i18n_keys(card):
    i18n = open(I18N, encoding='utf-8').read()
    keys = set(re.findall(r'data-i18n="(tk\.cr[A-Za-z]+)"', card))
    keys |= set(re.findall(r'data-i18n-ph="(tk\.cr[A-Za-z]+)"', card))
    assert keys, '卡内未找到 tk.cr* i18n 锚点'
    for k in sorted(keys):
        assert ("'%s'" % k) in i18n, k + ' 未在 i18n.js 定义'
        # zh / en 两个字典各一份（同名键应出现两次）
        assert i18n.count("'%s'" % k) >= 2, k + ' 中英不成对'


def test_card_result_zones(card):
    """结果三区：轨道容器 + 表格容器 + 详情行插入目标。"""
    assert 'id="crTracks"' in card
    assert 'id="crTable"' in card
    assert 'id="crResult"' in card
    # 卡片运行按钮防呆：空输入先 alert
    assert "if (!val('cr_input'))" in open(
        TOOLS, encoding='utf-8').read()
