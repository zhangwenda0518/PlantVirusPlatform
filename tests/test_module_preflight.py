# -*- coding: utf-8 -*-
"""运行前依赖自愈（module_deps）单测。

三层保证：
  ① 映射完整性——每个工具页 / 管道阶段都必须在 MODULE_TOOLS 里有条目，
     否则新增页面时 preflight 会静默漏检（这正是要防的回归）
  ② 检测逻辑——缺失判定、插件匹配、装不上时给指引
  ③ 自愈落盘——record_autofix / read_autofix_log 供 selfcheck 展示
"""
import json
import os

import pytest

from Virus_Platform_Core import module_deps as md


# ── ① 映射完整性（防止新增页面漏登记）────────────────────
def test_all_tool_pages_declared():
    from Virus_Platform_Core.web.tools_api import TOOL_REGISTRY
    missing = sorted(set(TOOL_REGISTRY) - set(md.MODULE_TOOLS))
    assert not missing, f'工具页未在 MODULE_TOOLS 声明依赖: {missing}'


def test_all_pipeline_stages_declared():
    from Virus_Platform_Core.pipeline import STAGE_REGISTRY
    missing = sorted(set(STAGE_REGISTRY) - set(md.MODULE_TOOLS))
    assert not missing, f'管道阶段未在 MODULE_TOOLS 声明依赖: {missing}'


def test_no_orphan_entries():
    """MODULE_TOOLS 里的键必须真实存在（防归档后留死条目）。"""
    from Virus_Platform_Core.web.tools_api import TOOL_REGISTRY
    from Virus_Platform_Core.pipeline import STAGE_REGISTRY
    valid = set(TOOL_REGISTRY) | set(STAGE_REGISTRY)
    orphan = sorted(set(md.MODULE_TOOLS) - valid)
    assert not orphan, f'MODULE_TOOLS 有已不存在的模块: {orphan}'


def test_declared_keys_are_real_tool_keys():
    """声明的工具键必须能在 detect_tools 里找到（防拼写错误）。

    特例：这几类不参与 detect_tools 探测，但在 _NO_PLUGIN_HINT 里有指引，
    或由插件以 install_to 方式还原到固定路径（模块自己拼路径找）：
      - rscript       ：phylodyn 自己找
      - gbdraw/orfipy ：随 3rd/python 运行时，detect_tools 走 Scripts 探测
      - srna-target   ：mirna_target 硬编码 3rd/tools/srna-target
    """
    from Virus_Platform_Core.config import detect_tools
    EXTRA = {'rscript', 'srna-target'}
    known = set(detect_tools()) | set(md._NO_PLUGIN_HINT) | EXTRA
    bad = []
    for mod, keys in md.MODULE_TOOLS.items():
        for k in keys:
            if k not in known:
                bad.append((mod, k))
    assert not bad, f'未知工具键: {bad}'


# ── ② 检测逻辑 ─────────────────────────────────────────
@pytest.fixture(autouse=True)
def _isolate_autofix_log(tmp_path, monkeypatch):
    """所有用例都不写真实 run/plugin_autofix.json（避免污染平台状态）。"""
    monkeypatch.setattr(md, '_AUTOFIX_LOG', str(tmp_path / '_af.json'))


class _Cfg:
    """最小 Cfg 替身：只有 tool_status / refresh_plugins。"""
    def __init__(self, have):
        self.have = have
        self.refreshed = 0

    def tool_status(self):
        return {k: (f'/x/{k}' if k in self.have else None)
                for k in self.have}

    def refresh_plugins(self):
        self.refreshed += 1


def test_missing_tools_detects_gap():
    cfg = _Cfg({'seqkit'})
    assert md.missing_tools(cfg, 'convert') == []
    # hostremoval 需要 seqkit/minibwa/kunpeng，只有 seqkit
    assert md.missing_tools(cfg, 'hostremoval') == ['minibwa', 'kunpeng']


def test_missing_tools_unknown_module():
    assert md.missing_tools(_Cfg({'seqkit'}), 'no-such-module') == []


def test_ensure_tools_noop_when_complete():
    cfg = _Cfg({'seqkit', 'fastp', 'minibwa', 'kunpeng'})
    r = md.ensure_tools(cfg, 'hostremoval')
    assert r['missing'] == [] and r['installed'] == []
    assert r['still_missing'] == []


def test_ensure_tools_reports_hint_when_no_plugin(monkeypatch):
    """插件源里没有的键（rscript/gbdraw/orfipy）要给出人工安装指引。"""
    cfg = _Cfg(set())
    monkeypatch.setattr(md.pm, 'available_from_index',
                        lambda *a, **kw: ([], None))
    r = md.ensure_tools(cfg, 'genoplot')          # 只需 gbdraw
    assert r['still_missing'] == ['gbdraw']
    assert 'pip install' in r['hints']['gbdraw']
    assert 'gbdraw' in md.format_error('genoplot', r)


def test_ensure_tools_installs_from_index(monkeypatch, tmp_path):
    """索引里有能提供缺失键的插件 → 调用安装并重扫。"""
    installed = []

    def fake_index(cfg, force=False):
        return ([{'name': 'blast', 'provides': {'blastn': 'bin/blastn.exe',
                                                'makeblastdb': 'bin/mk.exe'},
                  'urls': ['http://x/blast.zip'], 'sha256': 'abc',
                  'size': 1000, 'installed': False, 'needed': True}], None)

    def fake_install(urls, expected_sha256=None):
        installed.append((urls, expected_sha256))
        return {'name': 'blast'}

    monkeypatch.setattr(md.pm, 'available_from_index', fake_index)
    monkeypatch.setattr(md.pm, 'install_urls', fake_install)
    cfg = _Cfg({'orfipy'})
    # orfa 需要 blastn + makeblastdb，两个都由 blast 插件提供
    r = md.ensure_tools(cfg, 'orfa')
    assert r['installed'] == ['blast']
    assert installed and installed[0][1] == 'abc'   # sha256 透传
    assert cfg.refreshed == 1                        # 装完重扫


def test_ensure_tools_install_failure_is_reported(monkeypatch):
    """安装失败不抛异常，如实返回 still_missing + 指引。"""
    def fake_index(cfg, force=False):
        return ([{'name': 'blast', 'provides': {'blastn': 'b'},
                  'urls': ['http://x/b.zip'], 'sha256': 'x', 'size': 1,
                  'installed': False, 'needed': True}], None)

    def boom(urls, expected_sha256=None):
        raise OSError('网络不通')

    monkeypatch.setattr(md.pm, 'available_from_index', fake_index)
    monkeypatch.setattr(md.pm, 'install_urls', boom)
    r = md.ensure_tools(_Cfg({'makeblastdb'}), 'orfa')
    assert r['installed'] == []
    assert r['still_missing'] == ['blastn']
    assert 'blastn' in r['hints']


def test_ensure_tools_auto_install_false_only_checks(monkeypatch):
    called = []
    monkeypatch.setattr(md.pm, 'available_from_index',
                        lambda *a, **kw: called.append(1) or ([], None))
    r = md.ensure_tools(_Cfg(set()), 'convert', auto_install=False)
    assert r['still_missing'] == ['seqkit']
    assert called == []            # 只检查不下载


def test_format_error_lists_all_parts():
    r = {'missing': ['a', 'b'], 'installed': ['p1'],
         'still_missing': ['a', 'b'],
         'hints': {'a': '装 A', 'b': '装 B'}}
    txt = md.format_error('demo', r)
    assert 'demo' in txt and '装 A' in txt and '装 B' in txt
    assert 'p1' in txt


# ── ③ 自愈记录落盘 ─────────────────────────────────────
def test_autofix_log_roundtrip(tmp_path, monkeypatch):
    monkeypatch.setattr(md, '_AUTOFIX_LOG', str(tmp_path / 'af.json'))
    assert md.read_autofix_log() == {}
    md.record_autofix('orfa', ['blast'], [])
    log = md.read_autofix_log()
    assert log['orfa']['installed'] == ['blast']
    assert 'at' in log['orfa']
    # 空结果不写盘
    md.record_autofix('convert', [], [])
    assert 'convert' not in md.read_autofix_log()


def test_autofix_log_survives_corrupt_file(tmp_path, monkeypatch):
    p = tmp_path / 'af.json'
    p.write_text('not json', encoding='utf-8')
    monkeypatch.setattr(md, '_AUTOFIX_LOG', str(p))
    assert md.read_autofix_log() == {}
    md.record_autofix('x', ['y'], [])
    assert md.read_autofix_log()['x']['installed'] == ['y']
