# -*- coding: utf-8 -*-
"""结果图整理（内置 svguest）回归守护。

锁四件事：
1. 内置单文件存在、来源/许可齐全、`?open=` 补丁仍在且只认同源；
2. `/svguest/` 路由可访问并返回 HTML；
3. 结果中心导航有入口；
4. 主要绘图助手确实产出 png/svg/pdf（SVG 是二次排版的输入）。
"""
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

VENDOR = os.path.join(ROOT, 'webapp', 'static', 'vendor', 'svguest')


def _read(name):
    with open(os.path.join(VENDOR, name), encoding='utf-8') as f:
        return f.read()


def test_vendored_files_and_patch():
    app_html = _read('app.html')
    assert len(app_html) > 400_000
    # 上游标记
    assert '<title>svguest</title>' in app_html
    # 平台补丁：?open= + 同源限制 + 复用上游 readFile
    assert 'VirusPlatform 集成补丁' in app_html
    assert "get('open')" in app_html
    assert 'u.origin !== location.origin' in app_html
    assert 'readFile(f, false)' in app_html
    # 许可与来源说明必须随目录保留（GPL-3.0）
    lic = _read('LICENSE')
    assert 'GNU GENERAL PUBLIC LICENSE' in lic
    assert 'Version 3' in lic
    src = _read('SOURCE.md')
    assert 'github.com/JieqiangHe/svguest' in src
    assert 'GPL-3.0' in src


def test_svguest_route_serves_patched_html():
    import app as platform_app
    client = platform_app.app.test_client()
    r = client.get('/svguest/')
    assert r.status_code == 200
    assert 'text/html' in r.headers.get('Content-Type', '')
    assert b'VirusPlatform' in r.data
    assert b"get('open')" in r.data
    assert r.headers.get('Cache-Control') == 'no-cache'


def test_result_center_nav_has_entry():
    from Virus_Platform_Core.web.pages import NAV_GROUPS, _PATH_TO_GROUP
    result = next(g for g in NAV_GROUPS if g['id'] == 'result')
    assert any(it.get('href') == '/svguest/' for it in result['items'])
    assert _PATH_TO_GROUP['/svguest/'] == 'result'


def test_figure_helper_writes_three_formats(tmp_path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    from Virus_Platform_Core.known_virus_suite import kv_fig

    fig, ax = plt.subplots(figsize=(2, 2))
    ax.plot([0, 1], [0, 1])
    kv_fig._save(fig, str(tmp_path), 'demo_fig')

    names = set(os.listdir(tmp_path))
    assert {'demo_fig.png', 'demo_fig.svg', 'demo_fig.pdf'} <= names


def test_variant_default_formats_include_svg():
    from Virus_Platform_Core.known_virus_suite import (kv_variant_evo,
                                                       kv_variant_plot)
    for mod in (kv_variant_plot, kv_variant_evo):
        assert 'svg' in mod.FIG_FORMATS
        assert mod._save.__defaults__[-1] == mod.FIG_FORMATS


def test_frontend_helper_present():
    with open(os.path.join(ROOT, 'webapp', 'static', 'app.js'),
              encoding='utf-8') as f:
        js = f.read()
    assert 'function svgEditHref' in js
    assert "return '/svguest/?open=' + encodeURIComponent(href)" in js
