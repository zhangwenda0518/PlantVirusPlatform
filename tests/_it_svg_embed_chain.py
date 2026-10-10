# -*- coding: utf-8 -*-
"""P1-4 端到端链路回归：残缺 SVG 不得进入 report.html。

背景
----
`viz.build_report` 会把 `09_genome_plots/*.svg`（gbdraw 外部工具产物）读进来
内嵌进 `report.html`。gbdraw 曾产出 `</path` 缺 `>` 的残缺 SVG
（2026-10-04 审查实测 18 处），浏览器容错渲染成残图 —— 肉眼几乎看不出，
但矢量交付物已损坏。出口侧（`gbdraw_plot` / `genome_diag` / `dfv_plot`）挂了
`verify_svg` 后，这里再补消费侧兜底：**先修后读，修不动就跳过嵌入**。

三条臂
------
A. **修得动**：闭合标签缺 `>`（开标签合法）→ 被就地修复、以良构片段内嵌。
B. **修不动**：开标签也没闭合 → 不得内嵌（对照臂：直接读原文证明它确实非良构，
   否则本测试测不出东西）。
C. **本来就良构**：原样内嵌，且报告内**每一个** `<svg>...</svg>` 片段都可解析。

用法: python tests/_it_svg_embed_chain.py
退出码: 0 = 通过；1 = 失败。
"""
import json
import os
import shutil
import sys
import xml.etree.ElementTree as ET

for _s in (sys.stdout, sys.stderr):
    try:
        _s.reconfigure(errors='replace')
    except (AttributeError, OSError):
        pass

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from Virus_Platform_Core import viz  # noqa: E402

PROBE = os.path.join(ROOT, 'run', '_it_svg_chain')

# A. 与历史缺陷同形：开标签合法，仅闭合 </path 缺 '>'
FIXABLE = ('<?xml version="1.0"?>\n'
           '<svg xmlns="http://www.w3.org/2000/svg" width="120" height="80">'
           '<g id="gfix"><path d="M 0 0 L 10 10"></path\n</g></svg>\n')
# B. 开标签也没闭合 —— verify_svg 的正则补不回来
UNFIXABLE = ('<?xml version="1.0"?>\n'
             '<svg xmlns="http://www.w3.org/2000/svg" width="120" height="80">'
             '<g id="gbad"><path d="M 0 0 L 10 10"</path\n</g></svg>\n')
# C. 本来就良构
GOOD = ('<?xml version="1.0"?>\n'
        '<svg xmlns="http://www.w3.org/2000/svg" width="120" height="80">'
        '<g id="gok"><path d="M 0 0 L 10 10"/></g></svg>\n')

FAILS = []


def _check(cond, msg):
    print(('  ok   ' if cond else '  FAIL ') + msg)
    if not cond:
        FAILS.append(msg)


def _wellformed(text):
    try:
        ET.fromstring(text)
        return True, ''
    except ET.ParseError as e:
        return False, str(e)


def _svg_fragments(html):
    """从报告 HTML 里抽出每个 <svg>…</svg> 片段。"""
    out, i = [], 0
    while True:
        a = html.find('<svg', i)
        if a < 0:
            return out
        b = html.find('</svg>', a)
        if b < 0:
            return out
        out.append(html[a:b + 6])
        i = b + 6


def main():
    print('P1-4 链路回归：残缺 SVG 不得进入 report.html')
    print('=' * 62)
    shutil.rmtree(PROBE, ignore_errors=True)
    gdir = os.path.join(PROBE, '09_genome_plots')
    os.makedirs(gdir)
    try:
        for name, text in (('fixable.svg', FIXABLE), ('unfixable.svg', UNFIXABLE),
                           ('good.svg', GOOD)):
            with open(os.path.join(gdir, name), 'w', encoding='utf-8',
                      newline='\n') as f:
                f.write(text)
        with open(os.path.join(gdir, 'summary.json'), 'w', encoding='utf-8') as f:
            json.dump({'stage': 'gbdraw', 'n_seqs': 1,
                       'plots': ['fixable.svg', 'unfixable.svg', 'good.svg']}, f)

        # ---------- 对照臂前置：两个坏文件必须确实非良构 ----------
        print('[对照] 两个坏样本读原文确实非良构（否则本测试恒过）')
        for name, text in (('fixable.svg', FIXABLE),
                           ('unfixable.svg', UNFIXABLE)):
            ok, err = _wellformed(text)
            _check(not ok, f'{name} 读原文非良构（{err or "可解析"}）')

        # ---------- 被测 ----------
        rp = viz.build_report(PROBE)
        html = open(rp, encoding='utf-8').read()
        print(f'\n[被测] build_report → {os.path.relpath(rp, ROOT)}')

        # A/B: 残缺字节串不得出现在报告里
        _check('</path\n' not in html,
               '报告内不含残缺字节串 "</path\\n"')

        frags = _svg_fragments(html)
        n_bad = sum(1 for fr in frags if not _wellformed(fr)[0])
        print(f'         报告内 SVG 片段 {len(frags)} 个，非良构 {n_bad} 个')
        _check(n_bad == 0, f'非良构片段 = {n_bad} 个（应 0）')

        # A: 修得动的那张必须进来了（且带 gfix 标记 = 确实是它的内容）
        _check(any('gfix' in fr for fr in frags),
               '修得动的 SVG 已修复并内嵌（标记 gfix 可见）')
        # C: 本来就好的那张也在
        _check(any('gok' in fr for fr in frags),
               '本来良构的 SVG 正常内嵌（标记 gok 可见）')
        # B: 修不动的那张不得进来
        _check(not any('gbad' in fr for fr in frags),
               '修不动的 SVG 已被跳过、未内嵌（标记 gbad 不可见）')

        # A: 落盘文件被就地修复
        p_fix = os.path.join(gdir, 'fixable.svg')
        ok_fix, err_fix = _wellformed(open(p_fix, encoding='utf-8').read())
        _check(ok_fix, f'fixable.svg 落盘文件已被就地修复（{err_fix or "可解析"}）')
        # B: 落盘文件保持原样（不写入半成品）
        p_un = os.path.join(gdir, 'unfixable.svg')
        _check(open(p_un, encoding='utf-8').read() == UNFIXABLE,
               'unfixable.svg 落盘内容未被改动（不留半成品）')
    finally:
        shutil.rmtree(PROBE, ignore_errors=True)

    print()
    if FAILS:
        print(f'✘ P1-4 链路回归失败（{len(FAILS)} 条）')
        for f in FAILS:
            print(f'   · {f}')
        return 1
    print('✓ P1-4 链路回归通过：报告内 SVG 全部良构，修不动的被跳过')
    return 0


if __name__ == '__main__':
    sys.exit(main())
