# -*- coding: utf-8 -*-
"""前端 JS 语法静态校验（内联 script 块 + static/*.js）。

为什么需要这个测试：`webapp/templates/tools.html` 的单页 JS 有 4000+ 行，
后端测试（`_audit_contract.py` / `_it_platform.py`）只看路由与 i18n 键，
**一行都不看 JS 语法**；真浏览器测试只在跑得到那几条路径时才暴露问题。
2026-09-18 实测：一次批量替换把 `async function pgGifExport(...)` 的函数头
连删 5 行，语法已坏（后续函数体成了顶层语句），而当时**没有任何测试报错** ——
只有真打开页面点那个按钮才会发现。本测试把这个缺口补上。

做法：把每个模板的内联 `<script>` 主体（以及 `static/*.js` 全文）抽出来，
用 node `--check` 只做**语法**校验，不执行。
  · Jinja 占位符（`{{ … }}` / `{% … %}`）在检查前替换成 `null` / 空串 ——
    它们是模板变量，不是 JS 语法的一部分；
  · 带 `src=` 的 `<script>` 跳过（外部文件，另有各卡的资源存在性检查覆盖）；
  · 本机没有 node 时**跳过**而非失败（与 playwright 缺失时的处理一致）：
    这是加强项，不该让没装 node 的环境整批测试变红。

用法：
    python tests/_check_js_syntax.py
"""
from __future__ import annotations

import io
import os
import re
import shutil
import subprocess
import sys
import tempfile
from collections import Counter

for _stream in (sys.stdout, sys.stderr):
    try:
        _stream.reconfigure(errors='replace')
    except (AttributeError, OSError):
        pass

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# 托管 node 优先（隔离环境），其次系统 PATH —— 与 binary_context 的口径一致
NODE_CANDIDATES = [
    os.path.join(os.path.expanduser('~'), '.workbuddy', 'binaries', 'node',
                 'versions', '22.22.2-3', 'node.exe'),
    r'C:\Program Files\nodejs\node.exe',
]

SCRIPT_RE = re.compile(r'<script([^>]*)>(.*?)</script>', re.S)

FAILS: list = []


def check(cond, msg):
    print(('  ok   ' if cond else '  FAIL ') + msg, flush=True)
    if not cond:
        FAILS.append(msg)


def find_node():
    for p in NODE_CANDIDATES:
        if os.path.isfile(p):
            return p
    return shutil.which('node')


NON_JS_TYPES = ('text/plain', 'text/template', 'text/html',
                'application/json', 'application/ld+json')


def blocks_of(path):
    """→ [(起始行号, JS 源码, 标签)]；跳过带 src 的空块与非 JS 的 script。

    必须按 type 过滤：help.html 的英文手册用
    `<script type="text/plain" id="vp-manual-src">` 内嵌 markdown（惰性文本，
    由页面渲染器读取），把它当 JS 交给 node 会报"语法错误"——这是误报。
    只有缺省 type 或 JS MIME 的块才是真正要检查的脚本。
    """
    src = io.open(path, encoding='utf-8', errors='replace').read()
    out = []
    if path.endswith('.js'):
        return [(1, src, os.path.basename(path))]
    for m in SCRIPT_RE.finditer(src):
        attrs, body = m.group(1), m.group(2)
        if 'src=' in attrs or not body.strip():
            continue
        low = attrs.lower()
        if any(t in low for t in NON_JS_TYPES):
            continue
        line = src[:m.start()].count('\n') + 1
        out.append((line, body, 'inline script @%d' % line))
    return out


def main():
    node = find_node()
    if not node:
        print('[SKIP] 未找到 node，跳过前端 JS 语法校验')
        return 0
    print('node:', node)

    targets = []
    tdir = os.path.join(ROOT, 'webapp', 'templates')
    for name in sorted(os.listdir(tdir)):
        if name.endswith('.html'):
            targets.append(os.path.join(tdir, name))

    def _ours(name):
        """只校验**我们自己写的** JS。

        第三方压缩包（`*.min.js`、`vendor/` 下的库）不属本测试范围：
        node --check 对它们会误报（plotly.min.js 实测直接报语法错），
        而我们改不到、也不该改它们。把它们排除掉，测试的失败才等于
        「我们的代码坏了」。"""
        low = name.lower()
        return (low.endswith('.js') and not low.endswith('.min.js')
                and 'vendor' not in low)

    sdir = os.path.join(ROOT, 'webapp', 'static')
    for name in sorted(os.listdir(sdir)):
        if _ours(name):
            targets.append(os.path.join(sdir, name))

    tmp = os.path.join(tempfile.gettempdir(), '_vp_js_syntax_check.js')
    n_blocks = n_ok = 0
    print('\n[1] 内联 script 块 / static/*.js 语法')
    for path in targets:
        rel = os.path.relpath(path, ROOT).replace('\\', '/')
        for line, body, tag in blocks_of(path):
            n_blocks += 1
            js = re.sub(r'\{\{.*?\}\}', 'null', body, flags=re.S)
            js = re.sub(r'\{%.*?%\}', '', js, flags=re.S)
            # 模板里刻意留的未闭合注释/行尾续行会误报，这里只要求
            # 「能通过 node 的语法解析」，不做语义检查。
            with io.open(tmp, 'w', encoding='utf-8') as f:
                f.write(js)
            r = subprocess.run([node, '--check', tmp],
                               capture_output=True, text=True, timeout=60)
            if r.returncode == 0:
                n_ok += 1
            else:
                first = (r.stderr or '').strip().split('\n')
                detail = ' | '.join(first[:3])[:300]
                check(False, '%s 的 %s 语法错误：%s' % (rel, tag, detail))
    check(n_blocks > 0, '抽到 %d 个 JS 块（0 个说明抽取逻辑失效了）' % n_blocks)
    print('  · %d/%d 个块语法通过' % (n_ok, n_blocks))

    print('\n[2] 关键函数仍然存在（防「批量替换把函数头删了」这类事故）')
    # 这 4 个都是被跨卡复用的函数：删掉任何一个，另一张卡会静默失去功能。
    tools = io.open(os.path.join(ROOT, 'webapp', 'templates', 'tools.html'),
                    encoding='utf-8', errors='replace').read()
    for fn, why in (('function _pgArc(', '迁移弧线的唯一几何实现（网页/GIF/静态图共用）'),
                    ('function drawPgGeo(', '弧线地图渲染（系统地理卡 + 时间与地理推断卡共用）'),
                    ('async function pgGifExport(', 'GIF 导出入口'),
                    ('async function loadPhylogeoResult(', '系统地理卡结果装配')):
        check(fn in tools, '%s 仍在（%s）' % (fn, why))

    print('\n[3] tools.html 无重复 id（getElementById 只返回靠前那个）')
    # 只查 tools.html：它是唯一一个把 20+ 张卡的 DOM 塞进同一页的模板，
    # 重复 id 会表现为「后画的图跑到前一个容器里 / 另一张卡看起来没出图」。
    # 2026-09-18 实测踩过一次（drawPgGeo 通用化后内层 id 仍写死 'pgGeoPlot'）。
    # 其余模板不查：explorer.html 有 1 处历史重复（该页已停用），
    # logan.html 的重复是 JS 模板字面量（`${name}`）造成的误报。
    ids = re.findall(r'\bid="([^"]+)"', tools)
    dup = {k: v for k, v in Counter(ids).items()
           if v > 1 and '${' not in k}
    check(not dup, 'tools.html 无重复 id（重复的：%s）' % (dup or '无'))
    print('  · 共 %d 个静态 id' % len(ids))

    print('\n' + '=' * 52)
    if FAILS:
        print('FAILED: %d 项' % len(FAILS))
        for f in FAILS:
            print('  - ' + f)
        return 1
    print('前端 JS 语法检查通过 ✔')
    return 0


if __name__ == '__main__':
    sys.exit(main())
