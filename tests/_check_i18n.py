# -*- coding: utf-8 -*-
"""i18n 缺键校验（Python 包装）：把 tests/_verify_i18n.js 挂进批量回归。

为什么要这层包装
----------------
`_verify_i18n.js` 一直不在 `tests/_run_all.py` 的清单里（批量入口只跑 .py），
所以它**从来没被执行过** —— 2026-09-15 记的缺陷后来修了、fixture 没跟着改，
脚本就变成永久红灯也没人发现。这个包装把「模板锚点 ↔ 字典键」的检查
纳入常态回归。

校验内容：
  1) 静态扫描（纯 Python，不依赖 node）：webapp/static/*.js（跳过 min/bak）与
     webapp/templates/*.html 内嵌 <script> 中，alert / confirm / prompt 的
     第一参数不得含「裸 CJK 字面量」——必须走 t('key')。
     · console.* 调用整体剔除（调试文案不算回归）；
     · t('key', '中文回退') 的参数会被剔除，故「已走 i18n、仅带中文回退」合法。
  2) 实现见 _verify_i18n.js（需 node）：扫 webapp/templates/*.html 里所有
     data-i18n / -ph / -title 的值，逐个喂给真实的 applyI18n，断言不会被写成
     key 字面量；断言中英成对（只缺一种语言同样是坏文案）；并扫源码文本断言
     字典**无重复键**（JS 对象字面量后者静默覆盖前者，解析后看不见）。

本机没有 node 时只跳过第 2 项（第 1 项是纯 Python，照跑）：
这是加强项，不该让没装 node 的环境整批测试变红。

用法: python tests/_check_i18n.py
退出码: 0 = 通过或（node 缺失时）部分跳过；1 = 裸 CJK 对话框文案 / 缺键 / 中英不成对 / 重复键。
"""
import glob
import os
import re
import shutil
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# 托管 node 优先（隔离环境），其次系统 PATH —— 与 _check_js_syntax.py 口径一致
NODE_CANDIDATES = [
    os.path.join(os.path.expanduser('~'), '.workbuddy-ai', 'binaries', 'node',
                 'versions', '22.22.2-3', 'node.exe'),
    os.path.join(os.path.expanduser('~'), '.workbuddy', 'binaries', 'node',
                 'versions', '22.22.2-3', 'node.exe'),
    r'C:\Program Files\nodejs\node.exe',
]


def find_node():
    for p in NODE_CANDIDATES:
        if os.path.isfile(p):
            return p
    return shutil.which('node')


# ---------------------------------------------------------------- 静态扫描
# 2026-10-03 新增：alert/confirm/prompt 文案不得再出现「裸 CJK 字面量」。
# 扫描 webapp/static/*.js（跳过 min/bak）与 webapp/templates/*.html 的内嵌
# <script>（跳过 src= 外链）；console.* 调用整体剔除后再查——调试用中文不算
# 文案回归。t('key', '中文回退') 的字符串参数会被剔除，因此「已走 i18n、
# 仅带中文回退」的写法合法，只有裸字符串才是 FAIL。
CJK_RE = re.compile(r'[\u4e00-\u9fff]')
DIALOG_RE = re.compile(r'\b(alert|confirm|prompt)\s*\(')
CONSOLE_RE = re.compile(r'\bconsole\s*\.\s*\w+\s*\(')
T_CALL_RE = re.compile(r"\bt\(\s*'((?:[^'\\]|\\.)*)'")
T_ARG2_RE = re.compile(r"\s*,\s*'((?:[^'\\]|\\.)*)'")
SCRIPT_RE = re.compile(r'<script(?![^>]*\bsrc=)[^>]*>(.*?)</script>', re.S)


def _arg_end(src, i):
    """i 指向 '(' 之后首个字符；返回与之配对的 ')' 的下标（无则 len(src)）。"""
    depth, inq, esc = 1, None, False
    while i < len(src):
        ch = src[i]
        if inq:
            if esc:
                esc = False
            elif ch == '\\':
                esc = True
            elif ch == inq:
                inq = None
        else:
            if ch in '\'"`':
                inq = ch
            elif ch == '(':
                depth += 1
            elif ch == ')':
                depth -= 1
                if depth == 0:
                    return i
        i += 1
    return len(src)


def _strip_t_calls(s):
    """剔除 t('key') / t('key', '中文回退') 子调用，只留其余实参文本。"""
    out, i = [], 0
    while i < len(s):
        m = T_CALL_RE.match(s, i)
        if m:
            j = m.end()
            m2 = T_ARG2_RE.match(s, j)
            if m2:
                j = m2.end()
            if j < len(s) and s[j] == ')':
                i = j + 1
                continue
        out.append(s[i])
        i += 1
    return ''.join(out)


def _strip_console_calls(src):
    out, i = [], 0
    while i < len(src):
        m = CONSOLE_RE.match(src, i)
        if m:
            i = _arg_end(src, m.end()) + 1
            continue
        out.append(src[i])
        i += 1
    return ''.join(out)


def _scan_js(rel, src, base_line, problems):
    src = _strip_console_calls(src)
    for m in DIALOG_RE.finditer(src):
        arg = src[m.end():_arg_end(src, m.end())]
        if CJK_RE.search(_strip_t_calls(arg)):
            ln = base_line + src[:m.start()].count('\n')
            problems.append((rel, ln, m.group(1),
                             arg.strip()[:60].replace('\n', ' ')))


def scan_cjk_dialogs(root):
    """返回 [(相对路径, 行号, 调用, 实参片段), ...]；模板只扫内嵌 script。"""
    problems = []
    for path in sorted(glob.glob(os.path.join(root, 'webapp', 'static', '*.js'))):
        name = os.path.basename(path)
        if name.endswith('.min.js') or '.bak' in name:
            continue
        with open(path, encoding='utf-8', errors='replace') as fh:
            _scan_js(os.path.relpath(path, root), fh.read(), 1, problems)
    for path in sorted(glob.glob(os.path.join(root, 'webapp', 'templates', '*.html'))):
        if '.bak' in os.path.basename(path):
            continue
        with open(path, encoding='utf-8', errors='replace') as fh:
            src = fh.read()
        rel = os.path.relpath(path, root)
        for bm in SCRIPT_RE.finditer(src):
            base = src[:bm.start(1)].count('\n') + 1
            _scan_js(rel, bm.group(1), base, problems)
    return problems


def main():
    # 第 1 项：裸 CJK 对话框文案静态扫描（纯 Python，node 缺失也照跑）
    problems = scan_cjk_dialogs(ROOT)
    if problems:
        for rel, ln, fn, seg in problems:
            print(f'  ✘ {rel}:{ln}  {fn}({seg}…) 含裸 CJK 字面量')
        print('\n[FAIL] alert/confirm/prompt 第一参数出现裸中文字面量——'
              '请抽成 i18n key 并改用 t(\'key\')')
        return 1
    print('裸 CJK 对话框文案 0 处 ✔（alert/confirm/prompt 已全部走 t()）\n')

    node = find_node()
    if not node:
        print('[SKIP] 未找到 node，跳过 i18n 缺键校验')
        return 0
    print('node:', node)

    script = os.path.join(ROOT, 'tests', '_verify_i18n.js')
    if not os.path.isfile(script):
        print(f'[FAIL] 找不到 {script}')
        return 1

    r = subprocess.run([node, script], cwd=ROOT, capture_output=True,
                       text=True, encoding='utf-8', errors='replace')
    out = (r.stdout or '') + (r.stderr or '')
    print(out.rstrip())

    if r.returncode != 0:
        print('\n[FAIL] i18n 校验未通过 —— 锚点缺字典键 / 中英不成对 / 字典有重复键')
        print('       看上面的 ✘ 清单；加 --verbose 可看完整（本包装固定不带）')
        return 1
    print('\ni18n 校验通过 ✔')
    return 0


if __name__ == '__main__':
    for _s in (sys.stdout, sys.stderr):
        try:
            _s.reconfigure(errors='replace')
        except (AttributeError, OSError):
            pass
    sys.exit(main())
