# -*- coding: utf-8 -*-
"""SVG 良构守卫（2026-10-04 审查 M3 立；2026-10-07 审查 P1-4 扩面）。

背景
----
2026-10-04 审查实测：`examples/results/synteny/compare.svg` 曾有 18 处
`</path` 闭合缺 `>`，矢量交付物损坏（浏览器容错渲染成残图，肉眼不易察觉）。
`utils.verify_svg` 已实现幂等自修，但**只有 `dfv_plot.py` 挂了调用**，
`gbdraw_plot.py`（外部工具产物）与 `genome_diag.py`（手写 SVG）两个出口
都没有校验，而 `viz.py` 会把 gbdraw 的 SVG 静默读入并内嵌进 `report.html`
—— 残图就这样进了交付报告。

判据（五条）
-----------
A. **覆盖范围**：`examples/` 与 `run/` 下全部 `*.svg` 必须良构（只读解析，
   **不修改任何文件**）；`KNOWN_BAD` 白名单内除外（故意残缺的回归样本）。
B. **灵敏度自测（对照臂）**：临时目录造一张残缺 SVG（`</path` 缺 `>`）与一张
   良构 SVG，用**同一个** `_is_wellformed()` 判据 —— 残缺必须被判坏、
   良构必须不被判坏。否则本守卫可能恒过（判据本身失灵）。
C. **修复路径可用**：残缺 fixture 的副本交给 `utils.verify_svg`，必须返回
   True 且修复后文件可解析（钉住幂等自修机制，D3 类回归）。
D. **白名单不腐烂**：`KNOWN_BAD` 每条必须**仍然**非良构；若已变成良构，
   说明条目过期 —— 报 FAIL 提示清理，防止白名单无限吞掉真问题。
E. **出口覆盖（AST 判据）**：`dfv_plot.py` / `gbdraw_plot.py` / `genome_diag.py`
   的产物出口、以及 `viz.py` 的内嵌前必须调用 `verify_svg`；且 `viz.py` 中
   对同一变量（`svg_p`）的 `verify_svg` 调用必须**早于** `safe_open` 读取。
   用 AST 而非文本搜索 —— 注释里就写着 `verify_svg` 这几个字。

用法: python tests/_check_svg_wellformed.py
退出码: 0 = 通过；1 = 失败。
"""
import ast
import os
import shutil
import sys
import tempfile
import xml.etree.ElementTree as ET

for _s in (sys.stdout, sys.stderr):
    try:
        _s.reconfigure(errors='replace')
    except (AttributeError, OSError):
        pass

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

# 故意残缺的样本（非缺陷）：留着是为了随时能回归 verify_svg 的修复路径。
# 精确到相对路径（正斜杠，与 _rel() 口径一致）—— 用目录通配会吞掉真正的新残缺。
KNOWN_BAD = {
    'run/_tmp/bad.svg':          '最小残缺 fixture（判据灵敏度样本）',
    'run/_audit_figs/synteny.svg': '2026-10-04 审查复现图（</path 缺 >）',
}

SCAN_ROOTS = ('examples', 'run')

# E. 出口覆盖：模块 → 该模块必须出现 verify_svg 调用
EXPORT_MODULES = (
    'Virus_Platform_Core/dfv_plot.py',
    'Virus_Platform_Core/gbdraw_plot.py',
    'Virus_Platform_Core/genome_diag.py',
    'Virus_Platform_Core/viz.py',
)

# 与历史缺陷同形：**开标签合法**、闭合标签 `</path` 缺 `>`（verify_svg 可自修）。
# 注意不能写成 `<path d="M0"</path` —— 那样开标签也没闭合，属不可修，
# 测不出"修复路径可用"（本守卫第一版就踩了这个坑）。
BROKEN_SVG = '<svg><g><path d="M0 0 L10 10"></path\n</g></svg>'
GOOD_SVG = ('<?xml version="1.0"?>\n'
            '<svg xmlns="http://www.w3.org/2000/svg" width="10" height="10">'
            '<g><path d="M0 0 L10 10"/></g></svg>\n')

FAILS = []


def _check(cond, msg):
    print(('  ok   ' if cond else '  FAIL ') + msg)
    if not cond:
        FAILS.append(msg)


def _is_wellformed(path):
    """返回 (是否良构, 错误串)。只读，不改文件。"""
    try:
        ET.parse(path)
        return True, ''
    except ET.ParseError as e:
        return False, str(e)


def _iter_svgs():
    for base in SCAN_ROOTS:
        root = os.path.join(ROOT, base)
        if not os.path.isdir(root):
            continue
        for dp, _dirs, files in os.walk(root):
            for f in files:
                if f.endswith('.svg'):
                    yield os.path.join(dp, f)


def _rel(p):
    return os.path.relpath(p, ROOT).replace('\\', '/')


def _calls(tree):
    """产出 (函数名, 首个实参名, lineno) —— 只取 `f(x)` 形式。"""
    out = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        f = node.func
        name = f.attr if isinstance(f, ast.Attribute) else getattr(f, 'id', None)
        if not name:
            continue
        arg = None
        if node.args and isinstance(node.args[0], ast.Name):
            arg = node.args[0].id
        out.append((name, arg, node.lineno))
    return out


def part_a_scan():
    print('[A] 覆盖范围：examples/ + run/ 全部 SVG 良构（只读，不改文件）')
    n = 0
    bad = []
    for p in _iter_svgs():
        n += 1
        ok, err = _is_wellformed(p)
        if ok:
            continue
        if _rel(p) in KNOWN_BAD:
            continue
        bad.append((_rel(p), err))
    print(f'         扫描 {n} 个 SVG（白名单 {len(KNOWN_BAD)} 个故意残缺）')
    for p, e in bad:
        _check(False, f'{p}: {e}')
    _check(not bad, f'非白名单残缺 SVG = {len(bad)} 个（应 0）')
    return n


def part_b_sensitivity(tmpdir):
    print('\n[B] 灵敏度自测（对照臂）：判据必须能区分残缺与良构')
    pb = os.path.join(tmpdir, 'broken.svg')
    pg = os.path.join(tmpdir, 'good.svg')
    with open(pb, 'w', encoding='utf-8', newline='\n') as f:
        f.write(BROKEN_SVG)
    with open(pg, 'w', encoding='utf-8', newline='\n') as f:
        f.write(GOOD_SVG)
    ok_b, err_b = _is_wellformed(pb)
    ok_g, _err_g = _is_wellformed(pg)
    _check(not ok_b, f'残缺 fixture 被判坏（err={err_b or "无"}）')
    _check(ok_g, '良构 fixture 未被误判')
    # 若把残缺 fixture 也判成良构 → 判据失灵，本守卫会恒过
    return not ok_b and ok_g


def part_c_repair(tmpdir):
    print('\n[C] 修复路径：verify_svg 幂等自修仍然可用')
    from Virus_Platform_Core.utils import verify_svg
    pc = os.path.join(tmpdir, 'tofix.svg')
    with open(pc, 'w', encoding='utf-8', newline='\n') as f:
        f.write(BROKEN_SVG)
    try:
        r = verify_svg(pc)
    except Exception as e:                            # noqa: BLE001
        _check(False, f'verify_svg 抛异常：{type(e).__name__}: {e}')
        return
    _check(r is True, f'verify_svg 返回 True（实际 {r!r}）')
    ok, err = _is_wellformed(pc)
    _check(ok, f'修复后文件可解析（err={err or "无"}）')
    # 幂等：再修一次仍是良构
    verify_svg(pc)
    ok2, _ = _is_wellformed(pc)
    _check(ok2, '幂等：二次调用后仍良构')


def part_d_whitelist():
    print('\n[D] 白名单不腐烂：KNOWN_BAD 每条必须仍然残缺')
    stale = []
    n_gone = 0
    for rel, why in sorted(KNOWN_BAD.items()):
        p = os.path.join(ROOT, rel)
        if not os.path.isfile(p):
            # run/_tmp/ 是会被清理的临时目录，样本不在属正常，不算失败
            n_gone += 1
            print(f'  note  {rel} 不存在（{why}）—— 无待抑制项，跳过')
            continue
        ok, _err = _is_wellformed(p)
        if ok:
            stale.append(f'{rel}（已变良构，条目过期 · {why}）')
    for s in stale:
        _check(False, f'白名单过期：{s}')
    _check(not stale,
           f'{len(KNOWN_BAD) - n_gone}/{len(KNOWN_BAD)} 条白名单仍然残缺'
           f'（{n_gone} 条样本缺失已跳过）')


def part_e_exports():
    print('\n[E] 出口覆盖（AST 判据）：产物出口与内嵌前必须调 verify_svg')
    for rel in EXPORT_MODULES:
        p = os.path.join(ROOT, rel)
        if not os.path.isfile(p):
            _check(False, f'{rel} 不存在')
            continue
        src = open(p, encoding='utf-8').read()
        tree = ast.parse(src)
        cs = _calls(tree)
        n_v = sum(1 for n, _a, _l in cs if n == 'verify_svg')
        _check(n_v > 0, f'{rel}: verify_svg 调用 {n_v} 处（应 >0）')

    # viz.py 顺序：verify_svg(svg_p) 必须早于 safe_open(svg_p)
    rel = 'Virus_Platform_Core/viz.py'
    p = os.path.join(ROOT, rel)
    if os.path.isfile(p):
        src = open(p, encoding='utf-8').read()
        cs = _calls(ast.parse(src))
        v = [l for n, a, l in cs if n == 'verify_svg' and a == 'svg_p']
        s = [l for n, a, l in cs if n == 'safe_open' and a == 'svg_p']
        if not v:
            _check(False, f'{rel}: 找不到 verify_svg(svg_p) 调用')
        elif not s:
            _check(False, f'{rel}: 找不到 safe_open(svg_p) 调用')
        else:
            _check(min(v) < min(s),
                   f'{rel}: verify_svg(svg_p)@{min(v)} 早于 '
                   f'safe_open(svg_p)@{min(s)}（先修后读）')


def main():
    print('SVG 良构守卫（2026-10-04 M3 / 2026-10-07 P1-4）')
    print('=' * 62)
    part_a_scan()
    tmpdir = tempfile.mkdtemp(prefix='vp_svgwf_')
    try:
        part_b_sensitivity(tmpdir)
        part_c_repair(tmpdir)
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)
    part_d_whitelist()
    part_e_exports()

    print()
    if FAILS:
        print(f'✘ SVG 良构守卫失败（{len(FAILS)} 条）')
        for f in FAILS:
            print(f'   · {f}')
        return 1
    print('✓ SVG 良构守卫通过（范围含 run/ 产物；判据经灵敏度自测）')
    return 0


if __name__ == '__main__':
    sys.exit(main())
