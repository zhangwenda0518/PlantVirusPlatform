# -*- coding: utf-8 -*-
"""SDT 分布图边界守卫（2026-10-07 审查 P1-5）。

背景
----
`sdt_exact.plot_distribution` 在**全同输入**（所有两两 identity 都恰好 100，
例如输入含 ≥2 条完全相同的序列）下：
`min_v = floor(nanmin(vals)) = 100` → `np.arange(100, 101, 1)` 只有 1 个边界
→ `histogram` 产出 0 个 bin → `max(空数组)` 抛
`ValueError: max() iterable argument is empty`。
调用处（`run_sdt_exact`）无 try/except，而此刻 CSV 与热图**已落盘**、
`sdt_matrix.json` 未写 —— 表现为「产物在目录里、任务却报失败」。

判据（三条）
-----------
A. **新实现不崩**：全同矩阵（全 100）能正常出图（png + pdf 均生成）。
B. **灵敏度（对照臂）**：从 `archive/_bak_20261007/` **原样抽取**的旧实现
   在同一输入上**必须**抛异常 —— 否则本测试恒过。
C. **旧修复未回退**：常规矩阵上 Proportion 合计 == 1.0（钉住 2026-10-04 的 D3 修复）。

用法: python tests/_check_sdt_distribution.py
退出码: 0 = 通过；1 = 失败。
"""
import ast
import os
import sys
import tempfile

for _s in (sys.stdout, sys.stderr):
    try:
        _s.reconfigure(errors='replace')
    except (AttributeError, OSError):
        pass

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import matplotlib  # noqa: E402
matplotlib.use('Agg')

import numpy as np  # noqa: E402

from Virus_Platform_Core import sdt_exact as se  # noqa: E402

BACKUP = os.path.join(ROOT, 'archive', '_bak_20261007',
                      'Virus_Platform_Core', 'sdt_exact.py')
FAILS = []


def _check(cond, msg):
    print(('  ok   ' if cond else '  FAIL ') + msg)
    if not cond:
        FAILS.append(msg)


def extract_fn(path, name):
    src = open(path, encoding='utf-8').read()
    for n in ast.parse(src).body:
        if isinstance(n, ast.FunctionDef) and n.name == name:
            return ast.get_source_segment(src, n)
    raise RuntimeError(f'{path} 中找不到 {name}')


def identical_matrix(n=4):
    """全同输入：任意两两 identity 都是 100（含 NaN 对角线语义由掩码处理）。"""
    m = np.full((n, n), 100.0)
    np.fill_diagonal(m, np.nan)
    return m


def normal_matrix():
    """常规矩阵：15 对组合，覆盖 71–99 的整段值域（含低值，钉住 D3）。"""
    vals = [71, 74, 78, 80, 83, 85, 88, 89, 91, 93, 95, 96, 97, 98, 99]
    n = 6
    m = np.full((n, n), np.nan)
    k = 0
    for i in range(n):
        for j in range(i + 1, n):
            m[i][j] = m[j][i] = float(vals[k % len(vals)])
            k += 1
    return m, vals


def main():
    if not os.path.isfile(BACKUP):
        print(f'✘ 找不到修复前备份：{BACKUP}')
        return 1
    ns = dict(vars(se))
    exec(compile(extract_fn(BACKUP, 'plot_distribution'),
                 '<old plot_distribution>', 'exec'), ns)
    old_plot = ns['plot_distribution']
    print(f'对照臂：旧实现抽自 {os.path.relpath(BACKUP, ROOT)}')

    tmpdir = tempfile.mkdtemp(prefix='vp_sdtdist_')
    try:
        # ---------- B. 对照臂：旧实现必须崩 ==========
        print('\n[B] 灵敏度（对照臂）：旧实现在全同矩阵上必须抛异常')
        base_old = os.path.join(tmpdir, 'old_all100')
        raised = None
        try:
            old_plot(identical_matrix(), 'all-identical', base_old)
        except Exception as e:                        # noqa: BLE001
            raised = e
        _check(raised is not None,
               f'旧实现抛异常 = {type(raised).__name__ if raised else None}'
               f'（应非 None）')
        if raised is not None:
            print(f'         {type(raised).__name__}: {raised}')

        # ---------- A. 新实现不崩 ==========
        print('\n[A] 新实现：全同矩阵能正常出图')
        base_new = os.path.join(tmpdir, 'new_all100')
        try:
            se.plot_distribution(identical_matrix(), 'all-identical', base_new)
            ok = True
            err = ''
        except Exception as e:                        # noqa: BLE001
            ok, err = False, f'{type(e).__name__}: {e}'
        _check(ok, f'不抛异常{"" if ok else "（" + err + "）"}')
        _check(os.path.isfile(base_new + '.png'), '生成 .png')
        _check(os.path.isfile(base_new + '.pdf'), '生成 .pdf')

        # ---------- C. 常规矩阵：D3 修复未回退 ==========
        print('\n[C] 常规矩阵：Proportion 合计 == 1.0（钉住 2026-10-04 D3 修复）')
        mat, vals = normal_matrix()
        min_v = max(0.0, float(np.floor(np.nanmin(mat[~np.isnan(mat)]))))
        edges = np.arange(min_v, 101, 1)
        counts, _ = np.histogram(np.asarray(vals, dtype=float), bins=edges)
        total = counts.sum() / len(vals)
        _check(abs(total - 1.0) < 1e-12,
               f'bins 覆盖全域：合计 = {total:.6f}（min_v={min_v:g}）')
        base_n = os.path.join(tmpdir, 'new_normal')
        se.plot_distribution(mat, 'normal', base_n)
        _check(os.path.isfile(base_n + '.png'), '常规矩阵出图正常')
    finally:
        for f in os.listdir(tmpdir):
            os.remove(os.path.join(tmpdir, f))
        os.rmdir(tmpdir)

    print()
    if FAILS:
        print(f'✘ SDT 分布图守卫失败（{len(FAILS)} 条）')
        return 1
    print('✓ SDT 分布图守卫通过（全同输入不再崩，常规输入行为不变）')
    return 0


if __name__ == '__main__':
    sys.exit(main())
