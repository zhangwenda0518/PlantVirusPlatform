# -*- coding: utf-8 -*-
"""CRISPR 引物交付序列纯度守卫（2026-10-07 审查 P0-2）。

背景
----
`crispr_design._enumerate_primers`（路线 B，primer3 兜底）内的 `_ok_primer`
原来只挡字面 `'N'`，于是 W/S/M/K/R/Y 等 IUPAC 兼并码全部放行 —— 而路线 A
（primer3 分支）与探针设计都用 `set(seq) - _PURE_ACGT` 守卫，模块注释也明确
声明"引物序列的纯 ACGT 约束由…引物后置守卫另行把关"。模板（consensus）在
2 碱基平票处会记兼并码（`build_consensus`），因此含兼并码的引物会一路流到
`ranked_sets.tsv` / `order_sheet.xlsx`（用户无法合成）。

判据（三条，缺一不可）
--------------------
A. **新实现**：含兼并码的窗口上，路线 B 产出的引物**全部为纯 ACGT**。
B. **灵敏度（对照臂）**：从 `archive/_bak_20261007/` **原样抽取**的旧实现
   在同一输入上**必须**产出含兼并码的引物 —— 否则本测试恒过。
C. **不过度过滤**：新实现在**纯 ACGT 窗口**上仍能产出引物（证明守卫没有把
   正常候选一并毙掉）。

用法: python tests/_check_crispr_primer_purity.py
退出码: 0 = 通过；1 = 失败。
"""
import ast
import os
import random
import sys

for _s in (sys.stdout, sys.stderr):
    try:
        _s.reconfigure(errors='replace')
    except (AttributeError, OSError):
        pass

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from Virus_Platform_Core import crispr_design as cd  # noqa: E402

BACKUP = os.path.join(ROOT, 'archive', '_bak_20261007',
                      'Virus_Platform_Core', 'crispr_design.py')
ACGT = set('ACGT')
FAILS = []


def _check(cond, msg):
    print(('  ok   ' if cond else '  FAIL ') + msg)
    if not cond:
        FAILS.append(msg)


def extract_fn(path, name):
    """从备份**原样抽取**函数源码（不重写、不复刻）。"""
    src = open(path, encoding='utf-8').read()
    for n in ast.parse(src).body:
        if isinstance(n, ast.FunctionDef) and n.name == name:
            return ast.get_source_segment(src, n)
    raise RuntimeError(f'{path} 中找不到 {name}')


PROFILE = {'product_min': 100, 'product_max': 200, 'min_cons_primer': 0.0}


def make_window(with_ambiguity):
    """600bp 窗口；with_ambiguity=True 时在侧翼埋 IUPAC 兼并码 'W'。"""
    random.seed(0)
    w = list(''.join(random.choice('ACGT') for _ in range(600)))
    if with_ambiguity:
        # 侧翼区（正/反向引物都可能落在这里）
        for pos in (10, 11, 12, 20, 21, 30, 31, 32,
                    300, 301, 302, 310, 311, 320, 321):
            w[pos] = 'W'
    return ''.join(w)


def impure(pairs):
    return [p for p in pairs
            if (set(p['f_seq']) | set(p['r_seq'])) - ACGT]


def main():
    if not os.path.isfile(BACKUP):
        print(f'✘ 找不到修复前备份：{BACKUP}')
        return 1

    # ---- 对照臂：旧实现（整函数从备份原样抽出执行，命名空间取模块全局）----
    ns = dict(vars(cd))
    exec(compile(extract_fn(BACKUP, '_enumerate_primers'),
                 '<old _enumerate_primers>', 'exec'), ns)
    old_enum = ns['_enumerate_primers']
    print(f'对照臂：旧实现抽自 {os.path.relpath(BACKUP, ROOT)}')

    args = lambda win: (win, [1.0] * len(win), 250, 270, 'RPA', PROFILE, 25, [])  # noqa: E731

    # ================= A. 新实现：含兼并码窗口 → 产出必须纯 ACGT =============
    print('\n[A] 新实现：含 "W" 的窗口上，路线 B 产出的引物必须全部纯 ACGT')
    win = make_window(True)
    new_pairs = cd._enumerate_primers(*args(win))
    bad = impure(new_pairs)
    _check(len(new_pairs) > 0, f'仍有引物产出（{len(new_pairs)} 对，证明未把候选全毙）')
    _check(len(bad) == 0, f'含非 ACGT 碱基的引物对 = {len(bad)}（应 0）')

    # ================= B. 对照臂：旧实现必须漏出兼并码 ======================
    print('\n[B] 灵敏度（对照臂）：旧实现必须漏出含兼并码的引物，否则本测试恒过')
    old_pairs = old_enum(*args(win))
    old_bad = impure(old_pairs)
    _check(len(old_bad) > 0,
           f'旧实现漏出的含兼并码引物对 = {len(old_bad)}（应 > 0）')
    for p in old_bad[:3]:
        print(f'         旧 f={p["f_seq"]}  r={p["r_seq"]}')

    # ================= C. 不过度过滤：纯 ACGT 窗口仍能产出 ==================
    print('\n[C] 不过度过滤：纯 ACGT 窗口上新实现仍产出引物')
    win2 = make_window(False)
    _check(all(c in 'ACGT' for c in win2), '构造窗口本身是纯 ACGT')
    new2 = cd._enumerate_primers(*args(win2))
    _check(len(new2) > 0, f'纯 ACGT 窗口上产出 {len(new2)} 对（应 > 0）')
    _check(len(impure(new2)) == 0, '纯 ACGT 窗口上产出亦全为纯 ACGT')

    print()
    if FAILS:
        print(f'✘ CRISPR 引物纯度守卫失败（{len(FAILS)} 条）')
        return 1
    print('✓ CRISPR 引物纯度守卫通过（路线 B 交付序列已限纯 ACGT）')
    return 0


if __name__ == '__main__':
    sys.exit(main())
