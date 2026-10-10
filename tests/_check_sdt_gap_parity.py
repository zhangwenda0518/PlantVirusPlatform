# -*- coding: utf-8 -*-
"""SDT gap 口径一致性守卫（2026-10-07 审查 P0-1）。

背景
----
`sdt_exact._sdt_formula` 原来只把 `'-'` 当 gap，而 `phylo.pairwise_identity_matrix`
早已统一为 `'-.'`（其注释："gap 判定统一 `in '-.'` …旧实现只认 '-'，把 '.'
当成碱基位点参与比较"）。后果：同一份含 `.` 的已比对 MSA，在 /sdt 页与
/compare 页会给出两个答案（实测一侧为 `.`、另一侧为碱基时差 9.72 个百分点），
足以跨过 80% 的物种划界阈值。

判据（三条，缺一不可）
--------------------
A. **一致性**：对含 `.` / `-` 的 MSA，`_sdt_formula` 逐对 == `pairwise_identity_matrix`。
B. **灵敏度（对照臂）**：从 `archive/_bak_20261007/` **原样抽取**的旧实现在同一输入上
   必须给出**不同**的值。若旧实现也"通过"，说明本测试抓不到该缺陷（恒过）。
C. **改动面受限**：不含 `.` 的输入上，新旧实现必须**零差异** —— 证明本次只动了
   含 `.` 的那一档，没有顺手改到别的语义。

用法: python tests/_check_sdt_gap_parity.py
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

from Virus_Platform_Core.sdt_exact import _sdt_formula as new_fn  # noqa: E402
from Virus_Platform_Core.phylo import pairwise_identity_matrix  # noqa: E402

BACKUP = os.path.join(ROOT, 'archive', '_bak_20261007',
                      'Virus_Platform_Core', 'sdt_exact.py')

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


def write_msa(path, rows):
    with open(path, 'w', encoding='utf-8', newline='\n') as f:
        for name, seq in rows:
            f.write(f'>{name}\n{seq}\n')


def main():
    # ---- 对照臂：旧实现（从备份原样抽出执行）----
    if not os.path.isfile(BACKUP):
        print(f'✘ 找不到修复前备份：{BACKUP}')
        return 1
    ns = {}
    exec(compile(extract_fn(BACKUP, '_sdt_formula'),
                 '<old sdt_exact._sdt_formula>', 'exec'), ns)
    old_fn = ns['_sdt_formula']
    print(f'对照臂：旧实现抽自 {os.path.relpath(BACKUP, ROOT)}')

    tmpdir = tempfile.mkdtemp(prefix='vp_sdtgap_')
    try:
        # ================= 用例 1：含 '.' 的 MSA =================
        # s2 在 idx4 是 '.'（应视为 gap）、idx11 是真错配
        s1 = 'ACGTACGTACGT'
        s2 = 'ACGT.CGTACGA'
        s3 = 'ACGT-CGTACGT'
        rows = [('s1', s1), ('s2', s2), ('s3', s3)]
        p1 = os.path.join(tmpdir, 'with_dot.fasta')
        write_msa(p1, rows)

        names, mat = pairwise_identity_matrix(p1)
        idx = {n: i for i, n in enumerate(names)}
        seqs = dict(rows)

        print('\n[A] 一致性：_sdt_formula == pairwise_identity_matrix（含 "." 输入）')
        worst = 0.0
        for i in range(len(names)):
            for j in range(i + 1, len(names)):
                a, b = seqs[names[i]], seqs[names[j]]
                got = new_fn(a, b)
                exp = mat[i][j]
                worst = max(worst, abs(got - exp))
                _check(abs(got - exp) < 1e-9,
                       f'{names[i]} vs {names[j]}: sdt={got:.4f} phylo={exp:.4f}')
        _check(worst < 1e-9, f'全部对的最大偏差 {worst:.2e} < 1e-9')

        print('\n[B] 灵敏度（对照臂）：旧实现必须给出不同值，否则本测试恒过')
        diffs = []
        for i in range(len(names)):
            for j in range(i + 1, len(names)):
                a, b = seqs[names[i]], seqs[names[j]]
                o, n = old_fn(a, b), new_fn(a, b)
                if abs(o - n) > 1e-9:
                    diffs.append((names[i], names[j], o, n))
        _check(len(diffs) > 0, f'旧实现与新实现存在差异的对数 = {len(diffs)}（应 > 0）')
        for nm1, nm2, o, n in diffs:
            print(f'         {nm1} vs {nm2}: 旧 {o:.4f} → 新 {n:.4f}'
                  f'（差 {n - o:+.4f} 个百分点）')

        # ================= 用例 2：不含 '.' 的 MSA（改动面受限）=================
        print('\n[C] 改动面受限：不含 "." 的输入上新旧必须零差异')
        r2 = [('t1', 'ACGT-ACGTACGT'), ('t2', 'ACGT-ACGTACGA'),
              ('t3', 'ACGTACGTACGT')]
        p2 = os.path.join(tmpdir, 'no_dot.fasta')
        write_msa(p2, r2)
        mx = 0.0
        for i in range(len(r2)):
            for j in range(i + 1, len(r2)):
                a, b = r2[i][1], r2[j][1]
                mx = max(mx, abs(old_fn(a, b) - new_fn(a, b)))
        _check(mx < 1e-12, f'无 "." 输入上的最大新旧差 {mx:.2e} == 0')
    finally:
        for f in os.listdir(tmpdir):
            os.remove(os.path.join(tmpdir, f))
        os.rmdir(tmpdir)

    print()
    if FAILS:
        print(f'✘ SDT gap 口径守卫失败（{len(FAILS)} 条）')
        return 1
    print('✓ SDT gap 口径守卫通过（_sdt_formula 与 pairwise_identity_matrix 同口径）')
    return 0


if __name__ == '__main__':
    sys.exit(main())
