# -*- coding: utf-8 -*-
"""AA 同一性回归守卫：锚定框翻译不得退回「独立最长 ORF」伪值口径。

背景（2026-10-04 审查 M1）：旧实现逐条独立取最长 ORF 翻译，分离物/变体
序列里的一个无义突变即可打断共享长 ORF，使该条的最长 ORF 落到不同源
短框，identity_table 出现 NT 88.5% / AA 33.5% 这类生物学不可能值。
修复 = 分簇锚定翻译（_anchored_aa_seqs）。本守卫用内置示例集复跑一次，
断言高 NT 对不再出现 AA 塌陷。

用法: python tests/_check_identity_aa.py
"""
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

FAILS = []


def check(cond, msg):
    tag = 'ok  ' if cond else 'FAIL'
    print(f'  {tag} {msg}')
    if not cond:
        FAILS.append(msg)


def main():
    from Virus_Platform_Core.sdt_exact import run_identity_table
    from Virus_Platform_Core.config import get_config

    mafft = get_config().tool('mafft')
    if not mafft or not os.path.exists(mafft):
        print('SKIP: mafft 不可用')
        return 0

    import tempfile
    td = os.path.join(ROOT, 'run', '_tmp', 'id_aa_check')
    os.makedirs(td, exist_ok=True)
    try:
        class _L:
            def log(self, m):
                pass
        run_identity_table(
            os.path.join(ROOT, 'examples', 'example_virus_set.fasta'),
            td, mafft, logger=_L())
        table = os.path.join(td, 'identity_table.csv')
        check(os.path.exists(table), f'identity_table.csv 生成: {table}')

        bad, n_high = [], 0
        import csv as _csv
        with open(table, encoding='utf-8') as f:
            for row in _csv.DictReader(f):
                if not row.get('NT_identity(%)') or not row.get('AA_identity(%)'):
                    continue
                nt = float(row['NT_identity(%)'])
                aa = float(row['AA_identity(%)'])
                if nt >= 75.0:
                    n_high += 1
                    if aa < 55.0:
                        bad.append((row['Comparison'][:40], nt, aa))
        check(n_high >= 3, f'高 NT 对样本量足够（NT≥75% 共 {n_high} 对）')
        check(not bad,
              '高 NT 对无 AA 塌陷（NT≥75% 且 AA<55% 的对 = 0）'
              + (f'，违例: {bad[:3]}' if bad else ''))
    finally:
        import shutil
        shutil.rmtree(td, ignore_errors=True)

    if FAILS:
        print(f'✘ {len(FAILS)} 项未通过')
        return 1
    print('✓ AA 同一性回归守卫全部通过')
    return 0


if __name__ == '__main__':
    sys.exit(main())
