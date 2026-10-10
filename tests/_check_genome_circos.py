# -*- coding: utf-8 -*-
"""圈图回归守卫：plot_genome_circos 必须产出真圈图（返回 True），
不得静默回退线性图。

背景（2026-10-04 审查 M2）：pycirclize>=1.9 移除 Sector.text(deg=) 形参，
旧调用 deg=90 经 **kwargs 传给 matplotlib Text 抛 AttributeError，
report 重生成时 3/3 contig 圈图静默回退线性图。修复后以本守卫钉住。

用法: python tests/_check_genome_circos.py
"""
import os
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
    from Bio import SeqIO
    from Virus_Platform_Core.viz import plot_genome_circos, _pretty_contig_name

    gb = os.path.join(ROOT, 'examples', 'example_genome.gb')
    if not os.path.isfile(gb):
        print('SKIP: 示例 GenBank 不存在')
        return 0
    rec = next(SeqIO.parse(gb, 'genbank'))
    seq = str(rec.seq)
    bed = [
        {'start': 100, 'end': 1500, 'strand': '+'},
        {'start': 2000, 'end': 4200, 'strand': '-'},
        {'start': 4400, 'end': 5900, 'strand': '+'},
    ]
    out = os.path.join(ROOT, 'run', '_tmp', 'circos_check.png')
    os.makedirs(os.path.dirname(out), exist_ok=True)
    ok = plot_genome_circos('NODE_1_length_6240_cov_29.664667', seq, bed, out)
    check(ok is True, '返回 True（真圈图，未回退线性图）')
    check(os.path.isfile(out) and os.path.getsize(out) > 20000,
          f'PNG 产出且非空（{os.path.getsize(out) if os.path.exists(out) else 0} B）')

    check(_pretty_contig_name('NODE_1_length_6240_cov_29.664667')
          == 'NODE_1（6,240 bp）', 'SPAdes 节点名 → 可读标题')
    check(_pretty_contig_name('ctg_abc') == 'ctg_abc', '非节点名原样返回')

    if os.path.exists(out):
        os.remove(out)
    if FAILS:
        print(f'✘ {len(FAILS)} 项未通过')
        return 1
    print('✓ 圈图回归守卫全部通过')
    return 0


if __name__ == '__main__':
    sys.exit(main())
