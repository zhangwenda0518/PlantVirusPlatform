# -*- coding: utf-8 -*-
"""校验下载的 Plant_Virus_Ref 数据集：表-序列对齐、TSWV 线索。"""
import io
import os

BASE = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                    'databases', '_dl_staging')
info_p = os.path.join(BASE, 'Plant_Virus_Ref.Info.tsv')
fa_p = os.path.join(BASE, 'Plant_Virus_Ref.fasta')

info_acc = set()
with io.open(info_p, encoding='utf-8') as f:
    head = f.readline().rstrip('\n').split('\t')
    ai = head.index('Accession')
    for ln in f:
        c = ln.split('\t')
        if len(c) > ai and c[ai].strip():
            info_acc.add(c[ai].strip())

fa_acc = set()
with io.open(fa_p, encoding='utf-8') as f:
    for ln in f:
        if ln.startswith('>'):
            fa_acc.add(ln[1:].strip().split()[0])

print('info accessions:', len(info_acc), '| fasta accessions:', len(fa_acc))
only_info = sorted(info_acc - fa_acc)
only_fa = sorted(fa_acc - info_acc)
print('info 有 序列无:', len(only_info), only_info[:5])
print('序列有 info 无:', len(only_fa), only_fa[:5])

# TSWV 各类线索计数
with io.open(info_p, encoding='utf-8') as f:
    rows = [ln.rstrip('\n').split('\t') for ln in f]
def count(*terms):
    n = 0
    for r in rows:
        line = '\t'.join(r).lower()
        if all(t.lower() in line for t in terms):
            n += 1
    return n
print('含 tomato spotted wilt:', count('tomato spotted wilt'))
print('含 tomatomaculae:', count('tomatomaculae'))
print('含 tswv:', count('tswv'))
