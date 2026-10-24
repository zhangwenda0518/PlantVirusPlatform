# -*- coding: utf-8 -*-
"""GCVA-32 (Cytorhabdovirus sp. lycii OR489165.1) RDP5 断点分布图
数据 = 平台 rdp 工具真实产物（rdp5 breakpoint distribution.csv + ORFCoords）
图式 = RDP5 RGrapher 风格（红色柱 + CI 带 + ORF 盒）"""
import csv, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

RD = r'run/tool_runs/rdp_20261012_051924_db318d/rdp'
ALN_LEN = 16078

pos, num, hi99 = [], [], []
with open(f'{RD}/rdp5 breakpoint distribution.csv', encoding='utf-8') as f:
    for row in csv.reader(f):
        if not row or row[0].startswith('Position'):
            continue
        try:
            pos.append(float(row[0])); num.append(float(row[1])); hi99.append(float(row[2]))
        except (ValueError, IndexError):
            continue
pos, num, hi99 = np.array(pos), np.array(num), np.array(hi99)

orf = []
with open(r'D:/桌面/重组分析/bdp_debug/breakpoint_data.csvORFCoords.csv', encoding='utf-8') as f:
    for row in csv.reader(f):
        if not row or row[0] == 'Start':
            continue
        orf.append((int(row[0]), int(row[1]), row[2]))

# 基因组坐标 → 比对坐标（参考行去 gap 映射）
g2a = [0]
name = None
ref_seq = []
for line in open(r'D:/桌面/重组分析/gcva_32.fasta', encoding='utf-8'):
    line = line.strip()
    if line.startswith('>'):
        name = line[1:].split()[0]
        if name != 'OR489165.1' and ref_seq:
            break
    elif name == 'OR489165.1':
        ref_seq.append(line)
g = 0
for base in ''.join(ref_seq):
    if base != '-':
        g += 1
    g2a.append(g)  # g2a[g] = 该基因组位置后的比对列

fig, ax = plt.subplots(figsize=(15, 3.2), dpi=200)
# 99% CI 带
ax.fill_between(pos, 0, hi99, color='#f5c6c6', alpha=0.55, linewidth=0,
                label='Upper 99% CI')
# 断点柱（200nt 窗）
ax.bar(pos, num, width=200, color='#F08080', edgecolor='#a03c3c',
       linewidth=0.6, align='center', label='Recombination breakpoints')
# ORF 盒
ybox = max(hi99.max(), num.max()) * 1.12
box_colors = ['#5f9ea0', '#48d1cc', '#7fb2c9', '#5f9ea0', '#48d1cc', '#7fb2c9']
for i, (s, e, gsym) in enumerate(orf):
    a0, a1 = g2a[min(s, len(g2a) - 1)], g2a[min(e, len(g2a) - 1)]
    c = box_colors[i % len(box_colors)]
    ax.fill([a0, a1, a1 - 120, a0], [ybox, ybox, ybox - ybox*0.32, ybox - ybox*0.32],
            color=c, alpha=0.9, linewidth=0)
    ax.text((a0 + a1) / 2, ybox * 0.875, gsym, ha='center', va='center',
            fontsize=9, fontweight='bold', color='white' if c != '#48d1cc' else '#0b3c3c')

ax.set_xlim(-150, ALN_LEN + 150)
ax.set_ylim(0, ybox * 1.05)
ax.set_title('Cytorhabdovirus sp. lycii OR489165.1 (GCVA-32) — 20 isolates, 6 recombination events',
             fontsize=13, fontweight='bold')
ax.set_xlabel('Nucleotide position in the Cytorhabdovirus sp. lycii genome alignment', fontsize=11)
ax.set_ylabel('Recombination breakpoint number\n(200nt window)', fontsize=10)
ax.legend(loc='upper left', fontsize=9, framealpha=0.9)
ax.spines[['top', 'right']].set_visible(False)
fig.tight_layout()
fig.savefig(r'figures/assets/gcva_rdp_breakpoints.png', dpi=200)
print('saved gcva_rdp_breakpoints.png', 'events=6 aln=%d' % ALN_LEN)
