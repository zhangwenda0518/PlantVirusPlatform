# -*- coding: utf-8 -*-
r"""fig_b9：min contig 门槛扫描——覆盖率/一致率 vs 最短 contig 阈值

输入：run_contig_sweep.py 的扫描日志（外部产物）。
路径解析：argv[1] > 环境变量 VP_BENCH_SWEEP_LOG > 仓库内
results/sweep.log（已归档的本次运行日志）。输入缺失时打印说明并跳过该图，
不报堆栈（2026-10-03 参数化，此前硬编码 C:\vp_srna_test\bench\sweep.log）。
"""
import json
import os
import re
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams['font.sans-serif'] = ['Microsoft YaHei']
plt.rcParams['axes.unicode_minus'] = False

HERE = os.path.dirname(os.path.abspath(__file__))
LOG = (sys.argv[1] if len(sys.argv) > 1
       else os.environ.get('VP_BENCH_SWEEP_LOG')
       or os.path.join(HERE, 'results', 'sweep.log'))
if not os.path.isfile(LOG):
    print(f'[跳过] 未找到门槛扫描日志：{LOG}\n'
          f'       该图为外部运行产物（run_contig_sweep.py 的输出日志）。\n'
          f'       请用 argv[1] 或环境变量 VP_BENCH_SWEEP_LOG 指定日志路径，'
          f'或先跑 run_contig_sweep.py 生成。')
    sys.exit(0)
txt = open(LOG, encoding='utf-8', errors='replace').read()
m = re.search(r'RESULT:(\{.*)', txt)
data = json.loads(m.group(1).split('\n')[0])
json.dump(data, open(os.path.join(HERE, 'results', 'contig_sweep.json'), 'w',
                     encoding='utf-8'), ensure_ascii=False, indent=1)
FIGS = os.path.join(HERE, 'figs')
sweep = data['sweep']
cuts = [r['cutoff'] for r in sweep]
MID, WARN, OK = '#4C72B0', '#D9A441', '#5B9A57'

fig, (axl, axr) = plt.subplots(1, 2, figsize=(11.6, 4.4))
r1 = [r['targets']['NC_043447.1']['coverage_pct'] for r in sweep]
r2 = [r['targets']['NC_043448.1']['coverage_pct'] for r in sweep]
ns = [r['n_contigs'] for r in sweep]

axl.plot(cuts, r1, 'o-', color=MID, lw=1.8, label='PVB RNA1 (7,148bp)')
axl.plot(cuts, r2, 's-', color=WARN, lw=1.8, label='PVB RNA2 (4,527bp)')
for c, v in zip(cuts, r2):
    axl.annotate(f'{v:.1f}%', (c, v), textcoords='offset points',
                 xytext=(0, -14), ha='center', fontsize=8, color=WARN)
for c, v in zip(cuts, r1):
    axl.annotate(f'{v:.1f}%', (c, v), textcoords='offset points',
                 xytext=(0, 7), ha='center', fontsize=8, color=MID)
axl.set_xlabel('min contig 门槛 (bp)')
axl.set_ylabel('PVB 基因组覆盖率 (%)')
axl.set_title('放宽 min contig 的覆盖收益（本数据集有限）', fontsize=11)
axl.set_xticks(cuts)
axl.set_ylim(0, 22)
axl.legend(fontsize=8.5)
axl.spines[['top', 'right']].set_visible(False)

axr.plot(cuts, ns, 'o-', color=OK, lw=1.8)
for c, n in zip(cuts, ns):
    axr.annotate(str(n), (c, n), textcoords='offset points', xytext=(0, 7),
                 ha='center', fontsize=8, color=OK)
axr.set_xlabel('min contig 门槛 (bp)')
axr.set_ylabel('保留 contig 数')
axr.set_title('保留 contig 数（一致率全程 100%）', fontsize=11)
axr.set_xticks(cuts)
axr.spines[['top', 'right']].set_visible(False)

fig.suptitle('min contig 门槛扫描：Bari-S3 srna 组装 → PVB RNA1/RNA2', fontsize=11.5)
fig.tight_layout(rect=(0, 0, 1, 0.93))
fig.savefig(os.path.join(FIGS, 'fig_b9_sweep.png'), dpi=200)
fig.savefig(os.path.join(FIGS, 'fig_b9_sweep.svg'))
print('fig_b9 ->', FIGS)
