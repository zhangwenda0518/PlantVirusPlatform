# -*- coding: utf-8 -*-
"""smoke_e2e 最小样品夹具：AI 助手相关测试（_it_ai_* / _check_ai_*）共用。

背景：这批测试历史上依赖 run/results/smoke_e2e 这个**遗留样品**——
它一旦被清理（tool-runs 归档清理 / 手工删 run/ / 换机器首跑），5 个
AI 测试一起红，且报错点与 AI 代码无关、极难排查（2026-10-13 实测：
digest 返回 None → 广播层静默跳过 → '关键阶段触发 ai_read' 失败）。

改为按需自建、进程退出时只删**自己新建**的那份（幂等：样品已存在时
只补缺失的表文件、退出时不删，不碰别人的产物）。
"""
import atexit
import os
import shutil

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# all_viruses.summary.tsv 的表头列名与 kvsuite 真实产物一致
# （read_result_table 断言返回含 'header'，digest 的 kvsuite 段读它取物种数）
HEADER = ('Species\tReads\tDNA\tGenome\tRefSeq_acc\tCoverage\t'
          'Reads_pct\tDepth\tTaxID')
ROWS = [
    ('Cucumber mosaic virus', '1204', '0.982', 'NC_002035.1',
     'NC_002035.1', '0.91', '16.7', '48.2', '3052486'),
    ('Tomato mosaic virus', '412', '0.977', 'NC_011807.1',
     'NC_011807.1', '0.44', '5.7', '21.0', '12236'),
    ('Potato virus Y', '218', '0.991', 'NC_001616.1',
     'NC_001616.1', '0.12', '3.0', '9.6', '5730'),
]


def ensure():
    """确保 smoke_e2e 样品与 ②b 鉴定表存在（幂等），返回样品目录路径。"""
    sample_dir = os.path.join(ROOT, 'run', 'results', 'smoke_e2e')
    tsv = os.path.join(sample_dir, '02b_kvsuite', 'identify',
                       'all_viruses.summary.tsv')
    created = not os.path.isdir(sample_dir)
    os.makedirs(os.path.dirname(tsv), exist_ok=True)
    if not os.path.isfile(tsv):
        with open(tsv, 'w', encoding='utf-8', newline='') as f:
            f.write(HEADER + '\n')
            for row in ROWS:
                f.write('\t'.join(row) + '\n')
    if created:
        atexit.register(_cleanup, sample_dir)
    return sample_dir


def _cleanup(sample_dir):
    shutil.rmtree(sample_dir, ignore_errors=True)
