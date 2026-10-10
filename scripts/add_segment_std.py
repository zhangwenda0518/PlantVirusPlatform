# -*- coding: utf-8 -*-
"""给 docs/data 的 Info 表追加 Segment_std 列（用平台规范化规则）。

⚠️ 仅数据库维护者使用（维护者一次性工具，非平台运行链路）：
   默认 --data-dir 指向维护者机器上的 plant_virus_db_pipeline/docs/data，
   脚本会**原位改写**该目录下的 Info 表（先备份 .bak_segmentstd）。
   普通用户运行平台不需要本脚本。

- 原始文件先备份为 <名>.bak_segmentstd；
- 列追加在末尾，其余列与行序原样保留（逐行文本处理，不重排引号）；
- 规范化用 norm_segment_ctx（Segment + Family + Molecule_type 上下文）。

用法: python scripts/add_segment_std.py [--data-dir <docs/data 目录>]
"""

import argparse
import os
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
from Virus_Platform_Core.segment_norm import norm_segment_ctx  # noqa: E402

# 默认 = 维护者机器上的路径；换机/换目录请用 --data-dir 覆盖
DEFAULT_DATA_DIR = r'D:\桌面\C-host_classify\plant_virus_db_pipeline\docs\data'
TARGET_NAMES = [
    'Plant_Virus_Full.Info.tsv',
    'Plant_Virus.complete_ref_info.tsv',
    'Plant_Virus_Ref.Info.tsv',
]


def col_index(header):
    low = [h.strip().lstrip('#').lower() for h in header]
    def find(*names):
        for n in names:
            if n in low:
                return low.index(n)
        return None
    return (find('segment'), find('family'), find('molecule_type'),
            find('moleculetype'))


def main(data_dir: str) -> int:
    # 入口先校验输入目录：默认路径只存在于维护者机器上，普通用户误运行时
    # 应给出明确指引，而不是静默"全部跳过"后返回 0。
    if not os.path.isdir(data_dir):
        print(f'✘ 数据目录不存在: {data_dir}')
        print('  本脚本仅供数据库维护者使用；请用 --data-dir 指向'
              ' plant_virus_db_pipeline 的 docs/data 目录')
        return 1
    targets = [os.path.join(data_dir, n) for n in TARGET_NAMES]
    for path in targets:
        if not os.path.isfile(path):
            print('跳过（不存在）:', path)
            continue
        bak = path + '.bak_segmentstd'
        if not os.path.isfile(bak):
            shutil.copy2(path, bak)
        else:
            # 原位改写型脚本的二跑陷阱：备份存在≠可沿用——当前文件可能
            # 已是上一轮的改写产物，再拿它当原稿会把改写结果固化进备份
            import filecmp
            if filecmp.cmp(path, bak, shallow=False):
                print('  备份已是原稿，沿用:', bak)
            else:
                print('  备份保留首跑原稿（当前文件与备份不同，'
                      '视为已改写产物）:', bak)
        with open(path, encoding='utf-8', errors='replace') as f:
            lines = f.read().split('\n')
        # 去掉末尾空行后处理，再统一补回
        trail = 0
        while lines and lines[-1] == '':
            lines.pop()
            trail += 1
        header = lines[0].split('\t')
        if col_index(header)[0] is None:
            print('跳过（无 Segment 列）:', path)
            continue
        i_seg, i_fam, i_mol = col_index(header)[:3]
        if i_seg is None:
            print('跳过（无 Segment 列）:', path)
            continue
        out = [lines[0] + '\tSegment_std']
        n = 0
        for line in lines[1:]:
            if not line:
                continue
            cols = line.split('\t')

            def at(idx):
                if idx is not None and idx < len(cols):
                    return cols[idx].strip()
                return ''

            std = norm_segment_ctx(at(i_seg), at(i_fam), at(i_mol))
            out.append(line + '\t' + std)
            n += 1
        with open(path, 'w', encoding='utf-8', newline='') as f:
            f.write('\n'.join(out) + '\n')
        name = os.path.basename(path)
        print(f'✔ {name}: {n} 行已追加 Segment_std 列')
    print('DONE')
    return 0


if __name__ == '__main__':
    _p = argparse.ArgumentParser(
        description='给数据库 Info 表追加 Segment_std 列（仅数据库维护者使用）')
    _p.add_argument('--data-dir', default=DEFAULT_DATA_DIR,
                    help='维护者机器上的 docs/data 目录'
                         f'（默认: {DEFAULT_DATA_DIR}）')
    _args = _p.parse_args()
    sys.exit(main(_args.data_dir))
