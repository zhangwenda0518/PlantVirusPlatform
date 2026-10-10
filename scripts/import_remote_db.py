# -*- coding: utf-8 -*-
"""把远端 39.106.101.94「植物病毒数据库平台」的全库导出导入为本地鉴定库。

用法（先在服务器侧导出两个文件）：
  python scripts/import_remote_db.py \
      --csv plant_virus_db.csv --fasta reference.fasta --name full_db

CSV 要求：含 accession 列（或第一列），可含 taxid / species / segment 列
（列名大小写不敏感、自动识别）。FASTA 头第一个词 = accession。
产物：databases/virusref_db/<name>/  自包含鉴定库布局
  reference.fasta + reference.ref_info.tsv
（之后在「已知病毒识别与定量」卡的「鉴定库」下拉即可选用；salmon 索引
 由引擎按需自建，或到「数据库构建」页预建。）
"""

import argparse
import csv
import os
import re
import shutil
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from Virus_Platform_Core.config import DIRS  # noqa: E402
from Virus_Platform_Core.utils import check_path  # noqa: E402

COL_MAP = {
    'accession': 'accession', 'acc': 'accession', 'accn': 'accession',
    'taxid': 'taxid', 'tax_id': 'taxid', 'tax': 'taxid',
    'species': 'species', 'sp': 'species',
    'segment': 'segment', 'seg': 'segment',
}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument('--csv', required=True, help='远端导出的全库 CSV/TSV')
    ap.add_argument('--fasta', required=True, help='全库序列 FASTA（头首词=accession）')
    ap.add_argument('--name', required=True, help='本地库目录名（如 full_db）')
    args = ap.parse_args()

    # 先校验输入（must_exist），再建输出目录：输入写错时不应留下一个空库目录，
    # 否则下次运行会撞上"库已存在"而看不出真正原因（2026-10-03）。
    try:
        csv_path = check_path(args.csv, must_exist=True)
        fasta_path = check_path(args.fasta, must_exist=True)
    except (FileNotFoundError, ValueError) as e:
        print(f'输入校验失败: {e}')
        return 2

    name = re.sub(r'[^\w\-.]+', '_', args.name).strip('_-')
    if not name:
        print('无效库名'); return 2
    lib_dir = os.path.join(DIRS['virus_src'], name)
    if os.path.isdir(lib_dir):
        print(f'库已存在: {lib_dir}（请换 --name 或先删除）'); return 2
    os.makedirs(lib_dir, exist_ok=False)
    try:
        return _build_lib(csv_path, fasta_path, lib_dir, name)
    except Exception:
        # 建库中途失败：回滚本次新建的目录，避免留下半成品让下次误报"库已存在"
        shutil.rmtree(lib_dir, ignore_errors=True)
        raise


def _build_lib(csv_path, fasta_path, lib_dir, name) -> int:
    """读 CSV/FASTA → 落盘 reference.fasta / reference.ref_info.tsv / 原始元数据。"""
    # ---- 读 CSV：自动识别列 ----
    with open(csv_path, encoding='utf-8-sig', errors='replace', newline='') as f:
        sample = f.read(65536)
        f.seek(0)
        delim = '\t' if sample.count('\t') > sample.count(',') else ','
        rows = list(csv.reader(f, delimiter=delim))
    if not rows:
        print('CSV 为空'); return 2
    header = [h.strip().lstrip('#').lower() for h in rows[0]]
    col_of = {}
    species_cols = []
    for i, h in enumerate(header):
        key = COL_MAP.get(h)
        if not key:
            # 前缀匹配：Species_ICTV / Species_NCBI → species 等
            for cand in COL_MAP:
                if h.startswith(cand):
                    key = cand
                    break
        if key == 'species':
            species_cols.append(i)      # 两列都记，行取值时 ICTV 优先
            continue
        if key and key not in col_of:
            col_of[key] = i
    if 'accession' not in col_of:
        col_of['accession'] = 0        # 第一列兜底
    print('列映射:', {k: header[v] for k, v in col_of.items()},
          '| 其余列忽略')

    want = set()
    ref_info_rows = []
    for row in rows[1:]:
        if not row or not any(c.strip() for c in row):
            continue
        acc = row[col_of['accession']].strip()
        if not acc:
            continue
        want.add(acc.split('.')[0])
        sp = ''
        for i in species_cols:
            if i < len(row) and row[i].strip():
                sp = row[i].strip()       # 文件列序在后者优先（本管线 ICTV 列在后）
        ref_info_rows.append((acc,
                              row[col_of['taxid']].strip() if 'taxid' in col_of else '',
                              sp,
                              row[col_of['segment']].strip() if 'segment' in col_of else ''))

    # ---- 过滤 FASTA ----
    kept, total = 0, 0
    with open(fasta_path, encoding='utf-8', errors='replace') as fin, \
         open(os.path.join(lib_dir, 'reference.fasta'), 'w',
              encoding='utf-8', newline='\n') as fout:
        keep = False
        for line in fin:
            if line.startswith('>'):
                total += 1
                acc = line[1:].split()[0] if line[1:].split() else ''
                keep = acc.split('.')[0] in want or acc in want
                if keep:
                    kept += 1
            if keep:
                fout.write(line)
    print(f'FASTA 共 {total} 条，命中库记录 {kept} 条')

    # 原始元数据整份随库保存（宿主/地理/年份等列是系统地理、RDP 选株的原料）
    shutil.copyfile(csv_path, os.path.join(lib_dir, 'reference.meta.original'
                    + os.path.splitext(csv_path)[1]))
    # ---- ref_info.tsv（平台口径列）----
    with open(os.path.join(lib_dir, 'reference.ref_info.tsv'), 'w',
              encoding='utf-8', newline='') as f:
        f.write('Accession\tTaxid\tSpecies\tSegment\n')
        for acc, taxid, sp, seg in ref_info_rows:
            f.write(f'{acc}\t{taxid or "Unannotated"}\t{sp or acc}\t{seg}\n')

    print(f'完成：{lib_dir}')
    print('下一步：「已知病毒识别与定量」卡 → 鉴定库下拉选', name,
          '（或到「数据库构建」页预建 salmon 索引）')
    return 0


if __name__ == '__main__':
    sys.exit(main())
