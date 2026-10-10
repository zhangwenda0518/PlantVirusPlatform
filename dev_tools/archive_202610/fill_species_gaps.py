# -*- coding: utf-8 -*-
"""从 Plant_Virus_Full 库中找出平台库缺失的物种，每物种取最优代表补入。

⚠️ 历史一次性脚本（2026-10-03 归档至 dev_tools/archive_202610/）：直接写
平台参考库（databases/virusref_db/kv_index，该目录已被 plant_ref 取代而
不存在），并依赖外部全量库路径 D:\\桌面\\C-host_classify\\...。**默认 dry-run
（只打印清单，不写任何文件）**；确需写入时必须显式加 --write，且先自行核对
库路径与备份。

用法: python fill_species_gaps.py [--write]
"""
import csv
import os
import sys
import time

import polars as pl

PLATFORM = r'D:\桌面\植物病毒分析平台'
FULL_INFO = r'D:\桌面\C-host_classify\plant_virus_db_pipeline\docs\data\Plant_Virus_Full.Info.tsv'
FULL_FASTA = r'D:\桌面\C-host_classify\plant_virus_db_pipeline\docs\data\Plant_Virus_Full.fasta'
LIB_FASTA = os.path.join(PLATFORM, 'databases', 'virusref_db', 'kv_index', 'reference.fasta')
LIB_INFO = os.path.join(PLATFORM, 'databases', 'virusref_db', 'kv_index', 'reference.ref_info.tsv')

# 默认 dry-run；只有显式 --write 才会写库（2026-10-03 改：此前默认即写）
DRY_RUN = '--write' not in sys.argv
MIN_LEN = 200


def norm_sp(s):
    import re
    return re.sub(r'[^a-z0-9 ]', '', (s or '').lower()).strip()


def main():
    t0 = time.time()
    print('=' * 72)
    print('[归档脚本 · 危险操作提示] fill_species_gaps.py 会直接追加写入平台'
          '参考库，\n  且绑定 2026-09 的旧库路径与外部全量库 D:\\桌面\\'
          'C-host_classify\\...。')
    print(f'  当前模式：{"DRY-RUN（只列清单，不写文件）" if DRY_RUN else "★ 写入模式（--write）★"}')
    print('=' * 72)

    # 1) 库内已有物种集合
    lib_sp = set()
    lib_info_cols = []
    with open(LIB_INFO, encoding='utf-8') as f:
        reader = csv.reader(f, delimiter='\t')
        lib_info_cols = next(reader)
        for row in reader:
            if len(row) > 2:
                lib_sp.add(norm_sp(row[2]))
    print(f'库内已有物种: {len(lib_sp)}')

    # 2) 加载全量库信息，找出缺失物种并选最优代表
    print('加载全量库信息 (199k)...')
    full = pl.scan_csv(FULL_INFO, separator='\t', infer_schema_length=0)

    # 只保留 Complete 且长度达标的
    full = full.filter(
        (pl.col('Nuc_Completeness') == 'complete') &
        (pl.col('Length').cast(pl.Int64, strict=False) >= MIN_LEN)
    )

    # 排序优先级: RefSeq > GenBank, 长度降序
    full = full.sort('Sequence_Type', descending=True).sort('Length', descending=True)

    # 加载到内存（筛选后应该少很多）
    full_df = full.collect()
    print(f'全量库 Complete ≥{MIN_LEN}bp: {full_df.height} 条')

    # 找缺失物种
    missing_accs = {}   # {norm_species: best_row}
    for row in full_df.iter_rows(named=True):
        sp_key = norm_sp(row.get('Species_NCBI') or row.get('Species_ICTV'))
        if not sp_key or sp_key in lib_sp:
            continue
        if sp_key not in missing_accs:
            missing_accs[sp_key] = row

    n_missing = len(missing_accs)
    print(f'缺失物种: {n_missing} 种')

    if n_missing == 0:
        print('无缺口，退出')
        return

    # dry-run 只输出清单
    if DRY_RUN:
        for sp_key, row in sorted(missing_accs.items()):
            print(f"  {row['Accession']:14s} {row['Species_NCBI'][:50]:52s} "
                  f"{row['Length']:>7}bp {row['Sequence_Type']}")
        print(f'\n共 {n_missing} 种缺失。'
              f'{"dry-run 结束——确需写入请显式加 --write。" if DRY_RUN else "去掉 --write 即为只读。"}')
        return

    # 3) 需要提取的 accession 集合
    target_accs = set()
    for row in missing_accs.values():
        acc = row['Accession'].split('.')[0]
        target_accs.add(acc)
    print(f'待提取 accession: {len(target_accs)}')

    # 4) 流式提取 FASTA 序列
    extracted = {}  # {base_acc: (header, seq)}
    nm, hdr, cur = None, '', []
    with open(FULL_FASTA, encoding='utf-8', errors='replace') as f:
        for ln in f:
            if ln.startswith('>'):
                if nm is not None and nm in target_accs:
                    extracted[nm] = ('>' + hdr, ''.join(cur))
                base = ln[1:].split()[0].split('.')[0].upper()
                nm = base
                hdr = ln[1:].strip()
                cur = []
            elif nm is not None:
                cur.append(ln.strip())
    if nm is not None and nm in target_accs:
        extracted[nm] = ('>' + hdr, ''.join(cur))
    print(f'FASTA 提取: {len(extracted)} / {len(target_accs)}')

    # 5) 追加到平台库
    n_added_fa = n_added_info = 0
    with open(LIB_FASTA, 'a', encoding='utf-8') as fa:
        with open(LIB_INFO, 'a', encoding='utf-8', newline='') as fi:
            writer = csv.writer(fi, delimiter='\t')
            for sp_key, row in sorted(missing_accs.items()):
                acc_base = row['Accession'].split('.')[0]
                if acc_base not in extracted:
                    continue
                hdr, seq = extracted[acc_base]
                # 写 FASTA
                fa.write(f'>{row["Accession"]} {row.get("GenBank_Title", sp_key)[:120]}\n')
                for i in range(0, len(seq), 70):
                    fa.write(seq[i:i+70] + '\n')
                n_added_fa += 1
                # 写 ref_info（列对齐库内格式）
                writer.writerow([
                    row['Accession'], row.get('Taxid', ''),
                    row.get('Species_NCBI', ''), row.get('Species_ICTV', ''),
                    row.get('Segment', ''), row.get('Sequence_Type', ''),
                    row.get('Category', ''), row.get('Topology', ''),
                    row.get('Molecule_type', ''), row.get('Molecule_Type2', ''),
                    row.get('Length', ''), row.get('Nuc_Completeness', ''),
                    row.get('GenBank_Title', ''), row.get('Geo_Location', ''),
                    row.get('USA', ''), row.get('Host', ''),
                    row.get('Isolation_Source', ''), row.get('Collection_Date', ''),
                    row.get('Release_Date', ''), row.get('NCBI_Status', ''),
                    row.get('Virus name(s)', ''), row.get('VMR_Species', ''),
                    row.get('VMR_Family', ''), row.get('VMR_Genus', ''),
                    row.get('Virus GENBANK accession', ''),
                    row.get('Genome coverage', ''),
                ])
                n_added_info += 1

    # 统计 FASTA 总条数
    n_total = sum(1 for _ in open(LIB_FASTA, encoding='utf-8') if _.startswith('>'))
    elapsed = time.time() - t0
    print(f'\n补充完成: FASTA +{n_added_fa}, ref_info +{n_added_info}')
    print(f'库内总条数: {n_total}')
    print(f'耗时: {elapsed:.1f}s')


if __name__ == '__main__':
    main()
