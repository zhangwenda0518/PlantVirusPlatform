# -*- coding: utf-8 -*-
"""修复 Plant_Virus_Ref 数据集：表-序列对齐 + 补 TSWV 三段。

输入（databases/_dl_staging/）：
  Plant_Virus_Ref.fasta / Plant_Virus_Ref.Info.tsv   —— 服务器下载的破损版
  tswv_3seg.fasta / tswv_3seg.info.rows              —— 服务器 Full(09-18) 提取
输出：
  Plant_Virus_Ref.fixed.fasta / Plant_Virus_Ref.fixed.Info.tsv
修复规则：
  * accession 版本号归一（AJ238493.1 ↔ AJ238493）后表-序列配对；
  * 配对成功的记录：info 的 Accession 列改写为 fasta 头的精确形式；
  * 同 key 多版本只保留与 info 精确匹配的那条（代表集不留重复版本）；
  * 仍无配对（只有表/只有序列）的记录剔除并计数；
  * TSWV（Orthotospovirus tomatomaculae）M/S/L 三段注入，列格式对齐
    近缘种行；Taxid 用该表既有 tomatomaculae 行的 3052585。
"""
import io
import os
import sys

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
STG = os.path.join(BASE, 'databases', '_dl_staging')


def read_fasta(path):
    """→ [(header_accession, [lines])]（保持顺序）。"""
    out = []
    cur = None
    with io.open(path, encoding='utf-8') as f:
        for ln in f:
            ln = ln.rstrip('\n\r')
            if ln.startswith('>'):
                cur = (ln[1:].strip().split()[0], [])
                out.append(cur)
            elif cur is not None:
                cur[1].append(ln)
    return out


def base_key(acc):
    return acc.strip().split('.')[0].strip()


def main():
    fa = read_fasta(os.path.join(STG, 'Plant_Virus_Ref.fasta'))
    with io.open(os.path.join(STG, 'Plant_Virus_Ref.Info.tsv'),
                 encoding='utf-8') as f:
        head = f.readline().rstrip('\n').split('\t')
        info_rows = [ln.rstrip('\n').split('\t') for ln in f]

    ai = head.index('Accession')
    info_by_key = {}
    dup_info = 0
    for r in info_rows:
        k = base_key(r[ai])
        if k in info_by_key:
            dup_info += 1
            continue
        info_by_key[k] = r

    # fasta: base_key → [(acc, lines)]（保序）
    fa_by_key = {}
    for acc, lines in fa:
        fa_by_key.setdefault(base_key(acc), []).append((acc, lines))

    paired, dropped_fa_only, dropped_info_only, multi_ver = [], [], [], 0
    seen_keys = set()
    for acc, lines in fa:
        k = base_key(acc)
        if k in seen_keys:
            continue                     # 同 key 的第二个版本：跳过
        row = info_by_key.get(k)
        if row is None:
            dropped_fa_only.append(acc)
            continue
        seen_keys.add(k)
        # 有精确版本匹配就用它，否则用第一条同 key 序列
        cands = fa_by_key[k]
        pick = next((c for c in cands if c[0] == row[ai].strip()), cands[0])
        if len(cands) > 1:
            multi_ver += 1
        row = row[:]

        row[ai] = pick[0]                # Accession 改写为 fasta 精确头
        paired.append((pick, row))
    info_keys = {base_key(r[ai]) for r in info_rows}
    dropped_info_only = [k for k in info_keys
                         if k not in {base_key(a) for a, _ in fa}]
    dropped_info_only = [r[ai] for r in info_rows
                         if base_key(r[ai]) in set(dropped_info_only)]

    # ---- TSWV 三段 ----
    tswv_fa = read_fasta(os.path.join(STG, 'tswv_3seg.fasta'))
    with io.open(os.path.join(STG, 'tswv_3seg.info.rows'),
                 encoding='utf-8') as f:
        tswv_full = [ln.rstrip('\n').split('\t') for ln in f if ln.strip()]
    # Full(21列): 0 Acc,1 Taxid(空),2 Sp_NCBI,3 Sp_ICTV,4 Family,5 Segment,
    #             6 Topology,7 Molecule_type,8 Mol_Type2,9 Seq_Type,10 Length,
    #             11 Nuc_Comp,12 GenBank_Title,13..20 其余（Geo/USA/Host/
    #             Isol/CollDate/RelDate/NCBI_Status/Segment_std）
    TSWV_TAXID = '3052585'               # 表内既有 tomatomaculae 行同款
    tswv_rows = []
    for r in tswv_full:
        out = [''] * len(head)
        out[head.index('Accession')] = r[0]
        out[head.index('Taxid')] = TSWV_TAXID
        out[head.index('Species_NCBI')] = r[2]
        out[head.index('Species_ICTV')] = r[3]
        out[head.index('Segment')] = r[5]
        out[head.index('Sequence_Type')] = r[9]
        out[head.index('Category')] = 'Segmented_Complete'
        out[head.index('Topology')] = r[6]
        out[head.index('Molecule_type')] = r[7]
        out[head.index('Molecule_Type2')] = r[8]
        out[head.index('Length')] = r[10]
        out[head.index('Nuc_Completeness')] = r[11]
        out[head.index('GenBank_Title')] = r[12]
        out[head.index('Collection_Date')] = r[18] if len(r) > 18 else ''
        out[head.index('Virus name(s)')] = 'Tomato spotted wilt virus'
        out[head.index('VMR_Species')] = r[3]
        out[head.index('VMR_Family')] = r[4]
        out[head.index('VMR_Genus')] = 'Orthotospovirus'
        tswv_rows.append((r[0], out))

    # ---- 写出 ----
    fa_out = os.path.join(STG, 'Plant_Virus_Ref.fixed.fasta')
    info_out = os.path.join(STG, 'Plant_Virus_Ref.fixed.Info.tsv')
    with io.open(fa_out, 'w', encoding='utf-8', newline='\n') as f:
        for (acc, lines), _row in paired:
            f.write('>' + acc + '\n')
            for ln in lines:
                f.write(ln + '\n')
        for acc, lines in tswv_fa:
            f.write('>' + acc + '\n')
            for ln in lines:
                f.write(ln + '\n')
    with io.open(info_out, 'w', encoding='utf-8', newline='\n') as f:
        f.write('\t'.join(head) + '\n')
        for _, row in paired:
            f.write('\t'.join(row) + '\n')
        for _, row in tswv_rows:
            f.write('\t'.join(row) + '\n')

    n_fa = len(paired) + len(tswv_fa)
    n_info = len(paired) + len(tswv_rows)
    print(f'配对保留: {len(paired)} | 注入 TSWV: {len(tswv_fa)}')
    print(f'修复后 fasta={n_fa} info={n_info}（相等={n_fa == n_info}）')
    print(f'剔除-序列有表无: {len(dropped_fa_only)} 示例{dropped_fa_only[:3]}')
    print(f'剔除-表有序列无: {len(dropped_info_only)} '
          f'示例{dropped_info_only[:3]}')
    print(f'同 key 多版本: {multi_ver} | info 重复行: {dup_info}')

    # ---- 自检 ----
    chk = read_fasta(fa_out)
    keys = [base_key(a) for a, _ in chk]
    assert len(keys) == len(set(keys)), 'base key 仍有重复！'
    with io.open(info_out, encoding='utf-8') as f:
        h2 = f.readline().rstrip('\n').split('\t')
        rows2 = [ln.rstrip('\n').split('\t') for ln in f]
    a2 = h2.index('Accession')
    keys2 = [base_key(r[a2]) for r in rows2]
    assert set(keys) == set(keys2), '配对后仍有错位！'
    hit = [r for r in rows2 if 'tomatomaculae' in '\t'.join(r).lower()]
    print(f'自检通过: {len(chk)} 条序列，表序列一致，TSWV 行 {len(hit)} 条')


if __name__ == '__main__':
    sys.exit(main())
