# -*- coding: utf-8 -*-
"""检查参考库内各属/科的覆盖情况与宿主分布（基准评测配套脚本）。

用法:
    python check_coverage.py [--db-dir <库目录>] [--expect <期望科清单.json>]

库目录解析优先级：--db-dir > 环境变量 VP_BENCH_DB > 旧默认路径
`databases/virusref_db/kv_index`（2026-09-29 起被 plant_ref 取代，现已不存在；
缺失时报清晰中文错误——本脚本的快照库需先重建或显式指定）。

2026-10-03 修复：
  ① 库路径参数化（此前硬编码 D:\\桌面\\植物病毒分析平台 + 老库路径）；
  ② 原"落入 0 覆盖的科"硬编码名单（Reoviridae/... 7 个科）移除，改为从库
     数据动态生成两类结果：a) 库内无植物宿主物种的科（按 Host 列计算）；
     b) 关键科缺口对照——需要外部期望清单时用 --expect 传入 JSON
     （{"families": ["Reoviridae", ...]}），未提供时明确打印"跳过"，不猜。
"""
import argparse
import collections
import csv
import json
import os
import sys

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))
DB_DEFAULT = os.path.join(REPO, 'databases', 'virusref_db', 'kv_index')


def parse_args():
    ap = argparse.ArgumentParser(description='参考库科/属覆盖与宿主分布检查')
    ap.add_argument('--db-dir', default=None,
                    help='参考库目录（含 reference.ref_info.tsv）；'
                         '默认取环境变量 VP_BENCH_DB，否则为旧路径 kv_index')
    ap.add_argument('--expect', default=None,
                    help='可选：期望科清单 JSON（{"families": [...]}），'
                         '用于对照"关键科是否在库内"')
    return ap.parse_args()


def main():
    args = parse_args()
    db_dir = args.db_dir or os.environ.get('VP_BENCH_DB') or DB_DEFAULT
    info_path = os.path.join(db_dir, 'reference.ref_info.tsv')
    if not os.path.isfile(info_path):
        raise SystemExit(
            f'[基准库缺失] 找不到注释表：{info_path}\n'
            f'说明：基准快照绑定 2026-09 的 kv_index 库（manifest 2026-09-19），'
            f'2026-09-29 起平台默认库已改为 virusref_db/plant_ref/。\n'
            f'请先重建等价库，或用 --db-dir <库目录> / 环境变量 VP_BENCH_DB '
            f'显式指定（目录内需有 reference.ref_info.tsv）。')

    fam_sp = collections.defaultdict(set)
    fam_plant = collections.defaultdict(bool)
    host_cat = collections.Counter()
    all_species = set()

    for row in csv.reader(open(info_path, encoding='utf-8'), delimiter='\t'):
        if row and not row[0].startswith('#') and len(row) > 22:
            fam = row[22].strip() or row[3].strip()
            sp = row[2].strip()
            host = row[14].strip().lower()
            if fam:
                fam_sp[fam].add(sp)
            if sp:
                all_species.add(sp)
            if sp:
                is_plant = 'viridiplantae' in host or 'plant' in host
                if is_plant:
                    host_cat['plant'] += 1
                elif any(x in host for x in ('arthropod', 'insect', 'animal',
                                             'vertebrate')):
                    host_cat['animal'] += 1
                elif host:
                    host_cat['other'] += 1
                else:
                    host_cat['empty'] += 1
                if fam and is_plant:
                    fam_plant[fam] = True

    print("=== Host 列分布 ===")
    for k, v in host_cat.most_common():
        print(f"  {k:10s}: {v}")

    test_genera = [
        ('Potexvirus', 'Potexvirus '),
        ('Cucumovirus', 'Cucumovirus '),
        ('Ilarvirus', 'Ilarvirus '),
        ('Luteovirus', 'Luteovirus '),
        ('Polerovirus', 'Polerovirus '),
        ('Orthotospovirus', 'Orthotospovirus '),
        ('Begomovirus', 'Begomovirus '),
        ('Potyvirus', 'Potyvirus '),
        ('Tobamovirus', 'Tobamovirus '),
        ('Nepovirus', 'Nepovirus '),
        ('Ampelovirus', 'Ampelovirus '),
        'Trichovirus ',
        'Vitivirus ',
        'Crinivirus ',
        'Tospovirus ',
    ]

    print(f"\n=== 关键属覆盖（{len(all_species)} 物种中）===")
    for item in test_genera:
        if isinstance(item, tuple):
            genus, prefix = item
        else:
            genus, prefix = item.rstrip(), item
        n = sum(1 for sp in all_species if sp.startswith(prefix))
        print(f"  {genus:24s} {n:>4d} 种")

    print(f"\n=== 库内全部 Orthotospovirus/Tospovirus 物种 ===")
    ortho = sorted(sp for sp in all_species
                   if 'Orthotospovirus' in sp or 'tospovirus' in sp.lower())
    for sp in ortho:
        print(f"  {sp}")
    print(f"  (共 {len(ortho)} 条)")

    # ── 0 植物宿主覆盖的科：从库数据动态生成（2026-10-03 起不再硬编码名单）──
    zero_plant = sorted(f for f in fam_sp if not fam_plant[f])
    print(f"\n=== 库内无植物宿主物种的科（动态生成，共 {len(zero_plant)} 个）===")
    for fam in zero_plant:
        print(f"  {fam:28s} species={len(fam_sp[fam])}")

    # ── 关键科缺口对照（需外部期望清单；未提供时明确跳过）──
    if args.expect and os.path.isfile(args.expect):
        exp = json.load(open(args.expect, encoding='utf-8'))
        fams = exp.get('families') or []
        print(f"\n=== 关键科缺口对照（期望清单：{args.expect}，{len(fams)} 科）===")
        for fam in fams:
            n = len(fam_sp.get(fam, set()))
            flag = '存在' if n else '缺失(0 覆盖)'
            print(f"  {fam:28s} species={n:<4d} {flag}")
    else:
        print('\n[跳过] 关键科缺口对照未提供期望清单'
              '（--expect <json>，格式 {"families": ["Reoviridae", ...]}）。'
              '\n       原脚本内的 7 个科硬编码名单已于 2026-10-03 移除：'
              '缺失与否应由数据+期望清单判定，不再写死在代码里。')


if __name__ == '__main__':
    main()
