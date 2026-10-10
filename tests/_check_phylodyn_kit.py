# -*- coding: utf-8 -*-
"""进化动力学工具组（phylodyn 新卡）后端巡检：Example 数据逐模块验证 + 灵敏度自检。

覆盖 2026-09-17 新增的 VirPhyKit 对齐引擎：
    phylodyn_kit.py   —— SeqIDRenamer / SeqGrouper / GeoSubsampler / VirSpaceTime /
                         MJRM 生成器 / 数据接入三来源（zip、手动、GenBank）
    phylodyn_trees.py —— MCC 解析 / TempMig / RRT / RSPP / BSP / TreeTime / LTT

两个「精确对拍」契约（对齐外部 VirPhyKit 的 Example 产物，改动引擎必须保持）：
    [MJRM]  PVS_no_matrix.xml + 六状态 → 与 PVS_with_matrix.xml 逐字节一致
            （行尾归一后；CRLF 是参照产物的历史格式）
    [TempMig] MCC_Mascot_strict.tre → 年×方向矩阵与 Migration_matrix.txt 逐格一致
第三个契约（2026-09-18 补，源自一次真实故障）：
    [date 口径] 数据接入必须同时吃得下 ISO 与 explorer 导出包的**小数年**
                （`2008.58197` → `2008-08-01` 折算入库）。原实现只认严格 ISO，
                拿自家 explorer 的真导出包（209 条全小数年）导入会整包被拒、
                只抛一句「没有一条序列通过校验」——上下游当场断开。[1]/[10] 段钉住。
本脚本的可信度仿 tests/_check_virphykit_export.py：末段「灵敏度自检」把同一套
断言喂故意做坏的输入，必须逐条报错——证明这些断言「会红」，不是恒真。

VirPhyKit 示例目录不存在时（其他机器），Example 段整体 SKIP，自造 fixture 段
与灵敏度自检段照常执行。

用法: python tests/_check_phylodyn_kit.py   [--fast 跳过 TreeTime 全量跑]
退出码: 0 = 全部通过；1 = 有失败
"""
import io
import os
import re
import shutil
import sys
import tempfile

for _s in (sys.stdout, sys.stderr):
    try:
        _s.reconfigure(errors='replace')
    except (AttributeError, OSError):
        pass

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from Virus_Platform_Core import phylodyn_kit as pk  # noqa: E402
from Virus_Platform_Core import phylodyn_trees as pt  # noqa: E402

FAILS = []
# 外部参照仓（延伸基因组/MMPV-RNA 工作区）：默认仍是维护者本机路径，
# 可用环境变量 VP_EXT_ROOT 覆盖（指向 MMPV-RNA 工作区根）；
# 目录不存在时相关小节自动跳过。
EXT_ROOT = os.environ.get('VP_EXT_ROOT', r'D:\桌面\延伸基因组\MMPV-RNA')
EXAMPLE = os.path.join(EXT_ROOT, 'biosoft', 'VirPhyKit', 'Example')
WORK = os.path.join(ROOT, 'run', '_check_phylodyn_kit')
FAST = '--fast' in sys.argv


def chk(cond, msg):
    if cond:
        print('  \u2714 ' + msg)
    else:
        print('  \u2718 ' + msg)
        FAILS.append(msg)
    return bool(cond)


def out(name):
    d = os.path.join(WORK, name)
    os.makedirs(d, exist_ok=True)
    return d


# ════════════════════════════════════════════════════════════════
print('[1] 通用件：日期校验 / 小数年换算')
ok, v, err = pk.validate_date_str('2019-03-01')
chk(ok and v == '2019-03-01', '完整日期通过')
ok, v, _ = pk.validate_date_str('2019-07')
chk(ok and v == '2019-07', '年-月通过')
ok, v, _ = pk.validate_date_str('2019')
chk(ok and v == '2019', '年仅通过')
chk(not pk.validate_date_str('2019.5')[0], '小数年被拒（默认仍是严格 ISO 口径）')
ok, v, _ = pk.validate_date_str('2019.5', allow_decimal=True)
chk(ok and v == '2019.5', 'allow_decimal=True 时小数年被接受（explorer 导出包默认口径）')
chk(not pk.validate_date_str('03/2019')[0], '斜杠格式被拒')
chk(not pk.validate_date_str('2019-13-01')[0], '13 月被拒')
chk(not pk.validate_date_str('2019-02-30')[0], '2 月 30 被拒')
chk(not pk.validate_date_str('19.5', allow_decimal=True)[0], '非四位年份的小数仍被拒')
chk(abs(pk.date_to_decimal('2019-01-01') - 2019.0) < 0.003, '1 月 1 日 ≈ 2019.0')
chk(abs(pk.date_to_decimal('2019-12-31') - 2020.0) < 0.01, '12 月 31 日 ≈ 年末')
chk(abs(pk.date_to_decimal('2007-07-27') - 2007.5695) < 0.001, '小数年与 TreeTime 元数据口径一致')
# 逆换算：explorer 导出包的 209 条小数年靠它入库（2026-09-18 修的真实断点）
chk(pk.decimal_to_iso('2008.58197') == '2008-08-01', '小数年 → ISO（2008.58197 → 2008-08-01）')
chk(abs(pk.date_to_decimal(pk.decimal_to_iso('2007.569473')) - 2007.569473) < 1e-6,
    'decimal_to_iso 是 date_to_decimal 的精确逆（往返误差 < 1e-6 年）')
chk(pk.decimal_to_iso('2019.5') == '2019-07-02', '小数年 .5 折成年中（2019-07-02）')

# ════════════════════════════════════════════════════════════════
print('[2] SeqIDRenamer / GeoSubsampler / VirSpaceTime（Example）')
if os.path.isdir(EXAMPLE):
    o = out('rename')
    r = pk.rename_sequences(os.path.join(EXAMPLE, 'SeqIDRenamer', 'rename.fasta'),
                            os.path.join(EXAMPLE, 'SeqIDRenamer', 'SeqIDRenamer.txt'), o)
    seqs = pk.read_fasta(r['out_file'])
    chk(r['n_seq'] == 16 and r['n_renamed'] == 16 and r['n_unmapped'] == 0,
        f"重命名 16/16（实际 {r['n_seq']}/{r['n_renamed']}）")
    chk(list(seqs)[0] == 'Seq01_trait01_date01', '首个新 ID 与 VirPhyKit 契约一致')

    o = out('sub')
    fa = os.path.join(EXAMPLE, 'Geosubsampler', 'RSV', 'RSV_209CP.fasta')
    s1 = pk.subsample_fasta(fa, o, mode='random', n=20, seed=7)
    chk(s1['n_picked'] == 20, '随机抽取 20 条')
    chk(len(pk.read_fasta(s1['out_file'])) == 20, 'extract.fas 恰 20 条')
    s1b = pk.subsample_fasta(fa, o, mode='random', n=20, seed=7)
    chk(s1b['n_picked'] == 20 and
        list(pk.read_fasta(s1b['out_file'])) == list(pk.read_fasta(s1['out_file'])),
        '同 seed 结果可复现')
    s2 = pk.subsample_fasta(fa, o, mode='equal', seed=7)
    chk(min(s2['per_region_sampled'].values()) ==
        min(s2['per_region'].values()) and s2['n_picked'] == 5 * min(s2['per_region'].values()),
        f"均等采样 = 5 区 × 最小区 {min(s2['per_region'].values())} 条")
    s3 = pk.subsample_fasta(fa, o, mode='region', region='CCD', n=5, seed=7)
    chk(s3['n_removed'] == 5 and s3['n_remaining'] == 204, '定区剔除 5 条，余 204')
    try:
        pk.subsample_fasta(fa, o, mode='region', region='CCD', n=999, seed=7)
        chk(False, '超额剔除被拒')
    except ValueError:
        chk(True, '超额剔除被拒')

    t = pk.parse_temporal(os.path.join(EXAMPLE, 'VirSpaceTime', 'ToMV', 'Temporal.txt'))
    sp = pk.parse_spatial(os.path.join(EXAMPLE, 'VirSpaceTime', 'ToMV', 'Spatial.txt'))
    chk(t['years'][0] == 1975 and t['years'][-1] == 2018 and len(t['regions']) == 3,
        'Temporal：1975-2018，3 个区域')
    chk(sp['n_rows'] == 102 and sp['n_bad'] == 0 and sp['regions'] ==
        ['East Asia', 'Europe', 'Middle East'], 'Spatial：102 行落点，0 坏行')
else:
    print('  \u26a0 VirPhyKit Example 目录不存在，跳过 Example 段（自造 fixture 照常）')

# ════════════════════════════════════════════════════════════════
print('[3] SeqGrouper（自造 GenBank + 映射表，含 united states 特判）')
o = out('group')
gb_text = []
for i, (acc, iso, host, geo, date) in enumerate([
        ('LC001', 'IsoA', 'Solanum lycopersicum', 'China: Wuhan', '01-MAR-2019'),
        ('LC002', 'IsoB', 'Solanum lycopersicum', 'USA: Florida', '2019-07-15'),
        ('LC003', 'IsoC', 'Nicotiana tabacum', 'Kenya', '2020'),
], 1):
    gb_text.append(f"""LOCUS       {acc}                  8 bp    DNA     linear   VRL {date}
DEFINITION  Test virus isolate {iso}.
ACCESSION   {acc}
FEATURES             Location/Qualifiers
     source          1..8
                     /organism="Test virus"
                     /isolate="{iso}"
                     /host="{host}"
                     /country="{geo}"
                     /collection_date="{date}"
ORIGIN
        1 acgtacgt
//
""")
gb_path = os.path.join(o, 'in.gb')
with open(gb_path, 'w', encoding='utf-8', newline='\n') as f:
    f.write(''.join(gb_text))
mapping_path = os.path.join(o, 'map.tsv')
with open(mapping_path, 'w', encoding='utf-8', newline='\n') as f:
    f.write('Asia\tchina\nAmericas\tusa\nAmericas\tbrazil\n')

from Virus_Platform_Core import gb_collection as gbc  # noqa: E402
records = gbc.parse_flatfile(pk.read_text(gb_path))
rows = pk.gb_record_rows(records)
chk(len(rows) == 3 and rows[1]['Geo Location'] == 'USA', '记录表 3 行，USA 冒号前缀已剥')
chk(rows[0]['Collection Date'] == '2019-03-01' and rows[2]['Collection Date'] == '2020',
    'GenBank 日期归一（01-MAR-2019 → 2019-03-01）')
gmap = pk.load_group_mapping(mapping_path)
rows2, counts, unmatched = pk.apply_grouping(rows, gmap=gmap)
chk([r['Group'] for r in rows2] == ['Asia', 'Americas', 'Unknown'],
    '分组：china→Asia、united states 特判→Americas、Kenya→Unknown')
chk(list(unmatched) == ['Kenya'], '未匹配值如实回收')

# ════════════════════════════════════════════════════════════════
print('[4] MJRM 生成器（Example 逐字节对拍 —— 锁死契约）')
if os.path.isdir(EXAMPLE):
    o = out('mjrm')
    src = os.path.join(EXAMPLE, 'MJRM Generator', 'PVS', 'PVS_no_matrix.xml')
    exp = os.path.join(EXAMPLE, 'MJRM Generator', 'PVS', 'PVS_with_matrix.xml')
    dst = os.path.join(o, 'mine.xml')
    pk.mjrm_insert(['AS', 'OC', 'EU', 'ME', 'NAm', 'SAm'], src, dst)
    a = open(dst, 'rb').read().replace(b'\r\n', b'\n')
    b = open(exp, 'rb').read().replace(b'\r\n', b'\n')
    chk(a == b, '插入结果与 PVS_with_matrix.xml 逐字节一致（行尾归一后）')
    mat, rewards = pk.mjrm_blocks(['A', 'B', 'C'])
    chk(mat.count('<parameter id=') == 6 and '<rewards>' in rewards
        and rewards.count('_reward') == 3, '纯片段模式：6 方向 + 3 rewards')
try:
    pk.mjrm_blocks(['OnlyOne'])
    chk(False, '少于 2 个状态被拒')
except ValueError:
    chk(True, '少于 2 个状态被拒')
_bad = os.path.join(out('mjrm'), 'no_marker.xml')
with open(_bad, 'w', encoding='utf-8', newline='\n') as f:
    f.write('<beast><taxa id="taxa"/></beast>\n')
try:
    pk.mjrm_insert(['A', 'B'], _bad, os.path.join(out('mjrm'), 'x.xml'))
    chk(False, '无插入标记的 XML 被拒')
except ValueError:
    chk(True, '无插入标记的 XML 被拒')

# ════════════════════════════════════════════════════════════════
print('[5] TempMig 迁移矩阵（Example 逐格对拍 —— 锁死契约）')
if os.path.isdir(EXAMPLE):
    o = out('tempmig')
    s = pt.migration_matrix_from_mcc(
        os.path.join(EXAMPLE, 'TempMig', 'RSV', 'MCC_Mascot_strict.tre'), out_dir=o)
    exp_lines = open(os.path.join(EXAMPLE, 'TempMig', 'RSV',
                                  'Migration_matrix.txt')).read().splitlines()
    hdr = exp_lines[0].split('\t')
    exp = {}
    for ln in exp_lines[1:]:
        c = ln.split('\t')
        exp[int(c[0])] = {k: int(v) for k, v in zip(hdr[1:], c[1:])}
    chk(s['root_year'] == 1907 and s['mostcurrent_year'] == 2013,
        f"年轴 1907-2013（实际 {s['root_year']}-{s['mostcurrent_year']}）")
    mine = {y: {d: s['matrix'][d][i] for d in s['directions']}
            for i, y in enumerate(s['years'])}
    common = sorted(set(exp) & set(mine))
    diffs = [(y, d, exp[y].get(d), mine[y].get(d))
             for y in common for d in set(exp[y]) | set(mine[y])
             if exp[y].get(d, 0) != mine[y].get(d, 0)]
    chk(len(common) == len(exp) == 107 and not diffs,
        f'矩阵逐格一致（107 年 × {len(s["directions"])} 方向，差异格 {len(diffs)}）')
    chk(os.path.isfile(s['out_file']) and
        open(s['out_file']).read().splitlines()[0].startswith('Year\t'),
        'TransposedMatrix.txt 已写出（Year 首列）')
    chk(s['trait'] == 'max', 'BEAST2 max.set 性状自动探测')
    # Translate 块解析（叶名是数字编号的 NEXUS）
    root, T = pt.load_mcc_tree(os.path.join(EXAMPLE, 'TempMig', 'RSV',
                                            'MCC_Mascot_strict.tre'))
    tips = [n['name'] for n, _p, _d in pt._walk(root) if not n['children']]
    chk(any('_' in x for x in tips[:5]), 'Translate 块已还原真实叶名')

# ════════════════════════════════════════════════════════════════
print('[6] RRT / RSPP（Example）')
if os.path.isdir(EXAMPLE):
    import glob
    o = out('rrt')
    orig = os.path.join(EXAMPLE, 'RRT', 'PVS', 'Original', 'PVS_origin.tre')
    rands = sorted(glob.glob(os.path.join(EXAMPLE, 'RRT', 'PVS', 'Randomized', '*.tre')))
    r = pt.rrt_from_mcc(orig, rands)
    chk(r['n_random'] == 20, f"20 棵随机副本（实际 {r['n_random']}）")
    chk(r['target_state'] == 'SAm' and abs(r['real_prob'] - 0.3958) < 0.001,
        '真数据根状态 SAm 后验 ≈ 0.3958')
    chk(r['passed'] is True and '通过' in r['verdict'], 'RRT 判定 PASS（0.3958 > 0.3150）')
    for f_ in ('DAT_MCC.tre', 'Mascot_MCC.tre', 'MTT_MCC.tre'):
        s = pt.rspp_from_mcc(os.path.join(EXAMPLE, 'RSPP-Viz', f_))
        chk(abs(s['prob_sum'] - 1.0) < 1e-6, f'{f_}: 后验和 ≈ 1')
    s = pt.rspp_from_mcc(os.path.join(EXAMPLE, 'RSPP-Viz', 'DAT_MCC.tre'))
    chk(s['trait'] == 'Region' and s['map_state'] == 'SAm' and len(s['states']) == 6,
        'DAT_MCC：Region 六状态，MAP=SAm')

# ════════════════════════════════════════════════════════════════
print('[7] BSP-Viz（合成 TSV + 合成 BEAST log）')
o = out('bsp')
tsv = os.path.join(o, 'sky.tsv')
with open(tsv, 'w', encoding='utf-8', newline='\n') as f:
    f.write('Time\tMedian\tLower\tUpper\n')
    t_ = 1990.0
    while t_ <= 2020:
        med = 100 * 1.1 ** (2020 - t_)
        f.write('%.1f\t%.2f\t%.2f\t%.2f\n' % (t_, med, med * 0.6, med * 1.5))
        t_ += 0.5
b = pt.bsp_from_tsv(tsv)
chk(b['n'] == 61 and not b['approx'] and abs(b['data']['median'][0] - 1744.94) < 1,
    'TSV 精确口径：61 区间，中位数首点正确')
import random  # noqa: E402
random.seed(0)
lg = os.path.join(o, 'sky.log')
with open(lg, 'w', encoding='utf-8', newline='\n') as f:
    cols = (['Sample'] + [f'bayesianSkyline.popSizes{i + 1}' for i in range(3)]
            + [f'bayesianSkyline.groupSizes{i + 1}' for i in range(3)]
            + ['treeModel.rootHeight'])
    f.write('\t'.join(cols) + '\n')
    for i in range(300):
        ps = [random.uniform(20, 200) for _ in range(3)]
        f.write('\t'.join([str(i)] + ['%.2f' % x for x in ps]
                          + ['10', '10', '10', '2010.0']) + '\n')
b2 = pt.bsp_from_beast_log(lg, present_year=2020)
chk(b2['approx'] and b2['n_samples'] == 300 and b2['n'] >= 100,
    'log 近似口径：300 样本，网格 ≥100 点')
_bad = os.path.join(o, 'bad.tsv')
with open(_bad, 'w', encoding='utf-8', newline='\n') as f:
    f.write('A\tB\n1\t2\n')
try:
    pt.bsp_from_tsv(_bad)
    chk(False, '缺列 TSV 被拒')
except ValueError:
    chk(True, '缺列 TSV 被拒')

# ════════════════════════════════════════════════════════════════
print('[8] TreeTime-RTT（Example H3N2 全量；--fast 跳过）')
if os.path.isdir(EXAMPLE) and not FAST:
    o = out('treetime')
    s = pt.treetime_rtt(
        os.path.join(EXAMPLE, 'TreeTime-RTT', 'H3N2', 'fasta_file', 'h3n2_na_500.fasta'),
        os.path.join(EXAMPLE, 'TreeTime-RTT', 'H3N2', 'nwk_file', 'h3n2_na_500.nwk'),
        os.path.join(EXAMPLE, 'TreeTime-RTT', 'H3N2', 'metadata_file',
                     'h3n2_na_500.metadata.csv'),
        o, mapping_path=os.path.join(EXAMPLE, 'TreeTime-RTT', 'Default_mapping.txt'),
        n_perm=50)
    chk(400 <= s['n_tips'] <= 500, f"叶数 {s['n_tips']} 在 400-500")
    chk(1e-3 < (s['clock_rate'] or 0) < 1e-2,
        f"速率 {s['clock_rate']:.2e} 落在 H3N2 NA 合理带 (1e-3, 1e-2)")
    chk((s['r2'] or 0) > 0.9, f"R² {s['r2']:.3f} > 0.9（该示例时间信号极强）")
    chk(s['perm_p'] <= 0.05, f"date-permutation p={s['perm_p']:.3f} ≤ 0.05")
    chk(os.path.isfile(os.path.join(o, 'timetree_inferred.nwk')) and
        os.path.isfile(os.path.join(o, 'RootToTip.pdf')) and
        os.path.isfile(os.path.join(o, 'rtt_points.tsv')), '三件产物齐全')
    chk((s['n_regions'] or 0) >= 5, f"地区着色 {s['n_regions']} 组（Mapping 表生效）")
elif FAST:
    print('  \u26a0 --fast：跳过 TreeTime 全量跑')
else:
    print('  \u26a0 Example 不存在，跳过')

# ════════════════════════════════════════════════════════════════
print('[9] LTT（自研口径：时间树日历年轴 + 遗传距离轴检测）')
o = out('ltt')
tt_nwk = os.path.join(o, 'tt.nwk')
with open(tt_nwk, 'w', encoding='utf-8', newline='\n') as f:
    f.write('((A_2010.5:5,(B_2008.0:8,C_2008.0:4):6):7,(D_2012.0:2,E_2012.0:2):15);')
meta = os.path.join(o, 'm.csv')
with open(meta, 'w', encoding='utf-8', newline='\n') as f:
    f.write('name,date\nA_2010.5,2010.5\nB_2008.0,2008.0\nC_2008.0,2008.0\n'
            'D_2012.0,2012.0\nE_2012.0,2012.0\n')
l = pt.ltt_from_tree(tt_nwk, meta_path=meta, out_dir=o)
chk(l['axis'] == 'calendar_year' and abs(l['curve'][0]['x'] - 1991.0) < 1e-6,
    '时间树 → 日历年轴（根 1991 = 2012 - 21）')
chk([(p['x'], p['n']) for p in l['curve'][:4]] ==
    [(1991.0, 2), (1998.0, 3), (2004.0, 4), (2006.0, 5)], 'LTT 阶梯事件正确')
gen_nwk = os.path.join(o, 'gen.nwk')
if os.path.isdir(EXAMPLE):
    l2 = pt.ltt_from_tree(os.path.join(EXAMPLE, 'TreeDater-LTT', 'nwk_file',
                                       'h3n2_na_500.nwk'),
                          meta_path=os.path.join(EXAMPLE, 'TreeDater-LTT',
                                                 'metadata_file',
                                                 'h3n2_na_500.metadata.csv'))
    chk(l2['axis'] == 'subst_distance' and l2['n_tips'] == 476,
        '遗传距离树被识别（不冒充日历年轴），476 叶')

# ════════════════════════════════════════════════════════════════
print('[10] 数据接入三来源 + 时间地点同检')
o = out('prep')
fasta = '>SeqA\nACGTACGTACGT\n>SeqB\nTTTTGGGGCCCC\n>SeqC\nAAAACCCCGGGG\n'
# 2026-09-18 口径变更：数据接入**接受 explorer 导出包的小数年**（折算成 ISO 入库），
# 所以「坏行」只剩真正非法的日期；SeqB 的小数年改由下方正向断言验证折算结果。
meta_mixed = ('name\tdate\tlocation\tlat\tlon\n'
              'SeqA\t2019-03-01\tWuhan\t30.5\t114.3\n'
              'SeqB\t2019.5\tWuhan\t30.5\t114.3\n'
              'SeqC\t2020-13-01\tNairobi\n')
try:
    pk.import_manual(fasta, meta_mixed, o, strict=True)
    chk(False, '严格模式对坏行中止')
except ValueError as e:
    chk('第 4 行' in str(e) and '月份 13' in str(e) and '共 1 行不合格' in str(e),
        f'严格模式只认非法日期那 1 行（SeqC 13 月），小数年不再算坏行（{str(e)[:42]}…）')
r = pk.import_manual(fasta, meta_mixed, o, strict=False)
chk(r['n_kept'] == 2 and [i['name'] for i in r['issues']] == ['SeqC'],
    '宽松模式保留 SeqA+SeqB（小数年被折算而非剔除），只剔 SeqC')
with open(r['metadata'], encoding='utf-8-sig') as fh:
    _rows_txt = [ln.strip() for ln in fh if ln.strip()]
chk(any(x.startswith('SeqB,2019-07-02,') for x in _rows_txt),
    f'小数年 2019.5 落盘已折成 ISO（末行 = {_rows_txt[-1]}）')
chk('location' in meta_mixed and r['summary']['locations'] == ['Wuhan'],
    '保留行汇总正确（地点 Wuhan）')
# 缺地点（时间地点必须同时给出）
meta_noloc = 'name\tdate\tlocation\nSeqA\t2019-03-01\t\n'
try:
    pk.import_manual(fasta, meta_noloc, o, strict=True)
    chk(False, '缺地点被拒')
except ValueError as e:
    chk('缺少地点' in str(e), '缺地点被拒且提示明确')
# FASTA ↔ 元数据互查
meta_orphan = 'name\tdate\tlocation\nSeqA\t2019-03-01\tWuhan\nSeqZ\t2019-03-02\tRome\n'
r2 = pk.import_manual(fasta, meta_orphan, o, strict=False)
chk(r2['n_kept'] == 1 and any('SeqZ' in (i.get('name') or '') for i in r2['issues']),
    '元数据孤儿行被点名，FASTA↔元数据配对收紧')
# 文件对来源：explorer 导出的「FASTA + 元数据」两份文件
# （zip 那条 2026-09-18 已整条下线 —— 卡上的来源与 `/api/phylodyn/import_export_zip` 都删了）
# ⚠️ explorer 的 date 列**可能是小数年**（2008.58197 / 2019.5），所以必须混入小数年
#    ——2026-09-18 的真实故障就是拿 209 条小数年的导出整包被拒（当时只认 ISO）。
#    见 [1] 段的逆换算断言。
fa3 = os.path.join(o, 'pair.fasta')
mt3 = os.path.join(o, 'pair.csv')
with open(fa3, 'w', encoding='utf-8', newline='\n') as f:
    f.write('>SeqA\nACGTACGTACGT\n>SeqB\nTTTTGGGGCCCC\n>SeqC\nAAAACCCCGGGG\n')
with open(mt3, 'w', encoding='utf-8', newline='\n') as f:
    f.write('name,date,location\nSeqA,2019-03-01,Wuhan\n'
            'SeqB,2019-07,Nairobi\nSeqC,1986.76813,Kunming\n')
r3 = pk.import_from_files(fa3, mt3, o)
chk(r3['n_kept'] == 3 and not r3['issues'],
    '文件对来源：YYYY / YYYY-MM / 小数年 三类 date 全收，issues 为空')
with open(r3['metadata'], encoding='utf-8-sig') as fh:
    _z_rows = [ln.strip() for ln in fh if ln.strip()]
chk(any(x.startswith('SeqC,1986-10-08,') for x in _z_rows),
    f'小数年落盘后统一为 ISO（{_z_rows[-1]}）')
chk(os.path.basename(r3['fasta']) == 'sequences.fasta'
    and os.path.basename(r3['metadata']) == 'metadata.csv', '标准数据集两个文件名固定')
# GenBank 本地来源：一条完整 + 一条 collection_date 是区间值（不合规必须被剔）
gb1 = os.path.join(o, 'x.gb')
gb_nodate = gb_text[1].replace('LC002', 'LC002x').replace(
    '/collection_date="2019-07-15"', '/collection_date="2019-2020"')
with open(gb1, 'w', encoding='utf-8', newline='\n') as f:
    f.write(gb_text[0] + gb_nodate)
r4 = pk.import_from_genbank(o, files=[gb1])
chk(r4['n_kept'] == 1 and any('LC002x' in (i.get('name') or '') for i in r4['issues']),
    f"GenBank 来源：不合规 collection_date 被剔（保留 {r4['n_kept']}）")
# 2026-09-18：zip 来源整条下线（explorer 不再导出 zip）→ 底层函数 `import_from_zip`
# 也按用户要求删掉了。这里改成**防回潮**断言（原来这段测的是"非 zip 被拒"）。
chk(not hasattr(pk, 'import_from_zip'),
    'import_from_zip 已从 kit 移除（zip 来源下线，2026-09-18）')

# ════════════════════════════════════════════════════════════════
print('[11] 灵敏度自检：坏输入必须让断言报错（防恒绿）')
o = out('sens')
bad_tree = os.path.join(o, 'noann.tre')
with open(bad_tree, 'w', encoding='utf-8', newline='\n') as f:
    f.write('((A:1,B:1):1,C:1);')
try:
    pt.rspp_from_mcc(bad_tree)
    chk(False, '无注释树 RSPP 被拒')
except ValueError:
    chk(True, '无注释树 RSPP 被拒')
try:
    pt.migration_matrix_from_mcc(bad_tree)
    chk(False, '无性状树 TempMig 被拒')
except ValueError:
    chk(True, '无性状树 TempMig 被拒')
bad_dates = os.path.join(o, 'bad_dates.csv')
with open(bad_dates, 'w', encoding='utf-8', newline='\n') as f:
    f.write('name,date\nA,not-a-date\nB,also-bad\nC,2019-99-99\n')
try:
    pt._read_dates_meta(bad_dates)
    chk(False, '全坏日期元数据被拒')
except ValueError:
    chk(True, '全坏日期元数据被拒')
try:
    pk._sniff_table_text('name,date\n', '空表')
    chk(False, '空数据表被拒')
except ValueError:
    chk(True, '空数据表被拒')

# ════════════════════════════════════════════════════════════════
print()
if FAILS:
    print(f'PHYLODYN KIT CHECKS FAILED: {len(FAILS)} 项')
    for m in FAILS:
        print('  ✘', m)
    sys.exit(1)
print('PHYLODYN KIT CHECKS PASSED')
