# -*- coding: utf-8 -*-
"""VirPhyKit Example 全目录验证：每一个数据文件都要被对应的工具真实消费并断言。

覆盖矩阵（Example 下 21 个文件，[10] 段做全目录覆盖检查，漏一个即失败）：

  Example 文件                                → 工具          → 断言
  ─────────────────────────────────────────────────────────────────────
  SeqIDRenamer/rename.fasta                   → pdrename      → 16 条头按映射逐条替换、序列不变
  SeqIDRenamer/SeqIDRenamer.txt                 （映射表）
  Geosubsampler/RSV/RSV_209CP.fasta           → pdsub         → 区域普查 + 三模式 + 可复现
  VirSpaceTime/ToMV/Temporal.txt              → pdspacetime   → 解析 + Total=各区行和自洽
  VirSpaceTime/ToMV/Spatial.txt                 （坐标）       → 102 落点全合法、区域与 Temporal 一致
  RRT/PVS/Original/PVS_origin.tre             → pdrrt         → 根状态 SAm=0.3958
  RRT/PVS/Randomized/PVS_random1..20.tre        （20 副本）    → 判定 PASS、Min/Max 与副本一致
  RSPP-Viz/DAT_MCC.tre                        → pdrspp        → Region 六态、MAP=SAm
  RSPP-Viz/Mascot_MCC.tre                       → 同上        → max → SEA=1.0
  RSPP-Viz/MTT_MCC.tre                          → 同上        → type 九态、MAP=Spain
  TempMig/RSV/MCC_Mascot_strict.tre           → pdtempmig     → 与期望矩阵逐格一致（锁死契约）
  TempMig/RSV/Migration_matrix.txt              （期望输出）
  MJRM Generator/PVS/PVS_no_matrix.xml        → pdmjrm        → 逐字节一致（锁死契约）；
  MJRM Generator/PVS/PVS_with_matrix.xml        （期望输出）    状态集合从 XML 自动抽取核对
  TreeTime-RTT/H3N2/fasta_file/h3n2_na_500.fasta → pdtreetime → na_500 全量跑通（见下）
  TreeTime-RTT/H3N2/nwk_file/h3n2_na_500.nwk  → 同上          → 476 叶、速率/R²/置换 p、地区着色
  TreeTime-RTT/H3N2/nwk_file/h3n2_na_200.nwk  → 同上（负例）  → 可解析 198 叶；但菌株不在
  TreeTime-RTT/H3N2/nwk_file/h3n2_na_20.nwk   → 同上（负例）  → na_500 比对/元数据里（2/19、
    13/198 重叠）→「缺日期」「缺序列」两级精确拒绝（Example 未随附其数据，属预期）
  TreeTime-RTT/H3N2/nwk_file/ebola.nwk        → 同上（负例）  → 0 重叠必须被「缺日期」明确拒绝
  TreeTime-RTT/H3N2/metadata_file/*.csv         （元数据）     → 476 行；na_20/200 配它会缺叶日期 → 拒绝
  TreeTime-RTT/Default_mapping.txt              （地点→区划）  → 覆盖 21/27 地点，未覆盖→Unknown 如实呈现
  TreeDater-LTT/nwk_file/h3n2_na_500.nwk      → pdltt         → 476 叶、遗传距离轴 + 诚实提示
  TreeDater-LTT/metadata_file/*.csv             （元数据）

用法: python tests/_check_phylodyn_example.py
退出码: 0 = 全部通过；1 = 有失败
"""
import csv
import glob
import io
import os
import re
import sys

for _s in (sys.stdout, sys.stderr):
    try:
        _s.reconfigure(errors='replace')
    except (AttributeError, OSError):
        pass

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from Virus_Platform_Core import phylodyn_kit as pk      # noqa: E402
from Virus_Platform_Core import phylodyn_trees as pt    # noqa: E402

# 外部参照仓（延伸基因组/MMPV-RNA 工作区）：默认仍是维护者本机路径，
# 可用环境变量 VP_EXT_ROOT 覆盖；EX 不存在时本脚本 [SKIP] 返回 0。
EXT_ROOT = os.environ.get('VP_EXT_ROOT', r'D:\桌面\延伸基因组\MMPV-RNA')
EX = os.path.join(EXT_ROOT, 'biosoft', 'VirPhyKit', 'Example')
WORK = os.path.join(ROOT, 'run', '_check_phylodyn_example')
FAILS = []
CONSUMED = set()


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


def use(path):
    """登记一个 Example 文件被本测试真实消费。"""
    CONSUMED.add(os.path.normpath(path))
    return path


if not os.path.isdir(EX):
    print(f'[SKIP] VirPhyKit Example 不存在: {EX}')
    sys.exit(0)

# ════════════════════════════════════════════════════════════════
print('[1] SeqIDRenamer：16 条头逐条替换、序列内容不变')
fa = use(os.path.join(EX, 'SeqIDRenamer', 'rename.fasta'))
mp = use(os.path.join(EX, 'SeqIDRenamer', 'SeqIDRenamer.txt'))
src_seqs = pk.read_fasta(fa)
ren_map = {}
for ln in pk.read_text(mp).splitlines():
    if ln.strip():
        a, b = ln.split('\t')
        ren_map[a] = b
r = pk.rename_sequences(fa, mp, out('renamed'))
new_seqs = pk.read_fasta(r['out_file'])
chk(list(new_seqs) == [ren_map[h] for h in src_seqs],
    '输出头序列 = 映射表新 ID 按原顺序逐条对应')
chk(all(new_seqs[ren_map[h]] == s for h, s in src_seqs.items()),
    '序列内容逐条未变')

# ════════════════════════════════════════════════════════════════
print('[2] GeoSubsampler：区域普查 + 三模式 + 可复现')
fa = use(os.path.join(EX, 'Geosubsampler', 'RSV', 'RSV_209CP.fasta'))
all_seqs = pk.read_fasta(fa)
census = {}
for h in all_seqs:
    census[pk.region_of_name(h)] = census.get(pk.region_of_name(h), 0) + 1
chk(sum(census.values()) == len(all_seqs) == 209, f'209 条序列全部分区（{census}）')
o = out('sub')
s1 = pk.subsample_fasta(fa, o, mode='random', n=30, seed=42)
s1b = pk.subsample_fasta(fa, o, mode='random', n=30, seed=42)
ex1 = list(pk.read_fasta(s1['out_file']))
chk(len(ex1) == 30 and set(ex1) <= set(all_seqs), '随机抽取 30 条且都来自原文件')
chk(ex1 == list(pk.read_fasta(s1b['out_file'])), '同 seed 两次结果逐条一致')
s2 = pk.subsample_fasta(fa, o, mode='region', region='SWM', n=8, seed=42)
ex2 = pk.read_fasta(s2['extract_file'])
chk(len(ex2) == 8 and all(pk.region_of_name(h) == 'SWM' for h in ex2),
    '定区剔除：抽出的 8 条全是 SWM')
chk(s2['n_remaining'] == 201, '剩余 201 条')
mn = min(census.values())
s3 = pk.subsample_fasta(fa, o, mode='equal', seed=42)
ex3 = pk.read_fasta(s3['out_file'])
cnt3 = {}
for h in ex3:
    cnt3[pk.region_of_name(h)] = cnt3.get(pk.region_of_name(h), 0) + 1
chk(set(cnt3) == set(census) and set(cnt3.values()) == {mn} and len(ex3) == mn * len(census),
    f'均等采样：{len(census)} 区 × 每区 {mn} 条 = {len(ex3)}')
for mode, kw in [('random', {'n': 999}), ('region', {'n': 1, 'region': 'ZZZ'})]:
    try:
        pk.subsample_fasta(fa, o, mode=mode, seed=1, **kw)
        chk(False, f'{mode} 非法参数被拒')
    except ValueError:
        chk(True, f'{mode} 非法参数被拒')

# ════════════════════════════════════════════════════════════════
print('[3] VirSpaceTime：两表解析 + 内部自洽')
tp = use(os.path.join(EX, 'VirSpaceTime', 'ToMV', 'Temporal.txt'))
sp = use(os.path.join(EX, 'VirSpaceTime', 'ToMV', 'Spatial.txt'))
t = pk.parse_temporal(tp)
s = pk.parse_spatial(sp)
chk(all(tot == sum(t['series'][rg][i] for rg in t['regions'])
        for i, tot in enumerate(t['totals'])),
    'Temporal：Total 列 = 各区域列行和（全部行自洽）')
chk(t['years'] == sorted(t['years']), 'Temporal：年份单调不减')
chk(s['n_rows'] == 102 and s['n_bad'] == 0, 'Spatial：102 行、0 坏行')
chk(set(s['regions']) == set(t['regions']),
    f'Spatial 区域集合与 Temporal 列一致（{sorted(s["regions"])}）')

# ════════════════════════════════════════════════════════════════
print('[4] RRT：1 真树 + 20 副本')
orig = use(os.path.join(EX, 'RRT', 'PVS', 'Original', 'PVS_origin.tre'))
rands = [use(p) for p in sorted(glob.glob(os.path.join(EX, 'RRT', 'PVS', 'Randomized', '*.tre')))]
chk(len(rands) == 20, f'20 棵随机副本全部入测（{len(rands)}）')
r = pt.rrt_from_mcc(orig, rands)
chk(abs(r['real_prob'] - 0.39577836411609496) < 1e-9 and r['target_state'] == 'SAm',
    '真树根状态 SAm 后验 = 0.39577836…（与文件注释逐位一致）')
mins, maxs = {}, {}
for rw in r['randomized']:
    for st, p in rw['probs'].items():
        mins[st] = min(mins.get(st, 1), p)
        maxs[st] = max(maxs.get(st, 0), p)
chk(len(mins) == len(r['states']) and len(r['randomized']) == 20,
    f'副本后验表 {len(mins)} 状态 × 20 副本齐全')
chk(abs(mins[r['target_state']] - r['min_rand_target']) < 1e-12
    and abs(maxs[r['target_state']] - r['max_rand_target']) < 1e-12,
    'Min/Max 与 20 棵副本逐棵核对一致')
chk(r['passed'] and r['real_prob'] > r['max_rand_target'],
    f"RRT 判定 PASS（{r['real_prob']:.4f} > {r['max_rand_target']:.4f}）")

# ════════════════════════════════════════════════════════════════
print('[5] RSPP-Viz：三棵 MCC 树')
expect = {
    'DAT_MCC.tre': ('Region', 'SAm', 6),
    'Mascot_MCC.tre': ('max', 'SEA', 1),
    'MTT_MCC.tre': ('type', 'Spain', 9),
}
for fn, (trait, map_state, n_states) in expect.items():
    p = use(os.path.join(EX, 'RSPP-Viz', fn))
    rs = pt.rspp_from_mcc(p)
    chk(rs['trait'] == trait and rs['map_state'] == map_state
        and len(rs['states']) == n_states and abs(rs['prob_sum'] - 1.0) < 1e-6,
        f'{fn}: trait={rs["trait"]} MAP={rs["map_state"]} {len(rs["states"])} 态 后验和={rs["prob_sum"]:.4f}')

# ════════════════════════════════════════════════════════════════
print('[6] TempMig：与期望矩阵逐格一致（锁死契约）')
tree = use(os.path.join(EX, 'TempMig', 'RSV', 'MCC_Mascot_strict.tre'))
expf = use(os.path.join(EX, 'TempMig', 'RSV', 'Migration_matrix.txt'))
s = pt.migration_matrix_from_mcc(tree, out_dir=out('tempmig'))
lines = pk.read_text(expf).splitlines()
hdr = lines[0].split('\t')
exp = {}
for ln in lines[1:]:
    c = ln.split('\t')
    exp[int(c[0])] = {k: int(v) for k, v in zip(hdr[1:], c[1:])}
mine = {y: {d: s['matrix'][d][i] for d in s['directions']}
        for i, y in enumerate(s['years'])}
common = sorted(set(exp) & set(mine))
diffs = sum(1 for y in common for d in set(exp[y]) | set(mine[y])
            if exp[y].get(d, 0) != mine[y].get(d, 0))
chk(len(common) == len(exp) == 107 and len(s['directions']) == 25 and diffs == 0,
    f'矩阵逐格一致：107 年 × 25 方向，差异格 {diffs}')
offd = sum(v for k, v in s['totals'].items() if '_to_' in k
           and k.split('_to_')[0] != k.split('_to_')[1])
diag = sum(v for k, v in s['totals'].items()
           if k.split('_to_')[0] == k.split('_to_')[1])
chk(diag > offd > 0, f'口径 sanity：对角线谱系-年 {diag} > 方向迁移分支-年 {offd} > 0')

# ════════════════════════════════════════════════════════════════
print('[7] MJRM：逐字节一致（锁死契约）+ 状态集合从 XML 自动抽取')
src = use(os.path.join(EX, 'MJRM Generator', 'PVS', 'PVS_no_matrix.xml'))
expf = use(os.path.join(EX, 'MJRM Generator', 'PVS', 'PVS_with_matrix.xml'))
xml_text = pk.read_text(src)
regions_in_xml = sorted(set(re.findall(r'<attr name="Region">\s*(\S+?)\s*</attr>', xml_text)))
chk(regions_in_xml == ['AS', 'EU', 'ME', 'NAm', 'OC', 'SAm'],
    f'从输入 XML 自动抽到 6 个状态（{regions_in_xml}）')
# BEAST 状态序 = Region.count 对角 0 的位置序；状态名按 XML 首现序对应行号
count_val = re.search(r'<parameter id="Region.count" value="([^"]+)"', xml_text)
states_order = None
if count_val:
    vals = count_val.group(1).split()
    zero_pos = [i for i, v in enumerate(vals) if float(v) == 0.0]
    by_first = []
    for mm in re.finditer(r'<attr name="Region">\s*(\S+?)\s*</attr>', xml_text):
        if mm.group(1) not in by_first:
            by_first.append(mm.group(1))
    if len(zero_pos) == 6 and len(by_first) == 6:
        states_order = [by_first[i] for i in range(6)]
chk(states_order == ['AS', 'OC', 'EU', 'ME', 'NAm', 'SAm'],
    f'状态序从 Region.count 对角线 + XML 首现序复原 = BEAST 状态序（{states_order}）')
dst = os.path.join(out('mjrm'), 'mine.xml')
pk.mjrm_insert(states_order, src, dst)
a = open(dst, 'rb').read().replace(b'\r\n', b'\n')
b = open(expf, 'rb').read().replace(b'\r\n', b'\n')
chk(a == b, '插入结果与 PVS_with_matrix.xml 逐字节一致（行尾归一后）')

# ════════════════════════════════════════════════════════════════
print('[8] TreeTime-RTT：na_500 全量 + na_200/na_20 派生元数据 + 两个负例')
EX_TT = os.path.join(EX, 'TreeTime-RTT', 'H3N2')
fasta500 = use(os.path.join(EX_TT, 'fasta_file', 'h3n2_na_500.fasta'))
meta500 = use(os.path.join(EX_TT, 'metadata_file', 'h3n2_na_500.metadata.csv'))
o = out('treetime')
# Default_mapping.txt 在 TreeTime-RTT\ 根目录（不在 H3N2\ 里——放错了会静默
# 回退成原始地点名着色，这里路径必须指对）
mapping = use(os.path.join(EX, 'TreeTime-RTT', 'Default_mapping.txt'))


def tips_of(nwk):
    text = pk.read_text(nwk)
    return set(re.findall(r'[\(,]([^():,\[\]]+):', text))


def derive_meta_from_names(nwk, path):
    """叶名 `...|MM/DD/YYYY|LOC...` / `...|YYYY|...` → name,date,location CSV。
    小数年用平台口径（year + doy/365.25），与 VirPhyKit 元数据同源。"""
    rows = ['name,date,location']
    for tip in sorted(tips_of(nwk)):
        parts = tip.split('|')
        d = parts[2] if len(parts) > 2 else ''
        m = re.fullmatch(r'(\d{2})/(\d{2})/(\d{4})', d)
        mm = re.fullmatch(r'(\d{2})/(\d{4})', d)
        if m:
            dec = pk.date_to_decimal(f'{m.group(3)}-{m.group(1)}-{m.group(2)}')
        elif mm:
            dec = float(mm.group(2)) + (int(mm.group(1)) - 0.5) / 12.0
        elif re.fullmatch(r'\d{4}', d):
            dec = float(d) + 0.5
        else:
            continue
        loc = parts[3] if len(parts) > 3 else ''
        rows.append(f'{tip},{dec},{loc}')
    with io.open(path, 'w', encoding='utf-8', newline='\n') as f:
        f.write('\n'.join(rows) + '\n')
    return path


o = out('treetime')
s = pt.treetime_rtt(fasta500,
                    use(os.path.join(EX_TT, 'nwk_file', 'h3n2_na_500.nwk')),
                    meta500, o, mapping_path=mapping, n_perm=50)
chk(s['n_tips'] == 476, f'na_500：476 叶全部配到日期（{s["n_tips"]}）')
chk(1e-3 < (s['clock_rate'] or 0) < 1e-2,
    f'速率 {s["clock_rate"]:.3e} 在 H3N2 NA 合理带')
chk((s['r2'] or 0) > 0.9 and s['perm_p'] <= 0.05,
    f'R²={s["r2"]:.3f}，date-permutation p={s["perm_p"]:.4f}')
chk(os.path.isfile(os.path.join(o, 'timetree_inferred.nwk'))
    and os.path.isfile(os.path.join(o, 'RootToTip.pdf')), '两件产物齐全')
locs = {p['region'] for p in s['points']}
chk('Unknown' in locs and 5 <= len(locs) <= 12,
    f'映射生效：{len(locs)} 组区划（含 Unknown；Guam/HongKong/Serbia 无映射如实归 Unknown）')

# na_20 / na_200：与 na_500 同命名规范的较小示例树，但它们的菌株**不在**
# na_500 的比对/元数据里（重叠 2/19、13/198）——VirPhyKit 没有随 Example 提供
# 这两棵树的数据。平台必须给出**精确拒绝**（先查日期、再查比对覆盖）。
fasta_seqs = pk.read_fasta(fasta500)
for nwk_name, n_expect in [('h3n2_na_200.nwk', 198), ('h3n2_na_20.nwk', 19)]:
    nwk = use(os.path.join(EX_TT, 'nwk_file', nwk_name))
    tips = tips_of(nwk)
    chk(len(tips) == n_expect, f'{nwk_name}：树可解析，{n_expect} 叶')
    meta_d = derive_meta_from_names(nwk, os.path.join(out('treetime'),
                                                      nwk_name + '.meta.csv'))
    sub = [(h, sq) for h, sq in fasta_seqs.items() if h in tips]
    sub_fa = os.path.join(out('treetime'), nwk_name + '.fasta')
    pk.write_fasta(sub_fa, sub)
    n_no_seq = len(tips) - len(sub)
    try:
        pt.treetime_rtt(sub_fa, nwk, meta_d, out('treetime'))
        chk(False, f'{nwk_name} 数据不全时被拒')
    except ValueError as e:
        chk(f'{n_no_seq} 个树叶的序列' in str(e),
            f'{nwk_name}：比对缺 {n_no_seq} 条序列 → 明确拒绝'
            '（Example 未随附这两棵树的比对，属预期）')

# 负例 1：na_20 配 na_500 的元数据（2/19 重叠）→ 必须精确拒绝
meta_names = {row['name'] for row in csv.DictReader(
    io.open(meta500, encoding='utf-8-sig'))}
bad = tips_of(os.path.join(EX_TT, 'nwk_file', 'h3n2_na_20.nwk')) - meta_names
try:
    pt.treetime_rtt(fasta500,
                    os.path.join(EX_TT, 'nwk_file', 'h3n2_na_20.nwk'), meta500,
                    out('treetime'))
    chk(False, 'na_20 配 500 元数据被拒')
except ValueError as e:
    chk(str(len(bad)) in str(e) and '缺日期' in str(e),
        f'na_20 配 500 元数据 →「{len(bad)} 个叶缺日期」明确拒绝')
# 负例 2：ebola.nwk（0 重叠）→ 同样拒绝
ebola = use(os.path.join(EX_TT, 'nwk_file', 'ebola.nwk'))
try:
    pt.treetime_rtt(fasta500, ebola, meta500, out('treetime'))
    chk(False, 'ebola.nwk 被拒')
except ValueError as e:
    chk('362 个叶缺日期' in str(e), f'ebola.nwk（362 叶 0 重叠）→ 明确拒绝')

# ════════════════════════════════════════════════════════════════
print('[9] TreeDater-LTT：同一棵 h3n2 树')
ltt_nwk = use(os.path.join(EX, 'TreeDater-LTT', 'nwk_file', 'h3n2_na_500.nwk'))
ltt_meta = use(os.path.join(EX, 'TreeDater-LTT', 'metadata_file',
                            'h3n2_na_500.metadata.csv'))
l = pt.ltt_from_tree(ltt_nwk, meta_path=ltt_meta, out_dir=out('ltt'))
chk(l['n_tips'] == 476 and l['axis'] == 'subst_distance',
    f'LTT：{l["n_tips"]} 叶；遗传距离树如实用距离轴（不冒充日历年）')
ns = [p['n'] for p in l['curve']]
chk(ns == sorted(ns) and ns[0] >= 2 and ns[-1] == 476,
    f'LTT 曲线单调不减：根 {ns[0]} 谱系（该树根为 {ns[0]} 叉）→ 现存 476')
chk(l['log_linear_slope'] is not None and os.path.isfile(
    os.path.join(out('ltt'), 'ltt_curve.tsv')), '曲线文件已写出，斜率有限')

# ════════════════════════════════════════════════════════════════
# ════════════════════════════════════════════════════════════════
print('[11] MCC 区域随机化器（RSV MCC，Fitch 口径）')
o = out('pdrand')
rr = pt.mcc_randomize_test(tree, n=20, out_dir=o, seed=7)
chk(rr['n_tips'] == 209 and rr['trait'] == 'max',
    '随机化：209 叶、trait=max（实际 trait=%s, tips=%d）' % (rr['trait'], rr['n_tips']))
chk(rr['steps_real'] < rr['steps_random_min'],
    '真实步数 %d < 随机最小 %d（结构信号强）' % (rr['steps_real'], rr['steps_random_min']))
chk(rr['p_value'] <= 0.05 and rr['significant'],
    'p=%.4f 显著' % rr['p_value'])
chk(len(rr['random_trees']) == 20, '20 棵随机树落盘')
r0 = pt.load_mcc_tree(os.path.join(o, 'random1.tre'))
tips_r = [pt.node_state(n, rr['trait']) for n, _p, _d in pt._walk(r0[0])
          if not n['children']]
import collections as _cl
chk(dict(_cl.Counter(tips_r)) == dict(_cl.Counter(rr['census'])),
    '随机树叶状态多重集与原树守恒（纯置换）')
rt = pt.rrt_from_mcc(tree, sorted(glob.glob(os.path.join(o, '*.tre'))))
chk(rt['n_random'] == 20 and isinstance(rt['passed'], bool),
    '随机树可直接喂 RRT 卡（n=20）')

print('[12] BEAST .trees 后验集（合成 12 棵）')
o = out('posterior')
root, _T = pt.load_mcc_tree(tree)
tree_str = pt.write_beast_newick(root)
trees12 = os.path.join(o, 'twelve.trees')
with io.open(trees12, 'w', encoding='utf-8', newline=chr(10)) as f:
    f.write('#NEXUS' + chr(10) + 'Begin trees;' + chr(10))
    for i in range(12):
        f.write(chr(9) + 'tree t%d = %s;' % (i, tree_str) + chr(10))
    f.write('End;' + chr(10))
rs = pt.rspp_posterior(trees12, burnin=0)
one = pt.rspp_from_mcc(tree)
_one_probs = dict(zip(one['states'], one['probs']))
chk(rs['n_trees'] == 12 and rs['map_consensus'] == one['map_state']
    and abs(rs['mean_probs'][one['map_state']] - _one_probs[one['map_state']]) < 1e-9,
    '后验 RSPP：单树 .trees 与直接解析一致；12 棵均值正确')
r6 = pt.rspp_posterior(trees12, burnin=0.5)
chk(r6['n_trees'] == 6, 'burn-in 0.5 → 只留 6 棵')
tp = pt.tempmig_posterior(trees12, burnin=0.05, out_dir=o)
tm1 = pt.migration_matrix_from_mcc(tree)
_y0, _y1 = tp['years'][0], tp['years'][-1]
_msg = '后验 TempMig：同树均值=单树、sd=0；交集年轴 %d-%d' % (_y0, _y1)
chk(tp['n_trees'] == 12
    and abs(tp['totals_mean']['JAP_to_JAP'] - tm1['totals']['JAP_to_JAP']) < 1e-6
    and tp['totals_sd']['JAP_to_JAP'] == 0, _msg)

print('[13] TreeDater 真引擎（便携 R + treedater；不可用时 SKIP）')
td_ok, td_info = pt.treedater_available()
if not td_ok:
    # 环境依赖段：本机未装便携 R / treedater 时按本段标题口径 SKIP，
    # 不计入失败（其余 [1]-[12]/[14] 段照常守护引擎契约）
    print('[SKIP] TreeDater 真引擎不可用（%s）' % td_info)
else:
    chk(True, '便携 R + treedater 可用')
if td_ok:
    o = out('treedater')
    lt = use(os.path.join(EX, 'TreeDater-LTT', 'nwk_file', 'h3n2_na_500.nwk'))
    lm = use(os.path.join(EX, 'TreeDater-LTT', 'metadata_file',
                          'h3n2_na_500.metadata.csv'))
    td = pt.treedater_ltt(lt, lm, o, seq_len=1400)
    chk(td['engine'] == 'treedater' and td['rate'] is not None,
        'treedater 速率 %s' % (('%.2e' % td['rate']) if td['rate'] else '未解析'))
    # 2026-09-18 修：ML 树多为**未定根**（根三分叉），根到端距离无意义 →
    # R 侧先判 is.rooted()，未定根时以最老样本为外类群定根再回归。
    # 修前此条报 −3.49e−3（定根前的假斜率），修后 +3.8e−3、R²≈0.99。
    chk(1e-3 < (td['rate'] or 0) < 1e-1, '速率在合理带（1e-3, 1e-1）且为正')
    chk((td.get('rtt_r2') or 0) > 0.9,
        'root-to-tip 回归 R²=%.3f > 0.9（时间信号强 ⇒ 定根正确）'
        % (td.get('rtt_r2') or 0))
    chk(td['axis'] == 'calendar_year', '定年树 LTT 横轴 = 日历年')
    chk(td['n_tips'] == 476 and td['curve'][-1]['n'] == 476, '曲线终点 = 476 叶')
    chk(os.path.isfile(os.path.join(o, 'dated_tree.nwk'))
        and os.path.isfile(os.path.join(o, 'Phylogeny.pdf')), '定年树 + PDF 齐全')

print('[10] 全目录覆盖：Example 下每个文件都被上面真实消费')
all_files = {os.path.normpath(p) for p in glob.glob(os.path.join(EX, '**', '*'), recursive=True)
             if os.path.isfile(p)}
missing = sorted(all_files - CONSUMED)
chk(not missing, f'覆盖 {len(CONSUMED)}/{len(all_files)} 个文件'
                 + (f'；未覆盖: {[os.path.relpath(m, EX) for m in missing]}' if missing else ''))

# ════════════════════════════════════════════════════════════════
print()
if FAILS:
    print(f'PHYLODYN EXAMPLE CHECKS FAILED: {len(FAILS)} 项')
    for m in FAILS:
        print('  ✘', m)
    sys.exit(1)
print(f'PHYLODYN EXAMPLE CHECKS PASSED（Example 全部 {len(all_files)} 个文件已验证）')
