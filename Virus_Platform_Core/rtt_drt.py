# -*- coding: utf-8 -*-
r"""DRT（Date-Randomization Test，日期随机化检验）—— 时间信号的**显著性检验**。

## 为什么必须有它（`#t-rtt` 的 AUDIT 缺陷）

平台原先只报 `phylogeo.root_to_tip()` 的 **最大** R²，再拿固定阈值 0.5 判断
"有没有时间信号"。但"取最大"等价于在 ~2n−3 个候选根位上取**极值** —— 零假设
（日期与拓扑无关）下这个最大值同样会被抬高。AUDIT 实测放大 4.4–5.4 倍，
按 0.5 判会误判：8 叶 **23.5%** / 12 叶 7.0% / 20 叶 0.5% / 30 叶 0%。
而 8–20 叶正是病毒数据集的家常便饭（GCVA 32 条、分区段更少）。

DRT 的解法**不是换阈值，而是造零分布**：把日期在 tip 之间随机置换 N 次，
每次用同一个统计量重算，看真实值在零分布中的**百分位**（≥95% 才通过）。
阈值是数据自己的分位数，不再依赖叶数。

## 来源与口径

移植自参照管道两版并存实现，口径逐条对齐：

* ``virome_phylo_pipeline/utils/temporal_signal.py::date_randomization_test``
  （Duchêne et al. 2015, *Mol Biol Evol* 32:1895–1906）—— 发布版管道调用它，
  内部走 TreeTime 且**每次 reroot(root='least-squares')**，即重根参与零分布。
* ``virome_phylo_pipeline/phylo_results/_drt_tt3.py`` —— 后来的刻意变体：
  *"定根: 在原始拓扑 + 原始根下做比较, 不做 optimal_reroot 重根 …… 固定根同时
  保证 DRT 只改变日期映射这一个变量, 语义干净"*。

两者都提供（`rooting=`），因为**统计含义不同**且都自洽：

  · ``'optimize'``（默认）：每轮（含真实数据）各自取 R² 最大的根位。
    零分布 = "随机日期下能达到的最好 R²" → 与真实统计量同构，比较公平。
    这是**发布版管道的行为**。
  · ``'asis'``：全部轮次固定用树**自带**的根位，不搜索。零分布 = "给定根位下
    随机日期的 R²"。语义最干净（只动日期一个变量，与 `_drt_tt3.py` 一致），
    但树是无根的 NJ/FastTree 产物时根位本身任意 → 可能低估（假阴性）。
  · ``'optimize_once'``：在真实数据上挑一次最优根位，**冻结**给所有轮次。
    ⚠️ **反保守** —— 真实值被最大化、零分布没有 → 通过率被抬高。
    仅供诊断/复现历史编号用，结果里会带 ``warnings``。

## 依赖

只用 ``numpy``（平台已依赖）做置换与分位数；**不引入 TreeTime / scipy**。
置换流用 ``np.random.default_rng(seed).permutation(n)`` —— 与参照管道**同一个
整数置换序列**，所以同一份数据、同一个种子下，第 i 次随机化的日期映射与管道
同构（便于逐轮对拍）。
"""

from __future__ import annotations

import math
import os

from .phylogeo import resolve_years

# 默认 999 次置换、种子 42 固定（审查 P2：原 20 次的 p 值分辨率下限
# 1/21≈0.048，"显著与否"几乎由最极端一次结果决定；999 次给到 1/1000
# 分辨率，DRT 置换为轻量计算，代价可接受）
DEFAULT_N_PERM = 999
DEFAULT_SEED = 42
# DRT 回归至少要有的带日期叶数（与 root_to_tip 的内部门槛一致）
DEFAULT_MIN_TIPS = 5

ROOTINGS = ('optimize', 'asis', 'optimize_once')

# 判据显著性水平（与 conclusion / passed 一致；写进结论字段供下游引用）
DRT_ALPHA = 0.05


# ---------------------------------------------------------------- 树几何

def tree_geometry(tree_path):
    """把树的**几何**一次算好：每个候选根位 → {叶名: 根到该叶距离}。

    返回 ``{'branches': [(key, {叶名: 距离}), …], 'root': {叶名: 距离},
    'names': [叶名…], 'n_branches': int}``

    * ``branches`` 的每个候选＝把**某条枝的中点**当根（与原 ``root_to_tip``
      的搜索空间逐字一致，key 是该枝两端节点 id 排序后的二元组）。
    * ``root`` ＝树**自带**根位（不重根），供 ``rooting='asis'``。

    ⚠️ 耗时部分是遍历求距离（O(枝数 × 叶数)）。DRT 要跑 N 次置换，而**日期置换
    只改年份映射、不动几何** —— 所以这里一次算好，置换循环里只做回归算术，
    把 N 次遍历降成 1 次。这是本模块与"照抄原实现"最大的工程差别。
    """
    if not os.path.isfile(tree_path):
        raise RuntimeError(f'树文件不存在：{tree_path}')
    with open(tree_path, encoding='utf-8', errors='replace') as f:
        text = f.read()
    # ⚠️ 必须用 `treetime_ml.parse_annotated_newick` 而非 `phylogeo.parse_newick`：
    #    后者不认 `[&…]` 注释，`(A:1[&height=0.5],B:2)` 会把 A 的枝长吞成 0.0
    #    （实测）→ 全部根到叶距离错、DRT 结果无意义。TreeTime 的
    #    divergence/timetree.nexus 树**总是带注释**，是这里的常规输入。
    #    两者的节点 schema 同构（name/length/children），可直接替换。
    from .treetime_ml import parse_annotated_newick
    root = parse_annotated_newick(text)

    adj = {}
    names = {}

    def _reg(n):
        i = id(n)
        names[i] = n['name'] or ''
        adj.setdefault(i, [])
        return i

    def _walk(n):
        i = _reg(n)
        for c in n['children']:
            ci = _walk(c)
            ln = c['length'] or 0.0
            adj[i].append((ci, ln))
            adj[ci].append((i, ln))
        return i

    ri = _walk(root)

    def _tip_dists(a, block):
        out = {}
        stack = [(a, 0.0, block)]
        while stack:
            cur, d, frm = stack.pop()
            nb = [(x, l) for x, l in adj[cur] if x != frm]
            if not nb:
                out[cur] = d
                continue
            for x, l in nb:
                stack.append((x, d + l, cur))
        return out

    leaves = [i for i, nb in adj.items() if len(nb) == 1]
    if not leaves:
        raise RuntimeError('树里没有叶节点（度数＝1 的节点）：无法做根到尾回归')
    name_of = {i: names.get(i, '') for i in adj}
    leaf_names = [name_of[l] for l in leaves]

    branches = []
    seen = set()
    for a in list(adj):
        for b, ln in adj[a]:
            key = (min(a, b), max(a, b))
            if key in seen:
                continue
            seen.add(key)
            half = (ln or 0.0) / 2.0
            da = _tip_dists(a, b)
            db = _tip_dists(b, a)
            dists = {}
            for leaf in leaves:
                nm = name_of.get(leaf, '')
                if not nm:
                    continue
                d0 = da.get(leaf)
                if d0 is None:
                    d0 = db.get(leaf)
                    if d0 is None:
                        continue
                dists[nm] = d0 + half
            branches.append((key, dists))

    root_dists = _tip_dists(ri, None)
    root_pos = {name_of[l]: d for l, d in root_dists.items() if name_of.get(l)}
    return {'branches': branches, 'root': root_pos, 'names': leaf_names,
            'n_branches': len(branches)}


# ---------------------------------------------------------------- 回归

def rtt_fit(dists, year_by_leaf, min_tips=DEFAULT_MIN_TIPS):
    """在**固定**根位下做根到尾 OLS。返回 dict 或 None（不满足门槛）。

    ``dists`` = {叶名: 根到该叶距离}；``year_by_leaf`` = {叶名: 4 位年份}。
    返回 ``{'r2','slope','intercept','n','points'}``，其中 ``points`` 是
    ``[(叶名, 距离, 年份), …]`` 原始点（供散点图/tsv）。

    与 ``phylogeo.root_to_tip`` 的回归公式**逐字相同**（x=年份、y=距离，
    slope=替换/位点/年），差别只是这里**不搜索根位** —— 根位由调用方给定。
    """
    pts = []
    for nm, d in dists.items():
        yr = year_by_leaf.get(nm)
        if yr is None:
            continue
        pts.append((nm, d, float(yr)))
    n = len(pts)
    if n < min_tips:
        return None
    mx = sum(p[2] for p in pts) / n
    my = sum(p[1] for p in pts) / n
    sxx = sum((p[2] - mx) ** 2 for p in pts)
    sxy = sum((p[2] - mx) * (p[1] - my) for p in pts)
    syy = sum((p[1] - my) ** 2 for p in pts)
    if sxx == 0 or syy == 0:
        return None          # 采样年份全同 or 所有距离全同 → R² 无定义
    slope = sxy / sxx
    r2 = (sxy ** 2) / (sxx * syy)
    return {'r2': float(r2), 'slope': float(slope),
            'intercept': float(my - slope * mx), 'n': n, 'points': pts}


def _best_over_branches(geo, year_by_leaf, min_tips):
    """逐候选根位回归，返回 (best_fit 或 None, 全部有效 fit 的 dict)。"""
    all_fit = {}
    for key, dists in geo['branches']:
        f = rtt_fit(dists, year_by_leaf, min_tips)
        if f is not None:
            all_fit[key] = f
    if not all_fit:
        return None, {}
    bk = max(all_fit, key=lambda k: all_fit[k]['r2'])
    return all_fit[bk], all_fit


def rtt_root_scan(tree_path, dates, date_trait='year',
                  min_tips=DEFAULT_MIN_TIPS):
    """把**所有**候选根位的 R² 摊开（诊断用，就是 AUDIT 缺陷的量化入口）。

    返回 ``{'n_branches','n_valid','best':{…},'asis':{…}|None,
    'median','q05','q95','hist'}``。

    用途：对比"取最大"（``best``）与"不搜索"（``asis``）的差距。差得越多，
    说明固定阈值 0.5 越不可信 —— 因为那个最大值很大程度上是**搜索出来的**。
    """
    import numpy as np

    geo = tree_geometry(tree_path)
    ybl = resolve_years(geo['names'], dates, date_trait)
    best, all_fit = _best_over_branches(geo, ybl, min_tips)
    if best is None:
        raise RuntimeError('没有任何候选根位能算出回归（带日期的叶不足 '
                           f'{min_tips} 条，或采样年份无变化）')
    r2s = np.array([f['r2'] for f in all_fit.values()], dtype=float)
    asis = rtt_fit(geo['root'], ybl, min_tips)
    hist = {'<0.1': 0, '0.1-0.25': 0, '0.25-0.5': 0, '0.5-0.75': 0, '>=0.75': 0}
    for v in r2s:
        if v < 0.1:
            hist['<0.1'] += 1
        elif v < 0.25:
            hist['0.1-0.25'] += 1
        elif v < 0.5:
            hist['0.25-0.5'] += 1
        elif v < 0.75:
            hist['0.5-0.75'] += 1
        else:
            hist['>=0.75'] += 1
    return {
        'n_branches': geo['n_branches'],
        'n_valid': int(len(r2s)),
        'n_dated_tips': len(ybl),
        'best': {'r2': round(best['r2'], 4), 'slope': best['slope'],
                 'n': best['n']},
        'asis': ({'r2': round(asis['r2'], 4), 'slope': asis['slope'],
                  'n': asis['n']} if asis else None),
        'median': round(float(np.median(r2s)), 4),
        'q05': round(float(np.percentile(r2s, 5)), 4),
        'q95': round(float(np.percentile(r2s, 95)), 4),
        'hist': hist,
    }


# ---------------------------------------------------------------- 统计

def _percentileofscore_rank(values, score):
    """复刻 ``scipy.stats.percentileofscore(values, score)``（默认 kind='rank'）。

    公式 ＝ ``(left + right + plus1) × 50 / n``，其中 ``left`` = 小于 score 的个数、
    ``right`` = 小于等于 score 的个数、``plus1`` = score 在样本里时取 1 否则 0。

    ⚠️ 旧实现写的是 ``100 × (n_< + 0.5 × n_=) / n`` —— 它等价于「并列值全取
    下秩」，而 scipy 的 rank 口径取**并列值的平均秩**，两者在有并列时差
    ``1/(2n)``：实测 ``values=[1,2,2,3,4], score=2`` scipy 给 50.0、旧实现 40.0。
    判据 ``passed_percentile95`` 直接依赖这个值，差 1/(2n) 会让边界轮次改判，
    所以这里必须与 scipy 逐位一致。自己写是为了**不引入 scipy** —— 本模块只依赖
    numpy，打包时少一个隐性依赖（平台的 `3rd/python` 就吃过"可选库静默降级"的亏）。
    """
    n = len(values)
    if n == 0:
        return None
    n_less = sum(1 for v in values if v < score)
    n_le = sum(1 for v in values if v <= score)
    return (n_le + n_less + (1 if n_le > n_less else 0)) * 50.0 / n


def _ci95(values):
    """中央 95% 区间 ``[2.5%, 97.5%]``（与 numpy.percentile 一致的线性插值）。"""
    import numpy as np
    return (float(np.percentile(values, 2.5)),
            float(np.percentile(values, 97.5)))


def date_randomization_test(tree_path, dates, date_trait='year',
                            n_perm=DEFAULT_N_PERM, seed=DEFAULT_SEED,
                            rooting='optimize', min_tips=DEFAULT_MIN_TIPS,
                            progress=None):
    """日期随机化检验（Duchêne et al. 2015）。

    参数
    ----
    tree_path : 无根/有根 Newick（平台 NJ / FastTree / RAxML-NG 产物都行）
    dates     : {叶名或原头首词: 日期}（口径同 ``phylogeo.resolve_years``）
    n_perm    : 置换次数（参照管道默认 20；建议 ≥20）
    seed      : 置换 RNG 种子（**固定**才能复现；参照管道用 42）
    rooting   : ``'optimize'`` / ``'asis'`` / ``'optimize_once'``（见模块 docstring）

    返回 dict（可 JSON 序列化）
    -------------------------
    ``real_r2`` / ``real_slope``      真实数据统计量
    ``perm_r2`` / ``perm_slope``      零分布（长度 = 有效次数）
    ``n_valid`` / ``n_invalid``       有效/失败轮次（**失败不再被静默吞掉**）
    ``r2_percentile``                 真实 R² 在零分布中的百分位（rank 口径）
    ``p_value`` / ``n_ge``            经验 p = ``(1 + #{置换 ≥ 真实}) / (1 + 有效数)``
    ``passed``                        **经验 p ≤ 0.05** ⇒ 有时间信号（默认判据）
    ``passed_percentile95``           参照管道的口径（百分位 ≥95），见下
    ``rate_overlap`` / ``rate_ci95``  真实速率是否落在零分布中央 95%（≤0 或 ≥0 同号
                                      也单独报，见 ``sign_consistent``）
    ``sign_consistent``               真实与**全部**有效轮次速率同号 ⇒ 方向一致
    ``r2_best_root``                  真实数据的"最大根位 R²"（旧判据），用于对比
    ``old_threshold_pass``            旧判据（``r2_best_root ≥ 0.5``）的结论
    ``conclusion`` / ``warnings``     人读结论与告警

    ⚠️ 为什么 ``rate_overlap`` 单独不够：零分布的速率可能全是**负**的，而真实速率
    是正的 —— 此时"落在中央 95% 区间内"为 False 也说明不了问题（区间本身就跨了
    正负）。所以另报 ``sign_consistent``：**真实速率方向与零分布必须不一致**，
    才算有信息。两个条件一起看。

    ⚠️ 为什么默认判据是经验 p 而不是参照管道的"百分位 ≥95"：两者在 ``n_perm=20``
    下**不等价**，而且后者的实际 α 随 ``n_perm`` 漂。
    ``percentileofscore`` 的 rank 口径为 ``100×(n_< + 0.5×n_=)/n``，于是"≥95"
    在 ``n_perm=20`` 时只要求真实值进**前 2 名**，在 ``n_perm=100`` 时进**前 6 名**
    → 零假设下的实际 α ≈ 2/21 = 9.5%（20 次）／6/101 = 5.9%（100 次）。
    经验 p 用 ``(1 + #{置换 ≥ 真实})/(1 + 有效数)`` 则**恒定**：
    "p ≤ 0.05" 就等于"真实值超过全部置换"，α = 1/(n_perm+1)，与 n_perm 无关。

    实测（真实 H3N2 拓扑+日期，随机取子集后**打乱日期**＝零假设成立；
    8/12 叶各 60 次、20 叶 200 次，seed=42）：

        叶数       旧判据 maxR²≥0.5   管道口径 百分位≥95   本模块 经验 p≤0.05
          8             23.3%              11.7%               6.7%
         12              3.3%               8.3%               3.3%
         20              1.5%             11.5%               6.0%
        20(n_perm=100)   1.5%              6.5%               5.5%

    同批数据的**灵敏度**（用真实日期，判"显著"的比例）：8/12/20 叶下
    旧判据与 DRT 都 ≈ 95–100% —— 即 DRT 收紧的是假阳性，没有牺牲检出力。

    两个口径都保留（``passed`` = 经验 p；``passed_percentile95`` = 管道口径），
    分歧时写 ``warnings`` 提示。**要更细的 p 就加大 ``n_perm``**：
    ``n_perm=20`` 时 p 的分辨率下限是 ``1/21 = 0.0476``，恰好卡在 0.05 边缘，
    想要"p<0.01"必须 ``n_perm ≥ 100``。
    """
    import numpy as np

    if rooting not in ROOTINGS:
        raise RuntimeError(f'rooting 只支持 {ROOTINGS}，收到 {rooting!r}')
    n_perm = max(1, int(n_perm))
    min_tips = max(3, int(min_tips))

    geo = tree_geometry(tree_path)
    ybl = resolve_years(geo['names'], dates, date_trait)
    if len(ybl) < min_tips:
        raise RuntimeError(
            f'带日期的叶只有 {len(ybl)} 条（需 ≥{min_tips}）—— 检查日期列名/'
            '日期格式，或改用更小的 min_tips（但 <5 的回归不可靠）')

    warnings = []
    real_best, real_all = _best_over_branches(geo, ybl, min_tips)
    if real_best is None:
        raise RuntimeError('没有任何候选根位能算出回归：采样年份无变化？')
    real_asis = rtt_fit(geo['root'], ybl, min_tips)
    r2_best_root = round(real_best['r2'], 4)

    # ── 冻结根位（仅 optimize_once 用）──
    frozen_key = max(real_all, key=lambda k: real_all[k]['r2'])
    frozen_dists = dict(geo['branches'])[frozen_key]

    def _stat(ybl_i):
        """按 rooting 口径算一轮统计量 → fit dict 或 None。"""
        if rooting == 'asis':
            return rtt_fit(geo['root'], ybl_i, min_tips)
        if rooting == 'optimize_once':
            return rtt_fit(frozen_dists, ybl_i, min_tips)
        return _best_over_branches(geo, ybl_i, min_tips)[0]

    real_fit = _stat(ybl)
    if real_fit is None:
        raise RuntimeError(
            f'rooting={rooting!r} 下真实数据算不出回归 —— '
            + ('树自带根位的带日期叶不足（改用 rooting="optimize"）'
               if rooting == 'asis' else '年份无变化或带日期叶不足'))

    # ── 置换：与参照管道同一个整数置换流 ──
    tip_names = list(ybl)
    vals = [ybl[n] for n in tip_names]
    perm_rng = np.random.default_rng(seed)
    perm_r2, perm_slope = [], []
    n_invalid = 0
    for i in range(n_perm):
        p = perm_rng.permutation(len(vals))
        yi = {tip_names[k]: vals[p[k]] for k in range(len(vals))}
        f = _stat(yi)
        if f is None:
            n_invalid += 1
            continue
        perm_r2.append(round(f['r2'], 6))
        perm_slope.append(f['slope'])
        if progress:
            progress('drt', 0.5 + 0.45 * (i + 1) / n_perm,
                     f'日期随机化 {i + 1}/{n_perm}')

    n_valid = len(perm_r2)
    out = {
        'engine': 'local-drt',
        'tree': os.path.basename(tree_path),
        'rooting': rooting,
        'seed': int(seed),
        'n_perm': n_perm,
        'n_valid': n_valid,
        'n_invalid': n_invalid,
        'n_dated_tips': len(ybl),
        'n_candidate_roots': geo['n_branches'],
        'real_r2': round(real_fit['r2'], 4),
        'real_slope': real_fit['slope'],
        'real_n': real_fit['n'],
        'real_r2_asis_root': round(real_asis['r2'], 4) if real_asis else None,
        'r2_best_root': r2_best_root,
        'perm_r2': perm_r2,
        'perm_slope': perm_slope,
        'r2_percentile': None,
        'p_value': None,
        'n_ge': None,
        'passed': None,
        'passed_percentile95': None,
        'rate_overlap': None,
        'rate_ci95': None,
        'sign_consistent': None,
        'r2_ci95': None,
        'old_threshold_pass': bool(r2_best_root >= 0.5),
        # 单一结论字段（动力学轮暂缓项 4）—— 下游别再加自己的组合逻辑
        'verdict': 'inconclusive',
        'verdict_basis': None,
        'verdict_dissent': [],
        'verdict_text': '',
        'alpha': DRT_ALPHA,
        'min_valid_for_alpha': None,
        'conclusion': '',
        'warnings': warnings,
    }

    if n_invalid:
        warnings.append(
            f'{n_invalid}/{n_perm} 次随机化算不出回归（带日期叶不足或年份无变化），'
            f'零分布实际只有 {n_valid} 个点 —— 名义次数不等于有效次数')
    # 经验 p 的最小可达值 = 1/(1+n_valid)：达不到 α 就说明这一轮**没有功效**，
    # 此时报 "未通过" 是把没有检出力的检验当阴性证据（旧行为）。结论必须给
    # 'inconclusive'，并说明要加多少次置换。
    min_valid = int(math.ceil(1.0 / DRT_ALPHA)) - 1
    out['min_valid_for_alpha'] = min_valid
    if n_valid == 0:
        out['verdict'] = 'inconclusive'
        out['verdict_basis'] = '零分布为空（所有随机化都失败）→ 无法判读'
        out['verdict_text'] = '⚠ 无结论：日期随机化全部失败，无法判断时间信号'
        out['conclusion'] = 'DRT 无结论：零分布为空（所有随机化都失败）'
        return out

    pct = _percentileofscore_rank(perm_r2, real_fit['r2'])
    out['r2_percentile'] = round(float(pct), 2)
    # 经验 p（North et al. 2002 的 +1 修正）：零分布里 ≥ 真实的轮数越多越不显著。
    # 用 +1 而不是裸比例，是因为真实数据本身就是零假设下的一次抽样，必须计入分母
    # （裸比例会在 n_ge=0 时给出 p=0，"不可能更极端"是假象）。
    n_ge = sum(1 for v in perm_r2 if v >= real_fit['r2'])
    p_value = (1.0 + n_ge) / (1.0 + n_valid)
    out['n_ge'] = int(n_ge)
    out['p_value'] = round(p_value, 4)
    out['passed'] = bool(p_value <= 0.05)
    out['passed_percentile95'] = bool(pct >= 95.0)
    out['r2_ci95'] = [round(v, 4) for v in _ci95(perm_r2)]

    # conclusion 由**结论分支**统一赋值（原先在这里按 passed 写死，功效不足时
    # 会出现「verdict=无结论」+「conclusion=✗ 未通过」自相矛盾的一屏两份结论）。

    if out['passed'] != out['passed_percentile95']:
        warnings.append(
            f'两种口径分歧：经验 p 判 {"通过" if out["passed"] else "不通过"}，'
            f'参照管道的"百分位≥95"判 '
            f'{"通过" if out["passed_percentile95"] else "不通过"}'
            + ('（本模块以经验 p 为准：n_perm 小的时候"百分位≥95"的名额被放宽 —— '
               'n_perm=20 进前 2 名、n_perm=100 进前 6 名都算过，实际 α 可到 9.5%）'
               if out['passed_percentile95'] else ''))

    # ── 速率层：区间重叠 + 方向一致性 ──
    if len(perm_slope) >= 5:
        lo, hi = _ci95(perm_slope)
        rs = float(real_fit['slope'])
        out['rate_ci95'] = [lo, hi]
        out['rate_overlap'] = bool(lo <= rs <= hi)
        out['sign_consistent'] = bool(
            (rs > 0 and all(v > 0 for v in perm_slope))
            or (rs < 0 and all(v < 0 for v in perm_slope)))
        if out['sign_consistent']:
            warnings.append(
                f'真实速率方向（{"正" if rs > 0 else "负"}）与零分布**全部**轮次一致 '
                '—— 方向本身不含时间信息（时间结构被随机日期复现了）')
    else:
        warnings.append(f'有效速率只有 {len(perm_slope)} 个（需 ≥5）→ '
                        'rate_overlap / sign_consistent 不计算')

    # ── 反保守口径的显式告警 ──
    if rooting == 'optimize_once':
        warnings.append(
            'rooting="optimize_once" 是**反保守**口径：真实值的根位是在真实数据上'
            '挑的最优，而零分布用的是同一个冻结根位（没有再做搜索）→ 真实侧多了一次'
            '"取最大"，通过率被抬高。仅用于诊断/复现，判定请用 optimize 或 asis。')
    if rooting == 'optimize':
        warnings.append(
            'rooting="optimize"：零分布同样逐轮取最大根位 R²（与参照发布版管道 '
            f'temporal_signal.py 一致）。真实数据"最大根位 R²" = {r2_best_root}，'
            '而树自带根位只有 '
            + (f'{out["real_r2_asis_root"]}' if out['real_r2_asis_root'] is not None
               else '算不出')
            + ' —— 两者差距就是"根位搜索"贡献的量，**单看 R² 数值会高估信号**。')

    # ── 单一结论字段（2026-09-16，动力学轮暂缓项 4）──────────────────
    # 此前三个判据（经验 p / 管道"百分位≥95" / 旧固定阈值 R²≥0.5）各自成字段，
    # 前端、summary.csv、报表、脚本各写一套组合逻辑 → 同一份结果在不同入口
    # 可能显示不同结论。这里给出**权威**字段，下游只读它：
    #   verdict ∈ {'pass','fail','inconclusive'}
    # · 'pass' 只认经验 p（α 恒定，见模块 docstring）；
    # · **功效不足时不给 'fail'** —— n_valid 太小时最小可达 p = 1/(1+n_valid)
    #   已经 > α，该检验不可能通过，此时"未通过"是没有检出力的检验被当成
    #   阴性证据；改成 'inconclusive' 并给出需要的置换次数。
    if out['passed']:
        out['verdict'] = 'pass'
        out['verdict_basis'] = f'经验 p={p_value:.4f} ≤ α={DRT_ALPHA}'
        out['verdict_text'] = (
            f'✓ 有时间信号（DRT 通过）：经验 p={p_value:.3f}，'
            f'真实 R²={real_fit["r2"]:.4f}，可继续做分子钟定年')
        out['conclusion'] = (
            f'✓ 通过：真实 R²={real_fit["r2"]:.4f}，经验 p={p_value:.3f}'
            f'（{n_ge}/{n_valid} 次日期随机化 ≥ 真实值，阈值 {DRT_ALPHA}）'
            f'—— 时间信号显著，适合做分子钟定年')
    elif n_valid < min_valid:
        out['verdict'] = 'inconclusive'
        out['verdict_basis'] = (
            f'有效置换只有 {n_valid} 次 < {min_valid} —— 最小可达 p='
            f'{1.0 / (1.0 + n_valid):.3f} 已 > α={DRT_ALPHA}，本次检验无功效')
        out['verdict_text'] = (
            f'⚠ 无结论（检验无功效）：有效置换 {n_valid} 次，最小可达 '
            f'p={1.0 / (1.0 + n_valid):.3f} > α={DRT_ALPHA}，'
            f'无论真实 R² 多高都不可能判通过 —— 请把置换次数提到 ≥'
            f'{min_valid + 1}（且确保每轮都能算出回归）后重跑；'
            f'本次真实 R²={real_fit["r2"]:.4f} 仅供参考')
        # conclusion 必须与 verdict 同口径：写「✗ 未通过」会把"没检出力"
        # 当成阴性证据（且与页面上的「⚠ 无结论」自相矛盾）。
        out['conclusion'] = (
            f'⚠ 无结论：真实 R²={real_fit["r2"]:.4f}，经验 p={p_value:.3f}'
            f'（{n_ge}/{n_valid} 次日期随机化 ≥ 真实值）—— 有效置换只有 '
            f'{n_valid} 次，最小可达 p={1.0 / (1.0 + n_valid):.3f} 已 > '
            f'{DRT_ALPHA}，本次检验无功效：既不能判"有信号"也不能判"无信号"，'
            f'请把置换次数提到 ≥{min_valid + 1} 后重跑')
        warnings.append(
            f'零分布只有 {n_valid} 次有效置换 → DRT 无功效（最小可达 p='
            f'{1.0 / (1.0 + n_valid):.3f} > {DRT_ALPHA}），'
            f'结论为「无结论」而非「无信号」：把 n_perm 提到 ≥{min_valid + 1}')
    else:
        out['verdict'] = 'fail'
        out['verdict_basis'] = f'经验 p={p_value:.4f} > α={DRT_ALPHA}（有效置换 {n_valid} 次）'
        out['verdict_text'] = (
            f'✗ 无时间信号（DRT 未通过）：经验 p={p_value:.3f}，'
            f'真实 R²={real_fit["r2"]:.4f} —— 分子钟定年不可靠')
        out['conclusion'] = (
            f'✗ 未通过：真实 R²={real_fit["r2"]:.4f}，经验 p={p_value:.3f}'
            f'（{n_ge}/{n_valid} 次日期随机化 ≥ 真实值，阈值 {DRT_ALPHA}，'
            f'百分位 {pct:.0f}%）—— 时间信号不足，分子钟定年不可靠')

    # ── 与旧判据的分歧提示（这条最有信息量；措辞随 verdict 走，否则
    #    「无结论」的那一轮会被写成"DRT 判**不通过**"）──
    if out['old_threshold_pass'] and not out['passed']:
        if out['verdict'] == 'inconclusive':
            warnings.append(
                f'⚠️ 新旧判据分歧：旧的固定阈值法（最大根位 R²≥0.5，此处 '
                f'{r2_best_root}）会判"有时间信号"，但 DRT 本次**无结论**'
                f'（有效置换 {n_valid} 次 < {min_valid}，检验无功效）—— '
                '既不能说"有信号"也不能说"无信号"，请把置换次数提到 ≥'
                f'{min_valid + 1} 后重跑再引用本条分歧。')
        else:
            warnings.append(
                f'⚠️ 新旧判据分歧：旧的固定阈值法（最大根位 R²≥0.5，此处 '
                f'{r2_best_root}）会判"有时间信号"，但 DRT 判**不通过**'
                f'（经验 p={p_value:.3f}，百分位 {pct:.0f}%）—— '
                '这正是 AUDIT 记的"取最大 + 固定阈值"在小样本上的高估，'
                '**以 DRT 为准**。实测零假设下：8 叶时旧判据误判 23.3%、DRT 5.0%。')

    # 分歧清单：其它判据是否与权威结论不一致（供前端一行说明，不必再解释规则）
    dissent = []
    if out['passed_percentile95'] != out['passed']:
        dissent.append('参照管道「百分位≥95」= '
                       + ('通过' if out['passed_percentile95'] else '不通过'))
    if out['old_threshold_pass'] != out['passed']:
        dissent.append('旧固定阈值 R²≥0.5 = '
                       + ('通过' if out['old_threshold_pass'] else '不通过'))
    if out['sign_consistent']:
        dissent.append('真实速率与全部置换同号（方向无信息）')
    out['verdict_dissent'] = dissent
    out['warnings'] = warnings
    return out


# ══════════════════════════════════════════════════════════════════════════
# DRT · TreeTime 统计量版（2026-09-18，`#t-rtt` 换 TreeTime 引擎时加）
# ══════════════════════════════════════════════════════════════════════════
#
# 为什么要有第二个 DRT：
#   `#t-rtt` 的根到尾回归引擎从「平台自研最小二乘」换成「TreeTime 定年」之后，
#   上面那个 `date_randomization_test()` 的统计量（`rtt_fit` / `rtt_root_scan`
#   的 R²）就与面板上展示的 R² **不是同一个量**了 —— 真实值与零分布必须同口径，
#   否则 p 值无意义。所以零分布也得用 TreeTime 的 R² 来造。
#
# ⚠️ 代价比老版本高一个数量级：**每一轮置换都要完整重跑一次 TreeTime**。
#   476 叶实测 87 s/次（本机，`--verbose 4`、`--time-marginal true`），
#   n_perm=10 → 约 16 分钟。调用方**必须**先把预估耗时写进日志告诉用户。

# ⚠️ 默认 **20**，不是 10。理由：DRT 的判据是经验 p ≤ α=0.05，而 n 次置换能
# 达到的最小 p 是 1/(n+1)。n=10 → 最小 p = 0.0909 **永远 > 0.05** ⇒ verdict
# 恒为「无结论（检验无功效）」，默认跑出来的判据等于没用（实测 80 叶、n=10：
# 真实 R²=0.9897、零分布最高才 0.0502，分离度极大，却只能判「无结论」）。
# n=20 → 最小 p = 1/21 = 0.0476 ≤ 0.05，才真的能判「通过 / 不通过」。
# 代价：每轮一次完整 TreeTime（476 叶实测 87 s/次 → 20 次约 30 分钟）。
DEFAULT_N_PERM_TT = 20


def estimate_treetime_perm_cost(sec_per_run, n_perm):
    """预估 DRT-TreeTime 的总耗时（秒）：真实值 1 次 + 置换 n_perm 次。"""
    try:
        return float(sec_per_run) * (1 + max(0, int(n_perm)))
    except (TypeError, ValueError):
        return None


def _stat_treetime(res):
    """从一次 TreeTime-RTT 结果里取统计量 → `(值, 口径说明)`。

    优先用 `r2_precise`（从 TreeTime 自己写的 `divergence_tree.nexus` 把**同一个
    回归**按全精度重算）；退化到 `r2`（`molecular_clock.txt` 的原值）。

    ⚠️ **必须优先用全精度那个**：TreeTime 把 `--r^2` 写死成 `%1.2f`
    （`treetime/utils.py:32`），零分布会塌成 {0.98, 0.99, 1.00} 几个格子，
    秩检验根本没有分辨率。
    """
    if not res or not res.get('ok'):
        return None, None
    if res.get('r2_precise') is not None:
        return float(res['r2_precise']), 'r2_precise（TreeTime 回归的全精度重算）'
    if res.get('r2') is not None:
        return float(res['r2']), 'r2（molecular_clock.txt，2 位小数）'
    return None, None


def date_randomization_test_treetime(aln_path, tree_path, dates, work_dir, *,
                                     date_trait='year', n_perm=DEFAULT_N_PERM_TT,
                                     seed=DEFAULT_SEED, seq_len=None,
                                     timeout=1800, rtt_kw=None, say=None,
                                     progress=None, cancel=None,
                                     sec_per_run_hint=None,
                                     drt_permutation=False):
    """日期随机化检验（DRT）—— 统计量 = **TreeTime 的根到尾回归 R²**。

    与同模块 `date_randomization_test()` 的唯一区别是统计量的来源：那个用平台
    自研的最小二乘回归，这个用 TreeTime 跑完整条定年流程后对分歧树的回归。
    RTT 引擎换成 TreeTime 之后，真实值与零分布必须同口径。

    零假设：采样日期与拓扑无关。把日期在 tip 之间随机置换 n_perm 次，每次
    **重跑一次 TreeTime**，得到零分布；真实 R² 在零分布中的位置给出经验 p。

    ⚠️ 代价：每轮一次完整 TreeTime。476 叶实测 87 s/次，n_perm=20 ≈ 30 分钟。
    **因此置换默认不执行**：`drt_permutation=False`（默认）时只做真实数据这一轮
    TreeTime，结果里 `drt_permutation=False`、`p_value=None`、verdict 恒为
    'inconclusive'；要真正的零分布检验必须显式传 `drt_permutation=True`
    （调用方/`run_clock` 侧透传）。结果 dict 的 `drt_permutation` 字段记录
    本次是否真的执行了置换，**引用 p 值前必须先看它**。`cancel` 在每轮之间检查。

    `dates`: `{叶名: 十进制年}`。`work_dir`: 每轮一个子目录（`real/`、`perm_000/`…），
    TreeTime 的全部产物都留着（可逐轮复核，不删中间态）。

    返回 dict，字段与 `date_randomization_test()` **同名同义**（`verdict` /
    `p_value` / `perm_r2` / `conclusion` …），另加 `engine='treetime-drt'`、
    `stat_basis`、`sec_per_run`、`elapsed`、`drt_permutation`。
    """
    import time as _time
    # 本模块的既有约定：numpy **只在函数内局部导入**（见 213 / 269 / 334 行），
    # 保持模块顶层只依赖 stdlib + phylogeo。
    import numpy as np
    from Virus_Platform_Core import treetime_rtt as trr

    say = say or (lambda *_a, **_k: None)
    rtt_kw = dict(rtt_kw or {})
    warnings = []
    work = os.path.join(work_dir, 'drt_tt')
    os.makedirs(work, exist_ok=True)

    names = [n for n in dates if n is not None]
    if len(names) < DEFAULT_MIN_TIPS:
        raise RuntimeError(f'带日期的叶只有 {len(names)} 个（< {DEFAULT_MIN_TIPS}），'
                           '跑不了 DRT')
    vals = [dates[n] for n in names]
    n_perm = max(1, int(n_perm))

    est = estimate_treetime_perm_cost(sec_per_run_hint, n_perm)
    say(f'[DRT-TreeTime] {len(names)} 叶，置换 {n_perm} 次 —— '
        '**每轮完整重跑一次 TreeTime**'
        + (f'，按 {sec_per_run_hint:.0f}s/次预估总耗时约 {est / 60.0:.1f} 分钟'
           if sec_per_run_hint and est else ''))

    def _one(d, tag):
        sub = os.path.join(work, tag)
        return trr.run_treetime_rtt(aln_path, tree_path, d, sub,
                                    date_trait=date_trait, seq_len=seq_len,
                                    timeout=timeout, say=None, cancel=cancel,
                                    **rtt_kw)

    def _cancelled():
        return bool(cancel is not None and cancel.is_set())

    t0 = _time.time()
    # ── 真实数据 ──
    if progress:
        progress('drt', 0.45, 'DRT-TreeTime：真实数据')
    real = _one(dict(zip(names, vals)), 'real')
    real_stat, stat_basis = _stat_treetime(real)
    if real_stat is None:
        raise RuntimeError('真实数据这一轮没跑出 R²（%s）—— DRT 无法进行'
                           % (real.get('error') or '未写出 molecular_clock.txt'))
    sec_per_run = _time.time() - t0
    say(f'[DRT-TreeTime] 真实 R²={real_stat:.6f}，统计量口径：{stat_basis}；'
        f'本轮 {sec_per_run:.1f}s；预计总耗时 '
        f'{estimate_treetime_perm_cost(sec_per_run, n_perm) / 60.0:.1f} 分钟')

    # ── 置换（可选：默认关闭，见 docstring 的代价说明）──
    loop_n = n_perm if drt_permutation else 0
    if not drt_permutation:
        say('[DRT-TreeTime] drt_permutation=False（默认）→ **跳过置换**，'
            '只保留真实数据这一轮的 TreeTime 结果；要零分布检验请显式打开。')
    rng = np.random.default_rng(seed)
    perm_r2, perm_rate, perm_ok, perm_fail = [], [], [], []
    n_invalid = 0
    for i in range(loop_n):
        if _cancelled():
            warnings.append(f'用户取消：置换做到 {i}/{n_perm} 轮就停了，'
                            f'零分布只有 {len(perm_r2)} 个点')
            break
        p = rng.permutation(len(vals))
        perm_dates = {names[k]: vals[p[k]] for k in range(len(vals))}
        r = _one(perm_dates, 'perm_%03d' % i)
        v, _b = _stat_treetime(r)
        if v is None:
            n_invalid += 1
            perm_fail.append({'i': i, 'error': r.get('error') or '无 R²'})
        else:
            perm_r2.append(round(v, 6))
            if r.get('rate') is not None:
                perm_rate.append(float(r['rate']))
            perm_ok.append(i)
        if progress:
            progress('drt', 0.45 + 0.5 * (i + 1) / n_perm,
                     f'DRT-TreeTime 置换 {i + 1}/{n_perm}'
                     + (f'（本轮失败，累计 {n_invalid} 次）' if v is None else ''))
        say(f'[DRT-TreeTime] 置换 {i + 1}/{n_perm}：'
            + (f'R²={v:.6f}' if v is not None
               else f'失败（{r.get("error") or "无 R²"}）')
            + f'（累计 {_time.time() - t0:.0f}s）')

    n_valid = len(perm_r2)
    out = {
        'engine': 'treetime-drt',
        'stat': 'r2',
        'stat_basis': stat_basis,
        'tree': os.path.basename(tree_path),
        'aln': os.path.basename(aln_path),
        'seed': int(seed),
        'n_perm': n_perm,
        'n_valid': n_valid,
        'n_invalid': n_invalid,
        'n_dated_tips': len(names),
        # 本次是否真的执行了置换零分布（False ＝ 只跑了真实数据一轮，
        # p_value/perm_r2 等字段不可用于判读）。
        'drt_permutation': bool(drt_permutation),
        'sec_per_run': round(sec_per_run, 2),
        'elapsed': round(_time.time() - t0, 1),
        'real_r2': round(real_stat, 6),
        'real_r2_treetime': real.get('r2'),
        'real_rate': real.get('rate'),
        'real_rate_precise': real.get('rate_precise'),
        'real_n': real.get('n_points'),
        'perm_r2': perm_r2,
        'perm_rate': perm_rate,
        'perm_ok': perm_ok,
        'perm_fail': perm_fail[:20],
        'r2_percentile': None,
        'p_value': None,
        'n_ge': None,
        'passed': None,
        'passed_percentile95': None,
        'rate_overlap': None,
        'rate_ci95': None,
        'sign_consistent': None,
        'r2_ci95': None,
        'old_threshold_pass': bool(real_stat >= 0.5),
        'verdict': 'inconclusive',
        'verdict_basis': None,
        'verdict_dissent': [],
        'verdict_text': '',
        'alpha': DRT_ALPHA,
        'min_valid_for_alpha': None,
        'conclusion': '',
        'work_dir': work,
        'warnings': warnings,
    }

    if n_invalid:
        warnings.append(
            f'{n_invalid}/{n_perm} 轮置换 TreeTime 没跑出 R²，'
            f'零分布实际只有 {n_valid} 个点 —— 名义次数不等于有效次数'
            + (f'；失败原因示例：{perm_fail[0]["error"]}' if perm_fail else ''))
    min_valid = int(math.ceil(1.0 / DRT_ALPHA)) - 1
    out['min_valid_for_alpha'] = min_valid
    if n_valid == 0:
        if not drt_permutation:
            out['verdict_basis'] = ('置换未执行（drt_permutation=False，默认）'
                                    '→ 无零分布，无法判读')
            out['verdict_text'] = (
                '⚠ 无结论：未执行日期随机化置换（drt_permutation=False）'
                f'—— 以上真实 R²={real_stat:.4f} 只是 TreeTime 的时钟拟合，'
                '不是 DRT 判据；需要检验请显式传 drt_permutation=True')
            out['conclusion'] = ('DRT（TreeTime 口径）无结论：置换未执行'
                                 '（drt_permutation=False）')
            warnings.append(
                'drt_permutation=False：只跑了真实数据一轮 TreeTime，'
                '**没有**零分布 → p 值/显著性判据本次不可用。')
            return out
        out['verdict_basis'] = '零分布为空（所有置换都没跑出 R²）→ 无法判读'
        out['verdict_text'] = '⚠ 无结论：日期随机化全部失败，无法判断时间信号'
        out['conclusion'] = 'DRT（TreeTime 口径）无结论：零分布为空'
        return out

    pct = _percentileofscore_rank(perm_r2, real_stat)
    out['r2_percentile'] = round(float(pct), 2)
    n_ge = sum(1 for v in perm_r2 if v >= real_stat)
    p_value = (1.0 + n_ge) / (1.0 + n_valid)
    out['n_ge'] = int(n_ge)
    out['p_value'] = round(p_value, 4)
    out['passed'] = bool(p_value <= 0.05)
    out['passed_percentile95'] = bool(pct >= 95.0)
    out['r2_ci95'] = [round(v, 4) for v in _ci95(perm_r2)]

    if out['passed'] != out['passed_percentile95']:
        warnings.append(
            f'两种口径分歧：经验 p 判 {"通过" if out["passed"] else "不通过"}，'
            f'"百分位≥95"判 {"通过" if out["passed_percentile95"] else "不通过"}'
            '（本模块以经验 p 为准）')

    if len(perm_rate) >= 5:
        lo, hi = _ci95(perm_rate)
        rs = float(real.get('rate') or 0.0)
        out['rate_ci95'] = [lo, hi]
        out['rate_overlap'] = bool(lo <= rs <= hi)
        out['sign_consistent'] = bool(
            (rs > 0 and all(v > 0 for v in perm_rate))
            or (rs < 0 and all(v < 0 for v in perm_rate)))
        if out['sign_consistent']:
            warnings.append(
                '真实速率方向与零分布**全部**轮次一致 —— 方向本身不含时间信息'
                '（时间结构被随机日期复现了）')
    else:
        warnings.append(f'有效速率只有 {len(perm_rate)} 个（需 ≥5）→ '
                        'rate_overlap / sign_consistent 不计算')

    if out['passed']:
        out['verdict'] = 'pass'
        out['verdict_basis'] = f'经验 p={p_value:.4f} ≤ α={DRT_ALPHA}'
        out['verdict_text'] = (
            f'✓ 有时间信号（DRT 通过）：经验 p={p_value:.3f}，'
            f'真实 R²={real_stat:.4f}（TreeTime 口径），可继续做分子钟定年')
        out['conclusion'] = (
            f'✓ 通过：真实 R²={real_stat:.4f}，经验 p={p_value:.3f}'
            f'（{n_ge}/{n_valid} 次日期随机化 ≥ 真实值，阈值 {DRT_ALPHA}）'
            '—— 时间信号显著（统计量＝TreeTime 根到尾回归 R²）')
    elif n_valid < min_valid:
        out['verdict'] = 'inconclusive'
        out['verdict_basis'] = (
            f'有效置换只有 {n_valid} 次 < {min_valid} —— 最小可达 p='
            f'{1.0 / (1.0 + n_valid):.3f} 已 > α={DRT_ALPHA}，本次检验无功效')
        out['verdict_text'] = (
            f'⚠ 无结论（检验无功效）：有效置换 {n_valid} 次，最小可达 '
            f'p={1.0 / (1.0 + n_valid):.3f} > α={DRT_ALPHA} —— 请把置换次数'
            f'提到 ≥{min_valid + 1} 后重跑；本次真实 R²={real_stat:.4f} 仅供参考')
        out['conclusion'] = (
            f'⚠ 无结论：真实 R²={real_stat:.4f}，经验 p={p_value:.3f}'
            f'（{n_ge}/{n_valid}）—— 有效置换只有 {n_valid} 次，本次检验无功效：'
            '既不能判"有信号"也不能判"无信号"')
        warnings.append(
            f'零分布只有 {n_valid} 次有效置换 → DRT 无功效（最小可达 p='
            f'{1.0 / (1.0 + n_valid):.3f} > {DRT_ALPHA}），'
            '结论为「无结论」而非「无信号」')
    else:
        out['verdict'] = 'fail'
        out['verdict_basis'] = (f'经验 p={p_value:.4f} > α={DRT_ALPHA}'
                                f'（有效置换 {n_valid} 次）')
        out['verdict_text'] = (
            f'✗ 无时间信号（DRT 未通过）：经验 p={p_value:.3f}，'
            f'真实 R²={real_stat:.4f} —— 分子钟定年不可靠')
        out['conclusion'] = (
            f'✗ 未通过：真实 R²={real_stat:.4f}，经验 p={p_value:.3f}'
            f'（{n_ge}/{n_valid} 次日期随机化 ≥ 真实值，阈值 {DRT_ALPHA}，'
            f'百分位 {pct:.0f}%）—— 时间信号不足，分子钟定年不可靠')

    dissent = []
    if out['passed_percentile95'] != out['passed']:
        dissent.append('「百分位≥95」= '
                       + ('通过' if out['passed_percentile95'] else '不通过'))
    if out['old_threshold_pass'] != out['passed']:
        dissent.append('旧固定阈值 R²≥0.5 = '
                       + ('通过' if out['old_threshold_pass'] else '不通过'))
    if out['sign_consistent']:
        dissent.append('真实速率与全部置换同号（方向无信息）')
    out['verdict_dissent'] = dissent
    return out
