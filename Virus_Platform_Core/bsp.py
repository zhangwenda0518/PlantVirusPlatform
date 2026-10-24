# -*- coding: utf-8 -*-
r"""BSP（Bayesian Skyline Plot）—— 贝叶斯天际线图的**读取/复算/绘图**。

## 来源与口径

对应 VirPhyKit 的 `src/BSP/`（菜单名 **BSP-Viz**，Quick_Guide 原文：
"Visualizes the export results of Bayesian skyline plots"）。已确认它**是纯绘图器，
不是估计器**：`function_bsp.py::generate_single_plot` 只做

    Skyline = read.table(tsv, header=T)          # 需要 Time / Lower / Median / Upper 四列
    plot(Skyline$Median ~ Skyline$Time, log="y", xlim=c(1974,2016), …)
    polygon(Time, Lower…Upper)                    # 95% 带
    → pdf(output)

即：**输入是一张已经算好的 TSV，输出是 PDF**。所以我们分两条路补齐：

* ``mode='table'``  —— 与上游**1:1**：读 ``Time/Lower/Median/Upper`` 的 tsv/csv，
  出图；并顺手加两项上游没有的（① x 轴范围**自适应**而不是写死 1974–2016；
  ② 时间轴反向时**不靠硬编码**，按下标的 min/max 反排）。
* ``mode='log'``    —— 上游没有的增值项：直接吃 **BEAST 的 BSP 日志**
  （``skyline.popSize1..N`` + ``skyline.groupSize1..N``），自己算分位数。
  这一条必须看清下面的口径警告。

## ⚠️ 两条路给的"Lower/Upper"含义不同，别混引用

* ``table``：Lower/Upper 是**别人算好的**（通常来自 Tracer 的 BSP 重建导出）。
  本模块**原样转发**，不改数、不重算 —— 引用时须写明"经 Tracer/上游导出"。
* ``log``：本模块取 **每个分组的 ``skyline.popSize_i`` 的后验分位数**。
  这部分是**精确**的：BSP 里第 *i* 段的有效种群大小就是一个独立参数，其后验样本
  就是日志里的那一列，分位数没有歧义。

## ⚠️⚠️ ``log`` 模式的时间轴是**近似**，务必连带报告

BSP 的第 *i* 段覆盖哪个时间区间，取决于该段里**实际发生的那几次溯祖事件在树上的时刻**
—— 严格重建要读 ``.trees`` 文件、逐样本重算（Tracer 的
"Bayesian skyline reconstruction" 就是这么做的）。**只有 ``.log`` 时拿不到这些时刻**，
本模块按"溯祖事件在时间上均匀分布"折算：

    第 i 段的边界 = tMRCA × cumsum(groupSize)[i] / sum(groupSize)      （逐样本算，再取中位）

这是常见的简化口径，**不是 BSP 的定义**。因此：

* 返回值里 ``time_axis`` 字段会标 ``'approx-uniform-coalescent'``，
  并在 ``warnings`` 里明确写出"时间轴为近似，正式引用请用 Tracer 重建结果或提供 .trees"；
* 页面上必须显示这条告警（别让图看起来像严格结果）；
* 若你手上有 Tracer 导出的 TSV，走 ``table`` 模式 —— 那才是严格口径。
"""

import math
import os
import re

import numpy as np

from Virus_Platform_Core.utils import check_path, safe_open
from Virus_Platform_Core.io_utils import read_rows as _io_read_rows

# 日志里 BSP 参数列的两种写法：多段 skyline.popSize1..N / 单段 skyline.popSize
_POP_RE = re.compile(r'^skyline\.popSize(\d*)$')
_GRP_RE = re.compile(r'^skyline\.groupSize(\d*)$')
# 可能承载"树高/tMRCA"的列名（不同 BEAST 版本/模板不一）
_HEIGHT_COLS = ('treeModel.rootHeight', 'tree.height', 'TreeHeight', 'rootHeight',
                'height', 'tMRCA', 'TreeHeight.t:Species')


# ---------------------------------------------------------------- 通用小工具

def _percentile(vals, q):
    """线性插值分位数（与 numpy.percentile 默认口径一致）。"""
    a = np.asarray(vals, dtype=float)
    a = a[np.isfinite(a)]
    if a.size == 0:
        return float('nan')
    return float(np.percentile(a, q * 100.0))


def effective_sample_size(vals):
    """ESS —— **初始正序列**估计量（Geyer 1992），Tracer 用的也是这一族。

        tau   = 1 + 2·Σ_{k=1..K} ρ_k ，取到 ρ_k + ρ_{k+1} < 0 处截断（首负对）
        ESS   = n / tau                  （钳在 [2, n]）

    为什么不用"取到自相关首次为负"：那会把噪声当信号、低估 ESS。
    白噪声自检：iid 样本应给出 ESS ≈ n（见 tests/_check_bsp.py）。
    """
    a = np.asarray(vals, dtype=float)
    a = a[np.isfinite(a)]
    n = a.size
    if n < 4:
        return float('nan')
    sd = a.std(ddof=1)
    if not (sd > 0):
        return float(n)
    x = a - a.mean()
    # 自相关（用 FFT 一次算全部 lag，n 大时比逐 lag 快得多）
    nfft = 1 << int(math.ceil(math.log2(2 * n)))
    f = np.fft.rfft(x, nfft)
    acov = np.fft.irfft(f * np.conjugate(f), nfft)[:n]
    acov /= np.arange(n, 0, -1)          # 无偏化：除以 (n-k)
    if acov[0] <= 0:
        return float(n)
    rho = acov / acov[0]
    tau = 1.0
    k = 1
    while k + 1 < n:
        pair = rho[k] + rho[k + 1]
        if pair < 0:                      # 首负对 → 截断
            break
        tau += 2.0 * pair
        k += 2
    if tau <= 0:
        return float(n)
    return float(min(max(n / tau, 2.0), float(n)))


def _read_table_rows(path):
    """读 TSV/CSV → (表头 list, 行 list, 分隔符)。跳过 ``#`` 注释行。"""
    hdr, rows, delim = _io_read_rows(path, skip_comments=True)
    if hdr is None:
        raise RuntimeError(f'文件里没有数据行：{path}')
    return hdr, rows, delim


def _pick_col(hdr, wanted, what):
    """按候选名找列（大小写/空格不敏感；再退化为子串匹配）。"""
    low = [h.strip().lower().replace(' ', '') for h in hdr]
    for w in wanted:
        if w.lower().replace(' ', '') in low:
            return low.index(w.lower().replace(' ', ''))
    for w in wanted:
        token = w.lower().replace(' ', '')
        for i, h in enumerate(low):
            if token in h:
                return i
    raise RuntimeError(f'表里找不到{what}列（试过 {list(wanted)}）；实际列：{hdr}')


# ---------------------------------------------------------------- 模式 1：读现成表

def read_skyline_table(path):
    """读 ``Time/Lower/Median/Upper`` 的表（= VirPhyKit BSP-Viz 的输入）。

    列名不敏感（``time`` / ``Time`` / ``TIME`` 都认），顺序任意；多出的列忽略。
    返回 ``{'time': [...], 'lower': [...], 'median': [...], 'upper': [...], 'n': n}``
    —— 原样搬运，不做任何数值改动（口径见模块头）。
    """
    path = check_path(path, must_exist=True)
    hdr, rows, _ = _read_table_rows(path)
    i_t = _pick_col(hdr, ('Time', 'time'), 'Time')
    i_m = _pick_col(hdr, ('Median', 'median'), 'Median')
    i_l = _pick_col(hdr, ('Lower', 'lower', 'Lower95', 'lower95'), 'Lower')
    i_u = _pick_col(hdr, ('Upper', 'upper', 'Upper95', 'upper95'), 'Upper')

    def conv(idx):
        out = []
        for r in rows:
            v = r[idx].strip()
            try:
                out.append(float(v))
            except ValueError:
                out.append(float('nan'))
        return out

    t, m, lo, up = conv(i_t), conv(i_m), conv(i_l), conv(i_u)
    keep = [k for k in range(len(t)) if all(np.isfinite([t[k], m[k], lo[k], up[k]]))]
    if not keep:
        raise RuntimeError('表里 Time/Lower/Median/Upper 四列没有一行是完整数值 —— '
                           f'首行={rows[0] if rows else None}')
    dropped = len(t) - len(keep)
    return {'time': [t[k] for k in keep], 'median': [m[k] for k in keep],
            'lower': [lo[k] for k in keep], 'upper': [up[k] for k in keep],
            'n': len(keep), 'n_dropped': dropped, 'columns': hdr,
            'time_axis': 'as-given', 'source': os.path.basename(path)}


# ---------------------------------------------------------------- 模式 2：读 BEAST 日志

def skyline_from_beast_log(log_path, burnin=0.1, ci=0.95, time_scale='auto'):
    """从 BEAST BSP 日志算天际线：**每段取该段 popSize 后验的分位数**。

    ``burnin``：丢掉的**比例**（0.1 = 前 10%），与平台其它 BEAST 日志解析一致。
    ``ci``：可信区间宽度（0.95 → Lower=2.5%、Upper=97.5%）。
    ``time_scale``：
      · ``'auto'``（默认）—— 有 groupSize 列就按**均匀溯祖**折算时间轴
        （**近似**，见模块头），否则退回"段序号"当横轴并告警；
      · ``'group'`` —— 强制用段序号（不猜时间）。

    返回 ``dict``，含 ``time/lower/median/upper``（与 ``read_skyline_table`` 同构）、
    每段的 ``groups``（median/lower/upper/ess）、``time_axis``、``warnings``。
    """
    log_path = check_path(log_path, must_exist=True)
    hdr, rows, _ = _read_table_rows(log_path)
    if len(rows) < 20:
        raise RuntimeError(f'日志只有 {len(rows)} 行，不足以做 burn-in 与分位数')

    warnings = []
    pop_idx = {}
    grp_idx = {}
    for i, h in enumerate(hdr):
        m = _POP_RE.match(h.strip())
        if m:
            pop_idx[int(m.group(1) or 1)] = i
            continue
        m = _GRP_RE.match(h.strip())
        if m:
            grp_idx[int(m.group(1) or 1)] = i
    if not pop_idx:
        cand = [h for h in hdr if 'pop' in h.lower() or 'skyline' in h.lower()]
        raise RuntimeError(
            '日志里没有 `skyline.popSize*` 列 —— 这不是一份 BSP（贝叶斯天际线）日志。'
            f'疑似相关列：{cand or "无"}；前 12 个列名：{hdr[:12]}')

    keys = sorted(pop_idx)
    if len(keys) == 1 and keys[0] == 1 and 1 not in grp_idx:
        # 只有单段（skyline.popSize，无编号）→ 这是"恒定种群"式输出，不是天际线
        warnings.append('日志只有单段 `skyline.popSize`（没有分组），天际线会退化成'
                        '一条水平线 —— 确认你跑的是 BSP 而不是 Coalescent:Constant Size')

    n0 = len(rows)
    cut = int(n0 * float(burnin))
    if cut >= n0 - 10:
        raise RuntimeError(f'burn-in {burnin:.0%} 之后只剩 {n0 - cut} 行 —— 比例过大')
    body = rows[cut:]
    if burnin > 0:
        warnings.append(f'已按 burn-in {burnin:.0%} 丢弃前 {cut} 行（{n0} → {len(body)}）')

    # 段序号 → 时间区间（近似口径）
    def col_float(idx):
        out = []
        for r in body:
            try:
                out.append(float(r[idx]))
            except (ValueError, IndexError):
                out.append(float('nan'))
        return np.asarray(out, dtype=float)

    heights = None
    for c in _HEIGHT_COLS:
        try:
            j = _pick_col(hdr, (c,), '树高')
            heights = col_float(j)
            warnings.append(f'时间轴用 `{hdr[j]}` 作 tMRCA')
            break
        except RuntimeError:
            continue

    groups = []
    times = []
    use_time = (time_scale != 'group') and bool(grp_idx) and heights is not None
    if time_scale == 'auto' and not (bool(grp_idx) and heights is not None):
        warnings.append(
            '缺 `skyline.groupSize*` 或树高列 → **无法换算真实时间轴**，'
            '横轴改用「第几段」（段序号）。要出真实年份轴，'
            '请提供 Tracer 导出的 Time/Lower/Median/Upper 表（`table` 模式）。')
    if use_time:
        warnings.append('⚠️ 时间轴为近似：按"溯祖事件在时间上均匀分布"折算'
                        '（tMRCA × cumsum(groupSize) / sum(groupSize)），'
                        '**不是 BSP 的定义**；正式引用请用 Tracer 重建结果或 .trees 逐样本重建。')

    gsum = None
    gcum = None
    if use_time:
        gvals = [col_float(grp_idx[k]) for k in sorted(grp_idx)]
        gsum = np.nansum(np.vstack(gvals), axis=0)
        gcum = np.cumsum(np.vstack(gvals), axis=0)

    for pos, k in enumerate(keys, start=1):
        v = col_float(pop_idx[k])
        seg = {'index': pos, 'name': hdr[pop_idx[k]],
               'median': _percentile(v, 0.5),
               'lower': _percentile(v, (1 - ci) / 2.0),
               'upper': _percentile(v, 1 - (1 - ci) / 2.0),
               'mean': float(np.nanmean(v)) if np.isfinite(v).any() else float('nan'),
               'ess': effective_sample_size(v)}
        if use_time and gcum is not None and k in sorted(grp_idx):
            gi = sorted(grp_idx).index(k)
            with np.errstate(invalid='ignore', divide='ignore'):
                end = heights * gcum[gi] / gsum
            prev = heights * (gcum[gi - 1] if gi > 0 else 0.0) / gsum
            seg['t_start'] = float(np.nanmedian(prev))
            seg['t_end'] = float(np.nanmedian(end))
            seg['t_mid'] = (seg['t_start'] + seg['t_end']) / 2.0
        else:
            seg['t_start'] = seg['t_end'] = seg['t_mid'] = float(pos)
        groups.append(seg)

    groups.sort(key=lambda s: s['t_mid'])
    for s in groups:
        times.append(s['t_mid'])
    out = {'time': times,
           'median': [s['median'] for s in groups],
           'lower': [s['lower'] for s in groups],
           'upper': [s['upper'] for s in groups],
           'groups': groups,
           'n': len(groups),
           'n_dropped': 0,
           'columns': hdr,
           'time_axis': 'approx-uniform-coalescent' if use_time else 'group-index',
           'source': os.path.basename(log_path),
           'n_samples': len(body), 'burnin': float(burnin), 'ci': float(ci),
           'warnings': warnings}
    _check_skyline(out, warnings)
    return out


def _check_skyline(res, warnings):
    """一致性自检：区间上下颠倒 / median 越界 / ESS 过低。"""
    bad = [g['name'] for g in res.get('groups', [])
           if not (g['lower'] <= g['median'] <= g['upper'])]
    if bad:
        warnings.append(f'{len(bad)} 段的 Lower/Median/Upper 顺序不成立（{bad[:3]}…）'
                        '—— 检查是不是把列认错了')
    low = [g['name'] for g in res.get('groups', [])
           if np.isfinite(g['ess']) and g['ess'] < 100]
    if low:
        warnings.append(f'{len(low)} 段的 ESS < 100（{low[:3]}…）→ 该段后验没跑够，'
                        '分位数不可信，建议加长链或调 ESS 阈值')
    if any(g['median'] <= 0 for g in res.get('groups', [])):
        warnings.append('有段的中位数 ≤ 0 —— 对数纵轴会画不出来，检查输入单位')
    return warnings


# ---------------------------------------------------------------- 写产物

def write_skyline_table(res, out_path):
    """把结果写成 ``Time/Lower/Median/Upper`` 表 —— **与上游 TSV 同格式**，
    所以可以直接丢给 VirPhyKit 的 BSP-Viz，也可以当作再分析输入。

    ⚠️ 用 ``%.10g`` 而不是 ``%.6g``：上游 R 的 `read.table` 默认按 double 读，
    6 位有效数字会把大数只保留到个位/十位（如 10005.9 → 10006），
    写出去再读回来就**对不上**了（测试里就是靠 round-trip 抓到这个的）。
    """
    out_path = os.path.abspath(out_path)
    with safe_open(out_path, 'w') as f:
        f.write('Time\tLower\tMedian\tUpper\n')
        for t, lo, md, up in zip(res['time'], res['lower'], res['median'], res['upper']):
            f.write(f'{t:.10g}\t{lo:.10g}\t{md:.10g}\t{up:.10g}\n')
    return out_path


def render_pdf(res, out_path, direction='forward', color=None,
               secondary_axis=False, xlim=None, title='Bayesian Skyline Plot',
               ylabel='Ne·τ (years)', xlabel='Sampling Year'):
    """matplotlib 出静态图（PDF/PNG，按扩展名）。Agg 后端，不弹窗。

    与上游 R 版的差别（都是有意的）：
    · x 轴**自适应** —— 上游写死 ``xlim=c(1974,2016)``，换个数据集就画不出来；
    · 纵轴按上游同样思路**对数 + 留白 1.6 倍**（上游把 log 范围乘 1.6 再转回）；
    · ``secondary_axis=True`` 时才画次刻度（上游同）。
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from Virus_Platform_Core.sci_palette import OKABE

    if color is None:          # 缺省天际线色：Okabe-Ito 蓝（原 #00A0E9 青不属于任何真源色板）
        color = OKABE['blue']

    t = list(res['time'])
    lo, md, up = list(res['lower']), list(res['median']), list(res['upper'])
    order = sorted(range(len(t)), key=lambda i: t[i])
    t = [t[i] for i in order]
    lo = [lo[i] for i in order]
    md = [md[i] for i in order]
    up = [up[i] for i in order]

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.fill_between(t, lo, up, color=color, alpha=0.35, linewidth=0)
    ax.plot(t, md, color='black', lw=2.6)
    ax.set_yscale('log')
    if xlim:
        ax.set_xlim(*xlim)
    elif direction == 'reverse':
        ax.set_xlim(max(t), min(t))
    else:
        ax.set_xlim(min(t), max(t))
    if ylabel:
        ax.set_ylabel(ylabel)
    if xlabel:
        ax.set_xlabel(xlabel)
    ax.set_title(title)
    if secondary_axis:
        ax.grid(True, which='major', axis='x', alpha=0.3)
    fig.tight_layout()
    out_path = os.path.abspath(out_path)
    fmt = 'pdf' if out_path.lower().endswith('.pdf') else 'png'
    buf = io.BytesIO()
    fig.savefig(buf, format=fmt, dpi=150, bbox_inches='tight')
    plt.close(fig)
    with safe_open(out_path, 'wb') as f:
        f.write(buf.getvalue())
    return out_path


def audit_skyline(res):
    """面向页面的体检（不抛异常，只回 ok/problems/warnings）。"""
    problems, warnings = [], list(res.get('warnings') or [])
    if res.get('n', 0) < 2:
        problems.append('有效段数 < 2 —— 天际线至少要有两段才画得出来')
    if res.get('time_axis') == 'group-index':
        problems.append('横轴是段序号而不是时间：无法与采样年对应，'
                        '请提供 Tracer 导出的表或 BSP 日志（含 groupSize + 树高）')
    if res.get('time_axis') == 'approx-uniform-coalescent':
        warnings.append('时间轴为**近似**（均匀溯祖假设）；'
                        '论文里若要用严格口径，请给 Tracer 重建的 TSV')
    if res.get('time_axis') == 'as-given' and res.get('n_dropped'):
        warnings.append(f'输入表里有 {res["n_dropped"]} 行四列不完整，已跳过')
    _check_skyline(res, warnings)
    return {'ok': not problems, 'problems': problems, 'warnings': warnings}
