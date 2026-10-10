# -*- coding: utf-8 -*-
r"""BEAST 接入层 —— 把服务器端 BEAST 跑出来的产物读进来（导入 + 核验 + 复用）。

定位（与平台其它"外部工具接入"同一纪律）：**本地平台不做 MCMC**。重后验
（长链 BEAST）在服务器/集群上跑，本层负责把它变成平台自己的数据结构：

1. **`.log`（BEAST tab 日志）→ 诊断汇总**：逐列后验均值/中位数/95% 分位 +
   ESS，并按列名分组（likelihood / prior / posterior / 树高 / 速率 /
   指示 / 跳转 / 停留时间）。ESS 用 bsp.effective_sample_size（Geyer 初始正
   序列，Tracer 同一族），阈值 200 与 Tracer 惯例一致。
2. **跳转计数列（`c_A-to-B[1]` / `A-to-B`）→ 走廊后验表**：逐走廊
   `P(跳转 > 0)`（即**后验使用概率**——走廊在链上被真实跳转过、非 bootstrap
   复本频率）、均值 /
   中位数 / 95% 区间 / 最大值。数值口径与 mjrm 工具**同一套解析**
   （mjrm.parse_beast_log + _classify_columns），同一份日志在两个面板里
   不会给出两个数。
3. **`.trees`（后验样本树组或 TreeAnnotator MCC 树）→ 走廊后验表**：
   * ≥2 棵树（后验样本组）：逐树数走廊 → `p_used` = 走廊出现在**样本树的
     比例**（真后验频率，非 bootstrap 复本比例），
     `mean_events` = 逐树计数均值，95% 区间取逐树计数分位；并给出
     `sample_events`（逐树的 [距根距离, 状态i, 状态j]）→ 后验时间轴动画。
   * 单棵树（MCC / TreeAnnotator 标注树）：没有"逐样本"，改用**枝端后验**：
     `support` = 承载该走廊的各枝两端节点后验的**下界**（最不确定的那个
     节点），`mean_events` = Σ P(父=a)·P(子=b)（端独立近似）。这时
     **不给**动画（没有逐样本事件，不伪造帧）。
4. **可选运行**：探测本机 BEAST（probe_beast），有就代跑一份 XML
   （run_beast，stdout 边跑边回吐、可取消），跑完自动回到第 1–3 步解析。

⚠️ 措辞纪律：
  · 本层读的**就是** BEAST 产物，可以明说来源，但只说列里/树里**有什么**，
    不引申；
  · **不计算贝叶斯因子**：`*_indicator` 列只原样报后验均值（BSSVS 的 BF
    是另一套计算，本层不做，也不拿它跟平台的分带权重比大小）；
  · MCC 单树模式下，`support` 是**枝端后验下界**（一个保守的可靠性指标），
    不是"该走廊的后验概率"—— 两句话在界面上不能混用。
"""

import io
import math
import os
import re
import shutil
import subprocess
import sys
import time

import numpy as np

from Virus_Platform_Core.utils import check_path

# BEAST1 CLI 是 beast / beast.exe；BEAST2 的 CLI 一般是 beast（同一个名字），
# 打包版另有 beast2.exe。逐个试，Path 顺序优先。
_BEAST_EXES = ('beast', 'beast.exe', 'beast2', 'beast2.exe', 'BEAST')

# 列名 → 分组的候选（大小写/空格不敏感，先精确后子串）
_GROUP_EXACT = {
    'posterior': 'posterior', 'likelihood': 'likelihood', 'prior': 'prior',
}
_RATE_HINT = re.compile(r'(^|[._-])rate(s)?(\d+)?$', re.I)
_IND_HINT = re.compile(r'indicator(s)?(\d+)?$', re.I)


# ---------------------------------------------------------------- 环境探测

def probe_beast(exe=None, timeout=30):
    """探测 BEAST 可执行文件 → {available, path, version, note, exe_tried}。

    `exe` 为空时按 `_BEAST_EXES` 顺序在 Path 里找（也接受显式给一个路径）。
    版本用 `-version` 取（BEAST1/BEAST2 都认这个开关）；取不到不影响
    `available=True`，只在 `note` 里说明 —— 探测失败**不静默**，也不假装
    跑得起来。整个函数不抛异常（探测本身不该污染调用方）。
    """
    tried = []
    cands = [exe] if exe else list(_BEAST_EXES)
    for c in cands:
        if not c:
            continue
        tried.append(c)
        p = c if os.path.isfile(c) else shutil.which(c)
        if not p:
            continue
        ver, err = None, None
        try:
            r = subprocess.run([p, '-version'], capture_output=True, text=True,
                               timeout=timeout)
            out = ((r.stdout or '') + '\n' + (r.stderr or '')).strip()
            for ln in out.splitlines():
                if ln.strip():
                    ver = ln.strip()[:200]
                    break
        except Exception as e:            # noqa: BLE001 — 探测期一律吞
            err = f'{type(e).__name__}: {e}'
        note = '已探测到 BEAST：%s%s' % (p, ('（%s）' % ver) if ver else '')
        if err:
            note += f'；版本探测失败（{err}）—— 仍可尝试运行'
        return {'available': True, 'path': p, 'version': ver,
                'exe_tried': tried, 'note': note}
    return {'available': False, 'path': None, 'version': None,
            'exe_tried': tried,
            'note': ('本机没探测到 BEAST（试过 %s）——导入产物仍然可用：'
                     '把服务器端跑出来的 .log / .trees 填进输入框即可；'
                     '要在这里直接跑，请把 beast 放进 Path 或给出完整路径。'
                     % '、'.join(tried))}


# ---------------------------------------------------------------- 表格读取

def read_table(path):
    """读 BEAST 表格日志 → ``(cols, rows, meta)``。

    优先 BEAST 原生口径：跳过 `#` 注释，**列头行以 `state` 开头**（BEAST1/2
    的 .log / .rates.log / .states.log 都是这样）。找不到 `state` 列头时退回
    通用表格口径（首条非注释行当列头，tab/空白分隔）并在 ``meta['header_kind']``
    里标成 ``'generic'`` —— 调用方可据此提醒"这不像 BEAST 原生日志"，
    但**不因此拒绝**（Tracer 导出的表、脚本重排过的表同样有分析价值）。
    """
    path = check_path(path, must_exist=True)
    with io.open(path, encoding='utf-8', errors='replace') as f:
        lines = [ln.rstrip('\n').rstrip('\r') for ln in f]
    body = [(i, ln) for i, ln in enumerate(lines)
            if ln.strip() and not ln.lstrip().startswith('#')]
    if not body:
        raise RuntimeError(f'{os.path.basename(path)} 里没有数据行'
                           '（只有注释/空行）')
    hdr_i, hdr_line = body[0]
    kind = 'beast-state'
    if hdr_line.strip().startswith('state\t') or \
            hdr_line.strip().startswith('state '):
        pass
    else:
        kind = 'generic'
    sep = '\t' if hdr_line.count('\t') >= max(1, hdr_line.count(' ') // 4) \
        else None
    cols = ([c.strip() for c in hdr_line.split('\t')] if sep
            else hdr_line.split())
    rows = []
    for _i, ln in body[1:]:
        parts = ([c.strip() for c in ln.split('\t')] if sep else ln.split())
        if len(parts) < len(cols):
            parts += [''] * (len(cols) - len(parts))
        rows.append(parts)
    meta = {'file': os.path.basename(path), 'path': path, 'header_kind': kind,
            'n_columns': len(cols), 'n_rows': len(rows),
            'sep': 'tab' if sep else 'space'}
    if kind == 'generic':
        meta['note'] = ('列头不以 `state` 开头 —— 不像是 BEAST 原生日志，'
                        '按通用表格读（列名/分隔符已自动判断）')
    return cols, rows, meta


def _rows_with_burnin(path, burnin):
    """(cols, rows, meta)，按 ``burnin`` 比例丢**开头**的采样点。

    走 mjrm 的严格口径（与「Markov jump 迁移计数」工具同一套解析，保证同一
    份日志在两处给出同一个数）；它认不出 `state` 列头时退回 read_table。
    """
    from Virus_Platform_Core.mjrm import parse_beast_log
    b = min(0.9, max(0.0, float(burnin or 0.0)))
    try:
        cols, rows, skipped = parse_beast_log(path, b)
        return cols, rows, {'file': os.path.basename(path), 'path': path,
                            'header_kind': 'beast-state', 'reader': 'mjrm',
                            'n_skipped': skipped,
                            'n_rows': len(rows) + skipped}
    except RuntimeError:
        cols, rows, meta = read_table(path)
        cut = int(len(rows) * b)
        meta = dict(meta, reader='generic', n_skipped=cut,
                    n_rows=len(rows))
        return cols, rows[cut:], meta


# ---------------------------------------------------------------- 列统计

def _col_values(rows, k):
    """第 k 列 → (float 数组, n_非数值)。BEAST 用 `-` 表示"该列这次没采到"。"""
    vals, bad = [], 0
    for r in rows:
        if k >= len(r):
            continue
        v = r[k].strip()
        try:
            vals.append(float(v))
        except ValueError:
            bad += 1
    return np.asarray(vals, dtype=float), bad


def _stats(vals, ci=0.95):
    a = np.asarray(vals, dtype=float)
    a = a[np.isfinite(a)]
    if a.size == 0:
        return None
    lo_q = (1.0 - float(ci)) / 2.0 * 100.0
    hi_q = (1.0 + float(ci)) / 2.0 * 100.0
    from Virus_Platform_Core.bsp import effective_sample_size
    return {'n': int(a.size),
            'mean': round(float(a.mean()), 6),
            'median': round(float(np.median(a)), 6),
            'min': round(float(a.min()), 6),
            'max': round(float(a.max()), 6),
            'ci95': [round(float(np.percentile(a, lo_q)), 6),
                     round(float(np.percentile(a, hi_q)), 6)],
            'ess': (None if a.size < 4
                    else round(float(effective_sample_size(a)), 1))}


def _group_of(name, route_keys, reward_keys, total_keys=()):
    """列名 → 分组 key（只按**名字**分组，不改数值）。

    jump/dwell/total 列要先剥掉 BEAST 自己加的外壳再比对走廊/性状名：
    `c_` 前缀、`[1]` 下标（`c_AS-to-EC[1]`）、`_reward` 后缀
    （`c_AS_reward[1]`）、`.count` 后缀（`c_AS.count[1]`）—— 原样比对会把
    实测最常见的这几种写法全漏到 `other` 里。
    """
    low = name.strip().lower().replace(' ', '')
    if low in _GROUP_EXACT:
        return _GROUP_EXACT[low]
    if low.endswith('indicator') or low.endswith('indicators') \
            or _IND_HINT.search(low):
        return 'indicator'
    core = re.sub(r'\[\d+\]$', '', name.strip())
    core = re.sub(r'^c_', '', core)
    for key in route_keys:
        if core == key or core.endswith(key):
            return 'jump'
    dwell_core = re.sub(r'_reward$', '', core)
    for key in reward_keys:
        if dwell_core == key or dwell_core.endswith(key):
            return 'dwell'
    total_core = re.sub(r'\.count$', '', core)
    if total_core in total_keys:
        return 'total'
    if _RATE_HINT.search(low) or 'clock.rate' in low:
        return 'rate'
    return 'other'


def log_summary(path, burnin=0.1, ci=0.95):
    """BEAST `.log` → 诊断汇总（逐列后验 + ESS + 分组 + 告警）。

    返回 dict：``columns``（数值列全量：name/group/mean/median/ci95/ess/n）、
    分组子集（``groups`` 里的 posterior / likelihood / prior / rate /
    indicator / jump / dwell / total / other）、
    ``tree_height``（tMRCA 列，供时间标度核对）、``ess_min`` / ``n_low_ess``、
    ``n_nonnumeric``（`-` 那样的缺采样次数，逐列非数值总数）、``warnings``。

    ⚠️ 只做**汇总**，不做模型比较：这里没有贝叶斯因子，`indicator` 列的
    均值就是"该指示变量被采样为 1 的比例"，到此为止。
    """
    cols, rows, meta = _rows_with_burnin(path, burnin)
    if not rows:
        raise RuntimeError(f'{os.path.basename(path)} burn-in {burnin:.0%} '
                           f'之后没有数据行（总行数 {meta.get("n_rows")}）')
    warnings = []
    if meta.get('note'):
        warnings.append(meta['note'])
    if meta.get('n_skipped'):
        warnings.append('已按 burn-in %.0f%% 丢弃前 %d 行（剩 %d 行）'
                        % (float(burnin) * 100.0, meta['n_skipped'], len(rows)))

    from Virus_Platform_Core import mjrm
    try:
        cls = mjrm._classify_columns(cols)
    except Exception:                          # noqa: BLE001 — 分组失败不致命
        cls = {'route': {}, 'reward': {}, 'total': {}}
    route_keys = {f'{a}-to-{b}' for (a, b) in cls.get('route', {})}
    reward_keys = set(cls.get('reward', {}))
    total_keys = set(cls.get('total', {}))

    columns, nonnum = [], 0
    for k, name in enumerate(cols):
        vals, bad = _col_values(rows, k)
        nonnum += bad
        st = _stats(vals, ci=ci)
        if st is None:
            continue
        st['name'] = name
        st['group'] = _group_of(name, route_keys, reward_keys, total_keys)
        st['n_nonnum'] = bad
        columns.append(st)
    if not columns:
        raise RuntimeError(f'{os.path.basename(path)} 里没有可解析的数值列'
                           f'（列名：{cols[:12]}）')

    low = [c for c in columns
           if c.get('ess') is not None and c['ess'] < 200]
    low.sort(key=lambda c: c['ess'])
    ess_min = (dict(name=low[0]['name'], ess=low[0]['ess']) if low else None)
    if low:
        warnings.append('ESS < 200 的列 %d 个（最低 `%s` = %.0f）—— '
                        '该链对这组量还没跑稳，均值和区间都可能不可靠'
                        % (len(low), low[0]['name'], low[0]['ess']))
    if nonnum:
        warnings.append('有 %d 个格子是非数值（BEAST 用 `-` 表示该列这次没采到）'
                        '——逐列统计已跳过这些格子' % nonnum)
    if len(rows) < 200:
        warnings.append('有效采样点只有 %d 个 —— 后验区间偏粗，'
                        '建议加长链或减小 logEvery' % len(rows))
    if not any(c['group'] == 'posterior' for c in columns):
        warnings.append('没有 `posterior` 列 —— 这不像标准 BEAST 日志'
                        '（或该 <log> 只记录了一部分参数），核对一下来源')

    groups = {}
    for key in ('posterior', 'likelihood', 'prior', 'rate', 'indicator',
                'jump', 'dwell', 'total', 'other'):
        groups[key] = [c for c in columns if c['group'] == key]
    from Virus_Platform_Core.bsp import _HEIGHT_COLS
    th = None
    want = [w.lower().replace(' ', '') for w in _HEIGHT_COLS]
    for c in columns:
        nm = c['name'].strip().lower().replace(' ', '')
        if nm in want or any(nm.endswith(w) for w in want):
            th = c
            break
    return {
        'file': meta['file'],
        'header_kind': meta['header_kind'],
        'reader': meta.get('reader'),
        'n_total': meta.get('n_rows'),
        'n_samples': len(rows),
        'n_skipped': meta.get('n_skipped'),
        'burnin': float(burnin),
        'ci': float(ci),
        'columns': columns,
        'groups': groups,
        'tree_height': th,
        'ess_min': ess_min,
        'n_low_ess': len(low),
        'n_nonnumeric': nonnum,
        'warnings': warnings,
    }


# ---------------------------------------------------------------- 跳转日志 → 走廊

_J_EPS = 1e-9


def jump_corridors(path, burnin=0.1, ci=0.95):
    """BEAST 跳转日志（`c_A-to-B[1]` 那些列）→ 走廊后验表。

    每条走廊给：
      ``p_used``     = **P(该走廊跳转数 > 0)**（后验样本里的比例）
      ``mean_events``/``median``/``ci95``/``max_events`` = 跳转计数的后验
      ``n_zero``     = 采样点里该走廊一次都没跳的次数
      ``n_obs``      = 该列的有效采样点数（BEAST 用 `-` 的格子不算）

    ⚠️ 跳转计数是**采样历史下的期望跳转次数**（不是树上边数、也不是
    BSSVS 的速率指示），与 Fitch 观测事件数是两条口径 —— 面板上分开摆。
    """
    cols, rows, meta = _rows_with_burnin(path, burnin)
    if not rows:
        raise RuntimeError(f'{os.path.basename(path)} burn-in 之后没有数据行')
    from Virus_Platform_Core import mjrm
    cls = mjrm._classify_columns(cols)
    route, reward, total = cls['route'], cls['reward'], cls['total']
    if not route:
        cand = [c for c in cols if 'indicator' in c.lower()]
        raise RuntimeError(
            f'{os.path.basename(path)} 里没有按走廊命名的跳转计数列'
            '（`c_A-to-B[1]` / `A-to-B`）—— 给不出走廊后验表。'
            + (f' 见到 `*indicator*` 列：{cand[:6]}；'
               'BEAST 把矩阵写成 `X.indicators1..N` 扁平向量时**不猜顺序**'
               '（猜错会把走廊全标反），请改用带走廊名的日志或 MCC 树。'
               if cand else f' 可用列（前 20）：{cols[:20]}'))
    warnings = []
    if meta.get('n_skipped'):
        warnings.append('已按 burn-in %.0f%% 丢弃前 %d 行（剩 %d 行）'
                        % (float(burnin) * 100.0, meta['n_skipped'], len(rows)))

    corridors = []
    for (src, dst), k in sorted(route.items()):
        vals, bad = _col_values(rows, k)
        if vals.size < 5:
            warnings.append(f'{src}→{dst} 有效采样仅 {vals.size} 个，已跳过')
            continue
        st = _stats(vals, ci=ci)
        p_used = float(np.mean(vals > _J_EPS))
        corridors.append({
            'from': src, 'to': dst,
            'p_used': round(p_used, 4),
            'mean_events': st['mean'], 'median': st['median'],
            'ci95': st['ci95'], 'max_events': st['max'],
            'n_zero': int(np.sum(vals <= _J_EPS)), 'n_obs': st['n'],
            'n_nonnum': bad,
        })
    corridors.sort(key=lambda c: (-c['p_used'], -c['mean_events'],
                                  str(c['from']), str(c['to'])))

    dwell = {}
    for name, k in sorted(reward.items()):
        st = _stats(_col_values(rows, k)[0], ci=ci)
        if st:
            dwell[name] = st
    tot = sum(v['mean'] for v in dwell.values())
    if tot > 0:
        for name in dwell:
            dwell[name]['fraction'] = round(dwell[name]['mean'] / tot, 4)

    total_mean = None
    total_source = 'c_<性状>.count 列'
    if total:
        nm = next(iter(total))
        st = _stats(_col_values(rows, total[nm])[0], ci=ci)
        if st:
            total_mean = st['mean']
    if total_mean is None:
        total_mean = round(sum(c['mean_events'] for c in corridors), 6)
        total_source = '逐走廊求和（该日志没有 `c_<性状>.count` 列）'
        warnings.append('日志里没有 `c_<性状>.count` 列，总跳转数由逐走廊求和得到')
    if len(rows) < 200:
        warnings.append('有效采样点只有 %d 个 —— 后验区间偏粗，'
                        '建议加长链或减小 logEvery' % len(rows))
    return {'corridors': corridors, 'n_samples': len(rows),
            'burnin': float(burnin), 'ci': float(ci),
            'total_mean': total_mean, 'total_source': total_source,
            'dwell': dwell,
            'route_cols': len(route), 'warnings': warnings,
            'file': os.path.basename(path)}


# ---------------------------------------------------------------- 树组 → 走廊

def _tree_fmt(head):
    return 'nexus' if head.lstrip().upper().startswith('#NEXUS') else 'newick'


# BEAST 1.x 的 TreeLog 把树级注释写在等号**前**：
#   tree STATE_1 [&lnP=…,posterior=…] = [&R] (…)
# Bio.Nexus 只认 MrBayes 式 `tree NAME = [&注释]`，会 NexusError（真例：
# git-repo WNV_BSG_tree_*.trees）。BEAST 2 / TreeAnnotator 的输出版式不受影响。
_TREE_PRE_EQUALS = re.compile(
    r'^([ \t]*(?:tree|utree)\s+\S+)[ \t]*(\[[^\]]*\])[ \t]*=', re.M | re.I)


def _load_trees(path):
    """读树文件 → (trees, fmt, layout_fix)。外部文件只读，改写只在内存里做。

    两个现实版式的坑：
      ① BEAST 1.x TreeLog 的树级注释在等号前（见 _TREE_PRE_EQUALS）——NEXUS
         解析直接报错；失败时把注释平移过等号重解析。
      ② Bio.Nexus 用 safename() 回填 Translate 标签，含 `-` 等标点的名字会被
         重新加引号（`'AF…-73.899889_…'`）——带引号的名字会让一切按名查找
         落空（PEDV 案的同类坑），统一剥掉。
    """
    from Bio import Phylo
    with io.open(path, encoding='utf-8', errors='replace') as f:
        fmt = _tree_fmt(f.read(4096))
    layout_fix = False
    try:
        trees = list(Phylo.parse(path, fmt))
    except Exception:
        if fmt != 'nexus':
            raise
        with io.open(path, encoding='utf-8', errors='replace') as f:
            text = f.read()
        text, n_fix = _TREE_PRE_EQUALS.subn(r'\1 = \2', text)
        if not n_fix:
            raise
        layout_fix = True
        trees = list(Phylo.parse(io.StringIO(text), fmt))
    for t in trees:
        for cl in t.find_clades():
            nm = (cl.name or '').strip()
            if len(nm) >= 2 and nm[0] == nm[-1] and nm[0] in ('"', "'"):
                cl.name = nm[1:-1].strip() or None
    return trees, fmt, layout_fix


def _known_states(tree, max_nodes=400):
    """从标注注释里收集状态全集（供 mot.parse_node_states 的回退层用）。

    两路：`max.set={…}`（BEAST1 逐状态格式的带头键）优先；没有再退一步扫
    任意 `<X>.set={…}`（TreeAnnotator 对离散性状也写这种）。两路都没有时
    返回空 —— 那说明树里只有 `max`/`location` 这类整体性状（mot 的 ①②③ 层
    照样能解析），不猜。
    """
    out = set()
    for i, cl in enumerate(tree.find_clades()):
        if i >= max_nodes:
            break
        c = getattr(cl, 'comment', '') or ''
        m = re.search(r'\bmax\.set=\{([^}]*)\}', c)
        if m:
            out |= {x.strip().strip('"') for x in m.group(1).split(',')
                    if x.strip()}
    if out:
        return sorted(out)
    for i, cl in enumerate(tree.find_clades()):
        if i >= max_nodes:
            break
        c = getattr(cl, 'comment', '') or ''
        for m in re.finditer(r'(?:^|[,&])\s*([A-Za-z][\w%.\-]*)\.set='
                             r'\{([^}]*)\}', c):
            out |= {x.strip().strip('"') for x in m.group(2).split(',')
                    if x.strip()}
        if out:
            break
    return sorted(out)


def _node_heights(tree, info):
    """节点高度（距末次采样的时间）→ (rh, dist_from_root 字典)。

    与 mot 同一套回退：先认 `height=` 注释（BEAST 的 .trees 每个节点都带），
    缺的用**累积枝长**从根推（rh − cum）。返回的 dist 是"距根距离"，
    也就是**从根往现在走**的坐标（根 = 0，最新采样 = rh）——后验时间轴
    动画的横轴就是这个距离（不是日历年）。
    """
    hs = [v['height'] for v in info.values() if v['height'] is not None]
    if hs:
        rh = max(hs)
        for v in info.values():
            if v['height'] is None:
                v['height'] = rh - v['cum']
    else:
        rh = max(v['cum'] for v in info.values()) if info else 0.0
        for v in info.values():
            v['height'] = rh - v['cum']
    dist = {k: max(0.0, rh - v['height']) for k, v in info.items()}
    return rh, dist


def tree_corridors(path, burnin=0.0, max_trees=2000, ci=0.95,
                   max_anim_trees=500, n_bins=40, dates=None,
                   min_prob=0.0, progress=None):
    """`.trees` / MCC 树 → 走廊后验表（两种模式，见模块头）。

    burnin      丢**开头**的树比例（后验样本组才有意义；单树时忽略）
    max_trees   最多用多少棵树（超出如实上报 n_trees_used < n_trees）
    max_anim_trees  进 sample_events（动画用）的树上限
    min_prob    只统计"父子两端后验都 ≥ 该阈值"的枝（可靠性敏感性；0=不过滤）

    返回 dict：
      mode='sample'（≥2 棵树）/ 'mcc'（单树）/ 'empty'
      corridors[]  from/to/p_used/mean_events/ci95/max_events/n_events
                   （mcc 模式另给 support=枝端后验下界、p_used=同一数值，
                    见模块头的措辞纪律）
      labels[]     状态序（动画的下标口径）
      n_trees / n_trees_used / n_branches / n_unresolved / n_same_state
      tree_height / time_hist / sample_events / anim_truncated
      warnings[]
    """
    from Virus_Platform_Core import mot
    path = check_path(path, must_exist=True)
    trees, fmt, layout_fix = _load_trees(path)
    if not trees:
        raise RuntimeError(f'无树可解析（{fmt}）：{path}')
    warnings = []
    if layout_fix:
        warnings.append('这是 BEAST 1.x 版式的树组（树级注释写在 `=` 前），'
                        '已按原样改写后解析（外部文件未动）')
    n_all = len(trees)
    b = min(0.9, max(0.0, float(burnin or 0.0)))
    if n_all > 1 and b > 0:
        cut = int(n_all * b)
        trees = trees[cut:]
        warnings.append('已按 burn-in %.0f%% 丢弃前 %d 棵树（剩 %d 棵）'
                        % (b * 100.0, cut, len(trees)))
    if len(trees) > int(max_trees or 0) > 0:
        # 等距抽样（审查 P1-10）：旧实现无说明地取**前** max_trees 棵——
        # MCMC 链前段恰是收敛最差的部分，且按顺序截取没有抽稀依据。改为
        # 覆盖 burn-in 后全段的等距抽样（确定性、可复现），并在警告里写明。
        _n_burn = len(trees)
        _mt = int(max_trees)
        import numpy as _np
        _idx = _np.linspace(0, len(trees) - 1, _mt).round().astype(int)
        trees = [trees[i] for i in _idx]
        warnings.append('burn-in 后树数 %d 超过上限 %d，已用**等距抽样**取'
                        '%d 棵（覆盖全段，非简单截取前段）'
                        % (_n_burn, _mt, len(trees)))
    n_used = len(trees)
    dates = dict(dates or {})

    known = _known_states(trees[0])
    if not known:
        for t in trees[1:6]:
            known = _known_states(t)
            if known:
                break
    if not known:
        warnings.append('树里没解析出状态集合（`max.set` / `<X>.set` 都没有）'
                        '—— 只统计 `max`/`location` 这类整体性状')

    # 逐树解析（同一套 walk：mot.parse_node_states 的四层回退）
    mode = 'mcc' if n_used == 1 else 'sample'
    per_tree = []          # [{'counts','events','rh','max_dist','branches'}]
    n_unresolved = 0
    n_same = 0
    n_branches = 0
    n_lowprob = 0
    label_order = []
    for ti, tree in enumerate(trees):
        info = {}

        def walk(cl, parent=None, cum=0.0):
            stt = mot.parse_node_states(getattr(cl, 'comment', None), known)
            top = (max(stt, key=stt.get) if stt else None)
            info[id(cl)] = {'states': stt, 'loc': top,
                            'top_prob': (stt[top] if stt else None),
                            'height': mot.node_height(getattr(cl, 'comment',
                                                              None)),
                            'parent': parent, 'clade': cl,
                            'length': cl.branch_length or 0.0, 'cum': cum,
                            'is_tip': len(cl.clades) == 0}
            for ch in cl.clades:
                walk(ch, cl, cum + (ch.branch_length or 0.0))
        walk(tree.root)
        rh, dist = _node_heights(tree, info)
        counts, events, branches = {}, [], []
        for v in info.values():
            if v['clade'] is tree.root:
                continue
            p = v['parent']
            if p is None:
                continue
            pv = info[id(p)]
            # ⚠️ 两端都要检查（2026-10-03 修复）：旧实现只数"父端缺状态"，
            #    子端缺状态在循环开头就被 `continue` 掉了（而下游告警原文写的
            #    是"某一端没解析出状态"）→ n_unresolved 少报。现在按"任一端缺"
            #    统一计数，与告警口径一致。
            if pv['loc'] is None or v['loc'] is None:
                n_unresolved += 1
                continue
            n_branches += 1
            if pv['loc'] == v['loc']:
                n_same += 1
                continue
            a, b = pv['loc'], v['loc']
            pa = float(pv['top_prob'] or 0.0)
            pb = float(v['top_prob'] or 0.0)
            # ⚠️ 这里必须用**节点对象**的 id：info 的键是 id(clade)，
            # 写成 id(v)（v 是值字典）会查不到/查错（实测：真 MCC 树上
            # 204 个节点全部 KeyError）。
            t_mid = (dist[id(p)] + dist[id(v['clade'])]) / 2.0
            if float(min_prob or 0) > 0 and not (pa >= float(min_prob)
                                                and pb >= float(min_prob)):
                n_lowprob += 1
                continue
            counts[(a, b)] = counts.get((a, b), 0) + 1
            events.append((t_mid, a, b))
            if mode == 'mcc':
                branches.append((a, b, pa, pb, t_mid))
        per_tree.append({'counts': counts, 'events': events,
                         'rh': rh, 'max_dist': (max(dist.values())
                                                if dist else 0.0),
                         'branches': branches})
        for ev in events:
            for s in (ev[1], ev[2]):
                if s not in label_order:
                    label_order.append(s)

    corridors = []
    keys = set()
    for pt in per_tree:
        keys |= set(pt['counts'])
    if mode == 'sample':
        for (a, b) in keys:
            arr = np.asarray([pt['counts'].get((a, b), 0) for pt in per_tree],
                             dtype=float)
            st = _stats(arr, ci=ci)
            corridors.append({
                'from': a, 'to': b,
                'p_used': round(float(np.mean(arr > 0)), 4),
                'support': round(float(np.mean(arr > 0)), 4),
                'mean_events': st['mean'], 'median': st['median'],
                'ci95': st['ci95'], 'max_events': int(arr.max()),
                'n_zero': int(np.sum(arr <= 0)), 'n_events': int(arr.sum()),
                'n_obs': st['n'],
            })
        corridors.sort(key=lambda c: (-c['n_events'], -c['p_used'],
                                      str(c['from']), str(c['to'])))
    else:
        acc = {}
        for (a, b, pa, pb, _t) in per_tree[0]['branches']:
            d = acc.setdefault((a, b), {'n': 0, 'exp': 0.0, 'worst': 1.0,
                                        'best': 0.0})
            d['n'] += 1
            d['exp'] += pa * pb
            d['worst'] = min(d['worst'], pa, pb)
            d['best'] = max(d['best'], pa * pb)
        for (a, b), d in acc.items():
            corridors.append({
                'from': a, 'to': b,
                # mcc 模式下 support / p_used 是**同一根线**：枝端后验下界
                # （该走廊各枝两端里最不确定的那个节点）。字段名沿用 p_used
                # 是为了让前端/分带复用同一个字段位；driver_label / exp_note
                # 会把这句讲清楚，不许简写成"该走廊的后验概率"。
                'support': round(d['worst'], 4),
                'p_used': round(d['worst'], 4),
                'mean_events': round(d['exp'], 4),
                'median': None, 'ci95': None,
                'max_events': round(d['best'], 4),
                'n_events': d['n'], 'n_zero': 0, 'n_obs': 1,
            })
        corridors.sort(key=lambda c: (-c['mean_events'], -c['support'],
                                      str(c['from']), str(c['to'])))
        warnings.append('这是**单棵树**（MCC / TreeAnnotator 标注树）：没有逐样本，'
                        '所以支持率 = 承载该走廊的各枝端节点后验的**下界**'
                        '（最不确定的那个节点），期望事件数 = Σ P(父=a)·P(子=b)'
                        '（端独立近似）—— 与样本树组模式的"使用概率"不是同一件事，'
                        '那句口径只对样本组模式成立')
    if n_lowprob:
        warnings.append('按 min_prob=%.2f 过滤掉 %d 条低置信枝（两端后验必须'
                        '都过线）' % (float(min_prob), n_lowprob))
    if n_unresolved:
        warnings.append('%d 条枝的某一端没解析出状态（性状名不在支持范围内？）'
                        '—— 未计入走廊表' % n_unresolved)

    # 时间轴直方图 + 动画样本（只对样本组模式）
    time_hist, sample_events, anim_truncated = None, [], False
    if mode == 'sample' and per_tree:
        hi = max(max(pt['max_dist'], pt['rh']) for pt in per_tree)
        nb = max(1, int(n_bins))
        bw = hi / nb if hi > 0 else 1.0
        hist = np.zeros((n_used, nb), dtype=float)
        for ti, pt in enumerate(per_tree):
            for (t, _a, _b) in pt['events']:
                bi = min(nb - 1, max(0, int(t / bw)))
                hist[ti, bi] += 1.0
        mean = hist.mean(axis=0)
        lo = (np.percentile(hist, 2.5, axis=0) if n_used >= 2
              else np.zeros(nb))
        hiq = (np.percentile(hist, 97.5, axis=0) if n_used >= 2
               else np.zeros(nb))
        time_hist = {'lo': 0.0, 'hi': round(hi, 4), 'n_bins': nb,
                     'bin_width': round(bw, 4),
                     'axis': '距根距离（枝长单位）',
                     'mean': [round(float(x), 4) for x in mean],
                     'p2_5': [round(float(x), 4) for x in lo],
                     'p97_5': [round(float(x), 4) for x in hiq],
                     'per_tree': '每窗计数 = 单棵样本树上落在该窗的迁移事件数；'
                                 '均值/区间为逐树统计'}
        idx = {s: i for i, s in enumerate(label_order)}
        cap = int(max_anim_trees or 0)
        # 均匀抽（不是只取前 N 棵）：后验样本的链序不该被"只取开头"截断
        keep = ([i for i in range(0, n_used, max(1, n_used // cap))][:cap]
                if cap > 0 else list(range(n_used)))
        for ti in keep:
            pt = per_tree[ti]
            sample_events.append([[round(t, 4), idx[a], idx[b]]
                                  for (t, a, b) in pt['events']])
        anim_truncated = len(sample_events) < n_used
        if anim_truncated:
            warnings.append('动画只用了 %d / %d 棵样本树（逐树事件按 '
                            'max_anim_trees 截断）' % (len(sample_events),
                                                       n_used))
    return {
        'mode': mode, 'file': os.path.basename(path),
        'corridors': corridors, 'labels': label_order,
        'n_trees': n_all, 'n_trees_used': n_used,
        'n_branches': n_branches,
        'n_branches_per_tree': (round(n_branches / float(n_used), 2)
                                if n_used else None),
        'n_unresolved': n_unresolved, 'n_same_state': n_same,
        'n_lowprob': n_lowprob,
        'states': label_order,
        'tree_height': (round(per_tree[0]['rh'], 4) if per_tree else None),
        'time_hist': time_hist, 'sample_events': sample_events,
        'anim_truncated': anim_truncated,
        'max_anim_trees': int(max_anim_trees or 0),
        'min_prob': float(min_prob),
        'warnings': warnings,
    }


# ---------------------------------------------------------------- 可选运行

def xml_outputs(xml_path):
    """从 BEAST XML 里读 `fileName="…"` → 绝对路径清单（相对 XML 所在目录）。"""
    d = os.path.dirname(os.path.abspath(xml_path))
    with io.open(xml_path, encoding='utf-8', errors='replace') as f:
        text = f.read()
    names = re.findall(r'fileName\s*=\s*"([^"]+)"', text)
    out = []
    for n in names:
        p = n if os.path.isabs(n) else os.path.join(d, n)
        if p not in out:
            out.append(p)
    return out


def run_beast(xml_path, exe=None, seed=None, extra_args=(), timeout=None,
              on_line=None, cancel=None):
    """跑一份 BEAST XML（可选能力：没探测到 BEAST 就明确报错，不静默降级）。

    stdout+stderr 合并，逐行回吐给 ``on_line``（任务日志里能实时看到 BEAST
    的进度），``cancel()`` 为真时终止子进程。返回 dict：``ok`` / ``cmd`` /
    ``seconds`` / ``returncode`` / ``log_files`` / ``tree_files`` /
    ``new_files`` / ``tail``（最后 20 行）。

    ⚠️ 输出文件靠**两路合并**判定：XML 里 `fileName=` 声明的 + XML 目录下
    mtime 晚于启动时刻的新文件 —— 后者兜住"模板里改了输出名/多写了几份
    日志"的情况，前者兜住"输出落在别的目录"。两路都空 = 真没产出。
    """
    xml_path = check_path(xml_path, must_exist=True)
    pr = probe_beast(exe)
    if not pr['available']:
        raise RuntimeError(pr['note'])
    t0 = time.time()
    cmd = [pr['path']]
    if seed:
        cmd += ['-seed', str(int(seed))]
    cmd += [str(a) for a in (extra_args or ())]
    cmd.append(xml_path)
    lines = []

    def _emit(ln):
        ln = ln.rstrip('\r\n')
        if not ln:
            return
        lines.append(ln)
        if on_line:
            on_line(ln)

    # Windows 下不给 CREATE_NO_WINDOW：BEAST 是 GUI 程序（beast.exe 带窗口），
    # 代跑时会在桌面弹黑框/窗口干扰用户（同仓 utils/treetime_ml 的惯例）。
    creationflags = 0
    if sys.platform == 'win32':
        creationflags = getattr(subprocess, 'CREATE_NO_WINDOW', 0)
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, text=True,
                            encoding='utf-8', errors='replace', bufsize=1,
                            creationflags=creationflags)
    killed = False
    try:
        for ln in iter(proc.stdout.readline, ''):
            _emit(ln)
            if cancel and cancel():
                killed = True
                proc.kill()
                _emit('（已按取消请求终止 BEAST 进程）')
                break
            if timeout and (time.time() - t0) > float(timeout):
                killed = True
                proc.kill()
                _emit('（超过时限 %.0fs，已终止 BEAST 进程）' % float(timeout))
                break
    finally:
        try:
            proc.stdout.close()
        except Exception:                       # noqa: BLE001
            pass
        rc = proc.wait()
    secs = round(time.time() - t0, 2)

    declared = xml_outputs(xml_path)
    d = os.path.dirname(os.path.abspath(xml_path))
    fresh = []
    try:
        for nm in os.listdir(d):
            p = os.path.join(d, nm)
            if os.path.isfile(p) and os.path.getmtime(p) >= t0 - 1:
                fresh.append(p)
    except OSError:
        pass
    files = []
    for p in declared + fresh:
        if os.path.isfile(p) and p not in files:
            files.append(p)
    logs = [p for p in files if p.lower().endswith('.log')]
    trees = [p for p in files if p.lower().endswith(('.trees', '.tree', '.tre'))]
    ok = (rc == 0) and not killed and bool(logs or trees)
    # cmd 给 **列表**（消费端要能直接 join / 复制命令复现）；不要在这里拼成
    # 字符串——前端与 job 都会再 join 一次，拼成字符串会被逐字符拆开。
    return {'ok': bool(ok), 'killed': killed, 'returncode': rc,
            'cmd': list(cmd), 'beast': pr['path'],
            'version': pr['version'], 'seconds': secs,
            'files': [os.path.basename(p) for p in files],
            'log_files': logs, 'tree_files': trees,
            'new_files': [os.path.basename(p) for p in fresh],
            'xml': os.path.basename(xml_path),
            'tail': lines[-20:]}
