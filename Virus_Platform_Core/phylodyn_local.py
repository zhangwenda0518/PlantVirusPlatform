# -*- coding: utf-8 -*-
r"""本地时间与地理推断流水线（轨道 A）—— A0 → A5 六阶段。

## 为什么要有这个模块

平台的进化动力学此前是 **8 张各自为战的卡**：定年结果要人工搬去地理卡，
地理卡又不知道时间树从哪来；`t-mot`/`t-bsp`/`t-mjrm`/`t-beast` 四张卡还都
必须交服务器跑 BEAST。四仓调研（phymap-workflow / Mugration-Analysis /
pathogen-phylo-genomics-pipeline / MiniMev）提炼出的六个共同模式里，最关键的
一条是：

> **时间与地理是分开的阶段，靠「时间树」这一个接口产物衔接。**

本模块把这条落地成一条本地端到端链，接口产物固定为
**`clock/timetree.nwk`**（枝长＝年）。

## 阶段与输入输出契约

| 阶段 | 输入 | 输出 |
|---|---|---|
| **A0 prep** | 比对 FASTA（必需）+ 元数据（可选）+ 树（可选，缺则 FastTree 现场建） | `prep/dates.csv`、`prep/states.csv`、`prep/prep.fasta`、`prep/prep.nwk`、`prep/prep_report.json` |
| **A1 clock** | A0 全部 | **`clock/timetree.nwk`**、`clock/clock_stats.tsv`（**带方法列**）、`clock/rtt_scatter.tsv`、`clock/clock_qc/drt.json`、`clock/divergence_tree.nwk` |
| **A2 phylogeo** | `clock/timetree.nwk` + `prep/states.csv` | `phylogeo/GTR.txt`、`confidence.csv`、`annotated_tree.nexus`、**`corridors.tsv`**（ML 与 Fitch 并列）、**`events_year.tsv`**、`phylogeo_report.json` |
| **A3 coalescent** | `clock/timetree.nwk` + `prep/prep.fasta` | `coalescent/skyline.tsv`、`coalescent_report.json` |
| **A4 ancestral** | `prep/prep.fasta` + `clock/divergence_tree.nwk` | `ancestral/ancestral_sequences.fasta`、`branch_mutations.txt`、`auspice_tree.json` |
| **A5 report** | A0–A4 的产物 | **`run_summary.json`**、**`config.resolved.json`**、`RUN_REPORT.md` |

## 五条写死在代码里的红线

1. **每个数字都带方法名**（`clock_stats.tsv` 有 `method` 列，走廊表分 `ml`/`fitch` 两列）
2. **点估计与后验分开列表**，不混排名
3. **ML 边际概率 ≠ 贝叶斯后验**，措辞不简化
4. **采样偏倚检查默认开启**（`sampling_bias_correction='auto'`）
5. **时间信号未过 DRT 时不输出定年结论**（RSV 209 实测 R²=0.02–0.04，
   DRT 不显著时速率/tMRCA 只是把噪声拟合出的斜率）

## A3 为什么要**单独再跑一遍** clock

A1 与 A3 是**两次独立的 TreeTime 调用**，不是同一次。原因：`--coalescent skyline`
会改变速率估计（RSV 209 实测：不带 3.07e-4，带 2.45e-4），把两个模型的数字
混在一行里就没法说清"这个速率是哪个模型的"。所以 `clock_stats.tsv` 只记 A1，
A3 的天际线连同它自己那次的速率单独落在 `coalescent/` 下。
"""

import csv
import datetime
import hashlib
import io
import json
import os
import re
import time

from Virus_Platform_Core import treetime_ml as tt
from Virus_Platform_Core.config import PLATFORM_ROOT
from Virus_Platform_Core.utils import check_path, open_write, safe_open

# 阶段顺序与显示名（前端阶段流按这个顺序渲染）
STAGE_ORDER = ('prep', 'clock', 'phylogeo', 'coalescent', 'ancestral', 'report')
STAGE_META = {
    'prep':       {'code': 'A0', 'name': '数据准备', 'en': 'Data prep'},
    'clock':      {'code': 'A1', 'name': '时间推断', 'en': 'Molecular clock'},
    'phylogeo':   {'code': 'A2', 'name': '地理推断', 'en': 'Phylogeography'},
    'coalescent': {'code': 'A3', 'name': '群体动态', 'en': 'Coalescent skyline'},
    'ancestral':  {'code': 'A4', 'name': '祖先序列', 'en': 'Ancestral sequences'},
    'report':     {'code': 'A5', 'name': '汇总溯源', 'en': 'Run report'},
}

# A1 的定年通路：**只留 TreeTime**（2026-09-17 定稿）。
# 原先并列的 LSD2（严格钟最小二乘）与 treedater（非相关松弛钟）已摘除：
# 三法速率互不可比、维护三套结果口径的成本高于收益，且下游
# （A2 地理 / A3 天际线 / A4 祖先序列）本来就只吃 TreeTime 的节点年代。
# 底层实现仍留在 phylogeo.py / treedater.py 里未接线，将来要恢复只需重新接线。
CLOCK_METHODS = ('treetime',)

_MONTH_DAYS = (31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)


def _say(log):
    return log or (lambda *_a, **_k: None)


def _prog(progress):
    return progress or (lambda *_a, **_k: None)


def _dump(obj, path, indent=1):
    """落盘 JSON：NaN/Inf → null（与 web/tool_jobs._dump_json 同口径）。

    ⚠️ 前端 `r.json()` 是浏览器 JSON.parse，裸 NaN 会直接抛错 →
    结果面板静默空白。这里自己再兜一层，模块被脚本直接调用时也不会写出非法 JSON。
    """
    def _clean(v):
        if isinstance(v, float):
            if v != v or v in (float('inf'), float('-inf')):
                return None
            return v
        if isinstance(v, dict):
            return {k: _clean(x) for k, x in v.items()}
        if isinstance(v, (list, tuple)):
            return [_clean(x) for x in v]
        return v
    with open_write(path) as f:
        json.dump(_clean(obj), f, ensure_ascii=False, indent=indent,
                  allow_nan=False)
    return path


def _sha256(path, limit=None):
    h = hashlib.sha256()
    try:
        with open(path, 'rb') as f:
            while True:
                b = f.read(1 << 20)
                if not b:
                    break
                h.update(b)
    except OSError:
        return None
    return h.hexdigest()[:16]


def _fnum(v, nd=6):
    """数值 → 稳定字符串（None/NaN → 空串；避免 1e-05 与 1E-5 两写法）。"""
    if v is None:
        return ''
    try:
        x = float(v)
    except (TypeError, ValueError):
        return str(v)
    if x != x:
        return ''
    return f'{x:.{nd}g}'


# ---------------------------------------------------------------- 日期归一

_ISO_RE = re.compile(r'^(\d{4})(?:-(\d{1,2}))?(?:-(\d{1,2}))?$')
_RANGE_RE = re.compile(r'^\[\s*([^:]+?)\s*:\s*([^]]+?)\s*\]$')


def to_decimal_year(v, mode='mid', stats=None):
    """日期 → 十进制年。认 `2013` / `2013.42` / `2013-05` / `2013-05-01` /
    `[2013.2:2013.7]`（区间取中点）。

    返回 `(decimal_year 或 None, 是否精确到日)`。**认不出就返回 None**，
    绝不猜 —— 猜错会让整条链的年份轴错位。

    `stats`（可选 dict）：给了就累加本函数遇到的**二义日期**计数
    （`stats['n_ambiguous_frac']`，形如 `2019.5` 这种"小数年 vs 年.月"写法）；
    不给就不统计（纯转换调用零副作用）。

    ⚠️ 小数年快路径（2026-10-03 修，与 `phylodyn_govern.parse_date` 对齐）
    ----------------------------------------------------------------
    原实现 `if '.' in s: return float(s)` 带**前导零**的两位小数会走错口径：
    `2014.11` 被当小数年 2014.11（≈2 月），而唯一权威 `parse_date` 判它是
    「年.月」= 2014-11（≈11 月中）—— 同一份数据两个答案，差 **0.76 年**。
    现在快路径**只放行小数位 ≥3** 的写法（`2014.056`，explorer 导出的
    小数年都是这个形状）；其余（`.01`–`.12` 的「年.月」、`.5` 这种二义短小数）
    一律交回 `parse_date`，口径由它定，二义写法经 `stats` 计数上报。

    ⚠️ 口径（2026-09-18 修 P0，唯一实现见 `phylodyn_govern.parse_date`）
    ----------------------------------------------------------------
    原来日精度分支写的是 `y + (d - 0.5) / dim` —— **月份项整个丢了**：
    `2013-01-15`、`2013-03-15`、`2013-05-15`、`2013-08-15`、`2013-12-15`
    **全都算出 2013.4677**（只按"当月第几天"定位）。对 A0 全链是致命的：
    `stage_prep` 把 `dates[n]` 直接交给 A1 定年 / A2 的 MOTP 时间分箱，
    而用户的数据里**日精度有 2,651 行** → 整条时间轴在年内被压扁、不同月份的样本撞车。

    现在统一到**服务器 / `beast_handoff._decimal_of` 的同一条口径**：:

        yr + (mo - 1 + (dy - 1) / days_in_month) / 12        # mode='mid'
        yr + (mo - 1) / 12                                    # mode='start'

    · 只给年月 → `mode='mid'` 取该月 15 日、`'start'` 取月初
    · 只给年   → `mode='mid'` 取 **6 月 15 日**（服务器/上游口径）、`'start'` 取 1 月 1 日
      ⚠️ 这一档与修复前**不同**（修复前是"年初"）：用户数据里只给年的有 210 行，
      会整体后移约 0.46 年。之所以跟服务器对齐，是因为**两条路径必须给出同一个答案** ——
      `explorer/virphykit_export.py` 早就写明"服务器 `to_decimal_year` 把只给年的
      当 6 月 15 日"，本地再按年初算就是同一份数据两个答案（上游建
      `utils/decimal_year.py` 要治的正是这个病）。
    """
    if v is None:
        return None, False
    s = str(v).strip().strip('"')
    if not s:
        return None, False
    m = _RANGE_RE.match(s)
    if m:
        # 区间写法 `[2013.2:2013.7]` 本身就声明了"小数年" → 端点**不**计二义
        # （不把 stats 往下传，否则 `[2013.2:…]` 的两个端点会被当成二义写法）
        a, _ = to_decimal_year(m.group(1), mode=mode)
        b, _ = to_decimal_year(m.group(2), mode=mode)
        if a is not None and b is not None:
            return (a + b) / 2.0, False
        return (a if a is not None else b), False
    # ⚠️ 只有**带小数点**的才算"已经算好的小数年"（`2013.42`）。
    #    光秃秃的 `2013` **不能**走这条路：服务器口径把它当"只给年"→ 6 月 15 日
    #    （2013.4556），而 `float('2013')=2013.0` 是"年初"，差 0.46 年。
    #    两条路径（本地 A0 vs 服务器 BEAST）必须给同一个答案，所以这里收紧。
    #
    #    ⚠️ 2026-10-03：快路径再收紧到**小数位 ≥3**（`2014.056`）。原来只要含
    #    '.' 就直通 float，于是 `2014.11` 被读成 2 月（2014.11），而
    #    `phylodyn_govern.parse_date` 判它是「年.月」= 11 月中（2014.8722）——
    #    同一份数据两个答案，差 **0.76 年**。`.01`–`.12` 这种"年.月"写法与
    #    `.5` 这种二义短小数现在一律交回 parse_date（唯一权威实现）定夺；
    #    小数位 ≥3 的（explorer 导出的小数年形状）才走快路径。
    _frac = s.rsplit('.', 1)
    if '.' in s and len(_frac) == 2 and _frac[0].isdigit() and _frac[1].isdigit() \
            and len(_frac[1]) >= 3:
        try:
            return float(s), False
        except ValueError:
            pass
    if not _ISO_RE.match(s):
        # 服务器口径接受 `YYYY/MM/DD`（`decimal_year.to_decimal_year` 会先 replace）→ 跟着收；
        # 小数年形状（`2014.5` / `2014.11`）同样放行去 parse_date（口径唯一权威）
        s = s.replace('/', '-')
        if not _ISO_RE.match(s) and not re.fullmatch(r'\d{4}\.\d+', s):
            return None, False
    # 交给唯一实现（stdlib-only，无循环导入；懒加载避免包初始化顺序问题）
    from Virus_Platform_Core.phylodyn_govern import parse_date
    info = parse_date(s, mode=mode)
    if info.decimal is None:
        return None, False
    if stats is not None and info.ambiguous_frac:
        # 二义写法（`2019.5`：既是 5 月也可能是小数年 0.5）按小数年处理，
        # 但**必须计数**：调用方要把这个数报出去，否则"选了哪一个"没人知道。
        stats['n_ambiguous_frac'] = stats.get('n_ambiguous_frac', 0) + 1
    return info.decimal, (info.precision == 'day')


def _read_fasta_names(path, limit=0):
    """读 FASTA 头首词（叶名口径与 `phylogeo.safe_first_token` 一致）。"""
    names = []
    with safe_open(path) as f:
        for ln in f:
            if isinstance(ln, bytes):
                ln = ln.decode('utf-8', 'replace')
            if ln.startswith('>'):
                tok = ln[1:].strip().split()
                names.append(tok[0] if tok else '')
                if limit and len(names) >= limit:
                    break
    return names


def _write_fasta(records, path, width=70):
    with open_write(path) as f:
        for nm, seq in records:
            f.write('>' + nm + '\n')
            for i in range(0, len(seq), width):
                f.write(seq[i:i + width] + '\n')
    return path


def _read_fasta_records(path):
    recs = []
    name, buf = None, []
    with safe_open(path) as f:
        for ln in f:
            if isinstance(ln, bytes):
                ln = ln.decode('utf-8', 'replace')
            s = ln.rstrip('\r\n')
            if s.startswith('>'):
                if name is not None:
                    recs.append((name, ''.join(buf)))
                tok = s[1:].strip().split()
                name = tok[0] if tok else ''
                buf = []
            elif name is not None:
                buf.append(s.strip())
    if name is not None:
        recs.append((name, ''.join(buf)))
    return recs


def _read_table(path):
    """读 CSV/TSV → (header, rows)，自动判分隔符，剥 BOM 与 `#`。"""
    with io.open(path, encoding='utf-8-sig', errors='replace', newline='') as f:
        sample = f.read(8192)
        f.seek(0)
        delim = '\t' if sample.count('\t') > sample.count(',') else ','
        rdr = csv.reader(f, delimiter=delim)
        rows = [r for r in rdr]
    if not rows:
        return [], []
    header = [h.strip().lstrip('#').strip() for h in rows[0]]
    return header, [r for r in rows[1:] if any((c or '').strip() for c in r)]


# ================================================================ A0 prep

# ============================================ A0 冗余剪枝（对齐上游 --prune）
# 上游参照：`phymap-workflow/scripts/prep_metadata.py::prune_sequences`
# （`--prune` 开关；默认关闭）。这里把它的**两种方法、阈值口径、代表选择规则**
# 原样搬过来，只在两处做了平台化：
#   ① 区划只有一个（平台 States 是一列），不像上游能把多列拼成一个 location 串；
#   ② 距离矩阵优先 numpy 向量化（上游 pairsnp 可选 / numpy 兜底），
#      两条路都做，报告里写明实际走了哪条。

#: 参与比较的碱基（只在**双方都是 ACGT** 的位上计数 —— N/简并一律不比）
_SNP_BASES = frozenset('ACGT')

#: 剪枝规模上限：超过就**明说跳过**（连矩阵都是 n²，硬跑只会拖垮整条链）
_PRUNE_MAX_N = 6000

#: clade（完全连锁聚类）的规模上限 —— 它最坏是 O(n³)，比 fps 贵得多
_PRUNE_CLADE_MAX_N = 1500

PRUNE_METHODS = ('fps', 'clade')
PRUNE_RESOLUTIONS = ('year', 'month', 'week', 'day')


def _date_group(raw, resolution='month'):
    """原始日期串 → 归组键（`2013` / `2013-05` / `2013-W19` / `2013-05-01`）。

    与上游 `get_date_group` 同口径，**包括那两个看起来别扭的默认值**：
    只写年份的日期补 `07` 月、`15` 日。这是上游的既定行为，它决定"年内无月份
    信息"的样品归到哪一组 —— 不改成"按年归组"是为了让两边结果可对照。

    ⚠️ **一处刻意偏离**（已在 `_check_t9_prune.py` 里断言并记录）：
    `YYYY[-MM[-DD]]` 之外的形态（最典型的是**十进制年** `2008.58197`）上游一律
    返回 `unknown` —— 结果是**所有样品挤进同一个日期组**，`fps` 会跨年份剪枝，
    时间分辨率被剪没。这里改成按**年份**归组（`2008`）。对 ISO 形态两边逐字相同。
    """
    if raw is None:
        return 'unknown'
    val = str(raw).strip()
    if not val or val.lower() in ('nan', 'null', 'none'):
        return 'unknown'
    val = val.replace('/', '-')
    parts = val.split('-')
    if parts and len(parts[0]) == 4 and parts[0].isdigit():
        year = parts[0]
        month = parts[1] if len(parts) >= 2 and parts[1] else '07'
        day = parts[2] if len(parts) >= 3 and parts[2] else '15'
        if resolution == 'year':
            return year
        if resolution == 'day':
            return '%s-%s-%s' % (year, str(month).zfill(2), str(day).zfill(2))
        if resolution == 'week':
            try:
                dt = datetime.date(int(year), int(month), int(day))
                return '%s-W%02d' % (year, dt.isocalendar()[1])
            except (ValueError, TypeError):
                return '%s-%s' % (year, str(month).zfill(2))
        return '%s-%s' % (year, str(month).zfill(2))
    # 非 ISO 形态（如十进制年 `2008.58197`）：上游给 'unknown'（＝不按时间分组），
    # 这里给**年份**。这是刻意的偏离，理由见 docstring。
    try:
        return str(int(float(val)))
    except (TypeError, ValueError):
        return str(val)


def _pair_snp(a, b):
    """两条（已 upper 的）序列的 SNP 数。"""
    n = min(len(a), len(b))
    return sum(1 for i in range(n)
               if a[i] != b[i] and a[i] in _SNP_BASES and b[i] in _SNP_BASES)


class _SnpDist(object):
    """成对 SNP 距离。numpy 可用 → 一次性向量化；否则逐对现算。

    `how` 暴露给报告 —— "这次用的是哪条路"必须能查到（两条路数值口径一致，
    但速度差一个量级，排障时要能分辨）。
    """

    def __init__(self, seqs):
        up = [str(s).upper() for s in seqs]
        self._seqs = up
        self.n = len(up)
        self._m = None
        self.how = '纯 Python 逐对（numpy 不可用）'
        try:
            import numpy as np
        except ImportError:
            np = None
        if np is None or not self.n:
            return
        L = max((len(s) for s in up), default=0)
        if not L:
            return
        enc = np.full((self.n, L), 45, dtype=np.uint8)      # 45 = '-'
        for i, s in enumerate(up):
            b = s.encode('ascii', 'replace')
            enc[i, :len(b)] = np.frombuffer(b, dtype=np.uint8)
        ok = (enc == 65) | (enc == 67) | (enc == 71) | (enc == 84)
        m = np.zeros((self.n, self.n), dtype=np.int32)
        for i in range(self.n):
            m[i] = ((enc != enc[i]) & ok & ok[i]).sum(axis=1)
        self._m = m
        self.how = 'numpy 向量化'

    def d(self, i, j):
        if self._m is not None:
            return int(self._m[i, j])
        return _pair_snp(self._seqs[i], self._seqs[j])

    @property
    def matrix(self):
        return self._m


def prune_sequences(names, seq_of, state_of, raw_date_of, method='fps',
                    resolution='month', max_reps=3, min_snp_diff=2,
                    clade_cutoff=5, say=None):
    """按 (区划, 日期组[, 支系]) 去冗余 → `(kept_ids, report)`。

    两种方法（与上游 `--prune` 同名同义）：
      - **fps**（默认）：按 (区划, 日期组) 分组。组内条数 ≤ `max_reps` 直接全留；
        否则从"N 最少"那条起步做**最远点采样**，每轮挑"到已选集合的最小距离"最大的
        候选，直到凑满 `max_reps` 条或该最小距离 < `min_snp_diff`。
      - **clade**：先算**全量**成对 SNP 距离做**完全连锁聚类**（每轮在所有簇对里
        取"两簇成员两两最大距离"最小的那一对合并，直到该最小值 > `clade_cutoff`），
        再按 (区划, 日期组, 支系) 每组留 1 条。合并顺序与并列取舍**逐字复刻上游**
        （含 `pop(j)` 造成的簇移位语义）—— 并列距离在实际数据里很常见（SNP 距离
        是整数），随手改成别的顺序就会得到另一组代表。

    代表一律取 **N 最少**的那条（上游口径）—— 按文件顺序或名字排序挑代表会让
    "剔谁"不可复现，而剪枝结果直接影响速率/tMRCA，必须可复现。

    返回的 `report` 逐条记录：走了哪个方法、阈值、每个分组留了几条、
    **被剔掉的是谁**（前 200 条随报告落盘，全量另存 `prune_dropped.tsv`）。
    """
    say = say or (lambda *_a, **_k: None)
    rep = {'enabled': True, 'method': method, 'resolution': resolution,
           'max_reps': max_reps, 'min_snp_diff': min_snp_diff,
           'clade_cutoff': clade_cutoff, 'n_in': len(names),
           'n_out': None, 'n_dropped': None, 'kept': [], 'dropped': [],
           'groups': [], 'how': None, 'warnings': []}

    ids = [n for n in names if n in seq_of]
    missing = [n for n in names if n not in seq_of]
    if missing:
        rep['warnings'].append(
            '%d 条序列在 FASTA 里找不到序列本体（无法算距离）→ 一律**保留**，'
            '不参与剪枝：%s' % (len(missing), '、'.join(missing[:10])))
        rep['kept'].extend(missing)
    if len(ids) < 3:
        rep['warnings'].append('可用于剪枝的序列不足 3 条 → 不剪')
        rep['kept'].extend(ids)
        rep['kept'] = sorted(set(rep['kept']))
        rep['n_out'] = len(rep['kept'])
        rep['n_dropped'] = len(names) - rep['n_out']
        return rep['kept'], rep

    def n_count(nid):
        return str(seq_of[nid]).upper().count('N')

    # ⚠️ 并列时的先后规则：**按输入 FASTA 的顺序**先到先得。
    # 上游用 DataFrame 顺序（也＝输入顺序），所以这样两边结果一致；而按名字
    # 排序看齐"确定性"会把结果改成另一组代表 —— 剪枝直接影响速率/tMRCA，
    # 不能"看着等价"就行。（同输入 → 同输出，一样可复现。）
    order = {n: i for i, n in enumerate(ids)}

    def pick(members):
        return min(members, key=lambda x: (n_count(x), order[x]))

    def loc(nid):
        return (state_of.get(nid) or '未知区划')

    def dgrp(nid):
        return _date_group(raw_date_of.get(nid), resolution)

    kept = set(rep['kept'])
    dist = _SnpDist([seq_of[n] for n in ids])
    idx = {n: i for i, n in enumerate(ids)}
    rep['how'] = dist.how
    say('[prune] %d 条序列，SNP 距离：%s' % (len(ids), dist.how))

    if method == 'fps':
        groups = {}
        for nid in ids:
            groups.setdefault((loc(nid), dgrp(nid)), []).append(nid)
        for (lo, dg), members in sorted(groups.items()):
            members = sorted(members, key=lambda x: order[x])   # 回到输入顺序
            if len(members) <= max_reps:
                kept.update(members)
                rep['groups'].append(
                    {'loc': lo, 'date_group': dg, 'n_in': len(members),
                     'n_kept': len(members), 'kept': members,
                     'why': '组内 %d ≤ max_reps=%d，整组保留'
                            % (len(members), max_reps)})
                continue
            # ---- 最远点采样：起点＝N 最少（并列按输入顺序先到先得）----
            start = pick(members)
            sel = [start]
            while len(sel) < max_reps:
                best, best_d = None, -1
                for cand in members:            # 输入顺序（与上游同序）
                    if cand in sel:
                        continue
                    md = min(dist.d(idx[cand], idx[s]) for s in sel)
                    if md > best_d:
                        best_d, best = md, cand
                if best is None or best_d < min_snp_diff:
                    break
                sel.append(best)
            sel = sorted(sel, key=lambda x: order[x])
            kept.update(sel)
            rep['groups'].append(
                {'loc': lo, 'date_group': dg, 'n_in': len(members),
                 'n_kept': len(sel), 'kept': sel,
                 'why': 'FPS：起点取 N 最少的 %s；此后每轮挑「到已选集合最小距离」'
                        '最大的候选，最小距离 < min_snp_diff=%d 即停（实收 %d 条）'
                        % (start, min_snp_diff, len(sel))})

    elif method == 'clade':
        # 完全连锁聚类：按距离**升序**处理候选对，只在"合并后簇内最大跨距
        # ≤ clade_cutoff"时才合并 —— 这与上游"每轮取簇间最小距离（完全连锁
        # 距离）再合并"是同一个划分，但少了 O(簇数²) 的重复扫描。
        # ⚠️ **逐字复刻上游的贪心**：每轮在**所有簇对**里取「完全连锁距离」
        # （两簇成员两两最大距离）**最小**的一对合并；并列时取簇下标最小的一对
        # （上游 `for i: for j in range(i+1,…)` + 严格 `<` 就是这个语义）；
        # 合并用 `clusters[i].extend(clusters[j]); clusters.pop(j)`——
        # **pop(j) 会让簇顺序整体移位**，而移位会改变下一轮的并列取舍，
        # 所以这里连 pop 的写法都保持一致。
        # 唯一没照抄的是内层的"逐成员求最大距离"：换成矩阵取块最大值
        # （同样的数，快一个量级），并做了版本号失效的缓存。
        cl = {}
        ver = {}
        members = {i: [i] for i in range(len(ids))}
        cl_order = list(range(len(ids)))       # 当前簇的顺序（pop 会让它移位）

        def _cl(ca, cb):
            """两簇的完全连锁距离（缓存；合并后靠版本号失效）。"""
            ent = cl.get((ca, cb))
            if ent is not None and ent[0] == ver.get(ca) and ent[1] == ver.get(cb):
                return ent[2]
            m = dist.matrix
            v = 0
            va, vb = members[ca], members[cb]
            if m is not None:
                import numpy as _np
                v = int(m[_np.ix_(va, vb)].max())
            else:
                for x in va:
                    for y in vb:
                        dd = dist.d(x, y)
                        if dd > v:
                            v = dd
            cl[(ca, cb)] = (ver.get(ca), ver.get(cb), v)
            return v

        n_cl = len(ids)
        skip_clade = n_cl > _PRUNE_CLADE_MAX_N
        if skip_clade:
            rep['warnings'].append(
                '序列数 %d 超过 clade 方法上限 %d（完全连锁聚类最坏 O(n³)）→ '
                '**本次不做剪枝**（保留全部样本）。请改用 fps 方法，'
                '或先用更粗的日期粒度。' % (n_cl, _PRUNE_CLADE_MAX_N))
            kept.update(ids)
        else:
            while cl_order:
                min_dist, to_merge = float('inf'), None
                for i in range(len(cl_order)):
                    for j in range(i + 1, len(cl_order)):
                        md = _cl(cl_order[i], cl_order[j])
                        if md < min_dist:
                            min_dist, to_merge = md, (i, j)
                if min_dist <= clade_cutoff and to_merge is not None:
                    i, j = to_merge
                    ca, cb = cl_order[i], cl_order[j]
                    members[ca].extend(members[cb])
                    ver[ca] = ver.get(ca, 0) + 1
                    cl_order.pop(j)           # 同上：与上游同写法（含移位语义）
                else:
                    break
        clade_of = {}
        for cid in cl_order:
            for pos in members[cid]:
                clade_of[ids[pos]] = 'C%05d' % cid
        rep['n_clades'] = len(cl_order) if cl_order else None

        if not skip_clade:
            groups = {}
            for nid in ids:
                groups.setdefault((loc(nid), dgrp(nid), clade_of[nid]),
                                  []).append(nid)
            for (lo, dg, cid), members_g in sorted(groups.items()):
                # 每个 (区划·日期组·支系) 留 1 条代表：N 最少优先，并列按输入顺序
                best = pick(members_g)
                kept.update([best])
                rep['groups'].append(
                    {'loc': lo, 'date_group': dg, 'clade': cid,
                     'n_in': len(members_g), 'n_kept': 1, 'kept': [best],
                     'why': 'clade：支系 %s 内按（区划·日期组）合并，'
                            '取 N 最少的代表' % cid})
    else:
        rep['warnings'].append('未知剪枝方法 %r → 不剪' % method)
        kept.update(ids)

    rep['kept'] = sorted(kept)
    rep['dropped'] = sorted(n for n in names if n not in kept)
    rep['n_out'] = len(rep['kept'])
    rep['n_dropped'] = len(rep['dropped'])
    return rep['kept'], rep


def apply_rename_map(aln, meta, tree, map_path, out_dir,
                     name_col='name', log=None):
    """可选前置：按「旧ID→新ID」映射批量重命名 FASTA / 元数据键 / 输入树 tip。

    返回 (aln_new, meta_new, tree_new, stats)。未映射的 ID 保持原名；
    tree 为 None 或文件不存在时原样返回（A1 TreeTime 自己建树的场景不受影响）。
    映射格式与 🏷 序列 ID 重命名工具一致：每行 `旧ID <TAB/逗号> 新ID`。
    """
    import csv as _csv
    import re as _re
    from Bio import SeqIO as _SeqIO
    m = {}
    # errors='replace'：与同模块其余读入口径一致 —— 非 UTF-8 字节不让整条链
    # 因 UnicodeDecodeError 中止（坏字节在 ID 上会明显到用户能看出来）
    for line in open(map_path, encoding='utf-8-sig', errors='replace'):
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        for sep in ('\t', '⇥', ','):
            if sep in line:
                old, _, new = line.partition(sep)
                old, new = old.strip(), new.strip()
                if old and new:
                    m[old] = new
                break
    if not m:
        raise ValueError('重命名映射为空或格式不对（每行应为 旧ID<TAB>新ID）')
    _new_names = list(m.values())
    dup_new = sorted({v for v in _new_names if _new_names.count(v) > 1})
    if dup_new:
        raise ValueError('映射里出现重复的新 ID: %s' % dup_new[:5])

    os.makedirs(out_dir, exist_ok=True)

    # ① FASTA
    recs = list(_SeqIO.parse(aln, 'fasta'))
    n_hit = 0
    for r in recs:
        new_id = m.get(r.id, m.get(r.name, r.id))
        if new_id != r.id:
            n_hit += 1
        r.id, r.name, r.description = new_id, new_id, new_id
    aln_new = os.path.join(out_dir, 'prep.renamed.fasta')
    _SeqIO.write(recs, aln_new, 'fasta')

    # ② 元数据（键列：优先 name，否则首列）
    meta_new, n_meta = meta, 0
    if meta and os.path.exists(meta):
        with open(meta, encoding='utf-8-sig', errors='replace', newline='') as f:
            sample = f.readline()
        delim = '\t' if sample.count('\t') >= sample.count(',') else ','
        with open(meta, encoding='utf-8-sig', errors='replace', newline='') as f:
            rows = list(_csv.DictReader(f, delimiter=delim))
        if rows:
            key = next((c for c in rows[0] if c.strip().lower() == name_col),
                       list(rows[0].keys())[0])
            for r in rows:
                old = str(r.get(key, '') or '').strip()
                if old in m:
                    r[key] = m[old]
                    n_meta += 1
            meta_new = os.path.join(out_dir, 'meta.renamed' +
                                    ('.tsv' if delim == '\t' else '.csv'))
            with open(meta_new, 'w', encoding='utf-8', newline='') as f:
                w = _csv.DictWriter(f, fieldnames=list(rows[0].keys()),
                                    delimiter=delim)
                w.writeheader()
                w.writerows(rows)

    # ③ 输入树 tip（文本级替换：newick 里 tip 名只出现在 ( , ) 定界之间）
    tree_new, n_tree = tree, 0
    if tree and os.path.exists(tree):
        txt = open(tree, encoding='utf-8', errors='replace').read()
        for old, new in m.items():
            pat = _re.compile(r'(?<=[(,])' + _re.escape(old) + r'(?=[,:)])')
            txt, k = pat.subn(new, txt)
            n_tree += k
        tree_new = os.path.join(out_dir, 'tree.renamed.nwk')
        open(tree_new, 'w', encoding='utf-8').write(txt)

    stats = {'n_map': len(m), 'n_fasta_renamed': n_hit,
             'n_meta_renamed': n_meta, 'n_tree_tips_renamed': n_tree}
    if log:
        log('重命名: 映射 %d 条, fasta 命中 %d, 元数据改键 %d, 树 tip 改名 %d'
            % (stats['n_map'], n_hit, n_meta, n_tree))
    return aln_new, meta_new, tree_new, stats


def stage_prep(aln, meta, tree, out_dir, date_col=None, trait='region',
               date_trait=None, log=None, progress=None, build_tree=True,
               prune=False, prune_method='fps', prune_date_resolution='month',
               prune_max_reps=3, prune_min_snp_diff=2, prune_clade_cutoff=5):
    """A0 数据准备：FASTA × 元数据 → 日期表 / 状态表 / 规范化比对 / 树。

    ⚠️ 匹配门槛的**真实口径**（2026-10-03 校正：原 docstring 写"匹配率必须 100%
    才算通过"，与实现不符 —— 实际是半数放行，照它读会误判链路的严格程度）：
      · 缺日期：逐条**告警**（A1 定年时缺日期的叶会被剔除，树变小；**一条日期
        都解析不出**才 error 中止）；
      · 缺区划：告警；**超过半数**拿不到区划才 error（多半是列名/表头格式不对）；
    未匹配名单原样列出（红线 7）—— TreeTime 那边名字对不上是**静默算错**，
    不报就是埋雷。完成日志同时报出实际缺失比例（缺日期 x/y、缺区划 x/y）。

    `prune=True` 时先做**冗余剪枝**（对齐上游 phymap-workflow `--prune`）：
    按 (区划, 日期组[, SNP 支系]) 去冗余，默认关闭。**有损**，详见
    `prune_sequences`；剪枝记录落 `prune_report.json` + `prune_dropped.tsv`。

    返回 dict（`ok` / `dates` / `states` / `n_*` / `unmatched_*` / `outputs`）。
    """
    say, prog = _say(log), _prog(progress)
    os.makedirs(out_dir, exist_ok=True)
    say('A0 数据准备')
    res = {'stage': 'prep', 'ok': False, 'error': None, 'out_dir': out_dir,
           'aln': aln, 'meta': meta, 'tree_in': tree, 'date_col': date_col,
           'trait': trait, 'warnings': [], 'outputs': {}, 'notes': []}

    names = _read_fasta_names(aln)
    if len(names) < 3:
        res['error'] = f'比对里只有 {len(names)} 条序列 —— 时间/地理推断至少要 3 条'
        return res
    res['n_seqs'] = len(names)
    res['n_unique_names'] = len(set(names))
    if len(set(names)) != len(names):
        dup = sorted({n for n in names if names.count(n) > 1})[:10]
        res['warnings'].append(
            f'比对里有重名序列（{len(names) - len(set(names))} 条）—— '
            f'元数据匹配与走廊计数都会错位。样例：{"、".join(dup)}')
    say(f'比对：{len(names)} 条序列')

    # ---- 元数据（可选）----
    hdr, rows = ([], [])
    if meta:
        hdr, rows = _read_table(meta)
        if not hdr:
            res['warnings'].append('元数据表是空的 → 退回从 FASTA 头解析区划/年份')
    meta_map = {}
    if hdr:
        key_i = 0
        for i, h in enumerate(hdr):
            if h.lower() in ('seq_id', 'seqid', 'accession', 'sample', 'name'):
                key_i = i
                break
        for r in rows:
            if key_i < len(r) and r[key_i].strip():
                meta_map[r[key_i].strip()] = {
                    h: (r[i].strip() if i < len(r) else '')
                    for i, h in enumerate(hdr)}
    # ---- 列名归一化 + 自动择优（2026-09-18）----
    # explorer 新版导出的 25 列（Accession / Collection_Date / Release_Date /
    # Geo_Location / Country / Year …）在这里被翻译成规范列 name/date/location：
    # 日期采样类优先（Collection_Date，空则 Release_Date 并**会在 summary 里告警
    # 时间轴偏晚**），地点详细类优先（Geo_Location 'Russia: Stavropol' > Country）。
    _pick_notes = {}
    if meta_map:
        from Virus_Platform_Core.phylodyn_kit import normalize_meta_rows as _norm
        _keys = list(meta_map)
        _rows2, _pick_notes = _norm([meta_map[k] for k in _keys])
        for _k, _r2 in zip(_keys, _rows2):
            _r2.setdefault('name', _k)
            meta_map[_k] = _r2
        for _canon in ('name', 'date', 'location'):
            if _canon not in [h.strip().lower() for h in hdr]:
                hdr = list(hdr) + [_canon]
        res['meta_pick'] = _pick_notes
        for _w in _pick_notes.get('warnings') or []:
            res['warnings'].append(_w)
    res['meta_columns'] = hdr

    # ---- 日期列 / 区划列自动判定 ----
    if hdr and not date_col:
        low = {h.strip().lower(): h for h in hdr}
        for want in ('date', 'collection_date', 'sampling_date', 'release_date',
                     'decimal_date', 'time', '日期'):
            if want in low:
                date_col = low[want]
                break
    if hdr and not date_trait:
        low = {h.strip().lower(): h for h in hdr}
        for want in ('year', 'decimal_year', 'sampling_year', '年份'):
            if want in low:
                date_trait = low[want]
                break
    if hdr and trait not in hdr:
        low3 = {h.strip().lower(): h for h in hdr}
        for want in ('location', 'geo_loc_name', 'geo_location', 'country',
                     'region', 'state', 'province', '地点', '区域'):
            if want in low3:
                res['warnings'].append(
                    f'区划列自动识别为「{low3[want]}」'
                    + (f'（卡上给的是「{trait}」，表里没有）' if trait else '（未指定列名）'))
                trait = low3[want]
                break
    res['trait'] = trait          # 覆盖 res 初始化时写入的**卡上传入值**
    res['trait_used'] = trait
    res['date_col'] = date_col
    res['date_trait'] = date_trait

    from Virus_Platform_Core import phylogeo as pg

    # ---- 逐序列取日期与区划 ----
    dates, states, no_date, no_state, date_src, state_src = {}, {}, [], [], {}, {}
    raw_dates = {}                                  # 原始日期串（剪枝按它归组）
    _dstats = {}                                    # 二义日期计数（见 to_decimal_year）
    for n in names:
        row = meta_map.get(n) or {}
        if not row and meta_map:                    # 前缀/包含兜底（同 phylogeo 口径）
            for k, v in meta_map.items():
                if k and (n.startswith(k) or k in n):
                    row = v
                    break
        dv = None
        if row:
            if date_col and date_col in row:
                dv, _p = to_decimal_year(row[date_col], stats=_dstats)
                if dv is not None:
                    date_src[n] = f'元数据列 {date_col}'
                    raw_dates[n] = row[date_col]
            if dv is None and date_trait and date_trait in row:
                dv, _p = to_decimal_year(row[date_trait], stats=_dstats)
                if dv is not None:
                    date_src[n] = f'元数据列 {date_trait}'
                    raw_dates[n] = row[date_trait]
        if dv is None:                              # 退回 FASTA 头（REGION_ACC_YEAR）
            tr = pg._header_traits(n)
            dv, _p = to_decimal_year(tr.get('year'), stats=_dstats)
            if dv is not None:
                date_src[n] = 'FASTA 头（名称_登录号_年）'
                raw_dates[n] = tr.get('year')
        if dv is not None:
            dates[n] = dv
        else:
            no_date.append(n)

        sv = ''
        if row:
            sv = (row.get(trait) or '').strip()
            if sv:
                state_src[n] = f'元数据列 {trait}'
        if not sv:
            tr = pg._header_traits(n)
            sv = (tr.get('region') or '').strip()
            if sv:
                state_src[n] = 'FASTA 头第 2 段'
        if sv:
            states[n] = sv
        else:
            no_state.append(n)

    res['n_seqs_input'] = len(names)

    # ---- 冗余剪枝（可选；对齐上游 phymap-workflow `--prune`）----
    # ⚠️ 这是**有损**操作：它改变样本组成，进而改变速率 / tMRCA / 走廊计数。
    # 所以：① 默认关闭；② 逐条记录剔了谁、为什么；③ 剪完的样本数写进报告，
    # 引用结果时必须连"剪过"一起说。
    # 位置在「日期/区划解析之后、落盘之前」—— 分组要靠区划与日期组，
    # 而 dates.csv / states.csv / prep.fasta / 树都必须按留下的样本重写。
    keep_for_write = None               # 非 None 时落盘只写这些样本（剪枝后）
    drop_ctx = None                     # 剪枝前的 区划/原始日期 副本（写剔除清单用）
    if prune and len(names) > 3:
        prog('prep', 0.12, 'A0 冗余剪枝')
        recs_all = _read_fasta_records(aln)
        seq_of = {}
        for nm, sq in recs_all:
            if nm and nm not in seq_of:
                seq_of[nm] = sq
        if len(names) > _PRUNE_MAX_N:
            res['prune'] = {
                'enabled': True, 'skipped': True, 'n_in': len(names),
                'reason': '序列数 %d 超过剪枝上限 %d —— 成对距离矩阵是 n²，'
                          '这个量级硬跑会拖垮整条链' % (len(names), _PRUNE_MAX_N)}
            res['warnings'].append('**剪枝已跳过**：' + res['prune']['reason'])
            say('[WARN] 剪枝已跳过：' + res['prune']['reason'])
        else:
            kept_ids, prep_prune = prune_sequences(
                names, seq_of, states, raw_dates, method=prune_method,
                resolution=prune_date_resolution, max_reps=prune_max_reps,
                min_snp_diff=prune_min_snp_diff,
                clade_cutoff=prune_clade_cutoff, say=say)
            keep_set = set(kept_ids)
            res['prune'] = prep_prune
            for w in (prep_prune.get('warnings') or []):
                res['warnings'].append(w)
            if len(keep_set) < 3:
                res['prune']['abandoned'] = True
                msg = ('剪枝后只剩 %d 条（< 3）→ **放弃剪枝、保留全部样本**；'
                       '阈值（max_reps=%s / min_snp_diff=%s / clade_cutoff=%s）'
                       '对这份数据太激进' % (len(keep_set), prune_max_reps,
                                            prune_min_snp_diff,
                                            prune_clade_cutoff))
                res['warnings'].append(msg)
                say('[WARN] ' + msg)
            else:
                names = [n for n in names if n in keep_set]
                # 先留一份剪枝前的 区划/原始日期：`prune_dropped.tsv` 要给**被剔掉的**
                # 样本也写上它原来的区划与日期组，否则"为什么剔它"就没法核对。
                drop_ctx = {'states': dict(states), 'raw': dict(raw_dates)}
                dates = {k: v for k, v in dates.items() if k in keep_set}
                states = {k: v for k, v in states.items() if k in keep_set}
                no_date = [n for n in no_date if n in keep_set]
                no_state = [n for n in no_state if n in keep_set]
                date_src = {k: v for k, v in date_src.items() if k in keep_set}
                state_src = {k: v for k, v in state_src.items() if k in keep_set}
                keep_for_write = keep_set
                res['n_seqs'] = len(names)
                note = ('A0 冗余剪枝已启用（方法 %s，日期粒度 %s）：'
                        '%d → %d 条（剔 %d）。**样本组成已被改变**，'
                        '本轮的速率 / tMRCA / 走廊计数都在剪枝后的样本上得出，'
                        '引用时必须一起说明。'
                        % (prune_method, prune_date_resolution,
                           prep_prune['n_in'], prep_prune['n_out'],
                           prep_prune['n_dropped']))
                res['notes'].append(note)
                say('[prune] ' + note)
                say('[prune] 距离计算：%s；分组 %d 个；被剔名单见 '
                    'prune_dropped.tsv' % (prep_prune.get('how'),
                                           len(prep_prune.get('groups') or [])))

    res.update({'n_dates': len(dates), 'n_states': len(states),
                'n_no_date': len(no_date), 'n_no_state': len(no_state),
                'unmatched_date': no_date[:50],
                'unmatched_state': no_state[:50],
                'date_source_mix': _tally(date_src.values()),
                'state_source_mix': _tally(state_src.values())})

    if _dstats.get('n_ambiguous_frac'):
        # 二义日期必须报数（to_decimal_year 的 stats 累加）：`2019.5` 按小数年
        # 处理是平台口径（数据源是 explorer 的 date 列），但"选了哪一个"要让用户看见。
        res['n_date_ambiguous_frac'] = _dstats['n_ambiguous_frac']
        res['warnings'].append(
            '%d 条序列的日期形如 `2019.5`（既是 5 月、也可能是小数年 0.5）—— '
            '本平台按**小数年**处理（数据源是 explorer 的 date 列）；'
            '若这批数据来自 BioAider 那种"年.月"写法，请改用它导出的 ISO'
            % _dstats['n_ambiguous_frac'])

    if no_date:
        res['warnings'].append(
            f'{len(no_date)} 条序列没有日期 —— TreeTime 定年要求**每个叶都有日期**，'
            f'缺的会被剔除后定年（树会变小）。未匹配名单（前 50）：'
            + '、'.join(no_date[:50]))
    if no_state:
        res['warnings'].append(
            f'{len(no_state)} 条序列没有区划 —— 地理重构会把这些叶当"缺数据"处理'
            f'（不产生迁移，但会摊薄走廊计数）。未匹配名单（前 50）：'
            + '、'.join(no_state[:50]))
    if no_state and len(no_state) > len(names) // 2:
        res['error'] = ('超过半数的序列拿不到区划 —— 多半是区划列名不对或 FASTA 头'
                        '格式不符（约定 `名称|区划|年份` 或 `区划_登录号_年`）。'
                        f'元数据实际列：{hdr or "（未提供元数据）"}；'
                        f'当前区划列名：{trait}')
        return res
    if not dates:
        res['error'] = '一条日期都没解析出来 —— 请提供含日期列的元数据表'
        return res

    # ---- 落盘 ----
    out = res['outputs']
    dates_csv = os.path.join(out_dir, 'dates.csv')
    tt.write_dates_csv(sorted(dates.items()), dates_csv)
    out['dates'] = dates_csv
    states_csv = os.path.join(out_dir, 'states.csv')
    tt.write_states_csv(sorted(states.items()), states_csv, trait)
    out['states'] = states_csv

    # 规范化比对（叶名统一成 FASTA 首词，TreeTime 与树才对得上）
    # 剪枝后只写留下来的序列 —— prep.fasta 是 A1/A2 的输入，
    # 写全量会让"剪了但没生效"，报告里的样本数与实际算的根本不是一回事。
    prep_fa = os.path.join(out_dir, 'prep.fasta')
    recs = _read_fasta_records(aln)
    if keep_for_write is not None:
        recs = [(nm, sq) for nm, sq in recs if nm in keep_for_write]
    _write_fasta([(nm, seq) for nm, seq in recs if nm], prep_fa)
    out['prep_fasta'] = prep_fa
    seqlen = max((len(s) for _n, s in recs), default=0)
    res['aln_len'] = seqlen
    res['seq_len_modal'] = _modal([len(s) for _n, s in recs])
    if res['seq_len_modal'] and seqlen != res['seq_len_modal']:
        res['warnings'].append(
            f'比对不是等长的（最长 {seqlen}、众数 {res["seq_len_modal"]}）—— '
            'TreeTime 会按最长的算，序列长度参数请填众数')

    # 剪枝记录落盘：剔了谁、为什么，必须能逐条查（红线：有损操作不许只留个数字）
    if res.get('prune'):
        pr = res['prune']
        dropped = pr.get('dropped') or []
        _st = (drop_ctx or {}).get('states') or {}
        _rw = (drop_ctx or {}).get('raw') or {}
        if dropped or not pr.get('skipped'):
            dp = os.path.join(out_dir, 'prune_dropped.tsv')
            with open_write(dp) as f:
                f.write('name\tstate\tdate_group\tkept\n')
                for n in sorted(pr.get('kept') or []):
                    f.write('%s\t%s\t%s\t1\n'
                            % (n, _st.get(n) or states.get(n) or '',
                               _date_group(_rw.get(n), prune_date_resolution)))
                for n in dropped:
                    f.write('%s\t%s\t%s\t0\n'
                            % (n, _st.get(n) or '', _date_group(
                                _rw.get(n), prune_date_resolution)))
            out['prune_dropped'] = dp
        rp = os.path.join(out_dir, 'prune_report.json')
        _dump(pr, rp)
        out['prune_report'] = rp

    # ---- 树：给了就用，没给就现场建 ----
    prep_nwk = os.path.join(out_dir, 'prep.nwk')
    _tree_ok = True
    if tree:
        try:
            with safe_open(tree) as f:
                txt = f.read()
            if isinstance(txt, bytes):
                txt = txt.decode('utf-8', 'replace')
            if txt.lstrip().upper().startswith('#NEXUS'):
                nwk = tt.nexus_tree_newick(txt)
            else:
                nwk = txt.strip()
            # 剪枝后：输入树也要**跟着剪叶**。不剪的话 prep.nwk 上留着已经被剔掉的
            # 叶，A1 会在"树上有 N 个叶在日期表里找不到"那儿再剔一次 —— 最终树一样，
            # 但用户看到的是"剪枝没生效"，而 prep.nwk 也不再是 A0 的真实输出。
            if keep_for_write is not None:
                try:
                    _rt0 = tt.parse_annotated_newick(nwk)
                    _pt, _pstat = tt.prune_tree(_rt0, keep_for_write)
                    if _pt is None or _pstat['n_kept_tips'] < 3:
                        _m = ('输入树按剪枝结果剪叶后只剩 %d 个叶（< 3）→ '
                              '**弃用输入树，改用剪枝后的 prep.fasta 现场建树**'
                              % _pstat['n_kept_tips'])
                        res['warnings'].append(_m)
                        say('[WARN] ' + _m)
                        _tree_ok = False
                    else:
                        nwk = tt.to_newick(_pt) + ';'
                        res['tree_pruned'] = {
                            'n_removed': _pstat['n_removed'],
                            'n_kept_tips': _pstat['n_kept_tips'],
                            'n_collapsed': _pstat['n_collapsed'],
                            'removed_tips': _pstat['removed_tips'][:200]}
                        say('[prune] 输入树剪叶：- %d 叶，剩 %d 叶'
                            % (_pstat['n_removed'], _pstat['n_kept_tips']))
                except ValueError as e:
                    res['warnings'].append(
                        f'输入树按剪枝结果剪叶失败（{e}）→ 保留原树（未剪叶）')
            if _tree_ok:
                with open_write(prep_nwk) as f:
                    f.write(nwk.rstrip().rstrip(';') + ';\n')
                out['prep_nwk'] = prep_nwk
                res['tree_source'] = '用户提供'
                try:
                    rt = tt.parse_annotated_newick(nwk)
                    res['tree_tips'] = len([x for x in tt.iter_pre(rt)
                                            if not x['children']])
                except ValueError as e:
                    res['warnings'].append(f'输入的树解析异常：{e}')
        except OSError as e:
            res['error'] = f'读不到输入树：{e}'
            return res
    if (not tree or not _tree_ok) and build_tree and not out.get('prep_nwk'):
        say('没有可用的树 → 用当前 prep.fasta 现场建（FastTree / NJ）')
        prog('prep', 0.35, 'FastTree 建树')
        try:
            from Virus_Platform_Core.phylo import _run_fasttree
            # ⚠️ 关键字是 `logger`，不是 `log` —— 传 `log=` 会 TypeError，
            # 而它被下面的 except 吞成"现场建树失败"，于是"A0 不给树就自动建树"
            # 这条路一直没通（2026-09-17 由 tests/_check_t9_prune.py 的 F 段查出）。
            _run_fasttree(prep_fa, prep_nwk, logger=None)
            out['prep_nwk'] = prep_nwk
            res['tree_source'] = 'FastTree（本次现场建）'
            rt = tt.parse_annotated_newick(open(prep_nwk, encoding='utf-8').read())
            res['tree_tips'] = len([x for x in tt.iter_pre(rt)
                                    if not x['children']])
        except Exception as e:
            res['error'] = (f'现场建树失败：{type(e).__name__}: {e} —— '
                            '请直接提供一棵 newick 树（A1 定年必须有树）')
            return res
    elif not out.get('prep_nwk'):
        # 只有"确实没有可用树"才告警 —— 用户给了树且建好了就别再喊"没给树"
        res['warnings'].append('没给树且关闭了现场建树 —— A1 定年无法进行')

    # ---- 树上的叶与日期/状态表的交集（这是真正决定成败的匹配率）----
    if out.get('prep_nwk'):
        try:
            rt = tt.parse_annotated_newick(
                open(out['prep_nwk'], encoding='utf-8').read())
            tips = [x['name'] for x in tt.iter_pre(rt) if not x['children']]
            res['tree_tips'] = len(tips)
            miss_d = [t for t in tips if t not in dates]
            miss_s = [t for t in tips if t not in states]
            res['tree_tips_no_date'] = len(miss_d)
            res['tree_tips_no_state'] = len(miss_s)
            res['tree_unmatched_date'] = miss_d[:50]
            res['tree_unmatched_state'] = miss_s[:50]
            if miss_d:
                res['warnings'].append(
                    f'**树上有 {len(miss_d)}/{len(tips)} 个叶在日期表里找不到** —— '
                    'TreeTime 定年会把这些叶剪掉（树变小、tMRCA 受影响）。'
                    f'名单（前 50）：{"、".join(miss_d[:50])}')
        except ValueError as e:
            res['warnings'].append(f'规范化后的树解析失败：{e}')

    res['date_range'] = ([round(min(dates.values()), 4), round(max(dates.values()), 4)]
                         if dates else None)
    if res['date_range']:
        span = res['date_range'][1] - res['date_range'][0]
        res['date_span'] = round(span, 4)
        if span < 3:
            res['warnings'].append(
                f'采样年跨度只有 {span:.2f} 年 —— 分子钟定年需要足够的时间跨度，'
                '跨度太小时速率与 tMRCA 都不可信')
    res['state_counts'] = _tally(states.values())
    res['ok'] = not res['error']
    prog('prep', 1.0, f'prep 完成：{len(dates)} 日期 / {len(states)} 区划')
    _n_seq = max(1, len(names))
    say(f'A0 完成：日期 {len(dates)}/{len(names)}（缺 {len(no_date)} 条 = '
        f'{100.0 * len(no_date) / _n_seq:.1f}%），'
        f'区划 {len(states)}/{len(names)}（缺 {len(no_state)} 条 = '
        f'{100.0 * len(no_state) / _n_seq:.1f}%），'
        f'年份跨度 {res.get("date_span")} 年，'
        f'树 {res.get("tree_tips")} 叶（{res.get("tree_source")}）')
    _dump(res, os.path.join(out_dir, 'prep_report.json'))
    return res


def _tally(vals):
    out = {}
    for v in vals:
        out[v] = out.get(v, 0) + 1
    return dict(sorted(out.items(), key=lambda x: -x[1]))


def _modal(vals):
    if not vals:
        return None
    t = _tally(vals)
    return int(max(t.items(), key=lambda x: (x[1], x[0]))[0])


# ================================================================ A1 clock

def stage_clock(prep, out_dir, methods=('treetime',), aln=None, tree=None,
                dates=None, seq_len=None, relax=None, reroot='least-squares',
                keep_root=False, confidence=True, time_marginal='only-final',
                clock_filter=None, clock_rate=None,
                drt_perm=20, drt_seed=42, do_drt=True,
                rng_seed=0, timeout=3600, log=None, progress=None):
    """A1 时间推断：用 **TreeTime** 定年，产出接口产物 `clock/timetree.nwk`。

    `methods`: 只接受 `('treetime',)`（= `CLOCK_METHODS`）。
    ⚠️ 2026-09-17 起 LSD2 与 treedater 两条通路已摘除，传进来会被过滤掉并告警
    （见 `CLOCK_METHODS` 的说明：三法速率不可比，维护成本高于收益）。
    下游 A2/A3/A4 本来就只吃 TreeTime 的节点年代，因此收敛后不影响链路。

    `clock_rate`: 固定分子钟速率（见 `treetime_ml.run_clock` 的口径说明）。
    给了就**关掉速率估计**，`clock_stats.tsv` 的 `model` 列会写成
    「严格钟（速率固定，用户指定）」——**不许**让它看起来像 ML 估出来的。
    这是慢演化病原（自动估不出钟）时的官方兜底。

    接口产物固定为 `clock/timetree.nwk`，由 TreeTime 产出。它同时给时间树 +
    节点年代 + 地理/天际线/祖先序列，是**唯一**的通路，因此不再有"退用别家
    时间树"的分支 —— 失败就是失败，报告里如实写，不回退到已摘除的实现。

    另外把该接口产物转成**前端可画的直角时间树坐标对象**放进 `res['timetree']`
    （见 `build_clock_timetree`），与 RTT / phylogeo 卡同一个形状，
    前端复用同一个 `drawTimeTree` —— 不再让 phylodyn 只给数字表不给图。
    """
    say, prog = _say(log), _prog(progress)
    methods = tuple(str(m).strip().lower() for m in (methods or ()))
    dropped = [m for m in methods if m not in CLOCK_METHODS]
    methods = tuple(m for m in methods if m in CLOCK_METHODS) or ('treetime',)
    os.makedirs(out_dir, exist_ok=True)
    res = {'stage': 'clock', 'ok': False, 'error': None, 'out_dir': out_dir,
           'methods': list(methods), 'warnings': [], 'notes': [],
           'outputs': {}, 'stats': [], 'timetree_source': None,
           'tt_env': None, 'r2': None, 'drt': None, 'timetree': None}
    if dropped:
        # 不静默：老配置/老脚本传 'lsd2'/'treedater' 时要看得见它被忽略了
        res['warnings'].append(
            'A1 已收敛到 TreeTime：忽略不支持的通路 ' + '、'.join(dropped)
            + '（LSD2 / treedater 已于 2026-09-17 摘除）')
        say(f"[WARN] {res['warnings'][-1]}")
    say('A1 时间推断：' + '、'.join(methods))

    env = tt.check_tt_env()
    res['tt_env'] = {k: env.get(k) for k in
                     ('ok', 'exe', 'version', 'how', 'problems')}
    if not env['ok']:
        for x in env['problems']:
            say(f'[WARN] {x}')
    else:
        say(f"TreeTime {env['version']}（{env['exe']}）")
    argv = [env['exe']] if env['ok'] else None

    aln = aln or (prep.get('outputs') or {}).get('prep_fasta')
    tree = tree or (prep.get('outputs') or {}).get('prep_nwk')
    dates = dates or (prep.get('outputs') or {}).get('dates')
    if not tree:
        res['error'] = '没有可用的树（A0 未产出 prep.nwk）—— 定年必须有树'
        return res
    seq_len = int(seq_len or prep.get('seq_len_modal') or 0)

    out = res['outputs']
    stats = res['stats']

    # ---------- TreeTime ----------
    if 'treetime' in methods:
        if not argv:
            res['warnings'].append('TreeTime 不可用 → 跳过该通路')
        else:
            prog('clock', 0.15, 'TreeTime 分子钟定年')
            say(f'TreeTime 定年（seq_len={seq_len}, relax={relax or "严格钟默认"}, '
                f'reroot={reroot}'
                + (f', clock_rate={clock_rate}（固定，跳过速率估计）'
                   if clock_rate else '') + '）')
            c = tt.run_clock(argv, os.path.join(out_dir, 'treetime'),
                             tree=tree, aln=aln, dates=dates, seq_len=seq_len,
                             relax=relax, reroot=reroot, keep_root=keep_root,
                             confidence=confidence, time_marginal=time_marginal,
                             clock_filter=clock_filter, rng_seed=rng_seed,
                             clock_rate=clock_rate,
                             timeout=timeout, say=say)
            res['treetime'] = _slim_clock(c)
            res['warnings'] += c.get('warnings') or []
            if not c['ok']:
                res['warnings'].append('TreeTime 定年失败：' + str(c.get('error')))
                say('[WARN] TreeTime 定年失败：' + str(c.get('error')))
            else:
                ck = c['clock']
                # 口径：固定速率下**不许**把 model 写成"严格钟（默认）" ——
                # 那会让读者以为是数据估出来的速率。
                if clock_rate:
                    _model = '严格钟（速率固定，用户指定）'
                    _note = ('速率由用户固定为 %g，**不是 ML 估计**；'
                             'r² 为根到尾回归（与 DRT 不是一回事）'
                             % float(clock_rate))
                else:
                    _model = ('非相关松弛钟' if relax and not relax[1]
                              else ('松弛钟' if relax else '严格钟（默认）'))
                    _note = 'ML 分子钟；r² 为根到尾回归（与 DRT 不是一回事）'
                if clock_rate:
                    res['notes'].append(
                        '分子钟速率已**由用户固定**为 %g 每位点每年：'
                        'TreeTime 这次没有做速率估计，引用该速率时必须注明'
                        '「固定速率（慢演化兜底）」，不能写成数据推断出的结果。'
                        % float(clock_rate))
                stats.append({
                    'method': 'TreeTime', 'model': _model,
                    'rate': ck.get('rate'), 'rate_unit': '每位点每年',
                    'r2': ck.get('r2'), 'tmrca_calendar': None,
                    'tmrca_ybp': None, 'cov': None, 'loglik': None,
                    'seconds': None,
                    'note': _note})
                res['r2'] = ck.get('r2')
                for tag, fn in (('timetree_nexus', 'timetree.nexus'),
                                ('timetree_nwk', 'timetree.nwk'),
                                ('divergence_nexus', 'divergence_tree.nexus'),
                                ('molecular_clock', 'molecular_clock.txt'),
                                ('rtt_pdf', 'root_to_tip_regression.pdf'),
                                ('tree_pdf', 'timetree.pdf'),
                                ('auspice', 'auspice_tree.json'),
                                ('ancestral_fasta', 'ancestral_sequences.fasta')):
                    p = os.path.join(out_dir, 'treetime', fn)
                    if os.path.isfile(p):
                        out['tt_' + tag] = p
                # 接口产物：timetree.nwk（A2/A3 唯一入口）
                src = os.path.join(out_dir, 'treetime', 'timetree.nwk')
                if os.path.isfile(src):
                    dst = os.path.join(out_dir, 'timetree.nwk')
                    with safe_open(src) as f:
                        body = f.read()
                    if isinstance(body, bytes):
                        body = body.decode('utf-8', 'replace')
                    with open_write(dst) as f:
                        f.write(body)
                    out['timetree_nwk'] = dst
                    res['timetree_source'] = 'TreeTime'
                    # ---- tMRCA：TreeTime 的 clock 输出**没有** tMRCA 字段 ----
                    # 时间树的枝长严格＝年（自洽性实测 0 条边不自洽），所以
                    # `root_year = date_tip − dist(root→tip)` 直接给出 tMRCA 的
                    # 日历年；`max(date) − root_year` 就是距今年数。
                    try:
                        with safe_open(dst) as f:
                            troot = tt.parse_annotated_newick(f.read())
                        dmap0 = _load_dates_map(dates)
                        _nd, ry = tt.node_dates_from_branches(troot, dmap0)
                        if ry is not None and dmap0:
                            res['root_year'] = round(ry, 4)
                            stats[-1]['tmrca_calendar'] = round(ry, 4)
                            stats[-1]['tmrca_ybp'] = round(max(dmap0.values()) - ry, 4)
                            say(f"tMRCA（由时间树枝长反推）＝ {ry:.2f} 日历年"
                                f"（距今 {stats[-1]['tmrca_ybp']:.2f} 年）")
                        else:
                            res['warnings'].append(
                                '时间树算不出根年代（叶日期覆盖不足）→ '
                                'tMRCA 留空，不拿相对量冒充日历年')
                    except (OSError, ValueError) as e:
                        res['warnings'].append(
                            f'tMRCA 未算出：{type(e).__name__}: {e}')
                # divergence 树也规范化一份（A4 祖先序列按替换单位喂更规范）
                dsrc = os.path.join(out_dir, 'treetime', 'divergence_tree.nexus')
                if os.path.isfile(dsrc):
                    try:
                        with safe_open(dsrc) as f:
                            dnwk = tt.nexus_tree_newick(f.read())
                        ddst = os.path.join(out_dir, 'divergence_tree.nwk')
                        with open_write(ddst) as f:
                            f.write(dnwk.rstrip().rstrip(';') + ';\n')
                        out['divergence_nwk'] = ddst
                    except (OSError, ValueError) as e:
                        res['warnings'].append(
                            f'divergence 树未规范化：{type(e).__name__}: {e}')


    if not out.get('timetree_nwk'):
        res['error'] = ('TreeTime 没有产出时间树 —— 时间与地理推断无法继续。'
                        '（A1 已收敛到 TreeTime 单通路，不再回退到已摘除的'
                        ' LSD2/treedater）'
                        + '；'.join(res['warnings'][-3:]))
        _dump(res, os.path.join(out_dir, 'clock_report.json'))
        return res

    # ---------- 数值合理性闸门 ----------
    # ⚠️ 来历（保留这条教训，与具体实现无关）：实测（RSV 209 + FastTree 树）时
    # LSD2 曾给出 rate=1e-10、tMRCA=-2.51e8（−2.5 亿年）这种**拟合失败**的结果，
    # 而 `lsd2_dating` 不抛异常、也不带 `ok` 字段 —— 不看数值就会把垃圾写进
    # clock_stats.tsv。LSD2 通路虽已于 2026-09-17 摘除，但 TreeTime 同样可能
    # 在时间信号很弱时把速率推到荒谬区间，闸门照旧保留。
    # 这里按「速率必须为正 + tMRCA 必须落在采样年域附近」做闸门，
    # 可疑项标 suspect 并说明理由，**不删除**（删了就查不到发生了什么）。
    _gate_clock_stats(stats, _load_dates_map(dates), res)

    # ---------- 时间信号：RTT 散点 + DRT ----------
    qc_dir = os.path.join(out_dir, 'clock_qc')
    os.makedirs(qc_dir, exist_ok=True)
    res['qc_dir'] = qc_dir
    if do_drt:
        prog('clock', 0.8, f'时间信号检验（DRT {drt_perm} 次置换）')
        dmap = _load_dates_map(dates)
        try:
            from Virus_Platform_Core import rtt_drt
            dr = rtt_drt.date_randomization_test(
                tree, dmap, n_perm=int(drt_perm), seed=int(drt_seed))
            res['drt'] = {k: dr.get(k) for k in (
                'real_r2', 'real_slope', 'n_valid', 'n_invalid',
                'r2_percentile', 'p_value', 'n_ge', 'passed',
                'passed_percentile95', 'r2_best_root', 'old_threshold_pass',
                'sign_consistent', 'conclusion', 'verdict', 'verdict_dissent')}
            res['drt']['warnings'] = dr.get('warnings') or []
            _dump(dr, os.path.join(qc_dir, 'drt.json'))
            out['drt_json'] = os.path.join(qc_dir, 'drt.json')
            pv = dr.get('p_value')
            say(f"DRT：R²={dr.get('real_r2')}，p={pv}，判据={'通过' if dr.get('passed') else '不通过'}")
            if not dr.get('passed'):
                res['warnings'].append(
                    '**DRT 日期随机化检验未通过** —— 时间信号不成立，'
                    '本阶段的速率/tMRCA 只能作为记录，**不得**作为定年结论引用。'
                    '（红线 5：时间信号未过 DRT 不输出定年结论）')
            res['warnings'] += dr.get('warnings') or []
        except Exception as e:
            res['warnings'].append(f'DRT 检验未完成：{type(e).__name__}: {e}')
            say(f'[WARN] DRT 未完成：{e}')

    # RTT 散点（用 divergence 树，枝长＝替换/位点，与 R² 同口径）
    prog('clock', 0.9, '根到尾散点')
    try:
        scat_tree = out.get('divergence_nwk') or out.get('timetree_nwk')
        with safe_open(scat_tree) as f:
            sroot = tt.parse_annotated_newick(f.read())
        dist = tt.tree_distances(sroot)
        dmap = _load_dates_map(dates)
        pts = []
        for n in tt.iter_pre(sroot):
            if n['children']:
                continue
            if n['name'] in dmap:
                pts.append((n['name'], dmap[n['name']], dist[id(n)]))
        if len(pts) >= 3:
            p = os.path.join(qc_dir, 'rtt_scatter.tsv')
            with open_write(p) as f:
                f.write('name\tdate\tdistance\tsource\n')
                for nm, y, d in sorted(pts, key=lambda x: x[1]):
                    f.write(f'{nm}\t{y}\t{d}\t{scat_tree}\n')
            out['rtt_scatter'] = p
            try:
                from Virus_Platform_Core import rtt_drt
                fit = rtt_drt.rtt_fit({nm: d for nm, _y, d in pts},
                                      {nm: y for nm, y, _d in pts})
                if fit:
                    res['rtt_fit'] = {'r2': fit['r2'], 'slope': fit['slope'],
                                      'intercept': fit['intercept'], 'n': fit['n']}
                    res['rtt_points'] = [[round(y, 4), round(d, 8)]
                                         for _nm, y, d in sorted(pts,
                                                                 key=lambda x: x[1])]
                    say(f"根到尾回归（固定根位）：R²={fit['r2']:.4f}, "
                        f"slope={fit['slope']:.4g}")
            except Exception as e:
                res['warnings'].append(f'RTT 回归未完成：{type(e).__name__}: {e}')
    except Exception as e:
        res['warnings'].append(f'RTT 散点未产出：{type(e).__name__}: {e}')

    # ---------- 时间树坐标对象（前端画图；与 RTT/phylogeo 两卡同构）----------
    # 统一放在这里、而不是通路分支里：横轴口径只由**接口产物**
    # `clock/timetree.nwk` 决定。A1 收敛到 TreeTime 后产物只有一份来源，
    # 不会再出现"图画的是 TreeTime、表里写的是别家"这种对不上的情况。
    # 缺年代时返回只带 note 的对象 → 前端显示说明文字
    # 而不是留一个空 div（静默空白会让人分不清"没跑"和"跑坏了"）。
    try:
        tto = build_clock_timetree(
            out.get('timetree_nwk'), dates,
            source_label=res.get('timetree_source') or '')
        if tto:
            res['timetree'] = tto
            if tto.get('x'):
                say(f"时间树（{res.get('timetree_source')}）：{tto['n_tips']} 叶，"
                    f"年代区间 {tto['min_x']}–{tto['max_x']}"
                    + (f"（{tto['n_backfilled']} 个节点无年代）"
                       if tto.get('n_backfilled') else ''))
            else:
                res['warnings'].append('时间树坐标未产出：'
                                       + str(tto.get('note')))
    except Exception as e:      # 附加产物，失败不该拖垮已经算好的定年结果
        res['warnings'].append(f'时间树坐标未产出：{type(e).__name__}: {e}')
        say(f'[WARN] 时间树坐标未产出：{e}')

    # ---------- clock_stats.tsv（带方法列，红线 1）----------
    p = os.path.join(out_dir, 'clock_stats.tsv')
    cols = ['method', 'model', 'rate', 'rate_unit', 'r2', 'tmrca_calendar',
            'tmrca_ybp', 'cov', 'loglik', 'seconds', 'suspect',
            'suspect_reason', 'note']
    with open_write(p) as f:
        f.write('\t'.join(cols) + '\n')
        for s in stats:
            f.write('\t'.join(
                str(s.get(c) or '') if c in ('method', 'model', 'rate_unit',
                                             'note', 'suspect_reason')
                else ('1' if c == 'suspect' and s.get(c) else
                      ('0' if c == 'suspect' else _fnum(s.get(c))))
                for c in cols) + '\n')
    out['clock_stats'] = p
    res['n_methods_ok'] = len(stats)
    res['n_methods_usable'] = sum(1 for s in stats if not s.get('suspect'))
    # ⚠️ 2026-10-03 删：这里原有一条 `if len(stats) > 1:` 的"多通路速率不可比"
    #    告警，但 A1 自 2026-09-17 收敛到 CLOCK_METHODS=('treetime',) 之后
    #    stats 至多 1 条 → 分支**永不成立**（死代码），故删除。
    res['ok'] = not res['error']
    prog('clock', 1.0, 'A1 完成')
    say(f"A1 完成：时间树来源 {res['timetree_source']}；成功通路 {len(stats)} 条")
    _dump(res, os.path.join(out_dir, 'clock_report.json'))
    return res


def _gate_clock_stats(stats, dates_map, res):
    """给定年结果加数值合理性闸门（就地改 `stats`）。

    判据（都能从采样年域直接推出来，不需要外部知识）：
      ① 速率必须是**有限正数** —— 1e-10 这种是拟合塌缩，不是"很慢的钟"
      ② tMRCA 的日历年必须 ≤ 最新采样年（祖先不可能比后代晚）
      ③ tMRCA 不能比最老样品还老 10 倍采样跨度以上（病毒不是恐龙）
      ④ tMRCA 若给了但超出 1000–2100，基本是单位/数值事故

    可疑项标 `suspect=True` + `suspect_reason`，并**保留**在表里 ——
    删掉就没人知道发生了什么。同时把理由写进 `res['warnings']`。
    """
    if not stats:
        return stats
    ys = [v for v in (dates_map or {}).values() if v is not None]
    hi = max(ys) if ys else None
    lo = min(ys) if ys else None
    span = (hi - lo) if (hi is not None and lo is not None) else None
    for s in stats:
        reasons = []
        r = s.get('rate')
        try:
            rf = float(r) if r is not None else None
        except (TypeError, ValueError):
            rf = None
        if rf is None:
            reasons.append('速率缺失')
        elif rf != rf or rf in (float('inf'), float('-inf')):
            reasons.append('速率非有限值')
        elif rf <= 0:
            reasons.append(f'速率不是正数（{r}）')
        elif rf < 1e-9:
            reasons.append(f'速率 {r} 小到不合理（拟合塌缩到下界，不是"很慢的钟"）')
        tm = s.get('tmrca_calendar')
        if tm is not None:
            try:
                tf = float(tm)
            except (TypeError, ValueError):
                tf = None
            if tf is not None:
                if hi is not None and tf > hi + 1:
                    reasons.append(f'tMRCA（{tf:.6g}）晚于最新采样年（{hi:.4g}）—— 祖先比后代晚，不可能')
                elif lo is not None and span and tf < lo - 10 * max(span, 1.0):
                    reasons.append(f'tMRCA（{tf:.6g}）比最老样品（{lo:.4g}）还早 '
                                   f'>{10:.0f} 倍采样跨度 —— 数值事故，不是真实年代')
                elif not (1000 <= tf <= 2100):
                    reasons.append(f'tMRCA（{tf:.6g}）不在 1000–2100 之间 —— 单位或数值异常')
        if reasons:
            s['suspect'] = True
            s['suspect_reason'] = '；'.join(reasons)
            res['warnings'].append(
                f"{s.get('method')} 的结果**不可用**：" + '；'.join(reasons)
                + '。该行已标 suspect 保留在 clock_stats.tsv 里备查，'
                  '但**不要**引用它的数字（换树/换根位/换 seq_len 后重跑）')
    res['n_methods_usable'] = sum(1 for s in stats if not s.get('suspect'))
    return stats


def _slim_clock(c):
    # clock_filter/time_marginal/reroot/argv：定年的关键开关必须随结果落盘
    # 可追溯（clock_filter 缺省时 TreeTime CLI 默认 4.0 会剔叶，与 RTT 卡的
    # 显式 0 语义不同，丢失后无法区分两种口径）。
    return {k: c.get(k) for k in ('ok', 'error', 'returncode', 'clock', 'relax',
                                  'coalescent', 'skyline_meta', 'warnings',
                                  'stdout_tail', 'clock_rate_fixed',
                                  'clock_filter', 'time_marginal', 'reroot',
                                  'argv')}


def _lsd2_tree_path(lr):
    """【已停用，2026-09-17】从 lsd2_dating 的返回里找时间树路径。

    A1 定年收敛到 TreeTime 后不再有"退用 LSD2 时间树"的分支，此函数已无调用点。
    底层实现（`phylogeo.lsd2_dating`）按约定保留在仓库里未接线 —— 将来若要恢复
    LSD2 通路，这里也一并接回即可。刻意保留函数体而不是删掉：删了会让"曾经有过
    这条通路"的痕迹消失，恢复时得从备份里翻。
    """
    if not isinstance(lr, dict):
        return None
    for k in ('dated_tree', 'tree', 'out_tree', 'nwk'):
        v = lr.get(k)
        if isinstance(v, str) and os.path.isfile(v):
            return v
    outs = lr.get('outputs') or {}
    if isinstance(outs, dict):
        for v in outs.values():
            if isinstance(v, str) and v.lower().endswith(('.nwk', '.tree',
                                                          '.nexus', '.newick')):
                return v
    return None


def _load_dates_map(dates_path):
    """日期表 → {叶名: 十进制年}。

    ⚠️ 2026-10-03：改用 `csv.reader` 解析 —— 写端是 `treetime_ml.write_dates_csv`
    的 `_csv_cell`（含逗号/引号/换行的单元格会**加引号**），原来手撕
    `split(',')` + `strip('"')` 遇到 `"a,b"` 这种叶名会切错列、数据静默错位。
    读写两端统一走 csv 口径。
    """
    out = {}
    if not dates_path or not os.path.isfile(dates_path):
        return out
    with io.open(dates_path, encoding='utf-8-sig', errors='replace',
                 newline='') as f:
        for i, cells in enumerate(csv.reader(f)):
            if not cells:                       # 空行（csv 给出 []）
                continue
            if i == 0 and (cells[0] or '').strip().lower().lstrip('#') == 'name':
                continue
            name = (cells[0] or '').strip()
            if len(cells) >= 2 and name:
                v, _p = to_decimal_year(cells[1].strip())
                if v is not None:
                    out[name] = v
    return out


def build_clock_timetree(timetree_nwk, dates_path, source_label='',
                         tip_dates=None):
    """接口产物 `clock/timetree.nwk` → 前端可画的**直角时间树坐标对象**。

    返回 `phylogeo.time_tree_segments` 的那套 dict（x/y/tips/min_x/max_x/
    n_tips/n_internal/n_backfilled/axis_label/source/note），与 phylogeo 卡
    **同一个形状** —— 前端 `drawTimeTree` 直接吃，
    不需要为 phylodyn 单写一套画法。

    ⚠️ 横轴＝**日历年**，靠枝长反推（`treetime_ml.node_dates_from_branches`）：
    时间树的枝长单位＝年，所以 `date_node = date_root + dist(root→node)`。
    不读 `[&date=...]` 注释 —— mugration 会把标注树上的 date 覆盖成区划（实测），
    从枝长算既不依赖命名也不依赖注释。

    ⚠️ 叶日期匹配不到（树名与日期表对不上 / 定年结果没写节点年代）时**不产出坐标**，
    改回一个只带 `note` 的 dict —— 前端据此显示说明文字而不是留一个空 div，
    也不拿"相对量"冒充日历年。
    """
    p = timetree_nwk
    if not p or not os.path.isfile(p):
        return {'x': [], 'y': [], 'tips': [], 'n_tips': 0, 'n_internal': 0,
                'n_backfilled': 0, 'axis_label': '', 'source': '',
                'note': '接口产物 clock/timetree.nwk 不存在 → 不画时间树'}
    try:
        with safe_open(p) as f:
            body = f.read()
        if isinstance(body, bytes):
            body = body.decode('utf-8', 'replace')
    except OSError as e:
        return {'x': [], 'y': [], 'tips': [], 'n_tips': 0, 'n_internal': 0,
                'n_backfilled': 0, 'axis_label': '', 'source': '',
                'note': f'时间树读不回：{type(e).__name__}: {e}'}
    # 产物可能是 newick 也可能是 nexus（历史 run 里存在过 nexus 形态）——
    # `parse_annotated_newick` 自己认 `#NEXUS` 头，这里不用先转。
    try:
        root = tt.parse_annotated_newick(body)
    except (ValueError, AttributeError, KeyError) as e:
        return {'x': [], 'y': [], 'tips': [], 'n_tips': 0, 'n_internal': 0,
                'n_backfilled': 0, 'axis_label': '', 'source': '',
                'note': f'时间树解析失败：{type(e).__name__}: {e}'}

    dmap = tip_dates if tip_dates is not None else _load_dates_map(dates_path)
    x_of, root_year = tt.node_dates_from_branches(root, dmap or {})
    if not x_of:
        return {'x': [], 'y': [], 'tips': [], 'n_tips': 0, 'n_internal': 0,
                'n_backfilled': 0, 'axis_label': '', 'source': '',
                'note': ('时间树的叶名与日期表匹配不上（或带日期叶不足半数）'
                         '→ 算不出日历年坐标，不画。'
                         f'树 {len(dmap or {})} 个日期可用')}

    from Virus_Platform_Core import phylogeo as pg
    tto = pg.time_tree_segments(
        root, x_of, axis_label='采样年（日历年）',
        source=os.path.basename(p))
    tto['root_year'] = round(root_year, 4) if root_year is not None else None
    tto['n_dated'] = len(x_of)
    tto['n_date_map'] = len(dmap or {})
    # 口径写清楚：横轴是**枝长反推**出来的日历年（叶钉在采样年上），
    # 不是读 `[&date=...]` 注释 —— mugration 之后那份标注树上的 date 已被
    # 覆盖成区划，照读会得到一串区划名当"年份"。
    # 另：`time_tree_segments` 的 note 自带「时间树（label）：」前缀，这里剥掉，
    # 否则前端提示行会出现两次"时间树（…）"。
    _inner = re.sub(r'^时间树（[^）]*）：', '', str(tto.get('note') or ''))
    tto['note'] = ('横轴＝由时间树枝长反推的日历年'
                   '（date_node = date_root + dist(root→node)，叶钉在采样年上）；'
                   + _inner)
    if source_label:
        tto['dating_source'] = source_label
    return tto


# ================================================================ A2 phylogeo

def stage_phylogeo(clock, prep, out_dir, trait='region',
                   sampling_bias_correction='auto', pc=None,
                   cross_check=True, rssp_bs=0, seed=0, confidence=True,
                   coord_table=None, timeout=3600, log=None, progress=None):
    """A2 地理推断：TreeTime mugration（ML）**与**平台 Fitch 简约**并列**。

    - 接口输入：`clock/timetree.nwk` + `prep/states.csv`
    - 产物：`corridors.tsv`（`ml` / `fitch` 两列并列）、`events_year.tsv`
      （`Year, Origin, Destination`，借 Mugration-Analysis）、`GTR.txt`、
      `confidence.csv`、`annotated_tree.nexus`、`phylogeo_report.json`

    `coord_table`（可选）：区划坐标表（`区划/纬度/经度` 三列，兼容 VirPhyKit
    Spatial.txt）。**只影响「迁移弧线地图」的落点**，不参与任何统计量：
    迁移次数 / 方向 / 走廊排序全部只由树与状态表决定。没给坐标表时用
    `phylogeo.BUILTIN_REGION_COORDS` 的**近似质心兜底**（出处标 `builtin`，
    前端必须显式写明「近似」）；两者都没有的区划进 `coord_missing`，
    **绝不静默少画点**。

    ⚠️ 两条口径**不是同一个统计量**：`ml` 是 ML 点估计，`fitch` 是最少变化数下的
    一种悬挂（tie-break 依赖树序）。实测 RSV 209：总数 70 vs 71（差 1.4%），
    但**逐走廊方向差别很大**（NCS→KOR 24 vs 2、CCD→NCS 3 vs 17）——
    所以两列都要给，不许只留一个，也不许相加。
    """
    say, prog = _say(log), _prog(progress)
    os.makedirs(out_dir, exist_ok=True)
    res = {'stage': 'phylogeo', 'ok': False, 'error': None, 'out_dir': out_dir,
           'trait': trait, 'warnings': [], 'notes': [], 'outputs': {},
           'sampling_bias_correction': sampling_bias_correction}
    say('A2 地理推断（ML mugration + Fitch 并列）')

    out = res['outputs']
    cout = clock.get('outputs') or {}
    tree = cout.get('timetree_nwk')
    states = (prep.get('outputs') or {}).get('states')
    if not tree or not os.path.isfile(str(tree)):
        res['error'] = '缺少时间树（clock/timetree.nwk）—— A2 的唯一接口输入'
        return res
    if not states or not os.path.isfile(str(states)):
        res['error'] = '缺少状态表（prep/states.csv）'
        return res
    res['tree_in'] = tree
    res['states_in'] = states

    env = tt.check_tt_env()
    argv = [env['exe']] if env['ok'] else None

    # ---------- TreeTime mugration ----------
    if not argv:
        res['warnings'].append('TreeTime 不可用 → 只有 Fitch 一条口径')
    else:
        prog('phylogeo', 0.2, 'TreeTime mugration（ML 地理重构）')
        xc = cout.get('divergence_nwk') if cross_check else None
        m = tt.run_mugration(argv, tree, states, os.path.join(out_dir, 'ml'),
                             attribute=trait, confidence=confidence,
                             sampling_bias_correction=sampling_bias_correction,
                             pc=pc, also_tree=xc, rng_seed=seed,
                             timeout=timeout, say=say)
        res['ml'] = {k: m.get(k) for k in ('ok', 'error', 'returncode', 'gtr',
                                           'confidence_summary',
                                           'sampling_bias_correction',
                                           'sampling_bias_correction_suggested',
                                           'stdout_tail')}
        res['warnings'] += m.get('warnings') or []
        if not m['ok']:
            res['warnings'].append('mugration 失败：' + str(m.get('error')))
        else:
            for tag, fn in (('gtr', 'GTR.txt'), ('confidence_csv', 'confidence.csv'),
                            ('annotated_nexus', 'annotated_tree.nexus')):
                p = os.path.join(out_dir, 'ml', fn)
                if os.path.isfile(p):
                    out['ml_' + tag] = p
            # 节点年代（供 events_year）
            with safe_open(tree) as f:
                root = tt.parse_annotated_newick(f.read())
            nd, root_year = tt.node_dates_from_branches(
                root, _load_dates_map(
                    (prep.get('outputs') or {}).get('dates')))
            res['root_year'] = root_year
            cc = tt.count_corridors(root, m['state_of'], node_year=nd)
            res['ml_n_changes'] = cc['n_changes']
            res['ml_edges_missing'] = cc['n_edges_missing']
            res['ml_counts'] = {'%s->%s' % k: v for k, v in cc['counts'].items()}
            res['ml_events'] = cc['rows']
            # 对拍
            if m.get('cross_check'):
                xcc = tt.cross_check_corridors(
                    root, m['cross_check'],
                    cc_tree_path=(m['cross_check'] or {}).get('tree'),
                    primary={'n_changes': cc['n_changes']})
                res['cross_check'] = xcc
                if xcc and xcc.get('align_warning'):
                    res['warnings'].append(xcc['align_warning'])
                if xcc and xcc.get('divergent'):
                    w = ('枝长口径敏感性告警：换成 divergence 树跑同一次 mugration，'
                         '迁移次数从 %d 变成 %d（差 %.0f%%）—— '
                         '两棵树**不成比例**说明时间信号弱，地理重建对枝长口径敏感，'
                         '走廊计数只能当量级参考，方向性结论请交 BEAST 求后验'
                         % (cc['n_changes'], xcc['n_changes'],
                            100 * (xcc.get('rel_diff') or 0)))
                    res['warnings'].append(w)
                    say('[WARN] ' + w)

    # ---------- Fitch（同一棵树，平台口径） ----------
    prog('phylogeo', 0.6, 'Fitch 简约重构（平台口径，同树对拍）')
    try:
        from Virus_Platform_Core import phylogeo as pg
        with safe_open(tree) as f:
            froot = tt.parse_annotated_newick(f.read())
        leaf_states = _load_states_map(states)
        _, trans, nf = pg.fitch_mugration(froot, leaf_states)
        fc = {}
        frows = []
        nd2, _ry = tt.node_dates_from_branches(
            froot, _load_dates_map((prep.get('outputs') or {}).get('dates')))
        for pnm, cnm, a, b in trans:
            fc[(a, b)] = fc.get((a, b), 0) + 1
        # 逐条事件（带年份）：从 fitch 的 child 名找回节点
        by_name = {}
        for n in tt.iter_pre(froot):
            if n['name']:
                by_name.setdefault(n['name'], n)
        for pnm, cnm, a, b in trans:
            cn = by_name.get(cnm)
            frows.append({'from': a, 'to': b, 'parent': pnm, 'child': cnm,
                          'year': (round(nd2[id(cn)], 3)
                                   if cn is not None and nd2.get(id(cn)) is not None
                                   else None)})
        res['fitch_n_changes'] = nf
        res['fitch_counts'] = {'%s->%s' % k: v for k, v in fc.items()}
        res['fitch_events'] = frows
        say(f'Fitch（同树）迁移 {nf} 次')
        # 可选 RSPP bootstrap（复用平台实现，很贵，默认关）
        if rssp_bs:
            prog('phylogeo', 0.75, f'RSPP bootstrap {rssp_bs} 复本')
            try:
                aln = (prep.get('outputs') or {}).get('prep_fasta')
                # ⚠️ 2026-10-03 修：`observed` 参数要的是**迁移数（数值）**，
                #    不是 dict —— 原实现把 {'n_changes': …} 传进去，`lo <= observed <= hi`
                #    直接 TypeError，被下面的 except 吞成"RSPP 未完成"，于是 A2 的
                #    bootstrap 置信区间**永远拿不到**。这里传 Fitch 的真实迁移数 nf。
                rb = pg.rssp_bootstrap(aln, leaf_states,
                                       os.path.join(out_dir, '_rssp'),
                                       'fasttree', int(rssp_bs), seed,
                                       None, nf)
                # 返回键按 rssp_bootstrap 的实际返回结构取（见 phylogeo.py）。
                res['rssp'] = {k: rb.get(k) for k in (
                    'n_boot', 'n_ok', 'failed', 'observed', 'observed_in_range',
                    'transitions', 'corridor_support', 'note', 'clades')}
                say(f"RSPP 完成：{rb.get('n_ok')}/{rb.get('n_boot')} 复本有效"
                    + ('' if rb.get('observed_in_range') is not False
                       else f"；点估计 {nf} 落在 bootstrap 区间外"
                            '（迁移数对拓扑选择敏感，解读以区间为准）'))
            except Exception as e:
                res['warnings'].append(f'RSPP 未完成：{type(e).__name__}: {e}')
    except Exception as e:
        res['warnings'].append(f'Fitch 通路异常：{type(e).__name__}: {e}')
        say(f'[WARN] Fitch 异常：{e}')

    if not (res.get('ml_counts') or res.get('fitch_counts')):
        res['error'] = '两条口径都没算出走廊 —— 见 warnings'
        _dump(res, os.path.join(out_dir, 'phylogeo_report.json'))
        return res

    # ---------- corridors.tsv（两列并列） ----------
    ml_c = {tuple(k.split('->')): v for k, v in (res.get('ml_counts') or {}).items()}
    ft_c = {tuple(k.split('->')): v for k, v in (res.get('fitch_counts') or {}).items()}
    tab = tt.corridors_table(ml_c, ft_c)
    p = os.path.join(out_dir, 'corridors.tsv')
    with open_write(p) as f:
        f.write('from\tto\tml\tsum_ml\tfitch\tsum_fitch\tdelta\n')
        tot_ml = res.get('ml_n_changes')
        tot_ft = res.get('fitch_n_changes')
        for r in tab:
            f.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (
                r['from'], r['to'], _fnum(r['ml'], 8), _fnum(tot_ml, 8),
                _fnum(r['fitch'], 8), _fnum(tot_ft, 8), _fnum(r['delta'], 8)))
    out['corridors'] = p

    # ---------- events_year.tsv（Year, Origin, Destination） ----------
    p = os.path.join(out_dir, 'events_year.tsv')
    with open_write(p) as f:
        f.write('Year\tOrigin\tDestination\tMethod\tChildNode\tBranchLength\n')
        for r in sorted(res.get('ml_events') or [],
                        key=lambda x: (x['year'] is None, x['year'] or 0)):
            f.write('%s\t%s\t%s\tTreeTime-ML\t%s\t%s\n' % (
                _fnum(r['year'], 8), r['from'], r['to'], r['child'],
                _fnum(r['branch_len'], 10)))
        for r in sorted(res.get('fitch_events') or [],
                        key=lambda x: (x['year'] is None, x['year'] or 0)):
            f.write('%s\t%s\t%s\tFitch\t%s\t\n' % (
                _fnum(r['year'], 8), r['from'], r['to'], r['child']))
    out['events_year'] = p
    res['n_events_ml'] = len(res.get('ml_events') or [])
    res['n_events_fitch'] = len(res.get('fitch_events') or [])
    res['n_corridors'] = len(tab)

    # ---------- 采样偏倚检查（红线 4） ----------
    counts = _load_state_counts(states)
    tot = sum(counts.values()) or 1
    res['sampling_bias'] = {
        'state_counts': counts,
        'shares': {k: round(v / tot, 4) for k, v in counts.items()},
        'max_share': round(max(counts.values()) / tot, 4) if counts else None,
        'min_share': round(min(counts.values()) / tot, 4) if counts else None,
        'corrected': bool(res.get('ml', {}).get('sampling_bias_correction')),
        'factor': res.get('ml', {}).get('sampling_bias_correction'),
    }
    if counts and len(counts) >= 2:
        mx, mn = max(counts.values()), min(counts.values())
        if mn and mx / mn >= 10:
            res['warnings'].append(
                f'采样极不均衡（最多的区划 {mx} 条 vs 最少 {mn} 条，{mx/mn:.0f} 倍）——'
                '未校正时迁移方向的置信度会系统性虚高；本次'
                + ('已按系数 %.4g 校正' % res['sampling_bias']['factor']
                   if res['sampling_bias']['corrected'] else '**未校正**'))

    # ---------- 区域坐标（迁移弧线地图的落点；不影响任何统计量） ----------
    # 口径：用户坐标表优先 → 内置质心表兜底 → 都没有就进 coord_missing 并显式报出。
    # 这一段**只**决定弧线画在哪：迁移次数 / 方向 / 走廊排序全部只由树与状态表
    # 决定，所以缺坐标不会改变任何结论，只会让图少画弧线（且必须说明少画了）。
    try:
        from Virus_Platform_Core import phylogeo as pgxy
        regions = sorted(set(counts)
                         | {r['from'] for r in tab} | {r['to'] for r in tab})
        user_coords, user_skipped = {}, []
        has_table = bool(coord_table) and os.path.isfile(str(coord_table))
        if has_table:
            user_coords, user_skipped = pgxy.parse_coord_table(coord_table)
            if user_skipped:
                res['warnings'].append(
                    '坐标表有 %d 行没解析出坐标（已跳过）：%s'
                    % (len(user_skipped),
                       '；'.join(str(x) for x in user_skipped[:3])))
            say('坐标表 %s：解析出 %d 个区划'
                % (os.path.basename(str(coord_table)), len(user_coords)))
        elif coord_table:
            res['warnings'].append('坐标表不存在：%s → 退回内置质心表' % coord_table)
        coords, csrc, cmiss = pgxy.merge_coords(regions, user_coords)
        res['coords'] = coords
        res['coord_src'] = csrc
        res['coord_missing'] = cmiss
        res['coord_regions'] = {r: int(counts.get(r, 0)) for r in regions}
        res['coord_table'] = (os.path.basename(str(coord_table)) if has_table
                              else None)
        n_builtin = sum(1 for v in csrc.values() if v == 'builtin')
        if cmiss:
            res['warnings'].append(
                '区划缺坐标 %d 个（%s）—— 迁移弧线地图会少画这些区划相关的弧线；'
                '可在「区域坐标表」补 区划/纬度/经度 三列后重跑'
                % (len(cmiss), '、'.join(cmiss[:5])))
        if n_builtin:
            res['notes'].append(
                '坐标：%d 个区划来自用户坐标表，%d 个用内置**近似质心**兜底'
                '（不是精确地理中心，引用时须写明出处）'
                % (len(coords) - n_builtin, n_builtin))
        say('坐标：%d/%d 个区划可用（内置质心兜底 %d 个，缺 %d 个）'
            % (len(coords), len(regions), n_builtin, len(cmiss)))
    except Exception as e:  # noqa: BLE001
        res['warnings'].append(
            f'坐标解析失败（不影响统计量，只是没有弧线地图）：'
            f'{type(e).__name__}: {e}')

    res['ok'] = not res['error']
    prog('phylogeo', 1.0, 'A2 完成')
    say(f"A2 完成：ML {res.get('ml_n_changes')} 次 / Fitch {res.get('fitch_n_changes')} 次"
        f"，走廊 {len(tab)} 条，事件 {res['n_events_ml']}（ML）")
    _dump(res, os.path.join(out_dir, 'phylogeo_report.json'))
    return res


def _load_states_map(states_path):
    """州/区划表 → {叶名: 状态}。

    ⚠️ 2026-10-03：同 `_load_dates_map`，改用 `csv.reader` —— 写端
    `treetime_ml.write_states_csv` 会给含逗号/引号的单元格加引号，手撕
    `split(',')` 会切错列。读写两端同一口径。
    """
    out = {}
    if not states_path or not os.path.isfile(states_path):
        return out
    with io.open(states_path, encoding='utf-8-sig', errors='replace',
                 newline='') as f:
        for i, cells in enumerate(csv.reader(f)):
            if not cells:
                continue
            if i == 0 and (cells[0] or '').strip().lower().lstrip('#') == 'name':
                continue
            name = (cells[0] or '').strip()
            if len(cells) >= 2 and name:
                out[name] = (cells[1] or '').strip() or 'Unknown'
    return out


def _load_state_counts(states_path):
    c = {}
    for v in _load_states_map(states_path).values():
        if v and v not in tt.MISSING_STATES:
            c[v] = c.get(v, 0) + 1
    return dict(sorted(c.items(), key=lambda x: -x[1]))


# ================================================================ A3 coalescent

def stage_coalescent(clock, prep, out_dir, n_skyline=None, gen_per_year=50.0,
                     rng_seed=0, timeout=3600, log=None, progress=None):
    """A3 群体动态：TreeTime `--coalescent skyline`（BSP 的**免 MCMC 预览版**）。

    ⚠️ **独立再跑一遍 clock**，不是复用 A1 的结果：`--coalescent skyline` 会改变
    速率估计（RSV 209 实测 3.07e-4 → 2.45e-4）。混在一行里说不清"这个速率是哪个
    模型的"，所以本阶段自己那一份速率单独记在 `coalescent_report.json` 里。

    ⚠️ `skyline.tsv` 的上下界是 TreeTime 自己的**近似区间**（±2 SD of the LH），
    **不是**贝叶斯 95% HPD；`N_e` 的绝对标度取决于 `--gen-per-year`（默认 50）。
    正式 95% 带仍需 BEAST 的 BSP（轨道 B）。
    """
    say, prog = _say(log), _prog(progress)
    os.makedirs(out_dir, exist_ok=True)
    res = {'stage': 'coalescent', 'ok': False, 'error': None, 'out_dir': out_dir,
           'gen_per_year': gen_per_year, 'n_skyline': n_skyline,
           'warnings': [], 'outputs': {}, 'skyline': [], 'skyline_meta': {}}
    say('A3 群体动态（TreeTime skyline，免 MCMC 预览）')
    env = tt.check_tt_env()
    if not env['ok']:
        res['error'] = 'TreeTime 不可用：' + '；'.join(env['problems'])
        return res
    cout = clock.get('outputs') or {}
    tree = cout.get('divergence_nwk') or cout.get('timetree_nwk')
    aln = (prep.get('outputs') or {}).get('prep_fasta')
    dates = (prep.get('outputs') or {}).get('dates')
    if not tree:
        res['error'] = '缺少可用的树（divergence / timetree）'
        return res
    prog('coalescent', 0.3, 'TreeTime skyline 拟合')
    c = tt.run_clock([env['exe']], out_dir, tree=tree, aln=aln, dates=dates,
                     seq_len=prep.get('seq_len_modal'), coalescent='skyline',
                     n_skyline=n_skyline or 10, gen_per_year=gen_per_year,
                     confidence=False, time_marginal='only-final',
                     rng_seed=rng_seed, timeout=timeout, say=say)
    res['clock'] = {k: c.get(k) for k in ('ok', 'error', 'returncode', 'clock',
                                          'coalescent', 'skyline_meta',
                                          'stdout_tail')}
    res['warnings'] += c.get('warnings') or []
    if not c['ok']:
        res['error'] = 'skyline 未产出：' + str(c.get('error'))
        _dump(res, os.path.join(out_dir, 'coalescent_report.json'))
        return res
    res['skyline'] = c.get('skyline') or []
    res['skyline_meta'] = c.get('skyline_meta') or {}
    for tag, fn in (('skyline_tsv', 'skyline.tsv'), ('skyline_pdf', 'skyline.pdf'),
                    ('timetree_nexus', 'timetree.nexus'),
                    ('molecular_clock', 'molecular_clock.txt')):
        p = os.path.join(out_dir, fn)
        if os.path.isfile(p):
            res['outputs'][tag] = p
    res['rate_this_pass'] = (c.get('clock') or {}).get('rate')
    res['r2_this_pass'] = (c.get('clock') or {}).get('r2')
    if (clock.get('stats') or [{}])[0].get('rate') and res['rate_this_pass']:
        try:
            a = float(clock['stats'][0]['rate'])
            b = float(res['rate_this_pass'])
            if a and abs(b - a) / a > 0.02:
                res['warnings'].append(
                    'A1（无天际线）与 A3（带 skyline）的速率不同：%.4g vs %.4g'
                    '（差 %.1f%%）—— 这是**模型不同**导致的，不是算错；'
                    '引用时务必写明是哪一次。' % (a, b, 100 * abs(b - a) / a))
        except (TypeError, ValueError, IndexError):
            pass
    res['warnings'].append(
        'skyline 的上下界是 TreeTime 的**近似区间**（±%s SD of the LH），'
        '**不是**贝叶斯 95%% HPD；N_e 的绝对标度取决于 `--gen-per-year` 假设'
        '（本次 %s）。正式 95%% 带请交 BEAST 的 BSP（轨道 B）。'
        % (res['skyline_meta'].get('sd_multiple', 2),
           res['skyline_meta'].get('gen_per_year', gen_per_year)))
    res['ok'] = True
    prog('coalescent', 1.0, 'A3 完成')
    say(f"A3 完成：{len(res['skyline'])} 个网格点，"
        f"本次速率 {res['rate_this_pass']}")
    _dump(res, os.path.join(out_dir, 'coalescent_report.json'))
    return res


# ================================================================ A4 ancestral

def stage_ancestral(clock, prep, out_dir, marginal=False, method_anc=None,
                    rng_seed=0, timeout=3600, log=None, progress=None):
    """A4 祖先序列：TreeTime `ancestral`。

    喂**替换单位**的树（`clock/divergence_tree.nwk`），因为 `ancestral` 要推断的是
    核苷酸 GTR —— 枝长口径按替换/位点是 TreeTime 文档的假定；拿不到 divergence 树
    时退回时间树（重建对枝长整体缩放不变，结果等价，只是少了那层规范）。

    平台此前只有 BEAST 产物导入路径，**本地从无到有**。
    """
    say, prog = _say(log), _prog(progress)
    os.makedirs(out_dir, exist_ok=True)
    res = {'stage': 'ancestral', 'ok': False, 'error': None, 'out_dir': out_dir,
           'warnings': [], 'outputs': {}, 'n_mutations': 0, 'tree_used': None}
    say('A4 祖先序列（TreeTime ancestral）')
    env = tt.check_tt_env()
    if not env['ok']:
        res['error'] = 'TreeTime 不可用：' + '；'.join(env['problems'])
        return res
    cout = clock.get('outputs') or {}
    tree = cout.get('divergence_nwk') or cout.get('timetree_nwk')
    aln = (prep.get('outputs') or {}).get('prep_fasta')
    if not tree or not aln:
        res['error'] = '缺少树或比对'
        return res
    res['tree_used'] = os.path.basename(str(tree))
    if not cout.get('divergence_nwk'):
        res['warnings'].append(
            '没有 divergence 树 → 退回时间树（枝长＝年）。重建对枝长整体缩放不变，'
            '结果等价；但若要做速率相关解读，请改用替换单位的树')
    prog('ancestral', 0.3, 'TreeTime ancestral 重构')
    a = tt.run_ancestral([env['exe']], aln, tree, out_dir, marginal=marginal,
                         method_anc=method_anc, rng_seed=rng_seed,
                         timeout=timeout, say=say)
    res['ancestral'] = {k: a.get(k) for k in ('ok', 'error', 'returncode',
                                              'n_mutations', 'stdout_tail')}
    res['warnings'] += a.get('warnings') or []
    if not a['ok']:
        res['error'] = 'ancestral 未产出：' + str(a.get('error'))
        _dump(res, os.path.join(out_dir, 'ancestral_report.json'))
        return res
    res['n_mutations'] = a.get('n_mutations') or 0
    for tag, fn in (('ancestral_fasta', 'ancestral_sequences.fasta'),
                    ('branch_mutations', 'branch_mutations.txt'),
                    ('annotated_nexus', 'annotated_tree.nexus'),
                    ('auspice', 'auspice_tree.json')):
        p = os.path.join(out_dir, fn)
        if os.path.isfile(p):
            res['outputs'][tag] = p
    res['ok'] = True
    prog('ancestral', 1.0, 'A4 完成')
    say(f"A4 完成：{res['n_mutations']} 处逐枝替换")
    _dump(res, os.path.join(out_dir, 'ancestral_report.json'))
    return res


# ================================================================ A5 report

def stage_report(stages, out_dir, params=None, log=None, progress=None):
    """A5 汇总溯源：`run_summary.json` + `config.resolved.json` + `RUN_REPORT.md`。

    照 pathogen-phylo-genomics-pipeline 的溯源做法：每个产物记 **sha256 前 16 位**
    （便于事后确认"交付的就是这一份"），并把**实际生效的配置**单独落一份
    （与用户填的表单区分开 —— 默认值/钳制后的值只有 resolved 才准）。
    """
    say, prog = _say(log), _prog(progress)
    os.makedirs(out_dir, exist_ok=True)
    say('A5 汇总溯源')
    prog('report', 0.4, '汇总产物与校验和')

    files = []
    for stage, sr in (stages or {}).items():
        if not isinstance(sr, dict):
            continue
        for tag, p in (sr.get('outputs') or {}).items():
            if isinstance(p, str) and os.path.isfile(p):
                try:
                    size = os.path.getsize(p)
                except OSError:
                    size = None
                files.append({'stage': stage, 'tag': tag, 'path': p,
                              'rel': os.path.relpath(p, os.path.dirname(out_dir))
                              .replace(os.sep, '/'),
                              'size': size, 'sha256_16': _sha256(p)})
    files.sort(key=lambda x: (x['stage'], x['tag']))

    warns = []
    for stage in STAGE_ORDER:
        sr = (stages or {}).get(stage) or {}
        for w in (sr.get('warnings') or []):
            warns.append({'stage': stage, 'text': str(w)})
    errs = [{'stage': k, 'text': str(v.get('error'))}
            for k, v in (stages or {}).items()
            if isinstance(v, dict) and v.get('error')]

    summary = {'tool': 'phylodyn_local', 'created': time.strftime('%Y-%m-%d %H:%M:%S'),
               'out_dir': os.path.dirname(out_dir),
               'stages': {k: {'ok': (v or {}).get('ok'),
                              'error': (v or {}).get('error')}
                          for k, v in (stages or {}).items()},
               'n_files': len(files), 'files': files,
               'warnings': warns, 'errors': errs}
    _dump(summary, os.path.join(out_dir, 'run_summary.json'))
    _dump({'params': params or {}, 'stages': list((stages or {}).keys()),
           'resolved_at': time.strftime('%Y-%m-%d %H:%M:%S')},
          os.path.join(out_dir, 'config.resolved.json'))

    md = _render_report(stages, params, files, warns, errs)
    with open_write(os.path.join(out_dir, 'RUN_REPORT.md')) as f:
        f.write(md)

    prog('report', 1.0, 'A5 完成')
    say(f'A5 完成：{len(files)} 个产物已登记校验和，{len(warns)} 条告警')
    return {'stage': 'report', 'ok': not errs, 'error': None,
            'out_dir': out_dir, 'warnings': [], 'n_files': len(files),
            'n_warnings': len(warns), 'n_errors': len(errs),
            'outputs': {'run_summary': os.path.join(out_dir, 'run_summary.json'),
                        'config_resolved': os.path.join(out_dir,
                                                        'config.resolved.json'),
                        'run_report': os.path.join(out_dir, 'RUN_REPORT.md')}}


def _render_report(stages, params, files, warns, errs):
    """人读版报告：用了什么方法、口径提示、警告、产物清单。"""
    L = ['# 本地时间与地理推断 · 运行报告', '']
    L.append(f"生成时间：{time.strftime('%Y-%m-%d %H:%M:%S')}")
    L.append('')
    L.append('> 本报告由平台「本地时间与地理推断流水线」（轨道 A）自动生成。'
             '**全程本机完成，未跑 MCMC**；凡涉及后验显著性的结论仍需 BEAST（轨道 B）。')
    L.append('')

    L.append('## 一、阶段状态')
    L.append('')
    L.append('| 阶段 | 名称 | 状态 |')
    L.append('|---|---|---|')
    for k in STAGE_ORDER:
        sr = (stages or {}).get(k)
        if k == 'report' and not isinstance(sr, dict):
            # A5 就是**这份报告本身** —— 它当然在运行，只是还没把自己写进
            # `stages`（stage_report 是被 _finish 调用后才回填的）。
            # 不特判的话报告会自称"A5 未运行"，等于交付物自我否定。
            L.append(f"| {STAGE_META[k]['code']} | {STAGE_META[k]['name']} | ✅ 完成（本报告） |")
            continue
        if not isinstance(sr, dict):
            L.append(f"| {STAGE_META[k]['code']} | {STAGE_META[k]['name']} | 未运行 |")
            continue
        st = '✅ 完成' if sr.get('ok') else ('⛔ 失败' if sr.get('error') else '⚠ 有告警')
        L.append(f"| {STAGE_META[k]['code']} | {STAGE_META[k]['name']} | {st} |")
    L.append('')

    ck = (stages or {}).get('clock') or {}
    if ck.get('stats'):
        L.append('## 二、定年结果（**每条都带方法名，不可互相比较**）')
        L.append('')
        # ⚠️ 两个 "R²" **不是同一个量**，必须各自标注（2026-10-03 修复）：
        #    · 本表列 = 每条定年方法的**时钟拟合 R²**（TreeTime molecular_clock.txt，
        #      TreeTime 口径）；
        #    · 下面 DRT 行的 real_r2 = **DRT 最大根位 R²**（逐候选根位取最大，
        #      与时钟拟合 R² 不可直接比较，数值差异主要来自"根位搜索"）。
        L.append('| 方法 | 钟模型 | 速率（每位点每年） | R²（时钟拟合, TreeTime 口径） | tMRCA | CoV | 备注 |')
        L.append('|---|---|---|---|---|---|---|')
        for s in ck['stats']:
            L.append('| %s | %s | %s | %s | %s | %s | %s |' % (
                s.get('method'), s.get('model'), _fnum(s.get('rate'), 5),
                _fnum(s.get('r2'), 4),
                _fnum(s.get('tmrca_calendar'), 6) or _fnum(s.get('tmrca_ybp'), 6),
                _fnum(s.get('cov'), 4), s.get('note') or ''))
        L.append('')
        L.append(f"时间树接口产物来源：**{ck.get('timetree_source')}**"
                 '（A2/A3 用的就是这一棵）。')
        L.append('')
        if ck.get('drt'):
            d = ck['drt']
            L.append(f"时间信号检验（DRT，{d.get('n_valid')} 次有效置换）："
                     f"R²（DRT 最大根位）={_fnum(d.get('real_r2'), 4)}，"
                     f"经验 p={_fnum(d.get('p_value'), 4)}，"
                     f"判据 **{'通过' if d.get('passed') else '不通过'}**。"
                     '（此 R² 与上表"时钟拟合"列的 R² 不是同一个量，'
                     '不可互相比较。）')
            L.append('')
            if not d.get('passed'):
                L.append('> ⚠️ **DRT 未通过 → 上面的速率与 tMRCA 只是记录，'
                         '不得作为定年结论引用**（红线 5）。')
                L.append('')

    pg_ = (stages or {}).get('phylogeo') or {}
    if pg_.get('ml_counts') or pg_.get('fitch_counts'):
        L.append('## 三、地理迁移（两条口径**并列**）')
        L.append('')
        L.append(f"- TreeTime mugration（**ML 点估计**）：迁移 **{pg_.get('ml_n_changes')}** 次，"
                 f"走廊 {pg_.get('n_corridors')} 条，事件 {pg_.get('n_events_ml')} 条")
        L.append(f"- 平台 Fitch 简约（同树）：迁移 **{pg_.get('fitch_n_changes')}** 次")
        L.append('')
        L.append('> ⚠️ 两列**不是同一个统计量**，只能看趋势一致性：**不可相加、'
                 '不可平均、不可排名**。实测 RSV 209 总数 70 vs 71（差 1.4%），'
                 '但逐走廊方向差别很大（如 NCS→KOR 24 vs 2）—— 方向性结论务必谨慎。')
        L.append('')
        L.append('> ⚠️ mugration 的 `confidence.csv` 是 **ML 边际概率，不是贝叶斯后验**'
                 '（无先验、无拓扑不确定性），与 BEAST 后验概率**不可比大小**。')
        L.append('')
        cc = pg_.get('cross_check')
        if cc:
            L.append(f"- 枝长口径对拍（换 `{cc.get('tree')}` 再跑一次 mugration）："
                     f"迁移 {cc.get('n_changes')} 次，"
                     f"相对主口径 {100 * (cc.get('rel_diff') or 0):.1f}%"
                     f"{'，**分歧显著**' if cc.get('divergent') else ''}")
            L.append('')
        sb = pg_.get('sampling_bias') or {}
        if sb:
            L.append(f"- 采样偏倚：最多区划占 {sb.get('max_share')}，"
                     f"最少占 {sb.get('min_share')}；"
                     + ('已按系数 %s 校正' % sb.get('factor')
                        if sb.get('corrected') else '**未校正**'))
            L.append('')

    co = (stages or {}).get('coalescent') or {}
    if co.get('skyline'):
        L.append('## 四、群体动态（天际线，**免 MCMC 预览**）')
        L.append('')
        L.append(f"- 网格点 {len(co['skyline'])} 个；本次（带 skyline 的）速率 "
                 f"{_fnum(co.get('rate_this_pass'), 5)}")
        L.append(f"- 区间口径：±{co.get('skyline_meta', {}).get('sd_multiple')} "
                 f"SD of the LH（**近似**）；N_e 标度假设 "
                 f"{co.get('skyline_meta', {}).get('gen_per_year')} gen/year")
        L.append('')
        L.append('> ⚠️ 这**不是**贝叶斯 95% HPD，也不是 BSP。正式 95% 带请交 BEAST（轨道 B）。')
        L.append('')

    an = (stages or {}).get('ancestral') or {}
    if an.get('ok'):
        L.append('## 五、祖先序列')
        L.append('')
        L.append(f"- 逐枝替换 **{an.get('n_mutations')}** 处；"
                 f"树用的是 `{an.get('tree_used')}`")
        L.append('')

    if warns:
        L.append('## 六、告警（逐条看，别跳过）')
        L.append('')
        for w in warns:
            L.append(f"- **[{STAGE_META.get(w['stage'], {}).get('code', w['stage'])}]** "
                     f"{w['text']}")
        L.append('')
    if errs:
        L.append('## 七、错误')
        L.append('')
        for e in errs:
            L.append(f"- **[{STAGE_META.get(e['stage'], {}).get('code', e['stage'])}]** "
                     f"{e['text']}")
        L.append('')

    L.append('## 八、产物清单（含校验和）')
    L.append('')
    L.append('| 阶段 | 产物 | 大小 | sha256(16) |')
    L.append('|---|---|---|---|')
    for f in files:
        L.append(f"| {STAGE_META.get(f['stage'], {}).get('code', f['stage'])} | "
                 f"`{f['rel']}` | {f['size']} | `{f['sha256_16']}` |")
    L.append('')
    if params:
        L.append('## 九、实际生效的参数')
        L.append('')
        L.append('```json')
        L.append(json.dumps(params, ensure_ascii=False, indent=1))
        L.append('```')
        L.append('')
    return '\n'.join(L)


# ================================================================ 编排

def run_pipeline(aln, out_dir, meta=None, tree=None, stages=None, **kw):
    """跑轨道 A 全链（或指定阶段）。返回 `{'ok','stages':{…},'summary':…}`。

    阶段依赖：`clock` 需要 `prep`；`phylogeo`/`coalescent`/`ancestral` 需要 `clock`
    （其中 `ancestral` 只用 aln + 树）。缺依赖的阶段会被**跳过并记原因**，
    不静默。

    `prep_in` / `clock_in`：**续跑**用。给上一轮已经落盘的阶段结果字典（至少含
    `outputs`），本阶段就不再重算，直接用它的接口产物往下走。阶段间只认文件，
    所以续跑 = 把文件重新指一遍，不做任何反序列化猜测。

    `should_cancel`：可调用，返回 True 就中止。**只在阶段之间检查** ——
    单个阶段内部（比如一次 TreeTime 定年要跑 30 s）不打断，因为那些子进程
    各自有 timeout，硬中断只会留下半截产物、更难排查。
    """
    log = kw.pop('log', None)
    progress = kw.pop('progress', None)
    prep_in = kw.pop('prep_in', None)
    clock_in = kw.pop('clock_in', None)
    should_cancel = kw.pop('should_cancel', None)
    rename_map = kw.pop('rename_map', None)

    def _abort_if_cancelled():
        if should_cancel is not None and should_cancel():
            raise RuntimeError('已取消（阶段之间检查）')

    say, prog = _say(log), _prog(progress)
    stages = tuple(stages or STAGE_ORDER)
    bad = [s for s in stages if s not in STAGE_ORDER]
    if bad:
        raise ValueError('未知阶段：%s' % '、'.join(bad))
    os.makedirs(out_dir, exist_ok=True)
    res = {'ok': False, 'out_dir': out_dir, 'stages': {}, 'error': None,
           'aln': aln, 'meta': meta, 'tree': tree}

    # ---- 可选前置: 序列 ID 重命名 (A0 之前; fasta/元数据/树三处同改, 全链一致) ----
    if rename_map and 'prep' in stages:
        say('可选前置: 应用序列 ID 重命名映射')
        aln, meta, tree, _rn_stats = apply_rename_map(
            aln, meta, tree, rename_map, os.path.join(out_dir, 'prep'), log=say)
        res['rename'] = _rn_stats
    elif rename_map:
        say('rename_map 已给但本轮不含 prep 阶段 → 跳过（映射应在首轮 prep 时应用）')

    # ---- A0 ----
    prep = None
    if 'prep' in stages:
        _abort_if_cancelled()
        prog('prep', 0.05, 'A0 数据准备')
        prep = stage_prep(aln, meta, tree, os.path.join(out_dir, 'prep'),
                          date_col=kw.get('date_col'),
                          trait=kw.get('trait') or 'region',
                          date_trait=kw.get('date_trait'),
                          log=say, progress=prog,
                          build_tree=kw.get('build_tree', True),
                          prune=bool(kw.get('prune')),
                          prune_method=kw.get('prune_method') or 'fps',
                          prune_date_resolution=(kw.get('prune_date_resolution')
                                                 or 'month'),
                          # ⚠️ 阈值用显式 None 判空、不用 `or 默认值`：
                          # 「0」是合法阈值（0 表示不做最小差异约束），
                          # `or` 会把它悄悄换成默认值 2/5。
                          prune_max_reps=int(
                              kw.get('prune_max_reps')
                              if kw.get('prune_max_reps') is not None else 3),
                          prune_min_snp_diff=int(
                              kw.get('prune_min_snp_diff')
                              if kw.get('prune_min_snp_diff') is not None else 2),
                          prune_clade_cutoff=int(
                              kw.get('prune_clade_cutoff')
                              if kw.get('prune_clade_cutoff') is not None else 5))
        res['stages']['prep'] = prep
        if not prep['ok']:
            res['error'] = 'A0 失败：' + str(prep.get('error'))
            say('[ERROR] ' + res['error'])
            return _finish(res, out_dir, stages, kw, say, prog)
    elif prep_in:
        prep = prep_in
        say('A0 跳过（续跑：复用上一轮 prep 产物）')
        res['stages']['prep'] = {'stage': 'prep', 'ok': True, 'error': None,
                                 'resumed': True, 'warnings': [],
                                 'outputs': prep_in.get('outputs') or {},
                                 'seq_len_modal': prep_in.get('seq_len_modal')}
        if kw.get('seq_len') is None:
            kw['seq_len'] = prep_in.get('seq_len_modal')
    else:
        prep = {'outputs': {'prep_fasta': aln,
                            'prep_nwk': tree,
                            'dates': kw.get('dates'),
                            'states': kw.get('states')},
                'seq_len_modal': kw.get('seq_len'), 'ok': True}

    # ---- A1 ----
    clock = None
    if 'clock' in stages:
        _abort_if_cancelled()
        prog('clock', 0.1, 'A1 时间推断')
        clock = stage_clock(
            prep, os.path.join(out_dir, 'clock'),
            methods=kw.get('methods') or CLOCK_METHODS,
            seq_len=kw.get('seq_len'), relax=kw.get('relax'),
            reroot=kw.get('reroot') or 'least-squares',
            keep_root=bool(kw.get('keep_root')),
            confidence=kw.get('confidence', True),
            time_marginal=kw.get('time_marginal') or 'only-final',
            clock_filter=kw.get('clock_filter'),
            clock_rate=kw.get('clock_rate'),
            drt_perm=kw.get('drt_perm', 20), drt_seed=kw.get('drt_seed', 42),
            do_drt=kw.get('do_drt', True),
            rng_seed=kw.get('rng_seed', 0),
            timeout=kw.get('timeout', 3600), log=say, progress=prog)
        res['stages']['clock'] = clock
        if not clock['ok']:
            res['error'] = 'A1 失败：' + str(clock.get('error'))
            say('[ERROR] ' + res['error'])
            return _finish(res, out_dir, stages, kw, say, prog)
    elif clock_in:
        clock = dict(clock_in)
        clock.setdefault('outputs', {})
        clock.setdefault('stats', [])
        say('A1 跳过（续跑：复用上一轮 clock 产物）')
        # 阶段字典照抄上一轮（含 timetree_source / root_year / drt / rtt_fit），
        # 告警也带过来但打上「上一轮」前缀 —— 上一轮标过的数值事故在这一轮
        # 同样是事实，抹掉才是隐瞒。
        resumed = dict(clock)
        resumed.update({
            'stage': 'clock', 'ok': True, 'error': None, 'resumed': True,
            'warnings': ['[上一轮 A1] ' + str(w)
                         for w in (clock.get('warnings') or [])]})
        res['stages']['clock'] = resumed
    elif any(s in stages for s in ('phylogeo', 'coalescent', 'ancestral')):
        clock = {'outputs': {'timetree_nwk': tree,
                             'divergence_nwk': kw.get('divergence_tree')},
                 'stats': [], 'ok': True, 'warnings': []}

    # ---- A2 ----
    if 'phylogeo' in stages:
        _abort_if_cancelled()
        prog('phylogeo', 0.1, 'A2 地理推断')
        pg_res = stage_phylogeo(
            clock or {}, prep, os.path.join(out_dir, 'phylogeo'),
            trait=kw.get('trait') or 'region',
            sampling_bias_correction=kw.get('sampling_bias_correction', 'auto'),
            pc=kw.get('pc'), cross_check=kw.get('cross_check', True),
            rssp_bs=kw.get('rssp_bs', 0), seed=kw.get('seed', 0),
            confidence=kw.get('confidence', True),
            coord_table=kw.get('coord_table'),
            timeout=kw.get('timeout', 3600), log=say, progress=prog)
        res['stages']['phylogeo'] = pg_res
        if not pg_res['ok']:
            res['error'] = res['error'] or ('A2 失败：' + str(pg_res.get('error')))
            say('[ERROR] A2：' + str(pg_res.get('error')))

    # ---- A3 ----
    if 'coalescent' in stages:
        _abort_if_cancelled()
        prog('coalescent', 0.1, 'A3 群体动态')
        co = stage_coalescent(
            clock or {}, prep, os.path.join(out_dir, 'coalescent'),
            n_skyline=kw.get('n_skyline'),
            gen_per_year=kw.get('gen_per_year', 50.0),
            rng_seed=kw.get('rng_seed', 0),
            timeout=kw.get('timeout', 3600), log=say, progress=prog)
        res['stages']['coalescent'] = co
        if not co['ok']:
            res['error'] = res['error'] or ('A3 失败：' + str(co.get('error')))
            say('[ERROR] A3：' + str(co.get('error')))

    # ---- A4 ----
    if 'ancestral' in stages:
        _abort_if_cancelled()
        prog('ancestral', 0.1, 'A4 祖先序列')
        an = stage_ancestral(
            clock or {}, prep, os.path.join(out_dir, 'ancestral'),
            marginal=bool(kw.get('marginal')),
            method_anc=kw.get('method_anc'),
            rng_seed=kw.get('rng_seed', 0),
            timeout=kw.get('timeout', 3600), log=say, progress=prog)
        res['stages']['ancestral'] = an
        if not an['ok']:
            res['error'] = res['error'] or ('A4 失败：' + str(an.get('error')))
            say('[ERROR] A4：' + str(an.get('error')))

    _abort_if_cancelled()
    return _finish(res, out_dir, stages, kw, say, prog)


def _finish(res, out_dir, stages, kw, say, prog):
    if 'report' in stages:
        rep = stage_report(res['stages'], os.path.join(out_dir, 'report'),
                           params={k: v for k, v in kw.items()
                                   if isinstance(v, (str, int, float, bool,
                                                     list, tuple, type(None)))},
                           log=say, progress=prog)
        res['stages']['report'] = rep
    # 顶层 ok：只要有阶段跑成功且没有 A0/A1 级致命错误
    fatal = [k for k, v in res['stages'].items()
             if k in ('prep', 'clock') and isinstance(v, dict) and not v.get('ok')]
    res['ok'] = not fatal and any(
        isinstance(v, dict) and v.get('ok') for v in res['stages'].values())
    return res
