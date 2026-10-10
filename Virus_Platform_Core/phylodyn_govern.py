# -*- coding: utf-8 -*-
"""元数据入口治理（时间 + 地理）：归一 / 标注 / 体检 —— **不猜、不静默、不改语义**。

来源与定位
----------
本模块的口径**移植自上游** `MMPV-RNA/virome_phylo_pipeline/utils/metadata_governance.py`
（2026-09-16 建立，作者列的动机是管线入口踩过四类"静默出错"）。上游那四类，
我们平台上**同样会踩**（我们的 `normalize_meta_rows` 只做"别名识别 + 择优"，
不做占位符识别、不认 `DD-Mon-YYYY`、把"只剩国家名"当地点用）：

  A. **坐标占位符**  `XX.XX N XXX.XX E` —— 填了个"形似坐标"的模板，看着有值。
  B. **日期占位符**  字面量 `YYYY-MM-DD` / `YYYY-MM` / `YYYY` —— 形似合法日期，
     **任何格式校验都会放行**（我们的 `validate_date_str` 就会放行）。
  C. **地理层级 Unknown**  `China, Unknown, Yinchuan_AI` —— 若静默降级成"国家质心"，
     误差可达上千公里，却仍被标成"已解析"。
  D. **只给年 / GenBank 日期**  `YYYY` 与 `DD-Mon-YYYY` —— 后者我们的
     `validate_date_str` 会**直接判不合格并剔行**。

三条原则（与平台既有口径一致，此处写死）
---------------------------------------
  1. **不猜** —— 解析不出就标 `issue` + 置空，绝不套默认值；
  2. **不静默** —— 任何降精度 / 任何丢弃都进报告，可统计条数；
  3. **不改语义** —— 只归一 + 标注，**默认不丢行**（清空该字段而非删行），
     要丢行必须由调用方**显式**开开关。

⚠️ 小数年口径（本模块唯一权威，见 `decimal_year`）
------------------------------------------------
我们平台此前有**四处**各自换算，同一份数据能算出不同答案（实测）：

    | 数据        | phylodyn_local | phylodyn_kit | virphykit_export | 上游 decimal_year |
    |---|---|---|---|---|
    | `2013-05-01`| 2013.375       | 2013.3313    | 2013.3288        | 2013.3333 |
    | `2013-05`   | 2013.375       | 2013.3696    | 2013.3288        | 2013.3709 |
    | `2013`      | **2013.0**     | **2013.5366**| **2013.0**       | 2013.4556 |

差异最大的是「只给年」：**2013.0（年初） vs 2013.5366（年中）**，差 **0.54 年** ——
对以年为尺度的定年/时间信号而言这是"同一份数据两个答案"。
本模块给出唯一实现，并按上游口径把 `mode` **显式暴露**（`'mid'` 年中 / `'start'` 年初），
默认 `'mid'`（年只给到时取 6 月 15 日）—— 这正是上游 `decimal_year.py` 当初的立项目的。
"""

import calendar
import re
from collections import Counter, OrderedDict
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

# ══════════════════════════════════════════════════════════════════════
# 占位符
# ══════════════════════════════════════════════════════════════════════

PLACEHOLDER_LITERALS = {
    "", "na", "n/a", "nan", "none", "null", "missing", "unknown",
    "not applicable", "not_applicable", "not provided", "not_provided",
    "tbd", "unspecified", "未", "未知", "待定", "不详",
}

# ⚠️ 只精确匹配**已知**占位字面量，**不许**用泛化正则：
#    上游的教训 —— 曾写 `^[a-z]+\s*[:：]\s*[a-z]+$` 想兜住"任意 Word:Word"，
#    结果把真实地理 `China:Ningxia` / `China:Beijing` 全打成占位符（72 条）。
PLACEHOLDER_PATTERNS = [
    re.compile(r"^y{4}[-/.]m{2}[-/.]d{2}$", re.I),                       # YYYY-MM-DD
    re.compile(r"^y{4}[-/.]m{2}$", re.I),                                # YYYY-MM
    re.compile(r"^y{4}$", re.I),                                         # YYYY
    re.compile(r"^x+(\.[x]+)?\s*[ns]\s*x+(\.[x]+)?\s*[ew]$", re.I),     # XX.XX N XXX.XX E
    re.compile(r"^country\s*[:：]\s*region$", re.I),                     # Country:Region
]


def is_placeholder(v) -> bool:
    """是不是「未填值」占位符（含**形似格式的假日期/假坐标**）。"""
    s = str(v or "").strip()
    if s.lower() in PLACEHOLDER_LITERALS:
        return True
    return any(p.match(s) for p in PLACEHOLDER_PATTERNS)


# ══════════════════════════════════════════════════════════════════════
# 日期
# ══════════════════════════════════════════════════════════════════════

_MONTHS = {'jan': 1, 'feb': 2, 'mar': 3, 'apr': 4, 'may': 5, 'jun': 6,
           'jul': 7, 'aug': 8, 'sep': 9, 'sept': 9, 'oct': 10, 'nov': 11,
           'dec': 12}

# 认了哪些格式、归一后精度是什么（`DateInfo.fmt` 的取值表，报告 `date_format`
# 直方图里的键就出自这里 —— 逐个对应着读报告）。
# ⚠️ 2026-10-03：原为模块常量 `DATE_FORMAT_TABLE`，但全仓无代码引用（只被
#    注释/文档提到）→ 降级成注释表：文档信息保留，不再是一个没人用的死常量。
#     iso_day      YYYY-MM-DD                                day
#     iso_month    YYYY-MM                                  month
#     iso_year     YYYY                                     year
#     genbank      DD-Mon-YYYY                              day
#     euro_day     DD/MM/YYYY 或 DD.MM.YYYY                 day
#     mon_year     Mon-YYYY                                 month
#     decimal_year 已转好的小数年（小数位 ≥4，如 2007.569473） day


@dataclass
class DateInfo:
    """日期解析结果（只承载事实，不做取舍）。"""
    raw: str = ''
    iso: Optional[str] = None            # 归一后的 ISO（YYYY / YYYY-MM / YYYY-MM-DD）
    decimal: Optional[float] = None      # 小数年（按 mode）
    precision: str = 'none'              # day / month / year / none
    fmt: str = ''                        # 格式名（取值见上方注释表）
    issue: str = ''                      # placeholder / unparsable / out_of_range / ''
    ambiguous_frac: bool = False         # 形如 2019.5 / 2014.11 的"小数年 vs 年.月"二义写法
    mode: str = 'mid'

    @property
    def ok(self) -> bool:
        return self.iso is not None


def _days_in_month(y: int, mo: int) -> int:
    return calendar.monthrange(y, mo)[1]


def decimal_year(raw, *, mode: str = 'mid') -> Optional[float]:
    """日期 → 小数年（**唯一权威实现**）。认不出返回 None（绝不猜）。

    约定（与上游 `utils/decimal_year.py` 一致，也与 BEAST 侧的"当月实际天数"一致）::

        YYYY-MM-DD → yr + (mo - 1 + (dy - 1) / days_in_month) / 12
        YYYY-MM    → 该月 15 日（mode='mid'）/ 该月 1 日（mode='start'）
        YYYY       → 6 月 15 日（mode='mid'）/ 1 月 1 日（mode='start'）
        小数年（小数位 ≥4）→ 原样返回

    ⚠️ `'mid'` 与 `'start'` 只在「只给年 / 只给年月」两级有差异 ——
    这是**有意的口径选择**，由调用方显式给，不藏在实现里。
    """
    info = parse_date(raw, mode=mode)
    return info.decimal


def parse_date(raw, *, mode: str = 'mid') -> DateInfo:
    """解析日期字符串 → `DateInfo`。认 7 种格式，见上方格式注释表。

    · 占位符（字面量 `YYYY-MM-DD` / `YYYY` / `XX`…）→ `issue='placeholder'`
    · 小数年与「年.月」的**消歧判据**：小数位 **≥4** 才算小数年
      （实测 VirPhyKit 数据小数位 4~6；而 `2014.05` 这种 `年.月` 写法恒为 2 位）
    · 只给年/只给年月 → 按 `mode` 取年中或年初（**这是唯一影响数值的分支**）
    """
    assert mode in ('mid', 'start'), f"mode 必须是 'mid'|'start'，收到 {mode!r}"
    info = DateInfo(raw=str(raw or '').strip(), mode=mode)
    if not info.raw or is_placeholder(info.raw):
        # 空串也归 placeholder（与上游一致：少一个 issue 取值，报告更好读）
        info.issue = 'placeholder'
        return info
    s = info.raw.replace('/', '-').strip()

    # ① 小数年 vs「年.月」——消歧判据（**比上游多一条，因为我们的数据源不同**）
    #
    #   上游只按"小数位 ≥4 才算小数年"，那是因为他们面对的是 BioAider GUI 的
    #   `2014.05`（年.月）。但**本平台的数据源是 explorer 导出的 date 列**，
    #   那里的小数年是 `y + doy/days_in_year`，位数不定 —— `2019.5` 这种 1 位小数
    #   也是小数年（平台既有 `decimal_to_iso('2019.5') → 2019-07-02`，并有测试钉住）。
    #   若照抄上游规则，`2019.5` 会被读成"5 月"，时间偏半年。
    #
    #   实际判据：**带前导零的两位**（`.05` / `.11`）＝ 年.月；
    #   其余（`.5` / `.75` / `.569473`）＝ 小数年。
    #   ⚠️ `.11` 这种"既是 11 月、也可能是 0.11 年"的写法**本质二义** →
    #   记为 `ambiguous_frac` 计数并进报告，不静默选一个。
    m_frac = re.fullmatch(r'(\d{4})\.(\d+)', s)
    if m_frac:
        frac = m_frac.group(2)
        is_month_style = bool(re.fullmatch(r'0[1-9]|1[0-2]', frac))
        info.ambiguous_frac = (not is_month_style and len(frac) <= 2
                               and 1 <= int(frac) <= 12)
        if not is_month_style:
            y = int(m_frac.group(1))
            if not 1900 <= y <= 2100:
                info.issue = 'out_of_range'
                return info
            try:
                v = float(s)
            except ValueError:
                info.issue = 'unparsable'
                return info
            doy = max(1, min(366, round((v - y) * 365) + 1))
            import datetime as _dt
            info.iso = (_dt.date(y, 1, 1) + _dt.timedelta(days=doy - 1)).isoformat()
            info.decimal = v
            info.precision, info.fmt = 'day', 'decimal_year'
            return info
        # 带前导零的两位小数（`.05` / `.11`）＝「年.月」→ 改写成 ISO 写法继续走下面的分支
        s = '%s-%s' % (m_frac.group(1), frac)

    # ③ ISO 三档
    m = re.fullmatch(r'(\d{4})(?:-(\d{1,2})(?:-(\d{1,2}))?)?', s)
    if m:
        y = int(m.group(1))
        if not 1900 <= y <= 2100:
            info.issue = 'out_of_range'
            return info
        mo = int(m.group(2)) if m.group(2) else None
        dy = int(m.group(3)) if m.group(3) else None
        return _finish_iso(info, y, mo, dy, 'iso')

    # ④ DD-Mon-YYYY / Mon-YYYY（GenBank 口径）
    m = re.fullmatch(r'(\d{1,2})-([A-Za-z]{3,4})-(\d{4})', s)
    if m and m.group(2).lower() in _MONTHS:
        return _finish_iso(info, int(m.group(3)), _MONTHS[m.group(2).lower()],
                           int(m.group(1)), 'genbank')
    m = re.fullmatch(r'([A-Za-z]{3,4})-(\d{4})', s)
    if m and m.group(1).lower() in _MONTHS:
        info0 = _finish_iso(info, int(m.group(2)), _MONTHS[m.group(1).lower()],
                            None, 'mon_year')
        info0.fmt = 'mon_year'
        return info0

    # ⑤ DD/MM/YYYY 或 DD.MM.YYYY（欧洲写法；`/` 已在上面换成 `-` 吗？没有——
    #    这里用原始串判断，避免与 ISO 的 `-` 混淆）
    raw2 = info.raw
    m = re.fullmatch(r'(\d{1,2})[/.](\d{1,2})[/.](\d{4})', raw2)
    if m:
        d, mo, y = int(m.group(1)), int(m.group(2)), int(m.group(3))
        # ⚠️ **先认形状再验范围**：`13/13/2020` 要判「越界」而不是「认不出」——
        #    诊断信息差一档（上游同口径）。范围校验在 `_finish_iso` 里做。
        return _finish_iso(info, y, mo, d, 'euro_day')

    info.issue = 'unparsable'
    return info


def _finish_iso(info: 'DateInfo', y: int, mo: Optional[int], dy: Optional[int],
                fmt: str) -> 'DateInfo':
    """把 (y, mo, dy) 归一成 ISO + 小数年；`mode` 只在这一步起作用。

    ⚠️ **越界一律判 `out_of_range`，绝不夹取**（与上游一致）。
    曾写成 `min(max(mo,1),12)` 把 `2020-13-01` 静默夹成 `2020-12-01` —— 那样
    "13 月"这种明确错值会被当成好数据放行，而平台上原有的
    `validate_date_str` 与 UI 测试都要求它被点名剔除。
    """
    if mo is not None and not 1 <= mo <= 12:
        info.issue = 'out_of_range'
        return info
    if mo is not None and dy is not None and not 1 <= dy <= _days_in_month(y, mo):
        info.issue = 'out_of_range'
        return info
    if mo is None:
        iso = f'{y:04d}'
        prec = 'year'
        mo_eff = 6 if info.mode == 'mid' else 1
        dy_eff = 15 if info.mode == 'mid' else 1
    elif dy is None:
        iso = f'{y:04d}-{mo:02d}'
        prec = 'month'
        mo_eff = mo
        dy_eff = 15 if info.mode == 'mid' else 1
    else:
        iso = f'{y:04d}-{mo:02d}-{dy:02d}'
        prec = 'day'
        mo_eff, dy_eff = mo, dy
    info.iso = iso
    info.precision = prec
    info.fmt = fmt if prec != 'year' or fmt != 'iso' else 'iso_year'
    if info.fmt == 'iso' and prec == 'month':
        info.fmt = 'iso_month'
    if info.fmt == 'iso' and prec == 'year':
        info.fmt = 'iso_year'
    dim = _days_in_month(y, mo_eff)
    info.decimal = y + (mo_eff - 1.0 + (min(dy_eff, dim) - 1.0) / dim) / 12.0
    return info


# ══════════════════════════════════════════════════════════════════════
# 地理
# ══════════════════════════════════════════════════════════════════════

# "没有信息"的层级：删掉它们、但**留痕**（不是当地点用）
UNKNOWN_LEVEL_TOKENS = {'unknown', 'unk', 'n/a', 'na', 'none', 'null', 'unspecified',
                        'not applicable', 'not_provided', '未知', '未', '不详', '待定'}

# "纯国家"判据：只剩这些词时，该条已无分辨率 —— **不得**赋国家中心点
COUNTRY_ONLY_TOKENS = {
    'china', 'cn', 'russia', 'united states', 'usa', 'us', 'japan', 'india',
    'brazil', 'australia', 'canada', 'france', 'germany', 'italy', 'spain',
    'korea', 'south korea', 'republic of korea', 'mexico', 'argentina', 'chile',
    'south africa', 'egypt', 'iran', 'iraq', 'saudi arabia', 'turkey', 'poland',
    'ukraine', 'netherlands', 'belgium', 'sweden', 'norway', 'finland', 'denmark',
    'austria', 'switzerland', 'portugal', 'greece', 'hungary', 'romania',
    'bulgaria', 'serbia', 'croatia', 'slovenia', 'slovakia', 'czech republic',
    'new zealand', 'indonesia', 'malaysia', 'thailand', 'vietnam', 'philippines',
    'singapore', 'pakistan', 'bangladesh', 'sri lanka', 'nepal', 'kenya',
    'nigeria', 'ethiopia', 'morocco', 'algeria', 'tunisia', 'ghana', 'uganda',
    'tanzania', 'colombia', 'peru', 'venezuela', 'ecuador', 'bolivia', 'uruguay',
    'paraguay', 'cuba', 'guatemala', 'honduras', 'panama', 'costa rica', 'syria',
    'jordan', 'lebanon', 'israel', 'yemen', 'oman', 'kuwait', 'qatar', 'uae',
    'united arab emirates', 'afghanistan', 'kazakhstan', 'uzbekistan',
    'azerbaijan', 'georgia', 'armenia', 'belarus', 'lithuania', 'latvia',
    'estonia', 'moldova', 'albania', 'iceland', 'ireland', 'luxembourg', 'malta',
    'cyprus', 'mongolia', 'myanmar', 'cambodia', 'laos', 'brunei', 'fiji',
    'papua new guinea', '汉城',
}
# 中国及其特别行政区的各国别写法（归属表述要一致）
_CN_ALIASES = {'china', 'cn', 'prc', "people's republic of china",
               'mainland china', '中国'}


@dataclass
class GeoInfo:
    """地理解析结果。`usable` 是**唯一**判据：能不能拿去当分析用的地点。"""
    raw: str = ''
    components: List[str] = field(default_factory=list)   # 删掉 Unknown 层级后剩下的段
    dropped_levels: List[str] = field(default_factory=list)
    normalized: Optional[str] = None                      # 'China, Yinchuan'
    issue: str = ''                                       # placeholder / unknown_level /
                                                          # country_only / ''
    level: str = 'none'                                   # country / sub / none（粗略分级）

    @property
    def usable(self) -> bool:
        """能不能用作「地点」（**唯一判据**）。

        `country_only`（删掉 Unknown 层级后**只剩国名**）**不算** —— 那种情况
        误差可达上千公里，宁可算"没有地点"，也不给一个假精度很高的点
        （上游 2026-09-16 老师拍板）。粒度信息另放 `level`（country / sub），
        给报告统计用，不参与可用性判断。
        """
        return bool(self.normalized) and self.issue != 'country_only'


def split_location(raw) -> List[str]:
    """拆地理层级：**冒号与逗号都拆**。

    `China:Ningxia` → `['China','Ningxia']`；`China, Ningxia, Yinchuan_AI` → 三段。
    （上游实测：只拆 `[,;]` 会让 `China:Ningxia` 整串去查表、查不到就静默降级成国家质心。）
    """
    return [p.strip() for p in re.split(r'[,;:：]', str(raw or '')) if p.strip()]


def _is_unknown_level(p: str) -> bool:
    """是不是"没有信息"的层级。

    ⚠️ 不只看精确词表：`Unknown_AI` / `unknown_region` 这类**带后缀**的写法同样是
    占位（上游会把它们一起删掉）。判据 = 词表命中 **或** 以 unknown/未/不详 开头。
    """
    t = (p or '').strip().lower()
    if t in UNKNOWN_LEVEL_TOKENS:
        return True
    return bool(re.match(r'^(unknown|unk|unspecified|未|不详)([_\- ]|$)', t))


def _is_country_only_term(p: str) -> bool:
    return (p or '').strip().lower() in COUNTRY_ONLY_TOKENS


def parse_location(raw, *, drop_unknown: bool = True) -> GeoInfo:
    """解析地理字符串。

    `drop_unknown=True`（默认，上游口径）：**删除** `Unknown` 层级并在
    `dropped_levels` 留痕 —— 这样 `China, Unknown, Yinchuan_AI` → `China, Yinchuan`，
    而不是像"静默降级"那样返回一个国家质心却标记为"已解析"。

    删完只剩纯国家名时（`China, Unknown, Unknown_AI`）→ `issue='country_only'`、
    `normalized=None`、`usable=False`：该记录**没有地理分辨率**，
    由调用方决定是清空还是剔行（默认清空，不丢行）。
    """
    info = GeoInfo(raw=str(raw or '').strip())
    if not info.raw or is_placeholder(info.raw):
        info.issue = 'placeholder'
        return info

    kept, dropped = [], []
    for p in split_location(info.raw):
        if _is_unknown_level(p):
            (dropped if drop_unknown else kept).append(p)
        else:
            kept.append(p)
    info.components = kept
    info.dropped_levels = dropped
    if dropped:
        info.issue = 'unknown_level'

    if not kept:
        # 全被删光（如 `Unknown:Unknown`）→ 保留 `unknown_level` 这个**原因**，
        # 不要改写成 placeholder（那是"值本身是占位符"，不是同一回事）
        info.issue = info.issue or 'unknown_level'
        info.normalized = None
        return info

    info.normalized = ', '.join(kept)
    # 只剩"纯国家"**且删过层级** → 该条已无地理分辨率，**不得**赋国家中心点
    # （误差可达上千公里）。⚠️ 光秃秃一个 `Russia` **不**判 country_only ——
    # 它就是个（粗粒度但真实）的地点；我们的数据里 438 条都是这种，
    # 一律判死会把真地点全清掉。粒度另用 `level` 记录，供报告统计。
    if dropped and all(_is_country_only_term(p) for p in kept):
        info.issue = 'country_only'
        info.level = 'country'
        info.normalized = None      # 与上游一致：没分辨率的**不给** normalized
        return info                 # （`components` 里仍留着 'China' 便于回溯）
    info.level = 'country' if all(_is_country_only_term(p) for p in kept) else 'sub'
    return info


# ══════════════════════════════════════════════════════════════════════
# 表级治理
# ══════════════════════════════════════════════════════════════════════

def govern_rows(rows, *, date_key: str = 'date', loc_key: str = 'location',
                name_key: str = 'name', mode: str = 'mid',
                drop_unknown_geo: bool = True,
                blank_placeholders: bool = True) -> Tuple[List[Dict], Dict]:
    """对每行做「归一 + 标注」，返回 `(rows2, report)`。

    处置口径（默认，与上游一致）：

    · **只清空字段，不丢行** —— 占位符/不可用地点的格子清空，行保留，
      `row['_gov']` 里留痕（`date_issue` / `loc_issue` / `precision` / `dropped_levels`）。
    · `normalize_meta_rows` 已经做过的「别名识别 + 择优」在**之前**跑；
      本函数只处理**选出来的那个值**干不干净。
    · `report` 是给人看的账：每种 issue 多少条、精度分布、被删的 Unknown 层级、
      `country_only` 多少条（这些"看着有地点、其实没分辨率"的必须点出来）。

    ⚠️ 本函数**不改** 已有 `lat`/`lon`；坐标的补算/分级由 `attach_coords` 与
    调用方负责（补算约值的 `source` 必须与真实观测分开记，见 `phylodyn_kit.attach_coords`）。
    """
    out: List[Dict] = []
    date_issue, loc_issue, prec = Counter(), Counter(), Counter()
    fmt_hist, dropped_hist = Counter(), Counter()
    n_country_only = 0
    n_unknown_level_all = 0    # 层级全是 Unknown（清空该字段）的行数
    n_ambig = 0          # 形如 `2019.5`/`2014.11` 的"小数年 vs 年.月"二义写法计数
    for r in rows:
        r2 = dict(r)
        gov = {}

        di = parse_date(r2.get(date_key), mode=mode)
        if di.issue:
            date_issue[di.issue] += 1
            gov['date_issue'] = di.issue
            # ⚠️ **只清"占位符"**（字面量 `YYYY-MM-DD` 之类本来就是"没填"）；
            #    真错值（`2020-13-01` 越界、认不出的串）**原样留着**，
            #    让下游校验点名"日期不合规：13 月超范围" —— 清掉会退化成
            #    "缺少日期"，用户就看不到到底哪错了。
            if blank_placeholders and di.issue == 'placeholder':
                r2[date_key] = ''
        else:
            r2[date_key] = di.iso            # 归一成 ISO（下游只认 ISO）
            gov['date_precision'] = di.precision
            gov['date_fmt'] = di.fmt
            if di.ambiguous_frac:
                n_ambig += 1
                gov['date_ambiguous'] = True
            prec[di.precision] += 1
            fmt_hist[di.fmt] += 1

        gi = parse_location(r2.get(loc_key), drop_unknown=drop_unknown_geo)
        if gi.issue:
            loc_issue[gi.issue] += 1
            gov['loc_issue'] = gi.issue
        if gi.dropped_levels:
            dropped_hist['、'.join(gi.dropped_levels)] += 1
            gov['dropped_levels'] = gi.dropped_levels
        if gi.issue == 'country_only':
            n_country_only += 1
        if gi.issue == 'unknown_level' and not gi.usable:
            # 层级**全是** Unknown（`Unknown:Unknown` / `unknown_region`）：删完
            # 什么都不剩，等同于"没有地点"
            n_unknown_level_all += 1
        if gi.usable:
            r2[loc_key] = gi.normalized        # 归一：删掉 Unknown 层级
            gov['loc_level'] = gi.level
        elif blank_placeholders:
            # ⚠️ 清空名单补 `unknown_level`（2026-10-03 修）：`Unknown:Unknown`
            #    这类全 Unknown 层级的地点 normalized=None、usable=False，但原实现
            #    只清 placeholder/country_only —— 字面量 "Unknown" 会原样留下，
            #    继续参与等量降采样 / 迁移统计（等于把一个不存在的区划当数据用）。
            if gi.issue in ('placeholder', 'country_only', 'unknown_level'):
                r2[loc_key] = ''               # 没分辨率 → 清空，但不丢行
        # 保底：把原始值留在旁边，便于回溯（下游只读规范列）
        if gov:
            r2['_gov'] = gov
        out.append(r2)

    report = {
        'mode': mode,
        'n_rows': len(rows),
        'date_issue': dict(date_issue),
        'location_issue': dict(loc_issue),
        'date_precision': dict(prec),
        'date_format': dict(fmt_hist.most_common(8)),
        'dropped_unknown_levels': dict(dropped_hist.most_common(8)),
        'n_country_only': n_country_only,
        'n_unknown_level_all': n_unknown_level_all,
        'n_date_ambiguous_frac': n_ambig,
        'n_date_ok': sum(prec.values()),
    }
    warns = []
    if date_issue.get('placeholder'):
        warns.append('%d 行的日期是**占位符**（字面量 YYYY-MM-DD / YYYY 之类），'
                     '形似合法日期、格式校验会放行 —— 已清空并在报告留痕'
                     % date_issue['placeholder'])
    if n_country_only:
        warns.append('%d 行清理掉 Unknown 层级后**只剩国家名**（无地理分辨率）：'
                     '已就地清空，**没有**赋国家中心点（那样误差可达上千公里）'
                     % n_country_only)
    if loc_issue.get('placeholder'):
        warns.append('%d 行的地点是占位符/Unknown —— 已**清空该字段**（行保留）：'
                     '下游会把它当"没有地点"，这比拿 Unknown 当区划去参与'
                     '时空调度/降采样可靠' % loc_issue['placeholder'])
    if n_unknown_level_all:
        warns.append('%d 行的地点层级**全是 Unknown**（如 `Unknown:Unknown`、'
                     '`unknown_region`）—— 已就地清空，不赋任何坐标、也不再'
                     '参与等量降采样/迁移统计' % n_unknown_level_all)
    if dropped_hist:
        warns.append('删除过 Unknown 地理层级的行 %d 条（例：%s）'
                     % (sum(dropped_hist.values()),
                        '、'.join(list(dropped_hist)[:3])))
    if n_ambig:
        warns.append('%d 行的日期形如 （既是 11 月、也可能是小数年 0.11）——'
                     '本平台按**小数年**处理（数据源是 explorer 的 date 列）；'
                     '若这批数据来自 BioAider 那种"年.月"写法，请改用它导出的 ISO'
                     % n_ambig)
    if prec.get('year'):
        warns.append('%d 行的日期只精确到**年**，按 %s 口径折算（'
                     'mid＝6 月 15 日 / start＝1 月 1 日）'
                     % (prec['year'], mode))
    if warns:
        report['warnings'] = warns
    return out, report


def inspect_table(path, *, date_col: Optional[str] = None,
                  loc_col: Optional[str] = None, name_col: Optional[str] = None,
                  mode: str = 'mid', **kw) -> Dict:
    """读一张 metadata 表并体检（不改文件，只出报告）。列名不给就自动认。"""
    import csv
    import io as _io
    text = _io.open(path, encoding='utf-8-sig', errors='replace').read()
    sep = '\t' if text[:4096].count('\t') > text[:4096].count(',') else ','
    rows = [{ (k or '').strip(): (v or '').strip() for k, v in r.items() if k}
            for r in csv.DictReader(_io.StringIO(text), delimiter=sep)]
    if not rows:
        return {'n_rows': 0, 'error': '空表'}
    low = {k.lower(): k for k in rows[0]}
    date_key = date_col or next((low[k] for k in
                                 ('date', 'collection_date', 'release_date',
                                  'sampling_date') if k in low), None)
    loc_key = loc_col or next((low[k] for k in
                               ('location', 'geo_location', 'country', 'region')
                               if k in low), None)
    _, rep = govern_rows(rows, date_key=date_key or 'date',
                         loc_key=loc_key or 'location',
                         name_key=name_col or 'name', mode=mode, **kw)
    rep['columns'] = {'date': date_key, 'location': loc_key}
    return rep
