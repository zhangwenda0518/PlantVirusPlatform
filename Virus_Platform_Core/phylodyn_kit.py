# -*- coding: utf-8 -*-
"""进化动力学·序列/元数据工具箱（phylodyn 组新增卡片的后端引擎）。

与 VirPhyKit（Yin et al., 2025, Ecol Evol）的**方法学对齐**、代码全部自研
（VirPhyKit 为 GPL-3.0，不复制其代码，只对齐输入/输出契约与算法口径）。
本模块负责「序列/元数据」一侧：

    1. SeqIDRenamer 同款：按 TSV 映射表批量重命名 FASTA 序列 ID
    2. SeqGrouper 同款：GenBank 记录表 + 分组映射表 → Group 列 / 分组计数
    3. GeoSubsampler 同款：FASTA 时空降采样（均等采样 / 定区剔除 / 随机抽取）
    4. VirSpaceTime 同款：Temporal.txt / Spatial.txt 解析（前端 plotly 绘图）
    5. MJRM 生成器：BEAST Markov jumps/rewards 矩阵 XML 片段生成与插入
    6. 进化动力学数据接入：三来源（explorer 导出包 / 手动输入 / 在线下载）
       + 时间地点同检 + 标准化数据集 + 汇总（时间轴 / 地点分布 / 地图点位）

树类引擎（TempMig / RRT / RSPP / BSP / TreeTime / LTT）见 phylodyn_trees.py。
网页卡接线见 web/tool_jobs.py（_tool_job_pd*）。
"""
import csv
import io
import math
import os
import random
import re
from collections import OrderedDict

# ================================================================
# 通用小件
# ================================================================

_FASTA_WS = re.compile(r'\s+')


def read_text(path):
    """容错读文本：utf-8-sig 优先，失败退 gbk。"""
    with open(path, 'rb') as f:
        raw = f.read()
    for enc in ('utf-8-sig', 'utf-8', 'gbk'):
        try:
            return raw.decode(enc)
        except UnicodeDecodeError:
            continue
    return raw.decode('utf-8', errors='replace')


def read_fasta(path):
    """FASTA → OrderedDict[头, 序列]（头保留 > 后原文，序列去空白）。"""
    out = OrderedDict()
    name, chunks = None, []
    for line in read_text(path).splitlines():
        if line.startswith('>'):
            if name is not None:
                out[name] = ''.join(chunks)
            name = line[1:].strip()
            chunks = []
        elif name is not None:
            chunks.append(_FASTA_WS.sub('', line))
    if name is not None:
        out[name] = ''.join(chunks)
    return out


def write_fasta(path, pairs, width=80):
    with open(path, 'w', encoding='utf-8', newline='\n') as f:
        for name, seq in pairs:
            f.write(f'>{name}\n')
            for i in range(0, len(seq), width):
                f.write(seq[i:i + width] + '\n')


# ---------------------------------------------------------------- 日期口径
# 平台统一口径：YYYY / YYYY-MM / YYYY-MM-DD；小数年只在树/元数据
# 明确允许的地方（TreeTime 元数据、MCC 树叶名）使用。
# ⚠️ 2026-09-18 修正：**数据接入不再拒收小数年**，而是折算成 ISO 后入数据集。
#   原因：explorer「导出 VirPhyKit 输入包」的 date 列**默认就是小数年**
#   （`virphykit_export.build_package(..., date_mode=DATE_DECIMAL)`）——
#   原口径与自家上游产物直接冲突，实测拿真导出包（209 条小数年）导入
#   会整包被拒、只抛一句「没有一条序列通过校验」，上下游根本接不上。
#   折算（而非原样透传）是为了让数据集保持**单一 ISO 口径**：ISO 是
#   explorer 侧自己标注的「两边都读得对」的写法，而小数年在自家管线
#   （`to_decimal_year(strict=True)`）里会被错读成 +0.46 年。
# ----------------------------------------------------------------
_DOY_CUM = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334]
_DEC_RE = re.compile(r'\d{4}\.\d+')


def validate_date_str(s, allow_decimal=False):
    """校验 YYYY[-MM[-DD]]。返回 (ok, 规范值或 None, 原因)。

    allow_decimal=True 时**额外接受小数年**（如 2008.58197）——那是 explorer
    「导出 VirPhyKit 输入包」的默认 date 口径；此时返回值就是该小数的原样
    字符串，由调用方决定折算还是透传（数据接入走 `decimal_to_iso`）。
    默认 False：树/其它调用方仍按严格 ISO 口径。
    """
    s = (s or '').strip()
    if allow_decimal and _DEC_RE.fullmatch(s):
        return True, s, ''
    m = re.fullmatch(r'(\d{4})(?:-(\d{2})(?:-(\d{2}))?)?', s)
    if not m:
        if _DEC_RE.fullmatch(s):
            return False, None, ('小数年（如 2007.57）不被接受，请改用 '
                                 'YYYY / YYYY-MM / YYYY-MM-DD')
        if re.fullmatch(r'\d{1,2}[/-]\d{1,2}[/-]\d{2,4}', s):
            return False, None, '斜杠/连字符日-月-年格式不被接受，请改用 YYYY-MM-DD'
        return False, None, '日期必须是 YYYY / YYYY-MM / YYYY-MM-DD'

    y = int(m.group(1))
    mo = int(m.group(2)) if m.group(2) else None
    dy = int(m.group(3)) if m.group(3) else None
    if y < 1000 or y > 3000:
        return False, None, f'年份 {y} 超出合理范围 (1000-3000)'
    if mo is not None and not (1 <= mo <= 12):
        return False, None, f'月份 {mo:02d} 不在 01-12'
    if dy is not None:
        dim = {1: 31, 2: 28, 3: 31, 4: 30, 5: 31, 6: 30,
               7: 31, 8: 31, 9: 30, 10: 31, 11: 30, 12: 31}[mo]
        leap = (y % 4 == 0 and (y % 100 != 0 or y % 400 == 0))
        if mo == 2 and leap:
            dim = 29
        if not (1 <= dy <= dim):
            return False, None, f'日 {dy:02d} 不在 {y}-{mo:02d} 的 01-{dim:02d}'
    out = str(y) + (f'-{mo:02d}' if mo is not None else '') + (f'-{dy:02d}' if dy is not None else '')
    return True, out, ''


def date_to_decimal(s):
    """YYYY[-MM[-DD]] → 小数年（TreeTime 口径：year + day_of_year/365.25，
    Jan 1 = 1；实测对齐 VirPhyKit 元数据 2007-07-27 → 2007.569473）。"""
    ok, v, err = validate_date_str(s)
    if not ok:
        # 已经是小数年 / 年-月-日以外的合法树元数据格式，交给上层决定
        try:
            return float(s)
        except (TypeError, ValueError):
            raise ValueError(f'无法解析日期 [{s}]: {err}')
    parts = [int(x) for x in v.split('-')]
    y = parts[0]
    mo = parts[1] if len(parts) > 1 else 7     # 缺省按年中
    dy = parts[2] if len(parts) > 2 else 15
    doy = _DOY_CUM[mo - 1] + dy
    return y + doy / 365.25


def decimal_to_iso(x):
    """小数年 → YYYY-MM-DD：`date_to_decimal` 的**精确逆**（同一 365.25 口径）。

    用途：explorer「VirPhyKit 输入包」的 date 列是小数年，数据接入要先折算成
    ISO 再入数据集（保持数据集单一口径）。往返关系由测试钉住：
    `date_to_decimal(decimal_to_iso(v)) ≈ v`（误差 < 1e-6 年）。

    ⚠️ 与 explorer 侧 `virphykit_export.decimal_year()` 的口径略有差异
    （它用「当年实际天数」，我们用 365.25）→ 换算回来可能与原始采集日期差
    **≤1 天**。这在年/月粒度的采样日期上无影响，而且远小于小数年被
    `to_decimal_year(strict=True)` 错读时的 0.46 年。
    """
    v = float(x)
    y = int(v)
    doy = int(round((v - y) * 365.25))
    doy = max(1, min(365, doy))
    mo = 1
    for k in range(12, 0, -1):
        if _DOY_CUM[k - 1] < doy:
            mo = k
            break
    return '%04d-%02d-%02d' % (y, mo, doy - _DOY_CUM[mo - 1])



# ================================================================
# 1. SeqIDRenamer —— 按 TSV（旧ID\t新ID）批量重命名
# ================================================================

def rename_sequences(seq_path, map_path, out_dir):
    """返回 dict(out_file, n_seq, n_renamed, n_unmapped, n_map_unused, missing示例)。"""
    seqs = read_fasta(seq_path)
    ren = OrderedDict()
    bad_lines = []
    for i, line in enumerate(read_text(map_path).splitlines(), 1):
        line = line.rstrip('\r\n')
        if not line.strip():
            continue
        parts = line.split('\t')
        if len(parts) != 2 or not parts[0].strip():
            bad_lines.append((i, line[:80]))
            continue
        ren[parts[0].strip()] = parts[1].strip()
    if bad_lines:
        raise ValueError('映射表第 %s 行不是「旧ID<TAB>新ID」两列格式（示例: %r）'
                         % (bad_lines[0][0], bad_lines[0][1]))

    out_rows, n_renamed, missing = [], 0, []
    for name in seqs:
        new = ren.get(name)
        if new is None:
            head = name.split()[0] if name else ''
            new = ren.get(head)
            if new is None:
                missing.append(name)
                out_rows.append((name, seqs[name]))
                continue
        n_renamed += 1
        out_rows.append((new, seqs[name]))
    os.makedirs(out_dir, exist_ok=True)
    out_file = os.path.join(out_dir, 'Renamed_' + os.path.basename(seq_path))
    write_fasta(out_file, out_rows)
    used = set()
    for name in seqs:
        head = name.split()[0] if name else name
        for k in (name, head):
            if k in ren:
                used.add(k)
    return {
        'out_file': out_file,
        'n_seq': len(seqs),
        'n_renamed': n_renamed,
        'n_unmapped': len(missing),
        'missing_examples': missing[:10],
        'n_map_unused': len([k for k in ren if k not in used]),
        'n_map_bad': len(bad_lines),
    }


# ================================================================
# 2. SeqGrouper —— GenBank 记录表 + 分组映射
# ================================================================

GB_HEADERS = ['Isolate', 'ID', 'Organism', 'Length', 'Host',
              'Geo Location', 'Collection Date']


def gb_record_rows(records):
    """gb_collection.parse_flatfile 的元数据 dict → SeqGrouper 表行。"""
    rows = []
    for _chunk, m in records:
        geo = m.get('country', '') or ''
        if ':' in geo:
            geo = geo.split(':', 1)[0].strip()
        rows.append({
            'Isolate': m.get('isolate', '') or 'N/A',
            'ID': m.get('acc', ''),
            'Organism': m.get('organism', '') or 'N/A',
            'Length': m.get('length', 0),
            'Host': m.get('host', '') or 'N/A',
            'Geo Location': geo or 'N/A',
            'Collection Date': _norm_gb_date(m.get('date', '')),
        })
    return rows


def _norm_gb_date(s):
    """GenBank collection_date 常见格式 → YYYY-MM-DD / 原样。"""
    s = (s or '').strip()
    for fmt in ('%d-%b-%Y', '%Y-%m-%d', '%d-%b-%Y_%H:%M:%S'):
        try:
            import datetime as _dt
            return _dt.datetime.strptime(s.replace('_', ' '), fmt).strftime('%Y-%m-%d')
        except ValueError:
            continue
    return s or 'N/A'


def load_group_mapping(path):
    """映射表：每行 `组<TAB>值`（值小写匹配；united states 特判为 usa）。"""
    gmap = {}
    for i, line in enumerate(read_text(path).splitlines(), 1):
        line = line.rstrip('\r\n')
        if not line.strip():
            continue
        parts = line.split('\t')
        if len(parts) != 2 or not parts[0].strip() or not parts[1].strip():
            raise ValueError(f'分组映射表第 {i} 行不是「组<TAB>值」两列格式: {line[:80]!r}')
        gmap[parts[1].strip().lower()] = parts[0].strip()
    return gmap


def apply_grouping(rows, gmap=None, col='Geo Location'):
    """按列分组。返回 (带 Group 的行, 组计数, 未匹配值集合)。"""
    counts, unmatched = OrderedDict(), set()
    out = []
    for r in rows:
        v = str(r.get(col, 'N/A') or 'N/A')
        if gmap is None:
            g = v if v != 'N/A' else 'N/A'
        else:
            key = v.lower()
            if key == 'united states':
                key = 'usa'
            g = gmap.get(key)
            if g is None:
                g = 'Unknown'
                if v not in ('', 'N/A'):
                    unmatched.add(v)
        r2 = dict(r)
        r2['Group'] = g
        out.append(r2)
        if g != 'N/A':
            counts[g] = counts.get(g, 0) + 1
    return out, counts, unmatched


def save_group_csv(path, rows):
    headers = GB_HEADERS + ['Group']
    with open(path, 'w', encoding='utf-8-sig', newline='') as f:
        w = csv.writer(f)
        w.writerow(headers)
        for r in rows:
            w.writerow([r.get(h, '') for h in headers])


# ================================================================
# 3. GeoSubsampler —— FASTA 降采样（均等 / 定区剔除 / 随机抽取）
# ================================================================

# VirPhyKit 契约：序列名首段 token 即区域段（如 CCD_AB622861_2008.58 -> CCD）
NAME_REGION_RX = re.compile(r'(?:^|>|_)([A-Za-z0-9]+)(?:_|$)')


def region_of_name(name):
    m = NAME_REGION_RX.search(name)
    return m.group(1) if m else None


# ── 「无区划信息」判定 + 区域坐标表并入（2026-09-18，数据准备合并）──────────
# 为什么单列一个判定：实测 explorer 新版导出的 2,967 条里 `Unknown` 占 2,529（85%）。
# 它**不是一个区划**，只是一句"不知道"；把它当区划参与"各分区等量"会让结果
# 塌成 47 条（等于悄悄丢掉 98% 数据）。所以降采样与统计一律先判定、再排除、并如实报告。
#
# ⚠️ 2026-10-03：词表与 `phylodyn_govern._is_unknown_level` **对齐（唯一权威）** ——
#   治理把 `unknown_region` / `Unknown_AI` / `未` 这类带后缀写法清成空地点，判定
#   若只认精确的 'unknown' 就会与治理口径脱节（同一条数据两种结论）。这里只额外
#   收治理不处理的"符号型空值"（`?` / `-` / `--` / `nan` 等）。
_NO_REGION_SYMBOLS = {'', '?', '-', '--', 'nan'}
_NO_REGION_TOKENS = _NO_REGION_SYMBOLS | {'unknown', 'n/a', 'na', 'none', 'null'}


def is_no_region(v):
    """是不是"没有区划信息"（Unknown / 空 / ? / N/A / `Unknown:Unknown` …）。

    判定 = 符号型空值（本模块口径）**或** `phylodyn_govern._is_unknown_level`
    （治理口径，认 `unknown` 后缀写法）；两处共用同一套"这不是区划"的语义。
    复合写法按 `phylodyn_govern.split_location` 拆段：**每段都是 Unknown**
    （`Unknown:Unknown` / `unknown_region: Unknown`）才算"没有区划"；
    `China: Unknown` 这种还有真实国家的**不**在此列（治理另有 country_only 口径管）。
    """
    t = (v or '').strip()
    if t.lower() in _NO_REGION_SYMBOLS:
        return True
    from Virus_Platform_Core.phylodyn_govern import _is_unknown_level, split_location
    parts = split_location(t)
    return bool(parts) and all(_is_unknown_level(p) for p in parts)


def attach_coords(rows, points, region_key='location'):
    """把区域坐标表并进元数据行：给每行补 `lat`/`lon`（**已有坐标的行不动**）。

    `points` 取 `parse_spatial()` 的返回值（`{'points': [...]}`）或裸列表
    （每项 `{'region','lat','lon'}`）。

    匹配顺序（与平台其余"区域落点"口径一致）：区划名**精确** → **小写** →
    **截到国家**（`Russia: Stavropol` → `Russia`）。
    返回 `(rows2, stats)`；`stats` 含 `n_points` / `n_matched`（按行计）/
    `missing_regions`（坐标表里有、但这批数据没用到的）/ `unplaced_regions`
    （这批数据用到、但坐标表里没有的）。
    """
    if isinstance(points, dict):
        points = points.get('points') or []
    exact, lower, country = {}, {}, {}
    for pt in points:
        rag = (pt.get('region') or '').strip()
        if not rag:
            continue
        try:
            la, lo = float(pt['lat']), float(pt['lon'])
        except (KeyError, TypeError, ValueError):
            continue
        exact.setdefault(rag, (rag, la, lo))
        lower.setdefault(rag.lower(), (rag, la, lo))
        country.setdefault(rag.split(':')[0].strip().lower(), (rag, la, lo))

    out, used, unmatched = [], set(), {}
    n_matched = 0
    for r in rows:
        r2 = dict(r)
        has = str(r2.get('lat') or '').strip() and str(r2.get('lon') or '').strip()
        if not has:
            rag = (r2.get(region_key) or '').strip()
            hit = None
            if rag:
                hit = exact.get(rag) or lower.get(rag.lower())
                if not hit:
                    ck = rag.split(':')[0].strip().lower()
                    hit = country.get(ck)
                if hit:
                    used.add(hit[0])      # 记**坐标表里那一条的**区划名，才能和 points 对账
                else:
                    unmatched[rag] = unmatched.get(rag, 0) + 1
            if hit:
                r2['lat'], r2['lon'] = hit[1], hit[2]
                n_matched += 1
        out.append(r2)
    stats = {'n_points': len(points), 'n_matched': n_matched,
             'missing_regions': sorted(
                 {(p.get('region') or '').strip() for p in points}
                 - used - {''}),
             'unplaced_regions': sorted(unmatched.items(),
                                        key=lambda kv: -kv[1])[:20]}
    return out, stats


# ── 三元组去重（2026-09-18，移植自上游 data_collector.dedup_triplet 的思想）──
# 上游口径：本地 + 线上合并后按 **(date, location, sequence) 三元组**去重，
# 出 `dedup_report.csv`。比"按序列去重"更贴系统发育分析：序列相同但时空不同的
# 记录，对时空分层是有信息量的（不该当重复删掉）；而三者全同的记录，
# 就是同一份观察被导入了两次。
def dedup_rows(rows, seqs=None, *, name_key='name', date_key='date',
               loc_key='location', seq_of=None, require_loc=True):
    """按 (date, location, sequence) 三元组去重。

    返回 `(keep_rows, drop_rows, report)`；**保持原顺序**，保留**首次出现**的那条。
    · `seq_of`（名字→序列）或 `seqs`（路径）都没给时，退化成按 (date, location) 去重，
      并在 report 里写明 `by='date+location'`（**不静默**降级）。
    · ⚠️ 2026-10-03 修：序列索引**同时**以 FASTA 完整头与首词建键 ——
      FASTA 头常带描述（`>NC_001 黄瓜花叶病毒`），而元数据的 name 列可能写完整头、
      也可能只写首词；原实现只按 `nm.split()[0]` 查，头带描述时**恒空串**，
      序列分量等于失效（三元组退化成二元组还不自知）。report 的 `by` 如实标注
      实际生效的键（一条序列都没命中 → 不冒充"按序列去重"）。
    """
    seq_map = seq_of
    if seq_map is None and seqs is not None:
        seq_map = read_fasta(seqs)
    if seq_map is not None:
        # 两键指向同一条序列（完整头 + 首词，先到先得）；不覆盖已有键
        idx = {}
        for h, sq in seq_map.items():
            for k in (h, (h.split()[0] if h else '')):
                if k and k not in idx:
                    idx[k] = sq
        seq_map = idx
    by = 'date+location+sequence' if seq_map else 'date+location'
    seen, keep, drop = {}, [], []
    n_skipped = 0
    n_seq_hit = 0
    for r in rows:
        nm = (r.get(name_key) or '').strip()
        loc = (r.get(loc_key) or '').strip()
        # ⚠️ 地点为空**不参与去重**：治理会把 `Unknown` 清空，不同国家的记录于是
        #    都变成空地点，再按 (date,'',sequence) 去重会把**不同来源**的记录当重复删掉
        #    —— 那是"猜"，违反本模块第一条原则。空地点的行原样保留并单独报数。
        if require_loc and not loc:
            n_skipped += 1
            keep.append(r)
            continue
        sq = ''
        if seq_map is not None:
            sq = seq_map.get(nm) or (seq_map.get(nm.split()[0]) if nm else '')
            if sq:
                n_seq_hit += 1
        key = ((r.get(date_key) or '').strip(), loc, sq)
        if key in seen:
            drop.append({'name': nm, 'dup_of': seen[key], 'date': key[0],
                         'location': key[1]})
        else:
            seen[key] = nm
            keep.append(r)
    if seq_map is not None and not n_seq_hit:
        # 给了序列但从没命中 → 实际生效的是二元组，如实标注（别让报告骗人）
        by = 'date+location（序列索引 0 条命中，等同 date+location）'
    rep = {'by': by, 'n_in': len(rows), 'n_kept': len(keep), 'n_dropped': len(drop),
           'n_seq_matched': n_seq_hit,
           'dropped': drop[:200], 'n_dropped_listed': min(200, len(drop)),
           'groups': len(seen), 'n_skipped_no_loc': n_skipped}
    return keep, drop, rep


def subsample_fasta(fasta_path, out_dir, mode='random', n=50, region=None,
                    region_from='name', meta_path=None, meta_region_col='location',
                    seed=None):
    """三种模式（对齐 GeoSubsampler 口径，seed 可复现）：

    mode='random'      随机抽 n 条 → extract.fas
    mode='region'      从指定 region 剔除 n 条 → extract.fas + remaining.fas
    mode='equal'       各 region 等量（取最小区组内条数）→ extract.fas
    region_from='name' 区域取自序列名首 token（VirPhyKit 契约）
    region_from='meta' 区域取自元数据 CSV 列（平台扩展）
    """
    seqs = read_fasta(fasta_path)
    if not seqs:
        raise ValueError('输入 FASTA 没有序列')
    rng = random.Random(seed)

    regions = None
    if region_from == 'meta':
        if not meta_path:
            raise ValueError('region_from=meta 需要提供元数据 CSV')
        regions = _meta_region_lookup(meta_path, meta_region_col)

    def _region(name):
        if regions is not None:
            key = name.split()[0]
            if key in regions:
                return regions[key]
            head = os.path.splitext(key)[0]
            return regions.get(head, '?')
        return region_of_name(name) or '?'

    if mode == 'equal':
        buckets = OrderedDict()
        for name in seqs:
            buckets.setdefault(_region(name), []).append(name)
        # ⚠️ 2026-09-18：把"没有区划信息"的组（''/`?`/Unknown/N/A）**排除在区集之外**。
        # 不排除的话，实测一批 2,967 条里 Unknown 占 85% → equal 把结果塌成 47 条
        # （等于悄悄丢掉 98% 数据）。排除掉并**如实报数**，让用户自己决定。
        no_info = OrderedDict((r, v) for r, v in buckets.items() if is_no_region(r))
        informative = OrderedDict((r, v) for r, v in buckets.items()
                                  if not is_no_region(r))
        if not informative:
            raise ValueError('各样本都没有可用区划（全是 Unknown/空）—— '
                             '等量降采样无从下手，请先补区划或改用 random 模式')
        sizes = {r: len(v) for r, v in informative.items()}
        k = min(sizes.values())
        picked = []
        for r, names in informative.items():
            picked.extend(rng.sample(names, k))
        out = [(nm, seqs[nm]) for nm in picked]
        os.makedirs(out_dir, exist_ok=True)
        out_file = os.path.join(out_dir, 'extract.fas')
        write_fasta(out_file, out)
        return {'mode': mode, 'out_file': out_file, 'n_picked': len(picked),
                'per_region': sizes, 'per_region_sampled': {r: k for r in sizes},
                'n_excluded_no_region': sum(len(v) for v in no_info.values()),
                'excluded_no_region': {r: len(v) for r, v in no_info.items()},
                'seed': seed}

    if mode == 'region':
        if not region:
            raise ValueError('mode=region 需要指定 region')
        in_r = [nm for nm in seqs if _region(nm) == region]
        rest = [nm for nm in seqs if _region(nm) != region]
        if n > len(in_r):
            raise ValueError(f'要求从 {region} 剔除 {n} 条，但该区只有 {len(in_r)} 条')
        if n == len(in_r):
            raise ValueError('不能剔除该区全部序列，请先备份原文件再操作')
        drop = rng.sample(in_r, n)
        os.makedirs(out_dir, exist_ok=True)
        ex = os.path.join(out_dir, 'extract.fas')
        rm = os.path.join(out_dir, 'remaining.fas')
        write_fasta(ex, [(nm, seqs[nm]) for nm in drop])
        write_fasta(rm, [(nm, seqs[nm]) for nm in rest if nm not in set(drop)])
        return {'mode': mode, 'region': region, 'n_removed': n,
                'extract_file': ex, 'remaining_file': rm,
                'n_remaining': len(seqs) - n, 'seed': seed}

    # mode == 'random'
    if n > len(seqs):
        raise ValueError(f'要求抽取 {n} 条，但总共只有 {len(seqs)} 条')
    picked = rng.sample(list(seqs.keys()), n)
    os.makedirs(out_dir, exist_ok=True)
    out_file = os.path.join(out_dir, 'extract.fas')
    write_fasta(out_file, [(nm, seqs[nm]) for nm in picked])
    return {'mode': mode, 'out_file': out_file, 'n_picked': n,
            'n_total': len(seqs), 'seed': seed}


def _meta_region_lookup(meta_path, col):
    """元数据 CSV/TSV → {name: region}（表头含 name 与指定列，大小写不敏感）。"""
    rows = _read_table(meta_path)
    if 'name' not in rows[0] or col not in rows[0]:
        raise ValueError(f'元数据表需要包含 name 与 {col} 列；实际表头: {list(rows[0])}')
    out = {}
    for r in rows:
        nm = (r.get('name') or '').strip()
        if nm:
            out[nm] = (r.get(col) or '').strip() or '?'
    return out


def _read_table(path):
    """CSV/TSV 自动嗅探 → list[dict]（表头小写化）。统一走 io_utils 容错读取。"""
    from Virus_Platform_Core.io_utils import read_dicts  # 懒加载：本模块保持零包依赖
    if not os.path.isfile(path):
        raise FileNotFoundError(path)          # 与旧 read_text 行为一致
    rows = read_dicts(path, lower_header=True)
    if not rows:
        raise ValueError(f'{path} 没有数据行')
    return rows


# ================================================================
# 4. VirSpaceTime —— Temporal.txt / Spatial.txt 解析
# ================================================================

def parse_spatial(path):
    """Spatial.txt：Region<TAB>Latitude<TAB>Longitude（每行一个样本落点）。"""
    rows = _read_table(path)
    latc = next((c for c in ('latitude', 'lat') if c in rows[0]), None)
    lonc = next((c for c in ('longitude', 'lon', 'long') if c in rows[0]), None)
    regc = 'region' if 'region' in rows[0] else None
    if not (latc and lonc and regc):
        raise ValueError(f'Spatial 表需要 Region/Latitude/Longitude 列；实际表头: {list(rows[0])}')
    pts, bad = [], []
    for i, r in enumerate(rows, 2):
        try:
            la, lo = float(r[latc]), float(r[lonc])
            if not (-90 <= la <= 90 and -180 <= lo <= 180):
                raise ValueError('坐标越界')
        except ValueError:
            bad.append(i)
            continue
        pts.append({'region': r[regc], 'lat': la, 'lon': lo})
    if not pts:
        raise ValueError('Spatial 表没有可用的坐标行')
    return {'points': pts, 'n_rows': len(rows), 'n_bad': len(bad), 'bad_rows': bad[:20],
            'regions': sorted({p['region'] for p in pts})}


def parse_temporal(path):
    """Temporal.txt：Year<TAB>区1<TAB>区2...<TAB>Total（长表每区一列）。

    区域列名保留原始大小写（plotly 图例直接显示；表头匹配仍大小写不敏感）。
    """
    text = read_text(path)
    sample = text[:4096]
    sep = '\t' if sample.count('\t') > sample.count(',') else ','
    rdr = csv.DictReader(io.StringIO(text), delimiter=sep)
    rows = []
    for r in rdr:
        rows.append({(k or '').strip().lower(): (v if v is not None else '').strip()
                     for k, v in r.items() if k is not None})
    if not rows:
        raise ValueError(f'{path} 没有数据行')
    if 'year' not in rows[0]:
        raise ValueError(f'Temporal 表需要 Year 列；实际表头: {list(rows[0])}')
    # 原始表头（保大小写）：小写键 → 原列名
    raw_hdr = {}
    rdr2 = csv.reader(io.StringIO(text), delimiter=sep)
    for raw in rdr2:
        raw_hdr = {(h or '').strip().lower(): (h or '').strip() for h in raw}
        break
    cols = [c for c in rows[0] if c not in ('year', 'total')]
    if not cols:
        raise ValueError('Temporal 表除 Year/Total 外还需要至少一个区域列')
    disp = {c: raw_hdr.get(c, c) for c in cols}
    years, series, totals = [], {c: [] for c in cols}, []
    for r in rows:
        try:
            years.append(int(float(r['year'])))
        except (TypeError, ValueError):
            continue
        for c in cols:
            try:
                series[c].append(int(float(r.get(c) or 0)))
            except (TypeError, ValueError):
                series[c].append(0)
        try:
            totals.append(int(float(r.get('total') or 0)))
        except (TypeError, ValueError):
            totals.append(sum(series[c][-1] for c in cols))
    return {'years': years,
            'regions': [disp[c] for c in cols],
            'series': {disp[c]: series[c] for c in cols},
            'totals': totals}


# ================================================================
# 5. MJRM 生成器 —— BEAST Markov jumps/rewards 矩阵片段
# ================================================================

def mjrm_blocks(traits):
    """离散状态列表 → (矩阵 XML 片段, rewards XML 片段)。

    输出格式与 BEAST1 markovJumpsTreeLikelihood 的登记约定一致：
    每个方向一个 `A-to-B` 指示矩阵参数 + 每个状态一个 reward 参数。
    """
    t = [x.strip() for x in traits if x.strip()]
    if len(t) < 2:
        raise ValueError('至少需要 2 个离散状态（英文逗号分隔）')
    n = len(t)
    mat = []
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            row = ['0'] * n
            row[j] = '1'
            mat.append(f'<parameter id="{t[i]}-to-{t[j]}" value="\n')
            for k in range(n):
                mat.append((' '.join(row) if k == i else ' '.join(['0'] * n)) + '\n')
            mat.append('"/>\n')
    rewards = ['<rewards>\n']
    # rewards 块整体按状态字母序生成（名称与 1.0 位置都是字母序下标，
    # 与 Example PVS_with_matrix.xml 一致；矩阵方向块保持用户输入序。
    # rewards 各行自含，BEAST 语义与次序无关）
    for i, s in enumerate(sorted(t)):
        vals = ['0.0'] * n
        vals[i] = '1.0'
        # 复刻 Example（PVS_with_matrix.xml）的旧版行格式：无缩进 + 双空格
        rewards.append(f'<parameter id="{s}_reward"  value="{" ".join(vals)}" />\n')
    rewards.append('</rewards>\n')
    return ''.join(mat), ''.join(rewards)


_MJRM_MARKERS = ('<!--  END Ancestral state reconstruction',
                 '<!-- END Ancestral state reconstruction',
                 '</markovJumpsTreeLikelihood>')


def mjrm_insert(traits, in_xml_path, out_xml_path):
    """把矩阵与 rewards 片段插到 BEAST XML 的 markovJumpsTreeLikelihood 元素内。

    定位规则：先找第一个 `</markovJumpsTreeLikelihood>` 闭合行，再回溯紧邻其上的
    `<!-- END Ancestral state reconstruction -->` 注释行，插到**注释行之前**
    （片段落在元素内部、END 注释之上——BEAST 才能读到；与 VirPhyKit Example
    的 PVS_with_matrix.xml 逐字节一致）。找不到该元素时退回第一个 END 注释行。
    """
    with open(in_xml_path, 'r', encoding='utf-8') as f:
        lines = f.readlines()
    close_idx = next((i for i, ln in enumerate(lines)
                      if '</markovJumpsTreeLikelihood>' in ln), -1)
    idx = -1
    if close_idx != -1:
        j = close_idx - 1
        while j >= 0 and not lines[j].strip():
            j -= 1
        if j >= 0 and 'END Ancestral state reconstruction' in lines[j]:
            idx = j
        else:
            idx = close_idx
    if idx == -1:
        idx = next((i for i, ln in enumerate(lines)
                    if 'END Ancestral state reconstruction' in ln), -1)
    if idx == -1:
        raise ValueError('XML 中找不到插入标记：</markovJumpsTreeLikelihood> 或 '
                         '<!-- END Ancestral state reconstruction -->（仅支持 '
                         'BEAST1/BEAUti 风格 XML）')
    mat, rewards = mjrm_blocks(traits)
    out = lines[:idx] + [mat, rewards] + lines[idx:]
    with open(out_xml_path, 'w', encoding='utf-8', newline='') as f:
        f.writelines(out)
    return out_xml_path


# ================================================================
# 6. 进化动力学数据接入（三来源 + 时间地点同检 + 汇总）
# ================================================================

# ================================================================
# 元数据列名归一化 + 自动择优（2026-09-18）
#
# explorer 新版「导出 metadata」给的是 25 列（Accession / Collection_Date /
# Release_Date / Geo_Location / Country / Year / Host_Name / …），而本平台的
# 数据集契约只有 name,date,location。与其让用户手改表头，这里做**自动识别**：
# 列名按别名表匹配，同族多列非空时按下面写死的优先级择优，并**逐行记账**用了哪一列
# （口径纪律：系统替用户做的判断必须让人看得见）。
# ================================================================

# ⚠️ 键一律按**小写**匹配（`_sniff_table_text` 读表时已把表头小写化；
#    A0 那边是原样表头，所以这里统一 lower 后再查）。
_EMPTYISH = {'', 'unknown', 'n/a', 'na', 'null', 'none', 'nan', '--', '-'}

# 序列名：越靠前越优先
_META_NAME_KEYS = ('name', 'accession', 'seq_id', 'seqid', 'sequence_id',
                   'isolate', 'strain', 'sample', 'id', '序列名')

# 日期：① 采样类（语义正确且必然更早）② 发布/提交类（兜底）③ 只有年的兜底
_META_DATE_SAMPLING = ('date', 'collection_date', 'sampling_date', 'collect_date',
                       'decimal_date', 'time', '日期', '采集日期', '采样日期')
_META_DATE_RELEASE = ('release_date', 'submission_date', 'publication_date',
                      'release date', 'first_release_date', '发布日期', '提交日期')
_META_DATE_YEARONLY = ('year', 'decimal_year', 'sampling_year', 'collection_year',
                       '年份', '年')

# 地点：① 详细的（国家:省州）② 只到国家的
_META_LOC_DETAILED = ('geo_loc_name', 'geo_location', 'geolocation', 'location',
                      'region', 'state', 'province', 'locality', '地点', '区域')
_META_LOC_COARSE = ('country', '国家')


def _meta_emptyish(v):
    """'Unknown' / 'N/A' / '-' 这类"看着有值其实没信息"的写法。"""
    return (v or '').strip().lower() in _EMPTYISH


def _meta_date_key(v):
    """日期串 → (年, 月, 日) 可比元组（缺月=0、缺日=0）；认不出返回 (9999,0,0)。"""
    s = (v or '').strip()
    m = re.match(r'(\d{4})(?:-(\d{1,2}))?(?:-(\d{1,2}))?', s)
    if m:
        return (int(m.group(1)), int(m.group(2) or 0), int(m.group(3) or 0))
    m = re.match(r'(\d{4})\.\d+', s)
    if m:
        return (int(m.group(1)), 0, 0)
    return (9999, 0, 0)


def _meta_detail_score(v):
    """地点串的"详细度"：带冒号（国家:省州）> 带逗号 > 更长。"""
    s = (v or '').strip()
    n = 0
    if ':' in s or '：' in s:
        n += 3
    if ',' in s or '，' in s:
        n += 2
    if len(s) > 12:
        n += 1
    return n


def meta_view(row):
    """一行（任意大小写键）→ {小写键: 值} 视图。"""
    return {(k or '').strip().lower(): (v.strip() if isinstance(v, str) else v)
            for k, v in (row or {}).items() if isinstance(k, str)}


def meta_pick_date(view):
    """返回 (值, 列名)。采样类 > 发布类 > 只有年；组内多列非空则取**最早**。"""
    for group in (_META_DATE_SAMPLING, _META_DATE_RELEASE, _META_DATE_YEARONLY):
        cands = []
        for k in group:
            v = (view.get(k) or '')
            if isinstance(v, str) and v.strip() and not _meta_emptyish(v):
                cands.append((v.strip(), k))
        if cands:
            cands.sort(key=lambda t: _meta_date_key(t[0]))
            return cands[0]
    return '', ''


def meta_pick_location(view):
    """返回 (值, 列名)。详细类 > 国家类；`Unknown` 让位但不丢行。"""
    best = ('', '', -1)
    for rank, group in ((2, _META_LOC_DETAILED), (1, _META_LOC_COARSE)):
        for k in group:
            v = (view.get(k) or '')
            if not isinstance(v, str) or not v.strip() or _meta_emptyish(v):
                continue
            score = rank * 10 + _meta_detail_score(v)
            if score > best[2]:
                best = (v.strip(), k, score)
    if best[2] >= 0:
        return best[0], best[1]
    # 两列都空/都是 Unknown → 原样保留（历史行为：Unknown 算"有地点"，只是没信息）
    for k in _META_LOC_DETAILED + _META_LOC_COARSE:
        v = (view.get(k) or '')
        if isinstance(v, str) and v.strip():
            return v.strip(), k
    return '', ''


def meta_pick_name(view):
    for k in _META_NAME_KEYS:
        v = (view.get(k) or '')
        if isinstance(v, str) and v.strip():
            return v.strip(), k
    return '', ''


def normalize_meta_rows(rows):
    """给每行**补**上规范列 name / date / location（保留原列），返回 (rows2, notes)。

    - 已有同名规范列的**不覆盖**（显式优先级最高）；
    - 逐行记 `_date_src` / `_loc_src`（用了哪一列），notes 里做汇总；
    - notes['release_date_rows'] 是"用了发布日当采样日"的行数 ——
      这个数必须报给用户，否则时间轴偏晚没人知道。
    """
    from collections import Counter
    out = []
    date_src, loc_src, name_src = Counter(), Counter(), Counter()
    n_release = 0
    for r in rows:
        view = meta_view(r)
        r2 = dict(r)
        if not (view.get('name') or '').strip():
            v, c = meta_pick_name(view)
            if v:
                r2['name'] = v
                name_src[c] += 1
        if not (view.get('date') or '').strip():
            v, c = meta_pick_date(view)
            if v:
                r2['date'] = v
                r2['_date_src'] = c
                date_src[c] += 1
                if c in _META_DATE_RELEASE:
                    n_release += 1
        if not (view.get('location') or '').strip():
            v, c = meta_pick_location(view)
            if v:
                r2['location'] = v
                r2['_loc_src'] = c
                loc_src[c] += 1
        out.append(r2)
    notes = {'n_rows': len(rows),
             'date_src': dict(date_src),
             'location_src': dict(loc_src),
             'name_src': dict(name_src),
             'release_date_rows': n_release}
    if n_release:
        notes.setdefault('warnings', []).append(
            f'{n_release} 行没有采集日期，退用了**发布/提交日期**当采样时间：'
            '发布日通常晚于采样日 → 时间轴会整体偏晚，做定年结论前请留意。')
    return out, notes


_META_REQUIRED = ('name', 'date', 'location')


def _sniff_table_text(text, label, *, govern=True, date_mode='mid'):
    sample = text[:4096]
    sep = '\t' if sample.count('\t') > sample.count(',') else ','
    rdr = csv.DictReader(io.StringIO(text), delimiter=sep)
    rows = []
    for raw in rdr:
        r = {}
        for k, v in raw.items():
            if k is None:
                continue
            r[k.strip().lower()] = (v or '').strip() if isinstance(v, str) else v
        rows.append(r)
    if not rows:
        raise ValueError(f'{label} 没有数据行')
    hdr_in = [h for h in rows[0] if h]
    # ⚠️ 「必需列」要按**原始表头**判（对齐名表），不能用归一化之后的表头：
    #    `normalize_meta_rows` 会给每行**补上**规范列 name/date/location（哪怕是空的），
    #    于是"这张表根本没有日期列"会被放过、拖到 validate_dataset 才报
    #    "没有一条序列通过校验" —— 报错变晚且说不清原因。
    _low = {h.strip().lower() for h in hdr_in}

    def _has_alias(canon):
        groups = {'name': _META_NAME_KEYS, 'date': (_META_DATE_SAMPLING
                                                    + _META_DATE_RELEASE
                                                    + _META_DATE_YEARONLY),
                  'location': _META_LOC_DETAILED + _META_LOC_COARSE}[canon]
        return canon in _low or any(a in _low for a in groups)

    missing = [c for c in _META_REQUIRED if not _has_alias(c)]
    if missing:
        raise ValueError(
            f'{label} 缺少必需列：{",".join(missing)}（表头需为 name,date,location[,lat,lon]）。'
            f'已自动识别的别名：序列名 {list(_META_NAME_KEYS)[:4]}…；'
            f'日期 {list(_META_DATE_SAMPLING)[:3]}…/发布日 {list(_META_DATE_RELEASE)[:2]}…/'
            f'年 {list(_META_DATE_YEARONLY)[:2]}…；'
            f'地点 {list(_META_LOC_DETAILED)[:3]}…/{list(_META_LOC_COARSE)[:2]}…。'
            f'实际表头: {hdr_in}')
    # ---- 列名归一化 + 自动择优（2026-09-18）：explorer 新版 25 列也能直接吃 ----
    rows, _pick_notes = normalize_meta_rows(rows)
    # ---- 入口治理（2026-09-18，移植上游 metadata_governance 口径）----
    # ⚠️ **必须在校验之前**：`DD-Mon-YYYY` / 字面量 `Yyyy-Mm-Dd` / `XX.XX N XXX.XX E`
    #    这类值只有先归一，才能通过严格校验；放在校验之后等于"坏行先被剔掉了"。
    if govern:
        from Virus_Platform_Core import phylodyn_govern as _gv
        rows, _gov_notes = _gv.govern_rows(rows, date_key='date',
                                          loc_key='location', name_key='name',
                                          mode=date_mode)
        _pick_notes['govern'] = _gov_notes
    hdr = [h for h in rows[0] if h and not h.startswith('_')]
    missing = [c for c in _META_REQUIRED if c not in hdr]
    if missing:
        raise ValueError(
            f'{label} 缺少必需列：{",".join(missing)}（表头需为 name,date,location[,lat,lon]）。'
            f'已自动识别的别名：序列名 {list(_META_NAME_KEYS)[:4]}…；'
            f'日期 {list(_META_DATE_SAMPLING)[:3]}…/发布日 {list(_META_DATE_RELEASE)[:2]}…/'
            f'年 {list(_META_DATE_YEARONLY)[:2]}…；'
            f'地点 {list(_META_LOC_DETAILED)[:3]}…/{list(_META_LOC_COARSE)[:2]}…。'
            f'实际表头: {hdr_in}')
    return rows, _pick_notes


def validate_dataset(seqs, rows, strict=False, require_location=True):
    """FASTA + 元数据行 → 校验并配对。

    规则（时间地点必须同时给出，格式必须正确）：
      * name 非空、在 FASTA 中存在；FASTA 每条序列都必须有元数据行
      * date 接受 YYYY / YYYY-MM / YYYY-MM-DD，**以及小数年**（explorer
        「VirPhyKit 输入包」的默认口径，会被折算成 ISO 入库）；日/月/年 仍拒
      * location 非空（除非 `require_location=False`：只做**定年**时可以只要
        日期 —— 对应上游 `drop_no_location=False` 的口径）；可选 lat/lon 为数值且在范围内
      * strict=True：任何一行不合格即抛错；否则剔除不合格行并在报告里逐条说明
    返回 (keep_seqs, keep_rows, issues)
    """
    issues, keep_rows = [], []
    seen = set()
    for i, r in enumerate(rows, 2):
        nm = (r.get('name') or '').strip()
        why = []
        if not nm:
            why.append('name 为空')
        elif nm in seen:
            why.append('name 重复')
        else:
            seen.add(nm)
        raw_date = (r.get('date') or '').strip()
        ok, _v, err = validate_date_str(raw_date, allow_decimal=True)
        if not raw_date:
            why.append('缺少日期（时间地点必须同时输入）')
        elif not ok:
            why.append(f'日期不合规：{err}')
        elif _DEC_RE.fullmatch(raw_date):
            # 小数年（explorer 导出包默认口径）→ 折算成 ISO，保持数据集单一口径。
            # 折算而非透传：见文件头「日期口径」注释。
            r = dict(r, date=decimal_to_iso(raw_date))

        if require_location and not (r.get('location') or '').strip():
            why.append('缺少地点（时间地点必须同时输入）')
        lat, lon = r.get('lat', ''), r.get('lon', '')
        try:
            lav = float(lat) if lat else None
            if lav is not None and not (-90 <= lav <= 90):
                why.append(f'纬度 {lav} 超出 [-90,90]')
        except ValueError:
            why.append(f'纬度不是数字: {lat}')
        try:
            lov = float(lon) if lon else None
            if lov is not None and not (-180 <= lov <= 180):
                why.append(f'经度 {lov} 超出 [-180,180]')
        except ValueError:
            why.append(f'经度不是数字: {lon}')
        if why:
            issues.append({'row': i, 'name': nm or '(空)', 'problems': why})
        else:
            keep_rows.append(r)
    if issues and strict:
        first = issues[0]
        raise ValueError(f'严格模式：第 {first["row"]} 行不合格（{"; ".join(first["problems"])}）；'
                         f'共 {len(issues)} 行不合格，全部见报告')
    # FASTA ↔ 元数据互相配对（已在上面记过问题的行不再重复点名）
    reported = {i['name'] for i in issues if i.get('name')}
    headers = {h.split()[0] if h else h: h for h in seqs}
    meta_names = {r.get('name') for r in keep_rows}
    orphan_meta = [r.get('name') for r in keep_rows if r.get('name') not in headers]
    orphan_fasta = [h for h in seqs if h.split()[0] not in meta_names]
    keep_rows = [r for r in keep_rows if r.get('name') in headers]
    for nm in orphan_meta[:10]:
        if nm not in reported:
            issues.append({'row': '-', 'name': nm,
                           'problems': ['元数据行在 FASTA 中没有对应序列']})
    for h in orphan_fasta[:10]:
        key = h.split()[0] if h else h
        if key not in reported:
            issues.append({'row': '-', 'name': key,
                           'problems': ['FASTA 序列缺少元数据行（时间地点缺失视同缺行）']})
    keep = OrderedDict()
    for h, s in seqs.items():
        key = h.split()[0] if h else h
        if key in {r.get('name') for r in keep_rows}:
            keep[h] = s
    return keep, keep_rows, issues


def _write_dataset(out_dir, seqs, rows):
    os.makedirs(out_dir, exist_ok=True)
    fa = os.path.join(out_dir, 'sequences.fasta')
    write_fasta(fa, list(seqs.items()))
    cs = os.path.join(out_dir, 'metadata.csv')
    # 标准五列**永远**是这五列且这个顺序；只有确实带了 host 才追加一列
    # （没有 host 时 schema 一字不变 → 下游 validate_dataset/govern/A0 不受影响）
    cols = ['name', 'date', 'location', 'lat', 'lon']
    if any((r.get('host') or '').strip() for r in rows):
        cols.append('host')
    with open(cs, 'w', encoding='utf-8-sig', newline='') as f:
        w = csv.writer(f)
        w.writerow(cols)
        for r in rows:
            w.writerow([(r.get(c) or '') for c in cols])
    return fa, cs


def dataset_summary(seqs, rows):
    """导入后的可视化汇总：按年直方 / 按地点计数 / 地图点位。"""
    per_year, per_loc, points = OrderedDict(), OrderedDict(), []
    for r in rows:
        y = (r.get('date') or '')[:4]
        per_year[y] = per_year.get(y, 0) + 1
        loc = r.get('location') or 'Unknown'
        per_loc[loc] = per_loc.get(loc, 0) + 1
        try:
            la, lo = float(r.get('lat') or ''), float(r.get('lon') or '')
            points.append({'name': r.get('name', ''), 'location': loc,
                           'lat': la, 'lon': lo, 'year': y})
        except (TypeError, ValueError):
            pass
    years = sorted(per_year)
    loc_sorted = sorted(per_loc.items(), key=lambda kv: (-kv[1], kv[0]))
    return {
        'n_sequences': len(seqs),
        'n_rows': len(rows),
        'n_locations': len(per_loc),
        'year_span': [years[0], years[-1]] if years else None,
        'years': years,
        'per_year': [per_year[y] for y in years],
        'locations': [k for k, _ in loc_sorted[:50]],
        'location_counts': [v for _, v in loc_sorted[:50]],
        'n_locations_hidden': max(0, len(per_loc) - 50),
        'map_points': points,
        'has_coords': bool(points),
    }


def import_from_files(fasta_path, meta_path, out_dir, strict=False, govern=True, date_mode='mid', require_location=True):
    """来源①/②共用：explorer 导出包内的 fasta+csv 对，或手动整理的文件。"""
    seqs = read_fasta(fasta_path)
    if not seqs:
        raise ValueError(f'{fasta_path} 没有序列')
    rows, _pick_notes = _sniff_table_text(read_text(meta_path), '元数据表', govern=govern, date_mode=date_mode)
    keep, keep_rows, issues = validate_dataset(
        seqs, rows, strict=strict, require_location=require_location)
    if not keep:
        raise ValueError('没有一条序列通过校验（时间/地点必须同时给出且格式正确），详见报告')
    fa, cs = _write_dataset(out_dir, keep, keep_rows)
    summ = dataset_summary(keep, keep_rows)
    summ['meta_pick'] = _pick_notes          # 口径纪律：替用户做的选择要看得见
    summ['govern'] = (_pick_notes or {}).get('govern')
    return {'fasta': fa, 'metadata': cs, 'n_kept': len(keep),
            'n_input': len(seqs), 'n_rows_input': len(rows), 'issues': issues,
            'meta_pick': _pick_notes,
            'govern': (_pick_notes or {}).get('govern'), 'summary': summ}


def import_manual(fasta_text, meta_text, out_dir, strict=False, govern=True, date_mode='mid', require_location=True):
    """来源②：手动输入。FASTA 文本 + 元数据表文本（TSV/CSV，name,date,location[,lat,lon]）。

    ⚠️ 2026-10-03 修：`govern` / `date_mode` / `require_location` 三个参数原来
    在内部调用 `import_from_files` 时被**丢掉**（只传了 strict），于是卡上取消勾选
    「入口治理」、改日期口径（mid/start）、放宽缺地点要求，对手动粘贴这条路**全都不生效**
    —— 现在原样透传（与文件来源同口径）。
    """
    if not fasta_text.strip():
        raise ValueError('请粘贴 FASTA 序列')
    if not meta_text.strip():
        raise ValueError('请粘贴元数据表（name, date, location 三列必需；时间地点必须同时给出）')
    import tempfile
    with tempfile.TemporaryDirectory() as td:
        fa_p = os.path.join(td, 'in.fasta')
        cs_p = os.path.join(td, 'in.tsv')
        with open(fa_p, 'w', encoding='utf-8', newline='\n') as f:
            f.write(fasta_text if fasta_text.endswith('\n') else fasta_text + '\n')
        with open(cs_p, 'w', encoding='utf-8', newline='\n') as f:
            f.write(meta_text if meta_text.endswith('\n') else meta_text + '\n')
        return import_from_files(fa_p, cs_p, out_dir, strict=strict,
                                 govern=govern, date_mode=date_mode,
                                 require_location=require_location)


_GB_DATE_RX = re.compile(r'/collection_date="([^"]+)"')


def gb_chunk_collection_date(gb_chunk, fallback=''):
    """单条 GenBank 文本 → /collection_date 限定符（LOCUS 行日期是提交日期，
    不是采集日期，只作兜底）。"""
    m = _GB_DATE_RX.search(gb_chunk)
    return m.group(1).strip() if m else (fallback or '')


def import_from_genbank(out_dir, collection=None, files=None, accessions=None,
                        term=None, max_records=200, strict=False,
                        species=None, taxid=None, full_length=False,
                        date_from=None, date_to=None,
                        govern=True, date_mode='mid'):
    """来源③：在线下载 / 本地 GenBank → **标准数据集（与我们的元数据 schema 一致）**。

    四种给法（按优先级）：
      · `files` / `collection`  本地 GenBank 文件或已下载集合（离线）
      · `accessions`            登录号列表（在线）
      · `species` / `taxid`     **一键采集**：自动查 taxid + 拼检索式（在线）
      · `term`                  手写 Entrez 检索式（在线，最自由）

    `full_length=True` 只收完整基因组/序列（`complete genome/cds/sequence`），
    过滤片段 —— 与上游 SeqHarvester 同口径。`date_from/date_to` 拼 `[PDAT]`。

    **输出 = 我们的 metadata**：`name,date,location,lat,lon`；GenBank 里若带
    `host` 则**追加一列 host**（没有 host 时 schema 一字不变，向后兼容）——
    带 host 之后可直接把 A0 的性状列（trait）设成 `host` 做宿主分化。
    """
    from Virus_Platform_Core import gb_collection as gbc
    logger = _ListLogger()

    def _collection_text(name):
        cdir = gbc.gb_collection_dir(name)
        paths = [os.path.join(cdir, f) for f in sorted(os.listdir(cdir))
                 if f.lower().endswith(('.gb', '.gbk', '.genbank', '.gbff'))]
        if not paths:
            raise FileNotFoundError(f'集合 [{name}] 内没有 GenBank 文件')
        return '\n'.join(read_text(p) for p in paths)

    if files:
        if isinstance(files, str):
            files = [p.strip() for p in re.split(r'[,\n;]+', files) if p.strip()]
        text = '\n'.join(read_text(p) for p in files)
    elif collection:
        text = _collection_text(collection)
    elif accessions:
        accs = [a.strip() for a in re.split(r'[,;\s]+', accessions) if a.strip()]
        if not accs:
            raise ValueError('accession 列表为空')
        pairs, missing = gbc._fetch_accessions_gb(accs[:max_records], logger.log)
        if missing:
            logger.log(f'{len(missing)} 个 accession 未取到: {missing[:10]}')
        records = pairs
        text = None
    elif term or species or taxid:
        # ---- 一键采集：物种名/taxid → taxid → 检索式（2026-09-18）----
        q = term
        if not q:
            if species and not taxid:
                taxid = gbc.resolve_taxid(species=species, logger=logger.log)
                if taxid:
                    logger.log(f'物种「{species}」→ taxid {taxid}')
                else:
                    logger.log(f'⚠️ 按学名查不到 taxid，退回学名检索「{species}」')
            q = gbc.build_query(species=species, taxid=taxid, term=term,
                                full_length=full_length,
                                date_from=date_from, date_to=date_to)
        logger.log(f'Entrez 检索式：{q}')
        name = 'pdprep_' + _stamp()
        # ⚠️ **多取一些再挑**：RefSeq 参考序列常常**一个 source 限定符都没有**
        #    （实测 NC_002030.1 无 country/collection_date）→ 若只按 max_records
        #    取前 N 条，很可能一条可用的都没有（报"没有记录同时具备 time+地点"）。
        #    所以按 4 倍超取（上限 2000），挑够 max_records 条可用的为止。
        over = min(2000, max(int(max_records or 200) * 4, 40))
        gbc.download_gb_collection(name=name, term=q, max_records=over,
                                   logger=logger)
        logger.collections.append(name)
        text = _collection_text(name)
    else:
        raise ValueError('在线下载需要提供 accession 列表或检索式；也可给本地 .gb 文件/集合名')

    # ⚠️ 2026-09-18 修 P0：原来这里还有一句**无条件**的
    # `records = gbc.parse_flatfile(text)`，它把上面 accession 分支设好的
    # `records = pairs` 覆盖成空（`parse_flatfile(None)` 返回 []）→
    # 「在线下载 · 按 accession」一直静默返回 0 条。
    if text is not None:
        records = gbc.parse_flatfile(text)
    seqs, rows, issues = OrderedDict(), [], []
    for _chunk, m in records:
        acc = m.get('acc', '')
        if not acc or acc in seqs:
            continue
        iso = gbc_records_sequence(_chunk)
        if not iso:
            continue
        # 地点：`geo_loc_name` **优先于** `country`（上游 SeqHarvester 同口径；
        # 现代 GenBank 基本只用 geo_loc_name）。
        # ⚠️ **不再截断到国家**：`Japan:Chiba` 的信息量比 `Japan` 大得多，
        #    归一交给治理层（`phylodyn_govern.parse_location` 会把冒号拆成层级），
        #    截断等于把省级分辨率白扔。
        geo = (m.get('geo_loc_name') or m.get('country') or '').strip()
        # 采集日期优先取 source 的 /collection_date 限定符（LOCUS 行日期是
        # 提交日期，常与采集时间差好几年）
        date_raw = gb_chunk_collection_date(_chunk, m.get('date', '') or '')
        date_ok, date_v, _ = validate_date_str(_gb_date_iso(date_raw))
        keep = date_ok and bool(geo)
        if not keep:
            why = []
            if not date_ok:
                why.append(f'GenBank 记录缺/异 collection_date（原始值: {date_raw or "无"}）')
            if not geo:
                why.append('GenBank 记录缺 country/geo_loc_name')
            issues.append({'row': '-', 'name': acc, 'problems': why})
            continue
        # host：GenBank 有就带上（`_record_summary` 早就解析了，此前组装数据集时丢了）
        host_v = (m.get('host') or '').strip()
        row = {'name': acc, 'date': date_v, 'location': geo, 'lat': '', 'lon': ''}
        if host_v and host_v.upper() not in ('N/A', 'NA', 'UNKNOWN', 'NONE'):
            row['host'] = host_v
        rows.append(row)
        seqs[acc] = iso
    if strict and issues:
        raise ValueError(f'严格模式：{len(issues)} 条记录缺时间或地点，例如 {issues[0]["name"]}: '
                         f'{"; ".join(issues[0]["problems"])}')
    # 在线「一键采集」会超量下载 → 这里截到用户要的条数（只截**可用**的那些，
    # 顺序即 NCBI 命中顺序）；本地文件/显式 accession 列表不截（用户点名要的）。
    if (term or species or taxid) and max_records and len(rows) > int(max_records):
        keep_names = {r['name'] for r in rows[:int(max_records)]}
        logger.log(f'可用 {len(rows)} 条 → 取前 {max_records} 条')
        rows = [r for r in rows if r['name'] in keep_names]
        seqs = OrderedDict((k, v) for k, v in seqs.items() if k in keep_names)
    if not rows:
        raise ValueError('没有记录同时具备 collection_date 与地点，无法构建进化动力学数据集'
                         + (f'（{len(issues)} 条被剔，例如 {issues[0]["name"]}: '
                            f'{"; ".join(issues[0]["problems"])}）' if issues else ''))
    # ⚠️ **也要走治理**：另外两个来源（文件对 / 手动）在 `_sniff_table_text` 里已治理，
    #    GenBank 来源是自己组行的 —— 不治理就会出现"同一平台两种口径"
    #    （实测：`Japan:Chiba` 在这里不归一会原样留下，而文件对来源会变成
    #    `Japan, Chiba`）。
    gov_notes = None
    if govern:
        from Virus_Platform_Core import phylodyn_govern as _gv
        rows, gov_notes = _gv.govern_rows(rows, date_key='date',
                                          loc_key='location', name_key='name',
                                          mode=date_mode)
    fa, cs = _write_dataset(out_dir, seqs, rows)
    return {'fasta': fa, 'metadata': cs, 'n_kept': len(rows),
            'n_input': len(records), 'issues': issues,
            'govern': gov_notes,
            'summary': dataset_summary(seqs, rows)}


def gbc_records_sequence(gb_chunk):
    """单条 GenBank 文本 → 序列（ORIGIN 段）。"""
    m = re.search(r'\nORIGIN[^\n]*\n(.*)', gb_chunk, re.S)
    if not m:
        return ''
    return ''.join(re.findall(r'[acgtnACGTNKYRWSBDMHVrykwsbdmhv]', m.group(1))).upper()


def _gb_date_iso(s):
    """GenBank date（17-FEB-2009 / FEB-2009 不规范 / 2009）→ 尽量 ISO。"""
    s = (s or '').strip()
    months = {m: i + 1 for i, m in enumerate(
        ['JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN',
         'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC'])}
    m = re.fullmatch(r'(\d{1,2})-([A-Za-z]{3})-(\d{4})', s)
    if m and m.group(2).upper() in months:
        return f'{m.group(3)}-{months[m.group(2).upper()]:02d}-{int(m.group(1)):02d}'
    m = re.fullmatch(r'([A-Za-z]{3})-(\d{4})', s)
    if m and m.group(1).upper() in months:
        return f'{m.group(2)}-{months[m.group(1).upper()]:02d}'
    m = re.fullmatch(r'(\d{4})(?:-(\d{2}))?(?:-(\d{2}))?', s)
    if m:
        return '-'.join(x for x in m.groups() if x)
    return s


class _ListLogger:
    """gb_collection 需要的极简 logger（log(msg, level) 签名）/ 下载收集器。"""

    def __init__(self):
        self.lines = []
        self.collections = []

    def log(self, msg, level='INFO'):
        self.lines.append(f'[{level}] {msg}')

    def __call__(self, msg, level='INFO'):
        self.log(msg, level)


def _stamp():
    import time
    return time.strftime('%Y%m%d_%H%M%S')
