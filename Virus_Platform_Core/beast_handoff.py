# -*- coding: utf-8 -*-
"""本地 → 服务器 BEAST 作业包：与 virome_phylo_pipeline 的交接通道。

**两个工作面的分工**（用户的既定计划）
--------------------------------------------------------------
本机（本平台）做得动、且够快的：建树（NJ / FastTree / IQ-TREE / RAxML-NG /
MPBoot / DecentTree）、Fitch 迁移重构、RRT / RSPP、MOTP、地图与 GIF/PDF
—— 实测秒到分钟级（见 run/_bench/bench_summary.json）。
本机**不做**：BEAST 1.10 离散地理 BSSVS 的 MCMC。5 链 × 5000 万 state 是
小时~天级，还要 logcombiner / treeannotator 收尾，属服务器活。
分工的接缝就是本模块：把本机数据打成服务器能直接吃的包（并把两边的命令、
回传清单写死在 README 里），跑完把产物拿回来导入 #t-phylogeo。

**本模块做什么**
  ① 从比对 FASTA + 元数据导出服务器要的两份输入：
       aln.fasta     默认只留有「可解析日期 + 地点」的序列（不静默丢，逐条报）
       metadata.csv  表头固定 name,date,location（其余列原样带上，服务器忽略）
  ② 逐条核对并统计：ID 精确匹配、日期可解析、地点非空、重复 ID、长度是否齐整
  ③ 写出 README.txt（粘贴即用的服务器命令 + 回传步骤）与 job.json（机器可读）
  ④ 预演服务器侧的地点离散化（同样的逗号规则 + ≤12 组），把"服务器会用哪几个
     地点"提前摊开，避免跑完几小时才发现地点名没匹配上

**不做什么**：不跑 MCMC、不联网、不改用户自己的服务器仓库、不下载任何东西。

**服务器侧契约**（逐条读过实现；改服务器脚本时要回来核对这里）
  · utils/beast1_bridge.py
      CSV 表头 name / date / location（location_column 默认 'location'）；
      只有 name+date+location 齐全的行进 meta；FASTA 中 id 不在 meta 的序列被丢；
      ≥4 条有效序列；日期解析 = utils/decimal_year.to_decimal_year(strict=True)
      —— **只认 YYYY[-MM[-DD]]**（也收 '/'）：`2008.58197` 这种小数年会被当成
      "yr=2008.58 + 6 月 15 日"算成 2009.04。所以本模块一律先把小数年折成
      YYYY-MM-DD 再写出（见 _decimal_to_ymd），往返误差 <1 天。
      地点经 _simplify_location 降到省/州级（逗号分隔取倒数第二段、两段取末段），
      超过 12 组时尾部并入 Other。
      替换模型：substitution_model 默认 'auto' → 在 FASTA 同目录 + 上 2 级 +
      phylogeny/tree 子目录 glob *.iqtree / *.log，取第一个含
      "Best-fit model" 行的文件 → **一律走 HKY**（+G 类数取自同一行）；
      一个都没找到 → 保持 auto → 落 **GTR** 分支（其源码注释自述"GTR 易
      特征分解不收敛"，是历史坑）。注意：**本机引擎（RAxML-NG / FastTree）
      的日志服务器不读**——要 HKY 就必须随包给一份 IQ-TREE 格式日志
      （见 export_job 的 iqtree_log，以及 _find_iqtree_log 的自动探测）。
  · run_phylogeo.py
      python run_phylogeo.py --virus X --prior bdsky --chains 5 --threads 8
          --chain-length 50000000 --fasta aln.fasta --metadata metadata.csv
          --work-dir DIR
  · utils/merge_results.py
      python -m utils.merge_results --work-dir DIR [--status|--background]
      → DIR/merged/{merged.log, merged.trees, mcc.tree}（+ 收敛诊断 PDF）
  产物回本机后：本模块不负责导入，走 #t-phylogeo 的「BEAST 日志 / BEAST 树」两栏
  （phylogeo.analyze 的 beast_log / beast_trees 参数）。
"""
import calendar
import csv
import glob
import io
import json
import os
import re
import shutil
import time

__all__ = ['export_job', 'parse_date', 'decimal_to_ymd',
           'simplify_location', 'SERVER_DIR', 'SERVER_HOST', 'WORK_PREFIX']

# 服务器侧路径（读自用户自己的管线：datasets.yaml / phylo_results/_ess_final.py
# 的绝对路径，均为 /home/zhangwenda/MMPV-RNA/virome_phylo_pipeline）
SERVER_HOST = 'zhangwenda@202.119.189.246'
SERVER_DIR = '$HOME/MMPV-RNA/virome_phylo_pipeline'
WORK_PREFIX = 'work_'

_DATE_STRICT = re.compile(r'^(\d{4})(?:[-/](\d{1,2})(?:[-/](\d{1,2}))?)?$')
_DECIMAL_YEAR = re.compile(r'^(\d{4})\.(\d+)$')
_HDR_YEAR = re.compile(r'^\d{4}(?:\.\d+)?$')
# 登录号样子的字符串：SRR12805583 / OR489165.1 —— 头兜底解析的已知误判产物
_ACC_LIKE = re.compile(r'^[A-Za-z]{1,4}\d{5,}(?:\.\d+)?$')


def simplify_location(loc):
    """与服务器 _simplify_location 同规则：逗号分隔 → 省/州级。

    `"China, Ningxia, Yinchuan"` → `Ningxia`；`"China, Ningxia"` → `Ningxia`；
    `"SWM"` → `SWM`（本机区划名原样保留）。
    """
    loc = (loc or '').replace('"', '').strip()
    parts = [p.strip() for p in loc.split(',')]
    if len(parts) >= 3:
        return parts[-2]
    if len(parts) == 2:
        return parts[-1]
    return parts[0]


def _days_in_year(yr):
    return 366 if calendar.isleap(yr) else 365


def _day_to_ymd(yr, doy):
    """年内日序（0-based，可越界进位/退位）→ (年, 月, 日)。"""
    while doy < 0:
        yr -= 1
        doy += _days_in_year(yr)
    while doy >= _days_in_year(yr):
        doy -= _days_in_year(yr)
        yr += 1
    mo = 1
    while doy >= calendar.monthrange(yr, mo)[1]:
        doy -= calendar.monthrange(yr, mo)[1]
        mo += 1
    return yr, mo, doy + 1


def decimal_to_ymd(dec):
    """小数年 → (年, 月, 日)：在服务器可达的"天网格"上取误差最小的那天。

    服务器口径 `yr + (mo-1 + (dy-1)/days_in_month)/12` 的值只能落在**整天**上
    （网格步长 1 天 ≈ 0.0027 年），所以不可能完全精确；这里在相邻三天里挑
    误差最小的那天写出去 —— 年末的小数年（如 1999.99931）会进到次年 1/1 而不是
    被夹到 12/31，避免年末累积出整天以上的偏差。
    """
    yr = int(dec // 1)
    frac = dec - yr
    if frac < 0:
        yr -= 1
        frac += 1.0
    doy = frac * _days_in_year(yr)
    best, best_err = None, None
    for cand in range(int(doy) - 1, int(doy) + 3):
        y2, m2, d2 = _day_to_ymd(yr, cand)
        err = abs(_decimal_of(y2, m2, d2) - dec)
        if best_err is None or err < best_err:
            best, best_err = (y2, m2, d2), err
    return best


def _decimal_of(yr, mo, dy):
    """(年,月,日) → 小数年，口径与服务器 to_decimal_year 完全一致。"""
    dim = calendar.monthrange(yr, mo)[1]
    return yr + (mo - 1 + (min(dy, dim) - 1) / dim) / 12.0


def parse_date(v):
    """日期串 → (小数年, 规范化字符串, 错误原因)。

    认三种：
      · `YYYY-MM-DD` / `YYYY/MM/DD` / `YYYY-MM` / `YYYY` → 原样规范化（服务器直吃）
      · 小数年 `2008.58197`（VirPhyKit FASTA 头那种）→ **折成 YYYY-MM-DD**
        （不能直接给服务器：它的解析会再加约 0.46 年，见模块 docstring）
    其它一律判脏（服务器也判脏并剔除；这里提前拦，避免到服务器才发现）。
    """
    s = str(v or '').strip()
    if not s:
        return None, '', 'empty'
    m = _DATE_STRICT.match(s)
    if m:
        yr = int(m.group(1))
        if not (1900 <= yr <= 2100):
            return None, '', 'year-out-of-range'
        mo = int(m.group(2)) if m.group(2) else None
        dy = int(m.group(3)) if m.group(3) else None
        if mo is not None and not (1 <= mo <= 12):
            return None, '', 'month-out-of-range'
        if dy is not None and not (1 <= dy <= 31):
            return None, '', 'day-out-of-range'
        if mo is None:
            return _decimal_of(yr, 6, 15), '%04d' % yr, ''
        if dy is None:
            return _decimal_of(yr, mo, 15), '%04d-%02d' % (yr, mo), ''
        return _decimal_of(yr, mo, dy), '%04d-%02d-%02d' % (yr, mo, dy), ''
    m = _DECIMAL_YEAR.match(s)
    if m:
        dec = float(s)
        if not (1900 <= dec <= 2100):
            return None, '', 'year-out-of-range'
        yr, mo, dy = decimal_to_ymd(dec)
        # 回写用逆运算后的小数年（保证服务器侧解析出的时间与本机一致）
        return _decimal_of(yr, mo, dy), '%04d-%02d-%02d' % (yr, mo, dy), ''
    return None, '', 'unparsable'


def _read_fasta(path):
    """→ [(id, header, seq)]；id = 头第一个空白符前的串（与服务器 SeqIO 一致）。"""
    out = []
    name = header = None
    buf = []
    with io.open(path, encoding='utf-8-sig', errors='replace') as f:
        for ln in f:
            if ln.startswith('>'):
                if name is not None:
                    out.append((name, header, ''.join(buf)))
                header = ln[1:].strip()
                name = header.split()[0] if header.split() else ''
                buf = []
            elif name is not None:
                buf.append(ln.strip())
    if name is not None:
        out.append((name, header, ''.join(buf)))
    return out


def _read_table(path):
    """CSV/TSV → (表头, [字典行])；编码/分隔符容忍（同 phylogeo.load_metadata）。"""
    with io.open(path, encoding='utf-8-sig', errors='replace', newline='') as f:
        sample = f.read(8192)
        f.seek(0)
        delim = '\t' if sample.count('\t') > sample.count(',') else ','
        rows = list(csv.reader(f, delimiter=delim))
    rows = [r for r in rows if r and any(c.strip() for c in r)]
    if not rows:
        return [], []
    header = [h.strip().lstrip('#') for h in rows[0]]
    out = []
    for r in rows[1:]:
        out.append({header[i]: (r[i] if i < len(r) else '')
                    for i in range(len(header))})
    return header, out


def _pick_key(header):
    for h in header:
        if h.lower() in ('seq_id', 'seqid', 'accession', 'sample', 'name'):
            return h
    return header[0] if header else ''


def _resolve_meta_key(name, meta):
    """FASTA id → 元数据键：先逐字符精确匹配，再退一步试平台头的裸名前缀。

    平台自己的比对文件常把头写成 `名称|区划|年份`（phylogeo 的 `_header_traits`
    就是这么解析的），而服务器的 `r.id` 是**整段**到首个空白 —— 不退这一步，
    用户按平台习惯准备的两个文件在服务器上 0 命中（本机看不出来）。
    返回命中的键，都不中返回 None。
    """
    if not meta:
        return None
    if name in meta:
        return name
    bare = name.split('|', 1)[0].strip()
    if bare and bare != name and bare in meta:
        return bare
    return None


def _header_traits(header):
    """FASTA 头回退：`>acc|区域|年份` 与 VirPhyKit 的 `区域_登录号_小数年`。

    与 phylogeo._header_traits 同规则（② 是兜底，数据集命名不规律时请直接给
    元数据表——本函数认错会把登录号当区划，不静默：命中的会进 job.json 的
    source 字段，README 里也写明每条来自哪里）。
    """
    parts = header.split('|')
    if len(parts) >= 2:
        out = {'region': parts[1].strip()}
        if len(parts) >= 3:
            out['year'] = parts[2].strip()
        return out
    seg = header.split('_')
    if len(seg) == 3 and _HDR_YEAR.match(seg[2].strip()):
        return {'region': seg[0].strip(), 'year': seg[2].strip()}
    return {}


def _find_iqtree_log(aln_path):
    """比对 FASTA 旁边找服务器能识别的 IQ-TREE 格式日志（口径同服务器）。

    beast1_bridge 在 FASTA 同目录 + 上 2 级 + phylogeny/tree 子目录 glob
    *.iqtree / *.log，取第一个含 "Best-fit model" 行的文件。这里照抄同一
    口径 —— 保证「本机自动找到的 = 服务器会用的那一份」。
    """
    cands = []
    d = os.path.dirname(os.path.abspath(aln_path))
    for _ in range(3):
        cands += sorted(glob.glob(os.path.join(d, '*.iqtree')))
        cands += sorted(glob.glob(os.path.join(d, '*.log')))
        for sub in ('phylogeny', 'tree'):
            sd = os.path.join(d, sub)
            if os.path.isdir(sd):
                cands += sorted(glob.glob(os.path.join(sd, '*.iqtree')))
                cands += sorted(glob.glob(os.path.join(sd, '*.log')))
        d = os.path.dirname(d)
    for c in cands:
        try:
            with io.open(c, encoding='utf-8', errors='replace') as f:
                if 'Best-fit model' in f.read():
                    return c
        except OSError:
            continue
    return None


def export_job(aln_path, meta_path=None, out_dir=None, virus_tag=None,
               trait='region', date_trait='year',
               prior='bdsky', chains=5, threads=8, chain_length=50_000_000,
               server_dir=None, server_host=None, subset=True,
               iqtree_log=None, coords_path=None):
    """打一个服务器作业包，返回统计字典（同时落 README.txt / job.json）。

    aln_path    比对 FASTA（必给）
    meta_path   元数据 CSV/TSV（可选；不给她就从 FASTA 头兜底解析）
    out_dir     输出目录（默认 `<fasta 同目录>/beast_job_<病毒名>`）
    virus_tag   作业名（默认取 FASTA 文件名主干）——服务器命令里 --virus 的值
    trait       区划列名（本地约定，默认 region）
    date_trait  日期列名（本地约定，默认 year）
    subset      True（默认）= 只导出有效序列；False = 全导（服务器仍会自己剔）
    iqtree_log  可选的 IQ-TREE 格式日志（.log/.iqtree，含 "Best-fit model" 行）
                —— 服务器**只认这种**日志（读到就走 HKY+G，读不到走 GTR，
                GTR 特征分解易不收敛）。不传时自动在比对 FASTA 旁边找存量
                日志（口径同服务器，见 _find_iqtree_log）；本机引擎
                RAxML-NG/FastTree 的日志服务器不读，不会自动带上。
    coords_path 可选坐标表（区划/纬度/经度）——带上供服务器出图用
    """
    t0 = time.time()
    warnings = []
    if not os.path.isfile(aln_path):
        return {'ok': False, 'error': '比对 FASTA 不存在: %s' % aln_path}

    recs = _read_fasta(aln_path)
    if len(recs) < 4:
        return {'ok': False, 'error': 'FASTA 只有 %d 条序列（服务器要求 ≥4）' % len(recs)}

    # ── 元数据来源：CSV 优先，FASTA 头兜底 ──
    src = 'meta' if meta_path else 'header'
    meta = {}
    extra_cols = []
    key_col = date_col = loc_col = None
    if meta_path:
        if not os.path.isfile(meta_path):
            return {'ok': False, 'error': '元数据不存在: %s' % meta_path}
        header, rows = _read_table(meta_path)
        if not header:
            return {'ok': False, 'error': '元数据表读不出表头: %s' % meta_path}
        key_col = _pick_key(header)
        for cand in (date_trait, 'date', 'collection_date', 'year'):
            if cand in header:
                date_col = cand
                break
        for cand in (trait, 'location', 'region', 'locality'):
            if cand in header:
                loc_col = cand
                break
        if date_col is None or loc_col is None:
            return {'ok': False, 'error': '元数据缺列：日期列（%s/date）或区划列'
                                          '（%s/location）没找到；表头=%s'
                                          % (date_trait, trait, ','.join(header))}
        for h in header:
            if h not in (key_col, date_col, loc_col):
                extra_cols.append(h)
        for r in rows:
            k = (r.get(key_col) or '').strip()
            if k:
                meta[k] = r

    kept, dropped, seen = [], [], set()
    n_dup = 0
    n_id_fix = 0
    id_fix_example = None
    for name, hdr, seq in recs:
        if name in seen:
            n_dup += 1
            dropped.append((name, 'dup-id'))
            continue
        seen.add(name)
        key = _resolve_meta_key(name, meta)
        if meta and key is None:
            dropped.append((name, 'no-meta'))
            continue
        row = meta.get(key) if key else None
        export_id = key or name          # 服务器只认整段 id → 必须落成元数据的键
        if key and key != name:
            n_id_fix += 1
            if id_fix_example is None:
                id_fix_example = (name, key)
        if row is not None:
            date_raw = (row.get(date_col) or '').strip()
            loc_raw = (row.get(loc_col) or '').strip()
            row_src = 'meta'
        else:
            tr = _header_traits(hdr)
            date_raw = (tr.get('year') or '').strip()
            loc_raw = (tr.get(trait) or tr.get('region') or '').strip()
            row_src = 'header' if tr else 'none'
        if not loc_raw:
            dropped.append((name, 'no-location'))
            continue
        dec, canon, err = parse_date(date_raw)
        if dec is None:
            dropped.append((name, 'bad-date' if err != 'empty' else 'no-date'))
            continue
        kept.append({'id': export_id, 'src_id': name, 'seq': seq, 'len': len(seq),
                     'date': canon, 'decimal': round(dec, 6), 'loc': loc_raw,
                     'loc_simple': simplify_location(loc_raw), 'src': row_src})

    # 最常见的踩坑：ID 不匹配（服务器按 id 精确匹配，差一个字符就 0 命中）
    if meta and not kept:
        return {'ok': False,
                'error': '元数据里的序列标识与 FASTA **一条都没匹配上**（服务器按'
                         ' id 精确匹配；已试过裸名前缀，如 `%s` → `%s`）。'
                         '元数据键列=%s，FASTA 头示例=%s'
                         % (str(recs[0][1]).split('|')[0], recs[0][1][:60],
                            key_col, recs[0][1][:60]),
                'dropped': dropped[:10]}
    if len(kept) < 4:
        return {'ok': False,
                'error': '只有 %d 条序列同时有可解析日期 + 地点（服务器要求 ≥4）'
                         % len(kept),
                'dropped': dropped}

    if meta and len(kept) < len(recs) * 0.5:
        warnings.append('元数据只匹配上 %d/%d 条：确认键列（%s）里的标识与 FASTA 头'
                        '逐字符一致（服务器也是精确匹配）' % (len(kept), len(recs), key_col))
    if n_id_fix:
        warnings.append('%d 条 FASTA 头带 `名称|区划|年份` 后缀：导出包里已把 id 截成'
                        '裸名（如 `%s` → `%s`），否则服务器按整段 id 匹配会 0 命中'
                        % (n_id_fix, id_fix_example[0], id_fix_example[1]))
    acc_like = sorted({r['loc_simple'] for r in kept
                       if r['src'] == 'header' and _ACC_LIKE.match(r['loc_simple'])})
    if acc_like:
        warnings.append('区划名疑似登录号（%s）：这是 FASTA 头兜底解析的已知误判'
                        '（`run_登录号_年` 认成区划）——请直接给元数据 CSV'
                        % ', '.join(acc_like[:4]))

    # ── 长度齐整性（服务器按 aln[0] 长度写 XML，参差会被 BEAST 判错） ──
    lens = sorted({r['len'] for r in kept})
    if len(lens) > 1:
        warnings.append('序列长度不齐（%d 种，%d–%d nt）：BEAST 需要等长比对，'
                        '请先用 MAFFT/trimAl 对齐再导出' % (len(lens), lens[0], lens[-1]))
    if n_dup:
        warnings.append('%d 条重复 ID（只保留首条）' % n_dup)

    # ── 服务器侧地点离散化预演（同 simplify + ≤12 组，尾部 Other） ──
    counts = {}
    for r in kept:
        counts[r['loc_simple']] = counts.get(r['loc_simple'], 0) + 1
    ordered = sorted(counts.items(), key=lambda x: (-x[1], x[0]))
    if len(ordered) > 12:
        top = dict(ordered[:11])
        top['Other'] = sum(c for _, c in ordered[11:])
        loc_preview = sorted(top)
        warnings.append('区划 %d 组 > 12：服务器会把尾部 %d 个小组并入 Other'
                        % (len(ordered), len(ordered) - 11))
    else:
        loc_preview = sorted(counts)
    if len(loc_preview) < 2:
        warnings.append('有效区划只有 1 组：服务器会退化成"只定年、不做 BSSVS 地理"')

    # ── 落盘 ──
    tag = virus_tag or os.path.splitext(os.path.basename(aln_path))[0] or 'MYVIRUS'
    tag = re.sub(r'[^0-9A-Za-z_.-]', '_', str(tag))
    if not out_dir:
        out_dir = os.path.join(os.path.dirname(os.path.abspath(aln_path)),
                               'beast_job_' + tag)
    os.makedirs(out_dir, exist_ok=True)
    aln_out = os.path.join(out_dir, 'aln.fasta')
    meta_out = os.path.join(out_dir, 'metadata.csv')
    readme_out = os.path.join(out_dir, 'README.txt')
    job_out = os.path.join(out_dir, 'job.json')

    # subset=False 全导时也只能用**导出的 id**（元数据键）：服务器对没进元数据的
    # 序列会自己剔，但剔之前是按 id 找行 —— 用原始头就一条都找不到
    norm = {r['src_id']: r['id'] for r in kept}
    rows_out = (kept if subset else
                [{'id': norm.get(n, n), 'seq': s} for n, _, s in recs])
    with io.open(aln_out, 'w', encoding='utf-8', newline='\n') as f:
        for r in rows_out:
            f.write('>%s\n' % r['id'])
            sq = r['seq']
            for i in range(0, len(sq), 60):
                f.write(sq[i:i + 60] + '\n')

    with io.open(meta_out, 'w', encoding='utf-8-sig', newline='\n') as f:
        w = csv.writer(f)
        w.writerow(['name', 'date', 'location'] + extra_cols)
        for r in kept:
            r0 = meta.get(r['id']) if meta else None
            w.writerow([r['id'], r['date'], r['loc']]
                       + [(r0.get(c) or '').strip() if r0 else '' for c in extra_cols])

    log_src, log_how = None, None
    if iqtree_log:
        if os.path.isfile(iqtree_log):
            log_src, log_how = iqtree_log, 'explicit'
        else:
            warnings.append('指定的建树日志不存在，未随包带走：%s' % iqtree_log)
    else:
        log_src = _find_iqtree_log(aln_path)
        log_how = 'near_aln' if log_src else None
    if log_src:
        shutil.copyfile(log_src, os.path.join(out_dir, 'iqtree.log'))
    else:
        warnings.append('未附 IQ-TREE 格式日志（服务器探测不到模型日志 → 走 GTR '
                        '分支，特征分解易不收敛）；如需 HKY，请在高级选项里指定'
                        '旧 IQ-TREE 产物（iqtree.iqtree / *.log）')
    if coords_path and os.path.isfile(coords_path):
        shutil.copyfile(coords_path, os.path.join(out_dir, 'coords.tsv'))

    server_dir = server_dir or SERVER_DIR
    host = server_host or SERVER_HOST
    work_dir = os.path.join('%s/%s' % (server_dir, 'jobs'), WORK_PREFIX + tag)
    cmd_run = ('python run_phylogeo.py --virus %s --prior %s --chains %d '
               '--threads %d --chain-length %d \\\n'
               '    --fasta aln.fasta --metadata metadata.csv --work-dir %s'
               % (tag, prior, chains, threads, chain_length, work_dir))
    cmd_status = 'python -m utils.merge_results --work-dir %s --status' % work_dir
    cmd_merge = 'python -m utils.merge_results --work-dir %s' % work_dir

    years = [r['decimal'] for r in kept]
    job = {
        'tool': 'beast_handoff', 'version': 1,
        'created': time.strftime('%Y-%m-%d %H:%M:%S'),
        'virus_tag': tag,
        'inputs': {'aln': os.path.abspath(aln_path),
                   'meta': os.path.abspath(meta_path) if meta_path else None,
                   'meta_source': src},
        'counts': {'in': len(recs), 'exported': len(kept), 'dropped': len(dropped),
                   'dup_id': n_dup, 'id_normalized': n_id_fix},
        'dates': {'min': min(years), 'max': max(years),
                  'span_years': round(max(years) - min(years), 3)},
        'locations': {'n': len(loc_preview), 'list': loc_preview,
                      'counts': dict(ordered[:11]) if len(ordered) > 11 else counts},
        'dropped': [{'id': i, 'reason': r} for i, r in dropped[:200]],
        'server': {'host': host, 'dir': server_dir, 'work_dir': work_dir},
        'commands': {'run': cmd_run, 'status': cmd_status, 'merge': cmd_merge},
        'beast': {'prior': prior, 'chains': chains, 'threads': threads,
                  'chain_length': chain_length,
                  'total_states': chains * chain_length},
        'model_hint': ({'source': os.path.abspath(log_src), 'how': log_how,
                        'server_model': 'HKY+G（服务器读到日志即走 HKY）'}
                       if log_src else
                       {'source': None, 'how': None,
                        'server_model': 'GTR（无日志，服务器 auto 分支）'}),
        'bring_back': ['merged/mcc.tree', 'merged/merged.log', 'merged/merged.trees'],
        'files': {'aln': aln_out, 'metadata': meta_out, 'readme': readme_out},
        'warnings': warnings,
    }
    with io.open(job_out, 'w', encoding='utf-8', newline='\n') as f:
        f.write(json.dumps(job, ensure_ascii=False, indent=1))

    readme = _readme_text(job, kept, extra_cols)
    with io.open(readme_out, 'w', encoding='utf-8-sig', newline='\n') as f:
        f.write(readme)

    out = dict(job)
    out.update({'ok': True, 'out_dir': out_dir, 'seconds': round(time.time() - t0, 2)})
    return out


def _readme_text(job, kept, extra_cols):
    c = job['counts']
    d = job['dates']
    loc = job['locations']
    # ⚠️ 这里**不读** job['beast']（2026-10-03 修复）：旧代码有 `b = job['beast']`
    #    但 b 从未被使用（F841），且 job 缺 'beast' 键时直接 KeyError —— 本函数
    #    的文案只用 counts/dates/locations/server/commands，索性删掉这个赋值。
    s = job['server']
    loc_lines = '\n'.join('    %-24s %d 条' % (k, v)
                          for k, v in sorted(loc['counts'].items(),
                                             key=lambda x: -x[1]))
    drop_lines = '\n'.join('    %s  (%s)' % (x['id'], x['reason'])
                           for x in job['dropped'][:15])
    if len(job['dropped']) > 15:
        drop_lines += '\n    … 余下 %d 条见 job.json 的 dropped' % (
            len(job['dropped']) - 15)
    warn = '\n'.join('  ⚠ ' + w for w in job['warnings']) if job['warnings'] else '  （无）'
    mh = job.get('model_hint') or {}
    if mh.get('source'):
        model_line = ('本包已附 iqtree.log（源：%s）→ 服务器读到 "Best-fit model" '
                      '即走 HKY。' % mh['source'])
    else:
        model_line = ('本包**未附** IQ-TREE 日志 → 服务器探测不到 → 走 GTR 分支'
                      '（其源码自述特征分解易不收敛）。')
    return '''\
=== BEAST 作业包：%(tag)s ===
导出时间 %(created)s   本机只做打包，不跑 MCMC（MCMC 是服务器活）

数据体检
  序列：原 %(n_in)d 条 → 导出 %(n_out)d 条（剔除 %(n_drop)d 条）
  年份：%(y0).3f – %(y1).3f（跨度 %(span).1f 年）
  区划：%(n_loc)d 组（服务器会按逗号取省/州级，>12 组时尾部并入 Other）
%(loc_lines)s
%(drop_lines)s

提醒
%(warn)s

──────── ① 传到服务器 ────────
  scp -r "%(out_dir)s" %(host)s:~/jobs/

──────── ② 服务器上跑（在管线目录里） ────────
  cd %(sdir)s
  %(cmd_run)s

  说明
    · 先小规模试跑（加 --chain-length 2000000）确认 taxa/地点数与日志无异常，
      再上正式链长；--chain-length 是**每条链**的 state 数。
    · --prior 三选一：skyline（文献默认）/ constant / bdsky（Bayesian skyline 之外
      的 birth-death skyline；换 prior 会自动换工作目录名）。
    · 命令跑完只保证"已提交并启动"，不等待；链没跑完别急着合并。
    · 比对模型：服务器选模型只看 IQ-TREE 格式日志里的 "Best-fit model" 行
      （本机引擎 RAxML-NG / FastTree 的日志服务器不读）。
      %(model_line)s

──────── ③ 看进度 / 合并（链全部跑完后） ────────
  %(cmd_status)s
  %(cmd_merge)s
  合并产物：%(work_dir)s/merged/
      mcc.tree / merged.log / merged.trees（+ merged 里还有收敛诊断 PDF）

──────── ④ 拿回本机（三个文件） ────────
  scp %(host)s:%(work_dir)s/merged/mcc.tree    .
  scp %(host)s:%(work_dir)s/merged/merged.log  .
  scp %(host)s:%(work_dir)s/merged/merged.trees .

──────── ⑤ 在本平台导入（系统地理卡片） ────────
  打开「系统地理分析」卡片：
    · BEAST 日志 .log  → merged.log       （逐列后验 + ESS 诊断 + 跳转计数走廊）
    · BEAST 树 .trees  → merged.trees（后验样本，出逐树动画）
                          或 mcc.tree（单棵 MCC，只有保守下界，没有区间/动画）
  两种输入都填好再点运行；弧线线宽会切到 BEAST 口径（样本树占比 /
  P(跳转>0) / MCC 枝端后验下界），三种口径的数字含义**不同**，面板里逐条写明。

  口径提醒：后验概率、BF、ESS 这些数字**全部来自服务器产物**；本机不算
  贝叶斯因子、不冒充后验（本机不跑 MCMC，本地口径从不叫 BEAST/BSSVS）。
''' % {
        'tag': job['virus_tag'], 'created': job['created'],
        'n_in': c['in'], 'n_out': c['exported'], 'n_drop': c['dropped'],
        'y0': d['min'], 'y1': d['max'], 'span': d['span_years'],
        'n_loc': loc['n'], 'loc_lines': loc_lines, 'drop_lines': drop_lines,
        'warn': warn, 'out_dir': os.path.dirname(job['files']['aln']),
        'model_line': model_line,
        'host': s['host'], 'sdir': s['dir'], 'cmd_run': job['commands']['run'],
        'cmd_status': job['commands']['status'], 'cmd_merge': job['commands']['merge'],
        'work_dir': s['work_dir'],
    }
