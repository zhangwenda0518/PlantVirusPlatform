# -*- coding: utf-8 -*-
"""结果摘要层（digest）：把样品各阶段产物压成固定结构的紧凑 JSON。

为什么需要：原始产物（BAM、大 TSV、树文件）远超 LLM 上下文，且多数字段
是路径/指纹等无语义内容。digest 由**确定性代码**读取落盘结果计算关键数字
（检出物种、reads 数、contig 数、一致性矩阵规模…），LLM 只做解读——
科研场景下"数字必须来自文件"是防幻觉的底线。

防御性约定：
  - 任何文件缺失/损坏 → 该节返回 None，绝不让整个 digest 失败；
  - 表格只取表头 + 前 N 行、单元格截断（控制 token 预算）；
  - 不读取任何序列本体（隐私口径：默认出本机的数据只有统计与表行）。
"""
import glob
import json
import os
import re
import time

from Virus_Platform_Core.config import DIRS
from Virus_Platform_Core.utils import safe_open

# 单表最多带给 LLM 的行/列/单元格字符数（token 预算控制）
MAX_ROWS = 12
MAX_COLS = 8
MAX_CELL = 80


def _results_root():
    return DIRS['results']


def _sample_dir(sample):
    """样品名 → results/ 真实目录（不存在返回 None）。"""
    from Virus_Platform_Core.utils import resolve_sample_name
    try:
        real = resolve_sample_name(sample, _results_root())
    except Exception:
        return None
    sd = os.path.join(_results_root(), real)
    return sd if os.path.isdir(sd) else None


def _mtime_iso(path):
    """目录/文件修改时间 → 'YYYY-MM-DD HH:MM'（拿不到返回空串）。"""
    try:
        return time.strftime('%Y-%m-%d %H:%M',
                             time.localtime(os.path.getmtime(path)))
    except OSError:
        return ''


def _dir_brief(root, limit=6):
    """目录概览：名字 + 修改时间 + 一级条目数（**不解析内容**）。

    用于扩展分析（LOGAN / 提交准备）：这些模块没有统一的产物契约，硬解会把
    digest 变成第二个产物解析器，且任何模块改版都会静默污染 AI 的输入。
    只报"有没有、什么时候、大概多大"，足够支撑"我跑过溯源吗"这类问题。
    """
    out = []
    if not os.path.isdir(root):
        return out
    items = []
    for name in os.listdir(root):
        if name.startswith(('_', '.')):
            continue
        p = os.path.join(root, name)
        if not os.path.isdir(p):
            continue
        try:
            items.append((os.path.getmtime(p), name))
        except OSError:
            continue
    items.sort(reverse=True)
    for _mt, name in items[:limit]:
        p = os.path.join(root, name)
        try:
            entries = len(os.listdir(p))
        except OSError:
            entries = 0
        out.append({'name': name[:60], 'mtime': _mtime_iso(p),
                    'entries': entries})
    return out


# 工具箱任务目录名约定：<tool>_<YYYYMMDD>_<HHMMSS>_<hex>
_TOOL_RUN_RE = re.compile(r'^(?P<tool>.+)_\d{8}_\d{6}_[0-9a-fA-F]+$')


def samples_brief(limit=50):
    """全部样品的一行式简况（digest 头部 / 工具 list_samples 用）。

    直接读 results/ 目录 + 各阶段 summary.json 是否存在，不跑
    pipeline_overview（那是 GUI 全量视图，含工具探测，逐样品调用偏重）。
    """
    from Virus_Platform_Core.pipeline import PIPELINE_STAGES
    out = []
    root = _results_root()
    if not os.path.isdir(root):
        return out
    for name in sorted(os.listdir(root))[:limit]:
        if name.startswith('_'):
            continue
        sd = os.path.join(root, name)
        if not os.path.isdir(sd):
            continue
        done, total = 0, 0
        last_done = ''
        for key, _name, dirname, sfile in PIPELINE_STAGES:
            total += 1
            if key == 'report':
                ok = bool(glob.glob(os.path.join(sd, dirname, 'report.html')))
            elif key == 'fastp':
                ok = os.path.isfile(os.path.join(sd, dirname,
                                                 'fastp_R1.fastq.gz'))
            else:
                ok = os.path.isfile(os.path.join(
                    sd, dirname, sfile or 'summary.json'))
            if ok:
                done += 1
                last_done = key
        out.append({'name': name, 'done': done, 'total': total,
                    'last_done_stage': last_done})
    return out


def _tsv_head(path, max_rows=MAX_ROWS, max_cols=MAX_COLS):
    """TSV/CSV 表头 + 前 N 行的通用安全读取（列/单元格截断）。"""
    if not os.path.isfile(path):
        return None
    rows = []
    try:
        with open(path, encoding='utf-8', errors='replace') as f:
            for i, line in enumerate(f):
                if i > max_rows:
                    break
                cols = line.rstrip('\n').split('\t')[:max_cols]
                rows.append([c[:MAX_CELL] for c in cols])
    except OSError:
        return None
    if not rows:
        return None
    return {'file': os.path.basename(path), 'header': rows[0],
            'rows': rows[1:], 'truncated': i >= max_rows}


def _json_at(sd, relpath):
    p = os.path.join(sd, relpath)
    if not os.path.isfile(p):
        return None
    try:
        with open(p, encoding='utf-8') as f:
            return json.load(f)
    except (OSError, ValueError):
        return None


def _pick(d, keys):
    """从 summary dict 里挑标量键（过滤空值与嵌套结构，控制体积）。"""
    out = {}
    for k in keys:
        v = d.get(k)
        if isinstance(v, (str, int, float, bool)) and v != '':
            out[k] = v
    return out


def _round(v, nd=3):
    """数值保留 nd 位小数；非数值/None 原样返回。"""
    try:
        return round(float(v), nd)
    except (TypeError, ValueError):
        return v


def _compact(pairs, nd=3):
    """[(key, value)] → dict：丢掉 None / 空串，浮点统一取 nd 位小数。

    与 `_pick` 的区别：`_pick` 是"按键名从 dict 里挑"，本函数是"按显式键值对
    构造"——扩展分析的 summary 层级深（clock.rtt_fit.r2），且要现算派生值
    （len(warnings)），用 `_pick` 表达不出来。
    """
    out = {}
    for k, v in pairs:
        if v is None or v == '':
            continue
        out[k] = _round(v, nd) if isinstance(v, float) else v
    return out


def _stage_kvsuite(sd):
    s = _json_at(sd, os.path.join('02b_kvsuite', 'summary.json'))
    if not s:
        return None
    out = _pick(s, ['engine', 'kv_lib', 'species_detected', 'viral_reads',
                    'reads_unit', 'viral_pairs', 'passed', 'rejected'])
    t = _tsv_head(os.path.join(sd, '02b_kvsuite', 'identify',
                               'all_viruses.summary.tsv'))
    if t:
        out['identify_table'] = t
    return out


def _stage_assembly(sd):
    s = _json_at(sd, os.path.join('03_assembly', 'summary.json'))
    if not s:
        return None
    out = _pick(s, ['assembly_mode', 'raw_contigs', 'contigs', 'total_bp',
                    'viral_contigs'])
    t = _tsv_head(os.path.join(sd, '03_assembly', 'virus_classification.tsv'),
                  max_rows=15)
    if t:
        out['classification'] = t
    return out


def _stage_hostana(sd):
    s = _json_at(sd, os.path.join('08_host_analysis', 'summary.json'))
    if not s:
        return None
    out = _pick(s, ['n_contigs', 'n_host_warn', 'class_source'])
    cats = s.get('categories')
    if isinstance(cats, dict):
        out['host_categories'] = {str(k)[:40]: v for k, v in
                                  list(cats.items())[:8]}
    prev = s.get('rows_preview')
    if isinstance(prev, list) and prev:
        out['rows_preview'] = prev[:6]
    return out


def _stage_orfa(sd):
    s = _json_at(sd, os.path.join('04b_orf_annot', 'summary.json'))
    if not s:
        return None
    out = _pick(s, ['model', 'engine', 'n_orfs', 'n_annotated', 'n_families',
                    'n_hmm_only', 'n_bfvd_only'])
    fams = s.get('families')
    if isinstance(fams, dict):
        out['top_families'] = dict(list(fams.items())[:8])
    return out


def _stage_phylo(sd):
    s = _json_at(sd, os.path.join('05_phylo', 'summary.json'))
    if not s:
        return None
    out = _pick(s, ['tree_tool', 'skipped', 'reason'])
    groups = s.get('groups')
    if isinstance(groups, list):
        gs = []
        for g in groups[:6]:
            if isinstance(g, dict):
                gs.append(_pick(g, ['name', 'contig', 'tree', 'tree_nwk',
                                    'n_refs', 'sdt', 'identity_csv',
                                    'aln_len', 'n_taxa']))
        out['groups'] = gs
    return out


def _stage_primer(sd):
    s = _json_at(sd, os.path.join('06_primer', 'summary.json'))
    if not s:
        return None
    out = _pick(s, ['n_primers', 'reason', 'mode'])
    t = _tsv_head(os.path.join(sd, '06_primer', 'primers.tsv'))
    if t:
        out['primers_table'] = t
    return out


def _stage_consensus(sd):
    s = _json_at(sd, os.path.join('03c_consensus', 'summary.json'))
    if not s:
        return None
    return _pick(s, ['n_sequences', 'n_variants', 'skipped', 'reason'])


def _stage_verify(sd):
    s = _json_at(sd, os.path.join('03b_verify', 'summary.json'))
    if not s:
        return None
    out = _pick(s, ['n_candidates', 'n_passed', 'n_viroid', 'skipped',
                    'reason'])
    return out


def _stage_subsample(sd):
    """子采样：只报了 N 对 reads，但决定后面所有数字的"分母"，必须让 AI 知道。"""
    s = _json_at(sd, os.path.join('00_prep', 'subsample.json'))
    if not s:
        return None
    return _pick(s, ['n_pairs'])


def _stage_fq2fa(sd):
    """⓪b FASTQ→FASTA 预转换（可选）：skipped=True 表示输入本就是 FASTA。"""
    s = _json_at(sd, os.path.join('00_prep', 'fq2fa.json'))
    if not s:
        return None
    out = _pick(s, ['skipped'])
    files = s.get('files')
    if isinstance(files, list):
        out['n_files'] = len(files)
    return out


def _stage_host(sd):
    """① 宿主去除：**宿主占比是首要质控指标** —— 它比 fastp 的 Q30 更早暴露
    "宿主库选错 / reads 不配对 / 输入为空"。2026-10-08 补：此前 digest 缺这
    一段，AI 被问到"宿主去除正常吗"只能答"数据不支持"。"""
    s = _json_at(sd, os.path.join('01_host_removal', 'stats.json'))
    if not s:
        return None
    out = _pick(s, ['ok', 'host_ratio', 'total_pairs', 'kept_pairs',
                    'dropped_pairs', 'single_end', 'error'])
    hdb = s.get('host_db')
    if isinstance(hdb, dict):
        out['host_db'] = _pick(hdb, ['name', 'species', 'taxid',
                                     'conflicted', 'has_manifest'])
    cp = s.get('classify_params')
    if isinstance(cp, dict):
        out['classify_params'] = _pick(cp, ['confidence', 'min_hit_groups'])
    return out


def _stage_gbdraw(sd):
    """⑨ 基因组图：n_seqs 有值而 n_plots=0 说明"有 contig 却一张图没出"
    （gbdraw 缺装或输入异常），是下游报告缺图的直接原因。"""
    s = _json_at(sd, os.path.join('09_genome_plots', 'summary.json'))
    if not s:
        return None
    out = _pick(s, ['n_seqs', 'standalone'])
    plots = s.get('plots')
    if isinstance(plots, list):
        out['n_plots'] = len(plots)
    return out


def sample_digest(sample):
    """样品级 digest：阶段状态总览 + 各阶段关键数字与表格头。"""
    sd = _sample_dir(sample)
    if not sd:
        return None
    from Virus_Platform_Core.pipeline import PIPELINE_STAGES
    stages = []
    for key, disp, dirname, sfile in PIPELINE_STAGES:
        if key == 'report':
            done = bool(glob.glob(os.path.join(sd, dirname, 'report.html')))
        elif key == 'fastp':
            done = os.path.isfile(os.path.join(sd, dirname,
                                               'fastp_R1.fastq.gz'))
        else:
            done = os.path.isfile(os.path.join(sd, dirname,
                                               sfile or 'summary.json'))
        stages.append({'key': key, 'dir': dirname, 'done': done})
    # fastp 质控数字（reads 前后/Q30）对解读很有价值，单独取
    fp = _json_at(sd, os.path.join('00_prep', 'fastp_report.json')) or {}
    # 顺序与流水线一致（subsample→fastp→fq2fa→host→…），便于 LLM 顺着读。
    # 2026-10-08 补齐主链：此前缺 subsample/fq2fa/host/gbdraw 四段，
    # 「全过程解读」在质控与出图两处断链。
    deep = {'subsample': _stage_subsample(sd),
            'fq2fa': _stage_fq2fa(sd),
            'host': _stage_host(sd),
            'kvsuite': _stage_kvsuite(sd),
            'assembly': _stage_assembly(sd),
            'verify': _stage_verify(sd),
            'consensus': _stage_consensus(sd),
            'hostana': _stage_hostana(sd),
            'orf_annot': _stage_orfa(sd),
            'phylo': _stage_phylo(sd),
            'primer': _stage_primer(sd),
            'gbdraw': _stage_gbdraw(sd)}
    return {'sample': os.path.basename(sd),
            'stages_done': [s['key'] for s in stages if s['done']],
            'fastp': _pick(fp, ['reads_before', 'reads_after', 'q30_rate']),
            'detail': {k: v for k, v in deep.items() if v}}


def task_overview(limit=12):
    """任务中心简况（运行中/最近完成），供"全过程分析"提问的上下文。"""
    try:
        from Virus_Platform_Core.web.tasks import tm
        items = tm.list_all(logs=False) or []
    except Exception:
        return []
    out = []
    for t in items[:limit]:
        out.append({'name': (t.get('name') or '')[:60],
                    'status': t.get('status') or '',
                    'stage': t.get('stage') or '',
                    'pct': t.get('pct')})
    return out


def tool_runs_brief(limit=8):
    """最近的工具箱/扩展分析任务（`run/tool_runs/<tool>_<ts>_<hex>/`）。

    **三条口径**：
      1. **按工具去重**：每个工具只保留最近一次 —— 否则高频工具（实测比对预处理
         `pdprep` 一次批量能连出 14 个目录）会把名额占满，时间线退化成单一工具
         的历史，AI 看不到别的模块；
      2. **契约型工具保证在列**：`_EXT_CONTRACTS` 登记过的工具（进化动力学 /
         系统地理 / RTT / RDP）即便很久没跑也保留其最近一次，这样 AI 能答
         "上次进化动力学结论是什么"——否则同样会被高频工具挤掉；
      3. **有稳定契约才解析数字**：登记过的工具附 `metrics`（关键数字）与
         `input`（跑的是哪份输入）；未登记的工具只报存在性与时间（各工具无统一
         产物契约，硬解会让 digest 变成第二个产物解析器）。
    要看具体表格，走 `read_tool_run_file` 工具（根白名单 = tool_runs/）。
    """
    root = DIRS.get('tool_runs')
    out = []
    if not root or not os.path.isdir(root):
        return out
    items = []
    for name in os.listdir(root):
        if name.startswith(('_', '.')):
            continue
        p = os.path.join(root, name)
        if not os.path.isdir(p):
            continue
        m = _TOOL_RUN_RE.match(name)
        tool = m.group('tool') if m else ''
        try:
            items.append((os.path.getmtime(p), name, tool))
        except OSError:
            continue
    items.sort(key=lambda t: -t[0])
    # 按工具去重（名字不合约定的按目录名当独立键，避免被塌缩成一条）
    seen, dedup = set(), []
    for mt, name, tool in items:
        key = tool or name
        if key in seen:
            continue
        seen.add(key)
        dedup.append((mt, name, tool))
    contract_tools = {pre for pre, _rel, _fn in _EXT_CONTRACTS}
    picked = [t for t in dedup if t[2] in contract_tools]
    rest = [t for t in dedup if t[2] not in contract_tools]
    chosen = (picked + rest)[:limit]
    chosen.sort(key=lambda t: -t[0])
    for _mt, name, tool in chosen:
        p = os.path.join(root, name)
        try:
            entries = len(os.listdir(p))
        except OSError:
            entries = 0
        row = {'run': name[:70], 'tool': tool[:40],
               'mtime': _mtime_iso(p), 'entries': entries}
        row.update(_ext_metrics(tool, p))       # 有稳定契约才加 metrics/input
        out.append(row)
    return out


def _ext_metrics(tool, run_dir):
    """按 `_EXT_CONTRACTS` 读某次扩展分析的关键数字（无契约/读不到 → {}）。"""
    for prefix, rel, fn in _EXT_CONTRACTS:
        if tool != prefix:
            continue
        s = _json_at(run_dir, rel)
        if not s:
            return {}
        m = {'metrics': fn(s)}
        inp = _ext_input(s)
        if inp:
            m['input'] = inp
        return m
    return {}


def _ext_input(s):
    """从 summary 里取"跑的是哪份输入"（文件名即可，路径对 LLM 无信息量）。"""
    for k in ('aln', 'fasta', 'input', 'input_fasta'):
        v = s.get(k)
        if isinstance(v, str) and v:
            return os.path.basename(v.replace('\\', '/'))[:60]
    return ''


def _ext_phylodyn(d):
    """进化动力学：时间树来源 / 定根年 / RTT 拟合 / 迁移事件数。"""
    cl = d.get('clock') or {}
    pg = d.get('phylogeo') or {}
    rtt = cl.get('rtt_fit') or {}
    return _compact([
        ('timetree', cl.get('timetree_source')),
        ('root_year', cl.get('root_year')),
        ('n_methods_ok', cl.get('n_methods_ok')),
        ('rtt_r2', rtt.get('r2')),
        ('n_changes_ml', pg.get('n_changes_ml')),
        ('n_corridors', pg.get('n_corridors')),
        ('n_warn', len(d.get('warnings') or [])),
        ('n_err', len(d.get('errors') or [])),
    ])


def _ext_phylogeo(d):
    """系统地理：序列数 / 区划列 / 迁移事件数 / 未解析数 / 年份跨度。

    另含 **BEAST 后验导入**知情层：`phylogeo/summary.json` 的 `beast` 段非空
    = 这次分析导入了外部 BEAST 产物（`.log` / `.trees`，见 `beast_handoff.py`
    的导出侧）。导入成功与否直接决定弧线线宽用的是**后验**还是 RSPP 复本频率，
    是"结果可信度"的关键分叉，故进 digest。
    """
    m = _compact([
        ('n_seqs', d.get('n_seqs')), ('trait', d.get('trait')),
        ('method', d.get('method')), ('n_transitions', d.get('n_transitions')),
        ('unresolved', d.get('unresolved')), ('n_years', d.get('n_years')),
        ('n_warn', len(d.get('warnings') or [])),
    ])
    b = d.get('beast')
    if isinstance(b, dict):
        m.update(_compact([
            ('beast_ok', not b.get('error')),
            ('beast_err', b.get('error')),
            ('beast_corr', b.get('n_corr')),
            ('beast_units', b.get('n_units')),
            ('beast_log', b.get('log_file')),
            ('beast_trees', b.get('trees_file')),
        ]))
    return m


def _ext_rtt(d):
    """RTT（根到尖端回归）：斜率 / R² / 点数 / 是否强信号。"""
    return _compact([
        ('slope', d.get('slope')), ('r2', d.get('r2')), ('n', d.get('n')),
        ('method', d.get('method')), ('strong_signal', d.get('strong_signal')),
        ('n_warn', len(d.get('warnings') or [])),
    ])


def _ext_rdp(d):
    """重组检测：引擎 / 序列数 / 比对长 / 变异位点 / 事件数 / 质控保留数。"""
    qc = d.get('qc') or {}
    return _compact([
        ('engine', d.get('engine')), ('n_seqs', d.get('n_seqs')),
        ('aln_len', d.get('aln_len')), ('variable', d.get('variable')),
        ('n_events', d.get('n_events')), ('n_keep', qc.get('n_keep')),
    ])


# 扩展分析「稳定契约」登记表：**只有登记在此的工具才解析数字**。
# 每项 = (run 目录名的 tool 前缀, summary 相对路径, 提取函数)。
# 新增工具时**先确认其 summary 是稳定契约**（键名固定、有版本演进意识），
# 再登记；无契约的工具一律走 tool_runs_brief 的存在性 + read_tool_run_file。
_EXT_CONTRACTS = (
    ('phylodyn', 'phylodyn_summary.json', _ext_phylodyn),
    ('phylogeo', os.path.join('phylogeo', 'summary.json'), _ext_phylogeo),
    ('rtt', os.path.join('rtt', 'summary.json'), _ext_rtt),
    ('rdp', os.path.join('rdp', 'events.json'), _ext_rdp),
)


def extensions_brief(limit=6):
    """项目级扩展分析的**存在性**概览：LOGAN 溯源 / NCBI 提交准备。

    与 tool_runs_brief 同一口径（只报有没有、什么时候、多少条目）。这些模块
    多数是"用户手动跑的一次性任务"，不是样品级流水线的一环，所以放在平台级
    上下文里而不是某个样品下。

    公共数据检索（`run/meta_search/`）**不在这里** —— 它另有 `prescreen_brief()`
    给出预筛关键计数（判定通过数等），见该函数。
    """
    out = {}
    for key in ('logan', 'submissions'):
        root = DIRS.get(key)
        if root:
            rows = _dir_brief(root, limit=limit)
            if rows:
                out[key] = rows
    return out


def prescreen_brief(limit=6):
    """病毒样本预筛（/prescreen）与病毒谱的关键计数。

    产物在 `run/meta_search/<项目>/prescreen/`：`stat_run_summary.tsv` 是**权威
    判定文件**（一行一个 run，`pass` 列 = 判定通过）。这里只统计计数，不搬明细：
      - `n_runs`：纳入统计的 run 数（表行数）；
      - `n_pass`：`pass` 列非空的 run 数；
      - `n_ok`：`status == 'OK'`（拿到了分类分析）的 run 数；
      - `virome`：病毒谱缓存 `.cache/virome.json` 是否存在 —— 病毒谱是**纯派生**
        （`public_meta/virome.py` 读预筛产物算），故只报缓存是否就绪、何时生成。

    明细表（哪些 run 通过、检出哪些病毒）走 `read_prescreen_file` 通道
    （2026-10-08 新增，根白名单 = `run/meta_search/<项目>/prescreen/`）。
    这里**只给计数**、不搬明细，是为了不把大表塞进每一轮上下文；模型需要
    明细时按需调那条工具拉前若干行。
    """
    root = DIRS.get('meta_search')
    out = []
    if not root or not os.path.isdir(root):
        return out
    items = []
    for name in os.listdir(root):
        if name.startswith(('_', '.')):
            continue
        pdir = os.path.join(root, name, 'prescreen')
        if not os.path.isdir(pdir):
            continue
        try:
            items.append((os.path.getmtime(pdir), name, pdir))
        except OSError:
            continue
    items.sort(reverse=True)
    for _mt, name, pdir in items[:limit]:
        row = {'project': name[:60], 'mtime': _mtime_iso(pdir)}
        srs = os.path.join(pdir, 'stat_run_summary.tsv')
        if os.path.isfile(srs):
            row.update(_prescreen_counts(srs))
        cache = os.path.join(pdir, '.cache', 'virome.json')
        if os.path.isfile(cache):
            row['virome'] = _mtime_iso(cache)
        out.append(row)
    return out


def _prescreen_counts(path):
    """数 `stat_run_summary.tsv` 的行数与 pass/OK 计数（只读表头 + 统计，不留明细）。"""
    import csv
    try:
        with open(path, encoding='utf-8', errors='replace', newline='') as f:
            rdr = csv.DictReader(f, delimiter='\t')
            fields = rdr.fieldnames or []
            if 'pass' not in fields:
                return {}
            n = n_pass = n_ok = 0
            for r in rdr:
                n += 1
                if (r.get('pass') or '').strip():
                    n_pass += 1
                if (r.get('status') or '').strip() == 'OK':
                    n_ok += 1
        return {'n_runs': n, 'n_pass': n_pass, 'n_ok': n_ok}
    except OSError:
        return {}


def meta_collections_brief(cap=10):
    """公共检索集合清单（名字 + run 数 + 是否含 Unified 元数据）。

    「AI 找数据」场景的关键上下文：没有这份清单，AI 不知道平台上已经
    检索过什么（实测案例：用户问"查宁夏枸杞转录组"，平台里明明躺着
    1171 条 run 的历史集合，AI 却答"当前上下文没有 collection 记录"）。
    """
    import os
    import time as _time
    from Virus_Platform_Core.config import DIRS
    root = DIRS.get('meta_search') or ''
    if not root or not os.path.isdir(root):
        return []
    out = []

    def _rows(path):
        try:
            with open(path, encoding='utf-8-sig', errors='replace') as f:
                return max(0, sum(1 for _ in f) - 1)
        except OSError:
            return None
    try:
        for name in sorted(os.listdir(root), reverse=True):
            d = os.path.join(root, name)
            if not os.path.isdir(d):
                continue
            sf = os.path.join(d, 'search', 'SRA_GSA_Merged_Final.csv')
            if not os.path.isfile(sf):
                continue
            out.append({'name': name, 'n_runs': _rows(sf),
                        'mtime': int(os.path.getmtime(d))})
            if len(out) >= cap:
                break
    except OSError:
        return []
    _ = _time
    return out


def platform_brief():
    """平台级一键简况：数据库就绪状态 + 样品清单 + 任务概览 + 扩展分析存在性
    + 预筛/病毒谱关键计数 + 公共检索集合。"""
    from Virus_Platform_Core.config import get_config
    cfg = get_config()
    try:
        from Virus_Platform_Core.kunpeng import db_ready
        from Virus_Platform_Core.taxonomy import taxonomy_ready
        dbs = {'host': db_ready(cfg.databases['host']),
               'virus': db_ready(cfg.databases['virus']),
               'taxonomy': taxonomy_ready()}
    except Exception:
        dbs = {}
    try:
        tool_runs = tool_runs_brief()
    except Exception:
        tool_runs = []
    try:
        extensions = extensions_brief()
    except Exception:
        extensions = {}
    try:
        prescreen = prescreen_brief()
    except Exception:
        prescreen = []
    try:
        meta_collections = meta_collections_brief()
    except Exception:
        meta_collections = []
    return {'databases': dbs,
            'samples': samples_brief(),
            'tasks': task_overview(),
            'tool_runs': tool_runs,
            'extensions': extensions,
            'prescreen': prescreen,
            'meta_collections': meta_collections}
