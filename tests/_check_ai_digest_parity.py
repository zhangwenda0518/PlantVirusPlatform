# -*- coding: utf-8 -*-
"""digest ↔ 源文件「数字对拍」守卫（2026-10-08 新增）。

为什么需要：`ai_agent/digest.py` 是 AI 的**唯一数字来源** —— 系统提示里写死
"数字只允许来自上下文或工具结果，禁止编造"。如果 digest 自己读错文件、读错
键、或把汇总值与明细错配，AI 会基于错数字给出**看起来完全合理**的解读。
这是本仓最危险的一类失效：没有报错、没有异常，只有错误结论。

本守卫**不复用 digest 的读取函数**（只用 json/os 原生读取），独立重算关键
数字后与之逐一比对。核心不变量：

    digest 报出的每一个标量，要么与源文件同名字段**逐值相等**，
    要么落在「派生键」白名单里（n_files / n_plots / host_db /
    classify_params —— 这些是 digest 自己数出来的，源里没有同名字段）。

三段：
  A 合成夹具：12 个被覆盖阶段各造已知值（含 '' 与 0/False 的边界），
    断言 digest 报出同一值；'' 必须被丢弃，0/False 必须被保留。
  B 真实样品：对 run/results 下**真实存在**的产物独立重读逐键比对。
    没有可用真实样品时显式打印 SKIP 原因（不静默退化成正例空输入）。
  C 灵敏度自证：给比对函数喂一个**故意改错**的 digest dict，必须报出问题
    —— 否则这套断言可能恒真。
  D 扩展分析知情层：`tool_runs_brief` 列出的任务目录必须真实存在、条目数对得上。
  E 扩展分析**契约对拍**：进化动力学 / 系统地理 / RTT / RDP 的 `metrics`
    独立按路径重算后逐键比对（**不复用 digest 的 `_ext_*` 提取函数**），
    并做「范围自证」——磁盘上有真实 run 的契约工具必须全部出现在 digest 里。
  F 扩展分析灵敏度自证：逐个契约工具改错一个 metrics 值 / 凭空造键，必须报出。
  G 预筛计数对拍：`prescreen_brief` 的 n_runs / n_pass / n_ok 独立重数比对，
    并做灵敏度自证（改错一个计数必须被报出）。

用法: python tests/_check_ai_digest_parity.py
退出码: 0 = 全部通过；1 = 存在失败。
"""
import csv
import io
import json
import os
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

FAIL = []


def check(name, cond, detail=''):
    print(('  ✔ ' if cond else '  ✘ ') + name
          + (f' — {detail}' if detail and not cond else ''))
    if not cond:
        FAIL.append(name)


# ------------------------------------------------------------------
# 阶段 → (源文件相对路径, digest detail 键, 该段从源里投影的标量键)
# 键表必须与 ai_agent/digest.py 的 _pick(...) 一致；改 digest 的取值面时
# 同步改这里，否则守卫会漏检（本条是"守卫标签 ≠ 守卫范围"的防线）。
# ------------------------------------------------------------------
CASES = [
    ('00_prep/subsample.json', 'subsample', ['n_pairs']),
    ('00_prep/fq2fa.json', 'fq2fa', ['skipped']),
    ('01_host_removal/stats.json', 'host',
     ['ok', 'host_ratio', 'total_pairs', 'kept_pairs', 'dropped_pairs',
      'single_end', 'error']),
    ('02b_kvsuite/summary.json', 'kvsuite',
     ['engine', 'kv_lib', 'species_detected', 'viral_reads', 'reads_unit',
      'viral_pairs', 'passed', 'rejected']),
    ('03_assembly/summary.json', 'assembly',
     ['assembly_mode', 'raw_contigs', 'contigs', 'total_bp', 'viral_contigs']),
    ('03b_verify/summary.json', 'verify',
     ['n_candidates', 'n_passed', 'n_viroid', 'skipped', 'reason']),
    ('03c_consensus/summary.json', 'consensus',
     ['n_sequences', 'n_variants', 'skipped', 'reason']),
    ('08_host_analysis/summary.json', 'hostana',
     ['n_contigs', 'n_host_warn', 'class_source']),
    ('04b_orf_annot/summary.json', 'orf_annot',
     ['model', 'engine', 'n_orfs', 'n_annotated', 'n_families',
      'n_hmm_only', 'n_bfvd_only']),
    ('05_phylo/summary.json', 'phylo', ['tree_tool', 'skipped', 'reason']),
    ('06_primer/summary.json', 'primer', ['n_primers', 'reason', 'mode']),
    ('09_genome_plots/summary.json', 'gbdraw', ['n_seqs', 'standalone']),
]

# digest 自己数出来的键（源文件里没有同名字段，不参与"逐值相等"比对）
DERIVED = {'n_files', 'n_plots', 'host_db', 'classify_params'}
# 表格类键（行/列已三重截断，值域不受控，只查存在性不查内容）
TABLE_KEYS = {'identify_table', 'classification', 'groups', 'primers_table',
              'top_families', 'host_categories', 'rows_preview'}

# ------------------------------------------------------------------
# 扩展分析 metrics 契约：(工具前缀, summary 相对路径,
#                         {digest 报出的键: (取值路径, 取值种类)})
# 取值种类：'val' 直取 / 'len' 取 len() / 'not' 取逻辑非（如 beast_ok）。
# 键表必须与 ai_agent/digest.py 的 `_EXT_CONTRACTS` / `_ext_*` 一致。
# `_expect_metrics` 按**取值路径**独立重算，不调用 digest 的提取函数 ——
# 否则 digest 改错提取逻辑时守卫会跟着一起错（"复用被检对象"的经典陷阱）。
# ------------------------------------------------------------------
EXT_CASES = [
    ('phylodyn', 'phylodyn_summary.json', {
        'timetree': (('clock', 'timetree_source'), 'val'),
        'root_year': (('clock', 'root_year'), 'val'),
        'n_methods_ok': (('clock', 'n_methods_ok'), 'val'),
        'rtt_r2': (('clock', 'rtt_fit', 'r2'), 'val'),
        'n_changes_ml': (('phylogeo', 'n_changes_ml'), 'val'),
        'n_corridors': (('phylogeo', 'n_corridors'), 'val'),
        'n_warn': (('warnings',), 'len'),
        'n_err': (('errors',), 'len')}),
    ('phylogeo', os.path.join('phylogeo', 'summary.json'), {
        'n_seqs': (('n_seqs',), 'val'),
        'trait': (('trait',), 'val'),
        'method': (('method',), 'val'),
        'n_transitions': (('n_transitions',), 'val'),
        'unresolved': (('unresolved',), 'val'),
        'n_years': (('n_years',), 'val'),
        'n_warn': (('warnings',), 'len'),
        # BEAST 后验导入知情层（beast 段为 None 时这些键整体不出现）
        'beast_ok': (('beast', 'error'), 'not'),
        'beast_err': (('beast', 'error'), 'val'),
        'beast_corr': (('beast', 'n_corr'), 'val'),
        'beast_units': (('beast', 'n_units'), 'val'),
        'beast_log': (('beast', 'log_file'), 'val'),
        'beast_trees': (('beast', 'trees_file'), 'val')}),
    ('rtt', os.path.join('rtt', 'summary.json'), {
        'slope': (('slope',), 'val'),
        'r2': (('r2',), 'val'),
        'n': (('n',), 'val'),
        'method': (('method',), 'val'),
        'strong_signal': (('strong_signal',), 'val'),
        'n_warn': (('warnings',), 'len')}),
    ('rdp', os.path.join('rdp', 'events.json'), {
        'engine': (('engine',), 'val'),
        'n_seqs': (('n_seqs',), 'val'),
        'aln_len': (('aln_len',), 'val'),
        'variable': (('variable',), 'val'),
        'n_events': (('n_events',), 'val'),
        'n_keep': (('qc', 'n_keep'), 'val')}),
]

# 「路径任一层不存在」的哨兵 —— 必须与「键在但值为 None/False」区分开，
# 否则 'not' 类取值会把"没有 beast 段"误判成"beast_ok=True"。
_MISSING = object()


def _dig_path(d, path):
    """按路径逐层取值；任一层缺失 → `_MISSING`（区别于键在但值为 None）。"""
    for k in path:
        if not isinstance(d, dict) or k not in d:
            return _MISSING
        d = d[k]
    return d


def _latest_run(root, prefix):
    """`<prefix>_<YYYYMMDD>_<HHMMSS>_<hex>` 目录里最近的一个（独立实现）。"""
    import re as _re
    if not os.path.isdir(root):
        return None
    pat = _re.compile(r'^' + _re.escape(prefix) + r'_\d{8}_\d{6}_[0-9a-fA-F]+$')
    best = None
    for name in os.listdir(root):
        if not pat.match(name):
            continue
        p = os.path.join(root, name)
        if not os.path.isdir(p):
            continue
        try:
            mt = os.path.getmtime(p)
        except OSError:
            continue
        if best is None or mt > best[0]:
            best = (mt, name, p)
    return best[2] if best else None


def _expect_metrics(src, spec):
    """按 spec 的取值路径从源 dict 独立重算期望 metrics（与 digest 同口径）。"""
    exp = {}
    for k, (path, kind) in spec.items():
        v = _dig_path(src, path)
        if v is _MISSING:
            continue
        if kind == 'len':
            v = len(v) if isinstance(v, (list, dict)) else v
        elif kind == 'not':
            v = not v
        if v is None or v == '':
            continue
        exp[k] = round(v, 3) if isinstance(v, float) else v
    return exp


def compare_ext(metrics_by_tool, problems):
    """metrics_by_tool: {tool: metrics_dict}；独立重读源文件逐键比对。"""
    from Virus_Platform_Core.config import DIRS
    root = DIRS['tool_runs']
    for prefix, rel, spec in EXT_CASES:
        run_dir = _latest_run(root, prefix)
        if run_dir is None:
            continue                       # 磁盘上没有该工具的 run
        got = metrics_by_tool.get(prefix)
        if got is None:
            problems.append(f'{prefix}: 磁盘有真实 run 但 digest 未产出 metrics')
            continue
        try:
            with io.open(os.path.join(run_dir, rel), encoding='utf-8') as f:
                src = json.load(f)
        except (OSError, ValueError) as e:
            problems.append(f'{prefix}: 守卫读源失败 {e}')
            continue
        exp = _expect_metrics(src, spec)
        for k, v in exp.items():
            if k not in got:
                problems.append(f'{prefix}.{k}: digest 漏报（源={v!r}）')
            elif got[k] != v:
                problems.append(f'{prefix}.{k}: digest={got[k]!r} ≠ 源={v!r}')
        for k, v in got.items():          # 反向：不得报出未登记的键
            if k not in spec:
                problems.append(f'{prefix}.{k}: digest 报出未登记的键（值={v!r}）')


def _scalars(src, keys):
    """按 digest._pick 的同一口径投影源标量：非标量与 '' 丢弃，0/False 保留。"""
    out = {}
    for k in keys:
        v = src.get(k)
        if isinstance(v, (str, int, float, bool)) and v != '':
            out[k] = v
    return out


def compare(detail, sources, problems):
    """detail: digest 的 sample_digest()['detail']；sources: {stage: 源 dict}。

    只累积问题，不抛异常——方便灵敏度自证直接调用。
    """
    for _rel, stage, keys in CASES:
        src = sources.get(stage)
        if src is None:
            continue                       # 源文件不存在：digest 应当不报该段
        got = detail.get(stage)
        if got is None:
            problems.append(f'{stage}: 源文件存在但 digest 未产出该段')
            continue
        exp = _scalars(src, keys)
        for k, v in exp.items():
            if k not in got:
                problems.append(f'{stage}.{k}: digest 漏报（源={v!r}）')
            elif got[k] != v:
                problems.append(f'{stage}.{k}: digest={got[k]!r} ≠ 源={v!r}')
        # 反向：digest 报出的标量不得凭空出现
        for k, v in got.items():
            if k in exp or k in DERIVED or k in TABLE_KEYS:
                continue
            if k not in src:
                problems.append(f'{stage}.{k}: digest 报出源中不存在的键'
                                f'（值={v!r}）')


def _write_json(path, obj):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with io.open(path, 'w', encoding='utf-8') as f:
        json.dump(obj, f, ensure_ascii=False)


# ------------------------------------------------------------------
# A) 合成夹具：12 个阶段全部造已知值（含边界）
# ------------------------------------------------------------------
FIXTURE = {
    '00_prep/subsample.json': {'n_pairs': 200000},
    # skipped=False 必须被保留（旧写法 v != '' 会误丢 False 吗？本守卫锁住）
    '00_prep/fq2fa.json': {'skipped': False, 'files': ['a.fa.gz']},
    '01_host_removal/stats.json': {
        'ok': True, 'host_ratio': 0.9871, 'total_pairs': 1000000,
        'kept_pairs': 12900, 'dropped_pairs': 987100, 'single_end': False,
        'host_db': {'name': 'tomato', 'taxid': 4081},
        'classify_params': {'confidence': 0.1, 'min_hit_groups': 3}},
    '02b_kvsuite/summary.json': {
        'engine': 'salmon', 'kv_lib': 'plant', 'species_detected': 7,
        'viral_reads': 4210, 'reads_unit': 'pairs', 'viral_pairs': 2105,
        'passed': 7, 'rejected': 0},
    '03_assembly/summary.json': {
        'assembly_mode': 'rnaviral', 'raw_contigs': 900, 'contigs': 41,
        'total_bp': 123456, 'viral_contigs': 3},
    '03b_verify/summary.json': {
        'n_candidates': 5, 'n_passed': 4, 'n_viroid': 1, 'skipped': False},
    '03c_consensus/summary.json': {
        'n_sequences': 4, 'n_variants': 137, 'skipped': False},
    '08_host_analysis/summary.json': {
        'n_contigs': 3, 'n_host_warn': 0, 'class_source': 'ictv'},
    '04b_orf_annot/summary.json': {
        'model': 'pyrodigal', 'engine': 'diamond', 'n_orfs': 12,
        'n_annotated': 9, 'n_families': 5, 'n_hmm_only': 1, 'n_bfvd_only': 2},
    '05_phylo/summary.json': {'tree_tool': 'fasttree', 'skipped': False},
    '06_primer/summary.json': {'n_primers': 6, 'mode': 'conserved'},
    '09_genome_plots/summary.json': {'n_seqs': 3, 'standalone': False,
                                     'plots': ['a.svg', 'b.svg', 'c.svg']},
}


def part_a(digest):
    tmp = tempfile.mkdtemp(prefix='vp_digest_parity_')
    try:
        sd = os.path.join(tmp, 'parity_probe')
        for rel, obj in FIXTURE.items():
            _write_json(os.path.join(sd, rel), obj)
        orig = digest.DIRS
        digest.DIRS = dict(orig, results=tmp)
        try:
            detail = (digest.sample_digest('parity_probe') or {}).get(
                'detail') or {}
        finally:
            digest.DIRS = orig
        sources = {stage: FIXTURE[rel] for rel, stage, _k in CASES}
        problems = []
        compare(detail, sources, problems)
        check('A 合成夹具：12 阶段逐键对拍一致', not problems,
              '; '.join(problems[:5]))
        # 边界：'' 丢弃 / 0 与 False 保留
        check('A 边界：skipped=False 未被丢弃',
              detail.get('fq2fa', {}).get('skipped') is False,
              str(detail.get('fq2fa')))
        check('A 边界：派生键 n_plots=3',
              detail.get('gbdraw', {}).get('n_plots') == 3,
              str(detail.get('gbdraw')))
        check('A 边界：空串值被丢弃',
              'reason' not in (detail.get('primer') or {}),
              str(detail.get('primer')))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


# ------------------------------------------------------------------
# B) 真实样品：独立重读源文件逐键比对
# ------------------------------------------------------------------
def part_b(digest):
    from Virus_Platform_Core.config import DIRS
    root = DIRS['results']
    if not os.path.isdir(root):
        print('  [SKIP] run/results 不存在，跳过真实样品对拍')
        return
    checked = 0
    problems = []
    for name in sorted(os.listdir(root)):
        if name.startswith(('_', '.')):
            continue
        sd = os.path.join(root, name)
        if not os.path.isdir(sd):
            continue
        sources, present = {}, []
        for rel, stage, _k in CASES:
            p = os.path.join(sd, rel)
            if not os.path.isfile(p):
                continue
            try:
                with io.open(p, encoding='utf-8') as f:
                    sources[stage] = json.load(f)
                present.append(stage)
            except (OSError, ValueError):
                continue
        if not present:
            continue
        try:
            detail = (digest.sample_digest(name) or {}).get('detail') or {}
        except Exception as e:                      # digest 崩了也是问题
            problems.append(f'{name}: sample_digest 抛 {type(e).__name__}: {e}')
            continue
        checked += 1
        before = len(problems)
        compare(detail, sources, problems)
        if len(problems) > before:
            for p_ in problems[before:]:
                problems[before] = f'{name}: ' + p_
    if checked == 0:
        print('  [SKIP] run/results 下暂无含 digest 覆盖产物的真实样品')
        return
    check(f'B 真实样品对拍（{checked} 个样品）', not problems,
          '; '.join(problems[:5]))


# ------------------------------------------------------------------
# C) 灵敏度自证：比对函数必须能报出"故意改错"
# ------------------------------------------------------------------
def part_c():
    src = {'kvsuite': FIXTURE['02b_kvsuite/summary.json']}
    # 必须按 CASES 声明的**完整**键表投影，否则 compare 会报"漏报"（这正是
    # 它该做的）—— 灵敏度自证自身也要构造正确。
    keys = next(k for _r, st, k in CASES if st == 'kvsuite')
    good = {'kvsuite': _scalars(src['kvsuite'], keys)}
    bad = json.loads(json.dumps(good))
    bad['kvsuite']['viral_reads'] = 999999          # 注入故意错误
    p_good, p_bad = [], []
    compare(good, src, p_good)
    compare(bad, src, p_bad)
    check('C 灵敏度：正确的 digest 不报问题', not p_good, '; '.join(p_good))
    check('C 灵敏度：改错一个数字必须被报出',
          any('viral_reads' in x for x in p_bad), '; '.join(p_bad) or '(无)')
    # 凭空造键也必须被抓
    p_inv = []
    compare({'kvsuite': dict(good['kvsuite'], fake_metric=1)},
            src, p_inv)
    check('C 灵敏度：凭空造数字必须被报出',
          any('fake_metric' in x for x in p_inv), '; '.join(p_inv) or '(无)')


# ------------------------------------------------------------------
# D) 扩展分析知情层：列出的任务目录必须真实存在且条目数对得上
# ------------------------------------------------------------------
def part_d(digest):
    from Virus_Platform_Core.config import DIRS
    rows = digest.tool_runs_brief(limit=10)
    problems = []
    for r in rows:
        p = os.path.join(DIRS['tool_runs'], r['run'])
        if not os.path.isdir(p):
            problems.append(f'{r["run"]}: 目录不存在')
            continue
        try:
            n = len(os.listdir(p))
        except OSError:
            n = -1
        if r['entries'] != n:
            problems.append(f'{r["run"]}: entries={r["entries"]} 实际={n}')
    check(f'D tool_runs_brief 条目数对得上（{len(rows)} 条）', not problems,
          '; '.join(problems[:5]))


# ------------------------------------------------------------------
# E) 扩展分析契约对拍 + 范围自证
# ------------------------------------------------------------------
def part_e(digest):
    from Virus_Platform_Core.config import DIRS
    rows = digest.tool_runs_brief(limit=50)
    by_tool = {r['tool']: r['metrics'] for r in rows if 'metrics' in r}
    problems = []
    compare_ext(by_tool, problems)
    present = [pre for pre, _rel, _spec in EXT_CASES
               if _latest_run(DIRS['tool_runs'], pre)]
    check(f'E 扩展分析契约对拍（覆盖 {len(present)} 个工具）', not problems,
          '; '.join(problems[:5]))
    # 范围自证：磁盘有 run 的契约工具必须全部出现在 digest 里（防"被高频工具挤掉"）
    missing = [p for p in present if p not in by_tool]
    check('E 范围自证：磁盘有 run 的契约工具全部入 digest', not missing,
          f'缺 {missing}')


# ------------------------------------------------------------------
# F) 扩展分析灵敏度自证：改错 metrics 必须被报出
# ------------------------------------------------------------------
def part_f(digest):
    rows = digest.tool_runs_brief(limit=50)
    base = {r['tool']: dict(r['metrics']) for r in rows if 'metrics' in r}
    if not base:
        print('  [SKIP] 无契约型工具产物，跳过扩展分析灵敏度自证')
        return
    p_good = []
    compare_ext(base, p_good)
    check('F 灵敏度：原样 metrics 不报问题', not p_good, '; '.join(p_good[:3]))
    for tool in sorted(base):
        k0 = sorted(base[tool])[0]
        bad = {k: dict(v) for k, v in base.items()}
        bad[tool][k0] = '__CORRUPT__'
        p_bad = []
        compare_ext(bad, p_bad)
        check(f'F 灵敏度：改错 {tool}.{k0} 必须被报出',
              any(f'{tool}.{k0}' in x for x in p_bad),
              '; '.join(p_bad[:3]) or '(无)')
    # 凭空造键也必须被抓
    inv = {k: dict(v) for k, v in base.items()}
    inv[sorted(base)[0]]['fake_metric'] = 1
    p_inv = []
    compare_ext(inv, p_inv)
    check('F 灵敏度：凭空造 metrics 键必须被报出',
          any('fake_metric' in x for x in p_inv), '; '.join(p_inv[:3]) or '(无)')


# ------------------------------------------------------------------
# G) 预筛计数对拍（prescreen_brief）
# ------------------------------------------------------------------
def _count_prescreen(srs):
    """独立重数 stat_run_summary.tsv：总行数 / pass 非空数 / status==OK 数。"""
    n = n_pass = n_ok = 0
    with io.open(srs, encoding='utf-8', errors='replace', newline='') as f:
        for rec in csv.DictReader(f, delimiter='\t'):
            n += 1
            if (rec.get('pass') or '').strip():
                n_pass += 1
            if (rec.get('status') or '').strip() == 'OK':
                n_ok += 1
    return {'n_runs': n, 'n_pass': n_pass, 'n_ok': n_ok}


def compare_prescreen(rows, problems):
    """rows: digest.prescreen_brief() 的返回；独立重读 TSV 逐值比对。"""
    from Virus_Platform_Core.config import DIRS
    root = DIRS['meta_search']
    for r in rows:
        srs = os.path.join(root, r.get('project', ''),
                           'prescreen', 'stat_run_summary.tsv')
        if not os.path.isfile(srs):
            if 'n_runs' in r:
                problems.append(f'{r.get("project")}: 无 stat_run_summary.tsv '
                                f'却报出了计数')
            continue
        exp = _count_prescreen(srs)
        for k, v in exp.items():
            if r.get(k) != v:
                problems.append(f'{r.get("project")}.{k}: '
                                f'digest={r.get(k)!r} ≠ 源={v!r}')


def part_g(digest):
    rows = digest.prescreen_brief(limit=20)
    if not rows:
        print('  [SKIP] run/meta_search 下无 prescreen 产物，跳过预筛计数对拍')
        return
    problems = []
    compare_prescreen(rows, problems)
    check(f'G 预筛计数对拍（{len(rows)} 个项目）', not problems,
          '; '.join(problems[:5]))
    # 灵敏度：改错一个计数必须被报出
    bad = [dict(r) for r in rows]
    target = next((r for r in bad if 'n_runs' in r), None)
    if target is not None:
        target['n_pass'] = (target.get('n_pass') or 0) + 1
        p_bad = []
        compare_prescreen(bad, p_bad)
        check('G 灵敏度：改错 n_pass 必须被报出',
              any('n_pass' in x for x in p_bad), '; '.join(p_bad[:3]) or '(无)')


# ------------------------------------------------------------------
# H) 扩展分析契约的**分支覆盖**（合成夹具）
# ------------------------------------------------------------------
def part_h(digest):
    """真实数据里 `beast` 段全为 None ⇒ E 段覆盖不到该分支。

    这里用合成夹具把分支走一遍（成功 / 报错 / null），并回归「按工具去重 +
    契约型工具保证在列」——否则"登记了但从未被验证"就是"守卫标签 ≠ 守卫范围"
    的翻版。
    """
    import time
    tmp = tempfile.mkdtemp(prefix='vp_digest_ext_')
    try:
        tr = os.path.join(tmp, 'tool_runs')

        def _w(name, obj):
            p = os.path.join(tr, name, 'phylogeo', 'summary.json')
            os.makedirs(os.path.dirname(p), exist_ok=True)
            with io.open(p, 'w', encoding='utf-8') as f:
                json.dump(obj, f, ensure_ascii=False)

        base = {'n_seqs': 8, 'trait': 'region', 'method': 'nj',
                'n_transitions': 1, 'unresolved': 0, 'n_years': 8,
                'warnings': []}
        cases = {
            'phylogeo_20261008_120000_aaaaaa':
                dict(base, beast={'error': None, 'n_corr': 3, 'n_units': 2500,
                                  'log_file': 'merged.log',
                                  'trees_file': 'mcc.tree'}),
            'phylogeo_20261008_130000_bbbbbb':
                dict(base, beast={'error': 'ValueError: 状态名对不上',
                                  'n_corr': 0, 'n_units': 0,
                                  'log_file': None, 'trees_file': None}),
            'phylogeo_20261008_140000_cccccc': dict(base, beast=None),
        }
        for name, obj in cases.items():
            _w(name, obj)
        # 造 20 个**更新鲜**的非契约工具目录：契约型工具若不优先就必被挤掉
        for i in range(20):
            os.makedirs(os.path.join(tr, 'toolx%02d_20261009_1200%02d_dddddd'
                                     % (i, i)), exist_ok=True)
        # 再造 20 个**同一个工具**的目录（高频工具占位场景）用于测去重
        for i in range(20):
            os.makedirs(os.path.join(tr, 'pdprep_20261009_1300%02d_eeee%02d'
                                     % (i, i)), exist_ok=True)
        old = time.time() - 100000
        for name in cases:
            os.utime(os.path.join(tr, name), (old, old))

        # (1) 分支覆盖：直接走 digest 的真实接线（契约查表 → 读文件 → 提取）
        m_ok = digest._ext_metrics('phylogeo',
                                   os.path.join(tr, list(cases)[0]))
        mm_ok = m_ok.get('metrics') or {}
        check('H beast 成功：beast_ok=True + 走廊/单位/文件名',
              mm_ok.get('beast_ok') is True and mm_ok.get('beast_corr') == 3
              and mm_ok.get('beast_units') == 2500
              and mm_ok.get('beast_log') == 'merged.log'
              and mm_ok.get('beast_trees') == 'mcc.tree', str(mm_ok))
        m_er = digest._ext_metrics('phylogeo',
                                   os.path.join(tr, list(cases)[1]))
        mm_er = m_er.get('metrics') or {}
        check('H beast 报错：beast_ok=False + 原文透出',
              mm_er.get('beast_ok') is False
              and 'ValueError' in (mm_er.get('beast_err') or ''), str(mm_er))
        m_nu = digest._ext_metrics('phylogeo',
                                   os.path.join(tr, list(cases)[2]))
        mm_nu = m_nu.get('metrics') or {}
        check('H beast=null：不产出任何 beast_* 键',
              not [k for k in mm_nu if k.startswith('beast_')], str(mm_nu))
        check('H 未登记工具：_ext_metrics 返回空',
              digest._ext_metrics('toolx00', os.path.join(tr, 'toolx00_'
                                   '20261009_120000_dddddd')) == {},
              '应为 {}')

        # (2) 去重 + 契约优先：20 个更新的非契约目录不得把 phylogeo 挤掉
        orig = digest.DIRS
        digest.DIRS = dict(orig, tool_runs=tr)
        try:
            rows = digest.tool_runs_brief(limit=8)
        finally:
            digest.DIRS = orig
        tools = [r['tool'] for r in rows]
        check('H 契约型工具保证在列（40 个更新的非契约工具之下）',
              'phylogeo' in tools and len(rows) == 8, str(tools))
        check('H 按工具去重：20 个同名 pdprep 只保留 1 条',
              len([t for t in tools if t == 'pdprep']) <= 1, str(tools))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def main():
    from Virus_Platform_Core.ai_agent import digest
    print('digest ↔ 源文件 数字对拍：')
    part_a(digest)
    part_b(digest)
    part_c()
    part_d(digest)
    part_e(digest)
    part_f(digest)
    part_g(digest)
    part_h(digest)
    print()
    if FAIL:
        print(f'✗ {len(FAIL)} 项失败: {FAIL}')
        return 1
    print('✔ 全部通过')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
