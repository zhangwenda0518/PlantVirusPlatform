# -*- coding: utf-8 -*-
"""病毒样本预筛（数据资源组独立页）。

基于 NCBI STAT / CNCB GSA 的**预计算**物种分类结果，在**不下载测序数据**
的前提下判断一批 Run 是否含病毒 —— 让「公共数据下载」只下值得下的样本。

引擎：Virus_Platform_Core/public_meta/prescreen_engine.py（CLI，子进程运行）
产物：run/meta_search/<项目>/prescreen/{stat_run_summary.tsv, stat_viral_taxa.tsv,
      prescreen_pass.tsv, unresolved.tsv, unresolved.list, STAT_PRESCREEN_REPORT.md}

设计要点
--------
* **只出清单，不干预下载**：预筛只产出通过/未决清单，下哪些由使用者决定
  （判定语义级别的门控不能替使用者定）。
* **输入三种来源**：粘贴 Run 列表 / 选一个「公共数据检索」项目（直接吃它
  的 Core14 元数据表）/ 选平台内已有的列表或 CSV 文件。
* **缓存复用**：JSON 缓存落在项目目录内，同一项目重跑近乎秒级返回；
  GSA 的 FAILED 属临时故障、不写缓存，重跑即补。
"""
import csv
import json
import os
import re

from flask import Blueprint, abort, jsonify, render_template, request, send_file

from Virus_Platform_Core.config import DIRS, engine_cmd
from Virus_Platform_Core.public_meta.virome import build_virome_all
from Virus_Platform_Core.utils import check_path
from Virus_Platform_Core.web.meta import _meta_stream_task, _meta_table_path

bp = Blueprint('prescreen', __name__)

# 产物文件名（与引擎输出一一对应；新增产物时同步这里与前端的下载清单）
ARTIFACTS = ('stat_run_summary.tsv', 'stat_viral_taxa.tsv', 'prescreen_pass.tsv',
             'unresolved.tsv', 'unresolved.list', 'STAT_PRESCREEN_REPORT.md')

_NAME_RE = re.compile(r'[A-Za-z0-9_\-]+')

# TSV 里的数值列：读回时转成数字，前端才能直接画图（否则全是字符串）
_FLOAT_COLS = ('virus_pct', 'virus_pct_ident', 'top_virus_pct',
               'identified_pct', 'top_organism_pct', 'fetch_s')
_INT_COLS = ('analysed', 'identified', 'n_taxa', 'virus_count')


# ------------------------------------------------------------------
# 路径
# ------------------------------------------------------------------
def _proj_dir(name):
    """项目目录（与「公共数据检索」同项目，便于归档与缓存复用）。"""
    return check_path(os.path.join(DIRS['meta_search'], name),
                      must_exist=False, in_platform=True)


def _out_dir(name):
    return os.path.join(_proj_dir(name), 'prescreen')


def _valid_name(name):
    return bool(name) and bool(_NAME_RE.fullmatch(name))


# ------------------------------------------------------------------
# 读结果
# ------------------------------------------------------------------
def _read_tsv(path):
    if not os.path.isfile(path):
        return []
    with open(path, 'r', encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def _coerce(row):
    """数值列转数字（空串保持 ''，不要变成 0 —— 「无数据」与「0」语义不同）。"""
    out = dict(row)
    for c in _FLOAT_COLS:
        v = out.get(c)
        if v not in (None, ''):
            try:
                out[c] = float(v)
            except ValueError:
                pass
    for c in _INT_COLS:
        v = out.get(c)
        if v not in (None, ''):
            try:
                out[c] = int(float(v))
            except ValueError:
                pass
    return out


def _stats(summary):
    ok = [r for r in summary if r.get('status') == 'OK']
    return {
        'total': len(summary),
        'ok': len(ok),
        'with_virus': sum(1 for r in ok if float(r.get('virus_pct') or 0) > 0),
        'pass': sum(1 for r in summary if r.get('pass') == 'PASS'),
        'unresolved': len(summary) - len(ok),
        'ncbi_total': sum(1 for r in summary if r.get('source') == 'NCBI_STAT'),
        'ncbi_ok': sum(1 for r in ok if r.get('source') == 'NCBI_STAT'),
        'gsa_total': sum(1 for r in summary if r.get('source') == 'GSA'),
        'gsa_ok': sum(1 for r in ok if r.get('source') == 'GSA'),
    }


# ------------------------------------------------------------------
# 页面
# ------------------------------------------------------------------
@bp.route('/prescreen')
def page_prescreen():
    return render_template('prescreen.html')


@bp.route('/api/prescreen/projects')
def api_prescreen_projects():
    """可作为输入来源的检索项目（含可用元数据表）。"""
    root = DIRS['meta_search']
    out = []
    if os.path.isdir(root):
        for n in sorted(os.listdir(root)):
            if not _valid_name(n):
                continue
            pdir = os.path.join(root, n)
            if not os.path.isdir(pdir):
                continue
            tables = {}
            for which in ('core14', 'search', 'full'):
                p = _meta_table_path(n, which)
                if p and os.path.isfile(p):
                    tables[which] = p
            has_result = os.path.isfile(
                os.path.join(pdir, 'prescreen', 'stat_run_summary.tsv'))
            out.append({
                'name': n,
                'tables': sorted(tables.keys()),
                'preferred': 'core14' if 'core14' in tables else (
                    'search' if 'search' in tables else (
                        'full' if 'full' in tables else '')),
                'has_result': has_result,
            })
    return jsonify({'projects': out})


def _detect_run_col(csv_path):
    """从表头探测 Run 列名（平台检索产出的表列名大小写不一）。"""
    try:
        with open(csv_path, 'r', encoding='utf-8-sig', errors='replace') as f:
            head = f.readline()
    except OSError:
        return None
    cols = [c.strip().strip('"') for c in head.rstrip('\n').split(',')]
    for want in ('Run', 'run', 'RUN', 'RunID', 'run_accession', 'Accession'):
        if want in cols:
            return want
    for c in cols:
        if 'run' in c.lower():
            return c
    return None


# ------------------------------------------------------------------
# 启动预筛
# ------------------------------------------------------------------
@bp.route('/api/prescreen/run', methods=['POST'])
def api_prescreen_run():
    body = request.get_json(force=True) or {}
    name = str(body.get('name') or '').strip()
    if not _valid_name(name):
        abort(400, '项目名只能含字母/数字/下划线/连字符')

    source = str(body.get('source') or 'paste').strip().lower()
    if source not in ('paste', 'project', 'file'):
        abort(400, '无效的输入来源')

    out_dir = _out_dir(name)
    os.makedirs(out_dir, exist_ok=True)

    # ---- 组装引擎输入（统一成 --runs <列表文件>）----
    cmd_args = []
    if source == 'paste':
        runs = body.get('runs') or []
        if isinstance(runs, str):
            runs = runs.splitlines()
        runs = [str(x).strip() for x in runs if str(x).strip()]
        if not runs:
            abort(400, '请至少粘贴一条 Run 编号')
        if len(runs) > 20000:
            abort(400, '单次最多 20000 条（更大批量请分批）')
        list_p = check_path(os.path.join(out_dir, 'input.list'),
                            must_exist=False, in_platform=True)
        with open(list_p, 'w', encoding='utf-8') as f:
            f.write('\n'.join(runs) + '\n')
        cmd_args += ['--runs', list_p]

    elif source == 'project':
        src_name = str(body.get('project') or '').strip()
        if not _valid_name(src_name):
            abort(400, '无效的来源项目名')
        which = str(body.get('table') or '').strip()
        if which not in ('core14', 'search', 'full'):
            which = ''
        cands = ([which] if which else []) + ['core14', 'search', 'full']
        csv_p, run_col = None, None
        for w in cands:
            p = _meta_table_path(src_name, w)
            if p and os.path.isfile(p):
                col = _detect_run_col(p)
                if col:
                    csv_p, run_col = p, col
                    break
        if not csv_p:
            abort(400, f'项目 {src_name} 下没有可用的元数据表（或表里找不到 Run 列）；'
                       f'请先在「公共数据检索」跑一次，或改用粘贴列表')
        cmd_args += ['--csv', csv_p, '--run-col', run_col]

    else:                                    # file：平台内的列表 / CSV 文件
        rel = str(body.get('path') or '').strip()
        if not rel:
            abort(400, '请指定文件路径')
        p = check_path(rel, must_exist=True, in_platform=True)
        if not os.path.isfile(p):
            abort(400, '指定路径不是文件')
        ext = os.path.splitext(p)[1].lower()
        if ext in ('.csv', '.tsv'):
            col = _detect_run_col(p) or 'Run'
            cmd_args += ['--csv', p, '--run-col', col]
        else:
            cmd_args += ['--runs', p]

    # ---- 参数校验放工厂期（提交时就报，而不是等调度起来在日志里看）----
    try:
        virus_pct = body.get('virus_pct')
        if virus_pct not in (None, ''):
            virus_pct = float(virus_pct)
            if not (0 <= virus_pct <= 100):
                abort(400, '病毒占比阈值应在 0–100 之间')
            cmd_args += ['--virus-pct', str(virus_pct)]
        jobs = int(body.get('jobs') or 6)
        if not (1 <= jobs <= 20):
            abort(400, 'NCBI 并发应在 1–20 之间（过高会被限流）')
        gsa_jobs = int(body.get('gsa_jobs') or 3)
        if not (1 <= gsa_jobs <= 6):
            abort(400, 'GSA 并发应在 1–6 之间（NGDC 有反爬）')
    except (TypeError, ValueError):
        abort(400, '参数格式错误')
    cmd_args += ['-j', str(jobs), '--gsa-jobs', str(gsa_jobs)]
    if bool(body.get('no_taxdump')):
        cmd_args.append('--no-taxdump')
    if bool(body.get('no_gsa')):
        cmd_args.append('--no-gsa')
    if bool(body.get('recheck_nodata')):
        cmd_args.append('--recheck-nodata')

    from Virus_Platform_Core.public_meta import prescreen_engine
    cmd = engine_cmd(prescreen_engine.__file__, '-o', out_dir, *cmd_args)

    tid = _meta_stream_task(
        f'病毒预筛 {name}', cmd, f'预筛完成: {out_dir}',
        log_file=os.path.join(out_dir, 'run.log'), link='/prescreen')
    return jsonify({'task': tid, 'name': name})


# ------------------------------------------------------------------
# 读结果
# ------------------------------------------------------------------
@bp.route('/api/prescreen/result')
def api_prescreen_result():
    name = str(request.args.get('name') or '').strip()
    if not _valid_name(name):
        abort(400, '无效的项目名')
    out_dir = _out_dir(name)
    sum_p = os.path.join(out_dir, 'stat_run_summary.tsv')
    if not os.path.isfile(sum_p):
        return jsonify({'exists': False, 'name': name})

    summary = [_coerce(r) for r in _read_tsv(sum_p)]
    unresolved = [_coerce(r) for r in
                  _read_tsv(os.path.join(out_dir, 'unresolved.tsv'))]
    passed = [_coerce(r) for r in
              _read_tsv(os.path.join(out_dir, 'prescreen_pass.tsv'))]

    files = []
    for fn in ARTIFACTS:
        fp = os.path.join(out_dir, fn)
        if os.path.isfile(fp):
            files.append({'name': fn, 'size': os.path.getsize(fp)})

    return jsonify({
        'exists': True,
        'name': name,
        'out_dir': out_dir,
        'stats': _stats(summary),
        'summary': summary,
        'unresolved': unresolved,
        'passed': passed,
        'files': files,
    })


@bp.route('/api/prescreen/viral')
def api_prescreen_viral():
    """某个 run 的病毒 taxon 明细（按需拉，不随结果整包返回 —— 明细可达万行）。"""
    name = str(request.args.get('name') or '').strip()
    acc = str(request.args.get('acc') or '').strip()
    if not _valid_name(name) or not re.fullmatch(r'[A-Za-z0-9_\-\.]+', acc or ''):
        abort(400, '无效参数')
    rows = _read_tsv(os.path.join(_out_dir(name), 'stat_viral_taxa.tsv'))
    return jsonify({'acc': acc, 'rows': [_coerce(r) for r in rows
                                         if r.get('acc') == acc]})


# ------------------------------------------------------------------
# 病毒谱（病毒视角）
# ------------------------------------------------------------------
# 进程内缓存：name -> {'key': ..., 'data': ...}
_VIR_CACHE = {}

# 缓存文件名（放 .cache 子目录，避免和交付产物混在一起）
_VIR_CACHE_FILE = os.path.join('.cache', 'virome.json')


def _virome_key(out_dir):
    """缓存指纹：两份输入产物 + taxdump 的 mtime/size。

    taxdump 也要进指纹 —— 平台更新分类库后，rank/谱系会变，
    不重算就会拿着旧分类展示新结果。
    """
    parts = []
    for p in (os.path.join(out_dir, 'stat_viral_taxa.tsv'),
              os.path.join(out_dir, 'stat_run_summary.tsv'),
              os.path.join(DIRS['taxonomy'], 'nodes.dmp'),
              os.path.join(DIRS['taxonomy'], 'names.dmp')):
        try:
            st = os.stat(p)
            parts.append(f'{st.st_mtime_ns}:{st.st_size}')
        except OSError:
            parts.append('-')
    return '|'.join(parts)


@bp.route('/api/prescreen/virome')
def api_prescreen_virome():
    """病毒谱：一次算完所有层级，返回 {level: 结果}。

    ⚠ 首次计算很慢（15~45 秒，瓶颈是读 300 万行 nodes.dmp + 470 万行
    names.dmp，与数据量无关）。故做两层缓存：
      * 进程内 `_VIR_CACHE`（同一项目的重复访问瞬时返回）
      * 项目目录内 `.cache/virome.json`（跨重启、跨会话）
    前端在结果加载后会自动预取一次，用户点开「病毒谱」时通常已就绪。
    """
    name = str(request.args.get('name') or '').strip()
    if not _valid_name(name):
        abort(400, '无效的项目名')
    out_dir = _out_dir(name)
    taxa_p = os.path.join(out_dir, 'stat_viral_taxa.tsv')
    sum_p = os.path.join(out_dir, 'stat_run_summary.tsv')
    if not os.path.isfile(taxa_p) or not os.path.isfile(sum_p):
        return jsonify({'exists': False, 'name': name})

    key = _virome_key(out_dir)

    mem = _VIR_CACHE.get(name)
    if mem and mem.get('key') == key:
        return jsonify(mem['data'])

    cache_p = os.path.join(out_dir, _VIR_CACHE_FILE)
    try:
        with open(cache_p, 'r', encoding='utf-8') as f:
            obj = json.load(f)
        if isinstance(obj, dict) and obj.get('key') == key:
            _VIR_CACHE[name] = obj
            return jsonify(obj['data'])
    except (OSError, ValueError):
        pass

    data = build_virome_all(taxa_p, sum_p)
    payload = {'key': key, 'data': data}
    try:
        os.makedirs(os.path.dirname(cache_p), exist_ok=True)
        with open(cache_p, 'w', encoding='utf-8') as f:
            json.dump(payload, f, ensure_ascii=False)
    except OSError:
        pass                                          # 缓存写不了不影响出结果
    _VIR_CACHE[name] = payload
    return jsonify(data)


@bp.route('/prescreen_files/<name>/<path:subpath>')
def page_prescreen_file(name, subpath):
    """产物文件下载/预览（严格限制在该项目 prescreen 目录内）。"""
    if not _valid_name(name) or '..' in subpath:
        abort(400, '无效路径')
    # 点开头的目录/文件是内部缓存（如 .cache/virome.json），不是交付产物
    if any(seg.startswith('.') for seg in subpath.replace('\\', '/').split('/')):
        abort(404, '无此产物')
    base = check_path(_out_dir(name), must_exist=True, in_platform=True)
    p = check_path(os.path.join(base, subpath), must_exist=True, in_platform=True)
    if not p.startswith(base + os.sep):
        abort(400, '路径越界')
    return send_file(p, as_attachment=request.args.get('dl') == '1')
