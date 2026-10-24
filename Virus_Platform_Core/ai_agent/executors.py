# -*- coding: utf-8 -*-
"""AI 执行层：分级执行权（有边界的"放开"）。

    第一级 执行单（默认，永远在）：AI 只能**起草**动作（create_action），
            产生待批卡片（对话坞里看参数 diff 与理由），用户批准后才经
            **与界面按钮完全相同的代码路径**执行；驳回即作废。
    第二级 自动重试（cfg.ai['exec_auto']，默认关）：阶段失败后原参数重跑，
            不经确认——因为"重试"就是平台断点续跑机制本身，AI 只是替你
            按了重试按钮。每样品每阶段限 1 次，杜绝无限循环。
    第三级 自动调参闭环（cfg.ai['exec_auto_tune']，默认关，需手动开启）：
            失败后由 AI 提议参数（限定键 + 数值边界钳制）并自动执行，
            每样品限 1 次。调参失败回落第二级（若其开启）。

永远不做的事（注册表里根本没有这些动作）：
    删除样品/数据、全局重置、修改 AI 配置自身、修改设置中心、
    平台目录之外的任何写操作。

安全闸（全部硬编码，不受配置影响）：
    - 同一样品同时只允许一个运行（复用 _sample_busy_register）；
    - AI 发起的运行全平台并发上限 MAX_AI_RUNS=1；
    - 执行单参数走白名单 + 类型钳制，未知键直接拒绝；
    - 全部动作（提案/批准/驳回/自动执行/结果）发布到事件流**并异步落盘**
      （`ai_agent/audit.py`，按天轮转 JSONL）；重启后 `list_actions` 从盘上
      重建历史，进行中的状态标为 `interrupted`，且 aid 计数器会推高防撞号。
"""
import json
import threading
import time

from Virus_Platform_Core.config import DIRS, get_config

# ------------------------------------------------------------------
# 动作注册表（白名单：没有登记的动作 AI 永远无法发起）
# ------------------------------------------------------------------
ACTIONS = {
    'rerun_stage': {
        'title': ('重跑指定阶段', 'Re-run stage(s)'),
        'risk': 'confirm',        # Tier1 确认；Tier3 开启时可自动（仍受预算）
    },
    'retry_failed': {
        'title': ('原参数重试失败阶段', 'Retry failed stage(s), same params'),
        'risk': 'auto_ok',        # 可进入第二级白名单
    },
    'download_runs': {
        'title': ('下载公共数据', 'Download public data'),
        'risk': 'confirm',        # 外网批量下载，永远人工确认（不走自动级）
    },
    'run_tool': {
        'title': ('运行工具箱工具', 'Run toolbox tool'),
        'risk': 'confirm',        # 对齐 SDT/RDP/引物/比对等 ~30 个工具箱能力
    },
    'create_sample': {
        'title': ('创建样品', 'Create sample'),
        'risk': 'confirm',        # 从已落盘的 FASTQ 建样品（不下载、不联网）
    },
}

# 允许 AI 调参的键与边界（第三级）；之外的键一律拒绝
TUNABLE_PARAMS = {
    'assembly_mode': ('rnaviral', 'metaviral', 'rna', 'meta', 'isolate', 'srna'),
    'min_contig_len': (50, 5000),
    'confidence': (0.0, 1.0),
    'subsample': (0, 5_000_000),
    'subsample_mode': ('prefix', 'random', 'proportion'),
    'subsample_prop': (0.0, 1.0),
    'subsample_seed': (1, 2_147_483_647),
    'min_orf_aa': (30, 1000),
    'top_n_refs': (5, 50),
}

# run_analysis 参数装配请求体里合法的键（白名单之外的提案直接拒）
BODY_PARAM_KEYS = {
    'assembly_mode', 'assembly_input', 'confidence', 'threads', 'memory',
    'subsample', 'subsample_mode', 'subsample_prop', 'subsample_seed',
    'min_contig_len', 'min_orf_aa', 'top_n_refs', 'tree_tool',
    'tree_sampling', 'primer_mode', 'do_trim', 'specificity', 'force',
    'fastp_dedup', 'do_fq2fa', 'gbdraw_max', 'plot_engine', 'kv_lib',
    'db_host', 'db_virus', 'chunk_dir', 'verify_host', 'verify_combine',
}

MAX_PENDING = 20                 # 待批执行单上限
MAX_AI_RUNS = 1                  # AI 发起的运行全平台并发上限
AUTO_RETRY_MAX_PER_STAGE = 1     # 第二级：每样品每阶段自动重试次数
AUTO_TUNE_MAX_PER_SAMPLE = 1     # 第三级：每样品自动调参次数

# ------------------------------------------------------------------
# 待决动作存储 + 预算
# ------------------------------------------------------------------
_lock = threading.Lock()
_actions = {}                    # aid -> rec（含全部历史状态）
_order_n = 0                     # 已发出的 aid 序号（重启后从审计盘推高）
_auto_retried = {}               # (sample, stage) -> 次数
_auto_tuned = {}                 # sample -> 次数
_hydrated = False                # 是否已从审计盘重建过历史

# 内存里保留的历史条数上限（审计盘才是长期留痕，内存只为界面列表）
MAX_KEEP_ACTIONS = 500


def _ai_cfg():
    return get_config().ai or {}


def _stage_disp(stage):
    try:
        from Virus_Platform_Core.pipeline import STAGE_NAMES, STAGE_NAMES_EN
        lang = get_config().lang
        return ((STAGE_NAMES_EN if lang == 'en' else STAGE_NAMES)
                .get(stage) or STAGE_NAMES.get(stage) or stage)
    except Exception:
        return stage


def _publish(rec, event='ai_action'):
    """执行单 → 事件总线。

    ⚠ 执行单的业务 id（'act_N' 字符串）以 **aid** 字段发布——绝不能叫
    'id'：事件总线自己的整型 id 是客户端去重/水位线的根基，被覆盖会把
    events_since 的 int 比较打成 TypeError，整条 SSE 流持续 500（实测踩过）。
    """
    from . import broadcast
    payload = {k: v for k, v in rec.items()
               if k in ('action', 'title', 'sample', 'stages', 'params',
                        'base_params', 'reason', 'status', 'source',
                        'error', 'tid', 'runs', 'tool', 'tool_params')}
    if isinstance(payload.get('runs'), list):
        payload['runs'] = payload['runs'][:10]   # 事件里只带前 10 条防膨胀
        payload['runs_n'] = rec.get('params', {}).get('runs_n')
    payload['aid'] = rec['id']
    broadcast.publish(event, **payload)


def _snapshot(rec):
    return {k: rec.get(k) for k in
            ('id', 'action', 'title', 'sample', 'stages', 'params',
             'base_params', 'reason', 'status', 'source', 'error', 'tid',
             'ts', 'tool', 'tool_params', 'runs')}


def _audit(rec):
    """执行单状态 → 审计盘。

    用 `record_sync`（**同步**写）而不是异步队列：状态变迁低频高价值，
    而 Windows 的 `TerminateProcess` 不可捕获，异步队列的 2s 窗口可能把
    "谁批的、执行没执行"整段丢掉。
    """
    try:
        from . import audit
        audit.record_sync('actions', rec)
    except Exception:
        pass


def _hydrate():
    """从审计盘重建历史执行单（**重启后 `list_actions` 仍能看到过去**）。

    只回灌内存里没有的 id；`pending/executing/launched` 这些"进行中"状态在
    重启后已失去执行线程，标为 `interrupted`（诚实反映"没跑完就没了"，而不是
    假装还能被批准）。同时把 aid 计数器推到已用最大值之后，**避免重启后
    `act_1` 撞号**（撞号会让新单覆盖旧记录）。

    ⚠ 必须在**取 `_lock` 之前**调用，否则与本函数内部的 `_lock` 自锁。
    """
    global _hydrated, _order_n
    if _hydrated:
        return
    with _lock:
        if _hydrated:
            return
        _hydrated = True
    try:
        from . import audit
        recs = audit.read('actions', limit=MAX_KEEP_ACTIONS * 2)
    except Exception:
        return
    if not recs:
        return
    mx = 0
    with _lock:
        for r in recs:                       # read() 返回新→旧
            aid = r.get('id')
            if not isinstance(aid, str) or not aid.startswith('act_'):
                continue
            try:
                mx = max(mx, int(aid.split('_', 1)[1]))
            except (IndexError, ValueError):
                continue
            if aid in _actions:
                continue
            r = dict(r)
            if r.get('status') in ('pending', 'executing', 'launched'):
                r['status'] = 'interrupted'
            _actions[aid] = r
        _order_n = max(_order_n, mx)
        # 内存只为界面列表服务：超量时按 ts 丢最旧
        if len(_actions) > MAX_KEEP_ACTIONS:
            keep = sorted(_actions.values(), key=lambda x: -(x.get('ts') or 0))
            _actions.clear()
            for r in keep[:MAX_KEEP_ACTIONS]:
                _actions[r['id']] = r


def list_actions(limit=30):
    _hydrate()
    with _lock:
        recs = sorted(_actions.values(), key=lambda r: -r['ts'])
        return [_snapshot(r) for r in recs[:limit]]


def get_action(aid):
    _hydrate()
    with _lock:
        rec = _actions.get(aid)
        return _snapshot(rec) if rec else None


# ------------------------------------------------------------------
# 参数校验（白名单 + 类型钳制；未知键直接拒绝）
# ------------------------------------------------------------------
def validate_params(action, params):
    """返回 (normalized_params, error)。error 非 None 表示提案非法。

    白名单 = BODY_PARAM_KEYS；TUNABLE_PARAMS 里的键额外做取值边界钳制
    （第三级自动调参的护栏，对对话提案同样生效）。
    """
    if params is None:
        return {}, None
    if not isinstance(params, dict):
        return {}, 'params 须为对象'
    bad = [k for k in params if k not in BODY_PARAM_KEYS]
    if bad:
        return {}, f'参数不在允许清单内: {", ".join(map(str, bad))}'
    out = {}
    for k, v in params.items():
        if v is None or v == '':
            continue
        if k in ('assembly_mode', 'subsample_mode'):
            if str(v) not in TUNABLE_PARAMS[k]:
                return {}, f'{k} 只允许: {", ".join(TUNABLE_PARAMS[k])}'
            out[k] = str(v)
            continue
        if k in ('do_trim', 'specificity', 'force', 'fastp_dedup',
                 'do_fq2fa'):
            out[k] = bool(v)
            continue
        if k in ('confidence', 'subsample_prop'):
            try:
                num = float(v)
            except (TypeError, ValueError):
                return {}, f'{k} 须为数字'
            if k in TUNABLE_PARAMS:
                lo, hi = TUNABLE_PARAMS[k]
                num = min(hi, max(lo, num))
            out[k] = round(num, 3)
            continue
        if k in ('primer_mode', 'tree_tool', 'tree_sampling',
                 'assembly_input', 'plot_engine', 'verify_host',
                 'verify_combine', 'kv_lib', 'db_host', 'db_virus',
                 'chunk_dir'):
            out[k] = str(v)
            continue
        # 其余全部是整数键
        try:
            num = int(v)
        except (TypeError, ValueError):
            return {}, f'{k} 须为整数'
        if k in TUNABLE_PARAMS:
            lo, hi = TUNABLE_PARAMS[k]
            num = int(min(hi, max(lo, num)))
        out[k] = num
    return out, None


# ------------------------------------------------------------------
# 创建 / 决策
# ------------------------------------------------------------------
def create_action(action, sample, stages=None, params=None, reason='',
                  source='chat', base_params=None):
    """AI（或代码）起草一张执行单。返回 (rec, error)。"""
    if action not in ACTIONS:
        return None, f'未知动作: {action}（允许: {", ".join(ACTIONS)}）'
    # 先重建历史：否则重启后第一张单会拿到 act_1，与盘上旧记录撞号
    _hydrate()
    if action == 'download_runs':
        return _create_download(params, reason, source)
    if action == 'run_tool':
        return _create_run_tool(params, reason, source)
    if action == 'create_sample':
        return _create_sample_order(params, reason, source)
    sample = str(sample or '').strip()
    if not sample:
        return None, '需要 sample'
    from Virus_Platform_Core.utils import resolve_sample_name
    try:
        real = resolve_sample_name(sample, DIRS['results'])
        sd = __import__('os').path.join(DIRS['results'], real)
        if not __import__('os').path.isdir(sd):
            return None, f'样品不存在: {sample}'
    except Exception as e:
        return None, f'样品解析失败: {e}'
    if stages is not None:
        from Virus_Platform_Core.pipeline import STAGE_ORDER
        stages = [str(s) for s in stages if str(s) in STAGE_ORDER]
        if not stages:
            return None, 'stages 为空或全部非法'
    params, err = validate_params(action, params)
    if err:
        return None, err
    with _lock:
        pending = sum(1 for r in _actions.values() if r['status'] == 'pending')
        if pending >= MAX_PENDING:
            return None, f'待批执行单已达上限（{MAX_PENDING}），请先处理'
        global _order_n
        _order_n += 1
        aid = f'act_{_order_n}'
        lang = get_config().lang
        rec = {'id': aid, 'action': action,
               'title': ACTIONS[action]['title'][0 if lang == 'zh' else 1],
               'sample': real, 'stages': stages, 'params': params,
               'base_params': base_params or {},
               'reason': str(reason or '')[:800], 'status': 'pending',
               'source': source, 'error': '', 'tid': '', 'ts': time.time()}
        _actions[aid] = rec
    _snapshot_publish(rec, 'pending')
    return rec, None


def _create_download(params, reason, source):
    """download_runs 执行单：runs 编号清单 → 待批卡片。

    校验：runs 必须是 1~500 条非空编号（SRR/ERR/DRR/CRR 或 URL——合法性
    由下载管理器在执行时逐条解析，这里只挡空单与超量）。
    """
    params = params if isinstance(params, dict) else {}
    runs = [str(x).strip() for x in (params.get('runs') or [])
            if str(x).strip()]
    if not runs:
        return None, "download_runs 需要 params.runs（run 编号列表，1~500 条）"
    if len(runs) > 500:
        return None, '单批最多 500 条，请拆分'
    name = str(params.get('name') or '').strip() or 'AI 下载'
    convert_sra = bool(params.get('convert_sra', True))
    with _lock:
        pending = sum(1 for r in _actions.values() if r['status'] == 'pending')
        if pending >= MAX_PENDING:
            return None, f'待批执行单已达上限（{MAX_PENDING}），请先处理'
        global _order_n
        _order_n += 1
        aid = f'act_{_order_n}'
        lang = get_config().lang
        rec = {'id': aid, 'action': 'download_runs',
               'title': ACTIONS['download_runs']['title'][0 if lang == 'zh'
                                                          else 1],
               'sample': '', 'stages': None,
               'params': {'runs_n': len(runs), 'name': name,
                          'convert_sra': convert_sra},
               'base_params': {}, 'runs': runs,
               'reason': str(reason or '')[:800], 'status': 'pending',
               'source': source, 'error': '', 'tid': '', 'ts': time.time()}
        _actions[aid] = rec
    _snapshot_publish(rec, 'pending')
    return rec, None


def _create_run_tool(params, reason, source):
    """run_tool 执行单：把工具箱（TOOL_REGISTRY）里的任一工具作为待批动作。

    参数即该工具在界面上的表单键值（字符串口径，与 /api/tool/run 的
    body.params 一致）；文件存在性等校验在执行期由 job 工厂做——失败会
    如实变成 failed 状态与错误信息，AI 可据此修正后重新提案。
    """
    params = params if isinstance(params, dict) else {}
    tool = str(params.get('tool') or '').strip()
    try:
        from Virus_Platform_Core.web.tools_api import TOOL_REGISTRY
        spec = TOOL_REGISTRY.get(tool)
    except Exception as e:
        return None, f'工具注册表不可用: {e}'
    if not spec:
        return None, (f'未知工具: {tool or "(空)"}；'
                      '先用 list_platform_tools 查看可用工具')
    tool_params = {str(k): v for k, v in params.items() if k != 'tool'}
    with _lock:
        pending = sum(1 for r in _actions.values() if r['status'] == 'pending')
        if pending >= MAX_PENDING:
            return None, f'待批执行单已达上限（{MAX_PENDING}），请先处理'
        global _order_n
        _order_n += 1
        aid = f'act_{_order_n}'
        lang = get_config().lang
        rec = {'id': aid, 'action': 'run_tool',
               'title': (spec['title'] if lang == 'zh' else
                         spec.get('title_en') or spec['title']),
               'sample': '', 'stages': None,
               'params': {'tool': tool, **tool_params},
               'tool': tool, 'tool_params': tool_params,
               'base_params': {},
               'reason': str(reason or '')[:800], 'status': 'pending',
               'source': source, 'error': '', 'tid': '', 'ts': time.time()}
        _actions[aid] = rec
    _snapshot_publish(rec, 'pending')
    return rec, None


def _create_sample_order(params, reason, source):
    """create_sample 执行单：从已落盘的 FASTQ 建样品（不联网、不下载）。"""
    params = params if isinstance(params, dict) else {}
    r1 = str(params.get('r1') or '').strip()
    r2 = str(params.get('r2') or '').strip()
    sample = str(params.get('sample') or '').strip()
    if not r1:
        return None, 'create_sample 需要 params.r1（FASTQ 绝对路径）'
    import os
    for k, v in (('r1', r1), ('r2', r2)):
        if v and not os.path.isfile(v):
            return None, f'params.{k} 文件不存在: {v}'
    if not sample:
        sample = os.path.splitext(os.path.basename(r1))[0]
    with _lock:
        pending = sum(1 for r in _actions.values() if r['status'] == 'pending')
        if pending >= MAX_PENDING:
            return None, f'待批执行单已达上限（{MAX_PENDING}），请先处理'
        global _order_n
        _order_n += 1
        aid = f'act_{_order_n}'
        lang = get_config().lang
        rec = {'id': aid, 'action': 'create_sample',
               'title': ACTIONS['create_sample']['title'][0 if lang == 'zh'
                                                          else 1],
               'sample': sample, 'stages': None,
               'params': {'sample': sample, 'r1': r1, 'r2': r2,
                          'project': str(params.get('project') or '')},
               'base_params': {},
               'reason': str(reason or '')[:800], 'status': 'pending',
               'source': source, 'error': '', 'tid': '', 'ts': time.time()}
        _actions[aid] = rec
    _snapshot_publish(rec, 'pending')
    return rec, None


def _snapshot_publish(rec, status):
    """状态迁移的唯一出口：改状态 → 发事件 → 写审计盘。"""
    rec['status'] = status
    _publish(rec)
    _audit(rec)


def resolve(aid, decision):
    """批准 / 驳回。返回 (rec, error)。decision: approve | reject。"""
    _hydrate()
    with _lock:
        rec = _actions.get(aid)
        if not rec:
            return None, '执行单不存在'
        if rec['status'] != 'pending':
            return None, f'执行单已处理（{rec["status"]}）'
    if decision == 'reject':
        _snapshot_publish(rec, 'rejected')
        return rec, None
    if decision != 'approve':
        return None, 'decision 只允许 approve / reject'
    # 并发闸：AI 发起的运行全平台同时最多 MAX_AI_RUNS 个
    if _ai_active_count() >= MAX_AI_RUNS:
        return None, (f'AI 发起的运行已达并发上限（{MAX_AI_RUNS}），'
                      f'请等当前运行结束')
    _snapshot_publish(rec, 'executing')
    threading.Thread(target=_execute, args=(rec,), daemon=True,
                     name=f'ai-exec-{aid}').start()
    return rec, None


def _execute(rec):
    """在后台线程把执行单落到真实任务（与界面按钮同一通道）。"""
    try:
        act = rec.get('action')
        if act == 'download_runs':
            tid = _launch_download(rec)
        elif act == 'run_tool':
            tid = _launch_tool(rec.get('tool'), rec.get('tool_params'))
        elif act == 'create_sample':
            tid = _launch_create_sample(rec.get('params') or {})
        else:
            tid = _launch_analysis(rec['sample'], rec.get('stages'),
                                   rec.get('params') or {})
        rec['tid'] = str(tid)
        _snapshot_publish(rec, 'launched')
    except Exception as e:
        rec['error'] = f'{type(e).__name__}: {e}'[:400]
        _snapshot_publish(rec, 'failed')


def _launch_download(rec):
    """启动一批公共数据下载（与「公共数据下载」页同一管理器/同一队列）。"""
    from Virus_Platform_Core.web.download import get_dl_manager
    params = rec.get('params') or {}
    return get_dl_manager().create(
        str(params.get('name') or 'AI 下载'),
        rec.get('runs') or [],
        concurrency=2,
        convert_sra=bool(params.get('convert_sra', True)))


def _launch_tool(tool, tool_params):
    """运行工具箱工具——与「工具箱」页 /api/tool/run 完全同一套
    「TOOL_REGISTRY job 工厂 + SimpleNamespace ctx + tm.start」路径。

    tool_params 是界面表单口径的键值（值统一转字符串）；缺必填参数时
    job 工厂 abort/raise → 执行单转 failed 并带可读错误，AI 可据此修正
    后重新提案。
    """
    import os
    import uuid
    from types import SimpleNamespace
    from Virus_Platform_Core.config import PLATFORM_ROOT
    from Virus_Platform_Core.utils import check_path, safe_open
    from Virus_Platform_Core.web.state import tool_runs_root as \
        _tool_runs_root, RUN_OUT_LINK
    from Virus_Platform_Core.web.tools_api import TOOL_REGISTRY, \
        _resolve_virus_db, LIGHT_TOOLS
    from Virus_Platform_Core.web.tasks import tm

    spec = TOOL_REGISTRY.get(tool)
    if not spec:
        raise RuntimeError(f'未知工具: {tool}')
    # 表单口径：值统一字符串（job 工厂按字符串读取/strip）
    p = {}
    for k, v in (tool_params or {}).items():
        if v is None:
            continue
        p[str(k)] = v if isinstance(v, str) else json.dumps(
            v, ensure_ascii=False) if isinstance(v, (list, dict)) else str(v)

    def _input_abs(v):
        return check_path(v if os.path.isabs(str(v))
                          else os.path.join(PLATFORM_ROOT, str(v)),
                          must_exist=True)

    def _req(key, what):
        v = str(p.get(key) or '').strip()
        if not v:
            raise RuntimeError(f'缺少必填参数: {what}（{key}）')
        return _input_abs(v)

    def _opt(key):
        v = str(p.get(key) or '').strip()
        return _input_abs(v) if v else None

    db_virus = _resolve_virus_db(p.get('db_virus'))
    if spec.get('need_db'):
        from Virus_Platform_Core.kunpeng import db_ready
        if not db_ready(db_virus):
            raise RuntimeError(f'病毒库未就绪: {db_virus}'
                               '（先到「数据库构建」页构建）')

    ts = time.strftime('%Y%m%d_%H%M%S')
    run_name = f'{tool}_{ts}_{uuid.uuid4().hex[:6]}'
    base_dir = check_path(os.path.join(_tool_runs_root(), run_name),
                          must_exist=False, in_platform=True)
    run_dir = base_dir
    os.makedirs(run_dir, exist_ok=True)
    threads = None
    if p.get('threads'):
        try:
            threads = max(1, min(int(p['threads']),
                                 (os.cpu_count() or 4) * 4))
        except (TypeError, ValueError):
            threads = None
    ctx = SimpleNamespace(p=p, run_dir=run_dir, threads=threads,
                          db_virus=db_virus, req=_req, opt=_opt)
    job = spec['job'](ctx)          # 缺参/文件缺失在这里抛错
    name = get_config().tr(f'AI执行·工具·{spec["title"]} {ts}',
                           f'AI tool·{spec.get("title_en", tool)} {ts}')

    def _job_with_run(log, prog, cancel):
        res = job(log, prog, cancel)
        if isinstance(res, dict):
            res.setdefault('run', run_name)
            res['out_dir'] = str(run_dir)
        try:
            os.utime(base_dir)
        except OSError:
            pass
        return res

    return tm.start(name, _job_with_run,
                    log_file=os.path.join(run_dir, 'run.log'),
                    weight=('light' if tool in LIGHT_TOOLS else 'heavy'),
                    out_dir=str(run_dir))


def _launch_create_sample(params):
    """创建样品（与 /api/pipeline/create 同一套校验与落盘）。"""
    import os as _os
    from Virus_Platform_Core.utils import check_path, safe_sample_name
    from Virus_Platform_Core.pipeline import save_sample_input
    from Virus_Platform_Core.web.samples import _default_sample_name
    r1 = str(params.get('r1') or '').strip()
    r2 = str(params.get('r2') or '').strip()
    typed = str(params.get('sample') or '').strip()
    check_path(r1, must_exist=True)
    if r2:
        check_path(r2, must_exist=True)
    sample = safe_sample_name(typed or _default_sample_name(r1))
    sd = check_path(_os.path.join(DIRS['results'], sample),
                    must_exist=False, in_platform=True)
    if _os.path.isdir(sd) and any(_os.scandir(sd)):
        raise RuntimeError(f'样品 {sample} 已存在，请换一个样品名')
    _os.makedirs(sd, exist_ok=True)
    save_sample_input(sd, r1, r2 or None, sample,
                      project=str(params.get('project') or '') or None)
    return sample


# ------------------------------------------------------------------
# 执行器：复用与界面按钮完全相同的启动路径
# ------------------------------------------------------------------
def _ai_active_count():
    try:
        from Virus_Platform_Core.web.tasks import tm
        return sum(1 for t in (tm.list_all(logs=False) or [])
                   if t.get('status') == 'running'
                   and str(t.get('name') or '').startswith('AI执行·'))
    except Exception:
        return 0


def _launch_analysis(sample, stages, body):
    """启动一次样品分析（run_analysis），返回 tm 任务 id。

    注意：这里与 /api/analyze、批处理队列走的是同一套
    「_sample_busy_register + tm.start + TaskLogger + run_analysis」；
    body 是 _analysis_kwargs 口径的参数请求体。
    """
    import os
    from Virus_Platform_Core.utils import (check_path, resolve_sample_name,
                                           TaskLogger)
    from Virus_Platform_Core.web.samples import (_analysis_kwargs,
                                                 _sample_busy_register,
                                                 _sample_busy_release)
    from Virus_Platform_Core.web.tasks import tm

    s = resolve_sample_name(sample, DIRS['results'])
    sd = check_path(os.path.join(DIRS['results'], s), must_exist=True,
                    in_platform=True)
    from Virus_Platform_Core.pipeline import (load_sample_input,
                                              DEFAULT_ANALYZE_STAGES)
    r1, r2, _proj = load_sample_input(sd)
    if not r1:
        raise RuntimeError(f'样品 {s} 缺少输入记录（input.json），无法运行')
    check_path(r1, must_exist=True)
    if not _sample_busy_register(s):
        raise RuntimeError(f'样品 {s} 已有任务在运行')
    st = [x for x in (stages or DEFAULT_ANALYZE_STAGES)]
    kwargs = _analysis_kwargs(dict(body or {}))

    def job(log, prog, cancel):
        from Virus_Platform_Core.pipeline import run_analysis
        logger = TaskLogger(callback=log)
        try:
            return run_analysis(sample=s, r1=r1, r2=r2, stages=st,
                                logger=logger, progress=prog, **kwargs)
        finally:
            _sample_busy_release(s)
            logger.close()

    log_dir = check_path(os.path.join(sd, 'logs'), must_exist=False,
                         in_platform=True)
    os.makedirs(log_dir, exist_ok=True)
    try:
        tid = tm.start(get_config().tr(f'AI执行·{s}', f'AI run·{s}'), job,
                       log_file=os.path.join(
                           log_dir,
                           f'ai_run_{time.strftime("%Y%m%d_%H%M%S")}.log'))
    except BaseException:
        _sample_busy_release(s)      # 任务没起来也要放掉占位
        raise
    return tid


# ------------------------------------------------------------------
# 第二级 / 第三级：失败后的自动决策
# ------------------------------------------------------------------
def after_failure(sample, stage, error=''):
    """broadcast.note_stage_failed 的执行钩子：按开关走 Tier2/Tier3。

    Tier3（自动调参）优先；未开则 Tier2（原参数重试）；都关则什么都不做。
    Tier3 需要 AI 调参（秒级~几十秒），放后台线程。
    """
    ai = _ai_cfg()
    if not ai.get('enabled'):
        return
    if ai.get('exec_auto_tune'):
        threading.Thread(target=_auto_tune_run, args=(sample, stage, error),
                         daemon=True,
                         name=f'ai-tune-{sample}-{stage}').start()
    elif ai.get('exec_auto'):
        _auto_retry(sample, stage)


def _auto_retry(sample, stage):
    """第二级：原参数重跑失败阶段（每样品每阶段限 1 次）。"""
    with _lock:
        n = _auto_retried.get((sample, stage), 0)
        if n >= AUTO_RETRY_MAX_PER_STAGE:
            return None
        _auto_retried[(sample, stage)] = n + 1
    try:
        base, stages = _last_run_params(sample)
        rec, err = create_action('retry_failed', sample, stages=stages,
                                 params=None, reason='阶段失败，自动重试'
                                 '（第二级，原参数）', source='auto2',
                                 base_params=base)
        if err or not rec:
            return None
        rec2, err2 = resolve(rec['id'], 'approve')
        if err2:
            # 并发满等瞬态原因：执行单留在 pending 等人处理，不强行改状态
            import sys
            sys.stderr.write(f'[ai_agent.executors] 自动重试未启动: {err2}\n')
            return None
        return rec2
    except Exception as e:
        import sys
        sys.stderr.write(f'[ai_agent.executors] 自动重试失败（静默）: '
                         f'{type(e).__name__}: {e}\n')
        return None


def _auto_tune_run(sample, stage, error):
    """第三级：AI 提议参数 → 自动执行（每样品限 1 次；失败回落第二级）。"""
    ai = _ai_cfg()
    with _lock:
        n = _auto_tuned.get(sample, 0)
        if n >= AUTO_TUNE_MAX_PER_SAMPLE:
            return
        _auto_tuned[sample] = n + 1
    try:
        base, stages = _last_run_params(sample)
        proposal = _ai_propose_params(sample, stage, error, base)
        if not proposal:
            if ai.get('exec_auto'):
                _auto_retry(sample, stage)       # 调参不成回落原参数重试
            return
        rec, err = create_action('rerun_stage', sample, stages=stages,
                                 params=proposal.get('params'),
                                 reason=proposal.get('reason', ''),
                                 source='auto3', base_params=base)
        if err or not rec:
            return
        resolve(rec['id'], 'approve')
    except Exception as e:
        import sys
        sys.stderr.write(f'[ai_agent.executors] 自动调参失败（静默）: '
                         f'{type(e).__name__}: {e}\n')


def _ai_propose_params(sample, stage, error, base):
    """让 AI 在 TUNABLE_PARAMS 边界内提参；输出必须是 JSON，非法即放弃。"""
    from . import digest
    from .provider import ai_ready
    ai = _ai_cfg()
    if not ai_ready(ai):
        return None
    d = digest.sample_digest(sample)
    if d is None:
        return None
    dj = json.dumps({'stage': stage, 'error': str(error)[:400],
                     'current_params': base, 'digest': d},
                    ensure_ascii=False)
    if len(dj) > 12000:
        dj = dj[:12000] + '…(truncated)'
    prompt = (
        '你是植物病毒分析平台的自动调参器。一个阶段失败了，请提出**重跑参数'
        '**。只允许这些键（含义与边界）：assembly_mode（rnaviral/metaviral/'
        'rna/meta/isolate/srna）、min_contig_len（50~5000）、confidence'
        '（0~1）、subsample（0~5000000）、subsample_mode（prefix/random/'
        'proportion：前 N 对/等概率随机 N 对/按比例随机）、subsample_prop'
        '（0~1，仅 proportion 用，如 0.5=抽 50%，条数精确）、subsample_seed'
        '（1~2147483647，random/proportion 复现用）、min_orf_aa（30~1000）、'
        'top_n_refs（5~50）。只改确有依据的键，其余不要出现。\n'
        '输出**纯 JSON**（无 markdown 围栏）：{"params": {...}, '
        '"reason": "≤120字"}。失败原因与当前参数、结果摘要：\n'
        + dj)
    from .provider import chat_stream
    out = []
    try:
        for ev in chat_stream(ai, [{'role': 'user', 'content': prompt}],
                              tools=None, max_tokens=400, timeout=(15, 60)):
            if ev.get('type') == 'delta':
                out.append(ev['content'])
    except Exception:
        return None
    text = ''.join(out).strip()
    # 剥掉可能的围栏
    if text.startswith('```'):
        text = text.strip('`').lstrip('json').strip()
    try:
        obj = json.loads(text[text.index('{'):text.rindex('}') + 1])
        params = obj.get('params')
        if not isinstance(params, dict):
            return None
        norm, err = validate_params('rerun_stage', params)
        if err or not norm:
            return None
        return {'params': norm,
                'reason': str(obj.get('reason') or '')[:200]}
    except (ValueError, json.JSONDecodeError):
        return None


def _last_run_params(sample):
    """取上一轮运行的参数（project manifest）与待跑阶段清单。

    返回 (body_dict_for_analysis_kwargs, stages)。manifest 键与
    _analysis_kwargs 请求体键存在历史差异（memory_gb→memory、
    do_specificity→specificity），这里做映射。
    """
    import os
    from Virus_Platform_Core.config import PLATFORM_ROOT
    from Virus_Platform_Core.utils import resolve_sample_name
    from Virus_Platform_Core.pipeline import (load_project_manifest,
                                              pipeline_overview)
    try:
        real = resolve_sample_name(sample, DIRS['results'])
        sd = os.path.join(DIRS['results'], real)
        man = load_project_manifest(sd) or {}
    except Exception:
        man = {}
    p = dict(man.get('params') or {})
    body = {}
    keymap = {'memory_gb': 'memory', 'do_specificity': 'specificity'}
    for k, v in p.items():
        if k in ('db_host', 'db_virus', 'chunk_dir'):
            continue                      # 绝对路径类：留给 _analysis_kwargs 缺省
        body[keymap.get(k, k)] = v
    try:
        ov = pipeline_overview(sd)['stages']
        stages = [x['stage'] for x in ov
                  if x['stage'] != 'subsample'
                  and x.get('status') not in ('done', 'unavailable')]
    except Exception:
        stages = None
    return body, (stages or None)
