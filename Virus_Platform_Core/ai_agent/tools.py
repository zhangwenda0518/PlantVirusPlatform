# -*- coding: utf-8 -*-
"""工具注册表：把平台能力暴露成 LLM 可调用的函数（OpenAI tools 协议）。

安全边界（工具是代理唯一的"动手"通道，必须收紧）：
  - 全部工具**只读**：读样品列表 / 读 digest / 读任务状态 / 读结果表 / 读预筛产物；
  - 读表工具各有一个**单一根白名单**，且共用同一份路径校验实现
    （`_safe_read_table_under`）：`read_result_table` → `run/results/<样品>/`，
    `read_tool_run_file` → `run/tool_runs/<任务>/`，
    `read_prescreen_file` → `run/meta_search/<项目>/prescreen/`（**只到
    prescreen/**，不放开整个 meta_search —— 那里还有原始下载件与检索缓存）。
    根与目标都过 check_path（in_platform）+ 规范化前缀比对（防符号链接/`..`
    逃逸）+ 64MB 体积闸，行/列/单元格三重截断，防上下文炸弹与路径穿越；
  - 不提供任何"重跑阶段 / 改配置 / 删数据"的工具——分析动作仍由用户在
    界面上亲手发起（顾问层边界，见包 docstring）。
"""
import json
import os
import re

# OpenAI tools 协议的 schema 声明（provider 层原样透传）
TOOLS_SPEC = [
    {'type': 'function', 'function': {
        'name': 'list_samples',
        'description': '列出平台全部样品及各阶段完成度（done/total）',
        'parameters': {'type': 'object', 'properties': {}, 'required': []}}},
    {'type': 'function', 'function': {
        'name': 'get_sample_digest',
        'description': '读取指定样品的结果摘要：各阶段状态、检出病毒、组装/注释/'
                       '宿主预测/进化/引物的关键数字与表格头',
        'parameters': {'type': 'object', 'properties': {
            'sample': {'type': 'string',
                       'description': '样品名（results/ 下的目录名）'}},
            'required': ['sample']}}},
    {'type': 'function', 'function': {
        'name': 'get_task_overview',
        'description': '查看任务中心：运行中/最近完成的分析任务及进度',
        'parameters': {'type': 'object', 'properties': {
            'limit': {'type': 'integer', 'description': '最多返回条数，默认 12'}},
            'required': []}}},
    {'type': 'function', 'function': {
        'name': 'read_result_table',
        'description': '读取样品结果目录下某个 TSV/CSV 表的表头与前若干行'
                       '（只限 results/ 内的文本表格，不适用于序列/BAM）',
        'parameters': {'type': 'object', 'properties': {
            'sample': {'type': 'string', 'description': '样品名'},
            'relpath': {'type': 'string',
                        'description': '相对样品目录的路径，如 '
                                       '02b_kvsuite/identify/all_viruses.summary.tsv'},
            'max_rows': {'type': 'integer',
                         'description': '最多数据行数（默认 12，上限 30）'}},
            'required': ['sample', 'relpath']}}},
    {'type': 'function', 'function': {
        'name': 'read_tool_run_file',
        'description': '读取某次工具箱/扩展分析任务目录下的 TSV/CSV 文本表，'
                       '返回表头与前若干行。任务名从平台上下文的 tool_runs 里'
                       '取（如 align_20261004_042206_4e2f3e）。只限该任务目录'
                       '内的文本表格，不适用于序列/BAM。',
        'parameters': {'type': 'object', 'properties': {
            'run': {'type': 'string',
                    'description': '任务目录名（run/tool_runs/ 下的一级目录）'},
            'relpath': {'type': 'string',
                        'description': '相对该任务目录的路径'},
            'max_rows': {'type': 'integer',
                         'description': '最多数据行数（默认 12，上限 30）'}},
            'required': ['run', 'relpath']}}},
    {'type': 'function', 'function': {
        'name': 'read_prescreen_file',
        'description': '读取病毒样本预筛（/prescreen）项目产物目录内的表：'
                       'run/meta_search/<项目>/prescreen/。常用文件：'
                       'stat_run_summary.tsv（权威判定，一行一个 run，'
                       'pass 列 = 判定通过）、stat_viral_taxa.tsv（检出的病毒分类）、'
                       'prescreen_pass.tsv（通过的 run）、unresolved.tsv（未决）。'
                       '项目名从平台上下文的 prescreen 段取。返回表头与前若干行。',
        'parameters': {'type': 'object', 'properties': {
            'project': {'type': 'string',
                        'description': '项目目录名（run/meta_search/ 下的一级目录）'},
            'relpath': {'type': 'string',
                        'description': '相对该项目 prescreen/ 目录的路径，如 '
                                       'stat_run_summary.tsv'},
            'max_rows': {'type': 'integer',
                         'description': '最多数据行数（默认 12，上限 30）'}},
            'required': ['project', 'relpath']}}},
    {'type': 'function', 'function': {
        'name': 'propose_action',
        'description': '发起「AI 执行单」：你**不能直接执行**任何动作，'
                       '只能起草待批执行单（用户在对话坞批准后才会运行）。'
                       '动作：rerun_stage=按给定参数重跑样品管道阶段；'
                       'retry_failed=按上一轮原参数重跑未完成阶段；'
                       'download_runs=把公共数据编号（SRR/ERR/DRR/CRR）下载'
                       '到平台；run_tool=运行工具箱任一工具（比对/SDT/RDP/'
                       '引物设计/共识/ORF/系统地理等 ~30 个，params.tool 必填，'
                       '其余键=该工具界面表单的参数名，文件参数传绝对路径）；'
                       'create_sample=从已落盘的 FASTQ 创建样品。'
                       '必须在 reason 里写清为什么。',
        'parameters': {'type': 'object', 'properties': {
            'action': {'type': 'string',
                       'enum': ['rerun_stage', 'retry_failed',
                                'download_runs', 'run_tool',
                                'create_sample']},
            'sample': {'type': 'string',
                       'description': '样品名（download_runs 不需要）'},
            'stages': {'type': 'array', 'items': {'type': 'string'},
                       'description': '要重跑的阶段键（如 assembly/kvsuite/'
                                      'primer），rerun_stage 必填'},
            'params': {'type': 'object',
                       'description': '按动作填写：rerun_stage/retry_failed='
                                      'run_analysis 白名单键；download_runs='
                                      '{"runs": ["SRR…","CRR…"], "name": "…", '
                                      '"convert_sra": true}；run_tool='
                                      '{"tool": "align", …该工具表单参数，'
                                      '文件键推荐平台根相对路径，如 '
                                      'run/results/<样品>/03_assembly/'
                                      'viral_contigs.fasta}；create_sample='
                                      '{"sample": "名", "r1": "FASTQ 路径", '
                                      '"r2": "可选", "project": "可选"}',
                       },
            'reason': {'type': 'string', 'description': '提案理由（必填）'}},
            'required': ['action', 'reason']}}},
    {'type': 'function', 'function': {
        'name': 'list_platform_tools',
        'description': '列出工具箱全部可运行工具（名字/标题/是否需病毒库）——'
                       '这些都能通过 propose_action(run_tool) 起草执行单运行。'
                       '用户想做比对、SDT、RDP、引物设计、共识、ORF、系统地理、'
                       'dsRNA、miRNA 等操作时，先查这里拿 tool 名',
        'parameters': {'type': 'object', 'properties': {}, 'required': []}}},
    {'type': 'function', 'function': {
        'name': 'search_public_data',
        'description': '在公共测序库检索数据（NCBI SRA + NGDC GSA，双引擎）。'
                       '转录组/RNA-seq、基因组、病毒组等测序策略都能查——'
                       '用户想"找某物种的测序数据/转录组数据"就用本工具，'
                       '**不要**说平台查不了。species 传拉丁名'
                       '（如 Lycium barbarum）。检索是后台任务（数分钟），'
                       '启动后立即返回集合名；完成后平台自动播报结果，'
                       '期间可用 get_meta_results 查询进度。',
        'parameters': {'type': 'object', 'properties': {
            'species': {'type': 'string', 'description': '物种拉丁名'},
            'source': {'type': 'string',
                       'enum': ['TRANSCRIPTOMIC', 'GENOMIC', 'ALL'],
                       'description': '测序策略，默认 TRANSCRIPTOMIC（转录组）'},
            'db': {'type': 'string', 'enum': ['both', 'sra', 'gsa'],
                   'description': '检索库，默认 both（SRA+GSA）'}},
            'required': ['species']}}},
    {'type': 'function', 'function': {
        'name': 'list_meta_collections',
        'description': '列出平台上已有的公共检索集合（名字 + run 总数）。'
                       '平台上下文的 meta_collections 段通常已含此信息；'
                       '本工具用于按需刷新',
        'parameters': {'type': 'object', 'properties': {}, 'required': []}}},
    {'type': 'function', 'function': {
        'name': 'get_meta_results',
        'description': '读取一次公共检索的结果表（SRA+GSA 合并清单：'
                       'run 编号、物种、测序策略、标题等），返回表头与前若干行',
        'parameters': {'type': 'object', 'properties': {
            'collection': {'type': 'string',
                           'description': '检索集合名（search_public_data '
                                          '返回的 collection）'},
            'limit': {'type': 'integer',
                      'description': '最多数据行数（默认12，上限30）'}},
            'required': ['collection'],
        }}},
    {'type': 'function', 'function': {
        'name': 'mcp_tool_catalog',
        'description': '列出 MCP 执行层全部 vp_* 工具（名称/必填参数/键清单）：'
                       '工具箱执行 vp_run_tool、流水线 vp_pipeline_start、'
                       'CLI vp_cli_run·vp_cli_submit（建库/公共数据/NCBI提交/'
                       'LOGAN/ICTV/维护）、预筛 vp_prescreen_run、任务与运行'
                       '查询、vp_read_run_file。填参数前先调本工具查口径',
        'parameters': {'type': 'object', 'properties': {}, 'required': []}}},
    {'type': 'function', 'function': {
        'name': 'mcp_call',
        'description': '执行一个 MCP 层 vp_* 工具。任务类(vp_run_tool/'
                       'vp_pipeline_start/vp_cli_submit/vp_prescreen_run)提交即'
                       '返回 task_id，随后用 mcp_call 调 vp_get_task 轮询到 '
                       'done/failed，vp_read_run_file 读产物；信封里有 '
                       'summary/interpretation/next。注意：vp_delete_tool_run '
                       '被禁用',
        'parameters': {'type': 'object', 'properties': {
            'tool': {'type': 'string', 'description': 'vp_* 工具名'},
            'arguments': {'type': 'object',
                          'description': '该工具的参数字典'}},
            'required': ['tool']}}},
]


def dispatch(name, args):
    """执行一次工具调用，返回（给 LLM 的）字符串结果。"""
    from . import digest
    args = args or {}
    try:
        if name == 'list_samples':
            return json.dumps(digest.samples_brief(), ensure_ascii=False)
        if name == 'get_sample_digest':
            d = digest.sample_digest(str(args.get('sample') or ''))
            if d is None:
                return json.dumps({'error': '样品不存在'}, ensure_ascii=False)
            return json.dumps(d, ensure_ascii=False)
        if name == 'get_task_overview':
            try:
                limit = max(1, min(30, int(args.get('limit') or 12)))
            except (TypeError, ValueError):
                limit = 12
            return json.dumps(digest.task_overview(limit), ensure_ascii=False)
        if name == 'read_result_table':
            return _read_result_table(args)
        if name == 'read_tool_run_file':
            return _read_tool_run_file(args)
        if name == 'read_prescreen_file':
            return _read_prescreen_file(args)
        if name == 'propose_action':
            from . import executors
            rec, err = executors.create_action(
                str(args.get('action') or ''),
                str(args.get('sample') or ''),
                stages=args.get('stages') or None,
                params=args.get('params') or None,
                reason=str(args.get('reason') or ''),
                source='chat')
            if err:
                return json.dumps({'error': err}, ensure_ascii=False)
            return json.dumps(
                {'created': True, 'id': rec['id'], 'status': rec['status'],
                 'note': '执行单已创建，等待用户在对话坞批准；批准前不会运行'},
                ensure_ascii=False)
        if name == 'list_platform_tools':
            return _list_platform_tools()
        if name == 'list_meta_collections':
            from . import digest
            return json.dumps(digest.meta_collections_brief(),
                              ensure_ascii=False)
        if name == 'search_public_data':
            return _search_public_data(args)
        if name == 'get_meta_results':
            return _get_meta_results(args)
        if name == 'mcp_tool_catalog':
            return _mcp_catalog()
        if name == 'mcp_call':
            return _mcp_call(args)
        return json.dumps({'error': f'未知工具 {name}'}, ensure_ascii=False)
    except Exception as e:                      # 工具失败应作为结果回给 LLM
        return json.dumps({'error': f'{type(e).__name__}: {e}'},
                          ensure_ascii=False)


# ------------------------------------------------------------------
# MCP 桥（2026-10-16）：让对话坞 LLM 直接驱动 virusplatform_mcp 的
# 19 个 vp_* 工具（执行/查询/读结果全量）。只开两座桥而非注册 19 个
# 重叠工具：mcp_tool_catalog 负责「发现+参数口径」，mcp_call 负责「执行」。
# 破坏面收敛：vp_delete_tool_run 在对话坞禁用（删除走网页任务中心）；
# 任务类工具提交即返回 task_id，LLM 用 mcp_call(vp_get_task) 轮询。
# ------------------------------------------------------------------
_MCP_BLOCKED = {'vp_delete_tool_run'}
_MCP_TEXT_CAP = 4000          # 给 LLM 的结果上限（防撑爆上下文）


def _mcp_catalog():
    from Virus_Platform_Core.web.mcp_console import _request
    try:
        res = _request('tools/list', timeout=150)
    except Exception as e:
        return json.dumps({'error': f'MCP 服务器不可用: {e}'},
                          ensure_ascii=False)
    out = [{'name': t.get('name'),
            'description': (t.get('description') or '').split(chr(10))[0][:160],
            'required': (t.get('inputSchema') or {}).get('required') or [],
            'properties': sorted((t.get('inputSchema') or {})
                                 .get('properties') or {})}
           for t in (res.get('tools') or [])]
    return json.dumps({'count': len(out), 'tools': out}, ensure_ascii=False)


def _mcp_call(args):
    from Virus_Platform_Core.web.mcp_console import _request
    name = str(args.get('tool') or '').strip()
    if not name.startswith('vp_'):
        return json.dumps({'error': '只允许调用 vp_* MCP 工具'},
                          ensure_ascii=False)
    if name in _MCP_BLOCKED:
        return json.dumps(
            {'error': name + ' 属破坏面工具，对话坞不开放；'
                         '删除运行请用网页端任务中心'},
            ensure_ascii=False)
    a = args.get('arguments')
    if not isinstance(a, dict):
        a = {}
    try:
        res = _request('tools/call', {'name': name, 'arguments': a},
                       timeout=330)
    except Exception as e:
        return json.dumps({'error': f'MCP 调用失败: {e}'}, ensure_ascii=False)
    text = chr(10).join(
        c.get('text', '') for c in (res.get('content') or [])
        if isinstance(c, dict) and c.get('type') == 'text')
    data = res.get('structuredContent')
    if data is None:
        try:
            data = json.loads(text)
        except ValueError:
            data = None
    body = data if data is not None else text
    out = json.dumps({'isError': bool(res.get('isError')), 'result': body},
                     ensure_ascii=False)
    if len(out) > _MCP_TEXT_CAP:
        out = out[:_MCP_TEXT_CAP] + '…[结果过长已截断；可换更具体的查询，'               '或用 vp_read_run_file 按文件读]'
    return out


def _safe_read_table_under(root, relpath, max_rows):
    """在**给定白名单根**下安全读一张文本表，返回给 LLM 的 JSON 串。

    三重防线（两个入口共用同一份实现，避免"只修一处"）：
      1. `check_path(..., in_platform=True)` 两次 —— 根与目标都必须在平台内；
      2. 规范化后前缀比对 —— 防符号链接 / `..` 相对段逃逸出根；
      3. 体积闸 64MB —— 防上下文炸弹。
    """
    from Virus_Platform_Core.utils import check_path
    from . import digest
    relpath = str(relpath or '').replace('\\', '/').lstrip('/')
    if not relpath:
        return json.dumps({'error': '需要 relpath'}, ensure_ascii=False)
    try:
        root = check_path(root, must_exist=True, in_platform=True)
        target = check_path(os.path.join(root, relpath),
                            must_exist=True, in_platform=True)
    except (ValueError, OSError) as e:
        return json.dumps({'error': f'路径不合法: {e}'}, ensure_ascii=False)
    if os.path.normcase(os.path.abspath(target)).find(
            os.path.normcase(os.path.abspath(root))) != 0:
        return json.dumps({'error': '路径越界'}, ensure_ascii=False)
    if os.path.getsize(target) > 64 * 1024 * 1024:
        return json.dumps({'error': '文件过大，不予读取'}, ensure_ascii=False)
    try:
        max_rows = max(1, min(30, int(max_rows or 12)))
    except (TypeError, ValueError):
        max_rows = 12
    t = digest._tsv_head(target, max_rows=max_rows)
    if t is None:
        return json.dumps({'error': '文件不存在或为空'}, ensure_ascii=False)
    return json.dumps(t, ensure_ascii=False)


def _read_result_table(args):
    """读样品结果目录（run/results/<sample>/）内的表。根 = results/<sample>。"""
    from Virus_Platform_Core.config import DIRS
    sample = str(args.get('sample') or '')
    if not sample:
        return json.dumps({'error': '需要 sample 与 relpath'},
                          ensure_ascii=False)
    return _safe_read_table_under(os.path.join(DIRS['results'], sample),
                                  args.get('relpath'), args.get('max_rows'))


def _read_tool_run_file(args):
    """读某次工具箱任务目录（run/tool_runs/<run>/）内的表。

    2026-10-08 新增：扩 digest 覆盖面时，主链之外的扩展分析（比对 / 组装 /
    ORF / dsRNA / miRNA / CRISPR…）产物都落在 tool_runs 下，没有这条通道，
    AI 只能知道"跑过"，答不了"结果是什么"。根白名单仍收紧在 tool_runs/。
    """
    from Virus_Platform_Core.config import DIRS
    run = str(args.get('run') or '')
    if not run:
        return json.dumps({'error': '需要 run（任务目录名）'},
                          ensure_ascii=False)
    return _safe_read_table_under(os.path.join(DIRS['tool_runs'], run),
                                  args.get('relpath'), args.get('max_rows'))


# 项目目录名（与 web/prescreen.py 的 _NAME_RE 同口径）：挡掉分隔符与 `..`
_PRESCREEN_NAME_RE = re.compile(r'[A-Za-z0-9_\-]+')


def _read_prescreen_file(args):
    """读病毒预筛项目目录（run/meta_search/<项目>/prescreen/）内的表。

    2026-10-08 新增：预筛产物**不在** run/tool_runs/ 下（`web/prescreen.py`
    只写 `DIRS['meta_search']`），故 `read_tool_run_file` 的白名单覆盖不到 ——
    此前 `prescreen_brief()` 只能给计数（n_runs/n_pass/n_ok），答不了
    「哪些 run 通过、检出哪些病毒」。根白名单**收紧到 prescreen/**：
    不放开整个 meta_search（那里还有原始下载件 `stat_raw/`、`gsa_raw/`
    与检索缓存 `.cache/`）。
    """
    from Virus_Platform_Core.config import DIRS
    project = str(args.get('project') or '').strip()
    if not project or not _PRESCREEN_NAME_RE.fullmatch(project):
        return json.dumps(
            {'error': '需要 project（项目目录名，只允许字母/数字/_/-）'},
            ensure_ascii=False)
    root = os.path.join(DIRS['meta_search'], project, 'prescreen')
    return _safe_read_table_under(root, args.get('relpath'),
                                  args.get('max_rows'))


# ------------------------------------------------------------------
# 公共数据检索（SRA + GSA）——"AI 驱动找数据"的核心工具
# 复用与「公共数据检索」页面完全相同的引擎与任务通道
# （web/meta.py 的 _meta_stream_task + public_meta/search_engine.py），
# 结果落 run/meta_search/<集合>/search/SRA_GSA_Merged_Final.csv。
# ------------------------------------------------------------------
def _meta_collection_name(species):
    """物种名 → 合法集合名（与 web/meta.py 同口径：<slug>_<时间戳>）。"""
    import time as _time
    slug = re.sub(r'[^A-Za-z0-9_\-]', '_', str(species))[:40].strip('_') \
        or 'search'
    ts = _time.strftime('%Y%m%d_%H%M')
    from Virus_Platform_Core.config import DIRS
    name = f'{slug}_{ts}'
    n = 2
    while os.path.isdir(os.path.join(DIRS['meta_search'], name)):
        name = f'{slug}_{ts}_{n}'
        n += 1
    return name


def _tool_param_hints(job):
    """从 job 工厂源码提取参数提示：_req('key','说明') = 必填，
    _opt('key') = 可选。免手写 30 份参数文档，且与界面表单同源。"""
    import inspect as _inspect
    import re as _re
    try:
        src = _inspect.getsource(job)
    except (TypeError, OSError):
        return {}
    # job 工厂取参有两种形态：ctx.req('key','what')（tools_api 闭包）与
    # _job_req(ctx,'key','what')（tool_jobs 助手）；可选取参同理。
    req = _re.findall(r"ctx\.req\('(\w+)',\s*'([^']+)'", src)
    req += _re.findall(r"_job_req\(ctx,\s*'(\w+)',\s*'([^']+)'", src)
    opt = _re.findall(r"ctx\.opt\('(\w+)'", src)
    opt += _re.findall(r"_job_opt\(ctx,\s*'(\w+)'", src)
    seen, req2 = set(), []
    for k, w in req:
        if k not in seen:
            seen.add(k)
            req2.append((k, w))
    req = req2[:8]
    opt = list(dict.fromkeys(opt))[:10]
    out = {}
    if req:
        out['required'] = [{'key': k, 'what': w} for k, w in req]
    if opt:
        out['optional'] = opt
    return out


def _list_platform_tools():
    try:
        from Virus_Platform_Core.web.tools_api import TOOL_REGISTRY,             LIGHT_TOOLS
        out = []
        for k, v in sorted(TOOL_REGISTRY.items()):
            entry = {'tool': k, 'title': v['title'],
                     'need_db': bool(v.get('need_db')),
                     'light': k in LIGHT_TOOLS}
            hints = _tool_param_hints(v.get('job'))
            if hints:
                entry['params_hint'] = hints
            out.append(entry)
        return json.dumps({'count': len(out), 'tools': out,
                           'note': '用 propose_action(action=run_tool, '
                                   'params={tool:…, …params_hint 里的参数键}) '
                                   '起草执行单；文件键优先用平台根相对路径'},
                          ensure_ascii=False)
    except Exception as e:
        return json.dumps({'error': f'{type(e).__name__}: {e}'},
                          ensure_ascii=False)


def _search_public_data(args):
    """启动一次公共检索（后台任务），返回任务号与集合名。

    完成播报由 _watch_search 守护线程负责：检索完 → ai_read 事件自动到账
    （带 Top 结果摘要），用户不用盯着任务中心。
    """
    import threading
    import time as _time
    from Virus_Platform_Core.config import DIRS, engine_cmd
    from Virus_Platform_Core.utils import check_path
    species = str(args.get('species') or '').strip()
    if not species:
        return json.dumps({'error': '需要 species（物种拉丁名）'},
                          ensure_ascii=False)
    source = str(args.get('source') or 'TRANSCRIPTOMIC')
    if source not in ('TRANSCRIPTOMIC', 'GENOMIC', 'ALL'):
        source = 'TRANSCRIPTOMIC'
    db = str(args.get('db') or 'both')
    if db not in ('both', 'sra', 'gsa'):
        db = 'both'
    name = _meta_collection_name(species)
    out_dir = check_path(os.path.join(DIRS['meta_search'], name, 'search'),
                         must_exist=False, in_platform=True)
    os.makedirs(out_dir, exist_ok=True)
    from Virus_Platform_Core.public_meta import search_engine
    cmd = engine_cmd(search_engine.__file__,
                     '-q', species, '-s', source, '--db', db,
                     '--workers', '5', '-o', out_dir)
    try:
        from Virus_Platform_Core.web.meta import _meta_stream_task
        tid = _meta_stream_task(
            f'公共检索 {species}', cmd, f'检索完成: {out_dir}',
            log_file=os.path.join(out_dir, 'run.log'), link='/meta')
    except Exception as e:
        return json.dumps({'error': f'检索任务启动失败: {e}'},
                          ensure_ascii=False)
    threading.Thread(target=_watch_search,
                     args=(tid, name, species), daemon=True,
                     name=f'ai-meta-watch-{name}').start()
    # 已有同物种历史集合：秒级可读，主动告诉 AI（不必干等网络检索）
    existing = _existing_collections(species)
    return json.dumps(
        {'started': True, 'task': tid, 'collection': name,
         'source': source, 'db': db,
         'existing_collections': existing,
         'note': '检索任务已启动（数分钟）；完成后平台自动播报 Top 结果，'
                 '也可随时用 get_meta_results 查询'
                 + ('。注意：平台已有同物种历史检索集合 '
                    + '、'.join(existing)
                    + '，可先用 get_meta_results 直接读取，无需等本次检索'
                    if existing else '')},
        ensure_ascii=False)


def _existing_collections(species, cap=5):
    """同物种的历史检索集合（名字以物种 slug 开头且结果表有数据）。"""
    from Virus_Platform_Core.config import DIRS
    slug = re.sub(r'[^A-Za-z0-9_\-]', '_', str(species))[:40].strip('_')
    if not slug:
        return []
    root = DIRS['meta_search']
    out = []
    try:
        for nm in sorted(os.listdir(root), reverse=True):
            if not nm.startswith(slug):
                continue
            if os.path.isfile(os.path.join(root, nm, 'search',
                                           'SRA_GSA_Merged_Final.csv')):
                out.append(nm)
            if len(out) >= cap:
                break
    except OSError:
        pass
    return out


def _watch_search(tid, collection, species):
    """检索任务完成 → 自动播报结果摘要（ai_read 事件）。最长守 30 分钟。"""
    import time
    from . import broadcast
    try:
        from Virus_Platform_Core.web.tasks import tm
    except Exception:
        return
    deadline = time.time() + 1800
    while time.time() < deadline:
        try:
            snap = tm.snapshot(tid, log_lines=0, include_result=False)
        except Exception:
            return
        if snap is None:
            return
        st = snap.get('status')
        if st != 'running':
            if st == 'done':
                res = _read_meta_csv(collection, 5)
                if res and res.get('rows'):
                    head = res['header'] or []
                    lines = []
                    for row in res['rows'][:5]:
                        cell = ' | '.join(
                            f'{head[i]}: {v}' for i, v in enumerate(row)
                            if v) [:200]
                        lines.append(f'- {cell}')
                    text = (f'公共检索完成：**{species}**（集合 `{collection}`）'
                            f'共 **{res.get("n_total", "?")}** 条 run。'
                            f'Top 预览：\n' + '\n'.join(lines))
                else:
                    text = (f'公共检索完成：**{species}**（集合 `{collection}`），'
                            f'结果表暂无法读取，请到「公共数据检索」页查看。')
                broadcast.publish('ai_read', sample='', stage='meta_search',
                                  name=f'公共检索 · {species}', text=text)
            else:
                broadcast.publish('stage_failed', sample='',
                                  stage='meta_search',
                                  name=f'公共检索 · {species}',
                                  error=f'检索任务状态: {st}')
            return
        time.sleep(3.0)


def _read_meta_csv(collection, max_rows=12):
    """读检索结果 CSV（逗号分隔，与 digest._tsv_head 的 TSV 口径不同）。

    返回 {file, header, rows, n_total, truncated}；无结果返回 None。
    """
    import csv as _csv
    from Virus_Platform_Core.config import DIRS
    from Virus_Platform_Core.utils import check_path
    collection = str(collection or '').strip()
    if not collection or not _PRESCREEN_NAME_RE.fullmatch(collection):
        return None
    try:
        root = check_path(DIRS['meta_search'], must_exist=True,
                          in_platform=True)
        path = check_path(os.path.join(root, collection, 'search',
                                       'SRA_GSA_Merged_Final.csv'),
                          must_exist=True, in_platform=True)
    except (ValueError, OSError):
        return None
    if os.path.getsize(path) > 64 * 1024 * 1024:
        return None
    try:
        max_rows = max(1, min(30, int(max_rows or 12)))
    except (TypeError, ValueError):
        max_rows = 12
    header, rows, n_total = None, [], 0
    # utf-8-sig：Excel 导出的 CSV 常带 BOM，否则表头首列会出现 ﻿
    with open(path, encoding='utf-8-sig', errors='replace', newline='') as f:
        reader = _csv.reader(f)
        for i, row in enumerate(reader):
            if i == 0:
                header = [c[:80] for c in row][:10]
                continue
            if row and any(c.strip() for c in row):
                n_total += 1
                if len(rows) < max_rows:
                    rows.append([c[:80] for c in row][:10])
    if header is None:
        return None
    return {'file': 'SRA_GSA_Merged_Final.csv', 'header': header,
            'rows': rows, 'n_total': n_total,
            'truncated': n_total > len(rows)}


def _get_meta_results(args):
    res = _read_meta_csv(str(args.get('collection') or ''),
                         args.get('max_rows') or args.get('limit'))
    if res is None:
        return json.dumps(
            {'error': '结果表不存在或尚无可读结果（检索可能仍在进行）'},
            ensure_ascii=False)
    return json.dumps(res, ensure_ascii=False)
