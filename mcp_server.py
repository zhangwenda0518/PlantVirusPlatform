# -*- coding: utf-8 -*-
"""virusplatform_mcp —— 植物病毒分析平台 MCP 服务器（stdio）。

把平台三大能力原样暴露为 MCP 工具，供 Claude Desktop / ZCode 等
客户端直接驱动：
  1. 工具箱全量可执行（web/tools_api.py 的 TOOL_REGISTRY）
  2. 样品流水线全流程驱动（pipeline.run_analysis，⓪~⑩ 15 阶段任选）
  3. 任务/运行/产物查询（TaskManager 快照、run 目录文件读取）

实现口径与网页端完全同源：工具运行走同一套 job 工厂 + TaskManager
（heavy/light 并发闸门、run.log、结果面板可视化补全都保留），流水线
走同一套 run_analysis 编排。长任务一律「提交返回 task_id +
vp_get_task 轮询」两个动作，绝不同步阻塞。

启动（与平台同一解释器，缺 mcp 包自动安装）：
    3rd\\python\\python.exe mcp_server.py              # stdio 服务（默认）
    3rd\\python\\python.exe mcp_server.py --http       # streamable-HTTP 模式
    3rd\\python\\python.exe mcp_server.py --http 9000  # HTTP 自定义端口
    3rd\\python\\python.exe mcp_server.py --http --token s3cret
                                                       # HTTP + Bearer 鉴权
    3rd\\python\\python.exe mcp_server.py --check      # 自检（不进服务）

下线指定工具（不删代码）：环境变量 VP_MCP_DISABLED_TOOLS=vp_a,vp_b 或
平台根 deferred.json {"disabled_tools": ["vp_a"]}；文件存在但解析失败
会拒绝启动（fail-closed）。

stdio：任务随客户端进程存活；HTTP（127.0.0.1:8765/mcp）：任务与客户端
解耦，断开照跑、重连可查，局域网内多客户端可共享。两种模式任务都在
本进程执行；与网页端同时运行时并发闸门各自计数。
"""
import json
import os
import re
import sys
import time
import uuid
from typing import Annotated, Optional

# 嵌入式 Python（3rd/python，_pth 布局）不自动注入脚本目录——与 app.py /
# main.py 同口径手动补平台根，Virus_Platform_Core 才可导入
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# ---- stdio 纪律：协议通道独占 stdout，一切输出走 stderr ----
for _s in (sys.stdout, sys.stderr):
    try:
        _s.reconfigure(encoding='utf-8', errors='replace')
    except (AttributeError, ValueError):
        pass


def _ensure_mcp_sdk():
    """缺 mcp 包自动安装一次（与启动平台-网页.bat 的依赖兜底同口径）。

    全程只写 stderr：stdio 传输下 stdout 是协议通道，pip 的任何输出
    混进去都会让客户端握手失败。
    """
    try:
        import mcp  # noqa: F401
        return
    except ImportError:
        pass
    import subprocess
    sys.stderr.write('[virusplatform_mcp] 首次运行：安装 mcp SDK v1（约 1 分钟）...\n')
    for mirror in ('https://pypi.tuna.tsinghua.edu.cn/simple', None):
        # 锁 mcp<2：2.x 把 FastMCP 改名为 MCPServer 且 API 不兼容（本服务器
        # 按 v1 官方 SDK 编写；已有 v2 环境会被 pip 自动降级）
        cmd = [sys.executable, '-m', 'pip', 'install', '--quiet', 'mcp<2']
        if mirror:
            cmd += ['-i', mirror]
        try:
            if subprocess.run(cmd).returncode == 0:
                break
        except OSError:
            continue
    import mcp  # noqa: F401  再失败就让 ImportError 原样抛给客户端


_ensure_mcp_sdk()

import anyio                                                     # noqa: E402
from mcp.server.fastmcp import FastMCP                           # noqa: E402
from mcp.types import ToolAnnotations                            # noqa: E402
from pydantic import Field                                       # noqa: E402

# ---- 平台模块（脚本位于平台根，sys.path[0] 即根）----
from Virus_Platform_Core.config import DIRS, PLATFORM_ROOT       # noqa: E402
from Virus_Platform_Core.utils import (TaskLogger, check_path,   # noqa: E402
                                       safe_open)
from Virus_Platform_Core.web.state import (cfg, run_real_dir,    # noqa: E402
                                           tool_runs_root, RUN_OUT_LINK)
from Virus_Platform_Core.web.tasks import tm                     # noqa: E402
from Virus_Platform_Core.web.tools_api import (LIGHT_TOOLS,      # noqa: E402
                                               TOOL_REGISTRY,
                                               VIRUS_DB_PRESETS,
                                               _resolve_virus_db)
from Virus_Platform_Core.pipeline import (DEFAULT_ANALYZE_STAGES,  # noqa: E402
                                          STAGE_GROUPS,
                                          STAGE_NAMES,
                                          STAGE_NAMES_EN,
                                          STAGE_ORDER,
                                          run_analysis)

# 工厂期缺参/选库失败走 flask.abort——统一转成带原文案的 ValueError
from werkzeug.exceptions import HTTPException                    # noqa: E402

_PLATFORM_OK = None

SERVER_VERSION = '1.2.0'


def _gated_tools():
    """启动下线清单（借鉴 wisp gate 的 fail-closed 思路）。

    两个来源并集：环境变量 VP_MCP_DISABLED_TOOLS（逗号分隔）与平台根
    deferred.json（{"disabled_tools": [...]}）。deferred.json 存在但解析
    失败 → 直接抛错拒启：写了闸门文件的运维预期是它必须生效，静默忽略
    会把该摘除的工具放出去。
    """
    names = set()
    env = os.environ.get('VP_MCP_DISABLED_TOOLS', '')
    names.update(t.strip() for t in env.replace(';', ',').split(',') if t.strip())
    dj = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                      'deferred.json')
    if os.path.isfile(dj):
        try:
            with open(dj, encoding='utf-8') as f:
                data = json.load(f)
        except (OSError, ValueError) as e:
            raise SystemExit(f'[virusplatform_mcp] deferred.json 解析失败，'
                             f'拒绝启动（fail-closed）: {e}')
        items = data.get('disabled_tools') if isinstance(data, dict) else data
        if isinstance(items, list):
            names.update(str(t).strip() for t in items if str(t).strip())
    return names


def _capability_fingerprint(enabled_tools):
    """能力指纹：版本 + 启用工具清单的 sha256 前 12 位。

    agent 会话内可据此感知「服务端能力变了」（升级/摘除工具后指纹必变）。
    """
    import hashlib
    payload = SERVER_VERSION + '|' + ','.join(sorted(enabled_tools))
    return hashlib.sha256(payload.encode('utf-8')).hexdigest()[:12]


def _platform_check():
    """平台自检一次并缓存结论（--check 与提交类工具共用）。"""
    global _PLATFORM_OK
    if _PLATFORM_OK is not None:
        return _PLATFORM_OK
    problems = []
    if not TOOL_REGISTRY:
        problems.append('TOOL_REGISTRY 为空')
    _PLATFORM_OK = problems
    return problems


mcp = FastMCP(
    'virusplatform_mcp',
    instructions=(
        '植物病毒分析平台（Windows 本地）的 MCP 入口。'
        '典型流程：vp_platform_info 看环境 → vp_list_samples / vp_list_tools '
        '选目标 → vp_run_tool（工具箱）或 vp_pipeline_start（全流程）提交 → '
        'vp_get_task 轮询到终态 → vp_read_run_file 读产物文本。'
        '数据库构建、公共数据检索（NCBI/GenBank/宿主基因组）、NCBI 提交、'
        'LOGAN 溯源、ICTV 库等 CLI 能力：vp_cli_help 查清单与参数，快查用 '
        'vp_cli_run，长任务用 vp_cli_submit。'
        '所有执行类工具提交即返回 task_id，长任务（组装/建树/定量/下载）'
        '必须轮询，不要重复提交。'
    ),
)

_RO = ToolAnnotations(readOnlyHint=True, idempotentHint=True)
_RO_LIVE = ToolAnnotations(readOnlyHint=True)
_DESTR = ToolAnnotations(destructiveHint=True)


# ======================================================================
# 内部实现（同步函数；工具层经 anyio.to_thread 调用，不卡事件循环）
# ======================================================================

def _input_abs(v):
    """输入数据文件：相对路径按平台根解析；绝对路径允许电脑任意位置（只读）。"""
    return check_path(v if os.path.isabs(v)
                      else os.path.join(PLATFORM_ROOT, v), must_exist=True)


def _clamp_threads(threads):
    """数字参数钳制，与 /api/tool/run 同口径。"""
    if threads is None:
        return None
    try:
        return max(1, min(int(threads), (os.cpu_count() or 4) * 4))
    except (TypeError, ValueError):
        return None


def _submit_tool_run(tool, params, threads, out_dir):
    """复刻 /api/tool/run 的完整提交链（tools_api.api_tool_run）。"""
    problems = _platform_check()
    if problems:
        raise ValueError('平台模块未就绪（' + '; '.join(problems) +
                         '）；先在平台根跑 python mcp_server.py --check 排查')
    spec = TOOL_REGISTRY.get(tool)
    if not spec:
        raise ValueError('未知工具: %s（可用: %s）'
                         % (tool, ', '.join(sorted(TOOL_REGISTRY))))
    p = dict(params or {})
    threads = _clamp_threads(threads)

    def _req(key, what):
        v = (p.get(key) or '').strip()
        if not v:
            raise ValueError('缺少必填参数 %s（%s）' % (key, what))
        return _input_abs(v)

    def _opt(key):
        v = (p.get(key) or '').strip()
        return _input_abs(v) if v else None

    try:
        db_virus = _resolve_virus_db(p.get('db_virus'))
    except HTTPException as e:
        raise ValueError(e.description or str(e)) from e
    if spec['need_db']:
        from Virus_Platform_Core.kunpeng import db_ready
        if not db_ready(db_virus):
            raise ValueError('病毒库未就绪: %s（先在平台「数据库构建」页构建，'
                             '或 db_virus 换选其他内置库: %s）'
                             % (db_virus, ', '.join(VIRUS_DB_PRESETS)))

    ts = time.strftime('%Y%m%d_%H%M%S')
    run_name = f'{tool}_{ts}_{uuid.uuid4().hex[:6]}'
    base_dir = check_path(os.path.join(tool_runs_root(), run_name),
                          must_exist=False, in_platform=True)
    run_dir = base_dir
    if out_dir:
        out_root = check_path(out_dir, must_exist=False, in_platform=True)
        _out_key = os.path.normcase(os.path.normpath(os.path.abspath(out_root)))
        _bad_roots = (
            ('样品结果目录', DIRS['results']),
            ('数据库目录', DIRS['databases']),
            ('任务状态目录', DIRS['tasks']),
            ('前端资源目录', os.path.join(PLATFORM_ROOT, 'webapp')),
            ('第三方依赖目录', os.path.join(PLATFORM_ROOT, '3rd')),
        )
        for _label, _bad in _bad_roots:
            _bad_key = os.path.normcase(
                os.path.normpath(os.path.abspath(str(_bad))))
            if _out_key == _bad_key or _out_key.startswith(_bad_key + os.sep):
                raise ValueError('输出目录不能位于%s内（%s）；请换平台 run/ '
                                 '下或其他目录' % (_label, _bad))
        run_dir = check_path(os.path.join(out_root, run_name),
                             must_exist=False, in_platform=True)
    os.makedirs(run_dir, exist_ok=True)
    if run_dir != base_dir:
        os.makedirs(base_dir, exist_ok=True)
        try:
            with safe_open(os.path.join(base_dir, RUN_OUT_LINK), 'wt') as f:
                json.dump({'out_dir': run_dir, 'tool': tool}, f,
                          ensure_ascii=False, indent=2)
        except OSError:
            pass

    from types import SimpleNamespace
    ctx = SimpleNamespace(p=p, run_dir=run_dir, threads=threads,
                          db_virus=db_virus, req=_req, opt=_opt)
    try:
        job = spec['job'](ctx)
    except HTTPException as e:
        # 工厂期报错不留下空运行目录（与 api_tool_run 同口径）
        import shutil as _sh
        _sh.rmtree(run_dir, ignore_errors=True)
        if run_dir != base_dir:
            _sh.rmtree(base_dir, ignore_errors=True)
        raise ValueError(e.description or str(e)) from e
    except Exception:
        import shutil as _sh
        _sh.rmtree(run_dir, ignore_errors=True)
        if run_dir != base_dir:
            _sh.rmtree(base_dir, ignore_errors=True)
        raise

    name = cfg.tr(f"工具·{spec['title']} {ts}",
                  f"Tool·{spec['title_en']} {ts}")

    def _job_with_run(log, prog, cancel):
        """与 tools_api._job_with_run 同口径：结果注入 run/out_dir、
        表格产物自动补图、外移运行 touch 指针目录。"""
        try:
            res = job(log, prog, cancel)
        except HTTPException as e:
            raise RuntimeError(str(e.description or e)) from e
        if isinstance(res, dict):
            res.setdefault('run', run_name)
            res['out_dir'] = str(run_dir)
        try:
            from Virus_Platform_Core.toolbox_plots import add_plots
            for _fp in add_plots(str(run_dir)):
                log(f"图 -> {_fp}")
        except Exception:
            pass
        try:
            os.utime(base_dir)
        except OSError:
            pass
        return res

    weight = 'light' if tool in LIGHT_TOOLS else 'heavy'
    tid = tm.start(name, _job_with_run,
                   log_file=os.path.join(run_dir, 'run.log'),
                   weight=weight, out_dir=str(run_dir))
    return _fmt(
        'submitted', f"工具 {tool} 已提交（{weight} 档），任务 {tid}",
        result={'task_id': tid, 'run': run_name, 'out_dir': str(run_dir),
                'tool': tool, 'weight': weight},
        nxt=[f"vp_get_task(task_id='{tid}') 轮询到 done/failed",
             f"vp_read_run_file(run='{run_name}', path=…) 读产物文本"],
        hint='用 vp_get_task(task_id=...) 轮询；产物文本用 '
             'vp_read_run_file(run=%s, path=...) 读取' % run_name)


def _pipeline_job(sample, r1, r2, stages, kw):
    """流水线任务体：run_analysis 全参数透传，日志/进度接 TaskManager。"""

    def job(log, prog, cancel):
        logger = TaskLogger(callback=log)
        try:
            return run_analysis(
                sample=sample, r1=r1, r2=(r2 or None),
                stages=stages, logger=logger,
                progress=lambda stage, pct, msg, *eta: prog(
                    str(stage), float(pct), str(msg)),
                **kw)
        finally:
            logger.close()
    return job


def _submit_pipeline(args):
    from Virus_Platform_Core.pipeline import _safe_sample_name
    sample = args['sample']
    r1 = _input_abs(args['r1'])
    r2 = _input_abs(args['r2']) if args.get('r2') else None
    stages = args.get('stages') or list(DEFAULT_ANALYZE_STAGES)
    _bad = [s for s in stages if s not in STAGE_ORDER]
    if _bad:
        raise ValueError('未知阶段: %s（可用: %s；默认全流程: %s）'
                         % (', '.join(_bad), ', '.join(STAGE_ORDER),
                            ', '.join(DEFAULT_ANALYZE_STAGES)))
    vm = str(args.get('verify_methods') or 'blastx,bfvd')
    kw = {'db_host': args.get('db_host'), 'db_virus': args.get('db_virus'),
          'kv_lib': args.get('kv_lib'),
          'threads': _clamp_threads(args.get('threads')),
          'confidence': float(args.get('confidence') or 0.0),
          'assembly_mode': args.get('assembly_mode') or 'rnaviral',
          'memory_gb': int(args.get('memory_gb') or 64),
          'min_contig_len': int(args.get('min_contig_len') or 200),
          # chunk_dir 为内部分块目录，不对 MCP 暴露
          'assembly_input': args.get('assembly_input') or 'virus',
          'min_viral_pairs': int(args.get('min_viral_pairs') or 500),
          'do_trim': args.get('do_trim', True) is not False,
          'subsample': int(args.get('subsample') or 0),
          'subsample_mode': args.get('subsample_mode') or 'prefix',
          'subsample_seed': int(args.get('subsample_seed') or 42),
          'subsample_prop': float(args.get('subsample_prop') or 0.0),
          'tree_tool': args.get('tree_tool') or 'fasttree',
          'tree_sampling': args.get('tree_sampling') or 'blast',
          'top_n_refs': int(args.get('top_n_refs') or 10),
          'min_pident': args.get('min_pident'),
          'min_qcov': args.get('min_qcov'),
          'primer_mode': args.get('primer_mode') or 'conserved',
          'min_orf_aa': int(args.get('min_orf_aa') or 100),
          'fastp_dedup': bool(args.get('fastp_dedup')),
          'do_fq2fa': args.get('do_fq2fa', True) is not False,
          'gbdraw_max': int(args.get('gbdraw_max') or 12),
          'gbdraw_fasta': args.get('gbdraw_fasta'),
          'gbdraw_ann': args.get('gbdraw_ann'),
          'plot_engine': args.get('plot_engine') or 'auto',
          'verify_host': args.get('verify_host') or 'all',
          'verify_methods': tuple(
              m.strip() for m in vm.split(',') if m.strip()),
          'verify_combine': args.get('verify_combine') or 'union',
          'force': bool(args.get('force')),
          'do_specificity': bool(args.get('specificity'))}
    if args.get('ncbi_refs'):
        kw['ncbi_refs'] = args['ncbi_refs']

    safe = _safe_sample_name(sample)
    sd = check_path(os.path.join(DIRS['results'], safe),
                    must_exist=False, in_platform=True)
    os.makedirs(sd, exist_ok=True)
    ts = time.strftime('%Y%m%d_%H%M%S')
    name = cfg.tr(f'流水线·{safe} {ts}', f'Pipeline·{safe} {ts}')
    tid = tm.start(name, _pipeline_job(sample, r1, r2, stages, kw),
                   log_file=os.path.join(sd, 'pipeline_mcp.log'),
                   weight='heavy', out_dir=str(sd))
    return _fmt(
        'submitted', f"样品 {safe} 流水线已提交（{len(stages)} 个阶段）",
        result={'task_id': tid, 'sample': safe, 'sample_dir': str(sd),
                'stages': stages},
        nxt=[f"vp_get_task(task_id='{tid}') 轮询各阶段进度",
             f"完成后样品结果在 {sd}（网页端结果中心同一份）"],
        hint='轮询 vp_get_task；阶段产物在 sample_dir 下'
             '（00_prep/01_host_removal/03_assembly/...）')


def _task_snapshot(tid, log_lines, result_max_kb=256):
    """任务快照的稳定投影（借鉴 wisp records 层：固定键序、不透原始 dict）。

    键集合固定，便于 agent 结构化消费；result 超 result_max_kb 时截断为
    {_truncated, _note, _head}，完整结果指回 run/tasks/<tid>.json。
    """
    snap = tm.snapshot(tid, log_lines=max(0, min(int(log_lines), 500)))
    if not snap:
        raise ValueError('任务不存在或已超出内存保留窗口: %s'
                         '（历史任务可查 run/tasks/%s.json）' % (tid, tid))
    snap.pop('restart', None)            # 网页端专用，MCP 无意义
    if snap.get('started'):
        snap['started_h'] = time.strftime(
            '%Y-%m-%d %H:%M:%S', time.localtime(snap['started']))
    if snap.get('finished'):
        snap['finished_h'] = time.strftime(
            '%Y-%m-%d %H:%M:%S', time.localtime(snap['finished']))
    kb = max(1, min(int(result_max_kb), 4096))
    res = snap.get('result')
    if res is not None:
        text = json.dumps(res, ensure_ascii=False)
        if len(text.encode('utf-8', 'replace')) > kb * 1024:
            snap['result'] = {'_truncated': True,
                              '_note': f'终态结果超 {kb} KB；完整结果在 '
                                       f'run/tasks/{tid}.json 或产物目录 '
                                       f'{snap.get("out_dir") or ""}',
                              '_head': text[:kb * 1024]}
    keys = ('id', 'name', 'status', 'stage', 'pct', 'msg', 'error', 'eta',
            'out_dir', 'link', 'started_h', 'finished_h', 'log', 'result')
    out = {k: snap.get(k) for k in keys}
    # ---- 统一信封：摘要 / 解读 / 下一步（规则生成）----
    st, stage = out['status'], out['stage'] or '—'
    pct = int(round((out['pct'] or 0) * 100))
    out['summary'] = f"{st} · {stage} {pct}% · {out['msg'] or ''}".rstrip()
    interp, nxt = [], []
    if st == 'running':
        interp.append(f"当前阶段「{stage}」进行到 {pct}%。")
        nxt.append('稍后再次调用 vp_get_task 轮询，直到 status 变为 '
                   'done/failed/cancelled')
    elif st == 'done':
        interp.append('任务成功完成；result 字段为各步摘要，产物在 out_dir。')
        nxt.append('用 vp_read_run_file(run, path) 读取关键文本产物')
        nxt.append('在网页端「任务中心/结果中心」可看同一份结果')
    elif st == 'failed':
        err = out['error'] or ''
        interp.append(f"失败原因：{err[:300]}")
        hits = _interpret_error(err)
        interp.extend(hits or ['按错误信息修正后重新提交即可'])
        nxt.append('修正参数/环境后重新提交（vp_run_tool / vp_pipeline_start /'
                   ' vp_cli_submit）')
    elif st == 'cancelled':
        interp.append('任务已被取消；产物目录保留已生成的部分。')
    out['interpretation'] = interp
    out['next'] = nxt
    return out


def _task_list(limit, include_archived=False):
    limit = max(1, min(int(limit), 100))
    outs = [x for x in tm.list_all(limit=limit, logs=False) if x]
    for o in outs:
        o.pop('restart', None)
        o.pop('log', None)
    archived, archived_count = [], 0
    if include_archived:
        # 内存之外的历史（平台重启后只剩 tasks/*.json）：投影成与在册
        # 任务同口径的稳定键，剔除 log_tail 大字段
        mem_ids = {o.get('id') for o in outs}
        for d in tm.list_archived(limit=limit):
            if d.get('id') in mem_ids:
                continue
            archived.append({k: d.get(k) for k in
                             ('id', 'name', 'status', 'stage', 'pct', 'msg',
                              'error', 'out_dir', 'started', 'finished')})
        archived_count = tm.count_archived()
    return {'count': len(outs), 'tasks': outs,
            'archived_count': archived_count, 'archived': archived}


# ======================================================================
# 工具卡参数 schema（键, 类型, 是否必填, 说明）——口径源自
# web/tool_jobs.py / web/chains.py 各工厂的取参与校验逻辑，与网页端同源。
# 未列出的键会被工厂忽略；缺必填项在提交时当场报中文原因。
# ======================================================================

_TOOL_PARAM_SCHEMAS = {
    'convert': [
        ('input', 'path', True, '输入文件（.sra 或 .fastq[.gz]）'),
        ('input2', 'path', False, 'R2 FASTQ（双端时给）'),
        ('target', 'str', False, '目标格式 fastq|fasta，默认 fastq（仅 FASTQ 输入可转 fasta）'),
    ],
    'fastp': [
        ('r1', 'path', True, 'R1 FASTQ'),
        ('r2', 'path', False, 'R2 FASTQ'),
        ('dedup', 'bool', False, '质控去重'),
    ],
    'hostremoval': [
        ('r1', 'path', True, 'R1 FASTQ'), ('r2', 'path', False, 'R2 FASTQ'),
        ('confidence', 'float', False, 'kunpeng 分类置信度，默认 0'),
        ('db', 'path', False, '宿主库目录，缺省平台主宿主库'),
        ('classify_r1', 'path', False, '预转换 R1 FASTA（跳过质量行解压，可选）'),
        ('classify_r2', 'path', False, '预转换 R2 FASTA'),
    ],
    'orf': [
        ('fasta', 'path', True, '输入 FASTA（核酸 contigs/基因组）'),
        ('min_aa', 'int', False, '最短 ORF 氨基酸数，默认 100（1~5000）'),
        ('annotate', 'bool', False, '连跑 ⑥b 功能注释'),
        ('orf_tool', 'str', False, 'ORF 预测引擎（透传引擎选择，缺省自动）'),
        ('model', 'str', False, 'CDS 模型 pyrodigal_rv|pyrodigal'),
        ('engine', 'str', False, '注释引擎（缺省自动）'),
        ('db', 'str', False, '注释库路径'),
        ('libs', 'str', False, '注释层 prot|pfam，空=全部'),
        ('min_pident', 'float', False, '层1 Identity% 门槛，0=不过滤'),
        ('min_qcov', 'float', False, '层1 覆盖% 门槛'),
    ],
    'orfa': [
        ('run', 'str', False, '已有 orf_ 运行名（与 fasta 二选一）'),
        ('fasta', 'path', False, 'FASTA 输入（预测+注释一步）'),
        ('model', 'str', False, 'pyrodigal_rv|pyrodigal'),
        ('engine', 'str', False, '注释引擎'), ('db', 'str', False, '注释库'),
        ('libs', 'str', False, 'prot|pfam，空=全部'),
        ('min_pident', 'float', False, '层1 Identity% 门槛'),
        ('min_qcov', 'float', False, '层1 覆盖% 门槛'),
    ],
    'genoplot': [
        ('fasta', 'path', False, 'FASTA（与 ann 至少给一）'),
        ('ann', 'path', False, 'GFF3（配 fasta）或 GenBank（自带注释，忽略 fasta）'),
        ('max_plots', 'int', False, '最多出图数，默认 12'),
        ('mode', 'str', False, 'circular|linear|both，默认 both'),
        ('gb_opts', 'dict', False, 'gbdraw CLI 透传参数字典'),
    ],
    'primer': [
        ('fasta', 'path', True, '输入 FASTA（plain=基因组全长分窗；conserved=已比对 MSA）'),
        ('mode', 'str', False, 'plain|conserved，默认 plain'),
        ('num_return', 'int', False, '每条返回引物对数，默认 3（1~20）'),
        ('specificity', 'bool', False, '宿主特异性 BLAST 检查（需先建宿主 BLAST 库）'),
    ],
    'identify': [
        ('input', 'path', True, '输入文件（FASTQ 或 contigs FASTA）'),
        ('input_type', 'str', False, 'pe|single|fasta，默认 pe'),
        ('input2', 'path', False, 'R2 FASTQ（pe 必填）'),
        ('confidence', 'float', False, '分类置信度，默认 0'),
        ('db_virus', 'str', False, '病毒库：内置键 virus|refvirus 或库目录路径'),
    ],
    'assemble': [
        ('r1', 'path', True, 'R1 FASTQ'),
        ('r2', 'path', False, 'R2 FASTQ（srna 模式可空，其余必填）'),
        ('mode', 'str', False, 'SPAdes 模式，默认 rnaviral（srna 单端组装）'),
        ('memory', 'int', False, '内存上限 GB，默认 64'),
        ('min_len', 'int', False, 'contig 最小长度，默认 200'),
    ],
    'contigs': [
        ('contigs', 'path', True, 'contigs FASTA'),
        ('min_len', 'int', False, '最小长度，默认 200'),
        ('confidence', 'float', False, '分类置信度'),
        ('db_virus', 'str', False, '病毒库：virus|refvirus 或库目录'),
    ],
    'structcmp': [
        ('seqs', 'path', True, '序列 FASTA（≥2 条）'),
        ('max_n', 'int', False, '最多条数，默认 30'),
    ],
    'verify': [
        ('run', 'str', False, 'contig 分类运行名（与 fasta 二选一）'),
        ('fasta', 'path', False, '独立候选 FASTA'),
        ('host', 'str', False, '宿主筛选 all|...，默认 all'),
        ('methods', 'str', False, '证据方法组合 blastx,bfvd，默认全开'),
        ('combine', 'str', False, 'union|intersection，默认 union'),
    ],
    'consensus': [
        ('run', 'str', False, 'contig 分类运行名（参考来源）'),
        ('fasta', 'path', False, '独立参考 FASTA'),
        ('kv_run', 'str', False, 'kvsuite 运行名（配 kv_refs 选参考）'),
        ('kv_refs', 'str', False, 'accession 列表，逗号分隔'),
        ('reads', 'str', False, '回贴 reads 文件/目录，逗号分隔；缺省自动取该运行输入'),
        ('ambig', 'str', False, '模糊碱基字符，默认 N'),
        ('min_qual', 'int', False, '碱基质量门槛，默认 20'),
        ('min_depth', 'int', False, '最小深度，默认 4'),
        ('min_freq', 'float', False, '共识频率门槛，默认 0.5'),
        ('min_mapq', 'int', False, '比对质量门槛，默认 10'),
        ('min_minor_freq', 'float', False, '次要变异频率门槛，默认 0.02'),
        ('min_cov_pct', 'float', False, '最小覆盖度%，默认 10'),
    ],
    'kvsuite': [
        ('stage', 'str', False, 'identify+filter（默认）|all|index|identify|filter|consensus|plot|variant'),
        ('samples', 'str', False, '样品名，逗号分隔（与 reads/sample_sheet 三选一）'),
        ('sample_sheet', 'path', False, '样本表 TSV（name,r1,r2）'),
        ('reads', 'str', False, 'FASTQ 目录/文件，逗号分隔（自动配对 R1/R2）'),
        ('kv_lib', 'str', False, '内置鉴定库键（选库优先，忽略零散 reference/index）'),
        ('reference', 'path', False, '参考序列 FASTA（未选库时）'),
        ('ref_info', 'path', False, '参考注释 TSV'),
        ('engine', 'str', False, 'salmon|minibwa，默认 salmon'),
        ('index_dir', 'path', False, '鉴定库目录（未选库时）'),
        ('min_cov', 'float', False, '最小覆盖%，默认 10'),
        ('min_depth', 'float', False, '最小深度，默认 0.5'),
        ('min_reads', 'int', False, '最小 reads 数，默认 10'),
        ('min_poisson', 'float', False, '泊松门槛，默认 0.3'),
        ('variant_qual', 'float', False, '变异质量门槛，默认 3.5'),
        ('min_freq', 'float', False, '变异频率门槛，默认 0.05'),
        ('aa_label_cutoff', 'float', False, 'AA 标签阈值，默认 0.5'),
        ('max_aa_labels', 'int', False, 'AA 标签上限，默认 40'),
        ('input_vcf', 'path', False, '变异段：外部 VCF（优先于 bam_dir）'),
        ('bam_dir', 'path', False, '变异段：BAM 目录'),
        ('ncbi_email', 'str', False, 'NCBI Entrez 邮箱'),
        ('ncbi_api_key', 'str', False, 'NCBI API key（仅环境变量传递，不进日志）'),
        ('no_genes', 'bool', False, '跳过基因注释'),
        ('no_variant_evo', 'bool', False, '跳过变异演化'),
        ('all_variants', 'bool', False, '输出全部变异'),
        ('no_snpgenie', 'bool', False, '跳过 SNPGenie'),
    ],
    'quicktree': [
        ('seqs', 'path', True, '序列 FASTA'),
        ('method', 'str', False, 'nj|fasttree，默认 nj'),
        ('aligned', 'bool', False, '输入已是比对结果（跳过 MAFFT）'),
        ('max_n', 'int', False, '最多条数，默认 100'),
    ],
    'align': [
        ('seqs', 'path', True, '序列 FASTA（≥2 条，自动判核酸/蛋白）'),
        ('seqs_extra', 'path', False, '第二 FASTA（与 seqs 合并比对）'),
        ('strategy', 'str', False, 'auto|linsi|fast，默认 auto'),
        ('trimal', 'str', False, 'automated1|gappyout|strict|none，默认 automated1'),
        ('max_n', 'int', False, '最多条数，默认 200'),
    ],
    'sdt': [
        ('seqs', 'path', True, '序列 FASTA'),
        ('max_n', 'int', False, '最多条数，默认 30'),
        ('orient', 'bool', False, '比对前定向，默认开'),
        ('seqtype', 'str', False, 'auto|nt|aa，默认 auto'),
        ('palette', 'str', False, '热图配色 sdt|viridis|cividis|coolwarm|magma 等'),
        ('aligned', 'bool', False, '输入为已比对 MSA（直接算公式）'),
    ],
    'identity': [
        ('nt_seqs', 'path', True, '核苷酸 FASTA'),
        ('aa_seqs', 'path', False, '氨基酸 FASTA（空=最长 ORF 翻译）'),
        ('max_n', 'int', False, '最多条数，默认 30'),
        ('palette', 'str', False, '复合热图配色'),
        ('aligned', 'bool', False, '输入为已比对 MSA'),
    ],
    'virchain': [
        ('input', 'path', True, '去宿主序列 FASTA'),
        ('input2', 'path', False, 'R2 FASTA'),
        ('run_identify', 'bool', False, '连跑 ②b 病毒识别段'),
        ('do_verify', 'bool', False, '连跑 ③b 候选验证'),
        ('mode', 'str', False, '组装模式，默认 rnaviral'),
        ('memory', 'int', False, '组装内存 GB，默认 64'),
        ('min_len', 'int', False, 'contig 最小长度，默认 200'),
        ('confidence', 'float', False, '分类置信度'),
        ('host', 'str', False, '验证宿主筛选，默认 all'),
        ('methods', 'str', False, '验证方法 blastx,bfvd'),
        ('combine', 'str', False, 'union|intersection'),
        ('db_virus', 'str', False, '病毒库：virus|refvirus 或库目录'),
    ],
    'kvchain': [],   # 运行时复用 kvsuite（stage 固定 all 全流程）
    'prepchain': [
        ('input', 'path', True, '输入文件（.sra 或 FASTQ[.gz]）'),
        ('input2', 'path', False, 'R2 FASTQ（.sra 输入不需要）'),
        ('confidence', 'float', False, '宿主去除分类置信度'),
        ('dedup', 'bool', False, '质控去重'),
        ('db', 'path', False, '宿主库目录，缺省平台主库'),
    ],
    'annochain': [
        ('fasta', 'path', True, '输入 FASTA（核酸 contigs/基因组）'),
        ('min_aa', 'int', False, '最短 ORF 氨基酸数，默认 100'),
    ],
    'rdp': [
        ('input', 'path', True, '比对 FASTA'),
        ('qc', 'bool', False, '比对质控开关，默认开'),
        ('qc_ref', 'str', False, '质控参考序列名'),
        ('qc_min_len', 'float', False, '最小覆盖长度比例，默认 0.90'),
        ('orf_csv', 'path', False, 'ORF 注释 CSV（可选）'),
    ],
    'phylogeo': [
        ('input', 'path', True, '比对 FASTA'),
        ('meta', 'str', False, '元数据 CSV 路径（性状列）'),
        ('trait', 'str', False, '离散性状列名'),
        ('method', 'str', False, 'nj|fasttree，默认 nj'),
        ('dated_tree', 'path', False, '定年树（可选）'),
        ('coords', 'path', False, '坐标表（样本点层，可选）'),
        ('beast_log', 'path', False, 'BEAST .log（导入）'),
        ('beast_trees', 'path', False, 'BEAST .trees（导入）'),
        ('date_trait', 'str', False, '日期列，默认 year'),
        ('motp_bin', 'float', False, '迁移概率分箱，默认 0'),
    ],
    'dsrna': [
        ('input', 'path', True, '目标序列 FASTA'),
        ('max_mm', 'int', False, 'siRNA 最大错配数'),
        ('off_targets', 'path', False, '脱靶评估库'),
        ('panel', 'str', False, '面板输出选项'),
        ('primer_scope', 'str', False, '引物范围选项'),
        ('isolates', 'str', False, '分离株比对集'),
    ],
    'crispr': [
        ('input', 'path', True, '目标序列 FASTA'),
        ('amp_mode', 'str', False, 'RPA|PCR，默认 RPA'),
        ('profile', 'str', False, '严格度 strict 等，默认 strict'),
        ('probe_type', 'str', False, '探针类型，默认 none'),
        ('nuclease', 'str', False, 'cas12a 默认'),
        ('use_host', 'bool', False, '宿主脱靶检查，默认开'),
        ('group_by_segment', 'bool', False, '按片段分组，默认开'),
        ('near_fasta', 'path', False, '近缘参考 FASTA'),
    ],
    # 2026-10-20：CONAN-SWIFT（Cas3 / I-E）crRNA 设计（t-conan 卡）。
    # 引擎 Virus_Platform_Core/conan_design.py；参数面与网页端表单同源。
    'conan': [
        ('input', 'path', True, '目标 FASTA（多分离物，可多分体）'),
        ('negative', 'path', False, '负组 FASTA（近缘非目标病毒，可选）'),
        ('pam', 'str', False, 'PAM（5′ 前置），默认 AAG'),
        ('spacer_len', 'int', False, 'spacer 长度，默认 32（24–36）'),
        ('seed_len', 'int', False, '种子长度（PAM+1 起），默认 8（4–12）'),
        ('ignored', 'str', False, '种子忽略位（逗号分隔 1 基位置），默认 6'),
        ('max_nonseed_mm', 'int', False, '目标非种子错配上限，默认 4'),
        ('exclude_terminal', 'bool', False, '末端排除开关，默认开'),
        ('terminal_margin', 'int', False, '末端排除 bp，默认 150'),
        ('select_gc', 'bool', False, 'GC 区间过滤开关，默认开'),
        ('gc_min', 'int', False, 'GC 下限 %，默认 35'),
        ('gc_max', 'int', False, 'GC 上限 %，默认 65'),
        ('exclude_homopolymer', 'bool', False, '种子同碱基 run 过滤，默认开'),
        ('homopolymer', 'int', False, '种子同碱基 run 阈值，默认 3'),
        ('group_by_segment', 'bool', False, '按片段标签分组，默认开'),
        ('top_n', 'int', False, '输出候选数，默认 20'),
    ],
    'rtt': [
        ('input', 'path', True, '比对 FASTA'),
        ('meta', 'path', True, '元数据 CSV'),
        ('date_trait', 'str', False, '日期列，默认 year'),
        ('method', 'str', False, 'nj 默认'),
    ],
    'phylodyn': [
        ('mode', 'str', False, 'pipeline|stage，默认 pipeline'),
        ('stage', 'str', False, 'stage 模式下的阶段名'),
        ('from_run', 'str', False, '上一轮专项运行名（断点续跑）'),
        ('input', 'path', False, '比对 FASTA（stage 模式）'),
        ('meta', 'path', False, '元数据 CSV'),
        ('tree', 'path', False, 'Newick 树'),
        ('rename_map', 'path', False, '重命名映射 TSV'),
        ('divergence_tree', 'path', False, '分歧时间树'),
        ('coord_table', 'path', False, '坐标表'),
        ('trait', 'str', False, '离散性状列，默认 region'),
        ('date_col', 'str', False, '日期列名'),
        ('date_trait', 'str', False, '日期解析口径'),
        ('seq_len', 'int', False, '序列长度'),
        ('methods', 'str', False, '方法列表，逗号分隔，默认 treetime'),
        ('reroot', 'str', False, '重新生根方式，默认 least-squares'),
        ('keep_root', 'bool', False, '保留原根'),
        ('relax_sigma', 'float', False, '松弛钟 sigma（0=关）'),
        ('relax_coupling', 'float', False, '松弛钟耦合'),
        ('drt', 'bool', False, 'DRT 回归，默认开'),
        ('drt_perm', 'int', False, 'DRT 置换次数，默认 200'),
        ('clock_rate', 'float', False, '固定速率 1e-12~1.0（与松弛钟互斥）'),
        ('skyline', 'bool', False, '天际线'),
        ('n_skyline', 'int', False, '天际线块数，默认 10'),
        ('gen_per_year', 'float', False, '每年世代数，默认 50'),
        ('prune', 'bool', False, '剪枝'),
        ('prune_method', 'str', False, 'fps 等'),
        ('prune_date_resolution', 'str', False, 'month 等'),
        ('prune_max_reps', 'int', False, '默认 3'),
        ('prune_min_snp_diff', 'int', False, '默认 2'),
        ('prune_clade_cutoff', 'int', False, '默认 5'),
        ('sampling_bias_correction', 'str', False, 'auto 等'),
        ('cross_check', 'bool', False, '交叉校验，默认开'),
        ('rssp_bs', 'int', False, 'RTT bootstrap，默认 0'),
        ('seed', 'int', False, '随机种子'),
        ('ancestral', 'bool', False, '祖先状态重构'),
        ('homoplasy', 'bool', False, '同塑性检查'),
        ('timeout', 'int', False, '超时秒，默认 3600'),
    ],
    'pdprep': [
        ('source', 'str', False, 'manual|gb|ncbi，默认 manual'),
        ('fasta_text', 'str', False, '粘贴 FASTA（manual）'),
        ('meta_text', 'str', False, '粘贴元数据（manual）'),
        ('species', 'str', False, '物种名（ncbi）'),
        ('taxid', 'str', False, 'TaxID（ncbi）'),
        ('full_only', 'bool', False, '仅完整基因组'),
        ('date_mode', 'str', False, '日期解析 mid 等'),
        ('date_from', 'str', False, '日期下限'),
        ('date_to', 'str', False, '日期上限'),
        ('collection', 'str', False, 'GenBank 集合名（gb）'),
        ('gb_files', 'path', False, 'GenBank 文件路径'),
        ('accessions', 'str', False, 'accession 列表'),
        ('term', 'str', False, 'Entrez 检索式'),
        ('max_records', 'int', False, '最多条数，默认 200'),
        ('dedup', 'bool', False, '去冗余'),
        ('strict', 'bool', False, '严格元数据校验'),
        ('govern', 'bool', False, '治理规则，默认开'),
        ('allow_noloc', 'bool', False, '允许无定位'),
        ('qc', 'str', False, '比对质控 off 等，默认 off'),
        ('qc_max_gap', 'float', False, '默认 0.10'),
        ('qc_max_n', 'float', False, '默认 0.05'),
        ('qc_min_len', 'float', False, '默认 0.90'),
        ('threads', 'int', False, '比对线程，默认 4'),
        ('coords', 'str', False, '坐标表'),
        ('sub', 'bool', False, '抽样'),
        ('sub_mode', 'str', False, 'random 等'),
        ('sub_seed', 'int', False, '抽样种子'),
        ('sub_n', 'int', False, '抽样数'),
        ('sub_region', 'str', False, '抽样区域'),
    ],
    'mirna': [
        ('mirna_text', 'str', False, '粘贴 miRNA 序列（与 mirna_fasta 二选一）'),
        ('virus_text', 'str', False, '粘贴病毒序列（与 virus_fasta 二选一）'),
        ('mirna_fasta', 'path', False, 'miRNA FASTA'),
        ('virus_fasta', 'path', False, '病毒 FASTA'),
        ('full_mode', 'bool', False, '全模式'),
        ('algorithms', 'str', False, '算法列表，逗号/空格分隔'),
    ],
    'beastexport': [
        ('input', 'path', True, '比对 FASTA'),
        ('meta', 'str', False, '元数据 CSV 路径'),
        ('trait', 'str', False, '离散性状列，默认 region'),
        ('date_trait', 'str', False, '日期列，默认 year'),
        ('virus_tag', 'str', False, '病毒名标签'),
        ('iqtree_log', 'path', False, 'IQ-TREE 日志（定年用，可选）'),
        ('coords', 'path', False, '坐标表'),
        ('prior', 'str', False, '先验 bdsky 等，默认 bdsky'),
    ],
}
_TOOL_PARAM_SCHEMAS['kvchain'] = [
    (k, t, r, ('stage 固定 all；' + d if k == 'stage' else d))
    for (k, t, r, d) in _TOOL_PARAM_SCHEMAS['kvsuite']]


def _tool_list_impl():
    out = []
    for name in sorted(TOOL_REGISTRY):
        spec = TOOL_REGISTRY[name]
        doc = (spec['job'].__doc__ or '').strip()
        params = [{'key': k, 'type': t, 'required': bool(r), 'desc': d}
                  for (k, t, r, d) in _TOOL_PARAM_SCHEMAS.get(name, ())]
        out.append({'tool': name, 'title': spec['title'],
                    'title_en': spec['title_en'],
                    'need_db': bool(spec['need_db']),
                    'weight': 'light' if name in LIGHT_TOOLS else 'heavy',
                    'gated': name in _GATED,
                    'required_params': [p['key'] for p in params
                                        if p['required']],
                    'params': params, 'doc': doc[:400]})
    return {'count': len(out),
            'note': 'params 键值口径与网页端同源；未列出的键被忽略，'
                    '缺必填项提交时当场报中文原因。gated=true 的卡已被'
                    'deferred.json/环境变量下线，调用会被拒，勿选',
            'tools': out}


def _mtime_sorted_dirs(root, skip_underscore=True):
    """目录名按 mtime 倒序（与 /api/tool/runs 同口径）。"""
    try:
        names = [n for n in os.listdir(root)
                 if not n.startswith('.')
                 and not (skip_underscore and n.startswith('_'))
                 and os.path.isdir(os.path.join(root, n))]
    except OSError:
        return []
    mt = []
    for n in names:
        try:
            mt.append((os.path.getmtime(os.path.join(root, n)), n))
        except OSError:
            mt.append((0, n))
    mt.sort(reverse=True)
    return [n for _m, n in mt]


def _sample_list_impl(limit, offset):
    root = DIRS['results']
    names = _mtime_sorted_dirs(root, skip_underscore=False)
    total = len(names)
    lo = max(0, int(offset))
    page = names[lo:lo + max(1, min(int(limit), 100))]
    from Virus_Platform_Core.pipeline import load_sample_input
    out = []
    for n in page:
        d = os.path.join(root, n)
        info = {'sample': n, 'dir': d}
        try:
            r1, r2, project = load_sample_input(d)
            info['r1'] = r1 or ''
            info['paired'] = bool(r2)
            if project:
                info['project'] = project
        except Exception:
            pass
        info['has_manifest'] = os.path.isfile(os.path.join(d, 'project.json'))
        out.append(info)
    return {'total': total, 'count': len(out), 'samples': out}


def _run_list_impl(limit, offset):
    root = tool_runs_root()
    names = _mtime_sorted_dirs(root)
    total = len(names)
    lo = max(0, int(offset))
    page = names[lo:lo + max(1, min(int(limit), 30))]
    out = []
    for n in page:
        try:
            d = run_real_dir(n)
        except Exception:
            d = os.path.join(root, n)
        if not os.path.isdir(d):
            continue
        files = []
        for cur, _sub, fns in os.walk(d):
            rel = os.path.relpath(cur, d)
            for fn in fns:
                files.append(fn if rel == '.' else
                             os.path.join(rel, fn).replace(os.sep, '/'))
        out.append({'run': n, 'out_dir': d, 'n_files': len(files),
                    'files': sorted(files)[:40],
                    'files_truncated': len(files) > 40})
    return {'total': total, 'count': len(out), 'runs': out}


_BINARY_EXT = ('.gz', '.zip', '.png', '.jpg', '.jpeg', '.pdf', '.bam',
               '.bai', '.ttf', '.exe', '.dll', '.woff', '.woff2',
               '.7z', '.tar', '.bin', '.h5')


def _read_run_file_impl(run, path, max_kb):
    if not re.fullmatch(r'[A-Za-z0-9_\-]+', run):
        raise ValueError('无效的运行名: %r' % run)
    base = run_real_dir(run)
    if not os.path.isdir(base):
        raise ValueError('运行不存在: %s（用 vp_list_tool_runs 查看现有运行）'
                         % run)
    rel = (path or '').strip().replace('\\', '/').lstrip('/')
    if not rel or '..' in rel.split('/'):
        raise ValueError('无效的产物相对路径: %r' % path)
    target = os.path.normpath(os.path.join(base, rel))
    base_n = os.path.normcase(os.path.abspath(base))
    if not os.path.normcase(os.path.abspath(target)).startswith(base_n + os.sep):
        raise ValueError('产物路径越界（必须在该运行目录内）')
    target = check_path(target, must_exist=True, in_platform=True)
    if str(target).lower().endswith(_BINARY_EXT):
        raise ValueError('%s 是二进制产物，本工具只读文本（tsv/json/log/'
                         'fasta/nwk/svg 等）；二进制请直接取用磁盘路径: %s'
                         % (os.path.basename(target), target))
    cap = max(1, min(int(max_kb), 4096)) * 1024
    size = os.path.getsize(target)
    with open(target, 'r', encoding='utf-8', errors='replace') as f:
        text = f.read(cap)
    truncated = size > cap
    if truncated:
        text += '\n...[已截断：文件 %d 字节，上限 %d KB]' % (size, cap // 1024)
    return _fmt(
        'ok', f"{rel}（{size} 字节{'，已截断' if truncated else ''}）",
        result={'run': run, 'path': rel, 'abs': target, 'size': size,
                'truncated': truncated, 'text': text})


def _platform_info_impl():
    info = {'server': 'virusplatform_mcp', 'version': SERVER_VERSION,
            'platform_root': PLATFORM_ROOT,
            'python': sys.executable, 'results_dir': DIRS['results'],
            'tool_runs_dir': tool_runs_root(),
            'databases_dir': DIRS['databases'],
            'threads_default': getattr(cfg, 'threads', None),
            'language': getattr(cfg, 'lang', 'zh')}
    enabled = _enabled_tool_names()
    info['tools_enabled'] = enabled
    info['tools_count'] = len(enabled)
    info['tools_gated'] = sorted(_GATED)
    info['fingerprint'] = _capability_fingerprint(enabled)
    try:
        from Virus_Platform_Core.kunpeng import db_ready
        info['host_db_ready'] = bool(db_ready(cfg.databases['host']))
        info['virus_db_ready'] = bool(db_ready(cfg.databases['virus']))
    except Exception as e:
        info['db_check_error'] = repr(e)
    problems = _platform_check()
    info['ready'] = not problems
    if problems:
        info['problems'] = problems
    if problems:
        info['summary'] = '平台模块异常：' + '; '.join(problems)
    else:
        db_ok = info.get('host_db_ready') and info.get('virus_db_ready')
        info['summary'] = (f"平台就绪（MCP 工具 {info['tools_count']} 个；"
                           f"{'库全部就绪' if db_ok else '注意：库未就绪'}）")
    return info


# ======================================================================
# MCP 工具层（签名即 schema；Annotated+Field 出约束，docstring 即描述）
# ======================================================================

@mcp.tool(name='vp_platform_info', annotations=_RO)
async def vp_platform_info() -> dict:
    """平台环境总览：根目录、结果/运行/数据库目录、宿主与病毒库就绪状态、
    默认线程数、就绪体检结论；并带版本号、启用工具清单与能力指纹
    （fingerprint 变了即服务端能力有增删）。

    返回 {server, version, platform_root, results_dir, tool_runs_dir,
    databases_dir, host_db_ready, virus_db_ready, threads_default,
    language, ready, tools_enabled, tools_count, tools_gated,
    fingerprint}。首次接入先调本工具确认环境，再选后续动作。
    """
    return await anyio.to_thread.run_sync(_platform_info_impl)

@mcp.tool(name='vp_list_tools', annotations=_RO)
async def vp_list_tools() -> dict:
    """工具箱全量清单：名称、中英文标题、是否依赖病毒库、heavy/light
    档位、参数口径说明（doc）。

    返回 {count, note, tools:[{tool, title, title_en, need_db, weight,
    doc}]}。提交运行用 vp_run_tool(tool=..., params={...})。
    """
    return await anyio.to_thread.run_sync(_tool_list_impl)

@mcp.tool(name='vp_run_tool')
async def vp_run_tool(
        tool: Annotated[str, Field(min_length=1, max_length=40,
                                   description='工具名（vp_list_tools 的 '
                                               'tool 键）')],
        params: Annotated[Optional[dict], Field(
            description='工具参数字典（键值口径见 vp_list_tools 的 doc；'
                        '文件路径可用绝对路径或平台根相对路径）')] = None,
        threads: Annotated[Optional[int], Field(
            ge=1, le=256, description='线程数（缺省用平台默认）')] = None,
        out_dir: Annotated[Optional[str], Field(
            max_length=500,
            description='可选自定义产物父目录（产物写 <out_dir>/<run_name>/）'
                        '；缺省写 run/tool_runs/<run_name>/')] = None) -> dict:
    """执行工具箱任一工具（格式转换/质控/宿主去除/组装/鉴定定量/建树/SDT/
    引物/基因组图/dsRNA/CRISPR/一键链……全量）。

    提交即返回 {task_id, run, out_dir, weight, hint}，不同步等待。
    必填参数缺失/库未就绪当场报中文原因。轮询用 vp_get_task，
    产物文本用 vp_read_run_file(run=...)。
    """
    return await anyio.to_thread.run_sync(
        _submit_tool_run, tool, params or {}, threads, out_dir)

@mcp.tool(name='vp_pipeline_stages', annotations=_RO)
async def vp_pipeline_stages() -> dict:
    """流水线阶段清单：15 个阶段键、中英文名、默认全流程序列、功能分组。

    返回 {stage_order, default_stages, names, names_en,
    groups:[{group, stages}]}。阶段键喂 vp_pipeline_start 的 stages。
    """
    return {'stage_order': list(STAGE_ORDER),
            'default_stages': list(DEFAULT_ANALYZE_STAGES),
            'names': dict(STAGE_NAMES), 'names_en': dict(STAGE_NAMES_EN),
            'groups': [{'group': g, 'stages': list(ss)}
                       for g, ss in STAGE_GROUPS]}


@mcp.tool(name='vp_pipeline_start')
async def vp_pipeline_start(
        sample: Annotated[str, Field(min_length=1, max_length=120,
                                     description='样品名（自动安全化；产物'
                                                 '目录 run/results/<样品名>/）')],
        r1: Annotated[str, Field(min_length=1, max_length=1000,
                                 description='R1 FASTQ(.gz) 路径（绝对或平台'
                                             '根相对）')],
        r2: Annotated[Optional[str], Field(
            max_length=1000, description='R2 FASTQ(.gz)（双端必填）')] = None,
        stages: Annotated[Optional[list[str]], Field(
            description='阶段键列表（vp_pipeline_stages 的键）；缺省跑默认 '
                        '13 步全流程')] = None,
        threads: Annotated[Optional[int], Field(ge=1, le=256)] = None,
        db_host: Annotated[Optional[str], Field(
            max_length=1000, description='宿主库目录（缺省平台主宿主库）')] = None,
        db_virus: Annotated[Optional[str], Field(
            max_length=1000, description='病毒库目录（缺省平台主病毒库）')] = None,
        confidence: Annotated[float, Field(
            ge=0, le=1, description='kunpeng 分类置信度')] = 0.0,
        assembly_mode: Annotated[str, Field(
            max_length=20, description='SPAdes 模式 '
                                       '(rnaviral/srna/metaviral/rna...)')] =
        'rnaviral',
        memory_gb: Annotated[int, Field(
            ge=1, le=1024, description='SPAdes 内存上限 GB')] = 64,
        min_contig_len: Annotated[int, Field(ge=1, le=100000)] = 200,
        assembly_input: Annotated[str, Field(
            max_length=20,
            description='③ 组装输入 virus=病毒筛选后 reads（默认）|'
                        'all=全部去宿主 reads')] = 'virus',
        min_viral_pairs: Annotated[int, Field(
            ge=0, le=10**9,
            description='③ 病毒 reads 低于该对数时回退全部去宿主 reads，'
                        '默认 500')] = 500,
        do_trim: Annotated[bool, Field(
            description='⑧ 引物设计前比对清剪，默认开')] = True,
        subsample: Annotated[int, Field(
            ge=0, description='子采样读对数（0=不子采样）')] = 0,
        tree_tool: Annotated[str, Field(max_length=20)] = 'fasttree',
        primer_mode: Annotated[str, Field(max_length=20)] = 'conserved',
        min_orf_aa: Annotated[int, Field(ge=1, le=5000)] = 100,
        ncbi_refs: Annotated[Optional[str], Field(
            max_length=100,
            description='NCBI 参考集合名（main.py ncbi-dl 先下载；进化树'
                        '扩充用）')] = None,
        kv_lib: Annotated[Optional[str], Field(
            max_length=300, description='已知病毒鉴定库键或库目录（②b 段；'
                                        '缺省平台默认库）')] = None,
        subsample_mode: Annotated[str, Field(max_length=20)] = 'prefix',
        subsample_seed: Annotated[int, Field(ge=0)] = 42,
        subsample_prop: Annotated[float, Field(ge=0.0, le=1.0)] = 0.0,
        tree_sampling: Annotated[str, Field(
            max_length=20,
            description='建树序列取样 blast 等，默认 blast')] = 'blast',
        top_n_refs: Annotated[int, Field(
            ge=1, le=100, description='进化分析取近缘参考数，默认 10')] = 10,
        min_pident: Annotated[Optional[float], Field(
            ge=0, le=100, description='⑥b 层1 Identity% 门槛')] = None,
        min_qcov: Annotated[Optional[float], Field(
            ge=0, le=100, description='⑥b 层1 覆盖% 门槛')] = None,
        fastp_dedup: Annotated[bool, Field(description='⓪ 质控去重')] = False,
        do_fq2fa: Annotated[bool, Field(
            description='⓪b FASTQ→FASTA 转换，默认开')] = True,
        gbdraw_max: Annotated[int, Field(
            ge=1, le=200, description='⑨ 最多出图数，默认 12')] = 12,
        gbdraw_fasta: Annotated[Optional[str], Field(
            max_length=1000, description='⑨ 覆盖出图 FASTA（缺省用 ③ 产物）')] = None,
        gbdraw_ann: Annotated[Optional[str], Field(
            max_length=1000, description='⑨ 覆盖注释 GFF3/GenBank')] = None,
        plot_engine: Annotated[str, Field(
            max_length=20, description='绘图引擎 auto/gbdraw/dfv')] = 'auto',
        verify_host: Annotated[str, Field(
            max_length=60, description='③b 宿主筛选 all 等，默认 all')] = 'all',
        verify_methods: Annotated[str, Field(
            max_length=60, description='③b 证据方法 blastx,bfvd（逗号分隔）')] =
        'blastx,bfvd',
        verify_combine: Annotated[str, Field(
            max_length=20, description='③b 合并方式 union/intersection')] =
        'union',
        force: Annotated[bool, Field(description='强制重跑全部阶段')] = False,
        specificity: Annotated[bool, Field(
            description='引物宿主特异性 BLAST 检查')] = False) -> dict:
    """启动样品全流程分析（⓪质控→①宿主去除→②b识别定量→③组装→…→⑩报告）。

    与网页端/CLI 同一编排引擎（断点续跑、可选阶段自动降级都保留）。
    提交即返回 {task_id, sample, sample_dir, stages, hint}；轮询
    vp_get_task 看阶段进度（stage/pct/msg），终态 result 为各阶段摘要。
    """
    args = {k: v for k, v in dict(
        sample=sample, r1=r1, r2=r2, stages=stages, threads=threads,
        db_host=db_host, db_virus=db_virus, confidence=confidence,
        assembly_mode=assembly_mode, memory_gb=memory_gb,
        min_contig_len=min_contig_len, assembly_input=assembly_input,
        min_viral_pairs=min_viral_pairs, do_trim=do_trim,
        subsample=subsample,
        tree_tool=tree_tool, primer_mode=primer_mode,
        min_orf_aa=min_orf_aa, ncbi_refs=ncbi_refs, force=force,
        specificity=specificity).items() if v is not None}
    return await anyio.to_thread.run_sync(_submit_pipeline, args)

@mcp.tool(name='vp_get_task', annotations=_RO_LIVE)
async def vp_get_task(
        task_id: Annotated[str, Field(
            min_length=4, max_length=32,
            description='vp_run_tool / vp_pipeline_start 返回的 task_id')],
        log_lines: Annotated[int, Field(
            ge=0, le=500, description='随快照返回的日志尾行数')] = 30,
        result_max_kb: Annotated[int, Field(
            ge=1, le=4096,
            description='终态 result 的文本上限 KB，超出截断为 '
                        '{_truncated,_note,_head}')] = 256) -> dict:
    """查单个任务快照：状态（queued/running→终态 done/failed/cancelled）、
    阶段、百分比、消息、日志尾、终态结果。

    返回固定键集合 {id, name, status, stage, pct, msg, error, eta,
    out_dir, link, started_h, finished_h, log, result}——稳定投影，
    键恒在；error 非空即失败原因。result 超限时截断并标注完整位置。
    """
    return await anyio.to_thread.run_sync(
        _task_snapshot, task_id, log_lines, result_max_kb)

@mcp.tool(name='vp_list_tasks', annotations=_RO)
async def vp_list_tasks(
        limit: Annotated[int, Field(ge=1, le=100)] = 20,
        include_archived: Annotated[bool, Field(
            description='并入磁盘归档任务（平台重启前的历史，'
                        'run/tasks/*.json）')] = False) -> dict:
    """近期任务列表（新到旧，分页）：状态/进度/产物目录，不带日志。

    返回 {count, tasks:[...]}；include_archived=true 时另带 archived
    （稳定键投影）与 archived_count。逐任务日志用 vp_get_task /
    vp_get_task_log（两者对归档任务均可用）。
    """
    return await anyio.to_thread.run_sync(
        _task_list, limit, include_archived)

@mcp.tool(name='vp_cancel_task', annotations=_DESTR)
async def vp_cancel_task(
        task_id: Annotated[str, Field(min_length=4, max_length=32)]) -> dict:
    """取消运行中/排队中的任务（对当前外部步骤硬停止进程树）。

    返回 {task_id, accepted, status, msg}。accepted=false 表示任务
    不存在或已结束；取消后用 vp_get_task 确认终态。
    """

    def _impl():
        ok = tm.cancel(task_id)
        snap = tm.snapshot(task_id, log_lines=0, include_result=False) or {}
        return {'task_id': task_id, 'accepted': bool(ok),
                'status': snap.get('status') or 'unknown',
                'msg': snap.get('msg') or ''}
    return await anyio.to_thread.run_sync(_impl)

@mcp.tool(name='vp_get_task_log', annotations=_RO)
async def vp_get_task_log(
        task_id: Annotated[str, Field(min_length=4, max_length=32)],
        lines: Annotated[int, Field(ge=1, le=2000)] = 400) -> dict:
    """任务全量日志尾部（含已完成/归档任务：内存→磁盘 log→归档 json
    三级回源）。

    返回 {task_id, lines, log:[...]}。
    """

    def _impl():
        got = tm.full_log(task_id, lines=max(1, min(int(lines), 2000)))
        return {'task_id': task_id, 'lines': len(got or []),
                'log': got or []}
    return await anyio.to_thread.run_sync(_impl)

@mcp.tool(name='vp_list_samples', annotations=_RO)
async def vp_list_samples(
        limit: Annotated[int, Field(ge=1, le=100)] = 20,
        offset: Annotated[int, Field(ge=0)] = 0) -> dict:
    """已接入样品清单（run/results/，按最近变动倒序，分页）：输入文件、
    是否双端、所属项目、是否有 project.json 清单。

    返回 {total, count, samples:[{sample, dir, r1, paired, project,
    has_manifest}]}。样品输入即 vp_pipeline_start 所需的 r1/r2 参考。
    """
    return await anyio.to_thread.run_sync(
        _sample_list_impl, limit, offset)

@mcp.tool(name='vp_list_tool_runs', annotations=_RO)
async def vp_list_tool_runs(
        limit: Annotated[int, Field(ge=1, le=30)] = 10,
        offset: Annotated[int, Field(ge=0)] = 0) -> dict:
    """近期工具/链式运行清单（run/tool_runs/，按最近变动倒序，分页）：
    产物目录与文件清单（最多 40 个文件名）。

    返回 {total, count, runs:[{run, out_dir, n_files, files,
    files_truncated}]}。读具体文件内容用 vp_read_run_file。
    """
    return await anyio.to_thread.run_sync(
        _run_list_impl, limit, offset)

@mcp.tool(name='vp_read_run_file', annotations=_RO)
async def vp_read_run_file(
        run: Annotated[str, Field(min_length=1, max_length=80,
                                  description='运行名（vp_list_tool_runs '
                                              '返回的 run）')],
        path: Annotated[str, Field(min_length=1, max_length=500,
                                   description='运行目录内相对路径，如 '
                                               'kvsuite/filter/filtered.tsv')],
        max_kb: Annotated[int, Field(
            ge=1, le=4096, description='返回文本上限 KB（超出截断并标注）')] =
        512) -> dict:
    """读某次运行的文本产物（tsv/json/log/fasta/nwk/svg 等）内容。

    路径限定在该运行目录内（防越界）；二进制产物明确拒收并提示磁盘
    路径。返回 {run, path, abs, size, truncated, text}。
    """
    return await anyio.to_thread.run_sync(
        _read_run_file_impl, run, path, max_kb)

@mcp.tool(name='vp_delete_tool_run', annotations=_DESTR)
async def vp_delete_tool_run(
        run: Annotated[str, Field(min_length=1, max_length=80,
                                  description='运行名（vp_list_tool_runs 的 '
                                              'run）')],
        confirm: Annotated[bool, Field(
            description='必须显式 true 才执行删除（不可逆，连同全部产物）')]) -> dict:
    """删除一次工具/链式运行目录（含全部产物，不可逆）。

    守卫：confirm 不为 true 拒删；该运行上有任务在跑时拒删（先
    vp_cancel_task）。返回 {deleted, run, out_dir}。
    """
    def _impl():
        if not confirm:
            raise ValueError('删除不可逆：请把 confirm 设为 true 再执行'
                             '（将删除该运行目录下全部产物）')
        if not re.fullmatch(r'[A-Za-z0-9_\-]+', run):
            raise ValueError('无效的运行名: %r' % run)
        p = check_path(run_real_dir(run), must_exist=True, in_platform=True)
        _key = os.path.normcase(os.path.abspath(p))
        with tm.lock:
            _busy = [r for r in tm.tasks.values()
                     if r.get('status') in ('running', 'queued')
                     and r.get('out_dir')
                     and os.path.normcase(os.path.abspath(str(r['out_dir'])))
                     == _key]
        if _busy:
            raise ValueError('该运行目录上有任务正在跑（%s %s）；'
                             '请先 vp_cancel_task 再删除'
                             % (_busy[0].get('id'), _busy[0].get('name')))
        import gc as _gc
        import shutil as _sh
        for _ in range(2):          # Windows 句柄回收重试（同网页端口径）
            _sh.rmtree(p, ignore_errors=True)
            if not os.path.isdir(p):
                break
            _gc.collect()
            time.sleep(0.3)
        if os.path.isdir(p):
            raise RuntimeError('删除失败（产物文件被占用）：%s' % p)
        base = os.path.join(tool_runs_root(), run)
        if os.path.abspath(base) != os.path.abspath(p):
            _sh.rmtree(base, ignore_errors=True)
        return {'deleted': True, 'run': run, 'out_dir': str(p)}
    return await anyio.to_thread.run_sync(_impl)


# ======================================================================
# CLI 能力面（数据库构建 / 公共数据 / NCBI 提交 / LOGAN / ICTV / 维护）：
# 两个通用工具 vp_cli_run（快查，同步）+ vp_cli_submit（长任务进 tm），
# 逐命令参数 schema 与 main.py argparse 同源。规格：
#   (键, 类型, 是否必填, 说明)；键以下划线书写，转 argv 时转 --中划线。
# 类型：str/int/float/bool(store_true)/list(重复旗标)/raw(整段透传)
# ======================================================================

CLI_SENSITIVE = {'ncbi_api', 'deepseek_api', 'ncbi_api_key'}


def _C(desc, mode, args=(), positional=(), weight='heavy', confirm=False,
       timeout=180):
    return {'desc': desc, 'mode': mode, 'args': list(args),
            'positional': list(positional), 'weight': weight,
            'confirm': confirm, 'timeout': timeout}


COMMAND_SPECS = {
    # ---- 数据库构建 ----
    'init-taxonomy': _C('下载/准备 NCBI taxonomy（首次必做，耗时长）', 'task', [
        ('force', 'bool', False, '强制重新下载')]),
    'build-host-db': _C('kunpeng 宿主数据库构建（GB 级 I/O，耗时长）', 'task', [
        ('genome', 'path', True, '宿主基因组 FASTA'),
        ('taxid', 'int', True, '宿主 NCBI TaxID'),
        ('hash_capacity', 'str', False, '哈希容量，默认 256M'),
        ('threads', 'int', False, '线程数（默认自动）'),
        ('out_dir', 'path', False, '输出库目录（默认 host-db/<taxid>_<名>_host_db）'),
        ('clean_mid', 'bool', False, '建库后清理中间文件（约省 3GB）'),
        ('rebuild', 'bool', False, '强制重建')]),
    'build-virus-db': _C('kunpeng 病毒数据库构建', 'task', [
        ('fasta', 'path', True, '病毒参考 FASTA'),
        ('info', 'path', True, 'Accession/Taxid 信息表 TSV'),
        ('hash_capacity', 'str', False, '哈希容量，默认 64M'),
        ('threads', 'int', False, '线程数'),
        ('rebuild', 'bool', False, '强制重建')]),
    # ---- 参考数据采集 ----
    'ncbi-dl': _C('NCBI Entrez 批量下载参考序列（进化树扩充；集合名供 '
                  'vp_pipeline_start 的 ncbi_refs 引用）', 'task', [
        ('name', 'str', True, '集合名'),
        ('db', 'str', False, 'nucleotide|protein，默认 nucleotide'),
        ('max', 'int', False, '最多下载条数，默认 100')],
        positional=[('term', 'Entrez 检索式，如 "Tobamovirus[ORGN] AND '
                             'complete genome[TITL]"')]),
    'ncbi-list': _C('列出已下载的 NCBI 参考集合', 'sync', [], timeout=60),
    'gb-dl': _C('下载 GenBank 集合（同属共线性比较输入）', 'task', [
        ('name', 'str', True, '集合名'),
        ('term', 'str', False, 'Entrez 检索式（与 acc 二选一）'),
        ('acc', 'str', False, 'accession 列表（逗号分隔，支持版本号）'),
        ('max', 'int', False, '最多下载数，默认 50')]),
    'gb-import': _C('导入本机 GenBank 文件为集合（多记录自动拆分）', 'task', [
        ('name', 'str', True, '集合名'),
        ('files', 'path', True, '.gb/.gbk 文件路径（逗号分隔，支持 .gz）')]),
    'gb-list': _C('列出 GenBank 集合', 'sync', [], timeout=60),
    'gb-check': _C('巡检 GenBank 集合（记录数/CDS 数/警告）', 'sync', [
        ('name', 'str', True, '集合名')], timeout=120),
    'host-genome': _C('宿主参考基因组下载（NCBI / gget·Ensembl，供建宿主库）',
                      'task', [
        ('species', 'str', False, '物种拉丁名（与 accession 二选一）'),
        ('accession', 'str', False, '组装 accession（GCF_/GCA_）'),
        ('chromosomes', 'str', False, '染色体过滤（如 X,Y；仅 accession 模式）'),
        ('source', 'str', False, 'ncbi|gget，默认 ncbi'),
        ('out', 'path', False, '输出目录'),
        ('release', 'int', False, 'gget: Ensembl release'),
        ('taxid', 'int', False, '指定 TaxID（省去在线查询）'),
        ('include_organelles', 'bool', False, '含叶绿体/线粒体（仅 NCBI 源）'),
        ('ncbi_api', 'str', False, 'NCBI API Key（日志中打码）'),
        ('min_length', 'int', False, '最短序列过滤'),
        ('verify_only', 'bool', False, '仅检查已有下载'),
        ('skip_datasets', 'bool', False, '跳过 datasets CLI 直接 E-utilities')]),
    'ref-status': _C('查看病毒参考库版本与统计', 'sync', [], timeout=60),
    'ictv-update': _C('更新 ICTV VMR 参考库（在线下载或本地 xlsx）', 'task', [
        ('xlsx', 'path', False, '本地 VMR xlsx（缺省在线下载当前版）')]),
    'ictv-status': _C('查看 ICTV VMR 参考库状态', 'sync', [], timeout=60),
    'ictv-refs': _C('ICTV VMR 选参：按属/科/种挑参考（download=true 时'
                    '联网按需下载，耗时更长）', 'sync', [
        ('genus', 'str', False, '属名（如 Tobamovirus）'),
        ('family', 'str', False, '科名（如 Potyviridae）'),
        ('db', 'str', False, 'plant|ictv，默认 plant'),
        ('species', 'str', False, '种名'),
        ('limit', 'int', False, '最多挑多少条，默认 20'),
        ('any', 'bool', False, '不限 Complete genome'),
        ('download', 'bool', False, '本地没有的 accession 联网下载')], timeout=300),
    'datasets': _C('NCBI datasets CLI 下载（genome|gene|taxonomy|virus，'
                   '参数整段透传）', 'task', [], positional=[
        ('raw', '透传参数串，如 "genome --accession GCF_000001405.40 '
                '--chromosomes X,Y"')]),
    'datasets-setup': _C('安装/更新 NCBI datasets CLI（官方 exe）', 'task', []),
    # ---- 公共数据检索 ----
    'meta-search': _C('公共数据检索：SRA+GSA 双引擎按物种检索 Run', 'task', [
        ('species', 'str', True, '物种拉丁名（如 "Lycium barbarum"）'),
        ('source', 'str', False, 'TRANSCRIPTOMIC|GENOMIC|All，默认前者'),
        ('db', 'str', False, 'sra|gsa|both，默认 both'),
        ('out', 'path', False, '输出目录（默认 meta_search/<物种>/search）'),
        ('no_detailed', 'bool', False, '关闭详细模式（Tissue/Location 深提取）'),
        ('ncbi_api', 'str', False, 'NCBI API Key（日志中打码）'),
        ('deepseek_api', 'str', False, 'DeepSeek API Key（AI 元数据清洗；日志打码）')]),
    'meta-info': _C('Run 列表 → Core14/Full 统一元数据', 'task', [
        ('runs', 'str', True, 'Run 编号列表文件或单个编号（SRR/CRR…）'),
        ('out', 'path', False, '输出目录'),
        ('mode', 'str', False, 'local|api|both，默认 both'),
        ('threads', 'int', False, '默认 4'),
        ('fill_date', 'bool', False, '采集日期缺失用发布日期兜底'),
        ('ncbi_api', 'str', False, 'NCBI API Key（打码）'),
        ('deepseek_api', 'str', False, 'DeepSeek API Key（打码）')], weight='light'),
    'meta-plot': _C('元数据 SCI 可视化（时间/机构/组织/地理）', 'task', [
        ('input', 'path', True, '检索结果或 Core14/Full 元数据 CSV'),
        ('out', 'path', False, '图表输出目录，默认 SCI_Figures_Output')],
        weight='light'),
    # ---- NCBI 提交 ----
    'submit-list': _C('提交准备：列出提交项目', 'sync', [], timeout=60),
    'submit-init': _C('提交准备：新建提交项目（unified_metadata.csv）', 'sync', [
        ('name', 'str', True, '项目名（字母/数字/_/-）'),
        ('taxonomy', 'path', False, 'taxonomy.tsv → 自动生成初表'),
        ('metadata', 'path', False, '公共元数据表（Core14/Full）'),
        ('import_csv', 'path', False, '从已有 unified CSV/TSV/Excel 导入'),
        ('demo', 'bool', False, '用示例数据建表'),
        ('authors', 'str', False, '作者 "Last, First; ..."'),
        ('title', 'str', False, '提交标题'),
        ('bioproject', 'str', False, 'BioProject ID'),
        ('host', 'str', False, '宿主物种名'),
        ('lat_lon', 'str', False, '经纬度 "38.47 N 106.27 E"'),
        ('sequencer', 'str', False, '测序平台'),
        ('assembler', 'str', False, '组装方法串'),
        ('coverage', 'str', False, '覆盖度（如 42.5x）')], timeout=120),
    'submit-validate': _C('提交准备：必填字段校验（存在未解决问题时退出码 1）',
                          'sync', [
        ('name', 'str', True, '项目名')], timeout=120),
    'submit-fill': _C('提交准备：批量填充某列占位符/替换旧值', 'sync', [
        ('name', 'str', True, '项目名'),
        ('column', 'str', True, '列名'),
        ('value', 'str', True, '新值'),
        ('old', 'str', False, '旧值（缺省=填充所有占位符）')], timeout=120),
    'submit-export': _C('提交准备：生成 source.src/miuvig/assembly/BioSample/'
                        'report', 'task', [
        ('name', 'str', True, '项目名'),
        ('assembler', 'str', False, '默认 SPAdes;4.3.0;rnaviral'),
        ('sequencer', 'str', False, '默认 Illumina NovaSeq 6000'),
        ('enrichment', 'str', False, '默认 rRNA depletion')], weight='light'),
    'submit-adopt': _C('从 MMPV-RNA discovery 产物一键建提交项目', 'task', [
        ('name', 'str', True, '提交项目名'),
        ('out_dir', 'path', False, 'discovery 输出根'),
        ('from_server', 'str', False, '从服务器拉取（数据集名，需免密 ssh）'),
        ('server_host', 'str', False, 'ssh 目标'),
        ('server_data_root', 'str', False, '服务器数据根'),
        ('dataset', 'str', False, '数据集名（Isolate 前缀）'),
        ('samples', 'path', False, '样本元数据表（需 match 列）'),
        ('run_name', 'str', False, '平台分类运行名'),
        ('force', 'bool', False, '项目已存在时清空重建'),
        ('dry_run', 'bool', False, '只报会生成什么'),
        ('no_export', 'bool', False, '只出 unified_metadata.csv'),
        ('min_length', 'int', False, '过滤短序列'),
        ('isolate_prefix', 'str', False, 'Isolate 前缀'),
        ('public_metadata', 'path', False, '公共元数据表'),
        ('geo_loc', 'str', False, '采集地点 "China:Ningxia"'),
        ('lat_lon', 'str', False, '经纬度'),
        ('collection_date', 'str', False, '采集日期 YYYY[-MM[-DD]]'),
        ('host', 'str', False, '宿主物种名'),
        ('tissue', 'str', False, '组织'), ('cultivar', 'str', False, '品种'),
        ('dev_stage', 'str', False, '发育阶段'),
        ('collected_by', 'str', False, '采集人'),
        ('isolation_source', 'str', False, '分离来源'),
        ('authors', 'str', False, '作者'), ('title', 'str', False, '标题'),
        ('bioproject', 'str', False, 'BioProject'),
        ('biosample_prefix', 'str', False, 'BioSample 前缀'),
        ('assembler', 'str', False, '组装方法'), ('sequencer', 'str', False, '平台'),
        ('coverage', 'str', False, '覆盖度'),
        ('annotation_pipeline', 'str', False, '注释流程串'),
        ('enrichment', 'str', False, '富集方式')]),
    'submit-sbt': _C('提交准备：生成 template.sbt（作者/机构）', 'sync', [
        ('name', 'str', True, '项目名'), ('last', 'str', True, '姓'),
        ('first', 'str', True, '名'), ('middle', 'str', False, '中间名'),
        ('affil', 'str', True, '机构'), ('div', 'str', False, '院系'),
        ('city', 'str', True, '城市'), ('sub', 'str', False, '省/州'),
        ('country', 'str', True, '国家'), ('street', 'str', False, '街道'),
        ('email', 'str', True, '邮箱'), ('postal', 'str', False, '邮编'),
        ('title', 'str', False, '提交标题'),
        ('extra_authors', 'str', False, '其他作者（\\n 分隔 "Last, First"）')],
        timeout=60),
    # ---- LOGAN 溯源 ----
    'logan-create': _C('LOGAN 溯源：生成待提交 Logan-Search 的查询片段', 'task', [
        ('name', 'str', True, '查询名称'),
        ('sample', 'str', False, '来源样品（需已有 ③ 组装 contigs）'),
        ('contigs', 'str', False, 'contig 名，逗号分隔（默认全部）'),
        ('fasta', 'path', False, '或直接给 FASTA 文件'),
        ('segments', 'int', False, '每条 contig 切片段数 1-4，默认 2')]),
    'logan-import': _C('LOGAN 溯源：导入结果表并生成报告', 'task', [
        ('name', 'str', True, '查询名称'),
        ('segment', 'int', False, '片段编号，从 1 起'),
        ('result', 'path', True, '结果表 CSV/TSV')], weight='light'),
    'logan-jobs': _C('LOGAN 溯源：列出查询任务', 'sync', [], timeout=60),
    'logan-batch': _C('LOGAN 溯源：Selenium 批量提交并自动收结果（耗时长）',
                      'task', [
        ('name', 'str', True, '查询名称'),
        ('email', 'str', True, '通知邮箱（逗号分隔轮换）'),
        ('group', 'str', False, 'Fast_No_human 等 9 档，默认 Fast_No_human'),
        ('show_browser', 'bool', False, '显示浏览器（默认 headless）'),
        ('first_wait', 'int', False, '默认 300'),
        ('max_wait', 'int', False, '默认 1800')]),
    # ---- 专项重跑（对已有样品结果，免流水线）----
    'report': _C('对已有样品结果重新生成可视化报告', 'task', [
        ('sample', 'str', True, '样品名'), ('force', 'bool', False, '强制重跑')],
        weight='light'),
    'host-analysis': _C('对已有 ③ 组装结果做 ICTV 宿主预测', 'task', [
        ('sample', 'str', True, '样品名'),
        ('force', 'bool', False, '强制重跑'),
        ('min_confidence', 'float', False, 'ICTV 置信度门槛 0-1（默认不过滤）')],
        weight='light'),
    'orfa': _C('ORF 功能注释重跑（⑥b，对已有样品）', 'task', [
        ('sample', 'str', True, '样品名'), ('threads', 'int', False, '线程'),
        ('libs', 'str', False, 'prot|pfam|bfvd 组合，逗号分隔'),
        ('min_pident', 'float', False, '层1 Identity% 门槛'),
        ('min_qcov', 'float', False, '层1 覆盖% 门槛'),
        ('force', 'bool', False, '强制重跑')], weight='light'),
    # ---- 维护 ----
    'db-migrate': _C('数据库迁移/对接（软件与数据库分离部署；copy/move 落盘'
                     '都会切换平台库根，必须 confirm=true）', 'task', [
        ('to', 'path', True, '目标数据库根目录（平台目录之外）'),
        ('mode', 'str', False, 'copy|move，默认 copy'),
        ('check', 'bool', False, '只校验目标是否已有完整库'),
        ('dry_run', 'bool', False, '只预检（空间/路径）')],
        confirm='apply'),
    'tool-runs': _C('输出目录管理：status/organize/archive/clean（clean 不可逆，'
                    '无 dry_run 时必须 confirm=true）', 'task', [
        ('before', 'str', False, 'archive: 归档该日期 YYYYMMDD 之前的运行'),
        ('older_days', 'int', False, 'archive: 归档早于 N 天的运行'),
        ('active_days', 'int', False, 'clean: 删除活动区早于 N 天的运行'),
        ('archive_days', 'int', False, 'clean: 删除归档区早于 N 天的归档'),
        ('dry_run', 'bool', False, '只预演，不动盘')],
        positional=[('action', 'status|organize|archive|clean')],
        confirm='clean'),
    'samples-backfill': _C('给缺 project.json 的老样品回填清单', 'task', [
        ('sample', 'list', False, '只处理指定样品（可重复）'),
        ('dry_run', 'bool', False, '只报告不落盘')], weight='light'),
    'selfcheck': _C('平台环境自检（模块/工具/数据库/磁盘/依赖）', 'sync', [],
                    timeout=300),
}

def _fmt(status, summary, result=None, interpret=None, nxt=None, **extra):
    """统一输出信封：一句话摘要 + 结构化结果 + 规则化解读 + 下一步指引。

    解读全部由真实字段**规则生成**（不经 LLM，杜绝编造）；extra 原样并入，
    保证既有键不丢。"""
    out = {'status': status, 'summary': summary}
    if result is not None:
        out['result'] = result
    if interpret:
        out['interpretation'] = [x for x in interpret if x]
    if nxt:
        out['next'] = [x for x in nxt if x]
    out.update(extra)
    return out


_ERR_RULES = (
    ('未就绪', '对应库未构建——先在「数据库构建」页构建，或换选其他内置库'),
    ('缺少必填', '补齐必填参数（口径见 vp_list_tools / vp_cli_help）'),
    ('路径不存在', '确认输入路径存在（绝对路径或平台根相对路径均可）'),
    ('未检测到', '外部依赖缺失——到「设置 → 工具探测」检查'),
    ('超时', '加大超时或分批；长任务改用异步提交（vp_cli_submit / vp_run_tool）'),
    ('已存在', '样品重名——换名提交，或按提示续跑剩余阶段'),
    ('限流', '降低并发（jobs/gsa_jobs）后重试'),
)


def _interpret_error(err):
    """失败原因 → 可操作建议（关键词命中规则，命中几条给几条）。"""
    err = str(err or '')
    return [msg for kw, msg in _ERR_RULES if kw in err]


def _cli_spec(command):
    spec = COMMAND_SPECS.get(command)
    if not spec:
        raise ValueError('未知 CLI 命令: %s（可用: %s；analyze/tools/orfa/'
                         'report/host-analysis/tool-runs 的日常路径已由'
                         ' vp_run_tool / vp_pipeline_start / vp_list_tool_runs '
                         '覆盖）' % (command, ', '.join(sorted(COMMAND_SPECS))))
    return spec


def _cli_argv(command, args):
    """args 字典 → main.py argv（含必填校验与未知名拒绝）。"""
    spec = _cli_spec(command)
    a = dict(args or {})
    unknown = set(a) - {k for k, *_ in spec['args']} - \
        {p for p, _ in spec['positional']} - {'confirm'}
    if unknown:
        raise ValueError('未知参数: %s（%s 可用参数见 vp_cli_help）'
                         % (', '.join(sorted(unknown)), command))
    argv = [command]
    for pname, pdesc in spec['positional']:
        v = a.get(pname)
        if v in (None, ''):
            raise ValueError('缺少必填参数 %s（%s）' % (pname, pdesc))
        if command == 'datasets' and pname == 'raw':
            import shlex
            argv += shlex.split(str(v))
        else:
            argv.append(str(v))
    for k, t, req, d in spec['args']:
        v = a.get(k)
        if v is None or v is False or v == '':
            if req:
                raise ValueError('缺少必填参数 %s（%s）' % (k, d))
            continue
        flag = '--' + k.replace('_', '-')
        if t == 'bool':
            argv.append(flag)
        elif t == 'list':
            items = v if isinstance(v, (list, tuple)) else [v]
            for x in items:
                argv += [flag, str(x)]
        else:
            argv += [flag, str(v)]
    return argv


def _cli_display(argv):
    out = []
    redact = False
    for tok in argv:
        if redact:
            out.append('***')
            redact = False
            continue
        out.append(str(tok))
        if tok.replace('-', '').replace('_', '').lower() in \
                {s.lower() for s in CLI_SENSITIVE}:
            redact = True
    return ' '.join(out)


def _cli_env():
    env = dict(os.environ)
    env.setdefault('PYTHONIOENCODING', 'utf-8')
    return env


def _cli_check_confirm(command, spec, args):
    if not spec['confirm']:
        return
    a = args or {}
    if spec['confirm'] == 'apply':        # db-migrate：预检类免确认
        if a.get('check') or a.get('dry_run'):
            return
    elif spec['confirm'] == 'clean':      # tool-runs：仅 clean 落盘要确认
        if str(a.get('action') or '') != 'clean' or a.get('dry_run'):
            return
    if a.get('confirm') is not True:
        raise ValueError('%s 属于有破坏面的操作（%s）：请把 confirm 设为 '
                         'true 再执行' % (command, spec['desc'][:40]))


def _cli_run_impl(command, args, timeout_sec):
    spec = _cli_spec(command)
    if spec['mode'] != 'sync':
        raise ValueError('%s 是长任务命令（%s）：请用 vp_cli_submit 提交，'
                         '再用 vp_get_task 轮询' % (command, spec['desc']))
    _cli_check_confirm(command, spec, args)
    argv = [sys.executable, os.path.join(PLATFORM_ROOT, 'main.py')] + \
        _cli_argv(command, args)
    disp = _cli_display(argv)
    timeout = max(10, min(int(timeout_sec), 3600))
    import subprocess as _sp
    try:
        r = _sp.run(argv, cwd=PLATFORM_ROOT, capture_output=True, text=True,
                    encoding='utf-8', errors='replace', timeout=timeout,
                    stdin=_sp.DEVNULL,      # 不继承 MCP 协议管道（Windows 句柄继承会卡死子进程）
                    env=_cli_env())
        out = (r.stdout or '') + (r.stderr or '')
        rc = r.returncode
        timed_out = False
    except _sp.TimeoutExpired as e:
        out = ((e.stdout or b'').decode('utf-8', 'replace')
               if isinstance(e.stdout, bytes) else (e.stdout or '')) + \
              ((e.stderr or b'').decode('utf-8', 'replace')
               if isinstance(e.stderr, bytes) else (e.stderr or ''))
        rc, timed_out = None, True
    cap = 64 * 1024
    out = out[-cap:]
    ok = (rc == 0) and not timed_out
    return _fmt(
        'ok' if ok else 'failed',
        f"{command} {'完成' if ok else f'退出码 {rc}'}",
        result={'command': command, 'argv_display': disp, 'exit_code': rc,
                'timed_out': timed_out,
                'output': out if not timed_out else
                          out + f'\n[超时] 超过 {timeout}s 被终止；长耗时请'
                                f'改用 vp_cli_submit'},
        interpret=([] if ok else
                   ['非零退出/超时：看 output 尾部错误行；参数口径见 '
                    'vp_cli_help']),
        exit_code=rc, timed_out=timed_out,
        output=out if not timed_out else
               out + f'\n[超时] 超过 {timeout}s 被终止；长耗时请改用 '
                     f'vp_cli_submit')



def _cli_task_job(command, argv, disp):
    def job(log, prog, cancel):
        logger = TaskLogger(callback=log)
        log('$ ' + disp)
        prog('run', 0.3, f'{command} 运行中')
        import threading as _th
        import subprocess as _sp
        from Virus_Platform_Core.utils import (DEFAULT_CMD_TIMEOUT_SEC,
                                               task_register_proc)
        env = _cli_env()
        proc = _sp.Popen(argv, cwd=PLATFORM_ROOT, stdout=_sp.PIPE,
                         stderr=_sp.STDOUT, text=True, encoding='utf-8',
                         errors='replace', bufsize=1, env=env,
                         stdin=_sp.DEVNULL)   # 不继承 MCP 协议管道
        task_register_proc(proc)
        tail = []

        def _kill():
            _sp.run(['taskkill', '/F', '/T', '/PID', str(proc.pid)],
                    capture_output=True)
        try:
            for line in proc.stdout:
                line = line.rstrip('\n')
                if line:
                    log(line)
                    tail.append(line)
                    if len(tail) > 30:
                        tail.pop(0)
                if cancel is not None and cancel.is_set():
                    _kill()
                    proc.wait()
                    raise RuntimeError('用户取消')
        finally:
            rc = proc.wait()
        if rc != 0:
            raise RuntimeError(f'{command} 退出码 {rc}；末几行：\n'
                               + '\n'.join(tail[-8:]))
        prog('done', 1.0, f'{command} 完成')
        logger.close()
        return {'command': command, 'exit_code': 0,
                'output_tail': tail[-10:]}
    return job


def _cli_submit_impl(command, args, confirm):
    spec = _cli_spec(command)
    if spec['mode'] != 'task':
        raise ValueError('%s 是快查命令（%s）：请用 vp_cli_run 直接执行'
                         % (command, spec['desc']))
    a = dict(args or {})
    a['confirm'] = confirm
    _cli_check_confirm(command, spec, a)
    argv = [sys.executable, os.path.join(PLATFORM_ROOT, 'main.py')] + \
        _cli_argv(command, a)
    disp = _cli_display(argv)
    ts = time.strftime('%Y%m%d_%H%M%S')
    tid = tm.start(f'CLI·{command} {ts}', _cli_task_job(command, argv, disp),
                   weight=spec['weight'], out_dir=None)
    return _fmt(
        'submitted', f"CLI {command} 已提交",
        result={'task_id': tid, 'command': command, 'argv_display': disp},
        nxt=[f"vp_get_task(task_id='{tid}') 轮询",
             f"vp_get_task_log(task_id='{tid}', lines=200) 看实时输出"],
        hint='vp_get_task 轮询；长输出在任务日志（vp_get_task_log）')


@mcp.tool(name='vp_cli_help', annotations=_RO)
async def vp_cli_help() -> dict:
    """CLI 能力清单：37 个 main.py 子命令的用途、模式（sync=vp_cli_run /
    task=vp_cli_submit）、逐参数 schema（键/类型/必填/说明）。

    返回 {count, note, commands:{命令: {desc, mode, weight, positional,
    args:[{key,type,required,desc}], confirm}}}。
    """
    def _impl():
        out = {}
        for name in sorted(COMMAND_SPECS):
            s = COMMAND_SPECS[name]
            out[name] = {
                'desc': s['desc'], 'mode': s['mode'], 'weight': s['weight'],
                'positional': [{'key': p, 'desc': d} for p, d in
                               s['positional']],
                'args': [{'key': k, 'type': t, 'required': bool(r), 'desc': d}
                         for k, t, r, d in s['args']],
                'confirm': s['confirm'] or None}
        return {'count': len(out),
                'note': '参数键下划线书写自动转 --中划线；bool 传 true 即'
                        '加旗标；ncbi_api/deepseek_api 等键值在任务日志中'
                        '打码',
                'commands': out}
    return await anyio.to_thread.run_sync(_impl)


@mcp.tool(name='vp_cli_run', annotations=_RO)
async def vp_cli_run(
        command: Annotated[str, Field(min_length=1, max_length=40,
                                      description='vp_cli_help 里的 sync 命令')],
        args: Annotated[Optional[dict], Field(
            description='参数字典（键见 vp_cli_help；下划线自动转 --中划线）')] =
        None,
        timeout_sec: Annotated[int, Field(
            ge=10, le=3600, description='超时秒（默认 180；超时提示改用 '
                                        'vp_cli_submit）')] = 180) -> dict:
    """执行 CLI 快查命令（同步，秒~分钟级）：ncbi-list/gb-list/gb-check/
    ictv-status/ictv-refs/ref-status/submit-list/submit-init/submit-validate/
    submit-fill/submit-sbt/logan-jobs/selfcheck。

    返回 {ok, exit_code, output（尾部 ≤64KB）, argv_display}。
    长任务命令（下载/建库/提交批量）会拒收并提示改用 vp_cli_submit。
    """
    return await anyio.to_thread.run_sync(
        _cli_run_impl, command, args or {}, timeout_sec)


@mcp.tool(name='vp_cli_submit')
async def vp_cli_submit(
        command: Annotated[str, Field(min_length=1, max_length=40,
                                      description='vp_cli_help 里的 task 命令')],
        args: Annotated[Optional[dict], Field(
            description='参数字典（键见 vp_cli_help）')] = None,
        confirm: Annotated[bool, Field(
            description='db-migrate 落盘迁移 / tool-runs clean 等破坏面操作'
                        '必须 true')] = False) -> dict:
    """提交 CLI 长任务（异步）：数据库构建 / ncbi-dl / gb-dl / 公共数据检索
    (meta-search) / host-genome / ICTV 更新 / LOGAN 批量 / submit-export·adopt /
    db-migrate / tool-runs / samples-backfill 等全部耗命令。

    提交即返回 {task_id, command, argv_display, hint}；vp_get_task 轮询，
    实时日志 vp_get_task_log。敏感键（API key）在日志中打码。
    """
    return await anyio.to_thread.run_sync(
        _cli_submit_impl, command, args or {}, confirm)




# ======================================================================
# 病毒样本预筛（prescreen）：与 /api/prescreen/* 同源——提交（tm 流式
# 任务，日志实时落 out_dir/run.log）/ 项目列表 / 结果读取。
# ======================================================================

def _prescreen_submit_impl(a):
    from Virus_Platform_Core.web.prescreen import (_valid_name,
                                                   _meta_table_path,
                                                   _detect_run_col, _out_dir)
    name = str(a.get('name') or '').strip()
    if not _valid_name(name):
        raise ValueError('项目名只能含字母/数字/下划线/连字符')
    source = str(a.get('source') or 'paste').strip().lower()
    if source not in ('paste', 'project', 'file'):
        raise ValueError('无效的输入来源（paste|project|file）')
    out_dir = _out_dir(name)
    os.makedirs(out_dir, exist_ok=True)
    cmd_args = []
    if source == 'paste':
        runs = a.get('runs') or []
        if isinstance(runs, str):
            runs = runs.splitlines()
        runs = [str(x).strip() for x in runs if str(x).strip()]
        if not runs:
            raise ValueError('请至少给一条 Run 编号（runs）')
        if len(runs) > 20000:
            raise ValueError('单次最多 20000 条（更大批量请分批）')
        list_p = check_path(os.path.join(out_dir, 'input.list'),
                            must_exist=False, in_platform=True)
        with open(list_p, 'w', encoding='utf-8') as f:
            f.write('\n'.join(runs) + '\n')
        cmd_args += ['--runs', list_p]
    elif source == 'project':
        src_name = str(a.get('project') or '').strip()
        if not _valid_name(src_name):
            raise ValueError('无效的来源项目名')
        which = str(a.get('table') or '').strip()
        cands = ([which] if which in ('core14', 'search', 'full') else []) + \
            ['core14', 'search', 'full']
        csv_p = run_col = None
        for w in cands:
            mp = _meta_table_path(src_name, w)
            if mp and os.path.isfile(mp):
                col = _detect_run_col(mp)
                if col:
                    csv_p, run_col = mp, col
                    break
        if not csv_p:
            raise ValueError(f'项目 {src_name} 下没有可用的元数据表；'
                             '先跑公共数据检索，或改用 runs 粘贴')
        cmd_args += ['--csv', csv_p, '--run-col', run_col]
    else:
        rel = str(a.get('path') or '').strip()
        if not rel:
            raise ValueError('source=file 时须给平台内列表/CSV 路径 path')
        fp = check_path(rel, must_exist=True, in_platform=True)
        if not os.path.isfile(fp):
            raise ValueError('指定路径不是文件')
        ext = os.path.splitext(fp)[1].lower()
        if ext in ('.csv', '.tsv'):
            col = _detect_run_col(fp) or 'Run'
            cmd_args += ['--csv', fp, '--run-col', col]
        else:
            cmd_args += ['--runs', fp]

    vp = a.get('virus_pct')
    if vp not in (None, ''):
        try:
            vp = float(vp)
        except (TypeError, ValueError):
            raise ValueError('virus_pct 须为 0–100 的数字')
        if not (0 <= vp <= 100):
            raise ValueError('病毒占比阈值应在 0–100 之间')
        cmd_args += ['--virus-pct', str(vp)]

    def _int(k, default, lo, hi, label):
        try:
            v = int(a.get(k) or default)
        except (TypeError, ValueError):
            raise ValueError(f'{label} 须为 {lo}–{hi} 的整数')
        if not (lo <= v <= hi):
            raise ValueError(f'{label} 应在 {lo}–{hi} 之间（过高会被限流）')
        return v

    cmd_args += ['-j', str(_int('jobs', 6, 1, 20, 'NCBI 并发')),
                 '--gsa-jobs', str(_int('gsa_jobs', 3, 1, 6, 'GSA 并发'))]
    if a.get('no_taxdump'):
        cmd_args.append('--no-taxdump')
    if a.get('no_gsa'):
        cmd_args.append('--no-gsa')
    if a.get('recheck_nodata'):
        cmd_args.append('--recheck-nodata')

    from Virus_Platform_Core.config import engine_cmd
    from Virus_Platform_Core.public_meta import prescreen_engine
    cmd = engine_cmd(prescreen_engine.__file__, '-o', out_dir, *cmd_args)
    disp = ' '.join(str(x) for x in cmd)
    ts = time.strftime('%Y%m%d_%H%M%S')
    tid = tm.start(f'病毒预筛 {name} (MCP {ts})',
                   _cli_task_job('prescreen', cmd, disp),
                   log_file=os.path.join(out_dir, 'run.log'),
                   weight='light', out_dir=str(out_dir), link='/prescreen')
    return _fmt(
        'submitted', f"预筛项目 {name} 已提交",
        result={'task_id': tid, 'name': name, 'out_dir': str(out_dir)},
        nxt=[f"vp_get_task(task_id='{tid}') 轮询",
             f"完成后 vp_prescreen_result(name='{name}') 读 PASS 统计"],
        hint='vp_get_task 轮询；完成后 vp_prescreen_result 读结果')


def _prescreen_projects_impl():
    from Virus_Platform_Core.web.prescreen import _proj_dir
    root = DIRS['meta_search']
    try:
        names = sorted(os.listdir(root))
    except OSError:
        names = []
    out = []
    for n in names:
        if n.startswith('.') or n.startswith('_'):
            continue
        out.append({'project': n,
                    'prescreen_done': os.path.isfile(os.path.join(
                        _proj_dir(n), 'prescreen', 'stat_run_summary.tsv'))})
    return {'count': len(out), 'projects': out}


def _prescreen_result_impl(name):
    from Virus_Platform_Core.web.prescreen import (_valid_name, _out_dir,
                                                   _read_tsv, _coerce,
                                                   _stats, ARTIFACTS)
    if not _valid_name(name):
        raise ValueError('无效的项目名')
    out_dir = _out_dir(name)
    sum_p = os.path.join(out_dir, 'stat_run_summary.tsv')
    if not os.path.isfile(sum_p):
        return {'exists': False, 'name': name,
                'hint': '尚无结果；先用 vp_prescreen_run 提交'}
    summary = [_coerce(r) for r in _read_tsv(sum_p)]
    passed = [_coerce(r) for r in
              _read_tsv(os.path.join(out_dir, 'prescreen_pass.tsv'))]
    unresolved = [_coerce(r) for r in
                  _read_tsv(os.path.join(out_dir, 'unresolved.tsv'))]
    files = [{'name': fn, 'size': os.path.getsize(os.path.join(out_dir, fn))}
             for fn in ARTIFACTS
             if os.path.isfile(os.path.join(out_dir, fn))]
    st = _stats(summary)
    interp = [f"共 {st['total']} 条 Run：成功分析 {st['ok']}，"
              f"PASS {st['pass']}，检出病毒 {st['with_virus']}，"
              f"未解析 {st['unresolved']}（NCBI {st['ncbi_ok']}/{st['ncbi_total']}"
              f" · GSA {st['gsa_ok']}/{st['gsa_total']}）。"]
    if st['pass']:
        interp.append(f"PASS 清单 {st['pass']} 条已就绪，可作下载/流水线输入。")
    if st['unresolved']:
        interp.append('未解析多为源站无数据或超时，可用 recheck_nodata 复查。')
    return _fmt(
        'ok', f"预筛 {name}：PASS {st['pass']}/{st['total']}",
        result={'exists': True, 'name': name, 'out_dir': out_dir,
                'stats': st,
                'summary_rows': len(summary), 'summary_head': summary[:50],
                'passed_rows': len(passed), 'passed_head': passed[:50],
                'unresolved_rows': len(unresolved),
                'unresolved_head': unresolved[:30],
                'files': files},
        interpret=interp,
        nxt=(["PASS 明细按 passed_head / files 里的 prescreen_pass.tsv 取用"]
             if st['pass'] else []),
        hint='全量 TSV 按 abs 路径直接读取（prescreen 不在 tool_runs）')


@mcp.tool(name='vp_prescreen_run')
async def vp_prescreen_run(
        name: Annotated[str, Field(min_length=1, max_length=40,
                                   description='预筛项目名（字母/数字/_/-）')],
        source: Annotated[str, Field(
            max_length=10, description='paste|project|file，默认 paste')] =
        'paste',
        runs: Annotated[Optional[list[str]], Field(
            description='source=paste：Run 编号列表（SRR/CRR…，≤20000 条）')] =
        None,
        project: Annotated[Optional[str], Field(
            max_length=80,
            description='source=project：公共数据检索项目名')] = None,
        table: Annotated[Optional[str], Field(
            max_length=10,
            description='project 模式的元数据表 core14|search|full')] = None,
        path: Annotated[Optional[str], Field(
            max_length=1000,
            description='source=file：平台内 Run 列表或 CSV 路径')] = None,
        virus_pct: Annotated[Optional[float], Field(
            ge=0, le=100, description='病毒占比阈值 %（默认不过滤）')] = None,
        jobs: Annotated[int, Field(ge=1, le=20,
                                   description='NCBI 并发，默认 6')] = 6,
        gsa_jobs: Annotated[int, Field(ge=1, le=6,
                                       description='GSA 并发，默认 3')] = 3,
        no_taxdump: Annotated[bool, Field(description='跳过 taxdump 下载')] = False,
        no_gsa: Annotated[bool, Field(description='跳过 GSA（只查 NCBI）')] = False,
        recheck_nodata: Annotated[bool, Field(description='重查历史无数据记录')] = False,
) -> dict:
    """提交病毒样本预筛（SRA+GSA 在线核查病毒占比 → PASS 清单 + 汇总报告）。

    三种输入二选一：runs（粘贴 Run 列表）/ project（复用公共数据检索的
    元数据表）/ path（平台内列表或 CSV）。提交即返回 {task_id, name,
    out_dir}；vp_get_task 轮询、完成后 vp_prescreen_result 读结果。
    """
    a = {k: v for k, v in dict(
        name=name, source=source, runs=runs, project=project, table=table,
        path=path, virus_pct=virus_pct, jobs=jobs, gsa_jobs=gsa_jobs,
        no_taxdump=no_taxdump, no_gsa=no_gsa,
        recheck_nodata=recheck_nodata).items() if v is not None}
    return await anyio.to_thread.run_sync(_prescreen_submit_impl, a)


@mcp.tool(name='vp_prescreen_projects', annotations=_RO)
async def vp_prescreen_projects() -> dict:
    """列出可复用的公共数据检索项目（含预筛完成标记）。

    返回 {count, projects:[{project, prescreen_done}]}；project 名可作
    vp_prescreen_run(source='project') 的输入。
    """
    return await anyio.to_thread.run_sync(_prescreen_projects_impl)


@mcp.tool(name='vp_prescreen_result', annotations=_RO)
async def vp_prescreen_result(
        name: Annotated[str, Field(min_length=1, max_length=40,
                                   description='预筛项目名')]) -> dict:
    """读预筛结果：PASS/未解析统计、汇总与清单头部（各 ≤50 行）、产物
    文件清单与磁盘路径。

    返回 {exists, stats:{total,ok,with_virus,pass,unresolved,ncbi_ok,gsa_ok},
    summary_head, passed_head, files, ...}；不存在时 exists=false。
    """
    return await anyio.to_thread.run_sync(_prescreen_result_impl, name)



# ======================================================================
# OpenVirome（Serratus 云端病毒组数据）：公开 API 直连，无需凭据。
# 价值：按病毒物种/组织/地理等过滤，直接查出 **palmprint 阳性** 的
# SRA Run（Serratus 已在云端跑过 RdRp palmprint 的阳性记录）——比按
# 物种元数据检索精准一个量级。查到的 Run 列表可喂平台下载/预筛/流水线。
# 参考：github.com/serratus-bio/open-virome-mcp（我们只取其公开 API 层；
# 其 psql/Neo4j/Azure 工具需私有凭据，不纳入）。
# ======================================================================

_OV_API = 'https://zrdbegawce.execute-api.us-east-1.amazonaws.com/prod'
_OV_ROUTES = ('/identifiers', '/counts', '/results', '/mwas')
_OV_FILTER_HINT = ('organism, tissue, disease, sex, geo_attribute, biome, '
                   'virus_family, stat_organism')


def _ov_post(route, data, timeout=120):
    import httpx
    if route not in _OV_ROUTES:
        raise ValueError(f'非法 OpenVirome 路由: {route}')
    headers = {'Content-Type': 'application/json',
               'Accept': 'application/json',
               'User-Agent': 'VirusPlatform-MCP/1.2 (local virome platform; '
                             'github.com/zhangwenda0518/VirusPlatform)'}
    r = None
    for attempt in (1, 2):          # 上游网关偶发 5xx，重试一次
        r = httpx.post(_OV_API + route, json=data, timeout=timeout,
                       headers=headers)
        if r.status_code < 500:
            break
        if attempt == 1:
            time.sleep(2)
    if r.status_code >= 500:
        raise RuntimeError(
            'OpenVirome 上游返回 ' + str(r.status_code)
            + '——该过滤维度/取值可能不被 /identifiers 支持'
              '（已验证可用: organism、tissue）；换维度或取值重试')
    try:
        return r.json()
    except ValueError:
        import base64
        return json.loads(base64.b64decode(r.text).decode('utf-8'))


@mcp.tool(name='vp_ov_search_runs')
async def vp_ov_search_runs(
        filter_type: Annotated[str, Field(
            min_length=2, max_length=40,
            description=f'过滤维度，常见: {_OV_FILTER_HINT}')],
        filter_value: Annotated[str, Field(
            min_length=1, max_length=200,
            description='过滤值，如 "Tobacco mosaic virus" / "rhizosphere"')],
        palmprint_only: Annotated[bool, Field(
            description='只返回 palmprint 阳性 Run（推荐，默认开）')] = True,
        max_runs: Annotated[int, Field(
            ge=1, le=2000, description='返回 Run 编号上限，默认 200')] = 200,
        timeout_sec: Annotated[int, Field(
            ge=10, le=600, description='API 超时秒，默认 120')] = 120,
) -> dict:
    """OpenVirome：按元数据过滤查 **palmprint 阳性** 的 SRA Run 列表。

    这是比按物种元数据检索更强的入口——Serratus 已在云端对千万级 SRA
    跑过 RdRp palmprint，命中的 Run 是**确认有 RNA 病毒信号**的样本。
    返回 {total, runs(截断), truncated}；Run 列表可喂网页「数据下载」
    批次、或 meta-info 生成元数据、或逐条进预筛/流水线。
    """
    def _impl():
        resp = _ov_post('/identifiers', {
            'filters': [{'filterType': filter_type,
                         'filterValue': filter_value,
                         'groupByKey': filter_type}],
            'palmprintOnly': bool(palmprint_only)}, timeout=timeout_sec)
        run = (resp.get('run') or {}) if isinstance(resp, dict) else {}
        total = run.get('totalCount', 0)
        ids = run.get('single') or []
        cap = max(1, min(int(max_runs), 2000))
        return _fmt(
            'ok', f"OpenVirome 检出 palmprint 阳性 Run {total} 条"
                  f"（维度 {filter_type}={filter_value}）",
            result={'total': total, 'returned': min(len(ids), cap),
                    'runs': ids[:cap],
                    'truncated': len(ids) > cap,
                    'filter': f'{filter_type}={filter_value}',
                    'palmprint_only': bool(palmprint_only)},
            interpret=[f"这些 Run 在 Serratus 云端已确认含 RNA 病毒信号"
                       f"（palmprint 阳性={'是' if palmprint_only else '否'}）。"]
                   if total else ['该维度下没有命中的 Run——换过滤值或关闭 '
                                  'palmprint_only 试试。'],
            nxt=['把 runs 列表粘到网页「数据下载」建批次下载',
                 '或逐条 vp_run_tool(tool="kvsuite") 前先跑预筛确认',
                 f"用 vp_ov_counts 按 organism/tissue 等维度看这批 Run 的"
                 f"元数据分布"])
    return await anyio.to_thread.run_sync(_impl)


@mcp.tool(name='vp_ov_counts')
async def vp_ov_counts(
        run_ids: Annotated[list[str], Field(
            min_items=1, max_items=5000,
            description='SRA Run 编号列表（如 vp_ov_search_runs 的结果）')],
        group_by: Annotated[str, Field(
            min_length=2, max_length=40,
            description=f'分组维度，常见: {_OV_FILTER_HINT} | sra')],
        max_rows: Annotated[int, Field(
            ge=1, le=500, description='返回分组行数上限，默认 50')] = 50,
        timeout_sec: Annotated[int, Field(ge=10, le=600)] = 120,
) -> dict:
    """OpenVirome：对一批 Run 做元数据分组统计（宿主物种/组织/疾病/
    地理/生态位/病毒科…的分布）。

    返回 {rows:[{key,count}...]}；用于研判这批样本的价值与多样性。
    """
    def _impl():
        resp = _ov_post('/counts', {
            'table': group_by if group_by != 'sra' else 'sra',
            'groupBy': group_by, 'idColumn': 'run',
            'ids': [str(x).strip() for x in run_ids if str(x).strip()],
            'sortByColumn': 'count', 'sortByDirection': 'desc',
            'palmprintOnly': False,
            'pageStart': 0, 'pageEnd': 1000}, timeout=timeout_sec)
        rows = resp if isinstance(resp, list) else \
            (resp.get('results') or resp.get('data') or [])
        clean = []
        for r in rows[:max_rows]:
            if isinstance(r, dict):
                k = r.get('key') or r.get(group_by) or r.get('name') or ''
                v = r.get('count', '')
                clean.append({'key': k, 'count': v})
        return _fmt(
            'ok', f"元数据分组（{group_by}）共 {len(rows)} 组，返回前 "
                  f"{len(clean)} 组",
            result={'group_by': group_by, 'rows': clean,
                    'rows_total': len(rows)},
            nxt=['结合分布挑高价值 Run：vp_ov_search_runs 精筛，'
                 '或直接进预筛/流水线'])
    return await anyio.to_thread.run_sync(_impl)


@mcp.tool(name='vp_ov_mwas')
async def vp_ov_mwas(
        run_ids: Annotated[list[str], Field(
            min_items=1, max_items=5000,
            description='SRA Run 编号列表')],
        virus_families: Annotated[Optional[list[str]], Field(
            description='按病毒科过滤（如 ["Virgaviridae"]）；缺省=全部')] = None,
        max_rows: Annotated[int, Field(ge=1, le=500, description='默认 50')] = 50,
        timeout_sec: Annotated[int, Field(ge=10, le=600)] = 120,
) -> dict:
    """OpenVirome：对一批 Run 做 MWAS（元数据 × 病毒丰度关联分析）。

    返回 {rows:[{bioproject,family,metadata_field,metadata_value,fold_change,
    p_value}...]}——哪些宿主/组织/地理因素与病毒丰度显著相关。
    """
    def _impl():
        resp = _ov_post('/mwas', {
            'idColumn': 'run',
            'ids': [str(x).strip() for x in run_ids if str(x).strip()],
            'virusFamilies': virus_families,
            'pageStart': 0, 'pageEnd': 1000}, timeout=timeout_sec)
        rows = resp if isinstance(resp, list) else \
            (resp.get('results') or resp.get('data') or [])
        return _fmt(
            'ok', f"MWAS 关联分析返回 {len(rows)} 行",
            result={'rows_total': len(rows), 'rows_head': rows[:max_rows]},
            interpret=(['行按 fold_change/p_value 研判：|fold_change| 越大、'
                        'p_value 越小，关联越强。'] if rows else
                       ['该批 Run 无显著关联结果。']))
    return await anyio.to_thread.run_sync(_impl)


# ======================================================================
# 启动 gate：deferred.json / VP_MCP_DISABLED_TOOLS 摘除工具（fail-closed）
# ======================================================================
# ======================================================================
# 启动 gate：deferred.json / VP_MCP_DISABLED_TOOLS 摘除工具（fail-closed）
# ======================================================================

def _apply_gate():
    """从注册表摘除下线工具；返回 (已摘除, 未知名列表)。"""
    known = {t.name for t in mcp._tool_manager.list_tools()}
    wanted = {t for t in _gated_tools() if t.startswith('vp_')}
    unknown = sorted(wanted - known)
    gated = wanted & known
    for t in gated:
        mcp._tool_manager._tools.pop(t, None)
    if gated:
        sys.stderr.write('[virusplatform_mcp] 已下线工具: '
                         + ', '.join(sorted(gated)) + '\n')
    if unknown:
        sys.stderr.write('[virusplatform_mcp] 忽略下线清单里的未知工具名: '
                         + ', '.join(unknown) + '\n')
    return gated, unknown


_GATED, _GATE_UNKNOWN = _apply_gate()


def _enabled_tool_names():
    return sorted(t.name for t in mcp._tool_manager.list_tools())


# ======================================================================
# 入口
# ======================================================================

def _run_http(port, token=None):
    """--http [port] [--token X]：streamable-HTTP 模式（默认 127.0.0.1:8765）。

    参考 awesome-mcp-fastapi 的核心形态：MCP 挂到 HTTP 服务上——
    任务存活与客户端解耦（客户端断开任务照跑、重连可查），
    局域网内多个客户端可共享同一服务端。stdio 模式仍是默认。
    给 --token 时所有请求须带 Authorization: Bearer <token>（401 拦截）。
    """
    try:
        import uvicorn  # noqa: F401
    except ImportError:
        import subprocess
        subprocess.run([sys.executable, '-m', 'pip', 'install', '--quiet',
                        'uvicorn', '-i',
                        'https://pypi.tuna.tsinghua.edu.cn/simple'])
    app_inner = mcp.streamable_http_app()

    if token:
        async def _bearer_shield(scope, receive, send):
            if scope.get('type') == 'http':
                headers = {k.lower(): v for k, v in scope.get('headers') or []}
                got = headers.get(b'authorization', b'').decode('latin-1')
                if got != f'Bearer {token}':
                    await send({'type': 'http.response.start', 'status': 401,
                                'headers': [(b'content-type',
                                             b'text/plain; charset=utf-8')]})
                    await send({'type': 'http.response.body',
                                'body': b'unauthorized'})
                    return
            await app_inner(scope, receive, send)
        app_to_run = _bearer_shield
    else:
        app_to_run = app_inner

    import uvicorn
    sys.stderr.write(f'[virusplatform_mcp] streamable-HTTP: '
                     f'http://127.0.0.1:{port}/mcp'
                     f'{"（Bearer 鉴权已开）" if token else ""}\n')
    uvicorn.run(app_to_run, host='127.0.0.1', port=int(port), log_level='warning')


def _selfcheck():
    """--check：不进服务，打印环境体检结论（stdout）。"""
    print('virusplatform_mcp 自检')
    print('  python   :', sys.executable)
    print('  platform :', PLATFORM_ROOT)
    import mcp as _m
    print('  mcp sdk  :', getattr(_m, '__version__', 'installed'))
    info = _platform_info_impl()
    print('  host_db  :', info.get('host_db_ready'))
    print('  virus_db :', info.get('virus_db_ready'))
    print('  toolboxes:', len(TOOL_REGISTRY), '张卡注册')
    print('  mcp工具  :', info['tools_count'], '个启用',
          f"（指纹 {info['fingerprint']}）")
    if info['tools_gated']:
        print('  已下线   :', ', '.join(info['tools_gated']))
    try:
        import uvicorn
        print('  http模式  : 可用（uvicorn', uvicorn.__version__ + '；'
              '启动加 --http [--token X]）')
    except ImportError:
        print('  http模式  : 需先 pip install uvicorn（--http 首跑自动装）')
    problems = _platform_check()
    if problems:
        print('  !! 问题:', '; '.join(problems))
        return 1
    print('  ready    : yes（可直接挂载 stdio 服务）')
    return 0


if __name__ == '__main__':
    _argv = sys.argv[1:]
    if '--check' in _argv:
        sys.exit(_selfcheck())
    if '--http' in _argv:
        _i = _argv.index('--http')
        _port = int(_argv[_i + 1]) if len(_argv) > _i + 1 \
            and _argv[_i + 1].isdigit() else 8765
        _token = _argv[_argv.index('--token') + 1] \
            if '--token' in _argv else None
        _run_http(_port, _token)
    else:
        mcp.run()          # stdio 传输（默认）
