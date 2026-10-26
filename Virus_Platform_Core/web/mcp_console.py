# -*- coding: utf-8 -*-
"""MCP 控制台 API：网页端直接驱动 virusplatform_mcp 的全部工具。

后端以最小 stdio 客户端（JSON-RPC，换行分帧）拉起平台根的 mcp_server.py
—— 与外部客户端（Claude Desktop / ZCode）走完全相同的 MCP 协议；工具
执行与任务管理归 MCP 服务器进程（heavy/light 并发闸门与网页端独立计数）。

路由：
  GET  /api/mcp/status   → {running, pid, script}
  POST /api/mcp/restart  → 重拉服务器（长同步调用卡住时的手动出口）
  GET  /api/mcp/tools    → {count, tools:[{name, description, inputSchema}]}
  POST /api/mcp/call     → {name, arguments} → {isError, text, data}

注：工具全部「提交/快查即返回」，正常调用秒级；vp_cli_run 若把
timeout_sec 调到 300s 以上，本客户端（330s）会先行超时——超时不杀
服务器（可能有任务在跑），用 /api/mcp/restart 手动重置。
工具总数不写死：以 mcp_server.py 的 @mcp.tool 注册数为准（/api/mcp/tools
返回 count），避免文档与实际能力漂移。
"""
import atexit
import json
import os
import subprocess
import sys
import threading
import time

from flask import Blueprint, abort, jsonify, request

from Virus_Platform_Core.config import PLATFORM_ROOT

bp = Blueprint('mcpconsole', __name__)

_LOCK = threading.RLock()
_PROC = None
_INITED = False
_REQ_ID = 0
_SERVER_SCRIPT = os.path.join(PLATFORM_ROOT, 'mcp_server.py')
_BOOT_TIMEOUT = 150          # 首次拉起含平台模块导入
_CALL_TIMEOUT = 330          # 工具全部「提交即返回」；覆盖 vp_cli_run 默认 180s

# 破坏面工具：网页控制台不开放（与 ai_agent/tools.py 的 _MCP_BLOCKED 同口径）。
# 不写死总数：工具清单由 mcp_server.py 的 @mcp.tool 决定，任何一处写死数字
# 都会随增删漂移（曾出现系统提示写 19、实际 22 的口径不一致）。
_BLOCKED_TOOLS = {'vp_delete_tool_run'}


def _alive():
    return _PROC is not None and _PROC.poll() is None


def _kill():
    global _PROC, _INITED
    if _PROC is not None:
        try:
            _PROC.kill()
        except OSError:
            pass
    _PROC = None
    _INITED = False


atexit.register(_kill)


def _read_msg(proc, want_id, deadline):
    """逐行读 JSON；通知/无主消息跳过，直到拿到 want_id 的响应。"""
    while True:
        if time.monotonic() > deadline:
            raise TimeoutError('MCP 服务器响应超时')
        line = proc.stdout.readline()
        if not line:
            raise ConnectionError(
                f'MCP 服务器已退出（退出码 {proc.poll()}）')
        line = line.strip()
        if not line:
            continue
        try:
            msg = json.loads(line)
        except ValueError:
            continue
        if isinstance(msg, dict) and msg.get('id') == want_id:
            return msg


def _ensure():
    """拉起服务器并完成 initialize 握手。"""
    global _PROC, _INITED, _REQ_ID
    if _alive() and _INITED:
        return
    _kill()
    env = dict(os.environ)
    env.setdefault('PYTHONIOENCODING', 'utf-8')
    if getattr(sys, 'frozen', False):
        # 冻结分发没有独立 python，不能 `exe mcp_server.py`（exe 会把参数
        # 当网页参数）——用 app.py 顶部的 --mcp-stdio 入口在本 exe 进程内
        # 跑 mcp_server（模块与 mcp SDK 已随 PyInstaller 入包，2026-10-24）。
        cmd = [sys.executable, '--mcp-stdio']
    else:
        cmd = [sys.executable, _SERVER_SCRIPT]
    _PROC = subprocess.Popen(
        cmd, cwd=PLATFORM_ROOT,
        stdin=subprocess.PIPE, stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL, text=True, encoding='utf-8',
        errors='replace', bufsize=1, env=env)
    deadline = time.monotonic() + _BOOT_TIMEOUT
    _REQ_ID += 1
    rid = _REQ_ID
    _PROC.stdin.write(json.dumps({
        'jsonrpc': '2.0', 'id': rid, 'method': 'initialize',
        'params': {'protocolVersion': '2025-03-26', 'capabilities': {},
                   'clientInfo': {'name': 'virusplatform-web', 'version': '1.0'}}
    }) + '\n')
    _PROC.stdin.flush()
    resp = _read_msg(_PROC, rid, deadline)
    if resp.get('error'):
        raise RuntimeError('initialize 失败: '
                           + json.dumps(resp['error'], ensure_ascii=False))
    _PROC.stdin.write(json.dumps(
        {'jsonrpc': '2.0', 'method': 'notifications/initialized'}) + '\n')
    _PROC.stdin.flush()
    _INITED = True


def _request(method, params=None, timeout=_CALL_TIMEOUT):
    """发一次请求并等响应；进程半死（连接类错误）自动重拉并重试一次。

    超时**不**杀服务器——可能正有合法长调用在执行；用 /api/mcp/restart
    手动重置。"""
    global _REQ_ID
    last = None
    for attempt in (1, 2):
        try:
            with _LOCK:
                _ensure()
                deadline = time.monotonic() + timeout
                _REQ_ID += 1
                rid = _REQ_ID
                _PROC.stdin.write(json.dumps(
                    {'jsonrpc': '2.0', 'id': rid, 'method': method,
                     'params': params or {}}) + '\n')
                _PROC.stdin.flush()
                resp = _read_msg(_PROC, rid, deadline)
                if resp.get('error'):
                    err = resp['error'] or {}
                    raise RuntimeError(
                        str(err.get('message')
                            or json.dumps(err, ensure_ascii=False)))
                return resp.get('result') or {}
        except (ConnectionError, BrokenPipeError, OSError) as e:
            last = e
            _kill()               # 连接坏了才重拉；超时保留现场
            if attempt == 2:
                break
        except TimeoutError as e:
            raise RuntimeError(
                f'{e}（服务器可能仍在执行该调用——稍候重试，或 POST '
                f'/api/mcp/restart 重置）') from e
    raise RuntimeError(f'MCP 服务器连接失败（已尝试重启一次）: {last}')


@bp.route('/api/mcp/status')
def api_mcp_status():
    with _LOCK:
        return jsonify({'running': bool(_alive() and _INITED),
                        'pid': _PROC.pid if _PROC else None,
                        'script': _SERVER_SCRIPT})


@bp.route('/api/mcp/restart', methods=['POST'])
def api_mcp_restart():
    with _LOCK:
        _kill()
    _request('tools/list', timeout=_BOOT_TIMEOUT)
    return jsonify({'ok': True})


@bp.route('/api/mcp/tools')
def api_mcp_tools():
    res = _request('tools/list')
    tools = [{'name': t.get('name'),
              'description': t.get('description') or '',
              'inputSchema': t.get('inputSchema') or {}}
             for t in (res.get('tools') or [])]
    return jsonify({'count': len(tools), 'tools': tools})


@bp.route('/api/mcp/call', methods=['POST'])
def api_mcp_call():
    body = request.get_json(force=True) or {}
    name = str(body.get('name') or '').strip()
    if not name:
        abort(400, '缺少工具名 name')
    # 破坏面工具与 AI 对话坞同口径：一律不开放。
    # vp_delete_tool_run 会 rm 整个运行目录（不可逆），此前网页控制台是
    # 裸透传 —— 普通用户点开就能删结果。删除请走任务中心/运行历史的界面。
    if name in _BLOCKED_TOOLS:
        abort(403, f'{name} 属破坏面工具：网页控制台不开放。'
                   '删除运行请到「任务中心 / 运行历史」操作')
    args = body.get('arguments')
    if not isinstance(args, dict):
        args = {}
    res = _request('tools/call', {'name': name, 'arguments': args})
    text = '\n'.join(
        c.get('text', '') for c in (res.get('content') or [])
        if isinstance(c, dict) and c.get('type') == 'text')
    data = res.get('structuredContent')
    if data is None:
        try:
            data = json.loads(text)
        except ValueError:
            data = None
    return jsonify({'isError': bool(res.get('isError')), 'text': text,
                    'data': data})
