# -*- coding: utf-8 -*-
"""AI 助手 API（对话驱动智能代理的 HTTP 边界）。

路由：
    GET  /api/ai/config   读取 AI 配置（api_key 掩码）+ 提供方预设 + 就绪态
    POST /api/ai/config   保存 AI 配置到 platform.json（逐项校验；空体 400）
    POST /api/ai/test     连通性测试（可用请求体临时覆盖，不必先保存）
    GET  /api/ai/context  平台简况 + 可选样品 digest（对话坞开场上下文）
    GET  /api/ai/audit    AI 代理审计留痕（跨重启读盘：执行单状态变迁 / 事件流）
    POST /api/ai/chat     对话（SSE 流式；事件见 ai_agent/chat.py）

审计安全（tests/_audit_routes.py 口径）：两条 POST 在空请求体下一律 400
且无任何副作用（不写配置、不发起网络请求），故不进 DESTRUCTIVE_POST。
"""
import json
import os
import re
import time

from flask import (Blueprint, Response, abort, jsonify, request,
                   stream_with_context)

from Virus_Platform_Core.ai_agent import audit
from Virus_Platform_Core.ai_agent import broadcast
from Virus_Platform_Core.ai_agent import chat as ai_chat
from Virus_Platform_Core.ai_agent import digest
from Virus_Platform_Core.ai_agent.provider import (PROVIDER_PRESETS,
                                                   ai_ready, test_connection)
from Virus_Platform_Core.config import _ai_local_path
from Virus_Platform_Core.web.state import cfg

bp = Blueprint('ai_assist', __name__)

_SECRET = 'api_key'
_DAY_RE = re.compile(r'^\d{4}-\d{2}-\d{2}$')


def _masked(ai):
    out = dict(ai or {})
    key = out.get(_SECRET) or ''
    out['api_key_set'] = bool(key)
    # 只露尾 4 位便于核对是哪把 Key；本体绝不出后端
    out[_SECRET] = ('••••' + key[-4:]) if len(key) > 4 else ('••••' if key else '')
    return out




# ------------------------------------------------------------------
# 供应商档案（ZCode 式多端点并存）：列表/保存/切换/删除。Key 分离存
# run/ai_local.json provider_keys[name]；platform.json 永不落 Key。
# 「切换」= 把档案值灌进平铺 cfg.ai（base_url/model/温度/长度/Key）并
# 落盘——对话坞/执行单/播报全链路立即使用新端点，零重启。
# ------------------------------------------------------------------

def _pf_local_read():
    try:
        with open(_ai_local_path(), encoding='utf-8') as f:
            d = json.load(f)
        return d if isinstance(d, dict) else {}
    except (OSError, ValueError):
        return {}


def _pf_local_write(d):
    lp = _ai_local_path()
    os.makedirs(os.path.dirname(lp), exist_ok=True)
    with open(lp, 'w', encoding='utf-8') as f:
        json.dump(d, f, ensure_ascii=False, indent=1)


@bp.route('/api/ai/profiles', methods=['GET'])
def api_ai_profiles_get():
    keys = _pf_local_read().get('provider_keys') or {}
    active = str(cfg.ai.get('active_profile') or '')
    out = [{'name': p.get('name'), 'base_url': p.get('base_url'),
            'model': p.get('model'), 'temperature': p.get('temperature'),
            'max_tokens': p.get('max_tokens'),
            'key_set': bool(keys.get(p.get('name'))),
            'active': p.get('name') == active}
           for p in (cfg.ai.get('providers') or [])]
    return jsonify({'count': len(out), 'profiles': out, 'active': active})


@bp.route('/api/ai/profiles/save', methods=['POST'])
def api_ai_profiles_save():
    body = request.get_json(silent=True) or {}
    name = str(body.get('name') or '').strip()
    if not name or len(name) > 40 or '/' in name or chr(92) in name:
        abort(400, '档案名须为 1~40 字符，且不含路径分隔符')
    base_url = str(body.get('base_url') or '').strip()
    model = str(body.get('model') or '').strip()
    if not base_url or not model:
        abort(400, 'base_url 与 model 必填')
    try:
        temp = min(2.0, max(0.0, float(body.get('temperature', 0.3))))
    except (TypeError, ValueError):
        abort(400, 'temperature 须为 0~2 的数字')
    try:
        mtok = max(64, min(32768, int(body.get('max_tokens', 2048))))
    except (TypeError, ValueError):
        abort(400, 'max_tokens 须为 64~32768 的整数')
    provs = [p for p in (cfg.ai.get('providers') or [])
             if p.get('name') != name]
    provs.append({'name': name, 'base_url': base_url, 'model': model,
                  'temperature': temp, 'max_tokens': mtok})
    cfg.ai = dict(cfg.ai, providers=provs)
    key = str(body.get('api_key') or '').strip()
    if key and not key.startswith('••••'):     # 掩码回传 = 未修改，保留旧 Key
        local = _pf_local_read()
        pk = local.get('provider_keys') or {}
        pk[name] = key
        local['provider_keys'] = pk
        _pf_local_write(local)
    try:
        cfg.save()
    except OSError as e:
        abort(500, f'配置保存失败: {e}')
    return jsonify({'ok': True, 'name': name})


@bp.route('/api/ai/profiles/activate', methods=['POST'])
def api_ai_profiles_activate():
    body = request.get_json(silent=True) or {}
    name = str(body.get('name') or '').strip()
    prof = next((p for p in (cfg.ai.get('providers') or [])
                 if p.get('name') == name), None)
    if not prof:
        abort(400, f'档案不存在: {name}')
    key = _pf_local_read().get('provider_keys', {}).get(name, '')
    cfg.ai = dict(cfg.ai, provider='custom', base_url=prof['base_url'],
                  model=prof['model'],
                  temperature=prof.get('temperature', 0.3),
                  max_tokens=prof.get('max_tokens', 2048),
                  api_key=key, active_profile=name)
    try:
        cfg.save()
    except OSError as e:
        abort(500, f'配置保存失败: {e}')
    return jsonify({'ok': True, 'active': name, 'ready': ai_ready(cfg.ai),
                    'config': _masked(cfg.ai)})


@bp.route('/api/ai/profiles/delete', methods=['POST'])
def api_ai_profiles_delete():
    body = request.get_json(silent=True) or {}
    name = str(body.get('name') or '').strip()
    provs = [p for p in (cfg.ai.get('providers') or [])
             if p.get('name') != name]
    if len(provs) == len(cfg.ai.get('providers') or []):
        abort(400, f'档案不存在: {name}')
    cfg.ai = dict(cfg.ai, providers=provs,
                  active_profile=('' if cfg.ai.get('active_profile') == name
                                  else cfg.ai.get('active_profile')))
    local = _pf_local_read()
    pk = local.get('provider_keys') or {}
    if name in pk:
        pk.pop(name)
        local['provider_keys'] = pk
        _pf_local_write(local)
    try:
        cfg.save()
    except OSError as e:
        abort(500, f'配置保存失败: {e}')
    return jsonify({'ok': True})

@bp.route('/api/ai/config', methods=['GET'])
def api_ai_config_get():
    # stage_options：流水线全部阶段（key + 当前语言显示名），供设置页
    # 把「播报阶段 / 闸门阶段」渲染成勾选组而不是手填代码
    try:
        from Virus_Platform_Core.pipeline import (STAGE_ORDER, STAGE_NAMES,
                                                  STAGE_NAMES_EN)
        lang = cfg.lang
        stage_options = [{'key': k,
                          'name': (STAGE_NAMES_EN if lang == 'en'
                                   else STAGE_NAMES).get(k) or k}
                         for k in STAGE_ORDER]
    except Exception:
        stage_options = []
    return jsonify({'config': _masked(cfg.ai),
                    'presets': PROVIDER_PRESETS,
                    'stage_options': stage_options,
                    'ready': ai_ready(cfg.ai)})


def _validate_ai_body(body):
    """校验并归一 AI 配置请求体 → 可合并进 cfg.ai 的 dict（非法项 400）。"""
    if not isinstance(body, dict) or not body:
        abort(400, '请求体须为非空 JSON 对象')
    out = {}
    if 'enabled' in body:
        out['enabled'] = bool(body['enabled'])
    if 'provider' in body:
        p = str(body['provider'] or '').strip()
        if p not in {x['id'] for x in PROVIDER_PRESETS}:
            abort(400, f'未知提供方: {p}')
        out['provider'] = p
    for k in ('base_url', 'model'):
        if k in body:
            out[k] = str(body[k] or '').strip()
    if _SECRET in body:
        v = str(body[_SECRET] or '').strip()
        if v.startswith('••••'):
            # 前端把掩码原样传回 = 未修改，保留已存 Key 不覆盖
            pass
        else:
            out[_SECRET] = v
    if 'temperature' in body:
        try:
            out['temperature'] = min(2.0, max(0.0, float(body['temperature'])))
        except (TypeError, ValueError):
            abort(400, 'temperature 须为 0~2 的数字')
    if 'max_tokens' in body:
        try:
            out['max_tokens'] = max(64, min(32768, int(body['max_tokens'])))
        except (TypeError, ValueError):
            abort(400, 'max_tokens 须为 64~32768 的整数')
    # 主动驱动三层（讲解/解读/决策闸门）——见 ai_agent/broadcast.py
    for k in ('explain', 'broadcast', 'gate', 'exec_auto', 'exec_auto_tune'):
        if k in body:
            out[k] = bool(body[k])
    if 'gate_timeout_min' in body:
        try:
            out['gate_timeout_min'] = max(1, min(120,
                                                 int(body['gate_timeout_min'])))
        except (TypeError, ValueError):
            abort(400, 'gate_timeout_min 须为 1~120 的整数（分钟）')
    for k in ('broadcast_stages', 'gate_stages'):
        if k in body:
            v = body[k]
            if not isinstance(v, list) or not all(
                    isinstance(x, str) for x in v):
                abort(400, f'{k} 须为字符串数组')
            out[k] = [x.strip() for x in v if x.strip()][:15]
    return out


@bp.route('/api/ai/config', methods=['POST'])
def api_ai_config_set():
    body = request.get_json(silent=True)
    merged = dict(cfg.ai)
    merged.update(_validate_ai_body(body))
    cfg.ai = merged
    try:
        cfg.save()
    except OSError as e:
        abort(500, f'配置保存失败: {e}')
    return jsonify({'ok': True, 'config': _masked(cfg.ai),
                    'ready': ai_ready(cfg.ai)})


@bp.route('/api/ai/test', methods=['POST'])
def api_ai_test():
    """连通性测试。请求体给要测试的配置（前端"先测再存"），**必须非空**——
    空体一律 400：否则全路由巡检的空 POST 会拿用户已存配置发真实网络请求。"""
    body = request.get_json(silent=True) or {}
    if not body:
        abort(400, '需要请求体（要测试的 AI 配置；留空字段即测试已存值）')
    ai = dict(cfg.ai)
    ai.update(_validate_ai_body(body))
    if not ai_ready(ai):
        abort(400, '配置不完整：需启用并填写 Base URL / 模型'
                   '（除本地 Ollama 外还需 API Key）')
    return jsonify(test_connection(ai))


@bp.route('/api/ai/context', methods=['GET'])
def api_ai_context():
    """对话坞开场上下文：样品清单 + 可选的选中样品 digest + 任务概览。"""
    sample = (request.args.get('sample') or '').strip()
    try:
        ctx = {'samples': digest.samples_brief(),
               'tasks': digest.task_overview(8),
               'ready': ai_ready(cfg.ai)}
        if sample:
            d = digest.sample_digest(sample)
            if d is not None:
                ctx['sample'] = d
        return jsonify(ctx)
    except Exception as e:
        return jsonify({'samples': [], 'tasks': [], 'ready': ai_ready(cfg.ai),
                        'error': f'{type(e).__name__}: {e}'})


@bp.route('/api/ai/gate', methods=['POST'])
def api_ai_gate():
    """第三层决策闸门确认：{sample, stage, action: 'continue'|'stop'}。

    空参数 POST 一律 400 且无副作用（全路由巡检安全）。
    """
    body = request.get_json(silent=True) or {}
    sample = str(body.get('sample') or '').strip()
    stage = str(body.get('stage') or '').strip()
    action = str(body.get('action') or '').strip()
    if not sample or not stage or not action:
        abort(400, '需要 sample / stage / action 三项')
    if not broadcast.resolve_gate(sample, stage, action):
        abort(404, '该阶段没有待决的 AI 闸门（可能已超时或已处理）')
    return jsonify({'ok': True})


@bp.route('/api/ai/actions')
def api_ai_actions():
    """AI 执行单清单（有限 JSON，审计/轮询安全）。"""
    from Virus_Platform_Core.ai_agent import executors
    recs = executors.list_actions()
    return jsonify({'pending': [r for r in recs if r['status'] == 'pending'],
                    'recent': recs})


@bp.route('/api/ai/actions/resolve', methods=['POST'])
def api_ai_action_resolve():
    """执行单决策：{id, decision: 'approve'|'reject'}。

    批准后动作在后台线程经与界面按钮相同的通道执行（executors）。
    空参数 POST 一律 400 且无副作用（全路由巡检安全）。
    """
    from Virus_Platform_Core.ai_agent import executors
    body = request.get_json(silent=True) or {}
    aid = str(body.get('id') or '').strip()
    decision = str(body.get('decision') or '').strip()
    if not aid or not decision:
        abort(400, '需要 id 与 decision（approve/reject）')
    rec, err = executors.resolve(aid, decision)
    if err:
        abort(404, err) if '不存在' in err else abort(400, err)
    return jsonify({'ok': True, 'action': executors.get_action(aid)})

@bp.route('/api/ai/events')
def api_ai_events():
    """主动播报事件流（双层协议，按客户端类型自动适配）：

    - 带 `Accept: text/event-stream`（前端 EventSource）→ SSE 常驻流：
      先补发 Last-Event-ID 之后的事件，随后实时推送 + 心跳；
    - 普通 GET（巡检/轮询/测试）→ **有限** JSON 快照
      `{events: [...], pending_gates: [...]}`——审计脚本会等完整响应体，
      无限流会把它挂死，这条边界必须守住。?since=<id> 增量取。
    """
    since = request.args.get('since')
    if not since:
        # EventSource 自动重连带的是 Last-Event-ID 请求头（不是查询参数）；
        # 不认它的话每次重连都会从头回放整个环形缓冲（前端虽去重，但浪费）。
        since = request.headers.get('Last-Event-ID')
    try:
        since_id = int(since) if since else 0
    except ValueError:
        since_id = 0
    wants_sse = 'text/event-stream' in (request.headers.get('Accept') or '')
    if not wants_sse:
        return jsonify({'events': broadcast.events_since(since_id),
                        'pending_gates': [g for g in broadcast.pending_gate()
                                          if not g.get('action')]})

    def gen():
        last = since_id
        idle = 0
        started = time.time()
        # 立即冲一行 SSE 注释：werkzeug 在首个非空块产生前不发响应头，
        # 不冲的话浏览器 EventSource 会一直停在 CONNECTING（实测踩过）
        cutoff = broadcast.current_max_id()
        yield ': hello\n\n'
        while True:
            # 客户端断连在这里**检测不到**：浏览器半关后（CLOSE_WAIT）写心跳
            # 依旧"成功"，生成器会永生，每页浏览泄漏一个线程 + 一条套接字。
            # 所以给连接设寿命上限，到期主动收尾——EventSource 会自动带
            # Last-Event-ID 重连续传，前端无感（2026-10-11 连接泄漏修复）。
            if time.time() - started > 300:
                yield ': rotate\n\n'
                return
            evs = broadcast.events_since(last)
            for ev in evs:
                last = ev['id']
                ev = dict(ev)
                # 连接瞬间之前已发生的事件 = 历史回放（重连/刷新补发），
                # 客户端据此抑制"自动弹坞"这类强打扰动作
                if ev['id'] <= cutoff:
                    ev['backlog'] = True
                yield (f"id: {ev['id']}\n"
                       f"data: {json.dumps(ev, ensure_ascii=False)}\n\n")
            if evs:
                idle = 0
            else:
                idle += 1
                if idle >= 20:                     # ~20s 心跳注释行
                    idle = 0
                    yield ': ping\n\n'
            time.sleep(1.0)

    return Response(stream_with_context(gen()),
                    mimetype='text/event-stream',
                    headers={'Cache-Control': 'no-cache',
                             'X-Accel-Buffering': 'no'})


@bp.route('/api/ai/audit')
def api_ai_audit():
    """AI 代理审计留痕（**跨重启读盘**）。

    与 `/api/ai/events` 的分工：后者是进程内环形缓冲的实时流（id 单调、
    重启归零）；本接口读的是 `run/ai_audit/*.jsonl` 的长期留痕 ——
    "AI 提过什么动作、谁批的、闸门何时弹的、用户选了继续还是终止"。

    参数：kind=actions|events（默认 actions）、limit=1..1000（默认 200）、
    day=YYYY-MM-DD（省略 = 从最近有记录的日期往前取）。只读。
    """
    kind = (request.args.get('kind') or 'actions').strip()
    if kind not in audit.KINDS:
        return jsonify({'error': f'kind 只允许: {", ".join(audit.KINDS)}'}), 400
    try:
        limit = int(request.args.get('limit') or 200)
    except ValueError:
        return jsonify({'error': 'limit 需为整数'}), 400
    limit = max(1, min(1000, limit))
    day = (request.args.get('day') or '').strip() or None
    if day is not None and not _DAY_RE.match(day):
        return jsonify({'error': 'day 需为 YYYY-MM-DD'}), 400
    try:
        records = audit.read(kind, limit=limit, day=day)
    except Exception as e:                 # 审计盘损坏也不该 500
        return jsonify({'kind': kind, 'records': [], 'stats': {},
                        'error': f'{type(e).__name__}: {e}'})
    return jsonify({'kind': kind, 'records': records,
                    'stats': audit.stats(), 'kinds': list(audit.KINDS)})


@bp.route('/api/ai/chat', methods=['POST'])
def api_ai_chat():
    body = request.get_json(silent=True) or {}
    messages = body.get('messages')
    if not isinstance(messages, list) or not messages:
        abort(400, '请求体须含非空 messages 数组')
    sample = str(body.get('sample') or '').strip() or None

    def gen():
        try:
            for ev in ai_chat.run_chat(cfg.ai, messages, sample):
                yield ev
        except Exception as e:                 # 编排层兜底（chat.py 已兜一层）
            yield f"data: {json.dumps({'type': 'error', 'message': f'{type(e).__name__}: {e}'}, ensure_ascii=False)}\n\n"
            yield "data: " + json.dumps({'type': 'done'}) + "\n\n"

    return Response(stream_with_context(gen()),
                    mimetype='text/event-stream',
                    headers={'Cache-Control': 'no-cache',
                             'X-Accel-Buffering': 'no'})
