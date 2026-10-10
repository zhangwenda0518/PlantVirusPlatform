# -*- coding: utf-8 -*-
"""AI 提供方：OpenAI 兼容 /v1/chat/completions 客户端（流式 + 工具调用）。

为什么不用 openai 官方包：requirements.txt 里它长期处于注释态（平台依赖
策略是"能不新增就不新增"），且各家兼容端点的差异只在 base_url 与鉴权头，
requests 足以覆盖。Ollama 本地模型暴露的也是同一套 /v1/chat/completions，
因此一个客户端通吃云端与离线两类提供方。

流式协议：上游返回 SSE（`data: {...}` 行，`[DONE]` 结束）。本模块把碎片
增量组装成两类事件交出去：
    {'type': 'delta',      'content': '...'}            正文增量
    {'type': 'tool_calls', 'tool_calls': [...]}          完整组装后的工具调用
    {'type': 'message',    'message': {...}}             非流式/末帧的完整消息

工具调用的流式组装：OpenAI 协议下 tool_calls 以碎片到达（首片带 id/name，
后续片只有 arguments 的增量），按 index 归并拼接。
"""
import json
import time

import requests

# (connect, read) 超时；read 是**块间**超时——流式生成再长也不会触发，
# 只防上游卡死。本地 Ollama 首次加载模型可能要几十秒，connect 放宽到 15s。
TIMEOUT = (15, 300)

# 提供方预设：选预设后前端自动回填 base_url 与推荐模型；custom 由用户全填。
PROVIDER_PRESETS = [
    {'id': 'deepseek', 'name': 'DeepSeek（深度求索）',
     'base_url': 'https://api.deepseek.com/v1', 'model': 'deepseek-chat',
     'key_url': 'https://platform.deepseek.com'},
    {'id': 'zhipu', 'name': '智谱 GLM',
     'base_url': 'https://open.bigmodel.cn/api/paas/v4', 'model': 'glm-4-flash',
     'key_url': 'https://open.bigmodel.cn'},
    {'id': 'dashscope', 'name': '阿里云百炼（通义 Qwen）',
     'base_url': 'https://dashscope.aliyuncs.com/compatible-mode/v1',
     'model': 'qwen-plus', 'key_url': 'https://bailian.console.aliyun.com'},
    {'id': 'moonshot', 'name': '月之暗面 Kimi',
     'base_url': 'https://api.moonshot.cn/v1', 'model': 'moonshot-v1-8k',
     'key_url': 'https://platform.moonshot.cn'},
    {'id': 'openai', 'name': 'OpenAI',
     'base_url': 'https://api.openai.com/v1', 'model': 'gpt-4o-mini',
     'key_url': 'https://platform.openai.com'},
    {'id': 'ollama', 'name': 'Ollama（本地，离线可用）',
     'base_url': 'http://127.0.0.1:11434/v1', 'model': 'qwen2.5:7b',
     'key_url': ''},
    {'id': 'custom', 'name': '自定义（OpenAI 兼容端点）',
     'base_url': '', 'model': '', 'key_url': ''},
]


class AIConfigError(Exception):
    """AI 配置不完整（未启用 / 缺 Key / 缺模型），需要用户去设置页补全。"""


class AIProviderError(Exception):
    """上游端点错误（网络 / 鉴权 / 参数），message 面向用户可读。"""

    def __init__(self, message, status=0):
        super().__init__(message)
        self.status = status


def ai_ready(ai_cfg):
    """配置是否可用（启用 + base_url + model；Ollama 允许无 Key）。"""
    return bool(ai_cfg and ai_cfg.get('enabled') and ai_cfg.get('base_url')
                and ai_cfg.get('model')
                and (ai_cfg.get('api_key') or 'ollama' == ai_cfg.get('provider')))


def _headers(ai_cfg):
    h = {'Content-Type': 'application/json'}
    if ai_cfg.get('api_key'):
        h['Authorization'] = f"Bearer {ai_cfg['api_key']}"
    return h


def _endpoint(ai_cfg, path):
    base = str(ai_cfg.get('base_url') or '').rstrip('/')
    return f'{base}{path}'


def _raise_provider_error(resp):
    """把非 2xx 响应翻译成可读错误（各家错误体结构不一，尽量剥出 message）。"""
    detail = ''
    try:
        body = resp.json()
        err = body.get('error')
        if isinstance(err, dict):
            detail = err.get('message') or json.dumps(err, ensure_ascii=False)
        elif err:
            detail = str(err)
        else:
            detail = json.dumps(body, ensure_ascii=False)[:300]
    except Exception:
        detail = (resp.text or '')[:300]
    hint = ''
    if resp.status_code in (401, 403):
        hint = '（API Key 无效或无权限）'
    elif resp.status_code == 404:
        hint = '（端点路径不存在，请检查 Base URL 是否以 /v1 结尾）'
    elif resp.status_code == 429:
        hint = '（额度不足或限流）'
    raise AIProviderError(f'上游返回 {resp.status_code}{hint}: {detail}',
                          status=resp.status_code)


def _assemble_tool_calls(frags):
    """流式工具调用碎片按 index 归并 → OpenAI 消息格式的 tool_calls 列表。"""
    out = {}
    for f in frags:
        i = f.get('index') or 0
        slot = out.setdefault(i, {'id': '', 'type': 'function',
                                  'function': {'name': '', 'arguments': ''}})
        if f.get('id'):
            slot['id'] = f['id']
        fn = f.get('function') or {}
        if fn.get('name'):
            slot['function']['name'] += fn['name']
        if fn.get('arguments'):
            slot['function']['arguments'] += fn['arguments']
    return [out[k] for k in sorted(out)]


def chat_stream(ai_cfg, messages, tools=None, temperature=None, max_tokens=None,
                timeout=None):
    """流式对话：yield 组装后的事件 dict（见模块 docstring）。

    timeout: (connect, read) 覆盖；后台播报/闸门评估应传紧超时
    （见 broadcast._AI_CALL_TIMEOUT），别让一次挂起的上游把流水线
    闸门拖住几分钟。

    网络与协议错误统一抛 AIProviderError / AIConfigError；
    requests 层异常（DNS/超时/拒绝连接）也归一为 AIProviderError。
    """
    if not ai_ready(ai_cfg):
        raise AIConfigError('AI 助手未配置完整（设置页需启用并填写 Base URL / '
                            'API Key / 模型）')
    payload = {
        'messages': messages,
        'model': ai_cfg['model'],
        'stream': True,
        'temperature': (ai_cfg.get('temperature') if temperature is None
                        else temperature),
        'max_tokens': ai_cfg.get('max_tokens') or 2048,
    }
    if tools:
        payload['tools'] = tools
        payload['tool_choice'] = 'auto'
    try:
        resp = requests.post(_endpoint(ai_cfg, '/chat/completions'),
                             headers=_headers(ai_cfg), json=payload,
                             stream=True, timeout=timeout or TIMEOUT)
    except requests.RequestException as e:
        raise AIProviderError(f'无法连接 AI 端点: {type(e).__name__}: {e}')
    with resp:
        if resp.status_code != 200:
            _raise_provider_error(resp)
        # 兼容个别网关对 stream:true 返回整体 JSON（非 SSE）的情况
        ctype = (resp.headers.get('Content-Type') or '').lower()
        if 'text/event-stream' not in ctype:
            try:
                data = resp.json()
            except ValueError:
                raise AIProviderError('上游响应既非 SSE 也非 JSON')
            msg = ((data.get('choices') or [{}])[0].get('message')) or {}
            if msg.get('tool_calls'):
                yield {'type': 'tool_calls', 'tool_calls': msg['tool_calls']}
            if msg.get('content'):
                yield {'type': 'delta', 'content': msg['content']}
            yield {'type': 'message', 'message': msg}
            return
        content_started = False
        frags = []
        finish = ''
        for raw in resp.iter_lines(decode_unicode=True):
            if not raw or not raw.startswith('data:'):
                continue
            body = raw[5:].strip()
            if body == '[DONE]':
                break
            try:
                chunk = json.loads(body)
            except ValueError:
                continue
            choice = (chunk.get('choices') or [{}])[0]
            delta = choice.get('delta') or {}
            if delta.get('content'):
                content_started = True
                yield {'type': 'delta', 'content': delta['content']}
            if delta.get('tool_calls'):
                frags.extend(delta['tool_calls'])
            if choice.get('finish_reason'):
                finish = choice['finish_reason']
        if frags:
            tc = _assemble_tool_calls(frags)
            yield {'type': 'tool_calls', 'tool_calls': tc}
            yield {'type': 'message',
                   'message': {'role': 'assistant', 'content': '',
                               'tool_calls': tc}}
        elif content_started or finish != 'tool_calls':
            # 把已流出的正文再以完整消息收尾（工具循环追加历史时需要）
            yield {'type': 'message',
                   'message': {'role': 'assistant', 'content': None}}


def test_connection(ai_cfg):
    """连通性测试：发一个 1 token 的最小请求，返回时延与回复。"""
    t0 = time.time()
    payload = {'model': ai_cfg.get('model') or '', 'stream': False,
               'max_tokens': 8,
               'messages': [{'role': 'user', 'content': '回复"OK"两个字母即可'}]}
    try:
        resp = requests.post(_endpoint(ai_cfg, '/chat/completions'),
                             headers=_headers(ai_cfg), json=payload,
                             timeout=(10, 60))
    except requests.RequestException as e:
        return {'ok': False, 'error': f'{type(e).__name__}: {e}',
                'latency_ms': int((time.time() - t0) * 1000)}
    if resp.status_code != 200:
        try:
            _raise_provider_error(resp)
        except AIProviderError as e:
            return {'ok': False, 'error': str(e),
                    'latency_ms': int((time.time() - t0) * 1000)}
    try:
        data = resp.json()
        reply = (((data.get('choices') or [{}])[0].get('message') or {})
                 .get('content') or '').strip()
    except ValueError:
        reply = ''
    return {'ok': True, 'reply': reply[:80],
            'latency_ms': int((time.time() - t0) * 1000)}
