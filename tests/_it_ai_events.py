# -*- coding: utf-8 -*-
"""AI 主动驱动三层（讲解/解读/决策闸门）集成测试。

不依赖真实 AI：起线程内模拟 OpenAI 兼容端点（纯内容流式），闸门评估直接
打桩。覆盖：事件总线与双协议端点、三层开关行为、闸门 继续/终止/超时
三路、gate API 校验与审计安全。

用法: python tests/_it_ai_events.py
退出码: 0 = 全部通过；1 = 存在失败。
"""
import json
import os
import sys
import threading
import time
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

for _s in (sys.stdout, sys.stderr):
    try:
        _s.reconfigure(errors='replace')
    except (AttributeError, OSError):
        pass

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
os.environ.setdefault('VP_NO_RECOVER', '1')
# 本守卫会发布大量事件 / 建执行单；关掉审计落盘，避免往 run/ai_audit/
# 写测试垃圾（审计本身的正确性由 tests/_check_ai_audit.py 在临时目录专测）。
os.environ.setdefault('VP_AI_AUDIT', '0')

FAIL = []


def check(name, cond, detail=''):
    print(('  ✔ ' if cond else '  ✘ ') + name + (f' — {detail}' if detail and not cond else ''))
    if not cond:
        FAIL.append(name)


class _MockAI(BaseHTTPRequestHandler):
    """纯内容流式端点（供第二层播报/诊断消费）。"""
    def log_message(self, *a):
        pass

    def do_POST(self):
        n = int(self.headers.get('Content-Length') or 0)
        try:
            json.loads(self.rfile.read(n) or b'{}')
        except ValueError:
            pass
        self.send_response(200)
        self.send_header('Content-Type', 'text/event-stream')
        self.end_headers()
        for piece in ('模拟播报：', '检出 3 种病毒，质量正常。'):
            chunk = {'choices': [{'delta': {'content': piece}}]}
            self.wfile.write(f"data: {json.dumps(chunk)}\n\n".encode())
        self.wfile.write(b'data: [DONE]\n\n')
        self.wfile.flush()


def _wait(fn, timeout=8.0, step=0.05):
    """轮询直到 fn 为真（返回值真值）或超时；返回 fn 的最终返回值。"""
    deadline = time.time() + timeout
    out = None
    while time.time() < deadline:
        out = fn()
        if out:
            return out
        time.sleep(step)
    return out


def main():
    from app import app as flask_app
    from Virus_Platform_Core.ai_agent import broadcast
    from Virus_Platform_Core.web.state import cfg
    c = flask_app.test_client()
    orig_ai = dict(cfg.ai)

    base_ai = {'enabled': True, 'provider': 'custom', 'base_url': 'x',
               'api_key': 'k', 'model': 'm', 'temperature': 0.3,
               'max_tokens': 512, 'explain': True, 'broadcast': True,
               'gate': False, 'gate_timeout_min': 15,
               'broadcast_stages': ['kvsuite', 'assembly', 'report'],
               'gate_stages': ['kvsuite', 'assembly']}

    print('== 0) 事件总线与双协议端点 ==')
    ev = broadcast.publish('_test', hello=1)
    r = c.get('/api/ai/events?since=0')
    d = r.get_json()
    check('JSON 快照 200', r.status_code == 200)
    check('发布事件可读', any(e['id'] == ev['id'] for e in d['events']))
    r = c.get(f'/api/ai/events?since={ev["id"]}')
    check('since 增量过滤', all(e['id'] > ev['id'] for e in r.get_json()['events']))
    r = c.get('/api/ai/events', headers={'Accept': 'text/event-stream'})
    check('SSE 模式 200 + text/event-stream',
          r.status_code == 200 and 'text/event-stream'
          in (r.headers.get('Content-Type') or ''))
    try:
        first = next(r.response)
        # 首块是连接即冲的 ': hello' 注释（冲响应头用），或是事件/心跳
        check('SSE 流可读', b'hello' in first or b'ping' in first
              or b'data:' in first, str(first[:60]))
    finally:
        r.response.close()

    print('== 1) 第一层 讲解 ==')
    cfg.ai = dict(base_ai)
    n0 = broadcast.events_since(0)[-1]['id']
    broadcast.note_stage_start('smoke_e2e', 'kvsuite')
    evs = _wait(lambda: [e for e in broadcast.events_since(n0)
                         if e['type'] == 'stage_start'])
    check('stage_start 事件发布', bool(evs))
    check('讲解文案非空且含阶段名',
          bool(evs and evs[0].get('note') and evs[0].get('name')))
    cfg.ai['explain'] = False
    m = broadcast.events_since(0)[-1]['id']
    broadcast.note_stage_start('smoke_e2e', 'assembly')
    time.sleep(0.3)
    check('开关关闭 → 不发布',
          all(e['type'] != 'stage_start' for e in broadcast.events_since(m)))

    print('== 2) 第二层 解读 ==')
    srv = ThreadingHTTPServer(('127.0.0.1', 0), _MockAI)
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    cfg.ai = dict(base_ai)
    cfg.ai['base_url'] = f'http://127.0.0.1:{srv.server_address[1]}/v1'
    m = broadcast.events_since(0)[-1]['id']
    broadcast.note_stage_done('smoke_e2e', 'primer')       # 非播报阶段
    time.sleep(0.6)
    evs = broadcast.events_since(m)
    check('非关键阶段只有轻量完成条',
          any(e['type'] == 'stage_done' for e in evs)
          and not any(e['type'] == 'ai_read' for e in evs), str([e['type'] for e in evs]))
    m = broadcast.events_since(0)[-1]['id']
    broadcast.note_stage_done('smoke_e2e', 'kvsuite')      # 关键阶段 → AI 播报
    evs = _wait(lambda: [e for e in broadcast.events_since(m)
                         if e['type'] == 'ai_read'], timeout=10)
    check('关键阶段触发 ai_read', bool(evs))
    check('播报文本来自模拟端点',
          bool(evs and '模拟播报' in (evs[0].get('text') or '')),
          str(evs[:1]))
    cfg.ai['broadcast'] = False
    m = broadcast.events_since(0)[-1]['id']
    broadcast.note_stage_failed('smoke_e2e', 'assembly', RuntimeError('boom'))
    evs = _wait(lambda: [e for e in broadcast.events_since(m)
                         if e['type'] == 'stage_failed'])
    check('失败事件必发（与开关无关）', bool(evs))
    time.sleep(0.6)
    check('播报关 → 失败无 AI 诊断',
          all(e['type'] != 'ai_read' for e in broadcast.events_since(m)))

    print('== 3) 第三层 决策闸门 ==')
    cfg.ai = dict(base_ai)
    cfg.ai['gate'] = False
    check('闸门关 → 直接继续', broadcast.gate_point('smoke_e2e', 'kvsuite')
          == 'continue')
    cfg.ai['gate'] = True
    cfg.ai['gate_stages'] = ['kvsuite']
    broadcast._ai_assess = (lambda s, st: '建议继续：检出 3 种病毒，'
                            'reads 占比正常。')      # 打桩评估
    # —— 继续路 ——
    out = {}
    th = threading.Thread(target=lambda: out.update(
        r=broadcast.gate_point('smoke_e2e', 'kvsuite')), daemon=True)
    th.start()
    pend = _wait(lambda: broadcast.pending_gate('smoke_e2e', 'kvsuite'))
    check('闸门进入待决', bool(pend))
    gate_ev = _wait(lambda: [e for e in broadcast.events_since(0)
                             if e['type'] == 'ai_gate'
                             and e.get('status') == 'pending'])
    check('pending 闸门事件发布', bool(gate_ev))
    r = c.post('/api/ai/gate', json={'sample': 'smoke_e2e',
                                     'stage': 'kvsuite', 'action': 'continue'})
    check('确认 API 200', r.status_code == 200)
    th.join(timeout=5)
    check('闸门放行 → 线程返回 continue', out.get('r') == 'continue', str(out))
    resolved = _wait(lambda: [e for e in broadcast.events_since(0)
                              if e['type'] == 'ai_gate'
                              and e.get('status') == 'continue'])
    check('继续状态事件发布', bool(resolved))
    # —— 终止路 ——
    out = {}
    th = threading.Thread(target=lambda: out.update(
        err=_capture(lambda: broadcast.gate_point('smoke_e2e', 'kvsuite'))),
        daemon=True)
    th.start()
    _wait(lambda: broadcast.pending_gate('smoke_e2e', 'kvsuite'))
    r = c.post('/api/ai/gate', json={'sample': 'smoke_e2e',
                                     'stage': 'kvsuite', 'action': 'stop'})
    check('终止 API 200', r.status_code == 200)
    th.join(timeout=5)
    check('终止 → GateAborted', isinstance(out.get('err'),
                                           broadcast.GateAborted), str(out))
    # —— 超时路 ——
    out = {}
    th = threading.Thread(target=lambda: out.update(
        r=broadcast.gate_point('smoke_e2e', 'kvsuite', timeout_s=1)),
        daemon=True)
    th.start()
    th.join(timeout=5)
    check('超时自动继续', out.get('r') == 'timeout', str(out))
    # —— API 校验（审计安全）——
    check('空体 gate → 400',
          c.post('/api/ai/gate', json={}).status_code == 400)
    check('无待决闸门 → 404',
          c.post('/api/ai/gate', json={'sample': 'x', 'stage': 'y',
                                       'action': 'continue'}).status_code == 404)
    check('非法 action → 404',
          c.post('/api/ai/gate', json={'sample': 'smoke_e2e',
                                       'stage': 'kvsuite',
                                       'action': 'nuke'}).status_code == 404)

    print('== 4) 配置校验 ==')
    r = c.post('/api/ai/config', json={'gate': True, 'gate_timeout_min': 0,
                                       'broadcast_stages': ['kvsuite', 'bogus']})
    d = r.get_json()
    check('三层开关可保存', r.status_code == 200 and d['config']['gate'] is True)
    check('超时钳到下限 1', d['config']['gate_timeout_min'] == 1)
    check('阶段清单按字符串保存', 'kvsuite' in d['config']['broadcast_stages']
          and 'bogus' in d['config']['broadcast_stages'])
    check('非字符串数组 → 400',
          c.post('/api/ai/config', json={'broadcast_stages': [3]}).status_code == 400)
    r = c.post('/api/ai/config', json={'gate_timeout_min': 'abc'})
    check('非法超时 → 400', r.status_code == 400)

    srv.shutdown()
    cfg.ai = orig_ai
    cfg.save()

    print()
    if FAIL:
        print(f'✗ {len(FAIL)} 项失败: {FAIL}')
        return 1
    print('✔ 主动驱动三层全部通过')
    return 0


def _capture(fn):
    try:
        fn()
        return None
    except Exception as e:
        return e


if __name__ == '__main__':
    raise SystemExit(main())
