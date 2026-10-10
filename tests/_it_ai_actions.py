# -*- coding: utf-8 -*-
"""AI 执行层（三级执行权）集成测试。

用法: python tests/_it_ai_actions.py
退出码: 0 = 全部通过；1 = 存在失败。

执行安全设计验证点：白名单拒绝未知动作/未知参数、样品存在性校验、
执行单生命周期（pending→approved→launched / rejected）、并发上限、
第二级自动重试预算（每样品每阶段 1 次）、第三级自动调参预算与回落。
真实 pipeline 不启动——执行器 monkeypatch 记账。
"""
import json
import os
import sys
import time

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


def _wait(fn, timeout=8.0, step=0.05):
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
    from Virus_Platform_Core.ai_agent import executors, tools
    from Virus_Platform_Core.web.state import cfg
    c = flask_app.test_client()
    orig_ai = dict(cfg.ai)

    base_ai = {'enabled': True, 'provider': 'custom', 'base_url': 'x',
               'api_key': 'k', 'model': 'm', 'temperature': 0.3,
               'max_tokens': 512, 'explain': True, 'broadcast': True,
               'gate': False, 'gate_timeout_min': 15,
               'exec_auto': False, 'exec_auto_tune': False,
               'broadcast_stages': ['kvsuite'], 'gate_stages': []}
    cfg.ai = dict(base_ai)

    print('== 1) 执行单创建与校验 ==')
    out = tools.dispatch('propose_action', {
        'action': 'retry_failed', 'sample': 'smoke_e2e', 'reason': 'r1'})
    d = json.loads(out)
    check('合法提案创建 pending', d.get('created') and d.get('status') == 'pending',
          out[:120])
    aid1 = d.get('id')

    out = tools.dispatch('propose_action', {
        'action': 'delete_everything', 'sample': 'smoke_e2e', 'reason': 'x'})
    check('注册表外动作拒绝', '"error"' in out, out[:120])
    out = tools.dispatch('propose_action', {
        'action': 'rerun_stage', 'sample': 'no_such_sample_xx',
        'stages': ['assembly'], 'reason': 'x'})
    check('不存在样品拒绝', '"error"' in out, out[:120])
    out = tools.dispatch('propose_action', {
        'action': 'rerun_stage', 'sample': 'smoke_e2e',
        'stages': ['assembly'],
        'params': {'bogus_key': 1, 'min_contig_len': 50}, 'reason': 'x'})
    check('白名单外参数拒绝', 'bogus_key' in out, out[:150])
    out = tools.dispatch('propose_action', {
        'action': 'rerun_stage', 'sample': 'smoke_e2e',
        'stages': ['assembly'], 'params': {'min_contig_len': 999999},
        'reason': 'x'})
    d = json.loads(out)
    check('TUNABLE 边界钳制', d.get('created')
          and d.get('id') != aid1, out[:150])

    print('== 2) 决策端点 ==')
    r = c.get('/api/ai/actions')
    d = r.get_json()
    check('actions 清单 200', r.status_code == 200 and len(d['pending']) >= 2)
    r = c.post('/api/ai/actions/resolve', json={})
    check('空体 → 400（审计安全）', r.status_code == 400)
    r = c.post('/api/ai/actions/resolve',
               json={'id': 'act_nope', 'decision': 'approve'})
    check('未知执行单 → 404', r.status_code == 404)
    r = c.post('/api/ai/actions/resolve',
               json={'id': aid1, 'decision': 'nuke'})
    check('非法 decision → 4xx', r.status_code in (400, 404))

    print('== 3) 批准执行（执行器打桩） ==')
    calls = []
    executors._launch_analysis = (
        lambda sample, stages, body: calls.append((sample, stages, body))
        or 'fake_tid')
    r = c.post('/api/ai/actions/resolve',
               json={'id': aid1, 'decision': 'approve'})
    check('批准 200', r.status_code == 200)
    d = _wait(lambda: executors.get_action(aid1))
    check('状态翻转为 launched', d and d['status'] == 'launched', str(d))
    check('执行器收到样品', bool(calls and calls[-1][0] == 'smoke_e2e'),
          str(calls))

    print('== 4) 第二级 自动重试（预算 1 次/样品·阶段） ==')
    cfg.ai = dict(base_ai)
    cfg.ai['exec_auto'] = True
    calls.clear()
    executors.after_failure('smoke_e2e', 'assembly', 'boom')
    ok = _wait(lambda: any(c[0] == 'smoke_e2e' for c in calls))
    check('失败 → 自动重试启动', bool(ok), str(calls))
    n_after_1 = len(calls)
    executors.after_failure('smoke_e2e', 'assembly', 'boom again')
    time.sleep(1.0)
    check('同阶段第二次失败 → 预算拦截',
          len(calls) == n_after_1, f'{n_after_1} -> {len(calls)}')

    print('== 5) 第三级 自动调参（默认关；开后有预算与参数钳制） ==')
    cfg.ai = dict(base_ai)                     # 重置：exec_auto/tune 全关
    calls.clear()
    executors.after_failure('smoke_e2e', 'kvsuite', 'boom')
    time.sleep(1.0)
    check('默认关 → 不动作', not calls, str(calls))
    cfg.ai['exec_auto_tune'] = True
    executors._ai_propose_params = (lambda s, st, err, base: {
        'params': {'min_contig_len': 500, 'assembly_mode': 'metaviral'},
        'reason': 'contig 偏短，放宽阈值并换 metaviral 模式'})
    executors.after_failure('smoke_e2e', 'kvsuite', 'boom')
    # tier3 的启动在后台线程：等最后一个条目带上提议参数（而不是被
    # 并发到达的旧条目抢先）
    ok = _wait(lambda: calls and calls[-1][2].get('min_contig_len') == 500,
               timeout=10)
    check('开启后自动调参执行', bool(ok), str(calls))
    if ok:
        body = calls[-1][2]
        check('提议参数进入执行体',
              body.get('min_contig_len') == 500
              and body.get('assembly_mode') == 'metaviral', str(body))
    n = len(calls)
    executors.after_failure('smoke_e2e', 'orf', 'boom')
    time.sleep(1.5)
    check('每样品调参预算 1 次', len(calls) == n, f'{n} -> {len(calls)}')

    print('== 6) 配置持久化 ==')
    r = c.post('/api/ai/config', json={'exec_auto': True,
                                       'exec_auto_tune': False})
    d = r.get_json()
    check('两级开关可保存', r.status_code == 200
          and d['config']['exec_auto'] is True
          and d['config']['exec_auto_tune'] is False, str(d)[:200])

    cfg.ai = orig_ai
    cfg.save()

    print()
    if FAIL:
        print(f'✗ {len(FAIL)} 项失败: {FAIL}')
        return 1
    print('✔ 执行层三级全部通过')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
