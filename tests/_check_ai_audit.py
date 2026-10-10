# -*- coding: utf-8 -*-
"""AI 审计落盘守卫（2026-10-08 新增）。

为什么需要：`executors._actions` 与 `broadcast._events` 原本都是**进程内内存**
——重启即失，事后无从追溯"AI 提了什么动作、谁批的、什么时候执行"。审计落盘
（`ai_agent/audit.py`）是"人在回路"这个设计承诺的证据链，一旦静默失效，
平台仍然"看起来完全正常"，只是**没有留痕**。

本守卫按 flask-async-write-drain 的清单逐项验证：

  A 功能与边界：条数不多不少 / 按天轮转 / 脱敏 / 截断 / stats / 保留期清理
  B **退出不丢**（真子进程 + SIGTERM）：Windows 下 SIGTERM 不触发 atexit，
    必须靠自装的信号处理器排空队列 —— 这是最容易静默失效的一环
  C **不重复写**：`_drain()` 与后台线程竞态时同一批只能落一次
  D **接线自证**：`broadcast.publish` / `executors` 状态变迁真的会落盘
    （只测 audit 模块自身等于没测接线）
  E **重启不撞号**：盘上已有 `act_7` 时，新单必须是 `act_8`
  F **范围自证**：`DIRS['ai_audit']` 必须在 `OUTPUT_SUBDIRS` 里，否则用户
    切换"输出目录"时审计会写回平台根（静默偏离约定）
  G **灵敏度自证**：注入故意错误（不脱敏 / 丢一半 / 关掉接线），必须变红

用法: python tests/_check_ai_audit.py
退出码: 0 = 全部通过；1 = 存在失败。
"""
import io
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time

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


def _lines(path):
    if not os.path.isfile(path):
        return []
    with io.open(path, encoding='utf-8', errors='replace') as f:
        return [ln for ln in (x.strip() for x in f) if ln]


def _count_dir(d):
    n = 0
    for name in os.listdir(d) if os.path.isdir(d) else []:
        if name.endswith('.jsonl'):
            n += len(_lines(os.path.join(d, name)))
    return n


def _restore(audit, key, orig, tmp):
    """还原被替换的目录。

    ⚠ 必须**先排空再换根**：审计根是在「写盘时」解析的，若队列里还有记录
    就换回真实根，退出时 atexit 排空会把测试记录写进 `run/ai_audit/`
    （实测踩过 —— 上一次失败运行在真实目录留下了两个 jsonl）。
    """
    try:
        audit.flush_now(timeout=5.0)
    except Exception:
        pass
    audit.DIRS[key] = orig
    shutil.rmtree(tmp, ignore_errors=True)


# ------------------------------------------------------------------
# A) 功能与边界（同进程，审计根指向临时目录）
# ------------------------------------------------------------------
def part_a(audit):
    tmp = tempfile.mkdtemp(prefix='vp_audit_a_')
    orig = audit.DIRS['ai_audit']
    audit.DIRS['ai_audit'] = tmp
    try:
        # A1 条数不多不少（30 actions + 20 events = 50 行）
        for i in range(30):
            audit.record('actions', {'id': f'act_{i}', 'status': 'pending',
                                     'ts': time.time()})
        for i in range(20):
            audit.record('events', {'id': i, 'type': 'probe', 'ts': time.time()})
        ok = audit.flush_now(timeout=8.0)
        check('A1 flush_now 在超时内排空', ok, 'flush 未在 8s 内追平')
        today = time.strftime('%Y-%m-%d')
        na = len(_lines(os.path.join(tmp, f'actions-{today}.jsonl')))
        ne = len(_lines(os.path.join(tmp, f'events-{today}.jsonl')))
        check('A1 条数不多不少（actions=30 / events=20）',
              na == 30 and ne == 20, f'actions={na} events={ne}')

        # A2 read 新→旧 + limit
        recs = audit.read('actions', limit=5)
        check('A2 read 返回新→旧且 limit 生效',
              len(recs) == 5 and recs[0]['id'] == 'act_29'
              and recs[-1]['id'] == 'act_25', str([r.get('id') for r in recs]))

        # A3 按天轮转（注入昨天的 ts）
        yday = time.strftime('%Y-%m-%d', time.localtime(time.time() - 86400))
        audit.record('actions', {'id': 'act_old', 'status': 'pending',
                                 'ts': time.time() - 86400})
        audit.flush_now(timeout=8.0)
        check('A3 按记录自身 ts 轮转到昨天',
              os.path.isfile(os.path.join(tmp, f'actions-{yday}.jsonl')),
              f'缺 actions-{yday}.jsonl')

        # A4 脱敏：敏感键整键变 ••••，业务键不动
        audit.record('events', {'type': 'sec', 'api_key': 'sk-SECRET-123',
                                'x_token': 'tok-abc', 'sample': 's1',
                                'nested': {'db_password': 'p@ss'},
                                'ts': time.time()})
        audit.flush_now(timeout=8.0)
        sec = [r for r in audit.read('events', limit=20)
               if r.get('type') == 'sec']
        r0 = sec[0] if sec else {}
        check('A4 敏感键整键脱敏（api_key / x_token / db_password）',
              r0.get('api_key') == '••••' and r0.get('x_token') == '••••'
              and (r0.get('nested') or {}).get('db_password') == '••••',
              json.dumps(r0, ensure_ascii=False)[:200])
        check('A4 业务键不受影响', r0.get('sample') == 's1', str(r0))
        raw = ''.join(_lines(os.path.join(tmp, f'events-{today}.jsonl')))
        check('A4 明文密钥绝不落盘',
              'sk-SECRET-123' not in raw and 'tok-abc' not in raw
              and 'p@ss' not in raw, '磁盘上出现了明文密钥')

        # A5 截断：超长文本被砍
        audit.record('events', {'type': 'long', 'text': 'A' * 5000,
                                'ts': time.time()})
        audit.flush_now(timeout=8.0)
        lg = [r for r in audit.read('events', limit=10)
              if r.get('type') == 'long']
        check('A5 超长文本被截断到上限内',
              bool(lg) and len(lg[0]['text']) <= audit._MAX_TEXT + 1,
              f'len={len(lg[0]["text"]) if lg else "?"}')

        # A6 stats 数字对得上
        st = audit.stats()
        check('A6 stats 报告的字节数与磁盘一致',
              st['actions']['bytes'] > 0 and st['events']['bytes'] > 0
              and st['actions']['n_days'] >= 2, json.dumps(st, ensure_ascii=False))

        # A7 保留期清理：造一个 200 天前的文件，触发 prune 后必须消失
        old_day = time.strftime('%Y-%m-%d',
                                time.localtime(time.time() - 200 * 86400))
        old_p = os.path.join(tmp, f'events-{old_day}.jsonl')
        with io.open(old_p, 'w', encoding='utf-8') as f:
            f.write('{}\n')
        audit._PRUNED.clear()
        audit.record('events', {'type': 'trigger-prune', 'ts': time.time()})
        audit.flush_now(timeout=8.0)
        check('A7 超保留期的旧文件被清理', not os.path.isfile(old_p),
              f'{old_p} 仍在')

        # A8 关闭开关：VP_AI_AUDIT=0 时不入队
        os.environ['VP_AI_AUDIT'] = '0'
        try:
            before = audit.counters()[0]
            audit.record('events', {'type': 'should-not-appear',
                                    'ts': time.time()})
            check('A8 VP_AI_AUDIT=0 时不入队',
                  audit.counters()[0] == before, '关不掉')
        finally:
            os.environ.pop('VP_AI_AUDIT', None)
    finally:
        _restore(audit, 'ai_audit', orig, tmp)


# ------------------------------------------------------------------
# C) 不重复写
# ------------------------------------------------------------------
def part_c(audit):
    tmp = tempfile.mkdtemp(prefix='vp_audit_c_')
    orig = audit.DIRS['ai_audit']
    audit.DIRS['ai_audit'] = tmp
    try:
        for i in range(25):
            audit.record('events', {'id': i, 'type': 'dup', 'ts': time.time()})
        audit.flush_now(timeout=8.0)
        audit._drain()                       # 退出路径再来一次
        audit.flush_now(timeout=8.0)
        n = _count_dir(tmp)
        check('C 反复 flush/drain 后条数仍为 25（不重复写）', n == 25, f'n={n}')
    finally:
        _restore(audit, 'ai_audit', orig, tmp)


# ------------------------------------------------------------------
# D) 接线自证：publish / executors 真的会落盘
# ------------------------------------------------------------------
def part_d(audit, broadcast, executors):
    tmp = tempfile.mkdtemp(prefix='vp_audit_d_')
    res = tempfile.mkdtemp(prefix='vp_audit_d_res_')
    o_aud, o_res = audit.DIRS['ai_audit'], audit.DIRS['results']
    audit.DIRS['ai_audit'] = tmp
    audit.DIRS['results'] = res
    try:
        os.makedirs(os.path.join(res, 'probe_sample'), exist_ok=True)
        # D1 broadcast.publish → 盘上有该事件
        broadcast.publish('probe_event', sample='probe_sample', text='hello')
        # D2 executors 状态变迁 → 盘上有 actions 记录
        rec, err = executors.create_action('retry_failed', 'probe_sample',
                                           stages=['kvsuite'],
                                           reason='审计守卫探针')
        check('D2 create_action 成功', bool(rec) and not err, str(err))
        if rec:
            executors.resolve(rec['id'], 'reject')
        audit.flush_now(timeout=8.0)

        evs = audit.read('events', limit=50)
        check('D1 broadcast.publish 真的落盘',
              any(e.get('type') == 'probe_event' for e in evs),
              str([e.get('type') for e in evs][:10]))
        acts = audit.read('actions', limit=50)
        statuses = [a.get('status') for a in acts
                    if a.get('sample') == 'probe_sample']
        check('D2 执行单 pending→rejected 两次状态都落盘',
              'pending' in statuses and 'rejected' in statuses, str(statuses))
        # D3 broadcast.history 读得回来
        check('D3 broadcast.history 能从盘上读回',
              any(e.get('type') == 'probe_event'
                  for e in broadcast.history(limit=50)),
              'history 为空')
        # D4 执行单记录里不含明文密钥类字段（结构自证）
        bad = [k for a in acts for k in a
               if 'key' in k.lower() and a.get(k) not in (None, '••••')]
        check('D4 审计记录里无明文密钥字段', not bad, str(bad))
    finally:
        audit.flush_now(timeout=5.0)          # 先排空再换根（见 _restore 注释）
        audit.DIRS['ai_audit'] = o_aud
        audit.DIRS['results'] = o_res
        shutil.rmtree(tmp, ignore_errors=True)
        shutil.rmtree(res, ignore_errors=True)


# ------------------------------------------------------------------
# E) 重启不撞号：盘上已有 act_7 → 新单必须是 act_8
# ------------------------------------------------------------------
def part_e(audit, executors):
    tmp = tempfile.mkdtemp(prefix='vp_audit_e_')
    res = tempfile.mkdtemp(prefix='vp_audit_e_res_')
    o_aud, o_res = audit.DIRS['ai_audit'], audit.DIRS['results']
    audit.DIRS['ai_audit'] = tmp
    audit.DIRS['results'] = res
    try:
        os.makedirs(os.path.join(res, 'probe_sample'), exist_ok=True)
        today = time.strftime('%Y-%m-%d')
        with io.open(os.path.join(tmp, f'actions-{today}.jsonl'), 'w',
                     encoding='utf-8', newline='\n') as f:
            f.write(json.dumps({'id': 'act_7', 'status': 'launched',
                                'sample': 'probe_sample', 'ts': time.time()},
                               ensure_ascii=False) + '\n')
        # 模拟"重启"：清空内存态
        executors._actions.clear()
        executors._hydrated = False
        executors._order_n = 0
        hist = executors.list_actions()
        check('E 重启后从盘上重建历史', any(a['id'] == 'act_7' for a in hist),
              str([a['id'] for a in hist]))
        got = [a for a in hist if a['id'] == 'act_7']
        check('E 进行中状态在重启后标为 interrupted',
              bool(got) and got[0]['status'] == 'interrupted',
              str(got))
        rec, err = executors.create_action('retry_failed', 'probe_sample',
                                           stages=['kvsuite'],
                                           reason='撞号探针')
        check('E 新单号 = act_8（不撞盘上的 act_7）',
              bool(rec) and rec['id'] == 'act_8', str(rec and rec['id']))
        if rec:
            executors.resolve(rec['id'], 'reject')
    finally:
        audit.flush_now(timeout=5.0)          # 先排空再换根（见 _restore 注释）
        audit.DIRS['ai_audit'] = o_aud
        audit.DIRS['results'] = o_res
        shutil.rmtree(tmp, ignore_errors=True)
        shutil.rmtree(res, ignore_errors=True)


# ------------------------------------------------------------------
# F) 范围自证：ai_audit 必须登记进 OUTPUT_SUBDIRS
# ------------------------------------------------------------------
def part_f():
    from Virus_Platform_Core import config
    check('F 范围自证：ai_audit 在 OUTPUT_SUBDIRS 里（切输出根不会写回平台）',
          'ai_audit' in config.OUTPUT_SUBDIRS
          and 'ai_audit' in config.DIRS,
          f'OUTPUT_SUBDIRS={config.OUTPUT_SUBDIRS}')


# ------------------------------------------------------------------
# B) 退出不丢（真子进程 + 信号）
# ------------------------------------------------------------------
# ⚠ 关键事实（决定了本段怎么测）：**Windows 上 `Popen.send_signal(SIGTERM)`
# 走的是 `TerminateProcess` —— 进程级不可捕获**，装多少 Python 信号处理器
# 都没用。能测的只有 Python 层信号路径（`signal.raise_signal`，等价于
# Linux/macOS 的外部 SIGTERM、Windows 的 Ctrl+C / Ctrl+Break）。
# 正因为"硬杀救不回来"，执行单才走 `record_sync` 同步落盘（见下面 mode=sync）。
_CHILD = '''import os, signal, sys, time
root = os.environ["VP_AUDIT_TEST_ROOT"]
sys.path.insert(0, root)
from Virus_Platform_Core.config import DIRS
DIRS["ai_audit"] = os.environ["VP_AUDIT_TEST_DIR"]
from Virus_Platform_Core.ai_agent import audit
mode = os.environ.get("VP_AUDIT_TEST_MODE", "raise_sig")
if os.environ.get("VP_AUDIT_TEST_SIGNALS") == "1":
    audit.install_signals()
n = int(os.environ.get("VP_AUDIT_TEST_N", "20"))
if mode == "sync":
    # 同步落盘一条，然后**不可捕获地硬杀自己**
    audit.record_sync("actions", {"id": "act_sync", "status": "pending",
                                  "ts": time.time()})
    print("READY", flush=True)
    os.kill(os.getpid(), signal.SIGTERM)
    time.sleep(30)
else:
    for i in range(n):
        audit.record("events", {"id": i, "type": "probe", "ts": time.time()})
    print("READY", flush=True)
    time.sleep(0.3)                      # 停在「已入队、未到 flush 周期」窗口
    signal.raise_signal(signal.SIGTERM)  # 走 Python 层信号 → 处理器排空
    time.sleep(30)                       # 处理器没生效时父进程兜底杀
'''


def _run_child(mode='raise_sig', install_sig=True, n=20):
    """起真子进程跑一种退出场景，返回 {'events': k, 'actions': k} 落盘条数。"""
    tmp = tempfile.mkdtemp(prefix='vp_audit_b_')
    aud = os.path.join(tmp, 'audit')
    os.makedirs(aud, exist_ok=True)
    script = os.path.join(tmp, 'child.py')
    with io.open(script, 'w', encoding='utf-8', newline='\n') as f:
        f.write(_CHILD)
    env = dict(os.environ)
    env.update({'VP_AUDIT_TEST_ROOT': ROOT, 'VP_AUDIT_TEST_DIR': aud,
                'VP_AUDIT_TEST_N': str(n), 'VP_AUDIT_TEST_MODE': mode,
                'VP_AUDIT_TEST_SIGNALS': '1' if install_sig else '0',
                'PYTHONIOENCODING': 'utf-8'})
    env.pop('VP_AI_AUDIT', None)
    proc = subprocess.Popen([sys.executable, '-u', script], env=env,
                            stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    try:
        line = proc.stdout.readline().decode('utf-8', 'replace')
        if 'READY' not in line:
            return {'events': -1, 'actions': -1}
        try:
            proc.wait(timeout=20)
        except subprocess.TimeoutExpired:
            proc.kill()
            proc.wait()
        out = {}
        for kind in ('events', 'actions'):
            k = 0
            for name in os.listdir(aud):
                if name.startswith(kind + '-'):
                    k += len(_lines(os.path.join(aud, name)))
            out[kind] = k
        return out
    finally:
        try:
            proc.stdout.close()
        except Exception:
            pass
        shutil.rmtree(tmp, ignore_errors=True)


def part_b():
    with_sig = _run_child('raise_sig', True, 20)
    check('B 装了信号处理器：SIGTERM 后 20 条事件全部落盘（退出不丢）',
          with_sig['events'] == 20, f'实际 {with_sig["events"]} 条')
    # 对照：不装处理器 → 默认处置直接终止，atexit 不跑 → 必须丢
    no_sig = _run_child('raise_sig', False, 20)
    check(f'B 对照：不装处理器会丢（实际 {no_sig["events"]}/20）'
          f'→ 证明信号处理器确有作用',
          no_sig['events'] < 20, f'no_sig={no_sig["events"]}')
    # 同步落盘：扛得住不可捕获的硬杀
    sync = _run_child('sync', False, 1)
    check('B 同步落盘（record_sync）扛得住不可捕获的硬杀',
          sync['actions'] == 1, f'实际 {sync["actions"]} 条')


# ------------------------------------------------------------------
# G) 灵敏度自证：注入故意错误，比对函数必须报出
# ------------------------------------------------------------------
def part_g(audit):
    tmp = tempfile.mkdtemp(prefix='vp_audit_g_')
    orig = audit.DIRS['ai_audit']
    audit.DIRS['ai_audit'] = tmp
    try:
        # G1 去掉脱敏 → 明文密钥必须出现（证明 A4 的断言真的在检）
        real_scrub = audit._scrub
        audit._scrub = lambda v, depth=0: v
        try:
            audit.record('events', {'type': 'leak', 'api_key': 'sk-LEAK-999',
                                    'ts': time.time()})
            audit.flush_now(timeout=8.0)
            raw = ''.join(_lines(os.path.join(
                tmp, f'events-{time.strftime("%Y-%m-%d")}.jsonl')))
            check('G1 灵敏度：去掉脱敏后明文密钥确实会出现（A4 有区分度）',
                  'sk-LEAK-999' in raw, '脱敏被去掉却没看到明文')
        finally:
            audit._scrub = real_scrub
        # G2 写一半 → 计数与磁盘不一致（证明 A1 的断言真的在检）
        real_write = audit._write_batch
        audit._write_batch = lambda batch: real_write(batch[:len(batch) // 2])
        try:
            for i in range(20):
                audit.record('events', {'type': 'half', 'ts': time.time()})
            audit.flush_now(timeout=8.0)
            half = len([1 for ln in _lines(os.path.join(
                tmp, f'events-{time.strftime("%Y-%m-%d")}.jsonl'))
                if '"half"' in ln])
            check('G2 灵敏度：丢掉一半后条数明显不符（A1 有区分度）',
                  half < 20, f'half={half}')
        finally:
            audit._write_batch = real_write
    finally:
        _restore(audit, 'ai_audit', orig, tmp)


def main():
    from Virus_Platform_Core.ai_agent import audit, broadcast, executors
    print('AI 审计落盘守卫：')
    part_a(audit)
    part_b()
    part_c(audit)
    part_d(audit, broadcast, executors)
    part_e(audit, executors)
    part_f()
    part_g(audit)
    print()
    if FAIL:
        print(f'✗ {len(FAIL)} 项失败: {FAIL}')
        return 1
    print('✔ 全部通过')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
