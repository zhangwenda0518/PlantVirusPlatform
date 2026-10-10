# -*- coding: utf-8 -*-
"""任务终态回调 on_finish 守卫（2026-10-07 审查 P1-1）。

背景
----
`TaskManager._run` 在「排队途中被取消」分支里**直接 return，从不调用 fn**，
所以 `fn` 的 `finally` 不会执行。样品分析路径的清理（`_sample_busy_release`）
原本只写在 job 的 finally 里 → 该样品在 `_RUNNING_SAMPLES` 里**永久残留**，
此后重新入队/删除/清空结果全部被拒，直到重启进程。

修复：`TaskManager.start(..., on_finish=cb)` —— 终态必调一次（成功/失败/取消，
含排队取消）。

判据（四条）
-----------
A. **灵敏度（对照臂）**：从 `archive/_bak_20261007/` 加载的**旧 TaskManager**
   在「排队取消」时**不会**执行 fn 的 finally —— 证明旧行为确实会泄漏。
B. **新行为**：同一场景下 `on_finish` 被调用**恰一次**。
C. **正常路径**：未被取消的任务跑完后 `on_finish` 也被调用恰一次（不能只在
   取消分支生效）。
D. **不重复**：排队取消场景下 on_finish 只调一次（取消分支与 finally 互斥）。

用法: python tests/_check_task_on_finish.py
退出码: 0 = 通过；1 = 失败。
"""
import importlib.util
import os
import sys
import tempfile
import threading
import time

for _s in (sys.stdout, sys.stderr):
    try:
        _s.reconfigure(errors='replace')
    except (AttributeError, OSError):
        pass

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

# ⚠️ 必须在 import tasks 之前：TaskManager() 是模块级单例，构造时会
# _recover_interrupted() 改写 run/tasks/。VP_NO_RECOVER 是平台为此提供的
# 只读开关（tasks.py 文档串里点名给只读巡检脚本用）。
os.environ['VP_NO_RECOVER'] = '1'

BACKUP = os.path.join(ROOT, 'archive', '_bak_20261007',
                      'Virus_Platform_Core', 'web', 'tasks.py')

FAILS = []


def _check(cond, msg):
    print(('  ok   ' if cond else '  FAIL ') + msg)
    if not cond:
        FAILS.append(msg)


def load_backup_tasks():
    """把修复前的 tasks.py 当独立模块加载（其 import 全为绝对路径）。"""
    spec = importlib.util.spec_from_file_location('_vp_tasks_bak', BACKUP)
    m = importlib.util.module_from_spec(spec)
    sys.modules['_vp_tasks_bak'] = m
    spec.loader.exec_module(m)
    return m


def make_mgr(mod, tmpdir):
    """构造一个不与真实 run/tasks/ 交互的 TaskManager（_persist 打桩）。"""
    mgr = mod.TaskManager()
    mgr.heavy_slots = threading.Semaphore(1)      # 压到 1，便于制造排队
    mgr.light_slots = threading.Semaphore(1)
    mgr._persist = lambda rec: None               # 不写 run/tasks/*.json
    mgr._trim_memory = lambda: None
    return mgr


def run_scenario(mod, tmpdir, tag, use_on_finish):
    """返回 (job_finally_called, on_finish_calls)。"""
    mgr = make_mgr(mod, tmpdir)
    hold = threading.Event()
    state = {'finally': 0, 'on_finish': 0}

    def blocker(emit, prog, cancel):
        hold.wait(10)
        return 'ok'

    def job(emit, prog, cancel):
        try:
            return 'ok'
        finally:
            state['finally'] += 1

    cb = (lambda: state.__setitem__('on_finish', state['on_finish'] + 1)) \
        if use_on_finish else None

    t1 = mgr.start(f'{tag}-blocker', blocker,
                   log_file=os.path.join(tmpdir, f'{tag}_1.log'))
    time.sleep(0.3)
    kw = {'log_file': os.path.join(tmpdir, f'{tag}_2.log')}
    if cb is not None:
        kw['on_finish'] = cb
    t2 = mgr.start(f'{tag}-queued', job, **kw)
    time.sleep(0.3)
    mgr.cancel(t2)                     # 排队途中取消
    time.sleep(1.2)
    snap = mgr.snapshot(t2, log_lines=0, include_result=False)
    hold.set()
    mgr.cancel(t1)
    time.sleep(0.3)
    return state, (snap or {})


def main():
    if not os.path.isfile(BACKUP):
        print(f'✘ 找不到修复前备份：{BACKUP}')
        return 1

    from Virus_Platform_Core.web import tasks as live
    bak = load_backup_tasks()
    print(f'对照臂：旧 TaskManager 抽自 {os.path.relpath(BACKUP, ROOT)}')

    tmpdir = tempfile.mkdtemp(prefix='vp_onfinish_')
    try:
        # ---------- A. 对照臂：旧实现在排队取消时不执行 fn 的 finally ----------
        print('\n[A] 灵敏度（对照臂）：旧实现在「排队取消」时 fn 的 finally 不执行')
        st_old, snap_old = run_scenario(bak, tmpdir, 'old', use_on_finish=False)
        _check(snap_old.get('status') == 'cancelled'
               and '排队' in (snap_old.get('error') or ''),
               f'旧：任务确为排队取消（status={snap_old.get("status")}, '
               f'error={snap_old.get("error")}）')
        _check(st_old['finally'] == 0,
               f'旧：fn 的 finally 执行次数 = {st_old["finally"]}（应 0，即会泄漏）')

        # ---------- B. 新行为：on_finish 恰一次 ----------
        print('\n[B] 新实现：同一场景下 on_finish 被调用恰一次')
        st_new, snap_new = run_scenario(live, tmpdir, 'new', use_on_finish=True)
        _check(snap_new.get('status') == 'cancelled'
               and '排队' in (snap_new.get('error') or ''),
               f'新：任务确为排队取消（status={snap_new.get("status")}）')
        _check(st_new['finally'] == 0,
               f'新：fn 的 finally 仍为 {st_new["finally"]}（证明泄漏窗口真实存在）')
        _check(st_new['on_finish'] == 1,
               f'新：on_finish 调用次数 = {st_new["on_finish"]}（应恰 1）')

        # ---------- C. 正常路径：on_finish 也要被调用 ----------
        print('\n[C] 正常路径：未取消的任务跑完后 on_finish 也被调用恰一次')
        mgr = make_mgr(live, tmpdir)
        calls = {'n': 0, 'fin': 0}

        def ok_job(emit, prog, cancel):
            try:
                return 'ok'
            finally:
                calls['fin'] += 1

        tid = mgr.start('normal', ok_job,
                        log_file=os.path.join(tmpdir, 'normal.log'),
                        on_finish=lambda: calls.__setitem__('n', calls['n'] + 1))
        for _ in range(60):
            s = mgr.snapshot(tid, log_lines=0, include_result=False)
            if s and s.get('status') != 'running':
                break
            time.sleep(0.2)
        _check(calls['fin'] == 1, f'正常路径 fn 的 finally = {calls["fin"]}（应 1）')
        _check(calls['n'] == 1, f'正常路径 on_finish = {calls["n"]}（应恰 1）')

        # ---------- D. 失败路径：on_finish 仍恰一次 ----------
        print('\n[D] 失败路径：抛异常的任务 on_finish 仍恰一次')
        mgr2 = make_mgr(live, tmpdir)
        calls2 = {'n': 0}

        def bad_job(emit, prog, cancel):
            raise RuntimeError('故意失败')

        tid2 = mgr2.start('bad', bad_job,
                          log_file=os.path.join(tmpdir, 'bad.log'),
                          on_finish=lambda: calls2.__setitem__('n', calls2['n'] + 1))
        for _ in range(60):
            s = mgr2.snapshot(tid2, log_lines=0, include_result=False)
            if s and s.get('status') != 'running':
                break
            time.sleep(0.2)
        _check(calls2['n'] == 1, f'失败路径 on_finish = {calls2["n"]}（应恰 1）')
    finally:
        for f in os.listdir(tmpdir):
            os.remove(os.path.join(tmpdir, f))
        os.rmdir(tmpdir)

    print()
    if FAILS:
        print(f'✘ on_finish 守卫失败（{len(FAILS)} 条）')
        return 1
    print('✓ on_finish 守卫通过（排队取消/正常/失败三种终态都恰调一次）')
    return 0


if __name__ == '__main__':
    sys.exit(main())
