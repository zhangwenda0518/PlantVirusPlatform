# -*- coding: utf-8 -*-
"""全局重置（/api/global/reset）与任务/队列清理的集成验证。

覆盖：
  1. _recover_interrupted 不被 tasks/queue.json 带崩、也不把它当旧任务修剪掉
     （queue.json 与任务记录同住 tasks/，但 12 位十六进制以外的名字都不是任务）
  2. cancel_all 取消运行中/排队中任务
  3. SampleQueue.reset_all 连 queued/running 一起清空并落盘
  4. purge_all 只清任务记录，results/ tool_runs/ 里的结果文件分毫不动
  5. /api/global/reset 端到端：计数正确、queue.json 存活、结果文件存活；
     include_archive=false 保留磁盘归档
  6. /results 页面带出 project 字段（全局项目筛选的数据来源）

沙箱约束：safe_open / check_path(in_platform=True) 只认平台内的路径，
所以临时数据必须建在 PLATFORM_ROOT/run/ 下（脚本结束整体删除），
不能放系统 temp。全程不触碰真实 run/results、run/tasks。
用法：python tests/_it_global_reset.py
"""
import json
import os
import shutil
import sys
import threading

for _s in (sys.stdout, sys.stderr):
    try:
        _s.reconfigure(errors='replace')
    except (AttributeError, OSError):
        pass

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.environ.setdefault('VP_NO_RECOVER', '1')

FAIL = []


def check(name, cond, extra=''):
    print(('  ✔ ' if cond else '  ✘ ') + name + (f'  {extra}' if extra else ''))
    if not cond:
        FAIL.append(name)


def _rec(tid, status):
    """造一条与 TaskManager.start() 结构一致的任务记录。"""
    return {'id': tid, 'name': f'task-{tid}', 'status': status, 'stage': '',
            'pct': 0.0, 'msg': '', 'error': '', 'log': [],
            'started': 0.0, 'finished': None, 'result': None,
            'result_preview': None, 'eta': None,
            'cancel': threading.Event(), 'procs': [],
            'weight': 'light', 'log_file': None, 'link': None,
            'restart': None, 'thread': None}


def _write_json(path, obj):
    with open(path, 'w', encoding='utf-8') as f:
        json.dump(obj, f, ensure_ascii=False)


def _task_files(tdir):
    return sorted(n for n in os.listdir(tdir)
                  if n.endswith('.json') and n != 'queue.json')


def main():
    from Virus_Platform_Core import config
    from Virus_Platform_Core.web import state as S
    from Virus_Platform_Core.web import tasks as T
    from Virus_Platform_Core.web import tasks_api
    from Virus_Platform_Core.web.samples import SampleQueue

    root = os.path.join(config.PLATFORM_ROOT, 'run', '_it_global_reset_tmp')
    real_tasks_dir = config.DIRS['tasks']
    real_results_dir = config.DIRS['results']
    real_queue_file = SampleQueue.QUEUE_FILE
    real_tm = tasks_api.tm
    real_q = S.get_sample_queue()
    try:
        shutil.rmtree(root, ignore_errors=True)
        tdir = os.path.join(root, 'tasks')
        rdir = os.path.join(root, 'results')
        rundir = os.path.join(root, 'tool_runs')
        for d in (tdir, rdir, rundir):
            os.makedirs(d, exist_ok=True)
        config.DIRS['tasks'] = tdir
        config.DIRS['results'] = rdir
        SampleQueue.QUEUE_FILE = os.path.join(tdir, 'queue.json')
        qfile = SampleQueue.QUEUE_FILE
        # _recover_interrupted 受 VP_NO_RECOVER 短路；本测试正是要验证它，
        # 而此刻 DIRS['tasks'] 已指向临时目录，放开不会碰到真实 run/tasks。
        saved_norecov = os.environ.pop('VP_NO_RECOVER', None)

        # 必须被保活的结果文件
        keep_file = os.path.join(rdir, 'DEMO', '07_report', 'report.html')
        os.makedirs(os.path.dirname(keep_file), exist_ok=True)
        with open(keep_file, 'w', encoding='utf-8') as f:
            f.write('<html>keep me</html>')
        keep_run = os.path.join(rundir, 'fastp_20260101_000000')
        os.makedirs(keep_run, exist_ok=True)
        with open(os.path.join(keep_run, 'run.log'), 'w', encoding='utf-8') as f:
            f.write('tool output')

        # ---- 1a. queue.json 最旧时不被当任务修剪，且解析不炸恢复循环 ----
        print('[1a] _recover_interrupted：queue.json 不参与任务恢复')
        ghost = 'aaaaaaaaaaaa'
        _write_json(os.path.join(tdir, ghost + '.json'),
                    {'id': ghost, 'name': 'ghost', 'status': 'running'})
        os.utime(os.path.join(tdir, ghost + '.json'), (5000, 5000))
        # queue.json 最新：旧代码会在恢复循环里先读到它（list 无 .get）而抛
        # AttributeError，被外层 except 吞掉 → 幽灵任务永远留在 running
        _write_json(qfile, [])
        os.utime(qfile, (9000, 9000))
        tm = T.TaskManager()
        check('queue.json 存活', os.path.isfile(qfile))
        with open(os.path.join(tdir, ghost + '.json'), encoding='utf-8') as f:
            g = json.load(f)
        check('遗留 running 任务被标记为 failed（防幽灵生效）',
              g.get('status') == 'failed', g.get('status'))

        # ---- 1b. 超过 100 条时的修剪不会误删 queue.json ----
        print('[1b] _recover_interrupted：修剪只针对任务记录')
        for i in range(120):
            tid = f'{i:012x}'
            p = os.path.join(tdir, tid + '.json')
            _write_json(p, {'id': tid, 'name': f't{i}', 'status': 'done'})
            os.utime(p, (1000 + i, 1000 + i))
        # queue.json 变成最旧的一个：旧代码会把它排进 files[100:] 直接删掉
        os.utime(qfile, (1, 1))
        tm = T.TaskManager()
        check('修剪后 queue.json 仍存活', os.path.isfile(qfile))
        left = _task_files(tdir)
        check('任务记录修剪到 100 条', len(left) == 100, f'实际 {len(left)}')
        check('最老的任务记录被优先修剪',
              not os.path.isfile(os.path.join(tdir, '000000000000.json')))

        # ---- 2. cancel_all ----
        print('[2] cancel_all')
        tm.tasks, tm.order = {}, []
        for tid, st in (('111111111111', 'running'), ('222222222222', 'queued'),
                        ('333333333333', 'done'), ('444444444444', 'failed')):
            tm.tasks[tid] = _rec(tid, st)
            tm.order.insert(0, tid)
        check('active_ids 只认 running/queued',
              sorted(tm.active_ids()) == ['111111111111', '222222222222'],
              str(sorted(tm.active_ids())))
        n = tm.cancel_all()
        check('cancel_all 受理 2 条', n == 2, str(n))
        check('取消后 cancel 事件已置位',
              tm.tasks['111111111111']['cancel'].is_set()
              and tm.tasks['222222222222']['cancel'].is_set())
        check('终态任务未被 cancel_all 计入',
              not tm.tasks['333333333333']['cancel'].is_set())

        # ---- 3. SampleQueue.reset_all ----
        print('[3] SampleQueue.reset_all')
        q = SampleQueue()
        q.items = [{'id': 'a', 'sample': 'S1', 'status': 'running'},
                   {'id': 'b', 'sample': 'S2', 'status': 'queued'},
                   {'id': 'c', 'sample': 'S3', 'status': 'done'}]
        q._persist()
        check('reset_all 返回被清条数 3', q.reset_all() == 3)
        check('队列内存已空', q.items == [])
        with open(qfile, encoding='utf-8') as f:
            _qd = json.load(f)
        # _persist 现写 {'items': […], 'paused': …}（暂停标记后加）；兼容旧裸列表
        _qitems = _qd.get('items') if isinstance(_qd, dict) else _qd
        check('队列落盘为空列表', _qitems == [], str(_qd))
        check('reset_all 后 queue.json 仍在', os.path.isfile(qfile))

        # ---- 4. purge_all ----
        print('[4] purge_all 只清任务记录，不碰结果文件')
        for tid in ('111111111111', '222222222222'):
            tm.tasks[tid]['status'] = 'cancelled'
        removed, active_left = tm.purge_all(include_archive=True)
        check('active_left=0', active_left == 0, str(active_left))
        check('内存任务记录已清空', tm.tasks == {}, str(list(tm.tasks)))
        left = _task_files(tdir)
        check('磁盘任务记录已清空', left == [], str(left))
        check('queue.json 未被 purge_all 删除', os.path.isfile(qfile))
        check('results/ 结果文件未被触碰', os.path.isfile(keep_file))
        check('tool_runs/ 产物未被触碰',
              os.path.isfile(os.path.join(keep_run, 'run.log')))

        print('[4b] purge_all(include_archive=False) 保留磁盘归档')
        arch = '555555555555'
        _write_json(os.path.join(tdir, arch + '.json'),
                    {'id': arch, 'name': 'arch', 'status': 'done'})
        tm.tasks[arch] = _rec(arch, 'done')
        tm.order.insert(0, arch)
        removed, _ = tm.purge_all(include_archive=False)
        check('内存记录被清', arch not in tm.tasks)
        check('磁盘归档保留', os.path.isfile(os.path.join(tdir, arch + '.json')))
        check('removed 只计内存 1 条', removed == 1, str(removed))

        # ---- 5. /api/global/reset 端到端 ----
        print('[5] /api/global/reset 端到端')
        tm.tasks, tm.order = {}, []
        for tid in ('666666666666', '777777777777'):
            tm.tasks[tid] = _rec(tid, 'done')
            tm.order.insert(0, tid)
            _write_json(os.path.join(tdir, tid + '.json'),
                        {'id': tid, 'name': 'x', 'status': 'done'})
        tasks_api.tm = tm
        S.set_sample_queue(q)
        q.items = [{'id': 'z', 'sample': 'S9', 'status': 'queued'}]
        q._persist()

        import app as appmod
        cli = appmod.app.test_client()
        # 缺 confirm 令牌必须被拒（防误触：巡检/脚本对路由发空 POST 时
        # 不能真的执行破坏性清理）
        r = cli.post('/api/global/reset', json={'include_archive': True})
        check('无 confirm 令牌 → 400（防误触）', r.status_code == 400,
              str(r.status_code))
        check('400 后未清任何任务记录（磁盘记录仍在）',
              len(_task_files(tdir)) >= 2, str(_task_files(tdir)))
        _qd = json.load(open(qfile, encoding='utf-8'))
        _qitems = _qd.get('items') if isinstance(_qd, dict) else _qd
        check('400 后队列未被清空', len(_qitems) == 1, str(_qitems))

        r = cli.post('/api/global/reset',
                     json={'confirm': 'reset', 'include_archive': True})
        check('HTTP 200', r.status_code == 200, str(r.status_code))
        d = r.get_json() or {}
        check('cleared_queue=1', d.get('cleared_queue') == 1, str(d))
        check('active_left=0', d.get('active_left') == 0, str(d))
        check('removed_task_records>=2',
              (d.get('removed_task_records') or 0) >= 2, str(d))
        check('kept_results=True', d.get('kept_results') is True)
        check('接口调用后 queue.json 仍在', os.path.isfile(qfile))
        check('接口调用后结果文件仍在', os.path.isfile(keep_file))
        check('接口调用后磁盘任务记录清空', _task_files(tdir) == [],
              str(_task_files(tdir)))

        print('[5b] include_archive=false 保留归档')
        arch2 = '888888888888'
        _write_json(os.path.join(tdir, arch2 + '.json'),
                    {'id': arch2, 'name': 'arch2', 'status': 'done'})
        r = cli.post('/api/global/reset',
                     json={'confirm': 'reset', 'include_archive': False})
        check('HTTP 200', r.status_code == 200, str(r.status_code))
        check('归档文件保留',
              os.path.isfile(os.path.join(tdir, arch2 + '.json')))

        # ---- 6. /results 带出 project ----
        print('[6] /results 页面带出 project 字段')
        prep = os.path.join(rdir, 'DEMO', '00_prep')
        os.makedirs(prep, exist_ok=True)
        _write_json(os.path.join(prep, 'input.json'),
                    {'sample': 'DEMO', 'r1': 'x.fq', 'r2': '',
                     'project': '2026-烟草'})
        r = cli.get('/results')
        body = r.get_data(as_text=True)
        check('HTTP 200', r.status_code == 200, str(r.status_code))
        check('/results 渲染出项目标签', '2026-烟草' in body)
        check('/results 渲染出样品行', 'DEMO' in body)
    finally:
        if saved_norecov is not None:
            os.environ['VP_NO_RECOVER'] = saved_norecov
        config.DIRS['tasks'] = real_tasks_dir
        config.DIRS['results'] = real_results_dir
        SampleQueue.QUEUE_FILE = real_queue_file
        tasks_api.tm = real_tm
        S.set_sample_queue(real_q)
        shutil.rmtree(root, ignore_errors=True)

    print()
    if FAIL:
        print(f'✘ {len(FAIL)} 项未通过: ' + ', '.join(FAIL))
        return 1
    print('✔ 全部通过')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
