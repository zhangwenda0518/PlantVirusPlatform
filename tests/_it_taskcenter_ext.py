# -*- coding: utf-8 -*-
"""阶段二：另外两处「外部任务」挂载 + attach 重绑定语义（离线，秒级）。

阶段二把三样「本来完全不在任务中心」的东西挂成了外部任务：
  1. 公共数据下载批次 —— DownloadManager（真浏览器测试在
     tests/_it_taskcenter_dl_ui.py）
  2. **批处理队列里「排队中」的样品** —— 原先只有真正开跑才 tm.start
     出「队列·<样品>」，排队期间用户点了入队却在任务中心看不见
  3. **存储占用扫描** —— 后台线程走整棵目录树（含 databases/），
     原先只有设置页一个 scanning 布尔

外部任务的契约（本测试锁死）：
  - provider 返回 dict → 卡片按该 dict 渲染；返回 None → 真身没了、卡片移除
  - provider 抛异常 → 卡片**降级保留**，不凭空消失
  - 不占 heavy/light 并发名额、不落 tasks/<tid>.json（否则重启后会被
    _recover_interrupted 当成幽灵任务标失败）
  - 同一 key 重复 attach：tid 稳定不变，但 provider/on_cancel 换成最新闭包
    （执行对象被重建时卡片必须跟着走，否则快照永远 None → 卡片被误删）
  - 「停止」转发到真身的取消函数

全程离线、只动内存里的假队列；不写任何真实批次/队列文件。
用法: python tests/_it_taskcenter_ext.py
"""
import json
import os
import sys
import threading
import time

for _s in (sys.stdout, sys.stderr):
    try:
        _s.reconfigure(errors='replace')
    except (AttributeError, OSError):
        pass

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
os.chdir(ROOT)
os.environ.setdefault('VP_NO_RECOVER', '1')

FAIL = []


def check(label, ok, extra=''):
    print(f'  [{"PASS" if ok else "FAIL"}] {label}'
          + (f'  {extra}' if extra else ''), flush=True)
    if not ok:
        FAIL.append(label)


import app as appmod                                              # noqa: E402
from Virus_Platform_Core.config import DIRS                       # noqa: E402
from Virus_Platform_Core.web import samples as samples_mod        # noqa: E402
from Virus_Platform_Core.web import settings as settings_mod      # noqa: E402
from Virus_Platform_Core.web.tasks import tm                      # noqa: E402

client = appmod.app.test_client()


def names(logs=False):
    return [x['name'] for x in tm.list_all(logs=logs) if x]


def tid_of(name):
    for x in tm.list_all(logs=False):
        if x and x['name'] == name:
            return x['id']
    return None


def api_raw():
    """/api/tasks?full=1 的形状是 {'active': [...], 'archived': ...}。"""
    r = client.get('/api/tasks?full=1&logs=0')
    assert r.status_code == 200, r.status_code
    j = r.get_json()
    if isinstance(j, dict):
        return j.get('active') or []
    return j or []


def api_names():
    return [x['name'] for x in api_raw() if x]


# ══════════════ 1) 排队中的样品 ══════════════
print('\n== 1) 批处理队列「排队中」的样品挂进任务中心 ==', flush=True)


class FakeQueue:
    """只实现 bind_queue_entry 用到的那点接口（lock/items/remove）。"""

    def __init__(self, items):
        self.lock = threading.RLock()
        self.items = items
        self.removed = []

    def remove(self, eid):
        with self.lock:
            self.removed.append(eid)
            n = len(self.items)
            self.items[:] = [x for x in self.items if x['id'] != eid]
            return len(self.items) < n


_real_q = samples_mod.sample_queue
fake_q = FakeQueue([{'id': 'e1', 'sample': 'ITEXT1', 'status': 'queued'},
                    {'id': 'e2', 'sample': 'ITEXT2', 'status': 'running'}])
samples_mod.sample_queue = fake_q
try:
    samples_mod.bind_queue_entry({'id': 'e1', 'sample': 'ITEXT1'})
    samples_mod.bind_queue_entry({'id': 'e2', 'sample': 'ITEXT2'})

    check('（1）排队中的样品出现在任务中心', '排队·ITEXT1' in names(),
          str([n for n in names() if '排队' in n]))
    check('（1）已开跑的样品**不**重复挂排队卡（避免两张）',
          '排队·ITEXT2' not in names(), str(names()))
    check('（1）/api/tasks 也能看到', '排队·ITEXT1' in api_names())

    tid_q = tid_of('排队·ITEXT1')
    s = tm.snapshot(tid_q, log_lines=10) if tid_q else None
    check('（1）阶段文案为「排队中」', bool(s) and s.get('stage') == '排队中',
          (s or {}).get('stage'))
    check('（1）进度为 0（不是 None，前端要能画空进度条）',
          bool(s) and s.get('pct') == 0.0, str((s or {}).get('pct')))
    # 2026-09-19 起扫描导入改为「入队后待手动启动」（33a4bda）：
    # queued 且未点 ▶ 启动 时文案是「已入队，待手动启动」；
    # 「等待批处理名额」只出现在手动启动后的排队期（samples.py:388）。
    check('（1）文案说明已入队待手动启动', bool(s) and '已入队' in s.get('msg', ''),
          (s or {}).get('msg', '')[:80])
    check('（1）「前往页面」指向结果中心（2026-10-11 /pipeline 归档后改指）',
          bool(s) and s.get('link') == '/results',
          str((s or {}).get('link')))
    check('（1）外部任务不落盘（无 tasks/<tid>.json）',
          not os.path.exists(os.path.join(DIRS['tasks'], f'{tid_q}.json')),
          os.path.join('tasks', f'{tid_q}.json'))

    # 开跑 → provider 返回 None → 卡片让位给 tm 的真任务卡
    with fake_q.lock:
        fake_q.items[0]['status'] = 'running'
    check('（1）开跑后排队卡自动让位（下次列取即消失）',
          '排队·ITEXT1' not in names(), str(names()))
    check('（1）让位后 /api/tasks 也没有残留（不是 [null]）',
          '排队·ITEXT1' not in api_names()
          and all(x is not None for x in api_raw()), str(api_names()))

    # 又回到排队（重排）→ 重新 attach，tid 稳定
    with fake_q.lock:
        fake_q.items[0]['status'] = 'queued'
    samples_mod.bind_queue_entry({'id': 'e1', 'sample': 'ITEXT1'})
    check('（1）重新排队后卡片回来且 tid 稳定', tid_of('排队·ITEXT1') == tid_q,
          f'{tid_of("排队·ITEXT1")} vs {tid_q}')

    # 「停止」→ 转发到真身（把队列条目真的删掉）
    ok = tm.cancel(tid_q)
    check('（1）点停止转发到队列（真的移除了条目）',
          ok is True and fake_q.removed == ['e1'], f'ok={ok} removed={fake_q.removed}')
    check('（1）停止后卡片消失', '排队·ITEXT1' not in names(), str(names()))
finally:
    samples_mod.sample_queue = _real_q
    tm.detach('queue:e1')
    tm.detach('queue:e2')

# ══════════════ 2) 存储占用扫描 ══════════════
print('\n== 2) 存储占用扫描挂进任务中心 ==', flush=True)
_tid_scan = None
try:
    with settings_mod._storage_lock:
        _scan_prev = settings_mod._storage['scanning']
        settings_mod._storage['scanning'] = True
    settings_mod._bind_storage_scan()
    _tid_scan = tid_of('存储占用扫描')
    check('（2）扫描中出现在任务中心', bool(_tid_scan), str(_tid_scan))
    s2 = tm.snapshot(_tid_scan, log_lines=10) if _tid_scan else None
    check('（2）阶段文案为「统计目录占用」',
          bool(s2) and s2.get('stage') == '统计目录占用', (s2 or {}).get('stage'))
    check('（2）「前往页面」指向设置页',
          bool(s2) and s2.get('link') == '/settings', str((s2 or {}).get('link')))
    check('（2）扫描中 /api/tasks 可见', '存储占用扫描' in api_names())

    with settings_mod._storage_lock:
        settings_mod._storage['scanning'] = False
    check('（2）扫描结束后卡片自动消失', '存储占用扫描' not in names(), str(names()))
finally:
    with settings_mod._storage_lock:
        settings_mod._storage['scanning'] = _scan_prev
    tm.detach('settings:storage-scan')

# ══════════════ 3) attach 重绑定 / 异常降级 ══════════════
print('\n== 3) attach 重绑定与 provider 异常降级 ==', flush=True)
KEY = 'ext:ittest-rebind'
try:
    t1 = tm.attach('EXT·A', lambda: {'status': 'running', 'msg': 'A',
                                     'pct': 0.1, 'stage': 'sa'},
                   key=KEY, link='/x', weight='light')
    check('（3）第一次 attach 登记成功', bool(t1), str(t1))
    check('（3）快照取的是 A',
          (tm.snapshot(t1, log_lines=5) or {}).get('msg') == 'A',
          str((tm.snapshot(t1, log_lines=5) or {}).get('msg')))

    # 模拟「执行对象被重建」：同 key 再 attach 一个全新 provider
    t2 = tm.attach('EXT·A', lambda: {'status': 'running', 'msg': 'B',
                                     'pct': 0.9, 'stage': 'sb'},
                   key=KEY, link='/y', weight='light')
    check('（3）重复 attach 的 tid 稳定不变', t2 == t1, f'{t2} vs {t1}')
    s3 = tm.snapshot(t1, log_lines=5) or {}
    check('（3）重复 attach 后 provider 已换成最新闭包（msg=B）',
          s3.get('msg') == 'B', str(s3.get('msg')))
    check('（3）重复 attach 后 link 也刷新', s3.get('link') == '/y',
          str(s3.get('link')))

    # provider 报错 → 降级卡片保留（不是凭空消失）
    tm.tasks[t1]['provider'] = lambda: (_ for _ in ()).throw(RuntimeError('x'))
    s4 = tm.snapshot(t1, log_lines=5)
    check('（3）provider 报错时降级保留而非消失', bool(s4), str(s4))
    check('（3）降级文案提示会自动恢复',
          bool(s4) and '状态读取失败' in s4.get('msg', ''), (s4 or {}).get('msg'))
    check('（3）降级卡片仍带 started（不会算成 49 万小时）',
          bool(s4) and bool(s4.get('started')), str((s4 or {}).get('started')))
    check('（3）降级卡片仍出现在 /api/tasks',
          'EXT·A' in api_names(), str(api_names()))

    # provider 干净返回 None → 卡片移除
    tm.tasks[t1]['provider'] = lambda: None
    check('（3）provider 干净返回 None 时卡片移除', tm.snapshot(t1) is None)
    check('（3）移除后 /api/tasks 无残留', 'EXT·A' not in api_names(),
          str(api_names()))
finally:
    tm.detach(KEY)

print('\n== 4) 外部任务不占并发名额 ==', flush=True)
_h0, _l0 = tm.heavy_slots._value, tm.light_slots._value
_probe = tm.attach('EXT·PROBE', lambda: {'status': 'running', 'msg': 'p'},
                   key='ext:ittest-probe')
check('（4）挂载期间 heavy/light 名额一个都没被占',
      (tm.heavy_slots._value, tm.light_slots._value) == (_h0, _l0),
      f'heavy {_h0}→{tm.heavy_slots._value} '
      f'light {_l0}→{tm.light_slots._value}')
check('（4）外部任务不写 tasks/<tid>.json',
      not os.path.exists(os.path.join(DIRS['tasks'], f'{_probe}.json')),
      _probe)
check('（4）外部任务也在「运行中」列表里',
      _probe in tm.active_ids(), str(tm.active_ids()[:5]))
tm.detach('ext:ittest-probe')
check('（4）detach 后彻底消失', tm.snapshot(_probe) is None)

print('\n== 5) 跑完的外部任务必须退出 active_ids（否则全局重置会误取消） ==',
      flush=True)
# 注意：这里**故意不调用 tm.cancel_all()** —— 它会把取消转发给所有活跃
# 任务，其中包括用户真实在跑的下载批次。测试绝不允许有杀掉用户下载的
# 可能，所以只断言 cancel_all / 全局重置唯一依据的那个谓词
# （active_ids），再单独验证「转发取消」这条机制本身。
_KEY5 = 'ext:ittest-active'
try:
    _box = {'status': 'running', 'msg': 'x'}
    _t5 = tm.attach('EXT·ACT', lambda: dict(_box), key=_KEY5, weight='light')
    tm.list_all(logs=False)
    check('（5）运行中的外部任务在 active_ids 里', _t5 in tm.active_ids(),
          str(tm.active_ids()[:5]))

    _box['status'] = 'done'
    got = [x for x in tm.list_all(logs=False) if x and x['id'] == _t5]
    check('（5）provider 报 done 后仍出现在任务列表（卡片要留在那）',
          bool(got),
          str([x['name'] for x in tm.list_all(logs=False) if x][:5]))
    check('（5）同时退出 active_ids（不再被当成「运行中」）',
          _t5 not in tm.active_ids(), str(tm.active_ids()[:5]))
    check('（5）快照状态也如实是 done',
          (tm.snapshot(_t5, log_lines=5) or {}).get('status') == 'done',
          str((tm.snapshot(_t5, log_lines=5) or {}).get('status')))

    # 复刻 cancel_all()/全局重置的唯一遍历依据（只对 active_ids 转发取消）。
    # 不直接调 cancel_all()，因为那会连用户真实在跑的下载一起停掉。
    _hits = []
    tm.tasks[_t5]['on_cancel'] = lambda: _hits.append(1)
    for _tid in tm.active_ids():
        if _tid == _t5:
            tm.cancel(_tid)
    check('（5）跑完后：按全局重置的遍历方式碰不到它（不会被误取消）',
          not _hits, str(_hits))

    _box['status'] = 'running'
    tm.list_all(logs=False)
    check('（5）回到运行中后 active_ids 又能收录它（条件真实生效）',
          _t5 in tm.active_ids(), str(tm.active_ids()[:5]))
    for _tid in tm.active_ids():
        if _tid == _t5:
            tm.cancel(_tid)
    check('（5）此时按同一遍历方式会真的转发取消',
          _hits == [1], f'hits={_hits}')
finally:
    tm.detach(_KEY5)

print('\n' + ('ALL PASS' if not FAIL else f'{len(FAIL)} FAIL: {FAIL}'), flush=True)
sys.exit(0 if not FAIL else 1)
