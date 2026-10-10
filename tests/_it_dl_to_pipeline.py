# -*- coding: utf-8 -*-
"""「一键转入分析流程」（/api/dl/to_pipeline）的回归测试（离线，秒级）。

这是用户下载完成后**实际点击**的那条桥：把批次里已完成的 FASTQ 变成
平台样品（写 00_prep/input.json），可选直接入批处理队列。它的分支此前
完全没有测试覆盖，而这里的边界最容易出「下载完了却进不了流程」这类问题：

  1. 双端 r1+r2 必须成对写入；单端只写 r1
  2. 只有 _2 没有 _1 的条目必须**跳过**，且**不能在 results/ 下留空样品目录**
     （留下空目录会在样品列表里显示成一个坏样品）
  3. 失败 / 已删除 / 未完成的条目不参与
  4. 同名样品已存在（且非空）时跳过，不覆盖
  5. `runs` 过滤只转指定条目；没有任何可用文件时 400
  6. `enqueue=true` 时创建的样品要进批处理队列

`DIRS['results']` / `DIRS['downloads']` 都指向 run/ 下的临时目录，
SampleQueue 用假对象替换，绝不触碰真实样品与真实队列。
用法: python tests/_it_dl_to_pipeline.py
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

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
os.chdir(ROOT)
os.environ['VP_NO_RECOVER'] = '1'

FAIL = []


def check(label, ok, extra=''):
    print(f'  [{"PASS" if ok else "FAIL"}] {label}' + (f'  {extra}' if extra else ''),
          flush=True)
    if not ok:
        FAIL.append(label)


WORK = os.path.join(ROOT, 'run', '_it_dl_to_pipeline')
shutil.rmtree(WORK, ignore_errors=True)
DL = os.path.join(WORK, 'downloads')
RES = os.path.join(WORK, 'results')
TASKS = os.path.join(WORK, 'tasks')
BDIR = os.path.join(DL, 'itpipe_batch')
for d in (DL, RES, TASKS, BDIR):
    os.makedirs(d, exist_ok=True)

from Virus_Platform_Core import public_data                       # noqa: E402
from Virus_Platform_Core.config import DIRS                        # noqa: E402
from Virus_Platform_Core.web import state as _state                # noqa: E402
from Virus_Platform_Core.pipeline import load_sample_input         # noqa: E402

_prev_dl, _prev_res = DIRS.get('downloads'), DIRS.get('results')
_prev_tasks = DIRS.get('tasks')
DIRS['downloads'] = DL
DIRS['results'] = RES
# ⚠ 必须在 import app **之前**改 DIRS['tasks']：SampleQueue.QUEUE_FILE 是在
# 导入 samples.py 时求值的。否则即便是"假队列"没装成功的意外情况，
# 真实 run/tasks/queue.json 也会被写进测试样品（本轮实测踩过这个坑，
# 已清理）。双重保险：路径隔离 + 下面装假队列。
DIRS['tasks'] = TASKS
public_data._manager = None
mgr = public_data.get_manager()


class FakeQueue:
    """假批处理队列：只记录 add() 调用，绝不真的跑分析。"""
    def __init__(self):
        self.added = []

    def add(self, samples, stages, params, project=None):
        # 返回契约与真实 SampleQueue.add 一致：条目 dict 列表（含 id/
        # sample）。此前返回字符串列表，download.to_pipeline 的
        # bind_queue_entry(it) 拿字符串调 .get('id') 直接 500（2026-10-06
        # 批量回归定位）
        import uuid as _uuid
        self.added.append({'samples': list(samples), 'project': project})
        return [{'id': _uuid.uuid4().hex[:10], 'sample': s, 'status': 'queued',
                 'project': project or ''} for s in samples]


fakeq = FakeQueue()


def fentry(acc, name, status='done'):
    return {'acc': acc, 'url': '', 'out': os.path.join(BDIR, name),
            'size': 10, 'expect_bytes': 10, 'done_bytes': 10, 'progress': 1.0,
            'status': status, 'md5': '', 'verified': True, 'error': '',
            'db': 'ENA', 'organism': ''}


mgr.batches['itpipe'] = {
    'id': 'itpipe', 'name': 'itpipe', 'status': 'completed', 'created': '',
    'dir': BDIR, 'runs': {}, 'error': '', 'converting': False,
    'convert_msg': '', 'files': [
        fentry('PAIR1', 'PAIR1_1.fastq.gz'), fentry('PAIR1', 'PAIR1_2.fastq.gz'),
        fentry('SINGLE1', 'SINGLE1.fastq.gz'),
        fentry('ONLYR2', 'ONLYR2_2.fastq.gz'),          # 只有 _2，无 _1
        fentry('FAILED1', 'FAILED1_1.fastq.gz', 'failed'),
        fentry('PEND1', 'PEND1_1.fastq.gz', 'pending'),
        fentry('DEL1', 'DEL1_1.fastq.gz', 'deleted'),
    ]}
# 已存在的非空样品目录 → 必须跳过而不是覆盖
os.makedirs(os.path.join(RES, 'EXIST1'), exist_ok=True)
with open(os.path.join(RES, 'EXIST1', 'keep.txt'), 'w') as f:
    f.write('x')
mgr.batches['itpipe']['files'].append(fentry('EXIST1', 'EXIST1.fastq.gz'))

import app as appmod                                               # noqa: E402
# 装假队列必须在 import app 之后：app → samples.py 在导入时会把**真实**
# SampleQueue 注册进 state，先装会被覆盖。
_state.set_sample_queue(fakeq)
c = appmod.app.test_client()

print('== 1) 全部转入 ==', flush=True)
r = c.post('/api/dl/to_pipeline', json={'batch': 'itpipe'})
j = r.get_json() or {}
check('（1）返回 200', r.status_code == 200, f'status={r.status_code} {j}')
created = sorted(j.get('created') or [])
skipped = sorted(j.get('skipped') or [])
check('（1）双端与单端都被创建',
      created == ['PAIR1', 'SINGLE1'], str(created))
check('（1）只 _2 / 已存在 被跳过',
      set(skipped) >= {'ONLYR2', 'EXIST1'}, str(skipped))
check('（1）失败/未完成/已删除条目不参与',
      not ({'FAILED1', 'PEND1', 'DEL1'} & set(created)), str(created))

print('\n== 2) input.json 内容 ==', flush=True)
r1, r2, _proj = load_sample_input(os.path.join(RES, 'PAIR1'))
check('（2）双端样品 r1/r2 成对写入',
      r1 and r1.endswith('PAIR1_1.fastq.gz')
      and r2 and r2.endswith('PAIR1_2.fastq.gz'), f'r1={r1} r2={r2}')
r1s, r2s, _p = load_sample_input(os.path.join(RES, 'SINGLE1'))
check('（2）单端样品只写 r1',
      r1s and r1s.endswith('SINGLE1.fastq.gz') and not r2s,
      f'r1={r1s} r2={r2s}')
with open(os.path.join(RES, 'PAIR1', '00_prep', 'input.json'),
          encoding='utf-8') as f:
    raw = json.load(f)
check('（2）input.json 结构完整（sample/r1/r2/project）',
      {'sample', 'r1', 'r2', 'project'} <= set(raw), str(sorted(raw)))

print('\n== 3) 跳过的条目不留下空样品目录 ==', flush=True)
d_only = os.path.join(RES, 'ONLYR2')
leftover = os.path.isdir(d_only) and not any(os.scandir(d_only))
check('（3）只有 _2 的条目不留空样品目录（空目录会被当成坏样品）',
      not os.path.isdir(d_only) or not leftover,
      f'exists={os.path.isdir(d_only)} empty={leftover}')
check('（3）失败/未完成条目也不建目录',
      not any(os.path.isdir(os.path.join(RES, n))
              for n in ('FAILED1', 'PEND1', 'DEL1')))
check('（3）已存在样品的内容未被清空',
      os.path.isfile(os.path.join(RES, 'EXIST1', 'keep.txt')))

print('\n== 4) runs 过滤 ==', flush=True)
r = c.post('/api/dl/to_pipeline', json={'batch': 'itpipe',
                                       'runs': ['SINGLE1']})
j = r.get_json() or {}
check('（4）只转指定条目（PAIR1 已存在 → 跳过）',
      (j.get('created') or []) == [] and 'SINGLE1' in (j.get('skipped') or []),
      f"created={j.get('created')} skipped={j.get('skipped')}")

print('\n== 5) 无可用文件 → 400 ==', flush=True)
r = c.post('/api/dl/to_pipeline', json={'batch': 'itpipe',
                                       'runs': ['FAILED1']})
check('（5）指定不可用条目时 400', r.status_code == 400,
      f'status={r.status_code} {str(r.get_json())[:80]}')
r = c.post('/api/dl/to_pipeline', json={'batch': 'no_such_batch'})
check('（5）批次不存在时 400', r.status_code == 400,
      f'status={r.status_code}')
r = c.post('/api/dl/to_pipeline', json={'batch': 'itpipe..evil'})
check('（5）非法批次名被拒', r.status_code == 400, f'status={r.status_code}')

print('\n== 6) enqueue 入队 ==', flush=True)
shutil.rmtree(os.path.join(RES, 'PAIR1'), ignore_errors=True)
shutil.rmtree(os.path.join(RES, 'SINGLE1'), ignore_errors=True)
r = c.post('/api/dl/to_pipeline',
           json={'batch': 'itpipe', 'enqueue': True, 'project': 'projX'})
j = r.get_json() or {}
check('（6）enqueue 时返回 queued 数量',
      j.get('queued') == len(j.get('created') or []) > 0,
      f"created={j.get('created')} queued={j.get('queued')}")
check('（6）队列收到新建样品（含 project）',
      bool(fakeq.added) and set(fakeq.added[-1]['samples'])
      == set(j.get('created') or [])
      and fakeq.added[-1]['project'] == 'projX',
      str(fakeq.added[-1] if fakeq.added else None))
r = c.post('/api/dl/to_pipeline', json={'batch': 'itpipe'})
check('（6）不 enqueue 时不入队',
      (r.get_json() or {}).get('queued') == 0,
      str((r.get_json() or {}).get('queued')))

DIRS['downloads'], DIRS['results'] = _prev_dl, _prev_res
DIRS['tasks'] = _prev_tasks
public_data._manager = None
mgr.batches.pop('itpipe', None)
# 断言真实队列文件没有被本次测试碰过
_q = os.path.join(ROOT, 'run', 'tasks', 'queue.json')
check('（7）真实 run/tasks/queue.json 未被测试写入',
      (not os.path.isfile(_q))
      or all(s not in open(_q, encoding='utf-8').read()
             for s in ('PAIR1', 'SINGLE1')),
      f'{_q}')
shutil.rmtree(WORK, ignore_errors=True)
print('\n' + ('ALL PASS' if not FAIL else f'{len(FAIL)} FAIL: {FAIL}'), flush=True)
sys.exit(0 if not FAIL else 1)
