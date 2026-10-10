# -*- coding: utf-8 -*-
"""后台任务引擎（自 app.py 拆出）。

TaskManager：后台线程执行 + 状态快照（内存为主，tasks/ 留一份 JSON）。
并发控制：heavy / light 两档 Semaphore，任务先排队（status='queued'）
再等名额；日志经 SSE（/api/task/<id>/stream）实时推给前端。

内存保留策略：最近 _KEEP_FULL 个任务保留完整信息，更早的清空日志/结果
引用只留状态摘要；超过 _KEEP_TOTAL 的移出内存（磁盘 tasks/*.json 仍可回溯）。

注：原实现用 `app.logger.warning`，拆出后改用本模块 logger，避免对 Flask
应用对象的模块级依赖（后台线程内 app 对象不可靠）。
"""
import hashlib
import json
import logging
import os
import re
import subprocess
import sys
import threading
import time
import uuid
from collections import deque

from Virus_Platform_Core.config import DIRS
from Virus_Platform_Core.utils import (check_path, decode_output, fmt_size,
                                       safe_open, safe_remove)
from Virus_Platform_Core.web.state import cfg, tool_runs_root as _tool_runs_root

_log = logging.getLogger('Virus_Platform_Core.web.tasks')


def _slot_limits():
    """并发闸门上限（platform.json defaults 可覆盖：max_heavy_tasks /
    max_light_tasks）；未配置时用缺省值。

    heavy = 吃满多核/GB 级内存的任务（kunpeng 分类、SPAdes 组装、
    DIAMOND 注释、建库、建树、SDT）；light = 网络 IO 或轻量计算
    （下载、格式转换、绘图、检索）。
    """
    d = getattr(cfg, 'defaults', None) or {}
    try:
        heavy = max(1, min(int(d.get('max_heavy_tasks', 2)), 8))
    except (TypeError, ValueError):
        heavy = 2
    try:
        light = max(1, min(int(d.get('max_light_tasks', 4)), 16))
    except (TypeError, ValueError):
        light = 4
    return heavy, light


_KEEP_FULL = 60


_KEEP_TOTAL = 300


class TaskManager:
    """后台线程任务 + 状态快照（内存为主，tasks/ 目录留一份 JSON）。

    并发控制：heavy / light 两档 Semaphore。任务先排队（status='queued'，
    对外仍呈现为 running，仅 msg 标注排队位次），拿到名额后才真正开跑。
    """

    def __init__(self):
        self.tasks = {}
        self.order = []
        self.lock = threading.Lock()
        self._ext_keys = {}          # 外部任务：稳定 key -> task id
        heavy, light = _slot_limits()
        self.heavy_slots = threading.Semaphore(heavy)
        self.light_slots = threading.Semaphore(light)
        self.limits = {'heavy': heavy, 'light': light}
        # heavy 名额实时占用数（闸门-线程联动：pipeline 开跑时按
        # 物理核 ÷ 占用数 钳制本任务线程，防 2 任务 × 各自 19 线程
        # 超订 20 核——2026-10-04 审查遗留项）
        self._heavy_used = 0
        self._recover_interrupted()

    def heavy_active(self):
        """当前占用的 heavy 名额数（≥1，保底供除法）。"""
        with self.lock:
            return max(1, self._heavy_used)

    def _recover_interrupted(self):
        """服务启动时把上次进程遗留的 running 任务标记为中断（防幽灵），
        并把任务状态文件修剪到最近 100 个（防 tasks/ 无限膨胀）。

        设 `VP_NO_RECOVER=1` 可跳过（只读巡检脚本 `import app` 时不应改动
        磁盘：tests/_route_inventory.py 自称只读，却会触发这里的修剪与改写）。
        """
        if os.environ.get('VP_NO_RECOVER'):
            return
        try:
            # 只认合法任务 id（12 位十六进制）：tasks/ 下还住着
            # SampleQueue 的 queue.json，它不是任务记录——不排除掉的话
            # ① 修剪到 100 个时会把 queue.json 当旧任务删掉（整条批处理队列丢失）
            # ② json.load 得到 list，下面 d.get() 抛 AttributeError，
            #    被外层 except 吞掉后整个中断恢复循环失效（幽灵 running 任务）。
            files = [n for n in os.listdir(DIRS['tasks'])
                     if n.endswith('.json') and self._tid_ok(n[:-5])]
            files.sort(key=lambda n: os.path.getmtime(
                os.path.join(DIRS['tasks'], n)), reverse=True)
            for fn in files[100:]:                       # 只保留最近 100 个
                try:
                    safe_remove(os.path.join(DIRS['tasks'], fn))
                except (OSError, ValueError):
                    pass
            files = files[:100]
            for name in files:
                p = check_path(os.path.join(DIRS['tasks'], name),
                               must_exist=False, in_platform=True)
                try:
                    with safe_open(p) as f:
                        d = json.load(f)
                except (OSError, ValueError):
                    continue
                if d.get('status') in ('running', 'queued'):
                    d['status'] = 'failed'
                    d['error'] = '服务重启导致任务中断，请重新运行'
                    d['finished'] = time.time()
                    with safe_open(p, 'wt') as f:
                        json.dump(d, f, ensure_ascii=False)
        except Exception as e:
            # 恢复失败的后果：幽灵 running 任务照旧、旧任务文件无痕堆积，
            # 必须留痕而非静默
            print(f'[tasks] 启动时中断任务恢复失败: {e}', file=sys.stderr)

    def start(self, name, fn, log_file=None, weight='heavy', link=None,
              out_dir=None, on_finish=None):
        """创建任务。log_file: 任务日志持久化文件（缺省落 logs/tasks/）。

        weight: 'heavy'（吃多核/大内存，受 heavy 闸门限流）或 'light'
        （网络 IO / 轻量计算，受 light 闸门限流）。缺省 heavy——宁可多排队，
        也不要把不认识的任务放过去打满机器。

        link: 结果/操作页跳转地址（任务卡「前往」按钮）。

        out_dir: 本任务的输出目录（绝对路径）。任务卡与结果面板据此把
        「结果输出在哪」直接显示给用户（含打开/复制入口）；不传则无此行。

        on_finish: **任务进入终态时必调一次**的回调（成功/失败/取消，含
        「排队途中被取消」）。⚠️ 排队中取消时 `fn` **根本不会被调用**，所以
        清理逻辑**不能**只写在 `fn` 的 finally 里 —— 那会永久泄漏（2026-10-07
        审查 P1-1：样品占位表 `_RUNNING_SAMPLES` 因此永不摘除）。需要
        「无论任务是否真正开跑都要执行的收尾」请用本参数，而非 fn 的 finally。
        回调异常只打 stderr，不影响任务终态。

        重启上下文自动捕获：任务在 POST API handler 里发起时，记录
        request.path + JSON body，任务卡可一键重新提交同一请求。
        仅存内存（不落盘），避免 API key 等敏感字段写入 tasks/*.json。

        log() 同时写内存 deque（前端实时读）与磁盘文件（重启后可查）。
        """
        tid = uuid.uuid4().hex[:12]
        restart = None
        try:
            from flask import has_request_context, request as _req
            if has_request_context() and _req.method == 'POST':
                body = _req.get_json(silent=True)
                restart = {'url': _req.path, 'body': body} if body else None
        except Exception:
            restart = None
        if not log_file:
            safe = re.sub(r'[^\w\-.]+', '_', name, flags=re.UNICODE)
            safe = safe[:40].strip('_ ') or 'task'
            log_dir = check_path(os.path.join(DIRS['logs'], 'tasks'),
                                 must_exist=False, in_platform=True)
            os.makedirs(log_dir, exist_ok=True)
            log_file = os.path.join(
                log_dir, f'{safe}_{time.strftime("%Y%m%d_%H%M%S")}_{tid}.log')
        rec = {'id': tid, 'name': name, 'status': 'queued', 'stage': '',
               'pct': 0.0, 'msg': '', 'error': '', 'log': deque(maxlen=500),
               'started': time.time(), 'finished': None, 'result': None,
               'result_preview': None, 'eta': None,
               'cancel': threading.Event(),
               'procs': [],          # 本任务启动的子进程（硬停止用）
               'weight': weight if weight in ('heavy', 'light') else 'heavy',
               'log_file': log_file, 'link': link, 'restart': restart,
               'out_dir': out_dir}
        with self.lock:
            self.tasks[tid] = rec
            self.order.insert(0, tid)

        # 命名为 _emit 而非 _log：本闭包（追加一行到任务日志）与模块级
        # logger `_log` 同名会把它遮蔽掉 —— 曾导致任务失败分支的
        # `_log.debug(堆栈)` 抛 AttributeError，异常堆栈永远进不了日志。
        def _emit(line):
            with self.lock:
                rec['log'].append(line)
            try:
                with safe_open(rec['log_file'], 'at') as f:
                    f.write(time.strftime('[%H:%M:%S] ') + line + '\n')
            except Exception:
                pass

        def _progress(stage, pct, msg, *eta):
            # 加锁写入：snapshot() 在锁内组装字段，不锁会出现
            # 「status=done 但 result 仍为空」这类半更新帧（前端据此显示
            # 已完成却空结果面板）。
            with self.lock:
                rec['stage'], rec['pct'], rec['msg'] = stage, pct, msg
                if eta:
                    try:
                        rec['eta'] = max(float(eta[0]), 0)
                    except (TypeError, ValueError):
                        rec['eta'] = None
                # chunk 级进度回调可达每秒数十次：每次全量重写 tasks/<id>.json
                # （含 log_tail）是写盘风暴——进度帧 ≥1.2s 才落一次盘，
                # 状态/终态的直接 _persist 调用不受影响
                _now = time.time()
                _flush = (_now - rec.get('_last_prog_flush', 0)) >= 1.2
                if _flush:
                    rec['_last_prog_flush'] = _now
            if _flush:
                self._persist(rec)

        def _run():
            from Virus_Platform_Core.utils import task_bind
            slot = (self.light_slots if rec['weight'] == 'light'
                    else self.heavy_slots)
            # on_finish 必须「恰一次」：排队取消分支与 finally 两条终态路径
            # 都会调它（2026-10-07 审查 P1-1）
            _fin = {'called': False}

            def _call_on_finish():
                if _fin['called']:
                    return
                _fin['called'] = True
                if on_finish is None:
                    return
                try:
                    on_finish()
                except Exception as e:                    # noqa: BLE001
                    print(f'[tasks] on_finish 回调失败（任务 {tid}）: '
                          f'{type(e).__name__}: {e}', file=sys.stderr)
            # 排队等名额：每 0.5s 探一次，期间可被取消（不会被闸门永久卡住）
            got = False
            while not rec['cancel'].is_set():
                if slot.acquire(timeout=0.5):
                    if rec['cancel'].is_set():
                        # 竞态窗口：拿到名额的瞬间被取消——归还名额走
                        # 取消分支，不能带着名额开跑（2026-10-04 审查 P3）
                        slot.release()
                        break
                    got = True
                    break
                rec['msg'] = self._queue_msg(rec)
            if not got:                       # 排队途中被取消
                rec['status'] = 'cancelled'
                rec['error'] = '任务已停止（排队中取消）'
                rec['finished'] = time.time()
                _emit('[CANCEL] ' + cfg.tr('排队中已取消',
                                          'Cancelled while queued'))
                self._persist(rec)
                # fn 从未被调用 → 其 finally 不会执行；终态收尾必须在这里做
                _call_on_finish()
                return
            if rec['weight'] != 'light':
                with self.lock:
                    self._heavy_used += 1
            rec['status'] = 'running'
            rec['started'] = time.time()      # 已运行时长从真正开跑算起
            self._persist(rec)
            task_bind(rec['cancel'], rec['procs'])
            try:
                result = fn(_emit, _progress, rec['cancel'])
                # preview 组装可能整树 os.walk（大运行目录数千文件）：必须在
                # 锁外先算好再进锁。原来在锁内做，一次收尾会把 snapshot()、
                # _progress() 等所有轮询请求全部卡住。
                preview = build_result_preview(result, rec['name'])
                # 终态与结果必须一次性原子写入（snapshot 在锁内读）
                with self.lock:
                    rec['status'] = ('cancelled' if rec['cancel'].is_set()
                                     else 'done')
                    rec['result'] = result
                    rec['result_preview'] = preview
                if rec['status'] == 'done':
                    _emit('[OK] ' + cfg.tr('任务完成', 'Task finished'))
            except Exception as e:
                cancelled = rec['cancel'].is_set()
                # 带异常类型：只留裸 message 时 FileNotFoundError/KeyError
                # 这类信息量极低的失败无从排查（平台日志铁律）
                err = ('任务已停止（用户取消）' if cancelled
                       else f'{type(e).__name__}: {e}')
                with self.lock:
                    rec['status'] = 'cancelled' if cancelled else 'failed'
                    rec['error'] = err
                _emit(f"[ERROR] {err}")
                if not cancelled:
                    _log.debug('任务 %s 异常堆栈', rec.get('name'),
                               exc_info=True)
            finally:
                task_bind(None, None)
                with self.lock:
                    rec['procs'].clear()
                    rec['finished'] = time.time()
                    rec['msg'] = ''
                if rec['weight'] != 'light':
                    with self.lock:
                        self._heavy_used = max(0, self._heavy_used - 1)
                slot.release()          # 归还名额
                self._trim_memory()     # 回收旧任务内存
                self._persist(rec)
                _call_on_finish()       # 终态收尾（与排队取消分支恰一次互斥）

        th = threading.Thread(target=_run, daemon=True, name=f'task-{tid}')
        rec['thread'] = th
        th.start()
        self._persist(rec)
        return tid

    def _persist(self, rec):
        # 外部任务不落盘：它们的真身在别的子系统（如 downloads/<批次>/
        # batch.json），再写一份会被 _recover_interrupted 当成幽灵任务。
        if rec.get('external'):
            return
        try:
            p = check_path(os.path.join(DIRS['tasks'], rec['id'] + '.json'),
                           must_exist=False, in_platform=True)
            data = {'id': rec['id'], 'name': rec['name'],
                    'status': rec['status'], 'stage': rec['stage'],
                    'pct': rec['pct'], 'msg': rec['msg'],
                    'log_file': rec.get('log_file'),
                    'error': rec['error'],
                    'log_tail': list(rec['log'])[-40:],
                    'started': rec['started'], 'finished': rec['finished'],
                    'eta': rec.get('eta'), 'link': rec.get('link'),
                    'out_dir': rec.get('out_dir')}
            if rec.get('result') is not None:
                try:
                    json.dumps(rec['result'], ensure_ascii=False)
                    data['result'] = rec['result']
                except (TypeError, ValueError):
                    data['result'] = str(rec['result'])
            with safe_open(p, 'wt') as f:
                json.dump(data, f, ensure_ascii=False)
        except Exception as e:
            # 状态落盘失败不影响内存态，但会让刷新后丢进度：留痕便于排查
            _log.warning('任务状态保存失败 %s: %s', rec.get('id'), e)

    def _queue_msg(self, rec):
        """排队提示：统计同档位里排在自己前面的任务数。"""
        ahead = 0
        with self.lock:
            for tid in self.order:
                if tid == rec['id']:
                    break
                r = self.tasks.get(tid)
                if (r and r['status'] == 'queued'
                        and r.get('weight') == rec.get('weight')):
                    ahead += 1
        if ahead:
            return cfg.tr(f'排队中 · 前面还有 {ahead} 个任务',
                          f'Queued · {ahead} task(s) ahead')
        return cfg.tr('排队中 · 等待空闲名额', 'Queued · waiting for a slot')

    def _trim_memory(self):
        """回收内存：超出 _KEEP_FULL 的终态任务清掉日志与结果引用（磁盘已
        持久化，仍可经 tasks/*.json 回溯）；超出 _KEEP_TOTAL 的移出内存。

        服务长跑时任务只增不减会持续吃内存，这里做稳态回收。
        外部任务（external）不参与淘汰：它们不落盘、数量有限（下载批次/
        扫描等），被淘汰后卡片会凭空消失却又被下次 attach 复活。
        """
        with self.lock:
            for tid in self.order[_KEEP_FULL:]:
                rec = self.tasks.get(tid)
                if (not rec or rec['status'] in ('running', 'queued')
                        or rec.get('external')):
                    continue
                if rec['log']:
                    rec['log'].clear()
                rec['result'] = None
                rec['result_preview'] = None
                rec['thread'] = None
                rec['procs'] = []
            if len(self.order) > _KEEP_TOTAL:
                keep = [t for t in self.order[_KEEP_TOTAL:]
                        if (self.tasks.get(t) or {}).get('external')]
                for tid in self.order[_KEEP_TOTAL:]:
                    if tid in keep:
                        continue
                    self.tasks.pop(tid, None)
                self.order = self.order[:_KEEP_TOTAL] + keep

    # ------------------------------------------------------------------
    # 外部任务：执行归别的子系统（DownloadManager / SampleQueue / 存储扫描），
    # 本管理器只做统一展示与「停止」转发。
    #
    # 为什么不反过来让那些子系统 import 本模块：它们跑在自建线程里，
    # 反向依赖 web 层既乱分层、又容易碰到 Flask 上下文问题。这里用
    # provider 回调把状态「拉」过来，是单向依赖。
    # ------------------------------------------------------------------
    def attach(self, name, provider, on_cancel=None, link=None, key=None,
               weight='light'):
        """登记一个「外部任务」，返回其稳定 task id。

        provider() -> dict | None，每次快照时调用（**在锁外**，避免与
        provider 内部可能的重入互锁）：
            {'status': 'running'|'done'|'failed'|'cancelled',
             'pct': float|None, 'msg': str, 'stage': str, 'error': str,
             'eta': int|None, 'log': [str], 'started': float|None,
             'finished': float|None}
          返回 None 表示该任务已消失（卡片随之移除）。

        key: 稳定标识（同一 key 重复 attach 只登记一次），缺省用 name。
        外部任务**不占** heavy/light 并发名额、**不落** tasks/*.json
        （否则重启后会被 _recover_interrupted 当成幽灵任务标失败）。
        """
        key = str(key or name)
        with self.lock:
            tid = self._ext_keys.get(key)
            if tid and tid in self.tasks:
                # 重复 attach 只登记一次（tid 不变、started 不重置），但要把
                # provider/on_cancel 换成最新的闭包 —— 否则绑定方重建了执行
                # 对象（如 DownloadManager 单例被重建）后，卡片仍指向旧对象，
                # 快照会一直返回 None 导致卡片被误判成「已消失」。
                rec = self.tasks[tid]
                rec['provider'] = provider
                if on_cancel is not None:
                    rec['on_cancel'] = on_cancel
                if link:
                    rec['link'] = link
                if name:
                    rec['name'] = name
                return tid
            # 用 key 的哈希做稳定 id：形状与 uuid4().hex[:12] 一致（12 位
            # 十六进制），_tid_ok / 前端 data-tid 都能直接用。
            tid = hashlib.sha1(key.encode('utf-8')).hexdigest()[:12]
            if tid in self.tasks:
                return tid
            rec = {'id': tid, 'name': name, 'status': 'running', 'stage': '',
                   'pct': None, 'msg': '', 'error': '',
                   'log': deque(maxlen=500), 'started': time.time(),
                   'finished': None, 'result': None, 'result_preview': None,
                   'eta': None, 'cancel': threading.Event(), 'procs': [],
                   'weight': weight if weight in ('heavy', 'light') else 'light',
                   'log_file': None, 'link': link, 'restart': None,
                   'external': True, 'provider': provider,
                   'on_cancel': on_cancel, 'hidden': False}
            self.tasks[tid] = rec
            self.order.insert(0, tid)
            self._ext_keys[key] = tid
        return tid

    def detach(self, key):
        """移除某个外部任务的登记（provider 返回 None 时由内部自动清理）。"""
        with self.lock:
            tid = self._ext_keys.pop(str(key), None)
            if tid:
                self.tasks.pop(tid, None)
                if tid in self.order:
                    self.order.remove(tid)
        return bool(tid)

    def _provider_state(self, tid):
        """在锁外调用 provider，取回外部任务状态（异常按「消失」处理）。"""
        return self._provider_state_ex(tid)[0]

    def _provider_state_ex(self, tid):
        """同 _provider_state，但额外告诉调用方「是干净地没有了，还是读失败」。

        返回 (state, ok)：provider 返回 None → (None, True)；
        provider 抛异常 → (None, False)。前端展示要区分这两者：
        前者该移除卡片，后者只该降级显示（否则一次瞬时错误就让用户的
        卡片永久消失）。
        """
        with self.lock:
            rec = self.tasks.get(tid)
            if not rec or not rec.get('external'):
                return None, False
            prov = rec.get('provider')
        if prov is None:
            return None, False
        try:
            return prov(), True
        except Exception:
            _log.warning('外部任务状态读取失败 %s', tid, exc_info=True)
            return None, False

    def _prune_external(self):
        """清掉已隐藏 / 真身已消失的外部任务（卡片与内存都不堆积）。

        区分两种 None：provider **干净地**返回 None = 真身没了（批次被删）
        → 丢弃；provider 抛异常 = 只是这次读不到 → **保留**卡片，否则一个
        瞬时错误会让用户的卡片永久消失。
        """
        with self.lock:
            hidden = [t for t, r in self.tasks.items()
                      if r.get('external') and r.get('hidden')]
            ext = [t for t, r in self.tasks.items() if r.get('external')]
        for tid in hidden:
            self._drop_external(tid)
        for tid in ext:
            if self._provider_missing(tid):
                self._drop_external(tid)

    def _provider_missing(self, tid):
        """provider 是否明确表示「该任务已不存在」（异常不算）。"""
        with self.lock:
            rec = self.tasks.get(tid)
            prov = rec.get('provider') if (rec and rec.get('external')) else None
        if prov is None:
            return False
        try:
            return prov() is None
        except Exception:
            return False

    def _drop_external(self, tid):
        with self.lock:
            self.tasks.pop(tid, None)
            if tid in self.order:
                self.order.remove(tid)
            for k, v in list(self._ext_keys.items()):
                if v == tid:
                    self._ext_keys.pop(k, None)

    @staticmethod
    def _ext_snapshot(tid, name, link, st, log_lines):
        """把 provider 的字段映射成本管理器统一的快照结构。"""
        status = str(st.get('status') or 'running')
        pct = st.get('pct')
        try:
            pct = None if pct is None else max(0.0, min(float(pct), 1.0))
        except (TypeError, ValueError):
            pct = None
        log = st.get('log') or []
        return {'id': tid, 'name': name, 'status': status,
                'stage': str(st.get('stage') or ''), 'pct': pct,
                'msg': str(st.get('msg') or ''),
                'error': str(st.get('error') or ''),
                'started': st.get('started'), 'finished': st.get('finished'),
                'eta': st.get('eta'), 'link': st.get('link') or link,
                'restart': None,
                'log': [str(x) for x in log][-log_lines:] if log_lines else []}

    def snapshot(self, tid, log_lines=80, include_result=True):
        with self.lock:
            rec = self.tasks.get(tid)
            if not rec:
                return None
            if rec.get('external'):
                name, link = rec['name'], rec.get('link')
                started0 = rec.get('started')
            else:
                name = None
        if name is not None:
            st, ok = self._provider_state_ex(tid)
            if st is None:
                if ok:
                    return None          # provider 明确说：真身已不存在
                # 读失败：给降级卡片，别让用户的卡片凭空消失。
                # started 取本记录的挂载时间 —— 否则前端算「耗时」会拿
                # 一个空 started 当 0，显示成 49 万小时这种荒唐值。
                st = {'status': 'running', 'pct': None, 'stage': '',
                      'msg': cfg.tr('状态读取失败（会自动恢复）',
                                    'State read failed (recovers automatically)'),
                      'error': '', 'log': [], 'started': started0}
            return self._ext_snapshot(tid, name, link, st, log_lines)
        with self.lock:
            rec = self.tasks.get(tid)
            if not rec:
                return None
            # 'queued' 对外一律呈现为 running：前端据此显示取消按钮、进度条
            # 与自动滚动；排队状态本身通过 msg 文案告知用户。
            st = 'running' if rec['status'] == 'queued' else rec['status']
            out = {'id': rec['id'], 'name': rec['name'],
                   'status': st, 'stage': rec['stage'],
                   'pct': rec['pct'], 'msg': rec['msg'],
                   'error': rec['error'], 'started': rec['started'],
                   'finished': rec['finished'], 'eta': rec.get('eta'),
                   'link': rec.get('link'), 'restart': rec.get('restart'),
                   'out_dir': rec.get('out_dir'),
                   'log': list(rec['log'])[-log_lines:]}
            if include_result:
                out['result'] = rec.get('result_preview')
        return out

    def list_all(self, limit=_KEEP_FULL * 2, logs=True):
        """任务列表。

        前端每 2.5s 轮询一次本接口，历史上对所有任务都返回 30 行日志 +
        结果预览，任务累积后 payload 持续变大。这里做瘦身：
        - 运行中/排队中的任务：完整 30 行日志（实时看进度）
        - 已结束的任务：只带最后 5 行（日志已落盘，展开任务卡可单独拉取）
        - logs=False：完全不带日志（任务中心/徒轮询只需元信息，
          日志走 /api/task/<tid>/log 单独拉取）
        - 结果预览一律保留（任务完成后的结果面板依赖它）

        外部任务（external）一并列出：provider 返回 None（真身已删除）
        或已 hidden 的先剪掉，避免卡片凭空堆积。
        """
        self._prune_external()
        self._refresh_ext_status()      # 外部任务的活跃态由 provider 决定
        with self.lock:
            ids = list(self.order[:limit])
            active = {i for i in ids
                      if (self.tasks.get(i) or {}).get('status')
                      in ('running', 'queued')}
        if not logs:
            outs = [self.snapshot(i, log_lines=0, include_result=True)
                    for i in ids if i in self.tasks]
        else:
            outs = [self.snapshot(i, log_lines=(30 if i in active else 5),
                                  include_result=True)
                    for i in ids if i in self.tasks]
        # 外部任务的真身可能在两次取快照之间被删除（provider 变 None）：
        # 必须滤掉，否则 /api/tasks 会吐出 [null]，前端 list.filter(...)
        # 直接 TypeError。
        return [x for x in outs if x]

    @staticmethod
    def _kill_tree(proc):
        """杀进程树（Windows 用 taskkill /T 连子进程；其它平台 kill）。"""
        if proc.poll() is not None:
            return
        try:
            if sys.platform == 'win32':
                subprocess.run(['taskkill', '/F', '/T', '/PID', str(proc.pid)],
                               capture_output=True, timeout=15)
            else:
                proc.kill()
        except Exception:
            try:
                proc.kill()
            except Exception:
                pass

    def cancel(self, tid):
        with self.lock:
            rec = self.tasks.get(tid)
            ext_cancel = rec.get('on_cancel') if (rec or {}).get('external') \
                else None
        if ext_cancel is not None:
            # 外部任务的停止由真身负责（如 DownloadManager.cancel 会
            # terminate aria2c/sracha 子进程并置批次为 cancelled）。
            try:
                ok = ext_cancel() is not False
            except Exception:
                _log.warning('外部任务取消失败 %s', tid, exc_info=True)
                ok = False
            if ok:
                rec['msg'] = cfg.tr('已请求停止', 'Stop requested')
            return bool(ok)
        if rec and rec['status'] in ('running', 'queued'):
            rec['cancel'].set()
            if rec['status'] == 'queued':
                # 排队中：闸门等待循环会自行退出并归还名额
                rec['msg'] = '已取消排队'
                return True
            # 硬停止：杀掉本任务已注册的全部子进程（kunpeng/spades/mafft 等）
            for proc in list(rec.get('procs') or ()):
                self._kill_tree(proc)
            rec['msg'] = '已请求取消（当前步骤结束后停止）'
            return True
        return False


    def _refresh_ext_status(self):
        """把外部任务的活跃态刷新一遍（provider 的结论写回 rec['status']）。

        active_ids() / cancel_all() / 全局重置的「还有几个在收尾」都只看
        rec['status']，而外部任务的执行方是别人（下载线程、队列、扫描线程）。
        不在读取前刷新的话，一个早就跑完的下载批次会永远被当成「运行中」：
        全局重置会对它转发一次取消（把已完成的批次误标成「已取消」），
        收尾等待也永远等不到干净。provider 一律在锁外调用，避免互锁。
        """
        with self.lock:
            ext = [t for t, r in self.tasks.items() if r.get('external')]
        for tid in ext:
            st, ok = self._provider_state_ex(tid)
            if st is None:
                continue          # 真身没了/读失败：交给 _prune_external 处置
            with self.lock:
                rec = self.tasks.get(tid)
                if rec is not None:
                    rec['status'] = str(st.get('status') or 'running')

    def active_ids(self):
        """运行中/排队中的任务 id（全局重置与状态查询用）。"""
        self._refresh_ext_status()
        with self.lock:
            return [t for t, r in self.tasks.items()
                    if r['status'] in ('running', 'queued')]

    def cancel_all(self):
        """请求取消全部运行中/排队中的任务，返回被受理的条数。

        只发取消信号（含已注册子进程的硬停止），不删除任务记录，
        更不碰 results/ tool_runs/ 下任何结果文件。
        """
        n = 0
        for tid in self.active_ids():
            if self.cancel(tid):
                n += 1
        return n

    def purge_all(self, include_archive=True):
        """清空全部任务记录，返回 (removed, active_left)。

        全局重置用，边界写死：
        - 只动「任务状态」——内存记录 + tasks/<tid>.json；
          results/ tool_runs/ downloads/ submissions/ 等结果文件一律不碰。
        - 仍在运行/排队中的任务不强删（会撞上正在写盘的句柄），
          只统计条数返回；调用方应先 cancel_all() 等它们收敛。
        - include_archive=False 只清内存记录，磁盘历史归档原样保留
          （所以这里**不能**复用 delete()——它会顺手 unlink tasks/<tid>.json）。
        - removed 是「被清掉的不同任务 id 数」，内存与磁盘同一 id 只计一次。
        """
        with self.lock:
            active = {t for t, r in self.tasks.items()
                      if r['status'] in ('running', 'queued')}
            mem = [t for t in self.tasks if t not in active]
            for t in mem:
                self.tasks.pop(t, None)
                if t in self.order:
                    self.order.remove(t)
        ids = set(mem)
        if include_archive:
            try:
                files = [f for f in os.listdir(DIRS['tasks'])
                         if f.endswith('.json')]
            except OSError:
                files = []
            for f in files:
                tid = f[:-5]
                # _tid_ok 顺带把 queue.json 这类非任务文件挡在外面
                if not self._tid_ok(tid) or tid in active:
                    continue
                try:
                    p = check_path(os.path.join(DIRS['tasks'], f),
                                   must_exist=False, in_platform=True)
                    if os.path.isfile(p):
                        os.remove(p)
                        ids.add(tid)
                except (OSError, ValueError):
                    continue
        return len(ids), len(active)

    @staticmethod
    def _tid_ok(tid):
        return bool(tid) and bool(re.fullmatch(r'[0-9a-f]{6,16}', str(tid)))

    def delete(self, tid):
        """删除任务记录（仅终态）：移出内存 + 删 tasks/<tid>.json。

        外部任务只「隐藏卡片」：真身（下载批次等）归各自子系统管理，
        在这里 pop 掉之后下次 attach 又会复活，反而像是删不掉。
        """
        if not self._tid_ok(tid):
            return False, '非法任务 ID'
        with self.lock:
            rec = self.tasks.get(tid)
            if rec and rec.get('external'):
                rec['hidden'] = True
                return True, 'ok'
            if rec and rec['status'] in ('running', 'queued'):
                return False, '运行中的任务请先停止再删除'
            self.tasks.pop(tid, None)
            if tid in self.order:
                self.order.remove(tid)
        try:
            p = check_path(os.path.join(DIRS['tasks'], tid + '.json'),
                           must_exist=False, in_platform=True)
            if os.path.isfile(p):
                os.remove(p)
        except (OSError, ValueError):
            pass
        return True, 'ok'

    def clear_finished(self):
        """清空内存中的终态任务（含各自的 tasks/<tid>.json）。

        注意：不碰磁盘上的历史归档（平台重启后加载的 tasks/*.json）。
        清空归档是破坏性操作，只能由用户逐个删除。

        外部任务按「隐藏」处理（见 delete），并顺带剪掉 provider 已消失的。
        """
        self._prune_external()
        with self.lock:
            tids = [t for t, r in self.tasks.items()
                    if r['status'] not in ('running', 'queued')
                    or r.get('external')]
            # 外部任务的 rec.status 只是登记初值，一律按「可清」处理：
            # 真身是否还在由 provider 说话，隐藏掉即可。
            tids = [t for t in tids
                    if (self.tasks.get(t) or {}).get('external')
                    or (self.tasks.get(t) or {}).get('status')
                    not in ('running', 'queued')]
        n = 0
        for t in tids:
            ok, _ = self.delete(t)
            if ok:
                n += 1
        self._prune_external()
        return n

    def list_archived(self, limit=200):
        """归档任务（磁盘 tasks/*.json，不在内存）：平台重启后的历史回溯。"""
        items, total = self._read_archived(limit=limit, logs=True)
        return items

    def count_archived(self):
        """磁盘归档任务总数（不受 list_archived 的 limit 截断影响）。"""
        try:
            return len([f for f in os.listdir(DIRS['tasks'])
                        if f.endswith('.json') and self._tid_ok(f[:-5])])
        except OSError:
            return 0

    def _read_archived(self, limit=200, logs=True):
        items = []
        try:
            files = [f for f in os.listdir(DIRS['tasks']) if f.endswith('.json')]
        except OSError:
            return items, 0
        total = 0
        for f in files:
            tid = f[:-5]
            if tid in self.tasks or not self._tid_ok(tid):
                continue
            total += 1
            if len(items) >= limit:
                continue
            try:
                with safe_open(os.path.join(DIRS['tasks'], f)) as fh:
                    d = json.load(fh)
                if not logs:
                    d.pop('log_tail', None)
                items.append(d)
            except (OSError, ValueError):
                continue
        items.sort(key=lambda d: d.get('started') or 0, reverse=True)
        return items[:limit], total

    def full_log(self, tid, lines=400):
        """任务全量日志（尾部 lines 行）：优先内存，内存被回收/归档时读
        磁盘 log_file（json 里记录了路径），最后回退 json 内的 log_tail。"""
        if not self._tid_ok(tid):
            return None
        with self.lock:
            rec = self.tasks.get(tid)
            is_ext = bool(rec and rec.get('external'))
        if is_ext:
            # 外部任务没有内存 deque，日志由 provider 提供（如批次
            # batch.log 的尾部）。
            st = self._provider_state(tid) or {}
            return [str(x) for x in (st.get('log') or [])][-lines:]
        rec = self.tasks.get(tid)
        if rec and rec['log']:
            return list(rec['log'])[-lines:]
        log_file = (rec or {}).get('log_file')
        js = None
        if not log_file:
            try:
                with safe_open(os.path.join(DIRS['tasks'], tid + '.json')) as f:
                    js = json.load(f)
                log_file = js.get('log_file')
            except (OSError, ValueError):
                js = None
        if log_file and os.path.isfile(log_file):
            # 防御深度：归档 json 的 log_file 必须在平台目录内
            try:
                check_path(log_file, must_exist=True, in_platform=True)
            except (OSError, ValueError):
                log_file = None
        if log_file and os.path.isfile(log_file):
            tail = self._tail_log_file(log_file, lines)
            if tail is not None:
                return tail
        if js is not None:
            return (js.get('log_tail') or [])[-lines:]
        try:
            with safe_open(os.path.join(DIRS['tasks'], tid + '.json')) as f:
                return (json.load(f).get('log_tail') or [])[-lines:]
        except (OSError, ValueError):
            return None

    _logcache = {}   # {(path, size, mtime_ns, lines): [行...]} 日志尾部缓存

    @classmethod
    def _tail_log_file(cls, path, lines):
        """从日志文件尾部读 lines 行。

        长任务日志可达 MB 级，而任务中心展开时会 2s 轮询一次，
        逐行读完整个文件会持续重复磁盘 IO。这里 seek 到尾部只读
        必要字节，并以 (size, mtime_ns, lines) 作键缓存：文件未变
        时直接复用，轮询场景下完全命中。
        """
        try:
            st = os.stat(path)
        except OSError:
            return None
        key = (path, st.st_size, st.st_mtime_ns, lines)
        hit = cls._logcache.get(key)
        if hit is not None:
            return list(hit)
        try:
            size = st.st_size
            chunk = min(size, max(65536, lines * 256))
            with safe_open(path, 'rb') as f:
                if size > chunk:
                    f.seek(size - chunk)
                data = f.read()
        except OSError:
            return None
        # utf-8 优先、GBK 兜底：历史日志（2026-09-27 编码统一前）里子进程
        # 按系统 ANSI 写的中文段按 utf-8+replace 会整体变 U+FFFD，读不回来
        rows = decode_output(data).splitlines()
        if size > chunk and rows:
            rows = rows[1:]           # 首行可能被字节截断，丢弃
        tail = rows[-lines:]
        if len(cls._logcache) > 32:
            cls._logcache.pop(next(iter(cls._logcache)))
        cls._logcache[key] = tuple(tail)
        return tail


tm = TaskManager()


def _sample_result_preview(sample):
    """样品任务的完成预览：报告链接 + 关键数字 + 产物文件。"""
    from Virus_Platform_Core.pipeline import pipeline_overview
    from Virus_Platform_Core.utils import resolve_sample_name
    # sample 来自 os.path.relpath(结果目录, results/)，本身就是真实目录名；
    # 用 resolve 而不是字符转换，避免手工/历史命名的样品被改名后 404。
    s = resolve_sample_name(sample, DIRS['results'])
    try:
        sd = check_path(os.path.join(DIRS['results'], s), must_exist=True,
                        in_platform=True)
    except (ValueError, FileNotFoundError):
        return {'kind': 'sample', 'sample': s}
    ov = pipeline_overview(sd)
    stages = [{'stage': x['stage'], 'name': x['name'], 'status': x['status'],
               'summary': x.get('summary') or ''} for x in ov.get('stages', [])]
    rpt = os.path.join(sd, '07_report', 'report.html')
    return {'kind': 'sample', 'sample': s,
            'report': f'/report/{s}/' if os.path.isfile(rpt) else None,
            'out_dir': sd,
            'stages': stages}


def _tool_result_files(run_name, limit=12):
    """工具运行目录的产物清单（按修改时间倒序，取前 limit 个）。"""
    root = check_path(os.path.join(_tool_runs_root(), run_name),
                      must_exist=False, in_platform=True)
    if not os.path.isdir(root):
        return []
    items = []
    for cur, _sub, fns in os.walk(root):
        for fn in fns:
            p = os.path.join(cur, fn)
            try:
                items.append((os.path.getmtime(p), os.path.relpath(p, root),
                              os.path.getsize(p)))
            except OSError:
                continue
    items.sort(reverse=True)
    return [{'path': rel.replace(os.sep, '/'), 'size': fmt_size(sz)}
            for _mt, rel, sz in items[:limit]]


def _flatten_stats(result, depth=2):
    """把任务返回 dict 展平成单层统计（嵌套 dict 用点号连接键名）。

    旧实现只取顶层标量，ORF 任务的 pyrodigal/orfipy/orfa 嵌套统计全被丢掉
    （与 04_orf/summary.json 的 n_orfs 等字段对不上）。列表只记条数
    （<key>_n），非有限浮点（NaN/inf）剔除——NaN 不能进 JSON 响应体。
    """
    import math
    out = {}

    def _ok(v):
        if isinstance(v, bool) or v is None or isinstance(v, (str, int)):
            return True
        return isinstance(v, float) and math.isfinite(v)

    def _walk(node, prefix, depth_left):
        for k, v in node.items():
            key = f'{prefix}{k}'
            if _ok(v):
                out[key] = v
            elif isinstance(v, dict) and depth_left > 0:
                _walk(v, key + '.', depth_left - 1)
            elif isinstance(v, (list, tuple)):
                out[key + '_n'] = len(v)
    _walk(result, '', depth)
    return out


def build_result_preview(result, task_name=''):
    """把任务返回值规范化成前端可渲染的预览 dict（失败安全）。"""
    try:
        if isinstance(result, str) and result:
            # 样品任务：返回 results/<样品> 目录
            rel = os.path.relpath(result, DIRS['results'])
            if rel.startswith('..') or os.path.isabs(rel):
                return {'kind': 'raw', 'text': str(result)}
            return _sample_result_preview(rel)
        if isinstance(result, dict):
            out = {'kind': 'tool', 'stats': _flatten_stats(result)}
            run = result.get('run')
            if run:
                out['run'] = str(run)          # 前端 toolrun 需顶层 run 名
                out['files'] = _tool_result_files(str(run))
            out_dir = result.get('out_dir')
            if isinstance(out_dir, str) and out_dir:
                out['out_dir'] = out_dir       # 输出目录随结果下发，前端明示
            return out
    except Exception:
        pass
    return None
