# -*- coding: utf-8 -*-
"""AI 代理审计落盘：执行单状态变迁与事件流 → 追加式 JSONL。

**为什么需要**：`executors._actions` 与 `broadcast._events` 都是**进程内内存**
——重启即失，事后无从追溯"AI 提了什么动作、谁批的、什么时候执行、结果如何"。
自动动作必须有可审计的痕迹，否则"人在回路"只剩一个界面承诺。

**产物**：`DIRS['ai_audit']/actions-YYYY-MM-DD.jsonl`、`events-YYYY-MM-DD.jsonl`
（按天轮转便于清理与按日检索；只追加、不修改）。

**三条硬约束**（对齐 flask-async-write-drain 的实测结论）
  1. **不阻塞调用方**：`record()` 只入队，后台线程批量落盘。事件流在流水线
     线程里发布，绝不能因为写盘把分析任务拖慢。
  2. **退出不丢**：Windows 下 `SIGTERM` **不触发 atexit**（实测），只靠 atexit
     会丢队列里最后几秒的记录 ⇒ 另装 SIGINT/SIGTERM/SIGBREAK 处理器，先排空
     再按默认语义退出。
  3. **不重复写**：atexit 与后台线程会竞态取同一批 ⇒ `_STOP` 用 `Event`、
     `_drain()` 严格 `set() → join() → 排空 → 写`，且**写盘逻辑只有
     `_write_batch` 一处**，两个调用方共用。

**铁律**
  - 落盘失败**绝不抛异常**（审计挂了不能影响分析任务）；
  - 敏感键（`api_key` / `*_token` / `*_secret` / `*_password` …）**整键脱敏**；
  - 单条记录有深度/列表/文本长度上限，防一条事件撑爆文件。

**关闭**：环境变量 `VP_AI_AUDIT=0`（测试与"不想留盘"的场景用）。
"""
import atexit
import json
import os
import re
import signal
import threading
import time

from Virus_Platform_Core.config import DIRS

KINDS = ('actions', 'events')

RETENTION_DAYS = 90             # 按天轮转后的保留天数（超期文件自动清理）
_FLUSH_EVERY = 2.0              # 后台线程批量落盘周期（秒）
_MAX_BATCH = 200                # 单批上限（写事务别太长）
_MAX_TEXT = 2000                # 单个字符串字段上限
_MAX_LIST = 50                  # 单个列表字段上限
_MAX_DEPTH = 4                  # 嵌套深度上限
_TAIL_BYTES = 2 * 1024 * 1024   # 读回时单文件最多从尾部读多少字节

# 敏感键：整键脱敏（值不落盘）。用「精确名 + 后缀」两路，避免误伤
# 业务键（如 'key' 只按精确名匹配，'apikey_id' 这类不会被误杀）。
_SECRET_EXACT = {'key', 'token', 'secret', 'password', 'passwd',
                 'authorization', 'auth', 'credential', 'credentials'}
_SECRET_SUFFIX = ('_key', '_token', '_secret', '_password', '_passwd')

_FILE_RE = re.compile(r'^(actions|events)-(\d{4}-\d{2}-\d{2})\.jsonl$')

_Q_LOCK = threading.Lock()
_Q = []                          # [(kind, rec)]，旧 → 新
_STATE_LOCK = threading.Lock()
_enqueued = 0                    # 累计入队条数
_written = 0                     # 累计已写盘条数（测试等它追平 _enqueued）
_STOP = threading.Event()        # 退出用：置位后后台线程结束
_WAKE = threading.Event()        # 立刻写一轮用（不影响线程存活）
_WORKER = None
_STARTED = False
_ATEXIT_DONE = False
_START_LOCK = threading.Lock()

_SIGNALS_INSTALLED = False
_ORIG_HANDLERS = {}

_PRUNE_LOCK = threading.Lock()
_PRUNED = set()                  # 已清理过的 (root, day)


# ------------------------------------------------------------------
# 脱敏与截断
# ------------------------------------------------------------------
def _secret_key(k):
    kl = str(k).lower()
    return kl in _SECRET_EXACT or kl.endswith(_SECRET_SUFFIX)


def _scrub(v, depth=0):
    """递归脱敏 + 截断（输出保证可 JSON 序列化）。"""
    if depth > _MAX_DEPTH:
        return '…'
    if isinstance(v, dict):
        out = {}
        for k, vv in v.items():
            out[str(k)] = '••••' if _secret_key(k) else _scrub(vv, depth + 1)
        return out
    if isinstance(v, (list, tuple)):
        items = [_scrub(x, depth + 1) for x in list(v)[:_MAX_LIST]]
        if len(v) > _MAX_LIST:
            items.append(f'…(+{len(v) - _MAX_LIST})')
        return items
    if isinstance(v, str):
        return v if len(v) <= _MAX_TEXT else v[:_MAX_TEXT] + '…'
    if v is None or isinstance(v, (int, float, bool)):
        return v
    return str(v)[:_MAX_TEXT]


def _enabled():
    return os.environ.get('VP_AI_AUDIT', '').strip().lower() not in (
        '0', 'false', 'off', 'no')


# ------------------------------------------------------------------
# 落盘根目录
# ------------------------------------------------------------------
def _audit_root():
    """审计根目录（**运行时**读 DIRS —— `apply_output_root` 会原地改写它，
    所以这里每次现取，不能缓存成模块常量）。"""
    root = DIRS.get('ai_audit')
    if not root:
        return None
    try:
        os.makedirs(root, exist_ok=True)
    except OSError:
        return None
    return root


def _day_of(rec):
    """按记录自身的 ts 定日期（跨午夜的批次不会归错天）。"""
    try:
        return time.strftime('%Y-%m-%d',
                             time.localtime(float(rec.get('ts'))))
    except (TypeError, ValueError, OSError):
        return time.strftime('%Y-%m-%d')


# ------------------------------------------------------------------
# 写入（唯一入口）
# ------------------------------------------------------------------
def _write_batch(batch):
    """唯一的写盘入口（后台线程与 `_drain` 共用）。失败静默。"""
    global _written
    if not batch:
        return
    root = _audit_root()
    if not root:
        return
    groups = {}
    for kind, rec in batch:
        groups.setdefault((kind, _day_of(rec)), []).append(rec)
    done = 0
    for (kind, day), recs in groups.items():
        path = os.path.join(root, f'{kind}-{day}.jsonl')
        try:
            with open(path, 'a', encoding='utf-8', newline='\n') as f:
                for rec in recs:
                    f.write(json.dumps(rec, ensure_ascii=False,
                                       default=str) + '\n')
        except (OSError, TypeError, ValueError):
            continue
        done += len(recs)
        _maybe_prune(root, day)
    with _STATE_LOCK:
        _written += done


def _maybe_prune(root, day):
    """每个 (root, 当天) 只清理一次：删掉超过 RETENTION_DAYS 的旧文件。"""
    key = (root, day)
    with _PRUNE_LOCK:
        if key in _PRUNED:
            return
        _PRUNED.add(key)
    cutoff = time.time() - RETENTION_DAYS * 86400
    try:
        names = os.listdir(root)
    except OSError:
        return
    for name in names:
        m = _FILE_RE.match(name)
        if not m:
            continue
        try:
            t = time.mktime(time.strptime(m.group(2), '%Y-%m-%d'))
        except ValueError:
            continue
        if t < cutoff:
            try:
                os.remove(os.path.join(root, name))
            except OSError:
                pass


def _flush():
    """后台线程：**用 Event.wait 而非 time.sleep**，`_WAKE`/`_STOP` 都能立刻唤醒。"""
    while not _STOP.is_set():
        _WAKE.wait(_FLUSH_EVERY)
        _WAKE.clear()
        with _Q_LOCK:
            batch = _Q[:_MAX_BATCH]
            del _Q[:_MAX_BATCH]
        _write_batch(batch)


def _drain():
    """退出路径：先停后台线程 → join → 再排空（顺序反了会重复写）。"""
    _STOP.set()
    _WAKE.set()
    w = _WORKER
    if w is not None and w.is_alive():
        try:
            w.join(timeout=3.0)
        except Exception:
            pass
    with _Q_LOCK:
        batch = _Q[:]
        _Q.clear()
    _write_batch(batch)


def _ensure_started():
    """启动后台写线程 + 注册 atexit（幂等；**不**在这里装信号处理器）。

    后台线程若已死亡（例如 `_drain()` 已经跑过、或线程异常退出），这里会
    **重启一个** —— 否则后续所有 `record()` 都会静默堆在队列里永不落盘。
    """
    global _WORKER, _STARTED, _ATEXIT_DONE
    w = _WORKER
    if _STARTED and w is not None and w.is_alive():
        return
    with _START_LOCK:
        w = _WORKER
        if _STARTED and w is not None and w.is_alive():
            return
        _STOP.clear()
        _STARTED = True
        if not _ATEXIT_DONE:
            _ATEXIT_DONE = True
            try:
                atexit.register(_drain)
            except Exception:
                pass
        try:
            _WORKER = threading.Thread(target=_flush, daemon=True,
                                       name='ai-audit-writer')
            _WORKER.start()
        except Exception:
            _WORKER = None


# ------------------------------------------------------------------
# 公开 API
# ------------------------------------------------------------------
def record(kind, rec):
    """把一条记录排入审计队列。**永不抛异常**，绝不阻塞调用方。"""
    global _enqueued
    try:
        if not _enabled() or kind not in KINDS or not isinstance(rec, dict):
            return
        _ensure_started()
        item = (kind, _scrub(rec))
        with _Q_LOCK:
            _Q.append(item)
        with _STATE_LOCK:
            _enqueued += 1
    except Exception:
        pass


def record_many(items):
    """批量入队：items = [(kind, rec), ...]（读盘重建时一次性回灌用）。"""
    for kind, rec in items or ():
        record(kind, rec)


def record_sync(kind, rec):
    """**同步**落盘一条（仍走唯一的 `_write_batch`）。永不抛异常。

    为什么需要：异步队列最多会丢 `_FLUSH_EVERY`(2s) 内的记录，而 Windows 上
    `taskkill` / `TerminateProcess` 是**进程级不可捕获**的 —— 装信号处理器
    也救不回来。执行单的状态变迁（提案/批准/驳回/执行/失败）是低频高价值
    审计对象，量级是"每轮几个"，同步写的代价可忽略，但丢不起。
    高频的事件流仍走 `record()` 异步通道。
    """
    try:
        if not _enabled() or kind not in KINDS or not isinstance(rec, dict):
            return
        _write_batch([(kind, _scrub(rec))])
    except Exception:
        pass


def install_signals():
    """接管 SIGINT/SIGTERM/SIGBREAK：退出前先排空审计队列。

    ⚠ Windows 下 `SIGTERM` **不触发 atexit**（实测：rc=1 且回调一个都不跑），
    所以不能只靠 atexit。只在**主线程**装（`signal.signal` 有线程限制）。
    返回是否安装成功。故意不在 `record()` 里自动调用 —— 否则测试进程会被
    悄悄接管信号。
    """
    global _SIGNALS_INSTALLED
    if _SIGNALS_INSTALLED:
        return True
    if threading.current_thread() is not threading.main_thread():
        return False
    _ensure_started()
    ok = False
    for nm in ('SIGINT', 'SIGTERM', 'SIGBREAK'):   # SIGBREAK 是 Windows 专有
        sig = getattr(signal, nm, None)
        if sig is None:
            continue
        try:
            _ORIG_HANDLERS[sig] = signal.getsignal(sig)
            signal.signal(sig, _on_signal)
            ok = True
        except Exception:
            pass
    _SIGNALS_INSTALLED = ok
    return ok


def _on_signal(signum, frame):
    """收到终止信号：先落盘，再按默认方式退出（保住退出码语义）。"""
    try:
        _drain()
    except Exception:
        pass
    try:
        prev = _ORIG_HANDLERS.get(signum)
        signal.signal(signum, prev if prev is not None else signal.SIG_DFL)
        if signum == getattr(signal, 'SIGINT', None):
            raise KeyboardInterrupt
        os.kill(os.getpid(), signum)     # 重发 → 正常退出语义
    except Exception:
        os._exit(0)                      # 兜底：绝不能让进程卡死


# ------------------------------------------------------------------
# 读回（供 UI / CLI / 测试）
# ------------------------------------------------------------------
def days_of(kind):
    """某种记录已有的日期（新 → 旧）。"""
    root = _audit_root()
    if not root or kind not in KINDS:
        return []
    out = set()
    try:
        for name in os.listdir(root):
            m = _FILE_RE.match(name)
            if m and m.group(1) == kind:
                out.add(m.group(2))
    except OSError:
        return []
    return sorted(out, reverse=True)


def _read_file(path, max_bytes=_TAIL_BYTES):
    """读一个 JSONL 文件（超长时只读尾部 max_bytes）。坏行跳过，不炸。"""
    try:
        size = os.path.getsize(path)
    except OSError:
        return []
    try:
        with open(path, 'rb') as f:
            if size > max_bytes:
                f.seek(size - max_bytes)
                f.readline()             # 丢掉可能被截断的首行
            raw = f.read()
    except OSError:
        return []
    out = []
    for line in raw.decode('utf-8', errors='replace').splitlines():
        line = line.strip()
        if not line:
            continue
        try:
            out.append(json.loads(line))
        except ValueError:
            continue
    return out


def read(kind, limit=200, day=None):
    """读回审计记录（**新 → 旧**，最多 limit 条）。

    day 省略时从最近有记录的日期往前取；`limit<=0` 表示不截断。
    """
    if kind not in KINDS:
        return []
    root = _audit_root()
    if not root:
        return []
    days = [day] if day else days_of(kind)
    out = []
    for d in days:
        out.extend(_read_file(os.path.join(root, f'{kind}-{d}.jsonl')))
        if limit and limit > 0 and len(out) >= limit:
            break
    out.reverse()
    return out[:limit] if (limit and limit > 0) else out


def stats():
    """审计概况：每种记录的可用日期与字节数。"""
    root = _audit_root()
    out = {}
    for kind in KINDS:
        days = days_of(kind)
        total = 0
        for d in days:
            try:
                total += os.path.getsize(
                    os.path.join(root, f'{kind}-{d}.jsonl'))
            except OSError:
                pass
        out[kind] = {'n_days': len(days), 'days': days[:7], 'bytes': total,
                     'root': root or ''}
    return out


def pending():
    """队列里尚未被取走的条数。"""
    with _Q_LOCK:
        return len(_Q)


def counters():
    """(已入队, 已写盘) —— 测试用「写盘数追平入队数」判落盘完成。"""
    with _STATE_LOCK:
        return _enqueued, _written


def flush_now(timeout=5.0):
    """把队列排空并等落盘完成。返回是否在超时内完成。

    **不依赖后台线程存活** —— 自己动手排空（批次在锁下切分，与后台线程
    拿到的批次互不相交，所以不会重复写）。判据是**计数追平**
    （`_written >= _enqueued` 且队列空），而不是"等一个固定时长"——固定
    时长既慢又不确定（对齐 flask-async-write-drain 的「比对请求数 vs
    记录数」判法，不用时间戳分组）。
    """
    _ensure_started()
    _WAKE.set()
    deadline = time.time() + max(0.1, float(timeout))
    while time.time() < deadline:
        with _Q_LOCK:
            batch = _Q[:_MAX_BATCH]
            del _Q[:_MAX_BATCH]
        if batch:
            _write_batch(batch)
            continue
        enq, wr = counters()
        if wr >= enq:
            return True
        time.sleep(0.02)
    return False


def shutdown():
    """显式关停（测试收尾用）：停线程 + 排空。"""
    _drain()
