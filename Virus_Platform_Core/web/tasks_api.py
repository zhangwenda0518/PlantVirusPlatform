# -*- coding: utf-8 -*-
"""任务状态 API（自 app.py 拆出）。

/api/tasks 列表与详情、日志、SSE 流、取消/删除/清理。"""
import json
import time

from flask import (Blueprint, Response, abort, jsonify, request)

from Virus_Platform_Core.web.tasks import tm

bp = Blueprint('tasks_api', __name__)

# 破坏性操作确认令牌（防误触，非鉴权）。改动此值须同步：
# webapp/static/app.js 的 gxGlobalReset()、tests/_it_global_reset.py、
# tests/_it_global_ctx_ui.py 的打桩响应、tests/_audit_routes.py 的跳过表。
RESET_CONFIRM = 'reset'


@bp.route('/api/tasks')
def api_tasks():
    # logs=0：只要元信息不要日志（任务中心徒轮询用，日志走单独端点）
    logs = request.args.get('logs') not in ('0', 'false', 'no')
    if request.args.get('full'):
        arch, arch_total = tm._read_archived(logs=logs)
        return jsonify({'active': tm.list_all(logs=logs),
                        'archived': arch,
                        'archived_total': arch_total})
    return jsonify(tm.list_all(logs=logs))


@bp.route('/api/task/<tid>')
def api_task(tid):
    try:
        lines = min(int(request.args.get('log_lines') or 80), 1000)
    except (TypeError, ValueError):
        lines = 80
    snap = tm.snapshot(tid, log_lines=lines)
    if not snap:
        abort(404)
    return jsonify(snap)


@bp.route('/api/task/<tid>/log')
def api_task_log(tid):
    """任务全量日志（尾部 N 行）：内存任务 / 归档任务统一入口。"""
    try:
        lines = min(int(request.args.get('lines') or 400), 2000)
    except (TypeError, ValueError):
        lines = 400
    log = tm.full_log(tid, lines=lines)
    if log is None:
        abort(404, '日志不存在')
    return jsonify({'log': log})


@bp.route('/api/task/<tid>/stream')
def api_task_stream(tid):
    """SSE 实时任务进度流：快照有变化立即推送，任务结束自动关流。

    前端 EventSource 订阅（LOGAN 批量面板用），替代定时轮询。

    多路复用：`?also=<tid2>,<tid3>`（逗号分隔，单连接最多 8 个任务）
    在一次连接里同时推送多个任务的快照，每条 data 事件带 tid 字段标明
    来源；所有任务都到终态（或消失）才关流。HTTP/1.1 浏览器对同源
    SSE 有 6 连接上限，批量任务逐条开流会触顶把轮询也堵死——前端将其
    合并成一条连接即可（用法：/api/task/<主tid>/stream?also=a,b,c）。
    不传 also 时事件负载与旧版逐字节一致（单任务订阅无需改动）。
    """
    tids = [tid]
    for extra in (request.args.get('also') or '').split(','):
        extra = extra.strip()
        if extra and extra not in tids:
            tids.append(extra)
    tids = tids[:8]                       # 单连接上限，防被拿来放大轮询
    multi = len(tids) > 1

    def gen():
        last = {}
        idle = 0
        done = set()                      # 已到终态 / 已消失的任务
        while len(done) < len(tids):
            for t in tids:
                if t in done:
                    continue
                snap = tm.snapshot(t, log_lines=10)
                if snap is None:
                    done.add(t)
                    gone = json.dumps({'tid': t}, ensure_ascii=False) \
                        if multi else '{}'
                    yield f'event: gone\ndata: {gone}\n\n'
                    continue
                item = dict(snap, tid=t) if multi else snap
                data = json.dumps(item, ensure_ascii=False)
                if data != last.get(t):
                    last[t] = data
                    idle = 0
                    yield f'data: {data}\n\n'
                if snap.get('status') != 'running':
                    done.add(t)
            if len(done) >= len(tids):
                break
            idle += 1
            if idle >= 30:                     # ~20s 心跳注释行
                idle = 0
                yield ': ping\n\n'
            time.sleep(0.7)

    return Response(gen(), mimetype='text/event-stream',
                    headers={'Cache-Control': 'no-cache',
                             'X-Accel-Buffering': 'no'})


@bp.route('/api/task/<tid>/cancel', methods=['POST'])
def api_task_cancel(tid):
    return jsonify({'ok': tm.cancel(tid)})


@bp.route('/api/task/<tid>/delete', methods=['POST'])
def api_task_delete(tid):
    ok, msg = tm.delete(tid)
    if not ok:
        abort(400, msg)
    return jsonify({'ok': True})


@bp.route('/api/tasks/clear_finished', methods=['POST'])
def api_tasks_clear_finished():
    return jsonify({'removed': tm.clear_finished()})


@bp.route('/api/tasks/stop_all', methods=['POST'])
def api_tasks_stop_all():
    """一键停止全部运行中/排队中的任务（不删任务记录、不碰结果文件，
    与 /api/global/reset 的破坏性清理划清界限）。active_left：取消是异步
    的，仍有任务在收尾时前端据此提示"N 个任务正在停止"。"""
    cancelled = tm.cancel_all()
    deadline = time.time() + 2.0
    while time.time() < deadline and tm.active_ids():
        time.sleep(0.1)
    return jsonify({'ok': True, 'cancelled': cancelled,
                    'active_left': len(tm.active_ids())})


@bp.route('/api/global/reset', methods=['POST'])
def api_global_reset():
    """全局重置：一次收拢所有模块的运行态与相互关联（**不删任何结果文件**）。

    为什么需要：样品、专项工具、队列、任务记录分散在 results/ tool_runs/
    tasks/ queue.json 五处，切模块时没有单一入口能把"还在跑的 + 排队的 +
    互相交接的"一次性归零。这里就是那个总闸。

    做三件事：
      1. `tm.cancel_all()`  取消全部运行中/排队中的任务（含子进程硬停止）
      2. `queue.reset_all()` 清空批处理队列全部条目（含未启动的）
      3. `tm.purge_all()`   清空任务记录（内存 + tasks/*.json 归档）

    **明确不碰**：results/（样品与全部产物）、tool_runs/（专项运行与历史）、
    downloads/、submissions/、gb_collections/、meta_search/、logan/。

    请求体：
      confirm         必填，必须等于 'reset'。**不是鉴权**（本地平台没有鉴权
                      概念），而是"防误触"：本接口在空请求体下就会做破坏性
                      清理，任何泛化的调用方——最典型的是
                      tests/_audit_routes.py 这种对每条路由发一次空 POST 的
                      健壮性巡检——都会顺手把用户的全部任务记录清空。
                      （2026-09-11 实测踩过：巡检确实清掉了 run/tasks/*.json。）
      include_archive 可选，默认 true。false 时只清内存记录、保留磁盘归档。

    返回各步计数 + `active_left`：取消是异步的，仍有任务在收尾时前端应提示
    "还有 N 个任务正在停止"，而不是假装已经干净。
    """
    import time as _time

    from Virus_Platform_Core.web import state as _state
    body = request.get_json(silent=True) or {}
    if str(body.get('confirm') or '') != RESET_CONFIRM:
        abort(400, '破坏性操作需显式确认：请求体须含 '
                   f'{{"confirm": "{RESET_CONFIRM}"}}')
    include_archive = bool(body.get('include_archive', True))

    cancelled = tm.cancel_all()
    queue = _state.get_sample_queue()
    cleared_queue = queue.reset_all() if queue is not None else 0

    # 给已受理的取消一点收敛时间：多数阶段在 check_cancel 处即刻退出，
    # 短暂等待能显著提高"一次点完就干净"的比例；等不到也不阻塞，
    # 把剩余条数如实报给前端。
    deadline = _time.time() + 2.0
    while _time.time() < deadline and tm.active_ids():
        _time.sleep(0.1)

    removed, active_left = tm.purge_all(include_archive=include_archive)
    return jsonify({
        'ok': True,
        'cancelled_tasks': cancelled,
        'cleared_queue': cleared_queue,
        'removed_task_records': removed,
        'active_left': active_left,
        'kept_results': True,
    })
