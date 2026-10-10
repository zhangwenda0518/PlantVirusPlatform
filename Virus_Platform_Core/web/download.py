# -*- coding: utf-8 -*-
"""公共数据下载中心（自 app.py 拆出）。

SRR/ERR/DRR/CRR 或 URL → aria2c 批量下载 → 一键转入分析流程。"""
import os
import re

from flask import (Blueprint, abort, jsonify, render_template, request,
                   send_file)

from Virus_Platform_Core.config import DIRS
from Virus_Platform_Core.utils import (check_path)
from Virus_Platform_Core.web import state as _state


def _queue():
    """批处理队列单例（app.py 启动时注册，避免循环导入）。"""
    return _state.get_sample_queue()


def _fmt_bytes(n):
    """字节数 → 人类可读（与下载页 fmtBytes 同口径，用于任务卡文案）。"""
    try:
        n = float(n or 0)
    except (TypeError, ValueError):
        return '—'
    if n <= 0:
        return '—'
    for div, unit in ((1 << 30, 'GB'), (1 << 20, 'MB'), (1 << 10, 'KB')):
        if n >= div:
            return f'{n / div:.1f} {unit}'
    return f'{int(n)} B'


def _fmt_speed(bps):
    try:
        bps = float(bps or 0)
    except (TypeError, ValueError):
        return '—'
    if bps < 1024:
        return '—'
    for div, unit in ((1 << 30, 'GB/s'), (1 << 20, 'MB/s'), (1 << 10, 'KB/s')):
        if bps >= div:
            return f'{bps / div:.1f} {unit}'
    return '—'


# 下载批次 → 任务中心「外部任务」的状态映射
_DL_STATUS = {'resolving': 'running', 'downloading': 'running',
              'converting': 'running', 'completed': 'done',
              'failed': 'failed', 'cancelled': 'cancelled',
              'interrupted': 'failed'}
_DL_STAGE = {'resolving': '解析编号', 'downloading': '下载中',
             'converting': '解码 .sra', 'completed': '已完成',
             'failed': '失败', 'cancelled': '已取消',
             'interrupted': '已中断（可续传）'}


def download_task_provider(bid):
    """下载批次 → 任务中心快照字段（provider 回调，返回 None 表示已删除）。

    这是「下载接入任务中心」的桥：执行仍归 DownloadManager（自建线程、
    自管状态、自管取消），这里只把它的状态翻译成任务中心的统一口径。
    """
    m = get_dl_manager()
    snap = m.snapshot(bid, include_files=False)
    if not snap:
        return None
    st = _DL_STATUS.get(snap.get('status'), 'running')
    parts = [f"已下载 {_fmt_bytes(snap.get('done_bytes'))}"
             f" / {_fmt_bytes(snap.get('total_bytes'))}"]
    spd = _fmt_speed(snap.get('speed'))
    if spd != '—':
        parts.append(f'速度 {spd}')
    eta = int(snap.get('eta') or 0)
    if eta > 0:
        parts.append(f'剩余 {eta}s')
    parts.append(f"{snap.get('done', 0)}/{snap.get('n_files', 0)} 文件")
    if snap.get('failed'):
        parts.append(f"{snap['failed']} 失败")
    if snap.get('convert_failed'):
        parts.append(f"{snap['convert_failed']} 个 .sra 未转 FASTQ")
    err = snap.get('error') or ''
    return {'status': st,
            'pct': snap.get('overall'),
            'stage': _DL_STAGE.get(snap.get('status'), snap.get('status') or ''),
            'msg': ' · '.join(parts),
            'error': err,
            'eta': eta or None,
            'log': m._read_log_tail(m.batches.get(bid) or {}, 40),
            'started': snap.get('started'),
            'finished': snap.get('finished')}


def bind_downloads():
    """把下载批次挂进任务中心（app 启动时调用一次）。

    已有批次在 add_listener 里补挂，之后新建的批次由 _notify 自动挂上。
    外部任务不占并发名额、不落 tasks/*.json；「停止」经 on_cancel 转发到
    DownloadManager.cancel（会 terminate aria2c/sracha 子进程）。
    """
    from Virus_Platform_Core.web.tasks import tm

    def _attach(bid):
        tm.attach(f'下载·{bid}', lambda b=bid: download_task_provider(b),
                  on_cancel=lambda b=bid: get_dl_manager().cancel(b),
                  link='/download?pop=1', key=f'dl:{bid}', weight='light')

    get_dl_manager().add_listener(_attach)


bp = Blueprint('download', __name__)


@bp.route('/download')
def page_download():
    return render_template('download.html')


@bp.route('/api/dl/create', methods=['POST'])
def api_dl_create():
    body = request.get_json(force=True) or {}
    name = str(body.get('name') or '').strip() or 'batch'
    items = body.get('items') or []
    items = [x.strip() for x in items if x and x.strip()]
    if not items:
        abort(400, '请提供要下载的编号/URL 列表')
    if len(items) > 500:
        abort(400, '单批最多 500 条（大批次请拆分）')
    try:
        bid = get_dl_manager().create(
            name, items,
            concurrency=int(body.get('concurrency') or 2),
            convert_sra=bool(body.get('convert_sra', True)),
            out_dir=str(body.get('out_dir') or '').strip() or None)
    except ValueError as e:
        # 自定义输出位置越界（仅允许平台内/已挂载输出根内）等校验失败：
        # 把原因原样讲给用户（400），而不是 500 哑错
        abort(400, str(e))
    return jsonify({'batch': bid})


@bp.route('/api/dl/batches')
def api_dl_batches():
    return jsonify(get_dl_manager().list_snapshots())


@bp.route('/api/dl/batch/<bid>')
def api_dl_batch(bid):
    snap = get_dl_manager().snapshot(bid)
    if not snap:
        abort(404, '批次不存在')
    return jsonify(snap)


@bp.route('/api/dl/batch/<bid>/<action>', methods=['POST'])
def api_dl_batch_action(bid, action):
    """批次操作：cancel / retry / retry_convert（只重跑 .sra 转换）/
    delete（记录+文件全删）/ delete_files（只删文件，留记录与日志）/
    delete_record（只删记录，留文件）。"""
    m = get_dl_manager()
    if action == 'cancel':
        ok = m.cancel(bid)
    elif action == 'retry':
        n = m.retry_failed(bid)
        return jsonify({'ok': bool(n), 'retried': n})
    elif action == 'retry_convert':
        n = m.retry_convert(bid)
        return jsonify({'ok': bool(n), 'retried': n})
    elif action == 'delete':
        ok = m.delete(bid)
    elif action == 'delete_files':
        ok = m.delete_files(bid)
    elif action == 'delete_record':
        ok = m.delete_record(bid)
    else:
        abort(400, '无效操作')
    if not ok:
        abort(400, '操作失败（批次可能不存在）')
    return jsonify({'ok': True})


def _dl_manager():
    from Virus_Platform_Core import public_data
    return public_data.get_manager()


get_dl_manager = _dl_manager


@bp.route('/downloads/<bid>/<path:filename>')
def page_dl_file(bid, filename):
    """下载批次内文件（严格限定在该批次目录内）。"""
    if not re.fullmatch(r'[A-Za-z0-9_\-.]+', bid) or '..' in filename:
        abort(400, '无效路径')
    from Virus_Platform_Core.public_data import _batch_root
    p = check_path(os.path.join(_batch_root(), bid, filename),
                   must_exist=True, in_platform=True)
    bdir = check_path(os.path.join(_batch_root(), bid),
                      must_exist=True, in_platform=True)
    if not p.startswith(bdir + os.sep):
        abort(400, '路径越界')
    return send_file(p, as_attachment=request.args.get('dl') == '1')


@bp.route('/api/dl/to_pipeline', methods=['POST'])
def api_dl_to_pipeline():
    """把下载完成的 run 转成本地样品（可选加入批处理队列）。"""
    body = request.get_json(force=True) or {}
    bid = body.get('batch') or ''
    runs = body.get('runs') or []
    project = str(body.get('project') or '') or None
    enqueue = bool(body.get('enqueue'))
    if not bid or not re.fullmatch(r'[A-Za-z0-9_\-.]+', bid):
        abort(400, '无效批次 ID')
    m = get_dl_manager()
    ready = m.ready_files(bid)
    if runs:
        ready = {k: v for k, v in ready.items() if k in runs}
    if not ready:
        abort(400, '该批次没有可用的 FASTQ 文件（可能仍在下载或需要 .sra 转换）')
    from Virus_Platform_Core.pipeline import _safe_sample_name, save_sample_input
    created, skipped = [], []
    for acc, files in ready.items():
        sample = _safe_sample_name(acc)
        r1 = files.get('r1') or files.get('single')
        r2 = files.get('r2')
        # 先判可用性、再建目录：原实现先 makedirs 再判 r1，于是只有 _2
        # （没有 _1/single）的条目会在 results/ 下留下一个**空样品目录**，
        # 样品列表里就多出一个打不开的坏样品。
        if not r1:
            skipped.append(sample)
            continue
        sd = check_path(os.path.join(DIRS['results'], sample),
                        must_exist=False, in_platform=True)
        if os.path.isdir(sd) and any(os.scandir(sd)):
            skipped.append(sample)
            continue
        os.makedirs(sd, exist_ok=True)
        save_sample_input(sd, r1, r2, sample, project=project)
        created.append(sample)
    queued = 0
    if enqueue and created:
        added = _queue().add(created, None, {}, project=project)
        # 与 samples.py 的 /api/queue/add 对齐：给「待启动」条目挂任务中心
        # 卡片。不挂的话条目只在分析流程页可见，任务中心里整批"消失"
        # （用户以为入队没生效）；卡片在真任务开跑后自动让位。
        from Virus_Platform_Core.web.samples import bind_queue_entry
        for it in (added or ()):
            bind_queue_entry(it)
        queued = len(created)
    return jsonify({'created': created, 'skipped': skipped, 'queued': queued})
