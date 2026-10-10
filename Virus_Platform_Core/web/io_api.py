# -*- coding: utf-8 -*-
"""输入输出杂项 API（自 app.py 拆出）。

粘贴序列、拖拽上传、序列查看器、打开平台目录。"""
import os
import re
import threading
import time
import uuid

from flask import (Blueprint, abort, jsonify, request)

from Virus_Platform_Core.config import DIRS, PLATFORM_ROOT
from Virus_Platform_Core.utils import (check_path, safe_open)

bp = Blueprint('io_api', __name__)

# 序列查看器最多扫描的字节数（防止误点一个多 GB 的组装产物把工作线程占满）
_SEQVIEW_MAX_BYTES = 64 << 20
# 粘贴序列上限：粘贴框是给人贴 contig/短序列的，几十 MB 的文本只会
# 撑爆 JSON 解析和内存，应走文件输入
_PASTE_MAX_BYTES = 20 << 20
# 单次上传上限（Kraken2 库包可到数十 GB，故给得很宽；仅拦手滑/恶意超大包）
_UPLOAD_MAX_BYTES = 64 << 30


@bp.route('/api/paste_input', methods=['POST'])
def api_paste_input():
    """粘贴序列文本 → 写入 uploads/paste_<ts>_<rand>.<ext>，返回平台相对路径。

    body: {text: "FASTA/FASTQ 文本", ext: ".fasta" | ".fastq" | ".tsv"}
    所有工具的文件输入框均可使用返回的 path。
    """
    body = request.get_json(force=True) or {}
    text = (body.get('text') or '').strip()
    ext = (body.get('ext') or '.fasta').strip().lower()
    if not text:
        abort(400, '粘贴内容为空')
    if len(text.encode('utf-8', errors='replace')) > _PASTE_MAX_BYTES:
        abort(413, '粘贴内容过大（>20MB），请保存为文件后用文件方式输入')
    if not re.fullmatch(r'\.(fasta|fa|fna|fas|fastq|fq|tsv|txt)', ext):
        abort(400, '不支持的文件类型')
    # 基本格式校验
    first_line = text.split('\n', 1)[0].strip()
    if ext in ('.fasta', '.fa', '.fna', '.fas') and not first_line.startswith('>'):
        abort(400, 'FASTA 格式应以 > 开头')
    if ext in ('.fastq', '.fq') and not first_line.startswith('@'):
        abort(400, 'FASTQ 格式应以 @ 开头')
    # 文件名带毫秒 + 随机段：秒级时间戳在同秒两次提交时会互相覆盖
    ts = time.strftime('%Y%m%d_%H%M%S')
    fname = f'paste_{ts}_{uuid.uuid4().hex[:6]}{ext}'
    fdir = os.path.join(PLATFORM_ROOT, 'run', 'uploads')
    os.makedirs(fdir, exist_ok=True)
    fp = os.path.join(fdir, fname)
    with safe_open(fp, 'wt') as f:
        f.write(text + '\n')
    return jsonify({'path': f'run/uploads/{fname}', 'size': len(text)})


@bp.route('/api/upload', methods=['POST'])
def api_upload():
    """拖拽/选择上传数据文件到平台 uploads/ 目录（流式落盘，支持大文件）。

    返回 {'path': 'uploads/<文件名>'} 供前端填入输入框。"""
    # Content-Length 预检（不读 body 即拒绝）：Kraken2 库包可到数十 GB，
    # 上限设得很宽，只拦手滑/恶意的超大包
    cl = request.content_length
    if cl and cl > _UPLOAD_MAX_BYTES:
        abort(413, f'上传文件过大（上限 {_UPLOAD_MAX_BYTES >> 30}GB）')
    # 无 Content-Length 且 Transfer-Encoding: chunked 的请求没法预检大小：
    # 直接读会一路落盘到写满磁盘。这里选择「直接拒绝」（411 Length Required）
    # 而不是"边读边数、超限即断"——本平台的上传方（页面 fetch/XHR +
    # FormData、拖拽代理）都会带 Content-Length，chunked 只可能是异常/恶意
    # 客户端；拒绝比读完 64GB 再断更省资源，错误信息也更清楚。
    if not cl and 'chunked' in (request.headers.get('Transfer-Encoding')
                                or '').lower():
        abort(411, '上传请求需带 Content-Length（不支持 chunked 无长度上传）')
    f = request.files.get('file')
    if not f:
        abort(400, '缺少上传文件')
    name = os.path.basename(f.filename or '')
    name = re.sub(r'[^\w.\-\u4e00-\u9fff]+', '_', name).strip('._')
    if not name:
        abort(400, '文件名无效')
    updir = check_path(DIRS.get('uploads')
                       or os.path.join(PLATFORM_ROOT, 'run', 'uploads'),
                       must_exist=False, in_platform=True)
    os.makedirs(updir, exist_ok=True)
    dst = check_path(os.path.join(updir, name), must_exist=False,
                     in_platform=True)
    if os.path.isfile(dst):                      # 重名不覆盖：追加序号
        base, ext = os.path.splitext(name)
        i = 1
        while os.path.isfile(check_path(
                os.path.join(updir, f'{base}_{i}{ext}'),
                must_exist=False, in_platform=True)):
            i += 1
        dst = check_path(os.path.join(updir, f'{base}_{i}{ext}'),
                         must_exist=False, in_platform=True)
    f.save(str(dst))
    up_disp = (os.path.abspath(dst) if DIRS.get('uploads')
               and not DIRS['uploads'].startswith(PLATFORM_ROOT)
               else 'run/uploads/' + os.path.basename(dst))
    return jsonify({'path': up_disp,
                    'size': os.path.getsize(dst)})


_DATA_EXTS = ('.fasta', '.fa', '.fna', '.fas', '.ffn')


@bp.route('/api/seqview')
def api_seqview():
    """轻量序列查看器：流式统计 + 分页预览（FASTA / FASTA.gz）。

    大文件按 _SEQVIEW_MAX_BYTES 截断扫描（返回 byte_capped=true）：
    统计是逐记录累加的，读到一个多 GB 的组装产物会白占用工作线程，
    而查看器本身只展示前若干条记录。
    """
    rel = request.args.get('path') or ''
    page = max(0, int(request.args.get('page', 0) or 0))
    per = 50
    # 平台内相对路径或任意绝对路径均可（只读流式统计，无写风险）；
    # 与文件浏览对话框（可选任意盘文件）行为对齐。
    if os.path.isabs(rel):
        p = check_path(rel, must_exist=True)
    else:
        p = check_path(os.path.join(PLATFORM_ROOT, rel), must_exist=True,
                       in_platform=True)
    if not p.lower().endswith(_DATA_EXTS +
                              tuple(e + '.gz' for e in _DATA_EXTS)):
        abort(400, '仅支持 FASTA / FASTA.gz')
    rows = []
    total_bp = 0
    truncated = False
    byte_capped = False
    rec_i = -1
    n_read = 0

    def _push(h, seq):
        nonlocal rec_i, total_bp, truncated
        rec_i += 1
        seq = seq.upper()
        if not seq:
            return
        n_gc = seq.count('G') + seq.count('C')
        n_deg = sum(seq.count(c) for c in 'RYKMSWBDHVN')
        total_bp += len(seq)
        if rec_i >= 20000:
            truncated = True
            return
        row = {'id': h.split()[0] if h.split() else h[:30],
               'len': len(seq),
               'gc': round(n_gc * 100.0 / len(seq), 1),
               'deg': round(n_deg * 100.0 / len(seq), 1)}
        if page * per <= rec_i < page * per + per:
            row['preview'] = seq[:300]
        rows.append(row)

    try:
        with safe_open(p) as f:
            h, buf = None, []
            for line in f:
                n_read += len(line)
                if n_read > _SEQVIEW_MAX_BYTES:
                    byte_capped = truncated = True
                    break
                line = line.strip()
                if line.startswith('>'):
                    if h is not None:
                        _push(h, ''.join(buf))
                    h, buf = line[1:], []
                elif h is not None and line:
                    buf.append(line)
            if h is not None and not byte_capped:
                _push(h, ''.join(buf))
    except (OSError, ValueError) as e:
        abort(400, f'读取失败: {e}')
    page_rows = [r for r in rows if 'preview' in r]
    pages = max(1, (min(rec_i + 1, 20000) + per - 1) // per)
    return jsonify({'total': rec_i + 1, 'total_bp': total_bp,
                    'page': page, 'pages': pages, 'per': per,
                    'truncated': truncated, 'byte_capped': byte_capped,
                    'rows': page_rows})


@bp.route('/api/open_platform_dir')
def api_open_platform_dir():
    from Virus_Platform_Core.utils import open_in_explorer
    open_in_explorer(check_path(PLATFORM_ROOT, must_exist=True, in_platform=True))
    return jsonify({'ok': True})


@bp.route('/api/open_dir', methods=['POST'])
def api_open_dir():
    """在资源管理器中打开一个平台内的目录（任务卡/结果面板的「打开」按钮）。

    body: {path: "绝对或平台相对路径"}。check_path 强制目录位于平台或
    自定义输出根之内，杜绝任意目录打开。
    """
    body = request.get_json(force=True) or {}
    p = (body.get('path') or '').strip()
    if not p:
        abort(400, '缺少路径')
    try:
        d = check_path(p, must_exist=True, in_platform=True)
    except (ValueError, FileNotFoundError) as e:
        abort(400, f'无法打开（{e}）')
    if not os.path.isdir(d):
        # 给的是文件路径时打开其所在目录，体验更顺（结果面板常持有文件路径）
        d = os.path.dirname(d)
        if not os.path.isdir(d):
            abort(400, '目录不存在')
    from Virus_Platform_Core.utils import open_in_explorer
    open_in_explorer(d)
    return jsonify({'ok': True})


# ── 拖拽代理（scripts/drag_proxy.py 悬浮窗）的路径接收 ──
# 代理窗口收到资源管理器拖入的真实路径后 POST 到这里；页面轮询取走
# 填入「最后点过的输入框」。内存态即可（本机 127.0.0.1，无持久化需求）。
_drop_state = {'ts': 0.0, 'paths': []}
_drop_lock = threading.Lock()


@bp.route('/api/dropped_paths', methods=['POST'])
def api_dropped_paths_post():
    body = request.get_json(silent=True) or {}
    paths = [str(x).strip() for x in (body.get('paths') or [])
             if str(x).strip()]
    if not paths:
        abort(400, '缺少 paths')
    with _drop_lock:
        _drop_state['ts'] = time.time()
        _drop_state['paths'] = paths
    return jsonify({'ok': True, 'n': len(paths)})


@bp.route('/api/dropped_paths/last')
def api_dropped_paths_last():
    try:
        since = float(request.args.get('since') or 0)
    except ValueError:
        since = 0.0
    with _drop_lock:
        ts = _drop_state['ts']
        paths = list(_drop_state['paths'])
    if ts <= since:
        return jsonify({'fresh': False, 'ts': ts})
    return jsonify({'fresh': True, 'ts': ts, 'paths': paths})
