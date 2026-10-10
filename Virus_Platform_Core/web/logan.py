# -*- coding: utf-8 -*-
"""LOGAN 溯源（自 app.py 拆出）。

病毒序列 → Logan-Search → 物种/样本溯源；含 Selenium 批量提交。"""
import os

from flask import (Blueprint, abort, jsonify, render_template, request,
                   send_file)

from Virus_Platform_Core.utils import (TaskLogger, check_path)
from Virus_Platform_Core.web.state import cfg
from Virus_Platform_Core.web.tasks import tm
from Virus_Platform_Core.web.common import _safe_sample

bp = Blueprint('logan', __name__)


@bp.route('/logan')
def page_logan():
    return render_template('logan.html')


@bp.route('/logan/report/<name>/')
def page_logan_report(name):
    from Virus_Platform_Core.logan_trace import _job_dir
    p = os.path.join(_job_dir(name), 'trace_report.html')
    if not os.path.isfile(p):
        abort(404, '溯源报告尚未生成（请先导入结果文件）')
    resp = send_file(check_path(p, must_exist=True, in_platform=True))
    # 报告内容由样本数据（contig 名等）拼成，属于半可信输入：加 CSP 兜底，
    # 即使未来某处转义遗漏，注入的脚本也无法读取平台 API 响应
    resp.headers['Content-Security-Policy'] = (
        "default-src 'self'; script-src 'self' 'unsafe-inline' 'unsafe-eval'; "
        "style-src 'self' 'unsafe-inline'; img-src 'self' data:; "
        "connect-src 'self'; frame-ancestors 'self'")
    return resp


@bp.route('/logan/report/<name>/<path:filename>')
def page_logan_report_file(name, filename):
    """溯源报告目录内静态文件（plotly.min.js / 原始结果表 / 查询 FASTA）。"""
    from Virus_Platform_Core.logan_trace import _job_dir
    # base 归一 + 前缀包含校验：filename 传绝对路径时 os.path.join 会整段
    # 替换 base，check_path(in_platform=True) 拦不住读平台内其它文件。
    base = os.path.abspath(_job_dir(name))
    p = check_path(os.path.join(base, filename),
                   must_exist=True, in_platform=True)
    if not p.startswith(base + os.sep):
        abort(400, '非法路径')
    return send_file(p)


@bp.route('/api/logan/samples')
def api_logan_samples():
    """可选来源样品：已有 ③组装 病毒 contigs 产出的样品。"""
    from Virus_Platform_Core.logan_trace import list_query_samples
    return jsonify(list_query_samples())


@bp.route('/api/logan/contigs/<sample>')
def api_logan_contigs(sample):
    from Virus_Platform_Core.logan_trace import list_virus_contigs
    return jsonify(list_virus_contigs(_safe_sample(sample)))


@bp.route('/api/logan/jobs')
def api_logan_jobs():
    from Virus_Platform_Core.logan_trace import list_jobs
    return jsonify(list_jobs())


@bp.route('/api/logan/job/<name>')
def api_logan_job(name):
    from Virus_Platform_Core.logan_trace import job_detail
    return jsonify(job_detail(name))


@bp.route('/api/logan/job/<name>/segment/<int:idx>')
def api_logan_segment(name, idx):
    """片段序列（前端复制提交用）。"""
    from Virus_Platform_Core.logan_trace import get_segment_fasta
    header, seq = get_segment_fasta(name, idx)
    return jsonify({'header': header, 'seq': seq})


@bp.route('/api/logan/job/<name>/query')
def api_logan_query_dl(name):
    """下载查询 FASTA（全部片段）。"""
    from Virus_Platform_Core.logan_trace import _job_dir
    p = os.path.join(_job_dir(name), 'query_all.fasta')
    if not os.path.isfile(p):
        abort(404, '查询 FASTA 不存在')
    return send_file(check_path(p, must_exist=True, in_platform=True),
                     as_attachment=True)


@bp.route('/api/logan/create', methods=['POST'])
def api_logan_create():
    from Virus_Platform_Core.logan_trace import create_job
    body = request.get_json(force=True) or {}
    name = (body.get('name') or '').strip()
    if not name:
        abort(400, '请填写查询名称')
    pasted = (body.get('fasta') or '').strip()
    sample = (body.get('sample') or '').strip()
    contigs = body.get('contigs') or []
    if not pasted and not (sample and contigs):
        abort(400, '请勾选样品的病毒 contigs，或粘贴序列')
    try:
        # ⚠️ sample 为空时必须给 None：_safe_sample('') 会回退成默认名
        # 'sample'（truthy），`or None` 兜不住，粘贴模式也会被拽去读
        # results/sample/03_assembly/virus_contigs.tsv → 必 400。
        return jsonify(create_job(name, sample=_safe_sample(sample) if sample else None,
                                  contig_ids=contigs or None,
                                  pasted=pasted or None,
                                  n_seg=int(body.get('n_seg') or 2)))
    except ValueError as e:
        abort(400, str(e))


@bp.route('/api/logan/job/<name>/import/<int:idx>', methods=['POST'])
def api_logan_import(name, idx):
    """导入某片段的 Logan-Search 结果表（CSV/TSV），解析并自动生成报告。"""
    from Virus_Platform_Core.logan_trace import import_result, MAX_UPLOAD_BYTES
    f = request.files.get('file')
    if not f:
        abort(400, '缺少结果文件')
    raw = f.read(MAX_UPLOAD_BYTES + 1)
    if len(raw) > MAX_UPLOAD_BYTES:
        abort(400, '结果文件超过 64MB 上限')
    if not raw.strip():
        abort(400, '结果文件为空')
    try:
        return jsonify(import_result(name, idx, f.filename or '', raw))
    except ValueError as e:
        abort(400, str(e))


@bp.route('/api/logan/job/<name>/delete', methods=['POST'])
def api_logan_delete(name):
    from Virus_Platform_Core.logan_trace import delete_job
    delete_job(name)
    return jsonify({'ok': True})


@bp.route('/api/logan/batch_ready')
def api_logan_batch_ready():
    """批量模式可用性（selenium 是否已安装）。"""
    from Virus_Platform_Core.logan_trace import selenium_ready, GROUP_HINTS
    return jsonify({'selenium': selenium_ready(), 'groups': GROUP_HINTS})


@bp.route('/api/logan/job/<name>/batch', methods=['POST'])
def api_logan_batch(name):
    """后台批量提交任务的全部未导入片段，完成后自动导入并生成报告。"""
    from Virus_Platform_Core.logan_trace import batch_submit, _load_job
    _jdir, meta, segs = _load_job(name)
    if not any(not s.get('imported') for s in segs):
        abort(400, '该任务所有片段均已导入结果')
    for t in tm.list_all():
        if t.get('status') == 'running' and t.get('name') in {
            f'LOGAN 批量 {meta["name"]}',
            f'LOGAN batch {meta["name"]}',
        }:
            abort(400, f'任务 {meta["name"]} 已有批量提交在运行，请等它结束')
    body = request.get_json(force=True) or {}
    emails = [e.strip() for e in (body.get('emails') or '').split(',') if e.strip()]
    if not emails:
        abort(400, '请至少填写一个通知邮箱（用于接收 Logan-Search 完成通知）')
    group = body.get('group') or 'Fast_No_human'

    def job(log, prog, cancel):
        logger = TaskLogger(callback=log)
        batch_submit(meta['name'], emails, group=group,
                     headless=bool(body.get('headless', True)),
                     first_wait=int(body.get('first_wait') or 300),
                     max_wait=int(body.get('max_wait') or 1800),
                     logger=logger, progress=prog, cancel=cancel)
        logger.close()
        return 'ok'

    tid = tm.start(cfg.tr(f'LOGAN 批量 {meta["name"]}',
                         f'LOGAN batch {meta["name"]}'), job,
                   weight='light')
    return jsonify({'task': tid})
