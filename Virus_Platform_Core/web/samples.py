# -*- coding: utf-8 -*-
"""样品与批处理队列（自 app.py 拆出）。

SampleQueue 顺序队列、样品 CRUD、管道逐级运行、分析任务提交。"""
import json
import logging
import os
import re
import shutil
import threading
import time
import uuid

from flask import (Blueprint, abort, jsonify, request,
                   send_file)

from Virus_Platform_Core.config import DIRS, PLATFORM_ROOT
from Virus_Platform_Core.utils import (TaskLogger, check_path, fmt_size,
                                      resolve_sample_name, safe_open,
                                      total_ram)
from Virus_Platform_Core.web.common import _safe_sample
from Virus_Platform_Core.web.state import cfg
from Virus_Platform_Core.web.tasks import tm

_log = logging.getLogger('Virus_Platform_Core.web.samples')

bp = Blueprint('samples', __name__)


class SampleQueue:
    """样品批处理队列：FIFO 顺序跑，每项复用与单样品运行完全相同的
    run_analysis 代码路径；状态持久化到 tasks/queue.json。"""

    QUEUE_FILE = os.path.join(DIRS['tasks'], 'queue.json')

    def __init__(self):
        self.lock = threading.RLock()
        self.items = []
        self.thread = None
        self._worker_alive = False
        self._paused = False
        self._load()

    def _load(self):
        try:
            with safe_open(self.QUEUE_FILE) as f:
                data = json.load(f)
            # 兼容旧格式（裸列表）：暂停标记是后来才加的字段
            if isinstance(data, list):
                data = {'items': data, 'paused': False}
            self.items = data.get('items') or []
            self._paused = bool(data.get('paused'))
            for it in self.items:
                if it.get('status') in ('queued', 'running'):
                    it['status'] = 'queued'
                    it['task'] = None
        except (OSError, ValueError):
            self.items = []

    def _persist(self):
        with self.lock:
            try:
                os.makedirs(os.path.dirname(self.QUEUE_FILE), exist_ok=True)
                # 原子写：断电/崩溃不损坏现有队列
                tmp = self.QUEUE_FILE + '.tmp'
                with safe_open(tmp, 'wt') as f:
                    json.dump({'items': self.items, 'paused': self._paused},
                              f, ensure_ascii=False, indent=1)
                os.replace(tmp, self.QUEUE_FILE)
            except OSError:
                pass

    def add(self, samples, stages, params, project=None):
        """samples: [样品名]（须已创建并含 input.json）。

        入队只登记为 queued（**待手动启动**），不拉起 worker —— 2026-09-20
        用户要求「入队后必须手动点开始」。点「启动」（单条 start）或
        「全部启动」（start_all）把条目转为 ready 后才由 worker 依序执行。
        """
        added = []
        with self.lock:
            existing = {it['sample'] for it in self.items
                        if it['status'] in ('queued', 'ready', 'running')}
            for s in samples:
                if s in existing:
                    continue
                it = {'id': uuid.uuid4().hex[:10], 'sample': s,
                      'stages': stages, 'params': params,
                      'project': project or '', 'status': 'queued',
                      'task': None, 'error': '', 'added': time.time()}
                self.items.append(it)
                added.append(it)
        self._persist()
        return added

    def start(self, entry_id):
        """手动启动一条：queued/ready → ready，并拉起 worker。

        对已 ready 条目幂等（服务重启后 worker 不在，再点一次启动即恢复）。
        返回 False 表示条目不存在或状态不可启动。
        """
        with self.lock:
            it = next((x for x in self.items if x['id'] == entry_id), None)
            if not it or it['status'] not in ('queued', 'ready'):
                return False
            it['status'] = 'ready'
        self._persist()
        self._ensure_worker()
        return True

    def start_all(self):
        """全部「待启动」→「排队中」，返回转掉的条数（没有可转的就不拉 worker）。"""
        with self.lock:
            n = 0
            for it in self.items:
                if it['status'] == 'queued':
                    it['status'] = 'ready'
                    n += 1
        self._persist()
        if n:
            self._ensure_worker()
        return n

    def _ensure_worker(self):
        with self.lock:
            if self._worker_alive:
                return
            self._worker_alive = True
        threading.Thread(target=self._worker, daemon=True,
                         name='sample-queue').start()

    def _worker(self):
        while True:
            with self.lock:
                # 暂停中：worker 退出（resume 会重新拉起），已就绪条目原地不动。
                # 正在跑的那条不受影响（worker 正阻塞在 _run_entry 里）。
                if self._paused:
                    self._worker_alive = False
                    self.thread = None
                    return
                # 只消费 ready（用户已点启动）的条目；queued = 待手动启动
                nxt = next((it for it in self.items
                            if it['status'] == 'ready'), None)
                if not nxt:
                    # 取件与退出标志清除必须在同一把锁内：否则「worker 判空
                    # 退出」与「start() 启动后发现线程还活着不再起新 worker」
                    # 交错时，新条目会永远没人消费（队列假死）
                    self._worker_alive = False
                    self.thread = None
                    return
                nxt['status'] = 'running'
            self._persist()
            try:
                self._run_entry(nxt)
            except Exception as e:
                nxt['status'] = 'failed'
                nxt['error'] = str(e)
            self._persist()

    def pause(self):
        """暂停队列：正在跑的条目继续跑完，已就绪（ready）条目暂不开跑，
        直到恢复或再次点启动。"""
        with self.lock:
            self._paused = True
        self._persist()

    def resume(self):
        with self.lock:
            self._paused = False
        self._persist()
        self._ensure_worker()

    def _run_entry(self, it):
        from Virus_Platform_Core.pipeline import (load_sample_input, run_analysis,
                                 pipeline_overview, STAGE_ORDER)
        # 队列项里的样品名是"已存在的样品"，按磁盘真实目录名解析
        s = resolve_sample_name(it['sample'], DIRS['results'])
        sd = check_path(os.path.join(DIRS['results'], s),
                        must_exist=True, in_platform=True)
        r1, r2, _proj = load_sample_input(sd)
        if not r1:
            raise RuntimeError(f'样品 {s} 缺少输入记录（input.json）')
        check_path(r1, must_exist=True)
        stages = it.get('stages') or None
        if stages:
            stages = [x for x in stages if x in STAGE_ORDER]
        if not stages:
            # 与「依次运行剩余步骤」一致：所有未完成且可用的阶段
            try:
                ov = pipeline_overview(sd)['stages']
                done = {x['stage'] for x in ov if x['status'] == 'done'}
                stages = [x['stage'] for x in ov
                          if x['stage'] in STAGE_ORDER
                          and x['stage'] not in done
                          and x['status'] != 'unavailable']
            except Exception:
                stages = []
        if not stages:
            raise RuntimeError(f'样品 {s} 没有需要运行的阶段（全部已完成）')

        if not _sample_busy_register(s):
            raise RuntimeError(
                f'样品 {s} 已有任务在运行，队列项跳过（请勿重复入队）')

        def job(log, prog, cancel):
            logger = TaskLogger(callback=log)
            try:
                run_analysis(sample=s, r1=r1, r2=r2, stages=stages,
                             **_analysis_kwargs(it.get('params') or {}),
                             logger=logger, progress=prog)
                return sd
            finally:
                _sample_busy_release(s)
                logger.close()

        q_log_dir = check_path(os.path.join(sd, 'logs'),
                               must_exist=False, in_platform=True)
        os.makedirs(q_log_dir, exist_ok=True)
        try:
            # on_finish：终态必调一次（含「排队途中被取消」——那种情况下 job
            # 根本不会执行，其 finally 也就不存在）。2026-10-07 审查 P1-1。
            tid = tm.start(cfg.tr(f'队列·{s}', f'Queue·{s}'), job,
                           log_file=os.path.join(
                               q_log_dir,
                               f'queue_{time.strftime("%Y%m%d_%H%M%S")}.log'),
                           on_finish=lambda: _sample_busy_release(s))
        except BaseException:
            _sample_busy_release(s)      # 任务没起来也要放掉占位
            raise
        it['task'] = tid
        self._persist()
        while True:
            time.sleep(2.0)
            snap = tm.snapshot(tid, log_lines=0, include_result=False)
            if snap is None or snap['status'] != 'running':
                # 三态显式映射：原先把 cancelled 并进 failed 且 error 常为空，
                # 队列卡片只显示"失败"、没有任何原因（用户以为出错了，其实
                # 是自己停的任务）。done/failed/cancelled 分别落位，
                # cancelled 写明"用户取消"。
                st = (snap or {}).get('status')
                if st == 'done':
                    it['status'] = 'done'
                    it['error'] = ''
                elif st == 'cancelled':
                    it['status'] = 'cancelled'
                    it['error'] = ((snap or {}).get('error')
                                   or f'用户取消（{s} 的任务已停止）')
                else:
                    it['status'] = 'failed'
                    it['error'] = ((snap or {}).get('error')
                                   or ('任务记录已消失（可能被中断或清理）'
                                       if snap is None else
                                       f'{s} 的任务执行失败（详见任务日志）'))
                return

    def snapshot(self):
        with self.lock:
            items = [dict(it) for it in self.items]
            paused = self._paused
        running = any(it['status'] == 'running' for it in items)
        return {'running': running, 'paused': paused, 'items': items}

    def remove(self, entry_id):
        with self.lock:
            it = next((x for x in self.items if x['id'] == entry_id), None)
            if not it or it['status'] == 'running':
                return False
            self.items.remove(it)
        self._persist()
        return True

    def clear_finished(self):
        with self.lock:
            self.items = [it for it in self.items
                          if it['status'] in ('queued', 'running')]
        self._persist()
        return True

    def reset_all(self):
        """清空批处理队列全部条目（全局重置用），返回被清掉的条数。

        与 clear_finished 的差别：连 queued/running 条目一起清。
        不在这里取消任务——任务取消由 TaskManager.cancel_all() 负责；
        清空后 worker 取不到 queued 项会自行退出，正在跑的那条其
        `_run_entry` 等待循环在任务进终态后正常返回（条目已不在列表里，
        收尾的 _persist 只是把空列表再写一遍）。
        不碰 results/ 下任何文件。
        """
        with self.lock:
            n = len(self.items)
            self.items = []
        self._persist()
        return n


sample_queue = SampleQueue()
# 注册到 state，供 download blueprint 的「转入分析流程」取用
from Virus_Platform_Core.web import state as _web_state

_web_state.set_sample_queue(sample_queue)


def _task_matches_sample(task_name, sample):
    name = str(task_name or '')
    sample = str(sample or '')
    if not name or not sample:
        return False
    if name.startswith(sample + ' · '):
        return True
    # 批处理队列的任务名为「队列·<样品>」/「Queue·<样品>」（分隔符无空格），
    # 与上面的「<样品> · 」不是同一套写法，必须单独匹配，
    # 否则队列运行中删除/清空样品时 _sample_busy 会漏判。
    if re.match(rf'^(?:队列|Queue)·{re.escape(sample)}$', name):
        return True
    return bool(re.match(rf'^(?:样品分析|Sample analysis)\s+{re.escape(sample)}(?:\b|$)',
                         name))


# ── 同样品并发防护（登记表） ─────────────────────────────────
# 原先只靠「遍历 tm.list_all() 按任务名正则匹配」判定样品是否在跑：
# ① 检查与 tm.start 之间有 TOCTOU 窗口（双击可同时通过）；
# ② 工具页直接跑某样品阶段目录的任务名匹配不上那几条正则。
# 两个任务并发写同一 results/<sample>/ 会互相覆盖 summary.json 与
# .done 断点标记。这里用进程内登记表做权威判定：请求/入队时占位，
# 任务 finally 释放（取消、异常都会走到）。
_RUNNING_SAMPLES = set()
_RUNNING_SAMPLES_LOCK = threading.Lock()


def _sample_busy_register(sample):
    """占位。返回 False 表示该样品已有任务占位（应拒绝新任务）。"""
    with _RUNNING_SAMPLES_LOCK:
        if sample in _RUNNING_SAMPLES:
            return False
        _RUNNING_SAMPLES.add(sample)
        return True


def _sample_busy_release(sample):
    with _RUNNING_SAMPLES_LOCK:
        _RUNNING_SAMPLES.discard(sample)


@bp.route('/api/queue')
def api_queue():
    return jsonify(sample_queue.snapshot())


@bp.route('/api/queue/add', methods=['POST'])
def api_queue_add():
    body = request.get_json(force=True) or {}
    samples = body.get('samples') or []
    if not samples:
        abort(400, '请提供样品列表')
    valid = []
    for s in samples:
        # 入队的都是"已有样品"：按磁盘真实目录名解析后再校验 input.json，
        # 否则手工/历史命名的样品会被字符转换改名而误判为"不存在"。
        name = resolve_sample_name(s, DIRS['results'])
        sd = check_path(os.path.join(DIRS['results'], name),
                        must_exist=False, in_platform=True)
        if not os.path.isfile(os.path.join(sd, '00_prep', 'input.json')):
            abort(400, f'样品 {s} 不存在或缺少输入记录（请先创建样品）')
        valid.append(name)
    stages = body.get('stages') or None       # None = 依次运行剩余步骤
    if stages:
        from Virus_Platform_Core.pipeline import STAGE_ORDER
        stages = [s for s in stages if s in STAGE_ORDER]
        if not stages:
            abort(400, 'stages 参数不合法')
    added = sample_queue.add(valid, stages,
                             body.get('params') or {},
                             project=str(body.get('project') or '') or None)
    for it in added:
        bind_queue_entry(it)
    return jsonify({'added': len(added)})


def bind_queue_entry(it):
    """把「待启动 / 排队中」的样品挂进任务中心（外部任务）。

    为什么需要：队列条目只有在**真正开跑**时才在 samples.py 里 tm.start
    出「队列·<样品>」任务，排队期间在任务中心完全看不见 —— 用户点了入队
    却以为没生效。这里给未开跑的条目挂一张卡（queued=待手动启动，
    ready=等待名额）；一旦开跑/结束，provider 返回 None，卡片自动让位给
    真任务卡（不会重复两张）。
    """
    from Virus_Platform_Core.web.tasks import tm
    eid = it.get('id')
    sample = it.get('sample') or ''
    if not eid:
        return

    def prov(eid=eid, sample=sample):
        with sample_queue.lock:
            cur = next((x for x in sample_queue.items if x['id'] == eid), None)
            if not cur or cur.get('status') not in ('queued', 'ready'):
                return None          # 已开跑或已结束 → 交给 tm 的真任务卡
            ahead = sum(1 for x in sample_queue.items
                        if x['status'] == 'running')
            if cur.get('status') == 'queued':
                # 2026-10-11 队列 UI 随 /pipeline 页归档：API 保留但无启动入口，
                # 文案不再指向「分析流程页」。
                msg = cfg.tr(
                    f'样品 {sample} 已入队，待手动启动'
                    '（批处理队列 UI 已随「自动化分析流程」页归档）',
                    f'Sample {sample} queued (manual start required; '
                    f'queue UI archived)')
            else:
                msg = cfg.tr(f'样品 {sample} 等待批处理名额'
                             + (f'（前面有 {ahead} 个在跑）' if ahead else ''),
                             f'Sample {sample} queued')
        return {'status': 'running', 'pct': 0.0, 'stage': '排队中',
                'msg': msg, 'log': []}

    tm.attach(f'排队·{sample}', prov,
              on_cancel=lambda eid=eid: sample_queue.remove(eid),
              link='/results', key=f'queue:{eid}',
              weight='light')


@bp.route('/api/queue/<entry_id>/remove', methods=['POST'])
def api_queue_remove(entry_id):
    if not sample_queue.remove(entry_id):
        abort(400, '条目不存在或正在运行（请先取消任务）')
    return jsonify({'ok': True})


@bp.route('/api/queue/clear_finished', methods=['POST'])
def api_queue_clear():
    sample_queue.clear_finished()
    return jsonify({'ok': True})


@bp.route('/api/queue/pause', methods=['POST'])
def api_queue_pause():
    """暂停队列：正在跑的条目跑完后，已就绪条目暂不开跑，等手动恢复。"""
    sample_queue.pause()
    return jsonify({'ok': True, 'paused': True})


@bp.route('/api/queue/resume', methods=['POST'])
def api_queue_resume():
    sample_queue.resume()
    return jsonify({'ok': True, 'paused': False})


@bp.route('/api/queue/<entry_id>/start', methods=['POST'])
def api_queue_start(entry_id):
    """手动启动一条队列项（待启动/已就绪 → 排队中，由 worker 依序执行）。"""
    if not sample_queue.start(entry_id):
        abort(400, '条目不存在或状态不可启动（仅「待启动」可启动）')
    return jsonify({'ok': True})


@bp.route('/api/queue/start_all', methods=['POST'])
def api_queue_start_all():
    """全部「待启动」→「排队中」。空参数 POST 就会真启动任务，已登记
    巡检跳过表（tests/_audit_routes.py DESTRUCTIVE_POST）。"""
    return jsonify({'ok': True, 'started': sample_queue.start_all()})


def _analysis_kwargs(body):
    """analyze / pipeline run 共用的 run_analysis 参数装配。

    参数优先级：请求 body > 设置页默认值(cfg.defaults) > 代码硬编码缺省。
    """
    d = cfg.defaults

    def _val(key, fallback):
        v = body.get(key)
        if v is None or v == '':
            return d.get(key, fallback)
        return v

    def _chk(key, fallback=False):
        v = body.get(key)
        if v is None:
            return bool(d.get(key, fallback))
        return bool(v)

    def _opt_float(key):
        """可选数值门槛：空/0/非法 → None（含义"不过滤"）。

        ⑥b 层1 的 min_pident/min_qcov 用（0 门槛等于不过滤，无区分意义）。
        """
        v = body.get(key)
        if v is None or v == '':
            v = d.get(key)
        if v in (None, '', 0, '0'):
            return None
        try:
            f = float(v)
        except (TypeError, ValueError):
            return None
        return f if f > 0 else None

    def _norm_db(v):
        if not v:
            return None
        return v if os.path.isabs(str(v)) else os.path.join(PLATFORM_ROOT, v)
    def _norm_methods(v):
        """验证证据方法：接受列表/逗号串，只保留 blastx/bfvd，空则默认双路。"""
        if v is None or v == '':
            v = d.get('verify_methods') or ['blastx', 'bfvd']
        if isinstance(v, str):
            v = [x for x in v.replace(';', ',').split(',') if x.strip()]
        keep = tuple(x for x in (str(i).strip() for i in v)
                     if x in ('blastx', 'bfvd'))
        return keep or ('blastx', 'bfvd')

    # 资源钳制（2026-10-04 审查 E3）：线程 ≤ 物理核、内存 ≤ 物理内存
    # 60%——此前 POST 任意大值直通 SPAdes（-m 64 在 16GB 机器换页抖动）
    _phys_cores = os.cpu_count() or 4
    _ram_cap_gb = max(2, int((total_ram() or 16 * 1024 ** 3)
                             * 0.6 / 1024 ** 3))

    return dict(
        db_host=_norm_db(body.get('db_host')) or cfg.databases['host'],
        db_virus=_norm_db(body.get('db_virus')) or cfg.databases['virus'],
        # ②b 鉴定库：平台内库名或库目录绝对路径，空 = 默认 kv_index。
        # 不能走 _norm_db（它会把库名误拼成平台根下的路径），原样透传，
        # 由 kv_stage.resolve_kv_lib 统一解析。
        kv_lib=(str(body.get('kv_lib') or '').strip() or None),
        chunk_dir=_norm_db(body.get('chunk_dir')),
        fastp_dedup=_chk('fastp_dedup'),
        do_fq2fa=_chk('do_fq2fa', True),
        gbdraw_max=int(_val('gbdraw_max', 12) or 12),
        gbdraw_fasta=body.get('gbdraw_fasta') or None,
        gbdraw_ann=body.get('gbdraw_ann') or None,
        plot_engine=_val('plot_engine', 'auto'),
        threads=max(1, min(int(body.get('threads') or cfg.threads),
                           _phys_cores)),
        confidence=float(_val('confidence', 0) or 0),
        assembly_mode=_val('assembly_mode', 'rnaviral'),
        assembly_input=_val('assembly_input', 'virus'),
        memory_gb=min(int(_val('memory', 64) or 64), _ram_cap_gb),
        subsample=int(_val('subsample', 0) or 0),
        min_contig_len=int(_val('min_contig_len', 200) or 200),
        min_orf_aa=int(_val('min_orf_aa', 100) or 100),
        # ⑥b 层1 命中过滤门槛（百分数口径；空/0 = 不过滤，保持历史结果）
        min_pident=_opt_float('min_pident'),
        min_qcov=_opt_float('min_qcov'),
        top_n_refs=int(_val('top_n_refs', 10) or 10),
        tree_tool=_val('tree_tool', 'fasttree'),
        tree_sampling=_val('tree_sampling', 'blast'),
        ncbi_refs=body.get('ncbi_refs') or None,
        primer_mode=_val('primer_mode', 'conserved'),
        do_trim=_chk('do_trim', True),
        do_specificity=_chk('specificity'),
        force=_chk('force'),
        # ③b 候选序列验证（宿主类群 / 证据方法 / 并集-交集）
        verify_host=_val('verify_host', 'all'),
        verify_methods=_norm_methods(body.get('verify_methods')),
        verify_combine=_val('verify_combine', 'union'),
    )


@bp.route('/api/analyze', methods=['POST'])
def api_analyze():
    body = request.get_json(force=True) or {}
    if not body.get('r1'):
        abort(400, '缺少 R1')
    check_path(body['r1'], must_exist=True)
    if body.get('r2'):
        check_path(body['r2'], must_exist=True)
    from Virus_Platform_Core.pipeline import DEFAULT_ANALYZE_STAGES
    stages = body.get('stages') or list(DEFAULT_ANALYZE_STAGES)
    # 同名样品已存在 → 认磁盘真实目录名（不然会另建一个规范化后的新目录，
    # 把同一次分析拆成两个样品）；新名字则退回规范名，等同新建。
    sample = resolve_sample_name(
        body.get('sample') or _default_sample_name(body['r1']), DIRS['results'])
    for t in tm.list_all():
        if t.get('status') == 'running' and _task_matches_sample(t.get('name'), sample):
            abort(400, f'样品 {sample} 已有任务在运行（{t.get("name")}），'
                       f'请等它结束或先取消后再试')
    if not _sample_busy_register(sample):
        abort(400, f'样品 {sample} 已有任务在运行，请等它结束或先取消后再试')

    def job(log, prog, cancel):
        from Virus_Platform_Core.pipeline import run_analysis
        logger = TaskLogger(callback=log)
        try:
            return run_analysis(
                sample=sample, r1=body['r1'], r2=body.get('r2'),
                stages=stages, **_analysis_kwargs(body),
                logger=logger, progress=prog)
        finally:
            _sample_busy_release(sample)
            logger.close()

    log_dir = check_path(os.path.join(DIRS['results'], sample, 'logs'),
                         must_exist=False, in_platform=True)
    os.makedirs(log_dir, exist_ok=True)
    try:
        # on_finish：终态必调一次（含排队途中被取消，见 tasks.TaskManager.start）
        tid = tm.start(cfg.tr(f'样品分析 {sample}', f'Sample analysis {sample}'),
                       job,
                       log_file=os.path.join(
                           log_dir,
                           f'analyze_{time.strftime("%Y%m%d_%H%M%S")}.log'),
                       on_finish=lambda: _sample_busy_release(sample))
    except BaseException:
        _sample_busy_release(sample)     # 任务没起来也要放掉占位
        raise
    return jsonify({'task': tid, 'sample': sample})


def _sample_dir(sample):
    """样品名 → (真实目录名, 绝对路径)；不存在则 404。

    用 resolve_sample_name 而非纯字符转换：样品名来自 URL，可能是
    /api/samples 返回的真实目录名（含手工创建/历史/中文命名），
    字符转换会把它改成另一个名字而读不到。
    """
    from Virus_Platform_Core.config import DIRS
    s = resolve_sample_name(sample, DIRS['results'])
    return s, check_path(os.path.join(DIRS['results'], s),
                         must_exist=True, in_platform=True)


@bp.route('/api/samples')
def api_samples():
    from Virus_Platform_Core.pipeline import (pipeline_overview, load_sample_input,
                             load_project_manifest)
    from Virus_Platform_Core.config import DIRS
    out = []
    res = check_path(DIRS['results'], must_exist=False, in_platform=True)
    if os.path.isdir(res):
        for name in sorted(os.listdir(res)):
            if name.startswith('_'):      # _archive / 内部测试样品不展示
                continue
            sd = check_path(os.path.join(res, name), must_exist=False,
                            in_platform=True)
            if not os.path.isdir(sd):
                continue
            try:
                stages = pipeline_overview(sd)['stages']
            except Exception:
                stages = []
            _r1 = _r2 = None
            try:
                _r1, _r2, project = load_sample_input(sd)
            except Exception:
                project = None
            # 清单总是要读：即便 input.json 已带项目名，
            # last_status / last_run 也只能从 project.json 取。
            try:
                man = load_project_manifest(sd)
            except Exception:
                man = {}
            out.append({'name': name,
                        'project': project or man.get('project') or '',
                        'last_status': man.get('last_status') or '',
                        'last_run': man.get('last_run') or '',
                        'done': sum(1 for s in stages if s['status'] == 'done'),
                        'total': len(stages) or 7,
                        # 输入档案：样品表格展示类型（PE/SE）与文件用
                        'r1': _r1 or '', 'r2': _r2 or ''})
    return jsonify(out)


@bp.route('/api/samples/archived')
def api_samples_archived():
    """归档样品列表（results/_archive/ 下的历史项目；结果中心展示）。"""
    arch = check_path(os.path.join(DIRS['results'], '_archive'),
                      must_exist=False, in_platform=True)
    out = []
    if os.path.isdir(arch):
        for name in sorted(os.listdir(arch)):
            sd = os.path.join(arch, name)
            if not os.path.isdir(sd):
                continue
            has_report = os.path.isfile(os.path.join(sd, '07_report',
                                                     'report.html'))
            try:
                mt = max((os.path.getmtime(os.path.join(dp, fn))
                          for dp, _, fns in os.walk(sd) for fn in fns),
                         default=0)
            except OSError:
                mt = 0
            out.append({'name': name, 'has_report': has_report,
                        'mtime': mt})
    return jsonify(out)


def _default_sample_name(path):
    """从文件名推导默认样品名：去扩展名与常见测序后缀。

    Lycium_barbarum_R1.fastq.gz -> Lycium_barbarum
    NX-5_1.fq.gz                -> NX-5
    sample.fq                   -> sample
    """
    base = os.path.basename(str(path))
    base = re.sub(r'\.(fastq|fq|fasta|fa|fna)(\.gz)+$', '', base, flags=re.I)
    base = re.sub(r'\.(fastq|fq|fasta|fa|fna)$', '', base, flags=re.I)
    # 循环剥离技术后缀（Illumina 常见 _L001_R1_001 之类为多层）
    pat = re.compile(
        r'(?:[_\-.]?R[12]|_\d{2,3}|[_\-.][12]|[_\-.]L00\d)$', re.I)
    while True:
        stripped = pat.sub('', base)
        if stripped == base or not stripped:
            break
        base = stripped
    return base or 'sample'


@bp.route('/api/samples/scan_folder', methods=['POST'])
def api_samples_scan_folder():
    """扫描 FASTQ 文件夹 → 自动识别单/双端样本（零手填批量导入）。

    - 双端：R1/R2 命名自动配对（`_R1`/`.R1`/`-R1` 三种分隔，大小写不敏感）；
    - 单端：无 R1 标记的 FASTQ 各成一个样本；孤儿 R2（有 R2 没 R1）不当样本；
    - 样本名：自动取 R1 文件名去掉 R1 标记与扩展名（GQMIX.R1.fq.gz → GQMIX）；
    - 项目名：前端用文件夹名预填；
    - 顶层没有 FASTQ 时自动下探一层子目录（按样本分文件夹的存放习惯）。
    """
    body = request.get_json(force=True) or {}
    folder = str(body.get('folder') or '').strip()
    if not folder:
        abort(400, '请提供文件夹路径')
    p = check_path(folder, must_exist=True)
    if not os.path.isdir(p):
        abort(400, '路径不是目录: %s' % folder)
    from Virus_Platform_Core.consensus import (R1_RX, _r1_mate_name,
                                               find_read_pairs)
    r2_rx = re.compile(r'^(.*?)([_.\-])R2(\.(?:fq|fastq)(?:\.gz)?)$', re.I)
    fq_ext = ('.fastq', '.fq', '.fastq.gz', '.fq.gz')

    def _stem(fn):
        m = R1_RX.match(fn) or r2_rx.match(fn)
        if m:
            return m.group(1)
        return re.sub(r'\.(fq|fastq)(\.gz)?$', '', fn, flags=re.I)

    out, seen = [], set()

    def _unique(nm):
        nm = re.sub(r'[^\w\-.]+', '_', nm).strip('_-.') or 'sample'
        base, k = nm, 2
        while nm in seen:
            nm = '%s_%d' % (base, k)
            k += 1
        seen.add(nm)
        return nm

    def _scan_dir(d):
        consumed = set()
        for r1, r2 in find_read_pairs(d):
            consumed.add(os.path.realpath(r1))
            if r2:
                consumed.add(os.path.realpath(r2))
            out.append({'name': _unique(_stem(os.path.basename(r1))),
                        'r1': r1, 'r2': r2 or ''})
        for fn in sorted(os.listdir(d)):
            fp = os.path.join(d, fn)
            if (not os.path.isfile(fp)
                    or not fn.lower().endswith(fq_ext)
                    or os.path.realpath(fp) in consumed
                    or _r1_mate_name(fn)          # R1 命名已被上面收（含 R2 缺失）
                    or r2_rx.match(fn)):          # 孤儿 R2：不当独立样本
                continue
            out.append({'name': _unique(_stem(fn)), 'r1': fp, 'r2': ''})

    _scan_dir(p)
    if not out:
        # 顶层没有 FASTQ：下探一层子目录（每子目录一个/多个样本）
        for name in sorted(os.listdir(p)):
            sub = os.path.join(p, name)
            if os.path.isdir(sub):
                _scan_dir(sub)
    if not out:
        abort(400, '该文件夹（含一层子目录）下未找到 FASTQ'
                   '（.fastq/.fq/.gz；双端按 R1/R2 命名自动配对）')
    return jsonify({'folder': p, 'project': os.path.basename(p),
                    'samples': out})


@bp.route('/api/pipeline/create', methods=['POST'])
def api_pipeline_create():
    body = request.get_json(force=True) or {}
    if not body.get('r1'):
        abort(400, '缺少 R1')
    check_path(body['r1'], must_exist=True)
    if body.get('r2'):
        check_path(body['r2'], must_exist=True)
    from Virus_Platform_Core.pipeline import save_sample_input, _safe_sample_name
    from Virus_Platform_Core.config import DIRS
    typed = str(body.get('sample') or '').strip()
    sample = _safe_sample_name(typed or _default_sample_name(body['r1']))
    sd = check_path(os.path.join(DIRS['results'], sample),
                    must_exist=False, in_platform=True)
    if os.path.isdir(sd) and any(os.scandir(sd)):
        # 撞名时把"输入名 → 实际目录名"讲清楚，否则用户填「样品A」却被告知
        # 「A 已存在」，完全不知道名字被规范化过（这类静默改名很难自查）。
        if typed and typed != sample:
            abort(400, f'样品名「{typed}」会被规范化为「{sample}」，'
                       f'而「{sample}」已存在；请换一个样品名')
        abort(400, f'样品 {sample} 已存在，请换一个样品名')
    os.makedirs(sd, exist_ok=True)
    # 可选：创建时立即截取子样本作为管道输入（大样品先小规模验证）
    sub = int(body.get('subsample', 0) or 0)
    in1, in2 = body['r1'], body.get('r2')
    if sub > 0:
        from Virus_Platform_Core.pipeline import subsample_fastq
        in1, in2 = subsample_fastq(in1, in2, os.path.join(sd, '00_prep'), sub)
    save_sample_input(sd, in1, in2, sample,
                      project=str(body.get('project') or '') or None)
    # note：名字被规范化过时明确回给前端提示（非 ASCII/空格会被折叠成 '_'，
    # 因为样品名要当目录名进 SPAdes/BLAST 等不支持中文路径的工具）。
    note = ''
    if typed and typed != sample:
        note = (f'样品名「{typed}」已按文件名安全规则登记为「{sample}」'
                f'（目录名不能含中文/空格等字符）')
    return jsonify({'sample': sample, 'typed': typed, 'note': note})


@bp.route('/api/samples/<sample>/rename', methods=['POST'])
def api_sample_rename(sample):
    """样品改名：results/<旧名> 整目录重命名，分析产物随目录一起带走。

    拦截：① 有任务在运行（_sample_busy）；② 在批处理队列排队/运行中。
    新名走与创建同一套 _safe_sample_name 规范化（非 ASCII/空格折叠成
    '_'、截断 50，因为目录名要进 SPAdes/BLAST 等不支持中文路径的工具），
    规范化结果与输入不同时在 note 里如实告知。input.json 的 sample
    字段同步改写，保证档案与目录一致。"""
    from Virus_Platform_Core.pipeline import _safe_sample_name
    s, sd = _sample_dir(sample)
    body = request.get_json(force=True) or {}
    typed = str(body.get('new_name') or '').strip()
    if not typed:
        abort(400, '新样品名不能为空')
    if _sample_busy(s):
        abort(400, f'样品 {s} 有任务在运行，请先取消再改名')
    for it in sample_queue.items:
        if (it.get('sample') == s
                and it.get('status') in ('queued', 'running')):
            abort(400, f'样品 {s} 在批处理队列中（排队/运行），'
                       f'请先在队列里移除后再改名')
    new = _safe_sample_name(typed)
    if new == s:
        abort(400, '新名字与当前相同（或规范化后相同），无需改名')
    if new.startswith('_'):
        abort(400, f'样品名「{typed}」规范化为「{new}」——下划线开头是'
                   f'内部保留前缀，请换一个名字')
    target = check_path(os.path.join(DIRS['results'], new),
                        must_exist=False, in_platform=True)
    if os.path.isdir(target) and any(os.scandir(target)):
        if typed != new:
            abort(400, f'样品名「{typed}」会被规范化为「{new}」，'
                       f'而「{new}」已存在；请换一个样品名')
        abort(400, f'样品 {new} 已存在，请换一个样品名')
    try:
        os.rename(sd, target)
    except OSError as e:
        abort(400, f'改名失败: {e}')
    # input.json 的 sample 字段与目录名保持一致（历史档案可读不误导）
    try:
        inp = os.path.join(target, '00_prep', 'input.json')
        if os.path.isfile(inp):
            with safe_open(inp) as f:
                d = json.load(f)
            if d.get('sample') == s:
                d['sample'] = new
                with safe_open(inp, 'wt') as f:
                    json.dump(d, f, ensure_ascii=False, indent=2)
    except (OSError, ValueError):
        pass
    note = '' if typed == new else (
        f'样品名「{typed}」已按文件名安全规则登记为「{new}」'
        f'（目录名不能含中文/空格等字符）')
    return jsonify({'sample': new, 'note': note})


@bp.route('/api/pipeline/<sample>')
def api_pipeline(sample):
    from Virus_Platform_Core.pipeline import pipeline_overview, load_sample_input
    s, sd = _sample_dir(sample)
    r1, r2, project = load_sample_input(sd)
    lang = request.args.get('lang') or None
    return jsonify({'sample': s, 'r1': r1, 'r2': r2, 'project': project or '',
                    **pipeline_overview(sd, lang)})


@bp.route('/api/pipeline/<sample>/run', methods=['POST'])
def api_pipeline_run(sample):
    body = request.get_json(force=True) or {}
    s, sd = _sample_dir(sample)
    # 防重复：该样品已有运行中任务时拒绝（避免并发分类写爆磁盘）
    for t in tm.list_all():
        if t.get('status') == 'running' and _task_matches_sample(t.get('name'), s):
            abort(400, f'样品 {s} 已有任务在运行（{t.get("name")}），'
                       f'请等它结束或先取消后再试')
    from Virus_Platform_Core.pipeline import load_sample_input, run_analysis, STAGE_ORDER
    r1, r2, _proj = load_sample_input(sd)
    if not r1:
        abort(400, '样品缺少输入记录（input.json），请重新创建样品')
    check_path(r1, must_exist=True)
    stages = body.get('stages') or ([body['stage']] if body.get('stage') else [])
    stages = [st for st in stages if st in STAGE_ORDER]
    # fastp 未安装时剔除（可选步骤）
    if 'fastp' in stages:
        from Virus_Platform_Core.preprocess import fastp_available
        if not fastp_available():
            stages = [st for st in stages if st != 'fastp']
    if not stages:
        # 未指定阶段 = 「依次运行剩余步骤」：所有未完成且可用的阶段
        try:
            from Virus_Platform_Core.pipeline import pipeline_overview
            ov = pipeline_overview(sd)['stages']
            done = {s['stage'] for s in ov if s['status'] == 'done'}
            stages = [s['stage'] for s in ov
                      if s['stage'] in STAGE_ORDER
                      and s['stage'] not in done
                      and s['status'] != 'unavailable']
        except Exception:
            stages = []
    if not stages:
        abort(400, '没有需要运行的阶段（全部已完成；重跑请勾选「强制重跑」'
                   '并指定具体步骤）')
    if not _sample_busy_register(s):
        abort(400, f'样品 {s} 已有任务在运行，请等它结束或先取消后再试')

    def job(log, prog, cancel):
        logger = TaskLogger(callback=log)
        try:
            run_analysis(sample=s, r1=r1, r2=r2, stages=stages,
                         **_analysis_kwargs(body), logger=logger, progress=prog)
            return sd
        finally:
            _sample_busy_release(s)
            logger.close()

    from Virus_Platform_Core.pipeline import stage_name
    label = ' → '.join(stage_name(st, cfg.lang) for st in stages)
    log_dir = check_path(os.path.join(sd, 'logs'), must_exist=False,
                         in_platform=True)
    os.makedirs(log_dir, exist_ok=True)
    try:
        # on_finish：终态必调一次（含排队途中被取消，见 tasks.TaskManager.start）
        tid = tm.start(f'{s} · {label}', job,
                       log_file=os.path.join(
                           log_dir,
                           f'run_{time.strftime("%Y%m%d_%H%M%S")}.log'),
                       on_finish=lambda: _sample_busy_release(s))
    except BaseException:
        _sample_busy_release(s)          # 任务没起来也要放掉占位
        raise
    return jsonify({'task': tid, 'sample': s})


@bp.route('/api/open_report_dir/<sample>')
def api_open_report_dir(sample):
    safe = _safe_sample(sample)
    d = check_path(os.path.join(DIRS['results'], safe), must_exist=True,
                   in_platform=True)
    os.startfile(check_path(d, must_exist=True, in_platform=True))
    return jsonify({'ok': True})


@bp.route('/api/samples/<sample>/<path:rel>')
def api_sample_file(sample, rel):
    """结果文件下载（严格限制在样品目录内）。"""
    safe = _safe_sample(sample)
    # rel 传绝对路径（如 C:/.../run/ai_local.json）时 os.path.join 会整段
    # 替换 base；check_path(in_platform=True) 只保证「平台内」，拦不住读
    # 平台内其它文件。归一后强制落在样品目录内（样板见 meta_files_download）。
    base = os.path.abspath(check_path(os.path.join(DIRS['results'], safe),
                                      must_exist=True, in_platform=True))
    p = check_path(os.path.join(base, rel), must_exist=True, in_platform=True)
    if not p.startswith(base + os.sep):
        abort(400, '非法路径')
    return send_file(p, as_attachment=request.args.get('dl') == '1')


def _sample_busy(safe):
    """样品是否有正在运行的任务（删除/清除结果前拦截）。

    登记表优先：任务名正则匹配不到的场景（工具页直跑某阶段）也能拦住。"""
    with _RUNNING_SAMPLES_LOCK:
        if safe in _RUNNING_SAMPLES:
            return True
    for t in tm.list_all():
        if t.get('status') == 'running' and _task_matches_sample(
                t.get('name'), safe):
            return True
    return False


@bp.route('/api/samples/<sample>/clear', methods=['POST'])
def api_sample_clear(sample):
    """清除样品的全部分析结果文件，保留样品登记（00_prep/input.json）。

    与 delete 的区别：样品仍留在样品列表中，输入档案不丢，
    之后再跑管道即从第一步重新分析。样品有运行中任务时拒绝。
    """
    safe = _safe_sample(sample)
    if _sample_busy(safe):
        abort(400, f'样品 {safe} 有任务在运行，请先取消再清除')
    d = check_path(os.path.join(DIRS['results'], safe), must_exist=True,
                   in_platform=True)
    # 双保险：必须是 results/ 的直接子目录
    if os.path.dirname(os.path.abspath(d)) != os.path.abspath(DIRS['results']):
        abort(400, '仅允许操作样品目录')
    prep = os.path.join(d, '00_prep')
    # 顶层两个档案文件与结果无关，清除结果时保留：
    #   input.json   —— 样品输入登记（在 00_prep 内）
    #   project.json —— 项目清单（项目身份/输入/参数）
    keep_top = {'project.json'}
    # 二次检查：上面的 _sample_busy 与这里之间有读档案/建目录等磁盘操作，
    # 队列或手动任务可能刚好在这段窗口内启动。Windows 上 rmtree 遇到
    # 正在写入的文件会抛 PermissionError 并留下半删目录，所以必须重新确认。
    if _sample_busy(safe):
        abort(400, f'样品 {safe} 刚刚启动了任务，请先取消再清除')
    for name in os.listdir(d):
        p = os.path.join(d, name)
        if name in keep_top:
            continue
        if os.path.abspath(p) == os.path.abspath(prep):
            # 输入档案目录：只保留 input.json（fastp/fq2fa/子采样产物一并清）
            for f in os.listdir(prep):
                if f == 'input.json':
                    continue
                fp = os.path.join(prep, f)
                try:
                    if os.path.isdir(fp):
                        shutil.rmtree(fp)
                    else:
                        os.remove(fp)
                except PermissionError:
                    abort(400, f'文件 {f} 正在被任务占用，无法删除，'
                               f'请等待任务结束或先取消')
        elif os.path.isdir(p):
            try:
                shutil.rmtree(p)
            except PermissionError:
                abort(400, f'目录 {name} 正在被任务占用，无法删除，'
                           f'请等待任务结束或先取消')
        else:
            try:
                os.remove(p)
            except PermissionError:
                abort(400, f'文件 {name} 正在被任务占用，无法删除，'
                           f'请等待任务结束或先取消')
    # 结果已清空，清单里的结果类字段必须一并归零，
    # 否则 project.json 会声称还有已完成的阶段（与卡片状态矛盾）。
    # 仅保留项目身份：project / input / params / schema / created。
    try:
        from Virus_Platform_Core.pipeline import load_project_manifest, save_project_manifest
        _m = load_project_manifest(d)
        if _m:
            for _k in ('stages_done', 'runs', 'last_status', 'last_stage',
                       'last_duration_s', 'last_run'):
                _m.pop(_k, None)
            save_project_manifest(d, _m)
    except Exception as e:
        # 清单归零失败会让 project.json 与已清空的结果目录矛盾（卡片仍显示阶段已完成）
        _log.warning('样品清单归零失败 %s: %s', safe, e)
    return jsonify({'ok': True})


@bp.route('/api/samples/<sample>/delete', methods=['POST'])
def api_sample_delete(sample):
    """删除样品：结果目录（含 00_prep 输入档案与全部产物）一并移除。

    仅允许平台 results/ 内的一级样品目录；样品有运行中任务时拒绝。
    只想清空结果、保留样品时用 /clear。
    """
    safe = _safe_sample(sample)
    if _sample_busy(safe):
        abort(400, f'样品 {safe} 有任务在运行，请先取消再删除')
    d = check_path(os.path.join(DIRS['results'], safe), must_exist=True,
                   in_platform=True)
    # 双保险：必须是 results/ 的直接子目录
    if os.path.dirname(os.path.abspath(d)) != os.path.abspath(DIRS['results']):
        abort(400, '仅允许删除样品目录')
    # 二次检查：上面的 _sample_busy 与这里之间有路径校验等操作，任务可能
    # 刚好在这段窗口内启动。Windows 上 rmtree 遇到正在写入的文件会抛
    # PermissionError 并留下半删目录。
    if _sample_busy(safe):
        abort(400, f'样品 {safe} 刚刚启动了任务，请先取消再删除')
    try:
        shutil.rmtree(d)
    except PermissionError:
        abort(400, f'样品目录 {safe} 正在被任务占用，无法删除，'
                   f'请等待任务结束或先取消')
    return jsonify({'ok': True})


@bp.route('/api/sample_files/<sample>')
def api_sample_files(sample):
    """列出样品结果文件树（相对路径）。"""
    safe = _safe_sample(sample)
    root = check_path(os.path.join(DIRS['results'], safe), must_exist=True,
                      in_platform=True)
    out = []
    for dirpath, _dirs, fnames in os.walk(root):
        for fn in fnames:
            if fn.startswith('.') or fn == 'plotly.min.js':
                continue
            full = check_path(os.path.join(dirpath, fn), must_exist=False,
                              in_platform=True)
            if not os.path.isfile(full):
                continue
            relp = os.path.relpath(full, root)
            out.append({'path': relp.replace(os.sep, '/'),
                        'size': fmt_size(os.path.getsize(full))})
    out.sort(key=lambda x: x['path'])
    return jsonify({'files': out})


# 平台启动时把「已在排队」的条目补挂到任务中心：SampleQueue 会从
# tasks/queue.json 恢复条目，但那时没有 attach 过 —— 不补挂的话，
# 重启后任务中心看不到这些排队项（下载批次由 DownloadManager.add_listener
# 走同一套补挂）。
for _it in list(sample_queue.items):
    if _it.get('status') == 'queued':
        try:
            bind_queue_entry(_it)
        except Exception:
            pass
