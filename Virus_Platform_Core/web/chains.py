# -*- coding: utf-8 -*-
"""链式一键分析：把多个任务工厂串成一次运行（自 app.py 拆出）。

设计：每步仍调用 Virus_Platform_Core/web/tool_jobs.py 的现有 _tool_job_*，各建
自己的 run_dir，一步不改；汇总目录只放「符号链接 + 清单」，
把一次多步运行的结果收拢成一份。
"""
import json
import os
import re
import time
import uuid
from pathlib import Path

from flask import abort

from Virus_Platform_Core.config import PLATFORM_ROOT
from Virus_Platform_Core.utils import (TaskLogger, check_path, safe_open)
from Virus_Platform_Core.web.state import cfg, tool_runs_root as _tool_runs_root
from Virus_Platform_Core.web.tool_jobs import (
    _tool_job_identify,
    _tool_job_assemble,
    _tool_job_contigs,
    _tool_job_verify,
    _tool_job_kvsuite,
    _tool_job_convert,
    _tool_job_fastp,
    _tool_job_hostremoval,
    _tool_job_orf,
    _tool_job_orfa,
    _tool_job_genoplot,
)


_CHAIN_LINKABLE = ('.fasta', '.fa', '.fna', '.tsv', '.json', '.txt',
                   '.vcf', '.gz', '.html', '.pdf', '.png')


def _chain_link(dst, src, logger=None, target_is_dir=None):
    """在汇总目录建链接指向 src。

    Windows 下 os.symlink 需开发者模式/管理员；不可用时按序降级：
    目录 → 硬链接不可用则跳过并记录；文件 → 硬链接；再不行则跳过。
    返回实际使用的方式（'symlink' / 'hardlink' / None），失败不抛异常。
    """
    if not os.path.exists(src):
        return None
    if target_is_dir is None:
        target_is_dir = os.path.isdir(src)
    try:
        os.symlink(src, dst, target_is_directory=bool(target_is_dir))
        return 'symlink'
    except (OSError, NotImplementedError, AttributeError):
        pass
    if not target_is_dir:
        try:
            os.link(src, dst)
            return 'hardlink'
        except OSError:
            pass
    if logger:
        logger.log(f'链接创建跳过（源或权限不可用）: {os.path.basename(dst)}')
    return None


def _chain_write_manifest(chain_dir, meta):
    p = os.path.join(chain_dir, '_chain.json')
    with safe_open(p, 'wt') as f:
        json.dump(meta, f, ensure_ascii=False, indent=2)
    return p


def _chain_run_step(chain_dir, root, step_key, step_name, job_factory,
                    p, threads, db_virus, logger, prog, prog_lo, prog_hi,
                    cancel=None, collect=None, optional=()):
    """跑链中一步：新建 run_dir、构造 ctx、调用 job，返回 (run_name, result)。

    prog 为阶段名；进度在 [prog_lo, prog_hi] 区间内线性映射。
    collect: 该步要链接进汇总目录的产物文件名列表（可选）。
    optional: 这些 key 就算下游 job 用 ctx.req 取也允许为空（返回 None），
              用于兼容现有 job 把可选输入写成必填的情况。
    """
    from types import SimpleNamespace
    ts = time.strftime('%Y%m%d_%H%M%S')
    # 秒级时间戳会撞车：同一秒内两次链（双击/多标签页）会复用同一步
    # 运行目录并互相覆盖产物。追加 4 位随机后缀（口径同 tools_api 的
    # /api/tool/run 建运行名）。
    run_name = f'{step_key}_{ts}_{uuid.uuid4().hex[:4]}'
    run_dir = check_path(os.path.join(root, run_name), must_exist=False,
                         in_platform=True)
    os.makedirs(run_dir, exist_ok=True)

    def _input_abs(v):
        return check_path(v if os.path.isabs(v)
                          else os.path.join(PLATFORM_ROOT, v),
                          must_exist=True)

    def _req(key, what):
        v = (p.get(key) or '').strip()
        if not v:
            if key in optional:
                return None
            abort(400, f'请选择{what}')
        return _input_abs(v)

    def _opt(key):
        v = (p.get(key) or '').strip()
        return _input_abs(v) if v else None

    ctx = SimpleNamespace(p=p, run_dir=run_dir, threads=threads,
                          db_virus=db_virus, req=_req, opt=_opt)
    job = job_factory(ctx)

    def _sub_prog(_stage, frac, msg):
        f = prog_lo + max(0.0, min(1.0, float(frac))) * (prog_hi - prog_lo)
        prog(step_key, f, msg)

    prog(step_key, prog_lo, f'{step_name}：开始')
    res = job(lambda m: logger.log(m), _sub_prog, cancel)
    logger.log(f'{step_name} 完成: run={run_name}')

    # 目录链接（整步产物）
    _chain_link(os.path.join(chain_dir, f'step_{step_key}'), run_dir,
                logger=logger, target_is_dir=True)
    # 关键文件链接
    linked = []
    for fn in (collect or ()):
        src = os.path.join(run_dir, fn)
        if os.path.isfile(src):
            how = _chain_link(os.path.join(chain_dir, f'{step_key}__{fn}'),
                              src, logger=logger, target_is_dir=False)
            if how:
                linked.append(fn)
    return run_name, res, linked


def _tool_job_virchain(ctx):
    """病毒识别和分类分析 · 一键：②提取序列组装 → ③组装结果再鉴定 → ④候选序列验证。

    输入固定为去宿主序列（FASTA，可单/双端）。①（识别分类与提取）为可选：
      跳过 ①：去宿主序列 ────────────────────→ ②组装 → ③再鉴定 → ④验证
      运行 ①：去宿主序列 → ①提取病毒序列 → ②组装 → ③再鉴定 → ④验证
    即两条路都到 ②，区别只在 ② 的输入（原始去宿主序列 vs ① 提取的病毒序列）。
    组装（SPAdes）接受 FASTA，或自动 --only-assembler 跳过纠错。
    每步各建 run_dir，汇总目录放符号链接 + _chain.json。
    """
    root = _tool_runs_root()
    inp = ctx.req('input', '去宿主序列（FASTA）')
    inp2 = ctx.opt('input2')
    # 可选：先跑 ① 提取病毒序列，作为 ② 的输入
    run_identify = bool(ctx.p.get('run_identify'))
    # ① 的输入类型：双端 / 单端 / FASTA
    if inp2:
        id_type = 'pe'
    else:
        id_type = ('fasta' if inp.lower().endswith(('.fa', '.fasta', '.fna'))
                   else 'single')

    do_verify = ctx.p.get('do_verify')
    do_verify = True if do_verify is None else bool(do_verify)
    asm_mode = (ctx.p.get('mode') or 'rnaviral').strip()
    mem = int(ctx.p.get('memory') or 64)
    min_len = int(ctx.p.get('min_len') or 200)
    conf = float(ctx.p.get('confidence') or 0)
    host = (ctx.p.get('host') or 'all').strip() or 'all'
    combine = (ctx.p.get('combine') or 'union').strip()
    methods = ctx.p.get('methods') or ['blastx', 'bfvd']
    if isinstance(methods, str):
        methods = [m.strip() for m in methods.split(',') if m.strip()]
    methods = [m for m in methods if m in ('blastx', 'bfvd')]

    def job(log, prog, cancel):
        logger = TaskLogger(callback=log)
        chain_dir = ctx.run_dir
        chain_name = os.path.basename(chain_dir)
        meta = {'chain': 'virchain', 'run': chain_name,
                'label': '病毒识别和分类分析 一键流程',
                'run_identify': run_identify,
                'created': time.strftime('%Y-%m-%d %H:%M:%S'), 'steps': []}

        # 进度权重：identify(若开) 0.20 / assembly 0.35 / contigs 0.25 / verify 0.20
        w_id = 0.20 if run_identify else 0.0
        w_asm = 0.35
        w_ctg = 0.25
        w_vfy = 0.20 if do_verify else 0.0
        tot = w_id + w_asm + w_ctg + w_vfy

        def _band(lo, w):
            return (lo, lo + (w / tot if tot else 0.0))

        cur = 0.0

        # ---- ①（可选）病毒识别分类与提取：产物作为 ② 的输入 ----
        asm_in1, asm_in2 = inp, inp2
        if run_identify:
            lo, hi = _band(cur, w_id)
            p1 = {'input': inp, 'input_type': id_type, 'confidence': conf}
            if inp2:
                p1['input2'] = inp2
            run1, res1, _ = _chain_run_step(
                chain_dir, root, 'identify', '① 病毒识别分类与提取',
                _tool_job_identify, p1, ctx.threads, ctx.db_virus,
                logger, prog, lo, hi, cancel=cancel,
                collect=('viral_ids.tsv',))
            meta['steps'].append({'key': 'identify', 'run': run1,
                                  'title': '① 病毒识别分类与提取',
                                  'summary': {k: v for k, v in (res1 or {}).items()
                                              if not isinstance(v, (list, dict))}})
            # 取 ① 实际产出的病毒序列作为 ② 的输入（配对输入会带 _1/_2 后缀）
            id_dir = os.path.join(root, run1)
            vseqs = []
            if os.path.isdir(id_dir):
                for fn in sorted(os.listdir(id_dir)):
                    if fn.startswith('viral_sequences') and fn.endswith('.fasta'):
                        fp = os.path.join(id_dir, fn)
                        if os.path.getsize(fp) > 0:
                            vseqs.append(fp)
            if not vseqs:
                raise RuntimeError('① 未提取到病毒序列（viral_sequences*.fasta 为空），'
                                   '可取消「先跑 ①」开关直接用去宿主序列组装')
            # 链接进汇总目录
            for fp in vseqs:
                _chain_link(os.path.join(chain_dir,
                                         'identify__' + os.path.basename(fp)),
                            fp, logger=logger, target_is_dir=False)
            asm_in1 = vseqs[0]
            asm_in2 = vseqs[1] if len(vseqs) > 1 else None
            meta['steps'][-1]['output'] = [os.path.basename(x) for x in vseqs]
            logger.log('② 组装输入改用 ① 提取的病毒序列: '
                       + ', '.join(os.path.basename(x) for x in vseqs))
            cur = hi

        # ---- ② 提取序列组装（输入 = 去宿主序列 或 ① 的提取序列） ----
        lo, hi = _band(cur, w_asm)
        p2 = {'r1': asm_in1, 'mode': asm_mode, 'memory': mem,
              'min_len': min_len}
        if asm_in2:
            p2['r2'] = asm_in2
        run2, res2, _ = _chain_run_step(
            chain_dir, root, 'assembly', '② 提取序列组装',
            _tool_job_assemble, p2, ctx.threads, ctx.db_virus,
            logger, prog, lo, hi, cancel=cancel,
            collect=('contigs.filtered.fasta',), optional=('r2',))
        contigs = (res2 or {}).get('contigs')
        meta['steps'].append({'key': 'assembly', 'run': run2,
                              'title': '② 提取序列组装',
                              'summary': {k: v for k, v in (res2 or {}).items()
                                          if not isinstance(v, (list, dict))}})
        if not contigs or not os.path.isfile(contigs):
            raise RuntimeError('② 组装未产出 contigs，后续步骤无法继续')
        cur = hi

        # ---- ③ 组装结果再鉴定 ----
        lo, hi = _band(cur, w_ctg)
        p3 = {'contigs': contigs, 'min_len': min_len, 'confidence': conf}
        run3, res3, _ = _chain_run_step(
            chain_dir, root, 'contigs', '③ 组装结果再鉴定',
            _tool_job_contigs, p3, ctx.threads, ctx.db_virus,
            logger, prog, lo, hi, cancel=cancel,
            collect=('virus_classification.tsv', 'viral_contigs.fasta'))
        viral_fa = (res3 or {}).get('viral_fasta')
        meta['steps'].append({'key': 'contigs', 'run': run3,
                              'title': '③ 组装结果再鉴定',
                              'summary': {k: v for k, v in (res3 or {}).items()
                                          if not isinstance(v, (list, dict))}})
        cur = hi

        # ---- ④ 候选序列验证（依赖 ③ 的产物，直接读其 run_dir） ----
        if do_verify:
            if not viral_fa or not os.path.isfile(viral_fa):
                logger.log('③ 未提取到病毒 contig，跳过 ④ 候选序列验证')
                meta['steps'].append({'key': 'verify', 'run': None,
                                      'title': '④ 候选序列验证',
                                      'skipped': '③ 未产出 viral_contigs.fasta'})
            else:
                lo, hi = _band(cur, w_vfy)
                # verify 的 run 模式需要源 run 名；直接传 ③ 的 run
                p4 = {'run': run3, 'host': host, 'combine': combine,
                      'methods': methods}
                run4, res4, _ = _chain_run_step(
                    chain_dir, root, 'verify', '④ 候选序列验证',
                    _tool_job_verify, p4, ctx.threads, ctx.db_virus,
                    logger, prog, lo, hi, cancel=cancel,
                    collect=('verify_summary.json',))
                meta['steps'].append({'key': 'verify', 'run': run4,
                                      'title': '④ 候选序列验证',
                                      'summary': {k: v for k, v in (res4 or {}).items()
                                                  if not isinstance(v, (list, dict))}})
                cur = hi

        prog('chain', 1.0, '一键流程完成')
        logger.close()
        meta['done'] = True
        _chain_write_manifest(chain_dir, meta)
        return {'chain': chain_name, 'run': chain_name,
                'n_steps': len(meta['steps']), 'steps': meta['steps'],
                'chain_json': os.path.join(chain_dir, '_chain.json')}
    return job


def _tool_job_kvchain(ctx):
    """病毒定量与共识 · 一键：调用 kvsuite 全流程（鉴定→过滤→共识→绘图→变异）。

    kvsuite 卡里的 stage 下拉默认已是 all，本 job 把「一次跑完并收拢产物」
    做成显式入口：建汇总目录 + 链接 kvsuite run_dir + 关键产物。
    """
    root = _tool_runs_root()

    def job(log, prog, cancel):
        logger = TaskLogger(callback=log)
        chain_dir = ctx.run_dir
        chain_name = os.path.basename(chain_dir)
        # kvsuite 的 stage 固定为 all（一键语义）
        p = dict(ctx.p)
        p['stage'] = 'all'
        meta = {'chain': 'kvchain', 'run': chain_name,
                'label': '病毒定量与共识 一键流程',
                'created': time.strftime('%Y-%m-%d %H:%M:%S'), 'steps': []}

        run1, res1, _ = _chain_run_step(
            chain_dir, root, 'kvsuite', '已知病毒识别与定量（全流程）',
            _tool_job_kvsuite, p, ctx.threads, ctx.db_virus,
            logger, prog, 0.0, 1.0, cancel=cancel,
            collect=('run.log',))
        meta['steps'].append({'key': 'kvsuite', 'run': run1,
                              'title': '已知病毒识别与定量（全流程）',
                              'summary': {k: v for k, v in (res1 or {}).items()
                                          if not isinstance(v, (list, dict))}})
        prog('chain', 1.0, '一键流程完成')
        logger.close()
        meta['done'] = True
        _chain_write_manifest(chain_dir, meta)
        return {'chain': chain_name, 'run': chain_name,
                'n_steps': 1, 'steps': meta['steps'],
                'kvsuite_run': run1,
                'chain_json': os.path.join(chain_dir, '_chain.json')}
    return job


def _tool_job_prepchain(ctx):
    """样本处理 · 一键前处理：转换 → 质控 → 转换(FASTA) → 宿主去除（输入自适应）。

    .sra    ：① SRA→FASTQ → ② fastp 质控 → ③ FASTQ→FASTA → ④ 宿主去除
    .fastq  ：① fastp 质控 → ② FASTQ→FASTA → ③ 宿主去除（跳过 SRA 解码）
    FASTA   ：不接受——宿主去除要从 FASTQ 提取 kept reads，FASTA 直进会把
              宿主 read 一起带到下游组装/鉴定。

    自动降段（其余步骤照常，结果与 _chain.json 里注明原因）：
      - 宿主库未就绪 → 跳过宿主去除（「数据库构建」页构建后可单独跑）；
      - fastp 未安装 → 跳过质控（转换/宿主去除改用原始 FASTQ）。
    转换出的 FASTA 作为宿主去除的分类加速输入（remove_host.classify_r1/r2，
    与样品管道 ⓪b→① 同语义）；kept reads 仍从 FASTQ 提取，kept_*.fastq.gz
    可直接喂 virchain（②→③→④）或样品管道。
    """
    root = _tool_runs_root()
    inp = ctx.req('input', '输入文件（.sra 或 FASTQ）')
    inp2 = ctx.opt('input2')
    is_sra = str(inp).lower().endswith('.sra')
    if is_sra:
        if inp2:
            abort(400, '.sra 输入不需要 R2（单/双端由 sracha 解码结果决定）')
    elif not str(inp).lower().endswith(('.fastq', '.fq', '.fastq.gz', '.fq.gz')):
        abort(400, '输入需为 .sra 或 FASTQ(.fastq/.fq[.gz])；'
                   'FASTA 不能直接宿主去除（kept reads 必须从 FASTQ 提取）')
    conf = float(ctx.p.get('confidence') or 0)
    dedup = bool(ctx.p.get('dedup'))
    db_host = ctx.opt('db') or cfg.databases['host']

    def job(log, prog, cancel):
        from Virus_Platform_Core.kunpeng import db_ready
        from Virus_Platform_Core.preprocess import fastp_available
        logger = TaskLogger(callback=log)
        chain_dir = ctx.run_dir
        chain_name = os.path.basename(chain_dir)
        host_ok = bool(db_ready(db_host))
        fastp_ok = bool(fastp_available())
        meta = {'chain': 'prepchain', 'run': chain_name,
                'label': '样本处理 一键前处理',
                'input': inp, 'input2': inp2 or '', 'is_sra': is_sra,
                'host_db_ready': host_ok, 'fastp_available': fastp_ok,
                'created': time.strftime('%Y-%m-%d %H:%M:%S'), 'steps': []}
        if not host_ok:
            logger.log('宿主库未就绪 → 自动跳过宿主去除'
                       '（「数据库构建」页构建宿主库后可单独运行）')
        if not fastp_ok:
            logger.log('fastp 未安装 → 自动跳过质控'
                       '（可选步骤；后续转换/宿主去除改用原始 FASTQ）')

        # 计划与进度权重（缺哪个环节就把它剔出计划，权重重新归一）
        plan = []
        if is_sra:
            plan.append(('sra2fq', 'SRA→FASTQ 转换', 0.15))
        if fastp_ok:
            plan.append(('fastp', 'fastp 质控', 0.25))
        plan.append(('fq2fa', 'FASTQ→FASTA 转换', 0.20))
        if host_ok:
            plan.append(('host', '宿主去除（kunpeng 分类）', 0.40))
        tot = sum(w for _, _, w in plan) or 1.0
        cur = 0.0

        def _band(w):
            nonlocal cur
            lo, hi = cur, cur + w / tot
            cur = hi
            return lo, hi

        fq1, fq2 = (inp, inp2) if not is_sra else (None, None)
        kept, fasta = [], []

        # ---- SRA→FASTQ（仅 .sra 输入）----
        if is_sra:
            lo, hi = _band(0.15)
            run1, res1, _ = _chain_run_step(
                chain_dir, root, 'sra2fq', 'SRA→FASTQ 转换',
                _tool_job_convert, {'input': inp, 'target': 'fastq'},
                ctx.threads, ctx.db_virus, logger, prog, lo, hi, cancel=cancel)
            outs = [os.path.join(root, run1, f)
                    for f in (res1 or {}).get('files') or []]
            pair2 = [f for f in outs
                     if re.search(r'_2\.(fastq|fq)\.gz$', os.path.basename(f),
                                  re.I)]
            rest = [f for f in outs if f not in pair2]
            fq1 = rest[0] if rest else None
            fq2 = pair2[0] if pair2 else None
            if not fq1:
                raise RuntimeError('SRA→FASTQ 未产出 FASTQ（见运行日志）')
            meta['steps'].append({
                'key': 'sra2fq', 'run': run1, 'title': 'SRA→FASTQ 转换',
                'output': [os.path.basename(x) for x in outs],
                'pe': bool(fq2)})
            logger.log(f'sracha 解码为{"双端" if fq2 else "单端"} FASTQ')

        # ---- fastp 质控（未安装则跳过，沿用原始 FASTQ）----
        if fastp_ok:
            lo, hi = _band(0.25)
            p2 = {'r1': fq1, 'dedup': dedup}
            if fq2:
                p2['r2'] = fq2
            run2, res2, _ = _chain_run_step(
                chain_dir, root, 'fastp', 'fastp 质控', _tool_job_fastp, p2,
                ctx.threads, ctx.db_virus, logger, prog, lo, hi,
                cancel=cancel, optional=('r2',))
            prep2 = os.path.join(root, run2, '00_prep')
            qc1 = os.path.join(prep2, 'fastp_R1.fastq.gz')
            qc2 = os.path.join(prep2, 'fastp_R2.fastq.gz') if fq2 else None
            missing = [os.path.basename(x) for x in (qc1, qc2)
                       if x and not os.path.isfile(x)]
            if missing:
                raise RuntimeError(f'质控未产出 {", ".join(missing)}')
            fq1, fq2 = qc1, (qc2 if fq2 else None)
            meta['steps'].append({
                'key': 'fastp', 'run': run2, 'title': 'fastp 质控',
                'summary': {k: v for k, v in (res2 or {}).items()
                            if not isinstance(v, (list, dict))}})

        # ---- FASTQ→FASTA（分类加速输入；双端逐端转换）----
        lo, hi = _band(0.20)
        mates = [('fq2fa', 'FASTQ→FASTA（R1）', fq1)]
        if fq2:
            mates.append(('fq2fa_r2', 'FASTQ→FASTA（R2）', fq2))
        fa1 = fa2 = None
        for j, (key, nm, src) in enumerate(mates):
            seg_lo = lo + (hi - lo) * j / len(mates)
            seg_hi = lo + (hi - lo) * (j + 1) / len(mates)
            rr, rs, _ = _chain_run_step(
                chain_dir, root, key, nm, _tool_job_convert,
                {'input': src, 'target': 'fasta'},
                ctx.threads, ctx.db_virus, logger, prog, seg_lo, seg_hi,
                cancel=cancel)
            rels = (rs or {}).get('files') or []
            if not rels:
                raise RuntimeError(f'{nm} 未产出 FASTA（见运行日志）')
            f = os.path.join(root, rr, rels[0])
            if key == 'fq2fa':
                fa1 = f
            else:
                fa2 = f
            fasta.append(f)
            _chain_link(os.path.join(chain_dir, os.path.basename(f)), f,
                        logger=logger, target_is_dir=False)
            meta['steps'].append({'key': key, 'run': rr, 'title': nm,
                                  'output': [os.path.basename(f)]})

        # ---- 宿主去除（库未就绪则跳过）----
        if host_ok:
            lo, hi = _band(0.40)
            p4 = {'r1': fq1, 'classify_r1': fa1, 'confidence': conf}
            if fq2:
                p4['r2'] = fq2
            if fa2:
                p4['classify_r2'] = fa2
            if ctx.p.get('db'):
                p4['db'] = ctx.p['db']
            run4, res4, _ = _chain_run_step(
                chain_dir, root, 'host', '宿主去除（kunpeng 分类）',
                _tool_job_hostremoval, p4, ctx.threads, ctx.db_virus,
                logger, prog, lo, hi, cancel=cancel,
                optional=('r2', 'classify_r1', 'classify_r2', 'db'))
            hr = os.path.join(root, run4, '01_host_removal')
            if os.path.isdir(hr):
                for fn in sorted(os.listdir(hr)):
                    if fn.startswith('kept_') and fn.endswith('.fastq.gz'):
                        _chain_link(os.path.join(chain_dir, fn),
                                    os.path.join(hr, fn), logger=logger,
                                    target_is_dir=False)
                        kept.append(os.path.join(hr, fn))
            meta['steps'].append({
                'key': 'host', 'run': run4,
                'title': '宿主去除（kunpeng 分类）',
                'summary': {k: v for k, v in (res4 or {}).items()
                            if not isinstance(v, (list, dict))}})
        else:
            meta['steps'].append({
                'key': 'host', 'run': None, 'title': '宿主去除（kunpeng 分类）',
                'skipped': '宿主库未就绪（「数据库构建」页构建后可单独运行）'})

        prog('chain', 1.0, '一键前处理完成')
        logger.close()
        meta['done'] = True
        _chain_write_manifest(chain_dir, meta)
        out = {'chain': chain_name, 'run': chain_name,
               'n_steps': len(meta['steps']), 'steps': meta['steps'],
               'kept_fastq': kept, 'fasta': fasta,
               'chain_json': os.path.join(chain_dir, '_chain.json')}
        if not host_ok:
            out['note'] = '宿主库未就绪，宿主去除已跳过'
        if not fastp_ok:
            out['note'] = (out.get('note', '') +
                           ('；' if out.get('note') else '') +
                           'fastp 未安装，质控已跳过').strip()
        return out
    return job


def _tool_job_annochain(ctx):
    """病毒注释分析 · 一键：注释三连 ① ORF 预测 → ② 功能注释 → ③ 基因组图谱。

    输入 FASTA（核酸 contigs / 基因组）。三步各建 run_dir：
      ① 步骤键固定为 'orf'（run 名 orf_*，布局 03_assembly + 04_orf）——
        ② 复用 _tool_job_orfa 的 run 模式直接对该运行注释（04b_orf_annot
        落在 ① 的运行目录，与样品管道 / 独立注释同一套产物口径），
        ③ 从中取 GFF3 给 gbdraw 画带注释的图谱。
    提交衔接：① 的运行自带 03_assembly + 04_orf + 04b 布局，跑完即出现在
    「数据提交」页「关联注释」下拉里；项目关联后 FASTA 推断 / featuretable
    / 提交包全部走既有机制，无需额外对接。
    """
    root = _tool_runs_root()
    fasta = ctx.req('fasta', '输入 FASTA（核酸 contigs / 基因组）')
    min_aa = max(1, min(int(ctx.p.get('min_aa') or 100), 5000))

    def job(log, prog, cancel):
        logger = TaskLogger(callback=log)
        chain_dir = ctx.run_dir
        chain_name = os.path.basename(chain_dir)
        meta = {'chain': 'annochain', 'run': chain_name,
                'label': '病毒注释分析 注释三连',
                'created': time.strftime('%Y-%m-%d %H:%M:%S'), 'steps': []}
        cur = 0.0

        def _band(w):
            nonlocal cur
            lo, hi = cur, cur + w
            cur = hi
            return lo, hi

        # ---- ① ORF 预测（键名固定 'orf'：run 名 orf_*，② 依赖该命名）----
        lo, hi = _band(0.30)
        run1, res1, _ = _chain_run_step(
            chain_dir, root, 'orf', '① ORF 预测', _tool_job_orf,
            {'fasta': fasta, 'min_aa': min_aa, 'annotate': False},
            ctx.threads, ctx.db_virus, logger, prog, lo, hi, cancel=cancel)
        n_ctg = len((res1 or {}).get('viral_contigs') or
                    (res1 or {}).get('contigs') or []) or None
        meta['steps'].append({
            'key': 'orf', 'run': run1, 'title': '① ORF 预测',
            'summary': {k: v for k, v in (res1 or {}).items()
                        if not isinstance(v, (list, dict))}})
        logger.log(f'① ORF 预测完成: run={run1}'
                   + (f'（{n_ctg} 条 contigs）' if n_ctg else ''))

        # ---- ② ORF 功能注释（复用 orfa 的 run 模式，产物写回 ① 运行目录）----
        lo, hi = _band(0.45)
        run2, res2, _ = _chain_run_step(
            chain_dir, root, 'orfa', '② ORF 功能注释', _tool_job_orfa,
            {'run': run1},
            ctx.threads, ctx.db_virus, logger, prog, lo, hi, cancel=cancel)
        ann_tsv = os.path.join(root, run1, '04b_orf_annot',
                               'orf_annotation.tsv')
        ann_linked = None
        if os.path.isfile(ann_tsv):
            ann_linked = os.path.basename(ann_tsv)
            _chain_link(os.path.join(chain_dir, ann_linked), ann_tsv,
                        logger=logger, target_is_dir=False)
        meta['steps'].append({
            'key': 'orfa', 'run': run2, 'title': '② ORF 功能注释',
            'target_run': run1,
            'output': [ann_linked] if ann_linked else [],
            'summary': {k: v for k, v in (res2 or {}).items()
                        if not isinstance(v, (list, dict))}})

        # ---- ③ 基因组图谱（FASTA + ①/② 的 GFF3）----
        gff = None
        d1 = os.path.join(root, run1)
        for pat in ('04b_orf_annot/*.gff3', '04b_orf_annot/*.gff',
                    '04_orf/*.gff'):
            hits = sorted(Path(d1).glob(pat))
            if hits:
                gff = str(hits[0])
                break
        lo, hi = _band(0.25)
        p3 = {'fasta': fasta, 'max_plots': 12, 'mode': 'both'}
        if gff:
            p3['ann'] = gff
            logger.log(f'③ 图谱注释来源: {os.path.basename(gff)}')
        run3, res3, _ = _chain_run_step(
            chain_dir, root, 'genoplot', '③ 基因组图谱', _tool_job_genoplot,
            p3, ctx.threads, ctx.db_virus, logger, prog, lo, hi,
            cancel=cancel, optional=('ann',))
        plots = []
        for pth in (res3 or {}).get('plots') or []:
            pth = str(pth)
            if os.path.isfile(pth):
                _chain_link(os.path.join(chain_dir, os.path.basename(pth)),
                            pth, logger=logger, target_is_dir=False)
                plots.append(os.path.basename(pth))
        meta['steps'].append({
            'key': 'genoplot', 'run': run3, 'title': '③ 基因组图谱',
            'output': plots,
            'summary': {k: v for k, v in (res3 or {}).items()
                        if not isinstance(v, (list, dict))}})

        prog('chain', 1.0, '注释三连完成')
        logger.close()
        meta['done'] = True
        _chain_write_manifest(chain_dir, meta)
        return {'chain': chain_name, 'run': chain_name,
                'n_steps': len(meta['steps']), 'steps': meta['steps'],
                'orf_run': run1,
                'annotation_tsv': ann_tsv if os.path.isfile(ann_tsv) else None,
                'plots': plots,
                'submit_hint': f'「数据提交」页建/打开项目后，'
                               f'「关联注释」选择 {run1} 即可对接 featuretable',
                'chain_json': os.path.join(chain_dir, '_chain.json')}
    return job


def _tool_job_chainoverview(ctx):
    """占位：链式汇总目录的读取走 /api/tool/chain_result，不注册为工具。"""
    raise NotImplementedError
