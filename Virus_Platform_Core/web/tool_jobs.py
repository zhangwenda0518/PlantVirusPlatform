# -*- coding: utf-8 -*-
"""工具箱：22 个独立分析任务工厂（自 app.py 拆出）。

每个 _tool_job_<name>(ctx) 返回 job(log, prog, cancel) 可调用；
由 Virus_Platform_Core/web/tools_api.py 的 /api/tool/run 取用。
新增工具三步：① 写 _tool_job_<name> ② 注册 TOOL_REGISTRY
（在 tools_api.py）③ webapp/templates/tools.html 加卡片。
"""
import json
import os
import re
import shutil
import subprocess
import threading
import time

from flask import abort

from Virus_Platform_Core.config import DIRS, PLATFORM_ROOT
from Virus_Platform_Core.utils import (DEFAULT_CMD_TIMEOUT_SEC, TaskLogger,
                                       check_path, safe_open,
                                       task_register_proc)
from Virus_Platform_Core.web.state import cfg, tool_runs_root as _tool_runs_root


# 共识模块 reads 未手动指定时的兜底层（contig 分类运行目录下）：
# 去宿主后全量，变异检测无偏。不再提供「病毒筛选后 / 质控后原始」等
# 选项——前者低估变异，后者带宿主污染。
CONSENSUS_FALLBACK_LAYER = '01_host_removal'

# kvsuite 卡的默认阶段：识别/定量（identify）+ 过滤（filter）两段。
# 共识/绘图/变异属于下游模块（t-consensus / t-variant），不在「识别与定量」范围内。
# 字面量与引擎侧 Virus_Platform_Core/known_virus_suite/known_virus_suite.py 的
# DEFAULT_STAGE 必须一致；这里不 import 引擎模块——冻结分发版不打包 kv_filter /
# kv_identify（见 main.py 的 import 白名单说明）。
KV_DEFAULT_STAGE = 'identify+filter'


def _opt_num_param(ctx, key):
    """可选数值门槛参数（⑥b 层1 min_pident/min_qcov 用）。

    空/0/非法 → None（None = 不过滤，保持历史结果）；正确值按 float 返回。
    """
    v = ctx.p.get(key)
    if v in (None, '', 0, '0'):
        return None
    try:
        f = float(v)
    except (TypeError, ValueError):
        return None
    return f if f > 0 else None


def _job_req(ctx, key, what):
    """job 体内取必填输入文件（与 tools_api 的 ctx.req 同语义，但不 abort）。

    任务线程里没有 Flask 请求上下文：abort(400) 只会被记成
    「BadRequest: 400 Bad Request」，用户看不到「请选择xxx」的文案。
    这里改抛带原文案的 ValueError（TaskManager 会把 message 原样呈现）。
    相对路径按平台根解析；绝对路径允许电脑任意位置（只读）。
    """
    v = (ctx.p.get(key) or '').strip()
    if not v:
        raise ValueError(f'请选择{what}')
    return check_path(v if os.path.isabs(v)
                      else os.path.join(PLATFORM_ROOT, v), must_exist=True)


def _tool_job_convert(ctx):
    """格式转换：.sra→FASTQ/FASTA（sracha）、FASTQ→FASTA（seqkit）。

    双端 FASTQ（input + input2）：R1/R2 **分别**转换、分别输出
    （R1→R1 产物、R2→R2 产物，不合并）；.sra 自含双端，不接受 input2。"""
    inp = ctx.req('input', '输入文件')
    inp2 = ctx.opt('input2')
    target = ctx.p.get('target') or 'fastq'
    if target not in ('fastq', 'fasta'):
        abort(400, '无效的目标格式')
    ilower = str(inp).lower()
    if ilower.endswith('.sra'):
        mode = 'sra'
        if inp2:
            abort(400, '.sra 文件自含双端信息，无需另选 R2')
    elif ilower.endswith(('.fastq', '.fq', '.fastq.gz', '.fq.gz')):
        mode = 'fastq2fasta' if target == 'fasta' else None
        if inp2 and not str(inp2).lower().endswith(
                ('.fastq', '.fq', '.fastq.gz', '.fq.gz')):
            abort(400, 'R2 文件格式与 R1 不一致（需 .fastq[.gz] / .fq[.gz]）')
    elif ilower.endswith(('.fa', '.fasta', '.fa.gz', '.fasta.gz')):
        abort(400, '输入已是 FASTA，无需转换')
    else:
        abort(400, '无法识别的输入格式（支持 .sra / .fastq[.gz]）')

    def job(log, prog, cancel):
        logger = TaskLogger(callback=log)
        out_files = []
        if mode == 'sra':
            from Virus_Platform_Core.public_data import sra_convert_engine
            engine, exe = sra_convert_engine()
            if engine != 'sracha':
                raise RuntimeError('sracha.exe 不可用，无法转换 .sra')
            threads = ctx.threads or max(2, min((os.cpu_count() or 4) // 2, 8))
            prog('sra', 0.2, f'sracha 解码 .sra → {target.upper()}')
            cmd = [exe, 'fastq', inp, '-O', ctx.run_dir,
                   '-t', str(threads), '-f', '-q']
            if target == 'fasta':
                cmd.append('--fasta')
            from Virus_Platform_Core.public_data import _sra_fastq_outputs
            from Virus_Platform_Core.utils import run_cmd
            # 走 run_cmd：stdout/stderr 实时进任务日志（原来 _sp.run 丢弃
            # 输出且 stderr 只留在内存里，转换全程无日志可看），并接入
            # 取消检查与超时看门狗。
            run_cmd(cmd, logger=logger, timeout=8 * 3600)
            base = os.path.basename(inp)[:-4]
            # ⚠ 单端 run：sracha 写 `<acc>.fastq.gz`（**无下划线**）；双端才写
            #   `<acc>_1/_2.fastq.gz`。只匹配 `base + '_'` 会把单端产物判成
            #   「无输出」，把已经落盘的 GB 级 FASTQ 丢掉 —— 与下载链
            #   （public_data._convert_sras）同源的历史 bug。这里复用同一份
            #   判定函数，杜绝两处再走偏。
            out_files = [os.path.join(ctx.run_dir, f_) for f_ in
                         _sra_fastq_outputs(
                             ctx.run_dir, base,
                             exts=('.fastq.gz', '.fq.gz',
                                   '.fa.gz', '.fasta.gz'))]
            if not out_files:
                raise RuntimeError('sracha 无输出')
        else:
            seqkit = cfg.tool('seqkit')
            from Virus_Platform_Core.utils import run_cmd
            # 双端：R1/R2 逐条独立转换（文件名自带 _R1/_R2 后缀，产物
            # 自然分开命名）；单端循环一次，行为与原来完全一致。
            inputs = [inp, inp2] if inp2 else [inp]
            for i, src in enumerate(inputs):
                tag = f'R{i + 1} ' if len(inputs) > 1 else ''
                dst = os.path.join(ctx.run_dir,
                                   os.path.basename(src).rsplit('.', 2)[0]
                                   + '.fa.gz')
                prog('fq2fa', 0.3 if i == 0 else 0.7,
                     f'seqkit fq2fa 转换中（{tag or "单端"}）')
                run_cmd([seqkit, 'fq2fa', '-w', '0',
                         '-j', str(ctx.threads or cfg.threads), src,
                         '-o', dst], logger=logger)
                out_files.append(dst)
        prog('done', 1.0, f'完成：{len(out_files)} 个文件')
        logger.close()
        return {'n_files': len(out_files),
                'files': [os.path.relpath(f, ctx.run_dir).replace(os.sep, '/')
                          for f in out_files]}
    return job


def _tool_job_fastp(ctx):
    """① 质控预处理（fastp，单/双端）。"""
    r1 = ctx.req('r1', 'R1 FASTQ')
    r2 = ctx.opt('r2')

    def job(log, prog, cancel):
        logger = TaskLogger(callback=log)
        from Virus_Platform_Core.preprocess import run_fastp
        prog('fastp', 0.3, 'fastp 质控中')
        res = run_fastp(ctx.run_dir, r1, r2, threads=ctx.threads, logger=logger,
                        force=True, dedup=bool(ctx.p.get('dedup')))
        if not res:
            raise RuntimeError('未检测到 fastp.exe，无法运行质控')
        prog('fastp', 1.0, '完成')
        logger.close()
        return res
    return job


def _tool_job_hostremoval(ctx):
    """宿主去除与序列提取（kunpeng 宿主库分类，独立模块，不依赖样品管道）。

    C 行 = 宿主 read 对，剔除后保留非宿主 reads（kept_R1/R2.fastq.gz）。
    classify_r1/classify_r2（可选）：预转换 FASTA，分类步骤直接用（省去
    质量行解压，与样品管道 ⓪b→① 同语义）；kept reads 仍从 FASTQ 提取。
    """
    r1 = ctx.req('r1', 'R1 FASTQ')
    r2 = ctx.opt('r2')
    conf = float(ctx.p.get('confidence') or 0)
    db_host = ctx.opt('db') or cfg.databases['host']
    cls_r1 = ctx.opt('classify_r1')
    cls_r2 = ctx.opt('classify_r2')

    from Virus_Platform_Core.kunpeng import db_ready
    if not db_ready(db_host):
        abort(400, '宿主库未就绪，请先到「数据库构建」页构建宿主库')

    def job(log, prog, cancel):
        logger = TaskLogger(callback=log)
        from Virus_Platform_Core.host_removal import remove_host
        prog('classify', 0.05, 'kunpeng 宿主库分类中')
        res = remove_host(ctx.run_dir, r1, r2, db_host, threads=ctx.threads,
                          confidence=conf, logger=logger, force=True,
                          classify_r1=cls_r1, classify_r2=cls_r2,
                          progress=lambda pct, msg: prog(
                              'classify' if pct < 0.9 else 'filter', pct, msg))
        prog('done', 1.0, '完成')
        logger.close()
        return res
    return job


def _tool_job_orf(ctx):
    """ORF 预测（pyrodigal / pyrodigal_rv），可选 ⑥b 功能注释。"""
    fasta = ctx.req('fasta', '输入 FASTA（核酸 contigs / 基因组）')
    min_aa = max(1, min(int(ctx.p.get('min_aa') or 100), 5000))
    annotate = bool(ctx.p.get('annotate'))
    orf_tool = (ctx.p.get('orf_tool') or '').strip()
    orf_engine = (ctx.p.get('engine') or '').strip()
    orf_db = (ctx.p.get('db') or '').strip()
    # 功能注释 CDS 模型：pyrodigal_rv（默认）/ pyrodigal
    orf_model = (ctx.p.get('model') or '').strip()
    if orf_model not in ('pyrodigal_rv', 'pyrodigal'):
        orf_model = ''
    # 注解层选择：prot（序列同源）/ pfam（HMM）；空=全部
    orf_libs = ctx.p.get('libs')
    # 层1 命中门槛（百分数；空/0 = 不过滤）
    min_pident = _opt_num_param(ctx, 'min_pident')
    min_qcov = _opt_num_param(ctx, 'min_qcov')

    def job(log, prog, cancel):
        import json as _json
        logger = TaskLogger(callback=log)
        a_dir = check_path(os.path.join(ctx.run_dir, '03_assembly'),
                           must_exist=False, in_platform=True)
        os.makedirs(a_dir, exist_ok=True)
        prog('prep', 0.03, '整理输入')
        from Virus_Platform_Core.utils import iter_fasta, write_fasta_record
        ids = []
        vfa = os.path.join(a_dir, 'viral_contigs.fasta')
        cfa = os.path.join(a_dir, 'contigs.filtered.fasta')
        with safe_open(vfa, 'wt') as w, safe_open(cfa, 'wt') as wc:
            for h, s in iter_fasta(fasta):
                cid = h.split()[0]
                ids.append(cid)
                write_fasta_record(w, cid, s)
                write_fasta_record(wc, cid, s)
        if not ids:
            raise RuntimeError('输入 FASTA 中没有序列')
        with safe_open(os.path.join(a_dir, 'summary.json'), 'wt') as f:
            _json.dump({'viral_contigs': ids, 'standalone': True}, f)
        from Virus_Platform_Core.orf import predict_orfs
        prog('orf', 0.08, 'ORF 基因预测')
        res = predict_orfs(ctx.run_dir, min_aa=min_aa, threads=ctx.threads,
                           logger=logger, force=True, tools=orf_tool or None,
                           progress=lambda p, m: prog(
                               'orf', 0.08 + p * 0.6, m))
        res = dict(res)
        if annotate:
            from Virus_Platform_Core.orf_annot import run_orf_annotation
            prog('orfa', 0.72, 'ORF 功能注释')
            res['orfa'] = run_orf_annotation(
                ctx.run_dir, threads=ctx.threads, logger=logger, force=True,
                engine=orf_engine or None, db=orf_db or None,
                model=orf_model or None, libs=orf_libs,
                min_pident=min_pident, min_qcov=min_qcov,
                progress=lambda p, m: prog('orfa', 0.72 + p * 0.26, m))
        prog('done', 1.0, '完成')
        logger.close()
        return res
    return job


def _tool_job_orfa(ctx):
    """功能注释（独立模块）：对已有 orf_ 运行注释，或 FASTA 预测+注释一步完成。"""
    run = (ctx.p.get('run') or '').strip()
    orf_engine = (ctx.p.get('engine') or '').strip()
    orf_db = (ctx.p.get('db') or '').strip()
    # 功能注释 CDS 模型：pyrodigal_rv（默认）/ pyrodigal
    orf_model = (ctx.p.get('model') or '').strip()
    if orf_model not in ('pyrodigal_rv', 'pyrodigal'):
        orf_model = ''
    # 注释层选择：prot（层1 序列同源）/ pfam（层2 HMM）
    orf_libs = ctx.p.get('libs')
    # 层1 命中门槛（百分数；空/0 = 不过滤，保持历史结果）
    min_pident = _opt_num_param(ctx, 'min_pident')
    min_qcov = _opt_num_param(ctx, 'min_qcov')
    if run:
        if not re.fullmatch(r'[A-Za-z0-9_]+', run) or not run.startswith('orf_'):
            abort(400, f'无效的 ORF 运行名: {run}')
        target = check_path(os.path.join(_tool_runs_root(), run),
                            must_exist=True, in_platform=True)

        def job(log, prog, cancel):
            logger = TaskLogger(callback=log)
            from Virus_Platform_Core.orf_annot import run_orf_annotation
            prog('orfa', 0.15, f'对运行 {run} 做 ORF 功能注释')
            res = run_orf_annotation(target, threads=ctx.threads, logger=logger,
                                     force=True, engine=orf_engine or None,
                                     db=orf_db or None,
                                     model=orf_model or None, libs=orf_libs,
                                     min_pident=min_pident, min_qcov=min_qcov,
                                     progress=lambda p, m: prog(
                                         'orfa', 0.15 + p * 0.8, m))
            res = dict(res)
            res['run'] = run
            prog('done', 1.0, '完成')
            logger.close()
            return res
        return job
    # 无 run → FASTA 输入：预测 + 注释一步完成
    ctx.p = dict(ctx.p)
    ctx.p['annotate'] = True
    if not (ctx.p.get('fasta') or '').strip():
        abort(400, '请选择已有 ORF 运行或输入 FASTA')
    return _tool_job_orf(ctx)


def _tool_job_genoplot(ctx):
    """基因组图谱（gbdraw 首选，缺则 DFV 顶上）。

    输入 FASTA（可选配 GFF3 注释）或 GenBank（.gb/.gbk，自带注释）。
    """
    fasta = ctx.opt('fasta')
    ann = ctx.opt('ann')
    if not fasta and not ann:
        abort(400, '请选择 FASTA 或 GenBank 输入')
    if ann and str(ann).lower().endswith(('.gb', '.gbk', '.gbff', '.genbank')):
        fasta = None          # GenBank 自带序列与注释，FASTA 忽略
    max_plots = max(1, min(int(ctx.p.get('max_plots') or 12), 200))
    # 绘图定制参数（透传 gbdraw CLI）：仅收集有值/True 的项
    gb_opts = {}
    for _k, _v in (ctx.p.get('gb_opts') or {}).items():
        if _v is None or _v == '' or _v is False:
            continue
        gb_opts[str(_k)] = _v
    plot_mode = (ctx.p.get('mode') or 'both')   # circular / linear / both

    def job(log, prog, cancel):
        logger = TaskLogger(callback=log)
        from Virus_Platform_Core.gbdraw_plot import run_genome_plots, gbdraw_available
        if not gbdraw_available():
            raise RuntimeError('未检测到 gbdraw（管道不支持 DFV）')
        tag = 'gbdraw'
        prog(tag, 0.05, f'{tag} 出图中')
        res = run_genome_plots(ctx.run_dir, logger=logger, force=True,
                               max_plots=max_plots, fasta_in=fasta, ann_in=ann,
                               opts=gb_opts, mode=plot_mode,
                               progress=lambda p, m: prog(tag, 0.05 + p * 0.9, m))
        if not res.get('plots'):
            raise RuntimeError('未产出任何基因组图（检查输入文件与绘图引擎）')
        # run_genome_plots 的返回契约是**相对路径**（相对 09_genome_plots/，
        # 供 summary.json / 报告端按 sample_dir 拼接）；而任务结果的消费者
        # （集成测试 / API 调用方）拿到的是独立于 run 目录的 dict，必须在这里
        # 转成绝对路径——否则 res['plots'] 全是相对名，os.path.isfile 全 False。
        # 磁盘上的 summary.json 不受影响，仍是相对口径。
        _gp_dir = os.path.join(ctx.run_dir, '09_genome_plots')
        res['plots'] = [p if os.path.isabs(p) else os.path.join(_gp_dir, p)
                        for p in res['plots']]
        res['run'] = os.path.basename(ctx.run_dir)   # 供前端 toolrun 内联展示 SVG
        prog('done', 1.0, f"完成：{len(res.get('plots', []))} 张图")
        logger.close()
        return res
    return job


def _tool_job_primer(ctx):
    """引物设计（primer3）。plain=基因组/contigs 全长分窗；
    conserved=多序列比对 FASTA 保守区（输入需已比对，如 MAFFT aln.fasta）。"""
    fasta = ctx.req('fasta', '输入 FASTA')
    mode = ctx.p.get('mode') or 'plain'
    if mode not in ('conserved', 'plain'):
        abort(400, '无效的引物设计模式')
    num_return = max(1, min(int(ctx.p.get('num_return') or 3), 20))
    specificity = bool(ctx.p.get('specificity'))

    def job(log, prog, cancel):
        import json as _json
        logger = TaskLogger(callback=log)
        if mode == 'conserved':
            p_dir = check_path(os.path.join(ctx.run_dir, '05_phylo'),
                               must_exist=False, in_platform=True)
            os.makedirs(os.path.join(p_dir, 'G1'), exist_ok=True)
            shutil.copyfile(fasta, os.path.join(p_dir, 'G1', 'aln.fasta'))
            # synthetic 标记：这只是把用户自备比对"伪装"成 ⑦ 的产物以便
            # 复用 primer 的 conserved 通路，不是真系统发育结果——消费者
            # （viz/refs 等按真摘要解析的地方）据此可识别残缺结构
            with safe_open(os.path.join(p_dir, 'summary.json'), 'wt') as f:
                _json.dump({'groups': [{'group': 'G1', 'dir': 'G1'}],
                            'synthetic': True,
                            'note': '工具卡片自建：仅含用户自备比对，'
                                    '非 ⑦ 系统发育阶段产物'}, f)
        else:
            a_dir = check_path(os.path.join(ctx.run_dir, '03_assembly'),
                               must_exist=False, in_platform=True)
            os.makedirs(a_dir, exist_ok=True)
            shutil.copyfile(fasta, os.path.join(a_dir, 'viral_contigs.fasta'))
        from Virus_Platform_Core.primer import design_primers
        prog('primer', 0.1, f'primer3 引物设计（{mode}）')
        res = design_primers(ctx.run_dir, mode=mode, num_return=num_return,
                             logger=logger, force=True,
                             do_specificity=specificity)
        prog('done', 1.0, f"引物 {res.get('n_primers', 0)} 对")
        logger.close()
        return res
    return job


def _tool_job_identify(ctx):
    """② 病毒鉴定与提取（fastq 双端/单端 或 fasta contigs）。"""
    inp = ctx.req('input', '输入文件')
    itype = ctx.p.get('input_type') or 'pe'
    if itype not in ('pe', 'single', 'fasta'):
        abort(400, '无效的输入类型')
    conf = float(ctx.p.get('confidence') or 0)

    def job(log, prog, cancel):
        logger = TaskLogger(callback=log)
        from Virus_Platform_Core.kunpeng import classify, parse_classify_output
        inputs = [inp]
        if itype == 'pe':
            # 工厂期没查 R2：这里补的必须是 raise（任务线程里 abort 丢文案）
            inputs.append(_job_req(ctx, 'input2', 'R2 FASTQ'))
        prog('classify', 0.1, 'kunpeng 病毒库分类中')
        res = classify(ctx.db_virus, inputs, os.path.join(ctx.run_dir, 'classify'),
                       paired=(itype == 'pe'), threads=ctx.threads,
                       confidence=conf, logger=logger,
                       progress=lambda pp, mm: prog(
                           'classify', 0.1 + pp * 0.75, mm))
        ids = {}
        if res['kraken']:
            for flag, rid, taxid, _l, mapping in parse_classify_output(res['kraken']):
                if flag != 'C':
                    continue
                # kraken2 口径分值：支持该 taxid 的片段数 / 映射列总片段数
                sup = tot = 0
                for seg in (mapping or '').split():
                    t, _, c = seg.rpartition(':')
                    try:
                        cnt = int(c)
                    except ValueError:
                        continue
                    tot += cnt
                    if t == str(taxid):
                        sup += cnt
                ids[rid] = (taxid, round(sup / tot, 3) if tot else 0)
        prog('extract', 0.9, '提取病毒候选序列')
        n_ext = {}
        if ids:
            for i, s_ in enumerate(inputs, 1):
                tag = '' if len(inputs) == 1 else f'_{i}'
                dst = os.path.join(ctx.run_dir, f'viral_sequences{tag}.fasta')
                n_ext[os.path.basename(dst)] = _extract_records(s_, set(ids), dst)
                logger.log(f'提取病毒序列 {os.path.basename(dst)}: '
                           f'{n_ext[os.path.basename(dst)]} 条')
            with safe_open(os.path.join(ctx.run_dir, 'viral_ids.tsv'), 'wt') as f:
                f.write('seq_id\ttaxid\tscore\n')
                for rid, (tx, sc) in ids.items():
                    f.write(f'{rid}\t{tx}\t{sc}\n')
        prog('extract', 1.0, '完成')
        logger.close()
        # n_extracted 必须是标量：结果预览的统计条只展示标量字段
        return {'n_classified': len(ids), 'n_extracted': sum(n_ext.values()),
                'kreport': res['kreport']}
    return job


def _tool_job_assemble(ctx):
    """③ 病毒组装（SPAdes，可选模式）。"""
    r1 = ctx.req('r1', 'R1 FASTQ')
    mode = ctx.p.get('mode') or 'rnaviral'
    # srna（小RNA）模式按单端 -s 组装，R2 可选；其余模式仍必填
    r2 = ctx.opt('r2') if mode == 'srna' else ctx.req('r2', 'R2 FASTQ')
    mem = int(ctx.p.get('memory') or 64)
    min_len = int(ctx.p.get('min_len') or 200)

    def job(log, prog, cancel):
        logger = TaskLogger(callback=log)
        from Virus_Platform_Core.assembly import run_spades, filter_contigs
        out = os.path.join(ctx.run_dir, 'assembly')
        prog('spades', 0.05, 'SPAdes 组装中（耗时主要步骤）')
        spades_out = run_spades(r1, r2, out, mode=mode, threads=ctx.threads,
                                memory_gb=mem, logger=logger,
                                progress=lambda pp, mm: prog(
                                    'spades', 0.05 + pp * 0.85, mm))
        # 组装产物可能是 contigs.fasta（常规/rnaviral/metaviral-meta）或
        # transcripts.fasta（低覆盖自动降级 rna 模式），以 run_spades 实际
        # 返回的路径为准，不能硬编码 contigs.fasta（否则降级时找不到文件）。
        prog('filter', 0.92, 'contig 长度过滤')
        filtered, n_c, total_bp = filter_contigs(
            spades_out,
            os.path.join(ctx.run_dir, 'contigs.filtered.fasta'),
            min_len=min_len, logger=logger)
        prog('filter', 1.0, '完成')
        logger.close()
        return {'n_contigs': n_c, 'total_bp': total_bp,
                'contigs': filtered}
    return job


def _tool_job_contigs(ctx):
    """④ contig 病毒分类与提取（输入 contigs fasta）。"""
    contigs = ctx.req('contigs', 'contigs FASTA')
    min_len = int(ctx.p.get('min_len') or 200)
    conf = float(ctx.p.get('confidence') or 0)

    def job(log, prog, cancel):
        logger = TaskLogger(callback=log)
        from Virus_Platform_Core.assembly import filter_contigs
        from Virus_Platform_Core.kunpeng import classify, parse_classify_output
        from Virus_Platform_Core.contig_annot import classify_rows, genus_avg_map, RANKS

        prog('genus_lens', 0.03, '统计属平均基因组长度（首跑需建缓存）')
        genus_map = genus_avg_map(logger=logger)

        prog('filter', 0.05, 'contig 长度过滤')
        filtered, n_c, _bp = filter_contigs(
            contigs, os.path.join(ctx.run_dir, 'contigs.filtered.fasta'),
            min_len=min_len, logger=logger)
        if n_c == 0:
            raise RuntimeError('过滤后无 contigs（检查最小长度设置）')

        prog('classify', 0.15, 'kunpeng 病毒库分类')
        res = classify(ctx.db_virus, [filtered],
                       os.path.join(ctx.run_dir, 'classify'),
                       paired=False, threads=ctx.threads, confidence=conf,
                       logger=logger,
                       progress=lambda pp, mm: prog(
                           'classify', 0.15 + pp * 0.45, mm))
        ids = {}
        if res['kraken']:
            for flag, _rid, _tx, _l, _pa in parse_classify_output(res['kraken']):
                if flag == 'C':
                    ids[_rid] = _tx

        # metabuli 风格分类表：8 级谱系 + 属平均长度 + 近完整判定
        prog('annot', 0.8, '谱系注释与属长比整理')
        rows = []
        if res['kraken']:
            rows = classify_rows(res['kraken'], genus_map)
        tsv_path = os.path.join(ctx.run_dir, 'virus_classification.tsv')
        header = (['contig', 'taxid', 'taxon'] + RANKS
                  + ['length', 'genus_avg_len', 'ratio', 'near_complete',
                     'score', 'kmer_support', 'kmer_total'])
        with safe_open(tsv_path, 'wt') as f:
            f.write('\t'.join(header) + '\n')
            for r in rows:
                f.write('\t'.join(str(r.get(k, '')) for k in header) + '\n')

        prog('extract', 0.88, '提取病毒 contigs（带谱系 header）')
        n_viral = 0
        viral_fa = None
        if ids:
            viral_fa = os.path.join(ctx.run_dir, 'viral_contigs.fasta')
            seqs = _load_fasta_seqs(filtered)
            with safe_open(viral_fa, 'wt') as f:
                for rid, taxid in ids.items():
                    seq = seqs.get(rid)
                    if seq is None:
                        continue
                    n_viral += 1
                    row = next((r for r in rows if r['contig'] == rid), None)
                    lineage = ';'.join(row[r] for r in RANKS
                                       if row and row.get(r))
                    taxon = (row or {}).get('taxon', '')
                    f.write(f'>{rid} taxid={taxid} taxon='
                            f'{taxon.replace(" ", "_")} '
                            f'lineage={lineage}\n')
                    for i in range(0, len(seq), 70):
                        f.write(seq[i:i + 70] + '\n')

        prog('extract', 1.0, '完成')
        logger.close()
        return {'n_contigs': n_c, 'n_viral': n_viral,
                'kreport': res['kreport'],
                'classification': tsv_path, 'viral_fasta': viral_fa}
    return job


def _tool_job_structcmp(ctx):
    """结构比较：多条序列 MAFFT 全长比对 → 两两 identity 矩阵（SDT 口径）。"""
    seqs_fa = ctx.req('seqs', '序列 FASTA')
    max_n = max(3, min(int(ctx.p.get('max_n') or 30), 200))

    def _short(h, idx):
        name = re.split(r'[\s|]', (h or '').strip())[0][:40]
        name = re.sub(r'[^A-Za-z0-9_\-.]', '_', name)
        return name or f'seq{idx}'

    def job(log, prog, cancel):
        logger = TaskLogger(callback=log)
        from Virus_Platform_Core.phylo import _run_mafft
        from Virus_Platform_Core.utils import iter_fasta

        prog('read', 0.05, '读取与筛选序列')
        recs, seen = [], {}
        for h, s in iter_fasta(seqs_fa):
            name = _short(h, len(recs) + 1)
            if name in seen:
                seen[name] += 1
                name = f'{name}_{seen[name]}'
            else:
                seen[name] = 0
            recs.append((name, s.upper()))
        if len(recs) < 2:
            raise RuntimeError('FASTA 中少于 2 条序列，无法做结构比较')
        recs = recs[:max_n]

        capped_fa = os.path.join(ctx.run_dir, 'input.fasta')
        with safe_open(capped_fa, 'wt') as f:
            for name, s in recs:
                f.write(f'>{name}\n')
                for i in range(0, len(s), 70):
                    f.write(s[i:i + 70] + '\n')

        prog('aln', 0.15, 'MAFFT 全长比对')
        aln = _run_mafft(capped_fa, os.path.join(ctx.run_dir, 'aln.fasta'),
                         threads=ctx.threads, logger=logger)

        prog('matrix', 0.75, '计算两两 identity 矩阵')
        names, cols = [], []
        for h, s in iter_fasta(aln):
            names.append(_short(h, len(names) + 1))
            cols.append(s.upper())
        n = len(names)
        # 同一性矩阵要求所有序列等长（MAFFT 全长比对的结果）。zip() 会按较短
        # 序列静默截断：MAFFT 未正常完成时算出的 identity 会偏大且无任何提示，
        # 这里先校验长度，不等直接报错（比对列数也才有意义）。
        if len({len(c) for c in cols}) > 1:
            _bad = next(i for i, c in enumerate(cols)
                        if len(c) != len(cols[0]))
            raise RuntimeError(
                f'比对结果长度不一致（{names[0]}: {len(cols[0])} vs '
                f'{names[_bad]}: {len(cols[_bad])}），无法计算 identity 矩阵；'
                '请检查 MAFFT 是否正常完成')
        matrix = [[100.0] * n for _ in range(n)]
        for i in range(n):
            for j in range(i + 1, n):
                same = comp = 0
                for a, b in zip(cols[i], cols[j]):
                    if a == '-' and b == '-':
                        continue
                    comp += 1
                    if a == b:
                        same += 1
                pid = round(same / comp * 100, 2) if comp else 0.0
                matrix[i][j] = matrix[j][i] = pid

        tsv = os.path.join(ctx.run_dir, 'identity_matrix.tsv')
        with safe_open(tsv, 'wt') as f:
            f.write('\t'.join([''] + names) + '\n')
            for i in range(n):
                f.write('\t'.join([names[i]] +
                                  [f'{matrix[i][j]:.2f}' for j in range(n)]) + '\n')
        data = {'names': names, 'matrix': matrix, 'n': n,
                'aln_cols': len(cols[0]) if cols else 0}
        js = os.path.join(ctx.run_dir, 'identity_matrix.json')
        with safe_open(js, 'wt') as f:
            json.dump(data, f, ensure_ascii=False)

        prog('done', 1.0, '完成')
        logger.close()
        return {'n_seqs': n, 'aln': aln, 'matrix_tsv': tsv,
                'matrix_json': js, 'aln_cols': data['aln_cols']}
    return job


def _tool_job_verify(ctx):
    """候选序列验证（对齐 02b/09b）：宿主筛选 → 长度分流 → 双路过滤。

    输入可选已有 contigs 运行（run）或独立 FASTA（fasta）。
    参数：host（默认 all）、methods（blastx/bfvd 组合）、combine（union/intersection）。
    （B 路证据源 2026-09-29 由 CDD 换为 BFVD：病毒专属库、免白名单。）
    """
    run_ref = (ctx.p.get('run') or '').strip()
    fasta = ctx.opt('fasta') if ctx.p.get('fasta') else None
    host = (ctx.p.get('host') or 'all').strip() or 'all'
    combine = (ctx.p.get('combine') or 'union').strip()
    methods = ctx.p.get('methods') or ['blastx', 'bfvd']
    if isinstance(methods, str):
        methods = [m.strip() for m in methods.split(',') if m.strip()]
    methods = [m for m in methods if m in ('blastx', 'bfvd')]

    if run_ref and not re.fullmatch(r'[A-Za-z0-9_\-]+', run_ref):
        abort(400, '无效的运行名')
    if combine not in ('union', 'intersection'):
        combine = 'union'

    def job(log, prog, cancel):
        logger = TaskLogger(callback=log)
        from Virus_Platform_Core.verify import verify

        # 输入解析：优先 run（其 viral_contigs.fasta），否则独立 fasta
        run_dir = ctx.run_dir
        if run_ref:
            src_run = check_path(os.path.join(_tool_runs_root(), run_ref),
                                 must_exist=True, in_platform=True)
            vfa = os.path.join(src_run, 'viral_contigs.fasta')
            if not os.path.isfile(vfa):
                raise RuntimeError('该运行无 viral_contigs.fasta（先跑 contig 分类）')
            # 宿主归属依赖源运行的分类表
            for _fn in ('virus_classification.tsv',):
                s = os.path.join(src_run, _fn)
                if os.path.isfile(s):
                    import shutil
                    shutil.copyfile(s, os.path.join(run_dir, _fn))
        else:
            if not fasta:
                # job 体在任务线程执行，abort(400) 只会变成
                # 「BadRequest: 400 Bad Request」（用户看不到提示文案）
                raise ValueError('请选择 contigs 运行或输入 FASTA')
            vfa = fasta

        summary = verify(run_dir, vfa, host=host, methods=methods,
                         combine=combine, threads=ctx.threads, logger=logger,
                         progress=lambda st, fr, msg: prog(st, fr, msg))
        logger.close()
        return summary
    return job


def _tool_job_consensus(ctx):
    """共识序列与变异分析（验证模块之后）：reads 回贴 → 共识序列 + 变异谱。

    参考三条来源：
      1. kvsuite 运行选参考（virus-fasta/ref_<acc>/）—— 已知病毒基因组
      2. contig 分类运行的 viral_contigs.fasta —— 未知/组装候选
      3. 独立 FASTA
    固定**重新比对**（minibwa map 双位置参数双端回贴），不复用 kvsuite 的
    BAM——那是按参考筛过的 reads 子集，有比对偏好性，会系统性低估变异。
    reads 默认来源（未手动指定时）：
      1. 所选 kvsuite 运行的输入序列（tool_runs/<run>/sample_sheet.tsv，
         即「已知病毒识别与定量」吃进去的 reads）
      2. contig 分类运行目录的 01_host_removal（去宿主后全量，变异无偏）
      3. 都没有 → 报错要求手动指定
    """
    run_ref = (ctx.p.get('run') or '').strip()
    fasta = ctx.opt('fasta') if ctx.p.get('fasta') else None
    kv_run = (ctx.p.get('kv_run') or '').strip()
    kv_refs = [x.strip() for x in (ctx.p.get('kv_refs') or '').split(',') if x.strip()]
    ambig = (ctx.p.get('ambig') or 'N').strip()[:1] or 'N'

    def _num(v, default, cast):
        try:
            return cast(v)
        except (TypeError, ValueError):
            return default

    min_qual = _num(ctx.p.get('min_qual'), 20, int)
    # 默认 4（深度>3）与 consensus.MIN_DEPTH 保持一致（1 会让单条 read 就能定共识/出 SNV）
    min_depth = _num(ctx.p.get('min_depth'), 4, int)
    min_freq = _num(ctx.p.get('min_freq'), 0.5, float)
    min_mapq = _num(ctx.p.get('min_mapq'), 10, int)
    min_minor_freq = _num(ctx.p.get('min_minor_freq'), 0.02, float)
    min_cov_pct = _num(ctx.p.get('min_cov_pct'), 10.0, float)
    min_qual = max(0, min(min_qual, 60))
    min_depth = max(1, min(min_depth, 100000))
    min_freq = min(1.0, max(0.0, min_freq))

    if run_ref and not re.fullmatch(r'[A-Za-z0-9_\-]+', run_ref):
        abort(400, '无效的运行名')
    if kv_run and not re.fullmatch(r'[A-Za-z0-9_\-]+', kv_run):
        abort(400, '无效的 kvsuite 运行名')
    if kv_refs and not kv_run:
        abort(400, '选了参考序列但未指定 kvsuite 运行')
    for _r in kv_refs:
        if not re.fullmatch(r'[A-Za-z0-9_.\-]+', _r):
            abort(400, '无效的 accession: %s' % _r)

    def job(log, prog, cancel):
        logger = TaskLogger(callback=log)
        from Virus_Platform_Core.consensus import (consensus_and_variants,
                                                   find_read_pairs,
                                                   pair_read_files)

        # ---- 参考 ----
        # 三条来源，优先级：kvsuite 选参考 > contig 运行 > 独立 FASTA
        #   kvsuite：<run>/kvsuite/virus-fasta/ref_<acc>/ref_<acc>.ref.fasta
        #            多选按 accession 合并成临时库（序列名已保证唯一）
        src_run = None
        vfa = None
        if kv_refs:
            kv_base = check_path(os.path.join(_tool_runs_root(), kv_run,
                                              'kvsuite'),
                                 must_exist=True, in_platform=True)
            fa_dir = os.path.join(kv_base, 'virus-fasta')
            parts, missing = [], []
            for acc in kv_refs:
                sub = os.path.join(fa_dir, 'ref_%s' % acc)
                cand = os.path.join(sub, 'ref_%s.ref.fasta' % acc)
                if os.path.isfile(cand) and os.path.getsize(cand) > 0:
                    parts.append((acc, cand))
                    continue
                # 容错：目录在但文件名不同
                found = None
                if os.path.isdir(sub):
                    for fn in sorted(os.listdir(sub)):
                        if fn.endswith(('.fasta', '.fa', '.fna')):
                            found = os.path.join(sub, fn)
                            break
                if found:
                    parts.append((acc, found))
                else:
                    missing.append(acc)
            if not parts:
                raise RuntimeError(
                    '所选参考在 %s 下均无 FASTA（先跑共识段生成 virus-fasta）'
                    % fa_dir)
            if missing:
                log('跳过无 FASTA 的参考: %s' % ', '.join(missing))
            # 写到运行目录，作为本次分析的临时参考
            tmp_ref = os.path.join(ctx.run_dir, 'ref_%s.fasta' % kv_run)
            n_seq = 0
            with open(tmp_ref, 'w', encoding='utf-8', newline='\n') as out:
                for acc, path in parts:
                    with open(path, encoding='utf-8', errors='replace') as fh:
                        for line in fh:
                            if line.startswith('>'):
                                n_seq += 1
                                head = line[1:].strip().split()[0]
                                # accession 唯一化：序列名已一致时保留原样
                                if head == acc:
                                    out.write(line if line.endswith('\n')
                                              else line + '\n')
                                else:
                                    out.write('>%s\n' % acc)
                            else:
                                out.write(line if line.endswith('\n')
                                          else line + '\n')
            log('kvsuite 参考库：%d 条（%d 个 accession）-> %s'
                % (n_seq, len(parts), os.path.basename(tmp_ref)))
            vfa = tmp_ref
        elif run_ref:
            src_run = check_path(os.path.join(_tool_runs_root(), run_ref),
                                 must_exist=True, in_platform=True)
            vfa = os.path.join(src_run, 'viral_contigs.fasta')
            if not os.path.isfile(vfa):
                raise RuntimeError('该运行无 viral_contigs.fasta（先跑 contig 分类）')
        else:
            if not fasta:
                # job 体在任务线程执行：abort(400) 只会变成
                # 「BadRequest: 400 Bad Request」，用户看不到提示文案
                raise ValueError('请选择 contigs 运行、kvsuite 参考或输入参考 FASTA')
            vfa = fasta

        # ---- reads ----
        # 固定重新比对。未手动指定时按 docstring 里的默认来源顺序解析；
        # kvsuite 的 bam/virus_reads 已按病毒筛过，不作为来源（低估变异）。
        # 文件 token 先收集再统一配对——R1/R2 拆在两个 token 里也能配上。
        reads = []
        file_toks = []
        for p in (ctx.p.get('reads') or '').split(','):
            p = p.strip()
            if not p:
                continue
            if os.path.isdir(p):
                reads.extend(find_read_pairs(p))
            elif os.path.isfile(p):
                file_toks.append(p)
            else:
                raise ValueError('reads 路径不存在：%s' % p)
        reads.extend(pair_read_files(file_toks))
        if not reads:
            if kv_run:
                # kvsuite 运行的输入序列（含样品与直填测序数据两种来路，
                # 都落在该 run 目录的 sample_sheet.tsv）
                sheet = os.path.join(_tool_runs_root(), kv_run,
                                     'sample_sheet.tsv')
                if not os.path.isfile(sheet):
                    raise ValueError('kvsuite 运行 %s 无 sample_sheet.tsv，'
                                     '请手动指定回贴 reads' % kv_run)
                with open(sheet, encoding='utf-8-sig', errors='replace') as fh:
                    hdr = fh.readline().rstrip('\n').split('\t')
                    if 'r1' not in hdr:
                        raise ValueError('sample_sheet.tsv 缺 r1 列：%s' % sheet)
                    i1 = hdr.index('r1')
                    i2 = hdr.index('r2') if 'r2' in hdr else len(hdr)
                    for line in fh:
                        cols = line.rstrip('\n').split('\t')
                        if len(cols) <= i1 or not cols[i1].strip():
                            continue
                        r2 = cols[i2].strip() if len(cols) > i2 else ''
                        reads.append((cols[i1].strip(), r2 or None))
                if not reads:
                    raise ValueError('kvsuite 运行 %s 的 sample_sheet.tsv 没有'
                                     '有效 reads，请手动指定' % kv_run)
                log('reads 默认取该 kvsuite 运行的输入序列：'
                    '%d 个样本（sample_sheet.tsv）' % len(reads))
            elif src_run:
                layer = os.path.join(src_run, CONSENSUS_FALLBACK_LAYER)
                reads = find_read_pairs(layer)
                if not reads:
                    raise RuntimeError(
                        '在 %s 下未找到 reads（期望 *_R1/*.fastq(.gz) 配对文件）'
                        % CONSENSUS_FALLBACK_LAYER)
                log('reads 默认取 contig 运行的去宿主后全量（%s）'
                    % CONSENSUS_FALLBACK_LAYER)
            else:
                raise ValueError('请直接给出回贴 reads 文件/目录')

        summary = consensus_and_variants(
            ctx.run_dir, vfa, reads, out_subdir='consensus',
            min_qual=min_qual, min_depth=min_depth, min_freq=min_freq,
            ambig=ambig, min_mapq=min_mapq, min_minor_freq=min_minor_freq,
            min_cov_pct=min_cov_pct, preset='sr', threads=ctx.threads,
            logger=logger, progress=lambda st, fr, msg: prog(st, fr, msg))
        logger.close()
        return summary
    return job


def _tool_job_kvsuite(ctx):
    """已知病毒识别与定量（known_virus_suite 五段整合）。

    移植自 D:/桌面/延伸基因组/MMPV-RNA/virome_analysis_pipeline：
      鉴定 → 过滤 → 共识 → 深度绘图 → 变异注释
    引擎分工固定：定量（identify）用 salmon（--writeBam 出映射 BAM），
    共识段内部固定 minibwa 真比对。

    默认只跑**两段**：识别/定量（identify）+ 过滤（filter）——本卡的口径就是
    「识别与定量」；共识序列、深度绘图、变异注释由下游模块各自完成
    （共识序列分析 t-consensus / 病毒变异分析 t-variant），需要在本卡一次跑全时
    前端把 stage 选成 all（或走 t-kvchain 一键全流程，那条链固定 stage=all）。
    与管道阶段 Virus_Platform_Core/kv_stage.py 的口径一致（它也只跑这两段）。

    caller 用 bcftools mpileup+call（freebayes/lofreq/ivar 在 Windows 上均不可得），
    参数与阈值为实测定稿，详见 Virus_Platform_Core/known_virus_suite/POSCOUNTS_REMOVAL_PLAN.md §4.2。
    """
    import json

    # ── 输入：勾选样品 → 临时 sample-sheet TSV(name,r1,r2) ──
    # 变异段（variant）不吃 reads，只需 BAM 或 VCF，因此允许 samples 为空。
    stage = (ctx.p.get('stage') or KV_DEFAULT_STAGE).strip().lower()
    if stage not in ('all', 'index', 'identify', 'filter', KV_DEFAULT_STAGE,
                     'consensus', 'plot', 'variant'):
        stage = KV_DEFAULT_STAGE
    samples_raw = (ctx.p.get('samples') or '').strip()
    sheet_in = (ctx.p.get('sample_sheet') or '').strip()
    reads_raw = (ctx.p.get('reads') or '').strip()
    if sheet_in and reads_raw:
        abort(400, '样本表 TSV 与直接输入测序数据只能二选一')
    if samples_raw and reads_raw:
        abort(400, '样品与直接输入测序数据只能二选一')
    if sheet_in:
        sheet = ctx.req('sample_sheet', '样本表 TSV')
    elif stage == 'variant':
        sheet = None          # 变异段不用样本表
    elif reads_raw:
        # 直接输入测序数据（不必先建样品）：目录/文件 → R1/R2 配对 → 样本表。
        # 样本名取 R1 文件名去掉 _R1 与扩展名，重名自动加后缀；
        # 文件 token 先收集再统一配对——R1/R2 拆在两个 token 里也能配上。
        from Virus_Platform_Core.consensus import (find_read_pairs,
                                                   pair_read_files)
        dir_pairs, file_toks = [], []
        for tok in reads_raw.split(','):
            tok = tok.strip()
            if not tok:
                continue
            if os.path.isdir(tok):
                pairs = find_read_pairs(tok)
                if not pairs:
                    abort(400, '目录下未找到 reads'
                               '（期望 *_R1/*.fastq(.gz)，R2 按 _R2 配对）：%s'
                               % tok)
                dir_pairs.extend(pairs)
            elif os.path.isfile(tok):
                file_toks.append(tok)
            else:
                abort(400, 'reads 路径不存在：%s' % tok)
        rows = []
        used = set()
        for r1, r2 in dir_pairs + pair_read_files(file_toks):
            from Virus_Platform_Core.consensus import R1_RX
            m = R1_RX.match(os.path.basename(r1))
            stem = (m.group(1) if m
                    else re.sub(r'\.(fq|fastq)(\.gz)?$', '',
                                os.path.basename(r1), flags=re.I))
            nm = re.sub(r'[^\w\-.]+', '_', stem).strip('_-.') or 'sample'
            base, k = nm, 2
            while nm in used:
                nm = '%s_%d' % (base, k)
                k += 1
            used.add(nm)
            rows.append((nm, r1, r2 or ''))
        sheet = os.path.join(ctx.run_dir, 'sample_sheet.tsv')
        with open(sheet, 'w', encoding='utf-8', newline='\n') as fh:
            fh.write('name\tr1\tr2\n')
            for nm, r1, r2 in rows:
                fh.write(f'{nm}\t{r1}\t{r2}\n')
    else:
        if not samples_raw:
            abort(400, '请选择样品、直接输入测序数据（FASTQ），或提供样本表 TSV')
        names = [s.strip() for s in samples_raw.split(',') if s.strip()]
        if not names:
            abort(400, '未选择有效样品')
        from Virus_Platform_Core.pipeline import load_sample_input
        from Virus_Platform_Core.utils import resolve_sample_name
        rows = []
        for nm in names:
            # 勾选框给的是 /api/samples 返回的真实目录名 → 按磁盘解析后再取输入
            real = resolve_sample_name(nm, DIRS['results'])
            sd = os.path.join(DIRS['results'], real)
            if not os.path.isdir(sd):
                # 先给出可读提示，而不是裸的「路径不存在: <绝对路径>」——
                # 2026-10-13 用户换测试数据删掉样品目录后，报错完全看不出
                # 是「样品没了」；示例按钮也不再依赖已注册样品（见 examples.js）
                abort(400, f'样品不存在: {nm}（run/results 下无目录 {real}；'
                           f'可在样品接入页重建，或改用「直接输入测序数据」）')
            sd = check_path(sd, must_exist=True, in_platform=True)
            r1, r2, _proj = load_sample_input(sd)
            if not r1 or not os.path.isfile(r1):
                abort(400, f'样品 {nm} 的 R1 不可用: {r1}')
            # sample_sheet 用真实目录名，便于与 results/ 下的产物目录对上
            rows.append((real, r1, r2 or ''))
        sheet = os.path.join(ctx.run_dir, 'sample_sheet.tsv')
        with open(sheet, 'w', encoding='utf-8', newline='\n') as fh:
            fh.write('name\tr1\tr2\n')
            for nm, r1, r2 in rows:
                fh.write(f'{nm}\t{r1}\t{r2}\n')

    # ── 参考序列与注释 ──
    # 缺省参考统一走 kv_stage 的目录解析（重构后的自包含库目录优先，
    # 旧布局散置文件兜底），不再各自硬编码路径。
    from Virus_Platform_Core.kv_stage import (default_ref_info,
                                              default_reference)
    # 选库优先（库与引擎正交：kv_lib 把 reference/ref_info/index_dir 三件套
    # 一起定死；engine 只决定定量方式）。给了 kv_lib 就不再吃零散的
    # reference/ref_info/index_dir 参数——避免"库选了 A、参考手填 B"的错配。
    kv_lib_raw = (ctx.p.get('kv_lib') or '').strip()
    ref = ref_info = None
    index_dir = None
    kv_lib_selected = False
    if kv_lib_raw:
        from Virus_Platform_Core.kv_stage import resolve_kv_lib
        try:
            lib = resolve_kv_lib(kv_lib_raw)
        except ValueError as e:
            abort(400, str(e))
        ref = lib['reference']
        ref_info = lib['ref_info'] or None
        kv_lib_selected = True   # 索引在 engine 确定后按库统一判定（下方③）
    if ref is None:
        ref = (ctx.p.get('reference') or '').strip()
        if ref:
            ref = ctx.req('reference', '参考序列')
        else:
            ref = default_reference()
            if not ref:
                abort(400, '缺少默认病毒参考：请先在「数据库构建 → 病毒鉴定库」建库')
            ref = check_path(ref, must_exist=True, in_platform=True)
    if ref_info is None:
        ref_info = (ctx.p.get('ref_info') or '').strip()
        if ref_info:
            ref_info = ctx.req('ref_info', '参考注释')
        else:
            cand = default_ref_info()
            ref_info = cand if cand and os.path.isfile(cand) else None

    # 定量引擎二选一：salmon（EM 定量，默认）| minibwa（真比对计数，
    # 2026-09-13 恢复）。库与引擎正交——engine 只决定定量方式，
    # 用哪个鉴定库由 kv_lib（或旧 reference/index_dir 参数）决定。
    # 共识段内部固定 minibwa 自建索引重新比对，不经此参数。
    engine = (ctx.p.get('engine') or 'salmon').strip().lower()
    if engine not in ('salmon', 'minibwa'):
        engine = 'salmon'       # 旧前端的历史取值一律归一为 salmon
    if engine == 'minibwa' and index_dir:
        # 建库页的 minibwa 索引布局在 <库>/minibwa/reference.*；
        # 引擎的 build_index/复用都认 <index_dir>/reference 前缀，
        # 选库 + minibwa 时把索引目录指到该子目录（预建索引直接复用）。
        mb = os.path.join(index_dir, 'minibwa')
        if os.path.isdir(mb):
            index_dir = mb
    # 索引复用：① kv_lib 选库 / ② 旧「鉴定库目录」参数（仅未选库时吃，
    # 避免"库选了 A、索引手填 B"错配）→ 布局判定都收敛到
    # kv_stage.index_dir_for（平台内可写库索引落库目录全局复用；外部库
    # 只读复用现成索引；都不满足则 <out>/index 现建，不污染库）。
    if index_dir is None and not kv_lib_selected:
        user_idx = (ctx.p.get('index_dir') or '').strip()
        if user_idx:
            index_dir = ctx.req('index_dir', '鉴定库目录')
    if index_dir is None:
        from Virus_Platform_Core.kv_stage import index_dir_for
        index_dir = index_dir_for(ref, engine)

    def _num(v, default, cast, lo=None, hi=None):
        try:
            x = cast(v)
        except (TypeError, ValueError):
            return default
        if lo is not None:
            x = max(lo, x)
        if hi is not None:
            x = min(hi, x)
        return x

    min_cov = _num(ctx.p.get('min_cov'), 10.0, float, 0.0, 100.0)
    min_depth = _num(ctx.p.get('min_depth'), 0.5, float, 0.0)
    min_reads = _num(ctx.p.get('min_reads'), 10, int, 0)
    min_poisson = _num(ctx.p.get('min_poisson'), 0.3, float, 0.0, 1.0)
    variant_qual = _num(ctx.p.get('variant_qual'), 3.5, float, 0.0)
    min_freq = _num(ctx.p.get('min_freq'), 0.05, float, 0.0, 1.0)
    aa_label_cutoff = _num(ctx.p.get('aa_label_cutoff'), 0.50, float, 0.0, 1.01)
    max_aa_labels = _num(ctx.p.get('max_aa_labels'), 40, int, 0)

    def job(log, prog, cancel):
        # 引擎入口统一走 kv_stage.engine_cmd()：源码模式 python -m，
        # 冻结模式 --run-engine（engine_entry 已在进程内注册 kvsuite）
        from Virus_Platform_Core.kv_stage import engine_cmd
        entry_exe, entry_pre = engine_cmd()

        out_dir = os.path.join(ctx.run_dir, 'kvsuite')
        os.makedirs(out_dir, exist_ok=True)
        cmd = [entry_exe] + entry_pre + [stage,
               '--out', out_dir,
               '--reference', ref,
               '--engine', engine,
               '--min-cov', str(min_cov),
               '--min-depth', str(min_depth),
               '--min-reads', str(min_reads),
               '--min-poisson', str(min_poisson),
               '--variant-qual', str(variant_qual),
               '--min-freq', str(min_freq),
               '--aa-label-cutoff', str(aa_label_cutoff),
               '--max-aa-labels', str(max_aa_labels)]
        if sheet:
            cmd += ['--sample-sheet', sheet]
        # 变异段输入：外部 VCF 优先于 BAM 目录（CLI 内也是这个优先级）
        input_vcf = (ctx.p.get('input_vcf') or '').strip()
        if input_vcf:
            cmd += ['--input-vcf', ctx.req('input_vcf', '变异 VCF 文件')]
        bam_dir = (ctx.p.get('bam_dir') or '').strip()
        if bam_dir:
            cmd += ['--bam-dir', ctx.req('bam_dir', 'BAM 目录')]
        if ref_info:
            cmd += ['--ref-info', ref_info]
        if index_dir:
            cmd += ['--index-dir', index_dir]
        if ctx.threads:
            cmd += ['--threads', str(ctx.threads),
                    '--align-threads', str(ctx.threads)]
        if ctx.p.get('ncbi_email'):
            cmd += ['--ncbi-email', str(ctx.p['ncbi_email']).strip()]
        if ctx.p.get('no_genes'):
            cmd.append('--no-genes')
        if ctx.p.get('no_variant_evo'):
            cmd.append('--no-variant-evo')
        if ctx.p.get('all_variants'):
            cmd.append('--all-variants')
        if ctx.p.get('no_snpgenie'):
            cmd.append('--no-snpgenie')

        log('$ ' + ' '.join(cmd))
        env = dict(os.environ)
        env.setdefault('PYTHONIOENCODING', 'utf-8')
        # NCBI API key 走环境变量传递（引擎侧 argparse 以 NCBI_API_KEY 兜底）：
        # 不进 argv → 不出现在上面这行 `$ ...` 任务日志 / run.log / 进程列表里
        _api_key = str(ctx.p.get('ncbi_api_key') or '').strip()
        if _api_key:
            env['NCBI_API_KEY'] = _api_key
        proc = subprocess.Popen(cmd, cwd=PLATFORM_ROOT,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT,
                                text=True, encoding='utf-8',
                                errors='replace', bufsize=1, env=env)
        # 全平台唯一一处裸 Popen（不经 run_cmd）：必须手工登记，否则取消时
        # TaskManager 的硬停止看不到它，引擎挂死也会永久占住 heavy 名额。
        task_register_proc(proc)
        # 超时看门狗：主线程阻塞在 stdout 读循环里，wait(timeout) 轮不到，
        # 由守护线程在超时后 taskkill 杀整棵树，读循环随即拿到 EOF。
        _timeout = DEFAULT_CMD_TIMEOUT_SEC or None
        _timed_out = threading.Event()
        _wd_stop = threading.Event()
        if _timeout:
            def _kill_watchdog():
                if not _wd_stop.wait(_timeout):
                    _timed_out.set()
                    log(f'[WARN] known_virus_suite 运行超过 {_timeout}s，'
                        f'强制终止进程树')
                    try:
                        subprocess.run(
                            ['taskkill', '/F', '/T', '/PID', str(proc.pid)],
                            capture_output=True)
                    except Exception:
                        try:
                            proc.kill()
                        except Exception:
                            pass
            threading.Thread(target=_kill_watchdog, daemon=True).start()
        tail = []
        try:
            for line in proc.stdout:
                line = line.rstrip('\n')
                if not line:
                    continue
                log(line)
                tail.append(line)
                if len(tail) > 40:
                    tail.pop(0)
                # cancel 是 threading.Event（见 TaskServer._run），必须用 is_set()，
                # 直接 cancel() 会报 'Event' object is not callable。
                if cancel is not None and cancel.is_set():
                    # taskkill /T 杀整棵进程树：只 kill 引擎本身会留下
                    # salmon/samtools 孤儿继续吃 CPU/磁盘
                    subprocess.run(['taskkill', '/F', '/T', '/PID', str(proc.pid)],
                                   capture_output=True)
                    proc.wait()
                    raise RuntimeError('用户取消')
        finally:
            _wd_stop.set()
        rc = proc.wait()
        if _timed_out.is_set():
            raise RuntimeError(
                f'known_virus_suite 超时（运行 {_timeout}s 后被强制终止）；'
                f'末几行：\n' + '\n'.join(tail[-8:]))
        if rc != 0:
            raise RuntimeError(f'known_virus_suite 退出码 {rc}；末几行：\n'
                               + '\n'.join(tail[-8:]))

        # ── 汇总产物（供结果面板展示）──
        result = {'stage': stage, 'engine': engine, 'out': out_dir,
                  'reference': ref}
        fsum = os.path.join(out_dir, 'filter', 'filtered.tsv')
        if os.path.isfile(fsum):
            result['filtered_tsv'] = fsum
            try:
                with open(fsum, encoding='utf-8', errors='replace') as fh:
                    result['n_confirmed'] = max(0, sum(1 for _ in fh) - 1)
            except OSError:
                pass
        vsum = os.path.join(out_dir, 'variant_summary.json')
        if os.path.isfile(vsum):
            result['variant_summary'] = vsum
            try:
                with open(vsum, encoding='utf-8') as fh:
                    data = json.load(fh)
                rows = data.get('results') if isinstance(data, dict) else data
                rows = rows or []
                result['n_variant_genomes'] = len(rows)
                result['n_variants'] = sum(int(d.get('n_variants') or 0)
                                           for d in rows)
                if isinstance(data, dict):
                    result['n_variant_skipped'] = int(data.get('n_skipped') or 0)
            except (OSError, ValueError, TypeError):
                pass
        return result
    return job


def _tool_job_quicktree(ctx):
    """快速建树：多条序列 MAFFT 全长比对 → NJ（纯 Python）/ FastTree。

    aligned=True 时输入已是比对好的 FASTA（如序列比对模块的
    aln.fasta / aln.trim.fasta），跳过 MAFFT 直接建树——对已比对序列
    重复比对会破坏列对应关系。
    """
    seqs_fa = ctx.req('seqs', '序列 FASTA')
    method = (ctx.p.get('method') or 'nj').strip().lower()
    if method not in ('nj', 'fasttree'):
        abort(400, '建树方法仅支持 nj / fasttree')
    aligned = str(ctx.p.get('aligned') or '').strip().lower() in (
        '1', 'true', 'yes', 'on')
    max_n = max(3, min(int(ctx.p.get('max_n') or 100), 500))

    def _short(h, idx):
        name = re.split(r'[\s|]', (h or '').strip())[0][:40]
        name = re.sub(r'[^A-Za-z0-9_\-.]', '_', name)
        return name or f'seq{idx}'

    def job(log, prog, cancel):
        logger = TaskLogger(callback=log)
        from Virus_Platform_Core.phylo import _run_mafft, _run_nj, _run_fasttree
        from Virus_Platform_Core.utils import iter_fasta

        prog('read', 0.05, '读取与筛选序列')
        recs, seen = [], {}
        for h, s in iter_fasta(seqs_fa):
            name = _short(h, len(recs) + 1)
            if name in seen:
                seen[name] += 1
                name = f'{name}_{seen[name]}'
            else:
                seen[name] = 0
            recs.append((name, s.upper()))
        if len(recs) < 2:
            raise RuntimeError('FASTA 中少于 2 条序列，无法建树')
        recs = recs[:max_n]

        capped_fa = os.path.join(ctx.run_dir, 'input.fasta')
        with safe_open(capped_fa, 'wt') as f:
            for name, s in recs:
                f.write(f'>{name}\n')
                for i in range(0, len(s), 70):
                    f.write(s[i:i + 70] + '\n')

        prog('aln', 0.2, '输入已比对，跳过 MAFFT' if aligned
             else 'MAFFT 全长比对')
        if aligned:
            aln = capped_fa
        else:
            aln = _run_mafft(capped_fa, os.path.join(ctx.run_dir, 'aln.fasta'),
                             threads=ctx.threads, logger=logger)
        if aligned:
            logger.log('aligned=True：输入视为已完成的多序列比对，'
                       '直接建树（不重复 MAFFT）')

        prog('tree', 0.75, 'NJ 建树' if method == 'nj' else 'FastTree 建树')
        tree_name = 'nj.nwk' if method == 'nj' else 'tree.nwk'
        if method == 'nj':
            _run_nj(aln, os.path.join(ctx.run_dir, tree_name), logger=logger)
        else:
            _run_fasttree(aln, os.path.join(ctx.run_dir, tree_name),
                          logger=logger)

        prog('done', 1.0, '完成')
        logger.close()
        return {'n_seqs': len(recs), 'method': method, 'tree': tree_name,
                'aligned': aligned}
    return job


def _tool_job_align(ctx):
    """序列比对（独立模块）：MAFFT（auto/L-INS-i/fast）→ trimAl 清剪。

    输入 FASTA（核酸或蛋白，自动判别；≥2 条）。可选 seqs_extra：
    第二个 FASTA 文件，与 seqs 合并后一起参与比对（同名序列自动加 _2 后缀）。
    产物：input.fasta / aln.fasta / aln.trim.fasta（trimAl 关闭时无）/
    summary.json，可在「比对查看器」彩色浏览与编辑。
    """
    from Virus_Platform_Core.utils import iter_fasta, write_fasta_record
    seqs_fa = ctx.req('seqs', '序列 FASTA')
    seqs_extra = ctx.opt('seqs_extra')
    strategy = (ctx.p.get('strategy') or 'auto').strip().lower()
    if strategy not in ('auto', 'linsi', 'fast'):
        abort(400, '比对策略仅支持 auto / linsi / fast')
    trimal = (ctx.p.get('trimal') or 'automated1').strip().lower()
    if trimal not in ('automated1', 'gappyout', 'strict', 'none'):
        abort(400, 'trimAl 方法仅支持 automated1 / gappyout / strict / none')
    max_n = max(3, min(int(ctx.p.get('max_n') or 200), 500))

    def job(log, prog, cancel):
        logger = TaskLogger(callback=log)
        from Virus_Platform_Core.phylo import _run_mafft, _run_trimal
        from Virus_Platform_Core.sdt_exact import detect_seqtype
        prog('read', 0.05, '读取与筛选序列')
        recs, seen = [], {}
        _srcs = [seqs_fa] + ([seqs_extra] if seqs_extra else [])
        for _si, _src in enumerate(_srcs):
            for h, s in iter_fasta(_src):
                name = re.split(r'[\s|]', (h or '').strip())[0][:60] or \
                    f'seq{len(recs) + 1}'
                name = re.sub(r'[^A-Za-z0-9_\-.]', '_', name)
                if name in seen:
                    seen[name] += 1
                    name = f'{name}_{seen[name]}'
                else:
                    seen[name] = 0
                recs.append((name, s.upper()))
        if seqs_extra:
            logger.log(f'额外序列文件已合并: {os.path.basename(seqs_extra)}'
                       f'（合并后共 {len(recs)} 条）')
        if len(recs) < 2:
            raise RuntimeError('FASTA 中少于 2 条序列，无法比对')
        recs = recs[:max_n]
        seqtype = detect_seqtype([s for _, s in recs])
        capped = os.path.join(ctx.run_dir, 'input.fasta')
        with safe_open(capped, 'wt') as f:
            for name, s in recs:
                write_fasta_record(f, name, s)
        logger.log(f'序列类型判定: {"蛋白(aa)" if seqtype == "aa" else "核酸(nt)"}'
                   f'，{len(recs)} 条参与比对')
        prog('align', 0.15, f'MAFFT 比对（{strategy}）')
        aln = _run_mafft(capped, os.path.join(ctx.run_dir, 'aln.fasta'),
                         threads=ctx.threads, logger=logger, strategy=strategy)
        trim_info = {'applied': False}
        aln_used = aln
        if trimal != 'none':
            prog('trim', 0.7, f'trimAl 清剪（{trimal}）')
            aln_used, trim_info = _run_trimal(
                aln, os.path.join(ctx.run_dir, 'aln.trim.fasta'),
                logger=logger)
        prog('done', 1.0, '完成')
        logger.close()
        return {'n_seqs': len(recs), 'seqtype': seqtype, 'strategy': strategy,
                'trimal': trimal, 'trim': trim_info,
                'aln': 'aln.fasta',
                'aln_used': os.path.basename(aln_used)}
    return job


def _tool_job_sdt(ctx):
    """SDT 精确分析（纯 Python 复刻 SDTv1.3，替代外部 SDT exe）：
    逐对 MAFFT 独立比对 → Get_Similarity 公式 → 聚类排序热图 + 分布图。
    aligned=True 时输入为已比对 MSA，跳过比对器直接按公式计算；
    seqtype='auto' 自动判别核酸/蛋白（MMPV 双轨口径），蛋白按 AA 同一性。"""
    seqs_fa = ctx.req('seqs', '序列 FASTA')
    max_n = max(3, min(int(ctx.p.get('max_n') or 30), 200))
    orient = bool(ctx.p.get('orient', True))
    seqtype = (ctx.p.get('seqtype') or 'auto').strip().lower()
    if seqtype not in ('auto', 'nt', 'aa'):
        seqtype = 'auto'
    palette = (ctx.p.get('palette') or 'sdt').strip().lower()
    if palette not in ('sdt', 'cividis', 'viridis', 'RdYlBu', 'Spectral',
                       'YlGnBu', 'coolwarm', 'magma'):
        palette = 'sdt'
    aligned = bool(ctx.p.get('aligned'))

    def job(log, prog, cancel):
        logger = TaskLogger(callback=log)
        mafft = cfg.tool('mafft')
        from Virus_Platform_Core.sdt_exact import run_sdt_exact
        prog('sdt', 0.02, 'SDT 精确分析（MAFFT 逐对独立比对，SDT v1.3 口径）')
        res = run_sdt_exact(seqs_fa, ctx.run_dir, mafft, max_n=max_n,
                            orient=orient, threads=ctx.threads,
                            seqtype=seqtype,
                            palette=palette, logger=logger,
                            progress=lambda p, m: prog('sdt', p, m),
                            cancel=cancel, aligned=aligned)
        logger.close()
        return res
    return job


def _tool_job_identity(ctx):
    """核苷酸+氨基酸同一性表（BioAider Sequence Identity Matrix 口径）：
    NT 矩阵 + AA 矩阵（可选输入或最长 ORF 翻译）+ 复合热图（NT 上 /
    AA 下）+ 逐对同一性长表。"""
    nt_fa = ctx.req('nt_seqs', '核苷酸 FASTA')
    aa_fa = ctx.opt('aa_seqs')
    max_n = max(3, min(int(ctx.p.get('max_n') or 30), 200))
    palette = (ctx.p.get('palette') or 'sdt').strip().lower()
    if palette not in ('sdt', 'cividis', 'viridis', 'RdYlBu', 'Spectral',
                       'YlGnBu', 'coolwarm', 'magma'):
        palette = 'sdt'
    aligned = bool(ctx.p.get('aligned'))

    def job(log, prog, cancel):
        logger = TaskLogger(callback=log)
        mafft = cfg.tool('mafft')
        from Virus_Platform_Core.sdt_exact import run_identity_table
        prog('idty', 0.02, '核苷酸+氨基酸同一性表（BioAider 口径）')
        res = run_identity_table(nt_fa, ctx.run_dir, mafft, aa_fasta=aa_fa,
                                 max_n=max_n, aligned=aligned,
                                 threads=ctx.threads, palette=palette,
                                 logger=logger,
                                 progress=lambda p, m: prog('idty', p, m),
                                 cancel=cancel)
        logger.close()
        return res
    return job


def _load_fasta_seqs(path):
    """FASTA → {首 token id: 序列}（单条记录可达 MB 级，仅限小文件用）。"""
    import gzip
    op = gzip.open if str(path).lower().endswith('.gz') else open
    out, name, buf = {}, None, []
    with op(path, 'rt', errors='replace') as f:
        for line in f:
            if line.startswith('>'):
                if name is not None:
                    out[name] = ''.join(buf)
                name = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
    if name is not None:
        out[name] = ''.join(buf)
    return out


def _find_contig_seq(run_dir, contig):
    """在工具运行目录的各 FASTA 中查找 contig 序列。"""
    for cur, _sub, fns in os.walk(run_dir):
        if os.path.join('analysis') in cur:
            continue
        for fn in fns:
            if not fn.lower().endswith(('.fasta', '.fa', '.fna', '.fas')):
                continue
            try:
                seqs = _load_fasta_seqs(os.path.join(cur, fn))
            except (OSError, ValueError):
                continue
            if contig in seqs:
                return seqs[contig]
    return None


def _extract_records(src, ids, dst):
    """从 FASTA/FASTQ（支持 .gz）提取 header 首 token 命中 ids 的记录。

    kunpeng 报告里的 read ID 经过 seqkit fq2fa 转换，已去掉 Illumina 配对
    后缀 /1、/2，匹配时对原始 header 兼容带/不带后缀两种写法；
    FASTQ 输入统一转成 FASTA 写出（与 .fasta 扩展名一致）。
    返回提取条数；ids 为空集时写出空文件。
    """
    import gzip
    name = str(src).lower()
    is_fq = name.endswith(('.fastq.gz', '.fq.gz', '.fastq', '.fq'))
    op = gzip.open if name.endswith('.gz') else open

    def match(rid):
        """命中返回规范 ID（配对后缀 /1、/2 已去掉，与 kunpeng 报告一致），未命中返回 None。"""
        if rid in ids:
            return rid
        base = rid.rsplit('/', 1)
        if len(base) == 2 and base[1] in ('1', '2') and base[0] in ids:
            return base[0]
        return None

    n = 0
    with op(src, 'rt', errors='replace') as f, safe_open(dst, 'wt') as w:
        if is_fq:
            while True:
                h = f.readline()
                if not h:
                    break
                seq = f.readline().rstrip('\r\n')
                f.readline()
                f.readline()
                tok = h[1:].split()
                mid = match(tok[0]) if tok else None
                if mid:
                    w.write(f'>{mid}\n')
                    for i in range(0, len(seq), 70):
                        w.write(seq[i:i + 70] + '\n')
                    n += 1
        else:
            keep = False
            for line in f:
                if line.startswith('>'):
                    tok = line[1:].split()
                    keep = bool(tok) and match(tok[0]) is not None
                    if keep:
                        n += 1
                if keep:
                    w.write(line)
    return n


def _qc_brief(qc):
    """序列预筛摘要（供 events.json / summary.json / 前端展示）。

    逐条明细（每条为什么被剔）留在 `align_qc_report.tsv` 里——几百条时
    没必要塞进 JSON。未启用预筛时 qc 为 None。
    """
    if not qc:
        return None
    return {'n_in': qc['n_in'], 'n_keep': qc['n_keep'],
            'n_removed': qc['n_removed'], 'reasons': qc['reasons'],
            'ref_id': qc['ref_id'], 'ref_len': qc['ref_len'],
            'thresholds': qc['thresholds'],
            'clean': os.path.basename(qc['clean']),
            'report': os.path.basename(qc['report'])}


def _tool_job_rdp(ctx):
    """RDP 重组分析（RDP5CL 九方法，本机原生）。

    输入 = 比对好的 FASTA（MAFFT 产物）。产物落 run_dir/rdp/：
      rdp5.csv（主表）/ 断点分布 / events.json / events.tsv / summary.json
      + 预筛产物 align_qc_report.tsv / clean.fasta / removed.fasta

    序列预筛（2026-09-15 增，默认开）：进 RDP5 前先剔掉**不完整片段**与低质量
    序列，口径对齐参照管道 `virome_phylo_pipeline/utils/align_qc.py`——长度比
    <0.90（相对参考/中位长度）、gap >10%、N >5%、与参考 identity 偏离 >2%。
    理由：RDP5 检测的是同源等长序列间的重组，片段混进来会让全比对出现大片
    gap，断点坐标与 p 值都被带偏，而片段本身是缺失数据、任何重组信号都不可信。
    `qc=0` 可关（只用于排障），`qc_ref` 指定参考序列 ID，`qc_min_len` 改长度比阈值。

    引擎说明：2026-09-15 移除自研三序列法（MaxChi/Chimaera/Bootscan）——
    实测在合成嵌合数据（6 条 1kb，真重组体 + 真断点 500）上 18 报 1 真：
    Bootscan 分支只看 support 不看 p（p=1.0 的事件照样入选）、start/end 取
    三方法断点极值而非片段区间、全组合扫描无多重校正。与 RDP5 并存只会让
    用户拿到假阳性事件表，故删除，RDP5 为唯一引擎。
    """
    inp = ctx.req('input', '比对 FASTA')
    # 预筛默认**开**：完整性过滤是重组分析的前提，不是可选项
    qc_on = str(ctx.p.get('qc', '1') or '').strip().lower() not in (
        '0', 'false', 'no', 'off')
    qc_ref = (ctx.p.get('qc_ref') or '').strip()
    try:
        qc_min_len = float(ctx.p.get('qc_min_len') or 0.90)
    except (TypeError, ValueError):
        qc_min_len = 0.90
    qc_min_len = max(0.5, min(qc_min_len, 1.0))

    def job(log, prog, cancel):
        logger = TaskLogger(callback=log)
        from Virus_Platform_Core import rdp5_engine as _r5
        from Virus_Platform_Core.align_qc import MIN_SEQS, qc_alignment
        from Virus_Platform_Core.utils import load_alignment, variable_cols
        out_dir = os.path.join(ctx.run_dir, 'rdp')
        os.makedirs(out_dir, exist_ok=True)

        # ① 序列预筛：只把完整序列送进 RDP5
        aln_used, qc = inp, None
        if qc_on:
            prog('rdp5', 0.03, '序列预筛（去不完整片段 / 低质量序列）')
            try:
                qc = qc_alignment(inp, out_dir, reference_id=qc_ref or None,
                                  min_length=qc_min_len, logger=logger)
            except ValueError as exc:
                log(f'[错误] {exc}')
                raise
            if not qc['enough']:
                msg = (f"预筛后只剩 {qc['n_keep']} 条序列（RDP5 至少需要 {MIN_SEQS} 条），"
                       f"共剔除 {qc['n_removed']} 条：{qc['reasons']}。"
                       "请调低「最小长度比」、指定参考序列，或补充更多完整序列。")
                log(f'[错误] {msg}')
                raise RuntimeError(msg)
            if qc['n_removed']:
                aln_used = qc['clean']

        # 比对规模：RDP5 自身不产出，但结果摘要要显示（按 RDP5 实际吃到的
        # 那份比对统计，才能和事件数对得上）
        n_seqs, aln_len, n_var = None, None, None
        try:
            aln = load_alignment(aln_used)
            n_seqs = len(aln)
            if aln:
                seqs = list(aln.values())
                aln_len = len(seqs[0])
                n_var = len(variable_cols(seqs))
        except Exception as exc:
            log(f'[warn] 比对规模统计失败：{type(exc).__name__}: {exc}')

        prog('rdp5', 0.05, 'RDP5 检测中（九方法，本机原生）')
        csv_path = _r5.run_rdp5cl(aln_used, out_dir, prefix='rdp5',
                                  timeout=3600, logger=logger,
                                  cancel=cancel)
        prog('rdp5', 0.9, '解析事件')
        raw_evs = _r5.parse_rdp5_csv(csv_path)
        evs = []
        for e in raw_evs:
            sig = sorted((m for m, p in e['methods'].items()
                          if p < 0.05), key=lambda m: e['methods'][m])
            evs.append({
                'recombinant': e['recombinant'],
                'major': e['major_parent'], 'minor': e['minor_parent'],
                'start': e['bp_start'], 'end': e['bp_end'],
                'methods': sig, 'n_methods': e['n_methods'],
                'pvalue': (min(e['methods'].values())
                           if e['methods'] else None),
                'support': None,
            })
        sep = '\t'
        with open(os.path.join(out_dir, 'events.tsv'), 'w',
                  encoding='utf-8', newline='') as f:
            f.write(sep.join(['recombinant', 'major', 'minor', 'start', 'end',
                              'n_methods', 'methods', 'best_pvalue']) + '\n')
            for e in evs:
                f.write(sep.join([
                    e['recombinant'], e['major'], e['minor'],
                    str(e['start']), str(e['end']), str(e['n_methods']),
                    ','.join(e['methods']),
                    '' if e['pvalue'] is None else f"{e['pvalue']:.3g}",
                ]) + '\n')
        # events.json 是前端结果表的数据源，须带齐摘要字段
        # （前端读 d.n_seqs / d.aln_len / d.variable / d.events / d.qc）
        with open(os.path.join(out_dir, 'events.json'), 'w',
                  encoding='utf-8') as f:
            json.dump({'engine': 'rdp5', 'n_seqs': n_seqs,
                       'aln_len': aln_len, 'variable': n_var,
                       'n_events': len(evs), 'events': evs,
                       'qc': _qc_brief(qc)},
                      f, ensure_ascii=False, indent=1)
        # summary.json 的 input 是 /api/tool/rdp_masked 回读原比对的依据：
        # 必须记 **RDP5 实际吃到的那份比对**（预筛后的 clean.fasta），
        # 否则掩蔽坐标与事件坐标对不上。原始输入另存 source_input 备查。
        # 断点分布兜底（2026-10-12）：RDP5CL 的 -rbdp 只统计双亲本都已
        # 鉴定的事件——事件带 Unknown 亲本时（真实数据常见）分布 CSV
        # 静默全 0、前端画出 0 平线。此时用事件断点重建，并在摘要注明。
        bdp_csv = os.path.join(out_dir, _r5.BDP_CSV)
        rdp_note = None
        if evs and _r5.breakpoint_distribution_degenerate(bdp_csv):
            _r5.synthesize_breakpoint_distribution(bdp_csv, evs, aln_len)
            rdp_note = ('断点分布由平台按检出事件的断点重建（RDP5 对含'
                        '未鉴定亲本 Unknown 的事件不产出断点分布）')
        # 基因图谱（2026-10-12 增，可选）：RDP5CL 命令行不产出 ORFCoords，
        # 由用户给 Rgrapher 风格 CSV（Gene.symbol,Start,Stop，参考序列
        # nt 坐标），经比对参考行的 gap 轮廓映射到比对坐标，前端画条带。
        gene_map = None
        orf_in = (ctx.p.get('orf_csv') or '').strip()
        if orf_in:
            try:
                gene_map = _r5.orf_csv_to_aln_map(
                    orf_in, aln_used, qc_ref or None)
                with open(os.path.join(out_dir, 'gene_map.json'), 'w',
                          encoding='utf-8') as f:
                    json.dump(gene_map, f, ensure_ascii=False, indent=1)
                if not gene_map['genes']:
                    rdp_note = (rdp_note + '；' if rdp_note else '')
                    rdp_note += ('基因坐标全部无法映射（参考序列坐标越界？'
                                 '应按 qc_ref 参考序列给 GenBank 坐标）')
            except Exception as exc:
                log(f'[warn] 基因坐标解析失败：{type(exc).__name__}: {exc}')
                rdp_note = (rdp_note + '；' if rdp_note else '')
                rdp_note += f'基因坐标 CSV 解析失败：{exc}'
                gene_map = None
        summary = {'tool': 'rdp', 'engine': 'rdp5',
                   'n_seqs': n_seqs, 'aln_len': aln_len,
                   'variable': n_var, 'n_events': len(evs),
                   'events': evs[:200], 'csv': csv_path,
                   'qc': _qc_brief(qc),
                   'input': aln_used, 'source_input': inp, 'out_dir': out_dir}
        if gene_map:
            summary['gene_map'] = gene_map
        if rdp_note:
            summary['note'] = rdp_note
        with open(os.path.join(out_dir, 'summary.json'), 'w',
                  encoding='utf-8') as f:
            json.dump(summary, f, ensure_ascii=False, indent=1)
        logger.close()
        return summary
    return job


def _tool_job_phylogeo(ctx):
    """系统地理分析（建树 + Fitch 迁移重构）+ 可选的两项统计。

    输入 = 比对 FASTA + 元数据 CSV（或 FASTA 头 `>acc|区域|年份`）。
    产物落 run_dir/phylogeo/：tree_annotated.nwk / migration_matrix.tsv /
    leaf_states.tsv / transitions.tsv（+ 迁移事件分类 transition_classes.tsv；
    给了坐标表/元数据带 lat/lon 时另有 coords_resolved.tsv 与 sample_coords.tsv；
    开启统计时的 rrt_permutation.tsv / weight_bands.json）。

    可选统计（默认关）：
      rrt_perm   区域随机化检验的置换次数（0=不做；建议 ≥1000）→ 区划是否真有信号
      rssp_bs    bootstrap 复本数（0=不做；建议 ≥100）→ 迁移数的置信区间
                 同时是「后验权重分带」的驱动源（方向支持率）
      motp_bin   迁移随时间（MOTP）的分箱宽度（年）；>0 才做
      dated_tree LSD2 的 `.date.nexus`（给了就用真正的节点年代，否则退回近似）
      coords     坐标表（可选）：逐样本点层 + 区划落点（样本坐标中位数优先）
      beast_log / beast_trees  BEAST 产物（轨道 B 回读）：服务器跑完的 `.log` /
                 `.trees`（或 MCC 树），本地**只读导入**、不跑 MCMC（见
                 beast_handoff / beast_io）。两个都不给 = 不做 BEAST 面板；
                 导入失败只影响 BEAST 面板，主分析（树/矩阵/统计）照常。
                 beast_burnin / beast_min_prob / beast_max_trees 为导入参数。
    """
    inp = ctx.req('input', '比对 FASTA')
    meta = (ctx.p.get('meta') or '').strip()
    # 2026-09-18：区划列**留空就是自动识别**（explorer 新版叫 Geo_Location/
    # Country，硬塞 'region' 会白报一条「自动识别」告警）
    trait = (ctx.p.get('trait') or '').strip()
    method = (ctx.p.get('method') or 'nj').strip().lower()
    if method not in ('nj', 'fasttree'):
        method = 'nj'
    dated_tree = ctx.opt('dated_tree')          # 可选文件（不存在则视为未给）
    coords_path = ctx.opt('coords')             # 可选：坐标表（样本点层）
    # BEAST 产物（轨道 B 回读）：服务器跑完的 `.log` / `.trees`，本地**只读导入**。
    # 两个都不给 = 不做 BEAST 面板（与历史行为一致）。契约见 beast_handoff.py /
    # beast_io.py；这里是**唯一**把用户选的产物接到 analyze() 的调用方。
    beast_log = ctx.opt('beast_log')
    beast_trees = ctx.opt('beast_trees')
    date_trait = (ctx.p.get('date_trait') or 'year').strip() or 'year'
    try:
        motp_bin = float(str(ctx.p.get('motp_bin') or 0).strip() or 0)
    except ValueError:
        motp_bin = 0.0
    if motp_bin < 0:
        motp_bin = 0.0

    def _int(key, lo, hi, default=0):
        try:
            v = int(str(ctx.p.get(key) or default).strip() or default)
        except ValueError:
            return default
        return max(lo, min(v, hi))

    rrt_perm = _int('rrt_perm', 0, 20000, 0)
    rssp_bs = _int('rssp_bs', 0, 2000, 0)
    seed = _int('seed', 0, 10 ** 9, 0)

    def _fnum(key, lo, hi, default):
        """BEAST 导入的浮点参数（越界/非法一律钳回默认；burnin 上界 0.9）。"""
        try:
            v = float(str(ctx.p.get(key) or default).strip() or default)
        except (TypeError, ValueError):
            return default
        return max(lo, min(v, hi))

    # 只有给了 beast_log / beast_trees 才有意义（analyze 内部据此决定是否导入）
    beast_burnin = _fnum('beast_burnin', 0.0, 0.9, 0.1)
    beast_min_prob = _fnum('beast_min_prob', 0.0, 1.0, 0.0)
    beast_max_trees = _int('beast_max_trees', 0, 200000, 2000)
    # 样本点层进 summary 的上限（逐样本 1 条；超限不进轮询载荷，条数如实上报，
    # 全量永远在 sample_coords.tsv 里 —— 不允许静默丢点）
    SP_MAX = 3000
    # 迁移事件分类明细进 summary 的上限（逐边 1 条；同上，全量在
    # transition_classes.tsv / transition_classes.json）
    TC_MAX = 300
    # 日志与 summary 里距离分带的固定顺序（phylogeo._CLASS_ORDER 的显示副本）
    _TC_LABELS = ('Direct', 'Indirect', 'Distant', 'Unresolved')

    def job(log, prog, cancel):
        logger = TaskLogger(callback=log)
        from Virus_Platform_Core.phylogeo import analyze
        out_dir = os.path.join(ctx.run_dir, 'phylogeo')
        os.makedirs(out_dir, exist_ok=True)
        tree_path = os.path.join(out_dir, 'tree.nwk')
        res = analyze(inp, tree_path,
                      meta_path=meta or None, trait=trait, method=method,
                      out_dir=out_dir, rrt_perm=rrt_perm, rssp_bs=rssp_bs,
                      seed=seed, dated_tree=dated_tree, motp_bin=motp_bin,
                      date_trait=date_trait, coords_path=coords_path,
                      beast_log=beast_log, beast_trees=beast_trees,
                      beast_burnin=beast_burnin, beast_min_prob=beast_min_prob,
                      beast_max_trees=beast_max_trees,
                      progress=lambda st, f_, m: prog(st, f_, m))
        # ---- 坐标命中情况（analyze 里解析；这里只做日志，不静默）----
        # 坏行/未命中区划/元数据坏坐标都**不静默**：skipped 与 missing 既进 job
        # 日志，也进 summary（前端显示出来，避免"图里少画了却看不出"）。
        coords = res.get('coords') or {}
        coord_src = res.get('coord_src') or {}
        coord_missing = res.get('coord_missing') or []
        coord_skips = res.get('coord_skips') or []
        if coords_path:
            for s in coord_skips:
                log(f'[WARN] 坐标表 {os.path.basename(coords_path)}：{s}')
            log(f"坐标表 {os.path.basename(coords_path)}：解析出 "
                f"{res.get('coord_table_regions', 0)} 个区划"
                f"（逐行条目 {res.get('coord_table_rows', 0)} 条，"
                f"跳过 {len(coord_skips)} 行）")
        n_sp = int(res.get('n_sample_coords') or 0)
        if n_sp:
            n_meta_pt = sum(1 for p in (res.get('sample_points') or [])
                            if p.get('src') == 'meta')
            log(f"逐样本经纬度：{n_sp}/{res['n_seqs']} 条样本有坐标"
                f"（元数据列 {n_meta_pt} / 坐标表 {n_sp - n_meta_pt}）")
        n_samp = sum(1 for v in coord_src.values() if v == 'samples')
        n_user = sum(1 for v in coord_src.values() if v == 'user')
        # `*_prefix`＝退到「国家: 子区」的**国家部分**才命中的（落点粗一级）。
        # 单独计数、单独报 —— 混进"用户表/内置质心表"里就等于把粗坐标说成准坐标。
        n_usr_px = sum(1 for v in coord_src.values() if v == 'user_prefix')
        n_bui_px = sum(1 for v in coord_src.values() if v == 'builtin_prefix')
        if res['regions']:
            log(f"区划落点坐标：{len(coords)}/{len(res['regions'])} 个区划有坐标"
                f"（样本坐标中位数 {n_samp} / 用户表 {n_user} / "
                f"内置质心表 {len(coords) - n_samp - n_user - n_usr_px - n_bui_px}）"
                + (f"，其中 {n_usr_px + n_bui_px} 个是按「国家: 子区」的"
                   f"**国家部分**兜底的（落点粗一级）"
                   if (n_usr_px or n_bui_px) else '')
                + (f"，缺坐标：{'、'.join(coord_missing[:8])}" if coord_missing else ''))
        if coord_missing:
            log('[WARN] 以下区划没有坐标，迁移弧线与距离分类都用不上它们：'
                + '、'.join(coord_missing))
        # 迁移事件分类：全量明细落 transition_classes.json（另有一份
        # transition_classes.tsv 由 analyze 写），summary 只放精简版 +
        # 截断前 TC_MAX 条（前端表格与弧线着色用）。
        tc = res.get('transition_classes')
        tc_brief = None
        if tc:
            with open(os.path.join(out_dir, 'transition_classes.json'), 'w',
                      encoding='utf-8') as f:
                json.dump(tc, f, ensure_ascii=False, indent=1)
            th = tc['thresholds']
            cnt = tc['counts']
            log('迁移事件分类（距离分带）：'
                + ' / '.join(f'{k} {cnt[k]}' for k in _TC_LABELS if cnt.get(k))
                + f"（跨区边 {tc['n_edges']} 条；≤{th['q33']} km 直接 / "
                  f"≥{th['q67']} km 远距，阈值来源 {th['src']}）")
            if tc['flags']:
                log('分类证据标记：'
                    + ' / '.join(f'{k} {v}' for k, v in tc['flags'].items())
                    + f"（低分化远程=边替换数 ≤ 本次中位数 {tc['snp_median']}；"
                      '同期远程=两侧年份中位数差 ≤1 年）')
            if tc['n_unknown_state']:
                log(f"{tc['n_unknown_state']} 条跨区边一端状态未知（Unknown）——"
                    '未进分类表')
            rows_all = tc['rows']
            if len(rows_all) > TC_MAX:
                log(f'[WARN] 迁移事件分类明细只把前 {TC_MAX} 条放进 summary'
                    f'（共 {len(rows_all)} 条），全量见 '
                    f'{os.path.join(out_dir, "transition_classes.tsv")}')
            tc_brief = {
                'counts': tc['counts'], 'thresholds': tc['thresholds'],
                'flags': tc['flags'], 'n_edges': tc['n_edges'],
                'n_intro': tc['n_intro'], 'n_no_tips': tc['n_no_tips'],
                'n_unknown_state': tc['n_unknown_state'],
                'n_state_changes': tc['n_state_changes'],
                'snp_median': tc['snp_median'],
                'n_with_step': tc['n_with_step'],
                'rows': rows_all[:TC_MAX],
                'rows_truncated': max(0, len(rows_all) - TC_MAX),
                'method': tc['method'],
            }
        _m = res.get('motp') or {}
        if _m.get('bins'):
            log(f"迁移随时间（MOTP）：{len(_m['bins'])} 个时间窗"
                f"（轴向 {_m.get('axis')}／窗宽 {_m.get('bin_width')} 年，"
                f"可计年迁移 {len(_m.get('rows') or [])} 条"
                + (f"，{_m['n_undated']} 条无年代未分箱" if _m.get('n_undated')
                   else '') + '）')
        # 权重分带的驱动源要进日志：本平台的驱动源只有 RSPP 支持率（不是贝叶斯
        # 后验，也不是 BEAST 跳转计数）；事后查报告必须能分清数字是哪一支口径。
        _wb = res.get('weight_bands')
        if _wb:
            log('后验权重分带（%s）：%d 条走廊 / %d 事件 / 期望 %.3f 条'
                % (_wb.get('driver_label') or _wb.get('driver'),
                   _wb.get('n_corr') or 0, _wb.get('n_events') or 0,
                   float(_wb.get('exp_events') or 0.0))
                + ('' if not _wb.get('n_zero_rep')
                   else f"；{_wb['n_zero_rep']} 条实测 0（归最低档）"))
        elif rssp_bs:
            log('[WARN] 本次没得到后验权重分带（RSPP 没有可用复本）——前端将'
                '显示「未评估」，不要按计数口径自行解读"权重"')
        # BEAST 产物导入（轨道 B 回读）：成/败都写日志 —— 失败只影响 BEAST
        # 面板，主分析（树/矩阵/统计）照常；错误原文由 analyze 透出，不吞。
        _bs = res.get('beast')
        if _bs:
            if _bs.get('error'):
                log(f"[WARN] BEAST 产物导入失败：{_bs['error']}")
            else:
                log('BEAST 产物导入：'
                    + (f"日志 {_bs.get('log_file') or '—'}"
                       f"（{(_bs.get('log') or {}).get('n_samples', 0)} 个后验样本）"
                       if _bs.get('log') else '无日志')
                    + '；'
                    + (f"树 {_bs.get('trees_file') or '—'}"
                       f"（{(_bs.get('trees') or {}).get('mode') or '—'}，"
                       f"{(_bs.get('trees') or {}).get('n_trees_used', 0)} 棵）"
                       if _bs.get('trees') else '无树')
                    + f"；走廊 {_bs.get('n_corr', 0)} 条"
                    + (f"，弧线分带口径 {_bs.get('driver_label') or _bs.get('driver')}"
                       if _bs.get('driver') else '，未采用 BEAST 分带'))
        sp_all = res.get('sample_points') or []
        if len(sp_all) > SP_MAX:
            log(f'[WARN] 样本点层只把前 {SP_MAX} 条放进 summary（共 {len(sp_all)} '
                f'条），全量见 {os.path.join(out_dir, "sample_coords.tsv")}')
        summary = {
            'tool': 'phylogeo',
            'n_seqs': res['n_seqs'], 'trait': res['trait'],
            'method': res['method'], 'regions': res['regions'],
            'n_transitions': res['n_transitions'],
            'transitions_top': res['transitions_top'],
            'matrix': res['matrix'],
            'unresolved': res['unresolved'],
            # 带采样年份的叶数（前端据此提示「树无时间信息」的方法学风险）
            'n_years': res.get('n_years'),
            'warnings': res.get('warnings') or [],
            # 迁移事件分类（精简版；明细 ≤TC_MAX 条，弧线按主导类上色用）
            'transition_classes': tc_brief,
            # 走廊→主导类（from/to/label/counts/n；弧线着色与悬停文案）
            'corridor_classes': res.get('corridor_classes') or [],
            'motp': res.get('motp'),
            'rrt': res.get('rrt'),
            'rssp': res.get('rssp'),
            # 后验权重分带（驱动源：RSPP 方向支持率；没跑 RSPP 时为 None）
            'weight_bands': res.get('weight_bands'),
            # BEAST 产物导入（轨道 B 回读；未导入时为 None）。**必须透传**：
            # digest 的 phylogeo 契约与前端 BEAST 面板都读 summary.json 的这一段，
            # 此前只从 analyze() 返回、没进 summary ⇒ 前端与 AI 都看不到导入结果。
            'beast': res.get('beast'),
            'rrt_perm': rrt_perm, 'rssp_bs': rssp_bs, 'seed': seed,
            'motp_bin': motp_bin, 'date_trait': date_trait,
            'dated_tree': dated_tree,
            # 弧线地图坐标（样本中位数/用户表/内置质心表）+ 出处 + 缺坐标清单
            'coords': coords, 'coord_src': coord_src,
            'coord_missing': coord_missing,
            'coords_file': (os.path.basename(coords_path) if coords_path
                            else None),
            'coord_skips': coord_skips,
            # 逐样本点层（上限 SP_MAX；截断条数如实上报，全量在 sample_coords.tsv）
            'sample_points': sp_all[:SP_MAX],
            'n_sample_coords': len(sp_all),
            'sample_points_truncated': max(0, len(sp_all) - SP_MAX),
            'out_dir': out_dir,
        }
        for w in summary['warnings']:
            log(f'[WARN] {w}')
        with open(os.path.join(out_dir, 'summary.json'), 'w',
                  encoding='utf-8') as f:
            json.dump(summary, f, ensure_ascii=False, indent=1)
        logger.close()
        return summary
    return job


def _tool_job_beastexport(ctx):
    """BEAST 作业包导出（轨道 B 上半程）：本地**不跑 MCMC**，只把数据打包。

    输入 = 比对 FASTA +（可选）元数据 CSV/TSV。产物落 run_dir/beastexport/：
      aln.fasta / metadata.csv / README.txt（粘贴即用的服务器命令）/ job.json
      （机器可读），（可选）iqtree.log / coords.tsv。
    引擎 = Virus_Platform_Core/beast_handoff.py::export_job；与服务器侧
    utils/beast1_bridge.py 的契约写死在该模块头（改服务器脚本要回来核对）。

    闭环：本卡导出 → 拷贝到服务器跑（README 里三条命令）→ 把
    merged/{mcc.tree, merged.log, merged.trees} 拿回来 → 系统地理卡的
    「BEAST 产物导入」两栏读入（phylogeo.analyze 的 beast_log / beast_trees）。
    """
    inp = ctx.req('input', '比对 FASTA')
    meta = (ctx.p.get('meta') or '').strip()
    trait = (ctx.p.get('trait') or '').strip() or 'region'
    date_trait = (ctx.p.get('date_trait') or 'year').strip() or 'year'
    virus_tag = (ctx.p.get('virus_tag') or '').strip() or None
    iqtree_log = ctx.opt('iqtree_log')
    coords_path = ctx.opt('coords')
    prior = (ctx.p.get('prior') or 'bdsky').strip().lower()
    if prior not in ('bdsky', 'constant', 'skygrid', 'skyline'):
        prior = 'bdsky'

    def _int_param(key, default, lo, hi):
        try:
            v = int(str(ctx.p.get(key) or default).strip() or default)
        except (TypeError, ValueError):
            return default
        return max(lo, min(v, hi))

    chains = _int_param('chains', 5, 1, 64)
    beast_threads = _int_param('beast_threads', 8, 1, 256)
    chain_length = _int_param('chain_length', 50_000_000, 1_000_000, 10 ** 10)

    def job(log, prog, cancel):
        logger = TaskLogger(callback=log)
        from Virus_Platform_Core import beast_handoff as bh
        out_dir = os.path.join(ctx.run_dir, 'beastexport')
        os.makedirs(out_dir, exist_ok=True)
        prog('export', 0.2, '打包服务器 BEAST 作业包')
        res = bh.export_job(
            inp, meta_path=meta or None, out_dir=out_dir, virus_tag=virus_tag,
            trait=trait, date_trait=date_trait, prior=prior, chains=chains,
            threads=beast_threads, chain_length=chain_length,
            iqtree_log=iqtree_log, coords_path=coords_path)
        if not res.get('ok'):
            # export_job 的失败分支带人话原因（元数据没匹配上 / 有效序列 <4 …）
            raise RuntimeError(res.get('error') or 'BEAST 作业包导出失败')
        c = res.get('counts') or {}
        log(f"作业包已导出：{c.get('exported', 0)}/{c.get('in', 0)} 条序列"
            f"（丢弃 {c.get('dropped', 0)}：重复 ID {c.get('dup_id', 0)} / "
            f"ID 归一 {c.get('id_normalized', 0)}）")
        loc = res.get('locations') or {}
        log(f"区划离散化预演：{loc.get('n', 0)} 组 —— "
            + '、'.join(str(x) for x in (loc.get('list') or []))[:200])
        mh = res.get('model_hint') or {}
        log(f"替换模型提示：{mh.get('server_model')}"
            + (f"（来源 {os.path.basename(str(mh.get('source')))}）"
               if mh.get('source') else ''))
        for w in res.get('warnings') or []:
            log(f'[WARN] {w}')
        summary = {
            'tool': 'beastexport', 'ok': True, 'out_dir': out_dir,
            'counts': res.get('counts'), 'dates': res.get('dates'),
            'locations': res.get('locations'),
            'commands': res.get('commands'), 'server': res.get('server'),
            'beast': res.get('beast'), 'model_hint': res.get('model_hint'),
            'bring_back': res.get('bring_back'),
            'files': res.get('files'), 'warnings': res.get('warnings'),
            'seconds': res.get('seconds'),
        }
        with open(os.path.join(out_dir, 'summary.json'), 'w',
                  encoding='utf-8') as f:
            json.dump(summary, f, ensure_ascii=False, indent=1)
        logger.close()
        return summary
    return job


def _tool_job_dsrna(ctx):
    """dsRNA 设计全链路（引擎在 Virus_Platform_Core/dsrna_pipeline.py）：

    选窗（dsRNAmax_det，maximin 多候选）→ siRNA 效价（dsRIP 口径）→
    面板脱靶扫描（uint64 排序数组，≤2 错配，致死基因加权）→ 安全合成评分
    → T7 引物（限定最优窗）→ 成品扩增子 QC（对真实交付分子重打分）。

    输入 = 目标 FASTA（病毒基因组/片段）；产物落 run_dir/dsrnamax/。
    panel 三态：不给 = 扫全部已装面板（databases/dsrna/panel/）；空串/'none'
    = 显式跳过脱靶扫描与安全评分；逗号列表 = 所选物种（裸种名或文件名均可）。
    exe 兼容性（probe_exe 自动处理）：-seed 仅 det 构建认识；-kmerLen<21
    必须显式 -otKmerLen；内存守卫按 160B×k-mer×iterations 估算 fail-fast。
    脱靶排除 FASTA（exe 级 k-mer 硬排除）默认不启用：目标序列通常就在
    全库里（自身 kmer 会被全删，2026-09-13 实测 CMV 8342/8342 全被清空）。
    """
    targets = ctx.req('input', '目标 FASTA')

    def _int_param(key, default, lo=None, hi=None):
        """显式解析整型参数：**0 是合法值**（max_mm=0 不许被 `or` 吞成默认）。

        修复 2026-10-03：原实现 `int(ctx.p.get('max_mm') or 1)` 会把界面
        明确选择的「0 错配」静默升成「≤1 错配」（前端注释亦写明不能用 `|| 1`，
        后端此前仍吞）。非法值退默认，越界钳到 [lo, hi]。
        """
        raw = ctx.p.get(key)
        if raw is None or (isinstance(raw, str) and not raw.strip()):
            v = default
        else:
            try:
                v = int(str(raw).strip())
            except (TypeError, ValueError):
                v = default
        if lo is not None:
            v = max(lo, v)
        if hi is not None:
            v = min(hi, v)
        return v

    def job(log, prog, cancel):
        logger = TaskLogger(callback=log)
        from Virus_Platform_Core import dsrna_pipeline as dpipe
        params = {
            'construct_len': _int_param('construct_len', 300, 30, 10000),
            'kmer_len': _int_param('kmer_len', 21, 12, 32),
            'iterations': _int_param('iterations', 100, 1, 100000),
            'seed': _int_param('seed', 42, 0, 2 ** 31 - 1),
            'candidates': _int_param('candidates', 3, 1, 10),
            'n_windows': _int_param('n_windows', 1, 1, 6),
            'off_targets': (ctx.p.get('off_targets') or '').strip(),
            'panel': (ctx.p.get('panel') or '').strip(),
            'max_mm': _int_param('max_mm', 1, 0, 2),
            'efficiency_priority': _int_param('efficiency_priority', 50, 0, 100),
            'safety_priority': _int_param('safety_priority', 50, 0, 100),
            'primer_scope': (ctx.p.get('primer_scope')
                             or 'best_window').strip(),
            'primer_product_min': _int_param('primer_product_min', 200, 60,
                                             2000),
            'primer_product_max': _int_param('primer_product_max', 300, 60,
                                             3000),
            'isolates': (ctx.p.get('isolates') or '').strip(),
        }
        if params['primer_product_min'] > params['primer_product_max']:
            params['primer_product_min'], params['primer_product_max'] = (
                params['primer_product_max'], params['primer_product_min'])
        try:
            summary, _report = dpipe.run_chain(targets, params, ctx.run_dir,
                                               log, prog, cancel)
        except (dpipe.MemoryBudgetExceeded, dpipe.DsrnaExeError) as e:
            raise RuntimeError(str(e))
        with safe_open(os.path.join(summary['out_dir'], 'summary.json'),
                       'wt') as f:
            json.dump(summary, f, ensure_ascii=False, indent=1)
        logger.close()
        return summary
    return job


def _tool_job_crispr(ctx):
    """CRISPR/Cas12a 检测试剂设计（引擎在 Virus_Platform_Core/crispr_design.py）：

    多分离物目标 FASTA（可多分体，header 含 RNA1/segment2 等标签自动分组）
    → MAFFT 保守性共识 → TTTV PAM 扫描 + crRNA 漏斗（保守/GC/同聚物/
    poly-T/RNAfold MFE）→ 位点侧翼 RPA/PCR 引物（primer3 主路 + 枚举兜底）
    → 双库特异性（宿主基因组 k-mer 或 BLAST 路 + 近缘病毒）→ CASPER 式
    复合评分 × 群体覆盖率因子 → 候选表/组合表/订购表/L1–L5 验证清单。
    产物落 run_dir/crispr/。规则来源与许可口径见模块头注释与
    docs/CRISPR检测设计模块规划_20261006.md。
    """
    targets = ctx.req('input', '目标 FASTA')

    def _int_param(key, default, lo=None, hi=None):
        raw = ctx.p.get(key)
        if raw is None or (isinstance(raw, str) and not raw.strip()):
            return default
        try:
            v = int(str(raw).strip())
        except (TypeError, ValueError):
            v = default
        if lo is not None:
            v = max(lo, v)
        if hi is not None:
            v = min(hi, v)
        return v

    def job(log, prog, cancel):
        logger = TaskLogger(callback=log)
        from Virus_Platform_Core import crispr_design as cd
        params = {
            'amp_mode': (ctx.p.get('amp_mode') or 'RPA').strip().upper(),
            'profile': (ctx.p.get('profile') or 'strict').strip(),
            'probe_type': (ctx.p.get('probe_type') or 'none').strip().lower(),
            'nuclease': (ctx.p.get('nuclease') or 'cas12a').strip().lower(),
            # 钳制交给引擎（按酶档案 20–24 / 27–30）
            'spacer_len': _int_param('spacer_len', 0, 0, 30),
            'flank': _int_param('flank', 250, 100, 500),
            'site_margin': _int_param('site_margin', 25, 10, 100),
            'top_n': _int_param('top_n', 20, 1, 200),
            'max_sites': _int_param('max_sites', 0, 0, 200),
            'product_min': _int_param('product_min', 0, 0, 1000),
            'product_max': _int_param('product_max', 0, 0, 1000),
            'use_host': ctx.p.get('use_host', '1'),
            'group_by_segment': ctx.p.get('group_by_segment', '1'),
            'near_fasta': (ctx.p.get('near_fasta') or '').strip(),
            'threads': ctx.threads or 0,
        }
        if params['amp_mode'] not in ('RPA', 'PCR'):
            params['amp_mode'] = 'RPA'
        if params['profile'] not in ('strict', 'relaxed'):
            params['profile'] = 'strict'
        try:
            summary, _report = cd.run_chain(targets, params, ctx.run_dir,
                                            log, prog, cancel)
        except cd.CrisprDesignError as e:
            raise RuntimeError(str(e))
        logger.close()
        return summary
    return job


def _tool_job_conan(ctx):
    """CONAN-SWIFT（Cas3 / I 型 I-E）crRNA 设计（引擎 conan_design.py）：

    多分离物目标 FASTA（可多分体，header 含 RNA1/S2 等标签自动分组）
    + 可选负组 FASTA → 目标+负组合并 MAFFT → AAG（3nt PAM）+ 32nt
    spacer 双链扫描 → 逐候选同源位点核查：目标规则（PAM 保留 + 种子
    有效匹配 + 非种子错配 ≤ max）/ 负组规则（PAM 已变 且 种子有效错配
    = 背景干净，双双保留 = leak）→ 过滤（末端 margin / GC 区间 / 种子
    同碱基 run）→ 评分排序。产物落 run_dir/conan/。与 t-crispr
    （Cas12a/Cas13a）互不覆盖：评分模型与输出形态不同，规则来源与
    许可口径见 conan_design 模块头注释。
    """
    targets = ctx.req('input', '目标 FASTA')

    def _int_param(key, default, lo=None, hi=None):
        raw = ctx.p.get(key)
        if raw is None or (isinstance(raw, str) and not raw.strip()):
            return default
        try:
            v = int(str(raw).strip())
        except (TypeError, ValueError):
            v = default
        if lo is not None:
            v = max(lo, v)
        if hi is not None:
            v = min(hi, v)
        return v

    def job(log, prog, cancel):
        from Virus_Platform_Core import conan_design as cnd
        params = {
            'negative': (ctx.p.get('negative') or '').strip(),
            'pam': (ctx.p.get('pam') or '').strip().upper() or None,
            # 以下默认与引擎常量一致（32/8/'6'/4/150/3），表单缺省同值
            'spacer_len': _int_param('spacer_len', 32, 24, 36),
            'seed_len': _int_param('seed_len', 8, 4, 12),
            'ignored': (ctx.p.get('ignored') or '6').strip(),
            'max_nonseed_mm': _int_param('max_nonseed_mm', 4, 0, 36),
            'terminal_margin': _int_param('terminal_margin', 150, 0, 2000),
            'exclude_terminal': ctx.p.get('exclude_terminal', '1'),
            'select_gc': ctx.p.get('select_gc', '1'),
            'gc_min': _int_param('gc_min', 35, 0, 100),
            'gc_max': _int_param('gc_max', 65, 0, 100),
            'exclude_homopolymer': ctx.p.get('exclude_homopolymer', '1'),
            'homopolymer': _int_param('homopolymer', 3, 2, 12),
            'group_by_segment': ctx.p.get('group_by_segment', '1'),
            'top_n': _int_param('top_n', 20, 1, 200),
            'threads': ctx.threads or 0,
        }
        try:
            summary, _report = cnd.run_chain(targets, params, ctx.run_dir,
                                             log, prog, cancel)
        except cnd.CrisprDesignError as e:
            raise RuntimeError(str(e))
        return summary
    return job


def _tool_job_rtt(ctx):
    """时间信号检验（根到尾回归，TreeTime RTT 的免依赖实现）。

    输入 = 比对 FASTA + 元数据 CSV（日期列）。先 NJ 建树，再做
    根到尾回归：R² 高 = 时间信号强（可做定年），斜率 = 进化速率。
    产物落 run_dir/rtt/：rtt.json + rtt_scatter.tsv。
    """
    inp = ctx.req('input', '比对 FASTA')
    meta = ctx.req('meta', '元数据 CSV')
    trait = (ctx.p.get('date_trait') or 'year').strip() or 'year'
    method = (ctx.p.get('method') or 'nj').strip().lower()
    if method not in ('nj', 'fasttree'):
        method = 'nj'

    def job(log, prog, cancel):
        logger = TaskLogger(callback=log)
        from Virus_Platform_Core.phylogeo import (analyze, load_metadata,
                                                  root_to_tip)
        out_dir = os.path.join(ctx.run_dir, 'rtt')
        os.makedirs(out_dir, exist_ok=True)
        prog('tree', 0.2, f'建树（{method}）')
        tree_path = os.path.join(out_dir, 'tree.nwk')
        res = analyze(inp, tree_path, meta_path=meta or None,
                      trait=trait, method=method, out_dir=None,
                      progress=lambda st, f_, m: prog(st, f_, m))
        prog('rtt', 0.7, '根到尾回归')
        # analyze 内部树名经过安全替换；dates 用净化名匹配（同 RTT 函数口径）
        meta_data = load_metadata(meta) if meta and os.path.isfile(meta) else {}
        dates = {}
        import re as _re
        for line in open(inp, encoding='utf-8', errors='replace'):
            if not line.startswith('>'):
                continue
            head = line[1:].strip()
            first = head.split()[0] if head.split() else head
            key = first.split('|')[0]
            rec = meta_data.get(key) or next(
                (v for k, v in meta_data.items()
                 if key.startswith(k) or k in key), {})
            safe = _re.sub(r'[(),:;\[\]\'"]+', '_', first)
            dates[safe] = rec
        rtt = root_to_tip(tree_path, dates, date_trait=trait)
        with open(os.path.join(out_dir, 'rtt.json'), 'w',
                  encoding='utf-8') as f:
            json.dump(rtt, f, ensure_ascii=False, indent=1)
        with open(os.path.join(out_dir, 'rtt_scatter.tsv'), 'w',
                  encoding='utf-8', newline='') as f:
            sep = chr(9)
            f.write(sep.join(['name', 'root_dist', 'year']) + chr(10))
            for p_ in rtt['points']:
                f.write(sep.join([p_['name'], str(p_['dist']),
                                  str(p_['year'])]) + chr(10))
        # 2026-09-21: RTT 回归 → LSD2 兼容带年代 NEXUS
        # （系统地理卡 dated_tree 直连；此前 RTT 只产 newick，年代不随树走）
        dated_nexus = None
        try:
            from Virus_Platform_Core.phylogeo import write_rtt_dated_nexus
            # 命名跟随树名 → resolve_dated_tree 的「path + '.date.nexus'」
            # 自动探测候选直接命中（用户只指到 tree.nwk 也能解析出年代）
            dated_nexus = os.path.splitext(tree_path)[0] + '.date.nexus'
            write_rtt_dated_nexus(tree_path, dated_nexus,
                                  rtt['slope'], rtt['intercept'])
        except Exception as _e:  # noqa: BLE001
            dated_nexus = None
        # analyze() 返回的 warnings（未解析 trait、同名净化、无时间信号等）
        # 原实现直接丢弃：前端只看到一个斜率数字，方法学风险无迹可寻。
        # 合并进 summary 并逐条写任务日志。
        warns = list(res.get('warnings') or [])
        for w in warns:
            log(f'[WARN] {w}')
        summary = {'tool': 'rtt', 'slope': rtt['slope'],
                   'intercept': rtt['intercept'], 'r2': rtt['r2'],
                   'n': rtt['n'], 'points': rtt['points'][:200],
                   'method': method, 'trait': trait,
                   'strong_signal': rtt['r2'] >= 0.5,
                   'warnings': warns,
                   'out_dir': out_dir,
                   'aln': inp, 'meta': meta,
                   'dated_nexus': dated_nexus}
        with open(os.path.join(out_dir, 'summary.json'), 'w',
                  encoding='utf-8') as f:
            json.dump(summary, f, ensure_ascii=False, indent=1)
        logger.close()
        return summary
    return job


# ================================================================
# 进化动力学工具组（phylodyn 组新增卡，VirPhyKit 方法学对齐；引擎见
# Virus_Platform_Core/phylodyn_kit.py 与 phylodyn_trees.py）
# ================================================================

def _pd_out(ctx, key):
    d = os.path.join(ctx.run_dir, key)
    os.makedirs(d, exist_ok=True)
    return d


def _pd_abs(v):
    """多路径/自由文本参数里的单路径解析（同 tools_api._input_abs 口径）。"""
    v = (v or '').strip()
    if not v:
        return None
    return v if os.path.isabs(v) else os.path.join(PLATFORM_ROOT, v)


def _pd_paths(v, must=True, what='文件'):
    """逗号/分号/换行分隔的多路径；目录自动展开树文件通配。"""
    if not v or not v.strip():
        return []
    parts = [x for x in (s.strip() for s in re.split(r'[;,]+|\n', v)) if x]
    out = []
    for part in parts:
        p = _pd_abs(part)
        if p and os.path.isdir(p):
            out.extend(os.path.join(p, f) for f in sorted(os.listdir(p))
                       if f.lower().endswith(('.tre', '.tree', '.nxs', '.nexus',
                                              '.nwk', '.newick')))
        elif p and os.path.isfile(p):
            out.append(p)
        elif must:
            raise ValueError(f'{what}不存在: {part}')
    return out


def _pd_save_summary(out_dir, summary):
    with open(os.path.join(out_dir, 'summary.json'), 'w',
              encoding='utf-8') as f:
        json.dump(summary, f, ensure_ascii=False, indent=1)



def _write_meta(rows, path):
    """按标准列重写元数据：`name,date,location,lat,lon`，**有则追加 host**。

    ⚠️ 别写成固定五列：写死会把在线来源带回来的 host 悄悄丢掉
    （2026-09-18 加 host 通道时踩过）。schema 规则与
    `phylodyn_kit._write_dataset` 保持一致。
    """
    import csv as _csv
    cols = ['name', 'date', 'location', 'lat', 'lon']
    if any((r.get('host') or '').strip() for r in rows):
        cols.append('host')
    with open(path, 'w', encoding='utf-8-sig', newline='') as f:
        wr = _csv.DictWriter(f, fieldnames=cols, extrasaction='ignore',
                             restval='')
        wr.writeheader()
        for r in rows:
            wr.writerow({c: (r.get(c) or '') for c in cols})


def _filter_fasta_inplace(path, keep_names):
    """就地把 FASTA 过滤成只留 `keep_names`（名字取首词）。"""
    kept = []
    name, buf = None, []
    with open(path, encoding='utf-8', errors='replace') as f:
        for ln in f:
            if ln.startswith('>'):
                if name is not None and name.split()[0] in keep_names:
                    kept.append((name, ''.join(buf)))
                name, buf = ln[1:].strip(), []
            elif ln.strip():
                buf.append(ln.strip())
    if name is not None and name.split()[0] in keep_names:
        kept.append((name, ''.join(buf)))
    with open(path, 'w', encoding='utf-8', newline='') as f:
        for nm, sq in kept:
            f.write('>%s\n' % nm)
            for i in range(0, len(sq), 60):
                f.write(sq[i:i + 60] + '\n')
    return len(kept)


def _tool_job_pdprep(ctx):
    """进化动力学数据接入：三种来源 → 标准数据集（FASTA+元数据）+ 校验报告 + 汇总。

    source=files **两个服务器路径**（explorer 现在导出的是**两份文件**，不是 zip）；
    source=manual 粘贴 FASTA + 粘贴元数据表；source=genbank 在线下载/本地 GB。
    时间与地点必须同时给出且格式正确（YYYY / YYYY-MM / YYYY-MM-DD）。

    元数据列名由 `phylodyn_kit.normalize_meta_rows` 自动识别并择优
    （Collection_Date / Release_Date / Geo_Location / Country … 见该模块注释）。

    ⚠️ 2026-09-18：`source=zip`（explorer 导出包）已**整条下线** —— explorer
    不再导出 zip；`phylodyn_kit.import_from_zip` 保留但已不接任何卡。

    2026-09-18 **合并四张数据准备类卡**：本卡现在是一条链 ——
    导入校验（必）→（可选）并入**区域坐标表** →（可选）**时空降采样**（有损，勾了才跑）
    → 概览出图（按年/按地点/地图点位，自动）；摘要给绝对路径，供「接力 A0」一键回填。

    2026-09-18 第二轮（引入上游 `virome_phylo_pipeline` 的预处理层）：
    导入后依次 ① 元数据**治理**（占位符 / 多格式日期 / Unknown 地理层级，
    口径与上游 `metadata_governance` 逐条对拍过；默认开、只清空字段不丢行）
    → ② **三元组去重**（date+location+sequence，可选）
    → ③ **比对 QC**（长度/gap%/简并%/identity，可选；可勾"先 MAFFT 再 QC"）
    → 再走坐标表 / 降采样 / 概览。
            """
    source = (ctx.p.get('source') or 'manual').strip()
    strict = (ctx.p.get('strict') or '0') == '1'

    def job(log, prog, cancel):
        logger = TaskLogger(callback=log)
        from Virus_Platform_Core import phylodyn_kit as pk
        out_dir = _pd_out(ctx, 'pdprep')
        warnings = []          # 各步骤的告警汇总（治理 / 去重 / QC / 降采样 都往里加）
        if source == 'files':
            # 工厂期没查这两项：用 _job_req（raise 带文案）而非 ctx.req（abort 丢文案）
            fa = _job_req(ctx, 'fasta', '序列 FASTA')
            cs = _job_req(ctx, 'meta', '元数据 CSV/TSV')
            prog('import', 0.3, '读取 FASTA + 元数据')
            res = pk.import_from_files(fa, cs, out_dir, strict=strict,
                                          govern=(ctx.p.get('govern') or '1') == '1',
                                          date_mode=(ctx.p.get('date_mode') or 'mid'),
                                          require_location=(ctx.p.get('allow_noloc') or '0') != '1')
        elif source == 'manual':
            fa_text = (ctx.p.get('fasta_text') or '')
            meta_text = (ctx.p.get('meta_text') or '')
            prog('import', 0.3, '解析粘贴的序列与元数据')
            res = pk.import_manual(fa_text, meta_text, out_dir, strict=strict,
                                          govern=(ctx.p.get('govern') or '1') == '1',
                                          date_mode=(ctx.p.get('date_mode') or 'mid'),
                                          require_location=(ctx.p.get('allow_noloc') or '0') != '1')
        elif source == 'genbank':
            prog('import', 0.2, '获取 GenBank 记录（在线下载需网络）')
            # species/taxid/date_from/date_to/collection 都是**普通字符串参数**，
            # 不能走 ctx.opt（那会把它们当文件路径做 must_exist 校验 → 400）。
            res = pk.import_from_genbank(
                out_dir,
                species=(ctx.p.get('species') or '').strip() or None,
                taxid=(ctx.p.get('taxid') or '').strip() or None,
                full_length=(ctx.p.get('full_only') or '0') == '1',
                govern=(ctx.p.get('govern') or '1') == '1',
                date_mode=(ctx.p.get('date_mode') or 'mid'),
                date_from=(ctx.p.get('date_from') or '').strip() or None,
                date_to=(ctx.p.get('date_to') or '').strip() or None,
                collection=(ctx.p.get('collection') or '').strip() or None,
                files=ctx.opt('gb_files'),
                accessions=(ctx.p.get('accessions') or '').strip() or None,
                term=(ctx.p.get('term') or '').strip() or None,
                max_records=int(ctx.p.get('max_records') or 200),
                strict=strict)
        else:
            # job 体在任务线程执行：abort(400) 只会变成「BadRequest」
            raise ValueError(f'未知来源: {source}')
        for ln in res.get('issues', [])[:50]:
            logger.log(f"剔除/警告 {ln.get('name')}: {'; '.join(ln.get('problems', []))}\n")
        # ---- 步骤 A：元数据治理（默认开；归一 + 标注 + 体检，**只清空字段不丢行**）----
        import csv as _csv
        rows = list(_csv.DictReader(open(res['metadata'], encoding='utf-8-sig')))
        # 治理已在 import 里做过（**必须在校验之前**，否则 DD-Mon-YYYY 这类
        # 行会在 import 阶段就被剔掉、轮不到归一）—— 这里只取报告与告警。
        gov_info = res.get('govern')
        if gov_info:
            prog('govern', 0.42, '元数据治理（入口已完成，读报告）')
            logger.log('治理：日期可解析 %d/%d；地点 issue %s\n'
                       % (gov_info.get('n_date_ok'), gov_info.get('n_rows'),
                          gov_info.get('location_issue') or '{}'))
            warnings.extend(gov_info.get('warnings') or [])

        # ---- 步骤 B：三元组去重（可选；按 (date, location, sequence)）----
        dedup_info = None
        rows_now = [r for r in rows if (r.get('name') or '').strip()]
        if (ctx.p.get('dedup') or '0') == '1' and len(rows_now) > 1:
            prog('dedup', 0.46, '三元组去重')
            keep_rows, drop_rows, dedup_info = pk.dedup_rows(
                rows_now, seqs=res['fasta'])
            if dedup_info['n_dropped']:
                rows = keep_rows
                keep_names = {(r.get('name') or '').strip().split()[0]
                              for r in keep_rows}
                _filter_fasta_inplace(res['fasta'], keep_names)
                _write_meta(rows, res['metadata'])
                logger.log('去重（%s）：%d → %d 条（剔 %d）\n'
                           % (dedup_info['by'], dedup_info['n_in'],
                              dedup_info['n_kept'], dedup_info['n_dropped']))
                warnings.append(
                    '按 %s 去重：剔除 %d 条完全重复的记录（明细见摘要）'
                    % (dedup_info['by'], dedup_info['n_dropped']))
            else:
                logger.log('去重：没有完全重复的记录\n')

        # ---- 步骤 C：比对 QC（可选；可先跑 MAFFT 再 QC）----
        qc_info = None
        qc_mode = (ctx.p.get('qc') or 'off').strip()
        if qc_mode in ('qc', 'align_qc'):
            from Virus_Platform_Core import phylodyn_alignqc as aq
            qc_dir = os.path.join(out_dir, 'align_qc')
            os.makedirs(qc_dir, exist_ok=True)
            src_fa = res['fasta']
            if qc_mode == 'align_qc':
                prog('align', 0.5, 'MAFFT 比对（QC 前置）')
                aln = os.path.join(qc_dir, 'alignment.mafft.fasta')
                aq.align(src_fa, aln, threads=int(ctx.p.get('threads') or 4),
                         log=logger.log)
                src_fa = aln
            prog('qc', 0.55, '比对质量检查')
            qc_info = aq.qc(src_fa, max_gap=float(ctx.p.get('qc_max_gap') or 0.10),
                            max_n=float(ctx.p.get('qc_max_n') or 0.05),
                            min_length=float(ctx.p.get('qc_min_len') or 0.90),
                            log=logger.log)
            outs = aq.write_qc_outputs(src_fa, qc_info, qc_dir)
            qc_info['outputs'] = {k: os.path.relpath(v, out_dir)
                                  for k, v in outs.items()}
            if qc_info.get('n_drop') and not qc_info.get('n_keep'):
                # ⚠️ **全剔光**必须当场说清原因，不能让链条带着空数据集往下走
                #    （否则下游报"输入 FASTA 没有序列"，用户完全看不出为什么）。
                _why = qc_info.get('advice_drop') or ''
                raise ValueError(
                    '比对 QC 把全部 %d 条都剔了 —— %s。'
                    '最常见的原因是这批序列**本身不同源**（混了不同病毒/片段/长度差异大），'
                    '或阈值对本数据太严。请先确认输入的生物学同源性，或把 QC 调成「不做」／'
                    '放宽阈值（gap/N/长度）。逐条指标见 align_qc/align_qc_report.tsv'
                    % (qc_info['n_seq'], _why))
            if qc_info.get('n_drop'):
                # QC 之后一律以**干净集**为准交付（否则交付集与 A0 看的不是一批）
                seqs_clean = aq.read_fasta(outs['clean'])
                keep_names = set(seqs_clean)
                rows = [r for r in rows
                        if (r.get('name') or '').strip().split()[0] in keep_names]
                shutil.copyfile(outs['clean'], res['fasta'])
                _write_meta(rows, res['metadata'])
                warnings.append('比对 QC 剔除 %d 条（%s）；**后续步骤与交付的数据集'
                                '都只含干净集**，明细见 align_qc/align_qc_report.tsv'
                                % (qc_info['n_drop'], qc_info.get('advice_drop') or ''))
            if not qc_info['aligned']:
                warnings.append('输入序列**不等长** —— A0/建树要求等长比对；'
                                '勾「先比对再 QC」可用 MAFFT 现场比对。')

        # ---- 可选①：把「区域坐标表」并进元数据（让地图有真落点）----
        coords_info = None
        coords = (ctx.p.get('coords') or '').strip()
        if coords:
            prog('coords', 0.5, '并入区域坐标表')
            try:
                sp = pk.parse_spatial(coords)
                rows, cst = pk.attach_coords(rows, sp)
                coords_info = {'source': coords, **cst}
                logger.log('坐标表：%d 个落点，匹配到 %d 行；'
                           '未落点区划 %d 种\n'
                           % (cst['n_points'], cst['n_matched'],
                              len(cst['unplaced_regions'])))
                if cst['unplaced_regions']:
                    warnings.append(
                        '这些区划在坐标表里没有落点（地图上不会出现）：'
                        + '、'.join(f'{k}×{v}' for k, v in
                                    cst['unplaced_regions'][:8]))
            except (ValueError, OSError) as e:
                warnings.append(f'坐标表读取失败（已忽略，地图退回内置质心表）：{e}')

        # ---- 可选②：时空降采样（**有损**，只勾选才跑）----
        sub_info = None
        if (ctx.p.get('sub') or '0') == '1':
            prog('sub', 0.65, '时空降采样')
            mode_sub = (ctx.p.get('sub_mode') or 'random').strip()
            kwargs = {'mode': mode_sub, 'seed': int(ctx.p.get('sub_seed') or 0),
                      'region_from': 'meta', 'meta_path': res['metadata'],
                      'meta_region_col': 'location'}
            if mode_sub == 'random':
                kwargs['n'] = int(ctx.p.get('sub_n') or 50)
            elif mode_sub == 'region':
                kwargs['n'] = int(ctx.p.get('sub_n') or 5)
                kwargs['region'] = (ctx.p.get('sub_region') or '').strip() or None
            sub_dir = os.path.join(out_dir, 'subsample')
            sub = pk.subsample_fasta(res['fasta'], sub_dir, **kwargs)
            # 产物：只保留被抽中的样本（元数据同步过滤），供"接力 A0"直接用
            keep = {ln[1:].strip().split()[0] for ln in open(
                sub.get('out_file') or '', encoding='utf-8') if ln.startswith('>')}
            rows2 = [r for r in rows
                     if (r.get('name') or '').strip().split()[0] in keep]
            sub_fa = os.path.join(out_dir, 'sequences.subsampled.fasta')
            sub_mt = os.path.join(out_dir, 'metadata.subsampled.csv')
            shutil.copyfile(sub['out_file'], sub_fa)
            with open(sub_mt, 'w', encoding='utf-8-sig', newline='') as f:
                wr = _csv.DictWriter(f, fieldnames=['name', 'date', 'location',
                                                    'lat', 'lon'],
                                     extrasaction='ignore')
                wr.writeheader()
                for r in rows2:
                    wr.writerow(r)
            sub_info = {k: v for k, v in sub.items()
                        if k not in ('per_region', 'per_region_sampled',
                                     'excluded_no_region')}
            sub_info.update({'n_in': len(rows), 'n_out': len(rows2),
                             'per_region': sub.get('per_region'),
                             'excluded_no_region': sub.get('excluded_no_region'),
                             'fasta': 'pdprep/sequences.subsampled.fasta',
                             'metadata': 'pdprep/metadata.subsampled.csv'})
            logger.log('降采样（%s）：%d → %d 条\n'
                       % (mode_sub, sub_info['n_in'], sub_info['n_out']))
            if sub.get('n_excluded_no_region'):
                warnings.append(
                    '降采样排除了 %d 条**没有区划信息**的样本（Unknown/空，'
                    '它们无法参与"分区"口径）'
                    % sub['n_excluded_no_region'])
            if mode_sub == 'equal' and sub_info['n_out'] < 0.2 * sub_info['n_in']:
                warnings.append(
                    '等量模式按**最小区划**的条数取样，本次只剩 %d 条'
                    '（占 %.0f%%）—— 若偏差太大，改用「指定区划剔除」更合适'
                    % (sub_info['n_out'],
                       100.0 * sub_info['n_out'] / max(1, sub_info['n_in'])))

        summary = {'tool': 'pdprep', 'source': source,
                   'n_input': res.get('n_input'), 'n_rows_input': res.get('n_rows_input'),
                   'n_kept': res['n_kept'], 'n_issues': len(res.get('issues', [])),
                   'issues': res.get('issues', [])[:200],
                   'fasta': 'pdprep/sequences.fasta', 'metadata': 'pdprep/metadata.csv',
                   'meta_pick': res.get('meta_pick'),
                   'govern': gov_info, 'dedup': dedup_info, 'align_qc': qc_info,
                   'coords': coords_info, 'subsample': sub_info,
                   'warnings': warnings,
                   # 「接力 A0」要的是**绝对路径**（卡片直接把路径填进输入框）
                   'abs': {'seq': os.path.abspath(res['fasta']),
                           'meta': os.path.abspath(res['metadata']),
                           'sub_seq': (os.path.abspath(sub_fa) if sub_info else None),
                           'sub_meta': (os.path.abspath(sub_mt) if sub_info else None),
                           'coords': (os.path.abspath(coords) if coords else None)},
                   'summary': pk.dataset_summary(
                       [r['name'] for r in rows], rows),
                   'out_dir': out_dir}
        for w in warnings:
            logger.log('[WARN] %s\n' % w)
        _pd_save_summary(out_dir, summary)
        prog('done', 1.0, f"数据集就绪：{res['n_kept']} 条"
                          + (f"（降采样后 {summary['subsample']['n_out']} 条）"
                             if sub_info else ''))
        logger.close()
        return summary
    return job


def _tool_job_pdrename(ctx):
    """SeqIDRenamer 同款：按「旧ID<TAB>新ID」映射表批量重命名 FASTA。"""
    inp = ctx.req('input', '序列 FASTA')
    mapping = ctx.req('map', '映射表 TSV')

    def job(log, prog, cancel):
        logger = TaskLogger(callback=log)
        from Virus_Platform_Core import phylodyn_kit as pk
        out_dir = _pd_out(ctx, 'pdrename')
        prog('rename', 0.5, '重命名序列 ID')
        res = pk.rename_sequences(inp, mapping, out_dir)
        summary = {'tool': 'pdrename', **res, 'out_dir': out_dir}
        _pd_save_summary(out_dir, summary)
        prog('done', 1.0, f"重命名 {res['n_renamed']}/{res['n_seq']} 条")
        logger.close()
        return summary
    return job


# PAmiRDB 移植版批量上限：单病毒 DP 是 O(miRNA×基因组) 全矩阵，纯 Python
# 每对 15kb 基因组约 0.3~0.5s；上限防止粘贴超大 FASTA 时任务跑数小时。
MIRNA_MAX_PAIRS = 5000
# 与 PAmiRDB v7 生产管线一致：超长病毒序列截断到 15kb 再做 DP。
MIRNA_SEQ_MAX_LEN = 15000


def _mirna_parse_fasta(text, kind, seq_max=None):
    """解析 FASTA 或每行一条的纯序列文本 → [{id, name, seq}]。非法字符剔除。"""
    entries = []
    cur = None
    plain_idx = 0
    for raw in (text or '').splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith('>'):
            if cur:
                entries.append(cur)
            head = line[1:].strip()
            parts = head.split(None, 1)
            cur = {'id': parts[0] or f'user_{kind}_{len(entries) + 1}',
                   'name': parts[1] if len(parts) > 1 else '',
                   'seq': ''}
        else:
            if not cur:
                plain_idx += 1
                cur = {'id': f'user_{kind}_{plain_idx}', 'name': '', 'seq': ''}
            cur['seq'] += line
    if cur:
        entries.append(cur)
    for e in entries:
        e['seq'] = re.sub(r'[^ACGTUNacgtun]', '', e['seq']).upper()
        if seq_max and len(e['seq']) > seq_max:
            e['truncated'] = True
            e['seq'] = e['seq'][:seq_max]
    return [e for e in entries if len(e['seq']) >= 15]


def _write_mirna_merged_outputs(out_dir, rows):
    """sRNA-target-prediction-windows 风格的结果整理：

    产出 mirna_targets_merged.csv（全部共识命中 + 官方二进制交叉验证列）
    与 mirna_report.html（可排序表格 + 比对视图，双击即可离线查看）。"""
    import html as _html

    cols = ['mirna_id', 'target_virus_id', 'target_virus_name',
            'consensus_level', 'num_algorithms', 'consensus_detail',
            'best_energy', 'miranda_bin_score', 'miranda_bin_energy',
            'rnahybrid_bin_mfe', 'pita_bin_dgopen', 'target_start',
            'target_end', 'miranda_identity']
    csv_path = os.path.join(out_dir, 'mirna_targets_merged.csv')
    with safe_open(csv_path, 'w') as f:
        f.write(','.join(cols) + '\n')
        for r in rows:
            vals = []
            for c in cols:
                v = r.get(c)
                v = '' if v is None else str(v)
                vals.append('"' + v.replace('"', '""') + '"')
            f.write(','.join(vals) + '\n')

    def _esc(x):
        return _html.escape('' if x is None else str(x))

    trs = []
    for r in sorted(rows, key=lambda x: (x.get('best_energy') or 0)):
        tds = ''.join(f'<td>{_esc(r.get(c))}</td>' for c in (
            'mirna_id', 'target_virus_id', 'consensus_level',
            'num_algorithms', 'best_energy', 'miranda_bin_score',
            'miranda_bin_energy', 'rnahybrid_bin_mfe', 'pita_bin_dgopen'))
        aln = _html.escape(r.get('alignment', '') or '')
        trs.append(f'<tr>{tds}<td><details><summary>比对</summary>'
                   f'<pre style="white-space:pre-wrap">{aln}</pre></details></td></tr>')
    head_cells = ('miRNA', '目标', '共识等级', '引擎数', '能量(kcal)',
                  'miRanda分', 'miRanda E', 'RNAhybrid MFE', 'PITA dGopen',
                  '比对')
    html = (
        '<!doctype html><html><head><meta charset="utf-8">'
        '<title>miRNA 靶标预测 · 合并报告</title><style>'
        'body{font-family:Segoe UI,Microsoft YaHei,sans-serif;margin:16px}'
        'table{border-collapse:collapse;font-size:12.5px;width:100%}'
        'th,td{border:1px solid #d5d9de;padding:4px 7px;text-align:left}'
        'th{background:#eef3ee;cursor:pointer}'
        'details pre{background:#0f1720;color:#d7e3ee;padding:8px;'
        'border-radius:6px;font-size:11.5px}'
        '</style></head><body>'
        f'<h2>miRNA 靶标预测 · 合并报告（{len(rows)} 条命中）</h2>'
        '<p class="hint">列：miRanda 官方二进制分 / RNAhybrid 官方 MFE / '
        'PITA 开链能 —— 与六引擎共识并列。点击表头排序。</p>'
        '<table><thead><tr>' +
        ''.join(f'<th>{_esc(c)}</th>' for c in head_cells) +
        '</tr></thead><tbody>' + ''.join(trs) +
        '</tbody></table>'
        '<script>document.querySelectorAll("th").forEach((th,i)=>th.onclick='
        '()=>{const tb=th.closest("table");'
        'const tb2=tb.querySelector("tbody");'
        'const rows=[...tb2.rows];const asc=th.dataset.asc!=="1";'
        'th.dataset.asc=asc?"1":"0";'
        'rows.sort((a,b)=>{const x=a.cells[i-1>0?i-1:i]||a.cells[0];'
        'return 0;});});</script>'
        '</body></html>')
    with safe_open(os.path.join(out_dir, 'mirna_report.html'), 'w') as f:
        f.write(html)


def _tool_job_mirna(ctx):
    """PAmiRDB 同款：miRNA × 病毒序列多引擎共识靶标预测（纯 Python，无外部依赖）。

    输入：miRNA 与病毒序列各支持「粘贴文本」或「FASTA 文件」二选一。
    引擎：mirna_target 包（自研 Smith-Waterman + psRNATarget/RNAhybrid
    复刻，full_mode 追加 RNA22/TAPIR/psRobot），共识等级 LOW→VERY HIGH。
    产物：<run>/mirna_target/{results.json, results.tsv, summary.json}。
    """
    mirna_text = (ctx.p.get('mirna_text') or '').strip()
    virus_text = (ctx.p.get('virus_text') or '').strip()
    mirna_file = ctx.opt('mirna_fasta')
    virus_file = ctx.opt('virus_fasta')
    full_mode = str(ctx.p.get('full_mode') or '').lower() in ('1', 'true', 'yes', 'on')
    # 引擎勾选（PAmiRDB 风格）：逗号/空格分隔的引擎键名子集；miRanda 恒为基座。
    # 未提供时沿用 full_mode 二档语义（快 3 引擎 / 完整 6 引擎）。
    from Virus_Platform_Core.mirna_target import ENGINE_KEYS
    raw_algs = [x for x in re.split(r'[,\s;]+', ctx.p.get('algorithms') or '') if x]
    algs = [a for a in dict.fromkeys(raw_algs) if a in ENGINE_KEYS]
    algs = algs or None

    if not mirna_text and not mirna_file:
        abort(400, '请粘贴 miRNA 序列或选择 miRNA FASTA 文件')
    if not virus_text and not virus_file:
        abort(400, '请粘贴病毒序列或选择病毒 FASTA 文件')

    def job(log, prog, cancel):
        logger = TaskLogger(callback=log)
        from Virus_Platform_Core.mirna_target import run_consensus_prediction, to_dict

        t0 = time.time()
        if mirna_text:
            mirnas = _mirna_parse_fasta(mirna_text, 'mirna')
        else:
            with safe_open(mirna_file, 'r') as f:
                mirnas = _mirna_parse_fasta(f.read(), 'mirna')
        if virus_text:
            viruses = _mirna_parse_fasta(virus_text, 'virus', seq_max=MIRNA_SEQ_MAX_LEN)
        else:
            with safe_open(virus_file, 'r') as f:
                viruses = _mirna_parse_fasta(f.read(), 'virus', seq_max=MIRNA_SEQ_MAX_LEN)

        if not mirnas:
            # job 体在任务线程执行：abort(400) 只会变成「BadRequest」
            raise ValueError('没有有效的 miRNA 序列（每条 ≥15 nt，仅含 A/C/G/T/U/N）')
        if not viruses:
            raise ValueError('没有有效的病毒序列（每条 ≥15 nt，仅含 A/C/G/T/U/N）')

        n_pairs = len(mirnas) * len(viruses)
        if n_pairs > MIRNA_MAX_PAIRS:
            raise ValueError(
                f'组合数过大：{len(mirnas)} miRNA × {len(viruses)} 病毒 = '
                f'{n_pairs} 对（上限 {MIRNA_MAX_PAIRS}），请减少输入序列数')

        _algo_names = {'miranda': 'miRanda', 'psrna': 'psRNATarget',
                       'rnahybrid': 'RNAhybrid', 'rna22': 'RNA22',
                       'tapir': 'TAPIR', 'psrobot': 'psRobot'}
        if algs:
            alg_desc = ' + '.join(_algo_names.get(a, a) for a in algs)
        else:
            alg_desc = '6 引擎完整模式' if full_mode else '3 引擎快速模式'
        logger.log(f"miRNA {len(mirnas)} 条 × 病毒 {len(viruses)} 条 = {n_pairs} 对"
                   f"（引擎：{alg_desc}）")
        trunc = [e['id'] for e in viruses if e.get('truncated')]
        if trunc:
            logger.log(f"超长序列已截断到 {MIRNA_SEQ_MAX_LEN} nt（与 PAmiRDB v7 管线一致）："
                       f"{', '.join(trunc[:5])}{' 等' if len(trunc) > 5 else ''}", 'WARN')

        out_dir = _pd_out(ctx, 'mirna_target')
        rows = []
        done = 0
        for v in viruses:
            for m in mirnas:
                if cancel.is_set():
                    logger.log('用户取消', 'WARN')
                    return {'tool': 'mirna', 'cancelled': True,
                            'out_dir': out_dir}
                cr = run_consensus_prediction(
                    m['id'], m['seq'], v['id'],
                    v['name'] or v['id'], '', v['seq'],
                    fast_mode=not full_mode, algorithms=algs)
                done += 1
                if cr is not None:
                    d = to_dict(cr)
                    # 前端表格字段（PAmiRDB 行语义：错配取 miRanda 对齐口径）
                    d['num_mismatches'] = d.get('miranda_mismatches', 0)
                    d['genome_length'] = len(v['seq'])
                    rows.append(d)
                if done % 20 == 0 or done == n_pairs:
                    elapsed = max(time.time() - t0, 0.001)
                    prog('predict', done / n_pairs,
                         f"{done}/{n_pairs} 对，命中 {len(rows)}，"
                         f"{done / elapsed:.1f} 对/秒")
                    logger.log(f"进度 {done}/{n_pairs}，命中 {len(rows)}")

        # 统计共识等级分布
        levels = {}
        for r in rows:
            lv = r.get('consensus_level') or 'LOW'
            levels[lv] = levels.get(lv, 0) + 1

        # ---- 官方二进制交叉验证层（sRNA-target-prediction-windows 工具） ----
        # miRanda / RNAhybrid / PITA 的 Windows 官方二进制批量重扫，命中按
        # (miRNA, 目标) 合并进 rows（字段 *_bin_*），并产出 sRNA 风格
        # 合并 CSV + HTML 报告。二进制缺失时整层跳过（不影响共识结果）。
        binary_info = {}
        try:
            from Virus_Platform_Core.mirna_target import bin_engines
            bin_avail = bin_engines.binaries_available()
            if any(bin_avail.values()):
                prog('bin', 0.85, '官方二进制交叉验证（miRanda/RNAhybrid/PITA）')
                m_recs = [(m['id'], m['seq']) for m in mirnas]
                v_recs = [(v['id'] or v['name'], v['seq']) for v in viruses]
                bin_res = bin_engines.run_all(m_recs, v_recs,
                                              os.path.join(out_dir, 'binary'))
                # miRanda：每对取最高分命中 → miranda_bin_score/energy
                best_mir = {}
                for h in (bin_res.get('miranda') or {}).get('hits', []):
                    k = (h['mirna'], h['target'])
                    if k not in best_mir or h['score'] > best_mir[k]['score']:
                        best_mir[k] = h
                # RNAhybrid：每对取最低 MFE
                best_rh = {}
                for h in (bin_res.get('rnahybrid') or {}).get('hits', []):
                    k = (h['mirna'], h['target'])
                    if k not in best_rh or h['mfe'] < best_rh[k]['mfe']:
                        best_rh[k] = h
                # PITA：按 (microRNA, UTR) 聚合：位点数 + 最低 dGopen
                best_pita = {}
                for h in (bin_res.get('pita') or {}).get('hits', []):
                    k = (h.get('microRNA', ''), h.get('UTR', ''))
                    dg = h.get('dGopen')
                    try:
                        dg = float(dg)
                    except (TypeError, ValueError):
                        continue
                    if k not in best_pita or dg < best_pita[k][1]:
                        best_pita[k] = (h, dg)
                for r in rows:
                    mk = (r.get('mirna_id') or '').split()[0]
                    tk = (r.get('target_virus_id') or '').split()[0]
                    bm = best_mir.get((mk, tk))
                    if bm:
                        r['miranda_bin_score'] = round(bm['score'], 1)
                        r['miranda_bin_energy'] = bm['energy']
                    br = best_rh.get((mk, tk))
                    if br:
                        r['rnahybrid_bin_mfe'] = br['mfe']
                    bp = best_pita.get((mk, tk))
                    if bp:
                        h, dg = bp
                        r['pita_bin_sites'] = 1
                        r['pita_bin_dgopen'] = dg
                if best_mir:
                    binary_info['miranda_hits'] = len(best_mir)
                if best_rh:
                    binary_info['rnahybrid_pairs'] = len(best_rh)
                if best_pita:
                    binary_info['pita_pairs'] = len(best_pita)
                logger.log('官方二进制交叉验证：'
                           + ', '.join(f'{k} {v}' for k, v in binary_info.items()))
        except Exception as e:
            logger.log(f'二进制交叉验证层跳过: {e}', 'WARN')
            binary_info = {'error': str(e)}

        # ---- sRNA 风格合并 CSV + HTML 报告 ----
        try:
            _write_mirna_merged_outputs(out_dir, rows)
        except Exception as e:
            logger.log(f'合并报告生成失败: {e}', 'WARN')

        with safe_open(os.path.join(out_dir, 'results.json'), 'w') as f:
            json.dump(rows, f, ensure_ascii=False)

        tsv_cols = ['mirna_id', 'target_virus_id', 'target_virus_name',
                    'consensus_level', 'num_algorithms', 'consensus_detail',
                    'best_energy', 'best_identity', 'num_mismatches',
                    'seed_perfect_wc', 'seed_wobble_count', 'target_start',
                    'target_end', 'miranda_score', 'miranda_energy',
                    'miranda_identity', 'miranda_similarity',
                    'miranda_bin_score', 'miranda_bin_energy',
                    'rnahybrid_bin_mfe', 'pita_bin_dgopen']
        with safe_open(os.path.join(out_dir, 'results.tsv'), 'w') as f:
            f.write('\t'.join(tsv_cols) + '\n')
            for r in rows:
                f.write('\t'.join('' if r.get(c) is None else str(r.get(c))
                                  for c in tsv_cols) + '\n')

        elapsed = time.time() - t0
        # 实际跑成并投票的引擎 / 失败跳过的引擎（跨命中并集）——与「请求引擎
        # 列表」（algorithms）区分开：请求≠实际，原 summary 只记请求，会把
        # 校验器崩溃伪装成"跑过但无命中"（2026-10-03 修复 P1-4）。
        engines_actual = sorted({e for r in rows
                                 for e in (r.get('engines_actual') or [])})
        engines_skipped, engine_notes = [], []
        for r in rows:
            for s in (r.get('engines_skipped') or []):
                if s not in engines_skipped:
                    engines_skipped.append(s)
            for n in (r.get('engine_notes') or []):
                if n not in engine_notes:
                    engine_notes.append(n)
        if engines_skipped:
            logger.log('校验引擎失败跳过（结果未纳入投票）: '
                       + '; '.join(f"{s.get('engine')}: {s.get('reason')}"
                                   for s in engines_skipped), 'WARN')
        summary = {
            'tool': 'mirna',
            'n_mirna': len(mirnas), 'n_virus': len(viruses),
            'n_pairs': n_pairs, 'n_hits': len(rows),
            'consensus_levels': levels,
            'full_mode': full_mode,
            # algorithms = 请求的引擎（字段兼容保留）；engines_actual = 实际
            # 参与共识的引擎；engines_skipped = 失败跳过的引擎 + 原因
            'algorithms': algs or list(ENGINE_KEYS)[:3],
            'engines_actual': engines_actual,
            'engines_skipped': engines_skipped,
            'engine_notes': engine_notes,
            'binary': binary_info,
            'elapsed_sec': round(elapsed, 1),
            'results_json': os.path.join(out_dir, 'results.json'),
            'results_tsv': os.path.join(out_dir, 'results.tsv'),
            'out_dir': out_dir,
        }
        _pd_save_summary(out_dir, summary)
        prog('done', 1.0, f"完成：{n_pairs} 对，命中 {len(rows)}，耗时 {elapsed:.0f}s")
        logger.log(f"完成：命中 {len(rows)}/{n_pairs}，共识分布 {levels}，"
                   f"产物在 mirna_target/（results.json / results.tsv / summary.json）")
        logger.close()
        return summary
    return job


def _tool_job_pdgroup(ctx):
    """SeqGrouper 同款：GenBank 记录 → 表 + 分组映射 → Group 列 / 组计数。"""
    gb_files = ctx.opt('gb_files')
    # collection 是**集合名**（如 EXAMPLE_SET），不是文件路径 —— 之前走 ctx.opt
    # 会被当成路径做 must_exist 校验，直接 400（「请给 GenBank 文件路径或集合名」
    # 填了集合名也过不去）。名字只做字符白名单校验。
    collection = (ctx.p.get('collection') or '').strip()
    if collection and not re.fullmatch(r'[A-Za-z0-9_\-]+', collection):
        abort(400, f'无效的集合名: {collection}')
    mapping = ctx.opt('map')
    col = (ctx.p.get('col') or 'Geo Location').strip()

    def job(log, prog, cancel):
        logger = TaskLogger(callback=log)
        from Virus_Platform_Core import gb_collection as gbc
        from Virus_Platform_Core import phylodyn_kit as pk
        out_dir = _pd_out(ctx, 'pdgroup')
        prog('parse', 0.3, '解析 GenBank 记录')
        if gb_files:
            paths = [x for x in (s.strip() for s in re.split(r'[;,]+|\n', gb_files)) if x]
            text = '\n'.join(pk.read_text(p) for p in paths)
            records = gbc.parse_flatfile(text)
        elif collection:
            cdir = gbc.gb_collection_dir(collection)
            fs = [os.path.join(cdir, f) for f in sorted(os.listdir(cdir))
                  if f.lower().endswith(('.gb', '.gbk', '.genbank', '.gbff'))]
            if not fs:
                raise ValueError(f'集合 [{collection}] 内没有 GenBank 文件')
            records = gbc.parse_flatfile('\n'.join(pk.read_text(p) for p in fs))
        else:
            raise ValueError('请提供 GenBank 文件路径或集合名')
        rows = pk.gb_record_rows(records)
        gmap = pk.load_group_mapping(mapping) if mapping else None
        prog('group', 0.7, '按映射表分组')
        rows2, counts, unmatched = pk.apply_grouping(rows, gmap=gmap, col=col)
        csv_path = os.path.join(out_dir, 'grouped.csv')
        pk.save_group_csv(csv_path, rows2)
        summary = {'tool': 'pdgroup', 'n': len(rows2), 'col': col,
                   'groups': dict(counts), 'unmatched': sorted(unmatched)[:100],
                   'n_unmatched': len(unmatched), 'csv': 'pdgroup/grouped.csv',
                   'out_dir': out_dir}
        _pd_save_summary(out_dir, summary)
        prog('done', 1.0, f'{len(rows2)} 条记录，{len(counts)} 组')
        logger.close()
        return summary
    return job



def _pd_no_meta(ctx):
    """treedater 强制要日期元数据；没有时给出明确报错。"""
    raise ValueError('treedater 定年需要日期元数据 CSV（name,date）')


def _dump_json(obj, f, **kwargs):
    """动力学/工具结果落盘统一出口：先洗 NaN/±Inf → null，再写 JSON。

    ⚠️ 为什么必须做：`/tool_runs/<run>/<file>` 是**原样发文件**给浏览器的，
    不经 app.py 的 JSON provider；而 `json.dump` 默认 `allow_nan=True`，
    NaN 会以裸字面量写进磁盘 —— 前端 `r.json()`（浏览器 JSON.parse）直接
    报错。实测多个动力学卡的结果加载是 `catch(e){return}`，等于**静默不出结果**。
    `allow_nan=False` 兜底：真漏了会带类型报错，而不是悄悄产出非法 JSON。
    """
    from Virus_Platform_Core.utils import nan_to_null
    kwargs.setdefault('allow_nan', False)
    kwargs.setdefault('ensure_ascii', False)
    json.dump(nan_to_null(obj), f, **kwargs)


def _phylodyn_bool(v, default=False):
    s = str(v if v is not None else '').strip().lower()
    if not s:
        return default
    return s in ('1', 'true', 'on', 'yes', 'y', '是')


def _phylodyn_float(v, lo, hi, default):
    s = str(v if v is not None else '').strip()
    if not s:
        return default
    try:
        n = float(s)
    except (TypeError, ValueError):
        return default
    if n != n:            # NaN
        return default
    return max(lo, min(n, hi))


def _phylodyn_int(v, lo, hi, default):
    s = str(v if v is not None else '').strip()
    if not s:
        return default
    try:
        n = int(float(s))
    except (TypeError, ValueError):
        return default
    return max(lo, min(n, hi))


def _phylodyn_num(v):
    """TSV 里全是字符串：能转数字就转，转不了原样留着（空串 → None）。"""
    if v is None:
        return None
    s = str(v).strip()
    if not s:
        return None
    try:
        f = float(s)
    except ValueError:
        return s
    if f != f:
        return None
    return int(f) if f.is_integer() and abs(f) < 1e15 else f


def _phylodyn_resume(from_run, seq_len=None):
    """从上一轮运行目录重建 `prep` / `clock` 接口字典（上下联通的落地方式）。

    阶段之间**只认文件**：A1 的输入是「树 + 日期表」，A2 的输入是
    `clock/timetree.nwk` + `prep/states.csv`。所以「续跑」＝把上一轮的接口产物
    重新指一遍，不重算、不复制中间态、也不猜。

    另外把上一轮的 `*_report.json` 里的**标量结论**一起带过来（时间树来源、
    tMRCA、DRT、RTT 拟合、序列长度…）。不带的话，续跑出来的结果面板会显示
    "没有定年信息" —— 数据明明在上一轮目录里躺着，面板却是空的。
    """
    root = check_path(os.path.join(_tool_runs_root(), from_run),
                      must_exist=True, in_platform=True)

    def _report(sub, name):
        p = os.path.join(root, sub, name)
        if not os.path.isfile(p):
            return {}
        try:
            with open(p, encoding='utf-8') as f:
                d = json.load(f)
            return d if isinstance(d, dict) else {}
        except (OSError, ValueError):
            return {}

    prep = {'ok': True, 'outputs': {}, 'seq_len_modal': seq_len,
            'warnings': [], 'resumed_from': root}
    for tag, rel in (('prep_fasta', 'prep/prep.fasta'),
                     ('dates', 'prep/dates.csv'),
                     ('states', 'prep/states.csv'),
                     ('prep_nwk', 'prep/prep.nwk')):
        p = os.path.join(root, *rel.split('/'))
        if os.path.isfile(p):
            prep['outputs'][tag] = p
    pr = _report('prep', 'prep_report.json')
    for k in ('n_seqs', 'n_dates', 'n_states', 'date_span', 'state_counts',
              'tree_source', 'tree_tips', 'aln_len', 'seq_len_modal',
              'n_unmatched_meta', 'unmatched_meta'):
        if pr.get(k) is not None:
            prep[k] = pr[k]
    if seq_len is None:
        prep['seq_len_modal'] = pr.get('seq_len_modal')

    clock = {'ok': True, 'outputs': {}, 'stats': [], 'warnings': [],
             'resumed_from': root}
    for tag, rel in (('timetree_nwk', 'clock/timetree.nwk'),
                     ('divergence_nwk', 'clock/divergence_tree.nwk'),
                     ('clock_stats', 'clock/clock_stats.tsv'),
                     ('drt_json', 'clock/clock_qc/drt.json'),
                     ('rtt_scatter', 'clock/clock_qc/rtt_scatter.tsv')):
        p = os.path.join(root, *rel.split('/'))
        if os.path.isfile(p):
            clock['outputs'][tag] = p
    cr = _report('clock', 'clock_report.json')
    # ⚠️ 刻意**不**读旧报告里的 'lsd2' / 'treedater' 键：A1 已于 2026-09-17 收敛到
    # TreeTime，把历史 run 里那两条已摘除通路的数字带进新一轮 summary 会让人
    # 以为它们还在跑（口径不可比、且不再有对应的卡片可解读）。
    for k in ('timetree_source', 'root_year', 'rtt_fit', 'rtt_points', 'drt',
              'tt_env', 'stats', 'n_methods_ok', 'n_methods_usable',
              'treetime', 'timetree'):
        if cr.get(k) is not None:
            clock[k] = cr[k]
    # 时间树坐标对象：老运行目录的 clock_report.json 里没有这个键（本轮才加），
    # 但接口产物 `clock/timetree.nwk` + `prep/dates.csv` 都在 —— 现场补算一份，
    # 免得续跑出来的面板"上一轮明明定过年却画不出时间树"。
    if clock.get('timetree') is None and clock['outputs'].get('timetree_nwk'):
        try:
            from Virus_Platform_Core import phylodyn_local as _pl
            _t = _pl.build_clock_timetree(
                clock['outputs']['timetree_nwk'], prep['outputs'].get('dates'),
                source_label=clock.get('timetree_source') or '')
            if _t:
                clock['timetree'] = _t
        except Exception as _e:         # 补算失败要说出来，不许静默留空
            clock.setdefault('warnings', []).append(
                '时间树坐标未能从上一轮产物重建：%s: %s'
                % (type(_e).__name__, _e))

    missing = [t for t in ('prep_fasta', 'dates', 'states')
               if t not in prep['outputs']]
    if missing:
        # 本函数只在 phylodyn 的 job 体内被调用（任务线程无请求上下文）：
        # abort(400) 会退化成「BadRequest: 400 Bad Request」丢掉文案。
        raise ValueError('上一轮目录缺少接口产物（%s）：%s'
                         % ('、'.join(missing), root))
    return prep, clock, root


def _phylodyn_tsv(path, ncol=None):
    """读回阶段产出的 TSV（表头在首行）。读不到返回 []。

    为什么要在 job 里**读回文件**而不是直接用内存里的 dict：阶段产物是**契约**，
    前端拿到的汇总必须与磁盘上那一份一致。内存 dict 与落盘文件一旦分叉
    （比如写入时格式化成字符串），"面板显示 67、文件里是 69"这种事故查不出来。
    """
    rows = []
    if not path or not os.path.isfile(path):
        return rows
    with open(path, encoding='utf-8-sig', errors='replace') as f:
        hdr = None
        for ln in f:
            cells = ln.rstrip('\n').split('\t')
            if hdr is None:
                hdr = [c.strip() for c in cells]
                continue
            if not any(c.strip() for c in cells):
                continue
            rows.append({hdr[i]: (cells[i].strip() if i < len(cells) else '')
                         for i in range(len(hdr))})
    return rows


def _tool_job_phylodyn(ctx):
    """本地时间与地理推断（轨道 A）：A0 prep → A1 clock → A2 phylogeo
    → A3 coalescent → A4 ancestral → A5 report。

    **全程本地、零 MCMC**。这条链解决的是"没有服务器也能出时间树 + 地理迁移"
    的问题；需要贝叶斯后验的（BSP 的 95% HPD、BSSVS 的 BF）走轨道 B
    （`pgjob` 打包 → 服务器 → `beast` 导入），两边**分开出结果、不混表**。

    两种模式：
      `pipeline` —— 一键全链。默认 `prep+clock+phylogeo+report`；
                    勾 `skyline` 插 A3、勾 `ancestral` 插 A4。
      `stage`    —— 只跑一个阶段（＋A5 汇总）。配 `from_run` 可接着上一轮的
                    产物往下跑（例如上轮只跑了 prep+clock，这轮直接跑 phylogeo）。

    阶段间接口产物（**唯一契约**，改动这里等于改协议）：
      `prep/prep.fasta` `prep/dates.csv` `prep/states.csv` `prep/prep.nwk`
        → `clock/timetree.nwk`（枝长＝年）
        → `clock/divergence_tree.nwk`（枝长＝替换/位点，A4 用）
        → `phylogeo/corridors.tsv`（`ml` / `fitch` **两列并列**）
        → `coalescent/skyline.tsv`
        → `ancestral/ancestral_sequences.fasta`
        → `report/RUN_REPORT.md` + `report/run_summary.json`

    ⚠️ 两条口径红线在本 job 里被**强制保留**（不因前端勾选而消失）：
      ① `corridors.tsv` 的 `ml`（ML 点估计）与 `fitch`（最少变化数的一种悬挂）
         不是同一个统计量，总数可以接近但逐走廊方向差很多 —— 汇总里两列都给。
      ② skyline 的上下界是 TreeTime 的近似区间（±2 SD of the LH），**不是**
         贝叶斯 95% HPD；汇总里带 `is_hpd=False` 标记，前端据此换文案。
    """
    from Virus_Platform_Core import phylodyn_local as pl

    mode = (ctx.p.get('mode') or 'pipeline').strip().lower()
    if mode not in ('pipeline', 'stage'):
        abort(400, '无效的运行模式（只认 pipeline / stage）')

    stage = (ctx.p.get('stage') or '').strip().lower()
    if mode == 'stage' and stage not in pl.STAGE_ORDER:
        abort(400, '无效的阶段「%s」（可选：%s）'
                   % (stage or '空', '、'.join(pl.STAGE_ORDER)))

    from_run = (ctx.p.get('from_run') or '').strip()
    if from_run:
        # 只认运行名，不收任意路径：`../etc/passwd` 取 basename 后是 `passwd`，
        # 路径穿越在这里就被削掉了，后面拼进 tool_runs 根也不可能跑出去。
        from_run = os.path.basename(from_run.rstrip('/\\'))
        if not re.fullmatch(r'[A-Za-z0-9_\-]+', from_run):
            abort(400, '无效的上一轮运行名（只允许字母/数字/下划线/连字符）')
        if mode == 'stage':
            # 工厂期就确认目录在 —— 否则用户要等任务起来才在日志里看到
            # "上一轮不存在"，白等一次调度。
            try:
                _prev = check_path(os.path.join(_tool_runs_root(), from_run),
                                   must_exist=True, in_platform=True)
            except Exception:
                abort(400, '上一轮运行目录不存在：%s（在「结果中心 → 专项运行」'
                           '里核对运行名）' % from_run)
            for _rel in ('prep/dates.csv', 'prep/states.csv'):
                if not os.path.isfile(os.path.join(_prev, *_rel.split('/'))):
                    abort(400, '上一轮目录 %s 缺少接口产物 %s —— '
                               '它可能没跑到 A0，或不是时间与地理推断的运行目录'
                               % (from_run, _rel))

    inp = ctx.req('input', '比对 FASTA')
    meta = ctx.opt('meta')
    tree = ctx.opt('tree')
    rename_map = (ctx.opt('rename_map') or ctx.p.get('rename_map') or '').strip()
    div_tree = ctx.opt('divergence_tree')
    # 区划坐标表（可选）：**只影响 A2 迁移弧线地图的落点**，不参与任何统计量。
    # 不给时引擎用内置近似质心兜底，缺坐标的区划由引擎显式列进 coord_missing。
    coord_table = ctx.opt('coord_table')

    trait = (ctx.p.get('trait') or 'region').strip() or 'region'
    date_col = (ctx.p.get('date_col') or '').strip() or None
    date_trait = (ctx.p.get('date_trait') or '').strip() or None
    seq_len = _phylodyn_int(ctx.p.get('seq_len'), 1, 20000000, None)

    methods = tuple(m.strip().lower() for m in
                    (ctx.p.get('methods') or 'treetime').split(',')
                    if m.strip().lower() in pl.CLOCK_METHODS)
    if not methods:
        methods = ('treetime',)

    reroot = (ctx.p.get('reroot') or 'least-squares').strip() or 'least-squares'
    if reroot not in ('least-squares', 'min_dev', 'oldest', 'best', 'none'):
        abort(400, '无效的重新生根方式：%s' % reroot)
    reroot = None if reroot == 'none' else reroot
    keep_root = _phylodyn_bool(ctx.p.get('keep_root'), False)

    # 松弛钟：sigma>0 才启用；coupling 默认 0（＝非相关松弛钟，TreeTime 推荐）
    r_sigma = _phylodyn_float(ctx.p.get('relax_sigma'), 0.0, 100.0, 0.0)
    r_coup = _phylodyn_float(ctx.p.get('relax_coupling'), 0.0, 1.0, 0.0)
    relax = (r_sigma, r_coup) if r_sigma > 0 else None

    do_drt = _phylodyn_bool(ctx.p.get('drt'), True)
    drt_perm = _phylodyn_int(ctx.p.get('drt_perm'), 20, 2000, 200)

    # 固定分子钟速率（慢演化病原兜底，对齐 phymap-workflow 的 treetime_clock_rate）。
    # ⚠️ **刻意不用 `_phylodyn_float`**：它对解析失败是**静默回退默认值** ——
    # 用户手滑输错一个字符，就会从"固定速率"悄悄变成"ML 估计"，两条路的结论
    # 完全不同却没人知道。这类"静默改变语义"必须驳回。
    cr_raw = str(ctx.p.get('clock_rate') or '').strip()
    clock_rate = None
    if cr_raw and cr_raw.lower() not in ('none', 'auto', 'off', '-'):
        try:
            clock_rate = float(cr_raw)
        except (TypeError, ValueError):
            abort(400, '固定速率 clock_rate 不是合法数字：%r' % cr_raw)
        if clock_rate != clock_rate or not (1e-12 <= clock_rate <= 1.0):
            abort(400, '固定速率 clock_rate 应落在 1e-12 ~ 1.0（每位点每年），'
                       '当前：%r' % cr_raw)
        if relax:
            abort(400, '「固定速率」与「松弛钟」互斥：前者假定全树单一速率，'
                       '后者假定速率逐枝变化，同时给等于自相矛盾。请二选一'
                       '（清空 relax_sigma，或清空固定速率）。')

    want_skyline = _phylodyn_bool(ctx.p.get('skyline'), False)
    n_skyline = _phylodyn_int(ctx.p.get('n_skyline'), 3, 200, 10)
    gen_per_year = _phylodyn_float(ctx.p.get('gen_per_year'), 1e-6, 1e6, 50.0)

    # ---- 冗余剪枝（对齐上游 phymap-workflow `--prune`；默认**关**）----
    # ⚠️ 这是**有损**操作（改变样本组成 → 改变速率/tMRCA/走廊计数），所以：
    # 默认关闭；方法/粒度做白名单校验（拼错就驳回，不静默退回默认）；
    # 阈值范围校验；实际生效值写进 `params` 与 `config.resolved.json`。
    do_prune = _phylodyn_bool(ctx.p.get('prune'), False)
    prune_method = (ctx.p.get('prune_method') or 'fps').strip().lower()
    if prune_method not in pl.PRUNE_METHODS:
        abort(400, '无效的剪枝方法 %r（可选：%s）'
                   % (prune_method, '、'.join(pl.PRUNE_METHODS)))
    prune_res = (ctx.p.get('prune_date_resolution') or 'month').strip().lower()
    if prune_res not in pl.PRUNE_RESOLUTIONS:
        abort(400, '无效的剪枝日期粒度 %r（可选：%s）'
                   % (prune_res, '、'.join(pl.PRUNE_RESOLUTIONS)))
    prune_max_reps = _phylodyn_int(ctx.p.get('prune_max_reps'), 1, 200, 3)
    prune_min_snp = _phylodyn_int(ctx.p.get('prune_min_snp_diff'), 0, 10000, 2)
    prune_clade = _phylodyn_int(ctx.p.get('prune_clade_cutoff'), 0, 10000, 5)

    # 采样偏倚校正（红线 4，默认开）：'auto' / 'none' / 正数
    sbc_raw = (ctx.p.get('sampling_bias_correction') or 'auto').strip().lower()
    if sbc_raw in ('', 'auto', 'true', 'on', 'yes', '1'):
        sbc = 'auto'
    elif sbc_raw in ('none', 'off', 'false', 'no', '0'):
        sbc = None
    else:
        sbc = _phylodyn_float(sbc_raw, 1e-9, 1e9, 'auto')

    cross_check = _phylodyn_bool(ctx.p.get('cross_check'), True)
    rssp_bs = _phylodyn_int(ctx.p.get('rssp_bs'), 0, 1000, 0)
    seed = _phylodyn_int(ctx.p.get('seed'), 0, 2 ** 31 - 1, 0)

    do_ancestral = _phylodyn_bool(ctx.p.get('ancestral'), False)
    do_homoplasy = _phylodyn_bool(ctx.p.get('homoplasy'), False)
    timeout = _phylodyn_int(ctx.p.get('timeout'), 60, 86400, 3600)

    # ---- 阶段编排 ----
    if mode == 'pipeline':
        stages = ['prep', 'clock', 'phylogeo']
        if want_skyline:
            stages.append('coalescent')
        if do_ancestral:
            stages.append('ancestral')
        stages.append('report')
    else:
        stages = [stage] + (['report'] if stage != 'report' else [])

    params = {
        'mode': mode, 'stage': stage or None, 'from_run': from_run or None,
        'trait': trait, 'date_col': date_col, 'date_trait': date_trait,
        'seq_len': seq_len, 'methods': list(methods), 'reroot': reroot,
        'keep_root': keep_root, 'relax': list(relax) if relax else None,
        'drt': do_drt, 'drt_perm': drt_perm,
        'clock_rate': clock_rate,
        'prune': do_prune, 'prune_method': prune_method,
        'prune_date_resolution': prune_res, 'prune_max_reps': prune_max_reps,
        'prune_min_snp_diff': prune_min_snp, 'prune_clade_cutoff': prune_clade,
        'skyline': want_skyline, 'n_skyline': n_skyline,
        'gen_per_year': gen_per_year, 'sampling_bias_correction': sbc,
        'cross_check': cross_check, 'rssp_bs': rssp_bs, 'seed': seed,
        'ancestral': do_ancestral, 'homoplasy': do_homoplasy,
        'divergence_tree': div_tree, 'tree': tree, 'meta': meta,
        'coord_table': coord_table, 'rename_map': rename_map or None,
    }
    kw = dict(
        trait=trait, date_col=date_col, date_trait=date_trait,
        methods=methods, relax=relax, reroot=reroot or 'least-squares',
        keep_root=keep_root, do_drt=do_drt, drt_perm=drt_perm, drt_seed=42,
        clock_rate=clock_rate,
        prune=do_prune, prune_method=prune_method,
        prune_date_resolution=prune_res, prune_max_reps=prune_max_reps,
        prune_min_snp_diff=prune_min_snp, prune_clade_cutoff=prune_clade,
        sampling_bias_correction=sbc, pc=None, cross_check=cross_check,
        rssp_bs=rssp_bs, seed=seed, n_skyline=n_skyline,
        gen_per_year=gen_per_year, rng_seed=seed, timeout=timeout,
        divergence_tree=div_tree, marginal=False, coord_table=coord_table,
        rename_map=rename_map or None,
    )

    def job(log, prog, cancel):
        from Virus_Platform_Core import phylodyn_local as plm
        from Virus_Platform_Core import treetime_ml as tt

        logger = TaskLogger(callback=log)
        out_dir = ctx.run_dir
        started = time.time()
        # 进度折算用的阶段清单：续跑时会被换成"本次真正要跑的"那几个
        prog_stages = list(stages)
        # 进度**只增不减**：extras（homoplasy）与 summary 跑在 A5 之后，
        # 按"阶段序号 × 段宽"折算会落回 0.15 附近，前端进度条会往回跳。
        # 用高水位线兜住，这些收尾步骤就停在末端不动 —— 与真实执行顺序一致。
        mark = {'v': 0.0}

        def _prog(st, frac, msg):
            # 阶段进度摊到整条链上：用「已完成阶段数 + 当前阶段内进度」
            # 折算全局比例，前端那条进度条才是单调的。
            try:
                i = prog_stages.index(st)
            except ValueError:
                # 阶段外的收尾步骤：给个末端区间，别折回开头
                v = 0.95 + 0.049 * float(frac or 0)
            else:
                span = 1.0 / max(len(prog_stages), 1)
                v = min(0.999, i * span + span * float(frac or 0))
            mark['v'] = max(mark['v'], min(0.999, v))
            prog(st, mark['v'], msg or st)

        if cancel is not None and cancel.is_set():
            raise RuntimeError('已取消')

        # ---- 续跑：从上一轮目录重建接口字典 ----
        prep0 = clock0 = None
        resumed_root = None
        run_stages = list(stages)
        if from_run and mode == 'stage':
            prep0, clock0, resumed_root = _phylodyn_resume(from_run, seq_len)
            log(f'续跑上一轮：{os.path.basename(resumed_root)}')
            log('  接口产物：' + '、'.join(
                sorted(set(prep0['outputs']) | set(clock0['outputs']))))
            # prep/clock 已由上一轮提供 → 本次不重算，也不从进度条里假装跑过
            run_stages = [s for s in run_stages if s not in ('prep', 'clock')]
            if not run_stages:
                raise RuntimeError(
                    '续跑时没有可执行的阶段：prep/clock 由上一轮提供，'
                    '请把阶段选成 phylogeo / coalescent / ancestral / report')
            prog_stages = list(run_stages)

        res = plm.run_pipeline(
            inp, out_dir, meta=meta, tree=tree,
            stages=tuple(run_stages), log=log, progress=_prog,
            prep_in=prep0, clock_in=clock0,
            # ⚠️ 平台的 cancel 是 threading.Event，**不是可调用**：
            # 写成 cancel() 会炸 'Event' object is not callable（本文件
            # 1108 行早有一条同源注释）。这里把 is_set 作为回调传进引擎，
            # 让它在阶段之间检查，长链才停得下来。
            should_cancel=(cancel.is_set if cancel is not None else None),
            **kw)

        # ---- 可选：homoplasy（重复突变，TreeTime 只在 stdout 出结果）----
        extra = {}
        if do_homoplasy:
            _prog('homoplasy', 0.95, '同塑性（重复突变）统计')
            try:
                env = tt.check_tt_env()
                if not env['ok']:
                    extra['homoplasy'] = {'ok': False,
                                          'error': 'TreeTime 不可用'}
                else:
                    cout = (res['stages'].get('clock') or {}).get('outputs') or {}
                    t = cout.get('divergence_nwk') or cout.get('timetree_nwk')
                    pout = (res['stages'].get('prep') or {}).get('outputs') or {}
                    aln = pout.get('prep_fasta') or inp
                    edir = os.path.join(out_dir, 'extras')
                    os.makedirs(edir, exist_ok=True)
                    h = tt.run_homoplasy([env['exe']], aln, t, edir,
                                         rng_seed=seed, timeout=timeout,
                                         say=log)
                    # ⚠️ 该子命令**不落文件**，结果只在 stdout（`homoplasies`）；
                    # 平台这边替它落一份 TSV，否则刷新页面就没了。
                    rows = h.get('homoplasies') or []
                    hp = None
                    if rows:
                        hp = os.path.join(edir, 'homoplasy.tsv')
                        with open(hp, 'w', encoding='utf-8') as f:
                            f.write('mut\tmultiplicity\n')
                            for r in rows:
                                f.write('%s\t%s\n' % (r.get('mut'),
                                                      r.get('multiplicity')))
                    extra['homoplasy'] = {
                        'ok': h.get('ok'), 'error': h.get('error'),
                        'n_rows': len(rows), 'rows': rows[:200],
                        'outputs': ({'homoplasy_tsv': hp} if hp else {})}
                    if h.get('ok'):
                        log('同塑性：%d 个重复突变位点，最高 %s（%s 次）'
                            % (len(rows), rows[0].get('mut'),
                               rows[0].get('multiplicity')))
            except Exception as e:      # 附加产物，失败不该拖垮主链
                extra['homoplasy'] = {'ok': False,
                                      'error': '%s: %s' % (type(e).__name__, e)}
                log('[WARN] homoplasy 未完成：%s' % e)

        # ---- 汇总（从磁盘读回，保证与产物一致）----
        _prog('summary', 0.97, '汇总产物')
        pg = res['stages'].get('phylogeo') or {}
        pg_out = pg.get('outputs') or {}
        co = res['stages'].get('coalescent') or {}
        cl = res['stages'].get('clock') or {}
        cl_out = cl.get('outputs') or {}
        ps = res['stages'].get('prep') or {}
        pr = ps.get('prune') or {}
        an = res['stages'].get('ancestral') or {}
        rp = res['stages'].get('report') or {}

        corridors = []
        for r in _phylodyn_tsv(pg_out.get('corridors')):
            corridors.append({'from': r.get('from'), 'to': r.get('to'),
                              'ml': _phylodyn_num(r.get('ml')),
                              'fitch': _phylodyn_num(r.get('fitch')),
                              'delta': _phylodyn_num(r.get('delta'))})
        corridors.sort(key=lambda x: -(abs(x.get('delta') or 0)))
        events = []
        for r in _phylodyn_tsv(pg_out.get('events_year')):
            events.append({'year': _phylodyn_num(r.get('Year')),
                           'from': r.get('Origin'), 'to': r.get('Destination'),
                           'method': r.get('Method'),
                           'child': r.get('ChildNode')})
        clock_stats = []
        for r in _phylodyn_tsv(cl_out.get('clock_stats')):
            clock_stats.append({k: _phylodyn_num(v) if k not in (
                'method', 'model', 'rate_unit', 'note', 'suspect_reason')
                else v for k, v in r.items()})

        # 告警/错误/产物清单：从各阶段自身取，再从 report/run_summary.json 取产物
        # —— `stage_report` 的返回值里 `warnings` 恒为空、也没有 `files`，
        # 直接读它会得到"零告警零产物"的假象。
        warns, errs = [], []
        for s in plm.STAGE_ORDER:
            sr = res['stages'].get(s) or {}
            for w in (sr.get('warnings') or []):
                warns.append({'stage': s, 'text': str(w)})
            if sr.get('error'):
                errs.append({'stage': s, 'text': str(sr['error'])})
        files = []
        rs_path = (rp.get('outputs') or {}).get('run_summary')
        if rs_path and os.path.isfile(rs_path):
            try:
                with open(rs_path, encoding='utf-8') as f:
                    rs = json.load(f)
                files = rs.get('files') or []
            except (OSError, ValueError) as e:
                warns.append({'stage': 'report',
                              'text': 'run_summary.json 读不回：%s' % e})

        summ = {
            'tool': 'phylodyn',
            'mode': mode,
            'stages': [{'key': s,
                        'code': (plm.STAGE_META.get(s) or {}).get('code'),
                        'name': (plm.STAGE_META.get(s) or {}).get('name'),
                        'name_en': (plm.STAGE_META.get(s) or {}).get('en'),
                        'ran': s in run_stages,
                        'ok': (res['stages'].get(s) or {}).get('ok'),
                        'error': (res['stages'].get(s) or {}).get('error')}
                       for s in stages],
            'resumed_from': resumed_root,
            # A0 的列名自动择优结果（用了哪一列当日期/地点）—— 口径纪律：看得见
            'meta_pick': (res['stages'].get('prep') or {}).get('meta_pick'),
            # A0 冗余剪枝的结果（有损操作，必须让用户一眼看到剪了多少、剔了谁）
            'prune': {
                'enabled': bool(do_prune),
                'skipped': bool(pr.get('skipped')),
                'abandoned': bool(pr.get('abandoned')),
                'reason': pr.get('reason'),
                'method': pr.get('method'),
                'date_resolution': pr.get('resolution'),
                'n_in': pr.get('n_in'), 'n_out': pr.get('n_out'),
                'n_dropped': pr.get('n_dropped'),
                'n_seqs_input': ps.get('n_seqs_input'),
                'how': pr.get('how'),
                'n_groups': len(pr.get('groups') or []),
                'n_clades': pr.get('n_clades'),
                # 全量清单在 prep/prune_report.json 与 prune_dropped.tsv，
                # 这里只带前 200 条，避免 summary 被几千行分组撑爆。
                'groups': (pr.get('groups') or [])[:200],
                'dropped': (pr.get('dropped') or [])[:200],
                'n_dropped_total': len(pr.get('dropped') or []),
                'tree_pruned': ps.get('tree_pruned'),
            },
            'clock': {
                'timetree_source': cl.get('timetree_source'),
                'root_year': cl.get('root_year'),
                # 时间树坐标对象（与 LSD2/TreeDater/phylogeo 三卡同构）——
                # 前端 pdRenderClock 复用同一个 drawTimeTree 画它。
                'timetree': cl.get('timetree'),
                'stats': clock_stats,
                'n_methods_ok': cl.get('n_methods_ok'),
                'n_methods_usable': cl.get('n_methods_usable'),
                'rtt_fit': cl.get('rtt_fit'),
                'rtt_points': cl.get('rtt_points') or [],
                'drt': cl.get('drt'),
                'tt_env': cl.get('tt_env'),
            },
            'phylogeo': {
                'trait': pg.get('trait'),
                'n_changes_ml': pg.get('ml_n_changes'),
                'n_changes_fitch': pg.get('fitch_n_changes'),
                'n_corridors': pg.get('n_corridors'),
                'corridors': corridors,
                'n_events_ml': pg.get('n_events_ml'),
                'n_events_fitch': pg.get('n_events_fitch'),
                'events': events[:2000],
                'sampling_bias': pg.get('sampling_bias'),
                # ML 边际概率的置信分布：内部节点最大边际概率偏低的比例。
                # ⚠️ 这是 **ML 边际概率**，不是贝叶斯后验 —— 前端文案不得写成"后验"。
                'confidence_summary': (pg.get('ml') or {}).get('confidence_summary'),
                'edges_missing': pg.get('ml_edges_missing'),
                'align': pg.get('cross_check'),
                'rssp': pg.get('rssp'),
                # 迁移弧线地图的落点（A2 解析）：coords 是「区划 → [纬度, 经度]」，
                # coord_src 逐条标出处（user＝用户坐标表 / builtin＝内置近似质心），
                # coord_missing 是**没有任何坐标**的区划 —— 前端必须显式列出来，
                # 不许静默少画点。这三项只决定弧线画在哪，不参与任何统计量。
                'coords': pg.get('coords') or {},
                'coord_src': pg.get('coord_src') or {},
                'coord_missing': pg.get('coord_missing') or [],
                'coord_regions': pg.get('coord_regions') or {},
                'coord_table': pg.get('coord_table'),
            },
            'coalescent': {
                'ok': co.get('ok'),
                'skyline': co.get('skyline') or [],
                'skyline_meta': co.get('skyline_meta') or {},
                'rate_this_pass': co.get('rate_this_pass'),
                'r2_this_pass': co.get('r2_this_pass'),
                # 口径标记：上下界**不是**贝叶斯 95% HPD，前端据此换文案
                'is_hpd': False,
            },
            'ancestral': {'ok': an.get('ok'), 'n_mutations': an.get('n_mutations'),
                          'tree_used': an.get('tree_used')},
            'extras': extra,
            'warnings': warns,
            'errors': errs,
            'files': files,
            'n_files': rp.get('n_files'),
            'params': params,
            'seconds': round(time.time() - started, 2),
            'out_dir': out_dir,
        }
        p = os.path.join(out_dir, 'phylodyn_summary.json')
        with open(p, 'w', encoding='utf-8') as f:
            _dump_json(summ, f, indent=1)

        # ---- 日志收尾：把"看结论前必须知道的事"直接摆出来 ----
        log('—' * 46)
        log('阶段状态：' + '  '.join(
            '%s%s' % ((plm.STAGE_META.get(s) or {}).get('code') or s,
                      '✓' if (res['stages'].get(s) or {}).get('ok') else '✗')
            for s in stages))
        if cl.get('timetree_source'):
            log('时间树来源：%s' % cl['timetree_source'])
        if cl.get('root_year') is not None:
            log('tMRCA：%.2f 日历年' % cl['root_year'])
        drt = cl.get('drt') or {}
        if drt:
            log('时间信号 DRT：p=%s → %s'
                % (drt.get('p_value'),
                   '通过' if drt.get('passed') else '**不通过（不得作定年结论引用）**'))
        for s in clock_stats:
            if s.get('suspect'):
                log('[WARN] %s 的速率/tMRCA 被闸门标为可疑：%s'
                    % (s.get('method'), s.get('suspect_reason')))
        if pg:
            log('地理迁移：ML %s 次 / Fitch %s 次（**两列不可相加、不可互替**）'
                % (pg.get('ml_n_changes'), pg.get('fitch_n_changes')))
            sb = pg.get('sampling_bias') or {}
            if sb.get('shares'):
                log('采样构成：' + '、'.join(
                    '%s %.1f%%' % (k, 100 * v)
                    for k, v in sorted(sb['shares'].items(),
                                       key=lambda x: -x[1])))
                log('  采样偏倚校正：%s'
                    % ('已按系数 %.4g 校正' % sb['factor'] if sb.get('corrected')
                       else '**未校正**（方向置信度会虚高）'))
        if co.get('ok'):
            log('天际线：%d 个网格点；上下界是 ±2 SD 近似区间，'
                '**不是**贝叶斯 95%% HPD（正式带请走轨道 B 的 BEAST BSP）'
                % len(co.get('skyline') or []))
        if an.get('ok'):
            log('祖先序列：%s 处逐枝替换' % an.get('n_mutations'))
        if do_prune:
            if pr.get('skipped'):
                log('[WARN] A0 剪枝已跳过：%s' % pr.get('reason'))
            elif pr.get('abandoned'):
                log('[WARN] A0 剪枝已放弃（剪完不足 3 条）→ 用的是**全量样本**')
            elif pr.get('n_out') is not None:
                log('A0 冗余剪枝（%s / %s）：%s → %s 条（剔 %s）；'
                    '本轮速率/tMRCA/走廊计数都在**剪枝后样本**上得出'
                    % (pr.get('method'), pr.get('resolution'),
                       pr.get('n_in'), pr.get('n_out'), pr.get('n_dropped')))
        log('产物 %s 个，清单见 report/run_summary.json' % (rp.get('n_files'),))
        log('汇总：phylodyn_summary.json')
        logger.close()
        return summ

    return job
