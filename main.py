# -*- coding: utf-8 -*-
"""
植物病毒分析平台 - 命令行入口
用法:
  python main.py init-taxonomy
  python main.py build-host-db --genome host-db/112863_Lycium_barbarum/genome.fa --taxid 112863
  python main.py build-virus-db --fasta databases/virusref_db/final.cluster.ref.fasta --info databases/virusref_db/final.cluster.ref_info.tsv
  python main.py analyze --r1 a_R1.fq.gz --r2 a_R2.fq.gz --sample NX-5
  python main.py report --sample NX-5        # 重新生成可视化报告
  python main.py logan-create --name 查询名 --sample NX-5   # LOGAN 溯源查询
  python main.py gb-dl --acc "NC_001367,NC_002692" -n tobamo  # 下载 GenBank 集合
图形界面: 双击 启动平台-桌面窗口.bat
"""
import os
import re
import sys
import argparse

if sys.platform == 'win32':
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
    sys.stderr.reconfigure(encoding='utf-8', errors='replace')

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from Virus_Platform_Core.config import get_config, DIRS, PLATFORM_ROOT, engine_cmd
from Virus_Platform_Core.selfcheck import (FALLBACK_MODULES,
                                           iter_platform_modules
                                           as _iter_platform_modules)
from Virus_Platform_Core.utils import TaskLogger, check_path
from Virus_Platform_Core.pipeline import DEFAULT_ANALYZE_STAGES


def _bootstrap_taxonomy():
    """First-run extraction of the bundled NCBI taxonomy archive.

    The database package ships NCBI's official new_taxdump.tar.gz (155 MB) and
    app.py unpacks it in a background thread when the web UI starts.  Buyers are
    told to run `VirusPlatform.exe --cli selfcheck` as their very first command,
    so without this bootstrap that command reported "NCBI Taxonomy: not ready"
    (exit 1) on a cold product until the web UI had been opened once.
    Extraction is pure local work - no network access, about a minute.
    """
    try:
        from Virus_Platform_Core.config import DIRS
        from Virus_Platform_Core.taxonomy import prepare_taxonomy, taxonomy_ready
        if taxonomy_ready():
            return
        tar = os.path.join(DIRS['databases'], 'new_taxdump.tar.gz')
        if not os.path.isfile(tar):
            return
        _say('[init] extracting the bundled NCBI taxonomy archive '
              '(first run only, about a minute)...')
        prepare_taxonomy()
        _say('[init] NCBI taxonomy is ready.')
    except Exception as e:                              # noqa: BLE001
        _say(f'[init] taxonomy bootstrap skipped: {e}')


def _self_cmd(*parts):
    """提示语里的自调用命令：分发版（冻结 exe）用 VirusPlatform.exe --cli，
    源码树用 python main.py —— 否则买家照着提示跑会找不到 python。"""
    base = ('VirusPlatform.exe --cli' if getattr(sys, 'frozen', False)
            else 'python main.py')
    return ' '.join((base,) + parts)


def _say(*args, **kw):
    """CLI 输出统一出口：en 模式下把消息翻成英文（zh 原样）。"""
    try:
        from Virus_Platform_Core.i18n_cli import tr as _tr
        if args:
            args = tuple(_tr(a) if isinstance(a, str) else a for a in args)
    except Exception:
        pass
    print(*args, **kw)


def make_logger(name, echo=True):
    log_file = os.path.join(DIRS['logs'], f'{name}_{__import__("time").strftime("%Y%m%d_%H%M%S")}.log')
    return TaskLogger(log_file, echo=echo)


def cmd_init_taxonomy(args):
    from Virus_Platform_Core.taxonomy import prepare_taxonomy, taxonomy_ready
    logger = make_logger('taxonomy')
    prepare_taxonomy(logger=logger, force=args.force)
    _say("\n✔ Taxonomy 就绪:", taxonomy_ready())
    logger.close()


def cmd_build_host_db(args):
    from Virus_Platform_Core.kunpeng import build_host_db, db_ready
    from Virus_Platform_Core.config import host_db_info, host_db_name
    logger = make_logger('build_host_db')
    check_path(args.genome, must_exist=True)
    # 建库内存峰值与线程数成正比（约 1GB/线程），默认限 8
    threads = args.threads if args.threads else min(8, get_config().threads)
    # 输出目录：--out-dir 未给则按物种命名 host-db/<taxid>_<源目录名>_host_db
    out = build_host_db(args.genome, args.taxid, db_dir=args.out_dir,
                        hash_capacity=args.hash_capacity, threads=threads,
                        logger=logger, rebuild=args.rebuild,
                        clean_mid=args.clean_mid)
    cfg = get_config()
    if db_ready(out):
        cfg.set_active_host_db(out)          # 建完即设为当前，避免建了不用
        info = host_db_info(out)
        _say(f"\n✔ 宿主库就绪: {out}")
        _say(f"  物种 = {info.get('species') or '(未查到学名)'} "
              f"(taxid={info.get('taxid')})")
        if info.get('conflicted'):
            _say(f"  ⚠ seqid2taxid.map 含多个 taxid {info.get('taxids_in_map')}"
                  f" —— 元数据被历史建库污染，建议换新目录重建")
        _say(f"  已设为当前宿主库（platform.json active_host_db = "
              f"{cfg.active_host_db}）")
    else:
        _say("\n✘ 建库失败")
    logger.close()


def cmd_build_virus_db(args):
    from Virus_Platform_Core.kunpeng import build_virus_db, db_ready
    logger = make_logger('build_virus_db')
    build_virus_db(args.fasta, args.info,
                   hash_capacity=args.hash_capacity, threads=args.threads,
                   logger=logger, rebuild=args.rebuild)
    _say("\n✔ 病毒库就绪" if db_ready(get_config().databases['virus']) else "\n✘ 建库失败")
    logger.close()


def cmd_analyze(args):
    from Virus_Platform_Core.pipeline import run_analysis
    check_path(args.r1, must_exist=True)
    if args.r2:
        check_path(args.r2, must_exist=True)
    logger = make_logger(f'analyze_{args.sample}')
    run_analysis(
        sample=args.sample, r1=args.r1, r2=args.r2,
        stages=args.stages.split(','),
        db_host=args.db_host, db_virus=args.db_virus,
        threads=args.threads, confidence=args.confidence,
        assembly_mode=args.assembly_mode, memory_gb=args.memory,
        subsample=args.subsample, min_contig_len=args.min_contig_len,
        subsample_mode=args.subsample_mode, subsample_seed=args.subsample_seed,
        subsample_prop=args.subsample_p,
        top_n_refs=args.top_n_refs, tree_tool=args.tree_tool,
        tree_sampling=args.tree_sampling,
        ncbi_refs=args.ncbi_refs,
        primer_mode=args.primer_mode, do_specificity=args.specificity,
        min_orf_aa=args.min_orf_aa, force=args.force,
        do_fq2fa=not args.no_fq2fa,
        gbdraw_max=args.gbdraw_max,
        gbdraw_fasta=args.gbdraw_fasta, gbdraw_ann=args.gbdraw_ann,
        plot_engine=args.plot_engine,
        logger=logger)
    logger.close()


def cmd_ncbi_dl(args):
    from Virus_Platform_Core.ncbi_download import download_collection, list_collections
    logger = make_logger(f'ncbi_dl_{args.name}')
    res = download_collection(args.term, args.name, db=args.db,
                              max_records=args.max, logger=logger)
    _say(f"\n✔ 集合 [{res['name']}]: 命中 {res['total_hits']}，"
          f"新增 {res['downloaded']}，跳过 {res['skipped']}")
    _say("已下载集合：")
    for c in list_collections():
        _say(f"  {c['name']:20s} {c['n_seqs']:6d} 条  {c['query'][:50]}")
    logger.close()


def cmd_ncbi_list(args):
    from Virus_Platform_Core.ncbi_download import list_collections
    cols = list_collections()
    if not cols:
        _say("（无已下载集合，先用 ncbi-dl 下载）")
        return
    for c in cols:
        _say(f"  {c['name']:20s} {c['n_seqs']:6d} 条  db={c['db']}  {c['date']}  "
              f"{c['query'][:60]}")


def cmd_gb_dl(args):
    from Virus_Platform_Core.gb_collection import download_gb_collection, list_gb_collections
    logger = make_logger(f'gb_dl_{args.name}')
    res = download_gb_collection(args.name, term=args.term,
                                 accessions=args.acc, max_records=args.max,
                                 logger=logger)
    tail = f"，未找到 {len(res['missing'])}" if res['missing'] else ''
    _say(f"\n✔ 集合 [{res['name']}]: 命中 {res['total_hits']}，"
          f"新增 {res['downloaded']}，跳过 {res['skipped']}{tail}")
    for c in list_gb_collections():
        _say(f"  {c['name']:22s} {c['n_records']:4d} 条  {c['date']}  "
              f"{c['source']}")
    logger.close()


def cmd_gb_import(args):
    from Virus_Platform_Core.gb_collection import import_local_gb
    logger = make_logger(f'gb_import_{args.name}')
    res = import_local_gb(args.name, args.files, logger=logger)
    _say(f"\n✔ 集合 [{res['name']}]: 导入 {res['imported']} 条，"
          f"跳过 {res['skipped']} 条")
    logger.close()


def cmd_gb_list(args):
    from Virus_Platform_Core.gb_collection import list_gb_collections
    cols = list_gb_collections()
    if not cols:
        _say('（无 GenBank 集合，先用 gb-dl 下载或 gb-import 导入）')
        return
    for c in cols:
        if c['source'] == 'query':
            src = c['term']
        elif c['source'] == 'accessions':
            src = f"{c['accessions']} 个 accession"
        else:
            src = '本机导入'
        _say(f"  {c['name']:22s} {c['n_records']:4d} 条  {src[:52]:52s} "
              f"{c['date']}")


def cmd_gb_check(args):
    from Virus_Platform_Core.gb_collection import inspect_collection
    st = inspect_collection(args.name)
    _say(f"集合 [{st['name']}]（{st['dir']}）: {len(st['records'])} 条记录\n")
    _say(f"{'accession':18s}{'长度':>10s}{'CDS':>5s}{'成熟肽':>6s}  物种")
    for r in st['records']:
        _say(f"{r['acc'][:17]:18s}{r['length']:10d}{r['cds_count']:5d}"
              f"{r['mat_peptides']:6d}  {(r['organism'] or r['title'])[:44]}")
    if st['warnings']:
        _say('\n⚠ 警告:')
        for w in st['warnings']:
            _say(f'  - {w}')


def cmd_report(args):
    from Virus_Platform_Core.pipeline import run_report_only
    logger = make_logger(f'report_{args.sample}')
    run_report_only(args.sample, logger=logger, force=args.force)
    logger.close()


def cmd_host_analysis(args):
    import os
    from Virus_Platform_Core.config import DIRS
    from Virus_Platform_Core.pipeline import _safe_sample_name
    from Virus_Platform_Core.host_analysis import predict_hosts
    sample = _safe_sample_name(args.sample)
    sample_dir = os.path.join(DIRS['results'], sample)
    logger = make_logger(f'hostana_{sample}')
    s = predict_hosts(sample_dir, logger=logger, force=args.force,
                      min_confidence=args.min_confidence)
    logger.close()
    _say(f"宿主预测完成: {s['n_contigs']} 条 contigs")
    if s.get('min_confidence') is not None:
        _say(f"置信度门槛 {s['min_confidence']}: "
              f"{s.get('n_low_confidence', 0)} 条标注 Low_confidence")
    for cat, n in s['categories'].items():
        _say(f"  {cat:15s}: {n}")


def cmd_orfa(args):
    import os
    from Virus_Platform_Core.config import DIRS
    from Virus_Platform_Core.pipeline import _safe_sample_name
    from Virus_Platform_Core.orf_annot import run_orf_annotation, LIB_LABELS
    sample = _safe_sample_name(args.sample)
    sample_dir = os.path.join(DIRS['results'], sample)
    logger = make_logger(f'orfa_{sample}')
    s = run_orf_annotation(sample_dir, threads=args.threads, logger=logger,
                           force=args.force, libs=args.libs,
                           min_pident=args.min_pident,
                           min_qcov=args.min_qcov)
    logger.close()
    used = s.get('libs') or []
    _say("启用注释层: " + '、'.join(LIB_LABELS.get(l, l) for l in used))
    eng = s.get('engine') or '（层1 未启用）'
    _say(f"ORF 功能注释完成: {s['n_annotated']}/{s['n_orfs']} 个 ORF 获得"
          f"有效注释，覆盖 {s['n_families']} 个病毒科（引擎 {eng}）")
    for c, n in list(s['categories'].items())[:10]:
        _say(f"  {c:20s}: {n}")


def cmd_logan_create(args):
    from Virus_Platform_Core.logan_trace import create_job, get_segment_fasta
    d = create_job(args.name, sample=args.sample or None,
                   contig_ids=(args.contigs.split(',') if args.contigs else None),
                   pasted=args.fasta, n_seg=args.segments)
    _say(f"查询任务已创建: {d['name']}（{d['n_segments']} 个片段）")
    _say("请把每个片段的序列提交到 https://logan-search.org/dashboard"
          "（Groups 建议 All_No_viral_human，建议填邮箱拿下载链接）：")
    for s in d['segments']:
        _h, seq = get_segment_fasta(d['name'], s['index'])
        fa = os.path.join(DIRS['logan'], d['name'], s['file'])
        _say(f"  s{s['index']}  {s['contig']} {s['start']}-{s['end']}  → {fa}")
    _say('拿到结果表后：' + _self_cmd('logan-import --name '
                                         f"{d['name']}")
          + ' --segment 1 --result 结果.tsv')


def cmd_logan_import(args):
    from Virus_Platform_Core.logan_trace import import_result
    raw = open(check_path(args.result, must_exist=True), 'rb').read()
    d = import_result(args.name, args.segment, args.result, raw)
    _say(f"已导入片段 s{args.segment}（{args.result}）；"
          f"任务状态: {d['status']}（{d['n_imported']}/{d['n_segments']} 片段）")
    _say(f"溯源报告: {os.path.join(DIRS['logan'], d['name'], 'trace_report.html')}")


def cmd_logan_jobs(args):
    from Virus_Platform_Core.logan_trace import list_jobs
    jobs = list_jobs()
    if not jobs:
        _say("暂无 LOGAN 溯源查询任务")
        return
    for j in jobs:
        _say(f"{j['name']:30s} 样品={j['sample'] or '—':12s} "
              f"片段 {j['n_imported']}/{j['n_segments']} 已导入  "
              f"{'报告✔' if j['has_report'] else ''}")


def cmd_logan_batch(args):
    from Virus_Platform_Core.logan_trace import batch_submit
    emails = [e.strip() for e in args.email.split(',') if e.strip()]
    logger = make_logger(f'logan_batch_{args.name}')
    r = batch_submit(args.name, emails, group=args.group,
                     headless=not args.show_browser,
                     first_wait=args.first_wait, max_wait=args.max_wait,
                     logger=logger)
    logger.close()
    _say(f"批量完成: 待提交 {r['n_pending']} 片段, 回收导入 {r['n_imported']} 个"
          f"（未落定的片段可再次运行补漏）")


def cmd_submit_list(args):
    from Virus_Platform_Core.ncbi_submit import store
    tabs = store.list_tables()
    if not tabs:
        _say('暂无提交项目（submissions/ 为空）。新建：'
          + _self_cmd('submit-init --name demo --demo'))
        return
    for t in tabs:
        extra = [f for f in t['files'] if f != 'unified_metadata.csv']
        _say(f"{t['name']:30s} {t['rows']:4d} 行  产物: "
              f"{', '.join(extra) if extra else '—（submit-export 生成）'}")


def cmd_submit_init(args):
    import subprocess as _sp
    from Virus_Platform_Core.ncbi_submit import store
    from Virus_Platform_Core.ncbi_submit import unified_metadata as _um_path
    if args.demo:
        store.create_table(args.name, sample='demo')
        _say(f"提交项目已创建（示例数据）: submissions/{args.name}/")
        return
    if args.import_csv:
        _, n = store.import_table(args.name, check_path(args.import_csv,
                                                        must_exist=True))
        _say(f"提交项目已创建（导入 {n} 行）: submissions/{args.name}/")
        return
    if not args.taxonomy:
        raise SystemExit("需要 --taxonomy <taxonomy.tsv>（contig<TAB>taxonomy），"
                         "或 --demo / --import-csv")
    d = store.table_dir(args.name)
    os.makedirs(d, exist_ok=True)
    cmd = engine_cmd(_um_path.__file__,
           '--taxonomy', check_path(args.taxonomy, must_exist=True),
           '--run-title', args.name, '-o', d)
    for flag, val in (('--metadata', args.metadata), ('--authors', args.authors),
                      ('--title', args.title), ('--bioproject', args.bioproject),
                      ('--host', args.host), ('--lat-lon', args.lat_lon),
                      ('--sequencer', args.sequencer), ('--assembler', args.assembler),
                      ('--coverage', args.coverage)):
        if val:
            cmd += [flag, val]
    r = _sp.run(cmd)
    if r.returncode == 0:
        _say(f"提交项目已创建: submissions/{args.name}/（unified_metadata.csv + source.src 等）")
    _say('继续: ' + _self_cmd(f'submit-validate --name {args.name}')
          + ' && ' + _self_cmd(f'submit-export --name {args.name}'))


def cmd_submit_validate(args):
    from Virus_Platform_Core.ncbi_submit import store
    issues = store.validate_table(args.name)
    if not issues:
        _say("✓ 所有必填字段已填写，可以提交")
        return
    # validate_table 的 issue 形状：{column, kind, count, desc, examples}
    # kind ∈ missing_col / placeholder / format / duplicate / blank
    _KIND_LABEL = {'missing_col': '缺少列', 'placeholder': '必填未填',
                   'format': '格式', 'duplicate': '重复'}
    _say(f"校验：{len(issues)} 个问题")
    for it in issues:
        kind = it.get('kind', '')
        n = it.get('count')
        label = _KIND_LABEL.get(kind, kind or '问题')
        if kind == 'missing_col' or n == -1:
            _say(f"  [{label}] {it['column']} — {it['desc']}")
            continue
        ex = ', '.join(repr(e) for e in (it.get('examples') or [])[:3])
        _say(f"  [{label} {n} 条] {it['column']} — {it['desc']}")
        if ex:
            _say(f"      例: {ex}")
    _say("修复占位符后再提交 GenBank（Web 端「提交准备」页可在线编辑）")
    # 有未解决问题时以退出码 1 结束（main() 对返回 False 统一 sys.exit(1)），
    # 供脚本/CI 串联 `submit-validate && submit-export` 时正确短路。
    return False


def cmd_submit_fill(args):
    from Virus_Platform_Core.ncbi_submit import store
    n = store.batch_fill(args.name, args.column, args.value, old_value=args.old)
    _say(f"已更新 '{args.column}' 列 {n} 个单元格")


def cmd_submit_export(args):
    from Virus_Platform_Core.ncbi_submit import store
    out = store.export_files(args.name, assembler=args.assembler,
                             sequencer=args.sequencer,
                             enrichment=args.enrichment)
    _say("产物已生成：")
    for k, v in out.items():
        _say(f"  {k:20s} {v}")
    _say("提交：NCBI BankIt 上传 source.src / .fsa+.tbl，或 Sequin 导入；"
          "BioSample 用 biosample_template.tsv 批量注册。")


def cmd_submit_adopt(args):
    """MMPV-RNA discovery 产物一键收养为提交项目（class_KEEP.fasta + ref_info.tsv）。"""
    from Virus_Platform_Core.ncbi_submit import adopt_mmpv
    argv = []
    for flag in ('out_dir', 'name', 'dataset', 'samples', 'run_name',
                 'from_server', 'server_host', 'server_data_root',
                 'isolate_prefix', 'public_metadata', 'min_length',
                 'geo_loc', 'lat_lon', 'collection_date', 'host', 'tissue',
                 'cultivar', 'dev_stage', 'collected_by', 'isolation_source',
                 'authors', 'title', 'bioproject', 'biosample_prefix',
                 'assembler', 'sequencer', 'coverage', 'annotation_pipeline',
                 'enrichment'):
        val = getattr(args, flag, None)
        if val not in (None, ''):
            argv += ['--' + flag.replace('_', '-'), str(val)]
    for flag in ('force', 'dry_run', 'no_export'):
        if getattr(args, flag, False):
            argv.append('--' + flag.replace('_', '-'))
    return adopt_mmpv.main(argv)


def cmd_submit_sbt(args):
    from Virus_Platform_Core.ncbi_submit import store
    fields = {'last': args.last, 'first': args.first, 'middle': args.middle or '',
              'affil': args.affil, 'div': args.div or '', 'city': args.city,
              'sub': args.sub or '', 'country': args.country,
              'street': args.street or '', 'email': args.email,
              'postal': args.postal or ''}
    extra = [ln for ln in (args.extra_authors or '').splitlines() if ln.strip()] \
        if args.extra_authors else []
    tpl, path = store.generate_sbt(fields, extra_authors=extra,
                                   title=args.title or '', out_name=args.name)
    _say(f"template.sbt 已生成: {path}（{1 + len(extra)} 位作者）")


def _subproc_env(deepseek_api=None, ncbi_api=None):
    """子进程环境：API Key 经环境变量注入，不落命令行（进程表/日志可见）。"""
    env = dict(os.environ)
    if deepseek_api:
        env['DEEPSEEK_API_KEY'] = str(deepseek_api)
    if ncbi_api:
        env['NCBI_API_KEY'] = str(ncbi_api)
    return env


def cmd_meta_search(args):
    import subprocess as _sp
    from Virus_Platform_Core.public_meta import search_engine
    cmd = engine_cmd(search_engine.__file__,
           '-q', args.species, '-s', args.source, '--db', args.db)
    if args.out:
        cmd += ['-o', args.out]
    if args.no_detailed:
        cmd.append('--no-detailed')
    # key 走环境变量（search_engine 侧优先读 env）
    env = _subproc_env(args.deepseek_api,
                       args.ncbi_api or os.environ.get('NCBI_API_KEY'))
    r = _sp.run(cmd, env=env)
    if r.returncode == 0:
        _say('下一步: ' + _self_cmd('meta-info --runs <输出目录>/sra.list'))


def cmd_meta_info(args):
    import subprocess as _sp
    from Virus_Platform_Core.public_meta import info_engine
    # --runs 既可以是编号列表文件，也可以是单个 SRR/CRR 编号
    if os.path.isfile(args.runs):
        runs_arg = check_path(args.runs, must_exist=True)
    elif re.fullmatch(r'[A-Za-z]\w+', args.runs or ''):
        runs_arg = args.runs
    else:
        raise SystemExit(f"--runs 需为编号列表文件或 Run 编号: {args.runs}")
    cmd = engine_cmd(info_engine.__file__,
           '-i', runs_arg, '-m', args.mode, '-t', str(args.threads))
    if args.out:
        cmd += ['-o', args.out]
    if args.fill_date:
        cmd.append('--fill-date')
    # key 走环境变量（info_engine 侧优先读 env）
    env = _subproc_env(args.deepseek_api, args.ncbi_api)
    _sp.run(cmd, env=env)


def cmd_meta_plot(args):
    from Virus_Platform_Core.public_meta.landscape_plot import plot_sci_landscape
    plot_sci_landscape(check_path(args.input, must_exist=True), args.out)


def cmd_host_genome(args):
    from Virus_Platform_Core.public_meta.host_genome import main as _hg_main
    sys.argv = ['host_genome.py',
                '--species', args.species or '',
                '--source', args.source]
    if args.accession:
        sys.argv += ['--accession', args.accession]
    if args.chromosomes:
        sys.argv += ['--chromosomes', args.chromosomes]
    if args.out:
        sys.argv += ['--outdir', args.out]
    if args.release:
        sys.argv += ['--release', str(args.release)]
    if args.taxid:
        sys.argv += ['--taxid', str(args.taxid)]
    if args.include_organelles:
        sys.argv.append('--include-organelles')
    if args.ncbi_api:
        sys.argv += ['--ncbi-api', args.ncbi_api]
    if args.min_length:
        sys.argv += ['--min-length', str(args.min_length)]
    if args.verify_only:
        sys.argv.append('--verify-only')
    if args.skip_datasets:
        sys.argv.append('--skip-datasets')
    _hg_main()


def cmd_datasets(args):
    """透传给 datasets_cli 子命令（genome/gene/taxonomy/virus）。"""
    from Virus_Platform_Core.public_meta.datasets_cli import main as _dm
    sys.argv = ['datasets.py'] + [a for a in args.ds_args if a != '--']
    _dm()


def cmd_datasets_setup(args):
    from Virus_Platform_Core.public_meta.datasets_cli import setup
    setup(log=print)


def cmd_ref_status(args):
    from Virus_Platform_Core import virus_ref
    st = virus_ref.status()
    if not st['dir']:
        _say("✘ 参考库未部署（databases/virus_ref 已移除）")
        return
    _say(f"参考库目录 : {st['dir']}")
    _say(f"版本       : {st['version']}（最后增量 {st['last_incremental']}）")
    _say(f"数据来源   : ICTV {st['source_ictv']} / NCBI {st['source_ncbi']}")
    _say(f"非冗余参考 : {st['n_ref']} 条")
    _say(f"完整基因组 : {st['n_complete']} 条")
    _say(f"规范化元数据: {'✔' if st['meta_ready'] else '✘（首次使用时自动构建）'}")


def cmd_ictv_update(args):
    from Virus_Platform_Core.ictv_db import update
    logger = make_logger('ictv_update')
    update(logger=logger, xlsx=args.xlsx)
    cmd_ictv_status(args)
    logger.close()


def cmd_ictv_status(args):
    from Virus_Platform_Core import ictv_db
    st = ictv_db.status()
    if not st['dir']:
        _say("✘ 未找到 databases/tree_db/ictv_tree.db/（ICTV 参考库未部署）")
        return
    _say(f"ICTV 库目录 : {st['dir']}")
    if not st['taxa_ready']:
        _say("  ✘ taxa.txt 未生成（用 ictv-update 更新 VMR 后自动解析）")
        return
    _say(f"MSL 版本    : {st['msl'] or '?'}（VMR: {st['xlsx'] or '?'}）")
    _say(f"解析时间    : {st['parsed'] or '?'}")
    _say(f"accession   : {st['n_rows']} 条（{st['n_genus']} 属 / {st['n_family']} 科）")
    _say(f"本地覆盖    : 本地已有 {st['n_in_acvirus']} 条，"
          f"gb 缓存已下 {st['n_gb_cache']} 条")


def cmd_ictv_refs(args):
    from Virus_Platform_Core import ictv_db
    from Virus_Platform_Core import acvirus
    db = getattr(args, 'db', None) or acvirus.DEFAULT_DB
    if db not in acvirus._LIBS:
        raise SystemExit(f'未知建树库: {db}（可选 {"/".join(acvirus._LIBS)}）')
    if not acvirus.available(db):
        raise SystemExit(f'{db} 建树库未就绪（缺少谱系表 taxa，'
                         f'见 tree_db/{db}_tree.db）')
    if not any((args.genus, args.family, args.species)):
        raise SystemExit('至少指定 --genus / --family / --species 之一')
    rows, total = ictv_db.select_refs(
        genus=args.genus, family=args.family, species=args.species,
        limit=args.limit, genome='any' if args.any else 'complete', db=db)
    scope = args.genus or args.family or args.species
    _say(f"[{scope}] 命中 {total} 条 库={db}"
          + ('' if args.any else '（仅 Complete genome）')
          + f"，按本地优先显示前 {len(rows)} 条:\n")
    _say(f"{'accession':12s}{'来源':11s}{'基因组':22s}{'物种':40s}宿主组")
    for r in rows:
        _say(f"{r['Genbank']:12s}{r['Source']:11s}"
              f"{(r['Genome_Coverage'] or '')[:20]:22s}"
              f"{(r['Species'] or '')[:38]:40s}{(r['Host_Source'] or '')[:16]}")
    need = [r['Genbank'] for r in rows if r['Source'] == 'ncbi']
    if not need:
        return
    _say(f"\n其中 {len(need)} 条本地没有（本地缓存未命中）")
    if not args.download:
        _say("预览模式未下载；加 --download 执行按需下载（需联网）")
        return
    logger = make_logger('ictv_refs')
    res = ictv_db.ensure_gb(need, logger=logger)
    tail = f"，未找到 {len(res['missing'])}" if res['missing'] else ''
    _say(f"✔ 下载 {res['downloaded']} 条，本地已有跳过 {res['skipped']} 条{tail}")
    _say(f"缓存 fasta: {res['fasta']}")
    logger.close()


def cmd_tools(args):
    cfg = get_config()
    _say(f"线程数: {cfg.threads}")
    for name, path in sorted(cfg.tool_status().items()):
        _say(f"  {name:12s} {'✔ ' + path if path else '✘ 未找到'}")
    _say(f"宿主库: {cfg.databases['host']}")
    _say(f"病毒库: {cfg.databases['virus']}")


def cmd_db_migrate(args):
    """数据库迁移/对接：把 databases（含 virusref_db/）与 host-db 复制或搬移到
    外置位置，并切换 database_root，实现软件与数据库分离部署。"""
    from Virus_Platform_Core.db_migrate import migrate, check
    if args.check:
        r = check(args.to)
        _say(('\n✔ ' if r.get('ok') else '\n✘ ') + (r.get('detail')
              or r.get('error', '')))
        return None if r.get('ok') else False
    r = migrate(args.to, mode=args.mode, dry_run=args.dry_run)
    if r.get('dry_run'):
        _say('\n' + r['detail'])
        return None
    _say(('\n✔ ' if r.get('ok') else '\n✘ 迁移失败: ') + (r.get('detail')
          or r.get('error', '')))
    return None if r.get('ok') else False


def cmd_samples_backfill(args):
    """给缺 project.json 的老样品从既有产物回填清单。

    project.json 由跑流程时写；在此之前的样品（或手工拷进 results/ 的目录）
    没有清单，项目名/输入/stages_done 就缺失。这里把
    pipeline.backfill_manifest_from_products 的反推能力暴露出来——
    该函数此前没有任何调用点。
    只填 project / input / stages_done；不猜 last_run / last_status
    （它们描述"一次运行"，凭产物推断属臆测）。
    """
    get_config()  # 确保 DIRS 已按 platform.json（含自定义输出根）就绪
    from Virus_Platform_Core.pipeline import backfill_all_manifests
    r = backfill_all_manifests(dry_run=args.dry_run, only=args.sample)
    _say(f"扫描样品 {r['total']} 个，需回填 {r['backfilled']} 个，"
          f"已有清单跳过 {r['skipped']} 个"
          + ('（预演，未写盘）' if r['dry_run'] else ''))
    for it in r['items']:
        flag = '将写入' if r['dry_run'] else ('已写入' if it['wrote'] else '无需写')
        _say(f"  [{flag}] {it['sample']:32s} 项目={it['project'] or '-'} "
              f"已完成阶段={it['stages_done']}")
    if not r['items']:
        _say('  没有需要回填的样品')
    return None


def cmd_tool_runs(args):
    """输出目录管理：status / organize / archive / clean 四动作。"""
    get_config()  # 确保 DIRS 已按 platform.json（含自定义输出根）就绪
    from Virus_Platform_Core import tool_runs_admin as tra
    act = args.action
    if act == 'status':
        info = tra.status()
        _say(f"根目录: {info['root']}")
        act_sz = sum(a['size'] for a in info['active'])
        _say(f"\n活动区（动态运行，最新在前，共 {len(info['active'])} 个，"
              f"{tra._fmt(act_sz)}）:")
        for a in info['active'][:15]:
            _say(f"  {a['name']:44s} {a['size_h']:>10s}  {a['files']} 文件")
        if len(info['active']) > 15:
            _say(f"  … 另有 {len(info['active']) - 15} 个，详见 GUI 运行目录页")
        _say('\n固定区:')
        for d, f in info['fixed'].items():
            _say(f"  {d:12s} {f['dirs']:>3d} 子目录  {f['files']:>4d} 文件  "
                  f"{f['size_h']}")
        return None
    if act == 'organize':
        r = tra.organize(dry_run=args.dry_run)
        _say(f"\n{'预演' if args.dry_run else '完成'}: 归位 {r['moved']} 项，"
              f"跳过 {r['skipped']} 项")
        return None
    if act == 'archive':
        r = tra.archive(before=args.before, older_days=args.older_days,
                        dry_run=args.dry_run)
        _say(f"\n{'预演' if args.dry_run else '完成'}: 归档 {r['archived']} 个运行")
        return None
    if act == 'clean':
        if args.active_days is None and args.archive_days is None:
            _say('✘ clean 需至少给一个期限: --active-days N 和/或 --archive-days N')
            return False
        r = tra.clean(active_days=args.active_days,
                      archive_days=args.archive_days, dry_run=args.dry_run)
        _say(f"\n{'预演' if args.dry_run else '完成'}: 删除 {r['removed']} 项，"
              f"释放 {r['freed_h']}")
        return None


def cmd_crispr_design(args):
    """CRISPR/Cas12a 检测试剂设计（与 GUI 工具卡 t-crispr 同引擎同产物）。"""
    import time as _time
    from Virus_Platform_Core import crispr_design as cd
    get_config()                      # 确保 DIRS / 工具探测就绪
    out_dir = args.out or os.path.join(
        DIRS['tool_runs'],
        'crispr_cli_' + _time.strftime('%Y%m%d_%H%M%S'))
    logger = make_logger('crispr_design')
    params = {
        'amp_mode': args.amp, 'profile': args.profile,
        'nuclease': args.nuclease, 'probe_type': args.probe,
        'spacer_len': args.spacer_len, 'top_n': args.top_n,
        'flank': args.flank, 'site_margin': args.margin,
        'use_host': '0' if args.no_host else '1',
        'group_by_segment': '0' if args.no_group else '1',
        'near_fasta': args.near or '',
        'threads': args.threads or 0,
    }
    if args.product_min:
        params['product_min'] = args.product_min
    if args.product_max:
        params['product_max'] = args.product_max

    def _prog(stage, frac, msg=''):
        _say(f'  [{stage:8s}] {int(frac * 100):3d}%  {msg}')

    summary, _report = cd.run_chain(args.input, params, out_dir,
                                    logger, _prog)
    logger.close()
    _say(f"CRISPR 设计完成: {summary['n_sets_output']} 组合"
          f"（候选 {summary['n_sets_total']}，最优评分 "
          f"{summary['best_score']}，耗时 {summary['elapsed_sec']}s）")
    for seg, d in summary['segments'].items():
        _say(f"  {seg:12s}: {d['n_isolates']} 株 / 共识 {d['consensus_len']}bp"
              f" / PAM {d['n_pam_sites']} → 漏斗后 {d['n_after_funnel']}"
              f" → 组合 {d['n_designed']}")
    _say(f"产物目录: {summary['out_dir']}")
    _say("⚠ " + summary['evidence_note'])
    return None


def cmd_importprobe(args):
    """内部命令：在**子进程**里逐个导入模块并打印协议行。

    供 selfcheck 在冻结分发下探测「import 直接段错误」的模块 —— 段错误
    会杀死子进程，但父进程靠已收到的行数就能定位肇事模块，平台自身
    （web 服务）不再被一个坏模块拖死。协议：
        P <模块名>          开始导入
        F <模块名> <repr>   导入抛了异常
        E <模块名>          枚举期就导不进来的子包
        DONE <n> <n_fail>   全部跑完
    """
    import importlib
    from Virus_Platform_Core.selfcheck import iter_platform_modules, FALLBACK_MODULES
    mods, walk_errs = iter_platform_modules()
    if not mods:
        mods = FALLBACK_MODULES
    wanted = [m.strip() for m in (args.modules or '').split(',') if m.strip()]
    if wanted:
        # 指定名单=按名单原样导入（不与平台模块表求交）：primer3_runtime
        # 的子进程探测会传 `--cli importprobe primer3` 这类**第三方**模块名，
        # 平台模块表里没有它，求交会得到空名单 → 什么都没导、exit 0，
        # 探测就假阳性通过了。
        mods = wanted
        walk_errs = [w for w in walk_errs if w in wanted]
    out = sys.stdout
    for n in walk_errs:
        out.write(f'E {n}.*\n'); out.flush()
    n_fail = 0
    for m in mods:
        out.write(f'P {m}\n'); out.flush()
        try:
            importlib.import_module(m)
        except Exception as e:
            n_fail += 1
            out.write(f'F {m} {e!r}\n'); out.flush()
    out.write(f'DONE {len(mods)} {n_fail}\n'); out.flush()
    return True


def cmd_selfcheck(args):
    """环境自检：模块导入 / 外部工具 / 数据库 / 磁盘空间，一次摸底。"""
    import importlib
    import time as _t
    ok_all = True

    # 1) 模块导入（遍历 Virus_Platform_Core 下**全部**子模块；
    #    枚举实现在 Virus_Platform_Core.selfcheck，web API 共用同一份）
    _say('── ① 模块导入 ──')
    t0 = _t.time()
    mods, walk_errs = _iter_platform_modules()
    if not mods:                       # 遍历不可用时的兜底，别让自检失去意义
        mods = FALLBACK_MODULES
    fails = []
    for m in mods:
        try:
            importlib.import_module(m)
        except Exception as e:
            fails.append((m, repr(e)))
    for n in walk_errs:                # 遍历阶段就导不进来的子包（如缺依赖）
        fails.append((n + '.*', '子包导入失败，其下模块未能枚举'))
    total = len(mods) + len(walk_errs)
    _say(f'  {"✔" if not fails else "✘"} {total - len(fails)}/{total} '
          f'({_t.time() - t0:.1f}s)')
    for m, e in fails:
        _say(f'    ✘ {m}: {e}')
    ok_all &= not fails

    # 2) 外部工具
    _say('── ② 外部工具 ──')
    cfg = get_config()
    missing = [k for k, v in sorted(cfg.tool_status().items()) if not v]
    total = len(cfg.tool_status())
    _say(f'  {"✔" if not missing else "△"} {total - len(missing)}/{total} 可用'
          + (f'（缺: {", ".join(missing)}）' if missing else ''))
    if missing:
        _say('    缺失工具只影响对应步骤（可选步骤自动降级）；'

              '必装件 kunpeng/SPAdes/BLAST 缺失时对应管道无法运行。')

    # 3) 数据库
    _say('── ③ 数据库 ──')
    from Virus_Platform_Core.kunpeng import db_ready
    from Virus_Platform_Core.taxonomy import taxonomy_ready
    from Virus_Platform_Core.config import host_db_info
    # 分发口径（2026-10-24）：宿主库按物种而异、任何安装包都不随包——
    # legacy 固定槽位 host/classify 上"没有库"是全新部署的**预期状态**，
    # 判 ✘ 会让分发版永远过不了自检。降级为 △ 提示；非 legacy 槽位的
    # 真实库不完整仍判失败。与 selfcheck._check_databases 同口径。
    _h0 = host_db_info(cfg.databases['host']) or {}
    _host_warn = (not db_ready(cfg.databases['host']) and _h0.get('legacy')
                  and not _h0.get('ready') and not _h0.get('has_manifest'))
    checks = [
        ('NCBI Taxonomy', taxonomy_ready(), False),
        ('宿主库 host_db', bool(db_ready(cfg.databases['host'])
                                or _host_warn), _host_warn),
        ('病毒库 virus_db', db_ready(cfg.databases['virus']), False),
    ]
    for name, ok, warn in checks:
        _say(f'  {"✔" if ok else "✘"} {name}'
              + ('（未配置——分发版默认不带，属预期；到「数据库构建」页自建，'
                 '或「📂 加载已构建宿主库」）' if warn else ''))
        if not warn:
            ok_all &= ok
        if name.startswith('NCBI Taxonomy') and not ok:
            # 单独版不含数据库包：✘ 之外必须给出口，否则买家只看到一个红叉
            _say('      提示: 「数据库构建」页可一键下载；离线时把 '
                  'new_taxdump.tar.gz 放到 databases/ 下，'
                  '任一条命令首次运行会自动解压')
        if name.startswith('宿主库'):
            # 宿主库"是谁"必须可见：库本体不记录物种，历史上出过
            # "目录名是枸杞、库却按 4081 番茄建"的事故，用错宿主是静默的。
            h = host_db_info(cfg.databases['host'])
            if h:
                who = (f"{h.get('species')} (taxid={h.get('taxid')})"
                       if h.get('taxid') else '身份未记录（无 host_db.json 清单）')
                _say(f'      目录 = {h.get("name")} | 物种 = {who}')
                if h.get('conflicted'):
                    _say(f'      ⚠ 元数据冲突: seqid2taxid.map 含多个 taxid '
                          f'{h.get("taxids_in_map")}（历史建库污染）')
                if h.get('legacy'):
                    _say('      提示: 这是旧布局固定槽位 host/classify；'
                          '重建后会自动迁到 host-db/<物种>_host_db/')

    # 4) 磁盘空间
    _say('── ④ 磁盘空间 ──')
    import shutil
    for key, d in (('平台根', PLATFORM_ROOT), ('结果', DIRS['results'])):
        try:
            free = shutil.disk_usage(d).free
            _say(f'  {"✔" if free > 10e9 else "△"} {key}: 剩余 {free / 1e9:.1f} GB'
                  + ('' if free > 10e9 else '（建议保留 >10GB）'))
        except OSError as e:
            _say(f'  ✘ {key}: {e}')

    # 5) 结果目录里的样品
    _say('── ⑤ 样品 ──')
    try:
        samples = sorted(os.listdir(DIRS['results']))
        samples = [s for s in samples
                   if os.path.isdir(os.path.join(DIRS['results'], s))]
        _lst = ", ".join(samples[:8]) + ('…' if len(samples) > 8 else '')
        _say(f'  共 {len(samples)} 个样品: {_lst}')
    except OSError:
        _say('  （results/ 目录尚不存在）')

    # 6) 第三方依赖（对照 requirements.txt；缺必需依赖会影响打包分发）
    _say('── ⑥ 第三方依赖 ──')
    try:
        from dev_tools.audit_deps import analyze as _dep_analyze
        dep = _dep_analyze()
        if not dep['third']:
            # 冻结分发：audit_deps 靠遍历源码树收集 import，而 exe 里没有 .py
            # 源码，扫描结果恒为 0。此处若照常打 ✔，就是一句没有信息量的假阳性
            # （2026-09-11 实测：打包版扫到 0 个却报"全部就绪"，而当时 polars
            #  正因漏声明被 spec 排除、整个「已知病毒识别与定量」静默失效）。
            # 如实说明不适用，别给假信心。
            _say('  – 不适用：冻结分发无源码树可扫描'
                  '（依赖已在打包环节固化，请以 ① 模块导入 为准）')
        elif not dep['hard_missing']:
            _say(f"  ✔ 全部依赖就绪（扫描 {len(dep['third'])} 个）")
            # uncovered = 代码在用、但 requirements.txt 没声明。原先这里无条件
            # 打"requirements 已覆盖"，而 uncovered 从不展示 —— 正是这个盲区
            # 让 polars 漏声明了很久，最终演变成打包版整段功能静默失效
            # （2026-09-11）。如实列出来。
            if dep['uncovered']:
                _unc = dep['uncovered']
                _say(f"  △ requirements.txt 未声明 {len(_unc)} 个: "
                      + ', '.join(m for m, _pip in _unc))
            if dep['missing']:
                _say(f"  △ 可选依赖缺失（对应功能降级）: "
                      f"{', '.join(dep['missing'])}")
        else:
            for m in dep['hard_missing']:
                _say(f'  ✘ 必需依赖未装: {m}')
            _say('  安装: python -m pip install -r requirements.txt')
        ok_all &= not dep['hard_missing']
    except Exception as e:
        _say(f'  ✘ 依赖审计失败: {e}')

    _say('════════════════════════')
    _say('✔ 自检通过' if ok_all else '△ 存在未就绪项（见上），按需处理后重跑本命令')
    return ok_all


def build_parser(cfg):
    """构建 CLI 参数解析器。

    2026-10-26：从 main() 内联迁出——CLI 中英双语需要在 parse_args 之前
    就地改写 description/help（i18n_cli.localize_parser），测试也需要
    直接拿到解析器对象做锚点比对。
    """
    _frozen = bool(getattr(sys, 'frozen', False))
    _epilog = __doc__.replace('python main.py', _self_cmd())
    if _frozen:
        _epilog = _epilog.replace('双击 启动平台-桌面窗口.bat',
                                  '双击 VirusPlatform.exe')
    p = argparse.ArgumentParser(
        description="植物病毒分析平台 CLI",
        prog='VirusPlatform.exe --cli' if _frozen else None,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=_epilog)
    p.add_argument('--lang', choices=['zh', 'en'], default=None,
                   help='输出语言：zh=中文，en=English（默认跟随「设置 → 语言」，也可用环境变量 VP_CLI_LANG）')
    sub = p.add_subparsers(dest='cmd', required=True)

    sp = sub.add_parser('init-taxonomy', help='下载/准备 NCBI taxonomy（首次必做）')
    sp.add_argument('--force', action='store_true')
    sp.set_defaults(func=cmd_init_taxonomy)

    sp = sub.add_parser('build-host-db', help='kunpeng 宿主数据库构建')
    sp.add_argument('--genome', required=True, help='宿主基因组 FASTA')
    sp.add_argument('--taxid', type=int, required=True, help='宿主 NCBI TaxID')
    sp.add_argument('--hash-capacity', default='256M')
    sp.add_argument('--threads', type=int, default=0,
                    help='线程数（默认自动=8；内存不足时用 4）')
    sp.add_argument('--out-dir', default=None,
                    help='输出库目录（默认按物种命名 '
                         'host-db/<taxid>_<源目录名>_host_db）')
    sp.add_argument('--clean-mid', action='store_true',
                    help='建库后清理中间文件（library/prep，约省 3GB；'
                         '重建时会自动再生成）')
    sp.add_argument('--rebuild', action='store_true')
    sp.set_defaults(func=cmd_build_host_db)

    sp = sub.add_parser('build-virus-db', help='kunpeng 病毒数据库构建')
    sp.add_argument('--fasta', required=True, help='病毒参考 FASTA')
    sp.add_argument('--info', required=True, help='Accession/Taxid 信息表 TSV')
    sp.add_argument('--hash-capacity', default='64M')
    sp.add_argument('--threads', type=int, default=cfg.threads)
    sp.add_argument('--rebuild', action='store_true')
    sp.set_defaults(func=cmd_build_virus_db)

    sp = sub.add_parser('analyze', help='样品全流程分析')
    sp.add_argument('--r1', required=True)
    sp.add_argument('--r2', default=None)
    sp.add_argument('--sample', required=True)
    sp.add_argument('--stages',
                    default=','.join(DEFAULT_ANALYZE_STAGES),
                    help='阶段（逗号分隔）: '
                         + ','.join(DEFAULT_ANALYZE_STAGES))
    sp.add_argument('--db-host', default=cfg.databases['host'])
    sp.add_argument('--db-virus', default=cfg.databases['virus'])
    sp.add_argument('--threads', type=int, default=cfg.threads)
    sp.add_argument('--confidence', type=float, default=0.0)
    sp.add_argument('--assembly-mode', default='rnaviral',
                    choices=['rnaviral', 'metaviral', 'meta', 'rna', 'isolate'])
    sp.add_argument('--subsample', type=int, default=0,
                    help='子采样 reads 对数，0=全部。mode=prefix 时取前 N 对'
                         '（快速冒烟，可能偏向 lane 开头）；mode=random 时'
                         '蓄水池等概率随机抽 N 对（代表性抽检）')
    sp.add_argument('--subsample-mode', default='prefix',
                    choices=['prefix', 'random', 'proportion'],
                    help='子采样模式: prefix=取前 N 对(默认); random=随机抽'
                         ' N 对(等概率,seed 可复现); proportion=按 --subsample-p'
                         ' 比例随机抽(两遍法,条数精确,对齐 seqkit sample2 -p)')
    sp.add_argument('--subsample-p', type=float, default=0.0,
                    dest='subsample_p',
                    help='比例子采样比例，仅 --subsample-mode proportion 时'
                         '生效（0~1，如 0.5=随机抽 50%%）')
    sp.add_argument('--subsample-seed', type=int, default=42,
                    dest='subsample_seed',
                    help='random/proportion 子采样的随机种子（同 seed 完全可复现）')
    sp.add_argument('--min-contig-len', type=int, default=200)
    sp.add_argument('--top-n-refs', type=int, default=10, help='进化分析取近缘参考数')
    sp.add_argument('--tree-tool', default='fasttree', choices=['fasttree', 'raxml-ng', 'nj'])
    sp.add_argument('--tree-sampling', default='blast',
                    choices=['blast', 'macro', 'genus', 'lineage'],
                    help='参考挑选: blast=按比对hits(默认); macro=同科建树'
                         '(目标属+同科各属背景); genus=属级树; lineage=种级树'
                         '（后三者需分类元数据；旧 acvirus_db/virus_ref 已移除，'
                         '现用 tree_db/plant_tree.db 植物病毒参考库）')
    sp.add_argument('--ncbi-refs', default=None,
                    help='NCBI 参考集合名（ncbi-dl 下载，逗号分隔可多个），'
                         '追加进进化树比对')
    sp.add_argument('--primer-mode', default='conserved', choices=['conserved', 'plain'])
    sp.add_argument('--specificity', action='store_true', help='引物宿主特异性 BLAST 检查(需建宿主BLAST库)')
    sp.add_argument('--min-orf-aa', type=int, default=100)
    sp.add_argument('--memory', type=int, default=64, help='SPAdes 内存上限 GB')
    sp.add_argument('--no-fq2fa', action='store_true',
                    help='跳过 FASTQ→FASTA 预转换（默认开启，推荐保留）')
    sp.add_argument('--gbdraw-max', type=int, default=12,
                    help='基因组图出图 contigs 上限（默认 12）')
    sp.add_argument('--gbdraw-fasta', default=None,
                    help='gbdraw 自备 FASTA（默认用③的病毒 contigs）')
    sp.add_argument('--gbdraw-ann', default=None,
                    help='gbdraw 自备注释（GFF3 或 GenBank .gb/.gbk）')
    sp.add_argument('--plot-engine', default='auto',
                    choices=['auto', 'gbdraw', 'dfv'],
                    help='⑨基因组图引擎：auto=gbdraw 优先、缺则用 '
                         'dna_features_viewer（默认 auto）')
    sp.add_argument('--force', action='store_true')
    sp.set_defaults(func=cmd_analyze)

    sp = sub.add_parser('ncbi-dl', help='NCBI Entrez 批量下载参考序列（进化树扩充）')
    sp.add_argument('term', help='Entrez 检索式，如 "Tobamovirus[ORGN] AND '
                                 'complete genome[TITL]"')
    sp.add_argument('-n', '--name', required=True, help='集合名（ analyze --ncbi-refs 引用）')
    sp.add_argument('--db', default='nucleotide', choices=['nucleotide', 'protein'])
    sp.add_argument('--max', type=int, default=100, help='最多下载条数（默认 100）')
    sp.set_defaults(func=cmd_ncbi_dl)

    sp = sub.add_parser('ncbi-list', help='列出已下载的 NCBI 参考集合')
    sp.set_defaults(func=cmd_ncbi_list)

    sp = sub.add_parser(
        'crispr-design',
        help='CRISPR/Cas12a 检测试剂设计（crRNA × RPA/PCR 引物协同，'
             '与 GUI「检测与防治→CRISPR 检测试剂设计」同引擎）')
    sp.add_argument('-i', '--input', required=True,
                    help='目标 FASTA（多分离物；header 含 RNA1/S2 等标签自动分片段）')
    sp.add_argument('--near', default=None,
                    help='近缘病毒排除 FASTA（可选）')
    sp.add_argument('--no-host', action='store_true',
                    help='跳过宿主基因组筛查（默认启用 host-db）')
    sp.add_argument('--no-group', action='store_true',
                    help='不按片段标签分组（全部记录视为同一分子）')
    sp.add_argument('--amp', default='RPA', choices=['RPA', 'PCR'],
                    help='扩增方式（默认 RPA）')
    sp.add_argument('--profile', default='strict',
                    choices=['strict', 'relaxed'],
                    help='参数档（relaxed=高变异 RNA 病毒放松口径）')
    sp.add_argument('--nuclease', default='cas12a',
                    choices=['cas12a', 'cas13a'],
                    help='效应酶（cas12a=LbCas12a TTTV；cas13a=LwaCas13a '
                         '免 PAM/T7 体外转录，默认 cas12a）')
    sp.add_argument('--spacer-len', type=int, default=0,
                    help='crRNA spacer 长度（0=按酶档案默认：Cas12a 20 / '
                         'Cas13a 30；区间 20–24 / 27–30）')
    sp.add_argument('--probe', default='none', choices=['none', 'exo', 'nfo'],
                    help='RPA 探针（exo=荧光 FAM/THF/BHQ1；nfo=侧流 '
                         'FAM/THF/C3 + 配套生物素引物；默认 none）')
    sp.add_argument('--top-n', type=int, default=20, help='输出组合数（默认 20）')
    sp.add_argument('--flank', type=int, default=250, help='位点侧翼窗口 bp（默认 250）')
    sp.add_argument('--margin', type=int, default=25,
                    help='扩增子含位点的最小边距 bp（默认 25）')
    sp.add_argument('--product-min', type=int, default=0,
                    help='产物下限 bp（0=用档位默认）')
    sp.add_argument('--product-max', type=int, default=0,
                    help='产物上限 bp（0=用档位默认）')
    sp.add_argument('-o', '--out', default=None,
                    help='产物根目录（默认 run/tool_runs/crispr_cli_<时间戳>）')
    sp.add_argument('-t', '--threads', type=int, default=0)
    sp.set_defaults(func=cmd_crispr_design)

    sp = sub.add_parser('gb-dl',
                        help='下载 GenBank 集合（同属共线性比较输入）')
    sp.add_argument('-n', '--name', required=True, help='集合名')
    sp.add_argument('--term', default=None,
                    help='Entrez 检索式，如 "Potyvirus[ORGN] AND '
                         'complete genome[TITL]"（与 --acc 二选一）')
    sp.add_argument('--acc', default=None,
                    help='accession 列表（逗号/空格分隔，支持版本号）')
    sp.add_argument('--max', type=int, default=50, help='最多下载数（检索式时）')
    sp.set_defaults(func=cmd_gb_dl)

    sp = sub.add_parser('gb-import',
                        help='导入本机 GenBank 文件为集合（多记录自动拆分）')
    sp.add_argument('-n', '--name', required=True, help='集合名')
    sp.add_argument('--files', required=True,
                    help='.gb/.gbk 文件路径（逗号分隔，支持 .gz）')
    sp.set_defaults(func=cmd_gb_import)

    sp = sub.add_parser('gb-list', help='列出 GenBank 集合')
    sp.set_defaults(func=cmd_gb_list)

    sp = sub.add_parser('gb-check',
                        help='巡检 GenBank 集合（记录数/CDS 数/警告）')
    sp.add_argument('-n', '--name', required=True, help='集合名')
    sp.set_defaults(func=cmd_gb_check)

    sp = sub.add_parser('report', help='对已有结果重新生成可视化报告')
    sp.add_argument('--sample', required=True)
    sp.add_argument('--force', action='store_true')
    sp.set_defaults(func=cmd_report)

    sp = sub.add_parser('host-analysis',
                        help='对已有 ③组装 结果做 ICTV 宿主预测')
    sp.add_argument('--sample', required=True)
    sp.add_argument('--force', action='store_true')
    sp.add_argument('--min-confidence', type=float, default=None,
                    help='ICTV 级联置信度门槛（0-1 小数，默认不过滤）；'
                         'Integrated_Confidence/P_Max 低于门槛的预测在 '
                         'confidence_flag 列标注 Low_confidence')
    sp.set_defaults(func=cmd_host_analysis)


    sp = sub.add_parser('orfa', help='ORF 功能注释（⑥b，RefSeq 病毒蛋白搜索）')
    sp.add_argument('--sample', required=True)
    sp.add_argument('--threads', type=int, default=None)
    sp.add_argument('--libs', default=None,
                    help='启用的注释层，逗号分隔：prot（层1 序列同源）/ '
                         'pfam（层2 Pfam HMM）/ bfvd（层2 BFVD 病毒域，'
                         '与 NCBI 提交同源）；缺省=全部启用。例：--libs prot,pfam')
    sp.add_argument('--min-pident', type=float, default=None,
                    help='层1 命中 identity 门槛（百分数 0-100，默认不过滤）；'
                         'DIAMOND 侧默认只有 e≤1e-5 一道门槛，收紧可减少'
                         '低 identity 误配（经验起点 30）')
    sp.add_argument('--min-qcov', type=float, default=None,
                    help='层1 命中 query 覆盖度门槛（百分数 0-100，默认不过滤；'
                         '经验起点 50）')
    sp.add_argument('--force', action='store_true')
    sp.set_defaults(func=cmd_orfa)

    sp = sub.add_parser('logan-create',
                        help='LOGAN 溯源：生成待提交 Logan-Search 的查询片段')
    sp.add_argument('--name', required=True, help='查询名称')
    sp.add_argument('--sample', default=None,
                    help='来源样品（需已有 ③组装 病毒 contigs）')
    sp.add_argument('--contigs', default=None,
                    help='逗号分隔的 contig 名（默认全部）')
    sp.add_argument('--fasta', default=None, help='或直接粘贴 FASTA 文本')
    sp.add_argument('--segments', type=int, default=2,
                    help='每条 contig 切片段数 1-4（默认 2；≤2.5kb 恒为 1 段）')
    sp.set_defaults(func=cmd_logan_create)

    sp = sub.add_parser('logan-import',
                        help='LOGAN 溯源：导入 Logan-Search 结果表并生成报告')
    sp.add_argument('--name', required=True, help='查询名称')
    sp.add_argument('--segment', type=int, default=1, help='片段编号（从 1 起）')
    sp.add_argument('--result', required=True, help='结果表 CSV/TSV 文件')
    sp.set_defaults(func=cmd_logan_import)

    sp = sub.add_parser('logan-jobs', help='LOGAN 溯源：列出查询任务')
    sp.set_defaults(func=cmd_logan_jobs)

    sp = sub.add_parser('logan-batch',
                        help='LOGAN 溯源：Selenium 批量提交全部未导入片段并自动收结果')
    sp.add_argument('--name', required=True, help='查询名称')
    sp.add_argument('--email', required=True,
                    help='通知邮箱，逗号分隔多个则轮换')
    sp.add_argument('--group', default='Fast_No_human',
                    choices=['All', 'All_No_viral_human', 'Fast', 'Fast_No_human',
                             'Fast_No_RefSeq', 'Transcriptomic',
                             'Metatranscriptomic', 'Metagenomic', 'GenBank_RefSeq'])
    sp.add_argument('--show-browser', action='store_true',
                    help='显示浏览器窗口（默认后台 headless）')
    sp.add_argument('--first-wait', type=int, default=300)
    sp.add_argument('--max-wait', type=int, default=1800)
    sp.set_defaults(func=cmd_logan_batch)

    sp = sub.add_parser('submit-list', help='提交准备：列出提交项目')
    sp.set_defaults(func=cmd_submit_list)

    sp = sub.add_parser('submit-init', help='提交准备：新建提交项目（unified_metadata.csv）')
    sp.add_argument('--name', required=True, help='项目名（字母/数字/_/-）')
    sp.add_argument('--taxonomy', help='taxonomy.tsv（contig<TAB>taxonomy）→ 自动生成初表')
    sp.add_argument('--metadata', help='公共元数据表（Core14/Full，自动填日期/地点等）')
    sp.add_argument('--import-csv', help='从已有 unified CSV/TSV/Excel 导入')
    sp.add_argument('--demo', action='store_true', help='用示例数据建表')
    sp.add_argument('--authors', help='作者 "Last, First; ..."')
    sp.add_argument('--title', help='提交标题')
    sp.add_argument('--bioproject', help='BioProject ID')
    sp.add_argument('--host', help='宿主物种名')
    sp.add_argument('--lat-lon', help='经纬度 "38.47 N 106.27 E"')
    sp.add_argument('--sequencer', default='Illumina NovaSeq 6000')
    sp.add_argument('--assembler', default='SPAdes;4.3.0;rnaviral')
    sp.add_argument('--coverage', help='覆盖度（如 42.5x）')
    sp.set_defaults(func=cmd_submit_init)

    sp = sub.add_parser('submit-validate',
                        help='提交准备：必填字段校验（存在未解决问题时退出码 1）')
    sp.add_argument('--name', required=True)
    sp.set_defaults(func=cmd_submit_validate)

    sp = sub.add_parser('submit-fill', help='提交准备：批量填充某列占位符/替换旧值')
    sp.add_argument('--name', required=True)
    sp.add_argument('--column', required=True)
    sp.add_argument('--value', required=True, help='新值')
    sp.add_argument('--old', default=None, help='旧值（缺省=填充所有占位符）')
    sp.set_defaults(func=cmd_submit_fill)

    sp = sub.add_parser('submit-export', help='提交准备：生成 source.src/miuvig/assembly/BioSample/report')
    sp.add_argument('--name', required=True)
    sp.add_argument('--assembler', default='SPAdes;4.3.0;rnaviral')
    sp.add_argument('--sequencer', default='Illumina NovaSeq 6000')
    sp.add_argument('--enrichment', default='rRNA depletion')
    sp.set_defaults(func=cmd_submit_export)

    sp = sub.add_parser('submit-adopt',
                        help='提交准备：从 MMPV-RNA discovery 产物一键建提交项目')
    sp.add_argument('--out-dir',
                    help='discovery 输出根（含 09b_Analysis_Verify/ 或 08_Rescue/）')
    sp.add_argument('--from-server', metavar='DATASET',
                    help='从服务器拉取 discovery 产物（相对 server-data-root 的'
                         '数据集名，如 out10）；ssh/scp 免密需已配好')
    sp.add_argument('--server-host', help='ssh 目标（默认 zhangwenda@202.119.189.246）')
    sp.add_argument('--server-data-root',
                    help='服务器数据根（默认 /home/zhangwenda/data-test）')
    sp.add_argument('--name', required=True, help='提交项目名（字母/数字/_/-）')
    sp.add_argument('--dataset', help='数据集名（默认取输出目录名），作 Isolate 前缀')
    sp.add_argument('--samples',
                    help='样本元数据表（TSV/CSV，需 match 列=contig 前缀/glob）')
    sp.add_argument('--run-name', help='平台分类运行名（默认 contigs_mmpv_<数据集>）')
    sp.add_argument('--force', action='store_true', help='项目已存在时清空重建')
    sp.add_argument('--dry-run', action='store_true', help='只报会生成什么，不落盘')
    sp.add_argument('--no-export', action='store_true', help='只出 unified_metadata.csv')
    sp.add_argument('--min-length', type=int, default=0, help='过滤短于该长度的序列')
    sp.add_argument('--isolate-prefix', help='Isolate 前缀（默认用 --dataset）')
    sp.add_argument('--public-metadata',
                    help='公共元数据表（默认自动探测上游 metadata_output/）')

    g = sp.add_argument_group('提交元数据（单样本或作为多样本的兜底）')
    g.add_argument('--geo-loc', help='采集地点 "China:Ningxia"')
    g.add_argument('--lat-lon', help='经纬度 "37.48 N 105.68 E"')
    g.add_argument('--collection-date', help='采集日期 YYYY[-MM[-DD]]')
    g.add_argument('--host', help='宿主物种名')
    g.add_argument('--tissue', help='组织类型')
    g.add_argument('--cultivar', help='栽培品种')
    g.add_argument('--dev-stage', help='发育阶段')
    g.add_argument('--collected-by', help='采集人/机构')
    g.add_argument('--isolation-source', help='分离来源')
    g.add_argument('--authors', help='作者 "Last, First; ..."')
    g.add_argument('--title', help='提交标题')
    g.add_argument('--bioproject', help='你自己的 BioProject（PRJNA...）')
    g.add_argument('--biosample-prefix', help='BioSample 前缀（SAMN...）')
    g.add_argument('--assembler', help='组装方法（缺省从 pipeline_config 读）')
    g.add_argument('--sequencer', help='测序平台')
    g.add_argument('--coverage', help='基因组覆盖度（如 42.5x）')
    g.add_argument('--annotation-pipeline', help='注释流程串')
    g.add_argument('--enrichment', default='rRNA depletion', help='miuvig 富集方式')
    sp.set_defaults(func=cmd_submit_adopt)

    sp = sub.add_parser('submit-sbt', help='提交准备：生成 template.sbt（作者/机构）')
    sp.add_argument('--name', required=True)
    sp.add_argument('--last', required=True)
    sp.add_argument('--first', required=True)
    sp.add_argument('--middle', default='')
    sp.add_argument('--affil', required=True, help='机构')
    sp.add_argument('--div', default='', help='院系')
    sp.add_argument('--city', required=True)
    sp.add_argument('--sub', default='', help='省/州')
    sp.add_argument('--country', required=True)
    sp.add_argument('--street', default='')
    sp.add_argument('--email', required=True)
    sp.add_argument('--postal', default='')
    sp.add_argument('--title', default='', help='提交标题')
    sp.add_argument('--extra-authors', default='',
                    help='其他作者，每行一个 "Last, First"（\\n 分隔）')
    sp.set_defaults(func=cmd_submit_sbt)

    sp = sub.add_parser('meta-search', help='公共数据检索：SRA+GSA 双引擎按物种检索 Run')
    sp.add_argument('--species', required=True, help='物种拉丁名（如 "Lycium barbarum"）')
    sp.add_argument('--source', default='TRANSCRIPTOMIC',
                    help='测序类型过滤 TRANSCRIPTOMIC/GENOMIC/All（默认 TRANSCRIPTOMIC）')
    sp.add_argument('--db', default='both', choices=['sra', 'gsa', 'both'])
    sp.add_argument('--out', default=None, help='输出目录（默认 meta_search/<物种>/search）')
    sp.add_argument('--no-detailed', action='store_true', help='关闭详细模式（Tissue/Location 深提取）')
    sp.add_argument('--ncbi-api', default=None, help='NCBI API Key（也可用环境变量 NCBI_API_KEY）')
    sp.add_argument('--deepseek-api', default=None, help='DeepSeek API Key（AI 元数据清洗，可选）')
    sp.set_defaults(func=cmd_meta_search)

    sp = sub.add_parser('meta-info', help='公共数据检索：Run 列表 → Core14/Full 统一元数据')
    sp.add_argument('--runs', required=True, help='Run 编号列表文件（每行一个 SRR/CRR…）或单个编号')
    sp.add_argument('--out', default=None, help='输出目录（默认 meta_search/info/）')
    sp.add_argument('--mode', default='both', choices=['local', 'api', 'both'])
    sp.add_argument('-t', '--threads', type=int, default=4)
    sp.add_argument('--fill-date', action='store_true', help='采集日期缺失时用发布日期兜底')
    sp.add_argument('--ncbi-api', default=None)
    sp.add_argument('--deepseek-api', default=None)
    sp.set_defaults(func=cmd_meta_info)

    sp = sub.add_parser('meta-plot', help='公共数据检索：元数据 SCI 可视化（时间/机构/组织/地理）')
    sp.add_argument('--input', required=True, help='检索结果或 Core14/Full 元数据 CSV')
    sp.add_argument('--out', default='SCI_Figures_Output', help='图表输出目录')
    sp.set_defaults(func=cmd_meta_plot)

    sp = sub.add_parser('host-genome', help='宿主参考基因组下载（NCBI / gget·Ensembl，供建宿主库）')
    sp.add_argument('--species', default=None, help='物种拉丁名（如 "Lycium barbarum"；与 --accession 二选一）')
    sp.add_argument('--accession', default=None, help='组装 accession（GCF_/GCA_，datasets 精确下载）')
    sp.add_argument('--chromosomes', default=None, help='染色体过滤（如 X,Y；仅 --accession 模式）')
    sp.add_argument('--source', choices=('ncbi', 'gget'), default='ncbi',
                    help='下载来源：ncbi=RefSeq（默认）；gget=Ensembl/Ensembl Plants')
    sp.add_argument('--out', default=None, help='输出目录（默认 host-db/<taxid>_<物种>/）')
    sp.add_argument('--release', type=int, default=None, help='gget: Ensembl release（默认最新）')
    sp.add_argument('--taxid', type=int, default=None, help='指定 NCBI TaxID（省去在线查询）')
    sp.add_argument('--include-organelles', action='store_true', help='同时下载叶绿体/线粒体基因组（仅 NCBI 源）')
    sp.add_argument('--ncbi-api', default=None)
    sp.add_argument('--min-length', type=int, default=0)
    sp.add_argument('--verify-only', action='store_true', help='仅检查已有下载')
    sp.add_argument('--skip-datasets', action='store_true', help='跳过 datasets CLI，直接 E-utilities 回退')
    sp.set_defaults(func=cmd_host_genome)

    sp = sub.add_parser('datasets',
                        help='NCBI datasets CLI 下载：datasets genome|gene|taxonomy|virus …（--help 看各模式参数）')
    sp.add_argument('ds_args', nargs=argparse.REMAINDER,
                    help='透传参数，如 genome --accession GCF_000001405.40 --chromosomes X,Y')
    sp.set_defaults(func=cmd_datasets)

    sp = sub.add_parser('datasets-setup',
                        help='安装/更新 NCBI datasets CLI（官方 exe → 3rd/tools/datasets/）')
    sp.set_defaults(func=cmd_datasets_setup)

    sp = sub.add_parser('ref-status', help='查看病毒参考库版本与统计')
    sp.set_defaults(func=cmd_ref_status)

    sp = sub.add_parser('ictv-update',
                        help='更新 ICTV VMR 参考库（在线下载或本地 xlsx → taxa.txt）')
    sp.add_argument('--xlsx', default=None,
                    help='本地 VMR xlsx 路径（不给则从 ictv.global 下载当前版）')
    sp.set_defaults(func=cmd_ictv_update)

    sp = sub.add_parser('ictv-status', help='查看 ICTV VMR 参考库状态')
    sp.set_defaults(func=cmd_ictv_status)

    sp = sub.add_parser('ictv-refs',
                        help='ICTV VMR 选参：按属/科/种挑参考，缺的按需下载')
    sp.add_argument('--genus', default=None, help='属名（如 Tobamovirus）')
    sp.add_argument('--family', default=None, help='科名（如 Potyviridae）')
    sp.add_argument('--db', default='plant', choices=['plant', 'ictv'],
                    help='建树参考库：plant=植物口径 / ictv=全病毒界（默认 plant）')
    sp.add_argument('--species', default=None, help='种名（如 "Tobamovirus mosaic"）')
    sp.add_argument('--limit', type=int, default=20, help='最多挑多少条（默认 20）')
    sp.add_argument('--any', action='store_true',
                    help='不限定 Complete genome（含 Coding-complete 等）')
    sp.add_argument('--download', action='store_true',
                    help='对本地没有的 accession 执行 NCBI 按需下载')
    sp.set_defaults(func=cmd_ictv_refs)

    sp = sub.add_parser('tools', help='显示工具探测状态')
    sp.set_defaults(func=cmd_tools)

    sp = sub.add_parser('db-migrate',
                        help='数据库迁移/对接：软件与数据库分离部署')
    sp.add_argument('--to', required=True,
                    help='目标数据库根目录（绝对路径，平台目录之外）')
    sp.add_argument('--mode', default='copy', choices=['copy', 'move'],
                    help="copy=复制后切换（推荐，源保留）；move=复制校验后删源")
    sp.add_argument('--check', action='store_true',
                    help='只校验目标目录是否已有完整数据库（对接前检查）')
    sp.add_argument('--dry-run', action='store_true',
                    help='只做预检（空间/路径），不实际复制')
    sp.set_defaults(func=cmd_db_migrate)

    sp = sub.add_parser('tool-runs',
                        help='输出目录管理：status/organize/archive/clean')
    sp.add_argument('action', choices=['status', 'organize', 'archive', 'clean'],
                    help='status=现状统计；organize=杂项归位；archive=归档旧运行；'
                         'clean=删除超期运行/归档（不可逆，先 --dry-run）')
    sp.add_argument('--before', default=None,
                    help="archive: 归档该日期（YYYYMMDD）之前的运行")
    sp.add_argument('--older-days', type=int, default=None,
                    help='archive: 归档早于 N 天的运行')
    sp.add_argument('--active-days', type=int, default=None,
                    help='clean: 删除活动区早于 N 天的运行')
    sp.add_argument('--archive-days', type=int, default=None,
                    help='clean: 删除归档区早于 N 天的归档')
    sp.add_argument('--dry-run', action='store_true',
                    help='只预演，不实际移动/删除')
    sp.set_defaults(func=cmd_tool_runs)

    sp = sub.add_parser('samples-backfill',
                        help='给缺 project.json 的老样品从既有产物回填清单')
    sp.add_argument('--sample', action='append', default=None,
                    help='只处理指定样品（可重复；支持非规范命名），缺省处理全部')
    sp.add_argument('--dry-run', action='store_true',
                    help='只报告将要写入的内容，不落盘')
    sp.set_defaults(func=cmd_samples_backfill)

    sp = sub.add_parser('selfcheck', help='环境自检（模块/工具/数据库/磁盘/依赖）')
    sp.set_defaults(func=cmd_selfcheck)

    sp = sub.add_parser('importprobe', help='内部命令：逐模块导入探测（供 selfcheck 子进程调用）')
    sp.add_argument('modules', help='逗号分隔的模块名清单')
    sp.set_defaults(func=cmd_importprobe)
    return p


def main():
    cfg = get_config()
    # 退出前把 AI 审计队列排空：CLI 跑管道同样会发布事件（note_stage_*），
    # 而 Windows 下 SIGTERM 不触发 atexit（实测）。信号处理器只能在主线程装。
    try:
        from Virus_Platform_Core.ai_agent import audit as _ai_audit
        _ai_audit.install_signals()
    except Exception:
        pass
    # 用法横幅随运行形态变化（冻结分发版没有 python / main.py），见 build_parser。
    p = build_parser(cfg)
    try:
        from Virus_Platform_Core import i18n_cli as _i18n
        _i18n.set_lang(None)          # 先按环境变量/配置定语言
        _i18n.install_print_hook()   # 核心模块里的直接 print 也跟随语言
        _i18n.localize_parser(p)     # help/description 就地英文化（zh 时不动）
    except Exception:
        _i18n = None
    args = p.parse_args()
    if _i18n is not None and getattr(args, "lang", None):
        _i18n.set_lang(args.lang)    # --lang 覆盖本次运行
    # 首次运行把随包的 new_taxdump.tar.gz 解压出来（网页端是后台线程做的，
    # CLI 先跑 selfcheck/analyze 时也必须就绪，否则首条命令就报 Taxonomy ✘）
    if args.cmd != 'init-taxonomy':
        _bootstrap_taxonomy()
    # 顶层异常包装（2026-10-04 审查 E5）：此前 cmd_* 异常裸抛 traceback，
    # 脚本化批处理难以区分「参数/运行失败」与解释器崩溃。统一为
    # 一行可读错误 + 退出码 2（1 保留给 selfcheck 的「存在未就绪项」）。
    try:
        rc = args.func(args)
    except KeyboardInterrupt:
        _say('\n[中断] 用户取消', file=sys.stderr)
        sys.exit(130)
    except SystemExit:
        raise
    except Exception as e:
        _say(f'[错误] {type(e).__name__}: {e}\n'
              f'       （完整回溯: python -c "import traceback, '
              f'sys; sys.argv=...；或加环境变量 VP_CLI_TRACEBACK=1 重跑）',
              file=sys.stderr)
        if os.environ.get('VP_CLI_TRACEBACK'):
            raise
        sys.exit(2)
    # selfcheck 返回 False（存在未就绪项）时以非零码退出，供脚本判断
    if rc is False:
        sys.exit(1)


if __name__ == '__main__':
    main()
