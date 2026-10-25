# -*- coding: utf-8 -*-
"""数据库构建与库状态（自 app.py 拆出）。

Taxonomy / 宿主库 / 病毒库构建、Kraken2 转换、kv 索引、库状态、
文件浏览对话框。"""
import json
import logging
import os
import re
import time

from flask import (Blueprint, abort, jsonify, request)

from Virus_Platform_Core.config import DIRS, PLATFORM_ROOT, db_path
from Virus_Platform_Core.utils import (TaskLogger, check_path, fmt_size)
from Virus_Platform_Core.web.state import cfg
from Virus_Platform_Core.web.tasks import tm

bp = Blueprint('build', __name__)


def _selfcheck_invalidate():
    """库状态刚变化（建库/加载/切库/登记），丢掉自检的进程内缓存。

    不失效的话，总览页环境状态条会把变化前的"数据库 未就绪"继续展示到
    缓存过期（默认 5 分钟），看起来像"库加载了但总览没同步"。"""
    from Virus_Platform_Core.selfcheck import invalidate
    invalidate()


def _busy_build_task(target=None, name=None):
    """同目标（输出目录 / 任务名）的建库任务是否已在跑或排队。

    建库类接口原先无并发去重：对同一目录连点两次会起两个任务同时写
    同一个库目录（一个删旧索引另一个正在建 → 产物混合、互相覆盖）。
    heavy 闸门只管并发数量，不管目标是否重复。提交前用它拦一道 409。
    """
    def _norm(p):
        return os.path.normcase(os.path.normpath(os.path.abspath(str(p))))
    tgt = _norm(target) if target else None
    with tm.lock:
        recs = list(tm.tasks.values())
    for r in recs:
        if r.get('status') not in ('running', 'queued'):
            continue
        if tgt and r.get('out_dir') and _norm(r['out_dir']) == tgt:
            return r
        if name and r.get('name') == name:
            return r
    return None


def _reject_busy_build(target=None, name=None):
    """同目标建库任务在跑 → 409（busy 记录一并告知，便于用户去任务中心取消）。"""
    busy = _busy_build_task(target, name)
    if busy:
        abort(409, f'同一目标的建库任务已在运行（{busy.get("id")} '
                   f'{busy.get("name")}）；请等它完成或在任务中心取消后再提交')


@bp.route('/api/tools')
def api_tools():
    return jsonify({'threads': cfg.threads, 'tools': cfg.tool_status()})


@bp.route('/api/dbs')
def api_dbs():
    from Virus_Platform_Core.kunpeng import db_ready
    out = {}
    for key in ('host', 'virus'):
        d = cfg.databases[key]
        out[key] = {'path': d,
                    'ready': db_ready(d) if os.path.isdir(d) else False}
    # 通用参考库（refvirus / 其它 kunpeng 库目录自动探测）；RVDB 已于 2026-09-10 移除
    for key, mapping in (('refvirus', ('virus', 'ref')),
                         ('k2viral', None)):
        if mapping:
            d = db_path(*mapping)
        else:
            d = os.path.join(DIRS['databases'], 'k2viral_db')
        out[key] = {'path': d,
                    'ready': db_ready(d) if os.path.isdir(d) else False}
    from Virus_Platform_Core.taxonomy import taxonomy_ready
    out['taxonomy'] = {'ready': taxonomy_ready(), 'path': DIRS['taxonomy']}
    # 病毒分类库统一为「自备预构建」：扫描 databases/ 下的 kunpeng 库目录
    # （宿主库除外）。使用者把库目录放进来即自动识别，工具②/④的病毒库
    # 下拉与构建页列表都以此为准。准入只看 db_ready（结构完整可用），
    # 不看 .building 标记——平台意外退出会留下孤儿标记，不应隐藏可用库。
    libs = []
    base = DIRS['databases']
    for dp, dn, fn in os.walk(base):
        if not any(f.startswith('hash_') and f.endswith('.k2d') for f in fn):
            continue
        rel = os.path.relpath(dp, base).replace('\\', '/')
        if rel in ('.', ''):
            continue
        # 宿主分类库不列入病毒库下拉（现已移出 databases/，放 host-db/，
        # 此处保留 host/ 前缀跳过以兼容旧布局残留目录）
        if rel.startswith('host/'):
            continue
        try:
            if not db_ready(dp):
                continue
        except (ValueError, OSError):
            continue
        # path 必须是**绝对路径**：工具②/④ 把下拉值原样回传，后端据此解析。
        # 曾写 'databases/<rel>'（按平台根理解），在「数据库包与程序目录分离」
        # 的分发版下会被解析成 <程序目录>\databases\... —— 那里并不存在，
        # 于是选库即报「病毒库不存在」；而开发机上源码树恰好也有 databases/，
        # 会静默解析到**另一份**库，问题被完全掩盖。
        libs.append({'name': rel, 'path': os.path.normpath(dp)})
    # 外部登记的自备库（platform.json extra_virus_libs，可为平台外任意位置）
    seen = {os.path.normpath(l['path']) for l in libs}
    for p in cfg.extra_virus_libs:
        n = os.path.normpath(os.path.abspath(p))
        if n in seen or not _db_dir_ready(n):
            continue
        libs.append({'name': os.path.basename(n), 'path': p,
                     'external': True})
    out['virus_libs'] = libs
    return jsonify(out)


def _db_dir_ready(d):
    """kunpeng 库结构完整性（不限定平台目录，用于加载/登记自备库）。"""
    if not os.path.isdir(d):
        return False
    if not all(os.path.isfile(os.path.join(d, f))
               for f in ('opts.k2d', 'taxo.k2d', 'hash_config.k2d')):
        return False
    import glob as _glob
    return bool(_glob.glob(os.path.join(d, 'hash_*.k2d')))


def _abs_arg(raw):
    d = raw if os.path.isabs(raw) else os.path.join(PLATFORM_ROOT, raw)
    return os.path.normpath(os.path.abspath(d))


def _within_allowed_roots(d):
    from Virus_Platform_Core.config import write_roots
    d = os.path.normpath(os.path.abspath(d))
    return any(d == os.path.normpath(r)
               or d.startswith(os.path.normpath(r) + os.sep)
               for r in write_roots())


@bp.route('/api/load_taxonomy', methods=['POST'])
def api_load_taxonomy():
    """加载已有 NCBI Taxonomy（含 nodes.dmp/names.dmp 的目录）→
    复制 .dmp/.pkl 进平台标准位置 databases/tax_db（扁平布局）。"""
    body = request.get_json(force=True) or {}
    raw = (body.get('dir') or '').strip()
    if not raw:
        abort(400, '缺少 taxonomy 目录')
    d = _abs_arg(raw)
    if not os.path.isdir(d):
        abort(400, f'目录不存在: {raw}')
    missing = [f for f in ('nodes.dmp', 'names.dmp')
               if not os.path.isfile(os.path.join(d, f))]
    if missing:
        abort(400, f'缺少 {"、".join(missing)}（不是有效的 taxonomy 目录）')
    os.makedirs(DIRS['taxonomy'], exist_ok=True)
    import shutil as _shutil
    copied = []
    for fn in sorted(os.listdir(d)):
        src = os.path.join(d, fn)
        dst = os.path.join(DIRS['taxonomy'], fn)
        if not os.path.isfile(src) or not fn.lower().endswith(('.dmp', '.pkl')):
            continue
        if os.path.normpath(os.path.abspath(src)) == os.path.normpath(dst):
            copied.append(fn)                     # 已在标准位置，跳过
            continue
        _shutil.copyfile(src, dst)
        copied.append(fn)
    from Virus_Platform_Core.taxonomy import taxonomy_ready
    _selfcheck_invalidate()
    return jsonify({'ok': taxonomy_ready(), 'copied': copied,
                    'dir': DIRS['taxonomy']})


@bp.route('/api/load_host_db', methods=['POST'])
def api_load_host_db():
    """加载已构建宿主库（kunpeng）登记为**当前**宿主库（持久化 platform.json）。

    平台外的自备库不再一刀切拒绝：先把它所在目录登记为外部库根
    （config.register_extra_db_root，持久化），后续 db_ready / check_path
    等就绪检查即放行——与「设置 → 数据库目录」同一信任级别与机制。
    走 config.set_active_host_db：host-db 下的库记目录名（可移植），
    外部库记绝对路径；旧的读取方（host_removal / 自检）无需改。"""
    body = request.get_json(force=True) or {}
    raw = (body.get('dir') or '').strip()
    if not raw:
        abort(400, '缺少宿主库目录')
    d = _abs_arg(raw)
    if not _db_dir_ready(d):
        abort(400, '不是有效的 kunpeng 库目录（缺 hash_*.k2d / opts.k2d / '
                   'taxo.k2d / hash_config.k2d）: ' + raw)
    if not _within_allowed_roots(d):
        from Virus_Platform_Core.config import register_extra_db_root
        try:
            register_extra_db_root(os.path.dirname(d))
        except (ValueError, OSError) as e:
            abort(400, f'无法登记外部库根: {e}')
    from Virus_Platform_Core.config import host_db_info
    cfg.set_active_host_db(d)
    _selfcheck_invalidate()
    return jsonify({'ok': True, 'path': cfg.databases['host'],
                    'info': host_db_info(cfg.databases['host'])})


@bp.route('/api/host_dbs')
def api_host_dbs():
    """已建成的宿主库列表（按物种目录）+ 当前生效的那个。

    列表里带 taxid / 物种 / 是否元数据冲突（seqid2taxid.map 含多个 taxid）
    —— 历史上就是靠目录名猜库是谁、结果猜错过。"""
    from Virus_Platform_Core.config import (host_db_dirs, host_db_info,
                                            legacy_host_db_dir,
                                            get_active_host_db)
    items = [host_db_info(p) for p in host_db_dirs()]
    legacy = legacy_host_db_dir()
    if os.path.isdir(legacy) and not any(
            os.path.normcase(i['path']) == os.path.normcase(legacy)
            for i in items):
        items.append(host_db_info(legacy))
    active = get_active_host_db()
    for it in items:
        it['active'] = os.path.normcase(it['path']) == os.path.normcase(active or '')
    return jsonify({'active': active, 'items': items,
                    'configured': cfg.active_host_db})


@bp.route('/api/set_host_db', methods=['POST'])
def api_set_host_db():
    """切换当前宿主库。name = host-db 下的目录名（推荐）或绝对路径；'' = 自动。

    可用性校验在 config.set_active_host_db 里完成（落盘前校验、失败不留痕），
    这里只把 ValueError 转成 400。"""
    from Virus_Platform_Core.config import host_db_info
    body = request.get_json(force=True) or {}
    name = (body.get('name') or '').strip()
    try:
        path = cfg.set_active_host_db(name)
    except ValueError as e:
        abort(400, str(e))
    _selfcheck_invalidate()
    return jsonify({'ok': True, 'path': path, 'info': host_db_info(path)})


@bp.route('/api/register_virus_lib', methods=['POST'])
def api_register_virus_lib():
    """登记 / 移除外部已构建病毒库目录（持久化；分类只读，允许平台外任意位置）。"""
    body = request.get_json(force=True) or {}
    raw = (body.get('dir') or '').strip()
    if not raw:
        abort(400, '缺少病毒库目录')
    d = _abs_arg(raw)
    if body.get('remove'):
        cfg.extra_virus_libs = [p for p in cfg.extra_virus_libs
                                if os.path.normpath(os.path.abspath(p)) != d]
        cfg.save()
        _selfcheck_invalidate()
        return jsonify({'ok': True, 'libs': cfg.extra_virus_libs})
    if not _db_dir_ready(d):
        abort(400, '不是有效的 kunpeng 库目录（缺 hash_*.k2d / opts.k2d / '
                   'taxo.k2d / hash_config.k2d）: ' + raw)
    if d not in cfg.extra_virus_libs:
        cfg.extra_virus_libs.append(d)
        cfg.save()
    _selfcheck_invalidate()
    return jsonify({'ok': True, 'libs': cfg.extra_virus_libs})


@bp.route('/api/browse')
def api_browse():
    """目录浏览（只读）。平台内目录可用相对路径；也支持浏览整台电脑
    （绝对路径，含盘符列表）。选择的数据输入文件可为任意位置；
    分析输出仍严格限制在平台目录内。"""
    raw = (request.args.get('path', '') or '.').strip()
    # all=1（选目录模式）：文件不限扩展名，便于查看目录里有什么
    show_all = request.args.get('all') == '1'
    if raw in ('.', '', '/', '~', '此电脑'):
        drives = [f'{c}:\\' for c in 'CDEFGHIJKLMNOPQRSTUVWXYZ'
                  if os.path.exists(f'{c}:\\')]
        return jsonify({'cwd': '此电脑', 'parent': '', 'dirs': drives,
                        'files': []})
    p = os.path.normpath(os.path.abspath(raw))
    if not os.path.isdir(p):
        abort(400, '目录不存在')
    dirs, files = [], []
    try:
        for name in sorted(os.listdir(p)):
            child = os.path.join(p, name)
            try:
                if os.path.isdir(child):
                    dirs.append(name)
                elif os.path.isfile(child) and (
                        show_all or name.lower().endswith(
                        ('.fastq.gz', '.fq.gz', '.fastq', '.fq',
                         '.fasta', '.fa', '.fna', '.fas',
                         '.fa.gz', '.fasta.gz', '.tsv', '.sra',
                         '.dmp', '.pkl', '.map',
                         '.gb', '.gbk', '.gbff', '.genbank'))):
                    files.append({'name': name,
                                  'size': fmt_size(os.path.getsize(child))})
            except OSError:
                continue                                  # 跳过无权限项
    except PermissionError:
        abort(400, '无权限访问该目录')
    except OSError:
        abort(400, '无法访问该路径')
    parent = os.path.dirname(p)
    return jsonify({'cwd': p,
                    'parent': '' if parent == p else parent,
                    'dirs': dirs, 'files': files})


@bp.route('/api/path_info')
def api_path_info():
    """路径探测（对话框路径栏「粘贴 + 回车」跳转/选中用）。
    只报存在性与类型，不列目录——列目录走 /api/browse。"""
    raw = (request.args.get('path', '') or '').strip()
    if not raw:
        return jsonify({'exists': False})
    p = os.path.normpath(os.path.abspath(raw))
    parent = os.path.dirname(p)
    return jsonify({'path': p, 'exists': os.path.exists(p),
                    'is_dir': os.path.isdir(p),
                    'parent': '' if parent == p else parent})


def _job_build_taxonomy(log, prog, cancel):
    from Virus_Platform_Core.taxonomy import prepare_taxonomy
    logger = TaskLogger(callback=log)
    prog('taxonomy', 0.1, '准备 NCBI taxonomy')
    prepare_taxonomy(logger=logger)
    prog('taxonomy', 1.0, 'taxonomy 就绪')
    _selfcheck_invalidate()
    logger.close()
    return 'ok'


def _job_build_host_db(body):
    def job(log, prog, cancel):
        from Virus_Platform_Core.kunpeng import build_host_db
        logger = TaskLogger(callback=log)
        # 建库内存峰值与线程数成正比（每线程约 1GB 缓冲），默认限 8
        threads = int(body.get('threads') or min(8, cfg.threads))
        prog('build_host', 0.05, f"注入 taxid={body['taxid']} 并建库 (线程 {threads})")
        out = build_host_db(body['genome'], int(body['taxid']),
                            db_dir=(body.get('out_dir') or None),
                            hash_capacity=body.get('hash_capacity', '256M'),
                            threads=threads,
                            logger=logger, rebuild=bool(body.get('rebuild')),
                            clean_mid=bool(body.get('clean_mid')))
        # 建完即设为"当前宿主库"：否则会出现"库已建好、平台却还在用旧库"
        # 这种静默错配（旧行为正是如此——建完不切，得再手工加载一次）。
        try:
            cfg.set_active_host_db(out)
            logger.log(f"当前宿主库已切换: {out}")
        except ValueError as e:                      # 目录没建成才可能到这
            logger.log(f"宿主库已建但切换当前库失败: {e}", "WARN")
        _selfcheck_invalidate()
        prog('build_host', 1.0, f'宿主库就绪: {os.path.basename(out)}')
        logger.close()
        return 'ok'
    return job


def _job_build_virus_db(body):
    def job(log, prog, cancel):
        from Virus_Platform_Core.kunpeng import build_virus_db
        logger = TaskLogger(callback=log)
        prog('build_virus', 0.05, '解析 info 表并建库')
        build_virus_db(body['fasta'], body['info'],
                       db_dir=(body.get('db_dir') or None),
                       hash_capacity=body.get('hash_capacity', '64M'),
                       threads=int(body.get('threads', cfg.threads)),
                       logger=logger, rebuild=bool(body.get('rebuild')),
                       clean_mid=bool(body.get('clean_mid')))
        _selfcheck_invalidate()
        prog('build_virus', 1.0, '病毒库就绪')
        logger.close()
        return 'ok'
    return job


@bp.route('/api/build_taxonomy', methods=['POST'])
def api_build_taxonomy():
    tid = tm.start(cfg.tr('下载/准备 Taxonomy', 'Download/prepare Taxonomy'),
                   _job_build_taxonomy, weight='light')
    return jsonify({'task': tid})


@bp.route('/api/build_host_db', methods=['POST'])
def api_build_host_db():
    body = request.get_json(force=True) or {}
    for k in ('genome', 'taxid'):
        if not body.get(k):
            abort(400, f'缺少参数 {k}')
    check_path(body['genome'], must_exist=True)
    # 输出目录：留空 = 按物种自动命名 host-db/<TaxID>_<源目录名>_host_db/
    # （多宿主可共存；库内会写 host_db.json 清单记录它是谁）
    from Virus_Platform_Core.config import host_db_name
    raw = (body.get('out_dir') or '').strip()
    if raw:
        d = _abs_arg(raw)
        if not _within_allowed_roots(d):
            abort(400, '输出目录须位于平台目录内（或先在「设置 → 数据库目录」'
                       '配置自定义数据库根）')
    else:
        try:
            d = os.path.join(DIRS['host_src'],
                             host_db_name(int(body['taxid']), body['genome']))
        except (TypeError, ValueError):
            abort(400, f"TaxID 必须是数字，收到: {body['taxid']!r}")
    body['out_dir'] = d
    tid = tm.start(cfg.tr(f"宿主库构建 (taxid={body['taxid']})",
                         f"Host DB build (taxid={body['taxid']})"),
                   _job_build_host_db(body))
    return jsonify({'task': tid, 'out_dir': d})


@bp.route('/api/host_genome_download', methods=['POST'])
def api_host_genome_download():
    """宿主参考基因组下载（上游）→ 产物直接衔接「宿主库构建」（下游）。

    body: {species: 物种拉丁名（"Oryza sativa" 或 oryza_sativa）,
           accession: 可选组装 accession（GCF_/GCA_，给出后走 NCBI datasets
                      按号精确下载，species 可省）,
           chromosomes: 可选染色体过滤（仅 accession 模式，如 "X,Y"）,
           source: 'gget'（默认，Ensembl/Ensembl Plants，补 NCBI 覆盖不到
                   的植物宿主）| 'ncbi'（datasets taxon / E-utilities 回退）,
           release: 可选 Ensembl release（仅 gget；留空 = 最新）,
           taxid:   可选 NCBI TaxID（留空则在线查询）}

    产物：host-db/<taxid>_<物种>/genome.fa（平台建库锚点名，连同
    all.genome.uniq.fasta / genome_report.txt / gget_ref.json）。任务
    result 带回 genome_fa 与 taxid，前端（数据资源组「基因组查询下载」
    独立卡）据此显示结果并跳转建库页，实现 下载 → 建库 一条龙。
    """
    body = request.get_json(force=True) or {}
    species = (body.get('species') or '').strip()
    accession = (body.get('accession') or '').strip()
    chromosomes = (body.get('chromosomes') or '').strip()
    if not species and not accession:
        abort(400, '缺少物种拉丁名或 Assembly Accession')
    if accession and not re.fullmatch(r'[GC]C[AF]_\d+\.\d+', accession):
        abort(400, f'accession 格式应为 GCF_/GCA_数字.版本: {accession}')
    source = (body.get('source') or 'gget').strip()
    if source not in ('gget', 'ncbi'):
        abort(400, f'未知下载来源: {source}（可选 gget / ncbi）')
    if source == 'gget' and accession:
        abort(400, 'gget/Ensembl 源不支持 accession（那是 NCBI 组装号）；'
                   '请留空来源或改用 NCBI')
    if source == 'gget':
        import importlib.util
        if importlib.util.find_spec('gget') is None:
            abort(400, 'gget 未安装：请先在平台 Python 环境执行 '
                       'pip install "gget==0.30.5"，或改选 NCBI 源')
    elif accession:
        from Virus_Platform_Core.public_meta.datasets_cli import \
            datasets_available
        if not datasets_available():
            abort(400, 'accession 模式需要 NCBI datasets CLI：'
                       '执行 python main.py datasets-setup 自动安装')
    release = body.get('release') or None
    if release is not None:
        try:
            release = int(release)
        except (TypeError, ValueError):
            abort(400, 'Ensembl release 必须是数字')
    taxid = body.get('taxid') or None
    if taxid is not None:
        try:
            taxid = int(taxid)
        except (TypeError, ValueError):
            abort(400, 'TaxID 必须是数字')

    if source == 'ncbi' and not accession:
        # 先查后下保险（前端已有一道）：物种可能命中多个参考组装版本，
        # 盲下会把多个版本全部合并进同一个宿主库，互相污染。0 个/多个
        # 都不放行，让用户先查询并选定一个版本。
        from Virus_Platform_Core.public_meta.host_genome import \
            list_reference_assemblies
        try:
            hits = list_reference_assemblies(species)
        except Exception as e:                  # noqa: BLE001
            abort(502, f'NCBI 组装查询失败: {e}')
        if not hits:
            abort(400, f'NCBI 未找到 {species} 的参考组装：可改选 '
                       'gget·Ensembl 源，或填写 Assembly Accession '
                       '（GCF_/GCA_）精确下载')
        if len(hits) > 1:
            lst = '、'.join(h['accession'] or h['uid'] for h in hits[:8])
            abort(409, f'{species} 命中 {len(hits)} 个参考组装版本'
                       f'（{lst}）；请先点「🔍 查询 NCBI 组装」选用一个 '
                       'accession 再下载，不要把多个版本混建一个宿主库')

    def job(log, prog, cancel):
        from Virus_Platform_Core.public_meta import host_genome as hg
        from Virus_Platform_Core.config import DIRS
        logger = TaskLogger(callback=log)
        slug = re.sub(r'[^\w\-]+', '_', accession or species).strip('_') \
            or 'host'
        if accession:
            prog('download', 0.1, f'datasets accession {accession}')
            from Virus_Platform_Core.public_meta import datasets_cli as dsc
            res = hg.download_genome_accession(
                accession, os.path.join(DIRS['host_src'], f'genome_{slug}'),
                chromosomes=chromosomes, taxid=taxid, log=logger.log)
            _selfcheck_invalidate()
            prog('done', 1.0, f"基因组就绪: {res['genome_fa']}")
            logger.close()
            return res
        prog('lookup', 0.05, '解析 TaxID 与输出目录')
        t = taxid or hg._lookup_taxid(species)
        out_dir = os.path.join(
            DIRS['host_src'], f'{t}_{slug}' if t else f'genome_{slug}')
        if source == 'gget':
            prog('download', 0.15, f'gget/Ensembl 下载 {species}')
            res = hg.download_genome_gget(species, out_dir, release=release,
                                          taxid=t, log=logger.log)
        else:
            prog('download', 0.15, f'NCBI 下载 {species}')
            os.makedirs(out_dir, exist_ok=True)
            try:
                hg.download_genome_https(species, out_dir, log=logger.log)
            except Exception as e:                 # noqa: BLE001
                logger.log(f'NCBI 下载失败: {e}', 'ERROR')
            fasta_files = hg.collect_fasta_files(
                os.path.join(out_dir, 'extracted'))
            if not fasta_files:
                raise RuntimeError('NCBI 源未下载到任何基因组 FASTA'
                                   '（可改选 gget/Ensembl 源）')
            merged = os.path.join(out_dir, 'all.genome.uniq.fasta')
            hg.merge_and_deduplicate(fasta_files, merged, log=logger.log)
            res = {'out_dir': out_dir,
                   'genome_fa': hg.finalize_host_dir(out_dir, merged,
                                                     log=logger.log),
                   'taxid': t, 'species': hg.normalize_species(species),
                   'source': 'ncbi'}
        _selfcheck_invalidate()
        prog('done', 1.0, f"基因组就绪: {res.get('genome_fa', out_dir)}")
        logger.close()
        return res

    label = accession or species
    _task_name = cfg.tr(f'宿主基因组下载 {label}',
                       f'Host genome download {label}')
    # out_dir 在本接口内要到 job 里查完 taxid 才知道，故按任务名去重
    # （同名 = 同 accession/物种，落同一宿主目录）。
    _reject_busy_build(name=_task_name)
    tid = tm.start(_task_name, job)
    return jsonify({'task': tid})


@bp.route('/api/datasets_summary', methods=['POST'])
def api_datasets_summary():
    """先查后下：datasets summary 同步预览（查询基因组查询下载页「🔍 查询」）。

    按 kind 归一化成前端表格行，直接返回 JSON（不进任务中心——查询是
    秒级只读操作，CLI 侧 300s 超时兜底）。下载仍走各自下载接口，查询条件
    与下载参数一一对应，避免"查的和下的不是一批"。

    body: {kind: 'genome' | 'virus' | 'gene' | 'taxonomy',
           taxon / accession / host                    (genome, virus)
           gene-id / symbol / taxon / protein-accession (gene)
           taxons                                      (taxonomy)
           limit: 可选预览条数（默认 genome 200 / virus 500 / gene 200 /
                  taxonomy 50，超出的记 truncated=true）}
    """
    body = request.get_json(force=True) or {}
    kind = (body.get('kind') or '').strip()
    if kind not in ('genome', 'virus', 'gene', 'taxonomy'):
        abort(400, f'未知查询类型: {kind}（可选 genome / virus / gene / '
                   'taxonomy）')
    from Virus_Platform_Core.public_meta.datasets_cli import \
        datasets_available
    if not datasets_available():
        abort(400, 'NCBI datasets CLI 未安装：执行 python main.py '
                   'datasets-setup 自动安装（官方 exe，NCBI 域名白名单内）')

    try:
        limit = int(body.get('limit'))
    except (TypeError, ValueError):
        limit = None
    kwargs = {
        'taxon': (body.get('taxon') or '').strip() or None,
        'accession': (body.get('accession') or '').strip() or None,
        'host': (body.get('host') or '').strip() or None,
        'gene_id': (body.get('gene-id') or '').strip() or None,
        'symbol': (body.get('symbol') or '').strip() or None,
        'protein_acc': (body.get('protein-accession') or '').strip() or None,
        'taxons': (body.get('taxons') or '').strip() or None,
        'limit': limit,
    }
    from Virus_Platform_Core.public_meta.datasets_cli import query_records
    try:
        res = query_records(kind, **kwargs)
    except ValueError as e:
        abort(400, str(e))
    except RuntimeError as e:
        # CLI 缺失 / NCBI 查询失败 / 超时：502，前端把消息展示在结果区
        abort(502, str(e))
    return jsonify(res)


@bp.route('/api/datasets_virus', methods=['POST'])
def api_datasets_virus():
    """NCBI datasets 下载病毒数据集（taxon + 可选宿主过滤，或单条 accession）→
    直接落平台鉴定库数据集布局（reference.fasta + reference.ref_info.tsv，
    databases/virusref_db/<名>/），建库页「原始数据集」卡即时识别，可一键建库。

    body: {taxon: 病毒 taxon（如 "tomato brown rugose fruit virus"）,
           host: 可选宿主过滤（如 "Solanum lycopersicum"）,
           accession: 可选病毒核酸 accession（如 NC_045512.2，查询表里
                      挑单条下载时用；给定时忽略 taxon）,
           name: 可选数据集名（默认按 taxon/accession 命名）}"""
    body = request.get_json(force=True) or {}
    taxon = (body.get('taxon') or '').strip()
    host = (body.get('host') or '').strip()
    name = (body.get('name') or '').strip()
    accession = (body.get('accession') or '').strip()
    if not taxon and not accession:
        abort(400, '缺少病毒 taxon 或 accession')
    if accession and not re.fullmatch(r'[A-Za-z0-9_.\-]{4,32}', accession):
        abort(400, f'accession 格式可疑: {accession}')
    from Virus_Platform_Core.public_meta.datasets_cli import \
        datasets_available
    if not datasets_available():
        abort(400, 'NCBI datasets CLI 未安装：执行 python main.py '
                   'datasets-setup 自动安装（官方 exe，NCBI 域名白名单内）')
    if name and not re.fullmatch(r'[A-Za-z0-9_\-]+', name):
        abort(400, f'数据集名仅限字母数字-_: {name}')
    from Virus_Platform_Core.config import DIRS
    if name:
        tgt = os.path.join(DIRS['virus_src'], name)
        if os.path.isfile(os.path.join(tgt, 'reference.fasta')):
            abort(409, f'同名数据集已存在: {tgt}；请换一个名字')

    def job(log, prog, cancel):
        from Virus_Platform_Core.public_meta import datasets_cli as dsc
        logger = TaskLogger(callback=log)
        prog('download', 0.1,
             ('datasets virus taxon=' + taxon if taxon else
              'datasets virus accession=' + accession)
             + (f' host={host}' if host else ''))
        res = dsc.download_virus(taxon=taxon or None, host=host or None,
                                 name=name or None, include=('genome',),
                                 accession=accession or None,
                                 log=logger.log)
        prog('done', 1.0,
             f"数据集就绪: {res['reference']}（{res['n_seqs']} 条序列）")
        logger.close()
        return {'out_dir': res['out_dir'], 'reference': res['reference'],
                'ref_info': res['ref_info'], 'n_seqs': res['n_seqs'],
                'name': res['name']}

    _task_name = cfg.tr(f'病毒数据集下载 {taxon or accession}',
                        f'Virus dataset download {taxon or accession}')
    # 给了 name（目标目录已知）按目录去重，否则按任务名（同 taxon+host）。
    _reject_busy_build(tgt if name else None, _task_name)
    tid = tm.start(_task_name, job)
    return jsonify({'task': tid})


@bp.route('/api/datasets_download', methods=['POST'])
def api_datasets_download():
    """NCBI datasets 通用下载（gene / taxonomy / genome 通用模式）→
    datasets_dl/<模式>_<查询名>/，不做平台数据集布局（virus 走上面的
    /api/datasets_virus，宿主基因组走 /api/host_genome_download）。

    body: {mode: 'genome' | 'gene' | 'taxonomy',
           taxon / accession / chromosomes / include / exclude   (genome)
           gene-id / symbol / taxon / protein-accession          (gene)
           taxons                                                (taxonomy)}
    """
    body = request.get_json(force=True) or {}
    mode = (body.get('mode') or '').strip()
    if mode not in ('genome', 'gene', 'taxonomy'):
        abort(400, f'未知模式: {mode}（可选 genome / gene / taxonomy）')
    from Virus_Platform_Core.public_meta.datasets_cli import \
        datasets_available
    if not datasets_available():
        abort(400, 'NCBI datasets CLI 未安装：执行 python main.py '
                   'datasets-setup 自动安装（官方 exe，NCBI 域名白名单内）')

    kwargs = {}
    if mode == 'genome':
        taxon = (body.get('taxon') or '').strip()
        accession = (body.get('accession') or '').strip()
        if not taxon and not accession:
            abort(400, 'genome 模式需要 taxon 或 accession')
        if accession and not re.fullmatch(r'[GC]C[AF]_\d+\.\d+', accession):
            abort(400, f'accession 格式应为 GCF_/GCA_数字.版本: {accession}')
        kwargs = {'taxon': taxon or None, 'accession': accession or None,
                  'chromosomes': (body.get('chromosomes') or '').strip() or None,
                  'include': (body.get('include') or '').strip() or None}
        slug = re.sub(r'[^\w\-]+', '_', accession or taxon).strip('_') or 'ds'
    elif mode == 'gene':
        gene_id = str(body.get('gene-id') or '').strip()
        symbol = (body.get('symbol') or '').strip()
        gt = (body.get('taxon') or '').strip()
        acc = (body.get('protein-accession') or '').strip()
        if not gene_id and not acc and not (symbol and gt):
            abort(400, 'gene 模式需要 gene-id，或 symbol+taxon，'
                       '或 protein-accession')
        kwargs = {'gene_id': gene_id or None, 'symbol': symbol or None,
                  'taxon': gt or None, 'accession': acc or None}
        slug = re.sub(r'[^\w\-]+', '_', gene_id or acc or symbol).strip('_')
    else:
        taxons = (body.get('taxons') or '').strip()
        if not taxons:
            abort(400, 'taxonomy 模式需要 taxons（逗号分隔）')
        kwargs = {'taxons': taxons}
        slug = re.sub(r'[^\w\-]+', '_', taxons.split(',')[0]).strip('_')

    out_dir = os.path.join(PLATFORM_ROOT, 'datasets_dl',
                           f'{mode}_{slug or "ds"}')

    def job(log, prog, cancel):
        from Virus_Platform_Core.public_meta import datasets_cli as dsc
        logger = TaskLogger(callback=log)
        prog('download', 0.1, f'datasets {mode}（见日志）')
        res = dsc.download(mode, out_dir=out_dir, log=logger.log, **kwargs)
        prog('done', 1.0, f"下载完成: {res['out_dir']}")
        logger.close()
        return {'out_dir': res['out_dir'],
                'fasta': res.get('fasta') or [],
                'n_seqs': res.get('n_seqs'),
                'n_reports': len(res.get('reports') or [])}

    tid = tm.start(cfg.tr(f'datasets {mode} 下载 {slug}',
                         f'datasets {mode} download {slug}'), job)
    return jsonify({'task': tid, 'out_dir': out_dir})


def _universal_db_job(source):
    """通用病毒库建库 job 构造器。"""
    def job(log, prog, cancel):
        from Virus_Platform_Core.universal_ref import build_universal_db
        logger = TaskLogger(callback=log)
        prog('prep', 0.05, '解析 accession→taxid 映射')
        prog('build', 0.3, 'kunpeng add-library + build-db（见日志）')
        res = build_universal_db(
            source, hash_capacity='2G', logger=logger, rebuild=True)
        prog('build', 0.9, '建库完成，校验库文件')
        _selfcheck_invalidate()
        prog('done', 1.0, '完成')
        logger.close()
        return res
    return job


@bp.route('/api/convert_kraken2', methods=['POST'])
def api_convert_kraken2():
    """Kraken2 库包/目录 → kunpeng 分片库（kunpeng hashshard，方式 C）。

    body: {tar: "databases/k2_viral_20260626.tar.gz"（平台内路径）,
           name: "k2viral"（目标库目录名 databases/<name>）,
           hash_capacity: "1G",
           overwrite: 目标已是可用库时是否允许覆盖重建（缺省拒绝）}
    """
    body = request.get_json(force=True) or {}
    tar = (body.get('tar') or '').strip()
    name = (body.get('name') or 'k2viral').strip()
    if not tar:
        abort(400, '缺少 Kraken2 库包路径')
    check_path(tar, must_exist=True)
    if not re.fullmatch(r'[A-Za-z0-9_\-]+', name):
        abort(400, '库名仅限字母数字-_')
    db_dir = os.path.join(DIRS['databases'], name)
    # 防覆盖（与鉴定库建库同一口径）：目标已是结构完整的 kunpeng 库时，
    # 不显式勾 overwrite 就拒绝——转换写到一半的残目录不算可用库，放行。
    # 判定用 falsy（与下方 build_virus_db 的 `not body.get('overwrite')` 一致）：
    # 原先写 `is None`，前端传 overwrite=false / "" 时会被当成"要覆盖"放行。
    if not body.get('overwrite') and _db_dir_ready(db_dir):
        abort(409, f'目标库已存在且可用: {db_dir}；请换一个库名，'
                   '或请求里带 overwrite=true 覆盖重建')

    def job(log, prog, cancel):
        from Virus_Platform_Core.kunpeng import convert_kraken2
        logger = TaskLogger(callback=log)
        prog('convert', 0.2, '解包 + hashshard 转换（见日志）')
        res = convert_kraken2(tar, db_dir,
                              hash_capacity=body.get('hash_capacity') or '1G',
                              logger=logger)
        _selfcheck_invalidate()
        prog('done', 1.0, '完成')
        logger.close()
        return {'db_dir': res}

    _task_name = cfg.tr(f'Kraken2 库转换 {name}', f'Kraken2 convert {name}')
    _reject_busy_build(db_dir, _task_name)       # 同目标已有任务在跑 → 409
    tid = tm.start(_task_name, job, out_dir=db_dir)
    return jsonify({'task': tid})


def _walk_kv_libs(root):
    """root 下递归找鉴定库目录：含 reference.fasta 或 manifest.json 即命中。

    命中后**不再下钻**——库是自包含的（库内 backup_*/logs/salmon_k31 等
    子目录不是独立库），这样自动识别任意位置的库，又不把库内备份当库。
    隐藏目录（. 开头）跳过。
    """
    found = []
    if not os.path.isdir(root):
        return found
    stack = [root]
    while stack:
        cur = stack.pop()
        try:
            names = sorted(os.listdir(cur))
        except OSError:
            continue
        if 'reference.fasta' in names or 'manifest.json' in names:
            lib = _kv_index_probe(cur, os.path.basename(cur))
            if lib:
                lib['path'] = cur
                found.append(lib)
                continue                     # 剪枝：库内不再找"库"
        for n in names:
            p = os.path.join(cur, n)
            if not n.startswith('.') and os.path.isdir(p):
                stack.append(p)
    return found


def _discover_kv_libs():
    """自动识别全部鉴定库（不写死位置）：三级来源，路径去重。

    ① 数据库整根（DIRS['databases']，含自定义数据库根）递归——
       virusref_db/<库名> 只是其中最常见位置，任何子目录里的库都算；
    ② 外部库根（设置页登记的 extra_db_roots）整根递归；
    ③ 外部登记目录（extra_kv_indexes，逐个）。
    scope：platform = databases/virusref_db 下（可按库名引用，建库页落点）；
           internal = 数据库根下其他位置；external = 平台外（只读使用）。
    """
    from Virus_Platform_Core.config import write_roots
    out, seen = [], set()

    def _add(lib, scope):
        p = os.path.normcase(os.path.normpath(os.path.abspath(lib['path'])))
        if p in seen:
            return
        seen.add(p)
        lib['scope'] = scope
        out.append(lib)

    def _under(path, base):
        path = os.path.normcase(os.path.normpath(os.path.abspath(path)))
        base = os.path.normcase(os.path.normpath(os.path.abspath(base)))
        return path == base or path.startswith(base + os.sep)

    vsrc = DIRS['virus_src']
    for root in [DIRS['databases']] + [
            r for r in (getattr(cfg, 'extra_db_roots', []) or [])
            if not _under(r, DIRS['databases'])]:
        base_scope_root = DIRS['databases']
        for lib in _walk_kv_libs(root):
            scope = 'platform' if _under(lib['path'], vsrc) else (
                'internal' if _under(lib['path'], base_scope_root)
                else 'external')
            _add(lib, scope)
    for p in (getattr(cfg, 'extra_kv_indexes', []) or []):
        if not os.path.isdir(p):
            continue
        lib = _kv_index_probe(p)
        if lib:
            _add(lib, 'external')
    # 默认库 kv_index 即使缺参考/索引也保留条目（便于排查），但不翻转 scope
    if not any(l['name'] == 'kv_index' and l['scope'] == 'platform'
               for l in out):
        d = os.path.join(vsrc, 'kv_index')
        if os.path.isdir(d):
            lib = _kv_index_probe(d, 'kv_index') or {
                'name': 'kv_index', 'path': d, 'salmon': False,
                'minibwa': False, 'manifest': None}
            lib['path'] = d
            _add(lib, 'platform')
    return out


@bp.route('/api/datasets')
def api_datasets():
    """原始数据集与派生库的就绪态（建库页「原始数据集」卡）。

    口径（数据集与数据库区分开）：
      * 数据集 = 分发携带的原始参考（reference.fasta + ref_info.tsv）、
        RefSeq/RVDB 源 FASTA —— 随包分发、只读使用；
      * 派生库 = 引擎索引（salmon/minibwa，首跑自动现建）、kunpeng 分类库
        （构建页一键构建）—— 本机派生物，删了可随时由数据集重建。
    """
    from Virus_Platform_Core.kunpeng import db_ready
    out = {'virusref': [], 'universal': {}}
    for lib in _discover_kv_libs():
        d = lib['path']
        nm = lib['name']
        ref = os.path.join(d, 'reference.fasta')
        if not os.path.isfile(ref):
            continue                         # 无参考数据集的目录不进数据集卡
        # 同源 kunpeng 分类库：平台库（virusref_db 下）才有约定落点——
        # kv_index ↔ 主库 kunpeng_db/virus_db，其余同名对应 kunpeng_db/<name>；
        # 其他位置的库不预设分类库落点。
        if lib.get('scope') == 'platform':
            if nm in ('kv_index', 'plant_ref'):   # 平台默认库（2026-09-29 更名）
                kdir, kname = db_path('virus', 'virus_db'), 'virus_db'
            else:
                kdir = os.path.join(DIRS['databases'], 'kunpeng_db', nm)
                kname = nm
            kunpeng = {'name': kname, 'path': kdir,
                       'ready': bool(db_ready(kdir))}
        else:
            kunpeng = {'name': '', 'path': '', 'ready': False}
        ri = os.path.join(d, 'reference.ref_info.tsv')
        out['virusref'].append({
            'name': nm, 'path': d,
            'reference': ref,
            'ref_mb': round(os.path.getsize(ref) / 1e6, 1),
            'ref_info': ri if os.path.isfile(ri) else '',
            'salmon': bool(lib.get('salmon')),
            'minibwa': bool(lib.get('minibwa')),
            'kunpeng_name': kunpeng['name'],
            'kunpeng_path': kunpeng['path'],
            'kunpeng_ready': kunpeng['ready'],
        })
    uni = out['universal']
    try:
        from Virus_Platform_Core.universal_ref import _ref_fasta
        for src, key, libkey in (('refvirus', 'refvirus', 'ref'),
                                 ('rvdb', 'rvdb', 'rvdb')):
            try:
                fa = _ref_fasta(src)
            except (OSError, ValueError):
                fa = None
            kd = db_path('virus', libkey)
            uni[key] = {'fasta': fa or '',
                        'built': bool(db_ready(kd)),
                        'built_path': kd}
    except Exception as exc:                   # noqa: BLE001
        uni['error'] = str(exc)
    return jsonify(out)


@bp.route('/api/build_kv_index', methods=['POST'])
def api_build_kv_index():
    """构建「病毒鉴定库」：参考 FASTA → 各定量引擎比对索引。

    body: {fasta: 参考 FASTA（平台内外均可，只读）,
           name: 库名（默认 kv_index，限字母数字-_）,
           threads: 线程数,
           ref_info: 可选，参考注释 TSV，与 FASTA 一起归档进库目录,
           engines: 可选，要建的索引引擎列表（默认 ['salmon']；可勾
                    'minibwa'，两套索引互不影响、可分次补建）}

    库与引擎正交：库 = 参考真相（reference.fasta + 注释 + manifest）+
    各引擎派生索引子目录（salmon_k31/、minibwa/reference.*）。选库即选
    参考+已有索引整体；定量时选哪个引擎（--engine）与用哪个库互不干涉。
    产物：databases/virusref_db/<name>/ 下被勾选引擎的索引子目录与
    reference.fasta + reference.ref_info.tsv（自包含库目录：库带着自己的
    参考，默认参考解析、t-kvsuite / kvsuite 索引复用都以此为准），并写
    一份 manifest.json 记录来源与时间。与「病毒分类库」（kunpeng）无关，
    两者互不影响。旧布局的 virusref_db 根目录散置 final.cluster.ref.*
    仍被读取方兼容，但新库一律写自包含布局。
    """
    body = request.get_json(force=True) or {}
    fasta = (body.get('fasta') or '').strip()
    if not fasta:
        abort(400, '缺少参数 fasta')
    fasta = check_path(fasta if os.path.isabs(fasta)
                       else os.path.join(PLATFORM_ROOT, fasta),
                       must_exist=True)
    if not os.path.isfile(fasta):
        abort(400, f'参考 FASTA 不是文件：{fasta}')

    engines = body.get('engines') or ['salmon']
    if not isinstance(engines, list):
        abort(400, 'engines 必须是列表')
    engines = list(dict.fromkeys(
        str(e).strip().lower() for e in engines if str(e).strip()))
    bad = [e for e in engines if e not in ('salmon', 'minibwa')]
    if bad:
        abort(400, f'不支持的建库引擎: {", ".join(bad)}（可选 salmon / minibwa）')
    if not engines:
        engines = ['salmon']
    name = (body.get('name') or 'kv_index').strip()
    if not re.fullmatch(r'[A-Za-z0-9_\-]+', name):
        abort(400, '库名仅限字母数字-_')
    out_dir = check_path(os.path.join(DIRS['virus_src'], name),
                         must_exist=False, in_platform=True)

    # 同名库冲突保护（用户口径：别和现有库冲突）：
    #   * 目录不存在 / 空目录 / 只有日志 → 不算冲突，直接建；
    #   * 已是有效库（有 reference.fasta 或任一引擎索引）且未勾 overwrite
    #     → 409 拒绝，提示换名或勾「允许覆盖同名库」；
    #   * 勾了 overwrite → 建库前清掉旧参考与旧索引，绝不允许出现
    #     「新参考 + 旧索引」的错配库（引擎的"存在即复用"只会认文件存在，
    #     不校验内容一致，必须在这里把旧产物清干净）。
    overwrite = bool(body.get('overwrite'))
    _stale = [os.path.join(out_dir, p) for p in
              ('reference.fasta', 'salmon_k31', 'salmon', 'minibwa')]
    if any(os.path.exists(p) for p in _stale) and not overwrite:
        abort(409, f'同名鉴定库已存在: {out_dir}；请换一个库名，'
                   '或勾选「允许覆盖同名库」重建')

    ref_info = (body.get('ref_info') or '').strip()
    if ref_info:
        ref_info = check_path(ref_info if os.path.isabs(ref_info)
                              else os.path.join(PLATFORM_ROOT, ref_info),
                              must_exist=True)

    try:
        threads = max(1, min(int(body.get('threads') or cfg.threads),
                             (os.cpu_count() or 4) * 4))
    except (TypeError, ValueError):
        threads = cfg.threads

    def job(log, prog, cancel):
        from ..known_virus_suite.kv_common import ToolRegistry, setup_logger
        from ..known_virus_suite.kv_engines import make_engine
        import shutil as _shutil
        os.makedirs(out_dir, exist_ok=True)
        # 旧库清理放进 job 体内（提交时只做 409 判定）：原实现在排队前就
        # 删旧参考/索引——任务还在排队时用户已失去可用旧库，任务失败更
        # 是「旧的没了、新的没建成」。任务开始后才删，且只在勾了 overwrite
        # 时删；若已丢失（被别的任务/用户删掉）也不报错。
        if overwrite:
            for p in _stale:
                try:
                    if os.path.isdir(p):
                        _shutil.rmtree(p)
                    elif os.path.isfile(p):
                        os.remove(p)
                except OSError as exc:               # noqa: BLE001
                    raise RuntimeError(f'清理旧库失败（{p}）：{exc}')
            log('已清理旧库产物（overwrite=true）')
        # kv_engines 的 logger 需是标准 logging.Logger（.info/.warning/.error）。
        # tag 带时间戳避免多次建库复用同一 logger 时 handler 累积。
        logger, _logs = setup_logger(
            out_dir, tag='kv_index_build_' + time.strftime('%H%M%S'))

        class _PanelHandler(logging.Handler):
            """把 kv_engines 的日志推到任务面板"""
            def emit(self, rec):
                try:
                    log(f'{time.strftime("%H:%M:%S", time.localtime(rec.created))} '
                        f'{rec.levelname} - {rec.getMessage()}')
                except Exception:                # noqa: BLE001
                    pass

        logger.addHandler(_PanelHandler())
        reg = ToolRegistry(logger)
        reg.probe()
        logger.info(f"参考 FASTA: {fasta}")
        logger.info(f"引擎: {', '.join(engines)}  线程: {threads}")
        logger.info(f"产物目录: {out_dir}")

        made, failed = [], []
        for i, eng_name in enumerate(engines):
            if cancel is not None and cancel.is_set():
                raise RuntimeError('用户取消')
            prog(eng_name, 0.05 + 0.9 * i / len(engines), f'构建 {eng_name} 索引')
            if not reg.has(eng_name):
                failed.append(f'{eng_name}: 未找到可执行文件')
                logger.error(f'{eng_name} 不可用，跳过')
                continue
            try:
                engine = make_engine(eng_name, reg, threads=threads,
                                     logger=logger)
                # 各引擎的 build_index 落在 out_dir/<引擎子目录>
                idx = engine.build_index(fasta, out_dir, threads)
                made.append({'engine': eng_name, 'index': str(idx)})
                logger.info(f'{eng_name} 索引就绪 -> {idx}')
            except Exception as exc:                # noqa: BLE001
                failed.append(f'{eng_name}: {exc}')
                logger.error(f'{eng_name} 建索引失败：{exc}')

        if not made:
            raise RuntimeError('所有引擎建索引均失败：' + '; '.join(failed))

        # 参考与注释归档进库目录（自包含：换机器/整库拷走即含全部输入）
        try:
            _shutil.copyfile(fasta, os.path.join(out_dir, 'reference.fasta'))
            logger.info(f'参考归档 -> {os.path.join(out_dir, "reference.fasta")}')
        except OSError as exc:
            logger.warning(f'参考归档失败（默认参考解析将回退原位置）：{exc}')
        if ref_info:
            try:
                _shutil.copyfile(ref_info,
                                 os.path.join(out_dir, 'reference.ref_info.tsv'))
                logger.info('注释归档 -> reference.ref_info.tsv')
            except OSError as exc:
                logger.warning(f'注释归档失败：{exc}')

        manifest = {
            'name': name,
            'reference': fasta,
            'ref_info': ref_info or '',
            'engines': made,
            'failed': failed,
            'threads': threads,
            'built_at': time.strftime('%Y-%m-%d %H:%M:%S'),
        }
        try:
            with open(os.path.join(out_dir, 'manifest.json'), 'w',
                      encoding='utf-8') as fh:
                json.dump(manifest, fh, ensure_ascii=False, indent=2)
        except OSError as exc:
            logger.warning(f'manifest 写入失败：{exc}')
        # 关闭 FileHandler，否则日志文件句柄会泄漏并锁住产物目录
        for h in list(logger.handlers):
            try:
                h.close()
            except Exception:                    # noqa: BLE001
                pass
            logger.removeHandler(h)
        prog('done', 1.0, '建库完成')
        return {'out': out_dir, 'engines': [m['engine'] for m in made],
                'failed': failed}

    _task_name = cfg.tr(f'病毒鉴定库构建 {name}', f'Virus index build {name}')
    _reject_busy_build(out_dir, _task_name)      # 同目标已有任务在跑 → 409
    tid = tm.start(_task_name, job, out_dir=out_dir)
    return jsonify({'task': tid, 'out': out_dir})


def _kv_index_probe(d, name=None):
    """探测单个目录是否是一个可用的病毒鉴定库（自包含参考 + 任一引擎索引）。

    返回 lib 字典（含各引擎索引的有无与路径、库内参考），不可用时返回 None。
    仅探测结构，不判断目录所在位置，平台内与平台外目录共用。
    库与引擎正交：salmon/minibwa 索引各自独立探测，有一套即可用；
    旧布局里散置的 minibwa 索引文件不识别（索引只认 <库>/minibwa/ 子目录）。
    """
    lib = {'name': name or os.path.basename(d), 'path': d, 'salmon': False,
           'minibwa': False, 'manifest': None}
    for cand in (os.path.join(d, 'salmon_k31', 'info.json'),
                 os.path.join(d, 'salmon', 'salmon_k31', 'info.json')):
        if os.path.isfile(cand):
            lib['salmon'] = True
            lib['salmon_path'] = os.path.dirname(cand)
            break
    mbw = os.path.join(d, 'minibwa', 'reference.mbw')
    if os.path.isfile(mbw):
        lib['minibwa'] = True
        lib['minibwa_path'] = os.path.dirname(mbw)
    for cand in ('reference.fasta', 'manifest.json'):
        mf = os.path.join(d, cand)
        if os.path.isfile(mf):
            if cand == 'reference.fasta':
                lib['reference'] = mf
                continue
            try:
                with open(mf, encoding='utf-8') as fh:
                    lib['manifest'] = json.load(fh)
            except (OSError, ValueError):
                pass
    # 库有效性的锚点是参考真相（reference.fasta）；引擎索引是派生物，
    # 可随库分发，也可在首次运行时按需自建（见 kv_stage.index_dir_for），
    # 缺失不构成"无效库"——salmon k31 全库现建约 10s，k15 约 35s。
    return lib if (lib['salmon'] or lib['minibwa']
                   or lib.get('reference')) else None


@bp.route('/api/kv_index_list')
def api_kv_index_list():
    """列出可用的鉴定库（自动识别，不写死位置）与引擎可用性。

    发现范围（_discover_kv_libs）：数据库整根递归（含自定义数据库根，
    virusref_db/<库名> 只是常见位置之一）+ 外部库根递归 + 外部登记目录。
    libs 条目带 scope：platform（可按库名引用）/ internal（数据库根内其他
    位置，按绝对路径引用）/ external（平台外，只读，按绝对路径引用）。
    """
    out = {'base': DIRS['virus_src'], 'libs': [], 'engines': {}}
    try:
        from ..known_virus_suite.kv_common import ToolRegistry
        reg = ToolRegistry()
        reg.probe()
        out['engines']['salmon'] = {'available': reg.has('salmon'),
                                    'path': str(reg.paths.get('salmon') or '')}
        out['engines']['minibwa'] = {'available': reg.has('minibwa'),
                                     'path': str(reg.paths.get('minibwa') or '')}
    except Exception as exc:                    # noqa: BLE001
        out['engines_error'] = str(exc)
    for lib in _discover_kv_libs():
        lib = dict(lib)
        lib['external'] = lib['scope'] == 'external'
        lib.setdefault('reference', '')
        out['libs'].append(lib)
    out['extra'] = [p for p in (getattr(cfg, 'extra_kv_indexes', []) or [])]
    return jsonify(out)


@bp.route('/api/kv_index_register', methods=['POST'])
def api_kv_index_register():
    """登记 / 移除外部病毒鉴定库目录（持久化；只读使用，允许平台外任意位置）。

    与 /api/register_virus_lib 同语义：把平台外已建好的 salmon 索引目录
    挂进平台，之后在工具卡的「鉴定库」下拉里即可选择复用。
    """
    body = request.get_json(force=True) or {}
    raw = (body.get('dir') or '').strip()
    if not raw:
        abort(400, '缺少鉴定库目录')
    d = _abs_arg(raw)
    if body.get('remove'):
        cfg.extra_kv_indexes = [
            p for p in (getattr(cfg, 'extra_kv_indexes', []) or [])
            if os.path.normpath(os.path.abspath(p)) != d]
        cfg.save()
        return jsonify({'ok': True, 'dirs': cfg.extra_kv_indexes})
    if not os.path.isdir(d):
        abort(400, '目录不存在: ' + raw)
    lib = _kv_index_probe(d)
    if lib is None:
        abort(400, '不是有效的鉴定库目录（未找到 reference.fasta 或任何'
                   '引擎索引）: ' + raw)
    if d not in (getattr(cfg, 'extra_kv_indexes', []) or []):
        cfg.extra_kv_indexes = list(getattr(cfg, 'extra_kv_indexes', []) or []) + [d]
        cfg.save()
    return jsonify({'ok': True, 'dirs': cfg.extra_kv_indexes,
                    'lib': {'name': lib['name'], 'path': d,
                            'salmon': lib['salmon']}})


@bp.route('/api/build_universal_db', methods=['POST'])
def api_build_universal_db():
    """构建通用病毒参考库（refvirus=NCBI RefSeq Viral / rvdb=RVDB C-RVDB）。"""
    body = request.get_json(force=True) or {}
    source = body.get('source') or 'refvirus'
    if source not in ('refvirus', 'rvdb'):
        abort(400, 'source 需为 refvirus 或 rvdb')
    name = {'refvirus': '通用病毒库构建 (RefSeq Viral)',
            'rvdb': 'RVDB 库构建 (C-RVDB)'}[source]
    tid = tm.start(cfg.tr(name, name), _universal_db_job(source))
    return jsonify({'task': tid})


@bp.route('/api/build_virus_db', methods=['POST'])
def api_build_virus_db():
    """从参考 FASTA + info 表构建 kunpeng 病毒分类库（高级入口）。

    body: {fasta, info,
           name: 库名（默认 plant = 平台主库，兼容旧行为）,
           rebuild/clean_mid/hash_capacity/threads,
           overwrite: 目标已是可用库时是否允许覆盖重建（缺省拒绝）}

    防覆盖（用户口径「建库不要覆盖了原有的」）：build_db 本身是替换式
    重建，直接写主库会把旧库清掉；这里在提交前拦一道——目标已是结构
    完整的库且未显式 overwrite 时 409 拒绝，换名即可另建新库共存
    （databases/ 下的库会被自动扫描进工具②/④的病毒库下拉）。
    """
    body = request.get_json(force=True) or {}
    for k in ('fasta', 'info'):
        if not body.get(k):
            abort(400, f'缺少参数 {k}')
    check_path(body['fasta'], must_exist=True)
    check_path(body['info'], must_exist=True)
    name = (body.get('name') or 'plant').strip()
    if not re.fullmatch(r'[A-Za-z0-9_\-]+', name):
        abort(400, '库名仅限字母数字-_')
    if name == 'plant':
        db_dir = db_path('virus', 'plant')
    else:
        db_dir = os.path.join(DIRS['databases'], 'kunpeng_db', name)
    if _db_dir_ready(db_dir) and not body.get('overwrite'):
        abort(409, f'目标库已存在且可用: {db_dir}；请换一个库名，'
                   '或请求里带 overwrite=true 覆盖重建')
    body['db_dir'] = db_dir

    def job(log, prog, cancel):
        from Virus_Platform_Core.kunpeng import build_virus_db
        logger = TaskLogger(callback=log)
        prog('build_virus', 0.05, '解析 info 表并建库')
        build_virus_db(body['fasta'], body['info'], db_dir=db_dir,
                       hash_capacity=body.get('hash_capacity', '64M'),
                       threads=int(body.get('threads', cfg.threads)),
                       logger=logger, rebuild=bool(body.get('rebuild')),
                       clean_mid=bool(body.get('clean_mid')))
        _selfcheck_invalidate()
        prog('build_virus', 1.0, '病毒库就绪')
        logger.close()
        return 'ok'

    _task_name = cfg.tr(f'病毒库构建 {name}', f'Virus DB build {name}')
    _reject_busy_build(db_dir, _task_name)       # 同目标已有任务在跑 → 409
    tid = tm.start(_task_name, _job_build_virus_db(body), out_dir=db_dir)
    return jsonify({'task': tid, 'db_dir': db_dir})
