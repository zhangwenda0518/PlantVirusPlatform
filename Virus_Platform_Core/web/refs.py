# -*- coding: utf-8 -*-
"""参考序列与比较分析 API（自 app.py 拆出）。

NCBI Entrez 下载、ICTV 级联选参、GenBank 集合管理、共线性比较、
MSA / 进化树 / SDT 数据接口。"""
import json
import os
import re
import shutil

from flask import (Blueprint, abort, jsonify, request,
                   send_file)

from Virus_Platform_Core.config import DIRS, PLATFORM_ROOT
from Virus_Platform_Core.utils import (TaskLogger, check_path, safe_open)
from Virus_Platform_Core.web.common import tree_overlay_for
from Virus_Platform_Core.web.state import cfg
from Virus_Platform_Core.web.tasks import tm

bp = Blueprint('refs', __name__)


@bp.route('/api/ncbi/collections')
def api_ncbi_collections():
    from Virus_Platform_Core.ncbi_download import list_collections
    return jsonify([c for c in list_collections()
                    if not c.get('name', '').startswith('_')])


def _pick_db(value):
    """归一建树参考库名：plant / ictv（缺省 plant）。

    库本身决定宿主口径：plant = 植物口径，ictv = 全病毒界不过滤。
    """
    from Virus_Platform_Core import acvirus
    db = str(value or '').strip().lower() or acvirus.DEFAULT_DB
    if db not in acvirus._LIBS:
        abort(400, f'未知建树库: {db}（可选 {"/".join(acvirus._LIBS)}）')
    if not acvirus.available(db):
        abort(400, f'{db} 建树库未就绪（缺少 fasta 或 taxa）')
    return db


@bp.route('/api/ictv/cascade')
def api_ictv_cascade():
    """级联选参下拉数据：?Realm=..&Phylum=..&db=plant|ictv（已选等级作为过滤）→
    {selected, db, ranks:[{col, zh, options:[{name,n}]}]}。

    宿主口径由 db 决定：plant = 植物口径，ictv = 全病毒界不过滤。"""
    from Virus_Platform_Core import ictv_db
    if not ictv_db.available():
        # 没有数据库包的版本（01-Standalone）里没有 ICTV 谱系表。这不是错误：
        # 按 200 + available=false 正常返回，让前端显示中性提示。此前 abort(400)
        # 会让页面控制台出现 4xx 子资源（验收脚本记为缺陷），且提示文案里还留着
        # 过时的 `python main.py ictv-update` 命令。
        _db0 = str(request.args.get('db') or '').strip().lower() or 'plant'
        return jsonify({
            'available': False,
            'selected': {},
            'db': _db0,
            'dbs': [],
            'ranks': [{'col': col, 'zh': zh, 'options': []}
                      for col, zh in ictv_db.RANK_COLS],
            'hint': ('ICTV 分类库未就绪——可在「数据库构建」页建库，或运行 '
                     'VirusPlatform.exe --cli ictv-update'),
        })
    db = _pick_db(request.args.get('db'))
    levels = {col: (request.args.get(col) or '').strip()
              for col, _zh in ictv_db.RANK_COLS}
    levels = {k: v for k, v in levels.items() if v}
    try:
        opts = ictv_db.rank_options(levels, db=db)
    except Exception as e:
        abort(400, f'ICTV taxa 读取失败: {e}')
    return jsonify({
        'selected': levels,
        'db': db,
        'dbs': _db_options(),
        'ranks': [{'col': col, 'zh': zh, 'options': opts.get(col, [])}
                  for col, zh in ictv_db.RANK_COLS]})


def _db_options():
    """供前端「库」下拉渲染：[{key, label, n}]。

    不带本地序列（fasta 缺失，只有谱系表）的库追加「按需下载」标注——
    此时选参照常可用，只是 Source 列全为 `ncbi`，序列挑中后才联网取。
    """
    from Virus_Platform_Core import acvirus
    labels = {'plant': '植物病毒库 plant',
              'ictv': '全病毒界库 ictv'}
    out = []
    for k in acvirus._LIBS:
        try:
            ok = acvirus.available(k)
            seqs = acvirus.has_seqs(k)
        except Exception:
            ok, seqs = False, False
        if ok:
            label = labels.get(k, k) + ('' if seqs else '（按需下载）')
            out.append({'key': k, 'label': label})
    return out


def _ictv_levels(body):
    """从请求体取分类级选择 {列名: 值}（非空）+ 最深层（用于日志与定位）。"""
    from Virus_Platform_Core import ictv_db
    levels = {}
    for col, _zh in ictv_db.RANK_COLS:
        v = str((body.get(col) or '')).strip()
        if v:
            levels[col] = v
    if not levels:
        abort(400, '请至少选择一个分类等级（界/门/纲/目/科/属/种）')
    deepest_col = [c for c, _zh in ictv_db.RANK_COLS if c in levels][-1]
    return levels, deepest_col, levels[deepest_col]


@bp.route('/api/ictv/preview', methods=['POST'])
def api_ictv_preview():
    """body: {<分类级列名>: 值, ..., limit, genome, db} → 选参预览（不下载）。"""
    body = request.get_json(force=True) or {}
    levels, dcol, dval = _ictv_levels(body)
    db = _pick_db(body.get('db'))
    genome = body.get('genome') or 'complete'
    limit = max(1, min(int(body.get('limit', 200) or 200), 2000))
    per_genus = max(0, min(int(body.get('per_genus', 0) or 0), 50))
    per_species = max(0, min(int(body.get('per_species', 0) or 0), 20))
    from Virus_Platform_Core import ictv_db
    try:
        rows, total = ictv_db.select_refs(ranks=levels, genome=genome,
                                          db=db,
                                          limit=max(limit * 20, 2000))
    except Exception as e:
        abort(400, f'选参失败: {e}')
    # 抽样上限先于 limit 生效（先每属/每种封顶，再取前 limit 条）
    rows, dropped = ictv_db.cap_per_rank(rows, per_genus, per_species)
    rows = rows[:limit]
    cap_note = []
    if per_genus:
        cap_note.append(f'每属≤{per_genus}')
    if per_species:
        cap_note.append(f'每种≤{per_species}')
    return jsonify({'scope': f'{dcol}={dval}', 'total': total,
                    'n_shown': len(rows),
                    'sampling': '、'.join(cap_note) or '不限',
                    'dropped_by_cap': dropped,
                    'rows': [{'acc': r.get('Genbank', ''),
                              'source': r.get('Source', ''),
                              'species': r.get('Species', ''),
                              'coverage': r.get('Genome_Coverage', ''),
                              'genome': r.get('Genome', ''),
                              'host': r.get('Host_Source', '')}
                             for r in rows]})


@bp.route('/api/ictv/download', methods=['POST'])
def api_ictv_download():
    """body: {rank, name, coll, limit, genome} → 后台任务：
    ICTV 选参（本地优先）→ 下载整科/整属 GenBank 集合（含巡检清单）
    → 同步导出 FASTA 参考集（databases/ncbi_refs/<coll>/refs.fa，
    供样品流程建树追加参考）。"""
    body = request.get_json(force=True) or {}
    levels, dcol, dval = _ictv_levels(body)
    coll = (body.get('coll') or '').strip()
    genome = body.get('genome') or 'complete'
    db = _pick_db(body.get('db'))
    if not coll:
        abort(400, '参数不完整（请填集合名）')
    limit = max(1, min(int(body.get('limit', 200) or 200), 2000))
    per_genus = max(0, min(int(body.get('per_genus', 0) or 0), 50))
    per_species = max(0, min(int(body.get('per_species', 0) or 0), 20))

    def job(log, prog, cancel):
        logger = TaskLogger(callback=log)
        from Virus_Platform_Core import ictv_db
        from Virus_Platform_Core.gb_collection import (download_gb_collection, gb_collection_dir,
                                      collection_records)
        from Virus_Platform_Core.ncbi_download import collection_dir as fasta_coll_dir
        from Virus_Platform_Core.utils import write_fasta_record

        prog('ictv', 0.05, f'ICTV 选参：{dcol} = {dval}（{len(levels)} 级过滤，库={db}）')
        rows, total = ictv_db.select_refs(ranks=levels, genome=genome, db=db,
                                          limit=max(limit * 20, 2000))
        rows, _dropped = ictv_db.cap_per_rank(rows, per_genus, per_species)
        rows = rows[:limit]
        if not rows:
            raise RuntimeError(f'ICTV {dcol} "{dval}" 无匹配参考'
                               '（换一个分类或改 genome=any）')
        accs = [r['Genbank'] for r in rows]
        logger.log(f'ICTV 命中 {total} 条，取前 {len(accs)} 条'
                   f'（{"Complete genome" if genome == "complete" else "全部"}）')

        # 重名集合替换式重建（与建库口径一致，不累积重复记录）
        shutil.rmtree(gb_collection_dir(coll), ignore_errors=True)
        prog('ictv', 0.15, f'下载 GenBank 集合 {coll}（{len(accs)} 条 accession）')
        res = download_gb_collection(coll, accessions=','.join(accs),
                                     max_records=len(accs),
                                     logger=logger, prog=prog, cancel=cancel)

        # 同步 FASTA 参考集（样品流程「NCBI 参考集合」可直接填 coll）
        prog('ictv', 0.85, '同步 FASTA 参考集（样品建树用）')
        fdir = fasta_coll_dir(coll)
        os.makedirs(fdir, exist_ok=True)
        fa = os.path.join(fdir, 'refs.fa')
        n_fa = 0
        with safe_open(fa, 'wt') as f:
            for acc, _org, seq in collection_records(coll):
                if seq:
                    write_fasta_record(f, acc, seq)
                    n_fa += 1
        logger.log(f'FASTA 参考集: {fa}（{n_fa} 条）')
        prog('ictv', 1.0, f'完成：GenBank {res.get("n_records", n_fa)} 条 + '
                    f'FASTA {n_fa} 条')
        logger.close()
        return {'collection': coll, 'scope': f'{dcol}={dval}',
                'levels': levels, 'n_records': n_fa, 'n_fasta': n_fa}

    tid = tm.start(cfg.tr(f'ICTV下载 {coll}', f'ICTV download {coll}'), job,
                   weight='light')
    return jsonify({'task': tid})


@bp.route('/api/ncbi/search', methods=['POST'])
def api_ncbi_search():
    """body: {term, db, limit} → {count, rows:[{acc,title,organism,...}]}。"""
    body = request.get_json(force=True) or {}
    term = (body.get('term') or '').strip()
    db = body.get('db') or 'nucleotide'
    if not term or db not in ('nucleotide', 'protein'):
        abort(400, '参数不完整（term 必填，db=nucleotide/protein）')
    from Virus_Platform_Core.ncbi_download import search_preview
    try:
        return jsonify(search_preview(term, db=db,
                                      limit=int(body.get('limit', 50))))
    except RuntimeError as e:
        abort(502, str(e))


@bp.route('/api/ncbi/download', methods=['POST'])
def api_ncbi_download():
    """body: {term, name, db, max} → 后台任务下载集合。"""
    body = request.get_json(force=True) or {}
    term = (body.get('term') or '').strip()
    name = (body.get('name') or '').strip()
    db = body.get('db') or 'nucleotide'
    if not term or not name or db not in ('nucleotide', 'protein'):
        abort(400, '参数不完整（term/name 必填，db=nucleotide/protein）')
    max_records = max(1, min(int(body.get('max', 100) or 100), 10000))

    def job(log, prog, cancel):
        logger = TaskLogger(callback=log)
        from Virus_Platform_Core.ncbi_download import download_collection
        res = download_collection(term, name, db=db, max_records=max_records,
                                  logger=logger, prog=prog, cancel=cancel)
        logger.close()
        return res

    tid = tm.start(cfg.tr(f'NCBI下载 {name}', f'NCBI download {name}'), job,
                   weight='light')
    return jsonify({'task': tid})


@bp.route('/api/gb/collections')
def api_gb_collections():
    """集合列表（只读概览：读已生成的清单/警告文件，不重新解析 .gb）。

    下划线前缀的集合（_内部测试/_归档）不在列表展示。
    """
    from Virus_Platform_Core.gb_collection import list_gb_collections, read_manifest
    out = []
    for c in list_gb_collections():
        if c.get('name', '').startswith('_'):
            continue
        try:
            st = read_manifest(c['name'])
            c['records'] = [{'acc': r['acc'], 'organism': r['organism'],
                             'length': r['length'], 'cds': r['cds_count']}
                            for r in st['records']]
            c['warnings'] = st['warnings']
        except (OSError, FileNotFoundError, ValueError):
            c['records'], c['warnings'] = [], []
        c['has_phylo'] = os.path.isfile(
            os.path.join(c['dir'], 'phylo', 'aln.fasta'))
        out.append(c)
    return jsonify(out)


@bp.route('/api/gb/extract', methods=['POST'])
def api_gb_extract():
    """body: {name} → 后台任务：集合 .gb → extract/ 分类分目录提取
    （genome.fa / CDS.fa / PEP.fa / CDS/<基因>.fa / PEP/<基因>.fa / genes.tsv）。"""
    body = request.get_json(force=True) or {}
    name = (body.get('name') or '').strip()
    if not name:
        abort(400, '参数不完整（name 必填）')

    def job(log, prog, cancel):
        logger = TaskLogger(callback=log)
        from Virus_Platform_Core.gb_collection import extract_collection_features
        res = extract_collection_features(name, logger=logger, prog=prog)
        logger.close()
        return res

    tid = tm.start(cfg.tr(f'CDS/PEP 提取 {name}', f'Extract CDS/PEP {name}'),
                   job, weight='light')
    return jsonify({'task': tid})


@bp.route('/api/gb/inspect', methods=['POST'])
def api_gb_inspect():
    """body: {name} → 后台任务全量巡检（重解析全部 .gb，重建清单与警告）。"""
    body = request.get_json(force=True) or {}
    name = (body.get('name') or '').strip()
    if not name:
        abort(400, '参数不完整（name 必填）')

    def job(log, prog, cancel):
        logger = TaskLogger(callback=log)
        from Virus_Platform_Core.gb_collection import inspect_collection
        st = inspect_collection(name, logger=logger)
        logger.close()
        return {'name': name, 'n_records': len(st['records']),
                'n_warnings': len(st['warnings'])}

    tid = tm.start(cfg.tr(f'GenBank巡检 {name}', f'GenBank inspect {name}'), job,
                   weight='light')
    return jsonify({'task': tid})


@bp.route('/api/gb/phylo', methods=['POST'])
def api_gb_phylo():
    """body: {name, tree_tool} → 后台任务：集合全基因组比对 + 建树。

    只建全基因组树（molecule=genome，产物 phylo/，结果中心 gb:<name> 展示）。
    基因 CDS / 蛋白级建树已下线：解析 .gb 注释得到的基因序列用于建树
    不可靠，基因级建树请用「序列比对」模块的挑选序列集产物。
    """
    body = request.get_json(force=True) or {}
    name = (body.get('name') or '').strip()
    tree_tool = body.get('tree_tool') or 'fasttree'
    if not name:
        abort(400, '参数不完整（name 必填）')
    if tree_tool not in ('fasttree', 'raxml-ng', 'nj'):
        abort(400, 'tree_tool 仅支持 fasttree / raxml-ng / nj')

    def job(log, prog, cancel):
        logger = TaskLogger(callback=log)
        from Virus_Platform_Core.gb_collection import build_collection_phylo
        res = build_collection_phylo(name, tree_tool=tree_tool,
                                     molecule='genome',
                                     logger=logger, prog=prog, cancel=cancel)
        logger.close()
        return res

    tid = tm.start(cfg.tr(f'集合建树 {name}', f'Collection phylo {name}'), job)
    return jsonify({'task': tid})


@bp.route('/api/gb/download', methods=['POST'])
def api_gb_download():
    """body: {name, term?/accessions?, max} → 后台任务下载 GenBank 集合。"""
    body = request.get_json(force=True) or {}
    name = (body.get('name') or '').strip()
    term = (body.get('term') or '').strip()
    accs = (body.get('accessions') or '').strip()
    if not name or (not term and not accs):
        abort(400, '参数不完整（name 与 term/accessions 必填）')
    if term and accs:
        abort(400, '检索式与 accession 列表二选一')
    max_records = max(1, min(int(body.get('max', 50) or 50), 10000))

    def job(log, prog, cancel):
        logger = TaskLogger(callback=log)
        from Virus_Platform_Core.gb_collection import (download_gb_collection,
                                      collection_records)
        res = download_gb_collection(name, term=term or None,
                                     accessions=accs or None,
                                     max_records=max_records,
                                     logger=logger, prog=prog, cancel=cancel)
        # 同步 FASTA 参考集（样品流程「NCBI 参考集合」可直接填集合名）
        try:
            from Virus_Platform_Core.ncbi_download import collection_dir as fasta_coll_dir
            from Virus_Platform_Core.utils import write_fasta_record
            fdir = fasta_coll_dir(name)
            os.makedirs(fdir, exist_ok=True)
            fa = os.path.join(fdir, 'refs.fa')
            n_fa = 0
            with safe_open(fa, 'wt') as f:
                for acc, _org, seq in collection_records(name):
                    if seq:
                        write_fasta_record(f, acc, seq)
                        n_fa += 1
            logger.log(f'FASTA 参考集: {fa}（{n_fa} 条）')
            res['n_fasta'] = n_fa
        except Exception as e:
            logger.log(f'FASTA 参考集同步失败（不影响 GenBank 集合）: {e}',
                       'WARN')
        logger.close()
        return res

    tid = tm.start(cfg.tr(f'GenBank下载 {name}', f'GenBank download {name}'), job,
                   weight='light')
    return jsonify({'task': tid})


@bp.route('/api/gb/import', methods=['POST'])
def api_gb_import():
    """body: {name, files:[本机 .gb 路径]} → 后台任务导入。"""
    body = request.get_json(force=True) or {}
    name = (body.get('name') or '').strip()
    files = body.get('files') or []
    if not name or not files:
        abort(400, '参数不完整（name 与 files 必填）')

    def job(log, prog, cancel):
        logger = TaskLogger(callback=log)
        from Virus_Platform_Core.gb_collection import import_local_gb
        res = import_local_gb(name, files, logger=logger, prog=prog)
        logger.close()
        return res

    tid = tm.start(cfg.tr(f'GenBank导入 {name}', f'GenBank import {name}'), job,
                   weight='light')
    return jsonify({'task': tid})


def _msa_group_dir(sample, group):
    """sample+group → 组目录（含 gb: 集合伪样品解析）。"""
    if sample.startswith('gb:'):
        from Virus_Platform_Core.gb_collection import gb_collection_dir
        try:
            return check_path(os.path.join(gb_collection_dir(sample[3:]),
                                           'phylo'),
                              must_exist=True, in_platform=True)
        except ValueError as e:
            raise ValueError(f'集合名不合法: {e}')
    from Virus_Platform_Core.utils import resolve_sample_name
    from Virus_Platform_Core.phylo import _safe_name
    # 样品名来自 /api/msa/samples（真实目录名）或查询参数 → 按磁盘解析
    sdir = check_path(os.path.join(DIRS['results'],
                                   resolve_sample_name(sample, DIRS['results'])),
                      must_exist=True, in_platform=True)
    return check_path(os.path.join(sdir, '05_phylo', _safe_name(group)),
                      must_exist=True, in_platform=True)


def _list_phylo_groups(pdir):
    """目录 → [{group, aln, trees}]（含比对文件的分组，优先 aln.trim）。"""
    groups = []
    for g in sorted(os.listdir(pdir)):
        gdir = os.path.join(pdir, g)
        if not os.path.isdir(gdir):
            continue
        aln = os.path.join(gdir, 'aln.trim.fasta')
        kind = 'trim'
        if not os.path.isfile(aln):
            aln = os.path.join(gdir, 'aln.fasta')
            kind = 'full'
        if os.path.isfile(aln):
            groups.append({'group': g, 'aln': kind, 'trees': _tree_files(gdir)})
    return groups


@bp.route('/api/msa/samples')
def api_msa_samples():
    """列出含比对结果的样品与集合建树产物（伪样品 gb:<集合名>）。"""
    out = []
    res = check_path(DIRS['results'], must_exist=False, in_platform=True)
    if os.path.isdir(res):
        for name in sorted(os.listdir(res)):
            if name.startswith('_'):      # 归档/内部样品不进下拉
                continue
            pdir = os.path.join(res, name, '05_phylo')
            if not os.path.isdir(pdir):
                continue
            groups = _list_phylo_groups(pdir)
            if groups:
                out.append({'sample': name, 'groups': groups})
    # GenBank 集合建树产物（比较基因组：下载/导入集合 → 直接建树）
    from Virus_Platform_Core.gb_collection import list_gb_collections
    try:
        cols = list_gb_collections()
    except OSError:
        cols = []
    for c in cols:
        pdir = os.path.join(c['dir'], 'phylo')
        if not os.path.isdir(pdir):
            continue
        aln = os.path.join(pdir, 'aln.trim.fasta')
        if not os.path.isfile(aln):
            aln = os.path.join(pdir, 'aln.fasta')
        if os.path.isfile(aln):
            out.append({'sample': f'gb:{c["name"]}',
                        'label': f'🧬 {c["name"]}（集合）',
                        'groups': [{'group': c['name'], 'aln': 'full',
                                    'trees': _tree_files(pdir)}]})
    return jsonify(out)


_TREE_FILES = (('tree.nwk', 'FastTree'),
               ('nj.nwk', 'NJ 快速树（identity 距离）'),
               ('raxml.raxml.support', 'RAxML-NG ML 树（FBP 支持值）'),
               ('raxml.raxml.bestTree', 'RAxML-NG ML 树（无支持值）'),
               ('iqtree.treefile', 'IQ-TREE ML 树（旧产物）'),
               ('iqtree.contree', 'IQ-TREE 一致树（旧产物）'))


def _tree_files(gdir):
    """组目录下实际存在的树文件列表（供树查看器下拉）。"""
    return [{'file': fn, 'label': label}
            for fn, label in _TREE_FILES
            if os.path.isfile(os.path.join(gdir, fn))]


@bp.route('/api/tree/data')
def api_tree_data():
    """sample+group(+file) → 进化树 Newick 与建树摘要（前端 Archaeopteryx.js 渲染）。"""
    sample = request.args.get('sample') or ''
    group = request.args.get('group') or ''
    fname = request.args.get('file') or ''
    try:
        gdir = _msa_group_dir(sample, group)
    except ValueError as e:
        abort(400, f'参数不合法: {e}')
    # 只接受 _TREE_FILES 白名单内的文件名：原先直接 os.path.join(gdir, fname)，
    # 而 os.path.join 遇到绝对路径会丢弃 gdir → file=C:\Windows\win.ini 可读
    # 任意文件（safe_open 读模式不限定平台内）。
    allowed = {fn for fn, _ in _TREE_FILES}
    if fname and fname not in allowed:
        abort(400, '不支持的树文件名')
    choices = ([(fname, '')] if fname else list(_TREE_FILES))
    tree_file = None
    for fn, _label in choices:
        p = os.path.join(gdir, fn)
        if fn and os.path.isfile(p):
            tree_file = p
            break
    if not tree_file:
        abort(404, '该组没有树文件（tree.nwk / nj.nwk / raxml.raxml.support）')
    with safe_open(tree_file) as f:
        nwk = f.read().strip()
    if not nwk:
        abort(404, '树文件为空')
    # 从 summary.json 取建树工具/模型等摘要（无则忽略）
    meta = {}
    summary = os.path.join(gdir, 'summary.json')
    if os.path.isfile(summary):
        try:
            with safe_open(summary) as f:
                sj = json.load(f)
            g = (next((x for x in (sj.get('groups') or [])
                       if x.get('group') == group), None)
                 if isinstance(sj.get('groups'), list) else None) or {}
            base = os.path.basename(tree_file)
            if base == 'nj.nwk':
                tool = 'NJ'
            elif (base.startswith('raxml.')
                    or g.get('tree_tool') == 'raxml-ng'
                    or sj.get('tree_tool') == 'raxml-ng'):
                tool = 'RAxML-NG'
            elif (base.endswith('.treefile') or base.endswith('.contree')
                    or g.get('tree_tool') == 'iqtree'
                    or sj.get('tree_tool') == 'iqtree'):
                tool = 'IQ-TREE（旧产物）'
            else:
                tool = 'FastTree'
            meta = {'tool': tool,
                    'model': g.get('tree_model') or sj.get('model'),
                    'logl': g.get('tree_logl') or sj.get('logl')}
        except (OSError, ValueError):
            pass
    return jsonify({'newick': nwk, 'file': os.path.basename(tree_file),
                    # 「进化树 + 基因组叠加」：同名 <stem>.overlay.json 有就带上
                    'overlay': tree_overlay_for(tree_file), **meta})


@bp.route('/api/msa/data')
def api_msa_data():
    """sample+group → SNP-only 展示数据（变异列矩阵/共识/多样性）。"""
    sample = request.args.get('sample') or ''
    group = request.args.get('group') or ''
    try:
        gdir = _msa_group_dir(sample, group)
    except ValueError as e:
        abort(400, f'参数不合法: {e}')
    aln = os.path.join(gdir, 'aln.trim.fasta')
    if not os.path.isfile(aln):
        aln = os.path.join(gdir, 'aln.fasta')
    if not os.path.isfile(aln):
        abort(404, '该组无比对文件')
    from Virus_Platform_Core.msa_view import snp_view_data
    try:
        return jsonify(snp_view_data(aln))
    except ValueError as e:
        abort(400, str(e))


@bp.route('/api/msa/path')
def api_msa_path():
    """sample+group → 该组比对文件的平台相对路径（供比对查看器直接打开）。"""
    sample = request.args.get('sample') or ''
    group = request.args.get('group') or ''
    try:
        gdir = _msa_group_dir(sample, group)
    except ValueError as e:
        abort(400, f'参数不合法: {e}')
    aln = os.path.join(gdir, 'aln.trim.fasta')
    if not os.path.isfile(aln):
        aln = os.path.join(gdir, 'aln.fasta')
    if not os.path.isfile(aln):
        abort(404, '该组无比对文件')
    rel = os.path.relpath(aln, PLATFORM_ROOT).replace(os.sep, '/')
    return jsonify({'path': rel, 'name': os.path.basename(aln)})


@bp.route('/api/sdt/data')
def api_sdt_data():
    """sample+group → SDT 口径全长成对 identity 矩阵（结果中心在线热图）。

    优先读样品流程产物 sdt_matrix.csv；无则从比对 fasta 现算
    （与 Virus_Platform_Core/phylo.pairwise_identity_matrix 同口径）。
    """
    sample = request.args.get('sample') or ''
    group = request.args.get('group') or ''
    try:
        gdir = _msa_group_dir(sample, group)
    except ValueError as e:
        abort(400, f'参数不合法: {e}')
    csv = os.path.join(gdir, 'sdt_matrix.csv')
    if os.path.isfile(csv):
        names, mat = [], []
        with safe_open(csv) as f:
            rows = [line.rstrip('\n').split(',') for line in f if line.strip()]
        if len(rows) >= 2:
            names = rows[0][1:]
            for r in rows[1:]:
                mat.append([float(x) if x else 0.0 for x in r[1:]])
            return jsonify({'names': names, 'matrix': mat,
                            'source': 'sdt_matrix.csv'})
    aln = os.path.join(gdir, 'aln.trim.fasta')
    if not os.path.isfile(aln):
        aln = os.path.join(gdir, 'aln.fasta')
    if not os.path.isfile(aln):
        abort(404, '该组没有 SDT 矩阵，也没有比对文件')
    from Virus_Platform_Core.phylo import pairwise_identity_matrix
    names, mat = pairwise_identity_matrix(aln)
    return jsonify({'names': names, 'matrix': mat, 'source': 'aln'})
