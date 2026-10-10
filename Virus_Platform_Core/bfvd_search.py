# -*- coding: utf-8 -*-
"""BFVD 病毒域搜索统一入口（⑥b 注释层2' 与 ③b 验证 B 路共用）。

背景（2026-09-29）
------------------
CDD 层（mmseqs vs NCBI Cdd）因本地库 prefilter 索引一次性分配 ~5GB、
小内存机器反复 OOM 而移除；改用 **BFVD**（Steinegger lab 病毒蛋白域库，
NCBI 官方病毒提交工具链 suvtk 的功能注释证据源）顶替：

- 库已在平台内：databases/misc_db/suvtk/bfvd（212MB 全套，mmseqs 格式），
  与 SUVTK 提交链共用同一份文件，随数据库包分发；prefilter 索引百 MB 级
  ——整库直搜即可，无需 cdd_search 时代的 8 分片 workaround；
- 病毒专属库：命中即病毒证据，**无需** CDD 时代的 1,580 条病毒白名单/tier；
- 命名与提交链同款：product 取 bfvd_uniprot_names.tsv 的 Protein names
  （去尾部括号段），organism 取同表 Organism——与 suvtk features 写进
  .sqn featuretable 的 product 完全同源，⑥b 预览即提交所见。

suvtk features 的搜索口径（对齐它以保证结果一致）：
  mmseqs easy-search <faa> bfvd out.m8 tmp -s 7.5 --format-mode 0
    --format-output query,target,pident,alnlen,mismatch,gapopen,
                    qstart,qend,tstart,tend,evalue,bits
  evalue < 1e-3

名字解析链（见 load_names）：m8 target（AlphaFold 模型名，无 .pdb 后缀）
→ bfvd_metadata.tsv（model → Uniref_entry）
→ bfvd_uniprot_names.tsv（From → Protein names / Organism）。
"""
import csv
import os
import re

from .config import get_config, db_path
from .utils import run_cmd

BFVD_EVALUE_MAX = 1e-3
BFVD_MAX_HITS = 3

# 尾部括号段正则（与 suvtk features.py 同款）："Replicase polyprotein (polyprotein)"
# → "Replicase polyprotein"
_TRAIL_BRACKET = re.compile(r'[\(\[].*?[\)\]]$')

_NAMES_CACHE = None

# mmseqs 内存分配失败的报错特征（小库一般不会触发，降线程阶梯仅作保险）
_ALLOC_MSGS = ('Can not allocate', 'could not allocate', 'IndexTable')


def bfvd_db(logger=None):
    """mmseqs 格式 BFVD 库前缀。优先级：platform.json databases.bfvd 覆盖 >
    默认 databases/misc_db/suvtk/bfvd（与 SUVTK 提交链共用同一份库）。"""
    cfg = get_config()
    p = (getattr(cfg, 'databases', {}) or {}).get('bfvd')
    if p and os.path.isfile(p):
        return p
    prefix = os.path.join(db_path('misc', 'suvtk'), 'bfvd')
    return prefix if os.path.isfile(prefix) else None


def bfvd_db_dir():
    """BFVD 库所在目录（两张 TSV 名字表同目录）。"""
    d = (getattr(get_config(), 'databases', {}) or {}).get('bfvd')
    if d and os.path.isfile(d):
        return os.path.dirname(d)
    return db_path('misc', 'suvtk')


def db_trace():
    """BFVD 库可追溯信息 → (db_files, db_versions)。

    db_files 覆盖库目录下 bfvd* 全套（mmseqs 索引 + 两张名字表），
    db_versions 取 bfvd.version 首行内容；供 ⑥b summary.json 记录
    （2026-10-03 审计 P2）。
    """
    from .utils import db_trace as _trace
    d = bfvd_db_dir()
    files, versions = _trace([d], version_files=[os.path.join(d, 'bfvd.version')],
                             n_max=40)
    return ({k: v for k, v in files.items() if k.startswith('bfvd')},
            versions)


def refresh_paths():
    """切换数据库目录后失效名字缓存（config._refresh_module_paths 回调）。"""
    global _NAMES_CACHE
    _NAMES_CACHE = None


def load_names(db_dir=None):
    """{model(去 .pdb): {'product','organism','uniref'}}，全量载入并缓存。

    命名链与 suvtk features.py 同款（见模块 docstring）；两张 TSV 任一缺失
    时返回空 dict（搜索照常，只是 product 退化为模型名/Uniref 代号）。
    """
    global _NAMES_CACHE
    if _NAMES_CACHE is not None:
        return _NAMES_CACHE
    d = db_dir or bfvd_db_dir()
    model2uniref = {}
    meta_p = os.path.join(d, 'bfvd_metadata.tsv')
    if os.path.isfile(meta_p):
        with open(meta_p, encoding='utf-8', errors='replace') as f:
            for row in csv.reader(f, delimiter='\t'):
                # 列序: Uniref_entry  model(.pdb)  length  avg_pLDDT  pTM  splitted
                if len(row) >= 2 and row[0] and row[1]:
                    model2uniref[row[1].strip().removesuffix('.pdb')] = \
                        row[0].strip()
    out = {}
    names_p = os.path.join(d, 'bfvd_uniprot_names.tsv')
    if os.path.isfile(names_p) and model2uniref:
        need = set(model2uniref.values())
        with open(names_p, encoding='utf-8', errors='replace') as f:
            for row in csv.DictReader(f, delimiter='\t'):
                acc = (row.get('From') or '').strip()
                if acc not in need:
                    continue
                product = _TRAIL_BRACKET.sub(
                    '', (row.get('Protein names') or '').strip())
                out[acc] = {'product': product,
                            'organism': (row.get('Organism') or '').strip()}
    by_model = {}
    for model, uniref in model2uniref.items():
        info = out.get(uniref)
        by_model[model] = {
            'product': (info or {}).get('product', ''),
            'organism': (info or {}).get('organism', ''),
            'uniref': uniref,
        }
    _NAMES_CACHE = by_model
    return by_model


def easy_search(mmseqs, query, db, out_tsv, tmp, format_output, evalue,
                max_seqs=5, threads=4, extra=None, logger=None):
    """easy-search vs BFVD（整库直搜；内存分配失败按线程阶梯降档重试）。

    成功/失败都要清掉 mmseqs 的 tmp 目录：残留的非空 tmp 会让下次重试
    直接报错，还会被 /api/tool/runs 目录遍历扫到（CDD 时代的实测教训）。
    """
    import shutil
    ladder = []
    for t in (threads, min(threads, 8), min(threads, 4), min(threads, 2), 1):
        t = max(1, int(t))
        if t not in ladder:
            ladder.append(t)
    last = None
    try:
        for t in ladder:
            shutil.rmtree(tmp, ignore_errors=True)
            cmd = [mmseqs, 'easy-search', query, db, out_tsv, tmp,
                   '--format-output', format_output,
                   '-e', str(evalue), '--max-seqs', str(max_seqs),
                   '--threads', str(t)] + list(extra or [])
            try:
                run_cmd(cmd, logger=logger)
                if t != ladder[0] and logger:
                    logger.log(f'BFVD 搜索 {ladder[0]} 线程内存不足，'
                               f'降为 {t} 线程后成功', 'WARN')
                return
            except RuntimeError as e:
                if not any(m in str(e) for m in _ALLOC_MSGS):
                    raise
                last = e
                if logger and t != ladder[-1]:
                    logger.log(f'BFVD 搜索 {t} 线程内存分配失败，降线程重试',
                               'WARN')
        raise last
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def _mmseqs_bin():
    return get_config().tool('mmseqs')


def search_orfs(faa, out_dir, threads=None, logger=None):
    """⑥b 注释层：mmseqs easy-search ORF 蛋白 vs BFVD。

    口径与 suvtk features 一致（-s 7.5，evalue<1e-3）。返回
    {orf_id: [hit, ...]} 或 None（库/工具不可用）。
    hit: {'lib':'bfvd','target':模型名,'acc':Uniref,'desc':product,
          'organism':…,'evalue','cov'}
    """
    cfg = get_config()
    try:
        mmseqs = _mmseqs_bin()
    except FileNotFoundError:
        if logger:
            logger.log("mmseqs2 未安装，跳过 BFVD 病毒域层", "WARN")
        return None
    db = bfvd_db()
    if not db:
        if logger:
            logger.log("未找到 BFVD 库（databases/misc_db/suvtk/bfvd），"
                       "跳过病毒域层", "WARN")
        return None
    if logger:
        logger.log(f"  [bfvd] 使用库: {db} "
                   f"({os.path.getsize(db) / 1e6:.0f}MB)")
    out_tsv = os.path.join(out_dir, 'bfvd_hits.raw.tsv')
    tmp = out_tsv + '.mmseqs_tmp'
    easy_search(mmseqs, faa, db, out_tsv, tmp,
                format_output='query,target,evalue,pident,qlen,qstart,qend',
                evalue=BFVD_EVALUE_MAX, max_seqs=5,
                threads=threads or cfg.threads,
                extra=['-s', '7.5'], logger=logger)
    names = load_names()
    by_q = {}
    # 结果文件可能不存在（easy_search 异常路径/库为空时未产出）——
    # 对齐 search_contigs:242 的 isfile 保护，避免 FileNotFoundError
    if not os.path.isfile(out_tsv):
        if logger:
            logger.log(f"  [bfvd] 未生成结果文件 {os.path.basename(out_tsv)}，"
                       f"按零命中处理", "WARN")
        return by_q
    with open(out_tsv, encoding='utf-8', errors='replace') as f:
        for line in f:
            p = line.rstrip('\n').split('\t')
            if len(p) < 7:
                continue
            q, target = p[0].split()[0], p[1].split()[0]
            try:
                ev = float(p[2])
                qlen, qs, qe = int(p[4]), int(p[5]), int(p[6])
            except ValueError:
                continue
            if ev > BFVD_EVALUE_MAX:
                continue
            info = names.get(target, {})
            by_q.setdefault(q, []).append({
                'lib': 'bfvd', 'target': target,
                'acc': info.get('uniref', ''),
                'desc': info.get('product', ''),
                'organism': info.get('organism', ''),
                'evalue': ev,
                'cov': round((qe - qs + 1) / max(qlen, 1), 2),
            })
    for q in by_q:
        by_q[q] = sorted(by_q[q], key=lambda h: h['evalue'])[:BFVD_MAX_HITS]
    if logger:
        logger.log(f"  [bfvd] {len(by_q)} 个 ORF 获得病毒域命中")
    return by_q


def search_contigs(fa, out_tsv, threads=None, logger=None):
    """③b 验证 B 路：mmseqs easy-search 核酸 contig vs BFVD（自动翻译）。

    BFVD 病毒专属，命中即病毒域证据（无白名单步骤）。
    返回 {contig: [hit, ...]}；hit: {'target','acc','product','organism',
    'evalue','cov'}。
    """
    cfg = get_config()
    mmseqs = _mmseqs_bin()
    db = bfvd_db()
    if not db:
        raise FileNotFoundError('未找到 BFVD 库（databases/misc_db/suvtk/bfvd）')
    tmp = out_tsv + '.mmseqs_tmp'
    easy_search(mmseqs, fa, db, out_tsv, tmp,
                format_output='query,target,evalue,qlen,qstart,qend',
                evalue=BFVD_EVALUE_MAX, max_seqs=5,
                threads=threads or cfg.threads,
                extra=['-s', '7.5'], logger=logger)
    names = load_names()
    by_q = {}
    if not os.path.isfile(out_tsv):
        return by_q
    with open(out_tsv, encoding='utf-8', errors='replace') as f:
        for line in f:
            p = line.rstrip('\n').split('\t')
            if len(p) < 6:
                continue
            q, target = p[0].split()[0], p[1].split()[0]
            try:
                ev, qlen, qs, qe = float(p[2]), int(p[3]), int(p[4]), int(p[5])
            except ValueError:
                continue
            if ev > BFVD_EVALUE_MAX:
                continue
            info = names.get(target, {})
            by_q.setdefault(q, []).append({
                'target': target, 'acc': info.get('uniref', ''),
                'product': info.get('product', ''),
                'organism': info.get('organism', ''),
                'evalue': ev,
                'cov': round((qe - qs + 1) / max(qlen, 1), 3),
            })
    for q in by_q:
        by_q[q].sort(key=lambda h: h['evalue'])
    if logger:
        logger.log('  [bfvd] %d 条命中病毒域' % len(by_q))
    return by_q
