# -*- coding: utf-8 -*-
"""
阶段⑥b ORF 功能注释：
对 ⑥ORF 产出的病毒蛋白（pyrodigal_rv 默认优先，可选 pyrodigal meta 基因模型，
orfipy 仅作历史运行兜底）做蛋白级功能注释，回答"每个 ORF 是什么蛋白、属于
哪个科属、基因组什么类型"。

方法选型吸收自 LanderDC/annotation_benchmark（2026，对 1.1 万条 ICTV 参考
病毒蛋白的注释方法基准）：序列同源搜索是准确率最高、最可复现的主策略
（BLASTp 78.4% > 结构折叠法 72.2% > 蛋白语言模型 63.6%），嵌入/结构法仅对
大 dsDNA 病毒的残余假想蛋白有补充价值（植物 RNA 病毒收益低，暂不引入）。
参考库用 NCBI RefSeq release viral 蛋白（头部自带 "产物 [物种]"）；功能类别
归属采用其 E-value 加权类别投票规则（有效类别权重 ≥ 最佳无效类别的 50% 才
胜出）。物种 → 科/属经 NCBI taxonomy 谱系解析；科 → 宿主类别/基因组类型查
ICTV 科级对照表（数据源自 weinstein-bioinfo/LazypipeX 的 ICTV.virus.family.*
.tsv，265 个病毒科，覆盖全部植物病毒科）。

搜索引擎自动降级：DIAMOND（首选）→ MMseqs2 → BLAST+ blastp，三者输出均为
标准 tabular，下游解析统一。首次运行自动下载蛋白库（~107MB 压缩，NCBI
官方域名白名单）并构建搜索库与物种分类索引，之后全缓存。

输出（04b_orf_annot/）:
  orf_annotation.tsv           逐 ORF 注释表（最优命中 + 投票类别 + 科属 +
                               宿主/基因组类型）
  orf_function_summary.tsv     功能类别分布
  orf_family_summary.tsv       科级分布
  contig_function_profile.tsv  逐 contig 功能画像（ORF 数/已注释数/类别/主科）
  orf_annotation.gff3          pyrodigal GFF3 追加 product/organism/category
                               属性（供 ⑨基因组图按功能着色等下游使用）
  summary.json                 卡片/报告用汇总
"""
import os
import re
import json
import math
import glob
import socket
import ipaddress
import urllib.request
from urllib.parse import urlparse

from .config import DIRS, get_config, db_path
from .utils import (check_path, safe_open, safe_remove, run_cmd, iter_fasta,
                    count_fasta_seqs, is_step_done, mark_step_done, db_trace)
from .hmm_annot import (hmm_available, annotate_orfs_hmm,
                        merge_layer2, available_libs, lib_trace)
from .bfvd_search import (search_orfs as bfvd_search_orfs, bfvd_db,
                          db_trace as bfvd_db_trace)
from .genome_diag import render_sample_diagrams

DB_DIR = db_path('annot', 'prot')
REFSEQ_URL = ('https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/'
              'viral.{}.protein.faa.gz')
_DL_HOSTS = {'ftp.ncbi.nlm.nih.gov'}
ICTV_HOST_TSV = os.path.join(DB_DIR, 'ictv_family_host.tsv')
ICTV_GENCOMP_TSV = os.path.join(DB_DIR, 'ictv_family_gencomp.tsv')


def refresh_paths():
    """重算模块级库路径常量（切换「设置 → 数据库目录」后由 config 回调）。

    DB_DIR 在 import 时快照；不重算则 db-migrate 后仍指旧库，
    ICTV 宿主/基因组组成表会读不到（静默为空）。
    """
    global DB_DIR, ICTV_HOST_TSV, ICTV_GENCOMP_TSV
    DB_DIR = db_path('annot', 'prot')
    ICTV_HOST_TSV = os.path.join(DB_DIR, 'ictv_family_host.tsv')
    ICTV_GENCOMP_TSV = os.path.join(DB_DIR, 'ictv_family_gencomp.tsv')


def _rm_if_exists(path):
    if path and (os.path.isfile(path) or glob.glob(path + '*')):
        safe_remove(path)


def _preview_rank(r):
    """rows_preview 排序键：层1 命中按 bitscore 降序，层2 兜底行（bitscore
    为空串）统一排在最后。空串直接取负会 TypeError（此前层1 恒开、
    兜底行常同时有层1 命中，问题被掩盖）。"""
    try:
        return (0, -float(r.get('bitscore') or 0))
    except (TypeError, ValueError):
        return (1, 0.0)


def _ascii_blast_base():
    """BLAST LMDB 不支持中文路径：blastp 库固定建在 %TEMP%\\vp_blast。"""
    import tempfile
    return os.path.join(tempfile.gettempdir(), 'vp_blast')


def annotation_engine(prefer=None):
    """探测可用搜索引擎：DIAMOND → MMseqs2 → blastp。返回 (名称, 路径) 或 None。

    prefer 不为空时优先该引擎（如用户在前端指定 DIAMOND/MMseqs2/blastp）；
    指定引擎不可用时回退到自动探测。
    """
    cfg = get_config()
    order = (['diamond', 'mmseqs', 'blastp'] if prefer not in
             ('diamond', 'mmseqs', 'blastp') else [prefer, 'diamond', 'mmseqs', 'blastp'])
    for name in order:
        try:
            p = cfg.tool(name)
        except FileNotFoundError:
            continue
        if p and os.path.isfile(p):
            return name, p
    return None


# ------------------------------------------------------------------
# 参考库：下载 / 分类索引 / 搜索库构建
# ------------------------------------------------------------------
class _WhitelistRedirectHandler(urllib.request.HTTPRedirectHandler):
    """重定向逐跳复检白名单与解析 IP，防止经 30x 跳转到任意主机。"""

    def redirect_request(self, req, fp, code, msg, headers, newurl):
        u = urlparse(str(newurl))
        if u.scheme != 'https' or u.hostname not in _DL_HOSTS:
            raise ValueError(f"重定向域名不在白名单: {u.hostname}")
        for info in socket.getaddrinfo(u.hostname, 443,
                                       proto=socket.IPPROTO_TCP):
            ip = ipaddress.ip_address(info[4][0])
            if (ip.is_private or ip.is_loopback or ip.is_link_local
                    or ip.is_reserved or ip.is_multicast or ip.is_unspecified):
                raise ValueError(f"重定向解析到受限地址: {ip}")
        return super().redirect_request(req, fp, code, msg, headers, newurl)


_OPENER = urllib.request.build_opener(_WhitelistRedirectHandler)


def _download_whitelisted(url, dest, logger=None, progress=None,
                          label='RefSeq 病毒蛋白库'):
    """https 下载（域名白名单 + 拒绝解析到私网/环回/保留地址 + 重定向逐跳
    复检），带字节进度。dest 为模块内部固定路径（DB_DIR 下）。

    label 仅用于日志标题（universal_ref 下载 NCBI accession2taxid 时复用本
    函数，标题需区分，2026-10-03 审计 P3）。
    """
    u = urlparse(url)
    if u.scheme != 'https' or u.hostname not in _DL_HOSTS:
        raise ValueError(f"下载域名不在白名单: {u.hostname}")
    for info in socket.getaddrinfo(u.hostname, 443, proto=socket.IPPROTO_TCP):
        ip = ipaddress.ip_address(info[4][0])
        if (ip.is_private or ip.is_loopback or ip.is_link_local
                or ip.is_reserved or ip.is_multicast or ip.is_unspecified):
            raise ValueError(f"域名解析到受限地址: {ip}")
    if logger:
        logger.log(f"下载 {label}: {url}")
    req = urllib.request.Request(url, headers={'User-Agent': 'vp-orf-annot'})
    # 注意：不能直接 safe_open 写 .gz（会被 crabz 管道再压一层）——保存的是
    # NCBI 原生 gzip 字节流。先写 .part 临时名（普通二进制写），完成后原子改名。
    with _OPENER.open(req, timeout=300) as r:
        total = int(r.headers.get('Content-Length') or 0)
        done = 0
        dest = check_path(dest, must_exist=False, in_platform=True)
        part = dest + '.part'
        with safe_open(part, 'wb') as w:
            while True:
                chunk = r.read(1 << 20)
                if not chunk:
                    break
                w.write(chunk)
                done += len(chunk)
                if progress and total:
                    progress(done / total,
                             f"下载 {done // (1 << 20)}MB"
                             f"{' / ' + str(total // (1 << 20)) + 'MB' if total else ''}")
        os.replace(part, dest)
    if logger:
        logger.log(f"下载完成: {done / 1e6:.0f} MB")


def _collect_organisms(faa, logger=None):
    """扫描蛋白库头部：统计蛋白数、收集唯一物种名（头部末尾 [organism]）。
    返回 (n_proteins, set(organisms))。"""
    n, orgs = 0, set()
    tail = re.compile(r'\[([^\[\]]+)\]\s*$')
    with safe_open(faa) as f:
        for line in f:
            if not line.startswith('>'):
                continue
            n += 1
            m = tail.search(line)
            if m:
                orgs.add(m.group(1).strip())
    if logger:
        logger.log(f"参考库共 {n:,} 条蛋白，唯一物种 {len(orgs):,} 个")
    return n, orgs


def _build_organism_tax(faa, dest, logger=None):
    """经 NCBI taxonomy 解析每个物种的属/科：organism→taxid（names.dmp）
    →谱系（nodes.dmp）。输出 TSV: organism\ttaxid\tgenus\tfamily。"""
    n, orgs = _collect_organisms(faa, logger=logger)
    tax_dir = DIRS['taxonomy']
    names_dmp = os.path.join(tax_dir, 'names.dmp')
    nodes_dmp = os.path.join(tax_dir, 'nodes.dmp')
    if not (os.path.isfile(names_dmp) and os.path.isfile(nodes_dmp)):
        raise RuntimeError("NCBI taxonomy 未就绪（先在数据库页完成下载）")

    # 名称 → taxid（仅科学名；RefSeq 物种名为学名）。
    # dmp 行以 "\t|" 结尾，逐字段剥掉尾部空白与 '|' 再比对
    org2tid, remain = {}, set(orgs)
    fields = re.compile(r'\t\|\t')
    with safe_open(names_dmp) as f:
        for line in f:
            p = [x.strip(' \t|') for x in fields.split(line.rstrip('\r\n'))]
            if len(p) < 4 or p[3] != 'scientific name':
                continue
            name = p[1]
            if name in remain:
                org2tid[name] = p[0]
                remain.discard(name)
                if not remain:
                    break
    if remain and logger:
        logger.log(f"taxonomy 未收录物种 {len(remain)} 个（按空科属处理）", "WARN")

    # 谱系：parent 全量 + rank 稀疏（仅 genus/family 入表，控内存）
    parent, rank_map = {}, {}
    with safe_open(nodes_dmp) as f:
        for line in f:
            p = [x.strip(' \t|') for x in fields.split(line.rstrip('\r\n'))]
            if len(p) < 3:
                continue
            try:
                tid = int(p[0])
                parent[tid] = int(p[1])
            except ValueError:
                continue
            r = p[2]
            if r in ('genus', 'family'):
                rank_map[(tid, r)] = True

    def lineage(tid):
        # 沿 parent 链上溯，取首个 genus 与 family 节点（起始节点本身可能是
        # genus 级名称，如 RefSeq 未定种条目）
        genus = family = ''
        cur = tid
        seen = set()
        while cur and cur not in seen:
            seen.add(cur)
            if not genus and (cur, 'genus') in rank_map:
                genus = str(cur)
            if not family and (cur, 'family') in rank_map:
                family = str(cur)
            if genus and family:
                break
            cur = parent.get(cur, 0)
        return genus, family

    # 回填学名：genus/family 节点数量有限，二次流式查名
    need_names = set()
    t2gf = {}
    for org, tid_s in org2tid.items():
        try:
            g, fam = lineage(int(tid_s))
        except (ValueError, KeyError):
            g = fam = ''
        t2gf[org] = (g, fam)
        need_names.update(x for x in (g, fam) if x)
    name_of = {}
    if need_names:
        with safe_open(names_dmp) as f:
            for line in f:
                p = [x.strip(' \t|') for x in fields.split(line.rstrip('\r\n'))]
                if len(p) >= 4 and p[0] in need_names and \
                        p[3] == 'scientific name':
                    name_of[p[0]] = p[1]
                    if len(name_of) == len(need_names):
                        break
    with safe_open(dest, 'wt') as f:
        f.write('organism\ttaxid\tgenus\tfamily\n')
        for org in sorted(org2tid):
            g, fam = t2gf.get(org, ('', ''))
            f.write(f"{org}\t{org2tid[org]}\t"
                    f"{name_of.get(g, '')}\t{name_of.get(fam, '')}\n")
    return n, len(org2tid)


def _search_db_trace(db, engine, faa):
    """搜索库可追溯信息 → (db_files, db_versions)（含蛋白库 FASTA 与
    物种索引；blastp 库是 db + .p?? 一组文件）。"""
    paths = [db, faa, os.path.join(DB_DIR, 'organism_tax.tsv')]
    if engine == 'blastp':
        paths += sorted(glob.glob(db + '.p??'))
    return db_trace(paths, n_max=40)


def _ensure_search_db(faa, engine, exe, logger=None, force=False):
    """按引擎构建/复用搜索库。返回 (引擎侧库路径, db_files, db_versions)。

    库文件清单（文件名+size+mtime）与版本锚点一并维护在 db_info.json
    （可追溯：2026-10-03 审计 P2）。
    """
    meta_f = os.path.join(DB_DIR, 'db_info.json')
    meta = {}
    try:
        with safe_open(meta_f) as f:
            meta = json.load(f)
    except (OSError, ValueError):
        pass
    # 缓存必须校验路径仍然存在：数据库从 databases/viral_prot 迁到
    # databases/annot/prot 后，db_info.json 里的旧绝对路径会让 DIAMOND
    # 直接报「无法打开 -d」，且不重建（2026-09-09 实测踩到）。
    # 另：**随包分发的 db_info.json 里 db_path 是相对 DB_DIR 的文件名**
    # （scripts/package.py 打包时会把打包机的绝对路径净化掉），所以这里
    # 见到相对路径要按当前库目录解析——否则换个安装路径就命不中预置库，
    # 白白重建一次 RefSeq 病毒蛋白 DIAMOND 库。
    cached = meta.get('db_path')
    if cached and not os.path.isabs(cached):
        cached = os.path.join(DB_DIR, cached)
    if (not force and meta.get('engine') == engine
            and meta.get('db_ready') and cached and os.path.exists(cached)):
        # 缓存命中也要有库清单：预置库（随包分发）的 db_info.json 里可能
        # 没有 db_files——现算一次，与已记录不一致才回写，避免每次运行
        # 都改文件 mtime。
        db_files, db_versions = _search_db_trace(cached, engine, faa)
        if (meta.get('db_files') != db_files
                or meta.get('db_versions') != db_versions):
            meta.update({'db_files': db_files, 'db_versions': db_versions})
            try:
                with safe_open(meta_f, 'wt') as f:
                    json.dump(meta, f, ensure_ascii=False, indent=1)
            except OSError:
                pass
        return cached, db_files, db_versions

    if engine == 'diamond':
        db = os.path.join(DB_DIR, 'viral_prot.dmnd')
        if force:
            _rm_if_exists(db)
        if logger:
            logger.log("构建 DIAMOND 蛋白库（RefSeq viral ~72 万条）...")
        run_cmd([exe, 'makedb', '--in', faa, '-d', db, '--quiet'], logger=logger)
        ok = os.path.isfile(db)
    elif engine == 'mmseqs':
        db = os.path.join(DB_DIR, 'mmseqs_db')
        if force:
            _rm_if_exists(db)
        if logger:
            logger.log("构建 MMseqs2 蛋白库（RefSeq viral）...")
        env = os.environ.copy()
        env['PATH'] = os.path.dirname(exe) + os.pathsep + env.get('PATH', '')
        run_cmd([exe, 'createdb', faa, db], logger=logger, env=env)
        ok = os.path.isfile(db)
    else:
        db = os.path.join(_ascii_blast_base(), 'viral_prot')
        if logger:
            logger.log("构建 BLAST 蛋白库（RefSeq viral，库建在 ASCII 临时目录）...")
        run_cmd([exe.replace('blastp', 'makeblastdb') if os.path.isfile(
                     exe.replace('blastp', 'makeblastdb')) else
                 get_config().tool('makeblastdb'),
                 '-dbtype', 'prot', '-in', faa, '-out', db], logger=logger)
        ok = bool(glob.glob(db + '.p??'))
    if not ok:
        raise RuntimeError("搜索库构建失败")
    # db_path 能落在 DB_DIR 内就存**相对文件名**（可随包分发、与安装路径无关）；
    # 只有 blastp 那种建在 ASCII 临时目录的库才存绝对路径。
    try:
        _rel = os.path.relpath(db, DB_DIR)
        db_store = db if _rel.startswith('..') else _rel
    except ValueError:
        db_store = db
    # 库文件清单（文件名+size+mtime）与版本锚点写进 db_info.json：结果表/报告
    # 可回溯"这张注释是用哪个库、哪个版本的库跑出来的"（2026-10-03 审计 P2）。
    db_paths = [db, faa, os.path.join(DB_DIR, 'organism_tax.tsv')]
    if engine == 'blastp':
        db_paths += sorted(glob.glob(db + '.p??'))
    db_files, db_versions = db_trace(db_paths, n_max=40)
    meta.update({'engine': engine, 'db_ready': True, 'db_path': db_store,
                 'db_files': db_files, 'db_versions': db_versions})
    with safe_open(meta_f, 'wt') as f:
        json.dump(meta, f, ensure_ascii=False, indent=1)
    return db, db_files, db_versions


def ensure_viral_prot_db(logger=None, force=False, local_faa=None,
                         progress=None, engine=None):
    """确保 RefSeq 病毒蛋白参考库就绪（下载→解压→物种索引→搜索库）。

    progress(frac 0~1, msg) 反映整个建库过程。返回 db_info dict
    （含 db_files/db_versions 库文件清单与版本锚点，供 summary 记录）。
    engine: 指定搜索引擎；None 自动探测。
    """
    os.makedirs(DB_DIR, exist_ok=True)
    faa = check_path(os.path.join(DB_DIR, 'viral_prot.faa'),
                     must_exist=False, in_platform=True)
    tax_tsv = check_path(os.path.join(DB_DIR, 'organism_tax.tsv'),
                         must_exist=False, in_platform=True)
    eng = annotation_engine(engine)
    if not eng:
        raise RuntimeError("未找到 DIAMOND / MMseqs2 / blastp 任何一个搜索引擎")
    engine, exe = eng

    if local_faa:
        faa = check_path(local_faa, must_exist=True)
    elif not os.path.isfile(faa) or force:
        gz = os.path.join(DB_DIR, 'viral_prot.faa.gz')
        if force:
            _rm_if_exists(gz)
        if not os.path.isfile(gz):
            # RefSeq viral release 的蛋白分卷（viral.1 起，404 即止）
            urls, i = [], 1
            while i <= 8:
                u = REFSEQ_URL.format(i)
                req = urllib.request.Request(u, method='HEAD',
                                             headers={'User-Agent': 'vp'})
                try:
                    urllib.request.urlopen(req, timeout=60)
                    urls.append(u)
                except urllib.error.HTTPError as e:
                    if e.code == 404:
                        break
                    raise
                i += 1
            for k, u in enumerate(urls):
                _download_whitelisted(
                    u, gz, logger=logger,
                    progress=(lambda frac, msg, _k=k, _n=len(urls):
                              progress((_k + frac) / _n * 0.4, '下载蛋白库 ' + msg))
                    if progress else None)
        if logger:
            logger.log("解压蛋白库 ...")
        with safe_open(gz, 'rb') as fin, safe_open(faa, 'wb') as fout:
            while True:
                chunk = fin.read(1 << 22)
                if not chunk:
                    break
                fout.write(chunk)

    if not os.path.isfile(tax_tsv) or force:
        if logger:
            logger.log("构建物种 → 科/属 分类索引（NCBI taxonomy 谱系）...")
        _build_organism_tax(faa, tax_tsv, logger=logger)
    if progress:
        progress(0.6, '构建搜索引擎库')

    db, db_files, db_versions = _ensure_search_db(faa, engine, exe,
                                                  logger=logger, force=force)
    if progress:
        progress(1.0, '参考库就绪')
    return {'faa': faa, 'engine': engine, 'db': db,
            'tax_tsv': tax_tsv, 'n_proteins': None,
            'db_files': db_files, 'db_versions': db_versions}


# ------------------------------------------------------------------
# 功能类别（正则规则按优先级，吸收 annotation_benchmark 的类别口径）
# ------------------------------------------------------------------
_UNKNOWN = '假想蛋白（功能未定）'
_OTHER = '其他功能蛋白'
_CATEGORY_RULES = [
    (r'hypothetical|uncharacterized|unknown protein|duf\d+', _UNKNOWN),
    (r'rna[- ]dependent rna polymerase|rdrp|rna polymerase|replicase'
     r'|polymerase|\bl protein\b|\blarge protein\b|pb1|pb2|pb3|nsp\d+'
     r'|l protein', '聚合酶/复制相关'),
    (r'helicase|primase|methyltransferase|endoribonuclease|exoribonuclease'
     r'|guanylyltransferase|processivity|cap-snatch', '聚合酶辅因子'),
    (r'protease|proteinase|peptidase', '蛋白酶/加工'),
    (r'movement protein', '运动蛋白'),
    (r'capsid|coat|nucleocapsid|nucleoprotein|core protein|virion|envelope'
     r'|glycoprotein|spike|matrix|membrane protein|surface|structural'
     r'|tubule|head|portal|tail|fiber|penton|hexon', '结构蛋白'),
    (r'silenc|suppress|avirulence|virulence|effector|pathogenic|toxin', '宿主互作/致病'),
    (r'integrase|transposase|terminase|terminal protein|origin|rep protein'
     r'|replication|nonstructural|plasmid|ssb|single.stranded dna binding'
     r'|packaging|genome', '基因组复制/维持'),
    (r'transcription|rna binding|transactiv|translation|ribosomal|regulatory'
     r'|regulat', '转录/翻译调控'),
]


def classify_product(product):
    """产物名 → 功能类别（首个命中的规则生效）。"""
    t = (product or '').lower()
    if not t:
        return _UNKNOWN
    for pat, cat in _CATEGORY_RULES:
        if re.search(pat, t):
            return cat
    return _OTHER


def _vote_categories(hits):
    """E-value 加权类别投票（annotation_benchmark 口径）：
    每条命中权重 = -log10(E)（E=0 记 1000），按类别求和；有效类别需 ≥
    最佳无效类别（假想蛋白）权重的 50% 才胜出。"""
    w_inf, w_uninf = {}, {}
    for h in hits:
        try:
            ev = float(h['ev'])
        except (TypeError, ValueError):
            ev = 1.0
        w = 1000.0 if ev <= 0 else -math.log10(max(ev, 1e-180))
        product = _parse_title(h.get('title', ''))[1]
        cat = classify_product(product)
        d = w_uninf if cat == _UNKNOWN else w_inf
        d[cat] = d.get(cat, 0.0) + w
    best_inf = max(w_inf.items(), key=lambda x: x[1])[0] if w_inf else None
    best_un = max(w_uninf.items(), key=lambda x: x[1])[0] if w_uninf else None
    if best_inf is None:
        return best_un or _UNKNOWN
    if best_un is None or w_inf[best_inf] >= 0.5 * w_uninf[best_un]:
        return best_inf
    return best_un


# ------------------------------------------------------------------
# 搜索与解析
# ------------------------------------------------------------------
def _run_search(engine, exe, faa, db, out_tsv, threads, logger=None):
    """三引擎统一跑蛋白搜索。输出列序统一为
    q, acc, pident, qlen, qstart, qend, evalue, bits, title。"""
    if engine == 'diamond':
        run_cmd([exe, 'blastp', '-d', db, '-q', faa, '-o', out_tsv,
                 '--outfmt', '6', 'qseqid', 'sseqid', 'pident', 'qlen',
                 'qstart', 'qend', 'evalue', 'bitscore', 'stitle',
                 '--very-sensitive', '--max-target-seqs', '8', '-e', '1e-5',
                 '--threads', str(threads), '--quiet'], logger=logger)
    elif engine == 'mmseqs':
        import shutil
        env = os.environ.copy()
        env['PATH'] = os.path.dirname(exe) + os.pathsep + env.get('PATH', '')
        tmp = out_tsv + '.mmseqs_tmp'
        try:
            # identity 必须用 OUTFMT_PIDENT（0-100），不能用 fident（0-1）：
            # 下游 pident 列与 DIAMOND/blastp 同量纲比较，用 fident 会让
            # identity 整体缩小 100 倍（2026-10-03 审计 P1）。
            run_cmd([exe, 'easy-search', faa, db, out_tsv, tmp,
                     '--format-output', 'query,target,pident,qstart,qend,'
                                        'evalue,'
                                        'bits,qlen',
                     '-e', '1e-5', '--max-seqs', '8', '-s', '7.5',
                     '--threads', str(threads)], logger=logger, env=env)
        finally:
            # easy-search 的临时工作目录含 Windows 无法 stat 的特殊文件，
            # 残留会让 /api/tool/runs 等目录遍历出问题
            shutil.rmtree(tmp, ignore_errors=True)
    else:
        run_cmd([exe, '-query', faa, '-db', db, '-evalue', '1e-5',
                 '-max_target_seqs', '8', '-num_threads', str(threads),
                 '-outfmt', '6 qseqid sseqid pident qlen qstart qend '
                            'evalue bitscore stitle',
                 '-out', out_tsv], logger=logger)


def _parse_search(engine, path, logger=None):
    """解析搜索引擎输出 → {query: [hit, ...]}（文件已按分数降序）。

    pident 统一按 0-100 量纲（DIAMOND/blastp 原样；mmseqs 用 OUTFMT_PIDENT）。
    防御：mmseqs 命中里出现 ≤1.0 的 identity 视为旧格式串残留的 fident
    （0-1 量纲），×100 修正并告警——否则该命中的 identity 会小 100 倍。
    """
    by_q = {}
    n_scaled = 0
    with safe_open(path) as f:
        for line in f:
            p = line.rstrip('\n').split('\t')
            if engine == 'mmseqs':
                # 列序: query,target,pident,qstart,qend,evalue,bits,qlen
                # （query/target 均为完整头部）
                if len(p) < 8:
                    continue
                q, acc, pid = (p[0].split()[0], p[1].split()[0], p[2])
                qs, qe, ev, bits, qlen, title = \
                    p[3], p[4], p[5], p[6], p[7], p[1]
            else:
                if len(p) < 9:
                    continue
                q, acc, pid, qlen, qs, qe, ev, bits, title = p[:9]
            try:
                h = {'q': q, 'acc': acc, 'pid': float(pid) if pid else 0.0,
                     'qlen': int(float(qlen or 0)),
                     'qs': int(float(qs or 0)), 'qe': int(float(qe or 0)),
                     'ev': float(ev), 'bits': float(bits), 'title': title}
            except ValueError:
                continue
            if engine == 'mmseqs' and 0.0 < h['pid'] <= 1.0:
                h['pid'] *= 100.0
                n_scaled += 1
            by_q.setdefault(h['q'], []).append(h)
    if n_scaled and logger:
        logger.log(f"mmseqs 命中中 {n_scaled} 条 identity ≤1.0，疑似 fident "
                   f"量纲（0-1），已按 ×100 修正为 pident", "WARN")
    return by_q


def _parse_title(title):
    """RefSeq 蛋白头 'ACC 产物名 [物种]' → (acc, product, organism)。"""
    title = (title or '').strip()
    if not title:
        return '', '', ''
    acc = title.split(' ', 1)[0]
    body = title.split(' ', 1)[1] if ' ' in title else ''
    org = ''
    m = re.search(r'\[([^\[\]]+)\]\s*$', body)
    if m:
        org = m.group(1).strip()
        body = body[:m.start()].strip()
    return acc, body, org


def _parse_orf_header(q):
    """识别 ORF 坐标。pyrodigal: 'CONTIG #begin #end strand'（1 起始闭区间，
    与 Virus_Platform_Core/orf.py:182 的写入口径、pyrodigal.gff 一致）；
    orfipy: 'CONTIG_ORF.n [s-e](±)'（1 起始闭区间）。返回 (contig, s, e, strand)。"""
    m = re.match(r'^(.*?)\s+#(\d+)\s+#(\d+)\s+(-?1)\s*$', q)
    if m:
        # 不得再 +1：写入方已是 1-based，+1 会让每条 ORF 的 start 右移 1bp
        # 且区间长度少 1（与同基因的 pyrodigal.gff 对不上）。
        return (m.group(1), int(m.group(2)), int(m.group(3)),
                1 if int(m.group(4)) >= 0 else -1)
    m = re.match(r'^(.*?)_ORF\.\d+\s+\[(\d+)-(\d+)\]\(([+-])\)', q)
    if m:
        return (m.group(1), int(m.group(2)), int(m.group(3)),
                1 if m.group(4) == '+' else -1)
    cid = q.split()[0] if q else 'query'
    return cid, None, None, None


def _load_ictv_tables():
    """ICTV 科级表 → ({family: host_source}, {family: genome_composition})。"""
    def _load(path):
        d = {}
        if os.path.isfile(path):
            with safe_open(path) as f:
                for line in f:
                    p = line.rstrip('\n').split('\t')
                    if len(p) >= 2 and p[0] != 'family':
                        d[p[0].strip()] = p[1].strip()
        return d
    return _load(ICTV_HOST_TSV), _load(ICTV_GENCOMP_TSV)


_GFF_BAD = str.maketrans({c: ' ' for c in ';=,\t\r'})
#: GFF3 规范保留字符 → 百分号转义（% ; = & , 与制表/换行/回车）。
#: 注意未用 urllib.parse.quote(safe='')：那会把中文与空格也转成 %XX，
#: 而 gbdraw CLI 直接消费本 GFF3 渲染 feature 标签（gbdraw_plot.py 把
#: --gff 透传给命令行），标签会变成 %E5%A4%96... 乱码。按保留字符转义
#: 既满足 GFF3 规范，又保持标签可读（2026-10-03 审计 P3）。
_GFF_ESC = str.maketrans({c: '%%%02X' % ord(c) for c in '%;=&,\t\r\n'})


def _gff_clean(v):
    return (v or '').translate(_GFF_BAD).strip()


def _gff_val(v):
    """GFF3 属性值：先 _gff_clean 清洗（分号/等号/逗号/制表 → 空格），
    再对保留字符做百分号转义。"""
    return _gff_clean(v).translate(_GFF_ESC)


# ------------------------------------------------------------------
# 阶段入口
# ------------------------------------------------------------------
def _pick_orf_input(sample_dir, model=None):
    """选择待注释蛋白。

    model: 'pyrodigal_rv' / 'pyrodigal'（显式指定）；None/非法值 → 默认
    pyrodigal_rv 优先。两模型产物都缺时回退 orfipy（仅历史运行兜底，无 gff）。
    返回 (faa 路径, gff 路径或 None, 模型名) 或 (None, None, None)。
    """
    orf_dir = os.path.join(sample_dir, '04_orf')
    if model in ('pyrodigal_rv', 'pyrodigal'):
        cands = [(model, f'{model}.faa', f'{model}.gff')]
    else:
        cands = [
            ('pyrodigal_rv', 'pyrodigal_rv.faa', 'pyrodigal_rv.gff'),
            ('pyrodigal', 'pyrodigal.faa', 'pyrodigal.gff'),
        ]
    cands.append(('orfipy', 'orfipy_pep.fa', None))  # 历史运行兜底
    for model, faa_name, gff_name in cands:
        faa = os.path.join(orf_dir, faa_name)
        if os.path.isfile(faa):
            try:
                n = count_fasta_seqs(faa)
            except (OSError, ValueError):
                n = 0
            if n:
                gff = (os.path.join(orf_dir, gff_name)
                       if gff_name and os.path.isfile(os.path.join(orf_dir,
                                                                   gff_name))
                       else None)
                return faa, gff, model
    return None, None, None


# ------------------------------------------------------------------
# 注释分层选择（prot 序列同源 / pfam HMM / bfvd 病毒域）
# ------------------------------------------------------------------
#: 可选的注释层，顺序即执行顺序与展示顺序。
# 历史层 'cdd'（CDD 结构域，mmseqs2）因本地库 prefilter 索引 ~5GB 反复
# OOM，2026-09-29 移除；同日改接 'bfvd'（BFVD 病毒域，与 SUVTK 提交链
# 同库同命名，见 Virus_Platform_Core/bfvd_search.py）。旧调用传入 'cdd'
# 会按无法识别的层忽略并告警。
ANNOT_LIBS = ('prot', 'pfam', 'bfvd')
LIB_LABELS = {
    'prot': '层1 序列同源（RefSeq 病毒蛋白）',
    'pfam': '层2 Pfam HMM（Pfam-A-Viruses）',
    'bfvd': '层2 BFVD 病毒域（mmseqs2，与 NCBI 提交同源）',
}
#: 兼容别名 → 规范层名（大小写不敏感）。
_LIB_ALIASES = {
    'seq': 'prot', 'homolog': 'prot', 'homology': 'prot', 'blast': 'prot',
    'blastp': 'prot', 'diamond': 'prot', 'mmseqs': 'prot', 'mmseqs2': 'prot',
    'refseq': 'prot', 'protein': 'prot', 'layer1': 'prot', 'l1': 'prot',
    'hmm': 'pfam', 'hmmer': 'pfam', 'pyhmmer': 'pfam', 'pfam-a': 'pfam',
    'pfam_a': 'pfam', 'pfam-a-viruses': 'pfam', 'layer2': 'pfam',
    'domain': 'bfvd', 'domains': 'bfvd', 'bfvd': 'bfvd',
    'bfvd_domain': 'bfvd',
}


def normalize_libs(libs=None):
    """规范化注释层选择，返回 (规范层元组, 无法识别的原始项列表)。

    libs: None / '' / [] → 全部层；'prot,pfam' 字符串；['prot','bfvd'] 列表。
    已识别的层按 ANNOT_LIBS 顺序去重返回；全部无法识别时回退为全部层。
    """
    if libs is None or libs == '' or libs == []:
        return tuple(ANNOT_LIBS), []
    parts = (re.split(r'[,;+|/\s]+', libs) if isinstance(libs, str)
             else list(libs))
    picked, bad = [], []
    for raw in parts:
        k = str(raw).strip().lower()
        if not k:
            continue
        k = _LIB_ALIASES.get(k, k)
        if k in ANNOT_LIBS:
            if k not in picked:
                picked.append(k)
        else:
            bad.append(str(raw).strip())
    if not picked:
        return tuple(ANNOT_LIBS), bad
    return tuple(l for l in ANNOT_LIBS if l in picked), bad


def run_orf_annotation(sample_dir, threads=None, logger=None, force=False,
                       progress=None, min_pident=None, engine=None, db=None,
                       model=None, libs=None, min_qcov=None):
    """⑥b ORF 功能注释主入口。返回 summary dict。

    engine: 指定搜索引擎（'diamond'/'mmseqs'/'blastp'），None 自动探测。
    db: 参考蛋白库。None=默认 RefSeq 病毒蛋白；路径=自备本地蛋白 FASTA。
    model: ORF 基因模型。'pyrodigal_rv'（默认优先）/ 'pyrodigal'；None=rv 优先。
    libs: 启用的注释层，见 ANNOT_LIBS。None/空=全部层（prot+pfam+bfvd）。
          只选 'pfam'/'bfvd' 时不触发 DIAMOND/库下载，纯 HMM/域注释。
    min_pident / min_qcov: 层1 命中的 identity / query 覆盖度过滤门槛，
          百分数（0-100，与输出 pident/qcov 列同量纲）。默认 None=不过滤，
          与历史结果完全一致。搜索引擎侧默认只有 e≤1e-5 一道门槛（DIAMOND
          默认口径），短 ORF / 远端同源会带出低 identity、低覆盖的命中；
          要收紧到"高质量同源"可试验 min_pident=30、min_qcov=50（阈值为
          经验起点，建议按所用库与数据分布标定后再定值）。
    """
    step = 'orfa'
    out_dir = check_path(os.path.join(sample_dir, '04b_orf_annot'),
                         must_exist=False, in_platform=True)
    os.makedirs(out_dir, exist_ok=True)
    summary_file = os.path.join(out_dir, 'summary.json')
    if is_step_done(out_dir, step) and not force:
        if logger:
            logger.log("阶段⑥b ORF 功能注释已完成，跳过")
        with safe_open(summary_file) as f:
            return json.load(f)

    faa, gff, model = _pick_orf_input(sample_dir, model)
    if not faa:
        raise RuntimeError("④ORF 阶段未产出蛋白序列（请先运行 ⑥ORF 预测）")
    try:
        n_orfs = count_fasta_seqs(faa)
    except (OSError, ValueError):
        n_orfs = 0
    sel, bad_libs = normalize_libs(libs)
    do_prot, do_pfam, do_bfvd = ('prot' in sel, 'pfam' in sel, 'bfvd' in sel)
    if bad_libs and logger:
        logger.log("忽略无法识别的注释层: " + ', '.join(bad_libs), "WARN")
    if logger:
        logger.log(f"阶段⑥b ORF 功能注释: {n_orfs} 个 ORF（{model} 基因模型，"
                   f"启用层 {'+'.join(sel)}）")

    # ---- 层1：序列同源搜索（可选） ----
    engine_name = ''
    hits_tsv = ''
    info = None
    by_q = {}
    if do_prot:
        eng = annotation_engine(engine)
        if not eng:
            raise RuntimeError("未找到搜索引擎（DIAMOND/MMseqs2/blastp 任一即可）")
        engine_name, exe = eng
        db_path = check_path(db, must_exist=True, in_platform=False) if db else None

        def _db_prog(frac, msg):
            if progress:
                progress(0.03 + frac * 0.42, '参考库: ' + msg)

        info = ensure_viral_prot_db(logger=logger, progress=_db_prog if progress
                                    else None, local_faa=db_path)
        engine_name, db_ref = info['engine'], info['db']
        if progress:
            progress(0.5, f'{engine_name} 蛋白搜索'
                     + ('（自备库）' if db_path else '（RefSeq viral）'))

        # ---- 搜索 ----
        hits_tsv = check_path(os.path.join(out_dir, 'search_hits.raw.tsv'),
                              must_exist=False, in_platform=True)
        _run_search(engine_name, exe, faa, db_ref, hits_tsv,
                    threads or get_config().threads, logger=logger)
        if progress:
            progress(0.52, '解析命中 + 类别投票')

        # ---- 解析合并 ----
        by_q = _parse_search(engine_name, hits_tsv, logger=logger)
    elif logger:
        logger.log("层1 序列同源搜索按选择跳过")
    # ORF 坐标与所属 contig：搜索引擎 qseqid 只保留首词，pyrodigal/orfipy
    # 的坐标在 faa 头部空白后的描述里。contig 优先取 GFF 的 ID=→seqid
    # 权威映射，无 GFF 时按命名约定剥离基因编号后缀
    gff_map = {}
    if gff:
        id_re = re.compile(r'ID=([^;]+)')
        with safe_open(gff) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                p = line.rstrip('\n').split('\t')
                if len(p) < 9:
                    continue
                m = id_re.search(p[8])
                if m:
                    gff_map[m.group(1)] = p[0]
    orf_meta = {}
    for h, _seq in iter_fasta(faa):
        qid = h.split()[0]
        cid, s0, e, strand = _parse_orf_header(h)
        if cid == qid:
            cid = (re.sub(r'_ORF\.\d+$', '', qid) if '_ORF.' in qid
                   else re.sub(r'_\d+$', '', qid))
        orf_meta[qid] = (gff_map.get(qid) or cid, s0, e, strand)
    org2tax = {}
    if info and info.get('tax_tsv'):
        with safe_open(info['tax_tsv']) as f:
            f.readline()
            for line in f:
                p = line.rstrip('\n').split('\t')
                if len(p) >= 4:
                    org2tax[p[0]] = (p[1], p[2], p[3])
    host_map, gencomp_map = _load_ictv_tables()

    def _sp(frm, to):
        return (lambda p, m: progress(frm + p * (to - frm), m)
                if progress else None)

    rows = []
    # 层1 命中过滤（min_pident/min_qcov，默认 None=不过滤）。在"逐 ORF 建行"
    # 前统一过滤：命中列表被过滤后，类别投票与 top 命中都只看通过门槛的命中。
    def _hit_qcov(h):
        try:
            return (abs(int(h['qe']) - int(h['qs'])) + 1) \
                / max(int(h['qlen']), 1) * 100
        except (ValueError, KeyError, ZeroDivisionError):
            return 0.0

    def _pass_filters(h):
        if min_pident is not None and float(h.get('pid') or 0) < min_pident:
            return False
        if min_qcov is not None and _hit_qcov(h) < min_qcov:
            return False
        return True

    n_hits_filtered = 0
    if min_pident is not None or min_qcov is not None:
        n_before = sum(len(v) for v in by_q.values())
        by_q = {k: [h for h in v if _pass_filters(h)] for k, v in by_q.items()}
        n_hits_filtered = n_before - sum(len(v) for v in by_q.values())
        if logger:
            logger.log(f"层1 命中过滤: min_pident={min_pident} "
                       f"min_qcov={min_qcov} 百分数口径，剔除 "
                       f"{n_hits_filtered}/{n_before} 条命中")

    # 逐 ORF 建行（不是逐命中建行）：层2 HMM/BFVD 兜底只能更新已存在的行，
    # 若只为层1 命中的 ORF 建行，"无 RefSeq 同源命中"的 ORF 就永远得不到
    # 层2 注释——正是层2 存在的意义（实测 10 个 ORF 只出 8 行）。
    total_q = max(len(orf_meta), 1)
    for i, (qid, (cid, s, e, strand)) in enumerate(orf_meta.items()):
        if i % 200 == 0 and progress:
            progress(0.55 + i / total_q * 0.07,
                     f'注释 {i}/{total_q} 个 ORF')
        hits = by_q.get(qid) or []
        row = {
            'orf_id': qid, 'contig': cid, 'start': s, 'end': e,
            'strand': {1: '+', -1: '-', None: ''}[strand],
            'length_aa': 0,
            'hit_accession': '', 'product': '', 'organism': '',
            'pident': '', 'qcov': 0.0, 'bitscore': '', 'evalue': '',
            'genus': '', 'family': '', 'host_source': '',
            'genome_composition': '', 'category': _UNKNOWN,
            'informative': 'N', 'hmm_hits': '', 'bfvd_hits': '',
            'evidence': '',
        }
        if hits:
            top = hits[0]
            acc, product, organism = _parse_title(top['title'])
            tax = org2tax.get(organism, ('', '', ''))
            genus, family = tax[1], tax[2]
            category = _vote_categories(hits)
            try:
                qlen = int(top['qlen'])
                qcov = round((abs(int(top['qe']) - int(top['qs'])) + 1)
                             / max(qlen, 1) * 100, 1)
            except (ValueError, ZeroDivisionError):
                qlen, qcov = 0, 0.0
            row.update({
                'length_aa': qlen,
                'hit_accession': acc, 'product': product,
                'organism': organism, 'pident': top['pid'], 'qcov': qcov,
                'bitscore': top['bits'], 'evalue': top['ev'],
                'genus': genus, 'family': family,
                'host_source': host_map.get(family, ''),
                'genome_composition': gencomp_map.get(family, ''),
                'category': category,
                'informative': 'Y' if category != _UNKNOWN else 'N',
                'evidence': 'seq' if category != _UNKNOWN else '',
            })
        rows.append(row)
    rows.sort(key=lambda r: (r['contig'],
                             r['start'] if r['start'] is not None else 0))

    # ---- 层2：HMM（Pfam 病毒库，viral_fams 双门槛）+ BFVD 病毒域 ----
    # 过滤口径：域级 i-Evalue ≤ 1e-3 且 HMM 模型覆盖率 ≥ 0.5
    # 两层各自独立开关：BFVD 不依赖 HMM 库，反之亦然。
    hmm_by_q, bfvd_by_q = {}, {}
    pfam_libs = []
    if do_pfam:
        if hmm_available() and available_libs():
            if logger:
                logger.log("层2 HMM 扫描（Pfam-A-Viruses）...")
            hmm_by_q, pfam_libs = annotate_orfs_hmm(
                faa, threads=threads, logger=logger,
                progress=(lambda p, m: progress(0.62 + p * 0.2, '层2 ' + m))
                if progress else None)
            pfam_libs = list(pfam_libs or ())
        elif logger:
            logger.log("层2 Pfam HMM 跳过（pyhmmer 未安装或无 HMM 库）", "WARN")
    if do_bfvd:
        if bfvd_db():
            if progress:
                progress(0.84, 'BFVD 病毒域搜索（mmseqs2）')
            bfvd_by_q = bfvd_search_orfs(faa, out_dir, threads=threads,
                                         logger=logger) or {}
        elif logger:
            logger.log("层2 BFVD 病毒域跳过（mmseqs2 或 BFVD 库不可用）", "WARN")
    n_hmm_only, n_bfvd_only = merge_layer2(rows, hmm_by_q, bfvd_by_q,
                                           classify_product)
    # 层2 兜底行可能新给出科：补查 ICTV 宿主/基因组类型
    for r in rows:
        if r.get('family') and not r.get('host_source'):
            r['host_source'] = host_map.get(r['family'], '')
            r['genome_composition'] = gencomp_map.get(r['family'], '')
    if progress:
        progress(0.92, '生成注释表与示意图')

    # ---- 输出 ----
    fields = ['orf_id', 'contig', 'start', 'end', 'strand', 'length_aa',
              'hit_accession', 'product', 'organism', 'pident', 'qcov',
              'bitscore', 'evalue', 'genus', 'family', 'host_source',
              'genome_composition', 'category', 'informative',
              'evidence', 'hmm_hits', 'bfvd_hits']
    with safe_open(os.path.join(out_dir, 'orf_annotation.tsv'), 'wt') as f:
        f.write('\t'.join(fields) + '\n')
        for r in rows:
            f.write('\t'.join(str(r.get(k, '')) for k in fields) + '\n')

    cats, fams = {}, {}
    for r in rows:
        cats[r['category']] = cats.get(r['category'], 0) + 1
        if r['family']:
            fams[r['family']] = fams.get(r['family'], 0) + 1
    with safe_open(os.path.join(out_dir, 'orf_function_summary.tsv'), 'wt') as f:
        f.write('Category\tORFs\n')
        for k in sorted(cats, key=lambda x: -cats[x]):
            f.write(f'{k}\t{cats[k]}\n')
    with safe_open(os.path.join(out_dir, 'orf_family_summary.tsv'), 'wt') as f:
        f.write('Family\tORFs\tHost_source\tGenome_composition\n')
        for k in sorted(fams, key=lambda x: -fams[x]):
            f.write(f'{k}\t{fams[k]}\t{host_map.get(k, "")}\t'
                    f'{gencomp_map.get(k, "")}\n')

    # 逐 contig 功能画像
    by_contig = {}
    for r in rows:
        by_contig.setdefault(r['contig'], []).append(r)
    with safe_open(os.path.join(out_dir, 'contig_function_profile.tsv'),
                   'wt') as f:
        f.write('contig\tORFs\tAnnotated\tCategories\tDominant_family'
                '\tHost_source\tGenome_composition\n')
        for cid in sorted(by_contig, key=lambda c: -len(by_contig[c])):
            rs = by_contig[cid]
            cc = {}
            for r in rs:
                if r['informative'] == 'Y':
                    cc[r['category']] = cc.get(r['category'], 0) + 1
            ff = {}
            for r in rs:
                if r['family']:
                    ff[r['family']] = ff.get(r['family'], 0) + 1
            dom = max(ff.items(), key=lambda x: x[1])[0] if ff else ''
            cats_txt = '，'.join(f'{k} {v}' for k, v in
                                 sorted(cc.items(), key=lambda x: -x[1]))
            f.write(f"{cid}\t{len(rs)}\t{sum(cc.values())}\t{cats_txt}\t"
                    f"{dom}\t{host_map.get(dom, '')}\t"
                    f"{gencomp_map.get(dom, '')}\n")

    # 注释写回 GFF3（pyrodigal 模型时）
    annot_by_id = {r['orf_id']: r for r in rows}
    n_gff = 0
    if gff:
        with safe_open(gff) as fin, \
                safe_open(os.path.join(out_dir, 'orf_annotation.gff3'),
                          'wt') as fout:
            fout.write('##gff-version 3\n')
            for line in fin:
                if line.startswith('#'):
                    continue
                p = line.rstrip('\n').split('\t')
                if len(p) < 9:
                    continue
                m = re.search(r'ID=([^;]+)', p[8])
                r = annot_by_id.get(m.group(1)) if m else None
                if r:
                    attrs = (f"{p[8]};product={_gff_val(r['product'])}"
                             f";organism={_gff_val(r['organism'])}"
                             f";category={_gff_val(r['category'])}"
                             f";family={_gff_val(r['family'])}"
                             f";evidence={_gff_val(r.get('evidence', ''))}"
                             + (f";hmm={_gff_val(r['hmm_hits'])}"
                                if r.get('hmm_hits') else '')
                             + (f";bfvd={_gff_val(r['bfvd_hits'])}"
                                if r.get('bfvd_hits') else ''))
                    n_gff += 1
                else:
                    attrs = p[8]
                fout.write('\t'.join(p[:8] + [attrs]) + '\n')

    n_annot = sum(1 for r in rows if r['informative'] == 'Y')
    # 基因组示意图（littlegenomes 风格适配：ORF 类别着色箭头 + 结构域窄条）
    contig_lengths = {}
    vfa = os.path.join(sample_dir, '03_assembly', 'viral_contigs.fasta')
    if os.path.isfile(vfa):
        try:
            for h, s in iter_fasta(vfa):
                contig_lengths[h.split()[0]] = len(s)
        except (OSError, ValueError):
            pass
    diagrams = render_sample_diagrams(rows, out_dir, contig_lengths,
                                      max_plots=12, logger=logger)

    ev_counts = {}
    for r in rows:
        ev = r.get('evidence') or ('seq' if r['informative'] == 'Y' else 'none')
        ev_counts[ev] = ev_counts.get(ev, 0) + 1

    # 库版本可追溯（2026-10-03 审计 P2）：层1 搜索库（RefSeq 蛋白）的文件
    # 清单/版本锚点来自 db_info.json；层2 补 HMM（Pfam press 派生件）与
    # BFVD（bfvd.version）——结果表可回溯"用哪个库、哪版跑的"。
    db_files = dict((info or {}).get('db_files') or {})
    db_versions = dict((info or {}).get('db_versions') or {})
    if pfam_libs:
        try:
            _f, _v = lib_trace()
            db_files.update(_f)
            db_versions.update(_v)
        except Exception as e:
            if logger:
                logger.log(f"HMM 库版本记录失败: {e}", "WARN")
    if do_bfvd and bfvd_db():
        try:
            _f, _v = bfvd_db_trace()
            db_files.update(_f)
            db_versions.update(_v)
        except Exception as e:
            if logger:
                logger.log(f"BFVD 库版本记录失败: {e}", "WARN")

    summary = {
        'stage': step,
        'model': model,
        'engine': engine_name,
        'libs': list(sel),
        'libs_labels': [LIB_LABELS[l] for l in sel],
        'pfam_libs': pfam_libs,
        'filters': {'min_pident': min_pident, 'min_qcov': min_qcov,
                    'n_hits_filtered': n_hits_filtered},
        'db_files': db_files,
        'db_versions': db_versions,
        'n_orfs': len(rows),
        'n_annotated': n_annot,
        'n_families': len(fams),
        'gff_annotated': n_gff,
        'n_hmm_only': n_hmm_only,
        'n_bfvd_only': n_bfvd_only,
        'evidence': dict(sorted(ev_counts.items())),
        'categories': dict(sorted(cats.items(), key=lambda x: -x[1])),
        'families': dict(sorted(fams.items(), key=lambda x: -x[1])),
        'tsv': 'orf_annotation.tsv',
        'diagrams': diagrams,
        'rows_preview': [
            {k: r[k] for k in ('orf_id', 'product', 'organism', 'family',
                               'category', 'pident', 'qcov', 'bitscore',
                               'evidence')}
            for r in sorted((r for r in rows if r['informative'] == 'Y'),
                            key=_preview_rank)[:20]
        ],
    }
    with safe_open(summary_file, 'wt') as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)
    _rm_if_exists(hits_tsv)
    mark_step_done(out_dir, step)
    if progress:
        progress(1.0, 'ORF 功能注释完成')
    if logger:
        logger.log(f"阶段⑥b 完成: {n_annot}/{len(rows)} 个 ORF 获得有效功能注释，"
                   f"覆盖 {len(fams)} 个病毒科")
    return summary


# ------------------------------------------------------------------
# 图表（报告嵌入）
# ------------------------------------------------------------------
def fig_category_bar(summary, top_n=12):
    """功能类别分布条形图。"""
    import plotly.graph_objects as go
    cats = summary.get('categories') or {}
    items = [(k, v) for k, v in cats.items()][:top_n][::-1]
    if not items:
        return None
    palette = {c: '#95a5a6' for c in (_UNKNOWN,)}
    fig = go.Figure(go.Bar(
        x=[v for _, v in items],
        y=[k for k, _ in items],
        orientation='h',
        text=[f"{v}" for _, v in items],
        textposition='auto',
        marker_color=[palette.get(k, '#2e86c1') for k, _ in items]))
    fig.update_layout(title=f'ORF 功能类别分布（共 {summary.get("n_orfs", 0)} 个）',
                      height=420, margin=dict(l=200))
    return fig


def fig_family_bar(summary, top_n=12):
    """病毒科级分布条形图。"""
    import plotly.graph_objects as go
    fams = summary.get('families') or {}
    items = list(fams.items())[:top_n][::-1]
    if not items:
        return None
    fig = go.Figure(go.Bar(
        x=[v for _, v in items],
        y=[k for k, _ in items],
        orientation='h',
        text=[f"{v}" for _, v in items],
        textposition='auto',
        marker_color='#27ae60'))
    fig.update_layout(title=f'ORF 命中病毒科分布（Top {len(items)}）',
                      height=420, margin=dict(l=200))
    return fig
