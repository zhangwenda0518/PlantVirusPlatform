# -*- coding: utf-8 -*-
"""本地（离线）序列比对引擎：BLASTN / DIAMOND BLASTX。

与 Virus_Platform_Core/contig_annot.py 的在线 NCBI 版本输出**同一套命中结构**，前端「命中表」和
结构域图无需区分数据来源。差别只在：

  - 在线：NCBI nt / nr 全库，需联网，分钟级；
  - 本地：平台内置库，完全离线、秒级，可批量。

使用的库：
  - 核酸：databases/virusref_db/final.cluster.ref.fasta → **预置 v4 BLAST 库**
          `databases/virusref_db/blast/virus.*`（随包分发，安装后零构建；v4 格式
          在含中文的路径下也能用 —— v5 的 LMDB 后端不行。见
          Virus_Platform_Core/assembly.py:ensure_virus_blast_db）
  - 蛋白：databases/annot/prot/viral_prot.dmnd（DIAMOND，RefSeq 病毒蛋白）
  - 元数据：databases/virusref_db/final.cluster.ref_info.tsv（accession → 物种/科/属/宿主/GenBank 标题）

> 历史：曾有 mmseqs2 CDD 保守域搜索（本地库 databases/annot/cdd/cdd_db，
> translated search 等价 RPS-BLAST）。因 CDD prefilter 索引一次性分配
> ~5GB、小内存机器反复 OOM，2026-09-29 随库一并移除。
"""
from __future__ import annotations

import csv
import os

from .config import get_config, db_path, ref_annotation_path
from .utils import check_path, run_cmd, safe_open

BLAST_EVALUE_MAX = 1e-5
MAX_HITS = 50

REF_INFO_TSV = ref_annotation_path()
PROT_DMND = os.path.join(db_path('annot', 'prot'), 'viral_prot.dmnd')


def refresh_paths():
    """重算模块级路径常量（切换「设置 → 数据库目录」后由 config 回调）。"""
    global REF_INFO_TSV, PROT_DMND
    global _REF_INFO_CACHE
    REF_INFO_TSV = ref_annotation_path()
    PROT_DMND = os.path.join(db_path('annot', 'prot'), 'viral_prot.dmnd')
    # 缓存也要失效：它们按旧库内容构建
    _REF_INFO_CACHE = None
    return REF_INFO_TSV


_REF_INFO_CACHE = None


# ------------------------------------------------------------------ 元数据
def ref_info() -> dict:
    """{accession: {title, sciname, family, genus, host, length}}（无则空 dict）。"""
    global _REF_INFO_CACHE
    if _REF_INFO_CACHE is not None:
        return _REF_INFO_CACHE
    out = {}
    if os.path.isfile(REF_INFO_TSV):
        with safe_open(REF_INFO_TSV) as f:
            for row in csv.DictReader(f, delimiter='\t'):
                acc = (row.get('Accession') or '').strip()
                if not acc:
                    continue
                out[acc] = {
                    'title': (row.get('GenBank_Title') or '').strip(),
                    'sciname': ((row.get('Species_ICTV')
                                 or row.get('Species_NCBI') or '').strip()),
                    'family': (row.get('VMR_Family') or '').strip(),
                    'genus': (row.get('VMR_Genus') or '').strip(),
                    'host': (row.get('Host') or '').strip(),
                    'length': (row.get('Length') or '').strip(),
                }
    _REF_INFO_CACHE = out
    return out


def _short_desc(text, limit=160):
    text = ' '.join((text or '').split())
    return text if len(text) <= limit else text[:limit - 1] + '…'


# ------------------------------------------------------------------ 可用性
def engine_status() -> dict:
    """本地引擎可用性（供设置页/接口自检）。"""
    cfg = get_config()
    out = {'nuc_db': None, 'nuc_db_built': False, 'prot_db': PROT_DMND,
           'ref_info': REF_INFO_TSV}
    for key, tool in (('blastn', 'blastn'), ('blastx', 'diamond'),
                      ('makeblastdb', 'makeblastdb')):
        try:
            out[key] = cfg.tool(tool)
        except FileNotFoundError:
            out[key] = None
    from .assembly import find_virus_ref_fasta
    try:
        out['nuc_src'] = find_virus_ref_fasta()
    except FileNotFoundError:
        out['nuc_src'] = None
    out['ok'] = bool(out.get('blastn') and out.get('nuc_src'))
    out['ok_blastx'] = bool(out.get('blastx') and os.path.isfile(PROT_DMND))
    # BFVD 域库（③b 验证 B 路与 ⑥b 注释层共用，与 SUVTK 提交链同库）
    try:
        from .bfvd_search import bfvd_db
        out['bfvd_db'] = bfvd_db()
    except Exception:
        out['bfvd_db'] = None
    try:
        out['mmseqs'] = cfg.tool('mmseqs')
    except FileNotFoundError:
        out['mmseqs'] = None
    out['ok_bfvd'] = bool(out.get('mmseqs') and out.get('bfvd_db'))
    return out


def _nuc_db(logger=None) -> str:
    """病毒参考核酸 BLAST 库前缀（自动构建/复用）。"""
    from .assembly import ensure_virus_blast_db
    return ensure_virus_blast_db(logger=logger)


def _run(cmd, env=None, logger=None):
    run_cmd(cmd, env=env, logger=logger)


# ------------------------------------------------------------------ BLASTN
def blastn_local(query_fa, out_dir, max_hits=MAX_HITS, evalue=BLAST_EVALUE_MAX,
                 threads=None, logger=None, per_query_top1=False) -> list:
    """本地 blastn：核酸 query vs 病毒参考核酸库。返回命中列表。

    per_query_top1=True 时先按 query 分组取各自 bitscore 最高的 top1，再
    全局排序截断到 max_hits——批量"每 contig 一条最近参考"场景（见
    contig_blast_table）必须这样：直接全局截断时，命中密集的前几个 contig
    会吃满 max_hits，后面的 contig 一条命中都拿不到（2026-10-03 审计 P2）。
    默认 False 保留每 query 的多条命中（交互式单序列的本地命中表需要
    展示多个参考命中）。
    """
    cfg = get_config()
    blastn = cfg.tool('blastn')
    db = _nuc_db(logger)
    out6 = os.path.join(out_dir, 'blastn_local.m8')
    _run([blastn, '-query', query_fa, '-db', db, '-outfmt',
          '6 qseqid sseqid pident length qlen slen evalue bitscore',
          '-max_target_seqs', str(max_hits), '-evalue', str(evalue),
          '-num_threads', str(threads or cfg.threads), '-out', out6],
         logger=logger)
    info = ref_info()
    best = {}
    if os.path.isfile(out6):
        with safe_open(out6) as f:
            for line in f:
                p = line.rstrip('\n').split('\t')
                if len(p) < 8:
                    continue
                acc = p[1]
                try:
                    bs = float(p[7])
                except ValueError:
                    bs = 0.0
                # 去重键必须含 query：多序列查询时只按 accession 去重会把
                # 其它 contig 的命中整条丢掉（返回字段里的 query 就错了）。
                key = (p[0], acc)
                if key in best and best[key]['bit_score'] >= bs:
                    continue
                meta = info.get(acc, {})
                best[key] = {
                    'query': p[0], 'accession': acc,
                    'title': meta.get('title') or acc,
                    'sciname': meta.get('sciname', ''),
                    'family': meta.get('family', ''),
                    'genus': meta.get('genus', ''),
                    'host': meta.get('host', ''),
                    'identity': float(p[2]) if p[2] else 0,
                    'align_len': int(p[3]) if p[3] else 0,
                    'qlen': int(p[4]) if p[4] else 0,
                    'slen': int(p[5]) if p[5] else 0,
                    'evalue': p[6], 'bit_score': bs,
                }
    hits = _truncate_hits(best.values(), max_hits, per_query_top1)
    if logger:
        logger.log(f'  [blastn·本地] {len(hits)} 条命中（库 {os.path.basename(db)}）')
    return hits


def _truncate_hits(hits, max_hits, per_query_top1):
    """全局截断前的汇总步骤。

    per_query_top1=True：按 query 分组取各自 bitscore 最高的一条，再全局
    排序截断——保证 max_hits ≥ query 数时每个 query 都有命中（批量
    contig 场景）；False：原样全局排序截断（保留每 query 多命中）。
    """
    if per_query_top1:
        top1 = {}
        for h in hits:
            q = h.get('query') or ''
            cur = top1.get(q)
            if cur is None or h['bit_score'] > cur['bit_score']:
                top1[q] = h
        hits = top1.values()
    return sorted(hits, key=lambda h: -h['bit_score'])[:max_hits]


def contig_blast_table(contigs_fa, out_dir, threads=None, logger=None,
                       evalue=BLAST_EVALUE_MAX, write_tsv=True):
    """contig 级"最近参考"：本地 blastn vs 病毒参考核酸库。

    供 `03_assembly/virus_contigs.tsv` 的 blast_* 六列使用——这六列此前
    只有列名、没有任何生产写入方（host_analysis 的「BLAST 最后回退」
    因此从未真正生效；实测 212 条 host_prediction 里 class_source=blast 为 0）。
    同时落盘 `<out_dir>/contig_blast.tsv`，供报告与排查。

    注意：每个 contig 只取自己的 top1 命中（per_query_top1=True），因此
    max_hits 传 max(contig 数, 默认值) 才保证覆盖全部 contig；若只传默认
    50，第 51 个之后的 contig 仍会零命中（此前 docstring 误以为"max_hits
    覆盖 contig 数"就够——旧的全局截断下并不够，2026-10-03 审计 P2）。

    返回 (table, tsv_path|None)；引擎或参考库缺失时返回 ({}, None) 优雅降级。
    """
    from .utils import count_fasta_seqs
    table, tsv_path = {}, None
    if not contigs_fa or not os.path.isfile(contigs_fa):
        return table, None
    try:
        n_ctg = count_fasta_seqs(contigs_fa) or 0
    except (OSError, ValueError):
        n_ctg = 0
    if n_ctg <= 0:
        return table, None
    try:
        # 每个 contig 只留自己的 top1（per_query_top1=True），max_hits 至少
        # 覆盖 contig 数即保证无 contig 零命中
        hits = blastn_local(contigs_fa, out_dir, max_hits=max(n_ctg, MAX_HITS),
                            evalue=evalue, threads=threads, logger=logger,
                            per_query_top1=True)
    except FileNotFoundError as e:
        if logger:
            logger.log(f'  本地 BLASTN 不可用（跳过 blast 列）: {e}', 'WARN')
        return {}, None
    except Exception as e:                       # 库损坏/比对失败都不该中断主流程
        if logger:
            logger.log(f'  本地 BLASTN 失败（跳过 blast 列）: {e}', 'WARN')
        return {}, None

    rows = []
    for h in hits:
        cid = h.get('query') or ''
        if not cid or cid in table:              # blastn_local 已按 bitscore 降序
            continue
        qlen = h.get('qlen') or 0
        aln = h.get('align_len') or 0
        sp = (h.get('sciname') or '').strip()
        if not sp:
            # ref_info 无物种名时退回 GenBank 标题里的 [物种]
            _t, sp2 = _split_stitle(h.get('title') or '', h.get('accession') or '')
            sp = sp2
        table[cid] = {
            'accession': h.get('accession') or '',
            'identity': h.get('identity') or 0,
            'coverage_hsp': round(aln * 100.0 / qlen, 1) if qlen else 0,
            'aln_len': aln,
            'species': sp,
            'family': (h.get('family') or '').strip(),
        }
        rows.append((cid, table[cid]))

    if write_tsv and rows:
        try:
            tsv_path = check_path(os.path.join(out_dir, 'contig_blast.tsv'),
                                  must_exist=False, in_platform=True)
            with safe_open(tsv_path, 'wt') as f:
                f.write('contig\tblast_top_hit\tblast_identity(%)\t'
                        'blast_coverage_hsp(%)\tblast_aln_len\t'
                        'blast_species\tblast_family\n')
                for cid, v in rows:
                    f.write(f"{cid}\t{v['accession']}\t{v['identity']}\t"
                            f"{v['coverage_hsp']}\t{v['aln_len']}\t"
                            f"{v['species']}\t{v['family']}\n")
        except (OSError, ValueError):
            tsv_path = None
    if logger:
        logger.log(f'  [blastn·本地] contig 级最近参考 {len(table)}/{n_ctg} 条')
    return table, tsv_path


# ------------------------------------------------------------------ BLASTX
def _split_stitle(stitle: str, acc: str):
    """'YP_010806246.1 P2-RdRp [Pepper yellows virus]' → (标题, 物种)。"""
    s = (stitle or '').strip()
    if s.startswith(acc):
        s = s[len(acc):].strip()
    sciname = ''
    if s.endswith(']') and '[' in s:
        s, sciname = s.rsplit('[', 1)
        sciname = sciname[:-1].strip()
        s = s.strip()
    return s or acc, sciname


def blastx_local(query_fa, out_dir, max_hits=MAX_HITS, evalue=BLAST_EVALUE_MAX,
                 threads=None, logger=None, per_query_top1=False) -> list:
    """本地 DIAMOND blastx：核酸 query → 蛋白命中 vs RefSeq 病毒蛋白库。

    per_query_top1 语义同 blastn_local（按 query 取 top1 再全局截断）。
    """
    cfg = get_config()
    diamond = cfg.tool('diamond')
    if not os.path.isfile(PROT_DMND):
        raise FileNotFoundError(f'缺少 DIAMOND 病毒蛋白库: {PROT_DMND}')
    out6 = os.path.join(out_dir, 'blastx_local.m8')
    _run([diamond, 'blastx', '-q', query_fa, '--db', PROT_DMND,
          '--long-reads', '-o', out6, '-f', '6', 'qseqid', 'sseqid', 'pident',
          'length', 'qlen', 'slen', 'evalue', 'bitscore', 'stitle',
          '--max-target-seqs', str(max_hits), '-e', str(evalue),
          '--threads', str(threads or cfg.threads)], logger=logger)
    best = {}
    if os.path.isfile(out6):
        with safe_open(out6) as f:
            for line in f:
                p = line.rstrip('\n').split('\t')
                if len(p) < 8:
                    continue
                acc = p[1]
                try:
                    bs = float(p[7])
                except ValueError:
                    bs = 0.0
                key = (p[0], acc)
                if key in best and best[key]['bit_score'] >= bs:
                    continue
                title, sciname = _split_stitle(p[8] if len(p) > 8 else '', acc)
                best[key] = {
                    'query': p[0], 'accession': acc, 'title': title,
                    'sciname': sciname,
                    'identity': float(p[2]) if p[2] else 0,
                    'align_len': int(p[3]) if p[3] else 0,
                    'qlen': int(p[4]) if p[4] else 0,
                    'slen': int(p[5]) if p[5] else 0,
                    'evalue': p[6], 'bit_score': bs,
                }
    hits = _truncate_hits(best.values(), max_hits, per_query_top1)
    if logger:
        logger.log(f'  [blastx·本地] {len(hits)} 条命中（DIAMOND vs viral_prot）')
    return hits
