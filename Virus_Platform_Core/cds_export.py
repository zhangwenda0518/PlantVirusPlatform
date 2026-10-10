# -*- coding: utf-8 -*-
"""CDS / PEP 提取产物模块（独立）。

职责边界：
  - 本模块只负责「从 GenBank 集合里挑 CDS → 按基因名归组 → 导出 CDS/PEP」
    这一条业务线，包含明细表构建（build_cds_table）、人工挑选状态持久化
    （save_selection / load_selection）、按基因归组导出（export_selected）。
  - GenBank 集合的目录布局与解析基础设施仍归 Virus_Platform_Core.gb_collection（集合下载、
    清单、提取），本模块单向依赖它，不反向被依赖，无循环。

产物布局（落在集合的 extract/selected/ 下）：
    CDS.fa            全部选中 CDS 核酸（汇总）
    PEP.fa            全部选中 CDS 蛋白（汇总）
    CDS/<基因>.fa     按基因名归组的核酸（同源基因每基因组一条）
    PEP/<基因>.fa     按基因名归组的蛋白
    selected.tsv      选中明细

序列 id = <accession>|<基因名>；蛋白优先取 GenBank translation，
缺失则按 transl_table 翻译。
"""

import json
import os
import re
import time

from .gb_collection import (_GB_EXTS, _gene_tag_of, _taxonomy_ranks,
                           extract_dir, gb_collection_dir)
from .utils import check_path, safe_open

__all__ = ['build_cds_table', 'selection_path', 'save_selection',
           'load_selection', 'export_selected', 'normalize_gene']


def normalize_gene(name):
    """基因名归一：只留字母数字下划线点连字符，截断 40 字符，空则 gene。"""
    return re.sub(r'[^A-Za-z0-9_\-.]+', '_', str(name or '').strip())[:40]


def _gb_signature(name):
    """集合 GenBank 文件指纹：[(文件名, size, mtime_ns), ...] 排序后取摘要。

    缓存失效判定用「文件集合 + 大小 + mtime」：集合增删记录、重下集合、
    导入新 gb 都会改变它；只是重新打开页面则完全命中。
    """
    import hashlib
    cdir = gb_collection_dir(name)
    items = []
    try:
        for fn_ in sorted(os.listdir(cdir)):
            if not fn_.lower().endswith(_GB_EXTS):
                continue
            p = os.path.join(cdir, fn_)
            try:
                st = os.stat(p)
            except OSError:
                continue
            items.append((fn_, st.st_size, st.st_mtime_ns))
    except OSError:
        return ''
    h = hashlib.sha256(repr(items).encode('utf-8')).hexdigest()[:24]
    return h


def _table_cache_path(name):
    return check_path(os.path.join(gb_collection_dir(name),
                                   'cds_table_cache.json'),
                      must_exist=False, in_platform=True)


def _iter_gb_records(path, logger=None, warn=None):
    """逐 record 读 GenBank 文件（**逐文件容错**，参照 ictv_db._rebuild_gb_fasta）。

    SeqIO.parse 是生成器：单条坏记录会让整个集合的解析中断、该集合所有
    产物为空。这里包一层 try + WARN——坏文件跳过，其余文件继续
    （2026-10-03 审计 P3）。warn 给定时失败原因一并收集（无 logger 的
    调用方如 Web 接口也能把告警带回响应）。
    """
    from Bio import SeqIO
    try:
        for rec in SeqIO.parse(path, 'genbank'):
            yield rec
    except Exception as e:
        msg = (f'{os.path.basename(path)} 解析失败，跳过: '
               f'{type(e).__name__}: {e}')
        if logger:
            logger.log(msg, 'WARN')
        if warn is not None:
            warn.append(msg)


def build_cds_table(name, use_cache=True, logger=None):
    """集合 .gb → 供人工挑选的 CDS 明细表。

    返回 {'rows': [...], 'viruses': [...], 'families': [...], 'genera': [...]}；
    每行含：rid（行级稳定 id，勾选/存盘用）、acc、virus、family、genus、
    gene（归一后基因名）、product、start、end、strand、length、
    has_translation（False 则导出时按 transl_table 翻译）、protein_id。
    坐标一律 1-based 闭区间，与 GenBank 惯例一致。

    use_cache=True 时把结果缓存到集合目录 cds_table_cache.json，键为 GB 文件
    指纹。为什么必须缓存：本函数对集合里**每个** .gb 做 SeqIO.parse 并对
    每条 CDS 调 feat.extract()，一个科级集合可达数百 MB / 数千条记录，
    解析动辄数十秒；而前端「解析并挑选」按钮与页面重开都会调它，原实现
    每次都全量重解析（用户感受就是点一下卡半天）。
    """
    from Bio import SeqIO

    cdir = gb_collection_dir(name)
    sig = _gb_signature(name) if use_cache else ''
    cpath = None
    if use_cache:
        try:
            cpath = _table_cache_path(name)
            if sig and os.path.isfile(cpath):
                with safe_open(cpath) as f:
                    cached = json.load(f)
                if cached.get('sig') == sig and 'rows' in cached:
                    out = {k: v for k, v in cached.items() if k != 'sig'}
                    out['cached'] = True
                    return out
        except (OSError, ValueError):
            cpath = None

    gbs = [fn_ for fn_ in sorted(os.listdir(cdir))
           if fn_.lower().endswith(_GB_EXTS)]
    rows = []
    viruses = {}          # virus -> {'family':..., 'genus':..., 'n': int}
    warns = []            # 逐文件解析失败（返回给调用方，见 _iter_gb_records）
    for fn_ in gbs:
        for rec in _iter_gb_records(os.path.join(cdir, fn_), logger, warns):
            family, genus = _taxonomy_ranks(rec)
            virus = str(rec.annotations.get('organism') or '').strip() or rec.id
            acc = rec.id
            vst = viruses.setdefault(
                virus, {'virus': virus, 'family': family, 'genus': genus,
                        'acc': acc, 'n': 0})
            if not vst['family'] and family:
                vst['family'] = family
            if not vst['genus'] and genus:
                vst['genus'] = genus
            for fi, feat in enumerate(rec.features):
                if feat.type != 'CDS':
                    continue
                q = feat.qualifiers
                product = next((str(v) for v in q.get('product', [''])),
                               '').strip()
                tag = _gene_tag_of(q, product)
                try:
                    loc = feat.location
                    start = int(loc.start) + 1
                    end = int(loc.end)
                    strand = '+' if loc.strand and loc.strand > 0 else '-'
                    if loc.strand is None:
                        strand = '.'
                except Exception:
                    start = end = 0
                    strand = '.'
                nt = str(feat.extract(rec.seq)).upper()
                tr = next((str(v) for v in q.get('translation', [])), None)
                rows.append({
                    'rid': f'{acc}#{fi}',
                    'acc': acc,
                    'virus': virus,
                    'family': family,
                    'genus': genus,
                    'gene': tag,
                    'product': product,
                    'start': start,
                    'end': end,
                    'strand': strand,
                    'length': len(nt),
                    'has_translation': bool(tr),
                    'protein_id': next((str(v) for v in q.get('protein_id', [])),
                                       ''),
                })
                vst['n'] += 1
    vlist = sorted(viruses.values(), key=lambda v: v['virus'].lower())
    out = {
        'rows': rows,
        'viruses': vlist,
        'families': sorted({v['family'] for v in vlist if v['family']}),
        'genera': sorted({v['genus'] for v in vlist if v['genus']}),
        'n_cds': len(rows),
        'n_virus': len(vlist),
        'warnings': warns,
        'cached': False,
    }
    if use_cache and cpath and sig:
        try:
            with safe_open(cpath, 'wt') as f:
                json.dump({**out, 'sig': sig}, f, ensure_ascii=False)
        except (OSError, ValueError):
            pass
    return out


def selection_path(name):
    """集合的人工挑选状态文件 selection.json。"""
    return check_path(os.path.join(gb_collection_dir(name), 'selection.json'),
                      must_exist=False, in_platform=True)


def save_selection(name, items, logger=None):
    """保存人工挑选结果。

    items: [{'rid','gene',...}]，仅保留 rid + gene（gene 为人工确认/改名后的值）。
    返回落盘 dict。
    """
    def log(msg, level='INFO'):
        if logger:
            logger.log(msg, level)

    cdir = gb_collection_dir(name)
    if not os.path.isdir(cdir):
        raise FileNotFoundError(f'集合不存在: {name}')
    keep = []
    for it in items or []:
        rid = str((it or {}).get('rid') or '').strip()
        if not rid:
            continue
        gene = normalize_gene((it or {}).get('gene'))
        keep.append({'rid': rid, 'gene': gene or 'gene'})
    data = {'name': name, 'n': len(keep), 'items': keep,
            'date': time.strftime('%Y-%m-%d %H:%M:%S')}
    with safe_open(selection_path(name), 'wt') as f:
        json.dump(data, f, ensure_ascii=False, indent=2)
    log(f'挑选状态已保存：{len(keep)} 条 → selection.json')
    return data


def load_selection(name):
    """读取人工挑选状态，无则返回 {'items': []}。"""
    p = selection_path(name)
    if not os.path.isfile(p):
        return {'name': name, 'n': 0, 'items': []}
    try:
        with safe_open(p) as f:
            data = json.load(f)
        data.setdefault('items', [])
        data['n'] = len(data['items'])
        return data
    except (OSError, ValueError):
        return {'name': name, 'n': 0, 'items': []}


def export_selected(name, items, logger=None, prog=None, subdir='selected'):
    """按人工挑选导出 CDS + PEP，按「基因名」归组。

    items: [{'rid','gene'}]，gene 为人工确认后的基因名（决定归组）。
    产物：
      extract/selected/CDS.fa        全部选中 CDS 核酸（汇总）
      extract/selected/PEP.fa        全部选中 CDS 蛋白（汇总）
      extract/selected/CDS/<基因>.fa 按基因名分组核酸
      extract/selected/PEP/<基因>.fa 按基因名分组蛋白
      extract/selected/selected.tsv  选中明细（rid/acc/virus/gene/product/坐标）
    蛋白优先 GenBank translation，缺则按 transl_table 翻译。
    返回 summary dict。
    """
    from Bio import SeqIO
    from Bio.Seq import Seq as BioSeq
    from .utils import write_fasta_record

    def log(msg, level='INFO'):
        if logger:
            logger.log(msg, level)

    cdir = gb_collection_dir(name)
    if not os.path.isdir(cdir):
        raise FileNotFoundError(f'集合不存在: {name}')

    want = {}
    for it in items or []:
        rid = str((it or {}).get('rid') or '').strip()
        if not rid:
            continue
        gene = normalize_gene((it or {}).get('gene'))
        want[rid] = gene or 'gene'
    if not want:
        raise ValueError('未选择任何 CDS')

    out = os.path.join(extract_dir(name), subdir)
    cds_dir = os.path.join(out, 'CDS')
    pep_dir = os.path.join(out, 'PEP')
    for d in (out, cds_dir, pep_dir):
        os.makedirs(d, exist_ok=True)

    gbs = [fn_ for fn_ in sorted(os.listdir(cdir))
           if fn_.lower().endswith(_GB_EXTS)]
    if prog:
        prog('export', 0.05, f'扫描 {len(gbs)} 个 GenBank 记录')

    cds_by_gene, pep_by_gene = {}, {}
    rows, hit = [], set()
    n_cds = n_pep = 0
    for i, fn_ in enumerate(gbs, 1):
        for rec in _iter_gb_records(os.path.join(cdir, fn_), logger):
            family, genus = _taxonomy_ranks(rec)
            virus = str(rec.annotations.get('organism') or '').strip() or rec.id
            for fi, feat in enumerate(rec.features):
                if feat.type != 'CDS':
                    continue
                rid = f'{rec.id}#{fi}'
                if rid not in want:
                    continue
                hit.add(rid)
                gene = want[rid]
                q = feat.qualifiers
                product = next((str(v) for v in q.get('product', [''])),
                               '').strip()
                sid = f'{rec.id}|{gene}'
                nt = str(feat.extract(rec.seq)).upper()
                tr = next((str(v) for v in q.get('translation', [])), None)
                if tr:
                    aa = tr.upper()
                else:
                    table = next((str(v) for v in q.get('transl_table',
                                                       ['11'])), '11')
                    # codon_start（GenBank CDS qualifier，取值 1/2/3）：CDS
                    # 起点不在首个完整密码子时，须从 codon_start-1 处切片再
                    # 翻译，否则整条蛋白移码（2026-10-03 审计 P3）；
                    # feat.extract() 返回的是含偏移的全长核酸。
                    try:
                        cs = int(next((str(v) for v in
                                       q.get('codon_start', ['1'])), '1'))
                    except ValueError:
                        cs = 1
                    cs = max(1, cs)
                    aa = str(BioSeq(nt[cs - 1:]).translate(table=table,
                                                           cds=False)).rstrip('*')
                cds_by_gene.setdefault(gene, []).append((sid, nt))
                n_cds += 1
                if len(aa) >= 1:
                    pep_by_gene.setdefault(gene, []).append((sid, aa))
                    n_pep += 1
                rows.append({'rid': rid, 'acc': rec.id, 'virus': virus,
                             'family': family, 'genus': genus, 'gene': gene,
                             'product': product, 'length': len(nt),
                             'aa_len': len(aa)})
        if prog:
            prog('export', 0.05 + 0.75 * i / max(len(gbs), 1),
                 f'解析 {i}/{len(gbs)}')

    with safe_open(os.path.join(out, 'CDS.fa'), 'wt') as fc, \
            safe_open(os.path.join(out, 'PEP.fa'), 'wt') as fp:
        for gene in sorted(cds_by_gene):
            with safe_open(os.path.join(cds_dir, f'{gene}.fa'), 'wt') as f:
                for sid, nt in cds_by_gene[gene]:
                    write_fasta_record(fc, sid, nt)
                    write_fasta_record(f, sid, nt)
        for gene in sorted(pep_by_gene):
            with safe_open(os.path.join(pep_dir, f'{gene}.fa'), 'wt') as f:
                for sid, aa in pep_by_gene[gene]:
                    write_fasta_record(fp, sid, aa)
                    write_fasta_record(f, sid, aa)

    with safe_open(os.path.join(out, 'selected.tsv'), 'wt') as f:
        f.write('rid\taccession\tvirus\tfamily\tgenus\tgene\tproduct\t'
                'cds_length\taa_length\n')
        for r in rows:
            f.write('\t'.join(str(r[k]) for k in
                              ('rid', 'acc', 'virus', 'family', 'genus',
                               'gene', 'product', 'length', 'aa_len')) + '\n')

    miss = sorted(set(want) - hit)
    if miss:
        log(f'注意：{len(miss)} 条所选 CDS 未在集合中找到（记录可能已变动）: '
            f'{", ".join(miss[:5])}{" ..." if len(miss) > 5 else ""}', 'WARN')
    log(f'导出完成：{n_cds} 条 CDS / {n_pep} 条蛋白，'
        f'{len(cds_by_gene)} 个基因分组 → {out}')
    if prog:
        prog('done', 1.0, f'完成: {n_cds} CDS / {n_pep} PEP')
    return {'dir': out, 'n_cds': n_cds, 'n_pep': n_pep,
            'genes': sorted(cds_by_gene), 'missing': miss,
            'cds_all': os.path.join(out, 'CDS.fa'),
            'pep_all': os.path.join(out, 'PEP.fa'),
            'tsv': os.path.join(out, 'selected.tsv')}
