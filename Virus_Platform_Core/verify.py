# -*- coding: utf-8 -*-
"""候选序列验证（对齐 MMPV-RNA 02b_Filter / 09b_Analysis_Verify）。

设计要点
--------
1. **过滤器，不是裁决器**——不确定时倾向保留；**no-hit 必须保留**，
   那是新病毒 / 新类病毒候选池，绝不当作"排除"。
2. **不预测 ORF**——DIAMOND blastx 与 mmseqs 均直接吃核酸 contig，
   由工具自带翻译（与 02b 同构，全程不落 ORF 中间文件）。
3. **宿主分类筛选 → 长度分流 → 双路过滤（并集 / 交集）**。
   B 路证据源用 **BFVD**（Steinegger lab 病毒蛋白域库，NCBI 官方提交
   工具链 suvtk 的注释证据源，库随包在 databases/misc_db/suvtk/）：
   原本用 NCBI Cdd，因本地库 prefilter 索引 ~5GB 反复 OOM，
   2026-09-29 换为 BFVD——病毒专属库，命中即病毒证据，无需白名单。

对外接口
--------
    split_by_length(fasta)            按长度分流
    run_blastx(...)                   A 路：DIAMOND blastx vs viral_prot
    run_bfvd(...)                     B 路：mmseqs translated search vs BFVD
    run_viroid_blastn(...)            类病毒支：blastn vs viroids_v4
    verify(ctx_like)                  主编排，产出 filtered.fasta / calls.tsv / summary.json
"""
import os
import csv
import json
import glob

from .config import get_config, db_path
from .utils import run_cmd, safe_open

# ---------------------------------------------------------------- 常量
LEN_VIRUS_MIN = 1000      # ≥ 此长度走病毒分支
LEN_VIROID_MIN = 200      # [200, 1000) 走类病毒分支；< 200 不处理
                          # 2026-10-03 起下限与组装输出默认 min_contig_len=200 对齐：
                          # 类病毒全长 246–401 nt，>200bp 即可进入类病毒支体检

BLASTX_EVALUE = 1e-3      # 对齐 02b filter_virus.py
BFVD_EVALUE_MAX = 1e-3

# 类病毒判定阈值（2026-09-08 用户确认）
VIROID_MIN_IDENT = 80.0
VIROID_MIN_COV = 60.0
VIROID_MIN_COV_SHORT = 40.0     # 250–300 bp 放宽
VIROID_SHORT_MAX = 300

# 病毒支「已知 vs 新病毒」双阈值（对齐参考项目 virome_discovery_pipeline
# 的 integrate_rescue_evidence.py / gen_final_judgement.py 口径）：
#   known：identity≥85% 且 subject 覆盖≥85%（物种级，已知种）
#   novel：identity≥40% 且 subject 覆盖≥50%（跨种/属级，新种候选）
VIRUS_KNOWN_IDENT = 85.0
VIRUS_KNOWN_COV = 85.0
VIRUS_NOVEL_IDENT = 40.0
VIRUS_NOVEL_COV = 50.0

VIRAL_PROT_DIR = db_path('annot', 'prot')
VIROID_DIR = db_path('misc', 'viroids')


def refresh_paths():
    """重算模块级路径常量（切换「设置 → 数据库目录」后由 config 回调）。"""
    global VIRAL_PROT_DIR, VIROID_DIR
    VIRAL_PROT_DIR = db_path('annot', 'prot')
    VIROID_DIR = db_path('misc', 'viroids')
    return VIRAL_PROT_DIR

# 分类谱系列（与 virus_classification.tsv 对齐），随每条 call 透出
LINEAGE = ['taxid', 'taxon', 'realm', 'kingdom', 'phylum', 'class',
           'order', 'family', 'genus', 'species', 'near_complete', 'score']

# calls.tsv 全字段（病毒支 / 类病毒支取并集；缺失列写空）
CALL_FIELDS = (['contig', 'branch', 'call', 'len'] + LINEAGE
               + ['host_source']
               + ['blastx_target', 'blastx_pident', 'blastx_alen',
                  'blastx_scov', 'blastx_evalue', 'blastx_bitscore',
                  'blastx_stitle', 'blastx_nhits']
               + ['bfvd_domains', 'bfvd_descriptions',
                  'bfvd_best_evalue', 'bfvd_nhits']
               + ['viroid_target', 'viroid_pident', 'viroid_cov',
                  'viroid_evalue', 'viroid_nhits'])


# ---------------------------------------------------------------- 元数据
def load_ictv_host():
    """ICTV 科 → 宿主来源（plants / fungi / bacteria / invertebrates ...）。"""
    path = os.path.join(VIRAL_PROT_DIR, 'ictv_family_host.tsv')
    if not os.path.isfile(path):
        return {}
    out = {}
    with safe_open(path) as f:
        for r in csv.DictReader(f, delimiter='\t'):
            if r.get('family'):
                out[r['family'].strip()] = (r.get('host.source') or '').strip()
    return out


def load_viroid_taxonomy():
    """类病毒 accession → 分类（family / genus / species）。"""
    path = os.path.join(VIROID_DIR, 'viroids.taxonomy_info.tsv')
    if not os.path.isfile(path):
        return {}
    out = {}
    with safe_open(path) as f:
        for r in csv.DictReader(f, delimiter='\t'):
            acc = (r.get('accession') or '').strip()
            if acc:
                out[acc] = {'family': (r.get('family') or '-').strip(),
                            'genus': (r.get('genus') or '-').strip(),
                            'species': (r.get('species') or '-').strip(),
                            'description': (r.get('description') or '').strip()}
    return out


# ---------------------------------------------------------------- 分流
def split_by_length(fasta, taxid=None, logger=None):
    """按长度分流。返回 (virus, viroid, skipped)，元素为 {id: seq}。

    taxid: {contig: taxid}（kunpeng 分类结果，可空）。凡 taxid>0 的 contig
    即已被 kunpeng 判为病毒，**无论长度都走病毒支**——否则 500–1000bp 的
    病毒 contig 会被长度规则路由到类病毒支，丢掉 kunpeng 的分类身份
    （病毒支才有谱系/宿主归属与 blastx+BFVD 双路证据），甚至被当成类病毒
    误判。只有 kunpeng 未分类（taxid 缺省/为 0）的短序列才走类病毒支。

    重复 header 会追加 "_2"、"_3" 后缀去重，避免静默覆盖。
    """
    taxid = taxid or {}
    virus, viroid, skipped = {}, {}, []
    rescued = 0
    cur, buf, seen = None, [], {}
    with safe_open(fasta) as f:
        for line in f:
            if line.startswith('>'):
                if cur is not None:
                    rescued += _route(cur, ''.join(buf), virus, viroid,
                                      skipped, taxid)
                raw = line[1:].split()[0]
                n = seen.get(raw, 0) + 1
                seen[raw] = n
                cur = raw if n == 1 else f'{raw}_{n}'
                buf = []
            else:
                buf.append(line.strip())
    if cur is not None:
        rescued += _route(cur, ''.join(buf), virus, viroid, skipped, taxid)
    if rescued and logger:
        logger.log(f'长度分流：{rescued} 条短于 {LEN_VIRUS_MIN}bp 但已被 '
                   f'kunpeng 判定为病毒（taxid>0）的 contig 改走病毒支')
    return virus, viroid, skipped


def _contig_taxid(cid, taxid):
    """取 contig 的 kunpeng taxid；缺省/非法一律 0（未分类）。"""
    try:
        return int(taxid.get(cid) or 0)
    except (TypeError, ValueError):
        return 0


def _route(cid, seq, virus, viroid, skipped, taxid=None):
    n = len(seq)
    if _contig_taxid(cid, taxid or {}) > 0:
        # kunpeng 已判病毒：长度门槛不适用，短 contig 也走病毒支（返回 1
        # 表示这是一条被"救回"的短 contig，供调用方汇总日志）
        virus[cid] = seq
        return 1 if n < LEN_VIRUS_MIN else 0
    if n >= LEN_VIRUS_MIN:
        virus[cid] = seq
    elif n >= LEN_VIROID_MIN:
        viroid[cid] = seq
    else:
        skipped.append(cid)
    return 0


def _write_fasta(path, records):
    # newline='\n'：显式 LF，避免 Windows 下裸 open 的文本模式把 '\n'
    # 翻译成 CRLF（下游比对/建库工具与跨平台校验都按 LF 口径）
    with open(path, 'w', encoding='utf-8', newline='\n') as f:
        for cid, seq in records.items():
            f.write('>%s\n' % cid)
            for i in range(0, len(seq), 60):
                f.write(seq[i:i + 60] + '\n')
    return path


# ---------------------------------------------------------------- 引擎
def verify_engine():
    """探测验证所需引擎。返回 {'diamond': path|None, 'mmseqs': path|None}。

    两路都缺时阶段不可用（blastx 需 DIAMOND，BFVD 需 mmseqs2）；
    单缺一路仍可降级跑（methods 里缺的那路会被 verify() 自行跳过）。
    """
    cfg = get_config()
    out = {'diamond': None, 'mmseqs': None}
    for name in out:
        try:
            p = cfg.tool(name)
        except Exception:
            p = None
        if p and os.path.isfile(p):
            out[name] = p
    return out


def run_blastx(fa, out_tsv, threads=None, logger=None):
    """DIAMOND blastx：核酸 contig vs 病毒蛋白库（工具自带 6-frame 翻译）。

    返回 {query: [hit, ...]}。hit 仅含比对字段与 stitle，**不做宿主/科推断**
    ——宿主归属统一由 run_dir/virus_classification.tsv 提供（见 _load_cls_family）。
    """
    cfg = get_config()
    diamond = cfg.tool('diamond')
    dmnd = os.path.join(VIRAL_PROT_DIR, 'viral_prot.dmnd')
    if not os.path.isfile(dmnd):
        if logger:
            logger.log('未找到 viral_prot.dmnd，跳过 blastx 分支', 'WARN')
        return None
    # 不加 --long-reads：DIAMOND 帮助原文 "short for --range-culling --top 10
    # -F 15"，即强制只报 top score 10% 范围内的命中（--top 10；-F 15 为
    # frameshift penalty）——正好把"新病毒候选"赖以判定的远源/弱相似命中
    # 压掉（本模块的过滤器语义是不确定时倾向保留）。去掉后走默认口径。
    run_cmd([diamond, 'blastx', '-q', fa, '--db', dmnd,
             '-o', out_tsv, '-e', str(BLASTX_EVALUE),
             '--threads', str(threads or cfg.threads),
             '--max-target-seqs', '10',
             # 内存限幅：--block-size 默认 2.0（20 亿字母）会一次性吃数 GB；
             # 病毒蛋白库 229MB，0.4 分块 + 8 路索引足以在 ~1-2GB 内完成，
             # 系统 内存/页面文件 紧张时避免 0xC0000005 访问冲突崩溃。
             '--block-size', '0.4', '--index-chunks', '8',
             '--outfmt', '6', 'qseqid', 'sseqid', 'pident', 'length',
             'evalue', 'bitscore', 'stitle', 'qlen', 'slen'], logger=logger)
    by_q = {}
    if not os.path.isfile(out_tsv):
        return by_q
    with open(out_tsv, encoding='utf-8', errors='replace') as f:
        for line in f:
            p = line.rstrip('\n').split('\t')
            if len(p) < 9:
                continue
            q = p[0].split()[0]
            try:
                pid, alen, ev = float(p[2]), int(p[3]), float(p[4])
                slen = int(p[8])
            except ValueError:
                continue
            stitle = p[6] if len(p) > 6 else ''
            # subject 覆盖度（对齐参考项目 aa_qcov 口径：blastx 命中蛋白
            # 的比对覆盖比例，用于区分「已知病毒 vs 新病毒候选」双阈值）
            scov = round(alen / max(slen, 1) * 100, 1)
            by_q.setdefault(q, []).append({
                'target': p[1].split()[0], 'pident': pid, 'alen': alen,
                'evalue': ev, 'bitscore': float(p[5]) if p[5] else 0.0,
                'stitle': stitle, 'scov': scov,
            })
    # 按 (-bitscore, evalue) 排序：DIAMOND 输出顺序不保证，直接取 bh[0]
    # 当"最佳命中"会让 virus_call/calls.tsv 的最优命中随输出顺序漂移
    for q in by_q:
        by_q[q].sort(key=lambda h: (-h['bitscore'], h['evalue']))
    if logger:
        logger.log('  [blastx] %d 条获得命中' % len(by_q))
    return by_q


# ---------------------------------------------------------------- B 路
def run_bfvd(fa, out_tsv, threads=None, logger=None):
    """mmseqs easy-search：核酸 contig vs BFVD（自动 translated search）。

    BFVD 病毒专属库（NCBI 提交链 suvtk 的注释证据源，随包在
    databases/misc_db/suvtk/bfvd），命中即病毒域证据——无白名单/tier。
    返回 {query: [hit, ...]}；库缺失返回 None（降级单路 blastx）。
    """
    cfg = get_config()
    try:
        cfg.tool('mmseqs')
    except FileNotFoundError:
        if logger:
            logger.log('mmseqs2 未配置，跳过 BFVD 分支', 'WARN')
        return None
    from .bfvd_search import bfvd_db, search_contigs
    if not bfvd_db():
        if logger:
            logger.log('未找到 BFVD 库（databases/misc_db/suvtk/bfvd），'
                       '跳过 BFVD 分支', 'WARN')
        return None
    hits = search_contigs(fa, out_tsv, threads=threads, logger=logger)
    return hits


def _load_cls(run_dir, assembly_dir=None):
    """从病毒分类表读 contig → 完整分类谱系行。

    宿主归属唯一来源（不用 organism_tax.tsv）。返回 {contig: 原始行 dict}，
    含 taxid / taxon / realm..species / near_complete / score 等全部列。

    分类表位置由调用者明确指定（assembly_dir），因为两种入口语义不同：
      - 管道：③ 组装的产物在 <sample_dir>/03_assembly/
      - 工具：app.py 会把分类表复制到运行目录根
    assembly_dir 缺省时依次探测这两个位置，兼容旧调用。
    """
    cls = {}
    cands = []
    if assembly_dir:
        cands.append(os.path.join(assembly_dir, 'virus_classification.tsv'))
    else:
        cands.append(os.path.join(run_dir, 'virus_classification.tsv'))
        cands.append(os.path.join(run_dir, '03_assembly',
                                  'virus_classification.tsv'))
    cls_tsv = None
    for cand in cands:
        if os.path.isfile(cand):
            cls_tsv = cand
            break
    if not cls_tsv:
        return cls
    with safe_open(cls_tsv) as f:
        for r in csv.DictReader(f, delimiter='\t'):
            cid = (r.get('contig') or r.get('Contig')
                   or r.get('sequence') or '').strip()
            if cid:
                cls[cid] = r
    return cls


# ---------------------------------------------------------------- 类病毒支
def run_viroid_blastn(fa, out_tsv, threads=None, logger=None):
    """blastn：200–1000bp contig vs viroids_v4。返回 {query: [hit, ...]}。"""
    cfg = get_config()
    blastn = cfg.tool('blastn')
    db = os.path.join(VIROID_DIR, 'viroids_v4')
    # 版本无关判定：v4 库出 .nsq/.nhr…，v5 库出 .ndb/.njs…——只认 .nsq 会让
    # v5 库被静默跳过。按 `.n??` 后缀 glob（与 assembly._db_present 同口径）。
    if not glob.glob(db + '.n??'):
        if logger:
            logger.log('未找到 viroids_v4 库，跳过类病毒分支', 'WARN')
        return None
    run_cmd([blastn, '-query', fa, '-db', db,
             '-outfmt',
             '6 qseqid sseqid pident length qlen slen evalue bitscore stitle',
             '-evalue', '1e-5', '-max_target_seqs', '5',
             '-num_threads', str(threads or cfg.threads),
             '-out', out_tsv], logger=logger)
    tax = load_viroid_taxonomy()
    by_q = {}
    if not os.path.isfile(out_tsv):
        return by_q
    with open(out_tsv, encoding='utf-8', errors='replace') as f:
        for line in f:
            p = line.rstrip('\n').split('\t')
            if len(p) < 8:
                continue
            q = p[0].split()[0]
            try:
                pid, alen, qlen, slen = (float(p[2]), int(p[3]),
                                         int(p[4]), int(p[5]))
            except ValueError:
                continue
            # 覆盖度按较短序列计（类病毒序列短，比 identity 更关键）
            cov = round(alen / max(min(qlen, slen), 1) * 100, 1)
            acc = p[1].split()[0].replace('ref|', '').replace('gb|', '')
            acc = acc.rstrip('|')
            t = tax.get(acc, {})
            by_q.setdefault(q, []).append({
                'target': acc, 'pident': pid, 'alen': alen, 'cov': cov,
                'evalue': float(p[6]) if p[6] else 1.0,
                'family': t.get('family', '-'), 'genus': t.get('genus', '-'),
                'species': t.get('species', '-'),
                'description': t.get('description', p[8] if len(p) > 8 else ''),
            })
    for q in by_q:
        by_q[q].sort(key=lambda h: (-h['pident'], h['evalue']))
    if logger:
        logger.log('  [viroid-blastn] %d 条获得命中' % len(by_q))
    return by_q


def viroid_call(hits, qlen):
    """类病毒判定：identity ≥80% 且覆盖 ≥60%（250–300bp 放宽到 ≥40%）。"""
    if not hits:
        return 'unclassified'          # no-hit 保留：可能是新类病毒
    min_cov = VIROID_MIN_COV_SHORT if qlen < VIROID_SHORT_MAX else VIROID_MIN_COV
    for h in hits:
        if h['pident'] >= VIROID_MIN_IDENT and h['cov'] >= min_cov:
            return 'viroid'
    return 'unclassified'


# ---------------------------------------------------------------- 病毒支判定
def virus_call(blastx_hits, bfvd_hits):
    """病毒支细分判定（对齐参考项目「域证据判是不是病毒 + blast 判已知/新」）。

    返回四档：
      known         已知病毒（blastx identity≥85% 且 subject 覆盖≥85%）
      novel         新病毒候选（blastx identity≥40% 且 subject 覆盖≥50%）
      domain_only   仅 BFVD 病毒域证据（无达标 blastx，或仅有弱/远源命中）
      unclassified  no-hit（新病毒候选池，保留）

    viral_prot 是纯病毒蛋白库，不存在「命中非病毒需排除」；BFVD 亦为
    病毒专属库，命中直接视为病毒域证据。过滤器语义不变：no-hit 保留。
    """
    bh = blastx_hits or []
    ch = bfvd_hits or []
    has_viral_dom = bool(ch)
    if bh:
        best = bh[0]
        pid = best.get('pident', 0.0) or 0.0
        scov = best.get('scov', 0.0) or 0.0
        if pid >= VIRUS_KNOWN_IDENT and scov >= VIRUS_KNOWN_COV:
            return 'known'
        if pid >= VIRUS_NOVEL_IDENT and scov >= VIRUS_NOVEL_COV:
            return 'novel'
        # 弱 / 远源蛋白命中：有病毒域则归「仅域证据」，否则仍算远源新病毒候选
        return 'domain_only' if has_viral_dom else 'novel'
    if has_viral_dom:
        return 'domain_only'
    return 'unclassified'


def combine_pass(a_pass, b_pass, mode):
    """组合 A/B 两路结果。mode: union | intersection（别名 intersect）| a | b。"""
    if mode == 'intersect':
        # 历史前端值（旧版管道页发 'intersect'，界面语义=两路均通过）
        mode = 'intersection'
    if mode == 'union':
        return a_pass | b_pass
    if mode == 'intersection':
        return a_pass & b_pass
    if mode == 'a':
        return set(a_pass)
    if mode == 'b':
        return set(b_pass)
    return a_pass | b_pass


# ---------------------------------------------------------------- 主编排
def verify(run_dir, fasta, host='all', methods=('blastx', 'bfvd'),
           combine='union', threads=None, logger=None, progress=None,
           out_subdir='verify', assembly_dir=None):
    """候选序列验证。产物落在 run_dir/<out_subdir>/。

    out_subdir：工具调用用默认 'verify'；管道 03b_verify 阶段传 '03b_verify'。
    assembly_dir：分类谱系表（virus_classification.tsv）所在目录。管道传
    <sample_dir>/03_assembly；工具不传（分类表已被复制到 run_dir 根）。
    methods：'blastx' / 'bfvd' 组合；combine：union / intersection。
    返回摘要 dict（供任务引擎展示）。
    """
    def _prog(stage, frac, msg):
        if progress:
            progress(stage, frac, msg)
        if logger:
            logger.log(msg)

    out_dir = os.path.join(run_dir, out_subdir)
    os.makedirs(out_dir, exist_ok=True)

    # 缺工具的那一路要在这里剔除：run_blastx 直接 cfg.tool('diamond')
    # （缺则 FileNotFoundError），若不预筛会让整个 03b 阶段失败，
    # 而 verify_engine 的契约是"单缺一路仍可降级跑"。
    _eng = verify_engine()
    _methods = [m for m in (methods or ('blastx', 'bfvd'))
                if m in ('blastx', 'bfvd')]
    dropped = [m for m in (methods or ()) if m not in ('blastx', 'bfvd')]
    if dropped and logger:
        logger.log('未知证据方法 %s 忽略（仅支持 blastx/bfvd）'
                   % ','.join(dropped), 'WARN')
    if 'blastx' in _methods and not _eng.get('diamond'):
        _methods = [m for m in _methods if m != 'blastx']
        if logger:
            logger.log('未找到 DIAMOND，跳过 blastx 分支（降级为 BFVD）', 'WARN')
    if 'bfvd' in _methods and not _eng.get('mmseqs'):
        _methods = [m for m in _methods if m != 'bfvd']
        if logger:
            logger.log('未找到 mmseqs2，跳过 BFVD 分支', 'WARN')
    if not _methods:
        raise RuntimeError('候选序列验证无可用引擎（DIAMOND/mmseqs2 全缺）')
    methods = tuple(_methods)

    # 分类表先于分流加载：分流需要 kunpeng 的 taxid 判定（taxid>0 的短 contig
    # 必须走病毒支，见 split_by_length），不能等宿主筛选时才读。
    cls = _load_cls(run_dir, assembly_dir=assembly_dir)
    taxids = {}
    for cid, row in cls.items():
        try:
            t = int(row.get('taxid') or 0)
        except (TypeError, ValueError):
            t = 0
        if t > 0:
            taxids[cid] = t

    _prog('split', 0.05, '按长度分流（kunpeng taxid>0 的 contig 恒走病毒支）')
    virus, viroid, skipped = split_by_length(fasta, taxid=taxids, logger=logger)
    _prog('split', 0.10,
          '分流完成：病毒 %d / 类病毒 %d / 过短跳过 %d'
          % (len(virus), len(viroid), len(skipped)))

    calls = []          # 每行：contig, branch, call, detail...
    summary = {'virus': 0, 'viroid': 0, 'skipped': len(skipped),
               'calls': {}, 'combine': combine, 'methods': list(methods),
               'host_filter': host}

    # ---- 宿主筛选（作用于病毒支）----
    ictv_host = load_ictv_host()
    if host and host != 'all':
        virus, dropped = _filter_by_host(virus, host, cls, ictv_host)
        summary['host_dropped'] = len(dropped)
        _prog('host', 0.15, '宿主筛选 %s：保留 %d，剔除 %d'
              % (host, len(virus), len(dropped)))

    # ---- 病毒支 ----
    if virus:
        vfa = _write_fasta(os.path.join(out_dir, 'virus_input.fasta'), virus)
        a_hits = b_hits = None
        if 'blastx' in methods:
            _prog('blastx', 0.25, 'DIAMOND blastx vs 病毒蛋白库')
            a_hits = run_blastx(vfa, os.path.join(out_dir, 'blastx.tsv'),
                                threads=threads, logger=logger) or {}
        if 'bfvd' in methods:
            _prog('bfvd', 0.5, 'mmseqs2 病毒域搜索 vs BFVD')
            b_hits = run_bfvd(vfa, os.path.join(out_dir, 'bfvd_hits.tsv'),
                              threads=threads, logger=logger) or {}
        a_pass = set(a_hits) if a_hits is not None else set()
        b_pass = set(b_hits) if b_hits is not None else set()

        # 只用单路时，另一路不参与组合
        if a_hits is None:
            passed = b_pass
        elif b_hits is None:
            passed = a_pass
        else:
            passed = combine_pass(a_pass, b_pass, combine)

        # 过滤器语义：未通过双路者也保留为 unclassified（新病毒候选）
        for cid, seq in virus.items():
            bh = (a_hits or {}).get(cid)
            ch = (b_hits or {}).get(cid)
            call = virus_call(bh, ch) if cid in passed else 'unclassified'
            rec = {'contig': cid, 'branch': 'virus', 'call': call,
                   'len': len(seq)}
            # 完整分类谱系（唯一来源 virus_classification.tsv）
            cr = cls.get(cid, {})
            for k in LINEAGE:
                rec[k] = (cr.get(k) or '').strip()
            rec['host_source'] = ictv_host.get(rec['family'], '')
            # blastx 命中详情（最优命中 + 命中数）
            if bh:
                b0 = bh[0]
                rec.update({
                    'blastx_target': b0.get('target', ''),
                    'blastx_pident': b0.get('pident', ''),
                    'blastx_alen': b0.get('alen', ''),
                    'blastx_scov': b0.get('scov', ''),
                    'blastx_evalue': b0.get('evalue', ''),
                    'blastx_bitscore': b0.get('bitscore', ''),
                    'blastx_stitle': b0.get('stitle', ''),
                    'blastx_nhits': len(bh),
                })
            # BFVD 命中详情（病毒专属域库，命中即病毒证据）
            if ch:
                descs = [h.get('product', '') for h in ch if h.get('product')]
                rec.update({
                    'bfvd_domains': ';'.join(
                        h.get('acc') or h.get('target', '') for h in ch),
                    'bfvd_descriptions': '; '.join(descs) if descs else '',
                    'bfvd_best_evalue': min((h.get('evalue') for h in ch),
                                            default=''),
                    'bfvd_nhits': len(ch),
                })
            calls.append(rec)
        _prog('virus', 0.75, '病毒支：%d 条，其中 %d 条有证据'
              % (len(virus), len(passed)))

    # ---- 类病毒支 ----
    if viroid:
        _prog('viroid', 0.85, 'blastn vs 类病毒库（viroids_v4）')
        dfa = _write_fasta(os.path.join(out_dir, 'viroid_input.fasta'), viroid)
        hits = run_viroid_blastn(dfa, os.path.join(out_dir, 'viroid_blastn.tsv'),
                                 threads=threads, logger=logger) or {}
        for cid, seq in viroid.items():
            hh = hits.get(cid)
            call = viroid_call(hh, len(seq))
            best = hh[0] if hh else {}
            rec = {'contig': cid, 'branch': 'viroid', 'call': call,
                   'len': len(seq)}
            # 类病毒谱系来自 blastn 命中物种（taxonomy_info.tsv）
            for k in LINEAGE:
                rec[k] = ''
            rec['family'] = best.get('family', '')
            rec['genus'] = best.get('genus', '')
            rec['species'] = best.get('species', '')
            rec['taxon'] = best.get('species', '') or best.get('description', '')
            rec['host_source'] = ''
            # blastn 命中详情（最优命中 + 命中数）
            if best:
                rec.update({
                    'viroid_target': best.get('target', ''),
                    'viroid_pident': best.get('pident', ''),
                    'viroid_cov': best.get('cov', ''),
                    'viroid_evalue': best.get('evalue', ''),
                    'viroid_nhits': len(hh),
                })
            calls.append(rec)
        _prog('viroid', 0.92, '类病毒支：%d 条，判定为类病毒 %d 条'
              % (len(viroid),
                 sum(1 for c in calls
                     if c['branch'] == 'viroid' and c['call'] == 'viroid')))

    # ---- 产物 ----
    with open(os.path.join(out_dir, 'calls.tsv'), 'w', encoding='utf-8',
              newline='') as f:
        w = csv.DictWriter(f, fieldnames=CALL_FIELDS, delimiter='\t',
                           restval='')
        w.writeheader()
        for r in calls:
            w.writerow(r)

    for r in calls:
        summary['calls'][r['call']] = summary['calls'].get(r['call'], 0) + 1
    summary['virus'] = sum(1 for r in calls if r['branch'] == 'virus')
    summary['viroid'] = sum(1 for r in calls if r['branch'] == 'viroid')
    summary['n_calls'] = len(calls)
    summary['out_dir'] = out_dir
    with open(os.path.join(out_dir, 'summary.json'), 'w', encoding='utf-8') as f:
        json.dump(summary, f, ensure_ascii=False, indent=1)

    _prog('done', 1.0, '验证完成：%s' % json.dumps(summary['calls'],
                                                   ensure_ascii=False))
    return summary


def _filter_by_host(records, host, cls, ictv_host):
    """按 ICTV 科 → 宿主来源筛选。分类唯一来源 = virus_classification.tsv。

    host.source 是复合值（如 "algae+fungi+plants"），故用子串包含匹配；
    值为空（未标注）或匹配到所选宿主即保留，避免因缺信息误删。
    """
    if host == 'protozoa':
        # 历史前端值；ICTV host.source 的用词是 protists（复合值子串匹配）
        host = 'protists'
    keep, drop = {}, []
    for cid, seq in records.items():
        cr = cls.get(cid, {})
        fam = (cr.get('family') or '').strip()
        src = ictv_host.get(fam, '')
        if not src or host in src:
            keep[cid] = seq
        else:
            drop.append(cid)
    return keep, drop
