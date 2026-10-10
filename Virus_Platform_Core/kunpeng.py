# -*- coding: utf-8 -*-
"""
kunpeng 封装：数据库构建（add-library → build-db）与序列分类（classify）。
- build_host_db():  宿主基因组(单一 taxid)建库
- build_virus_db(): 病毒参考(按 info 表 accession→taxid)建库
- classify():       reads/contigs 分类，返回 Kraken 兼容输出与 kreport2
- parse_classify_output(): 解析 per-read 输出
"""
import os
import csv
import glob
import shutil
import time

from .config import get_config, DIRS, host_db_name, host_db_taxids, host_db_info
from .utils import (check_path, safe_open, run_cmd,
                    inject_taxid_map,
                    inject_taxid_chunked, est_decompressed, log_res_plan,
                    dir_size)

MANIFEST_NAME = 'host_db.json'

# Windows 异常退出码 → 常见原因（含本机 MSYS2 缺 DLL 的实测值）。
# 裸退出码对用户毫无意义，kv_consensus 早有同样的映射表，这里对齐口径。
_WIN_EXIT_HINTS = {
    3221225781: '0xC0000135 依赖 DLL 缺失（MSYS2 runtime/libgcc 未随工具拷贝）',
    3221225785: '0xC0000139 入口点找不到（DLL 版本不匹配）',
    3221225477: '0xC0000005 访问违例（被杀软拦截或二进制损坏）',
    3221226505: '0xC0000409 栈缓冲溢出/内存分配失败',
}


def write_host_db_manifest(db_dir, taxid, genome_fasta, n_seq, n_frag,
                           logger=None):
    """把"这个库是给谁建的"写成 <db_dir>/host_db.json。

    为什么需要：库本体（hash_*.k2d / taxo.k2d / seqid2taxid.map）不记录物种，
    历史上就出过事故——日志显示 2026-09-03 那次建库把 689 条序列（= 枸杞
    Lycium barbarum 基因组的序列数）注成了 taxid=4081（番茄 Solanum
    lycopersicum），此后该库的 map/taxo 长期带着两个 taxid，而用错宿主去做
    ①宿主去除是**静默**的。清单让库身份从"反推"变成"声明"：①宿主去除、自检、
    GUI 列表都能直接核对。

    写清单失败不阻断建库（库本身已建成），但一定留 WARN —— 沉默正是问题根源。
    """
    import json
    from .taxonomy import taxid_name
    db_dir = os.path.normpath(os.path.abspath(db_dir))
    genome = os.path.normpath(os.path.abspath(str(genome_fasta)))
    try:
        st = os.stat(genome)
        g_size, g_mtime = st.st_size, int(st.st_mtime)
    except OSError:
        g_size, g_mtime = None, None
    try:
        hash_bytes = sum(os.path.getsize(p)
                         for p in glob.glob(os.path.join(db_dir, 'hash_*.k2d')))
    except OSError:
        hash_bytes = None
    data = {
        'schema': 1,
        'host_db_name': os.path.basename(db_dir),
        'taxid': int(taxid),
        'species': taxid_name(taxid),
        'source_genome': genome,
        'source_genome_size': g_size,
        'source_genome_mtime': g_mtime,
        'n_seq': n_seq,
        'n_frag': n_frag,
        'taxids_in_map': sorted(host_db_taxids(db_dir)),
        'hash_bytes': hash_bytes,
        'built_at': time.strftime('%Y-%m-%d %H:%M:%S'),
        'kunpeng': get_config().tool('kunpeng'),
    }
    p = os.path.join(db_dir, MANIFEST_NAME)
    try:
        with open(p, 'w', encoding='utf-8') as f:
            json.dump(data, f, ensure_ascii=False, indent=2)
    except OSError as e:
        if logger:
            logger.log(f"宿主库清单写入失败（不影响库本身）: {p}: {e}", "WARN")
        return data
    if logger:
        bad = [t for t in data['taxids_in_map'] if t != int(taxid)]
        logger.log(f"宿主库清单: {p}")
        logger.log(f"  物种 = {data['species'] or '(未查到学名)'} (taxid={taxid})，"
                   f"源基因组 = {os.path.basename(genome)}，"
                   f"序列 {n_seq} 条 → 片段 {n_frag} 个")
        if bad:
            logger.log(f"⚠ 库内 seqid2taxid.map 仍含其它 taxid {bad} —— "
                       f"说明该目录混入了历史建库（元数据污染），"
                       f"建议换新目录重建", "WARN")
    return data


def db_ready(db_dir):
    """kunpeng 库是否可分类（hash 表 + 元数据齐全）。"""
    try:
        d = check_path(db_dir, must_exist=True, in_platform=True)
    except (FileNotFoundError, ValueError):
        return False
    has_hash = bool(glob.glob(os.path.join(d, 'hash_*.k2d')))
    has_meta = all(os.path.isfile(os.path.join(d, f))
                   for f in ('opts.k2d', 'taxo.k2d', 'hash_config.k2d'))
    return has_hash and has_meta


def db_building_marker(db_dir):
    return check_path(os.path.join(db_dir, '.building'), must_exist=False, in_platform=True)


def ensure_db_dirs(db_dir):
    """创建 library/ 与 taxonomy/（复制共享 taxonomy 的 dmp 文件）。"""
    d = check_path(db_dir, must_exist=False, in_platform=True)
    lib = os.path.join(d, 'library')
    tax = os.path.join(d, 'taxonomy')
    os.makedirs(lib, exist_ok=True)
    os.makedirs(tax, exist_ok=True)

    src_tax = check_path(DIRS['taxonomy'], must_exist=True)
    copied = []
    for name in ('nodes.dmp', 'names.dmp', 'merged.dmp'):
        src = check_path(os.path.join(src_tax, name), must_exist=False)
        dst = check_path(os.path.join(tax, name), must_exist=False, in_platform=True)
        if os.path.isfile(src) and not os.path.isfile(dst):
            shutil.copyfile(src, dst)
            copied.append(name)
    return d, copied


def _clean_mid_files(db_dir):
    """建库成功后清理中间产物（library/prep/taxonomy），只留 .k2d + 映射 + 清单。

    taxonomy/ 是 ensure_db_dirs 从 databases/tax_db **复制**进来的
    nodes/names/merged.dmp（实测 571 MB/库），只有 add-library/build-db 用得到：
    实测把它挪走后 classify 仍成功，kreport 与基线逐字段一致（物种名来自已生成的
    taxo.k2d）。所以它不是"运行期数据"，而是建库耗材——不需要为它做软链接
    （链接还带绝对路径，换盘/拷贝即失效），删掉即可；重建时 ensure_db_dirs 会
    自动再复制。
    """
    for sub in ('library', 'prep', 'taxonomy'):
        p = os.path.join(db_dir, sub)
        if os.path.isdir(p):
            shutil.rmtree(p, ignore_errors=True)
            if os.path.isdir(p):
                raise RuntimeError(f'中间文件清理失败: {p}')


def build_db(db_dir, tagged_fasta, hash_capacity='256M', threads=None,
             logger=None, rebuild=False, clean_mid=False):
    """add-library + build-db。tagged_fasta 头已含 kraken:taxid|N。"""
    cfg = get_config()
    kunpeng = cfg.tool('kunpeng')
    threads = threads or cfg.threads

    if db_ready(db_dir) and not rebuild:
        if logger:
            logger.log(f"数据库已存在，跳过建库: {db_dir}")
        return db_dir

    d, copied = ensure_db_dirs(db_dir)
    if logger:
        if copied:
            logger.log(f"taxonomy 文件复制: {', '.join(copied)}")
        logger.log(f"步骤1/2 add-library: {tagged_fasta}")
    marker = db_building_marker(db_dir)
    with safe_open(marker, 'wt') as f:
        f.write(time.strftime('%Y-%m-%d %H:%M:%S'))

    # 替换语义：tagged_fasta 总是全量文件，清掉历史 library 与失败残留，
    # 否则 add-library 累积追加会让重试时数据翻倍（曾导致 3.6GB 重复 fna）
    #
    # ⚠ 这里必须连**元数据**一起清：kunpeng 的 add-library 是往
    # seqid2taxid.map 追加的，早先只清 library/hash、漏了 map 与 taxo.k2d，
    # 于是每次重建都把陈旧的 seqid→taxid 并进 map。实测宿主库 map 里因此长期
    # 同时存在 4081（番茄）与 112863（枸杞）两个 taxid，taxo.k2d 也跟着带上
    # 了 Solanum —— 库"是谁"变得无法判定。hash_config/opts 一并清，避免旧
    # 容量/参数残留影响本次 build-db。
    lib_dir = os.path.join(d, 'library')
    stale = (glob.glob(os.path.join(lib_dir, 'library_*.fna'))
             + glob.glob(os.path.join(lib_dir, 'library_*.hllp_*.json'))
             + glob.glob(os.path.join(lib_dir, 'added.md5'))
             + glob.glob(os.path.join(d, 'hash_*.k2d'))   # 已覆盖 hash_config.k2d
             + glob.glob(os.path.join(d, 'chunk_*.k2'))
             + glob.glob(os.path.join(d, 'seqid2taxid.map'))
             + glob.glob(os.path.join(d, 'taxo.k2d'))
             + glob.glob(os.path.join(d, 'opts.k2d')))
    # 去重：glob 模式之间存在重叠（hash_*.k2d 同时命中 hash_1.k2d 与
    # hash_config.k2d），重复项第二次删除会抛 FileNotFoundError（实测踩过）。
    stale = list(dict.fromkeys(stale))
    if stale:
        if logger:
            logger.log(f"清理历史 library/残留文件 {len(stale)} 个（替换式重建）")
        for p in stale:
            # 已被并发/上一轮删掉就跳过：这里的目标是"清干净"，不是"必须删到"
            try:
                os.remove(check_path(p, must_exist=True, in_platform=True))
            except FileNotFoundError:
                pass

    try:
        run_cmd([kunpeng, 'add-library', '--db', d, '-i',
                 check_path(tagged_fasta, must_exist=True)], logger=logger)
        if logger:
            logger.log(f"步骤2/2 build-db (hash-capacity={hash_capacity}, threads={threads})")
        run_cmd([kunpeng, 'build-db', '--db', d,
                 '--hash-capacity', str(hash_capacity),
                 '-p', str(threads)], logger=logger)
        if not db_ready(d):
            raise RuntimeError(f"build-db 结束但库校验失败: {d}")
        if logger:
            sizes = sum(os.path.getsize(p) for p in glob.glob(os.path.join(d, 'hash_*.k2d')))
            logger.log(f"建库完成，hash 表总大小 {sizes / 1e9:.2f} GB")
    finally:
        if os.path.isfile(marker):
            os.remove(marker)
    if clean_mid:
        if logger:
            logger.log("清理中间产物（library/prep），只保留 .k2d + taxonomy + 映射")
        _clean_mid_files(d)
    return d


def build_host_db(genome_fasta, taxid, db_dir=None, hash_capacity='256M',
                  threads=None, logger=None, rebuild=False, chunk_bp=1000000,
                  clean_mid=False):
    """宿主库：任意物种基因组 FASTA + 该物种 NCBI TaxID。

    taxid 会先在 NCBI 分类(nodes/merged.dmp)中校验：已合并的提示改用新值、
    不存在的直接报错（避免建库到最后 resolve 阶段才失败）。
    序列名自动规范化：取首空白前 token 作 ID、清理历史标签、空 ID 自动编号。
    kunpeng convert 阶段按 60 条/批读取且无字节上限，多条大染色体同批会触发
    数十 GB 内存分配（详见 inject_taxid_chunked 注释），故基因组一律分块。

    db_dir=None 时**按物种建目录**：host-db/<taxid>_<源目录名>_host_db/
    （如 host-db/112863_Lycium_barbarum_host_db）。旧行为是全部塞进单一固定
    槽位 host-db/host/classify，导致多宿主无法共存、且库不记录物种（见
    write_host_db_manifest 的背景说明）。显式传 db_dir 仍可覆盖。
    建完在库目录写 host_db.json 清单。
    """
    if db_dir is None:
        db_dir = os.path.join(DIRS['host_src'], host_db_name(taxid, genome_fasta))
    try:
        taxid = int(taxid)
    except (TypeError, ValueError):
        raise ValueError(f"TaxID 必须是数字（NCBI 物种编号），收到: {taxid!r}")
    from .taxonomy import taxid_lookup
    try:
        status, extra = taxid_lookup(taxid)
    except (FileNotFoundError, NotADirectoryError):
        raise RuntimeError("NCBI taxonomy 未初始化：请先在「数据库构建」页"
                           "下载 NCBI 分类库（宿主/病毒建库都需要）")
    if status == 'missing':
        raise ValueError(
            f"TaxID {taxid} 不在 NCBI 分类库中。请核对编号"
            f"（https://www.ncbi.nlm.nih.gov/taxonomy 搜索物种），"
            f"或到「数据库构建」页更新 taxonomy 后重试")
    if status == 'merged':
        raise ValueError(
            f"TaxID {taxid} 已被 NCBI 合并到 {extra}，请改用 {extra} 重新构建")
    work = os.path.join(db_dir, 'prep')
    os.makedirs(work, exist_ok=True)
    tagged = check_path(os.path.join(work, 'host_tagged.fa'),
                        must_exist=False, in_platform=True)
    if logger:
        logger.log(f"宿主库构建开始: taxid={taxid}（目录 {db_dir}）")
    n_seq, n_frag = inject_taxid_chunked(genome_fasta, tagged, taxid,
                                        chunk_bp=chunk_bp, logger=logger)
    build_db(db_dir, tagged, hash_capacity=hash_capacity,
             threads=threads, logger=logger, rebuild=rebuild,
             clean_mid=clean_mid)
    write_host_db_manifest(db_dir, taxid, genome_fasta, n_seq, n_frag,
                           logger=logger)
    return db_dir


def parse_virus_info(info_tsv):
    """解析病毒 info 表 -> ({accession: taxid}, 所有列名)。"""
    info = check_path(info_tsv, must_exist=True)
    acc2taxid = {}
    with safe_open(info) as f:
        reader = csv.DictReader(f, delimiter='\t')
        fields = reader.fieldnames or []
        for row in reader:
            acc = (row.get('Accession') or '').strip()
            tax = (row.get('Taxid') or '').strip()
            if acc and tax.isdigit():
                acc2taxid[acc] = int(tax)
    return acc2taxid, fields


VIRUS_DOMAIN_TAXID = 10239        # NCBI「Viruses」域根 taxid


def vet_info_table(info_tsv, logger=None):
    """建库准入校验（**只告警、不阻断**）：info / ref_info 表脏行体检。

    为什么需要：reference.ref_info.tsv 由源表整合而来，偶有脏行——
    Taxid 为空/非数字的行会被 parse_virus_info **静默丢弃**（参考序列进不了
    库，表现为分类召回偏低且看不出原因）；Category 为空、Taxid 指向非病毒域
    （串表/串行）则让分类库元数据不完整。在建库入口扫一遍并打印 WARN 计数，
    把问题暴露在构建日志里，而不是静默入库。

    扫三类：① Taxid 缺失/非法；② Category 为空（表无该列则跳过）；
    ③ Taxid 不在病毒域（10239 谱系；含本地 taxonomy 未收录者）。
    ③ 需 TaxonomyGraph 整表加载一次，taxonomy 不可用时跳过并说明。
    返回计数 dict；任何异常都在函数内兜底（本地表脏不该让建库失败）。
    """
    rows = []
    cat_col = None
    try:
        with safe_open(check_path(info_tsv, must_exist=True)) as f:
            reader = csv.DictReader(f, delimiter='\t')
            cols = reader.fieldnames or []
            cat_col = next((c for c in cols
                            if (c or '').strip().lower() == 'category'), None)
            for row in reader:
                acc = (row.get('Accession') or '').strip()
                tax = (row.get('Taxid') or '').strip()
                cat = (row.get(cat_col) or '').strip() if cat_col else ''
                rows.append((acc, tax, cat))
    except Exception as e:
        if logger:
            logger.log(f"准入校验跳过（info 表不可读: {e}）", "WARN")
        return {}

    bad_taxid = [a for a, t, _ in rows if a and not t.isdigit()]
    blank_cat = [a for a, _, c in rows if a and cat_col is not None and not c]

    domain_checked = False
    non_n = 0
    non_samples = []
    try:
        from .taxonomy import TaxonomyGraph, taxonomy_ready
        if taxonomy_ready():
            g = TaxonomyGraph()
            g.load()
            domain_checked = True
            _cache = {}
            for acc, tax, _cat in rows:
                if not tax.isdigit():
                    continue
                tid = int(tax)
                if tid not in _cache:
                    _r = g.resolve(tid)     # 已合并 taxid 归一到现行值
                    _lin = g.lineage(_r) if _r in g.parent else []
                    _cache[tid] = any(_t == VIRUS_DOMAIN_TAXID
                                      for _t, _rk, _nm in _lin)
                if not _cache[tid]:
                    non_n += 1
                    if len(non_samples) < 10:
                        non_samples.append(f'{acc}(taxid={tax})')
    except Exception as e:
        if logger:
            logger.log(f"  准入校验：病毒域校验跳过（taxonomy 不可用: {e}）", "WARN")

    if logger:
        logger.log(f"准入校验（{len(rows)} 行）: Taxid 缺失/非法 {len(bad_taxid)}，"
                   f"Category 为空 {len(blank_cat) if cat_col else '表无该列'}，"
                   f"非病毒域 {non_n if domain_checked else '未校验（无 taxonomy）'}")
        for _label, _items in (('Taxid 缺失/非法', bad_taxid),
                               ('Category 为空', blank_cat)):
            if _items:
                logger.log(f"  ⚠ {_label}: {_items[:10]}"
                           f"{' …' if len(_items) > 10 else ''}", "WARN")
        if non_n:
            logger.log(f"  ⚠ 非病毒域/无法归入病毒域: {non_samples}"
                       f"{' …' if non_n > 10 else ''}", "WARN")
        if bad_taxid or blank_cat or non_n:
            logger.log("  以上为告警不阻断：建库继续，建议核对 ref_info 源表",
                       "WARN")
    return {'rows': len(rows), 'bad_taxid': len(bad_taxid),
            'blank_cat': len(blank_cat) if cat_col else None,
            'non_virus': non_n if domain_checked else None,
            'domain_checked': domain_checked}


def build_virus_db(virus_fasta, info_tsv, db_dir=None, hash_capacity='64M',
                   threads=None, logger=None, rebuild=False, clean_mid=False):
    """病毒库：按 info 表 accession→taxid 逐条注入。"""
    db_dir = db_dir or get_config().databases['virus']
    work = os.path.join(db_dir, 'prep')
    os.makedirs(work, exist_ok=True)
    tagged = check_path(os.path.join(work, 'virus_tagged.fa'),
                        must_exist=False, in_platform=True)
    acc2taxid, _ = parse_virus_info(info_tsv)
    if logger:
        logger.log(f"info 表解析: {len(acc2taxid)} 条 accession→taxid 映射")
    # 建库准入校验：脏行（Taxid 缺失/Category 空/非病毒域）只告警不阻断，
    # 避免"参考没进库、分类召回偏低"这类问题被静默带进产物。
    vet_info_table(info_tsv, logger=logger)
    n, missing = inject_taxid_map(virus_fasta, tagged, acc2taxid,
                                  logger=logger, missing_fatal=False)
    if missing:
        raise KeyError(
            f"病毒参考中有 {len(missing)} 条序列在 info 表无 Taxid，示例: {missing[:5]}；"
            f"请核对 info 表与 FASTA 是否配套")
    if logger:
        logger.log("病毒库构建开始")
    return build_db(db_dir, tagged, hash_capacity=hash_capacity,
                    threads=threads, logger=logger, rebuild=rebuild,
                    clean_mid=clean_mid)


# ------------------------------------------------------------------
# 分类
# ------------------------------------------------------------------
def _pick_chunk_dir(out_dir, inputs, logger=None, base=None):
    """选择 chunk 临时目录。

    base 给定（用户在界面选择）: 校验其所在盘剩余空间，不足则报错
    （需求 = 解压后数据的约 2 倍 ≈ gz 字节 ×16；7GB 双端实测 105.6GB）。
    未给定时: 平台盘余量 ≥ 2×需求 用平台盘；否则自动选剩余空间最大的
    盘根下的 vp_chunk 目录。

    kunpeng classify 的 splitr 把 reads 解压转成 k2 中间格式全部写盘，
    annotate/resolve 再写一份同量级中间结果——磁盘开销 ≈ 解压量 ×2。
    """
    need = sum(os.path.getsize(check_path(p, must_exist=True))
               for p in inputs) * 16

    def _new_chunk(root):
        return os.path.join(root, 'chunk_%d' % int(time.time() * 1000))

    if base:   # 用户指定的 chunk 目录：只校验空间，不擅自换位置
        base = os.path.abspath(str(base))
        os.makedirs(base, exist_ok=True)
        free = shutil.disk_usage(base).free
        if free < need * 1.1:
            raise RuntimeError(
                f"chunk 目录所在磁盘空间不足：分类需要约 {need / 1e9:.0f} GB"
                f"（≈解压后数据的 2 倍），该盘剩余仅 {free / 1e9:.0f} GB。"
                f"请在全局参数中改选其他磁盘上的目录，或用「子采样」减小数据量")
        chunk = _new_chunk(base)
        os.makedirs(chunk, exist_ok=True)
        if logger:
            logger.log(f"chunk 临时目录（用户指定）: {chunk}"
                       f"（需求约 {need / 1e9:.0f} GB，该盘剩余 {free / 1e9:.0f} GB）")
        return chunk

    # 优先：平台盘（余量 ≥ 2 倍需求）
    try:
        free_out = shutil.disk_usage(out_dir).free
        if free_out >= need * 2:
            chunk = os.path.join(out_dir, '_chunk_%d' % int(time.time() * 1000))
            os.makedirs(chunk, exist_ok=True)
            if logger and need > 20e9:
                logger.log(f"chunk 临时目录（平台盘）: {chunk}"
                           f"（需求约 {need / 1e9:.0f} GB，剩余 {free_out / 1e9:.0f} GB）")
            return chunk
    except OSError:
        pass
    # 回退：其余盘按剩余空间从大到小（枚举全部盘符；原 'CDEFGH' 漏掉 I 之后的盘）
    import string as _string
    best, best_free = None, 0
    for drv in _string.ascii_uppercase:
        root = '%s:%s' % (drv, os.sep)
        if not os.path.isdir(root):
            continue
        try:
            free = shutil.disk_usage(root).free
        except OSError:
            continue
        if free > best_free:
            best, best_free = os.path.join(root, 'vp_chunk'), free
    if best is None or best_free < need * 1.1:
        raise RuntimeError(
            f"磁盘空间不足：kunpeng 分类需要约 {need / 1e9:.0f} GB 临时空间"
            f"（≈解压后数据的 2 倍），可用盘最大剩余仅 "
            f"{best_free / 1e9:.0f} GB。建议：① 在全局参数中指定大盘上的"
            f"chunk 目录；② 用「子采样」减小数据量；③ 清理磁盘")
    os.makedirs(best, exist_ok=True)
    chunk = _new_chunk(best)
    if logger:
        logger.log(f"平台盘空间紧张，chunk 临时目录改用: {chunk}"
                   f"（需求约 {need / 1e9:.0f} GB，该盘剩余 {best_free / 1e9:.0f} GB）")
    os.makedirs(chunk, exist_ok=True)
    return chunk


_SPLIT_THRESHOLD_GB = 2.0   # 输入超过此 GB 数时自动分片分类


class _MemFail(RuntimeError):
    """kunpeng 内存分配失败（Windows abort 3221226505）。"""


def _free_commit_gb():
    """当前可用的虚拟内存提交空间（GB）。"""
    import ctypes

    class M(ctypes.Structure):
        _fields_ = [('dwLength', ctypes.c_ulong),
                    ('dwMemoryLoad', ctypes.c_ulong),
                    ('ullTotalPhys', ctypes.c_ulonglong),
                    ('ullAvailPhys', ctypes.c_ulonglong),
                    ('ullTotalPageFile', ctypes.c_ulonglong),
                    ('ullAvailPageFile', ctypes.c_ulonglong),
                    ('ullTotalVirtual', ctypes.c_ulonglong),
                    ('ullAvailVirtual', ctypes.c_ulonglong),
                    ('ullAvailExtendedVirtual', ctypes.c_ulonglong)]
    m = M()
    m.dwLength = ctypes.sizeof(M)
    try:
        ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(m))
        return m.ullAvailPageFile / 1e9
    except Exception:
        return 64.0


def _split_fasta_gz(src, n_parts, work_dir, tag):
    """把 seqkit fq2fa -w 0 产出的 FASTA.gz（记录恒为 2 行：头+单行序列）
    按记录轮流转拆成 n_parts 份 gz（配对两侧用同一轮转序号即保持同步）。"""
    work_dir = check_path(work_dir, must_exist=False, in_platform=True)
    os.makedirs(work_dir, exist_ok=True)
    src = check_path(src, must_exist=True)
    outs = [safe_open(os.path.join(work_dir, f'{tag}_part{j:02d}.fa.gz'),
                      'wb', compress_level=1) for j in range(n_parts)]
    n = 0
    try:
        with safe_open(src, 'rb') as f:
            for line in f:
                outs[(n // 2) % n_parts].write(line)
                n += 1
    finally:
        for w in outs:
            w.close()
    return [os.path.join(work_dir, f'{tag}_part{j:02d}.fa.gz')
            for j in range(n_parts)]


def _merge_kreports(part_kr_files, out_kr, logger=None):
    """合并各分片 kreport2：按 (rank, taxid, name) 累加 clade/direct 计数，
    保持父先于子的行序（首现位置），百分比按合并后总 reads 重算。"""
    order = {}
    rows = []
    for pf in part_kr_files:
        with safe_open(check_path(pf, must_exist=True)) as f:
            for line in f:
                p = line.rstrip('\n').split('\t')
                if len(p) < 6:
                    continue
                try:
                    clade = int(float(p[1]))
                    direct = int(float(p[2]))
                except ValueError:
                    continue
                key = (p[3], p[4], p[5])
                if key not in order:
                    order[key] = len(rows)
                    rows.append([clade, direct, p[3], p[4], p[5]])
                else:
                    i = order[key]
                    rows[i][0] += clade
                    rows[i][1] += direct
    total_u = sum(r[0] for r in rows if r[2] == 'U')
    total_c = sum(r[0] for r in rows if r[2] != 'U' and r[3] == '1')
    denom = max(total_u + total_c, 1)
    with safe_open(check_path(out_kr, must_exist=False, in_platform=True),
                   'wt') as f:
        for clade, direct, rk, taxid, name in rows:
            f.write(f'{clade / denom * 100:.2f}\t{clade}\t{direct}\t'
                    f'{rk}\t{taxid}\t{name}\n')
    if logger:
        logger.log(f"kreport 合并: {len(part_kr_files)} 片 → "
                   f"{len(rows)} 个分类单元（总 reads {denom:,}）")


def classify(db_dir, inputs, out_dir, paired=False, threads=None,
             confidence=0.0, min_hit_groups=2, logger=None, keep_chunk=False,
             chunk_dir=None, allow_convert=True, progress=None):
    """kunpeng classify。inputs: FASTA/FASTQ(.gz) 路径列表。

    allow_convert: 输入为 FASTQ 且装了 seqkit 时是否自动 fq2fa 预转换
    （管道已做预处理转换时传 False，避免同一文件转换两次）。
    progress: 可选回调 progress(pct 0~1, msg)，运行中按 chunk 目录增长
    上报进度（kunpeng 本体无进度输出，看门狗以中间盘占用近似）。
    返回 {'kraken': output_*.txt 路径或None, 'kreport': *.kreport2 路径或None,
          'out_dir': 输出目录}
    """
    cfg = get_config()
    kunpeng = cfg.tool('kunpeng')
    threads = threads or cfg.threads

    d = check_path(db_dir, must_exist=True, in_platform=True)
    if not db_ready(d):
        raise RuntimeError(f"kunpeng 库不可用: {d}")

    out = check_path(out_dir, must_exist=False, in_platform=True)
    os.makedirs(out, exist_ok=True)
    # 清掉上一轮的固定名产物：kunpeng 在 0 条 C 行时不写 output 文件，
    # 而末尾的 glob 只看"文件是否存在"——残留会被当成本次结果返回
    # （实测 results/NX-6/01_host_removal 留有 output_1.kreport2 且无 .done）。
    for _pat in ('output_*.txt', '*.kreport2', '_part*.txt',
                 '_part*.kreport2'):
        for _stale in glob.glob(os.path.join(out, _pat)):
            try:
                os.remove(_stale)
            except OSError:
                pass

    files = [check_path(p, must_exist=True) for p in inputs]
    need = sum(os.path.getsize(f) for f in files) * 16   # 经验系数（含余量）
    # 内存预估：库 hash 表（实测文件大小）+ resolve(≈每条记录 190B) + 基础开销
    _decomp = est_decompressed(files)
    _est_reads = max(1, int(_decomp / (316 if paired else 158)))
    _hash_gb = sum(os.path.getsize(p) for p in glob.glob(os.path.join(d, 'hash_*.k2d'))) / 1e9
    mem_est_gb = _hash_gb + _est_reads * 190 / 1e9 + 0.5
    if logger:
        log_res_plan(logger, 'kunpeng 分类', threads=threads,
                     mem_gb=mem_est_gb, disk_gb=need / 1e9,
                     note=f"内存=库hash表{_hash_gb:.1f}GB+resolve"
                          f"{_est_reads:,}条×190B；磁盘=输入gz×16（经验系数，"
                          f"结束记录实测值）；内存不足自动分片")
    common = ['-p', str(threads), '-g', str(int(min_hit_groups))]
    if paired:
        common.append('-P')
    if confidence and float(confidence) > 0:
        common += ['-T', str(confidence)]

    # chunk 模式（classify）：临时目录 = 用户指定优先，否则自动挑剩余空间足够的盘
    chunk = _pick_chunk_dir(out, inputs, logger=logger, base=chunk_dir)
    # FASTQ→FASTA（seqkit fq2fa，装了 seqkit 才启用）：分类不读质量值
    # （-Q 默认 0），预转换可砍掉质量行解压开销、加速读入。
    # 注意仅用于分类输入——下游组装（SPAdes 纠错）仍用原 FASTQ。
    try:
        seqkit = get_config().tool('seqkit')
    except (FileNotFoundError, RuntimeError):
        seqkit = None
    classify_inputs = list(files)
    conv_dir = None
    # FASTQ→FASTA（seqkit fq2fa，装了 seqkit 才启用）：分类不读质量值
    # （-Q 默认 0），预转换可砍掉 splitr 的质量行解压开销、加速读入。
    # 注意：chunk 中间盘占用由 k2 中间格式决定，与输入 FASTQ/FASTA 关系不大，
    # 空间占用基本不变（每次运行结束都会记录实测值供核对）。
    # 转换产物必须放 chunk 目录之外（kunpeng 要求 chunk-dir 为干净目录，
    # 启动时会清理其中文件）。
    if allow_convert and seqkit and any(str(f).lower().endswith(
            ('.fastq', '.fq', '.fastq.gz', '.fq.gz')) for f in files):
        conv_dir = check_path(
            os.path.join(out, f'_fq2fa_{int(time.time() * 1000)}'),
            must_exist=False, in_platform=True)
        os.makedirs(conv_dir, exist_ok=True)
        converted = []
        for i, f in enumerate(files, 1):
            fa = os.path.join(conv_dir, f'input_{i}.fa.gz')
            run_cmd([seqkit, 'fq2fa', '-w', '0', '-j', str(threads),
                     str(f), '-o', fa], logger=logger)
            converted.append(fa)
        if logger:
            logger.log(f"已用 seqkit fq2fa 预转换 {len(converted)} 个输入"
                       f"（分类不读质量值，砍掉质量行解压开销加速读入）")
        classify_inputs = converted
    # ---- 大输入自动分片：resolve 内存≈读数×190B（与线程数无关）。
    #      空闲提交内存不足时逐片分类（每片内存≈1/N），C 行输出按行合并、
    #      kreport 按 taxid 合并计数（分片读集合不相交，合并精确）。 ----
    total_in_gb = sum(os.path.getsize(str(f)) for f in classify_inputs) / 1e9
    parts = 1
    if total_in_gb > _SPLIT_THRESHOLD_GB:
        budget = max(_free_commit_gb() * 0.35, 4.0)   # 单片 resolve 预算 GB
        parts = min(6, max(2, int(total_in_gb * 1.6 / budget) + 1))
    attempt = parts
    chunk_base = os.path.dirname(chunk)
    chunk_name = os.path.basename(chunk)
    out_parts_txt, out_parts_kr = [], []
    _chunk_actual = 0          # 实测 chunk 磁盘占用（供日志核对预估）
    # conv_dir 是本次运行的 fq2fa 中间产物（可达数 GB）：无论成功、
    # 内存失败还是其他异常，都必须在函数退出前清掉，否则永久残留在
    # 阶段目录里（原实现只在"有 kreport"时删除）。
    _chunk_dirs = []          # 本次创建的 chunk 目录（失败时也要清）
    try:
        while True:
            try:
                part_inputs_list = [classify_inputs]
                if attempt > 1:
                    if conv_dir is None:
                        raise RuntimeError(
                            f"输入 {total_in_gb:.1f} GB 较大且未安装 seqkit，"
                            f"无法分片分类（resolve 需内存约 "
                            f"{total_in_gb * 1.6:.0f} GB）。请将 seqkit.exe "
                            f"放入平台目录后重试")
                    srcs = classify_inputs
                    if len(srcs) == 2:
                        p1 = _split_fasta_gz(srcs[0], attempt, conv_dir, 'mate1')
                        p2 = _split_fasta_gz(srcs[1], attempt, conv_dir, 'mate2')
                        part_inputs_list = [[p1[j], p2[j]] for j in range(attempt)]
                    else:
                        part_inputs_list = [[p] for p in
                                            _split_fasta_gz(srcs[0], attempt,
                                                            conv_dir, 'se')]
                    if logger:
                        logger.log(f"输入 {total_in_gb:.1f} GB，分 {attempt} 片"
                                   f"逐片分类（单片 resolve 内存≈1/{attempt}）")
                out_parts_txt, out_parts_kr = [], []
                for j, pin in enumerate(part_inputs_list):
                    # 分片 chunk 目录位于 _pick_chunk_dir 选定的基目录下
                    # （可能在平台外，如 E:\vp_chunk），无需再做平台内校验
                    chunk_j = (chunk if attempt == 1 else
                               os.path.join(chunk_base, f'{chunk_name}_p{j}'))
                    os.makedirs(chunk_j, exist_ok=True)
                    _chunk_dirs.append(chunk_j)

                    def _watch(chunk_j=chunk_j, j=j,
                               n_parts=len(part_inputs_list),
                               attempt_now=attempt):
                        """看门狗：kunpeng 无进度输出，用 chunk 中间盘增长近似。

                        n_parts / attempt_now 用默认参数显式绑定：本函数会被
                        run_cmd 的监控线程异步调用，若直接闭包读外层
                        part_inputs_list / attempt，下一轮重试（分片数变化）
                        时线程可能读到新值，进度分片序号错位（B023）。
                        """
                        if progress is None:
                            return
                        used = dir_size(chunk_j)
                        frac = min(used / need, 0.95) if need else 0
                        base = j / max(n_parts, 1)
                        span = 1 / max(n_parts, 1)
                        if progress:
                            progress(base + frac * span,
                                     f"分类中间数据已写 {used / 1e9:.1f} / "
                                     f"预估 {need / 1e9:.0f} GB"
                                     + (f"（分片 {j + 1}/{attempt_now}）"
                                        if attempt_now > 1 else ''))

                    cmd = [kunpeng, 'classify', '--db', d, '--chunk-dir', chunk_j,
                           '--output-dir', out] + common + \
                          ['--batch-size', '4'] + pin
                    if attempt > 1 and logger:
                        logger.log(f"分类分片 {j + 1}/{attempt} ...")
                    try:
                        run_cmd(cmd, logger=logger, monitor_fn=_watch if progress
                                else None, monitor_interval=5)
                    except RuntimeError as e:
                        if '3221226505' in str(e) or 'memory allocation' in str(e):
                            raise _MemFail(str(e)[:200])
                        msg = str(e)
                        for _rc, _why in _WIN_EXIT_HINTS.items():
                            if str(_rc) in msg:
                                raise RuntimeError(
                                    f"kunpeng 异常退出（{_why}）。原始错误: "
                                    + msg[:200])
                        raise
                    _chunk_actual += dir_size(chunk_j)
                    o1 = os.path.join(out, 'output_1.txt')
                    k1 = os.path.join(out, 'output_1.kreport2')
                    if os.path.isfile(o1):
                        dst = check_path(os.path.join(out, f'_part{j}_o.txt'),
                                         must_exist=False, in_platform=True)
                        shutil.move(o1, dst)
                        out_parts_txt.append(dst)
                    if os.path.isfile(k1):
                        dst = check_path(os.path.join(out, f'_part{j}.kreport2'),
                                         must_exist=False, in_platform=True)
                        shutil.move(k1, dst)
                        out_parts_kr.append(dst)
                    if not keep_chunk:
                        shutil.rmtree(chunk_j, ignore_errors=True)
                break
            except _MemFail as e:
                if attempt >= 6:
                    raise RuntimeError(
                        "kunpeng 分类内存不足（已自动分片+降 batch-size 重试仍"
                        "失败）。建议：① 关闭占内存的程序后重跑；② 用「子采样」"
                        "减小数据量。原始错误: " + str(e)[:200])
                attempt = min(6, attempt * 2)
                if logger:
                    logger.log(f"分片分类仍内存不足，加密到 {attempt} 片重试 ...")
    finally:
        if conv_dir and os.path.isdir(conv_dir):
            shutil.rmtree(conv_dir, ignore_errors=True)
        # 分片 chunk 目录可达数十 GB：正常路径已逐个删除，这里兜住
        # 内存失败/异常退出时残留的目录（keep_chunk 为调试用途，保留）。
        if not keep_chunk:
            for _cd in _chunk_dirs:
                if os.path.isdir(_cd):
                    shutil.rmtree(_cd, ignore_errors=True)

    # 合并分片输出
    if logger:
        logger.log(f"📊 磁盘核对 [kunpeng 分类] 预估 ~{need / 1e9:.0f} GB，"
                   f"实测 chunk 中间数据 {_chunk_actual / 1e9:.1f} GB"
                   f"（实测/预估 = {_chunk_actual / max(need, 1):.2f}）", "PLAN")
    if len(out_parts_txt) > 1 or (out_parts_txt and
                                  out_parts_txt[0] != os.path.join(
                                      out, 'output_1.txt')):
        with safe_open(os.path.join(out, 'output_1.txt'), 'wt') as w:
            for pf in out_parts_txt:
                with safe_open(pf) as f:
                    shutil.copyfileobj(f, w)
        for pf in out_parts_txt:
            try:
                os.remove(pf)
            except OSError:
                pass
    if out_parts_kr:
        _merge_kreports(out_parts_kr,
                        os.path.join(out, 'output_1.kreport2'), logger)
        for pf in out_parts_kr:
            try:
                os.remove(pf)
            except OSError:
                pass

    kraken = None
    kreport = None
    for p in sorted(glob.glob(os.path.join(out, 'output_*.txt'))):
        kraken = p          # 取最后一个（通常单输入只有一个）
    for p in sorted(glob.glob(os.path.join(out, '*.kreport2'))):
        kreport = p
    return {'kraken': kraken, 'kreport': kreport, 'out_dir': out}


def parse_classify_output(txt_path):
    """解析 Kraken per-read 输出。yield (flag, read_id, taxid, length, path_str)。

    flag: 'C'=classified 'U'=unclassified 'A'=ambiguous(部分工具)
    """
    p = check_path(txt_path, must_exist=True)
    with safe_open(p) as f:
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 4:
                continue
            flag, rid, taxid, length = parts[0], parts[1], parts[2], parts[3]
            path_str = parts[4] if len(parts) > 4 else ''
            try:
                taxid_i = int(taxid) if taxid.isdigit() else 0
            except ValueError:
                taxid_i = 0
            yield flag, rid, taxid_i, length, path_str


RANK_CODE = {
    'U': 'unclassified', 'R': 'root', 'R1': 'realm', 'R2': 'realm2',
    'D': 'domain', 'K': 'superkingdom', 'P': 'phylum', 'C': 'class',
    'O': 'order', 'F': 'family', 'G': 'genus', 'G1': 'subgenus',
    'S': 'species', 'S1': 'subspecies', '-': 'no rank',
}


def parse_kreport(kreport_path):
    """解析 kreport2 -> [ {percent, frags, rank, rank_code, taxid, name,
    depth, parent_taxid}, ... ]

    kreport 的 name 列保留缩进（每 2 空格 = 1 层），据此重建父子关系。
    rank 短码（S/G/F…）标准化为全名（species/genus/family…）。
    """
    p = check_path(kreport_path, must_exist=True)
    rows = []
    stack = []          # [(depth, taxid), ...] 当前路径
    with safe_open(p) as f:
        for line in f:
            line = line.rstrip('\n')
            parts = line.split('\t')
            if len(parts) >= 6:
                # 标准六列: percent, frags, cum_frags, rank, taxid, name(含缩进)
                try:
                    percent = float(parts[0])
                    frags = int(parts[1])
                    code = parts[3].strip()
                    taxid = int(parts[4])
                    raw_name = parts[5]
                except (ValueError, IndexError):
                    continue
            else:
                # 空格分隔的 kraken2 原生格式兜底
                sp = line.split()
                if len(sp) < 6:
                    continue
                try:
                    percent = float(sp[0])
                    frags = int(sp[1])
                    code = sp[3]
                    taxid = int(sp[4])
                    raw_name = line.split(sp[4], 1)[1] if sp[4] in line else sp[5]
                except (ValueError, IndexError):
                    continue
            name = raw_name.strip()
            depth = (len(raw_name) - len(raw_name.lstrip(' '))) // 2
            while stack and stack[-1][0] >= depth:
                stack.pop()
            parent_taxid = stack[-1][1] if stack else 0
            rows.append({'percent': percent, 'frags': frags,
                         'rank_code': code, 'rank': RANK_CODE.get(code, code),
                         'taxid': taxid, 'name': name, 'depth': depth,
                         'parent_taxid': parent_taxid})
            stack.append((depth, taxid))
    return rows


class KreportTree:
    """从 kreport 构建的轻量分类树（避免全量加载 NCBI taxonomy）。

    提供 name(taxid) / rank(taxid) / parent(taxid) / lineage_names(taxid) /
    get_ancestor_at_rank(taxid, rank)。
    """

    def __init__(self, rows):
        self.rows = rows
        self.by_taxid = {r['taxid']: r for r in rows}

    @classmethod
    def from_kreport(cls, kreport_path):
        return cls(parse_kreport(kreport_path))

    def name(self, taxid):
        r = self.by_taxid.get(int(taxid))
        return r['name'] if r else f'taxid:{taxid}'

    def rank_of(self, taxid):
        r = self.by_taxid.get(int(taxid))
        return r['rank'] if r else ''

    def parent(self, taxid):
        r = self.by_taxid.get(int(taxid))
        return r['parent_taxid'] if r else 0

    def lineage_names(self, taxid):
        chain, t, seen = [], int(taxid), set()
        for _ in range(40):
            r = self.by_taxid.get(t)
            if r is None or t in seen:
                break
            seen.add(t)
            chain.append(r['name'])
            t = r['parent_taxid']
            if t == 0:
                break
        return list(reversed(chain))

    def get_ancestor_at_rank(self, taxid, rank):
        t = int(taxid)
        for _ in range(40):
            r = self.by_taxid.get(t)
            if r is None:
                return None
            if r['rank'] == rank:
                return t
            t = r['parent_taxid']
            if t == 0:
                break
        return None


def convert_kraken2(k2_source, db_dir, hash_capacity='1G', logger=None):
    """Kraken2 库 → kunpeng 分片库（方式 C：kun_peng hashshard）。

    k2_source: Kraken2 库目录（含 hash.k2d/opts.k2d/taxo.k2d）
               或 .tar.gz/.tar 包（kraken2 官方预构建库下载格式）。
    db_dir:    目标 kunpeng 库目录（hashshard 输出 hash_*.k2d 等）。
    返回 db_dir。库就绪判定与自建库一致（db_ready）。
    """
    import tarfile
    cfg = get_config()
    kunpeng = cfg.tool('kunpeng')
    db_dir = check_path(db_dir, must_exist=False, in_platform=True)
    os.makedirs(db_dir, exist_ok=True)

    src = check_path(k2_source, must_exist=True, in_platform=True)
    tmp = None
    if os.path.isfile(src):
        # 包：解到临时目录
        if not src.lower().endswith(('.tar.gz', '.tgz', '.tar')):
            raise RuntimeError('Kraken2 库包仅支持 .tar.gz/.tar')
        tmp = check_path(os.path.join(db_dir, '_k2_extract'),
                         must_exist=False, in_platform=True)
        os.makedirs(tmp, exist_ok=True)
        if logger:
            logger.log(f"解包 Kraken2 库包: {os.path.basename(src)}")
        with tarfile.open(src, 'r:*') as tf:
            # Tar Slip 防护：库包可能来自第三方下载，先显式拒绝绝对路径、
            # .. 段与盘符成员（给出明确报错），再以 filter='data' 解包
            # （Python 3.12+ 拒绝绝对路径/外部链接成员，防写出平台目录之外）
            for m in tf.getmembers():
                name = m.name.replace('\\', '/')
                head = name.split('/')[0]
                if (name.startswith('/') or '..' in name.split('/')
                        or (len(head) == 2 and head[1] == ':')):
                    raise RuntimeError(f'Kraken2 库包含非法成员路径: {m.name}')
            tf.extractall(tmp, filter='data')
        # 找到含 hash.k2d 的目录（可能嵌套一层）
        k2_dir = None
        for cur, _sub, fns in os.walk(tmp):
            if 'hash.k2d' in fns and 'opts.k2d' in fns and 'taxo.k2d' in fns:
                k2_dir = cur
                break
        if not k2_dir:
            raise RuntimeError('解包后未找到 Kraken2 索引文件 '
                               '(hash.k2d/opts.k2d/taxo.k2d)')
    else:
        k2_dir = src
        for req in ('hash.k2d', 'opts.k2d', 'taxo.k2d'):
            if not os.path.isfile(os.path.join(k2_dir, req)):
                raise RuntimeError(f'Kraken2 库缺少 {req}（目录: {k2_dir}）')

    try:
        if logger:
            logger.log(f"kunpeng hashshard（就地转换）: {k2_dir} "
                       f"(hash-capacity {hash_capacity})")
        run_cmd([kunpeng, 'hashshard', '--db', k2_dir,
                 '--hash-capacity', hash_capacity], logger=logger)

        # hashshard 输出在 k2 目录内（hash_*.k2d / hash_config.k2d），
        # 复制（不是搬移）kunpeng 运行所需文件到目标 db_dir：源目录可能是
        # 用户的 Kraken2 官方库（opts/taxo.k2d 是它的原件），move 会把源库
        # 拆掉（下次还要用它做别的分析就没了）。复制后源库保持完整，
        # 多花的只是索引文件的磁盘副本（hash 分片与源库同量级）。
        import shutil
        need = [('hash_*.k2d', 'glob'), ('hash_config.k2d', 'file'),
                ('opts.k2d', 'file'), ('taxo.k2d', 'file')]
        copied = 0
        for pat, kind in need:
            if kind == 'glob':
                for f in glob.glob(os.path.join(k2_dir, pat)):
                    shutil.copyfile(f, os.path.join(db_dir,
                                                    os.path.basename(f)))
                    copied += 1
            else:
                f = os.path.join(k2_dir, pat)
                if os.path.isfile(f):
                    shutil.copyfile(f, os.path.join(db_dir, pat))
                    copied += 1
        if logger:
            logger.log(f"复制 kunpeng 库文件 {copied} 个 → {db_dir}"
                       f"（源库保持不变，确认不再需要后可手动清理）")

        # taxonomy dmp 也复制（classify 的 LCA 用，部分流程需要）
        for dmp in ('nodes.dmp', 'names.dmp'):
            f = os.path.join(k2_dir, dmp)
            if os.path.isfile(f) and not os.path.isfile(os.path.join(db_dir, dmp)):
                shutil.copyfile(f, os.path.join(db_dir, dmp))
    finally:
        # 清理解包临时目录（finally：hashshard/搬运中途失败也不残留数 GB 解包产物）
        if tmp:
            shutil.rmtree(tmp, ignore_errors=True)

    if not db_ready(db_dir):
        raise RuntimeError('hashshard 完成但库文件不齐全（db_ready 失败）')
    if logger:
        logger.log("Kraken2 库转换完成，kunpeng 可直接分类")
    return db_dir

