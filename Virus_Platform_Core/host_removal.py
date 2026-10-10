# -*- coding: utf-8 -*-
"""
阶段① 宿主序列去除：kunpeng 宿主库分类 reads → 剔除宿主 read 对 → 输出保留 reads。
宿主库仅含宿主单一物种，因此 classify 输出中 C 行（已分类）即宿主 read 对。
"""
import os
import json
import subprocess

from .utils import check_path, safe_open, iter_fastq_records, is_step_done, mark_step_done
from .kunpeng import classify, parse_classify_output


def _read_id(name_line):
    """@A00123:... 1:N:0:... -> A00123:...（取第一空白前 token，去 @）"""
    tok = name_line.split(None, 1)[0]
    return tok[1:] if tok.startswith('@') else tok


def _strip_pe_suffix(rid):
    """去掉 read 名尾部的 /1 或 /2 配对后缀。

    MGI/华大 FASTQ 头形如 @LH00278:...:26992/1（后缀与 ID 同 token，无空格），
    而 kunpeng/kraken 输出的 read 名不带该后缀。宿主去除与病毒提取都是按
    ID 精确比对，不剥后缀会导致两侧永不匹配（实测交集为 0）。
    Illumina 新格式（@A001:/1 空格 1:N:0:）的后缀在第二 token，_read_id
    取第一 token 时已天然丢掉，此处对无后缀 ID 是恒等变换。
    """
    if rid.endswith('/1') or rid.endswith('/2'):
        return rid[:-2]
    return rid


def _norm_read_id(name_line):
    """read 名 -> 与 kunpeng 输出对齐的规范 ID（去 @、去空白后段、去 /1 /2）。"""
    return _strip_pe_suffix(_read_id(name_line))


def _finfo(path):
    """文件指纹（供 stats.json 记录真实分类输入，便于事后复核）。"""
    p = str(path)
    try:
        st = os.stat(p)
        return {'name': os.path.basename(p), 'size': int(st.st_size),
                'mtime': int(st.st_mtime_ns)}
    except OSError:
        return {'name': os.path.basename(p), 'size': None, 'mtime': None}


def _id_patterns(ids, pe_suffix=''):
    """把 ID 集合转成 seqkit 模式文件行（精确匹配用，逐行一个）。

    pe_suffix: 追加到每个 ID 尾部的配对后缀，''/'/1'/'/2'。

    为什么不用 -r 正则：seqkit 精确匹配是 hash 查表，与模式数无关
    （452,400 条模式实测 2.4s）；-r 是逐 read 逐模式的正则扫描，
    5,000 条模式就要 197s，45 万条会直接卡死。故宁可先探测输入 FASTQ
    的真实后缀、把后缀写进模式文件，也要保住精确匹配。
    """
    return [f'{i}{pe_suffix}' for i in ids]


def detect_pe_suffix(fastq_path, logger=None, sample_n=100):
    """探测 FASTQ 记录头的配对后缀，返回 ''/'/1'/'/2'。

    MGI/华大：@LH00278:...:26992/1（单 token 带后缀）
    Illumina 新版：@A00123:... 1:N:0:...（后缀在第二 token，首 token 无）
    SRA 下载：@SRR39909438.9 A00459:.../1（首 token 无，后缀在第二 token）
    后两种取首 token 时天然无后缀，返回 ''即精确匹配原 ID。

    采样前 sample_n 条记录（只取每条记录的头行，跳过序列/质量行——质量串
    里也可能出现 "/1" 结尾）统计后缀分布：正常数据全一致；混杂时 WARN，
    但仍按首条结果返回（保持既有探测口径兼容，避免突然换行为）。
    """
    counts = {}
    first_suffix = ''
    try:
        with safe_open(fastq_path) as f:
            heads = 0
            idx = 0                      # 非空行序号：FASTQ 头落在 0,4,8…
            for line in f:
                s = line.rstrip('\n')
                if not s:
                    continue
                if idx % 4 == 0:
                    tok = s.split(None, 1)[0]
                    suf = tok[-2:] if (tok.endswith('/1') or
                                       tok.endswith('/2')) else ''
                    if heads == 0:
                        first_suffix = suf
                    counts[suf] = counts.get(suf, 0) + 1
                    heads += 1
                    if heads >= sample_n:
                        break
                idx += 1
    except Exception as e:
        # 读失败必须有痕迹：静默返回空后缀会退化成「无后缀精确匹配」，
        # 宿主过滤 0 条的最终告警才兜得住（2026-10-04 审查 P3）
        if logger:
            logger.log(f"配对后缀探测失败（按无后缀处理）: "
                       f"{type(e).__name__}: {e}", "WARN")
        return ''
    if logger and len(counts) > 1:
        dist = '，'.join(f'{k or "(无)"} ×{v}' for k, v in counts.items())
        logger.log(f"配对后缀探测：{os.path.basename(str(fastq_path))} 前 "
                   f"{sum(counts.values())} 条记录后缀不一致（{dist}），"
                   f"按首条取 {first_suffix or '(无)'}——请确认 R1/R2 文件"
                   f"配对正确", "WARN")
    return first_suffix


def collect_host_read_ids(kraken_out):
    """收集判定为宿主的 read id 集合（C 行）。

    统一剥掉 /1 /2 配对后缀：过滤侧（python 路径 `_norm_read_id`、seqkit 路径
    `_id_patterns` 追加探测后缀）都按"无后缀 ID"工作，若这里保留后缀，
    seqkit 模式会变成 `id/1/1`、python 侧则永不相等——两条路径同时失配且
    无异常（实测 dropped=0 静默漏掉全部宿主 reads）。
    """
    host_ids = set()
    total = 0
    for flag, rid, _taxid, _len, _path in parse_classify_output(kraken_out):
        total += 1
        if flag == 'C':
            host_ids.add(_strip_pe_suffix(rid))
    return host_ids, total


def _seqkit_filter(seqkit, r1, r2, out_r1, out_r2, ids, threads, work_dir,
                   logger=None, progress=None):
    """seqkit grep -v 反选加速：写出不在 ids 中的 reads（宿主去除的 kept）。

    kraken 输出只含 C 行（宿主 IDs），kept = 总数 - C 数；
    kunpeng 只输出 C 行，总数由 U/C 行 + kreport 未分类数还原不可行，
    故 total 由 kept + dropped 反推（dropped = C 行 id 数）。

    后缀处理：MGI 数据 FASTQ 头带 /1 /2 后缀（与 ID 同 token），而 kunpeng
    输出无后缀；先探测每个输入的实陒后缀，把它写进模式文件，从而保住
    seqkit 的精确匹配（hash 查表）。R1/R2 各用自身后缀，不会互相误匹。
    返回统计 dict；失败抛异常由调用方回退。
    """
    from .utils import run_cmd
    ids = list(ids)

    def _count_total(path):
        r = subprocess.run(
            [seqkit, 'stats', '-T', check_path(path, must_exist=True)],
            capture_output=True, text=True, encoding='utf-8',
            errors='replace')
        if r.returncode != 0:
            raise RuntimeError(f"seqkit stats 失败: {r.stderr[-200:]}")
        lines = r.stdout.strip().splitlines()
        hdr, val = lines[0].split('\t'), lines[1].split('\t')
        return int(dict(zip(hdr, val))['num_seqs'])

    srcs = ([(r1, out_r1), (r2, out_r2)] if r2 else [(r1, out_r1)])
    n_in1 = kept1 = dropped1 = kept2 = 0
    pats_made = []
    try:
        for i, (src, dst) in enumerate(srcs):
            suf = detect_pe_suffix(src, logger=logger)
            if progress:
                progress(0.15 + i * 0.4,
                         f"seqkit 反选 {os.path.basename(str(src))}")
            pat_i = check_path(
                os.path.join(work_dir, f'_host_ids_{i + 1}.txt'),
                must_exist=False, in_platform=True)
            with safe_open(pat_i, 'wt') as f:
                f.write('\n'.join(_id_patterns(ids, pe_suffix=suf)) + '\n')
            pats_made.append(pat_i)
            if logger:
                logger.log(f"后缀探测 {os.path.basename(str(src))}: "
                           f"{suf or '(无)'}")
            run_cmd([seqkit, 'grep', '-v', '-f', str(pat_i),
                     '-j', str(threads or 4),
                     '-o', check_path(dst, must_exist=False,
                                      in_platform=True),
                     check_path(src, must_exist=True)], logger=logger)
            # 保留数直接对**输出**统计：宿主 reads 通常占大头，输出远小于
            # 输入——比原先再对输入跑一遍 grep -v -C 计数省一整遍大文件扫描
            if i == 0:
                n_in1 = _count_total(src)
                kept1 = _count_total(dst)
                dropped1 = n_in1 - kept1
            else:
                kept2 = _count_total(dst)
        if r2 and kept1 != kept2:
            raise RuntimeError(
                f"R1/R2 kept 数不一致（{kept1} vs {kept2}），"
                f"输入可能不配对")
        if not r2:
            kept2 = kept1
        # dropped 用实测差值（不再用 len(ids) 充数）
        return {'total_pairs': n_in1, 'kept_pairs': kept1,
                'dropped_pairs': dropped1, 'kept_mate2': kept2}
    finally:
        for p in pats_made:
            try:
                os.remove(p)
            except OSError:
                pass


def _passthrough_fastq(r1, r2, out_r1, out_r2, logger=None):
    """0 条宿主 read（合法且常见）：整文件直通，避免逐条解压重压一遍。

    仍做 R1/R2 条数一致性校验（读取即校验完整性，截断文件会抛错）。
    """
    import shutil
    for src, dst in ((r1, out_r1), (r2, out_r2)):
        if src and dst:
            shutil.copyfile(check_path(src, must_exist=True),
                            check_path(dst, must_exist=False, in_platform=True))

    def _count(path):
        n = 0
        with safe_open(path) as f:
            for _rec in iter_fastq_records(f):
                n += 1
        return n

    n1 = _count(out_r1)
    kept2 = _count(out_r2) if r2 else n1
    if r2 and n1 != kept2:
        raise RuntimeError(f"R1/R2 read 数不一致（{n1} vs {kept2}），"
                           f"输入可能不配对")
    if logger:
        logger.log("0 条宿主 read：整文件直通（复制，不做逐条过滤）")
    return {'total_pairs': n1, 'kept_pairs': n1, 'dropped_pairs': 0,
            'kept_mate2': kept2, 'filter_mode': 'passthrough'}


def filter_host_fastq(r1, r2, out_r1, out_r2, drop_ids, logger=None,
                      threads=None, work_dir=None, progress=None):
    """宿主过滤：优先 seqkit -v 反选（多线程），失败回退纯 Python 流式。

    返回统计 dict（total_pairs/kept_pairs/dropped_pairs/filter_mode）。
    """
    if not drop_ids:
        return _passthrough_fastq(r1, r2, out_r1, out_r2, logger=logger)
    if work_dir is not None:
        try:
            from .config import get_config
            seqkit = get_config().tool('seqkit')
        except Exception:
            seqkit = None
        if seqkit:
            try:
                if logger:
                    logger.log(f"seqkit grep -v 加速过滤"
                               f"（宿主 IDs {len(drop_ids):,} 条）")
                res = _seqkit_filter(seqkit, r1, r2, out_r1, out_r2,
                                     drop_ids, threads, work_dir,
                                     logger=logger, progress=progress)
                res['filter_mode'] = 'seqkit'
                return res
            except Exception as e:
                if logger:
                    logger.log(f"seqkit 过滤异常，回退内置过滤 "
                               f"[{type(e).__name__}]: {e}", "WARN")
    if r2:
        res = filter_paired_fastq(r1, r2, out_r1, out_r2, drop_ids,
                                  logger=logger, progress=progress)
    else:
        res = filter_single_fastq(r1, out_r1, drop_ids, logger=logger,
                                  progress=progress)
    res['filter_mode'] = 'python'
    return res


def _est_total_pairs(paths):
    """按解压字节数粗估 reads 对数（进度条分母用，允许偏差）。"""
    from .utils import est_decompressed
    decomp = est_decompressed([p for p in paths if p])
    per_pair = 2 * (150 + 10)     # 2×(读长≈150 + 头行分摊)
    return max(int(decomp / per_pair), 1)


def filter_paired_fastq(r1, r2, out_r1, out_r2, ids, keep=False, logger=None,
                        progress_every=2_000_000, progress=None):
    """同步遍历双端 FASTQ 过滤 read 对。

    keep=False: 剔除 ids 中的 read 对（去宿主用，保留非宿主）；
    keep=True : 仅保留 ids 中的 read 对（提取病毒 reads 用）。
    返回统计。
    """
    kept = total = 0
    est_total = _est_total_pairs([r1]) if progress else 0
    with safe_open(r1) as f1, safe_open(r2) as f2, \
            safe_open(out_r1, 'wt') as w1, safe_open(out_r2, 'wt') as w2:
        it1, it2 = iter(f1), iter(f2)
        while True:
            rec1 = next(it1, None)
            rec2 = next(it2, None)
            if rec1 is None and rec2 is None:
                break
            if rec1 is None or rec2 is None:
                raise ValueError("R1/R2 read 数不一致（文件不配对）")
            # 各读 3 行补齐 4 行记录
            l1 = [rec1] + [next(it1, '') for _ in range(3)]
            l2 = [rec2] + [next(it2, '') for _ in range(3)]
            # 截断的 FASTQ（下载中断/尾部损坏）会补出空串并静默写出非法
            # 记录；utils.iter_fastq_records 对同样输入是抛错的，口径要一致。
            if not l1[3] or not l2[3]:
                raise ValueError(
                    f"FASTQ 记录不完整（{_norm_read_id(l1[0])}）："
                    f"文件可能被截断或 gz 尾部损坏")
            total += 1
            rid = _norm_read_id(l1[0])
            hit = rid in ids
            if hit == keep:
                w1.write(''.join(l1))
                w2.write(''.join(l2))
                kept += 1
            if progress and total % 200_000 == 0:
                progress(min(total / est_total, 0.98),
                         f"已处理 {total:,} 对（保留 {kept:,}）")
            if logger and progress_every and total % progress_every == 0:
                logger.log(f"  过滤进度: {total:,} 对 (保留 {kept:,})")
    return {'total_pairs': total, 'kept_pairs': kept,
            'dropped_pairs': total - kept}


def filter_single_fastq(r1, out_r1, drop_ids, logger=None, progress=None):
    """单端 FASTQ：剔除 drop_ids 中的 reads，返回统计。"""
    kept = total = 0
    est_total = _est_total_pairs([r1]) if progress else 0
    with safe_open(r1) as f1, safe_open(out_r1, 'wt') as w1:
        for rec in iter_fastq_records(f1):
            total += 1
            if _norm_read_id(rec[0]) not in drop_ids:
                w1.write(''.join(rec))
                kept += 1
            if progress and total % 200_000 == 0:
                progress(min(total / est_total, 0.98),
                         f"已处理 {total:,} 条（保留 {kept:,}）")
    return {'total_pairs': total, 'kept_pairs': kept,
            'dropped_pairs': total - kept}


def remove_host(sample_dir, r1, r2, db_host, threads=None, confidence=0.0,
                logger=None, force=False, chunk_dir=None,
                classify_r1=None, classify_r2=None, allow_convert=True,
                progress=None, min_hit_groups=2):
    """执行宿主去除（双端或单端）。返回阶段结果 dict（含 kept fastq 路径与统计）。

    classify_r1/classify_r2: 预处理阶段产出的 FASTA（可选）——分类步骤直接用
    FASTA（省去质量行解压，中间数据更小），过滤仍以 FASTQ 为准写出 kept reads。
    progress: 可选 progress(pct, msg) 供 GUI 进度条。
    """
    step = 'host_removal'
    out_dir = check_path(os.path.join(sample_dir, '01_host_removal'),
                         must_exist=False, in_platform=True)
    os.makedirs(out_dir, exist_ok=True)

    stats_file = os.path.join(out_dir, 'stats.json')
    kept_r1 = os.path.join(out_dir, 'kept_R1.fastq.gz')
    kept_r2 = os.path.join(out_dir, 'kept_R2.fastq.gz') if r2 else None
    if is_step_done(out_dir, step) and not force:
        # 复用校验（审查 P1-1）：旧实现只看 .done 就返回旧 stats——分类
        # 输入、宿主库、confidence/min_hit_groups 任一变化都必须重跑，
        # 否则换库/换输入后仍沿用上一轮的宿主占比与 kept reads。
        try:
            with safe_open(stats_file) as f:
                prev = json.load(f)
            _cur_fp = [_finfo(p) for p in
                       (classify_r1 or r1, classify_r2 or r2) if p]
            _cur_params = {'confidence': confidence,
                           'min_hit_groups': min_hit_groups}
            _prev_params = prev.get('classify_params') or {}
            _prev_db = (prev.get('host_db') or {}).get('path')
            # 口径对齐：stats 里记的是 host_db_info 的 abspath 规范化路径，
            # 本次入参可能是相对路径——两种写法都认，否则相对路径配置下
            # 永远判"库已切换"、每次都重跑（断点续跑被废）。
            _db_variants = {str(db_host),
                            os.path.normpath(os.path.abspath(str(db_host)))}
            _fp_ok = (prev.get('classify_inputs') == _cur_fp)
            _pm_ok = all(_prev_params.get(k) == v
                         for k, v in _cur_params.items())
            _db_ok = (_prev_db in (None, '',) or _prev_db in _db_variants)
        except (OSError, ValueError):
            _fp_ok = _pm_ok = _db_ok = False
        if _fp_ok and _pm_ok and _db_ok:
            if logger:
                logger.log("阶段①宿主去除已完成（输入/库/参数未变化），跳过")
            return prev
        if logger:
            _why = []
            if not _fp_ok:
                _why.append('分类输入已变化')
            if not _pm_ok:
                _why.append('confidence/min_hit_groups 已变化')
            if not _db_ok:
                _why.append(f'宿主库已切换（{_prev_db} → {db_host}）')
            logger.log('阶段①存在断点，但' + '、'.join(_why)
                       + '——重跑本阶段（不复用旧结果）', 'WARN')

    cls_in = [p for p in (classify_r1 or r1, classify_r2 or r2) if p]
    # 宿主库身份：记录 + 前置告警。
    # 为什么要在这里查：库本体不记录物种，用错宿主去做宿主去除是**静默**的
    # （历史事故：枸杞基因组建成了 taxid=4081 番茄库，map/taxo 长期带两个
    # taxid）。这里把"本次实际用的库是谁"写进 stats.json，让结果可追溯；
    # 身份不明或多 taxid 时明确告警，而不是照跑不误。
    from .config import host_db_info
    hdb = host_db_info(db_host)
    if logger:
        who = (f"{hdb.get('species') or '?'} (taxid={hdb.get('taxid')})"
               if hdb.get('taxid') else "身份未记录（无 host_db.json 清单）")
        logger.log(f"阶段① 宿主去除: kunpeng 宿主库分类"
                   f"（{'双端' if r2 else '单端'}，"
                   f"分类输入 {'FASTA' if classify_r1 else 'FASTQ'}）")
        logger.log(f"  宿主库 = {hdb.get('name') or db_host} | 物种 = {who}"
                   + (f" | ⚠ map 内含多个 taxid {hdb['taxids_in_map']}"
                      if hdb.get('conflicted') else ""))
        if hdb.get('conflicted'):
            logger.log("⚠ 该宿主库的 seqid2taxid.map 含多个 taxid，元数据已被"
                       "历史建库污染；若结果异常，请按物种目录重建宿主库", "WARN")
        elif not hdb.get('has_manifest'):
            logger.log("⚠ 该宿主库没有 host_db.json 清单，无法核对物种；"
                       "重建后会自动生成", "WARN")

    def _cls_prog(pct, msg):
        if progress:
            # 分类占本阶段 5%~88%，过滤占 88%~98%
            progress(0.05 + pct * 0.83, f'① 宿主去除：{msg}')

    res = classify(db_host, cls_in, out_dir, paired=bool(r2),
                   threads=threads, confidence=confidence,
                   min_hit_groups=min_hit_groups, logger=logger,
                   chunk_dir=chunk_dir, allow_convert=allow_convert,
                   progress=_cls_prog)

    # C 行为 0 时 kunpeng 不写 output 文件（全部未命中宿主库的合法场景）
    if res['kraken']:
        host_ids, n_records = collect_host_read_ids(res['kraken'])
    else:
        host_ids, n_records = set(), 0
    if logger:
        logger.log(f"分类记录 {n_records:,} 条；宿主 read 对 {len(host_ids):,} 个 "
                   f"({len(host_ids) / max(n_records, 1) * 100:.2f}%)")

    if logger:
        logger.log("过滤 FASTQ 保留非宿主 reads ...")
    if progress:
        progress(0.9, '① 宿主去除：过滤 FASTQ（seqkit 多线程反选）')
    filt = filter_host_fastq(r1, r2, kept_r1, kept_r2, host_ids,
                             logger=logger, threads=threads, work_dir=out_dir,
                             progress=(lambda p, m: progress(0.9 + p * 0.09,
                                                             f'① 过滤：{m}'))
                             if progress else None)
    if progress:
        progress(0.99, '① 宿主去除：写统计')

    stats = {
        'stage': step,
        'single_end': not bool(r2),
        'input_r1': os.path.basename(str(r1)),
        'input_r2': os.path.basename(str(r2)) if r2 else None,
        # 真实分类输入（可能是 ⓪b 的 conv_*.fa.gz，而非 input_r1/r2）：
        # 不记录会导致宿主占比不可复核（历史事故：分类用旧子集、过滤用全量，
        # 但 stats 只写 input_r1，事后无法看出分母分子错配）。
        'classify_inputs': [_finfo(p) for p in cls_in],
        'classify_params': {'confidence': confidence, 'threads': threads,
                            'min_hit_groups': min_hit_groups},
        'kraken_out': os.path.basename(res['kraken']) if res['kraken'] else None,
        'kreport': os.path.basename(res['kreport']) if res['kreport'] else None,
        # 宿主库溯源：谁被用来做宿主去除。原先 stats 里完全没有这项，
        # 出问题时无法判断"是不是用错宿主库了"。
        'host_db': {
            'name': hdb.get('name'),
            'path': hdb.get('path') or str(db_host),
            'taxid': hdb.get('taxid'),
            'species': hdb.get('species'),
            'taxids_in_map': hdb.get('taxids_in_map'),
            'conflicted': bool(hdb.get('conflicted')),
            'has_manifest': bool(hdb.get('has_manifest')),
        },
        'kept_r1': 'kept_R1.fastq.gz',
        'kept_r2': 'kept_R2.fastq.gz' if r2 else None,
        **filt,
        'host_ratio': round(filt['dropped_pairs'] / max(filt['total_pairs'], 1), 6),
    }
    if filt['kept_pairs'] == 0:
        # 失败态也落 stats（供排障），但显式标 ok=false：
        # 管道总览/清单据此显示「失败」而不是「已完成」。
        stats['ok'] = False
        stats['error'] = ('宿主去除后 0 条 reads 保留——输入可能为空、'
                          'R1/R2 不配对，或宿主库选错'
                          '（如用细菌宿主库处理植物数据）')
    with safe_open(stats_file, 'wt') as f:
        json.dump(stats, f, ensure_ascii=False, indent=2)
    if logger:
        logger.log(f"阶段① 完成: 保留 {filt['kept_pairs']:,}/{filt['total_pairs']:,} 对 "
                   f"(宿主占比 {stats['host_ratio'] * 100:.2f}%)")
        # 0 条被去除 = 宿主库与样品宿主可能不匹配（也可能是纯化病毒样品，
        # 所以只告警不报错）。原先这种情况完全静默，下游 ②b 会拿到满屏宿主 reads。
        if filt['dropped_pairs'] == 0 and filt['total_pairs'] > 0:
            if host_ids:
                logger.log(
                    f"⚠ kunpeng 判出 {len(host_ids):,} 个宿主 read 但过滤掉 0 条"
                    f"——疑似 read ID 口径不匹配（过滤模式 "
                    f"{filt.get('filter_mode')}）；请检查 FASTQ 头格式后重跑",
                    "WARN")
            else:
                logger.log("⚠ 没有任何 read 被判为宿主：请确认宿主库物种与样品宿主"
                           "一致（纯化病毒/无宿主污染样品也会如此）", "WARN")
    # 先校验再落断点标记：否则 0 reads 的失败态会留下 .done，重跑时
    # 281-285 行直接返回 stats.json，阶段被当成"成功"并让 ②③ 吃空输入。
    if filt['kept_pairs'] == 0:
        raise RuntimeError(stats['error'] + '。请核对输入文件与宿主库后重跑。')
    mark_step_done(out_dir, step)
    return stats
