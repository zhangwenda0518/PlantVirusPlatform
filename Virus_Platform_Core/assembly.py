# -*- coding: utf-8 -*-
"""
阶段③ 组装与 contig 分类：
- SPAdes 组装去宿主 reads（--rnaviral 默认，可选 --metaviral / --meta / --rna /
  --isolate / srna）
- srna（小RNA/siRNA）模式：--only-assembler + 奇数 k 序列 + 单端 -s
  （对齐 VirusDetect 的 sRNA 组装策略）；reads 最长 <32bp 时其他模式自动切换
- contig 长度过滤
- kunpeng 病毒库分类 contigs
- 汇总病毒 contig 表
"""
import os
import re
import csv
import glob
import shutil

from .config import get_config, DIRS, ref_annotation_path
from .utils import (check_path, safe_open, run_cmd, iter_fasta, write_fasta_record,
                    count_fasta_seqs, is_step_done, mark_step_done)
from .kunpeng import classify, parse_classify_output, parse_kreport, KreportTree

def _ascii_work_base(name='vp_blast'):
    """BLAST(LMDB)/mafft 不支持中文路径：返回纯 ASCII 持久工作目录。

    候选顺序：%TEMP%（ASCII 时）→ 系统盘根 vp_ascii_tmp → C:\\Windows\\Temp。
    C:\\Windows\\Temp 普通用户通常无写权限，只作最后兜底；全部失败时抛
    带指引的错误，而不是让 makedirs 在半路抛难懂的 PermissionError。
    """
    import tempfile
    cands = []
    t = os.path.abspath(tempfile.gettempdir())
    if all(ord(c) < 128 for c in t):
        cands.append(t)
    drv = (os.environ.get('SYSTEMDRIVE') or 'C:').rstrip('\\') or 'C:'
    cands.append(os.path.join(drv + os.sep, 'vp_ascii_tmp'))
    cands.append(r'C:\Windows\Temp')
    last_err = None
    for cand in cands:
        if not os.path.isdir(cand):
            continue
        d = os.path.join(cand, name)
        try:
            os.makedirs(d, exist_ok=True)
            return d
        except OSError as e:
            last_err = e
            continue
    raise RuntimeError(
        f'找不到可写的纯 ASCII 临时目录（外部工具不支持中文路径）。'
        f'请设置环境变量 TEMP 指向纯 ASCII 路径后重试（最后错误: {last_err}）')

def find_virus_ref_fasta():
    """定位病毒参考 FASTA。

    解析顺序：病毒源参考根 databases/virusref_db 下第一个 .fasta/.fa/.fna
    → 默认鉴定库 plant_ref/reference.fasta（库自包含布局后的常驻位置；
    2026-09-29 前默认库目录名为 kv_index）。
    """
    src = check_path(DIRS['virus_src'], must_exist=True)
    for pat in ('*.fasta', '*.fa', '*.fna', '*.fas'):
        hits = sorted(glob.glob(os.path.join(src, pat)))
        if hits:
            return check_path(hits[0], must_exist=True)
    # 库自包含布局（kv_index/reference.fasta）：顶层散置 FASTA 已不存在
    try:
        from .kv_stage import default_reference
        ref = default_reference()
        if ref and os.path.isfile(ref):
            return check_path(ref, must_exist=True)
    except ImportError:
        pass
    raise FileNotFoundError(f"{src} 下未找到病毒参考 FASTA")

# ------------------------------------------------------------------ BLAST 库
def _blast_db_sig(path):
    """库版本戳：源 FASTA 的 `大小|mtime`。取不到返回 ''（视为需重建）。"""
    try:
        st = os.stat(path)
        return f'{st.st_size}|{int(st.st_mtime)}'
    except OSError:
        return ''


def _db_present(prefix):
    """库文件是否齐（v4: .nhr/.nin/.nsd…；v5: .ndb/.njs/.not…，都是 `.n??`）。"""
    return bool(glob.glob(prefix + '.n??'))


def _drop_db(prefix):
    """删掉库的全部 .n* 文件（含 v4 遗留的 `.njs.lmdb-disabled`）。"""
    for p in glob.glob(prefix + '.n*'):
        try:
            os.remove(p)
        except OSError:
            pass


def _disable_lmdb_journal(prefix):
    """v4 库 makeblastdb 也会顺带生成 `.njs`；改名禁用它。

    留着不影响使用，但改名更保险（同 viroids_v4 的做法，见
    docs/DEVELOPMENT_NOTES.md 第五十九节）。
    """
    njs = prefix + '.njs'
    if not os.path.isfile(njs):
        return
    try:
        dst = njs + '.lmdb-disabled'
        if os.path.exists(dst):
            os.remove(dst)
        os.replace(njs, dst)
    except OSError:
        pass


def makeblastdb_nucl(prefix, ref, title, logger=None, v4=True):
    """建核酸 BLAST 库，返回 prefix。

    v4=True → 加 `-blastdb_version 4`，**任意路径可用**。v5 的 LMDB 后端在含
    中文的路径下打不开（`mdb_env_open: 系统找不到指定的路径`），历史上只能把库
    建到 ASCII 临时目录；v4 没这个限制，所以**随包分发的预置库一律用 v4**——
    安装到哪个盘都行，也不再依赖 %TEMP%。
    """
    cmd = [get_config().tool('makeblastdb'), '-in', ref, '-dbtype', 'nucl',
           '-out', prefix, '-title', title]
    if v4:
        cmd += ['-blastdb_version', '4']
    run_cmd(cmd, logger=logger)
    if v4:
        _disable_lmdb_journal(prefix)
    return prefix


def _rebuild_if_stale(prefix, ref, title, logger=None, v4=True, tag=''):
    """`prefix` 处库可用且戳未过期 → 直接返回；否则重建。返回 prefix。"""
    sig = _blast_db_sig(ref)
    stamp = prefix + '.stamp'
    if _db_present(prefix):
        old = ''
        try:
            with open(stamp, encoding='utf-8') as sf:
                old = sf.read().strip()
        except OSError:
            pass
        if old == sig:
            return prefix
        if logger:
            logger.log(f'{tag or prefix} 的源 FASTA 已更新，重建 BLAST 库')
        _drop_db(prefix)
    if logger:
        logger.log(f'构建 {tag or "病毒"} BLAST 库'
                   f'（{"v4" if v4 else "v5"}）: {ref} -> {prefix}')
    makeblastdb_nucl(prefix, ref, title, logger=logger, v4=v4)
    try:
        with open(stamp, 'w', encoding='utf-8') as sf:
            sf.write(sig)
    except OSError:
        pass
    return prefix


def ensure_virus_blast_db(logger=None):
    """病毒参考核酸 BLAST 库前缀（`local_search._nuc_db()` 用）。

    优先级：
      ① **随包预置的 v4 库** `databases/virusref_db/blast/virus.*`
         —— 安装后零构建、不依赖 %TEMP%、也不受安装路径是否含中文影响；
      ② 预置库缺失/过期 → **就地重建 v4**（含中文的路径照样能建、能用）；
      ③ 预置目录不可写（只读介质等）→ 回退 ASCII 临时目录建 v5（历史行为）。
    """
    ref = find_virus_ref_fasta()
    prefix = os.path.join(DIRS['virus_src'], 'blast', 'virus')
    try:
        os.makedirs(os.path.dirname(prefix), exist_ok=True)
        return _rebuild_if_stale(prefix, ref, 'virus_ref', logger=logger,
                                 v4=True, tag='病毒参考(virus)')
    except Exception as e:                       # 只读介质 / 权限不足
        if logger:
            logger.log(f'预置目录建库失败（{e}），回退 ASCII 临时目录（v5）',
                       'WARN')
    return _rebuild_if_stale(os.path.join(_ascii_work_base('vp_blast'), 'virus'),
                             ref, 'virus_ref', logger=logger, v4=False,
                             tag='病毒参考(virus)·临时')

# srna 模式的奇数 k 序列（VirusDetect 同款默认 9-23）。组装时按嗅探到的
# reads 最大长度截断（k 必须 < 最大读长，SPAdes 对 k ≥ 读长直接报错）。
SRNA_KMERS = (9, 11, 13, 15, 17, 19, 21, 23)

def _sniff_max_read_len(path, max_records=4000, logger=None):
    """嗅探 reads 头部若干条的最大长度；无法判断（空文件/读失败）返回 None。

    支持 FASTQ（含 .gz）与 FASTA（含 .gz）：自动 srna 切换与 k 截断都用它
    当参考值，不追求全量精确——超短 reads（小RNA）长度集中且均匀，头部
    采样足够判断。FASTA 场景真实存在：输入可直接是 FASTA（见
    `_seq_file_kind`），≤32bp 的 FASTA reads 同样需要自动切 srna。

    读失败（截断 gz 的 EOFError、权限/占用等 OSError）或格式不可识别也
    返回 None，但都会留 WARN：返回 None 意味着"无法判断是否超短 reads"，
    静默吞掉会让 sRNA 数据走错模式（平台日志铁律：跳过必须带类型和原因）。
    """
    try:
        with open(path, 'rb') as probe:
            gz = probe.read(2) == b'\x1f\x8b'
        if gz:
            import gzip
            fh = gzip.open(path, 'rt', encoding='utf-8', errors='replace')
        else:
            fh = open(path, 'r', encoding='utf-8', errors='replace')
        with fh:
            first = ''
            for ln in fh:
                s = ln.strip()
                if s:
                    first = s
                    break
            if first.startswith('>'):
                # FASTA：累计每条记录长度（兼容多行折行），只取头部 N 条
                m, cur, n = 0, 0, 0
                for ln in fh:
                    s = ln.strip()
                    if s.startswith('>'):
                        m = max(m, cur)
                        cur = 0
                        n += 1
                        if n >= max_records:
                            break
                    else:
                        cur += len(s)
                else:
                    m = max(m, cur)
                return m or None
            if not first.startswith('@'):
                if logger:
                    logger.log(f"读长嗅探：{os.path.basename(str(path))} 既非 "
                               f"FASTQ 也非 FASTA（首行 {first[:20]!r}），"
                               f"无法判断是否超短 reads", 'WARN')
                return None
            maxlen = 0
            for i, ln in enumerate(fh):
                if i >= max_records * 4:
                    break
                if i % 4 == 0:                  # 首行'@'已消费，第2/6/10…行为序列
                    maxlen = max(maxlen, len(ln.strip()))
            return maxlen or None
    except (OSError, EOFError) as e:
        if logger:
            logger.log(f"读长嗅探失败（{type(e).__name__}: {e}），"
                       f"无法判断是否超短 reads，"
                       f"{os.path.basename(str(path))} 将按常规模式处理", 'WARN')
        return None

def _is_ascii(p):
    return all(ord(c) < 128 for c in str(p))

def _seq_file_kind(path):
    """判断序列文件的真实类型，用于给 SPAdes 中转文件起正确扩展名。

    返回 (is_gzip, kind)，kind ∈ {'fastq', 'fasta', 'other'}。
    - is_gzip：由魔数（$\x1f\x8b）或 .gz 扩展名判断。
    - kind：优先读前几字节内容嗅探（'@'开头+质量行 → fastq；'>'开头 → fasta），
      嗅探不到时回退到扩展名。
    之所以必须做：旧版把任意输入强制复制成 in_R1.fastq.gz，当输入是 FASTA
    时 SPAdes 会用 gzip 解压 FASTA 文本 → BadGzipFile（用户实报的错误）。
    """
    is_gz = str(path).lower().endswith('.gz')
    try:
        with open(path, 'rb') as fh:
            magic = fh.read(2)
        if magic == b'\x1f\x8b':
            is_gz = True
    except OSError:
        pass

    kind = 'other'
    low = str(path).lower()
    # 先按扩展名给出候选
    if is_gz:
        base = low[:-3]
    else:
        base = low
    if base.endswith(('.fastq', '.fq')):
        kind = 'fastq'
    elif base.endswith(('.fasta', '.fa', '.fna', '.fas')):
        kind = 'fasta'

    # 内容嗅探：扩展名可能被改错（用户把 FASTA 存成 .fastq.gz），
    # 所以 gz 也要嗅探——否则扩展名直接决定 kind，SPAdes 会因
    # "reads should be in FASTQ format" 崩溃。
    try:
        if is_gz:
            import gzip as _gz
            with _gz.open(path, 'rt', encoding='utf-8', errors='replace') as fh:
                first = ''
                for ln in fh:
                    s = ln.strip()
                    if s:
                        first = s
                        break
        else:
            with open(path, 'r', encoding='utf-8', errors='replace') as fh:
                first = ''
                for ln in fh:
                    s = ln.strip()
                    if s:
                        first = s
                        break
        if first.startswith('>'):
            kind = 'fasta'
        elif first.startswith('@'):
            kind = 'fastq'
    except OSError:
        pass
    return is_gz, kind

def _shim_name(src, role):
    """给中转输入文件构成正确的目标名，保留真实格式与压缩状态。

    role: 'R1' / 'R2'。生成的扩展名与源文件实际格式一致，避免
    FASTA 被误命名为 .fastq.gz 而触发 SPAdes 的 gzip 解析崩溃。
    """
    is_gz, kind = _seq_file_kind(src)
    if kind == 'fastq':
        ext = '.fastq.gz' if is_gz else '.fastq'
    elif kind == 'fasta':
        ext = '.fasta.gz' if is_gz else '.fasta'
    else:
        # 未知类型：保留源扩展名，仅在确为 gz 时补 .gz
        low = str(src).lower()
        if is_gz and not low.endswith('.gz'):
            ext = low[low.rfind('.'):] + '.gz' if '.' in low else '.gz'
        else:
            ext = os.path.splitext(src)[1] or '.txt'
    return f'in_{role}{ext}'

def _spades_tmp_base():
    """SPAdes Windows 版不支持非 ASCII 路径：返回纯 ASCII 临时工作根目录。

    仅允许 %TEMP%（通常为 C:\\Users\\<user>\\AppData\\Local\\Temp）；
    若 TEMP 含非 ASCII 则回退系统盘根下的固定目录（原 C:\\vp_spades_tmp
    硬编码 C 盘，装在其它盘的机器上可能没有 C 盘写权限）。
    """
    import tempfile
    cand = os.path.abspath(tempfile.gettempdir())
    if os.path.isdir(cand) and _is_ascii(cand):
        return os.path.join(cand, 'vp_spades')
    drv = (os.environ.get('SYSTEMDRIVE') or 'C:').rstrip('\\') or 'C:'
    return os.path.join(drv + os.sep, 'vp_spades_tmp')

# SPAdes 各版本的 k 阶段行格式（实测 SPAdes 3.15 的 spades.log 为第一、
# 三种；旧版/摘要行为第二、三种）：
#   Assembling dataset ("dataset.info") with K=49
#   Assembling K=21
#   ===== K9 started. =====
# 旧实现只认第二种 → 主格式下 seen 永远为空，进度条卡在 5%（实测日志验证）。
_K_PROGRESS_RE = re.compile(
    r'with\s+K=(\d+)|Assembling\s+(?:k\s*=\s*)?(\d+)|K(\d+)\s+started',
    re.IGNORECASE)


def _spades_monitor(spades_dir, mode, progress):
    """SPAdes 看门狗：解析 spades.log 的 k 阶段行数估算进度。

    返回 monitor_fn（传给 run_cmd）；progress 为 None 时返回 None。
    """
    if progress is None:
        return None
    k_total = {'rnaviral': 3, 'rna': 3, 'metaviral': 4, 'isolate': 4,
               'meta': 5, 'srna': len(SRNA_KMERS)}.get(mode, 4)
    log_p = os.path.join(spades_dir, 'spades.log')
    seen = set()
    # 增量读：只解析新写入的字节。日志里 k 行出现在每个 k 阶段的开头，
    # 只读尾部 8KB 会整段错过（实测 110KB 日志的 K= 行都在前 380 行）。
    state = {'pos': 0, 'tail': ''}

    def _watch():
        try:
            if not os.path.isfile(log_p):
                return
            with open(log_p, 'rb') as f:
                if os.path.getsize(log_p) < state['pos']:
                    state['pos'] = 0        # 日志被清空重写（重试组装）
                f.seek(state['pos'])
                chunk = f.read()
                state['pos'] = f.tell()
            if not chunk:
                return
            text = state['tail'] + chunk.decode('utf-8', errors='replace')
            cut = text.rfind('\n') + 1
            state['tail'] = text[cut:]      # 末行可能写了一半，留待下次拼接
            for m in _K_PROGRESS_RE.finditer(text[:cut]):
                k = next((g for g in m.groups() if g), None)
                if k:
                    seen.add(k)
            k_done = min(len(seen), k_total)
            if k_done:
                progress(min(0.05 + k_done / k_total * 0.85, 0.92),
                         f'SPAdes 组装中：k 阶段 {k_done}/{k_total}')
        except OSError:
            pass
    return _watch

def run_spades(r1, r2, out_dir, mode='rnaviral', threads=None, memory_gb=64,
               logger=None, progress=None, info=None, extra_se=None):
    """运行 SPAdes。返回 contigs fasta 路径（rna 模式为 transcripts.fasta）。

    info：可选 dict，回填实际执行情况（requested_mode / mode / auto_srna /
    dropped_r2 / fallback_rna）——自动切换模式、丢弃 R2 这类会显著影响
    结果解释的行为必须进 summary.json，不能只留在日志里。

    srna（小RNA）模式：--only-assembler + 奇数 k 序列 + 单端 -s（R2 忽略），
    超出读长的 k 由 SPAdes 自动裁剪，故 k 序列全量直传。常规模式遇超短
    reads（最长 <32bp，小RNA/miRNA 建库的典型特征）自动切换 srna——这些
    reads 进不了默认 k=21+ 的图，实测 SPAdes 自动降 k 后还会在 rna 模式
    配置下崩溃。

    rnaviral/metaviral/meta 模式未拼出 contigs 时（低深度样品常见——病毒模式
    的"染色体外元件"筛选门槛较严），自动清空输出并用 rna 模式重试一次；
    srna 模式本身已是最后手段，失败不再重试。
    平台路径含中文时自动经纯 ASCII 临时目录中转（输入复制过去、结果拷回）。
    """
    max_read_len = None
    for _p in (r1, r2):
        if _p:
            _n = _sniff_max_read_len(_p, logger=logger)
            if _n and (max_read_len is None or _n > max_read_len):
                max_read_len = _n
    if info is not None:
        info['requested_mode'] = mode
        info['max_read_len_sampled'] = max_read_len
    auto_srna = False
    if max_read_len and max_read_len < 32 and mode != 'srna':
        auto_srna = True
        mode = 'srna'
        if info is not None:
            info['auto_srna'] = True
        if logger:
            logger.log(f"检测到超短 reads（头部采样最长 {max_read_len}bp，"
                       f"疑似小RNA/siRNA 数据），自动改用 srna 模式"
                       f"（--only-assembler + 奇数 k 序列 + 单端组装）", "WARN")
    try:
        return _run_spades_once(r1, r2, out_dir, mode, threads, memory_gb,
                                logger, progress, max_read_len=max_read_len,
                                info=info, extra_se=extra_se)
    except RuntimeError as e:
        if (not auto_srna and mode in ('rnaviral', 'metaviral', 'meta')
                and '未找到 contigs' in str(e)):
            if logger:
                logger.log(f"{mode} 模式未拼出 contigs（该模式对低深度病毒样品"
                           f"筛选较严），自动改用 rna 模式重试组装 ...", "WARN")
            if info is not None:
                info['fallback_rna'] = True
            spades_dir = check_path(os.path.join(out_dir, 'spades'),
                                    must_exist=False, in_platform=True)
            shutil.rmtree(spades_dir, ignore_errors=True)
            return _run_spades_once(r1, r2, out_dir, 'rna', threads,
                                    memory_gb, logger, progress,
                                    max_read_len=max_read_len, info=info,
                                    extra_se=extra_se)
        raise

def _run_spades_once(r1, r2, out_dir, mode, threads, memory_gb, logger,
                     progress=None, max_read_len=None, info=None,
                     extra_se=None):
    import uuid
    cfg = get_config()
    spades = cfg.tool('spades')
    threads = threads or cfg.threads
    # 内存钳制（2026-10-04 审查 E3）：-m 64 硬编码在 16GB 机器上换页抖动；
    # POST 侧已钳（web/samples.py），这里覆盖 CLI 直调路径
    from .utils import total_ram as _total_ram
    _ram_cap = max(2, int((_total_ram() or 16 * 1024 ** 3)
                          * 0.6 / 1024 ** 3))
    if int(memory_gb) > _ram_cap:
        if logger:
            logger.log(f"组装内存 {memory_gb}GB 超物理内存 60% 上限"
                       f"（{_ram_cap}GB），已自动钳制", "WARN")
        memory_gb = _ram_cap
    if mode == 'srna' and r2:
        if logger:
            logger.log("srna 模式按单端 -s 组装（VirusDetect 同款策略），忽略 R2",
                       "WARN")
        r2 = None
        if info is not None:
            info['dropped_r2'] = True
    if info is not None:
        info['mode'] = mode
    spades_dir = check_path(os.path.join(out_dir, 'spades'),
                            must_exist=False, in_platform=True)
    os.makedirs(spades_dir, exist_ok=True)

    if extra_se and not (os.path.isfile(extra_se)
                         and os.path.getsize(extra_se) > 0):
        extra_se = None
    use_shim = not (_is_ascii(spades_dir) and _is_ascii(r1)
                    and (not r2 or _is_ascii(r2)))
    spades_in1, spades_in2, spades_out = r1, r2, spades_dir
    shim_dir = None
    if use_shim:
        shim_dir = os.path.join(_spades_tmp_base(), uuid.uuid4().hex[:10])
        os.makedirs(shim_dir, exist_ok=True)
        spades_out = os.path.join(shim_dir, 'spades')
        spades_in1 = spades_in2 = None
        if logger:
            logger.log("平台路径含中文，SPAdes 经 ASCII 临时目录中转: " + shim_dir)
        _g1, _k1 = _seq_file_kind(r1)
        if _k1 == 'fasta' and logger:
            logger.log("输入为 FASTA（无质量值），SPAdes 将跳过纠错"
                       "仅做组装，建议改用 FASTQ 数据以获得更好结果")
        spades_in1 = os.path.join(shim_dir, _shim_name(r1, 'R1'))
        shutil.copyfile(check_path(r1, must_exist=True), spades_in1)
        if r2:
            spades_in2 = os.path.join(shim_dir, _shim_name(r2, 'R2'))
            shutil.copyfile(check_path(r2, must_exist=True), spades_in2)
        if extra_se:
            spades_inse = os.path.join(shim_dir, _shim_name(extra_se, 'SE'))
            shutil.copyfile(check_path(extra_se, must_exist=True), spades_inse)
            extra_se = spades_inse

    # FASTA reads 无质量值：SPAdes 默认纠错只认 FASTQ，会直接报错
    # （"to run read error correction, reads should be in FASTQ format"）。
    # 输入为 FASTA 时自动加 --only-assembler 跳过纠错，仅做组装。
    only_assembler = False
    try:
        _, _k1 = _seq_file_kind(r1)
        _k2 = None
        if r2:
            _, _k2 = _seq_file_kind(r2)
        if _k1 == 'fasta' or _k2 == 'fasta':
            only_assembler = True
    except Exception:
        only_assembler = False

    try:
        if spades_in2:
            in_args = ['-1', spades_in1, '-2', spades_in2]
        else:
            in_args = ['-s', spades_in1]          # 单端数据
        if extra_se:
            # ②b 抽取的孤儿单端 reads（只有一条 mate 比中病毒）：
            # SPAdes 支持重复 -s（各自成为独立单端库）
            in_args += ['-s', extra_se]
        if spades_in2 is None:
            if mode in ('metaviral', 'meta'):
                # SPAdes meta/metaviral 实际不使用单端 reads
                # （警告 "Single reads are not used in metagenomic mode"），
                # 单端自动降级 rna 模式（病毒转录组/低覆盖组装）
                if logger:
                    logger.log(f"单端数据：{mode} 模式不支持单端 reads，"
                               f"自动改用 rna 模式")
                mode = 'rna'
        if mode == 'srna':
            # 不按嗅探长度截 k：头部采样会低估（sRNA fastq 常按长度排序），
            # 而 SPAdes 对超出读长的 k 会自动裁剪（"K-mer sizes were set to
            # [19] because estimated read length (21) is less than 23"），
            # 全序列直传既安全又不丢高 k 阶段。
            if max_read_len and max_read_len <= SRNA_KMERS[0]:
                raise RuntimeError(
                    f"reads 过短（最长 {max_read_len}bp），srna 模式无法组装"
                    f"（最小 k={SRNA_KMERS[0]}，要求 k < 读长）")

        # 组装尝试序列：k-mer 覆盖模型（kmer_coverage_model）对偏斜数据/小子集
        # 可能报 "Invalid kmer coverage histogram"（退出码 64）——所有模式都
        # 按阶梯降级重试：srna 逐次去掉最高两个 k；rnaviral/rna 退到显式
        # -k 21,33 → 21；metaviral/meta/isolate 同理退到默认短梯。
        if mode == 'srna':
            ladder = list(SRNA_KMERS)
            attempt_extras = []
            while len(ladder) >= 3:
                attempt_extras.append(
                    ['--only-assembler', '-k', ','.join(str(k) for k in ladder)])
                ladder = ladder[:-2]
        else:
            base = ['--' + mode]
            if only_assembler:
                base.append('--only-assembler')
            attempt_extras = [list(base)]
            if mode in ('rnaviral', 'rna', 'metaviral', 'meta', 'isolate'):
                attempt_extras.append(base + ['-k', '21,33'])
                attempt_extras.append(base + ['-k', '21'])

        last_err = None
        for extra in attempt_extras:
            k_str = next((extra[i + 1] for i, x in enumerate(extra)
                          if x == '-k'), '默认')
            if info is not None and '-k' in extra:
                info['k'] = k_str
            cmd = ([spades] + extra + ['-t', str(threads), '-m', str(int(memory_gb)),
                                       '-o', spades_out] + in_args)
            if logger:
                extra_desc = f"k={k_str}" if '-k' in extra else "默认 k 阶梯"
                only_desc = "，仅组装不纠错" if (
                    '--only-assembler' in extra or only_assembler and mode == 'srna') else ""
                se_desc = "，单端 -s" if mode == 'srna' else ""
                logger.log(f"SPAdes 组装 (mode={mode}, {extra_desc}{se_desc}, "
                           f"threads={threads}, mem={memory_gb}GB{only_desc})")
            mon = _spades_monitor(spades_out, mode, progress)
            try:
                run_cmd(cmd, logger=logger, monitor_fn=mon, monitor_interval=15)
                last_err = None
                break
            except RuntimeError as e:
                last_err = e
                if '(退出码 64)' not in str(e) or extra is attempt_extras[-1]:
                    raise
                if logger:
                    logger.log("k-mer 覆盖模型失败（退出码 64），"
                               "降级 k 阶梯重试", "WARN")
                shutil.rmtree(spades_out, ignore_errors=True)
                os.makedirs(spades_out, exist_ok=True)

        result_dir = spades_out
        if use_shim:
            # 拷回关键产物（contigs/图/日志），丢弃体积大的中间 K* 目录
            keep_names = ['contigs.fasta', 'scaffolds.fasta', 'transcripts.fasta',
                          'hard_filtered_transcripts.fasta',
                          'soft_filtered_transcripts.fasta',
                          'assembly_graph_with_scaffolds.gfa', 'assembly_graph_after_simplification.gfa',
                          'spades.log', 'warnings.log', 'params.txt', 'input_dataset.yaml',
                          'before_rr.fasta', 'dataset.info', 'run_spades.sh', 'run_spades.yaml']
            moved = 0
            for name in keep_names:
                src = os.path.join(spades_out, name)
                if os.path.isfile(src):
                    shutil.copyfile(src, os.path.join(spades_dir, name))
                    moved += 1
            if logger:
                logger.log(f"SPAdes 产物拷回 {moved} 个文件 -> {spades_dir}")
            result_dir = spades_dir

        # 成功后清理体积大的 K* 中间目录（此前仅 shim 分支丢弃，非 shim
        # 全留在 03_assembly/spades/ 下可达数十 GB——2026-10-04 审查 E3）
        if not use_shim:
            import glob as _glob
            for _kd in _glob.glob(os.path.join(spades_dir, 'K*')):
                shutil.rmtree(_kd, ignore_errors=True)

        # 产物名随 SPAdes 版本/模式而异：常规/rnaviral/metaviral/meta 出
        # contigs.fasta，旧版 rnaSPAdes 出 transcripts.fasta，新版只出
        # hard_filtered_transcripts.fasta（主结果）/ soft_filtered_transcripts.fasta。
        candidates = ['contigs.fasta', 'transcripts.fasta',
                      'hard_filtered_transcripts.fasta',
                      'soft_filtered_transcripts.fasta', 'scaffolds.fasta']
        for name in candidates:
            p = os.path.join(result_dir, name)
            if os.path.isfile(p) and os.path.getsize(p) > 0:
                return check_path(p, must_exist=True)
        raise RuntimeError(f"SPAdes 结束但未找到 contigs 输出: {result_dir}")
    finally:
        if use_shim and shim_dir and os.path.isdir(shim_dir):
            shutil.rmtree(shim_dir, ignore_errors=True)

def filter_contigs(contigs_fasta, out_fasta, min_len=200, logger=None):
    """按长度过滤 contigs，返回 (过滤文件, 条数, 总长)。"""
    n, total_bp = 0, 0
    with safe_open(out_fasta, 'wt') as f:
        for header, seq in iter_fasta(contigs_fasta):
            if len(seq) >= min_len:
                cid = header.split()[0]
                write_fasta_record(f, cid, seq)
                n += 1
                total_bp += len(seq)
    if logger:
        logger.log(f"contig 过滤(≥{min_len}bp): 保留 {n} 条, 总长 {total_bp:,}bp")
    return out_fasta, n, total_bp

def _assembly_reuse_ok(fp_file, fp, logger=None):
    """断点复用判据：已有 SPAdes 产物对应的输入指纹 + 关键参数是否一致。

    指纹 = 输入文件（path+size+mtime_ns，对齐 preprocess.input_fingerprint）
    + 关键参数（assembly_mode / min_contig_len / assembly_input）。旧产物
    没有指纹记录、或记录与本次不一致 → 不可复用（重跑组装）并留日志：
    换输入/换模式/换长度阈值后沿用旧 contigs 会让"续跑"与"全量重跑"
    得到两套不同结果，且从产物上看不出来。
    """
    import json
    try:
        with safe_open(fp_file) as f:
            old = json.load(f)
    except (OSError, ValueError):
        if logger:
            logger.log("已有 SPAdes 产物但无输入指纹记录（旧版本产物），"
                       "为安全起见重新组装", "WARN")
        return False
    if old == fp:
        return True
    if logger:
        logger.log("已有 SPAdes 产物的输入指纹/关键参数与本次不一致，"
                   "重新组装（不复用旧 contigs）", "WARN")
    return False


def assemble_and_classify(sample_dir, r1, r2, db_virus, mode='rnaviral',
                          threads=None, memory_gb=64, min_contig_len=200,
                          logger=None, force=False, chunk_dir=None,
                          progress=None, assembly_input=None):
    step = 'assembly'
    out_dir = check_path(os.path.join(sample_dir, '03_assembly'),
                         must_exist=False, in_platform=True)
    os.makedirs(out_dir, exist_ok=True)
    summary_file = os.path.join(out_dir, 'summary.json')
    if is_step_done(out_dir, step) and not force:
        # .done 复用必须校验指纹（审查 P1-1）：下方 713-734 行已为 SPAdes
        # 产物实现了输入/参数指纹校验，但旧版 .done 快速返回根本走不到那
        # 里——换输入或改 min_contig_len/assembly_mode 后仍整段沿用旧结果。
        # 与 _fp 同口径：inputs + {assembly_mode, min_contig_len,
        # assembly_input}；旧版 summary 无 inputs_fp 键 → 不匹配 → 重跑
        # （安全方向：宁可重跑也不沿用来路不明的产物）。
        import json
        try:
            from .preprocess import input_fingerprint as _ifp_asm
            _cur_fp = {
                'inputs': [_ifp_asm(p) for p in (r1, r2) if p],
                'params': {'assembly_mode': mode,
                           'min_contig_len': min_contig_len,
                           'assembly_input': assembly_input},
            }
            with safe_open(summary_file) as f:
                _prev = json.load(f)
            if _prev.get('inputs_fp') == _cur_fp:
                if logger:
                    logger.log("阶段③组装已完成（输入与参数未变化），跳过")
                return _prev
            if logger:
                logger.log('阶段③存在断点，但输入 reads 或组装参数已变化'
                           '——重跑本阶段（不复用旧结果）', 'WARN')
        except (OSError, ValueError) as _e:
            if logger:
                logger.log(f'阶段③断点摘要不可读（{_e}），重跑本阶段', 'WARN')

    # 1. 组装（若上次运行已有 contigs 且未强制重跑则复用——断点续跑）
    # 复用判定要覆盖 SPAdes 全部可能产物名：rna 模式从不生成
    # contigs.fasta（只有 transcripts*/hard_filtered_transcripts.fasta），
    # 原判定只认 contigs.fasta → rna 模式每次都重跑数小时。
    _spades_dir = os.path.join(out_dir, 'spades')
    _fp_file = os.path.join(_spades_dir, '.input_fingerprint.json')
    # 本次组装输入指纹（复用判据，见 _assembly_reuse_ok；写/比都在这一个文件）
    from .preprocess import input_fingerprint
    _fp = {
        'inputs': [input_fingerprint(p) for p in (r1, r2) if p],
        'params': {'assembly_mode': mode, 'min_contig_len': min_contig_len,
                   'assembly_input': assembly_input},
    }
    _cands = ['contigs.fasta', 'transcripts.fasta',
              'hard_filtered_transcripts.fasta',
              'soft_filtered_transcripts.fasta', 'scaffolds.fasta']
    contigs_raw = None
    if not force:
        for _n in _cands:
            _p = os.path.join(_spades_dir, _n)
            if os.path.isfile(_p) and os.path.getsize(_p) > 0:
                contigs_raw = _p
                break
        if contigs_raw is not None and not _assembly_reuse_ok(
                _fp_file, _fp, logger=logger):
            contigs_raw = None            # 指纹/参数不匹配：不复用，重跑组装
    sp_info = {}
    if contigs_raw is None:
        contigs_raw = run_spades(r1, r2, out_dir, mode=mode, threads=threads,
                                 memory_gb=memory_gb, logger=logger,
                                 progress=(lambda p, m: progress(
                                     p * 0.55, f'③ SPAdes：{m}'))
                                 if progress else None, info=sp_info)
        # 记下本次组装输入的指纹与关键参数：下次断点复用前必须核对一致
        try:
            import json
            os.makedirs(_spades_dir, exist_ok=True)
            with safe_open(_fp_file, 'wt') as f:
                json.dump(_fp, f, ensure_ascii=False, indent=1)
        except OSError as e:
            if logger:
                logger.log(f'SPAdes 输入指纹落盘失败（不影响本次结果）: {e}',
                           'WARN')
    else:
        if logger:
            logger.log("复用已有 SPAdes 组装结果（输入指纹与关键参数一致，"
                       "跳过组装）")
        contigs_raw = check_path(contigs_raw, must_exist=True)
    if progress:
        progress(0.58, '③ contig 长度过滤')
    # 2. 过滤
    contigs_fa = os.path.join(out_dir, 'contigs.filtered.fasta')
    contigs_fa, n_contigs, total_bp = filter_contigs(
        contigs_raw, contigs_fa, min_len=min_contig_len, logger=logger)
    if n_contigs == 0:
        n_raw = count_fasta_seqs(contigs_raw)
        _hint = (f"SPAdes 产物 {os.path.basename(str(contigs_raw))} 共 {n_raw} "
                 f"条序列，无一条达到长度阈值 {min_contig_len}bp")
        _ml = _sniff_max_read_len(r1, logger=logger) if r1 else None
        if _ml and _ml < 32:
            _hint += (f"；输入 reads 超短（最长 {_ml}bp，小RNA/siRNA 数据），"
                      f"建议把「最小 contig 长度」调到 ≤100bp 并确认已自动"
                      f"切换 srna 模式")
        elif n_raw:
            _hint += "，建议调低「最小 contig 长度」或检查数据质量/深度"
        raise RuntimeError(f"组装后无满足长度要求的 contigs，无法继续（{_hint}）")

    # 3. kunpeng 分类 contigs
    if logger:
        logger.log("kunpeng 病毒库分类 contigs ...")
    if progress:
        progress(0.62, '③ kunpeng 分类 contigs')
    cls_out = os.path.join(out_dir, 'contig_classify')
    res = classify(db_virus, [contigs_fa], cls_out, paired=False,
                   threads=threads, logger=logger, chunk_dir=chunk_dir,
                   progress=(lambda p, m: progress(
                       0.62 + p * 0.13, f'③ 分类：{m}'))
                   if progress else None)
    contig_tax = {}
    if res['kraken']:
        for flag, cid, taxid, _len, _path in parse_classify_output(res['kraken']):
            contig_tax[cid] = (flag, taxid)
    # 轻量分类树（来自 contig kreport），供名称/谱系查询
    tree = KreportTree(parse_kreport(res['kreport'])) if res['kreport'] else KreportTree([])

    # 4. 汇总（kunpeng 分类 + contig 级 BLASTN 最近参考）
    #
    # 两件事在这里一起做：
    #  a) 落盘 viral_contigs.fasta —— 全平台有 9 个消费方按这个名字取（④宿主
    #     拆分、⑦/⑨出图、⑧引物、logan、suvtk），但此前只有 ⑥ORF 的
    #     extract_viral_contigs() 会写；默认阶段序里 hostana(④) 在 orf(⑥)
    #     之前，所以首次全量运行 ④ 拿不到它 → 08_host_analysis/*.classified.fasta
    #     被静默跳过（实测 ERR7586041 该产物 0 个）。
    #  b) 填 blast_* 六列 —— 此前只有列名没有写入方，host_analysis 的
    #     「BLAST 最后回退」从未生效（实测 212 条 host_prediction 中
    #     class_source=blast 为 0，14% 落到 Unknown）。
    viral_tsv = os.path.join(out_dir, 'virus_contigs.tsv')
    viral_fa = os.path.join(out_dir, 'viral_contigs.fasta')
    # 未分类候选池（审查 P1-3）：kunpeng 非 C / taxid≤0 的 contig 此前在
    # 生成 viral_contigs.fasta 时被直接丢弃，从未进入 blastx/BFVD 验证——
    # 对参考分类库差异大的新病毒会系统性漏检。这里把它们按长度合格为
    # 条件保留成独立文件，供 ③b 的独立验证通道消费（不进下游注释）。
    unclass_fa = os.path.join(out_dir, 'unclassified_contigs.fasta')
    viral_contigs = []
    n_unclass_pool = 0
    n_unclass_short = 0
    with safe_open(viral_fa, 'wt') as vf, safe_open(unclass_fa, 'wt') as uf:
        for header, seq in iter_fasta(contigs_fa):
            cid = header.split()[0]
            flag, taxid = contig_tax.get(cid, ('-', 0))
            if flag == 'C' and taxid > 0:
                viral_contigs.append(cid)
                write_fasta_record(vf, cid, seq)
            elif len(seq) >= min_contig_len:
                n_unclass_pool += 1
                write_fasta_record(uf, cid, seq)
            else:
                n_unclass_short += 1
    if progress:
        progress(0.755, '③ contig 级 BLASTN 最近参考')
    blast_tab, _blast_tsv = {}, None
    if viral_contigs:
        try:
            from .local_search import contig_blast_table
            blast_tab, _blast_tsv = contig_blast_table(
                viral_fa, out_dir, threads=threads, logger=logger)
        except Exception as e:              # 可选增强，失败不阻断组装
            if logger:
                logger.log(f"contig 级 BLASTN 跳过: {e}", "WARN")
    _BLAST_COLS = ('blast_top_hit', 'blast_identity(%)',
                   'blast_coverage_hsp(%)', 'blast_aln_len',
                   'blast_species', 'blast_family')
    with safe_open(viral_tsv, 'wt') as f:
        f.write('contig\tlength\tkunpeng_flag\tkunpeng_taxid\tkunpeng_species\t'
                + '\t'.join(_BLAST_COLS) + '\n')
        for header, seq in iter_fasta(contigs_fa):
            cid = header.split()[0]
            flag, taxid = contig_tax.get(cid, ('-', 0))
            if not (flag == 'C' and taxid > 0):
                continue
            species_k = tree.name(taxid) if taxid else ''
            b = blast_tab.get(cid) or {}
            f.write(f"{cid}\t{len(seq)}\t{flag}\t{taxid}\t{species_k}\t"
                    f"{b.get('accession', '')}\t{b.get('identity', '')}\t"
                    f"{b.get('coverage_hsp', '')}\t{b.get('aln_len', '')}\t"
                    f"{b.get('species', '')}\t{b.get('family', '')}\n")

    # 5. 完整分类谱系表（8 级 + 属长比）：03b_verify 的宿主归属来源。
    #    口径与工具④（组装结果再鉴定）一致——两者共用 classify_rows。
    _cls_warn = ''
    try:
        _write_classification_table(out_dir, res.get('kraken'), logger)
        if not os.path.isfile(os.path.join(out_dir,
                                          'virus_classification.tsv')):
            _cls_warn = '分类表未生成（无 kunpeng 分类输出）'
    except Exception as e:
        # 谱系表缺失不应阻断组装主流程（verify 会降级为无宿主归属），
        # 但必须让卡片/报告看得到降级原因，不能只留在日志里。
        _cls_warn = f'{type(e).__name__}: {e}'
    if _cls_warn:
        if logger:
            logger.log(f"警告：virus_classification.tsv 未生成：{_cls_warn}",
                       "WARN")

    summary = {
        'stage': step,
        'assembly_mode': mode,
        # 实际执行的 SPAdes 情况：requested_mode/mode/auto_srna/dropped_r2/
        # fallback_rna/max_read_len_sampled（含自动降级时 mode 与请求不同）
        'spades': (sp_info or None),
        'raw_contigs': count_fasta_seqs(contigs_raw),
        'contigs': n_contigs, 'total_bp': total_bp,
        'contigs_fasta': 'contigs.filtered.fasta',
        'viral_contigs': viral_contigs,
        'virus_contigs_tsv': 'virus_contigs.tsv',
        'classification_table': (None if _cls_warn else
                                 'virus_classification.tsv'),
        'classification_warning': _cls_warn,
        'kreport': os.path.basename(res['kreport']) if res['kreport'] else None,
        # 未分类候选池计数（P1-3）：进入池的数量 / 因长度不足被排除的数量。
        # 报告口径：viral_contigs = 分类器判定候选；unclassified_pool = 分类
        # 未定但长度合格的独立候选（03b 有独立验证通道）；short = 长度不足
        # 被排除，原因在此明示而不是静默丢弃。
        'unclassified_pool_n': n_unclass_pool,
        'unclassified_pool_fa': 'unclassified_contigs.fasta',
        'unclassified_excluded_short_n': n_unclass_short,
        # 复用判据（审查 P1-1）：.done 快速返回时据此校验——输入 reads /
        # 组装模式 / 长度阈值任一变化都必须重跑，不再只看 .done 存在。
        'inputs_fp': _fp,
    }
    import json
    with safe_open(summary_file, 'wt') as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)
    mark_step_done(out_dir, step)
    if logger:
        logger.log(f"阶段③ 完成: contigs {n_contigs} 条，病毒 contigs {len(viral_contigs)} 条")
    return summary


def _write_classification_table(out_dir, kraken_txt, logger=None):
    """kunpeng 原始输出 → virus_classification.tsv（8 级谱系 + 属长比）。

    与工具④（组装结果再鉴定）同口径，共用 contig_annot.classify_rows。
    下游 03b_verify 从本表取宿主归属（family → ICTV 宿主）。
    genus_avg_map 首跑建缓存，失败时退化为空 map（ratio/near_complete 为 0）。
    """
    from .contig_annot import classify_rows, genus_avg_map, RANKS
    tsv_path = os.path.join(out_dir, 'virus_classification.tsv')
    if not kraken_txt or not os.path.isfile(kraken_txt):
        if logger:
            logger.log("无 contig 分类输出，跳过谱系表")
        return None
    try:
        genus_map = genus_avg_map(logger=logger)
    except Exception as e:
        if logger:
            logger.log(f"属平均长度缓存不可用（{e}），ratio 置 0")
        genus_map = {}
    rows = classify_rows(kraken_txt, genus_map)
    header = (['contig', 'taxid', 'taxon'] + list(RANKS)
              + ['length', 'genus_avg_len', 'ratio', 'near_complete',
                 'score', 'kmer_support', 'kmer_total'])
    with safe_open(tsv_path, 'wt') as f:
        f.write('\t'.join(header) + '\n')
        for r in rows:
            f.write('\t'.join(str(r.get(k, '')) for k in header) + '\n')
    if logger:
        logger.log(f"谱系表 virus_classification.tsv: {len(rows)} 条")
    return tsv_path


def _load_virus_info_species():
    """病毒 info 表 -> {accession: {species, family, genus}}（存在才读）。"""
    for p in (ref_annotation_path(),):
        if os.path.isfile(p):
            out = {}
            with safe_open(p) as f:
                for row in csv.DictReader(f, delimiter='\t'):
                    acc = (row.get('Accession') or '').strip()
                    if acc:
                        out[acc] = {
                            'species': (row.get('Species_ICTV') or row.get('Species_NCBI') or '').strip(),
                            'family': (row.get('VMR_Family') or '').strip(),
                            'genus': (row.get('VMR_Genus') or '').strip(),
                        }
            return out
    return {}
