# -*- coding: utf-8 -*-
"""②b 已知病毒识别与定量（kvsuite）——管道阶段的封装。

把 reads 比对到**已知病毒参考索引**做识别与定量（engine 二选一：
salmon=EM 定量默认 / minibwa=真比对计数，见 kv_engines；salmon
--writeBam 与 minibwa 真比对都产出映射 BAM），产出「识别 + 定量」
（`identify/` + `filter/`），并用 samtools 把比对上的病毒 reads 抽成
`viral_R1/R2.fastq.gz` 供 ③组装（`assembly_input='virus'`）。
共识段（③c）内部固定 minibwa，与此处的定量引擎无关。

产物（<sample_dir>/02b_kvsuite/）：
    summary.json           本阶段摘要（卡片状态与 GUI 用）
    identify/              识别结果（all_viruses*.tsv）
    filter/filtered.tsv    二次过滤后的结果表（定量口径）
    align/*.sorted.bam     salmon 映射 BAM（含 .bai 与 pandepth 深度）
    viral_R1/R2.fastq.gz   比对上的病毒 reads（PE）；SE 只有 viral_R1
    kvsuite.log            引擎原始日志

为什么只跑 identify + filter：consensus / 变异注释属于 ③c 共识序列与变异
阶段（`consensus`），不在 ②b 的"识别与定量"范围内。
"""
import json
import os
import re
import subprocess
import sys
import time

from .config import DIRS, PLATFORM_ROOT, get_config
from .utils import check_path, run_cmd, run_cmd_redirect

STAGE_DIR = '02b_kvsuite'
SUMMARY = 'summary.json'
VIRAL_R1 = 'viral_R1.fastq.gz'
VIRAL_R2 = 'viral_R2.fastq.gz'
VIRAL_SE = 'viral_se.fastq.gz'
_ENGINE_SCRIPT = 'known_virus_suite.py'


def engine_cmd():
    """引擎入口：(可执行, 前置参数)。

    源码模式不能走 `python -m 包.模块`：内置绿色版 Python（3rd/python，
    带 python312._pth）不会把 CWD 加进 sys.path，模块必然找不到
    （任务日志表现为 "No module named 'Virus_Platform_Core'"）。
    改用 -c 引导：先注入平台根再进 main，argparse 吃后续参数没有问题。
    冻结分发（exe）走 `--run-engine <脚本名>`（由 engine_entry 进程内执行）。
    """
    if getattr(sys, 'frozen', False):
        return sys.executable, ['--run-engine', _ENGINE_SCRIPT]
    boot = (f"import sys; sys.path.insert(0, {PLATFORM_ROOT!r}); "
            "from Virus_Platform_Core.known_virus_suite.known_virus_suite "
            "import main; main()")
    return sys.executable, ['-c', boot]


def _kv_lib_dir():
    """默认鉴定库目录（参考病毒目录 databases/virusref_db/ 下的 plant_ref）。

    2026-09-29：默认库由 kv_index 更名 plant_ref（自包含布局不变：
    reference.fasta + reference.ref_info.tsv + manifest.json，引擎索引是
    库内派生子目录，首跑自建）。
    """
    return os.path.join(DIRS['virus_src'], 'plant_ref')


def default_reference():
    """平台默认病毒参考。

    参考病毒目录重构后每个库自包含：<库>/reference.fasta +
    <库>/reference.ref_info.tsv（引擎索引为派生物：可随库分发，也可在
    首次运行时按需自建，见 index_dir_for）。解析顺序：
    ① 默认库目录内的 reference.fasta；
    ② 库 manifest.json 登记的 reference（外部来源参考）；
    ③ 旧布局根目录散置的 final.cluster.ref.fasta（未迁移安装的兼容回退）。
    """
    p = os.path.join(_kv_lib_dir(), 'reference.fasta')
    if os.path.isfile(p):
        return p
    try:
        with open(os.path.join(_kv_lib_dir(), 'manifest.json'),
                  encoding='utf-8') as fh:
            ref = (json.load(fh) or {}).get('reference') or ''
        if ref and os.path.isfile(ref):
            return ref
    except (OSError, ValueError):
        pass
    p = os.path.join(DIRS['virus_src'], 'final.cluster.ref.fasta')
    return p if os.path.isfile(p) else None


def default_ref_info():
    """默认参考注释 TSV（新布局在库目录内，旧布局在 virusref_db 根下）。"""
    for p in (os.path.join(_kv_lib_dir(), 'reference.ref_info.tsv'),
              os.path.join(DIRS['virus_src'], 'final.cluster.ref_info.tsv')):
        if os.path.isfile(p):
            return p
    return None


def resolve_kv_lib(name_or_path):
    """把「鉴定库」解析成 (reference, ref_info, index_dir) 三件套。

    库与引擎正交的"选库"入口：一个库 = 参考真相 + 注释 + 各引擎索引集合，
    选库即三件套一起定，定量引擎（engine 参数）另行独立选择。
    接受：平台内库名（databases/virusref_db/<name>，推荐）或库目录绝对路径
    （含 platform.json 登记的外部库）。

    返回 dict；库不可用（目录不存在/无参考）时抛 ValueError，前端可直接展示。
    """
    import json as _json
    p = str(name_or_path or '').strip()
    if not p:
        raise ValueError('鉴定库为空')
    if os.path.isabs(p):
        lib = os.path.normpath(p)
    else:
        if not re_fullmatch_libname(p):
            raise ValueError(f'库名仅限字母数字-_: {p}')
        lib = os.path.normpath(os.path.join(DIRS['virus_src'], p))
    if not os.path.isdir(lib):
        raise ValueError(f'鉴定库目录不存在: {lib}')

    # 参考：① 库内 reference.fasta；② manifest.json 登记的路径（外部库）
    ref = os.path.join(lib, 'reference.fasta')
    manifest = {}
    mp = os.path.join(lib, 'manifest.json')
    if os.path.isfile(mp):
        try:
            with open(mp, encoding='utf-8') as fh:
                manifest = _json.load(fh) or {}
        except (OSError, ValueError):
            manifest = {}
    if not os.path.isfile(ref):
        ref = manifest.get('reference') or ''
        if ref and not os.path.isfile(ref):
            raise ValueError(f'鉴定库参考不存在（manifest 登记为 {ref}）: {lib}')
    if not ref:
        raise ValueError(f'鉴定库缺少 reference.fasta: {lib}')

    ref_info = os.path.join(lib, 'reference.ref_info.tsv')
    if not os.path.isfile(ref_info):
        ref_info = manifest.get('ref_info') or ''
        if ref_info and not os.path.isfile(ref_info):
            ref_info = ''
    return {'name': os.path.basename(lib), 'path': lib,
            'reference': ref, 'ref_info': ref_info or '', 'index_dir': lib}


def re_fullmatch_libname(name):
    """库名规则（与 web 建库 API 同一口径）：字母数字-_。"""
    import re as _re
    return bool(_re.fullmatch(r'[A-Za-z0-9_\-]+', name))


def index_dir_for(ref, engine='salmon'):
    """鉴定库索引的落点目录（命中即复用，缺失则引擎在该目录现建）。

    库与引擎正交：同一库目录下 salmon 索引在 <库>/salmon_k31/（小RNA 另有
    <库>/salmon_k15/），minibwa 索引在 <库>/minibwa/reference.mbw。

    平台内库（databases/virusref_db/<库名>/，2026-09-19 起从仅默认库
    kv_index 推广到所有平台内库）：直接指回库目录——不管引擎索引当前
    是否在位：已建则引擎「存在即复用」；未建（如分发包不带索引）则首次
    运行时在库目录内自建一次（k31 约 10s / k15 约 35s），全局复用。
    此前非默认库一律返回 None，引擎每次运行都在 <out>/index/ 现建索引，
    库稍大就明显拖慢每个任务。

    平台外库（外部登记，"只读使用"）：仅当库内已带该引擎的现成索引时
    返回库目录（纯只读复用）；否则返回 None，让引擎在 <out>/index/ 现建，
    绝不写外部目录。库目录只读（只读介质等）同理。
    （共识段的 minibwa 索引由共识段在运行目录内自建自用，不在这里复用。）
    """
    if engine not in ('salmon', 'minibwa'):
        return None
    lib = os.path.dirname(os.path.abspath(str(ref or '')))
    # 旧布局兜底：参考仍散置在 virusref_db 根（未迁移安装，且确为旧默认
    # 参考文件名）→ 索引照旧落默认库目录 kv_index（与 2026-09-12 参考
    # 目录重构前行为一致）；根下其它 FASTA 不认，避免错配 kv_index 索引
    if os.path.normcase(lib) == os.path.normcase(
            os.path.abspath(DIRS['virus_src'])) and \
            os.path.basename(os.path.abspath(str(ref or ''))) in \
            ('final.cluster.ref.fasta', 'final.cluster.ref.fa'):
        kv = os.path.join(DIRS['virus_src'], 'kv_index')
        return kv if os.access(kv, os.W_OK) else None
    # 平台内库（数据库根之下任意位置的 <库名>/reference.fasta，2026-09-19
    # 起不再限定 virusref_db 一层）：索引落库目录，缺失就地现建一次后全局复用
    dbroot = os.path.normcase(os.path.abspath(DIRS['databases']))
    libn = os.path.normcase(lib)
    if libn == dbroot or libn.startswith(dbroot + os.sep):
        # 只认库的"正身"参考（<库>/reference.fasta）：库目录里手放的
        # 其它 FASTA 不吃，否则会错配到该库的索引
        if os.path.basename(os.path.abspath(str(ref or ''))) != \
                'reference.fasta':
            return None
        if not os.path.isfile(os.path.join(lib, 'reference.fasta')):
            return None        # 参考未按库布局归档 → 索引位置无据可依
        if os.access(lib, os.W_OK):
            return lib         # 平台内可写库 → 索引落库目录，缺失就地现建
        return None            # 库目录只读 → 引擎在 <out>/index 现建
    # 平台外：已有现成索引才复用（只读），否则不写外部目录
    import os as _os
    if engine == 'salmon':
        ready = any(_os.path.isfile(_os.path.join(lib, *p)) for p in (
            ('salmon_k31', 'info.json'), ('salmon', 'salmon_k31', 'info.json')))
    else:
        ready = _os.path.isfile(_os.path.join(lib, 'minibwa', 'reference.mbw'))
    return lib if ready else None


_DTYPES = {'String', 'Int64', 'Int32', 'Float64', 'Float32', 'Boolean',
           'Utf8', 'Null', 'Object', 'Date', 'Datetime'}


def _read_tsv_rows(path):
    """读 TSV → (表头 list, 数据行 list)。跳过 polars 误写出的 dtype 行。"""
    if not os.path.isfile(path):
        return [], []
    with open(path, encoding='utf-8', errors='replace') as fh:
        lines = [ln.rstrip('\n') for ln in fh if ln.strip()]
    if not lines:
        return [], []
    head = lines[0].split('\t')
    rows = []
    for ln in lines[1:]:
        cols = ln.split('\t')
        if cols and all((c or '').strip() in _DTYPES for c in cols[:3]):
            continue                      # dtype 行（历史 bug 产物），不算数据
        rows.append(cols)
    return head, rows


def _col(head, *names):
    lut = {h.strip().lower(): i for i, h in enumerate(head)}
    for n in names:
        i = lut.get(n.strip().lower())
        if i is not None:
            return i
    return None


def _summarize(out_dir):
    """识别/定量结果 → {species, reads, passed, rejected, filtered_tsv}。

    识别与定量口径（audit 第 7 项修正）：
      * 物种名单（species_detected）取 **通过阈值**的结果：
        `identify/all_viruses.best.summary.tsv`（确诊表）→ 不存在则退
        `filter/filtered.tsv`（二次过滤表）→ 都不存在才退全量
        `identify/all_viruses.summary.tsv`，并在 species_basis 字段注明
        口径（passed_threshold / has_reads）。此前一律取未过滤全量表，
        任务卡「检出 N 个物种」把没通过阈值的也算了进去（虚高）。
      * species_with_reads：全量表里有 reads 的物种数（诊断用，与
        「通过阈值」区分）。
      * viral_reads：全量表 reads 求和（总数口径，不因阈值而变）。
      * 通过/剔除条数取 `filter/filtered.tsv` 与 `filter/discarded.tsv`
        （二次过滤 = 定量口径；多段物种会被节段完整性规则剔除）。
      * minibwa 的 summary 表对库内**每条参考**都出行（idxstats 口径），
        0 计数行若也计入物种，会把"库里有、样本里没有"的参考全算成检出
        （实测库内 8464 条参考时物种数虚高上百倍）。故 reads=0 的行只进
        reads 求和、不进物种名单。
    """
    full = os.path.join(out_dir, 'identify', 'all_viruses.summary.tsv')
    best = os.path.join(out_dir, 'identify', 'all_viruses.best.summary.tsv')
    filtered = os.path.join(out_dir, 'filter', 'filtered.tsv')
    discarded = os.path.join(out_dir, 'filter', 'discarded.tsv')

    def _sp_reads(path):
        head, rows = _read_tsv_rows(path)
        sp_i = _col(head, 'Species', 'Accession', 'Virus')
        rd_i = _col(head, 'Uniq_Reads', 'EM_Reads', 'Reads')
        species, reads = [], 0
        for cols in rows:
            n_reads = None
            if rd_i is not None and rd_i < len(cols):
                try:
                    n_reads = int(float(cols[rd_i]))
                except ValueError:
                    n_reads = None          # 读数列缺失/异常 → 不据此判定未检出
            if (n_reads is None or n_reads > 0) and sp_i is not None \
                    and sp_i < len(cols) and cols[sp_i].strip():
                name = cols[sp_i].strip()
                if name not in species:
                    species.append(name)
            if n_reads:
                reads += n_reads
        return species, reads

    # 物种名单：优先「通过阈值」（best.summary / filtered），
    # 都缺才退全量表（species_basis='has_reads' 注明口径）
    species_with_reads, reads = (_sp_reads(full) if os.path.isfile(full)
                                 else ([], 0))
    species, basis = species_with_reads, 'has_reads'
    for _p in (best, filtered):
        if os.path.isfile(_p):
            species, basis = _sp_reads(_p)[0], 'passed_threshold'
            break

    # 比对失败的样本：引擎在 identify/（旧布局在阶段根目录）落
    # failed_samples.json。失败样本与真阴性绝不可混同——不读这个文件，
    # salmon 崩掉的样本就会以「检出 0 物种」的假阴性面孔混进汇总。
    failed = []
    for cand in (os.path.join(out_dir, 'identify', 'failed_samples.json'),
                 os.path.join(out_dir, 'failed_samples.json')):
        if os.path.isfile(cand):
            try:
                with open(cand, encoding='utf-8') as fh:
                    failed = json.load(fh) or []
            except (OSError, ValueError):
                failed = []
            break
    return {
        'species_detected': species,
        'species_basis': basis,          # passed_threshold | has_reads
        'species_with_reads': species_with_reads,
        'viral_reads': reads,
        'passed': len(_read_tsv_rows(filtered)[1]),
        'rejected': len(_read_tsv_rows(discarded)[1]),
        'failed_samples': failed,
        'filtered_tsv': filtered if os.path.isfile(filtered) else '',
        'identify_tsv': full if os.path.isfile(full) else '',
    }


def _samtools():
    exe = get_config().tool('samtools')
    if not exe:
        raise FileNotFoundError('未找到 samtools（抽病毒 reads 需要）')
    return exe


def _input_fingerprint(path):
    """输入文件指纹（路径 + size + mtime）：断点复用前核对，防「换了输入
    重跑却拿到旧结果」。读不到 stat（文件刚被替换/网络抖动）时 size/mtime
    记 None——只有路径也一致才算命中。"""
    p = os.path.abspath(str(path))
    try:
        st = os.stat(p)
        return [p, st.st_size, st.st_mtime]
    except OSError:
        return [p, None, None]


# 可复现性字段（audit 第 15 项） ────────────────────────────────
_ENGINE_VERSION_CACHE = {}


def _engine_version(engine):
    """定量引擎版本（`<exe> --version` 探测一次并缓存；探测失败记 unknown）。"""
    if engine in _ENGINE_VERSION_CACHE:
        return _ENGINE_VERSION_CACHE[engine]
    ver = 'unknown'
    try:
        exe = get_config().tool(engine)
        if exe:
            r = subprocess.run([exe, '--version'], capture_output=True,
                               text=True, encoding='utf-8', errors='replace',
                               timeout=30)
            out = (r.stdout or '') or (r.stderr or '')
            ver = out.strip().splitlines()[0][:120] if out.strip() else 'unknown'
    except (OSError, subprocess.SubprocessError) as e:
        ver = f'unknown({type(e).__name__})'
    _ENGINE_VERSION_CACHE[engine] = ver
    return ver


def _index_k(index_dir, engine):
    """salmon 索引的 k（目录名 salmon_k<N> 或 index info.json 的 k 字段）。

    minibwa（BWA 0.7 seed-and-extend）索引无 k 概念，返回 None。
    """
    if engine != 'salmon' or not index_dir:
        return None
    cands = [str(index_dir), os.path.join(str(index_dir), 'salmon_k31')]
    for d in cands:
        m = re.match(r'salmon_k(\d+)$', os.path.basename(d.rstrip('\\/')))
        if m:
            return int(m.group(1))
    for d in cands:
        try:
            with open(os.path.join(d, 'info.json'), encoding='utf-8') as fh:
                k = int((json.load(fh) or {}).get('k') or 0)
            if k:
                return k
        except (OSError, ValueError):
            continue
    return None


def _default_thresholds():
    """识别/过滤段的默认阈值（本阶段不传阈值参数，实际生效的就是这些默认）。

    直接走 known_virus_suite.build_parser 的 argparse 默认值，避免手抄
    常量与 CLI 漂移；导入/解析失败返回 {}（summary 照写，不阻断）。
    """
    try:
        from .known_virus_suite.known_virus_suite import build_parser
        ns = build_parser().parse_known_args(['identify', '--out', '.'])[0]
        keys = ('min_cov', 'min_depth', 'min_reads', 'min_tpm',
                'min_poisson', 'ratio', 'ani_thresh')
        return {k: getattr(ns, k) for k in keys if hasattr(ns, k)}
    except Exception as e:  # noqa: BLE001 - 仅元信息，失败不阻断
        return {}


def _coverage_tool_used(out_dir, sample):
    """本样本实际生效的覆盖度后端（按产物判断，不猜参数）。"""
    align = os.path.join(out_dir, 'align')
    if os.path.isfile(os.path.join(
            align, f'{sample}_pandepth.pandepth.chr.stat.gz')):
        return 'pandepth'
    if os.path.isfile(os.path.join(align, f'{sample}_pandepth.samtools_cov.tsv')):
        return 'samtools'
    return 'none(无 BAM，覆盖度/深度为估算)'


def _bam_qual_missing(samtools_exe, bam, max_records=50):
    """BAM 是否无真实质量值（QUAL 列全为 '*'）。

    salmon --writeBam 的映射记录无 QUAL（0xff），samtools fastq 会把质量列
    写成占位 'B'（Phred 33，看着像高质量）——下游组装/变异若信任这些质量
    值会得到假的「高质量」结论。抽取后显式检测，返回 True/False/None（不确定）。
    """
    if not samtools_exe:
        return None
    try:
        proc = subprocess.Popen([samtools_exe, 'view', str(bam)],
                                stdout=subprocess.PIPE,
                                stderr=subprocess.DEVNULL,
                                text=True, encoding='utf-8', errors='replace')
    except OSError:
        return None
    n, saw = 0, False
    try:
        for line in proc.stdout:
            if not line or line.startswith('@'):
                continue
            p = line.split('\t')
            if len(p) < 11:
                continue
            saw = True
            if p[10] != '*':
                return False          # 任一记录有真实 QUAL 即视为有质量
            n += 1
            if n >= max_records:
                break
    finally:
        # 只看头部几十条，提前终止，不等整个 BAM 流完
        try:
            proc.stdout.close()
        except OSError:
            pass
        proc.terminate()
        proc.wait()
    return True if saw else None


def _extract_viral_reads(out_dir, bam, r2, threads, logger=None):
    """按名字排序后抽 fastq。

    PE：viral_R1/viral_R2（成对）+ viral_se（单端孤儿——只有一条 mate
    比中病毒的 reads；不单独落盘的话 samtools 会把它们按 flag 交错写进
    R1/R2，两文件行数错位，SPAdes -1/-2 按行配对直接崩溃）。
    SE：全部并入 viral_R1。
    返回 (r1, r2, se)。
    """
    st_exe = _samtools()
    name_sorted = os.path.join(out_dir, 'align', 'viral.namesort.bam')
    se_out = os.path.join(out_dir, VIRAL_SE)
    cmd1 = [st_exe, 'sort', '-n', f'-@{threads}', '-o', name_sorted, bam]
    cmd2 = [st_exe, 'fastq', f'-@{threads}',
            '-1', os.path.join(out_dir, VIRAL_R1)]
    if r2:
        cmd2 += ['-2', os.path.join(out_dir, VIRAL_R2),
                 '-s', se_out,
                 '-0', se_out]      # 两个 mate 都没比中的理论上不存在，兜底并入单端
    else:
        cmd2 += ['-0', os.path.join(out_dir, VIRAL_R1)]
    cmd2 += [name_sorted]
    # 走 run_cmd（非裸 subprocess）：进任务取消注册表（「停止」能杀掉进程树）
    # + 默认超时看门狗，samtools 挂死不再拖死整个任务
    for cmd in (cmd1, cmd2):
        run_cmd(cmd, logger=logger)
        if logger:
            logger.log(f'  抽病毒 reads: samtools {cmd[1]} 完成')
    try:
        os.remove(name_sorted)
    except OSError:
        pass
    r1_out = os.path.join(out_dir, VIRAL_R1)
    r2_out = os.path.join(out_dir, VIRAL_R2)
    # 空 gz（samtools 仍会写出 ~20 字节头）按不存在处理
    se_exists = False
    if r2 and os.path.isfile(se_out) and os.path.getsize(se_out) > 24:
        try:
            import gzip as _gz
            with _gz.open(se_out, 'rt', encoding='utf-8',
                          errors='replace') as f:
                se_exists = sum(1 for _ in zip(range(4), f)) >= 4
        except (OSError, EOFError):
            se_exists = False
    return (r1_out,
            r2_out if (r2 and os.path.isfile(r2_out)) else None,
            se_out if se_exists else None)


def run_kvsuite_stage(out_dir, r1, r2=None, sample=None, threads=8,
                      engine='salmon', force=False, need_viral_reads=True,
                      logger=None, progress=None, kv_lib=None):
    """跑 ②b 识别+定量，返回结果字典（同时落 summary.json）。

    engine 默认 salmon（EM 定量）；minibwa=真比对计数，两引擎共用同一
    覆盖度后端与 BAM 命名。kv_lib：鉴定库（平台内库名或库目录绝对路径，
    经 resolve_kv_lib 整体解析参考/注释/索引），缺省用平台默认库——
    一键分析管道由此与其他模块共用同一套「选库」入口，不再硬编码 kv_index。
    断点续跑：summary.json 已存在、不是 force、且输入指纹（r1/r2/参考的
    路径+size+mtime）与上次完全一致时才直接复用——换了 FASTQ 或换了
    鉴定库重跑会真正重跑，绝不拿旧结果冒充新结果。
    """
    out_dir = check_path(out_dir, must_exist=False, in_platform=True)
    os.makedirs(out_dir, exist_ok=True)
    summary_p = os.path.join(out_dir, SUMMARY)
    lib = None
    if kv_lib and str(kv_lib).strip():
        lib = resolve_kv_lib(kv_lib)
    ref = (lib or {}).get('reference') or default_reference()
    if not ref:
        raise FileNotFoundError(
            '缺少病毒参考：请先在「数据库构建 → 病毒鉴定库」建库'
            '（databases/virusref_db/plant_ref/）')
    if os.path.isfile(summary_p) and not force:
        done = None
        try:
            with open(summary_p, encoding='utf-8') as fh:
                done = json.load(fh)
        except (OSError, ValueError):
            pass
        # 复用指纹必须覆盖**实际影响结果的全部要素**（审查 P1-1）：旧实现
        # 只比 r1/r2/reference——engine 换了（salmon↔minibwa）、阈值调了、
        # ref_info 注释换了、need_viral_reads 开关变了，外层照样命中缓存，
        # 内层批次级的细指纹根本没机会参与判定。
        _ref_info_fp = _input_fingerprint(
            (lib or {}).get('ref_info') or default_ref_info()
        ) if (((lib or {}).get('ref_info')) or default_ref_info()) else None
        want_fp = {
            'r1': _input_fingerprint(r1),
            'r2': _input_fingerprint(r2) if r2 else None,
            'reference': _input_fingerprint(ref),
            'ref_info': _ref_info_fp,
            'engine': engine,
            'need_viral_reads': bool(need_viral_reads),
            'thresholds': _default_thresholds(),
        }
        prev_fp = dict(done.get('inputs_fp') or {})
        prev_repro = done.get('repro') or {}
        prev_fp.setdefault('engine', prev_repro.get('engine'))
        prev_fp.setdefault('thresholds', prev_repro.get('thresholds'))
        if done and prev_fp == want_fp:
            if logger:
                logger.log(f'②b 已知病毒识别与定量已完成且输入/引擎/阈值'
                           f'未变化，跳过（{summary_p}）')
            return done
        if logger:
            logger.log('②b 存在历史断点，但输入、引擎、参考注释或阈值已'
                       '变化——重跑本阶段（不复用旧结果）')

    sample = sample or os.path.basename(os.path.dirname(out_dir.rstrip(os.sep)))
    entry_exe, entry_pre = engine_cmd()
    log_p = os.path.join(out_dir, 'kvsuite.log')
    t0 = time.time()

    def _run(sub, extra):
        cmd = [entry_exe] + entry_pre + [sub, '--out', out_dir,
                                         '--reference', ref,
                                         '--engine', engine,
                                         '--threads', str(threads),
                                         '--align-threads', str(threads),
                                         '--sample-name', sample]
        if sub != 'filter':                       # filter 只用已有结果表
            cmd += ['-1', r1]
            if r2:
                cmd += ['-2', r2]
        ri = (lib or {}).get('ref_info') or default_ref_info()
        if ri and os.path.isfile(ri):
            cmd += ['--ref-info', ri]
        idx = index_dir_for(ref, engine)
        if idx:
            cmd += ['--index-dir', idx]
        cmd += list(extra)
        if logger:
            logger.log(f'  ②b kvsuite {sub}: {" ".join(cmd[1 + len(entry_pre):])}')
        with open(log_p, 'a', encoding='utf-8') as fl:
            fl.write(f'\n===== {sub} =====\n')
        # 走 run_cmd_redirect（非裸 subprocess.run）：引擎进程进任务取消
        # 注册表（「停止」可杀整棵进程树，salmon/samtools 不再成为孤儿），
        # 且带默认超时看门狗；输出继续追加进 kvsuite.log
        try:
            run_cmd_redirect(cmd, log_p, append=True)
        except RuntimeError as e:
            tail = ''
            try:
                with open(log_p, encoding='utf-8', errors='replace') as fl:
                    tail = ''.join(fl.readlines()[-8:])
            except OSError:
                pass
            raise RuntimeError(f'kvsuite {sub} 失败：{tail or e}') from e

    _run('identify', [])
    if progress:
        progress(0.7, '识别完成，二次过滤（定量）')
    _run('filter', [])

    stats = _summarize(out_dir)
    bam = os.path.join(out_dir, 'align', f'{sample}_pandepth.sorted.bam')
    if not os.path.isfile(bam):
        bam = os.path.join(out_dir, 'align', f'{sample}.sorted.bam')
    viral_r1 = viral_r2 = viral_se = None
    if need_viral_reads and os.path.isfile(bam):
        if progress:
            progress(0.9, '抽病毒 reads（供 ③组装）')
        viral_r1, viral_r2, viral_se = _extract_viral_reads(
            out_dir, bam, r2, threads, logger=logger)

    # 质量值声明（audit 第 14 项）：salmon 映射 BAM 无 QUAL，samtools fastq
    # 会把质量列写成占位 'B'（看着像高质量）。抽完 reads 显式检测+告警，
    # 不让下游组装/变异把占位质量当真实质量。
    viral_quality = ''
    if viral_r1:
        try:
            _noqual = _bam_qual_missing(get_config().tool('samtools'), bam)
        except Exception:  # noqa: BLE001 - 仅声明性字段，失败不阻断
            _noqual = None
        if _noqual:
            viral_quality = ('salmon_bam_no_qual(fastq 质量列为占位)'
                             if engine == 'salmon'
                             else 'bam_no_qual(fastq 质量列为占位)')
            if logger:
                logger.warning(
                    '抽出的病毒 reads 无真实质量值（BAM QUAL 为 *，samtools '
                    "fastq 写出占位质量 'B'）：下游组装/变异不应信任这些质量值")

    reads = stats['viral_reads']
    # 定量口径归一：salmon 的 NumReads 计**片段**（PE 一对 = 1），minibwa 的
    # 映射记录数计 **reads**（PE 一对 = 2）。viral_pairs 统一归一到「片段/
    # 对」——③组装 min_viral_pairs 门槛与前端展示都用它，此前不分引擎一律
    # `reads//2`，salmon+PE 样本会被腰斩一半（可能误触发"病毒 reads 不足"）。
    if engine == 'salmon':
        reads_unit, viral_pairs = 'fragments', int(reads)
    else:
        reads_unit, viral_pairs = 'reads', (int(reads) // 2 if r2 else int(reads))
    result = {
        'stage': 'kvsuite',
        'sample': sample,
        'engine': engine,
        'kv_lib': (lib or {}).get('name', '') if lib else '',
        'reference': ref,
        'index_dir': index_dir_for(ref, engine) or '',
        'identify_tsv': stats['identify_tsv'],
        'filtered_tsv': stats['filtered_tsv'],
        'bam': bam if os.path.isfile(bam) else '',
        'species_detected': stats['species_detected'],
        # 口径声明（audit 第 7 项）：passed_threshold=通过阈值表；
        # has_reads=通过表缺失、退全量表（物种数偏虚高，只作回退）
        'species_basis': stats.get('species_basis', ''),
        'species_with_reads': stats.get('species_with_reads', []),
        'viral_reads': reads,
        'reads_unit': reads_unit,
        'viral_pairs': viral_pairs,
        # 抽出 reads 的质量值声明（'' = 有真实质量；否则为占位质量说明）
        'viral_reads_quality': viral_quality,
        'passed': stats['passed'],
        'rejected': stats['rejected'],
        'viral_r1': viral_r1 or '',
        'viral_r2': viral_r2 or '',
        'viral_se': viral_se or '',
        'failed_samples': stats.get('failed_samples', []),
        'elapsed_sec': round(time.time() - t0, 1),
        # 可复现性字段（audit 第 15 项）：引擎版本/索引 k/libtype/阈值，
        # 缺这些字段时同一份结果无法复现（阈值与索引换了结果就变）。
        'repro': {
            'engine': engine,
            'engine_version': _engine_version(engine),
            'index_k': _index_k(index_dir_for(ref, engine), engine),
            'libtype': 'A' if engine == 'salmon' else '',
            'coverage_tool': _coverage_tool_used(out_dir, sample),
            'thresholds': _default_thresholds(),
        },
        # 输入指纹：下次断点续跑时核对（P1-1 扩展口径：含引擎/注释/阈值/
        # need_viral_reads，与复用判据完全一致），输入变了就重跑而不是复用
        'inputs_fp': {
            'r1': _input_fingerprint(r1),
            'r2': _input_fingerprint(r2) if r2 else None,
            'reference': _input_fingerprint(ref),
            'ref_info': (_input_fingerprint((lib or {}).get('ref_info')
                                            or default_ref_info())
                         if (((lib or {}).get('ref_info'))
                             or default_ref_info()) else None),
            'engine': engine,
            'need_viral_reads': bool(need_viral_reads),
            'thresholds': _default_thresholds(),
        },
    }
    with open(summary_p, 'w', encoding='utf-8') as fh:
        json.dump(result, fh, ensure_ascii=False, indent=2)
    if logger:
        _fail_tip = ''
        if stats.get('failed_samples'):
            _fail_tip = (f'，⚠ {len(stats["failed_samples"])} 个样本比对失败'
                         '（未计入结果，勿当作阴性！详见 failed_samples.json）')
        # 「通过阈值」与「有 reads」两个口径都报（audit 第 7 项）
        _basis_tip = ''
        if result['species_basis'] == 'passed_threshold':
            _basis_tip = (f'（通过阈值口径；有 reads 的物种 '
                          f'{len(result["species_with_reads"])} 个）')
        else:
            _basis_tip = ('（⚠ 未找到通过阈值结果表，物种数退为「有 reads」'
                          f'口径 {len(result["species_with_reads"])} 个，'
                          '偏虚高）')
        _qual_tip = (f'，⚠ 抽出 reads 质量值为占位（{viral_quality}）'
                     if viral_quality else '')
        logger.log(f'②b 已知病毒识别与定量完成：检出 '
                   f'{len(result["species_detected"])} 个物种{_basis_tip}'
                   f'（通过 {result["passed"]} / 剔除 {result["rejected"]}），'
                   f'病毒 {result["viral_reads"]:,} {reads_unit}'
                   f'（≈{result["viral_pairs"]:,} 片段/对）'
                   f'{_qual_tip}{_fail_tip}，'
                   f'用时 {result["elapsed_sec"]}s')
    return result
