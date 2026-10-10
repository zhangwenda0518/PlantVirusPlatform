# -*- coding: utf-8 -*-
"""共识序列与变异分析（验证模块之后）。

思路来源
--------
MultiVirusConsensus（Niema Lab）：把 reads 回贴到**多个**参考、逐参考统计
位置计数、据覆盖度判断哪些病毒真实存在。

为什么不是直接调用 MVC
----------------------
MVC 是纯编排层，硬性依赖 minimap2 + samtools + viral_consensus 三个可执行文件。
后两者在 Windows 上无官方构建（bioconda 的 minimap2/samtools 只有 linux-64 /
osx-64 / osx-arm64，无 win-64）。因此本模块：

- **比对**：minibwa（平台自带 0.7-r424 Windows 本地编译，静态单文件），
  短读用 `-x sr`。2026-09-08 由 minimap2 换为 minibwa（同一作者 lh3 的
  下一代比对器，可单文件静态分发）；两者在 88x 深度下经本模块
  pileup + call_consensus 后的共识序列完全一致，对拍详情见
  minibwa_win_build/README.md。注意 minibwa 映射前需先 `index` 建索引；
  双端 reads 用 `map 索引 R1 R2` 两个位置参数（实测输出带配对 FLAG
  99/147、RNEXT=`、PNEXT 互相指向），R1/R2 分次比对会丢配对信息，不采用
- **共识调用**：用 Python 重写 ViralConsensus 的算法（逐位置统计碱基，
  深度 ≥ min_depth 且频率 ≥ min_freq 才 call，否则给简并符）——
  与 sdt_exact.py 纯 Python 替代 SDTv1.3.exe 同一套路
- **samtools 不需要**：minibwa 直接输出 SAM，这里流式解析，省掉排序与格式转换

内存只与参考长度相关、与测序深度无关（ViralConsensus 的核心优点）得以保留。

产物
----
对齐 ViralConsensus（同名的三个产物语义一致）：
    consensus.fa       每条参考一条共识序列（与参考等长，低覆盖区为简并符）
    poscounts.tsv      逐位置碱基计数（每参考一个文件）
    inscounts.json     逐位置插入计数（每参考一个文件）
本模块附加（变异分析）：
    variants.tsv       变异位点（SNV / iSNV 准种 / 简并 / 低覆盖）
    coverage.tsv       每条参考的覆盖统计 + 是否判定为真实存在
    summary.json       摘要（供前端渲染）

注意：poscounts.tsv 比 ViralConsensus 多一列 N（未定碱基）。ViralConsensus
原格式为 `Pos A C G T - Total`，这里为 `Pos A C G T N - Total`。
"""
import os
import re
import json
import gzip
import shutil
from collections import Counter

from .config import get_config
from .utils import run_cmd

# ---------------------------------------------------------------- 常量
# 参数默认值对齐 ViralWasm-Consensus
MIN_BASE_QUALITY = 20      # -q / --min_qual：最小碱基质量
MIN_DEPTH = 4              # -d / --min_depth：最小覆盖深度（2026-10-03 定为 4，即"深度>3"：
                           # 原值 1 时单条 read（例如 1x 处的测序错误/污染）
                           # 即可定共识碱基并产出 SNV，低覆盖位点噪声全部入表；
                           # >3 与 ViralConsensus/samtools 常规口径一致，深度≤3
                           # 的位点共识仍给简并符、SNV 不再产出。要更宽松可显式
                           # 传小值（参数保留可配置）
MIN_FREQ = 0.5             # -f / --min_freq：最小等位基因频率
AMBIG = 'N'                # -a / --ambig：简并符号

MIN_MAPQ = 10              # 最小比对质量（低于此丢弃，避免多重比对污染）
MIN_MINOR_FREQ = 0.02      # 次等位基因频率阈值：≥ 此值记为 iSNV（准种/混合）
MIN_MINOR_COUNT = 2        # iSNV 的绝对 reads 数下限：1 条 read 在 10x 处就是 10%，
                           # 不设下限会把测序错误当成准种
MIN_COV_PCT = 10.0         # 判定「该病毒真实存在」的最低覆盖度百分比

# minibwa 的种子下限：实测 20-29bp reads 在任何种子参数下比对率都是 0
# （见 minibwa_win_build 对拍记录）。小RNA/siRNA 数据（21-24nt）落到这里时
# 覆盖/共识/变异全为 0，但那是工具不支持、不是"病毒不存在"——必须区分。
MINIBWA_MIN_READ_LEN = 30

BASES = 'ACGT'
IDX = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4, '-': 5}
NB = 6                     # A C G T N -

# 对齐 ViralConsensus common.h 的 FAIL_FLAGS
FAIL_FLAGS = 0x4 | 0x100 | 0x200 | 0x400   # UNMAP SECONDARY QCFAIL DUP

CIGAR_RE = re.compile(r'(\d+)([MIDNSHP=X])')

COVERAGE_FIELDS = [
    'contig', 'length', 'mapped_reads', 'mean_depth', 'median_depth',
    'min_depth', 'max_depth', 'covered_bp', 'coverage_pct', 'depth_ok_pct',
    'called_pct', 'n_variants', 'n_snv', 'n_isnv', 'present',
]

VARIANT_FIELDS = [
    'contig', 'pos', 'ref_base', 'cons_base', 'depth',
    'A', 'C', 'G', 'T', 'N', 'gap',
    'major_base', 'major_freq', 'minor_base', 'minor_freq', 'type',
]


def consensus_engine():
    """探测 minibwa。返回路径或 None（缺失时本模块不可用）。

    2026-09-08 起由 minibwa 承担回贴（minimap2 作者的下一代比对器）。
    实测对拍：88x 深度下两者经 pileup+call_consensus 后的共识序列完全一致。
    """
    cfg = get_config()
    for name in ('minibwa',):
        try:
            p = cfg.tool(name)
        except Exception:
            p = None
        if p and os.path.isfile(p):
            return p
    return None


def load_fasta(path):
    """读 FASTA。返回 OrderedDict-like {id: seq}（序列转大写）。"""
    out = {}
    name, buf = None, []
    op = gzip.open if str(path).endswith('.gz') else open
    with op(path, 'rt', encoding='utf-8', errors='replace') as f:
        for line in f:
            if line.startswith('>'):
                if name is not None:
                    out[name] = ''.join(buf).upper()
                name = line[1:].strip().split()[0]
                buf = []
            elif name is not None:
                buf.append(line.strip())
    if name is not None:
        out[name] = ''.join(buf).upper()
    return out


def _write_fasta(path, records, width=70):
    """写 FASTA。records: [(name, seq), ...] 或 {name: seq}。"""
    items = records.items() if isinstance(records, dict) else records
    with open(path, 'w', encoding='utf-8', newline='\n') as f:
        for name, seq in items:
            f.write('>%s\n' % name)
            for i in range(0, len(seq), width):
                f.write(seq[i:i + width] + '\n')
    return path


# R1/R2 命名约定：`_R1` / `.R1` / `-R1` 三种分隔（如 samp_R1.fq.gz、
# GQMIX.R1.fq.gz）——大小写不敏感，且必须紧跟 .fq/.fastq(.gz) 后缀，
# 避免把名字里恰好含 "R1" 的其他文件误判成 reads。
R1_RX = re.compile(r'^(.*?)([_.\-])R1(\.(?:fq|fastq)(?:\.gz)?)$', re.I)


def _r1_mate_name(fn):
    """R1 文件名 → 对应 R2 文件名；不是 R1 命名返回 None。"""
    m = R1_RX.match(fn)
    if not m:
        return None
    return m.group(1) + m.group(2) + 'R2' + m.group(3)


def find_read_pairs(d):
    """在目录下找 R1/R2 配对 reads，返回 [(r1, r2_or_None), ...]。

    工具入口与管道阶段共用：管道里 reads 固定在 01_host_removal 等层级目录，
    文件名前缀不统一（kept_R1 / viral_R1 / fastp_R1），按 R1 命名约定识别
    （见 R1_RX：`_R1` / `.R1` / `-R1` 均可）。只找到 R1 时按单端返回。
    """
    if not os.path.isdir(d):
        return []
    out = []
    for fn in sorted(os.listdir(d)):
        mate = _r1_mate_name(fn)
        if not mate:
            continue
        p1 = os.path.join(d, fn)
        p2 = os.path.join(d, mate)
        out.append((p1, p2 if os.path.isfile(p2) else None))
    return out


def pair_read_files(paths):
    """把手动指定的 reads 文件列表配成 [(r1, r2_or_None), ...]。

    R1 的伴侣只在给定列表内查找（R1/R2 两个框分别提交、R2 排在 R1 前
    都能配上）；不自动到磁盘上找 R2——用户只给 R1 就是单端。
    目录请用 find_read_pairs。
    """
    real = os.path.realpath
    by_path = {}
    for p in paths:
        if p and os.path.isfile(p):
            by_path.setdefault(real(p), p)
    out, consumed = [], set()
    # 第一遍：凡 R1 且列表内有 R2 伴侣的先配对
    for p in paths:
        if not p:
            continue
        rp = real(p)
        if rp not in by_path or rp in consumed:
            continue
        mate = _r1_mate_name(os.path.basename(p))
        if not mate:
            continue
        sib = os.path.join(os.path.dirname(p), mate)
        rs = real(sib)
        if rs in by_path and rs != rp and rs not in consumed:
            consumed.update((rp, rs))
            out.append((p, by_path[rs]))
    # 第二遍：落单的保持单端
    for p in paths:
        if not p:
            continue
        rp = real(p)
        if rp in by_path and rp not in consumed:
            consumed.add(rp)
            out.append((p, None))
    return out


# ---------------------------------------------------------------- 比对
# minibwa preset 语系与 minimap2 不同：minimap2 的 sr/map-ont/map-pb/map-hifi
# 四类在 minibwa 里归一为 sr（短读）/ lr（长读）。映射表写死，未知值报错而不静默回退。
PRESET_MAP = {
    'sr': 'sr',          # 短读
    'lr': 'lr',          # 长读
    'adap': 'adap',      # 自适应（长短混测）
    'map-ont': 'lr',
    'map-pb': 'lr',
    'map-hifi': 'lr',
}


def map_preset(preset):
    """把 minimap2 语系的 preset 映射到 minibwa 语系。

    minibwa 只有 sr/lr/adap 三种；minimap2 的长读 preset 统一走 lr。
    未登记的值直接报错，不静默回退——回退成错的 preset 会让长读被当短读处理。
    """
    p = str(preset or 'sr').strip().lower()
    if p not in PRESET_MAP:
        raise ValueError(
            '不支持的 preset: %r（minibwa 仅支持 %s）'
            % (preset, ', '.join(sorted(set(PRESET_MAP.values())))))
    return PRESET_MAP[p]


def build_index(ref_fa, logger=None):
    """minibwa index：生成 <ref_fa>.l2b + <ref_fa>.mbw，返回索引前缀。

    这两份索引是 minibwa map 的输入，缺一不可。
    """
    exe = consensus_engine()
    if not exe:
        raise RuntimeError('未找到 minibwa（bin/），无法回贴')
    run_cmd([exe, 'index', ref_fa], logger=logger)
    for ext in ('.l2b', '.mbw'):
        if not os.path.isfile(ref_fa + ext):
            raise RuntimeError('minibwa 索引未生成: %s%s' % (ref_fa, ext))
    return ref_fa


def _stage_ref_for_index(ref_fa, out_sam, logger=None):
    """给 minibwa 索引准备独占的 ASCII 安全临时目录（把参考复制进去）。

    为什么不在参考 FASTA 旁建索引（旧行为）：`.l2b/.mbw` 落在 run_dir 的
    参考旁边、比对完 finally 删除——同一参考的并发运行会互删对方的索引
    （A 建好、B 跑完删掉，A 的 map 中途失败）。改为每次调用复制参考到本次
    运行独占的临时目录，建索引→比对→整目录清理，彼此不相干。

    参考本身可能位于含中文的样品目录，而 minibwa 与 blast/mafft 同属对
    非 ASCII 路径敏感的 C 工具：索引目录优先走 %TEMP% 中转
    （assembly._ascii_work_base 的同款惯例），不可用时退回产物目录并 WARN。
    返回 (索引目录, 目录内参考路径)；调用方负责 rmtree 整个目录。
    """
    import tempfile
    try:
        from .assembly import _ascii_work_base
        base = _ascii_work_base('vp_mbw')
    except Exception as e:
        base = os.path.dirname(os.path.abspath(str(out_sam))) or '.'
        if logger:
            logger.log(f'ASCII 临时目录不可用（{e}），minibwa 索引目录改落 '
                       f'{base}', 'WARN')
    d = tempfile.mkdtemp(prefix='mbw_idx_', dir=base)
    local_ref = os.path.join(d, 'ref.fa')
    shutil.copyfile(ref_fa, local_ref)
    return d, local_ref


def run_minimap2(ref_fa, reads, out_sam, threads=None, preset='sr',
                 logger=None):
    """minibwa 回贴。reads 项为 str（单端）或 (r1, r2) 二元组（双端）。

    与 minimap2 的接口差异（实测确认，见 minibwa_win_build/README.md）：

    - **必须先用 index 建索引**，map 的输入是索引前缀而非 FASTA
    - **双端比对**：`map 索引 R1 R2` 两个位置参数，输出带配对 FLAG
      （99/147、RNEXT=`、PNEXT 互相指向）；R1/R2 分次比对会丢配对信息，
      旧版这里的「逐个比对后合并 SAM」写法已废弃
    - **preset 语系不同**：minimap2 的 map-ont/map-pb/map-hifi 均映射为 lr
    - `--secondary=no` 的等价物是 `-N 0`

    多个样本/多个 read 对时逐个比对到临时 SAM 再合并（保留一份 header），
    各样本 reads 汇入同一张覆盖统计——与旧版口径一致。

    索引不再写在参考旁（并发互删），改落本次运行独占的临时目录（见
    _stage_ref_for_index），函数退出前整目录清理。

    注意：minibwa 的 stderr 是 GBK 编码，平台 run_cmd 已用
    errors='replace' 兜底（报错文本可能含替换符，进程不受影响）。
    """
    exe = consensus_engine()
    if not exe:
        raise RuntimeError('未找到 minibwa（bin/），无法回贴')
    # 归一化成 [(r1, r2_or_None), ...]；容忍 None 项与缺失文件
    units = []
    for x in reads or []:
        if not x:
            continue
        if isinstance(x, str):
            r1, r2 = x, None
        else:
            r1, r2 = (list(x) + [None, None])[:2]
        if r1 and os.path.isfile(r1):
            units.append((r1, r2 if (r2 and os.path.isfile(r2)) else None))
    if not units:
        raise RuntimeError('没有可用的 reads 文件')

    nthread = str(threads or get_config().threads or 4)
    preset_mb = map_preset(preset)

    idx_dir, local_ref = _stage_ref_for_index(ref_fa, out_sam, logger=logger)
    try:
        idx = build_index(local_ref, logger=logger)

        def _map_cmd(dst, r1, r2):
            cmd = [exe, 'map', '-x', preset_mb, '-N', '0', '-t', nthread,
                   '-o', dst, idx, r1]
            if r2:
                cmd.append(r2)
            return cmd

        if len(units) == 1:
            r1, r2 = units[0]
            run_cmd(_map_cmd(out_sam, r1, r2), logger=logger)
            return out_sam

        parts = []
        try:
            for i, (r1, r2) in enumerate(units):
                part = '%s.part%d.sam' % (out_sam, i)
                run_cmd(_map_cmd(part, r1, r2), logger=logger)
                parts.append(part)

            seen_header = False
            with open(out_sam, 'w', encoding='utf-8', errors='replace') as fo:
                for part in parts:
                    with open(part, 'r', encoding='utf-8',
                              errors='replace') as fi:
                        for line in fi:
                            if line.startswith('@'):
                                if not seen_header:
                                    fo.write(line)
                            else:
                                fo.write(line)
                    seen_header = True
        finally:
            # 清理放在 finally：第 2 个样本起比对失败时，已产出的
            # part SAM（可达数 GB）不再残留 03c_consensus/
            for part in parts:
                try:
                    os.remove(part)
                except OSError:
                    pass
        return out_sam
    finally:
        # 本次运行独占的索引目录（参考副本 + .l2b/.mbw）用完即整目录清理
        shutil.rmtree(idx_dir, ignore_errors=True)


# ---------------------------------------------------------------- pileup
def pileup(sam_path, refs, min_qual=MIN_BASE_QUALITY, min_mapq=MIN_MAPQ,
           logger=None, progress=None):
    """逐位置统计碱基。返回 {ref: {'counts': [[A,C,G,T,N,-]...], 'ins': {pos: Counter}, 'nreads': int}}

    内存与测序深度无关，只随参考长度线性增长。
    插入语义对齐 ViralConsensus：记录「该位置之前」发生的插入序列；
    只登记**非空**插入（ViralConsensus 会为每个覆盖位置写 "" 空插入，
    高深度下会让 json 膨胀，这里省略——覆盖深度看 poscounts 的 Total 即可）。

    返回 dict 含 'qual_unavailable'：输入 QUAL 全为 '*'（无碱基质量）时为
    True，表示 min_qual 过滤实际未生效（minibwa 等工具的输出即如此）。
    """
    counts = {n: [[0] * NB for _ in range(len(s))] for n, s in refs.items()}
    ins = {n: {} for n in refs}
    nreads = Counter()
    # 是否实测到 QUAL=='*'（质量不可用）的比对记录：此时 min_qual 不过滤
    # 任何碱基，等于"质量过滤未生效"——必须让调用方/摘要可见（见下）。
    qual_unavailable = False

    is_gz = str(sam_path).endswith('.gz')
    total = 0 if is_gz else (os.path.getsize(sam_path) or 1)
    done = 0
    n_lines = 0
    last_pct = -1
    op = gzip.open if is_gz else open
    with op(sam_path, 'rt', encoding='utf-8', errors='replace') as f:
        for line in f:
            if line[0] == '@':
                continue
            done += len(line.encode('utf-8', 'replace')) if is_gz else len(line)
            n_lines += 1
            p = line.rstrip('\n').split('\t')
            if len(p) < 11:
                continue
            try:
                flag = int(p[1])
            except ValueError:
                continue
            if flag & FAIL_FLAGS:
                continue
            rname = p[2]
            if rname == '*' or rname not in refs:
                continue
            try:
                mapq = int(p[4])
            except ValueError:
                mapq = 0
            if mapq < min_mapq:
                continue

            cigar, seq, qual = p[5], p[9].upper(), p[10]
            if cigar == '*':
                continue
            # SAM 规范：QUAL 为 '*' 表示质量不可用。htslib（BAM）将其填充为
            # 0xFF(255)，ViralConsensus 直接读该值，因此这类记录一律通过
            # min_qual 过滤。若不过滤掉这里的 * 字符（ord('*')-33=9）会把
            # 所有缺质量记录误判为低质量而静默丢弃（实测 minibwa 输出即为此情形）。
            if qual == '*':
                qual = None
                if not qual_unavailable:
                    qual_unavailable = True
                    if logger:
                        logger.log(
                            '输入比对记录无碱基质量（QUAL 为 *，htslib 以 '
                            '255 填充）：min_qual 质量过滤不适用，'
                            '本结果等效于关闭质量过滤', 'WARN')
            try:
                rp = int(p[3]) - 1
            except ValueError:
                continue
            rows = counts[rname]
            rlen = len(rows)
            insmap = ins[rname]
            qp = 0
            insbuf = []
            nreads[rname] += 1

            for m in CIGAR_RE.finditer(cigar):
                ln = int(m.group(1))
                opc = m.group(2)
                if opc in ('M', '=', 'X'):
                    for i in range(ln):
                        qi = qp + i
                        ri = rp + i
                        if ri >= rlen:
                            break
                        # SEQ 比 CIGAR 短（异常/截断 SAM）时跳过该位点计数，
                        # 不让 IndexError 打断整阶段——此前只护了 QUAL 没护
                        # SEQ（2026-10-04 审查 E4）
                        if qi >= len(seq):
                            continue
                        # 碱基质量过滤（Phred+33）；qual=None 表示 SAM 的 '*'
                        # （质量不可用），按 htslib 语义视为 255 直接通过
                        if qual is None:
                            q = 255
                        else:
                            q = (ord(qual[qi]) - 33) if qi < len(qual) else 0
                        if q >= min_qual:
                            rows[ri][IDX.get(seq[qi], 4)] += 1
                        if insbuf:
                            s = ''.join(insbuf)
                            c = insmap.get(ri)
                            if c is None:
                                insmap[ri] = Counter({s: 1})
                            else:
                                c[s] += 1
                            insbuf = []
                    rp += ln
                    qp += ln
                elif opc == 'I':
                    insbuf.append(seq[qp:qp + ln])
                    qp += ln
                elif opc == 'D':
                    for i in range(ln):
                        ri = rp + i
                        if ri < rlen:
                            rows[ri][5] += 1
                    rp += ln
                elif opc == 'S':
                    qp += ln
                elif opc == 'N':
                    rp += ln
                # H：都不前进

            if progress:
                if total:
                    # 文本 SAM：分子分母同为字符量纲（2026-10-04 审查 P3：
                    # 此前 gz 的分母是压缩字节，pct 恒被截断在 99%）
                    pct = min(99, int(done * 100 / total))
                    if pct != last_pct and pct % 10 == 0:
                        last_pct = pct
                        progress('pileup', 0.4 + 0.4 * pct / 100,
                                 '统计位置计数 %d%%' % pct)
                else:
                    # gz：解压总量未知，按行数报进度（不给假百分比）
                    pct = min(99, int(n_lines / 100))
                    if pct != last_pct and pct % 10 == 0:
                        last_pct = pct
                        progress('pileup', 0.4 + 0.4 * pct / 100,
                                 f'统计位置计数（已读 {n_lines:,} 行）')
    return {'counts': counts, 'ins': ins, 'nreads': nreads,
            'qual_unavailable': qual_unavailable}


# ---------------------------------------------------------------- 共识
def call_consensus(rows, min_depth=MIN_DEPTH, min_freq=MIN_FREQ, ambig=AMBIG):
    """由位置计数调用共识序列。

    深度 < min_depth → 简并；最高频碱基频率 < min_freq → 简并（混合位点）。
    depth 计已通过质量过滤的碱基（A+C+G+T+N），不含缺失。
    """
    out = []
    called = 0
    amb = ambig or ''
    for row in rows:
        a, c, g, t, n, _gap = row
        depth = a + c + g + t + n
        if depth < min_depth:
            out.append(amb)
            continue
        acgt = (a, c, g, t)
        mx = max(acgt)
        if mx <= 0 or mx / depth < min_freq:
            out.append(amb)
            continue
        out.append(BASES[acgt.index(mx)])
        called += 1
    return ''.join(out), called


def coverage_stats(rows, min_depth=MIN_DEPTH):
    """覆盖统计。返回 dict。"""
    depths = [r[0] + r[1] + r[2] + r[3] + r[4] for r in rows]
    n = len(depths)
    if not n:
        return {'length': 0, 'mean_depth': 0, 'median_depth': 0, 'min_depth': 0,
                'max_depth': 0, 'covered_bp': 0, 'coverage_pct': 0.0,
                'depth_ok_pct': 0.0}
    cov = [d for d in depths if d > 0]
    ok = [d for d in depths if d >= min_depth]
    srt = sorted(depths)
    return {
        'length': n,
        'mean_depth': round(sum(depths) / n, 2),
        'median_depth': srt[n // 2],
        'min_depth': srt[0],
        'max_depth': srt[-1],
        'covered_bp': len(cov),
        'coverage_pct': round(len(cov) * 100.0 / n, 2),
        'depth_ok_pct': round(len(ok) * 100.0 / n, 2),
    }


def analyze_variants(ref_seq, cons, rows, min_minor_freq=MIN_MINOR_FREQ,
                     min_depth=MIN_DEPTH, ambig=AMBIG,
                     min_minor_count=MIN_MINOR_COUNT):
    """变异位点。

    type 可多选，分号分隔：
      SNV       共识 ≠ 参考（相对参考的置换，深度达标）
      iSNV      次等位基因频率 ≥ min_minor_freq 且 reads 数 ≥ min_minor_count
                （准种 / 混合感染信号）
      ambiguous 深度达标但共识仍为简并符（真实混合位点）
      lowcov    深度 < min_depth

    纯 lowcov（深度不足且无其他信号）**不入表**——覆盖缺口已由
    coverage.tsv 的 coverage_pct 表达，逐位点写出会产生上万行噪音。
    """
    out = []
    for i, row in enumerate(rows):
        a, c, g, t, n, gap = row
        depth = a + c + g + t + n
        ref_b = ref_seq[i] if i < len(ref_seq) else 'N'
        con_b = cons[i] if i < len(cons) else (ambig or 'N')
        cnt = {'A': a, 'C': c, 'G': g, 'T': t}
        rank = sorted(((v, k) for k, v in cnt.items()), reverse=True)
        major_n, major_b = rank[0]
        minor_n, minor_b = rank[1]
        major_f = round(major_n / depth, 4) if depth else 0.0
        minor_f = round(minor_n / depth, 4) if depth else 0.0

        types = []
        if depth < min_depth:
            types.append('lowcov')
            if con_b == (ambig or 'N') and depth > 0:
                types.append('ambiguous')
        else:
            if con_b == (ambig or 'N'):
                types.append('ambiguous')
            elif con_b in BASES and con_b != ref_b:
                types.append('SNV')
            if minor_n >= min_minor_count and minor_f >= min_minor_freq:
                types.append('iSNV')
        if not types or types == ['lowcov']:
            # 纯 lowcov（深度不足且无其他信号）不入表——覆盖缺口由
            # coverage.tsv 表达；否则 depth=0 的位点会把 variants.tsv
            # 灌满（实测 7235 行里 7154 行是 depth=0 的 lowcov）。
            continue
        out.append({
            'contig': '', 'pos': i + 1, 'ref_base': ref_b, 'cons_base': con_b,
            'depth': depth, 'A': a, 'C': c, 'G': g, 'T': t, 'N': n, 'gap': gap,
            'major_base': major_b, 'major_freq': major_f,
            'minor_base': minor_b, 'minor_freq': minor_f,
            'type': ';'.join(types),
        })
    return out


# ---------------------------------------------------------------- 主入口
def _probe_read_len(reads, max_records=2000):
    """头部采样 reads 的最大长度；读不到（空/异常）返回 None。

    reads 是 [(r1, r2), ...] 的配对列表（同 run_minimap2 的入参口径）；
    支持 FASTQ 与 FASTA（含 .gz）。取所有可读文件的最大读长——单端读长
    短不等于比对不了（另一端正常仍可映射），只有两端都短才判定为超短。
    只为「能否被 minibwa 比对」这一判断服务，不追求全量精确。
    """
    import gzip
    best = 0
    for pair in reads[:2]:
        for p in pair:
            if not p or not os.path.isfile(p):
                continue
            try:
                with open(p, 'rb') as probe:
                    gz = probe.read(2) == b'\x1f\x8b'
                if gz:
                    fh = gzip.open(p, 'rt', encoding='utf-8', errors='replace')
                else:
                    fh = open(p, 'r', encoding='utf-8', errors='replace')
                with fh:
                    first = ''
                    for ln in fh:
                        s = ln.strip()
                        if s:
                            first = s
                            break
                    m = 0
                    if first.startswith('>'):
                        # FASTA：累计每条记录长度（兼容多行折行）
                        cur, n = 0, 0
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
                    else:
                        for i, ln in enumerate(fh):
                            if i >= max_records * 4:
                                break
                            if i % 4 == 0:  # 首行'@'已消费，第 2/6/10… 行为序列
                                m = max(m, len(ln.strip()))
                    best = max(best, m)
            except OSError:
                continue
    return best or None


def consensus_and_variants(run_dir, ref_fasta, reads, out_subdir='consensus',
                           min_qual=MIN_BASE_QUALITY, min_depth=MIN_DEPTH,
                           min_freq=MIN_FREQ, ambig=AMBIG,
                           min_mapq=MIN_MAPQ, min_minor_freq=MIN_MINOR_FREQ,
                           min_minor_count=MIN_MINOR_COUNT,
                           min_cov_pct=MIN_COV_PCT, preset='sr',
                           threads=None, logger=None, progress=None,
                           keep_sam=False):
    """回贴 reads → 共识序列 + 变异谱。产物落在 run_dir/<out_subdir>/。

    ref_fasta  参考序列（验证后的候选 contig），每条参考独立统计
    reads      [(r1, r2_or_None), ...]（find_read_pairs / pair_read_files 的
               返回）；每对走 minibwa map 双位置参数双端比对
    返回摘要 dict（供任务引擎与前端展示）。
    """
    def _prog(stage, frac, msg):
        if progress:
            progress(stage, frac, msg)
        if logger:
            logger.log(msg)

    out_dir = os.path.join(run_dir, out_subdir)
    os.makedirs(out_dir, exist_ok=True)

    _prog('load', 0.05, '读取参考序列')
    refs = load_fasta(ref_fasta)
    if not refs:
        raise RuntimeError('参考序列为空：%s' % ref_fasta)
    total_len = sum(len(s) for s in refs.values())
    _prog('load', 0.10, '参考 %d 条，共 %d bp' % (len(refs), total_len))

    max_read_len = _probe_read_len(reads)
    short_reads = bool(max_read_len and max_read_len < MINIBWA_MIN_READ_LEN)
    if short_reads and logger:
        logger.log(
            f"⚠ 输入 reads 超短（头部采样最长 {max_read_len}bp）：minibwa 对 "
            f"<{MINIBWA_MIN_READ_LEN}bp reads 实测比对率为 0，本阶段覆盖/共识/"
            f"变异无法评估——相关 contig 记为 unknown（未知），不代表不存在；"
            f"如需共识请改用支持超短 reads 的比对器（如 bowtie 短读模式）",
            'WARN')

    sam = os.path.join(out_dir, 'aligned.sam')
    _prog('align', 0.15, 'minibwa 回贴 reads（preset=%s，双端配对）' % preset)
    run_minimap2(ref_fasta, reads, sam, threads=threads, preset=preset,
                 logger=logger)
    _prog('align', 0.40, '比对完成，统计位置计数')

    p = pileup(sam, refs, min_qual=min_qual, min_mapq=min_mapq,
               logger=logger, progress=progress)

    _prog('call', 0.85, '调用共识序列与变异位点')
    counts = p['counts']
    ins = p['ins']
    nreads = p['nreads']

    consensus = {}
    cov_rows = []
    all_variants = []
    for name, ref_seq in refs.items():
        rows = counts[name]
        seq, called = call_consensus(rows, min_depth, min_freq, ambig)
        consensus[name] = seq
        st = coverage_stats(rows, min_depth)
        st['contig'] = name
        st['mapped_reads'] = int(nreads.get(name, 0))
        st['called_pct'] = round(called * 100.0 / max(1, len(seq)), 2)
        vs = analyze_variants(ref_seq, seq, rows, min_minor_freq, min_depth,
                              ambig, min_minor_count)
        for v in vs:
            v['contig'] = name
        all_variants.extend(vs)
        st['n_variants'] = len(vs)
        # 精确匹配类型标签：'SNV' 是 'iSNV' 的子串，用 in 会误计
        st['n_snv'] = sum(1 for v in vs if 'SNV' in v['type'].split(';'))
        st['n_isnv'] = sum(1 for v in vs if 'iSNV' in v['type'].split(';'))
        if short_reads:
            # 短读比对不了 ≠ 病毒不存在：低覆盖结论不可采信，标为未知
            st['present'] = 'unknown'
        else:
            st['present'] = 'yes' if st['coverage_pct'] >= min_cov_pct else 'no'
        cov_rows.append(st)

    # ---- 产物 ----
    _prog('write', 0.92, '写出产物')
    safe_name = lambda s: re.sub(r'[^\w\-.]+', '_', s)[:120]

    _write_fasta(os.path.join(out_dir, 'consensus.fa'), consensus)

    for name, ref_seq in refs.items():
        tag = safe_name(name)
        with open(os.path.join(out_dir, '%s.poscounts.tsv' % tag), 'w',
                  encoding='utf-8', newline='') as f:
            f.write('Pos\tA\tC\tG\tT\tN\t-\tTotal\n')
            for i, row in enumerate(counts[name]):
                a, c, g, t, n, gap = row
                f.write('%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n'
                        % (i, a, c, g, t, n, gap, a + c + g + t + n + gap))
        with open(os.path.join(out_dir, '%s.inscounts.json' % tag), 'w',
                  encoding='utf-8') as f:
            json.dump({str(k): dict(v) for k, v in ins[name].items()}, f,
                      ensure_ascii=False)

    with open(os.path.join(out_dir, 'coverage.tsv'), 'w', encoding='utf-8',
              newline='') as f:
        w = __import__('csv').DictWriter(f, fieldnames=COVERAGE_FIELDS,
                                         delimiter='\t', restval='')
        w.writeheader()
        for r in cov_rows:
            w.writerow(r)

    with open(os.path.join(out_dir, 'variants.tsv'), 'w', encoding='utf-8',
              newline='') as f:
        w = __import__('csv').DictWriter(f, fieldnames=VARIANT_FIELDS,
                                         delimiter='\t', restval='')
        w.writeheader()
        for r in all_variants:
            w.writerow(r)

    if not keep_sam and os.path.isfile(sam):
        try:
            os.remove(sam)
        except OSError:
            pass

    present = [r['contig'] for r in cov_rows if r['present'] == 'yes']
    summary = {
        'n_refs': len(refs),
        'total_len': total_len,
        'n_mapped_reads': int(sum(nreads.values())),
        'n_variants': len(all_variants),
        'n_snv': sum(r['n_snv'] for r in cov_rows),
        'n_isnv': sum(r['n_isnv'] for r in cov_rows),
        'present': present,
        'low_support': [r['contig'] for r in cov_rows if r['present'] == 'no'],
        'unsupported_short_reads': bool(short_reads),
        'max_read_len_sampled': max_read_len,
        'unknown': [r['contig'] for r in cov_rows if r['present'] == 'unknown'],
        'min_cov_pct': min_cov_pct,
        # min_qual 是否真正生效：输入无碱基质量（QUAL 全为 '*'，minibwa 输出
        # 即如此）时为 false——下游/报告据此知道"这次没做质量过滤"，
        # 不能把低覆盖位点当作已按 Q20 把关的结果解读。
        'min_qual_effective': not bool(p.get('qual_unavailable')),
        'refs': cov_rows,
        'params': {'min_qual': min_qual, 'min_depth': min_depth,
                   'min_freq': min_freq, 'ambig': ambig,
                   'min_mapq': min_mapq, 'min_minor_freq': min_minor_freq,
                   'min_minor_count': min_minor_count, 'preset': preset},
        'out_dir': out_dir,
    }
    with open(os.path.join(out_dir, 'summary.json'), 'w', encoding='utf-8') as f:
        json.dump(summary, f, ensure_ascii=False, indent=1)

    _msg = ('完成：%d 条参考，%d 条 reads 比对上，变异位点 %d（SNV %d / iSNV %d），'
            '判定真实存在 %d 条'
            % (len(refs), summary['n_mapped_reads'], len(all_variants),
               summary['n_snv'], summary['n_isnv'], len(present)))
    if short_reads:
        _msg += ('；⚠ %d 条因超短 reads（最长 %dbp < %dbp）无法比对，'
                 '存在性记为 unknown'
                 % (len(summary['unknown']), max_read_len,
                    MINIBWA_MIN_READ_LEN))
    _prog('done', 1.0, _msg)
    return summary
