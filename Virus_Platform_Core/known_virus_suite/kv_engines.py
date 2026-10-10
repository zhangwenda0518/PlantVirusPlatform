#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
kv_engines.py — 引擎抽象层
=========================
定量（identify）引擎二选一（库与引擎正交：一个鉴定库可同时带两种索引）：
  salmon  —— 伪比对 EM 定量（quant.sf），--writeBam 出映射位点；多映射
             reads 按概率分摊，近缘参考多时更敏感（默认）。
  minibwa —— 真比对（mem）计数，每参考 Mapped_Reads 即真实比对数；
             CIGAR/NM/QUAL 为真。2026-09-13 恢复（原管线 traditional
             路径：bowtie2 → minibwa，覆盖度后端与命名完全一致）。

共识（consensus）仍固定由 kv_consensus 内部直接调 minibwa 自建索引重新
比对，不经本工厂（它要的是「对选出的单参考回贴」，与全库定量是两回事）。

索引布局（库 = 参考真相 + 各引擎派生索引，选库即选参考+索引整体）：
  <库>/salmon_k31/          salmon 索引
  <库>/minibwa/reference.*  minibwa 索引（.mbw/.l2b，前缀 minibwa/reference）

接口契约见 DESIGN.md 第三节。

设计来源: virome_analysis_pipeline/batch_virus_depth.py 的
build_index / process_sample_pseudo / process_sample_traditional
（复制后重构为引擎类，原文件不动）
"""

import gzip
import json
import os
import re
import shutil
import subprocess
import time
from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class AlignResult:
    """一次比对的统一返回结构"""
    sam_path: Path = None            # 比对记录（salmon 为排序后的映射 BAM）
    quant_path: Path = None          # 定量文件（salmon quant.sf）
    reads: dict = field(default_factory=dict)      # {accession: 定量 reads}
    refstats: dict = field(default_factory=dict)   # {accession: {Coverage(%), MeanDepth, Mapped_Reads, Covered_Bases}}
    mapped: int = 0                  # 实际比对上的记录数/片段数
    total: int = 0                   # 输入总片段数
    elapsed: float = 0.0
    has_positions: bool = False      # 是否含真实比对位点（决定共识段可用性）


# 小RNA（siRNA/miRNA）reads 的长度上界：21-24nt 建库为主，个别 27-30nt。
# 超过此长度的按常规 reads 处理（salmon k=31）。
SRNA_MAX_READ_LEN = 32
# 小RNA 场景的 salmon 索引 k：实测 250K×21-24nt 对 8465 参考库，
# k=31 映射率 0%，k=19 8.7%、k=15 10.9%（PVB 定量与 BWA 真值差 3-6%）。
SRNA_SALMON_K = 15


def sniff_read_len_stats(path, max_records=4000):
    """头部采样读长统计；FASTA/读失败返回 None。

    为什么不用「最长 read」判定小RNA：一条 40nt 的接头残留/污染 read 就会让
    整样本按常规 RNA 处理（k=31 对 22nt reads 映射率为 0 → 静默零命中）。
    改用 p90 + 短读占比判定（见 looks_like_srna）。
    """
    try:
        with open(path, 'rb') as probe:
            gz = probe.read(2) == b'\x1f\x8b'
        if gz:
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
                return None                     # FASTA：无 4 行周期，不嗅探
            lens = []
            for i, ln in enumerate(fh):
                if i >= max_records * 4:
                    break
                if i % 4 == 0:                  # 首行'@'已消费，第2/6/10…行为序列
                    s = ln.strip()
                    if s and s.startswith('@'):
                        continue                # 行序错位（空行/折叠序列）跳过
                    lens.append(len(s))
        if not lens:
            return None
        lens.sort()
        n = len(lens)
        return {'n': n, 'max': lens[-1], 'median': lens[n // 2],
                'p90': lens[min(n - 1, int(n * 0.9))],
                'short_frac': sum(1 for x in lens
                                  if x < SRNA_MAX_READ_LEN) / n}
    except OSError:
        return None


def looks_like_srna(stats, min_reads=20):
    """读长统计 → 是否小RNA 样本（p90 < 32bp；样本量太小时不判定）。"""
    return bool(stats and stats.get('n', 0) >= min_reads
                and stats.get('p90', 999) < SRNA_MAX_READ_LEN)


def sniff_max_read_len(path, max_records=4000):
    """嗅探 FASTQ 前 max_records 条 reads 的最大长度；FASTA/读失败返回 None。"""
    st = sniff_read_len_stats(path, max_records=max_records)
    return st['max'] if st else None


def is_srna_sample(sample):
    """样本是否小RNA reads（r1/r2 头部采样 p90 < SRNA_MAX_READ_LEN）。"""
    for p in (sample.get('r1'), sample.get('r2')):
        if p and looks_like_srna(sniff_read_len_stats(p)):
            return True
    return False


# ── 索引↔参考一致性指纹 ──────────────────────────────────────
# 复用判据原来只看"索引文件是否存在"，换库/更新参考/共享索引目录都会
# 静默用旧索引定量（实测：0 字节空参考也能"检出"旧库病毒）。建索引时
# 落 ref_fingerprint.json，复用时校验。
def ref_fingerprint(ref_fasta):
    """参考 FASTA 指纹；>256MB 只哈希头尾各 1MiB（避免每次跑等数十秒）。"""
    import hashlib
    p = Path(ref_fasta)
    st = p.stat()
    out = {'size': int(st.st_size), 'mtime': int(st.st_mtime_ns)}
    h = hashlib.sha1()
    with open(p, 'rb') as f:
        if st.st_size <= 256 * 1024 * 1024:
            for chunk in iter(lambda: f.read(1 << 20), b''):
                h.update(chunk)
            out['sha1'] = h.hexdigest()
        else:
            head = f.read(1 << 20)
            f.seek(max(0, st.st_size - (1 << 20)))
            h.update(head)
            h.update(f.read(1 << 20))
            out['sha1_head_tail'] = h.hexdigest()
    return out


def _write_index_fingerprint(fp_path, fp, logger=None):
    try:
        with open(fp_path, 'w', encoding='utf-8') as f:
            json.dump(fp, f, ensure_ascii=False, indent=1)
    except OSError as e:
        if logger:
            logger.warning(f"[index] 参考指纹写入失败 [{type(e).__name__}]: {e}"
                           f"——下次运行会重新校验（可能触发多余重建）")


def index_ref_ok(idx_dir, ref_fasta, probe_file, logger=None):
    """索引目录的参考指纹是否与当前参考一致；不一致调用方应重建。

    缺指纹（旧版构建）时按 mtime 启发式：索引比参考新 → 视为一致并补记；
    参考比索引新 → 视为不一致（保守：宁可重建一次）。
    """
    fp_path = Path(idx_dir) / 'ref_fingerprint.json'
    try:
        with open(fp_path, encoding='utf-8') as f:
            prev = json.load(f)
    except (OSError, ValueError):
        prev = None
    try:
        cur = ref_fingerprint(ref_fasta)
    except OSError as e:
        if logger:
            logger.warning(f"[index] 参考指纹计算失败 [{type(e).__name__}]: {e}"
                           f"——按沿用处理")
        return True
    if prev is None:
        try:
            if (os.path.getmtime(probe_file) + 5
                    >= os.path.getmtime(ref_fasta)):
                if logger:
                    logger.info("[index] 索引无参考指纹（旧版构建），"
                                "按 mtime 判定一致并补记")
                _write_index_fingerprint(fp_path, cur, logger=logger)
                return True
        except OSError:
            pass
        return False
    if prev.get('size') == cur.get('size') and \
            (prev.get('sha1') == cur.get('sha1')
             or prev.get('sha1_head_tail') == cur.get('sha1_head_tail')):
        return True
    return False


class VirusEngine:
    """引擎基类。子类实现 build_index / align。"""
    name = 'base'
    supports_consensus = False

    def __init__(self, exe, threads=8, logger=None, min_mapq=10):
        self.exe = exe
        self.threads = threads
        self.logger = logger
        self.min_mapq = min_mapq

    def log(self, msg, level='info'):
        if self.logger:
            getattr(self.logger, level)(f"[{self.name}] {msg}")

    def build_index(self, ref_fasta, index_dir, threads=None, kmer=31):
        raise NotImplementedError

    def align(self, index_path, sample, out_dir, threads=None):
        raise NotImplementedError


# ══════════════════════════════════════════════════════════
# salmon —— 定量引擎（伪比对定量 + --writeBam 映射位点）
# ══════════════════════════════════════════════════════════
class SalmonEngine(VirusEngine):
    name = 'salmon'
    supports_consensus = False   # 映射 CIGAR 是伪的（无缺口）、无 NM/QUAL——
                                 # 共识/变异的真实位点仍由共识段内部的 minibwa 提供

    def __init__(self, exe, threads=8, logger=None, min_mapq=10,
                 samtools=None, pandepth=None, coverage_tool='pandepth'):
        super().__init__(exe, threads=threads, logger=logger, min_mapq=min_mapq)
        # 覆盖度后端：默认 pandepth（与原 minibwa 路径同一口径：view -F 0x04
        # 建 BAM + pandepth -a 无 -q）。两者齐全才出真实覆盖度/深度与比对 BAM，
        # 缺失时退回「reads×读长」估算（引擎仍可用）。
        self.samtools = samtools
        self.pandepth = pandepth
        self.coverage_tool = (coverage_tool or 'pandepth').lower()

    def build_index(self, ref_fasta, index_dir, threads=None, kmer=31):
        """索引按 k 分目录（salmon_k31 / salmon_k15…），同库可并存。

        salmon_k31 是库有效性探针与常规数据的默认索引（kv_stage/build.py
        都按这个名字探测），不能动；小RNA 运行时按需在旁路建 salmon_k15。
        """
        threads = threads or self.threads
        index_dir = Path(index_dir)
        idx_path = index_dir / f'salmon_k{kmer}'
        if idx_path.is_dir() and (idx_path / 'info.json').exists():
            if index_ref_ok(idx_path, ref_fasta, idx_path / 'info.json',
                            logger=self.logger):
                self.log(f"复用已有索引 -> {idx_path}")
                # 旧索引兼容：salmon 默认折叠完全重复的参考（931 条重复被并掉，
                # 被折叠参考永远拿不到 reads）。新索引一律 --keepDuplicates；
                # 旧索引无法自动补救，这里只显式提示重建（不自动重建，避免
                # 每次跑都白等一轮建索引）。
                if self._index_collapsed_duplicates(idx_path):
                    self.log(
                        "旧索引折叠过完全重复参考（duplicate_clusters.tsv / "
                        "info.json:keep_duplicates=false）：被折叠的参考将永远"
                        "拿不到 reads，且无告警。请改名/删除该索引目录后重建"
                        "（本次不自动重建）", level='warning')
                return idx_path
            # 参考已变化（换库/更新参考）：旧索引定量结果是错的，必须重建。
            # 旧目录改名保留（不直接删），构建失败时还能人工回退。
            stale = idx_path.with_name(idx_path.name + '.stale')
            shutil.rmtree(stale, ignore_errors=True)
            try:
                idx_path.rename(stale)
            except OSError:
                pass
            self.log(f"索引与当前参考不一致（换库/更新参考后未重建），"
                     f"重建索引 -> {idx_path}", level='warning')

        index_dir.mkdir(parents=True, exist_ok=True)
        self.log(f"构建索引 (k={kmer}) -> {idx_path}")
        t0 = time.time()
        # --keepDuplicates（salmon 2.7 原生支持）：保留全部完全重复参考。
        # salmon 默认折叠完全重复序列（本平台库实测 931 条被折进 192 个
        # 保留簇），被折叠的 accession 在 quant.sf 里恒为 0 —— 换库/物种名
        # 不同的重复参考（如不同病毒名共用同一条序列）会被静默漏检。
        cmd = [self.exe, 'index', '-t', str(ref_fasta), '-i', str(idx_path),
               '-k', str(kmer), '-p', str(threads), '--keepDuplicates']
        with open(index_dir / f'salmon_index_k{kmer}.log', 'w',
                  encoding='utf-8') as fl:
            subprocess.run(cmd, stdout=fl, stderr=subprocess.STDOUT, check=True, timeout=7200)
        self.log(f"索引完成 {time.time()-t0:.1f}s"
                 f"（--keepDuplicates：保留全部重复参考，不折叠）")
        try:
            _write_index_fingerprint(idx_path / 'ref_fingerprint.json',
                                     ref_fingerprint(ref_fasta),
                                     logger=self.logger)
        except OSError as e:
            self.log(f"参考指纹写入失败 [{type(e).__name__}]: {e}",
                     level='warning')
        return idx_path

    @staticmethod
    def _index_collapsed_duplicates(idx_path):
        """已有索引是否由「折叠重复参考」的旧参数建成（即需重建才能保留重复）。

        判据优先读 salmon 的 info.json:keep_duplicates（2.x 明确记录）；
        更旧的索引没有该字段时，以 duplicate_clusters.tsv 是否存在为准
        （有该文件即说明建索引时发生了折叠）。
        """
        info = Path(idx_path) / 'info.json'
        try:
            with open(info, encoding='utf-8') as f:
                kd = json.load(f).get('keep_duplicates')
            if kd is not None:
                return kd is False
        except (OSError, ValueError):
            pass
        return (Path(idx_path) / 'duplicate_clusters.tsv').exists()

    def align(self, index_path, sample, out_dir, threads=None):
        threads = threads or self.threads
        out_dir = Path(out_dir)
        sname = sample['name']
        quant_dir = out_dir / 'quant' / sname
        quant_dir.mkdir(parents=True, exist_ok=True)

        # --writeBam：映射记录直接落 BAM（SEQ 为真、CIGAR 伪无缺口、无 QUAL）。
        # 只含 mapped 记录（salmon 只写 per-mapping 记录，unmapped 不落盘）。
        bam_raw = quant_dir / 'mappings.bam'
        cmd = [self.exe, 'quant', '-i', str(index_path), '-p', str(threads),
               '-l', 'A', '-o', str(quant_dir), '--writeBam', str(bam_raw)]
        if sample.get('r2'):
            cmd += ['-1', sample['r1'], '-2', sample['r2']]
        else:
            cmd += ['-r', sample['r1']]

        t0 = time.time()
        logf = out_dir / 'logs' / f'align_{sname}.log'
        logf.parent.mkdir(parents=True, exist_ok=True)
        with open(logf, 'w', encoding='utf-8') as fl:
            try:
                subprocess.run(cmd, stdout=fl, stderr=subprocess.STDOUT,
                               check=True)
            except subprocess.CalledProcessError as e:
                raise RuntimeError(
                    f"salmon quant 失败（exit {e.returncode}，详见 {logf}）。"
                    f"若输入为小RNA reads，请确认所用索引 k < 读长"
                    f"（超短 reads 需 k={SRNA_SALMON_K} 小索引）") from e
        elapsed = time.time() - t0

        res = AlignResult(quant_path=quant_dir / 'quant.sf', elapsed=elapsed,
                          has_positions=True)

        # quant.sf → EM 定量（salmon 的核心口径：多映射 reads 经 EM 分摊）
        qf = res.quant_path
        if qf.exists():
            with open(qf, 'r', encoding='utf-8') as f:
                header = f.readline().rstrip('\n').split('\t')
                ci = {h: i for i, h in enumerate(header)}
                for line in f:
                    p = line.rstrip('\n').split('\t')
                    if len(p) <= ci.get('NumReads', 1):
                        continue
                    acc = p[ci['Name']]
                    n = float(p[ci['NumReads']])
                    if n > 0:
                        res.reads[acc] = n

        # 覆盖度/深度/逐位点深度表：与原 minibwa 路径同一后端（pandepth），
        # 产物命名保持 <样本>_pandepth.* —— kv_plot 的深度表查找
        # （<样本>_pandepth.pandepth.SiteDepth.gz）与管道②b 抽病毒 reads 的
        # BAM 拾取（align/<样本>_pandepth.sorted.bam）都依赖这套名字。
        if not (self.samtools and self.pandepth):
            self.log('samtools/pandepth 不可用：不产出比对 BAM，覆盖度/深度退到'
                     '「reads×读长」估算', level='warning')
            res.total, res.mapped = self._meta_counts(quant_dir, logf, res.reads)
            res.has_positions = False
            try:
                bam_raw.unlink()     # 不会再用；GB 级原始 BAM 不留盘
            except OSError:
                pass
            return res

        try:
            backend = BamCoverageBackend(self.samtools, self.pandepth,
                                         threads=max(1, min(threads, 8)),
                                         logger=self.logger)
            align_dir = out_dir / 'align'
            name = f'{sname}_pandepth'
            sorted_bam = backend.to_bam(bam_raw, align_dir, name=name,
                                        strict_orig=True)
            if self.coverage_tool == 'pandepth':
                refstats = backend.coverage_pandepth(sorted_bam, align_dir,
                                                     name=name)
            else:
                refstats = backend.coverage_samtools(sorted_bam, align_dir,
                                                     name=name)
            mapped_reads = backend.count_mapped_by_ref(sorted_bam)
        finally:
            # 原始未排序映射 BAM 无论成败都不留盘：成功则已转换，
            # 失败则重跑本段会重新产出，残留在盘只会吃掉几十 GB
            try:
                bam_raw.unlink()
            except OSError:
                pass

        for acc, st in refstats.items():
            st['Mapped_Reads'] = mapped_reads.get(acc, 0)

        res.sam_path = sorted_bam
        res.refstats = refstats
        # 比对率口径修正：idxstats 的每参考计数含 secondary/supplementary
        # （多映射记录），求和的 mapped 会超过 -F 0x4 -F 0x100 的 total，
        # 比对率虚高甚至 >100%。统一改读 salmon 自己的 meta_info.json
        # （num_processed / num_mapped），与 BAM 统计无关。
        res.total, res.mapped = self._meta_counts(quant_dir, logf, res.reads)
        return res

    def _meta_counts(self, quant_dir, logf, reads):
        """salmon 的输入/比对片段数：优先 aux_info/meta_info.json。

        meta_info.json 是 salmon 自己写的权威统计（num_processed 输入片段、
        num_mapped 比对上片段），不依赖 BAM 里记录的 flag 口径；
        缺失/不可读时退回日志解析（_parse_counts）。
        """
        try:
            with open(Path(quant_dir) / 'aux_info' / 'meta_info.json',
                      encoding='utf-8') as f:
                meta = json.load(f)
            total = int(meta.get('num_processed') or 0)
            mapped = int(meta.get('num_mapped') or 0)
            if total > 0:
                return total, mapped
            self.log('meta_info.json 的 num_processed=0，比对统计退回日志解析',
                     level='warning')
        except (OSError, ValueError) as e:
            self.log(f'读取 meta_info.json 失败（{e}），比对统计退回日志解析',
                     level='warning')
        return self._parse_counts(logf, reads)

    def _parse_counts(self, logf, reads):
        """salmon 的比对统计口径：num_processed / num_mapped

        返回 (total, mapped)。日志里解析不到这两个数时退到 reads 总和，
        但记 warning——否则日志会显示成 100% 比对率，误导排查。
        （仅 samtools/pandepth 缺失、无 BAM 可统计时的兜底路径。）
        """
        total = mapped = 0
        try:
            txt = Path(logf).read_text(encoding='utf-8', errors='replace')
            m = re.search(r'(\d+)\s+fragments?\s+from\s+intermediate\s+RAD,\s+(\d+)\s+quantified', txt)
            if m:
                total, mapped = int(m.group(1)), int(m.group(2))
            else:
                m = re.search(r'mapped\s+(\d+)\s*/\s*(\d+)\s+fragments', txt)
                if m:
                    mapped, total = int(m.group(1)), int(m.group(2))
        except Exception as e:  # noqa: BLE001
            self.log(f'salmon 日志解析异常（{e}），比对统计改用 reads 总和',
                     level='warning')
        if total == 0:
            s = int(sum(reads.values()))
            self.log(f'未从 {Path(logf).name} 解析到比对统计，'
                     f'total/mapped 退到 reads 总和 {s}', level='warning')
            total = s or 1
            mapped = s
        return total, mapped


# ══════════════════════════════════════════════════════════
# minibwa —— 真比对定量引擎（原管线 traditional 路径，2026-09-13 恢复）
# ══════════════════════════════════════════════════════════
class BwaEngine(VirusEngine):
    """minibwa mem 真比对 → 排序 BAM → pandepth 覆盖度/深度。

    与 SalmonEngine 的差异只在「reads 怎么来」：
      salmon: quant.sf 的 EM NumReads（多映射分摊，非整数）
      minibwa: samtools idxstats 的每参考比对记录数（真实计数）
    覆盖度/深度/BAM 命名与 salmon 完全同一后端（<样本>_pandepth.*），
    下游 kv_plot 深度表与 ②b 抽病毒 reads 无需区分引擎。
    """
    name = 'minibwa'
    supports_consensus = False   # 共识段固定自建单参考索引重新比对（kv_consensus），
                                 # 不消费本引擎的 BAM——维持"定量/共识两条链路分离"

    def __init__(self, exe, threads=8, logger=None, min_mapq=10,
                 samtools=None, pandepth=None, coverage_tool='pandepth'):
        super().__init__(exe, threads=threads, logger=logger, min_mapq=min_mapq)
        if not samtools:
            raise ValueError('minibwa 定量需要 samtools（建排序 BAM 与计数）；'
                             '未找到 samtools 可执行文件')
        self.samtools = samtools
        self.pandepth = pandepth
        # pandepth 缺件时自动降级 samtools 后端（与 kv_config 口径一致），
        # 不像 salmon 那样能退到"无 BAM 估算"——真比对引擎没有 BAM 就没有定量
        if (coverage_tool or 'pandepth').lower() == 'pandepth' and not pandepth:
            self.log('pandepth 不可用，覆盖度/深度后端降级为 samtools coverage',
                     level='warning')
            coverage_tool = 'samtools'
        self.coverage_tool = (coverage_tool or 'pandepth').lower()

    def build_index(self, ref_fasta, index_dir, threads=None, kmer=31):
        """索引布局：<index_dir>/minibwa/reference.{mbw,l2b}。

        与库布局约定一致（build.py 建库、kv_stage.index_dir_for、
        _kv_index_probe 都以 reference.mbw 为存在判据）；复用判据与
        kv_consensus.build_minibwa_index 相同（.mbw 存在即复用）。
        kmer 参数仅为与 SalmonEngine.build_index 签名对齐，minibwa 的
        播种长度在 map/mem 时才起作用，索引本身与 k 无关。
        """
        index_dir = Path(index_dir)
        sub = index_dir / 'minibwa'
        prefix = sub / 'reference'
        if Path(str(prefix) + '.mbw').exists():
            if index_ref_ok(sub, ref_fasta, Path(str(prefix) + '.mbw'),
                            logger=self.logger):
                self.log(f"复用已有索引 -> {prefix}")
                return sub
            self.log(f"索引与当前参考不一致（换库/更新参考后未重建），"
                     f"重建索引 -> {prefix}", level='warning')
        sub.mkdir(parents=True, exist_ok=True)
        self.log(f"构建 minibwa 索引 -> {prefix}")
        t0 = time.time()
        cmd = [self.exe, 'index', str(ref_fasta), str(prefix)]
        with open(sub / 'index.log', 'w', encoding='utf-8') as fl:
            subprocess.run(cmd, stdout=fl, stderr=subprocess.STDOUT, check=True, timeout=7200)
        self.log(f"索引完成 {time.time()-t0:.1f}s")
        try:
            _write_index_fingerprint(sub / 'ref_fingerprint.json',
                                     ref_fingerprint(ref_fasta),
                                     logger=self.logger)
        except OSError as e:
            self.log(f"参考指纹写入失败 [{type(e).__name__}]: {e}",
                     level='warning')
        return sub

    def align(self, index_path, sample, out_dir, threads=None):
        threads = threads or self.threads
        out_dir = Path(out_dir)
        sname = sample['name']
        # 小RNA 硬拒绝：minibwa（BWA 0.7 seed-and-extend）对 21-24nt reads
        # 实测 0/250,000 映射（任何 seed 参数），静默跑完只会产出全零定量。
        _st = sniff_read_len_stats(sample.get('r1') or '')
        if looks_like_srna(_st):
            raise ValueError(
                f"[{sname}] minibwa 不支持小RNA reads（头部采样中位 "
                f"{_st['median']}bp / p90 {_st['p90']}bp，"
                f"实测 21-24nt 无法比对、定量全为零）——"
                f"小RNA 请改用 salmon 引擎（k={SRNA_SALMON_K} 索引自动适配）")
        prefix = Path(index_path) / 'reference'
        if not Path(str(prefix) + '.mbw').exists():
            raise FileNotFoundError(f"minibwa 索引缺失: {prefix}.mbw")

        # minibwa map → SAM(流) → samtools view -b 落未排序 BAM。
        # Python 桥接不经 shell（同 kv_consensus.align_sample 的管道写法）；
        # 不过滤记录，-F 0x04 的"只剔 unmapped"语义留给 BamCoverageBackend，
        # 两条引擎路径共用同一套 BAM 口径。
        align_dir = out_dir / 'align'
        align_dir.mkdir(parents=True, exist_ok=True)
        bam_raw = align_dir / f'{sname}.minibwa.bam'
        map_cmd = [self.exe, 'map', f'-t{threads}', '-o', '-',
                   str(prefix), str(sample['r1'])]
        if sample.get('r2'):
            map_cmd.append(str(sample['r2']))
        t0 = time.time()
        logf = out_dir / 'logs' / f'align_{sname}.log'
        logf.parent.mkdir(parents=True, exist_ok=True)
        # stderr 必须落文件而非 PIPE：大数据集比对时 stderr 产出可观，
        # PIPE 无人并发读取会在 64KB 缓冲写满后把 minibwa 整个写死
        # （表现：BAM 停止增长、minibwa/samtools 双进程 0% CPU 永久挂起）。
        with open(logf, 'ab') as ef:
            p_map = subprocess.Popen(map_cmd, stdout=subprocess.PIPE,
                                     stderr=ef)
            # 注册进任务取消表：任务硬停止/超时时 TaskManager 可杀到
            # 管道两进程（此前仅靠外层 6h 兜底——2026-10-04 审查 P3）
            from Virus_Platform_Core.utils import task_register_proc
            task_register_proc(p_map)
            p_bam = subprocess.Popen([self.samtools, 'view', '-b', '-o',
                                      str(bam_raw), '-'],
                                     stdin=p_map.stdout,
                                     stdout=subprocess.DEVNULL,
                                     stderr=ef)
            task_register_proc(p_bam)
            p_map.stdout.close()
            rc_bam = p_bam.wait()
            rc_map = p_map.wait()
        # minibwa/samtools 的 CRT 消息按系统 ANSI 代码页（中文 Windows = GBK）
        # 打印：捕获文件就地转成 UTF-8，与平台其余日志编码统一。
        try:
            raw = logf.read_bytes()
            if raw:
                from Virus_Platform_Core.utils import decode_output
                logf.write_bytes(decode_output(raw).encode('utf-8'))
        except (OSError, UnicodeError):
            pass
        if rc_map != 0 or rc_bam != 0:
            raise subprocess.CalledProcessError(
                max(rc_map, rc_bam), 'minibwa map | samtools view -b',
                stderr=f'详见 {logf}'.encode())
        with open(logf, 'a', encoding='utf-8') as fl:
            import shlex as _shlex
            fl.write(_shlex.join(str(x) for x in map_cmd)
                     + ' | samtools view -b')
        elapsed = time.time() - t0

        # 输入总量在过滤前数：主记录数（SE=read 数，PE=2×对数），
        # 才能得到真实的比对率（过滤后再数就恒为 100% 了）
        res = AlignResult(elapsed=elapsed, has_positions=True)
        res.total = self._count_primary(bam_raw)

        backend = BamCoverageBackend(self.samtools, self.pandepth,
                                     threads=max(1, min(threads, 8)),
                                     logger=self.logger)
        try:
            name = f'{sname}_pandepth'      # 命名与 salmon 路径一致（下游依赖）
            sorted_bam = backend.to_bam(bam_raw, align_dir, name=name,
                                        strict_orig=True)
            if self.coverage_tool == 'pandepth':
                refstats = backend.coverage_pandepth(sorted_bam, align_dir,
                                                     name=name)
            else:
                refstats = backend.coverage_samtools(sorted_bam, align_dir,
                                                     name=name)
            mapped_reads = backend.count_mapped_by_ref(sorted_bam)
        finally:
            # 原始未排序 BAM 无论成败不留盘（同 salmon 路径的口径）
            try:
                bam_raw.unlink()
            except OSError:
                pass

        for acc, st in refstats.items():
            st['Mapped_Reads'] = mapped_reads.get(acc, 0)

        # ANI（NM tag 口径，照原管线 batch_virus_depth38.py:_compute_ani_pi_worker）。
        # 与覆盖度同源：都读同一份 sorted BAM，不额外比对。
        # 失败不阻断鉴定（置 None，由分流逻辑按「未测定」处理）。
        try:
            ani_map = backend.ani_by_ref(sorted_bam)
        except Exception as e:  # noqa: BLE001
            self.log(f"ANI 计算失败（置空）: {e}", level='warning')
            ani_map = {}
        for acc, st in refstats.items():
            st['Avg_Read_ANI'] = ani_map.get(acc)

        res.sam_path = sorted_bam
        res.refstats = refstats
        # 真比对的"定量"就是每参考的真实比对计数（无 EM 分摊）。
        # 注意：idxstats 的每参考计数**含 secondary/supplementary 多映射记录**
        # （实测：一条 secondary + 一条 supplementary 都被计入参考行），
        # 与 identify 表的 Uniq_Reads 同口径保留；但样本级比对率不能拿它
        # 求和——会超过主记录总数（>100%）。统一按主记录口径另算：
        #   total  = -F 0x900（全部主记录，含未比对，见 _count_primary）
        #   mapped = -F 0x904（已比对主记录）
        res.reads = {acc: float(n) for acc, n in mapped_reads.items()}
        res.mapped = self._count_mapped_primary(sorted_bam)
        return res

    def _count_primary(self, bam):
        """主记录总数（-F 0x900：去 secondary/supplementary，保留未比对主记录）"""
        r = subprocess.run(
            [self.samtools, 'view', '-c', '-F', '0x900', str(bam)],
            capture_output=True, text=True, encoding='utf-8',
            errors='replace')
        return int((r.stdout or '0').strip() or 0)

    def _count_mapped_primary(self, bam):
        """已比对主记录数（-F 0x904：去 secondary/supplementary/unmapped）"""
        r = subprocess.run(
            [self.samtools, 'view', '-c', '-F', '0x904', str(bam)],
            capture_output=True, text=True, encoding='utf-8',
            errors='replace')
        return int((r.stdout or '0').strip() or 0)


# ══════════════════════════════════════════════════════════
# 覆盖度/深度：只走外部工具，严格照原管线
# ══════════════════════════════════════════════════════════
# 说明：模块内曾有一个自研的 SAM/CIGAR 解析器（parse_sam_coverage），已删除。
# 原管线（virome_analysis_pipeline/batch_virus_depth.py）的覆盖度与深度
# 只有一条来源：pandepth 的 <prefix>.chr.stat.gz（第 5/6 列），
# 全文没有任何 CIGAR 运算。自研实现属于平行口径，即使对拍一致也制造
# 静默分歧风险（尤其 D/N 处理、secondary/supplementary 取舍这类边界）。
# 工具缺失时直接 fail-fast，不降级。


# ══════════════════════════════════════════════════════════
# 外部覆盖度后端：pandepth（原管线口径）
# ══════════════════════════════════════════════════════════
class BamCoverageBackend:
    """
    SAM -> BAM -> sort -> index，然后交给 samtools coverage 或 pandepth 统计。

    ★ pandepth 路径严格照抄原管线（virome_analysis_pipeline/batch_virus_depth.py
      L372-404 pseudo 分支 / L493-503 traditional 分支）：
        samtools view -b -F 0x04   （只去 unmapped，不过滤 MAPQ）
        samtools sort -@ N
        pysam.index(bam)           （本地用 samtools index 等价替代）
        pandepth -a -i bam -o prefix -t N    ← **无 -q**（min MAPQ = 0）
        读 <prefix>.chr.stat.gz，第 5 列 Coverage(%)、第 6 列 MeanDepth

    口径要点：
      * 原管线不传 -q，所以 MAPQ 0 起算；本模块的 min_mapq 只作用于
        builtin / samtools 后端，pandepth 后端照原管线为 0。
      * pandepth -x 默认 1796，排除 0x4|0x100|0x200|0x800。
      * -a 语义是 "output all the site depth"，额外产逐位点文件，
        不改变统计口径。实测 truth30 / q06 均 5.5s。

    三条口径已在 truth30/truth3/q06 上对拍一致（差异仅小数位）：
      samtools coverage / pandepth / 内置 parse_sam_coverage
    """

    def __init__(self, samtools_exe, pandepth_exe=None, threads=4, logger=None):
        self.samtools = samtools_exe
        self.pandepth = pandepth_exe
        self.threads = threads
        self.logger = logger

    def log(self, msg, level='info'):
        if self.logger:
            getattr(self.logger, level)(f"[bam] {msg}")

    def to_bam(self, sam_path, out_dir, name='coverage', min_mapq=10,
               strict_orig=True):
        """SAM -> 排序后的 BAM + 索引。返回 bam_path

        strict_orig=True（默认）照原管线：samtools view -b -F 0x04
          只剔除 unmapped，不做 MAPQ 过滤。
        strict_orig=False 启用 min_mapq 过滤。
        """
        sam_path = Path(sam_path)
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        bam = out_dir / f'{name}.bam'
        logf = out_dir / f'{name}_samtools.log'

        if strict_orig:
            cmd = [self.samtools, 'view', '-b', '-F', '0x04',
                   '-o', str(bam), str(sam_path)]
        else:
            cmd = [self.samtools, 'view', '-b', '-q', str(min_mapq),
                   '-o', str(bam), str(sam_path)]
        with open(logf, 'w', encoding='utf-8') as fl:
            subprocess.run(cmd, stdout=fl, stderr=subprocess.STDOUT, check=True, timeout=7200)

        sorted_bam = out_dir / f'{name}.sorted.bam'
        # -m 512M：samtools sort 默认每线程 768MB 提交内存，病毒 BAM 通常
        # 很小；系统内存紧张时默认值会分配失败（couldn't allocate memory
        # for bam_mem），显式限幅更稳（超限部分自动走临时文件）。
        cmd = [self.samtools, 'sort', '-m', '256M', '-@', str(self.threads),
               '-o', str(sorted_bam), str(bam)]
        with open(logf, 'a', encoding='utf-8') as fl:
            subprocess.run(cmd, stdout=fl, stderr=subprocess.STDOUT, check=True, timeout=7200)

        cmd = [self.samtools, 'index', str(sorted_bam)]
        with open(logf, 'a', encoding='utf-8') as fl:
            subprocess.run(cmd, stdout=fl, stderr=subprocess.STDOUT, check=True, timeout=7200)

        try:
            bam.unlink()
        except OSError:
            pass
        return sorted_bam

    def count_mapped_by_ref(self, sorted_bam):
        """每个参考的比对记录数（samtools idxstats）

        ⚠️ 口径：idxstats 第 3 列**含 secondary/supplementary 多映射记录**
        （实测：secondary、supplementary 各一条时该参考计数为 3，主记录仅 1）——
        下游 identify 表的 Uniq_Reads 即此口径，用于「该参考拿到多少条 reads」；
        但样本级比对率（mapped/total）不可用它求和，BwaEngine.align 另按
        -F 0x904 的主记录口径统计。

        pandepth 的 .chr.stat.gz 只有 Coverage(%)/MeanDepth，没有读数；
        原管线在 traditional 分支用 samtools view -c 拿总记录数，
        这里按参考细分，同一工具、同一口径。

        实现：用 samtools idxstats 一次调用拿全部参考的映射数（第 3 列）。
        旧实现对 8464 条参考逐条启动 samtools（一次进程几十毫秒），
        单样本要多花几分钟；idxstats 是毫秒级。
        """
        r = subprocess.run([self.samtools, 'idxstats', str(sorted_bam)],
                           capture_output=True, text=True,
                           encoding='utf-8', errors='replace')
        counts = {}
        for line in r.stdout.splitlines():
            p = line.split('\t')
            if len(p) >= 4 and p[0] != '*':
                counts[p[0]] = int(p[2])
        return counts

    def count_records(self, sorted_bam):
        """比对上的记录总数（去 unmapped / secondary）"""
        r = subprocess.run(
            [self.samtools, 'view', '-F', '0x4', '-F', '0x100', '-c',
             str(sorted_bam)],
            capture_output=True, text=True,
            encoding='utf-8', errors='replace')
        return int((r.stdout or '0').strip() or 0)

    def ani_by_ref(self, sorted_bam, max_reads=10000):
        """逐参考的平均核酸一致性 ANI（NM tag 口径）

        口径严格照原管线 batch_virus_depth38.py:_compute_ani_pi_worker：
          ANI_read = (query_alignment_length - NM) / query_alignment_length
          对每个参考最多取 max_reads 条比对记录，取算术平均后 ×100
        区别只在实现方式：原管线用 pysam 逐记录读 NM，这里用
        `samtools view` 的 SAM 文本逐行解析（避免新增 pysam 依赖），
        两者的 query_alignment_length 与 NM 定义完全一致。

        采样策略：原管线用 step = max(1, total//10000) 等距抽样。这里从盘上
        流式读，按参考累计到 max_reads 就不再计入（前 N 条）。两种都只是
        采样，同一条 BAM 上结果差异在小数第二位上随机波动，不影响阈值判定；
        阈值 95% 与 99% ANI 的区分远大于此。

        返回 {reference: ani_percent(2 位小数)}。<br>缺 NM tag 的记录**不再按
        NM=0 计**（那会把 ANI 顶到 100% 且无告警）：逐参考统计无 NM 记录占比，
        >10% 的参考 ANI 置 None（未测定，由分流逻辑按「未测定不拦截」处理）
        并打 WARN；≤10% 的参考缺 NM 记录直接跳过分子/分母（不计入均值）。
        """
        # 流式逐行读：capture_output 会把整个 SAM 文本读进内存，
        # 大样品（数百万条比对，数 GB 文本）直接 MemoryError 且 ANI 静默置空
        acc: dict = {}
        proc = subprocess.Popen(
            [self.samtools, 'view', '-F', '0x904', str(sorted_bam)],
            stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
            text=True, encoding='utf-8', errors='replace')
        try:
            for line in proc.stdout:
                if not line or line.startswith('@'):
                    continue
                p = line.rstrip('\n').split('\t')
                if len(p) < 11:
                    continue
                ref = p[2]
                if ref == '*':
                    continue
                cur = acc.get(ref)
                if cur is not None and cur[2] >= max_reads:
                    continue
                aln_len = _cigar_query_aln_len(p[5])
                if aln_len <= 0:
                    continue
                nm = None
                for tag in p[11:]:
                    if tag.startswith('NM:i:'):
                        try:
                            nm = int(tag[5:])
                        except ValueError:
                            nm = None
                        break
                if cur is None:
                    acc[ref] = cur = [0, 0, 0, 0]   # [分子, 分母, 记录数, 缺NM数]
                cur[2] += 1
                if nm is None:
                    cur[3] += 1
                    continue
                cur[0] += aln_len - nm
                cur[1] += aln_len
        finally:
            proc.stdout.close()
            proc.wait()
        out = {}
        no_nm_refs = []
        for ref, (num, den, n, no_nm) in acc.items():
            if n > 0 and no_nm / n > 0.1:
                out[ref] = None            # 无 NM 占比过高 → 未测定，不给假 ANI
                no_nm_refs.append(ref)
            else:
                out[ref] = round(num / den * 100.0, 2) if den > 0 else None
        if no_nm_refs:
            shown = ', '.join(no_nm_refs[:5]) + ('…' if len(no_nm_refs) > 5 else '')
            self.log(f"ANI：{len(no_nm_refs)} 条参考的无 NM tag 记录占比 >10%，"
                     f"已置为未测定（None），不再按 NM=0 记成 100%：{shown}",
                     level='warning')
        return out

    def coverage_samtools(self, sorted_bam, out_dir, name='coverage'):
        """samtools coverage -> refstats 格式"""
        out_dir = Path(out_dir)
        raw = out_dir / f'{name}.samtools_cov.tsv'
        cmd = [self.samtools, 'coverage', str(sorted_bam)]
        with open(raw, 'w', encoding='utf-8') as fl:
            subprocess.run(cmd, stdout=fl, stderr=subprocess.STDOUT, check=True, timeout=7200)

        _reset_parse_failures()
        refstats = {}
        with open(raw, 'r', encoding='utf-8') as f:
            header = f.readline().rstrip('\n').split('\t')
            ci = {h.lstrip('#'): i for i, h in enumerate(header)}
            for line in f:
                p = line.rstrip('\n').split('\t')
                if len(p) < len(header):
                    continue
                rname = p[ci['rname']]
                refstats[rname] = {
                    'Length': _to_int(p[ci['endpos']]) - _to_int(p[ci['startpos']]) + 1,
                    'Mapped_Reads': _to_int(p[ci['numreads']]),
                    'Covered_Bases': _to_int(p[ci['covbases']]),
                    'Coverage(%)': _to_float(p[ci['coverage']]),
                    'MeanDepth': _to_float(p[ci['meandepth']]),
                }
        self._warn_parse_failures('samtools coverage')
        return refstats

    def _warn_parse_failures(self, tool):
        """解析失败显式告警（静默返 0 是假阴性来源）。"""
        n = _PARSE_FAILURES['int'] + _PARSE_FAILURES['float']
        if n:
            self.log(f"{tool} 输出解析失败 {n} 处（已按 0 处理）："
                     f"覆盖度/深度可能偏低，阈值判定受影响，请检查输出文件",
                     level='warning')

    def coverage_pandepth(self, sorted_bam, out_dir, name='coverage',
                          min_mapq=None):
        """pandepth -> refstats 格式（严格照原管线，无 -q）

        min_mapq=None 时不传 -q（原管线行为，MAPQ 0 起算）。

        pandepth 的 hts-3.dll 还有传递依赖（libcrypto/libcurl 等），
        这些 DLL 在 samtools\bin 里，所以调用时必须把该目录注入子进程 PATH。
        """
        out_dir = Path(out_dir)
        prefix = out_dir / f'{name}.pandepth'
        logf = out_dir / f'{name}_pandepth.log'
        cmd = [self.pandepth, '-a', '-i', str(sorted_bam), '-o', str(prefix),
               '-t', str(self.threads)]
        if min_mapq is not None:
            cmd += ['-q', str(min_mapq)]
        env = dict(os.environ)
        st_bin = str(Path(self.samtools).resolve().parent)
        pd_bin = str(Path(self.pandepth).resolve().parent)
        env['PATH'] = os.pathsep.join(
            [pd_bin, st_bin] + [env.get('PATH', '')])
        with open(logf, 'w', encoding='utf-8') as fl:
            subprocess.run(cmd, stdout=fl, stderr=subprocess.STDOUT, check=True,
                           env=env)

        stat_gz = Path(str(prefix) + '.chr.stat.gz')
        if not stat_gz.exists():
            raise RuntimeError(f"pandepth 未产生 {stat_gz}")

        # 原管线解析：跳过 # 行，取第 5 列 Cov / 第 6 列 Dep
        _reset_parse_failures()
        refstats = {}
        with gzip.open(stat_gz, 'rt', encoding='utf-8') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                p = line.rstrip('\n').split('\t')
                if len(p) >= 6:
                    refstats[p[0].strip()] = {
                        'Length': _to_int(p[1]),
                        'Covered_Bases': _to_int(p[2]),
                        'Coverage(%)': _to_float(p[4]),
                        'MeanDepth': _to_float(p[5]),
                    }
        self._warn_parse_failures('pandepth')
        return refstats


# pandepth/samtools 输出解析失败计数（模块级）：失败静默返 0 会让覆盖度/深度
# 偏低却照常参与阈值判定（假阴性来源）。每个覆盖度统计调用前清零、
# 结束后 _warn_parse_failures 显式告警（kv_engines 第 20 项）。
_PARSE_FAILURES = {'int': 0, 'float': 0}


def _reset_parse_failures():
    _PARSE_FAILURES['int'] = 0
    _PARSE_FAILURES['float'] = 0


def _to_int(v):
    try:
        return int(float(v))
    except (TypeError, ValueError):
        _PARSE_FAILURES['int'] += 1
        return 0


def _to_float(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        _PARSE_FAILURES['float'] += 1
        return 0.0


def _cigar_query_aln_len(cigar):
    """CIGAR 串中消耗 query 的碱基数（= pysam 的 query_alignment_length）

    消耗 query 的操作：M / I / = / X（samtools 的 S/H 是软/硬剪辑，
    pysam query_alignment_length 不计；D/N 消耗 reference 而非 query）。
    末段无长度数字（如 "*"）或解析失败时返回 0，由调用方跳过。
    """
    if not cigar or cigar == '*':
        return 0
    total = 0
    num = ''
    for ch in cigar:
        if ch.isdigit():
            num += ch
            continue
        if num and ch in 'MI=X':
            total += int(num)
        num = ''
    return total


def make_engine(name, tools, threads=8, logger=None, min_mapq=10,
                coverage_tool='pandepth'):
    """
    按名字构造引擎。tools 是 ToolRegistry。

    引擎与库正交：一个鉴定库可同时带 salmon_k31/ 与 minibwa/ 两套索引，
    选哪个引擎只决定「reads 怎么定量」，不决定「用哪个库」：
      salmon  —— EM 定量（quant.sf），--writeBam 出映射位点（默认）
      minibwa —— 真比对计数，排序 BAM 交给 BamCoverageBackend
    共识段不经本工厂（kv_consensus 内部固定 minibwa 自建单参考索引）。

    coverage_tool: pandepth | samtools（覆盖度/深度统计后端，默认 pandepth，
      严格照原管线：samtools view -F 0x04 建 BAM + pandepth -a，无 -q）。
      samtools 仅供快速对拍；两条路径共用同一份 BAM。
    """
    name = (name or '').lower()
    if name not in ('salmon', 'minibwa'):
        raise ValueError(
            f"不支持的引擎: {name}（可选 salmon / minibwa；"
            "共识段内部固定 minibwa，不经本工厂）")
    tools.require(name)
    tool = (coverage_tool or 'pandepth').lower()
    if tool == 'builtin':
        # builtin（自研 SAM/CIGAR 解析器）已删除：原管线覆盖度只有 pandepth
        # 一条口径，自研实现属平行口径。历史取值（旧前端/旧脚本）不再直接
        # ValueError 崩掉，告警后回退 pandepth（缺件时再由下面逻辑降级）。
        if logger:
            logger.warning("coverage_tool=builtin 已弃用（自研解析器已删除），"
                           "自动回退 pandepth")
        tool = 'pandepth'
    if tool not in ('pandepth', 'samtools'):
        raise ValueError(f"不支持的覆盖度后端: {tool}（可选 pandepth / samtools）")
    samtools = tools.paths.get('samtools')
    pandepth = tools.paths.get('pandepth') if tool == 'pandepth' else None
    if not samtools or (tool == 'pandepth' and not pandepth):
        if name == 'minibwa' and not samtools:
            # 真比对引擎没有 samtools 就产不出 BAM，也没有定量可言，fail-fast
            raise ValueError('minibwa 定量需要 samtools（建排序 BAM 与计数）')
        # 其余缺件不阻断：salmon 退"无 BAM 估算"；minibwa 的 pandepth
        # 缺件在 BwaEngine.__init__ 里降级 samtools 后端
        if logger:
            logger.warning(
                "覆盖度后端缺件（samtools=%s pandepth=%s），定量照跑，"
                "覆盖度/深度退到估算且无比对 BAM",
                '有' if samtools else '缺', '有' if pandepth else '缺')
    if name == 'minibwa':
        return BwaEngine(tools.paths['minibwa'], threads=threads,
                         logger=logger, min_mapq=min_mapq,
                         samtools=samtools, pandepth=pandepth,
                         coverage_tool=tool)
    return SalmonEngine(tools.paths['salmon'], threads=threads, logger=logger,
                        min_mapq=min_mapq, samtools=samtools,
                        pandepth=pandepth, coverage_tool=tool)
