#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
kv_identify.py — 鉴定段
======================
reads → 引擎比对/定量 → 覆盖度统计 → 泊松打假 → 双轨过滤 → 分类表

设计来源: virome_analysis_pipeline/batch_virus_depth.py 的
process_sample / summarize_results_polars 核心逻辑
（复制后重构，原文件不动）

保留的原管线逻辑:
  - 泊松打假 λ = Reads × AvgReadLen / Length, Predicted_Support = 1-exp(-λ)
  - Poisson_Ratio = (Coverage/100) / Predicted_Support
  - 双轨过滤: A 轨全基因组泊松 + B 轨活跃转录区
  - is_segmented 判定（多段物种须全段检出）
  - ANI 分物种阈值分流（Avg_Read_ANI：只有真比对引擎 minibwa 能算，
    salmon 定量无比对位点恒为 None → 视为「未测定」，不拦截）
"""

import json
import os
import subprocess
import time
from pathlib import Path

try:
    import polars as pl
except ImportError:
    raise SystemExit("需要 polars: pip install polars")

from .kv_common import fmt_time, load_ref_lengths, load_ref_info, ref_lookup
from .kv_engines import (SRNA_MAX_READ_LEN, SRNA_SALMON_K,
                         sniff_read_len_stats, looks_like_srna)


# ── 平均读长探测（对齐原管线 get_average_read_length）────────
def average_read_length(path, default=150.0):
    """流式读取前若干条记录估平均读长，避免全文件扫描"""
    try:
        buf = []
        total, n = 0, 0
        p_str = str(path)
        is_fasta = any(p_str.lower().endswith(e) for e in
                       ('.fa', '.fasta', '.fa.gz', '.fasta.gz'))
        with open(path, 'rt', encoding='utf-8', errors='replace') \
                if not p_str.lower().endswith('.gz') \
                else __import__('gzip').open(path, 'rt', encoding='utf-8',
                                             errors='replace') as f:
            for line in f:
                if not is_fasta:
                    if line.startswith('>'):      # FASTA 出现在 fastq 位置
                        is_fasta = True
                    else:
                        if len(buf) % 4 == 1:
                            total += len(line.strip())
                            n += 1
                        buf.append(line)
                        if n >= 2000:
                            break
                if is_fasta and not line.startswith('>'):
                    total += len(line.strip())
                    n += 1
                    if n >= 2000:
                        break
        return (total / n) if n else default
    except Exception as e:  # noqa: BLE001
        # 读长探测失败会让 salmon 分支的覆盖度估算偏掉，必须留痕
        import logging
        logging.getLogger('kv_identify').warning(
            '平均读长探测失败（%s: %s），退到默认 %.0f bp',
            path, e, default)
        return default


class IdentifyStage:
    """
    鉴定段。一次构造，可对多样本依次调用 process_sample()。
    """

    def __init__(self, args, engine, tools, logger, out_dir):
        self.args = args
        self.engine = engine
        self.tools = tools
        self.logger = logger
        self.out_dir = Path(out_dir)
        # 比对/处理失败的样本必须留痕：否则失败样本在汇总表里 0 行，
        # 与「真阴性」不可区分，是最危险的静默假阴性
        self.failed_samples = []

        self.ref_lengths = load_ref_lengths(args.reference, logger)
        self.ref_info = load_ref_info(getattr(args, 'ref_info', None), logger)
        self.avg_len_cache = {}
        # 引擎版本字符串缓存（批次输入指纹用，探测一次）
        self._engine_ver = None

    # ── 输入指纹（批次断点复用的安全判据）──────────────────
    @staticmethod
    def _file_fp(path):
        """文件指纹 path+size+mtime_ns（同 preprocess.input_fingerprint 口径）。"""
        if not path:
            return None
        p = Path(os.path.abspath(str(path)))
        try:
            st = p.stat()
            return [str(p), int(st.st_size), int(st.st_mtime_ns)]
        except OSError:
            return [str(p), None, None]

    def _engine_version(self):
        """引擎版本（`<exe> --version` 探测一次；失败记 unknown，不阻断）。"""
        if self._engine_ver is None:
            ver = ''
            try:
                r = subprocess.run([self.engine.exe, '--version'],
                                   capture_output=True, text=True,
                                   encoding='utf-8', errors='replace',
                                   timeout=30)
                out = (r.stdout or '') or (r.stderr or '')
                ver = out.strip().splitlines()[0] if out.strip() else ''
            except (OSError, subprocess.SubprocessError) as e:
                ver = f'unknown({type(e).__name__})'
            self._engine_ver = (ver or 'unknown')[:120]
        return self._engine_ver

    def _batch_fingerprint(self, batch, kmer):
        """批次输入指纹：r1/r2 文件 + 参考库/注释 + 引擎名/版本 + k + 关键阈值。

        断点复用批次缓存前必须核对——此前批次文件只存行数据，换了 FASTQ
        （或换库/换引擎/改阈值）重跑仍会复用旧批次结果，旧输入的行被当成
        新输入的结果（audit 第 2 项：identify 批次缓存无输入指纹）。
        """
        a = self.args
        return {
            'engine': self.engine.name,
            'engine_version': self._engine_version(),
            'coverage_tool': getattr(self.engine, 'coverage_tool', ''),
            'k': int(kmer),
            'reference': self._file_fp(getattr(a, 'reference', None)),
            'ref_info': self._file_fp(getattr(a, 'ref_info', None)),
            'samples': [
                {'name': s.get('name'),
                 'r1': self._file_fp(s.get('r1')),
                 'r2': self._file_fp(s.get('r2'))}
                for s in batch
            ],
            'thresholds': {
                'min_cov': float(getattr(a, 'min_cov', 0.0)),
                'min_depth': float(getattr(a, 'min_depth', 0.0)),
                'min_reads': int(getattr(a, 'min_reads', 0)),
                'min_tpm': float(getattr(a, 'min_tpm', 0.0)),
                'ratio': float(getattr(a, 'ratio', 0.0)),
                'ani_thresh': float(getattr(a, 'ani_thresh', 0.0)),
            },
        }

    # ── 单样本处理 ───────────────────────────────────────
    def process_sample(self, sample):
        t0 = time.time()
        sname = sample['name']
        try:
            ar = self.engine.align(self.index_path, sample, self.out_dir)
        except Exception as e:
            self.logger.error(f"[{sname}] 比对失败: {e}")
            self.failed_samples.append(
                {'sample': sname, 'stage': 'align', 'error': str(e)[:300]})
            return []

        # PE 片段长度（审查 P1-4）：salmon NumReads 计**片段**（PE 一对=1），
        # 无位点 fallback 的覆盖度/深度 = fragments × 片段长度 / L。旧实现
        # 只用 R1 平均读长，双端等长时约低估一半（min_depth/min_cov 阈值
        # 附近的参考可能被误筛）。片段长度 = mean(R1)+mean(R2)；来源在
        # Cov_Depth_Basis 列显式标注，靠近阈值的估算值不再冒充实测。
        _ck = (sample['r1'], sample.get('r2'))
        avg_len = self.avg_len_cache.get(_ck)
        if avg_len is None:
            avg_len = average_read_length(sample['r1'])
            if sample.get('r2'):
                avg_len += average_read_length(sample['r2'])
            self.avg_len_cache[_ck] = avg_len

        # refstats 来自映射 BAM + pandepth（salmon --writeBam 出位点后与
        # 真比对同源）；后端缺件时为空，覆盖度退到「reads×读长」估算
        refstats = ar.refstats or {}
        rows = []
        total_mapped = int(ar.mapped)

        for acc, qreads in ar.reads.items():
            st = refstats.get(acc, {})
            L = self.ref_lengths.get(acc) or st.get('Length') or 0
            cov = st.get('Coverage(%)', 0.0)
            depth = st.get('MeanDepth', 0.0)
            if not refstats:
                # salmon：无位点，覆盖度用「reads×读长/长度」估算，并标记来源
                cov = min(100.0, (qreads * avg_len / L * 100.0) if L else 0.0)
                depth = (qreads * avg_len / L) if L else 0.0

            info = ref_lookup(self.ref_info, acc)
            rows.append({
                'Sample': sname,
                'Accession': acc,
                'Length': int(L),
                'taxid': info.get('taxid', 'Unannotated'),
                'Species': info.get('species', acc),
                'Segment': info.get('segment', ''),
                'Molecule_Type2': info.get('molecule', '') or info.get('molecule_type', ''),
                # 宿主/地理/GenBank 标题：来自参考注释表（提取信息不再只有 4 列）
                'Host': info.get('host', ''),
                'Geo_Location': info.get('geo', ''),
                'GenBank_Title': info.get('genbank_title', ''),
                'Coverage(%)': float(cov),
                'MeanDepth': float(depth),
                # 覆盖/深度来源：真比对统计 or reads×读长估算（P1-4）
                'Cov_Depth_Basis': ('bam_refstats' if refstats
                                    else 'est_reads_x_len'),
                'EM_Reads': float(qreads),
                'Uniq_Reads': int(qreads),
                # ANI（只有真比对引擎能算；salmon 无比对位点，置 None）
                'Avg_Read_ANI': st.get('Avg_Read_ANI'),
                'Sample_Total_Mapped': int(total_mapped),
                'Avg_Read_Len': float(avg_len),
            })

        # 写每样本深度表
        stats_dir = self.out_dir / 'stats'
        stats_dir.mkdir(parents=True, exist_ok=True)
        if refstats:
            with open(stats_dir / f'{sname}.refstats.tsv', 'w', encoding='utf-8') as f:
                f.write("Accession\tLength\tMapped_Reads\tCovered_Bases\tCoverage(%)\tMeanDepth\n")
                for acc, st in refstats.items():
                    if st['Mapped_Reads'] > 0 or st['Covered_Bases'] > 0:
                        f.write(f"{acc}\t{st['Length']}\t{st['Mapped_Reads']}\t"
                                f"{st['Covered_Bases']}\t{st['Coverage(%)']:.4f}\t{st['MeanDepth']:.4f}\n")

        el = time.time() - t0
        self.logger.info(f"  [{sname}] {fmt_time(el)} | hits={len(rows)} | mapped={total_mapped}/{ar.total}")
        return rows

    # ── 汇总 ─────────────────────────────────────────────
    # 全量表列契约（零命中时也按此写**带表头的空表**——见下方空结果分支；
    # 列序与 process_sample 行键 + summarize 派生列一致，下游 kv_filter /
    # kv_stage._summarize / GUI / benchmark 解析器都按这套列名读）。
    _SUMMARY_COLS = (
        'Sample', 'Accession', 'Length', 'taxid', 'Species', 'Segment',
        'Molecule_Type2', 'Host', 'Geo_Location', 'GenBank_Title',
        'Coverage(%)', 'MeanDepth', 'Cov_Depth_Basis', 'EM_Reads',
        'Uniq_Reads', 'Avg_Read_ANI', 'Sample_Total_Mapped', 'Avg_Read_Len',
        'Predicted_Support', 'Poisson_Ratio', 'Asm_TPM', 'ani_screened',
        'ani_status',
    )

    def _write_empty_tables(self):
        """零命中也要产出**本轮的空结果**并清除上一轮阳性残留。

        审查 P0-1（2026-10-10 复核成立）：旧实现零命中时直接返回、不触碰
        磁盘，上一轮的 all_viruses.best.summary.tsv / unclassified.tsv 原样
        留在 identify/ 里；下游 filter 从磁盘读表、外层 _summarize 再把旧
        物种清单配上**本轮输入指纹**写进 summary.json——上次阳性被当成
        本次检出。铁律：每个阶段的磁盘产物必须反映"本轮"，空结果也是
        结果，一样要落盘覆盖。
        """
        out = self.out_dir / 'identify'
        out.mkdir(parents=True, exist_ok=True)
        empty = pl.DataFrame(schema={c: pl.Utf8 for c in self._SUMMARY_COLS})
        empty.write_csv(out / 'all_viruses.summary.tsv', separator='\t')
        for stale in ('all_viruses.best.summary.tsv',
                      'all_viruses.unclassified.tsv'):
            p = out / stale
            if p.exists():
                p.unlink()
                self.logger.warning(
                    f'本轮零命中：已清除上一轮残留产物 {stale}（防止旧阳性'
                    '被当作本次结果）')

    def summarize(self, all_rows):
        """
        全量结果 → 泊松打假 → 双轨过滤 → 分类输出。
        返回 (summary_df, confirmed_df, novel_df)
        """
        out = self.out_dir / 'identify'
        out.mkdir(parents=True, exist_ok=True)

        if not all_rows:
            self.logger.warning(
                '无任何比对命中：写入本轮空结果表并清除旧产物（P0-1 修复）')
            self._write_empty_tables()
            empty = pl.DataFrame()
            return empty, empty, empty

        # infer_schema_length=None 全量扫描定 schema：minibwa 的 Avg_Read_ANI
        # 是浮点且大量行仅前段为 None——默认只看前 100 行会推断成 Null 列，
        # 追加后面的浮点值直接 ComputeError（实测 D11 minibwa 全军覆没）
        df = pl.DataFrame(all_rows, infer_schema_length=None)

        # ── 泊松打假（对齐原管线公式）──
        df = df.with_columns([
            (pl.col('EM_Reads') * pl.col('Avg_Read_Len') / pl.col('Length')).alias('_lambda')
        ]).with_columns([
            (1.0 - (-pl.col('_lambda')).exp()).alias('Predicted_Support')
        ]).with_columns([
            pl.when(pl.col('Predicted_Support') > 0)
              .then((pl.col('Coverage(%)') / 100.0) / pl.col('Predicted_Support'))
              .otherwise(0.0).alias('Poisson_Ratio')
        ]).drop('_lambda')

        # ── TPM（识别表内产出，供 filter 的 --min-tpm 门槛消费）──
        # 口径与 kv_fig._add_derived_metrics 完全一致：EM_Reads 定量，
        # TPM = (EM_Reads / Length_kb) / Σ_同样本(EM_Reads / Length_kb) × 1e6，
        # 同一样本内 sum(TPM)=1e6。此前识别表没有任何 TPM 列，filter 的
        # tpm 列查找恒为 None，文档承诺的 --min-tpm 门槛从未生效（静默跳过）。
        df = df.with_columns([
            pl.when(pl.col('Length') > 0)
              .then(pl.col('EM_Reads') / (pl.col('Length') / 1000.0))
              .otherwise(0.0).alias('_tpk')
        ]).with_columns([
            pl.when(pl.col('_tpk').sum().over('Sample') > 0)
              .then(pl.col('_tpk') / pl.col('_tpk').sum().over('Sample') * 1e6)
              .otherwise(0.0).alias('Asm_TPM')
        ]).drop('_tpk')
        self.logger.info(
            "识别表已产出 TPM 列（Asm_TPM：EM_Reads/Length 归一化，"
            "同一样本内 sum(TPM)=1e6，口径同 kv_fig）——filter 的 --min-tpm "
            "门槛自此真实生效")

        # ── ANI 分流生效性标注 ──
        # salmon 无比对位点（NM/CIGAR 是伪的），Avg_Read_ANI 恒空 → 全部按
        # 已知种处理、unclassified 永不产出。给结果表加显式字段并在日志
        # 明说，避免把「没做 ANI 分流」误读成「ANI 全部合格」。
        ani_screened = bool(
            len(df) and df['Avg_Read_ANI'].is_not_null().any())
        df = df.with_columns(pl.lit(ani_screened).alias('ani_screened'))
        if not ani_screened:
            self.logger.warning(
                "ANI 分流未生效（Avg_Read_ANI 全空：默认 salmon 引擎无比对"
                "位点）：全部按已知种处理，疑似新种/unclassified 不会产出")

        df.write_csv(out / 'all_viruses.summary.tsv', separator='\t')
        self.logger.info(f"全量鉴定表: {len(df)} 行 -> {out / 'all_viruses.summary.tsv'}")

        # ── 基础门槛 ──
        base = (
            (pl.col('Sample_Total_Mapped') > 0) &
            (pl.col('Uniq_Reads') >= self.args.min_reads) &
            (pl.col('MeanDepth') >= self.args.min_depth)
        )

        # ── 双轨：A 轨全基因组 ──
        # is_rna 判定（RNA 病毒宽松/ DNA 病毒更严）：本库 Molecule_Type2 词表
        # 是 'RNA'/'DNA'/''ss-RNA'/'ds-RNA'（无 'ssRNA' 写法），旧正则
        # ssRNA|dsRNA|mRNA|cRNA 只命中 360/12159；改为能命中
        # RNA/ss-RNA/ds-RNA/cRNA/mRNA 的写法（实测命中 7400 行）。
        is_rna = pl.col('Molecule_Type2').str.contains(
            r'(?i)\b(ss|ds)?-?rna\b|mrna|crna').fill_null(False)
        n_rna = int(df.select(is_rna.sum()).item() or 0)
        self.logger.info(
            f"is_rna 判定（Molecule_Type2 词表 RNA/ss-RNA/ds-RNA/cRNA/mRNA）："
            f"命中 {n_rna}/{len(df)} 行")
        # DNA 更严分支：恢复原管线 batch_virus_depth.py:622-625 口径
        # ratio if is_rna else max(ratio, 0.5)（DNA 病毒泊松比门槛不低于 0.5）。
        ratio_eff = pl.when(is_rna).then(pl.lit(float(self.args.ratio))) \
            .otherwise(pl.lit(max(float(self.args.ratio), 0.5)))
        track_a = (
            (pl.col('Coverage(%)') >= self.args.min_cov) &
            (pl.col('Poisson_Ratio') >= ratio_eff)
        )

        # ── 双轨：B 轨活跃转录区（可选）──
        passed = df.filter(base & track_a)
        if not getattr(self.args, 'genes_cov', None):
            # 基因覆盖表在平台内没有生产者（只有本段消费），B 轨实际不可达。
            # 显式告警，避免「没有 B 轨结论」被当成「B 轨判定为阴性」。
            self.logger.warning(
                "B 轨跳过：无基因覆盖数据（--genes-cov 未提供；"
                "平台当前无基因覆盖表生产者，B 轨不可达）")
        else:
            genes = Path(self.args.genes_cov)
            if genes.exists():
                gdf = pl.read_csv(str(genes), separator='\t', ignore_errors=True)
                if {'seqid', 'gene_total_cov', 'gene_avr_cov'}.issubset(set(gdf.columns)):
                    gsub = gdf.select(['seqid', 'gene_total_cov', 'gene_avr_cov'])
                    merged = df.join(gsub, left_on='Accession', right_on='seqid', how='left')
                    n_rna_b = int(merged.select(is_rna.sum()).item() or 0)
                    self.logger.info(
                        f"B 轨 is_rna 判定命中 {n_rna_b}/{len(merged)} 行")
                    gene_pass = (
                        (pl.col('gene_total_cov').fill_null(0) >= self.args.min_gene_total_cov) &
                        (pl.col('gene_avr_cov').fill_null(0) >= self.args.min_gene_avr_cov)
                    )
                    rescue = merged.filter(base & is_rna & gene_pass).drop(
                        ['seqid', 'gene_total_cov', 'gene_avr_cov'])
                    passed = pl.concat([passed, rescue.select(passed.columns)], how='diagonal_relaxed')
                    self.logger.info(f"B 轨（转录区）额外放行: {len(rescue)} 行")
            else:
                self.logger.warning(f"genes_cov 文件不存在: {genes}")

        # 去重
        if len(passed):
            passed = passed.unique(subset=['Sample', 'Accession'], keep='first')

        # ── ANI 分流（minibwa 可算，salmon 置空＝未测定不拦截）──
        if 'Avg_Read_ANI' not in passed.columns:
            passed = passed.with_columns(
                pl.lit(None).cast(pl.Float64).alias('Avg_Read_ANI'))

        # ani_status（审查 P1-6）：把"筛选通过但 ANI 未核验"与"ANI 已核验
        # 达标/低于阈值"在表内显式分开——默认 salmon 引擎无比对位点，此前
        # 两类行混在同一 best.summary.tsv 里，"confirmed"名不副实（论文
        # 口径的 undetermined 无从体现）。文件布局保持不变（best/
        # unclassified 仍按惯例产出，GUI/解析器零改动），但每行自带核验
        # 状态：ani_confirmed / ani_below_threshold / ani_unavailable。
        sp_thresh = self.args.ani_thresh
        passed = passed.with_columns(
            pl.when(pl.col('Avg_Read_ANI').is_null())
              .then(pl.lit('ani_unavailable'))
              .when(pl.col('Avg_Read_ANI') >= sp_thresh)
              .then(pl.lit('ani_confirmed'))
              .otherwise(pl.lit('ani_below_threshold'))
              .alias('ani_status'))
        confirmed = passed.filter(
            pl.col('Avg_Read_ANI').is_null() | (pl.col('Avg_Read_ANI') >= sp_thresh)
        )
        novel = passed.filter(
            pl.col('Avg_Read_ANI').is_not_null() & (pl.col('Avg_Read_ANI') < sp_thresh)
        )

        # is_segmented 判定
        seg_counts = {}
        for acc, rec in self.ref_info.items():
            sp = rec.get('species')
            seg = (rec.get('segment') or '').strip()
            if sp and seg:
                seg_counts.setdefault(sp, set()).add(seg.upper())
        multi_seg = {sp for sp, s in seg_counts.items() if len(s) > 1}

        def _add_seg_flag(d):
            if not len(d):
                return d
            return d.with_columns(
                pl.col('Species').is_in(list(multi_seg)).alias('is_segmented'))

        confirmed = _add_seg_flag(confirmed)
        novel = _add_seg_flag(novel)

        if len(novel):
            novel.write_csv(out / 'all_viruses.unclassified.tsv', separator='\t')
        if len(confirmed):
            confirmed.write_csv(out / 'all_viruses.best.summary.tsv', separator='\t')

        self.logger.info(f"基础+双轨放行: {len(passed)} 行 -> 确诊 {len(confirmed)} / 疑似新种 {len(novel)}")
        return df, confirmed, novel

    # ── 主流程 ───────────────────────────────────────────
    def run(self, samples, resume=True):
        self.logger.info("=" * 60)
        self.logger.info(f"【鉴定段】引擎={self.engine.name} 样本数={len(samples)}")
        self.logger.info(f"阈值: cov>={self.args.min_cov}% depth>={self.args.min_depth} "
                         f"reads>={self.args.min_reads} poisson>={self.args.ratio}")
        self.logger.info("=" * 60)

        idx_dir = getattr(self.args, 'index_dir', None) or (self.out_dir / 'index')
        # 小RNA 自动适配：reads 超短（<32bp）时 salmon 索引改用小 k——
        # 实测 k=31 索引对 21-24nt reads 映射率为 0（reads < k 无种子可用），
        # k=15 对同一批数据映射率 ~11%，定量与 BWA 真值差 3-6%。
        # 索引按 k 分目录并存（salmon_k31 仍是库探针与常规默认），
        # 建一次后跨 run 复用。
        kmer = 31
        try:
            _stats = None
            for _s in samples[:4]:
                for _p in (_s.get('r1'), _s.get('r2')):
                    if _p:
                        _v = sniff_read_len_stats(_p)
                        # 任一样本像小RNA 即按小RNA 处理（p90 判定，
                        # 避免一条长 read 让整批走 k=31 → 映射率为 0）
                        if _v and looks_like_srna(_v):
                            _stats = _v
                            break
                if _stats:
                    break
            if _stats:
                kmer = SRNA_SALMON_K
                self.logger.warning(
                    f"检测到超短 reads（头部采样中位 {_stats['median']}bp / "
                    f"p90 {_stats['p90']}bp，疑似小RNA/siRNA 数据）："
                    f"salmon 索引改用 k={SRNA_SALMON_K}"
                    f"（默认 k=31 对 <32bp reads 映射率为 0）")
        except OSError:
            pass
        self.index_path = self.engine.build_index(self.args.reference, Path(idx_dir), self.args.align_threads, kmer=kmer)

        batches_dir = self.out_dir / 'batches'
        batches_dir.mkdir(parents=True, exist_ok=True)

        all_rows = []
        t0 = time.time()
        for i in range(0, len(samples), self.args.batch_size):
            batch = samples[i:i + self.args.batch_size]
            bi = i // self.args.batch_size + 1
            bfile = batches_dir / f'batch_{bi:04d}.json'
            # 批次缓存必须带输入指纹（旧 FASTQ/旧库/改阈值都不得复用）
            batch_fp = self._batch_fingerprint(batch, kmer)
            if resume and bfile.exists():
                payload = None
                try:
                    payload = json.loads(bfile.read_text(encoding='utf-8'))
                except (OSError, ValueError) as e:
                    self.logger.warning(f"批次 {bi} 缓存读取失败（{e}），重算该批次")
                if isinstance(payload, dict) and \
                        payload.get('input_fingerprint') == batch_fp:
                    self.logger.info(f"批次 {bi} 已存在且输入指纹一致，跳过")
                    all_rows.extend(payload.get('rows') or [])
                    continue
                self.logger.warning(
                    f"批次 {bi} 缓存与当前输入不匹配（旧 FASTQ/参数/库，"
                    f"或旧格式无指纹）：重算该批次，不复用旧结果")

            self.logger.info(f"批次 {bi}: {len(batch)} 个样本")
            rows = []
            for s in batch:
                try:
                    rows.extend(self.process_sample(s) or [])
                except Exception as e:
                    self.logger.error(f"  [{s['name']}] 异常: {e}")
                    self.failed_samples.append(
                        {'sample': s['name'], 'stage': 'process',
                         'error': str(e)[:300]})
            bfile.write_text(
                json.dumps({'input_fingerprint': batch_fp, 'rows': rows},
                           ensure_ascii=False), encoding='utf-8')
            all_rows.extend(rows)

        if self.failed_samples:
            try:
                fpath = self.out_dir / 'failed_samples.json'
                fpath.write_text(
                    json.dumps(self.failed_samples, ensure_ascii=False, indent=1),
                    encoding='utf-8')
                self.logger.warning(
                    f"⚠ {len(self.failed_samples)} 个样本处理失败（见 {fpath}）："
                    f"这些样本没有进入结果表，不能当作阴性解读！")
            except OSError as e:
                self.logger.warning(f"失败样本清单写盘失败: {e}")

        self.logger.info(f"鉴定段总耗时 {fmt_time(time.time()-t0)}，命中记录 {len(all_rows)}")
        return self.summarize(all_rows)
