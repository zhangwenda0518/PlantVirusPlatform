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

# 2026-10-20：由 polars 改为 pandas。原因：polars 的二进制运行时
# （polars-runtime-32，解压 167.7MB）是 site-packages 里最大的单项，
# 且无可替代的精简变体（compat/64/32 三个变体体积相同）。pandas 已是
# 平台硬依赖（图/表/统计处处在用），复用它能省下这 167.7MB。
# 语义对照见下方各处的注释；数值口径与 polars 版逐字对齐。
import numpy as np
import pandas as pd

from .kv_common import fmt_time, load_ref_lengths, load_ref_info, ref_lookup
from .kv_engines import (SRNA_MAX_READ_LEN, SRNA_SALMON_K,
                         sniff_read_len_stats, looks_like_srna)


def _write_tsv(df, path):
    """落盘 TSV，布尔列写成 true/false（与 polars 输出逐字节一致）。

    为什么较真：这些表被下游（kv_stage / kv_filter / kv_fig / GUI 解析器）
    按字符串读，历史产物里 ani_screened/is_segmented 一直是小写 true/false。
    pandas 默认写 True/False，虽然多数解析器都能容忍，但会让同一批数据在
    「改前 / 改后」两次运行之间出现 diff 噪声，也破坏与旧结果的逐字节可比性。
    """
    out = df
    bool_cols = [c for c in df.columns
                 if df[c].dtype == bool
                 or (df[c].dtype == object and len(df[c].dropna()) > 0
                     and df[c].dropna().map(
                         lambda v: isinstance(v, (bool, np.bool_))).all())]
    if bool_cols:
        out = df.copy()
        for c in bool_cols:
            out[c] = out[c].map(
                lambda v: '' if v is None or (isinstance(v, float) and np.isnan(v))
                else ('true' if v is True or v == True else
                      'false' if v is False or v == False else v))  # noqa: E712
    out.to_csv(path, sep='\t', index=False)


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
        # 记录配对性（定量双单位列用：EM_Reads=fragments，PE 的 Est_Reads=×2）
        if not hasattr(self, '_paired_samples'):
            self._paired_samples = {}
        self._paired_samples[sname] = bool(sample.get('r2'))
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
                # ICTV 种名：同种异名条目的聚合键（如 CGB 条目=HSVd 同种）
                'Ictv_Species': info.get('ictv', '') or info.get('species', ''),
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
        # 空表带 schema 落盘：pandas 的 to_csv 对空 DataFrame 只写表头，
        # 与 polars 行为一致（下游按列名读取，0 行即本轮空结果）。
        empty = pd.DataFrame({c: pd.Series(dtype='object')
                              for c in self._SUMMARY_COLS})
        _write_tsv(empty, out / 'all_viruses.summary.tsv')
        for stale in ('all_viruses.best.summary.tsv',
                      'all_viruses.unclassified.tsv',
                      'all_viruses.review.tsv'):
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
            # 零命中同样触发发现轨建议（HSVd 案例：比对轨不可见 ≠ 无病毒）
            (out / 'suggest_discovery_assembly.txt').write_text(
                "本样本零比对命中：建议运行组装发现轨（assemble_and_classify，"
                "srna/rnaviral 模式 + kunpeng contig 分类 + BLAST 定种）。\n"
                "实证：VIROMOCK D3 的 HSVd 零 reads 映射，组装发现轨以 "
                "NODE_24234（371bp，59.6×，94.5% to NC_001351.1）回收。\n",
                encoding='utf-8')
            self.logger.warning(
                "零命中 -> 已写 suggest_discovery_assembly.txt："
                "建议运行组装发现轨")
            empty = pd.DataFrame()
            return empty, empty, empty

        # polars 的 infer_schema_length=None 语义是「全量扫描定 schema」。
        # pandas 的 read 侧等价物不存在于 DataFrame(records) 构造里，但
        # pd.DataFrame(list[dict]) 会自动把混合 None/float 的列推成 float64，
        # 不会出现 polars 那个「前 100 行全 None → Null 列 → 后续追加浮点
        # 直接 ComputeError」的坑（实测 D11 minibwa 全军覆没的那个 bug）。
        # 数值列显式 to_numeric 兜底，字符串列保持 object。
        df = pd.DataFrame(all_rows)

        # ── 泊松打假（对齐原管线公式）──
        # λ = EM_Reads × Avg_Read_Len / Length
        # Predicted_Support = 1 - exp(-λ)
        # Poisson_Ratio = (Coverage/100) / Predicted_Support，分母 0 → 0
        # 与 polars 版逐字等价：polars 的 when(...>0).then(...).otherwise(0)
        # 对应 numpy 的 np.where；-(-λ).exp() 就是 1-exp(-λ) 的写法。
        _em = pd.to_numeric(df['EM_Reads'], errors='coerce')
        _arl = pd.to_numeric(df['Avg_Read_Len'], errors='coerce')
        _len = pd.to_numeric(df['Length'], errors='coerce')
        _cov = pd.to_numeric(df['Coverage(%)'], errors='coerce')
        _lambda = _em * _arl / _len
        _support = 1.0 - np.exp(-_lambda)
        df['Poisson_Ratio'] = np.where(
            _support > 0, (_cov / 100.0) / _support, 0.0)
        # 中间列不落表（polars 版 drop('_lambda')；_support 同理）

        # ── TPM（识别表内产出，供 filter 的 --min-tpm 门槛消费）──
        # 口径与 kv_fig._add_derived_metrics 完全一致：EM_Reads 定量，
        # TPM = (EM_Reads / Length_kb) / Σ_同样本(EM_Reads / Length_kb) × 1e6，
        # 同一样本内 sum(TPM)=1e6。
        # polars: col('_tpk').sum().over('Sample') → pandas: groupby().transform('sum')
        _tpk = np.where(_len > 0, _em / (_len / 1000.0), 0.0)
        _tpk_s = pd.Series(_tpk, index=df.index)
        _sum_by_sample = _tpk_s.groupby(df['Sample']).transform('sum')
        df['Asm_TPM'] = np.where(_sum_by_sample > 0,
                                 _tpk_s / _sum_by_sample * 1e6, 0.0)
        self.logger.info(
            "识别表已产出 TPM 列（Asm_TPM：EM_Reads/Length 归一化，"
            "同一样本内 sum(TPM)=1e6，口径同 kv_fig）——filter 的 --min-tpm "
            "门槛自此真实生效")

        # ── ANI 分流生效性标注 ──
        # salmon 无比对位点（NM/CIGAR 是伪的），Avg_Read_ANI 恒空 → 全部按
        # 已知种处理、unclassified 永不产出。给结果表加显式字段并在日志
        # 明说，避免把「没做 ANI 分流」误读成「ANI 全部合格」。
        ani_screened = bool(
            len(df) and pd.to_numeric(df['Avg_Read_ANI'],
                                      errors='coerce').notna().any())
        df['ani_screened'] = ani_screened
        if not ani_screened:
            self.logger.warning(
                "ANI 分流未生效（Avg_Read_ANI 全空：默认 salmon 引擎无比对"
                "位点）：全部按已知种处理，疑似新种/unclassified 不会产出")

        _write_tsv(df, out / 'all_viruses.summary.tsv')
        self.logger.info(f"全量鉴定表: {len(df)} 行 -> {out / 'all_viruses.summary.tsv'}")

        # ── 基础门槛 ──
        base = (
            (pd.to_numeric(df['Sample_Total_Mapped'], errors='coerce') > 0) &
            (pd.to_numeric(df['Uniq_Reads'], errors='coerce') >= self.args.min_reads) &
            (pd.to_numeric(df['MeanDepth'], errors='coerce') >= self.args.min_depth)
        )

        # ── 双轨：A 轨全基因组 ──
        # is_rna 判定（RNA 病毒宽松/ DNA 病毒更严）：本库 Molecule_Type2 词表
        # 是 'RNA'/'DNA'/''ss-RNA'/'ds-RNA'（无 'ssRNA' 写法），旧正则
        # ssRNA|dsRNA|mRNA|cRNA 只命中 360/12159；改为能命中
        # RNA/ss-RNA/ds-RNA/cRNA/mRNA 的写法（实测命中 7400 行）。
        # polars 的 str.contains(..., fill_null(False)) → pandas 的
        # str.contains(na=False)；regex=True 是默认，显式写出以免版本歧义。
        # 组用非捕获 (?:ss|ds)：pandas 对含捕获组的正则会发 UserWarning
        # （"has match groups"），匹配语义完全相同。
        is_rna = df['Molecule_Type2'].astype('object').str.contains(
            r'(?i)\b(?:ss|ds)?-?rna\b|mrna|crna', regex=True, na=False)
        n_rna = int(is_rna.sum())
        self.logger.info(
            f"is_rna 判定（Molecule_Type2 词表 RNA/ss-RNA/ds-RNA/cRNA/mRNA）："
            f"命中 {n_rna}/{len(df)} 行")
        # DNA 更严分支：恢复原管线 batch_virus_depth.py:622-625 口径
        # ratio if is_rna else max(ratio, 0.5)（DNA 病毒泊松比门槛不低于 0.5）。
        # polars 的 when/then/otherwise 逐行选择 → numpy 的 np.where。
        ratio_eff = np.where(is_rna, float(self.args.ratio),
                             max(float(self.args.ratio), 0.5))
        track_a = (
            (pd.to_numeric(df['Coverage(%)'], errors='coerce') >= self.args.min_cov) &
            (pd.to_numeric(df['Poisson_Ratio'], errors='coerce') >= ratio_eff)
        )

        # ── 双轨：B 轨活跃转录区（可选）──
        passed = df[base & track_a].copy()
        if not getattr(self.args, 'genes_cov', None):
            # 基因覆盖表在平台内没有生产者（只有本段消费），B 轨实际不可达。
            # 显式告警，避免「没有 B 轨结论」被当成「B 轨判定为阴性」。
            self.logger.warning(
                "B 轨跳过：无基因覆盖数据（--genes-cov 未提供；"
                "平台当前无基因覆盖表生产者，B 轨不可达）")
        else:
            genes = Path(self.args.genes_cov)
            if genes.exists():
                # ignore_errors=True 的 polars 语义 → pandas 的
                # engine='python' + on_bad_lines='skip'（容错跳坏行）。
                gdf = pd.read_csv(str(genes), sep='\t',
                                  on_bad_lines='skip', engine='python')
                if {'seqid', 'gene_total_cov', 'gene_avr_cov'}.issubset(set(gdf.columns)):
                    gsub = gdf[['seqid', 'gene_total_cov', 'gene_avr_cov']]
                    # polars join(left_on='Accession', right_on='seqid', how='left')
                    # → pandas merge(left_on='Accession', right_on='seqid', how='left')
                    merged = df.merge(gsub, left_on='Accession',
                                      right_on='seqid', how='left')
                    _is_rna_b = merged['Molecule_Type2'].astype('object').str.contains(
                        r'(?i)\b(?:ss|ds)?-?rna\b|mrna|crna', regex=True, na=False)
                    n_rna_b = int(_is_rna_b.sum())
                    self.logger.info(
                        f"B 轨 is_rna 判定命中 {n_rna_b}/{len(merged)} 行")
                    _base_b = (
                        (pd.to_numeric(merged['Sample_Total_Mapped'], errors='coerce') > 0) &
                        (pd.to_numeric(merged['Uniq_Reads'], errors='coerce') >= self.args.min_reads) &
                        (pd.to_numeric(merged['MeanDepth'], errors='coerce') >= self.args.min_depth)
                    )
                    gene_pass = (
                        (pd.to_numeric(merged['gene_total_cov'], errors='coerce').fillna(0)
                         >= self.args.min_gene_total_cov) &
                        (pd.to_numeric(merged['gene_avr_cov'], errors='coerce').fillna(0)
                         >= self.args.min_gene_avr_cov)
                    )
                    rescue = merged[_base_b & _is_rna_b & gene_pass].drop(
                        columns=['seqid', 'gene_total_cov', 'gene_avr_cov'])
                    # polars 的 concat(how='diagonal_relaxed') 允许列集不同；
                    # pandas 的 concat 默认同样按列名对齐、缺失填 NaN，
                    # 只需 reindex 到 passed 的列序保证写出稳定。
                    rescue = rescue.reindex(columns=passed.columns)
                    passed = pd.concat([passed, rescue], ignore_index=True)
                    self.logger.info(f"B 轨（转录区）额外放行: {len(rescue)} 行")
            else:
                self.logger.warning(f"genes_cov 文件不存在: {genes}")

        # 去重
        if len(passed):
            passed = passed.drop_duplicates(
                subset=['Sample', 'Accession'], keep='first').reset_index(drop=True)

        # ── ANI 分流（minibwa 可算，salmon 置空＝未测定不拦截）──
        if 'Avg_Read_ANI' not in passed.columns:
            passed['Avg_Read_ANI'] = np.nan

        # ani_status（审查 P1-6）：把"筛选通过但 ANI 未核验"与"ANI 已核验
        # 达标/低于阈值"在表内显式分开——默认 salmon 引擎无比对位点，此前
        # 两类行混在同一 best.summary.tsv 里，"confirmed"名不副实（论文
        # 口径的 undetermined 无从体现）。文件布局保持不变（best/
        # unclassified 仍按惯例产出，GUI/解析器零改动），但每行自带核验
        # 状态：ani_confirmed / ani_below_threshold / ani_unavailable。
        sp_thresh = self.args.ani_thresh
        _ani = pd.to_numeric(passed['Avg_Read_ANI'], errors='coerce')
        passed['ani_status'] = np.select(
            [_ani.isna(), _ani >= sp_thresh],
            ['ani_unavailable', 'ani_confirmed'],
            default='ani_below_threshold')
        confirmed = passed[_ani.isna() | (_ani >= sp_thresh)].copy()
        novel = passed[_ani.notna() & (_ani < sp_thresh)].copy()

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
            # polars is_in → pandas isin（返回布尔 Series，直接赋列）
            d = d.copy()
            d['is_segmented'] = d['Species'].isin(list(multi_seg))
            return d

        confirmed = _add_seg_flag(confirmed)
        novel = _add_seg_flag(novel)

        # ── 近缘溢出判别列（2026-10-18 基准优化项4）──
        # 同属近缘参考的 reads 分摊/溢出是该平台 FP 的最大来源（基准 D1–D10
        # 精确率 32.8% 的主因）。对每行确诊/疑似新种给出：同属其他种的最高
        # EM、其与本行 EM 的比值。Same_Genus_EM_Ratio > 0.1 即属内歧义信号
        #（近缘溢出或 EM 分摊嫌疑），供人工复核与基准 FP 归因；不改变判定
        # 本身（是否拦截留作后续可配开关）。属名取物种名首词（ICTV 双名法
        # 拉丁名首词即属名；俗名旧记录该列退化为跨"首词组"比较，仅作提示）。
        def _genus_of(sp):
            return (sp or '').strip().split(' ')[0] if (sp or '').strip() else ''

        _em_by = {}  # (Sample, 属) -> {物种: EM 合计}
        # polars iter_rows(named=True) → pandas to_dict('records')；
        # 用 itertuples 更快，但这里要保持 None/NaN 语义一致，走 records。
        for _r in df.to_dict('records'):
            _g = _genus_of(_r.get('Species'))
            if not _g:
                continue
            _key = (_r['Sample'], _g)
            _em_by.setdefault(_key, {})
            _em_by[_key][_r['Species']] = (
                _em_by[_key].get(_r['Species'], 0.0)
                + float(_r.get('EM_Reads') or 0.0))

        def _add_near_neighbor(d):
            if not len(d):
                return d
            nxt_sp, nxt_em, nxt_ratio = [], [], []
            for r in d.to_dict('records'):
                g = _genus_of(r.get('Species'))
                own = float(r.get('EM_Reads') or 0.0)
                others = {sp: em for sp, em
                          in _em_by.get((r['Sample'], g), {}).items()
                          if sp != r.get('Species')}
                if others:
                    sp2 = max(others, key=lambda k: others[k])
                    em2 = others[sp2]
                else:
                    sp2, em2 = '', 0.0
                nxt_sp.append(sp2)
                nxt_em.append(round(em2, 1))
                nxt_ratio.append(round(em2 / own, 3) if own > 0 else None)
            d = d.copy()
            d['Same_Genus_Next_Species'] = nxt_sp
            # None → NaN 才能落成 float 列（pandas 不会像 polars 那样
            # 保留 None 的 object 语义）
            d['Same_Genus_Next_EM'] = nxt_em
            d['Same_Genus_EM_Ratio'] = [
                np.nan if v is None else v for v in nxt_ratio]
            return d

        confirmed = _add_near_neighbor(confirmed)
        novel = _add_near_neighbor(novel)

        # ── 信号分级 + 疑似档（review tier，2026-10-18 基准优化）──
        # 依据（D:/vp_retest/sweep 16 组实测）：泊松比放宽是灵敏度主闸门
        # （0.3→0.1 在 S9 放行 GVA/GRVFV），但裸放宽会把噪声（EM≤76）与
        # 小基因组真信号（GSyV-1 EM=29）混在一起；EM 按参考长度归一
        # （EM/kb）后两者可分：真信号 ≥20/kb，噪声 <20/kb。故：
        #   confirmed  = 严格档（阈值不变，加 Signal_Grade 列）
        #   review     = 严格档未过、放宽档（cov2/dep0.1/reads1/泊松0.1）
        #                通过的行 -> all_viruses.review.tsv（疑似档独立成表，
        #                不污染确诊口径；含 Signal_Grade 供人工排序）
        def _add_signal_grade(d):
            if not len(d):
                return d
            d = d.copy()
            _l = pd.to_numeric(d['Length'], errors='coerce')
            _e = pd.to_numeric(d['EM_Reads'], errors='coerce')
            _emkb = np.where(_l > 0, _e / (_l / 1000.0), 0.0)
            d['Signal_Grade'] = np.select(
                [_emkb >= 100, _emkb >= 20], ['high', 'medium'], default='low')
            return d

        confirmed = _add_signal_grade(confirmed)
        novel = _add_signal_grade(novel)

        # ── 定量双单位（2026-10-18 基准优化）：EM_Reads 是 fragments
        #（PE 一对=1 片段），与按 reads 计的真值天生差 2×（第一轮基准
        # B3/B4 图口径假象的根因）。表内同时给出 Est_Reads（PE ×2 估算）
        # 与 Quant_Units 标注，报告侧零歧义。
        paired_map = getattr(self, '_paired_samples', {}) or {}
        if len(df) and paired_map:
            _df_paired = df['Sample'].map(
                lambda s: bool(paired_map.get(s, False)))
            df['Est_Reads'] = np.where(
                _df_paired,
                (pd.to_numeric(df['EM_Reads'], errors='coerce') * 2).round(1),
                pd.to_numeric(df['EM_Reads'], errors='coerce').round(1))
            df['Quant_Units'] = 'fragments'

            def _add_units(d):
                if not len(d):
                    return d
                d = d.copy()
                _p = d['Sample'].map(lambda s: bool(paired_map.get(s, False)))
                d['Quant_Units'] = 'fragments'
                d['Est_Reads'] = np.where(
                    _p,
                    (pd.to_numeric(d['EM_Reads'], errors='coerce') * 2).round(1),
                    pd.to_numeric(d['EM_Reads'], errors='coerce').round(1))
                return d
            confirmed = _add_units(confirmed)
            novel = _add_units(novel)

        review_cov = float(getattr(self.args, 'review_cov', 2.0))
        review_depth = float(getattr(self.args, 'review_depth', 0.1))
        review_reads = int(getattr(self.args, 'review_reads', 1))
        review_ratio = float(getattr(self.args, 'review_ratio', 0.1))
        base_r = (
            (pd.to_numeric(df['Sample_Total_Mapped'], errors='coerce') > 0) &
            (pd.to_numeric(df['Uniq_Reads'], errors='coerce') >= review_reads) &
            (pd.to_numeric(df['MeanDepth'], errors='coerce') >= review_depth)
        )
        ratio_r = np.where(is_rna, review_ratio, max(review_ratio, 0.5))
        track_a_r = (
            (pd.to_numeric(df['Coverage(%)'], errors='coerce') >= review_cov) &
            (pd.to_numeric(df['Poisson_Ratio'], errors='coerce') >= ratio_r)
        )
        review = df[base_r & track_a_r].copy()
        if len(passed):
            # polars join(how='anti') → pandas 的 indicator 反连接
            _pk = passed[['Sample', 'Accession']].drop_duplicates()
            _pk['_in_passed'] = True
            review = review.merge(_pk, on=['Sample', 'Accession'], how='left')
            review = review[review['_in_passed'].isna()].drop(
                columns=['_in_passed'])
        review = _add_seg_flag(review)
        review = _add_near_neighbor(review)
        review = _add_signal_grade(review)
        # 可选降噪：丢弃 low 分级行（--review-drop-low）。注意会连带丢掉
        # 低掺入小基因组的真信号（如 D3 的 GSyVd-1），默认关闭、仅供
        # 噪声极多的场景使用。
        if getattr(self.args, 'review_drop_low', False) and len(review):
            _n0 = len(review)
            review = review[review['Signal_Grade'] != 'low'].copy()
            self.logger.info(
                f"疑似档降噪（--review-drop-low）: 丢弃 low 信号 "
                f"{_n0 - len(review)} 行，剩余 {len(review)} 行")
        if len(review):
            review['ani_status'] = 'ani_unavailable'
            _write_tsv(review, out / 'all_viruses.review.tsv')
        self.logger.info(
            f"基础+双轨放行: {len(passed)} 行 -> 确诊 {len(confirmed)} / "
            f"疑似新种 {len(novel)} / 疑似档 {len(review)}"
            f"（review.tsv: cov≥{review_cov} depth≥{review_depth} "
            f"reads≥{review_reads} 泊松≥{review_ratio}）")

        if len(novel):
            _write_tsv(novel, out / 'all_viruses.unclassified.tsv')
        if len(confirmed):
            _write_tsv(confirmed, out / 'all_viruses.best.summary.tsv')

        # ── 发现轨自动联动标记（2026-10-18，HSVd 案例）：确诊为空说明
        # 比对/定量轨不可见（可能株系分歧或低丰度），提示/触发组装发现轨。
        marker = out / 'suggest_discovery_assembly.txt'
        if marker.exists():
            marker.unlink()
        if len(confirmed) == 0:
            marker.write_text(
                "本样本确诊为空：建议运行组装发现轨（assemble_and_classify，"
                "srna/rnaviral 模式 + kunpeng contig 分类 + BLAST 定种）。\n"
                "实证：VIROMOCK D3 的 HSVd 零 reads 映射，组装发现轨以 "
                "NODE_24234（371bp，59.6×，94.5% to NC_001351.1）回收。\n",
                encoding='utf-8')
            self.logger.warning(
                "确诊为空 -> 已写 suggest_discovery_assembly.txt："
                "建议对该样本运行组装发现轨（D3 HSVd 实证双轨互补）")

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
