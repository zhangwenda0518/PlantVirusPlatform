#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
kv_filter.py — 过滤段
====================
对鉴定段产出的表做二次阈值过滤，并施加节段病毒完整性规则。

设计来源: virome_analysis_pipeline/utils/filter_summary.py
（逻辑移植，原文件不动）

核心规则:
  1. 阈值门槛: cov / depth / reads / tpm / poisson
  2. 节段病毒（is_segmented 或多段物种）: 必须全段检出，否则整组移入 discarded
"""

from pathlib import Path

# 2026-10-20：由 polars 改为 pandas（与 kv_identify 同步，省掉 167.7MB
# 的 polars-runtime-32）。语义对照见各处注释。
import numpy as np
import pandas as pd

from .kv_identify import _write_tsv


def _find_col(df, *names):
    """按列名容错查找（去空格、大小写不敏感）"""
    lut = {str(c).strip().lower(): c for c in df.columns}
    for n in names:
        key = n.strip().lower()
        if key in lut:
            return lut[key]
    return None


def _stringify_objects(df):
    """写盘前把 object 列里残留的非字符串值统一成字符串。

    polars 版此处处理 Object dtype（其 write_csv 不支持 Object 直接抛错）。
    pandas 的 to_csv 能写任意对象，但**含 None 的 object 列**会写出空串而
    不是 'None'——这与 polars 版 map_elements(lambda v: '' if v is None
    else str(v)) 的语义一致，因此这里只需保持默认行为，函数保留为空实现
    以免调用点改动（未来若发现 dtype 问题可在此集中处理）。
    """
    return df


class FilterStage:
    def __init__(self, args, logger, out_dir, ref_info=None):
        self.args = args
        self.logger = logger
        self.out_dir = Path(out_dir)
        self.ref_info = ref_info or {}

    # ── 节段完整性规则 ────────────────────────────────────
    def _build_segment_map(self):
        """
        从 ref_info 构建 {species: {expected_segments}}。
        只有多段物种才参与全段检出规则。
        """
        seg = {}
        for acc, rec in self.ref_info.items():
            sp = rec.get('species')
            s = (rec.get('segment') or '').strip()
            if sp and s:
                seg.setdefault(sp, set()).add(s.upper())
        return {sp: s for sp, s in seg.items() if len(s) > 1}

    def apply(self, df):
        """
        返回 (kept_df, dropped_df)
        """
        out = self.out_dir / 'filter'
        out.mkdir(parents=True, exist_ok=True)

        if df is None or not len(df):
            self.logger.warning("过滤段输入为空")
            # 空输入也要落盘空表（P0-1）：否则上一轮 filtered.tsv 留在盘上
            # 被外层汇总当成本轮"通过"结果。schema 尽量沿用输入列。
            cols = list(df.columns) if df is not None else []
            empty = pd.DataFrame({c: pd.Series(dtype='object') for c in cols})
            _write_tsv(empty, out / 'filtered.tsv')
            _write_tsv(empty, out / 'discarded.tsv')
            return empty, empty.copy()

        a = self.args
        cov = _find_col(df, 'Coverage(%)', 'coverage', 'Rep_Coverage(%)')
        depth = _find_col(df, 'MeanDepth', 'depth', 'Rep_MeanDepth')
        reads = _find_col(df, 'Uniq_Reads', 'reads', 'Asm_EM_Reads')
        tpm = _find_col(df, 'TPM', 'tpm', 'Asm_TPM')
        poi = _find_col(df, 'Poisson_Ratio', 'poisson')

        missing = [n for n, c in [('Coverage', cov), ('MeanDepth', depth), ('Reads', reads)] if c is None]
        if missing:
            raise RuntimeError(f"过滤段缺少必要列: {missing}；实际列: {list(df.columns)}")

        # 逐行阈值掩码（polars 表达式 → pandas 布尔 Series）
        expr = ((pd.to_numeric(df[cov], errors='coerce') >= a.min_cov) &
                (pd.to_numeric(df[depth], errors='coerce') >= a.min_depth) &
                (pd.to_numeric(df[reads], errors='coerce') >= a.min_reads))
        if tpm and getattr(a, 'min_tpm', 0.0) > 0:
            expr = expr & (pd.to_numeric(df[tpm], errors='coerce').fillna(0)
                           >= a.min_tpm)
        if poi and getattr(a, 'min_poisson', 0.0) > 0:
            expr = expr & (pd.to_numeric(df[poi], errors='coerce').fillna(0)
                           >= a.min_poisson)

        self.logger.info(f"过滤阈值: cov>={a.min_cov} depth>={a.min_depth} reads>={a.min_reads} "
                         f"tpm>={getattr(a,'min_tpm',0)} "
                         f"poisson>={getattr(a,'min_poisson',0)}")
        pass_df = df[expr].copy()
        fail_df = df[~expr].copy()
        self.logger.info(f"阈值通过: {len(pass_df)} / {len(df)}")

        # ── 节段完整性 ──
        seg_map = self._build_segment_map()
        if not seg_map:
            self.logger.info("参考库无多段物种，跳过节段完整性规则")
        else:
            sp_col = _find_col(pass_df, 'Species')
            seg_col = _find_col(pass_df, 'Segment')
            if not sp_col or not seg_col:
                self.logger.warning(
                    "缺 Species/Segment 列，跳过节段完整性规则")
            elif 'Sample' not in pass_df.columns:
                # 没样本列时把不同样本的记录混为一组会误判“全段检出”，
                # 宁可不做这条规则。
                self.logger.warning(
                    "缺 Sample 列，无法逐样本判节段完整性，跳过该规则")
            else:
                keys = ['Sample', sp_col]
                group_rows, orphan_rows = [], []
                # polars group_by(keys) 含 null 分组 → pandas dropna=False
                for key_vals, grp in pass_df.groupby(keys, dropna=False,
                                                     sort=False):
                    sp = key_vals[1] if isinstance(key_vals, tuple) else key_vals
                    expected = seg_map.get(sp)
                    if not expected:
                        group_rows.append(grp)
                        continue
                    got = {(str(v).strip().upper() if v is not None else '')
                           for v in grp[seg_col].tolist()}
                    if expected.issubset(got):
                        group_rows.append(grp)
                    else:
                        dropped = grp.copy()
                        dropped['_drop_reason'] = (
                            'missing_segments:'
                            + ','.join(sorted(expected - got)))
                        orphan_rows.append(dropped)
                        self.logger.info(f"  节段不完整 [{sp}] 缺 {sorted(expected - got)}，剔除 {len(grp)} 段")
                # 保留原列序（polars 版用 diagonal_relaxed 合并，列集可能扩张）
                pass_df = (pd.concat(group_rows, ignore_index=True)
                           if group_rows else pass_df.iloc[0:0].copy())
                if orphan_rows:
                    _orph = pd.concat(orphan_rows, ignore_index=True)
                    _orph = _orph.drop(columns=['_drop_reason'], errors='ignore')
                    fail_df = pd.concat([fail_df, _orph], ignore_index=True)

        # 空结果也落盘（P0-1 后半）：旧实现只在非空时写 filtered/discarded，
        # 本轮全剔除/零命中时上一轮的旧表原样留在盘上，被当成"本轮通过"。
        # 现在无条件写——空表带 schema，下游读取得 0 行而不是旧数据。
        pass_df = _stringify_objects(pass_df)
        _write_tsv(pass_df, out / 'filtered.tsv')
        fail_df = _stringify_objects(fail_df)
        _write_tsv(fail_df, out / 'discarded.tsv')
        self.logger.info(f"过滤段完成: 保留 {len(pass_df)} 行 / 剔除 {len(fail_df)} 行 -> {out}")
        return pass_df, fail_df

    def run(self):
        """从磁盘读鉴定段产物再过滤"""
        src = self.out_dir / 'identify' / 'all_viruses.best.summary.tsv'
        if not src.exists():
            src = self.out_dir / 'identify' / 'all_viruses.summary.tsv'
        if not src.exists():
            self.logger.warning(f"未找到鉴定结果: {src}")
            return pd.DataFrame(), pd.DataFrame()
        self.logger.info(f"读取鉴定结果: {src}")
        # polars read_csv(ignore_errors=True) → pandas 的
        # on_bad_lines='skip'（容错跳坏行）；engine='python' 才支持该组合。
        df = pd.read_csv(str(src), sep='\t', on_bad_lines='skip',
                         engine='python')
        return self.apply(df)
