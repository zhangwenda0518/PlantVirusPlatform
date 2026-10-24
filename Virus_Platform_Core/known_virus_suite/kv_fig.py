#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
kv_fig.py — 识别/定量结果可视化段（丰度分布 / 共丰度）
================================================================
移植自 virome_analysis_pipeline/batch_plot_virus_depth.py 的
--mode freq / coabundance 两模块（--mode depth 已在 kv_plot.py；
sample 样品分布模块按需求裁掉，绘图对象 = 通过过滤（确认存在）的结果）。

与原管线的差异（平台适配）:
  - 输入 = 平台的 filter/filtered.tsv（通过过滤的 Sample×病毒记录，
    列继承自 identify/all_viruses.summary.tsv + Segment/is_segmented）；
  - 原管线的 FPKM/RPM/TPM 来自组装定量列（Asm_*），平台鉴定表没有这些列，
    按同一口径由 EM_Reads 现算：
      RPM  = EM_Reads / (Sample_Total_Mapped / 1e6)
      FPKM = RPM / (Length / 1e3)
      TPM  = (EM_Reads / Length_kb) / Σ_sample(EM_Reads / Length_kb) × 1e6
  - 产出 <run>/kvsuite/figures/*.png + *.pdf。过滤表是运行的不可变产物，
    生成一次即缓存（幂等：PNG 齐全则跳过），由 kvsuite_result 接口按需触发。
  - seaborn 0.13 box/stripplot 必须显式 dodge=False：hue=类别本身时
    auto-dodge 在有序 Categorical 下仍会把箱体挤成细条（实测）。
"""

import os
import textwrap
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from .. import sci_palette as sp
sp.apply_sci_style()   # 出版级全局样式：唯一真源

try:
    import seaborn as sns
    HAS_SEABORN = True
except ImportError:
    HAS_SEABORN = False

DPI = 200


# ── 原管线同名工具 ───────────────────────────────────────────
def _resolve_taxonomy_col(df):
    for c in ('taxonomy', 'Adjusted_Species', 'Species_NCBI',
              'Species_ICTV', 'Species'):
        if c in df.columns:
            return c
    return None


def _resolve_acc_col(df):
    for c in ('Virus', 'Rep_Accession', 'Accession'):
        if c in df.columns:
            return c
    return None


def _make_display_name(row, tax_col, acc_col):
    """tax\\n(accession) 两行显示名（同原管线）"""
    tax = str(row.get(tax_col, '')) if tax_col else ''
    acc = str(row.get(acc_col, '')) if acc_col else str(row.iloc[0])
    if tax and tax not in ('Unannotated', '-', 'nan', 'None', ''):
        return f"{tax}\n({acc})"
    return acc


def _load_summary(summary_path):
    """读鉴定/汇总表并归一列名（原管线别名规则的平台子集）。"""
    df = pd.read_csv(summary_path, sep='\t')
    df.columns = [c.strip().replace('\ufeff', '') for c in df.columns]
    for src, tgt in (('Rep_Accession', 'Accession'),
                     ('Adjusted_Species', 'Species'),
                     ('Species_NCBI', 'Species'),
                     ('Species_ICTV', 'Species')):
        if src in df.columns and tgt not in df.columns:
            df = df.rename(columns={src: tgt})
    return df


def _num(df, name):
    """取数值列；缺列返回全 NaN（老表兼容）。"""
    if name not in df.columns:
        return pd.Series(np.nan, index=df.index)
    return pd.to_numeric(df[name], errors='coerce')


def _add_derived_metrics(df):
    """补 RPM/FPKM/TPM（口径见模块 docstring；EM_Reads 定量口径）。"""
    df = df.copy()
    em = _num(df, 'EM_Reads').fillna(0.0)
    L = _num(df, 'Length')
    tot = _num(df, 'Sample_Total_Mapped')
    df['RPM'] = np.where((tot > 0) & em.notna(),
                         em / tot.where(tot > 0) * 1e6, np.nan)
    df['FPKM'] = np.where((L > 0) & df['RPM'].notna(),
                          df['RPM'] / L.where(L > 0) * 1e3, np.nan)
    tps = pd.Series(np.where(L > 0, em / L.where(L > 0) * 1e3, np.nan),
                    index=df.index)
    den = tps.groupby(df['Sample']).sum(min_count=1).rename('den')
    df = df.join(den, on='Sample')
    df['TPM'] = np.where((df['den'].fillna(0) > 0) & tps.notna(),
                         tps / df['den'].where(df['den'] > 0) * 1e6, np.nan)
    return df.drop(columns=['den'])


def _save(fig, out_dir, stem):
    """PNG + SVG + PDF 三格式落盘。

    SVG 是「结果图整理」的输入（svguest 可编辑文本与矢量元素）；
    某些 artist 不支持 SVG 时跳过该格式，不影响 PNG/PDF 主产物。
    """
    for ext in ('png', 'svg', 'pdf'):
        try:
            fig.savefig(os.path.join(out_dir, f'{stem}.{ext}'),
                        dpi=DPI, bbox_inches='tight')
        except Exception:
            if ext != 'svg':
                raise
    plt.close(fig)


# ── 模块 1: 丰度分布（综合多指标面板）─────────────────────────
def fig_freq(df, out_dir, metrics=None):
    """Comprehensive Multi-metric Virus Abundance 四联面板（箱线+散点，
    按 MeanDepth 中位数排序，log10 轴，标签带检出率 n=X/Y (Z%)）。"""
    tax_col = _resolve_taxonomy_col(df)
    acc_col = _resolve_acc_col(df)
    df = df.copy()
    df['Display_Name'] = df.apply(
        lambda r: _make_display_name(r, tax_col, acc_col), axis=1)
    df['Display_Name'] = df['Display_Name'].apply(
        lambda x: '\n'.join(textwrap.wrap(str(x), width=40)))
    total_samples = df['Sample'].nunique()
    prevalence = df.groupby('Display_Name')['Sample'].nunique()
    df['Display_Name'] = df['Display_Name'].apply(
        lambda x: f"{x}\nn={prevalence[x]}/{total_samples} "
                  f"({prevalence[x] / total_samples * 100:.0f}%)")

    metrics = [m for m in (metrics or ['MeanDepth', 'RPM', 'FPKM', 'TPM'])
               if m in df.columns]
    if not metrics:
        return 0
    unique_count = df['Display_Name'].nunique()

    if HAS_SEABORN:
        sns.set_style('ticks')
        sns.set_context('paper', font_scale=1.2)
    # 色盲安全定性色板：类目少用 colorblind，多用 tab20 循环（同原管线）
    if unique_count <= 8:
        palette = (sns.color_palette('colorblind', max(unique_count, 1))
                   if HAS_SEABORN else
                   plt.cm.tab10(np.linspace(0, 1, max(unique_count, 1))))
    else:
        palette = (sns.color_palette('tab20', max(unique_count, 1))
                   if HAS_SEABORN else
                   plt.cm.tab20(np.linspace(0, 1, max(unique_count, 1))))

    plot_data = df.melt(id_vars=['Display_Name'], value_vars=metrics,
                        var_name='Metric', value_name='Value').dropna()
    base_metric = metrics[0]
    medians = (plot_data[plot_data['Metric'] == base_metric]
               .groupby('Display_Name')['Value'].median().sort_values())
    plot_data['Display_Name'] = pd.Categorical(
        plot_data['Display_Name'], categories=medians.index, ordered=True)

    n_cols = len(metrics)
    fig, axes = plt.subplots(
        1, n_cols,
        figsize=(12 * 1.8, max(8, unique_count * 0.65)))
    if n_cols == 1:
        axes = [axes]
    for i, (ax, m) in enumerate(zip(axes, metrics)):
        mdf = plot_data[plot_data['Metric'] == m]
        if HAS_SEABORN:
            sns.boxplot(data=mdf, y='Display_Name', x='Value',
                        hue='Display_Name', ax=ax, palette=palette,
                        legend=False, width=0.6, linewidth=1.0,
                        fliersize=0, saturation=0.85, dodge=False)
            sns.stripplot(data=mdf, y='Display_Name', x='Value',
                          hue='Display_Name', ax=ax, palette=palette,
                          legend=False, size=2.5, alpha=0.55, jitter=0.2,
                          edgecolor='white', linewidth=0.2, dodge=False)
        else:
            positions = list(range(len(medians)))
            bp = ax.boxplot(
                [mdf[mdf['Display_Name'] == n]['Value'].values
                 for n in medians.index],
                positions=positions, patch_artist=True, widths=0.6,
                showfliers=False, vert=False)
            for pi, (patch, _) in enumerate(zip(bp['boxes'], medians.index)):
                patch.set_facecolor(palette[pi % len(palette)])
                patch.set_alpha(0.6)
        ax.set_title(m, fontweight='bold', fontsize=13)
        ax.set_xscale('log')
        # 数据常挤在不足 1 个数量级内（如 5.2e4~6.2e4）：默认刻度会密成一串
        # 互叠的 5.2×10⁴/5.3×10⁴/…（check_figures 体检实抓），只留主刻度又
        # 会一个数字都没有。窄跨度（<1.5 个数量级）改为首/中/尾 3 个
        # FixedLocator 数据刻度 + 普通数字格式：可读且绝无重叠。
        import math as _math
        import matplotlib.ticker as _mticker
        _vv = mdf['Value'].dropna()
        _vv = _vv[_vv > 0]
        if not _vv.empty and (_math.log10(_vv.max())
                              - _math.log10(_vv.min())) < 1.5:
            _lo, _hi = float(_vv.min()), float(_vv.max())
            _mid = 10 ** (_math.log10(_lo) / 2 + _math.log10(_hi) / 2)
            # 中位刻度取 2 位有效数字（333,333 这种原始值不体面），并夹在区间内
            _mid = max(_lo, min(_hi, round(_mid, 2 - 1
                                          - int(_math.floor(_math.log10(abs(_mid)))))))
            _ticks = [_lo, _mid, _hi]
            ax.xaxis.set_major_locator(_mticker.FixedLocator(_ticks))
            ax.xaxis.set_major_formatter(
                _mticker.FuncFormatter(lambda v, p: f'{v:,.0f}'))
        ax.xaxis.set_minor_formatter(_mticker.NullFormatter())
        if HAS_SEABORN:
            sns.despine(ax=ax)

    for i in range(1, n_cols):
        axes[i].set_ylabel('')
        axes[i].tick_params(left=False, labelleft=False)
    axes[0].set_ylabel('')
    fig.suptitle('Comprehensive Multi-metric Virus Abundance',
                 fontweight='bold', fontsize=18, y=1.01)
    _save(fig, out_dir, 'freq_multi_metrics')
    return 1


# ── 模块 2: 共丰度（TPM 热图 / 相对丰度 / 共现率 / Spearman）──
def fig_coabundance(df, out_dir, min_occurrence=2):
    tax_col = _resolve_taxonomy_col(df)
    acc_col = _resolve_acc_col(df)
    if not {'Sample', 'TPM'} <= set(df.columns):
        return 0

    label_key = [c for c in (tax_col, acc_col) if c]
    label_df = df[label_key].drop_duplicates().copy()
    label_df['Virus_Label'] = label_df.apply(
        lambda r: _make_display_name(r, tax_col, acc_col).replace('\n', ' '),
        axis=1)
    df = df.merge(label_df, on=label_key, how='left')

    tpm_pivot = df.pivot_table(index='Sample', columns='Virus_Label',
                               values='TPM', aggfunc='mean').fillna(0)
    virus_prevalence = (tpm_pivot > 0).sum(axis=0)
    keep = virus_prevalence[virus_prevalence >= min_occurrence].index.tolist()
    total_samples = tpm_pivot.shape[0]
    if len(keep) < 2 or total_samples < 2:
        return 0
    tpm_pivot = tpm_pivot[keep]
    n_viruses = len(keep)

    short = lambda cols: [  # noqa: E731
        (c.replace('\n', ' ').split('(')[-1].rstrip(')') if '(' in c
         else c[:15]) for c in cols]

    # 丰度热图 log10(TPM+1)
    tpm_display = tpm_pivot.copy()
    tpm_display.columns = short(tpm_display.columns)
    log_tpm = np.log10(tpm_display + 1)
    fig, ax = plt.subplots(
        figsize=(max(8, n_viruses * 1.2), max(5, total_samples * 0.6)))
    sns.heatmap(log_tpm, annot=True, fmt='.1f', cmap='YlGnBu', ax=ax,
                cbar_kws={'label': 'log10(TPM + 1)'})
    ax.set_title(f'Virus Abundance Heatmap\n({n_viruses} viruses, '
                 f'{total_samples} samples)', fontweight='bold', fontsize=13)
    ax.set_ylabel('Sample')
    ax.set_xlabel('')
    plt.xticks(rotation=45, ha='right', fontsize=8)
    plt.yticks(rotation=0, fontsize=9)
    _save(fig, out_dir, 'coabundance_abundance_heatmap')

    # 相对丰度堆叠柱（tab20 逐列取色：Set2 只有 8 色，>8 病毒必撞色）
    rel = tpm_pivot.div(tpm_pivot.sum(axis=1).replace(0, np.nan), axis=0) * 100
    rel = rel.fillna(0)
    rel.columns = short(rel.columns)
    fig, ax = plt.subplots(figsize=(max(10, total_samples * 0.8), 6))
    rel.plot(kind='bar', stacked=True,
             color=[plt.cm.tab20(i % 20) for i in range(len(rel.columns))],
             ax=ax, edgecolor='#333', linewidth=0.5)
    ax.set_title(f'Virus Relative Abundance Composition per Sample '
                 f'({total_samples} samples)', fontweight='bold', fontsize=13)
    ax.set_ylabel('Relative Abundance (%)', fontweight='bold')
    ax.set_xlabel('Sample')
    ax.legend(title='', bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=8)
    plt.xticks(rotation=0)
    _save(fig, out_dir, 'coabundance_relabund_stacked')

    # 共现率下三角热图（对角线 = 自身发生率）
    from itertools import combinations
    bin_mat = (tpm_pivot > 0).values.astype(int).T
    cooc = np.full((n_viruses, n_viruses), np.nan)
    annot = np.full((n_viruses, n_viruses), '', dtype=object)
    for i, j in combinations(range(n_viruses), 2):
        n_both = int(np.sum(bin_mat[i] * bin_mat[j]))
        rate = n_both / total_samples * 100
        cooc[i, j] = cooc[j, i] = rate
        annot[i, j] = annot[j, i] = f'{n_both}/{total_samples}\n{rate:.0f}%'
    for i in range(n_viruses):
        n_i = int(bin_mat[i].sum())
        cooc[i, i] = n_i / total_samples * 100
        annot[i, i] = f'{n_i}/{total_samples}\n{cooc[i, i]:.0f}%'
    s_names = short(tpm_pivot.columns)
    mask = np.triu(np.ones_like(cooc, dtype=bool), k=1)
    fig, ax = plt.subplots(
        figsize=(max(7, n_viruses * 1.0), max(5, n_viruses * 0.8)))
    sns.heatmap(cooc, annot=annot, fmt='', cmap='YlOrRd', vmin=0, vmax=100,
                mask=mask, ax=ax,
                cbar_kws={'label': 'Co-occurrence Rate (%)', 'shrink': 0.8},
                xticklabels=s_names, yticklabels=s_names,
                linewidths=0.5, linecolor='white')
    ax.set_title(f'Virus Co-occurrence Rate\n({n_viruses} viruses, '
                 f'{total_samples} samples)', fontweight='bold', fontsize=12)
    ax.set_xlabel('')
    ax.set_ylabel('')
    plt.xticks(rotation=45, ha='right', fontsize=8)
    plt.yticks(rotation=0, fontsize=8)
    _save(fig, out_dir, 'coabundance_cooccurrence_heatmap')

    # Spearman 丰度相关性热图（下三角）
    corr = tpm_pivot.corr(method='spearman')
    mask = np.triu(np.ones_like(corr, dtype=bool), k=1)
    fig, ax = plt.subplots(
        figsize=(max(6, n_viruses * 0.8), max(5, n_viruses * 0.7)))
    sns.heatmap(corr, annot=True, fmt='.2f', cmap='vlag', vmin=-1, vmax=1,
                center=0, mask=mask, ax=ax,
                cbar_kws={'label': 'Spearman ρ', 'shrink': 0.8},
                xticklabels=s_names, yticklabels=s_names)
    ax.set_title(f'Spearman Correlation between Virus Species\n'
                 f'({n_viruses} viruses, {total_samples} samples)',
                 fontweight='bold', fontsize=12)
    ax.set_xlabel('')
    ax.set_ylabel('')
    plt.xticks(rotation=45, ha='right', fontsize=8)
    plt.yticks(rotation=0, fontsize=8)
    _save(fig, out_dir, 'coabundance_spearman_heatmap')
    return 4


# ── 入口 ────────────────────────────────────────────────────
FREQ_PNG = 'freq_multi_metrics.png'
COAB_PNGS = ('coabundance_abundance_heatmap.png',
             'coabundance_relabund_stacked.png',
             'coabundance_cooccurrence_heatmap.png',
             'coabundance_spearman_heatmap.png')


def generate_figures(filtered_tsv, out_dir, logger=None):
    """按参考管线 freq / coabundance 两模块出图（绘图对象 = 通过过滤的
    结果），返回 PNG 文件名列表（含既有缓存）。

    幂等：全部 PNG + SVG 已在（此前生成过）直接返回；共丰度图在 <2 样品
    或 <2 个多样品病毒时按原管线语义跳过。SVG 一并作为幂等条件，旧运行
    首次访问时会补生成一份可编辑矢量图。
    """
    out_dir = Path(out_dir)
    have = {p.name for p in out_dir.glob('*.png')}
    have_svg = {p.name for p in out_dir.glob('*.svg')}
    need_svg = {f'{n[:-4]}.svg' for n in (FREQ_PNG, *COAB_PNGS)}
    if FREQ_PNG in have and set(COAB_PNGS) <= have and need_svg <= have_svg:
        return sorted(have)
    if not os.path.isfile(str(filtered_tsv)):
        return sorted(have)

    df = _load_summary(filtered_tsv)
    need = {'Sample', 'Accession'}
    if not need <= set(df.columns) or df.empty:
        return sorted(have)
    df = _add_derived_metrics(df)
    out_dir.mkdir(parents=True, exist_ok=True)

    n = 0
    try:
        n += fig_freq(df, out_dir)
    except Exception as e:  # noqa: BLE001
        if logger:
            logger.warning(f'  丰度分布图失败: {e}')
    try:
        n += fig_coabundance(df, out_dir)
    except Exception as e:  # noqa: BLE001
        if logger:
            logger.warning(f'  共丰度图失败: {e}')
    if logger:
        logger.info(f'  识别/定量图表: 新增 {n} 张 -> {out_dir}')
    return sorted(p.name for p in out_dir.glob('*.png'))
