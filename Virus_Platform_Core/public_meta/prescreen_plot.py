#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
prescreen_plot.py — SRA / GSA 病毒样本预筛结果可视化（四图组）

读 prescreen_engine.py 的产物 TSV，出 4 张 SCI 图（PNG 600dpi + PDF 矢量），
输出到 <预筛输出目录>/figures/：

    F1 prescreen_viral_pct_hist   病毒占比分布直方图（分来源 + 阈值线）
    F2 prescreen_depth_vs_viral   测序量 × 病毒占比散点（归一化口径 virus_pct_ident）
    F3 prescreen_top_viral_taxa   Top-N 检出病毒类群条形图
    F4 prescreen_funnel_status    预筛漏斗 + 未决状态构成

口径纪律（与引擎一致，不得混淆）：
    virus_pct        各来源自身口径（STAT=/analysed；GSA=/总 reads）——阈值
                     （PASS 判定）作用在这一列；
    virus_pct_ident  病毒命中/已鉴定命中 —— 跨来源可比，散点图用这一列，
                     阈值线不画在它上面（分母不同，画了就是错误引导）。

用法：
    python -m Virus_Platform_Core.public_meta.prescreen_plot <outdir> [--threshold 0.5]
    python -m Virus_Platform_Core.public_meta.prescreen_plot <outdir> --figdir D:/figs
"""

from __future__ import annotations

import argparse
import os
import sys

if sys.platform == "win32" and sys.stdout is not None:
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

import numpy as np
import pandas as pd
import textwrap

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator, FuncFormatter

try:                                    # 包内 import
    from .. import sci_palette as sp
except ImportError:                     # 独立脚本兜底
    _ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    if _ROOT not in sys.path:
        sys.path.insert(0, _ROOT)
    from Virus_Platform_Core import sci_palette as sp

sp.apply_sci_style()

_GREY = sp.OKABE["grey"]
_BLUE = sp.OKABE["blue"]
_VERM = sp.OKABE["vermillion"]
_THR_COLOR = sp.OKABE["vermillion"]
_WRAP_W = 28


# ------------------------------------------------------------------
# 读取与小工具
# ------------------------------------------------------------------
def _num(s: pd.Series) -> pd.Series:
    return pd.to_numeric(s, errors="coerce")


def load_summary(outdir: str) -> pd.DataFrame:
    p = os.path.join(outdir, "stat_run_summary.tsv")
    if not os.path.exists(p):
        raise FileNotFoundError(f"找不到 {p} —— 请先跑 prescreen_engine")
    df = pd.read_csv(p, sep="\t", dtype=str, keep_default_na=False)
    for c in ("virus_pct", "virus_pct_ident", "virus_count", "analysed",
              "identified", "identified_pct", "top_virus_pct", "fetch_s"):
        if c in df.columns:
            df[c] = _num(df[c])
    return df


def _wrap(s: str, width: int = _WRAP_W) -> str:
    return textwrap.fill(str(s), width=width)


def _save(fig, out_dir: str, stem: str) -> list:
    """PNG + SVG + PDF 三格式落盘（SVG 供结果图二次排版）。"""
    paths = []
    for ext, kw in (("png", dict(dpi=600)), ("svg", {}), ("pdf", {})):
        fp = os.path.join(out_dir, f"{stem}.{ext}")
        try:
            fig.savefig(fp, bbox_inches="tight", **kw)
            paths.append(fp)
        except Exception:
            if ext != "svg":
                raise
    plt.close(fig)
    return paths


def _note_empty(ax, msg: str = "No valid data") -> None:
    ax.text(0.5, 0.5, msg, ha="center", va="center",
            transform=ax.transAxes, color=_GREY, fontsize=11)
    ax.set_xticks([])
    ax.set_yticks([])


def _kfmt() -> FuncFormatter:
    def f(v, _):
        if v >= 1_000_000:
            return f"{v/1_000_000:g}M"
        if v >= 1_000:
            return f"{v/1_000:g}k"
        return f"{v:g}"
    return FuncFormatter(f)


# ------------------------------------------------------------------
# F1 病毒占比分布直方图（virus_pct，各来源自身口径 + 阈值线）
# ------------------------------------------------------------------
def plot_viral_pct_hist(df: pd.DataFrame, threshold: float | None,
                        out_dir: str) -> list:
    ok = df[df["status"] == "OK"].copy()
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    if ok.empty or not (ok["virus_pct"] > 0).any():
        _note_empty(ax, "No runs with taxonomy data (status=OK empty)")
        return _save(fig, out_dir, "prescreen_viral_pct_hist")

    vmax = float(np.nanpercentile(ok.loc[ok["virus_pct"] > 0, "virus_pct"], 99))
    hi = max(0.5, min(5.0, vmax * 1.15))
    bins = np.linspace(0.0, hi, 61)
    for src, color in (("NCBI_STAT", _BLUE), ("GSA", _VERM)):
        v = ok.loc[(ok["source"] == src) & (ok["virus_pct"] > 0), "virus_pct"]
        if len(v):
            ax.hist(v.clip(upper=hi), bins=bins, histtype="stepfilled",
                    alpha=0.55, color=color, label=f"{src} (n={len(v)})",
                    edgecolor="white", linewidth=0.4)
    n_clip = int((ok["virus_pct"] > hi).sum())
    if threshold is not None and 0 < threshold < hi:
        ax.axvline(threshold, color=_THR_COLOR, ls="--", lw=1.4,
                   label=f"Threshold {threshold:g}% (PASS)")
    ax.set_yscale("log")
    ax.set_xlabel(f"Viral read fraction virus_pct (%, view 0-{hi:.1f})")
    ax.set_ylabel("Runs (log scale)")
    ax.set_title("Prescreen: viral fraction distribution", fontsize=11)
    if n_clip:
        ax.text(0.99, 0.97, f"{n_clip} runs > {hi:.1f}% not shown",
                transform=ax.transAxes, ha="right", va="top",
                fontsize=8, color=_GREY)
    ax.legend(frameon=False, fontsize=8)
    sp.despine()
    return _save(fig, out_dir, "prescreen_viral_pct_hist")


# ------------------------------------------------------------------
# F2 测序量 × 病毒占比散点（virus_pct_ident 跨来源可比口径）
# ------------------------------------------------------------------
def plot_depth_vs_viral(df: pd.DataFrame, out_dir: str) -> list:
    ok = df[df["status"] == "OK"].copy()
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    pos = ok[(ok["virus_pct_ident"] > 0) & (ok["analysed"] > 0)]
    if ok.empty:
        _note_empty(ax)
        return _save(fig, out_dir, "prescreen_depth_vs_viral")

    ylo = 10 ** (np.floor(np.log10(max(pos["virus_pct_ident"].min(), 1e-4))) - 0.4) \
        if len(pos) else 1e-3
    zero = ok[(ok["virus_pct_ident"].fillna(0) <= 0) & (ok["analysed"] > 0)]
    if len(zero):
        ax.scatter(zero["analysed"], np.full(len(zero), ylo), s=12,
                   color=_GREY, alpha=0.25, linewidths=0, zorder=1,
                   label=f"No viral signal (n={len(zero)}, at floor)")
    for src, color, mark in (("NCBI_STAT", _BLUE, "o"), ("GSA", _VERM, "^")):
        sub = pos[pos["source"] == src]
        if len(sub):
            ax.scatter(sub["analysed"], sub["virus_pct_ident"], s=26,
                       color=color, alpha=0.75, marker=mark,
                       edgecolors="white", linewidths=0.4, zorder=3,
                       label=f"{src}（n={len(sub)}）")
    passed = pos[pos.get("pass", "") == "PASS"] if "pass" in pos.columns else pos.iloc[0:0]
    if len(passed):
        ax.scatter(passed["analysed"], passed["virus_pct_ident"], s=64,
                   facecolors="none", edgecolors="black", linewidths=0.9,
                   marker="o", zorder=4, label=f"PASS (n={len(passed)})")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylim(ylo * 0.8, None)
    ax.xaxis.set_major_locator(LogLocator(base=10))
    ax.yaxis.set_major_formatter(_kfmt())
    ax.set_xlabel("Analysed reads (log scale)")
    ax.set_ylabel("Viral / identified hits, virus_pct_ident (%)")
    # pad 给上方两行图例让位，避免标题与图例同带重叠
    ax.set_title("Sequencing depth vs viral signal (normalised)", fontsize=11, pad=34)
    if len(pos):
        # Top-5 标注：按位次交错垂直偏移，避免同簇标签互相压盖
        top = pos.nlargest(min(5, len(pos)), "virus_pct_ident").reset_index()
        for k, r in top.iterrows():
            ax.annotate(str(r["acc"]), (r["analysed"], r["virus_pct_ident"]),
                        textcoords="offset points",
                        xytext=(6, 4 + 8 * (k % 3) - (6 if k % 3 == 2 else 0)),
                        fontsize=6.5, color="black", alpha=0.85, zorder=5)
    # 图例固定画在坐标区上方外侧（两列），避免 best 自适应时压住数据点
    ax.legend(frameon=False, fontsize=8, loc="lower left",
              bbox_to_anchor=(0.0, 1.01), ncol=2, borderaxespad=0.0,
              columnspacing=1.2, handlelength=1.6)
    sp.despine()
    return _save(fig, out_dir, "prescreen_depth_vs_viral")


# ------------------------------------------------------------------
# F3 Top-N 检出病毒类群条形图（stat_viral_taxa 长表聚合）
# ------------------------------------------------------------------
def plot_top_viral_taxa(outdir: str, out_dir: str, topn: int = 20) -> list:
    p = os.path.join(outdir, "stat_viral_taxa.tsv")
    fig, ax = plt.subplots(figsize=(7.2, max(3.2, 0.32 * topn + 1.2)))
    if not os.path.exists(p):
        _note_empty(ax, "stat_viral_taxa.tsv not found")
        return _save(fig, out_dir, "prescreen_top_viral_taxa")
    vt = pd.read_csv(p, sep="\t", dtype=str, keep_default_na=False)
    if vt.empty:
        _note_empty(ax, "Long table empty (no virus detected)")
        return _save(fig, out_dir, "prescreen_top_viral_taxa")

    cnt = (vt.groupby("org")["acc"].nunique().sort_values(ascending=True)
             .tail(topn))
    labels = [_wrap(o) for o in cnt.index]
    colors = [_GREY if str(o).startswith("Viruses") else _BLUE for o in cnt.index]
    ax.barh(labels, cnt.values, color=colors, edgecolor="white", linewidth=0.4)
    for i, v in enumerate(cnt.values):
        ax.text(v, i, f" {v}", va="center", fontsize=8, color="black")
    ax.set_xlabel(f"Detected runs (top {len(cnt)}, unique)")
    ax.set_title("Prescreen: top detected viral taxa", fontsize=11)
    ax.tick_params(axis="y", labelsize=8)
    ax.xaxis.set_major_formatter(_kfmt())
    sp.despine()
    return _save(fig, out_dir, "prescreen_top_viral_taxa")


# ------------------------------------------------------------------
# F4 预筛漏斗 + 未决状态构成
# ------------------------------------------------------------------
def plot_funnel_status(df: pd.DataFrame, out_dir: str,
                       threshold: float | None) -> list:
    n_all = len(df)
    ok = df[df["status"] == "OK"]
    n_ok = len(ok)
    n_sig = int((ok["virus_pct"].fillna(0) > 0).sum())
    pass_col = df.get("pass", pd.Series(dtype=str))
    n_pass = int((pass_col == "PASS").sum())

    fig = plt.figure(figsize=(9.6, 4.2))
    gs = fig.add_gridspec(1, 2, width_ratios=[3, 2], wspace=0.42)
    axf = fig.add_subplot(gs[0, 0])
    axu = fig.add_subplot(gs[0, 1])

    if n_all == 0:
        _note_empty(axf)
        _note_empty(axu)
        return _save(fig, out_dir, "prescreen_funnel_status")

    stages = [("Candidate runs (input)", n_all, _GREY),
              ("With taxonomy data (OK)", n_ok, _BLUE),
              ("Virus detected (>0)", n_sig, _VERM)]
    if n_pass or threshold is not None:
        stages.append((f"PASS ({threshold:g}%)" if threshold else "PASS",
                       n_pass, sp.OKABE["green"]))
    y = np.arange(len(stages))[::-1]
    widths = [max(n, 0.35) for _, n, _ in stages]     # 0 层也给可见最小宽
    axf.barh(y, widths, color=[c for _, _, c in stages],
             edgecolor="white", linewidth=0.5, height=0.62)
    for yi, (_, n, _) in zip(y, stages):
        pct = f"（{n / n_all * 100:.1f}%）" if n_all else ""
        axf.text(max(n, 0.35), yi, f" {n}{pct}", va="center", fontsize=8.5)
    axf.set_yticks(y)
    axf.set_yticklabels([s for s, _, _ in stages], fontsize=9)
    axf.set_xscale("log")
    axf.set_xlabel("Runs (log scale)")
    axf.set_title("Prescreen funnel", fontsize=11)
    axf.xaxis.set_major_formatter(_kfmt())
    sp.despine()

    unres = df[df["status"] != "OK"]
    if len(unres):
        stat = (unres["source"] + "：" + unres["status"].str.split().str[0])
        vc = stat.value_counts().sort_values()
        ulabels = [_wrap(s, 22) for s in vc.index]      # 逐条换行，不能整列喂 textwrap
        axu.barh(ulabels, vc.values, color=_GREY, height=0.6)
        for i, v in enumerate(vc.values):
            axu.text(v, i, f" {v}", va="center", fontsize=8)
        axu.set_xlabel("Runs")
        axu.set_title("Unresolved (no taxonomy data)", fontsize=11)
        axu.tick_params(axis="y", labelsize=8)
        axu.xaxis.set_major_formatter(_kfmt())
    else:
        _note_empty(axu, "All runs have taxonomy data")
    sp.despine(ax=axu)
    return _save(fig, out_dir, "prescreen_funnel_status")


# ------------------------------------------------------------------
# 入口
# ------------------------------------------------------------------
def plot_prescreen(outdir: str, threshold: float | None = None,
                   figdir: str | None = None) -> list:
    """读 <outdir> 下的预筛 TSV，出 4 张图；返回图片路径列表。

    threshold：PASS 阈值（%，与 prescreen_engine --virus-pct 同口径）；
    不传时尝试从 stat_run_summary.tsv 的 pass 列推断是否画阈值层。
    """
    figdir = figdir or os.path.join(outdir, "figures")
    os.makedirs(figdir, exist_ok=True)
    df = load_summary(outdir)
    if threshold is None and "pass" in df.columns and (df["pass"] == "PASS").any():
        # 反推阈值：取 PASS 行里最小的 virus_pct（引擎判据即作用在该列）
        p = df.loc[df["pass"] == "PASS", "virus_pct"]
        threshold = float(np.nanmin(p)) if len(p) else None

    paths: list = []
    paths += plot_viral_pct_hist(df, threshold, figdir)
    paths += plot_depth_vs_viral(df, figdir)
    paths += plot_top_viral_taxa(outdir, figdir)
    paths += plot_funnel_status(df, figdir, threshold)
    return paths


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="SRA/GSA 病毒样本预筛结果可视化")
    ap.add_argument("outdir", help="prescreen_engine 的输出目录（含 stat_run_summary.tsv）")
    ap.add_argument("--threshold", type=float, default=None,
                    help="PASS 阈值 %%（与 --virus-pct 同口径；缺省从 pass 列反推）")
    ap.add_argument("--figdir", default=None, help="图片输出目录（缺省 <outdir>/figures）")
    a = ap.parse_args(argv)
    figs = plot_prescreen(a.outdir, threshold=a.threshold, figdir=a.figdir)
    for fp in figs:
        print(f"  图 -> {fp}")
    print(f"完成：{len(figs)} 个文件")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
