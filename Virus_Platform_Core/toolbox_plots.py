#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
toolbox_plots.py — 专项工具结果可视化补全（内容感知，纯附加层）

背景（2026-10-12 审计）：以下工具此前只落表格/JSON，结果面板没有配图——
引物设计、CRISPR 设计、dsRNA 设计、结构比较(structcmp)、RDP5 重组、
重组预筛(align_qc)、候选验证分类、miRNA 靶标、MOT 迁移矩阵、元数据治理。

本模块不改变任何工具的产物：run 目录扫到已知表格产物时，就地补出
配图（PNG 300dpi + PDF），工具结果面板按目录扫描列产物，图自动出现。
所有函数防御式实现：列名模糊匹配、缺列/空表/解析失败一律跳过，绝不抛错
拖垮工具任务本身（与 io_utils.safe_plot 同纪律）。

入口：add_plots(run_dir, logger=None) -> list[str]（写出的图片路径）
"""

from __future__ import annotations

import glob
from collections import Counter
import json
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

try:                                    # 包内
    from . import sci_palette as sp
except ImportError:                     # 独立兜底
    _ROOT = os.path.dirname(os.path.abspath(__file__))
    if _ROOT not in sys.path:
        sys.path.insert(0, _ROOT)
    from Virus_Platform_Core import sci_palette as sp

sp.apply_sci_style()

_GREY = sp.OKABE["grey"]
_BLUE = sp.OKABE["blue"]
_VERM = sp.OKABE["vermillion"]


# ------------------------------------------------------------------
# 公共小件
# ------------------------------------------------------------------
def _find(run_dir: str, *patterns: str) -> str | None:
    for pat in patterns:
        hits = sorted(glob.glob(os.path.join(run_dir, pat)))
        if hits:
            return hits[0]
    return None


def _read_tab(path: str) -> pd.DataFrame:
    sep = "\t" if path.lower().endswith((".tsv", ".txt")) else ","
    return pd.read_csv(path, sep=sep, dtype=str, keep_default_na=False)


def _num(df: pd.DataFrame, *keys: str) -> pd.Series | None:
    """按子串模糊找数值列（不区分大小写）。"""
    for c in df.columns:
        lc = c.lower()
        if any(k in lc for k in keys):
            v = pd.to_numeric(df[c], errors="coerce")
            if v.notna().sum() >= 2:
                return v
    return None


def _is_text(s: pd.Series) -> bool:
    """pandas ≥3 的字符串列是专用 str dtype，不再等于 object——按
    「非数值且非布尔」判文本，两个版本通吃。"""
    try:
        return not (pd.api.types.is_numeric_dtype(s)
                    or pd.api.types.is_bool_dtype(s))
    except (TypeError, ValueError):
        return False


def _text(df: pd.DataFrame, *keys: str, min_uniq: int = 2,
          max_uniq: int = 200) -> pd.Series | None:
    for c in df.columns:
        lc = c.lower()
        if any(k in lc for k in keys) and _is_text(df[c]):
            u = df[c].nunique()
            if min_uniq <= u <= max_uniq:
                return df[c]
    return None


def _save(fig, run_dir: str, stem: str) -> list:
    """PNG + SVG + PDF 三格式落盘（SVG 供结果图二次排版）。"""
    paths = []
    for ext, kw in (("png", dict(dpi=300)), ("svg", {}), ("pdf", {})):
        fp = os.path.join(run_dir, f"{stem}.{ext}")
        try:
            fig.savefig(fp, bbox_inches="tight", **kw)
            paths.append(fp)
        except Exception:
            pass
    plt.close(fig)
    return paths


def _wrap(s, width: int = 26) -> str:
    return textwrap.fill(str(s), width=width)


def _multi_hist(ax, series: pd.Series, title: str, xlabel: str,
                color: str) -> None:
    v = series.dropna()
    if len(v) == 0:
        ax.text(0.5, 0.5, "No data", ha="center", va="center",
                transform=ax.transAxes, color=_GREY)
        ax.set_title(title, fontsize=9)
        return
    ax.hist(v, bins=min(30, max(8, int(np.sqrt(len(v))) * 2)),
            color=color, edgecolor="white", linewidth=0.4)
    ax.axvline(float(v.median()), color="black", ls="--", lw=0.9,
               label=f"median {v.median():.2f}")
    ax.legend(frameon=True, framealpha=0.85, edgecolor="none", fontsize=7)
    ax.set_title(title, fontsize=9)
    ax.set_xlabel(xlabel, fontsize=8)
    # 剪掉两端最外侧刻度：多面板并排时，左面板最右刻度与右面板最左刻度
    # 都悬在面板边界外的同一行（tight_layout 的横向间距不算底部刻度的
    # 横向外溢），必相撞——体检实抓 "65"×"-0.5" 161px²（2026-10-13）。
    ax.xaxis.set_major_locator(
        matplotlib.ticker.MaxNLocator(nbins=7, prune='both'))


# ------------------------------------------------------------------
# 1. structcmp：同一性矩阵热图
# ------------------------------------------------------------------
def plot_identity_heatmap(run_dir: str) -> list:
    p = _find(run_dir, "identity_matrix.tsv")
    if not p:
        return []
    df = _read_tab(p)
    names = list(df.columns[1:])
    m = df.iloc[:, 1:].apply(pd.to_numeric, errors="coerce").to_numpy()
    if m.size == 0 or not np.isfinite(m).any():
        return []
    n = len(names)
    size = max(4.0, min(12.0, 0.42 * n + 2.4))
    fig, ax = plt.subplots(figsize=(size, size * 0.86))
    im = ax.imshow(m, cmap="viridis", vmin=np.nanmin(m), vmax=np.nanmax(m))
    ax.set_xticks(range(n)); ax.set_yticks(range(n))
    rot = 45 if n >= 4 else 0
    ax.set_xticklabels([_wrap(x, 14) for x in names],
                       rotation=rot, ha="right" if rot else "center", fontsize=7)
    ax.set_yticklabels([_wrap(x, 22) for x in names], fontsize=7)
    if n <= 14:
        for i in range(n):
            for j in range(n):
                ax.text(j, i, f"{m[i, j]:.1f}", ha="center", va="center",
                        fontsize=6, color="white" if m[i, j] < np.nanmean(m) else "black")
    fig.colorbar(im, ax=ax, shrink=0.8, label="Identity (%)")
    ax.set_title("Pairwise identity matrix (structcmp)", fontsize=11)
    return _save(fig, run_dir, "identity_matrix_heatmap")


# ------------------------------------------------------------------
# 2. RDP5 重组事件概览
# ------------------------------------------------------------------
def plot_rdp_events(run_dir: str) -> list:
    p = _find(run_dir, "events.tsv")
    if not p:
        return []
    df = _read_tab(p)
    if df.empty:
        return []
    startc = next((c for c in df.columns
                   if any(k in c.lower() for k in ("start", "begin", "left"))), None)
    endc = next((c for c in df.columns
                 if any(k in c.lower() for k in ("end", "stop", "right"))), None)
    namec = next((c for c in df.columns
                  if _is_text(df[c]) and df[c].nunique() >= 2
                  and not any(k in c.lower() for k in ("start", "end", "begin", "stop",
                                                       "left", "right", "pval", "p-"))), None)
    if not (startc and endc and namec):
        return []
    df = df.copy()
    df["_s"] = pd.to_numeric(df[startc], errors="coerce")
    df["_e"] = pd.to_numeric(df[endc], errors="coerce")
    df = df.dropna(subset=["_s", "_e"])
    if df.empty:
        return []
    counts = df.groupby(namec).size().sort_values(ascending=True).tail(30)
    fig, (ax1, ax2) = plt.subplots(
        1, 2, figsize=(11.5, max(3.0, 0.30 * len(counts) + 1.2)),
        gridspec_kw={"width_ratios": [1, 2.1], "wspace": 0.42})
    ax1.barh([_wrap(x, 20) for x in counts.index], counts.values,
             color=_BLUE, edgecolor="white", linewidth=0.4)
    for i, v in enumerate(counts.values):
        ax1.text(v, i, f" {v}", va="center", fontsize=8)
    ax1.set_xlabel("Recombination events")
    ax1.set_title("Events per sequence", fontsize=10)
    sp.despine()
    # 断点区间条带
    span = float(df[["_s", "_e"]].to_numpy().max()) or 1.0
    ymap = {name: i for i, name in enumerate(counts.index)}
    for _, r in df.iterrows():
        if r[namec] in ymap:
            ax2.broken_barh([(min(r["_s"], r["_e"]),
                              abs(r["_e"] - r["_s"]) or span * 0.002)],
                            (ymap[r[namec]] - 0.32, 0.64),
                            facecolors=_VERM, edgecolors="none", alpha=0.75)
    ax2.set_yticks(range(len(counts)))
    ax2.set_yticklabels([_wrap(x, 20) for x in counts.index], fontsize=7)
    ax2.set_xlim(0, span * 1.02)
    ax2.set_xlabel("Genome position (recombinant segment)")
    ax2.set_title("Breakpoint intervals", fontsize=10)
    sp.despine()
    return _save(fig, run_dir, "rdp_events_overview")


# ------------------------------------------------------------------
# 3. align_qc 重组预筛：QC 四联分布
# ------------------------------------------------------------------
def plot_align_qc(run_dir: str) -> list:
    p = _find(run_dir, "align_qc_report.tsv")
    if not p:
        return []
    df = _read_tab(p)
    if df.empty:
        return []
    cols = []
    for keys, title, xl in ((("length", "len"), "Sequence length", "bp"),
                            (("gap",), "Gap %", "%"),
                            (("degener", "ambig", "n_"), "Degenerate/ambiguous %", "%"),
                            (("ident", "pid"), "Pairwise identity", "%")):
        v = _num(df, *keys)
        if v is not None and v.notna().sum() >= 3:
            cols.append((v, title, xl))
    if not cols:
        return []
    removed = _text(df, "removed", "verdict", "pass", "status", "flag")
    fig, axes = plt.subplots(1, len(cols), figsize=(3.1 * len(cols), 3.4))
    if len(cols) == 1:
        axes = [axes]
    for ax, (v, title, xl) in zip(axes, cols):
        _multi_hist(ax, v, title, xl, _BLUE)
        if removed is not None:
            bad = v[removed.astype(str).str.strip().str.lower()
                    .isin(("true", "yes", "removed", "fail", "failed", "1"))]
            if len(bad):
                ax.hist(bad, bins=min(20, max(6, len(bad))), color=_VERM,
                        alpha=0.6, label=f"removed {len(bad)}")
                ax.legend(frameon=True, framealpha=0.85, edgecolor="none", fontsize=7)
    fig.suptitle("Alignment prescreen QC (align_qc)", fontsize=11)
    # w_pad 拉开面板间距：多面板直方图的相邻 x 刻度（如左图 "8" 与
    # 右图 "-1"）会贴到同一像素列（check_figures 实抓 80px² 重叠）
    fig.tight_layout(rect=(0, 0, 1, 0.94), w_pad=2.2)
    return _save(fig, run_dir, "align_qc_overview")


# ------------------------------------------------------------------
# 4. 引物设计：Tm / GC / 产物大小 + 每靶标对数
# ------------------------------------------------------------------
def plot_primers(run_dir: str) -> list:
    p = _find(run_dir, "primers*.tsv", "primer*.tsv")
    if not p:
        return []
    df = _read_tab(p)
    if df.empty:
        return []
    panels = []
    for keys, title, xl, color in ((("tm", "t_m"), "Melting temp Tm", "°C", _BLUE),
                                   (("gc",), "GC content", "%", _VERM),
                                   (("size", "length", "amplicon"), "Product/target length", "bp",
                                    sp.OKABE["green"]),
                                   (("penalty", "score", "specific"), "Penalty/score", "", _GREY)):
        v = _num(df, *keys)
        if v is not None and v.notna().sum() >= 3:
            panels.append((v, title, xl, color))
    target = _text(df, "target", "contig", "seq", "region", "group", "name",
                   min_uniq=1, max_uniq=60)
    if not panels and target is None:
        return []
    nrow = 2 if target is not None else 1
    ncol = max(1, len(panels))
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.1 * ncol, 3.0 * nrow),
                             squeeze=False)
    for j, (v, title, xl, color) in enumerate(panels):
        _multi_hist(axes[0][j], v, title, xl, color)
    if target is not None:
        tc = target.value_counts().sort_values(ascending=True).tail(15)
        # 第二行整行合并给「每靶标引物对数」，避免右侧留大片空白
        # （注意：第二行原有轴要全部移除，否则与新轴同槽叠加）
        gs = axes[1][0].get_gridspec()
        for jj in range(ncol):
            axes[1][jj].remove()
        ax = fig.add_subplot(gs[1, :])
        ax.barh([_wrap(x, 22) for x in tc.index], tc.values, color=_BLUE,
                edgecolor="white", linewidth=0.4)
        for i, v in enumerate(tc.values):
            ax.text(v, i, f" {v}", va="center", fontsize=8)
        ax.set_title("Primer pairs per target", fontsize=9)
        ax.set_xlabel("Pairs")
        sp.despine()
    fig.suptitle("Primer design overview", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.95), w_pad=2.2)
    return _save(fig, run_dir, "primer_design_overview")


# ------------------------------------------------------------------
# 5. CRISPR 设计：候选评分 + Top 组合
# ------------------------------------------------------------------
def plot_crispr(run_dir: str) -> list:
    p = _find(run_dir, "crrna_candidates.tsv")
    if not p:
        return []
    df = _read_tab(p)
    if df.empty:
        return []
    score = _num(df, "score", "efficiency", "rating", "rank_value")
    gc = _num(df, "gc")
    offt = _num(df, "offtarget", "off_target", "mismatch")
    panels = [(v, t, x, c) for v, t, x, c in (
        (score, "crRNA candidate score", "Score", _BLUE),
        (gc, "GC content", "%", _VERM),
        (offt, "Off-target/mismatch count", "count", _GREY)) if v is not None and v.notna().sum() >= 3]
    if not panels:
        return []
    fig, axes = plt.subplots(1, len(panels), figsize=(3.3 * len(panels), 3.4))
    if len(panels) == 1:
        axes = [axes]
    for ax, (v, t, x, c) in zip(axes, panels):
        _multi_hist(ax, v, t, x, c)
    fig.suptitle("CRISPR assay design overview", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.94), w_pad=2.2)
    paths = _save(fig, run_dir, "crispr_design_overview")
    # Top 组合条形图（ranked_sets.tsv）
    p2 = _find(run_dir, "ranked_sets.tsv")
    if p2:
        df2 = _read_tab(p2)
        s2 = _num(df2, "score", "rating") if not df2.empty else None
        id2 = _text(df2, "set", "id", "name", min_uniq=1, max_uniq=40) \
            if not df2.empty else None
        if s2 is not None and id2 is not None:
            d = pd.DataFrame({"id": id2, "s": s2}).dropna()
            d = d.reindex(d["s"].abs().sort_values(ascending=True).index).tail(15)
            fig2, ax = plt.subplots(figsize=(6.4, max(2.8, 0.28 * len(d) + 1.0)))
            ax.barh([_wrap(x, 26) for x in d["id"]], d["s"].abs(),
                    color=_BLUE, edgecolor="white", linewidth=0.4)
            for i, v in enumerate(d["s"].abs()):
                ax.text(v, i, f" {v:.3g}", va="center", fontsize=8)
            ax.set_xlabel("Set score (absolute)")
            ax.set_title("Top guide sets (ranked_sets)", fontsize=10)
            sp.despine()
            paths += _save(fig2, run_dir, "crispr_top_sets")
    return paths


# ------------------------------------------------------------------
# 6. miRNA 靶标：引擎共识 + Top 靶病毒
# ------------------------------------------------------------------
def plot_mirna(run_dir: str) -> list:
    p = _find(run_dir, "mirna_targets_merged.csv")
    if not p:
        return []
    df = _read_tab(p)
    if df.empty:
        return []
    nalg = _num(df, "num_algorithms")
    energy = _num(df, "best_energy")
    tv = _text(df, "target_virus_name", min_uniq=1, max_uniq=80)
    if tv is None:
        tv = _text(df, "target_virus_id", "target_virus", min_uniq=1, max_uniq=80)
    mi = _text(df, "mirna_id", min_uniq=1, max_uniq=80)
    npanel = 1 + (1 if tv is not None else 0) + (1 if energy is not None else 0)
    fig, axes = plt.subplots(1, npanel, figsize=(3.5 * npanel, 3.6))
    if npanel == 1:
        axes = [axes]
    k = 0
    if nalg is not None and nalg.notna().any():
        vc = nalg.dropna().astype(int).value_counts().sort_index()
        axes[k].bar(vc.index.astype(str), vc.values, color=_BLUE,
                    edgecolor="white", linewidth=0.4)
        axes[k].set_xlabel("Supporting engines")
        axes[k].set_ylabel("Hits")
        axes[k].set_title("Consensus strength", fontsize=9)
        k += 1
    if tv is not None:
        tc = tv.value_counts().sort_values(ascending=True).tail(12)
        axes[k].barh([_wrap(x, 24) for x in tc.index], tc.values,
                     color=_VERM, edgecolor="white", linewidth=0.4)
        for i, v in enumerate(tc.values):
            axes[k].text(v, i, f" {v}", va="center", fontsize=8)
        axes[k].set_xlabel("Hits")
        axes[k].set_title("Top target viruses", fontsize=9)
        k += 1
    if energy is not None and energy.notna().sum() >= 3:
        _multi_hist(axes[k], energy, "Best MFE energy", "kcal/mol", _GREY)
    fig.suptitle("miRNA target prediction overview", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.94), w_pad=2.2)
    return _save(fig, run_dir, "mirna_targets_overview")


# ------------------------------------------------------------------
# 6b. miRNA 关系图：Sankey + 二部网络（PAmiRDB 风格）
# ------------------------------------------------------------------
def _mirna_pairs(run_dir):
    """(miRNA, 靶) 计数对：优先 results.json 的 target_protein（PAmiRDB 口径），
    回退合并 CSV 的 target_virus_name。"""
    rj = _find(run_dir, "results.json", os.path.join("mirna_target", "results.json"))
    if rj:
        try:
            rows = json.load(open(rj, encoding="utf-8"))
            if (isinstance(rows, list) and rows and isinstance(rows[0], dict)
                    and rows[0].get("target_protein")):
                return (Counter((str(r.get("mirna_id")), str(r["target_protein"]))
                                for r in rows), True)
        except Exception:
            pass
    p = _find(run_dir, "mirna_targets_merged.csv")
    if not p:
        return None, False
    df = _read_tab(p)
    mi = _text(df, "mirna_id", min_uniq=1, max_uniq=60)
    tv = _text(df, "target_virus_name", min_uniq=1, max_uniq=60)
    if tv is None:
        tv = _text(df, "target_virus_id", "target_virus", min_uniq=1, max_uniq=60)
    if mi is None or tv is None:
        return None, False
    return Counter(zip(mi.astype(str), tv.astype(str))), False


def plot_mirna_relations(run_dir: str) -> list:
    pairs, protein_level = _mirna_pairs(run_dir)
    if not pairs:
        return []
    sm = Counter({m: sum(w for (a, _), w in pairs.items() if a == m)
                  for (m, _) in pairs})
    st = Counter({t: sum(w for (_, b), w in pairs.items() if b == t)
                  for (_, t) in pairs})
    top_m = [k for k, _ in sm.most_common(8)]
    top_t = [k for k, _ in st.most_common(6)]
    edges = [(m, t, w) for (m, t), w in pairs.items()
             if m in top_m and t in top_t]
    fig, (axs, axn) = plt.subplots(1, 2, figsize=(13.2, 4.6),
                                   gridspec_kw={"wspace": 0.18})
    GREEN, RED = "#2E7D5B", "#E15759"

    def layout(names):
        y, pos, pad = 0.0, {}, 0.4
        for n in sorted(names, key=lambda x: -sm.get(x, st.get(x, 0))):
            tot = sum(w for (a, b), w in pairs.items() if n in (a, b))
            pos[n] = (y, y + tot)
            y += tot + pad
        return pos, y

    lp, lh = layout(top_m)
    rp, rh = layout(top_t)
    H = max(lh, rh)
    for n, (y0, y1) in lp.items():
        axs.fill_betweenx([H - y1, H - y0], 0.06, 0.10, color=GREEN)
        axs.text(0.048, H - (y0 + y1) / 2, _wrap(n, 20), ha="right",
                 va="center", fontsize=8)
    for n, (y0, y1) in rp.items():
        axs.fill_betweenx([H - y1, H - y0], 0.90, 0.94, color=RED)
        axs.text(0.952, H - (y0 + y1) / 2, _wrap(n, 20), ha="left",
                 va="center", fontsize=8)
    from matplotlib.path import Path as _P
    from matplotlib.patches import PathPatch as _PP
    ycur = {n: lp[n][0] for n in top_m}
    rcur = {n: rp[n][0] for n in top_t}
    for m, t, w in sorted(edges, key=lambda e: -e[2]):
        y0 = ycur[m]
        y1 = rcur[t]
        ycur[m] += w
        rcur[t] += w
        xa, xb = 0.10, 0.90
        xm = (xa + xb) / 2
        yl0, yl1 = H - (y0 + w), H - y0
        yr0, yr1 = H - (y1 + w), H - y1
        poly = _P([(xa, yl0), (xm, yl0), (xm, yr0), (xb, yr0),
                   (xb, yr1), (xm, yr1), (xm, yl1), (xa, yl1), (xa, yl0)],
                  [_P.MOVETO, _P.CURVE4, _P.CURVE4, _P.CURVE4,
                   _P.LINETO, _P.CURVE4, _P.CURVE4, _P.CURVE4, _P.CLOSEPOLY])
        axs.add_patch(_PP(poly, facecolor=GREEN, alpha=0.28, edgecolor="none"))
    axs.set_xlim(-0.30, 1.32)
    axs.set_ylim(-0.4, H + 0.4)
    axs.axis("off")
    left_t = "Viral protein" if protein_level else "Target virus"
    axs.set_title(f"miRNA -> {left_t} Sankey", fontsize=10, loc="left")
    yl = np.linspace(0.85, 0.15, len(top_m))
    yr = np.linspace(0.85, 0.15, len(top_t))
    lpos = {n_: (0.22, y) for n_, y in zip(top_m, yl)}
    rpos = {n_: (0.78, y) for n_, y in zip(top_t, yr)}
    wmax = max(w for _, _, w in edges) or 1
    for m, t, w in edges:
        (x0, y0), (x1, y1) = lpos[m], rpos[t]
        axn.annotate("", xy=(x1, y1), xytext=(x0, y0),
                     arrowprops=dict(arrowstyle="-", color=GREEN, alpha=0.30,
                                     lw=0.8 + 2.6 * w / wmax,
                                     connectionstyle="arc3,rad=0.18"))
    for n_, (x, y) in lpos.items():
        axn.scatter([x], [y], s=240, color=GREEN, zorder=3)
        axn.text(x + 0.035, y, _wrap(n_, 18), fontsize=8, va="center")
    for n_, (x, y) in rpos.items():
        axn.scatter([x], [y], s=240, color=RED, zorder=3)
        axn.text(x - 0.035, y, _wrap(n_, 18), fontsize=8, va="center", ha="right")
    axn.set_xlim(0, 1)
    axn.set_ylim(0, 1)
    axn.axis("off")
    axn.set_title(f"miRNA <-> {left_t} network", fontsize=10, loc="left")
    axn.scatter([], [], s=60, color=GREEN, label="miRNA")
    axn.scatter([], [], s=60, color=RED, label=left_t)
    axn.legend(frameon=False, fontsize=8, loc="lower center",
               bbox_to_anchor=(0.5, -0.06), ncol=2)
    return _save(fig, run_dir, "mirna_relations")


# ------------------------------------------------------------------
# 6c. miRNA 双链配对结构图（PAmiRDB duplex：彩色碱基 + 配对连线 + seed 区）
# ------------------------------------------------------------------
_BASE_COLOR = {"A": "#D64545", "U": "#E69F00", "T": "#E69F00",
               "C": "#3B74C4", "G": "#3E9C5C"}


def plot_mirna_duplex(run_dir: str) -> list:
    rj = _find(run_dir, "results.json", os.path.join("mirna_target", "results.json"))
    if not rj:
        return []
    try:
        rows = json.load(open(rj, encoding="utf-8"))
    except Exception:
        return []
    hits = []
    for r in (rows if isinstance(rows, list) else []):
        if not isinstance(r, dict):
            continue
        mseq = str(r.get("mirna_seq") or r.get("mirna_sequence") or "")
        tseq = str(r.get("target_site_seq") or r.get("site_seq")
                   or r.get("target_seq") or "")
        pb = r.get("position_binding")
        if not (len(mseq) >= 12 and len(tseq) >= 12
                and isinstance(pb, dict) and pb):
            # 真实引擎形态：alignment 是三行比对串
            # miRNA  3' SEQ 5' / 配对行(|=WC, .=wobble, 空格=错配) / target 5' SEQ 3'
            aln = r.get("alignment")
            if isinstance(aln, str) and "|" in aln:
                try:
                    import re
                    ls = aln.splitlines()
                    mi = next(l for l in ls if "miRNA" in l)
                    ta = next(l for l in ls if "target" in l)
                    pl = next(l for l in ls
                              if l.strip() and set(l.strip()) <= set("|. "))

                    def _s(line):
                        m_ = re.search(r"[ACGUTNacgutn\-]{10,}", line)
                        return (m_.group(0).upper().replace("T", "U")
                                if m_ else "")
                    mrow, trow = _s(mi), _s(ta)          # 含 '-' 缺口
                    if mrow and trow:
                        mrow, prow, trow = (mrow[::-1], pl.rstrip()[::-1],
                                            trow[::-1])  # 统一成 miRNA 5'->3'
                        pb, ms, ts = {}, [], []
                        i_ = j_ = 0
                        for mb, pc, tb in zip(mrow, prow, trow):
                            if mb == "-" or tb == "-":
                                i_ += mb != "-"
                                j_ += tb != "-"
                                continue
                            i_ += 1
                            j_ += 1
                            ms.append(mb)
                            ts.append(tb)
                            pb[i_] = ("wc" if pc == "|" else
                                      "wobble" if pc == "." else "mismatch")
                        mseq, tseq = "".join(ms), "".join(ts)
                except Exception:
                    mseq = tseq = ""
        if len(mseq) >= 12 and len(tseq) >= 12 and isinstance(pb, dict) and pb:
            hits.append((r, mseq, tseq, pb))
    if not hits:
        return []
    try:
        hits.sort(key=lambda h: -abs(float(h[0].get("best_energy")
                                            or h[0].get("energy") or 0)))
    except Exception:
        pass
    hits = hits[:4]
    n = len(hits)
    W = max(7.0, 0.30 * max(len(h[2]) for h in hits) + 2.5)
    fig, axes = plt.subplots(n, 1, figsize=(W, 2.35 * n), squeeze=False)
    from matplotlib.patches import Rectangle
    for ax, (r, mseq, tseq, pb) in zip(axes[:, 0], hits):
        L = max(len(mseq), len(tseq))
        for x, b in enumerate(mseq, 1):
            ax.scatter([x], [1.15], s=230, color=_BASE_COLOR.get(b, _GREY), zorder=3)
            ax.text(x, 1.15, b, ha="center", va="center", fontsize=7,
                    color="white", zorder=4)
        for x, b in enumerate(tseq, 1):
            ax.scatter([x], [0.05], s=230, color=_BASE_COLOR.get(b, _GREY), zorder=3)
            ax.text(x, 0.05, b, ha="center", va="center", fontsize=7,
                    color="white", zorder=4)
        for k in range(1, L + 1):
            st = str(pb.get(k, pb.get(str(k), ""))).lower()
            if st in ("wc", "watson", "pair", "wc_pair"):
                ax.plot([k, k], [0.28, 0.92], color="#2E8B57", lw=1.8, zorder=2)
            elif st in ("wobble", "gu", "g_u"):
                ax.plot([k, k], [0.28, 0.92], color="#8C8C8C", lw=1.4,
                        ls=(0, (1, 1.6)), zorder=2)
            elif st in ("mismatch", "mm"):
                ax.text(k, 0.60, "X", ha="center", va="center", fontsize=9,
                        color="#D64545", fontweight="bold", zorder=4)
        ax.add_patch(Rectangle((1.5, -0.22), 7, 1.74, fill=True,
                               facecolor="#F5C542", alpha=0.16,
                               edgecolor="#B8860B", ls="--", lw=1.1))
        ax.text(5.0, 1.52, "seed (2-8)", ha="center", fontsize=7, color="#8A6D1A")
        ax.text(0.45, 1.15, "5' miRNA", ha="right", va="center", fontsize=8)
        ax.text(0.45, 0.05, "3' target", ha="right", va="center", fontsize=8)
        name = (f'{r.get("mirna_id", "?")} -> '
                f'{r.get("target_virus_name") or r.get("target_protein") or r.get("target_id", "?")}')
        ax.set_xlim(-0.2, L + 1.2)
        ax.set_ylim(-0.45, 1.72)
        ax.set_title(_wrap(name, 70), fontsize=9, loc="left")
        ax.axis("off")
    fig.suptitle("miRNA : target duplex structure (PAmiRDB style)", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    return _save(fig, run_dir, "mirna_duplex")

# ------------------------------------------------------------------
# 6d. miRNA 明细表（PAmiRDB 风格：可排序 HTML，徽章/seed n/7/靶基因/详情）
# ------------------------------------------------------------------
def plot_mirna_table(run_dir: str) -> list:
    rj = _find(run_dir, "results.json", os.path.join("mirna_target", "results.json"))
    rows = None
    if rj:
        try:
            rows = json.load(open(rj, encoding="utf-8"))
        except Exception:
            rows = None
    if not (isinstance(rows, list) and rows and isinstance(rows[0], dict)):
        p = _find(run_dir, "mirna_targets_merged.csv")
        if not p:
            return []
        df = _read_tab(p)
        rows = df.to_dict("records")
    rows = [r for r in rows if isinstance(r, dict) and r.get("mirna_id")]

    def fnum(r, *keys, default=None):
        for k in keys:
            try:
                v = r.get(k)
                if v not in (None, ""):
                    return float(v)
            except (TypeError, ValueError):
                pass
        return default

    def ftxt(r, *keys):
        for k in keys:
            v = str(r.get(k) or "").strip()
            if v:
                return v
        return ""

    recs = []
    for r in rows:
        lvl = ftxt(r, "consensus_level").upper() or (
            "HIGH" if fnum(r, "num_algorithms", default=0) >= 3 else "MEDIUM")
        recs.append({
            "miRNA": ftxt(r, "mirna_id"),
            "level": lvl,
            "energy": fnum(r, "best_energy", "energy", "miranda_energy"),
            "identity": fnum(r, "best_identity", "miranda_identity", "identity"),
            "mm": fnum(r, "num_mismatches", "mismatches"),
            "seed_wc": fnum(r, "seed_perfect_wc", "seed_wc"),
            "seed_wob": fnum(r, "seed_wobble_count", "seed_wob"),
            "algs": ftxt(r, "consensus_detail", "algorithm", "algorithms"),
            "gene": ftxt(r, "target_protein", "target_gene", "target_virus_name"),
            "mseq": ftxt(r, "mirna_seq", "mirna_sequence"),
            "tseq": ftxt(r, "target_site_seq", "site_seq", "target_seq"),
            "site": f'{ftxt(r, "target_virus_name")} @{int(fnum(r, "target_start", default=0))}'
                    f'-{int(fnum(r, "target_end", default=0))}',
        })
    recs.sort(key=lambda x: (x["energy"] if x["energy"] is not None else 0))

    def pairline(a, b):
        out = []
        for ca, cb in zip(a, b):
            comp = {"A": "U", "U": "A", "T": "A", "C": "G", "G": "C"}
            out.append("|" if comp.get(ca.upper()) == cb.upper() else
                       (":" if {ca.upper(), cb.upper()} in
                        ({"G", "U"}, {"G", "T"}) else " "))
        return "".join(out)

    trs, details = [], []
    for i, r in enumerate(recs):
        badge = {"HIGH": "#DFF0E3", "MEDIUM": "#DDF0F7", "LOW": "#F2F2F2"}.get(
            r["level"], "#F2F2F2")
        seed = (f'{int(r["seed_wc"])}/{int(r["seed_wc"] + (r["seed_wob"] or 0))}'
                if r["seed_wc"] is not None else "-")
        mseq, tseq = r["mseq"], r["tseq"]
        detail = ""
        if mseq and tseq:
            detail = (f"<pre>miRNA 5' {mseq} 3'\n"
                      f"       {pairline(mseq, tseq)}\n"
                      f"target 3' {tseq} 5'</pre>")
        ecell = f'{r["energy"]:.1f}' if r["energy"] is not None else "-"
        icell = f'{r["identity"]:.0f}%' if r["identity"] is not None else "-"
        mmcell = int(r["mm"]) if r["mm"] is not None else 0
        trs.append(
            f'<tr data-e="{r["energy"] or 0}">'
            f'<td>{r["miRNA"]}</td>'
            f'<td><span class="b" style="background:{badge}">{r["level"]}</span></td>'
            f'<td class="n">{ecell}</td>'
            f'<td class="n">{icell}</td>'
            f'<td class="n">{mmcell}</td>'
            f'<td class="n s">{seed}</td>'
            f'<td>{r["algs"]}</td>'
            f'<td class="g">{r["gene"]}</td>'
            f'<td><button class="v" data-i="{i}">View</button></td></tr>')
        details.append(f'<div class="d" id="d{i}" hidden><h4>{r["miRNA"]} &rarr; '
                       f'{r["site"]}</h4>{detail or "无序列数据"}</div>')
    html = """<!doctype html><meta charset="utf-8">
<title>miRNA targets (PAmiRDB style)</title>
<style>body{{font:14px/1.5 "Segoe UI",sans-serif;margin:24px;color:#222}}
table{{border-collapse:collapse;width:100%}}
th,td{{padding:9px 12px;border-bottom:1px solid #eee;text-align:left}}
th{{font-size:11px;letter-spacing:.08em;color:#666;cursor:pointer;user-select:none}}
th:hover{{color:#111}} td.n{{text-align:right;font-variant-numeric:tabular-nums}}
td.s{{color:#1c7c4c;font-weight:600}} td.g{{color:#c96a1b}}
.b{{padding:3px 10px;border-radius:10px;font-size:11px;font-weight:600}}
.v{{border:1px solid #1c7c4c;color:#1c7c4c;background:#fff;border-radius:6px;
padding:3px 12px;cursor:pointer}}
.d{{background:#181c24;color:#dfe6ee;border-radius:10px;padding:16px 20px;margin:8px 0}}
.d pre{{font-family:Consolas,monospace;font-size:13px;white-space:pre-wrap}}
h4{{margin:0 0 10px;color:#fff}}</style>
<h3>miRNA &harr; viral targets (__N__ hits) &mdash; click header to sort</h3>
<table><thead><tr>
<th data-k="miRNA">MIRNA</th><th data-k="level">CONSENSUS</th>
<th data-k="energy">ENERGY</th><th data-k="identity">IDENTITY</th>
<th data-k="mm">MM</th><th data-k="seed">SEED</th>
<th data-k="algs">ALGORITHM</th><th data-k="gene">TARGET GENE</th><th>DETAIL</th>
</tr></thead><tbody>__TRS__</tbody></table>
__DETAILS__
<script>
const tb=document.querySelector("tbody");let dir=1,col="";
document.querySelectorAll("th[data-k]").forEach(th=>th.onclick=()=>{
 const k=th.dataset.k;dir=(col===k)?-dir:1;col=k;
 const num=(k==="energy"||k==="identity"||k==="mm");
 [...tb.rows].sort((a,b)=>{const x=parseFloat(a.dataset.e),y=parseFloat(b.dataset.e);
  if(num)return (x-y)*dir;
  return a.cells[k==="miRNA"?0:1].textContent.localeCompare(
   b.cells[k==="miRNA"?0:1].textContent)*dir;})
 .forEach(r=>tb.appendChild(r));});
document.querySelectorAll(".v").forEach(btn=>btn.onclick=()=>{
 const d=document.getElementById("d"+btn.dataset.i);d.hidden=!d.hidden;});
</script>"""
    html = (html.replace("{{", "{").replace("}}", "}")
                .replace("__N__", str(len(recs)))
                .replace("__TRS__", "".join(trs))
                .replace("__DETAILS__", "".join(details)))
    out = os.path.join(run_dir, "mirna_targets_table.html")
    with open(out, "w", encoding="utf-8") as f:
        f.write(html)
    return [out]

# ------------------------------------------------------------------
# 7. MOT：route × year 迁移矩阵热图
# ------------------------------------------------------------------
def plot_mot(run_dir: str) -> list:
    p = _find(run_dir, "migration_over_time.csv", "*over_time.csv")
    if not p:
        return []
    df = _read_tab(p)
    if df.empty or df.shape[1] < 2:
        return []
    years = df.columns[1:]
    m = df.iloc[:, 1:].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
    if not np.isfinite(m).any() or m.shape[0] == 0:
        return []
    show = np.log1p(np.clip(m, 0, None))
    fig, ax = plt.subplots(figsize=(max(5.0, 0.55 * len(years) + 2.0),
                                    max(3.0, 0.38 * m.shape[0] + 1.4)))
    im = ax.imshow(show, cmap="viridis", aspect="auto")
    ax.set_xticks(range(len(years)))
    ax.set_xticklabels(years, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(m.shape[0]))
    ax.set_yticklabels([_wrap(str(r), 26) for r in df.iloc[:, 0]], fontsize=7)
    if m.shape[0] <= 15 and len(years) <= 20:
        for i in range(m.shape[0]):
            for j in range(len(years)):
                if m[i, j] > 0:
                    rel = (show[i, j] - np.nanmin(show)) / max(
                        np.nanmax(show) - np.nanmin(show), 1e-9)
                    ax.text(j, i, f"{int(m[i, j])}", ha="center", va="center",
                            fontsize=6, color="black" if rel > 0.55 else "white")
                else:
                    # 零值也要标——不标会被当成缺失数据
                    ax.text(j, i, "0", ha="center", va="center", fontsize=6,
                            color="0.55", alpha=0.9)
    fig.colorbar(im, ax=ax, shrink=0.8, label="log1p(migration count)")
    ax.set_title("Migration intensity, route x year (MOT)", fontsize=11)
    return _save(fig, run_dir, "mot_route_year_heatmap")


# ------------------------------------------------------------------
# 8. 候选验证 / contig 分类：判定构成
# ------------------------------------------------------------------
def plot_classification(run_dir: str) -> list:
    p = _find(run_dir, "virus_classification.tsv")
    if not p:
        return []
    df = _read_tab(p)
    if df.empty:
        return []
    verdict = _text(df, "verdict", "判定", "status", "result", min_uniq=1,
                    max_uniq=12)
    tax = _text(df, "kingdom", "phylum", "family", "科", "taxonomy", "lineage",
                "organism", "species", min_uniq=2, max_uniq=40)
    if verdict is None and tax is None:
        return []
    npanel = 1 + (1 if tax is not None else 0)
    fig, axes = plt.subplots(1, npanel, figsize=(4.2 * npanel, 3.6))
    if npanel == 1:
        axes = [axes]
    k = 0
    if verdict is not None:
        vc = verdict.fillna("(blank)").value_counts().sort_values(ascending=True)
        colors = [sp.verdict_color(str(i)) if str(i).lower() in
                  ("pass", "ok", "verified", "virus", "fail", "reject", "host")
                  else _BLUE for i in vc.index]
        axes[k].barh([_wrap(x, 22) for x in vc.index], vc.values,
                     color=colors, edgecolor="white", linewidth=0.4)
        for i, v in enumerate(vc.values):
            axes[k].text(v, i, f" {v}", va="center", fontsize=8)
        axes[k].set_xlabel("Count")
        axes[k].set_title("Verdict composition", fontsize=10)
        k += 1
    if tax is not None:
        tc = tax.value_counts().sort_values(ascending=True).tail(12)
        axes[k].barh([_wrap(x, 26) for x in tc.index], tc.values,
                     color=_GREY, edgecolor="white", linewidth=0.4)
        for i, v in enumerate(tc.values):
            axes[k].text(v, i, f" {v}", va="center", fontsize=8)
        axes[k].set_xlabel("Count")
        axes[k].set_title("Taxonomy composition (top N)", fontsize=10)
    fig.suptitle("Verification / classification overview", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    return _save(fig, run_dir, "classification_overview")


# ------------------------------------------------------------------
# 9. 元数据治理 / 降采样：时间-地理覆盖
# ------------------------------------------------------------------
def plot_metadata_coverage(run_dir: str) -> list:
    p = _find(run_dir, "*metadata*.csv", "*subsampled*.csv", "*govern*.csv")
    if not p:
        return []
    try:
        df = _read_tab(p)
    except Exception:
        return []
    if df.empty:
        return []
    # 先按列名过滤，再尝试日期解析（避免对每个文本列空跑 to_datetime）
    datec = next((c for c in df.columns
                  if any(k in c.lower() for k in ("date", "time", "year",
                                                  "采集", "日期"))
                  and pd.to_datetime(df[c], errors="coerce",
                                     format="mixed").notna().mean() > 0.7), None)
    geoc = _text(df, "location", "country", "region", "province", "geo",
                 "地点", "国家", min_uniq=2, max_uniq=60)
    if datec is None and geoc is None:
        return []
    npanel = 1 + (1 if geoc is not None else 0)
    fig, axes = plt.subplots(1, npanel, figsize=(4.6 * npanel, 3.5))
    if npanel == 1:
        axes = [axes]
    k = 0
    if datec is not None:
        yr = pd.to_datetime(df[datec], errors="coerce").dt.year.dropna().astype(int)
        vc = yr.value_counts().sort_index()
        if len(vc) >= 1:
            axes[k].bar(vc.index.astype(str), vc.values, color=_BLUE,
                        edgecolor="white", linewidth=0.4)
            axes[k].set_xlabel("Year")
            axes[k].set_ylabel("Sequences")
            axes[k].set_title("Sampling time coverage", fontsize=10)
            axes[k].tick_params(axis="x", rotation=45, labelsize=7)
            k += 1
    if geoc is not None and k < npanel:
        gc = geoc.value_counts().sort_values(ascending=True).tail(12)
        axes[k].barh([_wrap(x, 24) for x in gc.index], gc.values,
                     color=_VERM, edgecolor="white", linewidth=0.4)
        for i, v in enumerate(gc.values):
            axes[k].text(v, i, f" {v}", va="center", fontsize=8)
        axes[k].set_xlabel("Sequences")
        axes[k].set_title("Geographic coverage (top N)", fontsize=10)
    fig.suptitle("Metadata coverage overview", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    return _save(fig, run_dir, "metadata_coverage_overview")


# ------------------------------------------------------------------
# 10. 系统地理迁移概览（mugration：迁移矩阵热图 + 迁移路线计数）
# ------------------------------------------------------------------
def plot_phylogeo_migration(run_dir: str) -> list:
    mp = _find(run_dir, "migration_matrix.tsv")
    tp = _find(run_dir, "transitions.tsv")
    if not mp and not tp:
        return []
    fig, axes = plt.subplots(1, 2, figsize=(11.6, 4.4),
                             gridspec_kw={"wspace": 0.34})
    drew = False
    if mp:
        df = _read_tab(mp)
        if not df.empty and {"from", "to", "count"} <= {c.lower() for c in df.columns}:
            df.columns = [c.lower() for c in df.columns]
            df["count"] = pd.to_numeric(df["count"], errors="coerce")
            pv = df.pivot_table(index="from", columns="to",
                                values="count", aggfunc="sum", fill_value=0)
            show = pv.to_numpy(dtype=float)
            im = axes[0].imshow(show, cmap="viridis")
            axes[0].set_xticks(range(len(pv.columns)))
            axes[0].set_xticklabels(pv.columns, rotation=35, ha="right", fontsize=8)
            axes[0].set_yticks(range(len(pv.index)))
            axes[0].set_yticklabels(pv.index, fontsize=8)
            for i in range(show.shape[0]):
                for j in range(show.shape[1]):
                    if show[i, j] > 0:
                        rel = ((show[i, j] - show.min())
                               / max(show.max() - show.min(), 1e-9))
                        axes[0].text(j, i, f"{int(show[i, j])}", ha="center",
                                     va="center", fontsize=8,
                                     color="black" if rel > 0.55 else "white")
            fig.colorbar(im, ax=axes[0], shrink=0.8, label="Migration events")
            axes[0].set_title("Migration matrix (from x to)", fontsize=10)
            drew = True
    if tp:
        df = _read_tab(tp)
        if not df.empty and {"from", "to"} <= {c.lower() for c in df.columns}:
            df.columns = [c.lower() for c in df.columns]
            rc = (df["from"].astype(str) + " → " + df["to"].astype(str))
            vc = rc.value_counts().sort_values(ascending=True).tail(15)
            axes[1].barh([_wrap(x, 24) for x in vc.index], vc.values,
                         color=_VERM, edgecolor="white", linewidth=0.4)
            for i, v in enumerate(vc.values):
                axes[1].text(v, i, f" {v}", va="center", fontsize=8)
            axes[1].set_xlabel("Migration events")
            axes[1].set_title("Top migration routes (phylogeny)", fontsize=10)
            sp.despine()
            drew = True
    if not drew:
        plt.close(fig)
        return []
    fig.suptitle("Phylogeographic migration overview (mugration)", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    return _save(fig, run_dir, "phylogeo_migration_overview")

# ------------------------------------------------------------------
# 入口：内容感知分发
# ------------------------------------------------------------------
_RULES = (
    ("identity_matrix.tsv", plot_identity_heatmap),
    ("events.tsv", plot_rdp_events),
    ("align_qc_report.tsv", plot_align_qc),
    ("primers*.tsv", plot_primers),
    ("crrna_candidates.tsv", plot_crispr),
    ("mirna_targets_merged.csv", (plot_mirna, plot_mirna_relations,
                             plot_mirna_duplex, plot_mirna_table)),
    ("migration_over_time.csv", plot_mot),
    ("migration_matrix.tsv", plot_phylogeo_migration),
    ("virus_classification.tsv", plot_classification),
    ("*metadata*.csv", plot_metadata_coverage),
)


def add_plots(run_dir: str, logger=None) -> list:
    """扫 run 目录，给已知表格产物就地补配图；返回写出的图片路径。"""
    if not run_dir or not os.path.isdir(run_dir):
        return []
    paths: list = []
    for pat, fn in _RULES:
        base = run_dir
        # 真实工具产物常落在 run_dir/mirna_target/ 等子目录——规则未命中时下探一层
        if not _find(base, pat):
            for sub in sorted(os.listdir(base)):
                sd = os.path.join(base, sub)
                if os.path.isdir(sd) and _find(sd, pat):
                    base = sd
                    break
        hit = _find(base, pat)
        if not hit:
            continue
        try:
            out = []
            for f2 in (fn if isinstance(fn, tuple) else (fn,)):
                out += f2(base)
            for fp in out:
                if fp not in paths:
                    paths.append(fp)
        except Exception as e:                                  # noqa: BLE001
            if logger is not None:
                try:
                    logger.log(f"{fn.__name__} 补图失败（跳过）: {e}", "WARN")
                except Exception:
                    pass
    return paths


if __name__ == "__main__":
    d = sys.argv[1] if len(sys.argv) > 1 else "."
    for fp in add_plots(d):
        print(fp)
