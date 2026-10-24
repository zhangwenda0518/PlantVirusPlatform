#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
kv_variant_plot.py — 变异分析可视化段
=====================================
把变异段（kv_variant.py）的产物渲染成出版级图件，供平台页面内联展示。

输入（全部来自变异段与绘图段产物，不重新计算）:
  <out>/vcf/{accession}.vcf              每参考的过滤后变异（含 QUAL/DP/AD）
  <out>/annotated/{accession}.ann.tsv    SnpEff 功能注释（含 Gene_Name/Annotation/Impact）
  <out>/gbk_files/{accession}.gb         GenBank 注释（基因轨道坐标）
  <out>/variant_summary.json             变异段汇总（n_success/n_skipped/results/skipped）

输出:
  <out>/variant_plots/{accession}_variant_map.{png,svg,pdf}      位点分布棒棒糖图
  <out>/variant_plots/{accession}_variant_effect.{png,svg,pdf}   效应谱（Annotation × Impact）
  <out>/variant_plots/{accession}_variant_af_dp.{png,svg,pdf}   等位频率 vs 深度散点
  <out>/variant_plots/{accession}_variant_overview.{png,svg,pdf} 三联组合主图
  <out>/variant_plots/{accession}_snipit.{png,svg,pdf}           跨样本 SNP 分型图
  <out>/variant_plots/{accession}_snipit.snps.csv            snipit 位点明细表
  <out>/variant_plots/variant_plot_manifest.json             产物清单（含 snipit 段）

设计约束（大王 SCI 图偏好）:
  - 白底；Okabe-Ito 色盲友好配色；300 dpi；pdf.fonttype 42
  - 全英文标签；拉丁学名斜体（本模块只用 accession 与基因名，不含拉丁名）
  - 图例不遮挡数据
"""

import importlib.util
import json
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

# adjustText 为可选依赖：缺失时退化为无防重叠排版的纯文本标签
try:
    from adjustText import adjust_text as _adjust_text
except Exception:  # noqa: BLE001
    _adjust_text = None

from .. import sci_palette as sp
sp.apply_sci_style()   # 样式唯一真源（此前本模块手抄 rcParams，必漏键）

# 色板引用（Okabe-Ito 唯一真源；本模块历史名保留为别名）
OKABE_ITO = sp.OKABE
IMPACT_COLOR = sp.IMPACT_COLOR
IMPACT_ORDER = ['HIGH', 'MODERATE', 'LOW', 'MODIFIER']

# 每张变异图最多标注的氨基酸标签数（防止极端样本糊成一片）
MAX_AA_LABELS = 40
# 氨基酸标签默认显示的 AF 门槛（与 246 参考脚本一致的 50%）
AA_LABEL_AF_CUTOFF = 0.50

# 默认出图格式：PNG 供网页/汇报，SVG 供 svguest 二次排版，PDF 供期刊投稿
FIG_FORMATS = ('png', 'svg', 'pdf')


def _is_nonsynonymous(row: dict) -> bool:
    """判断注释行是否为非同义（missense / stop_gained / stop_lost / start_lost 等）。

    仅看 Annotation 字段；同义突变（synonymous_variant）与内含子等跳过。
    """
    eff = (row.get('Annotation') or '').lower()
    if not eff:
        return False
    if 'synonymous' in eff:
        return False
    # 只保留真正改变蛋白序列的效应
    return any(k in eff for k in (
        'missense', 'stop_gained', 'stop_lost', 'start_lost',
        'frameshift', 'inframe_insertion', 'inframe_deletion',
        'disruptive_inframe', 'conservative_inframe', 'initiator_codon',
    ))


def _clean_hgvsp(hgvsp: str) -> str:
    """把 SnpEff 的 HGVSp（如 'p.Gln11Leu'）整理成简洁形式 'Gln11Leu'。

    处理 'p.' 前缀、终止密码子 '*' → 'Ter'、以及多等位位点用逗号分隔的情况
    （只取第一个）。无法解析时返回空串。
    """
    if not hgvsp:
        return ''
    s = str(hgvsp).strip().split(',')[0].strip()
    if s in ('', 'p.'):
        return ''
    if s.startswith('p.'):
        s = s[2:]
    s = s.replace('*', 'Ter')
    # SnpEff 偶尔给出 'p.Gln11Leu' 之外的写法，只保留字母数字
    if not re.fullmatch(r'[A-Za-z0-9]+', s):
        return ''
    return s


# 效应类型显示名（SnpEff 的 snake_case → 阅读友好）
EFFECT_LABEL = {
    'missense_variant': 'Missense',
    'synonymous_variant': 'Synonymous',
    'stop_gained': 'Stop gained',
    'stop_lost': 'Stop lost',
    'start_lost': 'Start lost',
    'start_retained_variant': 'Start retained',
    'stop_retained_variant': 'Stop retained',
    'frameshift_variant': 'Frameshift',
    'inframe_insertion': 'Inframe ins',
    'inframe_deletion': 'Inframe del',
    'disruptive_inframe_insertion': 'Disruptive ins',
    'disruptive_inframe_deletion': 'Disruptive del',
    'conservative_inframe_insertion': 'Conservative ins',
    'conservative_inframe_deletion': 'Conservative del',
    'initiator_codon_variant': 'Initiator codon',
    'terminator_codon_variant': 'Terminator codon',
    'coding_sequence_variant': 'Coding sequence',
}


def safe_name(s):
    """文件名净化（与 kv_plot.safe_name 同口径）"""
    return "".join(c if c.isalnum() or c in "._-" else "_" for c in str(s))


def _decode_field(info: str, key: str):
    """从 VCF INFO 串取指定键的值（第一个匹配）。"""
    for item in info.split(';'):
        if item.startswith(key + '='):
            return item[len(key) + 1:]
        if item == key:
            return True
    return None


def parse_vcf_variants(vcf_path: Path) -> list[dict]:
    """读 VCF，取每条的 QUAL / DP / AD / AF / DP4。

    AF 口径与变异段过滤一致：AD[1] / INFO/DP（无 AD 时退回 AC/AN）。
    """
    rows = []
    with open(vcf_path, encoding='utf-8', errors='replace') as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 8:
                continue
            chrom, pos, _id, ref, alt, qual, filt, info = f[:8]
            try:
                pos_i = int(pos)
            except ValueError:
                continue
            dp = _decode_field(info, 'DP')
            ad = _decode_field(info, 'AD')
            ac = _decode_field(info, 'AC')
            an = _decode_field(info, 'AN')
            dp4 = _decode_field(info, 'DP4')
            try:
                dp_i = int(dp) if dp is not None else None
            except (TypeError, ValueError):
                dp_i = None
            ref_ad = alt_ad = None
            if ad:
                parts = str(ad).split(',')
                if len(parts) >= 2:
                    try:
                        ref_ad, alt_ad = int(parts[0]), int(parts[1])
                    except ValueError:
                        ref_ad = alt_ad = None
            # AF：优先 AD[1]/DP（与过滤口径一致）
            af = None
            if alt_ad is not None and dp_i:
                af = alt_ad / dp_i
            elif ac is not None and an:
                try:
                    af = float(ac) / float(an)
                except (TypeError, ValueError, ZeroDivisionError):
                    af = None
            try:
                qual_f = float(qual)
            except (TypeError, ValueError):
                qual_f = None
            # DP4: ref-fwd, ref-rev, alt-fwd, alt-rev（用于链偏倚）
            alt_fwd = alt_rev = None
            if dp4:
                p4 = str(dp4).split(',')
                if len(p4) == 4:
                    try:
                        alt_fwd, alt_rev = int(p4[2]), int(p4[3])
                    except ValueError:
                        alt_fwd = alt_rev = None
            rows.append({
                'CHROM': chrom, 'POS': pos_i, 'REF': ref, 'ALT': alt,
                'QUAL': qual_f, 'DP': dp_i, 'AF': af, 'FILTER': filt,
                'ALT_AD': alt_ad, 'REF_AD': ref_ad,
                'ALT_FWD': alt_fwd, 'ALT_REV': alt_rev,
                'is_indel': len(ref) != len(alt),
            })
    return rows


def load_ann(ann_path: Path) -> list[dict]:
    """读 SnpEff 注释 TSV（列名与 kv_variant.ANN_COLS 一致）。"""
    rows = []
    with open(ann_path, encoding='utf-8', errors='replace') as fh:
        header = fh.readline().rstrip('\n').split('\t')
        for line in fh:
            if not line.strip():
                continue
            vals = line.rstrip('\n').split('\t')
            vals += [''] * (len(header) - len(vals))
            rows.append(dict(zip(header, vals)))
    return rows


def parse_gb_genes(gb_path: Path, logger=None) -> list[dict]:
    """解析 GenBank 的 CDS/gene 坐标。

    直接复用 kv_plot.parse_gb_features，避免两份实现漂移；
    额外补一个 length 字段（记录全长，变异图坐标轴用）。
    """
    try:
        from Bio import SeqIO
    except ImportError:
        if logger:
            logger.warning('  Biopython 缺失，无法解析 GenBank 注释')
        return []
    from .kv_plot import parse_gb_features
    feats = parse_gb_features(str(gb_path), logger)
    if not feats:
        return []
    try:
        rec_len = len(next(SeqIO.parse(str(gb_path), 'genbank')).seq)
    except Exception as e:  # noqa: BLE001
        if logger:
            logger.warning(f'  {Path(gb_path).name} 记录长度读取失败: {e}')
        rec_len = 0
    for f in feats:
        f.setdefault('length', rec_len)
    # 去重：同坐标优先 CDS（原先这段在 return 之后，永不执行，
    # 且引用了不存在的 `out` —— 重复特征不会被去掉）
    seen = {}
    for f in feats:
        key = (f['start'], f['end'])
        if key not in seen or f['type'] == 'CDS':
            seen[key] = f
    return sorted(seen.values(), key=lambda x: x['start'])


def _style_axis(ax, hide_top_right=True):
    if hide_top_right:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
    ax.tick_params(labelsize=8, length=3)


def _save(fig, out_dir: Path, stem: str, formats=FIG_FORMATS) -> list[str]:
    out_dir.mkdir(parents=True, exist_ok=True)
    made = []
    try:
        for ext in formats:
            p = out_dir / f'{stem}.{ext}'
            try:
                fig.savefig(p, format=ext, dpi=300, bbox_inches='tight',
                            facecolor='white')
            except Exception:
                # SVG 遇到不支持的 artist 时跳过，不阻断 PNG/PDF 主产物
                if ext != 'svg':
                    raise
                continue
            made.append(str(p))
    finally:
        # 多个参考逐一出图时，savefig 抛异常不 close 会累积泄漏 figure
        plt.close(fig)
    return made


# ── 图一：位点分布棒棒糖图 ──────────────────────────────────────
def plot_variant_map(variants, ann_rows, genes, accession, genome_len,
                     out_dir, formats=FIG_FORMATS,
                     nonsyn_only=True, label_af_cutoff=AA_LABEL_AF_CUTOFF,
                     max_aa_labels=MAX_AA_LABELS):
    """沿参考坐标画变异位点；柱高=AF，按 Impact 着色，下方画基因轨道。

    nonsyn_only: 只画非同义突变（与 246 参考脚本一致，画面更聚焦）
    label_af_cutoff: AF 高于该值的非同义位点标注氨基酸变化（HGVS_p）
    max_aa_labels: 单张图标签数量上限
    """
    if not variants:
        return []
    # 注释按 POS 建索引（一个位点可能多条 ANN，取最严重的一条）
    ann_by_pos = {}
    for r in ann_rows:
        try:
            pos = int(r['POS'])
        except (TypeError, ValueError):
            continue
        prev = ann_by_pos.get(pos)
        if prev is None:
            ann_by_pos[pos] = r
        else:
            if (IMPACT_ORDER.index(r['Impact'])
                    if r['Impact'] in IMPACT_ORDER else 9) < (
                    IMPACT_ORDER.index(prev['Impact'])
                    if prev['Impact'] in IMPACT_ORDER else 9):
                ann_by_pos[pos] = r

    # 挑选要画的位点：柱高用 AF（缺 AF 时回退到固定高度）
    def _afeight(v):
        af = v.get('AF')
        if af is None:
            return None
        try:
            af = float(af)
        except (TypeError, ValueError):
            return None
        return max(0.0, min(1.0, af))

    drawn = []
    for v in variants:
        ann = ann_by_pos.get(v['POS'], {})
        if nonsyn_only:
            if not ann or not _is_nonsynonymous(ann):
                continue
        drawn.append((v, ann))
    if not drawn:  # 无非同义时退回全部变异，避免空图
        drawn = [(v, ann_by_pos.get(v['POS'], {})) for v in variants]

    has_genes = bool(genes)
    if has_genes:
        fig, (ax, axg) = plt.subplots(
            2, 1, figsize=(9.5, 4.0), sharex=True,
            gridspec_kw={'height_ratios': [3, 1], 'hspace': 0.12})
    else:
        fig, ax = plt.subplots(figsize=(9.5, 3.0))
        axg = None

    # 变异点：x=位置，y=AF（缺 AF 用 0.85 的固定高度）
    texts = []
    pending = []  # (af, pos, y, label)，用于按 AF 截断
    for v, ann in drawn:
        impact = ann.get('Impact') or 'MODIFIER'
        color = IMPACT_COLOR.get(impact, OKABE_ITO['grey'])
        af = _afeight(v)
        y = af if af is not None else 0.85
        ax.plot([v['POS'], v['POS']], [0.0, y], color=color, lw=0.8,
                alpha=0.8, zorder=2)
        ax.plot([v['POS']], [y], marker='o', ms=4.0, color=color,
                markeredgecolor='white', markeredgewidth=0.4, zorder=3,
                linestyle='none')
        if v['is_indel']:
            ax.plot([v['POS']], [min(1.0, y + 0.05)], marker='v', ms=3.4,
                    color=OKABE_ITO['black'], linestyle='none', zorder=4)
        # 氨基酸标签：仅非同义且 AF 达阈
        aa = _clean_hgvsp(ann.get('HGVS_p', ''))
        if aa and af is not None and af >= label_af_cutoff \
                and _is_nonsynonymous(ann):
            gene = (ann.get('Gene_Name') or '').strip()
            lab = f'{gene}: {aa} ({af * 100:.0f}%)' if gene \
                else f'{aa} ({af * 100:.0f}%)'
            pending.append((af, v['POS'], y, lab))

    # 标签数超上限时优先保留高 AF 位点（防拥挤）
    if max_aa_labels and len(pending) > max_aa_labels:
        pending = sorted(pending, key=lambda x: -x[0])[:max_aa_labels]
    for _af, _pos, y, lab in pending:
        texts.append(ax.text(_pos, y + 0.05, lab, fontsize=6.6,
                             color='#222222', zorder=5,
                             ha='center', va='bottom'))

    ax.set_ylim(0.0, 1.45)
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_ylabel('Allele frequency', fontsize=9)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    _style_axis(ax, hide_top_right=False)
    title = f'{accession} — variant landscape'
    # 轴上方三个文字带自下而上：统计行(1.01) → 图例(1.17) → 标题(y=1.30)。
    # 原版标题用 pad=8 会撞进图例带、统计行撞进标题带
    # （check_figures 体检实抓 108/160/69 px² 三对重叠）。
    ax.set_title(title, fontsize=10.5, fontweight='bold', y=1.30, pad=0)

    # 标签防重叠（adjustText 可用时）
    if texts and _adjust_text is not None:
        try:
            _adjust_text(texts, ax=ax,
                         arrowprops=dict(arrowstyle='-', color='#BBBBBB',
                                         lw=0.6, alpha=0.7, shrinkA=2,
                                         shrinkB=2),
                         expand_points=(1.6, 1.6),
                         expand_text=(1.2, 1.4),
                         force_text=(0.5, 1.2), force_points=(0.3, 0.6),
                         lim=400)
        except Exception:  # noqa: BLE001
            pass

    # 顶部标注统计
    n_high = sum(1 for _, a in drawn if a.get('Impact') == 'HIGH')
    n_mod = sum(1 for _, a in drawn if a.get('Impact') == 'MODERATE')
    n_total = len(drawn)
    stat = f'n = {n_total} nonsynonymous variants   HIGH {n_high}   MODERATE {n_mod}'
    if not nonsyn_only:
        stat = f'n = {n_total} variants   HIGH {n_high}   MODERATE {n_mod}'
    ax.text(0.0, 1.01, stat, transform=ax.transAxes, fontsize=8,
            color='#444444', va='bottom')

    # 图例（右上角外置，带位在统计行之上、标题之下，不遮挡）
    handles = [Line2D([], [], marker='o', color=c, linestyle='none',
                      ms=5, markeredgecolor='white', markeredgewidth=0.4,
                      label=k) for k, c in IMPACT_COLOR.items()]
    handles.append(Line2D([], [], marker='v', color=OKABE_ITO['black'],
                          linestyle='none', ms=4, label='Indel'))
    ax.legend(handles=handles, fontsize=7.5, frameon=False,
              loc='upper right', ncol=len(handles), handletextpad=0.4,
              columnspacing=1.0, bbox_to_anchor=(1.0, 1.17))

    # 基因轨道
    if axg is not None:
        _draw_gene_track(axg, genes, genome_len or max(v['POS'] for v in variants))
        axg.set_xlabel(f'Position along {accession} (nt)', fontsize=9)
    else:
        ax.set_xlabel(f'Position along {accession} (nt)', fontsize=9)
    return _save(fig, out_dir, f'{safe_name(accession)}_variant_map', formats)


def _draw_gene_track(ax, genes, genome_len):
    """基因轨道：横向矩形 + 链方向箭头，Okabe-Ito 循环配色。"""
    # 顶部留 0.75 给最高一条基因的名字（i+0.30 起 6.5pt 文字，
    # 原版 ylim 顶到 len-0.5 会把最上方基因名裁掉一半——目检实抓）
    ax.set_ylim(-0.5, len(genes) + 0.25 if len(genes) > 1 else 1.25)
    ax.set_yticks([])
    ax.spines['left'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    _style_axis(ax, hide_top_right=False)

    palette = [OKABE_ITO[k] for k in
               ('blue', 'orange', 'green', 'purple', 'vermillion', 'sky')]
    for i, g in enumerate(genes):
        c = palette[i % len(palette)]
        start, end = g['start'], g['end']
        ax.add_patch(Rectangle((start, i - 0.22), max(end - start, 1), 0.44,
                               facecolor=c, edgecolor='none', alpha=0.85))
        label = g['name'][:18]
        ax.text(start, i + 0.30, label, fontsize=6.5, color='#333333',
                va='bottom', ha='left', clip_on=True)
    ax.set_xlim(0, genome_len)


# ── 图二：效应谱 ───────────────────────────────────────────────────
def plot_effect_spectrum(ann_rows, accession, out_dir, formats=FIG_FORMATS):
    """Annotation 类型 × Impact 堆叠条形图。"""
    if not ann_rows:
        return []
    counts = {}
    for r in ann_rows:
        eff = r.get('Annotation') or 'unknown'
        imp = r.get('Impact') if r.get('Impact') in IMPACT_ORDER else 'MODIFIER'
        counts.setdefault(eff, {k: 0 for k in IMPACT_ORDER})[imp] += 1
    if not counts:
        return []
    # 排序：按总数降序，同级按最严重程度
    effects = sorted(counts.keys(),
                     key=lambda e: (-sum(counts[e].values()),
                                    min(IMPACT_ORDER.index(k)
                                        for k in IMPACT_ORDER
                                        if counts[e][k] > 0) if any(
                                        counts[e][k] > 0 for k in IMPACT_ORDER) else 9))
    labels = [EFFECT_LABEL.get(e, e.replace('_', ' ')) for e in effects]

    fig, ax = plt.subplots(figsize=(7.6, max(2.6, 0.42 * len(effects) + 1.5)))
    y = list(range(len(effects)))
    left = [0] * len(effects)
    for imp in IMPACT_ORDER:
        vals = [counts[e][imp] for e in effects]
        if not any(vals):
            continue
        ax.barh(y, vals, left=left, color=IMPACT_COLOR[imp], height=0.62,
                label=imp, edgecolor='white', linewidth=0.5)
        for i, (v, l0) in enumerate(zip(vals, left)):
            if v > 0:
                ax.text(l0 + v / 2, y[i], str(v), ha='center', va='center',
                        fontsize=7, color='white', fontweight='bold')
        left = [a + b for a, b in zip(left, vals)]

    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8.5)
    ax.invert_yaxis()
    ax.set_xlabel('Number of variants', fontsize=9)
    ax.set_title(f'{accession} — functional effect spectrum',
                 fontsize=10.5, fontweight='bold', pad=8)
    _style_axis(ax)
    ax.legend(fontsize=7.5, frameon=False, loc='lower right', title=None)
    ax.grid(axis='x', color='#DDDDDD', lw=0.5, alpha=0.7)
    ax.set_axisbelow(True)
    return _save(fig, out_dir, f'{safe_name(accession)}_variant_effect', formats)


# ── 图三：等位频率 vs 深度 ─────────────────────────────────────────
def plot_af_dp(variants, ann_rows, accession, out_dir, formats=FIG_FORMATS):
    """AF vs DP 散点；颜色=Impact，大小=QUAL，虚线标 min-freq 阈值。"""
    pts = [v for v in variants if v['AF'] is not None and v['DP']]
    if not pts:
        return []
    ann_by_pos = {}
    for r in ann_rows:
        try:
            pos = int(r['POS'])
        except (TypeError, ValueError):
            continue
        prev = ann_by_pos.get(pos)
        if prev is None or (IMPACT_ORDER.index(r['Impact'])
                            if r['Impact'] in IMPACT_ORDER else 9) < (
                            IMPACT_ORDER.index(prev['Impact'])
                            if prev['Impact'] in IMPACT_ORDER else 9):
            ann_by_pos[pos] = r

    fig, ax = plt.subplots(figsize=(6.2, 4.4))
    for v in pts:
        imp = ann_by_pos.get(v['POS'], {}).get('Impact') or 'MODIFIER'
        c = IMPACT_COLOR.get(imp, OKABE_ITO['grey'])
        size = 18 + 42 * min(1.0, math.log10(max(v['QUAL'] or 1, 1.0)) / 3.0)
        ax.scatter(v['DP'], v['AF'], s=size, color=c, alpha=0.75,
                   edgecolors='white', linewidths=0.5, zorder=3)
        if v['is_indel']:
            ax.scatter(v['DP'], v['AF'], s=size + 26, facecolors='none',
                       edgecolors=OKABE_ITO['black'], linewidths=0.6,
                       zorder=4)

    ax.axhline(0.05, color='#999999', ls='--', lw=0.8, zorder=1)
    ax.text(ax.get_xlim()[1], 0.052, ' min-freq 0.05', fontsize=7,
            color='#666666', va='bottom', ha='right')
    ax.set_xlabel('Read depth (DP)', fontsize=9)
    ax.set_ylabel('Allele frequency (AD[1]/DP)', fontsize=9)
    ax.set_title(f'{accession} — allele frequency vs depth',
                 fontsize=10.5, fontweight='bold', pad=8)
    ax.set_ylim(-0.03, 1.03)
    _style_axis(ax)
    ax.grid(color='#EEEEEE', lw=0.5, zorder=0)
    ax.set_axisbelow(True)

    handles = [Line2D([], [], marker='o', color=c, linestyle='none', ms=5.5,
                      markeredgecolor='white', label=k)
               for k, c in IMPACT_COLOR.items()]
    handles.append(Line2D([], [], marker='o', markerfacecolor='none',
                          markeredgecolor=OKABE_ITO['black'], linestyle='none',
                          ms=7, label='Indel'))
    ax.legend(handles=handles, fontsize=7.5, frameon=False,
              loc='lower right', ncol=2, handletextpad=0.4, columnspacing=1.0)
    return _save(fig, out_dir, f'{safe_name(accession)}_variant_af_dp', formats)


# ── 组合主图（三面板）─────────────────────────────────────────────
def plot_overview(variants, ann_rows, genes, accession, genome_len,
                  out_dir, formats=FIG_FORMATS):
    """三面板组合主图：上=位点分布，左下=效应谱，右下=AF/DP。"""
    if not variants and not ann_rows:
        return []
    fig = plt.figure(figsize=(11.5, 7.2))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.25],
                          width_ratios=[1.25, 1.0],
                          hspace=0.42, wspace=0.28)

    # 上：位点分布（跨两列）
    axt = fig.add_subplot(gs[0, :])
    _panel_variant_map(axt, variants, ann_rows, genes, accession, genome_len)

    # 左下：效应谱
    axl = fig.add_subplot(gs[1, 0])
    _panel_effect_bars(axl, ann_rows, accession)

    # 右下：AF/DP
    axr = fig.add_subplot(gs[1, 1])
    _panel_af_dp(axr, variants, ann_rows, accession)

    fig.suptitle(f'Variant analysis overview — {accession}',
                 fontsize=12, fontweight='bold', y=0.985)
    return _save(fig, out_dir, f'{safe_name(accession)}_variant_overview', formats)


def _panel_variant_map(ax, variants, ann_rows, genes, accession, genome_len):
    ann_by_pos = {}
    for r in ann_rows:
        try:
            pos = int(r['POS'])
        except (TypeError, ValueError):
            continue
        prev = ann_by_pos.get(pos)
        if prev is None or (IMPACT_ORDER.index(r['Impact'])
                            if r['Impact'] in IMPACT_ORDER else 9) < (
                            IMPACT_ORDER.index(prev['Impact'])
                            if prev['Impact'] in IMPACT_ORDER else 9):
            ann_by_pos[pos] = r
    for v in variants:
        imp = ann_by_pos.get(v['POS'], {}).get('Impact') or 'MODIFIER'
        c = IMPACT_COLOR.get(imp, OKABE_ITO['grey'])
        y = 1.0 + (min(0.9, math.log10(max(v['QUAL'] or 1, 1.0)) / 3.0)
                   if v['QUAL'] else 0.0)
        ax.plot([v['POS'], v['POS']], [0.45, y], color=c, lw=0.6,
                alpha=0.7, zorder=2)
        ax.plot([v['POS']], [y], marker='o', ms=3.6, color=c,
                markeredgecolor='white', markeredgewidth=0.35,
                linestyle='none', zorder=3)
        if v['is_indel']:
            ax.plot([v['POS']], [y + 0.13], marker='v', ms=3.0,
                    color=OKABE_ITO['black'], linestyle='none', zorder=4)
    ax.set_ylim(0.3, 2.15)
    ax.set_yticks([])
    ax.spines['left'].set_visible(False)
    _style_axis(ax)
    ax.set_xlabel(f'Position along {accession} (nt)', fontsize=8.5)
    ax.set_title('A  Variant landscape', fontsize=9.5, fontweight='bold',
                 loc='left', pad=6)
    handles = [Line2D([], [], marker='o', color=c, linestyle='none', ms=4.5,
                      markeredgecolor='white', label=k)
               for k, c in IMPACT_COLOR.items()]
    handles.append(Line2D([], [], marker='v', color=OKABE_ITO['black'],
                          linestyle='none', ms=3.6, label='Indel'))
    ax.legend(handles=handles, fontsize=6.8, frameon=False,
              loc='upper right', ncol=len(handles), handletextpad=0.35,
              columnspacing=0.8, bbox_to_anchor=(1.0, 1.20))


def _panel_effect_bars(ax, ann_rows, accession):
    counts = {}
    for r in ann_rows:
        eff = r.get('Annotation') or 'unknown'
        imp = r.get('Impact') if r.get('Impact') in IMPACT_ORDER else 'MODIFIER'
        counts.setdefault(eff, {k: 0 for k in IMPACT_ORDER})[imp] += 1
    if not counts:
        ax.axis('off')
        return
    effects = sorted(counts.keys(), key=lambda e: -sum(counts[e].values()))
    labels = [EFFECT_LABEL.get(e, e.replace('_', ' ')) for e in effects]
    y = list(range(len(effects)))
    left = [0] * len(effects)
    for imp in IMPACT_ORDER:
        vals = [counts[e][imp] for e in effects]
        if not any(vals):
            continue
        ax.barh(y, vals, left=left, color=IMPACT_COLOR[imp], height=0.6,
                label=imp, edgecolor='white', linewidth=0.4)
        left = [a + b for a, b in zip(left, vals)]
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=7.5)
    ax.invert_yaxis()
    ax.set_xlabel('Number of variants', fontsize=8.5)
    ax.set_title('B  Effect spectrum', fontsize=9.5, fontweight='bold',
                 loc='left', pad=6)
    _style_axis(ax)
    ax.grid(axis='x', color='#DDDDDD', lw=0.5, alpha=0.7)
    ax.set_axisbelow(True)


def _panel_af_dp(ax, variants, ann_rows, accession):
    pts = [v for v in variants if v['AF'] is not None and v['DP']]
    if not pts:
        ax.axis('off')
        return
    ann_by_pos = {}
    for r in ann_rows:
        try:
            pos = int(r['POS'])
        except (TypeError, ValueError):
            continue
        prev = ann_by_pos.get(pos)
        if prev is None or (IMPACT_ORDER.index(r['Impact'])
                            if r['Impact'] in IMPACT_ORDER else 9) < (
                            IMPACT_ORDER.index(prev['Impact'])
                            if prev['Impact'] in IMPACT_ORDER else 9):
            ann_by_pos[pos] = r
    for v in pts:
        imp = ann_by_pos.get(v['POS'], {}).get('Impact') or 'MODIFIER'
        c = IMPACT_COLOR.get(imp, OKABE_ITO['grey'])
        ax.scatter(v['DP'], v['AF'], s=16, color=c, alpha=0.75,
                   edgecolors='white', linewidths=0.4, zorder=3)
    ax.axhline(0.05, color='#999999', ls='--', lw=0.7, zorder=1)
    ax.set_xlabel('Read depth (DP)', fontsize=8.5)
    ax.set_ylabel('Allele frequency', fontsize=8.5)
    ax.set_ylim(-0.03, 1.03)
    ax.set_title('C  Allele frequency vs depth', fontsize=9.5,
                 fontweight='bold', loc='left', pad=6)
    _style_axis(ax)
    ax.grid(color='#EEEEEE', lw=0.5, zorder=0)
    ax.set_axisbelow(True)


# ── 跨样本 SNP 分型图（snipit 封装）───────────────────────────────
# snipit（GPL-3.0，pip install snipit，Bioinformatics 2024 btae510）把比对
# 序列相对参考的 SNP 画成「一行一样本」的分型图，回答现有三张图不覆盖的
# 问题：同一种病毒的 N 个样本，谁携带哪些共享/私有 SNP。
# 输入不重新比对——多样本 VCF 本身就在参考坐标系里（BcftoolsCaller 把全部
# BAM 一起 calling），按基因型把 REF/ALT 套回参考序列即得伪比对，snipit
# 只负责画。snipit 未安装时静默降级（与 adjustText 同一容错口径）。
SNIPIT_MIN_SAMPLES = 2
# vcf/ 里的 caller 中间产物（all.variants.raw/soft/filtered、demo 入口
# VCF 等）不构成参考，靠参考序列查找天然排除之外再按名字显式跳过，
# 避免 manifest 的 skipped 里堆噪音
_SNIPIT_VCF_SKIP = re.compile(r'^(all|demo)\.variants')


def _find_snipit_runner() -> list[str] | None:
    """snipit 可用性探测。以 `python -m snipit.command` 调用（console
    script 装在用户 Scripts 目录、不一定在 PATH 上）；未安装返回 None。"""
    if importlib.util.find_spec('snipit') is None:
        return None
    return [sys.executable, '-m', 'snipit.command']


def _read_vcf_genotypes(vcf_path: Path) -> tuple[list[str], list[dict]]:
    """读多样本 VCF → (样本名列表, 行列表)。

    行: {'CHROM', 'POS', 'REF', 'ALTS', 'GTS'}；GTS 与样本名对齐
    （'0/1'、'1|2'、'./.'）。sites-only VCF 返回空样本列表。
    """
    samples: list[str] = []
    rows: list[dict] = []
    with open(vcf_path, encoding='utf-8', errors='replace') as fh:
        for line in fh:
            if line.startswith('##'):
                continue
            if line.startswith('#CHROM'):
                samples = line.rstrip('\n').split('\t')[9:]
                continue
            if not line.strip():
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 8:
                continue
            gts: list = [None] * len(samples)
            if len(f) > 9 and len(samples) == len(f) - 9:
                fmt = f[8].split(':')
                gi = fmt.index('GT') if 'GT' in fmt else -1
                if gi >= 0:
                    for si, sc in enumerate(f[9:]):
                        parts = sc.split(':')
                        gts[si] = parts[gi] if gi < len(parts) else None
            try:
                pos = int(f[1])
            except ValueError:
                continue
            rows.append({'CHROM': f[0], 'POS': pos, 'REF': f[3].upper(),
                         'ALTS': [a.upper() for a in f[4].split(',')],
                         'GTS': gts})
    return samples, rows


def _load_ref_sequence(acc: str, out_dir: Path, gbk_dirs: list[Path]):
    """取参考序列 (seq_id, seq)：virus-fasta/ 优先（header 与 accession
    精确/裸版本号双向匹配，同 kv_variant.export_ref_fastas 的宽松口径），
    GenBank 兜底（SnpEff 建库必有 .gb）。都取不到返回 None。"""
    try:
        from Bio import SeqIO
    except ImportError:
        return None
    bare = str(acc).split('.')[0]
    vfa = out_dir / 'virus-fasta'
    if vfa.is_dir():
        for p in sorted(vfa.glob('ref_*/*.ref.fasta')):
            try:
                rec = next(SeqIO.parse(str(p), 'fasta'))
            except (StopIteration, ValueError, OSError):
                continue
            if rec.id == str(acc) or rec.id.split('.')[0] == bare:
                return rec.id, str(rec.seq).upper()
    for d in gbk_dirs:
        if not d.is_dir():
            continue
        for cand in (d / f'{acc}.gb', d / f'{acc}.gbk',
                     d / f'{bare}.gb', d / f'{bare}.gbk'):
            if not cand.is_file():
                continue
            try:
                rec = next(SeqIO.parse(str(cand), 'genbank'))
                return (rec.id or str(acc)), str(rec.seq).upper()
            except (StopIteration, ValueError, OSError):
                break
    return None


def build_snipit_alignment(vcf_path: Path, ref_id: str, ref_seq: str,
                           out_aln: Path, chrom: str | None = None) -> dict:
    """从多样本 VCF 重建参考坐标伪比对（参考 + 各样本基因型套回）。

    SNV/MNP 直接套换；纯缺失型（len(ALT)<len(REF)）套换后以 '-' 补齐
    REF 跨度，snipit 读作 deletion；行内任一等位含插入（len(ALT)>
    len(REF)）时整行跳过——定长伪比对表达不了插入（会错位），这些位点
    保留在 VCF 与其余图件里。无基因型（./.）的样本该位点记参考碱基：
    caller 无覆盖时等位本就未知，与共识序列口径一致。
    chrom: 只套该染色体的行（缺省 None 不筛，调用方传 ref_id 更稳）。
    返回统计 dict；样本数/可套位点数为 0 时调用方跳过。
    """
    samples, rows = _read_vcf_genotypes(vcf_path)
    stats = {'n_samples': len(samples), 'n_rows': len(rows), 'n_sites': 0,
             'n_rows_skipped_insertion': 0, 'samples': samples}
    if not samples:
        return stats
    ref_bare = str(ref_id).split('.')[0]
    seqs = {s: list(ref_seq) for s in samples}
    for r in rows:
        if chrom and not (r['CHROM'] == str(ref_id)
                          or r['CHROM'].split('.')[0] == ref_bare):
            continue
        span = len(r['REF'])
        end = r['POS'] - 1 + span
        if end > len(ref_seq):
            continue
        if any(len(a) > span for a in r['ALTS']):
            stats['n_rows_skipped_insertion'] += 1
            continue
        applied = False
        for si, s in enumerate(samples):
            gt = r['GTS'][si] if si < len(r['GTS']) else None
            if not gt:
                continue
            for tok in str(gt).replace('|', '/').split('/'):
                if tok in ('.', '0', ''):
                    continue
                try:
                    k = int(tok)
                except ValueError:
                    continue
                if not 1 <= k <= len(r['ALTS']):
                    continue
                alt = r['ALTS'][k - 1]
                # 缺失型：套换 ALT 后用 '-' 补齐跨度（snipit 读作 deletion）
                seqs[s][r['POS'] - 1:end] = list(alt) + ['-'] * (span - len(alt))
                applied = True
                break   # 一个样本一行只套第一个非参考等位
        if applied:
            stats['n_sites'] += 1

    out_aln.parent.mkdir(parents=True, exist_ok=True)

    # bcftools 多样本 calling 无 @RG/SM 时用 BAM 完整路径当样本列名，
    # 直接进图会是一长串路径——清洗成裸文件名（去目录、去 .bam 后缀），
    # 撞名时加序号兜底
    display: dict[str, str] = {}
    used: set[str] = set()
    for s in samples:
        name = Path(str(s)).name
        for suf in ('.sorted.bam', '.bam'):
            if name.endswith(suf):
                name = name[:-len(suf)]
                break
        base, k = name, 2
        while name in used:
            name = f'{base}_{k}'
            k += 1
        used.add(name)
        display[s] = name

    def _wrap(s):
        return '\n'.join(s[i:i + 80] for i in range(0, len(s), 80))

    with open(out_aln, 'w', encoding='utf-8') as fh:
        fh.write(f'>{ref_id}\n{_wrap(ref_seq)}\n')
        for s in samples:
            hdr = ''.join(c if not c.isspace() else '_' for c in display[s])
            fh.write(f'>{hdr}\n{_wrap("".join(seqs[s]))}\n')
    stats['samples'] = [display[s] for s in samples]
    return stats


def _snipit_one(acc, vcf_path, out_dir, gbk_dirs, plot_dir, runner,
                min_samples, formats, res, _log):
    """单参考的 snipit 出图（异常向上抛，由 run_snipit_plots 统一记 skipped）。"""
    samples, _ = _read_vcf_genotypes(vcf_path)
    if len(samples) < min_samples:
        reason = ('VCF 无样本列（sites-only）' if not samples
                  else f'样本数 {len(samples)} < {min_samples}（单样本 '
                       f'variant_map 已覆盖）')
        _log('info', f'  snipit {acc}: 跳过（{reason}）')
        res['skipped'].append({'accession': acc, 'reason': reason})
        return
    ref = _load_ref_sequence(acc, out_dir, gbk_dirs)
    if ref is None:
        _log('warning', f'  snipit {acc}: 参考序列不可得'
                        f'（virus-fasta/ 与 GenBank 均缺），跳过')
        res['skipped'].append({'accession': acc, 'reason': '参考序列不可得'})
        return
    ref_id, ref_seq = ref
    stem = f'{safe_name(acc)}_snipit'
    made = []
    # snipit 的 snps.csv 是固定文件名（snp_functions.write_out_snps 硬编码
    # output_dir/snps.csv），必须落暂存目录再改名归位，避免条目间互相覆盖
    with tempfile.TemporaryDirectory(prefix='snipit_') as td:
        aln = Path(td) / 'aln.fasta'
        st = build_snipit_alignment(vcf_path, ref_id, ref_seq, aln, chrom=ref_id)
        if st['n_sites'] <= 0:
            _log('info', f'  snipit {acc}: 无可套回变异，跳过')
            res['skipped'].append({'accession': acc, 'reason': '无可套回变异'})
            return
        for i, ext in enumerate(formats):
            cmd = runner + [str(aln), '-r', ref_id, '-o', stem,
                            '-d', td, '-f', ext, '-s',
                            '--sort-by-mutation-number', '--show-indels',
                            '--solid-background']
            p = subprocess.run(cmd, capture_output=True, text=True,
                               encoding='utf-8', errors='replace')
            src = Path(td) / f'{stem}.{ext}'
            if p.returncode != 0 or not src.is_file():
                tail = ((p.stderr or p.stdout or '').strip().splitlines()
                        or [''])[-3:]
                raise RuntimeError(
                    f'snipit -f {ext} 失败(rc={p.returncode}): '
                    + ' | '.join(tail))
            plot_dir.mkdir(parents=True, exist_ok=True)
            dst = plot_dir / f'{stem}.{ext}'
            shutil.move(str(src), str(dst))
            made.append(dst)
            if i == 0:
                snps_src = Path(td) / 'snps.csv'
                if snps_src.is_file():
                    dst_csv = plot_dir / f'{stem}.snps.csv'
                    shutil.move(str(snps_src), str(dst_csv))
                    made.append(dst_csv)
    for p in made:
        res['plots'].append(str(p))
    res['accessions'][acc] = {
        'n_samples': st['n_samples'], 'n_snps': st['n_sites'],
        'n_rows_skipped_insertion': st['n_rows_skipped_insertion'],
    }
    _log('info', f'  snipit {acc}: {st["n_samples"]} 样本 × '
                 f'{st["n_sites"]} 位点 → {len(made)} 个产物'
                 + (f'（{st["n_rows_skipped_insertion"]} 行含插入未套回）'
                    if st['n_rows_skipped_insertion'] else ''))


def run_snipit_plots(out_dir, logger=None, accession_filter=None,
                     gbk_dirs=None, min_samples=SNIPIT_MIN_SAMPLES,
                     formats=FIG_FORMATS) -> dict:
    """跨样本 SNP 分型图（snipit）编排入口。

    为 <out>/vcf/*.vcf 中样本列 ≥ min_samples 的参考出图；snipit 未安装、
    参考序列不可得、sites-only VCF、样本数不足均按条跳过，不抛异常。
    产物: variant_plots/{gid}_snipit.{png,pdf} + {gid}_snipit.snps.csv
    返回 {'available', 'n_accessions', 'plots', 'accessions', 'skipped'}
    """
    out_dir = Path(out_dir)
    vcf_dir = out_dir / 'vcf'
    plot_dir = out_dir / 'variant_plots'
    res = {'available': True, 'n_accessions': 0, 'plots': [],
           'accessions': {}, 'skipped': []}

    def _log(level, msg):
        if logger is not None:
            getattr(logger, level, logger.info)(msg)

    runner = _find_snipit_runner()
    if runner is None:
        _log('info', 'snipit 未安装，跳过跨样本 SNP 分型图'
                     '（pip install snipit 启用）')
        res['available'] = False
        return res
    if not vcf_dir.is_dir():
        res['available'] = False
        return res

    gbk_dirs = [Path(d) for d in (gbk_dirs or [])]
    accs = sorted(p.stem for p in vcf_dir.glob('*.vcf')
                  if not _SNIPIT_VCF_SKIP.match(p.stem))
    if accession_filter:
        keep = set(accession_filter)
        accs = [a for a in accs if a in keep]

    for acc in accs:
        try:
            _snipit_one(acc, vcf_dir / f'{acc}.vcf', out_dir, gbk_dirs,
                        plot_dir, runner, min_samples, formats, res, _log)
        except Exception as e:  # noqa: BLE001
            _log('warning', f'  snipit {acc}: 失败 {e}')
            res['skipped'].append({'accession': acc, 'reason': str(e)})
    res['n_accessions'] = len(res['accessions'])
    return res


# ── 编排入口 ───────────────────────────────────────────────────────
def run_variant_plots(out_dir, logger=None, formats=FIG_FORMATS,
                      accession_filter=None, nonsyn_only=True,
                      label_af_cutoff=AA_LABEL_AF_CUTOFF,
                      max_aa_labels=MAX_AA_LABELS, gbk_dir=None,
                      run_snipit=True) -> dict:
    """为 out_dir 下所有有变异+注释的参考出图。

    nonsyn_only: 变异图只画非同义突变（默认）
    label_af_cutoff: 标注氨基酸变化的 AF 下限（默认 0.50）
    max_aa_labels: 单张图标签数量上限（默认 40）
    gbk_dir: GenBank 目录；缺省时依次找 <out>/gbk_files、<out>/virus-annotations
    run_snipit: 追加跨样本 SNP 分型图（snipit，未安装时静默跳过）

    返回 {'n_accessions': n, 'plots': [...相对路径...], 'skipped': [...],
          'snipit': run_snipit_plots 的返回（run_snipit=False 时为空 dict）}
    """
    out_dir = Path(out_dir)
    vcf_dir = out_dir / 'vcf'
    ann_dir = out_dir / 'annotated'
    # 变异段把归一化注释导出到 <out>/virus-annotations/，而本段历史上只找
    # <out>/gbk_files/ —— 两处都找，否则基因轨道会静默为空。
    gbk_dirs: list[Path] = []
    if gbk_dir:
        gbk_dirs.append(Path(gbk_dir))
    for d in (out_dir / 'gbk_files', out_dir / 'virus-annotations'):
        if d not in gbk_dirs:
            gbk_dirs.append(d)
    plot_dir = out_dir / 'variant_plots'

    def _log(level, msg):
        if logger is None:
            return
        getattr(logger, level, logger.info)(msg)

    if not vcf_dir.is_dir():
        _log('warning', f'无 vcf 目录，跳过变异绘图: {vcf_dir}')
        return {'n_accessions': 0, 'plots': [], 'skipped': []}

    # 以注释目录为准（有注释才有基因/效应信息）
    accessions = []
    if ann_dir.is_dir():
        for fn in sorted(os.listdir(ann_dir)):
            if fn.endswith('.ann.tsv'):
                accessions.append(fn[:-len('.ann.tsv')])
    if accession_filter:
        keep = set(accession_filter)
        accessions = [a for a in accessions if a in keep]
    if not accessions:
        _log('warning', '无 .ann.tsv 注释文件，跳过变异绘图')
        return {'n_accessions': 0, 'plots': [], 'skipped': []}

    all_plots = []
    skipped = []
    for acc in accessions:
        vcf_path = vcf_dir / f'{acc}.vcf'
        ann_path = ann_dir / f'{acc}.ann.tsv'
        if not vcf_path.is_file():
            _log('warning', f'  {acc}: 缺 VCF，跳过')
            skipped.append({'accession': acc, 'reason': '缺 VCF'})
            continue
        try:
            variants = parse_vcf_variants(vcf_path)
            ann_rows = load_ann(ann_path)
        except Exception as e:
            _log('warning', f'  {acc}: 读取失败 {e}')
            skipped.append({'accession': acc, 'reason': f'读取失败: {e}'})
            continue
        if not variants and not ann_rows:
            _log('info', f'  {acc}: 无变异记录，跳过')
            skipped.append({'accession': acc, 'reason': '无变异记录'})
            continue

        gbk = None
        for _d in gbk_dirs:
            if not _d.is_dir():
                continue
            cand = _d / f'{acc}.gb'
            if cand.is_file():
                gbk = cand
                break
            # 注释文件名可能省略版本号（OR489165 vs OR489165.1），回退搜索
            found = sorted(_d.glob(f'{safe_name(acc)}*.gb'))
            if not found:
                found = sorted(_d.glob(f'{str(acc).split(".")[0]}*.gb'))
            if found:
                gbk = found[0]
                _log('debug', f'  {acc}: GenBank 匹配到 {gbk.name}（{_d.name}/）')
                break
        if gbk is None:
            gbk = gbk_dirs[0] / f'{acc}.gb'      # 仅用于日志/路径展示
        class _LogAdapter:
            """把 run_variant_plots 的 _log 包成 logger 对象，给复用函数用"""
            def __getattr__(self, lv):
                return lambda msg, *a: _log(lv if lv in (
                    'info', 'warning', 'error', 'debug') else 'info',
                    msg % a if a else msg)

        genes = parse_gb_genes(gbk, _LogAdapter()) if gbk.is_file() else []
        if not genes:
            _log('debug', f'  {acc}: 无基因轨道（{gbk.name if gbk.is_file() else "缺 GenBank"}）')
        genome_len = 0
        if genes:
            genome_len = max(g.get('length', 0) for g in genes)
        if not genome_len and variants:
            genome_len = max(v['POS'] for v in variants)

        made = []
        try:
            made += plot_variant_map(variants, ann_rows, genes, acc,
                                     genome_len, plot_dir, formats,
                                     nonsyn_only=nonsyn_only,
                                     label_af_cutoff=label_af_cutoff,
                                     max_aa_labels=max_aa_labels)
            made += plot_effect_spectrum(ann_rows, acc, plot_dir, formats)
            made += plot_vaf_dp_safe(variants, ann_rows, acc, plot_dir, formats)
            made += plot_overview(variants, ann_rows, genes, acc,
                                  genome_len, plot_dir, formats)
        except Exception as e:
            _log('error', f'  {acc}: 绘图失败 {e}')
            skipped.append({'accession': acc, 'reason': f'绘图失败: {e}'})
            continue

        for p in made:
            all_plots.append(os.path.relpath(p, out_dir).replace(os.sep, '/'))
        _log('info', f'  {acc}: {len(made)} 张图（变异 {len(variants)} 条，'
                     f'注释 {len(ann_rows)} 条）')

    manifest = {
        'n_accessions': len(accessions),
        'n_plots': len(all_plots),
        'plots': sorted(all_plots),
        'skipped': skipped,
    }

    # 跨样本 SNP 分型图（snipit）：样本列 ≥2 的参考逐条出图，独立于上面的
    # 注释驱动循环（snipit 不需要 SnpEff 注释，只要 VCF 带样本列 + 参考序列）
    if run_snipit:
        try:
            snipit_res = run_snipit_plots(out_dir, logger,
                                          accession_filter=accession_filter,
                                          gbk_dirs=gbk_dirs)
            if snipit_res.get('n_accessions'):
                manifest['n_snipit_plots'] = len(snipit_res.get('plots', []))
                _log('info', f'跨样本分型图（snipit）：'
                             f'{snipit_res["n_accessions"]} 条参考 → '
                             f'{len(snipit_res.get("plots", []))} 个产物')
        except Exception as e:  # noqa: BLE001
            _log('warning', f'跨样本分型图（snipit）失败（不影响其余图件）: {e}')
            snipit_res = {'available': True, 'error': str(e)}
        manifest['snipit'] = snipit_res

    plot_dir.mkdir(parents=True, exist_ok=True)
    (plot_dir / 'variant_plot_manifest.json').write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2), encoding='utf-8')
    _log('info', f'变异绘图完成：{len(accessions)} 条参考，{len(all_plots)} 张图')
    return manifest


def plot_vaf_dp_safe(variants, ann_rows, accession, out_dir, formats):
    """plot_af_dp 的容错包装（内部函数名历史原因保留）。"""
    return plot_af_dp(variants, ann_rows, accession, out_dir, formats)


if __name__ == '__main__':
    import argparse
    import logging
    ap = argparse.ArgumentParser(description='变异分析可视化')
    ap.add_argument('--out', required=True, help='variant 段产物目录')
    ap.add_argument('--acc', nargs='*', default=None, help='只画指定 accession')
    ap.add_argument('--verbose', action='store_true')
    a = ap.parse_args()
    logging.basicConfig(
        level=logging.DEBUG if a.verbose else logging.INFO,
        format='%(asctime)s %(levelname)s %(message)s', datefmt='%H:%M:%S')
    res = run_variant_plots(a.out, logging.getLogger(),
                            accession_filter=a.acc)
    print(json.dumps(res, ensure_ascii=False, indent=2))
