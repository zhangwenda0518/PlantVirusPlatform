# -*- coding: utf-8 -*-
"""样本管线流程漏斗三件套：TSV 汇总 + 梯形漏斗图 + 留存率条形图。

借鉴 MMPV-RNA 的 stage_funnel_summary（蓝本 docs/MMPV-RNA结果解析与绘图
学习笔记.md §2.1）：把「输入 reads → 质控 → 去宿主 → 组装 → 病毒 contigs
→ 验证 → 共识 → ORF」的层层衰减讲成一张图 + 一张表，供 ⑩ 报告与
结果中心复用。

设计约定（退化守卫）：
- 每个阶段一个提取函数，读对应 summary（路径与 pipeline.PIPELINE_STAGES
  一致）；产物缺失 → value=None，该层跳过而不是报 0；
- 全部层都缺 → build_funnel() 返回 None（调用方整节省略）；
- 漏斗层宽度按 log10(记录数) 归一——reads(10^6) 到 contigs(10^1) 有六个
  数量级，线性比例会让下层完全不可见；图上明确标注口径。

CLI 单跑：
    python -m Virus_Platform_Core.stage_funnel <样本目录> [--outdir DIR]
"""
from __future__ import annotations

import json
import os
import sys

from .io_utils import fnum, safe_plot
from . import sci_palette as sp

# ------------------------------------------------------------------
# 各阶段 summary 读取
# ------------------------------------------------------------------


def _load_summary(sample_dir, stage_dir, fname):
    p = os.path.join(sample_dir, stage_dir, fname)
    if not os.path.isfile(p):
        return None, p
    try:
        with open(p, encoding='utf-8', errors='replace') as fh:
            return json.load(fh), p
    except Exception:
        return None, p


def _fastp_reads(s, key_before='reads_before', key_after='reads_after'):
    """fastp 读数：优先平台回写字段，退回原生 report 的 summary 区。"""
    if s.get(key_before) is not None:
        return fnum(s.get(key_before))
    before = ((s.get('summary') or {}).get('before_filtering') or {})
    return fnum(before.get('total_reads'))


def _orf_count(s):
    n = fnum(s.get('n_orfs'))
    if n is None and isinstance(s.get('orfipy'), dict):
        n = fnum(s['orfipy'].get('n_orfs'))
    if n is None:
        n = 0.0
        for eng in ('pyrodigal', 'pyrodigal_rv'):
            n += fnum((s.get(eng) or {}).get('genes'), 0.0) or 0.0
        n = n or None
    return n


# ------------------------------------------------------------------
# 漏斗链：每层 (key, 中文名, 英文名, 提取函数)
# 提取函数 → (value, note, src_relpath)；value None = 产物缺失，整层跳过
# ------------------------------------------------------------------
def _x_input_reads(sd):
    s, src = _load_summary(sd, '00_prep', 'fastp_report.json')
    if s:
        v = _fastp_reads(s)
        if v is not None:
            return v, '', os.path.join('00_prep', 'fastp_report.json')
    s, src = _load_summary(sd, '01_host_removal', 'stats.json')
    if s and fnum(s.get('total_pairs')) is not None:
        # 无 fastp 时用宿主去除的对数回退。单位必须与下一层（去宿主后
        # reads = kept×2）一致：双端 ×2 折成 reads，单端 kept_pairs 本就是
        # 单条 reads（host_removal: kept_pairs=kept1=R1 记录数）不再乘。
        # 此前双端回退只记对数、下一层记条数——同一漏斗混用单位，
        # 留存率虚高 2 倍（2026-10-13 真实全流程实抓）。
        pairs = fnum(s.get('total_pairs'))
        se = bool(s.get('single_end'))
        # note 不要自带括号：漏斗图会把 note 再包一层「（…）」，
        # 自带就出「（（对×2 估算））」双括号（2026-10-13 实抓）。
        return (pairs if se else pairs * 2, '' if se else '对×2 估算',
                os.path.join('01_host_removal', 'stats.json'))
    return None, '', ''


def _x_qc_reads(sd):
    s, src = _load_summary(sd, '00_prep', 'fastp_report.json')
    if not s:
        return None, '', ''
    v = fnum(s.get('reads_after'))
    if v is None:
        v = fnum(((s.get('summary') or {}).get('after_filtering') or {})
                 .get('total_reads'))
    if v is None:
        return None, '', ''
    q30 = s.get('q30_rate')
    note = f"Q30 {q30}%" if q30 is not None else ''
    return v, note, os.path.join('00_prep', 'fastp_report.json')


def _x_host_free_reads(sd):
    s, src = _load_summary(sd, '01_host_removal', 'stats.json')
    if not s:
        return None, '', ''
    kept = fnum(s.get('kept_pairs'))
    if kept is None:
        return None, '', ''
    ratio = fnum(s.get('host_ratio'))
    note = (f"宿主占比 {ratio * 100:.2f}%" if ratio is not None else '')
    # 单位与输入层一致：双端 kept_pairs 是「对」→ ×2 折条；单端 kept_pairs
    # 本就是单条 reads（host_removal 的 kept1=R1 记录数）→ 不乘。
    se = bool(s.get('single_end'))
    return (kept if se else kept * 2, note,
            os.path.join('01_host_removal', 'stats.json'))


def _x_viral_reads(sd):
    s, src = _load_summary(sd, '02b_kvsuite', 'summary.json')
    if not s:
        return None, '', ''
    v = fnum(s.get('viral_reads'))
    if v is None:
        return None, '', ''
    n_sp = len(s.get('species_detected') or [])
    note = f"检出物种 {n_sp} 个" if n_sp else ''
    return v, note, os.path.join('02b_kvsuite', 'summary.json')


def _x_contigs(sd):
    s, src = _load_summary(sd, '03_assembly', 'summary.json')
    if not s:
        return None, '', ''
    v = fnum(s.get('contigs'))
    if v is None:
        return None, '', ''
    return v, '', os.path.join('03_assembly', 'summary.json')


def _x_viral_contigs(sd):
    s, src = _load_summary(sd, '03_assembly', 'summary.json')
    if not s:
        return None, '', ''
    v = s.get('viral_contigs')
    if v is None:
        return None, '', ''
    return float(len(v)), '', os.path.join('03_assembly', 'summary.json')


def _x_verified(sd):
    s, src = _load_summary(sd, '03b_verify', 'summary.json')
    if not s:
        return None, '', ''
    v = fnum(s.get('n_calls'))
    if v is None:
        return None, '', ''
    calls = s.get('calls') or {}
    order = ['known', 'novel', 'viroid', 'unclassified']
    note = '；'.join(f"{k} {calls[k]}" for k in order if k in calls)
    return v, note, os.path.join('03b_verify', 'summary.json')


def _x_consensus(sd):
    s, src = _load_summary(sd, '03c_consensus', 'summary.json')
    if not s:
        return None, '', ''
    v = fnum(s.get('n_refs'))
    if v is None:
        return None, '', ''
    variants = fnum(s.get('n_variants'))
    note = f"变异位点 {int(variants)}" if variants is not None else ''
    return v, note, os.path.join('03c_consensus', 'summary.json')


def _x_orf(sd):
    s, src = _load_summary(sd, '04_orf', 'summary.json')
    if not s:
        return None, '', ''
    v = _orf_count(s)
    if v is None:
        return None, '', ''
    return v, '', os.path.join('04_orf', 'summary.json')


FUNNEL_CHAIN = [
    ('input_reads', '输入 reads', 'Raw reads', _x_input_reads),
    ('qc_reads', '质控后 reads', 'QC-passed reads', _x_qc_reads),
    ('host_free_reads', '去宿主后 reads', 'Host-free reads', _x_host_free_reads),
    ('viral_reads', '病毒 reads', 'Viral reads', _x_viral_reads),
    ('contigs', '组装 contigs', 'Contigs', _x_contigs),
    ('viral_contigs', '病毒 contigs', 'Viral contigs', _x_viral_contigs),
    ('verified', '验证判定序列', 'Verified calls', _x_verified),
    ('consensus_refs', '共识序列', 'Consensus seqs', _x_consensus),
    ('orf_genes', 'ORF 基因', 'ORF genes', _x_orf),
]

# ------------------------------------------------------------------
# KPI 卡（非"记录衰减"性质的支路指标）
# ------------------------------------------------------------------


def kpi_stats(sample_dir):
    """支路 KPI → [(label, value_str, note)]；缺产物的条目自动省略。"""
    out = []
    s, _ = _load_summary(sample_dir, '02b_kvsuite', 'summary.json')
    if s:
        n_sp = len(s.get('species_detected') or [])
        if n_sp or fnum(s.get('viral_reads')) is not None:
            out.append(('已知病毒检出', f"{n_sp} 个物种",
                        f"通过 {s.get('passed', 0)} / 剔除 {s.get('rejected', 0)}"))
    s, _ = _load_summary(sample_dir, '03c_consensus', 'summary.json')
    if s and fnum(s.get('n_variants')) is not None:
        out.append(('变异位点', f"{int(fnum(s.get('n_variants')))}",
                    f"SNV {s.get('n_snv', 0)} / iSNV {s.get('n_isnv', 0)}"))
    s, _ = _load_summary(sample_dir, '08_host_analysis', 'summary.json')
    if s:
        cats = s.get('categories') or {}
        n = fnum(s.get('n_contigs'), len(cats)) or len(cats)
        if n:
            out.append(('宿主预测', f"{int(n)} 条 contigs",
                        '；'.join(f"{k} {v}" for k, v in
                                  list(cats.items())[:3])))
    s, _ = _load_summary(sample_dir, '04b_orf_annot', 'summary.json')
    if s and fnum(s.get('n_annotated')) is not None:
        out.append(('ORF 功能注释',
                    f"{int(fnum(s.get('n_annotated')))} / "
                    f"{int(fnum(s.get('n_orfs'), 0) or 0)}",
                    f"{s.get('n_families', 0)} 个病毒科"))
    s, _ = _load_summary(sample_dir, '05_phylo', 'summary.json')
    if s and s.get('groups') is not None:
        out.append(('进化分析组', f"{len(s.get('groups') or [])} 个", ''))
    s, _ = _load_summary(sample_dir, '06_primer', 'summary.json')
    if s and fnum(s.get('n_primers')) is not None:
        out.append(('引物', f"{int(fnum(s.get('n_primers')))} 对", ''))
    s, _ = _load_summary(sample_dir, '09_genome_plots', 'summary.json')
    if s and s.get('plots') is not None:
        eng = (s.get('engine') or 'gbdraw')
        eng = 'DFV' if 'dna_features_viewer' in eng else eng
        out.append(('基因组图', f"{len(s.get('plots') or [])} 张 SVG", eng))
    return out


# ------------------------------------------------------------------
# 汇总
# ------------------------------------------------------------------


def build_funnel(sample_dir, out_dir=None, logger=None):
    """汇总漏斗 → dict(rows, kpi, tsv)；全缺产物返回 None。

    rows: 每层 dict(key, stage_zh, stage_en, value, note, src)。
    out_dir 给定时写 stage_funnel.tsv（缺省不落盘，由调用方决定）。
    """
    sample_dir = str(sample_dir)
    rows = []
    for key, zh, en, fn in FUNNEL_CHAIN:
        try:
            value, note, src = fn(sample_dir)
        except Exception as e:  # 单层提取失败不拖垮漏斗
            if logger:
                logger.log(f"漏斗层 {key} 提取失败（跳过）: {e}", 'WARN')
            value, note, src = None, '', ''
        rows.append({'key': key, 'stage_zh': zh, 'stage_en': en,
                     'value': value, 'note': note, 'src': src})

    if not any(r['value'] is not None for r in rows):
        return None

    kpi = kpi_stats(sample_dir)

    tsv_path = None
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
        tsv_path = os.path.join(out_dir, 'stage_funnel.tsv')
        with open(tsv_path, 'w', encoding='utf-8') as fh:
            fh.write('stage_key\tstage\tvalue\tnote\tsource\n')
            for r in rows:
                v = '' if r['value'] is None else (
                    int(r['value']) if float(r['value']).is_integer()
                    else r['value'])
                fh.write(f"{r['key']}\t{r['stage_zh']}\t{v}\t"
                         f"{r['note']}\t{r['src']}\n")
    return {'rows': rows, 'kpi': kpi, 'tsv': tsv_path}


# ------------------------------------------------------------------
# 绘图：梯形漏斗 + 留存率
# ------------------------------------------------------------------


def _luminance(hex_or_rgba):
    """相对亮度（0 黑 1 白），用于标签黑白切换。"""
    if isinstance(hex_or_rgba, (tuple, list)):
        r, g, b = hex_or_rgba[:3]
    else:
        h = hex_or_rgba.lstrip('#')
        r, g, b = (int(h[i:i + 2], 16) / 255.0 for i in (0, 2, 4))
    return 0.2126 * r + 0.7152 * g + 0.0722 * b


def _clear_stale(out_base):
    """出图前清掉同基础名的陈旧产物——失败时不允许残留旧图冒充新结果。"""
    for ext in ('.png', '.svg', '.pdf'):
        p = f'{out_base}{ext}'
        if os.path.isfile(p):
            try:
                os.remove(p)
            except OSError:
                pass


def plot_funnel(result, out_base, logger=None):
    """梯形漏斗图（层宽按 log10 记录数归一）→ out_base.png/.svg/.pdf。"""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon
    sp.apply_sci_style()
    _clear_stale(out_base)

    layers = [r for r in result['rows'] if r['value'] is not None]
    n = len(layers)
    if n < 2:
        if logger:
            logger.log('漏斗有效层不足 2 层，跳过漏斗图', 'WARN')
        return None
    if n > sp.DENSE_ELEM_THRESHOLD:  # 防御性：当前链条最多 9 层
        layers = layers[:sp.DENSE_ELEM_THRESHOLD]
        n = len(layers)

    vals = [max(float(r['value']), 0.0) for r in layers]
    import math
    logs = [math.log10(v + 1.0) for v in vals]
    lo, hi = min(logs), max(logs)
    span = (hi - lo) or 1.0

    fig, ax = plt.subplots(figsize=(7.2, 0.62 * n + 1.2))
    colors = sp.grade_colors(n, cmap='Blues', lo=0.30, hi=0.92)
    total_h, layer_h, gap = 1.0, 0.82, 0.18
    y_top_full = n * (layer_h + gap)

    for i, (row, lg, c) in enumerate(zip(layers, logs, colors)):
        y1 = y_top_full - i * (layer_h + gap)          # 层顶
        y0 = y1 - layer_h                              # 层底
        w = 0.16 + 0.84 * (lg - lo) / span             # log10 归一宽度
        x0, x1 = -w / 2, w / 2
        ax.add_patch(plt.Rectangle((x0, y0), w, layer_h,
                                   facecolor=c, edgecolor='white', lw=1.2))
        if i + 1 < n:                                   # 到下一层的连接带
            lg2 = logs[i + 1]
            w2 = 0.16 + 0.84 * (lg2 - lo) / span
            ax.add_patch(Polygon([(x0, y0), (x1, y0),
                                  (w2 / 2, y0 - gap), (-w2 / 2, y0 - gap)],
                                 closed=True, facecolor=c, alpha=0.28,
                                 edgecolor='none'))
        label = f"{(row.get('stage_en') or row['stage_zh'])}  {int(row['value']):,}"
        if row['note']:
            label += f"（{row['note']}）"
        txt_color = 'white' if _luminance(c) < 0.55 else '#16233d'
        if w >= 0.5:                      # 层足够宽：文字居层中
            ax.text(0, (y0 + y1) / 2, label, ha='center', va='center',
                    fontsize=9.5, color=txt_color)
        else:                             # 窄层：文字放到层右侧，防裁剪
            ax.text(x1 + 0.02, (y0 + y1) / 2, label, ha='left', va='center',
                    fontsize=9.5, color='#16233d')

    ax.set_xlim(-0.72, 0.72)
    ax.set_ylim(-gap - 0.14, y_top_full + 0.04)
    ax.axis('off')
    ax.set_title('Pipeline funnel (width proportional to log10 records)',
                 fontsize=11, fontweight='bold', pad=10)
    fig.tight_layout()
    png, svg, pdf = (f'{out_base}.png', f'{out_base}.svg',
                     f'{out_base}.pdf')
    fig.savefig(png, dpi=300)
    try:
        fig.savefig(svg)
    except Exception:
        pass
    fig.savefig(pdf)
    plt.close(fig)
    return png


def plot_retention(result, out_base, logger=None):
    """留存率横条：各层相对输入（及相对上一层）的百分比 → png/svg/pdf。"""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    sp.apply_sci_style()
    _clear_stale(out_base)

    layers = [r for r in result['rows'] if r['value'] is not None]
    if len(layers) < 2:
        return None
    base = float(layers[0]['value'])
    if base <= 0:
        return None

    labels, rel_base, rel_prev, notes = [], [], [], []
    prev = None
    for r in layers:
        v = float(r['value'])
        labels.append(r['stage_zh'])
        rel_base.append(v / base * 100.0)
        rel_prev.append((v / prev * 100.0) if prev else 100.0)
        notes.append(r['note'])
        prev = v

    def _pct(x):
        """百分比人性化格式：1e+02% 这类科学计数绝对不能出现在图上；
        非零但四舍五入后为 0 的极小值（如 2/7200≈0.03%）显示 <0.1%，
        不显示误导性的 0%（2026-10-13 目检 fixture 实抓）。"""
        if 0 < x < 0.05:
            return '<0.1%'
        s = f'{x:.1f}'.rstrip('0').rstrip('.')
        return f'{s}%'

    y = list(range(len(labels)))[::-1]
    fig, ax = plt.subplots(figsize=(7.6, 0.52 * len(labels) + 1.4))
    ax.barh(y, rel_base, height=0.62, color=sp.OKABE['blue'],
            edgecolor='white', zorder=3)
    # 文字统一放条形右端外侧（条内放长文本必然溢出，2026-10-12 目检教训）
    for yi, rb, rp, r, nt in zip(y, rel_base, rel_prev, layers, notes):
        main = _pct(rb) + f"（环比 {_pct(rp)}）"
        if nt:
            main += f' · {nt}'
        ax.text(rb + 1.5, yi, main, va='center', ha='left',
                fontsize=8.3, color='#333333')
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.set_xlabel('Retention vs input records (%)')
    ax.set_xlim(0, max(115.0, max(rel_base) * 1.15))
    ax.xaxis.grid(True, ls='--', alpha=0.5, color='#d9d9d9', zorder=0)
    ax.set_axisbelow(True)
    sp.despine(ax)
    ax.set_title('Stage retention', fontsize=11, fontweight='bold', pad=10)
    fig.tight_layout()
    png, svg, pdf = (f'{out_base}.png', f'{out_base}.svg',
                     f'{out_base}.pdf')
    fig.savefig(png, dpi=300)
    try:
        fig.savefig(svg)
    except Exception:
        pass
    fig.savefig(pdf)
    plt.close(fig)
    return png


def render(sample_dir, out_dir=None, logger=None):
    """一键入口：汇总 + 双图。返回 dict(funnel=…, funnel_png=…, retention_png=…)
    或 None（全缺产物 / 有效层不足）。"""
    out_dir = out_dir or os.path.join(str(sample_dir), '07_report')
    result = build_funnel(sample_dir, out_dir=out_dir, logger=logger)
    if not result:
        return None
    base = os.path.join(out_dir, 'stage_funnel')
    # 图函数内部不再传 logger：safe_plot 已负责失败告警
    safe_plot('流程漏斗图', logger, plot_funnel, result, base)
    safe_plot('留存率图', logger, plot_retention, result,
              os.path.join(out_dir, 'stage_funnel_retention'))
    result['png'] = f'{base}.png' if os.path.isfile(f'{base}.png') else None
    return result


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser(
        description='样本管线流程漏斗：TSV + 梯形漏斗图 + 留存率图')
    ap.add_argument('sample_dir', help='样本结果目录（含各阶段阶段目录）')
    ap.add_argument('--outdir', default=None, help='输出目录（缺省 07_report/）')
    a = ap.parse_args()
    res = render(a.sample_dir, out_dir=a.outdir,
                 logger=type('L', (), {'log': staticmethod(
                     lambda m, lv='INFO': print(f'[{lv}] {m}'))})())
    if not res:
        print('没有可用的阶段产物，未生成漏斗。')
        sys.exit(1)
    print('stage_funnel.tsv:', res.get('tsv'))
    print('funnel png:', res.get('png'))
