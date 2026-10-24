# -*- coding: utf-8 -*-
"""科研配色与绘图样式唯一真源（single source of truth）。

借鉴 MMPV-RNA 平台的 sci_palette 体系（蓝本见
docs/MMPV-RNA结果解析与绘图学习笔记.md §1），本文件是全平台唯一的
配色/样式定义点。规则：

1. 业务色一律从 OKABE_ITO 取色相；派生色（深浅变体）只允许经
   lighten() 生成，模块里禁止新增游离 hex。此前 viz.py 用 NPG 十色、
   dfv_plot 与 genome_diag 各持一套 CATEGORY_COLORS（同语义不同色）、
   logan_trace 用 Set1、kv_* 各自手抄 rcParams——全部以本文件为准。
2. 同一语义跨图同色：判定、Impact、突变类型、功能类别等查本文件的
   语义字典（SEMANTIC 区），不要在模块里自建。
3. 样式只经 apply_sci_style() / sci_rcparams() 设置。不要在模块里
   手抄 rcParams——手抄必漏键（漏 spines/fonttype 会让同一份代码
   单独跑与跟跑产出 dpi、字体嵌入不一致的图）。
4. 连续数值用 SEQ_CMAP，发散数值用 DIV_CMAP，分级色带用 grade_colors()。
"""
from __future__ import annotations

import colorsys

# ------------------------------------------------------------------
# 主色板：Okabe-Ito（Wong 2011, Nature Methods，色盲安全），顺序固定
# ------------------------------------------------------------------
OKABE_ITO = ('#E69F00', '#56B4E9', '#009E73', '#F0E442',
             '#0072B2', '#D55E00', '#CC79A7', '#000000', '#999999')

OKABE = {
    'orange': '#E69F00',
    'sky': '#56B4E9',
    'green': '#009E73',
    'yellow': '#F0E442',
    'blue': '#0072B2',
    'vermillion': '#D55E00',
    'purple': '#CC79A7',
    'black': '#000000',
    'grey': '#999999',
}

# 分支/细线专用：去掉最浅的 yellow、sky，保证 0.5pt 细线白底可辨
OKABE_BRANCH = ('#0072B2', '#D55E00', '#009E73', '#E69F00', '#CC79A7')

SEQ_CMAP = 'viridis'    # 连续数值（丰度、深度、覆盖率……）
DIV_CMAP = 'RdBu_r'     # 发散数值（log2 倍数、相关性……）

# ------------------------------------------------------------------
# 栏宽常量：画布宽 = 最终排版宽度，印刷缩放≈1（否则字号被等比压小）。
# 本平台报告以网页呈现，这两个值用于论文级导出与"画布随元素放大"的上限参考。
# ------------------------------------------------------------------
COL_W_1 = 9.32   # 通栏（英寸）
COL_W_2 = 4.58   # 单栏（英寸）

# 元素超过该数量时自动关闭逐元素文字标注（防 23550px 宽的不可读巨图）
DENSE_ELEM_THRESHOLD = 60

# 画布尺寸上限（英寸）——超限说明该按 DENSE_ELEM_THRESHOLD 收敛而不是放大画布
MAX_FIG_W_IN = 20.0
MAX_FIG_H_IN = 28.0


# ------------------------------------------------------------------
# 颜色工具
# ------------------------------------------------------------------
def _clamp01(x):
    return 0.0 if x < 0.0 else (1.0 if x > 1.0 else x)


def lighten(hex_color: str, frac: float) -> str:
    """按比例调整明度：frac>0 变浅、frac<0 变深（Okabe-Ito 的唯一派生途径）。"""
    h = hex_color.lstrip('#')
    r, g, b = (int(h[i:i + 2], 16) / 255.0 for i in (0, 2, 4))
    hh, ll, ss = colorsys.rgb_to_hls(r, g, b)
    ll = _clamp01(ll + (1.0 - ll) * frac) if frac > 0 else _clamp01(ll * (1.0 + frac))
    r, g, b = colorsys.hls_to_rgb(hh, ll, ss)
    return '#{:02X}{:02X}{:02X}'.format(int(r * 255), int(g * 255), int(b * 255))


def seq_colors(n: int, start: int = 0) -> list:
    """从 OKABE_ITO 循环取 n 个离散色（分类系列着色的统一入口）。"""
    return [OKABE_ITO[(start + i) % len(OKABE_ITO)] for i in range(max(n, 0))]


def grade_colors(n: int, cmap: str = SEQ_CMAP, lo: float = 0.25, hi: float = 0.9) -> list:
    """连续色板取 n 档（漏斗层级、分级条带等），默认避开 viridis 过暗的底段。"""
    import numpy as np
    import matplotlib
    if n <= 0:
        return []
    cm = matplotlib.colormaps[cmap]   # mpl≥3.9 已移除 cm.get_cmap
    vals = np.linspace(lo, hi, n)
    return [cm(v) for v in vals]


def figsize_for(ncol: int, ratio: float = 0.62, width: float = None) -> tuple:
    """按栏数给画布尺寸：width 缺省时单栏/双栏/通栏自动对齐栏宽常量。"""
    if width is None:
        width = COL_W_2 if ncol <= 1 else (COL_W_1 if ncol <= 3 else COL_W_1 * 1.6)
    return (width, width * ratio)


# ------------------------------------------------------------------
# 语义色表（SEMANTIC）：同一含义跨图同色
# ------------------------------------------------------------------
# 判定/质检结论（verify、过滤、审计通用）
VERDICT_COLORS = {
    'PASS': OKABE['green'],
    'OK': OKABE['green'],
    'WARN': OKABE['orange'],
    'REVIEW': OKABE['orange'],
    'FAIL': OKABE['vermillion'],
    'SKIP': OKABE['grey'],
    'N/A': OKABE['grey'],
}

# SnpEff Impact 四级（沿用 kv_variant_plot 既定映射）
IMPACT_COLOR = {
    'HIGH': OKABE['vermillion'],
    'MODERATE': OKABE['orange'],
    'LOW': OKABE['sky'],
    'MODIFIER': OKABE['grey'],
}

# 突变类型（kv_variant_evo 分子突变谱）
MUTATION_CLASS_COLOR = {
    'Transition': OKABE['vermillion'],
    'Transversion': OKABE['blue'],
    'Indel': OKABE['orange'],
    'Complex': OKABE['grey'],
}

# 数据库来源（公共元数据六面板）
DB_COLORS = {'SRA': OKABE['blue'], 'GSA': OKABE['vermillion']}

# ORF 功能类别（genome_diag / dfv_plot / orf_annot 统一口径；
# 原两套字典对"结构蛋白"一个 #7ec8e3 一个 #e67e22，此处合并）。
# 派生色均由 lighten() 从 Okabe-Ito 生成，不引入外来色相。
CATEGORY_COLORS = {
    '聚合酶/复制相关': OKABE['blue'],
    '聚合酶辅因子': lighten(OKABE['blue'], 0.35),
    '蛋白酶/加工': OKABE['purple'],
    '结构蛋白': OKABE['orange'],
    '运动蛋白': OKABE['green'],
    '宿主互作/致病': OKABE['vermillion'],
    '基因组复制/维持': lighten(OKABE['green'], 0.30),
    '转录/翻译调控': lighten(OKABE['yellow'], -0.18),   # 黄加深，防白底发虚
    '其他功能蛋白': OKABE['grey'],
    '假想蛋白（功能未定）': lighten(OKABE['grey'], 0.40),
}
CATEGORY_FALLBACK = OKABE['grey']
CATEGORY_FALLBACK_LIGHT = lighten(OKABE['grey'], 0.40)

# 深度/均值等单指标线色（kv_plot 既定口径，原值即 Okabe-Ito，改为引用）
DEPTH_COLOR = OKABE['blue']
MEAN_COLOR = OKABE['vermillion']


def category_color(name: str) -> str:
    return CATEGORY_COLORS.get((name or '').strip(), CATEGORY_FALLBACK)


def verdict_color(name: str) -> str:
    return VERDICT_COLORS.get((name or '').strip().upper(), OKABE['grey'])


def impact_color(name: str) -> str:
    return IMPACT_COLOR.get((name or '').strip().upper(), OKABE['grey'])


# ------------------------------------------------------------------
# 样式真源：整份 rcParams（强制共享，模块不许手抄）
# ------------------------------------------------------------------
def sci_rcparams(base_size: float = 9.0) -> dict:
    """返回完整的 matplotlib rcParams 字典。

    font.family 直接给具体字体列表（不是 'sans-serif' 泛名）——matplotlib
    的逐字形回退只在 family 是字体列表时生效；这样英文用 Arial 字形、
    中文逐字回退到微软雅黑（平台报告图的中英混排两全）。
    """
    from cycler import cycler
    return {
        'font.family': ['Arial', 'Microsoft YaHei', 'SimHei', 'DejaVu Sans'],
        'font.sans-serif': ['Arial', 'Helvetica', 'Microsoft YaHei',
                            'SimHei', 'DejaVu Sans'],
        'font.size': base_size,
        'axes.titlesize': base_size + 1.5,
        'axes.labelsize': base_size,
        'xtick.labelsize': base_size - 0.5,
        'ytick.labelsize': base_size - 0.5,
        'legend.fontsize': base_size - 0.5,
        'axes.unicode_minus': False,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'axes.linewidth': 0.8,
        'axes.edgecolor': '#333333',
        'xtick.major.width': 0.8,
        'ytick.major.width': 0.8,
        'legend.frameon': False,
        'figure.facecolor': 'white',
        'savefig.facecolor': 'white',
        'savefig.dpi': 300,
        'pdf.fonttype': 42,      # TrueType 嵌入，PDF 文字可编辑（期刊要求）
        'ps.fonttype': 42,
        'svg.fonttype': 'none',  # SVG 保留可编辑文本
        'image.cmap': SEQ_CMAP,
        'axes.prop_cycle': cycler(color=list(OKABE_ITO)),
    }


def apply_sci_style(base_size: float = 9.0, **overrides) -> None:
    """把真源样式应用到全局 rcParams；overrides 可临时覆盖个别键。"""
    import matplotlib as mpl
    params = sci_rcparams(base_size)
    params.update(overrides)
    mpl.rcParams.update(params)


def despine(ax=None, top: bool = True, right: bool = True,
            left: bool = False, bottom: bool = False) -> None:
    """隐藏指定侧的脊线（默认上/右）。rc 已全局去顶右，此函数用于局部补刀。"""
    if ax is None:
        import matplotlib.pyplot as plt
        ax = plt.gca()
    for side, hide in (('top', top), ('right', right),
                       ('left', left), ('bottom', bottom)):
        if hide:
            ax.spines[side].set_visible(False)
