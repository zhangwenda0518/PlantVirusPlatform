# -*- coding: utf-8 -*-
"""sci_palette 配色/样式真源单测。

守护规则（docs/MMPV-RNA结果解析与绘图学习笔记.md §1）：
1. OKABE_ITO 九色顺序固定（业务代码按位置取色，顺序即契约）；
2. 语义色全部派生自 OKABE_ITO，不引入游离色相；
3. sci_rcparams() 是完整字典，fonttype 42 与 prop_cycle 必须在场；
4. 画布常量与分级色工具可用。
"""
import colorsys

from Virus_Platform_Core import sci_palette as sp


def test_okabe_ito_order_is_contract():
    assert sp.OKABE_ITO == ('#E69F00', '#56B4E9', '#009E73', '#F0E442',
                            '#0072B2', '#D55E00', '#CC79A7', '#000000',
                            '#999999')
    assert sp.OKABE['blue'] == '#0072B2'
    assert sp.OKABE['vermillion'] == '#D55E00'
    # 分支色板：去掉最浅的 yellow / sky，细线可辨
    assert '#F0E442' not in sp.OKABE_BRANCH
    assert '#56B4E9' not in sp.OKABE_BRANCH
    assert len(sp.OKABE_BRANCH) == 5


def test_semantic_colors_derived_from_okabe():
    allowed = set(sp.OKABE_ITO)

    def _hue(hex_color):
        h = hex_color.lstrip('#')
        r, g, b = (int(h[i:i + 2], 16) / 255.0 for i in (0, 2, 4))
        return colorsys.rgb_to_hls(r, g, b)[0]

    ok_hues = {_hue(c) for c in allowed}
    for dic in (sp.VERDICT_COLORS, sp.IMPACT_COLOR, sp.DB_COLORS,
                sp.MUTATION_CLASS_COLOR):
        for v in dic.values():
            assert v in allowed or min(
                abs(_hue(v) - h) for h in ok_hues) < 1e-9, (
                f'{dic.__name__ if hasattr(dic, "__name__") else dic} '
                f'出现游离色相: {v}')


def test_category_colors_cover_genome_diag_set():
    # genome_diag 的 10 类必须全部有着落（同语义跨图同色的前提）
    for cat in ('聚合酶/复制相关', '聚合酶辅因子', '蛋白酶/加工', '结构蛋白',
                '运动蛋白', '宿主互作/致病', '基因组复制/维持', '转录/翻译调控',
                '其他功能蛋白', '假想蛋白（功能未定）'):
        assert cat in sp.CATEGORY_COLORS
        assert sp.CATEGORY_COLORS[cat].startswith('#')
    # 未登录类别回退灰、结构蛋白 = orange（原 dfv/genome_diag 两套色并存的修复点）
    assert sp.category_color('不存在类别') == sp.CATEGORY_FALLBACK
    assert sp.category_color('结构蛋白') == sp.OKABE['orange']
    assert sp.impact_color('HIGH') == sp.OKABE['vermillion']


def test_sci_rcparams_is_complete():
    rc = sp.sci_rcparams()
    # 手抄必漏键——真源必须一次性给全（蓝本 §1.2 的教训）
    for key in ('pdf.fonttype', 'ps.fonttype', 'svg.fonttype',
                'axes.spines.top', 'axes.spines.right', 'axes.unicode_minus',
                'savefig.dpi', 'figure.facecolor', 'legend.frameon'):
        assert key in rc, f'sci_rcparams 缺键: {key}'
    assert rc['pdf.fonttype'] == 42
    # 逐字形回退要求 font.family 是具体字体列表（含中文回退），不是 'sans-serif' 泛名
    assert isinstance(rc['font.family'], (list, tuple))
    assert any('YaHei' in f for f in rc['font.family'])
    # prop_cycle 接管为 Okabe-Ito（防 matplotlib 默认 tab10 混入）
    assert [str(c) for c in rc['axes.prop_cycle'].by_key()['color']] == \
        list(sp.OKABE_ITO)


def test_helpers():
    assert len(sp.seq_colors(12)) == 12
    assert sp.seq_colors(3) == list(sp.OKABE_ITO[:3])
    g = sp.grade_colors(5, cmap='Blues')
    assert len(g) == 5
    assert sp.lighten('#000000', 1.0) == '#FFFFFF'
    assert sp.lighten('#FFFFFF', -1.0) == '#000000'
    w, h = sp.figsize_for(1)
    assert abs(w - sp.COL_W_2) < 1e-9 and h < w
