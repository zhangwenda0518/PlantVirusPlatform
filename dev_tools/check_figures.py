# -*- coding: utf-8 -*-
"""图形质量自动体检：文字重叠 + 越界检测（MMPV 蓝本 §1.6 的移植）。

把"图好不好"从人工目检变成可回归判据：劫持 matplotlib 的 Figure.savefig，
每次保存前 draw 画布、收集全部可见 Text 图元的窗口 bbox：
- 任意两两交叠面积 > OVERLAP_MIN_AREA → FAIL（文字重叠）；
- bbox 超出画布边界（TOL 0.5px）→ FAIL（文字越界）；
- 交叠面积 0.5~OVERLAP_MIN_AREA → WARN（贴边，不判失败）。

退出码：0 = 全部通过；2 = 存在 FAIL。

用法一：体检任一作图脚本（脚本收到的 argv 与直跑完全一致）
    python dev_tools/check_figures.py run -- path/to/script.py [args...]

用法二：库内使用（测试 / 作图模块自带守卫）
    import sys; sys.path.insert(0, os.path.join(_ROOT, 'dev_tools'))
    from check_figures import fig_quality_guard
    with fig_quality_guard('my_module') as rep:
        make_figure(...)
    assert not rep.fails, rep.summary()

用法三：内置演示（对平台漏斗三件套 fixture 跑一遍）
    python dev_tools/check_figures.py demo
"""
from __future__ import annotations

import os
import sys

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

import matplotlib
matplotlib.use('Agg')
import matplotlib.figure
import matplotlib.text

# 判据常量（px²）：>2 判 FAIL，>0.5 只记 WARN（亚像素噪声不抓）
OVERLAP_MIN_AREA = 2.0
WARN_MIN_AREA = 0.5
EDGE_TOL = 0.5


class FigureReport:
    """一次体检会话的违规收集器。"""

    def __init__(self):
        self.violations = []   # (level, figure_label, message)

    def add(self, level, label, message):
        self.violations.append((level, label, message))

    @property
    def fails(self):
        return [v for v in self.violations if v[0] == 'FAIL']

    @property
    def warns(self):
        return [v for v in self.violations if v[0] == 'WARN']

    def summary(self):
        if not self.violations:
            return '全部通过：0 重叠 / 0 越界'
        lines = [f'{len(self.fails)} FAIL / {len(self.warns)} WARN']
        for lv, label, msg in self.violations:
            lines.append(f'  [{lv}] {label}: {msg}')
        return '\n'.join(lines)


def _dedup_key(text_artist):
    """同字符串同位置的双绘 Text 去重（防止重复 draw 造成的假阳性）。"""
    try:
        return (text_artist.get_text(),
                tuple(round(v, 1) for v in text_artist.get_position()))
    except Exception:
        return None


def _axis_off_hidden_texts(fig):
    """收集"被 ax.axis('off') 隐藏"的 Text 实例集合。

    坑：ax.axis('off') 后刻度标签的 Text.get_visible() 仍为 True，但实际
    不渲染（刻度标签 Text 的 .axes 属性是 None，也无法从 Text 反查所属轴）
    ——只能正向遍历 fig.axes 收集。（tick_params(labelleft=False) 的情形
    在 mpl≥3.11 会正确落到 Text.get_visible()=False，无需在此处理。）
    """
    hidden = set()
    for ax in fig.axes:
        if getattr(ax, 'axison', True):
            continue
        for axis in (ax.xaxis, ax.yaxis):
            hidden.update(axis.get_ticklabels())
            hidden.update(axis.get_minorticklabels())
            hidden.add(axis.label)
            hidden.add(axis.get_offset_text())
    return hidden


def check_figure(fig: matplotlib.figure.Figure, label: str,
                 report: FigureReport, check_bounds: bool = True) -> None:
    """对单张 Figure 做"两两文字重叠 + 越界"体检。

    check_bounds=False 时跳过越界判定（bbox_inches='tight' 保存的图
    画布随内容生长，测量期越界会被裁剪救回）。
    """
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    fig_bb = fig.bbox
    hidden = _axis_off_hidden_texts(fig)

    texts = []
    seen = set()
    for t in fig.findobj(matplotlib.text.Text):
        s = ''
        try:
            s = (t.get_text() or '').strip()
        except Exception:
            continue
        if not s or not t.get_visible() or t in hidden:
            continue
        key = _dedup_key(t)
        if key in seen:
            continue
        seen.add(key)
        try:
            bb = t.get_window_extent(renderer=renderer)
        except Exception:
            continue
        if bb.width <= 0 or bb.height <= 0:
            continue
        texts.append((s, bb))

    # 越界：文字 bbox 超出画布（tight 保存时跳过——画布将随内容生长）
    if check_bounds:
        for s, bb in texts:
            if (bb.x0 < fig_bb.x0 - EDGE_TOL or bb.y0 < fig_bb.y0 - EDGE_TOL
                    or bb.x1 > fig_bb.x1 + EDGE_TOL
                    or bb.y1 > fig_bb.y1 + EDGE_TOL):
                report.add('FAIL', label,
                           f'文字越界: "{s[:40]}" bbox=({bb.x0:.0f},{bb.y0:.0f},'
                           f'{bb.x1:.0f},{bb.y1:.0f}) 画布=({fig_bb.x0:.0f},'
                           f'{fig_bb.y0:.0f},{fig_bb.x1:.0f},{fig_bb.y1:.0f})')

    # 两两重叠
    n = len(texts)
    for i in range(n):
        s1, ib = texts[i]
        for j in range(i + 1, n):
            s2, jb = texts[j]
            ix = min(ib.x1, jb.x1) - max(ib.x0, jb.x0)
            iy = min(ib.y1, jb.y1) - max(ib.y0, jb.y0)
            if ix <= 0 or iy <= 0:
                continue
            area = ix * iy
            if area > OVERLAP_MIN_AREA:
                report.add('FAIL', label,
                           f'文字重叠({area:.0f}px²): "{s1[:30]}" × "{s2[:30]}"')
            elif area > WARN_MIN_AREA:
                report.add('WARN', label,
                           f'文字贴边({area:.1f}px²): "{s1[:30]}" × "{s2[:30]}"')


class fig_quality_guard:
    """上下文管理器：会话内所有 savefig 先体检后保存。

    label 缺省用保存文件名。原 savefig 行为完全不变（体检不落盘任何东西）。
    """

    def __init__(self, label: str = ''):
        self.label = label
        self.report = FigureReport()
        self._orig = None

    def __enter__(self):
        self._orig = matplotlib.figure.Figure.savefig

        def _patched(fig_self, fname, *a, **kw):
            try:
                lab = self.label or os.path.basename(str(fname))
                # bbox_inches='tight' 会在保存时按内容重新定画布——测量期
                # 的"越界"由 tight 裁剪救回，不构成缺陷；重叠判定保留
                # （tight 只裁边不改图元相对位置）。
                tight = kw.get('bbox_inches', None) == 'tight'
                check_figure(fig_self, lab, self.report,
                             check_bounds=not tight)
            except Exception as e:  # 体检自身故障不得阻断出图
                self.report.add('WARN', self.label or str(fname),
                                f'体检器异常（跳过该图检查）: {e}')
            return self._orig(fig_self, fname, *a, **kw)

        matplotlib.figure.Figure.savefig = _patched
        return self.report

    def __exit__(self, *exc):
        matplotlib.figure.Figure.savefig = self._orig
        return False


def run_script(script_path, script_args):
    """以体检模式跑一个作图脚本（runpy，argv 原样透传）。"""
    import runpy
    old_argv = sys.argv
    sys.argv = [script_path] + list(script_args)
    try:
        with fig_quality_guard(os.path.basename(script_path)) as rep:
            try:
                runpy.run_path(script_path, run_name='__main__')
            except SystemExit as e:      # 脚本自带 sys.exit 属正常流程
                if e.code not in (0, None):
                    raise
    finally:
        sys.argv = old_argv
    return rep


def main(argv):
    if not argv or argv[0] in ('-h', '--help'):
        print(__doc__)
        return 0
    if argv[0] == 'demo':
        demo = os.path.join(_ROOT, 'tests', 'fixtures', 'sample_demo')
        from Virus_Platform_Core import stage_funnel
        with fig_quality_guard('stage_funnel(fixture)') as rep:
            stage_funnel.render(demo, out_dir=os.path.join(demo, '07_report'))
        print(rep.summary())
        return 2 if rep.fails else 0

    args = list(argv)
    if args[0] == 'run':
        args = args[1:]
    if args and args[0] == '--':     # 允许 run -- script 形式
        args = args[1:]
    if not args:
        print(__doc__)
        return 0
    script = args[0]
    if not os.path.isfile(script):
        print(f'脚本不存在: {script}')
        return 1
    rep = run_script(script, args[1:])
    print(rep.summary())
    return 2 if rep.fails else 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
