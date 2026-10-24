# -*- coding: utf-8 -*-
"""图形体检扫描守卫：对四个可复跑绘图模块跑 check_figures 质检。

覆盖（2026-10-12 起全部纳入体检回归）：
1. stage_funnel 漏斗三件套（tests/fixtures/sample_demo）
2. kv_fig 丰度图（examples/results/kvsuite 的 filtered.tsv）
3. kv_variant_plot 变异棒棒糖图（examples/results/variant/kvsuite，
   演示注释由 dev_tools/make_variant_demo.py 经真 snpEff 链生成；
   annotated/ 缺失时自动重建，重建失败按 SKIP 处理）
4. landscape_plot meta 六面板（figures/assets/meta_search 真实 CSV）

判据：每次 savefig 前的 Text 两两重叠/越界检查 0 FAIL（与
dev_tools/check_figures.py 同一引擎）。输出遵循 _check_* 脚本约定：
通过退出 0，失败退出 1，依赖缺失退出 0 并打 SKIP。
"""
import os
import sys
import tempfile
from pathlib import Path

_ROOT = Path(__file__).resolve().parents[1]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))
_DEV = str(_ROOT / 'dev_tools')
if _DEV not in sys.path:
    sys.path.insert(0, _DEV)

import logging

from check_figures import fig_quality_guard  # noqa: E402

FIXTURE_SAMPLE = _ROOT / 'tests' / 'fixtures' / 'sample_demo'
KVSUITE_EX = _ROOT / 'examples' / 'results' / 'kvsuite' / 'kvsuite'
VARIANT_EX = _ROOT / 'examples' / 'results' / 'variant' / 'kvsuite'
META_CSV = _ROOT / 'figures' / 'assets' / 'meta_search' / 'SRA_GSA_Merged_Final.csv'

FAILS = []


def _scan(name, fn):
    """跑一个扫描目标，FAIL 收进 FAILS，异常记 FAIL 不中断其余目标。"""
    try:
        with fig_quality_guard(name) as rep:
            fn()
        n_fail = len(rep.fails)
        if n_fail:
            FAILS.extend(f'{name}: {m}' for _lv, _lab, m in rep.fails)
            print(f'  ✘ {name}: {n_fail} FAIL')
            for _lv, _lab, m in rep.fails[:6]:
                print(f'      · {m}')
        else:
            print(f'  ✓ {name}: 0 重叠 / 0 越界')
    except Exception as e:  # noqa: BLE001
        FAILS.append(f'{name}: 体检运行异常 {e}')
        print(f'  ✘ {name}: 异常 {e}')


def scan_funnel():
    from Virus_Platform_Core import stage_funnel
    res = stage_funnel.build_funnel(str(FIXTURE_SAMPLE))
    if not res:
        raise RuntimeError('fixture 漏斗数据为空（fixtures 被动过？）')
    out = Path(tempfile.mkdtemp(prefix='_figqa_funnel_'))
    stage_funnel.plot_funnel(res, str(out / 'qa_funnel'))
    stage_funnel.plot_retention(res, str(out / 'qa_retention'))


def scan_kv_fig():
    from Virus_Platform_Core.known_virus_suite.kv_fig import generate_figures
    out = Path(tempfile.mkdtemp(prefix='_figqa_kv_'))
    generate_figures(str(KVSUITE_EX / 'filter' / 'filtered.tsv'), str(out))


def _variant_ann_ready():
    ann = VARIANT_EX / 'annotated'
    return ann.is_dir() and any(ann.glob('*.ann.tsv'))


def scan_variant_lollipop():
    if not _variant_ann_ready():
        print('  … annotated/ 缺失，自动重建演示注释（真 snpEff 链）…')
        import runpy
        try:
            runpy.run_path(str(_ROOT / 'dev_tools' / 'make_variant_demo.py'),
                           run_name='__main__')
        except SystemExit as e:
            if e.code not in (0, None):
                raise
    if not _variant_ann_ready():
        print('  SKIP 演示注释不可用（Java/snpEff 缺失？）——'
              '可手动运行 python dev_tools/make_variant_demo.py')
        return
    from Virus_Platform_Core.known_virus_suite.kv_variant_plot import (
        run_variant_plots)
    run_variant_plots(str(VARIANT_EX), logger=logging.getLogger('figqa'))


def scan_meta_landscape():
    if not META_CSV.is_file():
        print(f'  SKIP 缺输入 CSV: {META_CSV}')
        return
    import runpy
    out = Path(tempfile.mkdtemp(prefix='_figqa_meta_'))
    old_argv = sys.argv
    sys.argv = ['landscape_plot.py', '-i', str(META_CSV), '-o', str(out)]
    try:
        try:
            runpy.run_path(
                str(_ROOT / 'Virus_Platform_Core' / 'public_meta'
                    / 'landscape_plot.py'), run_name='__main__')
        except SystemExit as e:
            if e.code not in (0, None):
                raise
    finally:
        sys.argv = old_argv


def main():
    logging.basicConfig(level=logging.WARNING)
    print('== 图形体检扫描（check_figures 引擎）==')
    _scan('stage_funnel(fixture)', scan_funnel)
    _scan('kv_fig(example)', scan_kv_fig)
    _scan('kv_variant_plot(variant demo)', scan_variant_lollipop)
    _scan('landscape_plot(meta csv)', scan_meta_landscape)
    print()
    if FAILS:
        print(f'✘ 图形体检失败（{len(FAILS)} 条）')
        return 1
    print('✓ 图形体检通过：4 个扫描目标 0 重叠 / 0 越界')
    return 0


if __name__ == '__main__':
    sys.exit(main())
