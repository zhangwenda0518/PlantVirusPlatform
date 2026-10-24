# -*- coding: utf-8 -*-
"""表格容错读取与绘图容错三件套。

收编此前散落 7 处的重复实现（viz.py / kv_stage.py / bsp.py /
phylodyn_kit.py / phylodyn_local.py / ncbi_submit.adopt_mmpv.py 各自的
_read_tsv / _read_table），统一容错语义：

- 文件不存在 / 空文件 → 返回 default（不抛异常，调用方按"缺产物"降级）；
- BOM、非 UTF-8 尾巴、``#`` 注释行、空行一律容忍；
- 分隔符自动嗅探（tab 优先）；短行补空串、多余列丢弃（保留首现表头宽度）；
- fnum() 对 ``-``/NA/``12.3%``/千分位逗号等脏值返回 default 而不是抛错。

safe_plot() 是批量出图的刹车：单图失败只 WARN 不中断，与 MMPV 蓝本
"降级不崩溃"一致。
"""
from __future__ import annotations

import os
import io
import csv

__all__ = ['NA_TOKENS', 'sniff_delim', 'read_rows', 'read_dicts',
           'fnum', 'safe_plot']

# 视为缺失数值的 token（fnum 与 NA 单元格判定共用）
NA_TOKENS = {'', '-', '--', 'NA', 'N/A', 'na', 'nan', 'NaN', 'None',
             'null', 'NULL', '.'}


def sniff_delim(sample: str) -> str:
    """按前若干行统计 tab 与逗号数量判分隔符（tab 优先）。"""
    return '\t' if sample.count('\t') >= sample.count(',') else ','


def read_rows(path, sep='auto', skip_comments=False, skip_pred=None,
              encoding='utf-8-sig'):
    """读表格 → ``(header, rows, delim)``。

    header 为 list[str]（原始大小写、已 strip）；rows 为 list[list[str]]，
    全部单元格已 strip。文件缺失/空 → ``(None, [], '')``。

    skip_comments: 跳过行首 ``#`` 的行（BEAST/脚本重排日志常见）。
    skip_pred(cols) -> True 时跳过该数据行（如 kv_stage 的 dtype 行）。
    """
    if not os.path.isfile(path):
        return (None, [], '')
    header, rows = None, []
    delim = sep if sep != 'auto' else ','
    with io.open(path, 'r', encoding=encoding, errors='replace', newline='') as fh:
        if sep == 'auto':
            sample = fh.read(8192)
            fh.seek(0)
            if not sample.strip():
                return (None, [], '')      # 空文件：无分隔符可言
            delim = sniff_delim(sample)
        for raw in csv.reader(fh, delimiter=delim):
            cols = [(c or '').strip() for c in raw]
            if header is None:
                if not any(cols):
                    continue
                if skip_comments and cols[0].startswith('#'):
                    continue          # 表头前的注释行直接跳过（bsp 旧口径）
                header = [c.lstrip('#').strip() if skip_comments else c
                          for c in cols]
                continue
            if not any(cols):
                continue
            if skip_comments and cols[0].startswith('#'):
                continue
            if skip_pred is not None and skip_pred(cols):
                continue
            rows.append(cols)
    return (header, rows, delim)


def read_dicts(path, sep='auto', skip_comments=False, skip_pred=None,
               lower_header=False, encoding='utf-8-sig', default=None):
    """读表格 → list[dict]（缺列补 ``''``；文件缺失/空返回 default）。

    default 缺省为 None；老代码语义是"缺文件给空列表"的调用方传 ``[]``。
    lower_header=True 时表头小写化（phylodyn_kit 口径）。
    """
    header, rows, _ = read_rows(path, sep=sep, skip_comments=skip_comments,
                                skip_pred=skip_pred, encoding=encoding)
    if header is None:
        return default
    keys = [(h or '').strip().lower() if lower_header else (h or '').strip()
            for h in header]
    out = []
    for cols in rows:
        rec = {}
        for i, k in enumerate(keys):
            if not k:
                continue
            rec[k] = cols[i] if i < len(cols) else ''
        out.append(rec)
    return out


def fnum(v, default=None):
    """脏值容错转 float：``-``/NA/空 → default；容忍 ``12.3%`` 与千分位逗号。"""
    if v is None:
        return default
    if isinstance(v, (int, float)):
        return float(v)
    s = str(v).strip().rstrip('%').replace(',', '')
    if s in NA_TOKENS:
        return default
    try:
        return float(s)
    except ValueError:
        return default


def safe_plot(name, logger, fn, *args, **kwargs):
    """批量出图的容错包装：fn 失败只 WARN 不中断，返回是否成功。

    用法：``safe_plot('深度曲线', logger, plot_depth, sample, out_png)``
    """
    try:
        fn(*args, **kwargs)
        return True
    except Exception as e:  # noqa: BLE001 —— 任何单图失败都不得拖垮批次
        if logger:
            try:
                logger.log(f"{name} 出图失败（跳过）: {e}", 'WARN')
            except Exception:
                pass
        else:
            print(f"[WARN] {name} 出图失败（跳过）: {e}")
        return False
