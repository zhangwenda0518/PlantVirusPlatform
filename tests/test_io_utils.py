# -*- coding: utf-8 -*-
"""io_utils 容错读取三件套单测。

语义契约：文件缺失/空 → default 而非异常；BOM/# 注释/空行容忍；
分隔符自动嗅探；fnum 对脏值返回 default；safe_plot 失败不抛。
"""
import os

import pytest

from Virus_Platform_Core.io_utils import fnum, read_dicts, read_rows, safe_plot


def _write(path, text, encoding='utf-8-sig'):
    with open(path, 'w', encoding=encoding, newline='') as f:
        f.write(text)
    return str(path)


def test_missing_file_returns_default(tmp_path):
    nope = str(tmp_path / 'nope.tsv')
    assert read_dicts(nope) is None
    assert read_dicts(nope, default=[]) == []
    assert read_rows(nope) == (None, [], '')


def test_sniff_and_pad(tmp_path):
    p = _write(tmp_path / 't.csv', 'Name,Val\nA,1\nB\n')  # CSV + 短行
    rows = read_dicts(p)
    assert rows == [{'Name': 'A', 'Val': '1'}, {'Name': 'B', 'Val': ''}]


def test_comments_bom_blank_lines(tmp_path):
    # utf-8-sig 编码自带 BOM；文本里不再手工嵌 BOM 字符
    p = _write(tmp_path / 't.tsv', '# comment\nName\tVal\n\nA\t1.5\n',
               encoding='utf-8-sig')
    hdr, rows, delim = read_rows(p, skip_comments=True)
    assert hdr == ['Name', 'Val'] and delim == '\t'
    assert rows == [['A', '1.5']]


def test_skip_pred(tmp_path):
    p = _write(tmp_path / 't.tsv', 'H1\tH2\nInt64\tInt64\nA\t1\n')
    hdr, rows, _ = read_rows(
        p, sep='\t', skip_pred=lambda c: all(x in ('Int64', 'str') for x in c))
    assert hdr == ['H1', 'H2'] and rows == [['A', '1']]


def test_empty_file(tmp_path):
    p = _write(tmp_path / 't.tsv', '')
    assert read_dicts(p, default=[]) == []
    assert read_rows(p) == (None, [], '')


def test_fnum_dirty_values():
    assert fnum('-') is None
    assert fnum('NA') is None
    assert fnum('12.3%') == pytest.approx(12.3)
    assert fnum('1,234') == 1234.0
    assert fnum('abc', -1) == -1
    assert fnum(None, 0) == 0
    assert fnum(7) == 7.0


def test_safe_plot_swallow_and_log():
    class L:
        def __init__(self):
            self.msgs = []

        def log(self, m, lv='INFO'):
            self.msgs.append((lv, m))

    def boom(x):
        raise RuntimeError('炸了')

    lgr = L()
    assert safe_plot('测试图', lgr, boom, 1) is False
    assert lgr.msgs and '炸了' in lgr.msgs[0][1]

    def ok(x):
        return x + 1

    assert safe_plot('正常图', lgr, ok, 41) is True


def test_lower_header_option(tmp_path):
    p = _write(tmp_path / 't.tsv', 'Sample ID\tValue\nS1\t2\n')
    rows = read_dicts(p, sep='\t')
    assert 'Sample ID' in rows[0]
    rows2 = read_dicts(p, sep='\t', lower_header=True)
    assert 'sample id' in rows2[0] and os.path.isfile(p)
