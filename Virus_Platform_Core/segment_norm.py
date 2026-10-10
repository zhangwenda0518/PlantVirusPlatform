# -*- coding: utf-8 -*-
"""Segment（片段）名称规范化——把全库 299 种原始写法收敛到受控词表。

全库 Plant_Virus_Full.Info.tsv 的 Segment 列是录入自由的自由文本：
  DNA-A / DNA A / dna-a / DNA−A        → DNA-A
  RNA2 / RNA 2 / rna_2 / RNA-2         → RNA2
  Segment 1 / Segment_1 / Seg1 / seg 3 → Segment 1 / Segment 3
  S / M / L / S segment / M-RNA        → S / M / L（布尼亚病毒目口径）
  DNA / RNA / Genomic RNA（无分段号）   → 原样保留（单片段病毒的标注）
  其余（注释性文本）                    → 原样保留（长尾，不强行归类）

规范化只用于展示与筛选（Segment_std 派生列）；原始 Segment 列不动。
"""

import re

# (DNA|RNA) + 可选分隔 + 单个字母（Geminiviridae 等双组分为字母）
_RX_DNA_RNA_LETTER = re.compile(r'\b(DNA|RNA)\s*[-_ ]?\s*([A-Z])\b', re.I)
# (DNA|RNA) + 可选分隔 + 数字（多分节 RNA 病毒）
_RX_DNA_RNA_DIGIT = re.compile(r'\b(DNA|RNA)\s*[-_ ]?\s*([1-9])\b', re.I)
# Segment/Seg + 数字（允许两位数：呼肠孤病毒等有 10-12 个片段，
# 原 `([1-9])\b` 只认单数字，Segment 10/11/12 全部漏过规范化）
_RX_SEG_N = re.compile(r'\bSegment\s*[-_# ]?\s*(\d{1,2})\b', re.I)
_RX_SEG_AB = re.compile(r'\bSeg\s*[-_ ]?(\d{1,2})\b', re.I)
# 布尼亚病毒目 L/M/S 三段式
_RX_LMS = re.compile(r'^(?:segment\s*|seg\s*|strand\s*)?([LMS])(?:\s*(?:segment|rna))?$',
                     re.I)

_EMPTY_MARKERS = {'', '-', '--', 'n/a', 'na', 'none', 'unknown', 'null',
                  '未标注', '无'}


def norm_segment(raw):
    """原始 Segment 文本 → 规范化标签（'' 表示未标注/单片段）。"""
    s = (raw or '').strip()
    if not s or s.lower() in _EMPTY_MARKERS:
        return ''
    up = s.upper()

    m = _RX_DNA_RNA_LETTER.search(up)
    if m:
        return f'{m.group(1).upper()}-{m.group(2).upper()}'
    m = _RX_DNA_RNA_DIGIT.search(up)
    if m:
        return f'{m.group(1).upper()}{m.group(2)}'
    m = _RX_SEG_N.search(up) or _RX_SEG_AB.search(up)
    if m:
        return f'Segment {m.group(1)}'
    m = _RX_LMS.match(up)
    if m:
        return m.group(1).upper()
    return s


def norm_segment_ctx(raw, family='', molecule=''):
    """带上下文的规范化：利用科名与分子类型处理裸写法。

    在 norm_segment 的通用规则之后追加：
      - 裸字母 A/B + Gemininiviridae          → DNA-A / DNA-B
      - 裸字母 R/C/M/N/U…+ Nanoviridae        → DNA-R / DNA-C / …
      - 裸数字 + 分子类型含 RNA               → RNA1…RNA9
      - 裸数字 + 分子类型含 DNA               → DNA1…DNA9
    其余保持 norm_segment 的结果（长尾原样保留）。
    """
    up = (raw or '').strip().upper()
    fam = (family or '').upper()
    mol = (molecule or '').upper()

    # 裸字母（Geminiviridae A/B）与裸 U1-U4（Nanoviridae）：科可判定归属。
    # 必须是**单字符**完整值：子串判断会把 "AB"/"SM"/"ND" 这类多字符值
    # （如 "S segment" 之外的自由文本）误当成字母片段加前缀。
    if ((len(up) == 1 and up in 'ABCRSMDNU')
            or re.fullmatch(r'U[1-4]', up)):
        if 'NANO' in fam or 'BABU' in fam:
            return 'DNA-' + up
        if 'GEMINI' in fam and up in ('A', 'B'):
            return 'DNA-' + up
    # 裸数字：RNA 病毒（呼肠孤/双分/伴生等）的分节写法
    if len(up) == 1 and up in '123456789' and 'RNA' in mol and 'DNA' not in mol:
        return 'RNA' + up
    return norm_segment(raw)
