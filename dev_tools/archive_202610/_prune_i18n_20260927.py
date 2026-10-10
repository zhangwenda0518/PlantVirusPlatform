# -*- coding: utf-8 -*-
"""一次性脚本：从 i18n.js 的 zh/en 两块删除确认零引用的孤儿键（2026-09-27）。"""
import json
import re
import sys

for _s in (sys.stdout, sys.stderr):
    try:
        _s.reconfigure(errors='replace')
    except (AttributeError, OSError):
        pass

rem = json.load(open('run/_tmp_i18n_removable.json'))
p = 'webapp/static/i18n.js'
txt = open(p, encoding='utf-8').read()

# 单引号 JS 字符串（允许 \' 等转义）
STR = r"(?:'[^'\\]*(?:\\.[^'\\]*)*')"
removed, missed = 0, []
for k in rem:
    # 首选：'key': 'value', （连同尾逗号整体删除，前一项逗号保留 → 合法）
    pat = re.compile(r"'" + re.escape(k) + r"'\s*:\s*" + STR + r"\s*,\s*")
    new, n = pat.subn('', txt, count=2)          # zh / en 各 1 次
    if n == 0:
        # 兜底：行尾对（无尾逗号）→ 连同前导逗号删除
        pat2 = re.compile(r",\s*'" + re.escape(k) + r"'\s*:\s*" + STR)
        new, n = pat2.subn('', txt, count=2)
    if n == 0:
        missed.append(k)
        continue
    removed += n
    txt = new

open(p, 'w', encoding='utf-8', newline='\n').write(txt)
print(f'删除 {removed} 处（预期 {len(rem) * 2}），未命中 {len(missed)}: {missed[:8]}')
leftover = [k for k in rem if ("'" + k + "':") in txt]
print('残留键:', leftover or '无')
