# -*- coding: utf-8 -*-
"""生成 THIRD_PARTY_NOTICES.md：随平台分发的第三方组件清单。

为什么需要（2026-10-24）：平台随包分发 30+ 个第三方可执行工具（3rd/bin、
3rd/tools，其中含 GPL 系：RDP5 / SPAdes / snpEff…）与 130+ 个 Python 依赖
（3rd/python），公开发布必须给接收方一份「用了谁的什么东西、许可证是什么、
去哪看原文」的索引。

口径：本文件是**索引**而非法律意见——每个组件的权威许可证以其目录内随附的
LICENSE/COPYING 与上游声明为准；Python 包以 dist-info METADATA 的 License
字段为准（截断展示）。每次 3rd/ 或 requirements 变更后重跑一次。

2026-10-26：输出正文改为英文（面向英文分发），本模块注释仍保留中文。

用法:
    python dev_tools/third_party_notices.py           # 写 THIRD_PARTY_NOTICES.md
    python dev_tools/third_party_notices.py --check   # 只校验是否最新（旧则退出码 1）
"""
import argparse
import glob
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, 'THIRD_PARTY_NOTICES.md')

_LIC_PAT = re.compile(r'^(license|copying|notice|copyright|authors?|readme)'
                      r'[^/\\]*$', re.I)


def _license_files(tool_dir):
    """工具目录两层内随附的许可/说明文件（相对平台根的路径列表）。"""
    out = []
    for cur, dns, files in os.walk(tool_dir):
        depth = os.path.relpath(cur, tool_dir).count(os.sep)
        if depth > 1:
            dns[:] = []
            continue
        for fn in files:
            if _LIC_PAT.match(fn):
                out.append(os.path.relpath(os.path.join(cur, fn), ROOT)
                           .replace(os.sep, '/'))
    return sorted(out)[:4]


def _fmt_lic(meta_lic):
    """METADATA 的 License 字段常带整段许可证文本，取首行截断。"""
    lic = (meta_lic or '').strip().splitlines()
    lic = lic[0].strip() if lic else ''
    return (lic[:97] + '...') if len(lic) > 100 else (
        lic or '(not declared in metadata)')


def _tools_rows():
    rows = []
    # 3rd/bin：单文件工具
    for fn in sorted(glob.glob(os.path.join(ROOT, '3rd', 'bin', '*'))):
        if os.path.isfile(fn):
            rows.append((os.path.relpath(fn, ROOT).replace(os.sep, '/'),
                         ['3rd/bin/ (single-file tool; license in the '
                          'upstream Release notes)']))
    # 3rd/tools：工具套件目录
    tools_dir = os.path.join(ROOT, '3rd', 'tools')
    for name in sorted(os.listdir(tools_dir)):
        d = os.path.join(tools_dir, name)
        if not os.path.isdir(d):
            continue
        lics = _license_files(d)
        rows.append((f'3rd/tools/{name}/', lics or [
            '(no license file bundled in this directory; see the upstream '
            'project)']))
    return rows


def _python_rows():
    rows = []
    sp = os.path.join(ROOT, '3rd', 'python', 'Lib', 'site-packages')
    for di in sorted(glob.glob(os.path.join(sp, '*.dist-info'))):
        meta = os.path.join(di, 'METADATA')
        if not os.path.isfile(meta):
            continue
        fields = {}
        with open(meta, encoding='utf-8', errors='replace') as f:
            for line in f:
                if not line.strip():
                    break                      # 头部结束
                m = re.match(r'^([\w-]+):\s*(.*)$', line)
                if m:
                    fields[m.group(1)] = m.group(2)
        name = fields.get('Name') or os.path.basename(di)[:-len('.dist-info')]
        rows.append((name, fields.get('Version', '?'),
                     _fmt_lic(fields.get('License')
                              or fields.get('License-Expression'))))
    return rows


def build_md():
    t_rows = _tools_rows()
    py_rows = _python_rows()
    lines = [
        '# Third-Party Component Notices (THIRD PARTY NOTICES)',
        '',
        'This file is an **index**, not legal advice: the authoritative license',
        'for each component is the LICENSE/COPYING file bundled in its',
        'directory together with the upstream notice; for Python packages the',
        'License field of the installed metadata (dist-info METADATA) is',
        'authoritative (shown truncated). Regenerate with',
        '`python dev_tools/third_party_notices.py`.',
        '',
        f'## 1. Third-party tools shipped with the package '
        f'({len(t_rows)} items, 3rd/bin + 3rd/tools)',
        '',
        '| Component (path in package) | License/notice files bundled in the '
        'component directory |',
        '|---|---|',
    ]
    for path, lics in t_rows:
        lines.append(f'| {path} | {"; ".join(lics)} |')
    lines += [
        '',
        'Exception: `3rd/tools/rdp5/3seqTable` (105 MB, above the GitHub',
        'single-file limit) is distributed through Releases (tag:',
        'external-tools); RDP5 session runtime files (LastSave.rdp5 /',
        'tempfile) are not shipped in the repository or in plugin packages.',
        '',
        f'## 2. Third-party packages in the bundled Python runtime '
        f'({len(py_rows)} items, 3rd/python)',
        '',
        'The full version snapshot is in `requirements.lock`; the table below',
        'lists package name / version / license as declared in the metadata.',
        '',
        '| Package | Version | License (metadata) |',
        '|---|---|---|',
    ]
    for name, ver, lic in py_rows:
        lines.append(f'| {name} | {ver} | {lic} |')
    lines += [
        '',
        '## 3. Distribution notes',
        '',
        '- The platform code itself ships with the repository; the third-party',
        '  tool bodies are usable right after cloning (`3rd/` is tracked in',
        '  full, exceptions above), or can be installed on demand from the',
        '  plugin source in Settings -> Plugins (32 tool packages, per-package',
        '  sha256 verification).',
        '- GPL-family components (RDP5, SPAdes, snpEff, snipit, svguest and',
        '  others) ship as **separate processes / separate files**, not linked',
        '  against their sources; their license texts ship in their own',
        '  directories.',
        '- This list is generated by the build script and **does not constitute',
        '  legal advice**; before commercial redistribution, verify each',
        "  component's license yourself (especially the GPL copyleft and",
        '  trademark terms).',
        '',
    ]
    return '\n'.join(lines)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--check', action='store_true',
                    help='only validate that the existing file is up to date '
                         '(do not write)')
    a = ap.parse_args()
    md = build_md()
    if a.check:
        try:
            with open(OUT, encoding='utf-8') as f:
                cur = f.read()
        except OSError:
            cur = ''
        if cur != md:
            print('THIRD_PARTY_NOTICES.md is stale: re-run '
                  'python dev_tools/third_party_notices.py')
            return 1
        print('THIRD_PARTY_NOTICES.md is up to date')
        return 0
    with open(OUT, 'w', encoding='utf-8') as f:
        f.write(md)
    n_t = len(_tools_rows())
    n_py = len(_python_rows())
    print(f'wrote {OUT} ({n_t} tools / {n_py} Python packages)')
    return 0


if __name__ == '__main__':
    sys.exit(main())
