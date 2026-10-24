# -*- coding: utf-8 -*-
"""生成 THIRD_PARTY_NOTICES.md：随平台分发的第三方组件清单。

为什么需要（2026-10-24）：平台随包分发 30+ 个第三方可执行工具（3rd/bin、
3rd/tools，其中含 GPL 系：RDP5 / SPAdes / snpEff…）与 130+ 个 Python 依赖
（3rd/python），公开发布必须给接收方一份「用了谁的什么东西、许可证是什么、
去哪看原文」的索引。

口径：本文件是**索引**而非法律意见——每个组件的权威许可证以其目录内随附的
LICENSE/COPYING 与上游声明为准；Python 包以 dist-info METADATA 的 License
字段为准（截断展示）。每次 3rd/ 或 requirements 变更后重跑一次：

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
    return (lic[:97] + '…') if len(lic) > 100 else (lic or '（元数据未声明）')


def _tools_rows():
    rows = []
    # 3rd/bin：单文件工具
    for fn in sorted(glob.glob(os.path.join(ROOT, '3rd', 'bin', '*'))):
        if os.path.isfile(fn):
            rows.append((os.path.relpath(fn, ROOT).replace(os.sep, '/'),
                         ['3rd/bin/（单文件工具，许可见上游 Release 说明）']))
    # 3rd/tools：工具套件目录
    tools_dir = os.path.join(ROOT, '3rd', 'tools')
    for name in sorted(os.listdir(tools_dir)):
        d = os.path.join(tools_dir, name)
        if not os.path.isdir(d):
            continue
        lics = _license_files(d)
        rows.append((f'3rd/tools/{name}/', lics or ['（目录内未随附许可'
                                                    '文件，见上游项目）']))
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
        '# 第三方组件清单（THIRD PARTY NOTICES）',
        '',
        '本文件是**索引**而非法律意见：各组件的权威许可证以其目录内随附的',
        'LICENSE/COPYING 与上游声明为准；Python 包以安装元数据（dist-info',
        'METADATA）的 License 字段为准。重新生成：',
        '`python dev_tools/third_party_notices.py`。',
        '',
        f'## 一、随包分发的第三方工具（{len(t_rows)} 项，3rd/bin + 3rd/tools）',
        '',
        '| 组件（随包路径） | 目录内随附的许可/说明文件 |',
        '|---|---|',
    ]
    for path, lics in t_rows:
        lines.append(f'| {path} | {"; ".join(lics)} |')
    lines += [
        '',
        '例外：`3rd/tools/rdp5/3seqTable`（105MB，超 GitHub 单文件限制）经',
        'Release（tag: external-tools）分发；RDP5 会话运行态文件'
        '（LastSave.rdp5 / tempfile）',
        '不随版本库与插件包分发。',
        '',
        f'## 二、内置 Python 运行时的第三方包（{len(py_rows)} 项，3rd/python）',
        '',
        '完整版本快照见 `requirements.lock`；下表为包名 / 版本 / 元数据声明的许可证。',
        '',
        '| 包 | 版本 | 许可证（元数据） |',
        '|---|---|---|',
    ]
    for name, ver, lic in py_rows:
        lines.append(f'| {name} | {ver} | {lic} |')
    lines += [
        '',
        '## 三、分发口径备注',
        '',
        '- 平台自身代码随仓库分发；第三方工具本体**克隆即用**（3rd/ 全量纳管，',
        '  例外见上），或经「设置 → 插件管理」从插件源按需安装（32 个工具包，',
        '  逐包 sha256 校验）。',
        '- GPL 系组件（RDP5、SPAdes、snpEff、snipit、svguest 等）以**独立进程/'
        '  独立文件**形式随包，不与其源码链接；其许可证文本随各自目录分发。',
        '- 本清单由构建脚本扫描生成，**不构成法律意见**；商用再分发前请按',
        '  各组件许可证自行核对（尤其 GPL 系的传染条款与商标条款）。',
        '',
    ]
    return '\n'.join(lines)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--check', action='store_true',
                    help='只校验现有文件是否为最新（不写入）')
    a = ap.parse_args()
    md = build_md()
    if a.check:
        try:
            with open(OUT, encoding='utf-8') as f:
                cur = f.read()
        except OSError:
            cur = ''
        if cur != md:
            print('THIRD_PARTY_NOTICES.md 已过期：重跑 '
                  'python dev_tools/third_party_notices.py')
            return 1
        print('THIRD_PARTY_NOTICES.md 最新')
        return 0
    with open(OUT, 'w', encoding='utf-8') as f:
        f.write(md)
    n_t = len(_tools_rows())
    n_py = len(_python_rows())
    print(f'已写 {OUT}（工具 {n_t} 项 / Python 包 {n_py} 项）')
    return 0


if __name__ == '__main__':
    sys.exit(main())
