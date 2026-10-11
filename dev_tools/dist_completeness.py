# -*- coding: utf-8 -*-
"""清点比对：`package.py` 应拷贝的源码目录 vs 产物里实际存在的文件。

为什么需要：`package.py` 结尾只打印一句总数，**不校验"该拷的是否都拷了"**。
`copytree` 中途失败、`ignore_patterns` 写错、源目录被改名，都会表现为
"产物小了一点"，而 `selfcheck` 只看模块能否导入与工具能否探测到，
**看不见"少了一批示例数据 / 一批模板"** 这类缺失。

本工具做的是**集合差**（不比对内容哈希，快）：
  - 缺件  → 硬失败（退出码 1）
  - 多出  → 仅提示（正常，如 PyInstaller 生成的 `.dist-info`）

⚠️ 维护约定：`package.py` 里新增一处 `shutil.copytree(src, dst)` 时，
要在本文件的 `DEFAULT_PAIRS` 里加一条；否则那一块不受本工具保护。

用法:
  python dev_tools/dist_completeness.py
  python dev_tools/dist_completeness.py --pair 3rd/tools=3rd/tools
  python dev_tools/dist_completeness.py --app dist/VirusPlatform
"""
import argparse
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# (源码相对路径, 产物内相对路径) —— 与 dev_tools/package.py 的拷贝清单对齐
DEFAULT_PAIRS = [
    ('3rd/bin', '3rd/bin'),
    ('3rd/tools', '3rd/tools'),
    # 便携 R + treedater（2026-09-17 随包分发）；package.py 拷贝时排除
    # *.exe / doc / tests，见下方 PAIR_IGNORE 的同口径过滤
    ('3rd/R', '3rd/R'),
    ('3rd/open-virome/frontend/build', '3rd/open-virome/frontend/build'),
    ('webapp', '_internal/webapp'),
    # 非 .py 的包内资源（spec 的 datas 投放，代码本身在 exe 的 PYZ 里）
    ('Virus_Platform_Core/ncbi_submit/templates',
     '_internal/Virus_Platform_Core/ncbi_submit/templates'),
    ('Virus_Platform_Core/ncbi_submit/samples',
     '_internal/Virus_Platform_Core/ncbi_submit/samples'),
]

# 每个比对项在源码侧额外跳过的 文件模式/目录名（与 package.py 拷贝该项时
# 的 ignore_patterns 同口径；不在本表的项不受影响）。
PAIR_IGNORE = {
    '3rd/R': {'*.exe', 'doc', 'tests'},
}

# 开发期手工备份（webapp 下 20+ 个 *.bak_*）从不进分发版：spec._tree 与
# 本工具的源码侧遍历都跳过，避免"缺件"误报。
BAK_PATTERNS = ('*.bak', '*.bak_*')

IGNORE = {'.git', '.github', '__pycache__', '.mimosa'}


def rels(base, extra_ignore=()):
    import fnmatch
    out = set()
    for cur, dns, files in os.walk(base):
        dns[:] = [d for d in dns if d not in IGNORE
                  and d not in extra_ignore]
        for fn in files:
            if any(fnmatch.fnmatch(fn, p) for p in BAK_PATTERNS):
                continue
            # extra_ignore 同时承载 文件模式（如 3rd/R 的 *.exe）与目录名
            if any(fnmatch.fnmatch(fn, p) for p in extra_ignore):
                continue
            out.add(os.path.relpath(os.path.join(cur, fn), base)
                    .replace('\\', '/'))
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description='产物缺件清点')
    ap.add_argument('--app',
                    default=os.path.join(os.path.dirname(ROOT),
                                         '植物病毒分析平台-打包产物',
                                         'VirusPlatform'))
    ap.add_argument('--pair', action='append', default=[],
                    metavar='SRC=DST', help='追加/覆盖比对项（可重复）')
    ap.add_argument('--max-list', type=int, default=12,
                    help='每项最多列出的缺件数')
    a = ap.parse_args(argv)

    try:
        sys.stdout.reconfigure(encoding='utf-8', errors='replace')
    except Exception:
        pass

    pairs = list(DEFAULT_PAIRS)
    for p in a.pair:
        if '=' not in p:
            print(f'✘ --pair 格式应为 SRC=DST，收到 {p!r}')
            return 2
        s, d = p.split('=', 1)
        pairs = [x for x in pairs if x[0] != s] + [(s, d)]

    missing_total = 0
    for src_rel, dst_rel in pairs:
        src = os.path.join(ROOT, *src_rel.split('/'))
        dst = os.path.join(a.app, *dst_rel.split('/'))
        if not os.path.isdir(src):
            print(f'-- {src_rel}: 源码目录不存在，跳过')
            continue
        if not os.path.isdir(dst):
            print(f'!! {dst_rel}: 产物目录不存在 ← 整块缺失')
            missing_total += 1
            continue
        A, B = rels(src, extra_ignore=PAIR_IGNORE.get(src_rel, ())), rels(dst)
        miss = sorted(A - B)
        extra = sorted(B - A)
        status = '✔ 无缺件' if not miss else f'✘ 缺 {len(miss)} 个'
        print(f'-- {src_rel} ({len(A)}) -> {dst_rel} ({len(B)})  {status}')
        if miss:
            for m in miss[:a.max_list]:
                print(f'      - {m}')
            if len(miss) > a.max_list:
                print(f'      ... 另有 {len(miss) - a.max_list} 个')
            missing_total += len(miss)
        if extra:
            print(f'   ? 产物多出 {len(extra)} 个'
                  f'（正常情况：PyInstaller 生成物）: {extra[:3]}')

    print()
    if missing_total:
        print(f'✘ 共 {missing_total} 处缺件')
        return 1
    print('✔ 清点通过：应拷目录全部就位')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
