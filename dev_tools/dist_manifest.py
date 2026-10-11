# -*- coding: utf-8 -*-
"""生成 / 比对分发包的「构建清单」——让"文件数对不上"这类问题一眼可查。

背景：`package.py` 只在结尾打印一句「4743 文件 / 1.815 GB」。下次再打包时
如果这个数变了，**没有任何办法判断是正常的（新增依赖 / 目录调整）还是真缺件**。
本次（2026-09-15）就卡在这里：历史基线记的是 5555 文件 / 1.65 GB，
而实测 4743 / 1.815 GB，差 812 个文件，靠人肉翻 git log 也复原不出旧产物。

本工具把产物清点成一份 JSON 快照（按顶层目录汇总 + 全量文件清单），
于是"变了多少、变了哪些"变成一次 `--diff`。

用法:
  # 生成快照（默认写 run/dist_manifest_<YYYYmmdd-HHMMSS>.json）
  python dev_tools/dist_manifest.py

  # 指定输出
  python dev_tools/dist_manifest.py --out run/baseline.json

  # 比对两份快照
  python dev_tools/dist_manifest.py --diff run/baseline.json run/dist_manifest_xxx.json

  # 比对时只看汇总（不列具体文件）
  python dev_tools/dist_manifest.py --diff a.json b.json --summary

注意：`bytes` 是 **apparent size**（`os.path.getsize` 直接求和），
与 `package.py` 的口径一致。硬链接会被**重复计数**，所以它**不是磁盘占用**；
真实占用请看 `du`。这样设计是为了与历史记录可比。
"""
import argparse
import datetime
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SKIP_DIRS = {'__pycache__'}


def _fmt(n):
    for unit in ('B', 'KB', 'MB', 'GB', 'TB'):
        if abs(n) < 1024 or unit == 'TB':
            return f'{n:.1f} {unit}' if unit != 'B' else f'{n:.0f} B'
        n /= 1024


def scan(app):
    """返回 {rel_path: size} —— 全量文件清单。"""
    out = {}
    for cur, dns, files in os.walk(app):
        dns[:] = [d for d in dns if d not in SKIP_DIRS]
        for fn in files:
            p = os.path.join(cur, fn)
            try:
                out[os.path.relpath(p, app).replace('\\', '/')] = \
                    os.path.getsize(p)
            except OSError:
                pass
    return out


def top_of(rel):
    """文件相对路径 -> 顶层分组名（用于汇总表）。"""
    parts = rel.split('/')
    if len(parts) == 1:
        return '(根级文件)'
    return parts[0] + '/'


def build(app):
    files = scan(app)
    groups = {}
    for rel, sz in files.items():
        g = groups.setdefault(top_of(rel), [0, 0])
        g[0] += 1
        g[1] += sz
    return {
        'app': app,
        'generated': datetime.datetime.now().isoformat(timespec='seconds'),
        'n_files': len(files),
        'bytes': sum(files.values()),
        'groups': {k: {'n': v[0], 'bytes': v[1]}
                   for k, v in sorted(groups.items())},
        'files': dict(sorted(files.items())),
    }


def show(m):
    print(f'  文件数 {m["n_files"]}  apparent {_fmt(m["bytes"])}')
    for g, v in m['groups'].items():
        print(f'    {g:<24}{v["n"]:>7}  {_fmt(v["bytes"]):>10}')


def diff(a, b, summary_only):
    print(f'A(旧): {a["app"]}  {a["generated"]}')
    show(a)
    print(f'\nB(新): {b["app"]}  {b["generated"]}')
    show(b)

    fa, fb = a['files'], b['files']
    added = sorted(set(fb) - set(fa))
    removed = sorted(set(fa) - set(fb))
    changed = sorted(p for p in set(fa) & set(fb) if fa[p] != fb[p])

    print('\n== 差异 ==')
    print(f'  文件数 {a["n_files"]} -> {b["n_files"]} '
          f'({b["n_files"] - a["n_files"]:+d})')
    print(f'  体积   {_fmt(a["bytes"])} -> {_fmt(b["bytes"])} '
          f'({_fmt(b["bytes"] - a["bytes"])})')
    print(f'  新增 {len(added)} / 删除 {len(removed)} / 内容变化 {len(changed)}')

    # 按顶层分组统计增减
    ga, gb = a['groups'], b['groups']
    keys = sorted(set(ga) | set(gb))
    rows = []
    for k in keys:
        n0 = ga.get(k, {}).get('n', 0)
        n1 = gb.get(k, {}).get('n', 0)
        b0 = ga.get(k, {}).get('bytes', 0)
        b1 = gb.get(k, {}).get('bytes', 0)
        if (n0, n1, b0, b1) != (n0, n0, b0, b0):
            rows.append((k, n1 - n0, b1 - b0))
    if rows:
        print('\n  分组增减:')
        for k, dn, db in sorted(rows, key=lambda r: -abs(r[2])):
            print(f'    {k:<24}{dn:>+8}  {_fmt(db):>+10}')

    if not summary_only:
        for label, items in (('新增', added), ('删除', removed),
                             ('内容变化', changed)):
            if items:
                print(f'\n  {label}（{len(items)}）:')
                for p in items[:30]:
                    extra = ''
                    if label == '内容变化':
                        extra = f'  {_fmt(fa[p])} -> {_fmt(fb[p])}'
                    print(f'    {p}{extra}')
                if len(items) > 30:
                    print(f'    ... 另有 {len(items) - 30} 个')
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(description='分发包构建清单：生成 / 比对')
    ap.add_argument('--app',
                    default=os.path.join(os.path.dirname(ROOT),
                                         '植物病毒分析平台-打包产物',
                                         'VirusPlatform'),
                    help='产物目录（默认 仓库上一级/植物病毒分析平台-打包产物/VirusPlatform）')
    ap.add_argument('--out', help='快照输出路径')
    ap.add_argument('--diff', nargs=2, metavar=('A', 'B'), help='比对两份快照')
    ap.add_argument('--summary', action='store_true', help='只输出汇总')
    a = ap.parse_args(argv)

    try:
        sys.stdout.reconfigure(encoding='utf-8', errors='replace')
    except Exception:
        pass

    if a.diff:
        with open(a.diff[0], encoding='utf-8') as f:
            A = json.load(f)
        with open(a.diff[1], encoding='utf-8') as f:
            B = json.load(f)
        return diff(A, B, a.summary)

    if not os.path.isdir(a.app):
        print(f'✘ 产物目录不存在: {a.app}')
        return 1
    m = build(a.app)
    out = a.out or os.path.join(
        ROOT, 'run',
        'dist_manifest_' + datetime.datetime.now().strftime('%Y%m%d-%H%M%S')
        + '.json')
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, 'w', encoding='utf-8') as f:
        json.dump(m, f, ensure_ascii=False, indent=1)
    print(f'快照已写入: {out}')
    show(m)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
