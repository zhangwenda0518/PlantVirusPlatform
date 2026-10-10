#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""把当前工作目录整理成「程序 / 数据库 / 示例」分离布局。

目标布局（全部在 D: 同一卷，用改名实现，瞬时且不需要额外磁盘空间）：

  D:\\桌面\\植物病毒分析平台\\                ① 程序（软件）
      examples\\                               ② 示例（已在程序根目录，git 纳管不变）
      platform.json                            ← 指向 ③ ④

  D:\\桌面\\植物病毒分析平台-数据库\\          ③ 数据库
      databases\\（含 virusref_db\\）  host-db\\

  D:\\桌面\\植物病毒分析平台-数据\\            ④ 运行产物
      run\\（results / tool_runs / logs / tasks / uploads / submissions /
            meta_search / logan / fastq / downloads）

要点
----
**2026-09-10 布局整理后本脚本已同步**：运行期数据已收进 run/、外部依赖收进
3rd/、examples/ 已独立到平台根目录（不再位于 databases/ 内）。重跑前务必先
`--dry-run` 核对当前布局。

* examples/ 不在 databases/ 内了 → 「先移出 databases 再回迁示例」那一步已删除，
  示例天然留在程序目录；Virus_Platform_Core/config.py `_detect_examples_root` 第 1 候选
  (<平台根>/examples) 直接命中，示例零配置。
* 输出根必须整体切换：Virus_Platform_Core/config.py 的 OUTPUT_SUBDIRS 有 7 个目录，
  只搬 results/tool_runs 会让其余几个静默落到新根，造成数据分裂。
* 配置写入走平台自己的 Config.set_database_root / set_output_root
  （与 db_migrate.migrate 结尾同一调用），不手改 platform.json。
* 不用 db_migrate.migrate()：它不支持排除 examples，且要求目标盘
  剩余 ≥ 源+10GB；同卷改名不需要。

用法:
    python dev_tools/_reorg_split.py --dry-run     # 只看计划
    python dev_tools/_reorg_split.py               # 执行
"""
import argparse
import os
import shutil
import sys
# 控制台编码兜底：Windows 默认代码页是 GBK，本脚本的 ✔/✘/⚠ 等字符会让
# print 抛 UnicodeEncodeError（2026-09-11 实测多处踩过）。只改错误处理为
# replace（编码不动，中文照常可读），编不出的字符降级为 '?'。
for _stream in (sys.stdout, sys.stderr):
    try:
        _stream.reconfigure(errors='replace')
    except (AttributeError, OSError):
        pass


ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT)

DB_ROOT = r'D:\桌面\植物病毒分析平台-数据库'
OUT_ROOT = r'D:\桌面\植物病毒分析平台-数据'

# ③ 数据库：整目录迁移
DB_DIRS = ('databases', 'host-db')
# ④ 运行产物：2026-09-10 起已统一收在 run/ 下，整目录迁移即可
OUT_DIRS = ('run',)
# 迁移后要留在程序目录内的项（示例资产，git 纳管）。
# 2026-09-10 起 examples/ 已在平台根目录、不再位于 databases/ 内，
# 「先移出 databases 再回迁 examples」那一步已不再需要。
KEEP_IN_PLACE = ('examples',)


def stat_dir(path):
    """(文件数, 总字节)"""
    n = total = 0
    for cur, _d, files in os.walk(path):
        for fn in files:
            try:
                total += os.path.getsize(os.path.join(cur, fn))
            except OSError:
                continue
            n += 1
    return n, total


def fmt(n):
    for u in ('B', 'KB', 'MB', 'GB', 'TB'):
        if n < 1024 or u == 'TB':
            return f'{n:.2f} {u}' if u != 'B' else f'{n:.0f} B'
        n /= 1024


def plan(dry):
    print('=' * 74)
    print('  程序 / 数据库 / 示例 三分离整理' + ('  [DRY-RUN]' if dry else ''))
    print('=' * 74)
    print(f'  程序目录 : {ROOT}')
    print(f'  数据库   : {DB_ROOT}')
    print(f'  运行产物 : {OUT_ROOT}')
    print()

    steps = []          # (标签, 源, 目标)

    for d in DB_DIRS:
        src = os.path.join(ROOT, d)
        if os.path.isdir(src):
            steps.append((f'数据库/{d}', src, os.path.join(DB_ROOT, d)))
        else:
            print(f'  - 跳过（不存在）: {d}')

    for d in OUT_DIRS:
        src = os.path.join(ROOT, d)
        if os.path.isdir(src):
            steps.append((f'产物/{d}', src, os.path.join(OUT_ROOT, d)))

    tot = 0
    for label, src, dst in steps:
        n, b = stat_dir(src)
        tot += b
        print(f'  {label:22} {n:6d} 文件  {fmt(b):>10}   ->  {dst}')
    print(f'\n  合计迁移 {fmt(tot)}')

    # 示例：2026-09-10 起已独立在平台根目录（<ROOT>/examples），
    # 不在 databases/ 内，故无需回迁步骤。
    ex = os.path.join(ROOT, 'examples')
    if os.path.isdir(ex):
        n, b = stat_dir(ex)
        print(f'\n  示例保留在程序目录（不在 databases/ 内，无需回迁）: {ex}')
        print(f'    {n} 文件 / {fmt(b)}，git 纳管不变')
    return steps


def snapshot():
    snap = {}
    for d in list(DB_DIRS) + list(OUT_DIRS):
        p = os.path.join(ROOT, d)
        if os.path.isdir(p):
            snap[d] = stat_dir(p)
    return snap


def move(src, dst, label):
    if os.path.exists(dst):
        if os.listdir(dst) if os.path.isdir(dst) else True:
            sys.exit(f'✘ 目标已存在且非空，中止: {dst}')
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    shutil.move(src, dst)
    print(f'  ✔ {label}: {os.path.basename(src)} -> {dst}')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--dry-run', action='store_true')
    a = ap.parse_args()

    if (DB_ROOT == ROOT or DB_ROOT.startswith(ROOT + os.sep)
            or OUT_ROOT.startswith(ROOT + os.sep)):
        sys.exit('✘ 目标不能位于平台目录内')

    plan(a.dry_run)          # 只需副作用：打印/记录迁移计划
    print('\n  配置将写入:')
    print(f'    database_root = {DB_ROOT}')
    print(f'    output_root   = {OUT_ROOT}')
    print('    （走 Config.set_database_root / set_output_root）')
    if a.dry_run:
        print('\n  [DRY-RUN] 未做任何改动。')
        return

    before = snapshot()
    print('\n== 1/4 迁移数据库目录 ==')
    for d in DB_DIRS:
        src = os.path.join(ROOT, d)
        if not os.path.isdir(src):
            continue
        dst = os.path.join(DB_ROOT, d)
        move(src, dst, '数据库')

    print('\n== 2/4 示例目录 ==')
    print('  跳过：examples/ 已在平台根目录（不在 databases/ 内）')

    print('\n== 3/4 迁移运行产物 ==')
    for d in OUT_DIRS:
        src = os.path.join(ROOT, d)
        if os.path.isdir(src):
            move(src, os.path.join(OUT_ROOT, d), '产物')

    # ---- 校验：逐目录比对文件数与字节 ----
    print('\n== 校验迁移完整性 ==')
    bad = []
    for d, (n0, b0) in before.items():
        dst = os.path.join(DB_ROOT if d in DB_DIRS else OUT_ROOT, d)
        n1, b1 = stat_dir(dst)
        ok = (n0 == n1 and b0 == b1)
        print(f'  {"✔" if ok else "✘"} {d:14} {n0}文件/{fmt(b0)}'
              f'  ->  {n1}文件/{fmt(b1)}')
        if not ok:
            bad.append(d)
    if bad:
        sys.exit(f'✘ 校验不一致: {bad} —— 配置未改，可手工移回')

    # ---- 写配置（平台自己的 API）----
    print('\n== 4/4 切换配置 ==')
    from Virus_Platform_Core.config import get_config
    cfg = get_config()
    cfg.set_database_root(DB_ROOT)
    print(f'  ✔ database_root = {cfg.database_root}')
    cfg.set_output_root(OUT_ROOT)
    print(f'  ✔ output_root   = {cfg.output_root}')

    # 清理 platform.json 里指向旧位置的 databases 绝对路径
    print('\n== 清理 platform.json 中的陈旧绝对路径 ==')
    from Virus_Platform_Core.config import (get_database_root, get_examples_root, DIRS, db_path,
                                            host_db_dirs)
    import json
    with open('platform.json', encoding='utf-8') as f:
        raw = json.load(f)
    dbs = raw.get('databases') or {}
    changed = False
    for k, v in list(dbs.items()):
        if v and not os.path.isdir(v):
            print(f'  - 移除失效项 databases.{k} = {v}')
            dbs[k] = ''
            changed = True
    if changed:
        raw['databases'] = dbs
        with open('platform.json', 'w', encoding='utf-8', newline='\n') as f:
            json.dump(raw, f, ensure_ascii=False, indent=2)
        print('  platform.json 已更新')
    else:
        print('  无需清理')

    print('\n== 结果 ==')
    print(f'  database_root    = {get_database_root()}')
    print(f'  virus 库         = {db_path("virus", "plant")}')
    # 宿主库按物种分目录；db_path 会解析当前生效的那个（active_host_db）
    print(f'  host 库（当前）   = {db_path("host", "classify")}')
    print(f'  host 库（全部）   = {", ".join(host_db_dirs()) or "（无）"}')
    print(f'  taxonomy         = {db_path("tax", "core")}')
    print(f'  examples         = {get_examples_root()}')
    print(f'  results          = {DIRS["results"]}')
    print(f'  tool_runs        = {DIRS["tool_runs"]}')
    print('\n完成。请接着跑: python main.py selfcheck')


if __name__ == '__main__':
    main()
