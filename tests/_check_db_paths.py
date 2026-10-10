# -*- coding: utf-8 -*-
"""数据库/目录类 API 的路径绝对性检查。

背景（历史事故）：`/api/dbs` 曾对内置病毒库返回 `databases/<rel>` —— 相对
**数据库根**的路径。前端原样回传，后端却按 PLATFORM_ROOT 解析，于是分发版
（数据库包与程序目录分离）下选库直接 400「病毒库不存在」。源码树里
`PLATFORM_ROOT/databases` 恰好就是数据库根，所以只在分发版暴露；更阴的是
两边都存在时会**静默用错库**。

因此约定：凡是要回传给前端的「库/目录」路径，必须是绝对路径。

**第二条判据（2026-10-07 审查 P1-2 追加）**：后端也**不得**从
`PLATFORM_ROOT/databases` 读盘 —— 分发版数据库外置时该目录不是数据库根。
判据走 AST（见文末 `_scan_platform_databases`），白名单只有 `config.py`
（定义该根）与 `db_migrate.py`（迁移源）。含灵敏度自测。

用法: python tests/_check_db_paths.py
退出码: 0 = 全部绝对且无 PLATFORM_ROOT/databases 读盘；1 = 存在问题。
"""
import ast
import os
import sys

for _s in (sys.stdout, sys.stderr):
    try:
        _s.reconfigure(errors='replace')
    except (AttributeError, OSError):
        pass

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import app  # noqa: E402

SEP = ('/', '\\')

# 这些字段名装的是**给人看的标签/说明**，含斜杠也正常，不参与路径判定。
# 例：/api/dbs 的 virus_libs[].name = 'kunpeng_db/plant'（下拉显示名，
# option 的 value 走的是同一对象的 path 字段，那里必须是绝对路径）；
# /api/kv_index_list 的 manifest.layout 是一段中文说明；
# manifest.notes / fix_script / original_broken_backup 是数据集的**来源记录**
# （修复脚本名、破损备份位置等，全库 grep 无任何代码消费），同为说明字段。
LABEL_KEYS = ('name', 'label', 'title', 'note', 'notes', 'desc', 'layout',
              'hint', 'text', 'msg', 'message',
              'fix_script', 'original_broken_backup')


def _is_label_field(path):
    leaf = path.rsplit('.', 1)[-1].split('[')[0].lower()
    return leaf in LABEL_KEYS


def walk_paths(obj, path=''):
    """递归收集所有「看起来像路径」的字符串值及其字段路径。"""
    out = []
    if isinstance(obj, dict):
        for k, v in obj.items():
            out += walk_paths(v, f'{path}.{k}')
    elif isinstance(obj, list):
        for i, v in enumerate(obj):
            out += walk_paths(v, f'{path}[{i}]')
    elif isinstance(obj, str):
        if len(obj) > 3 and not obj.startswith('http') and any(s in obj for s in SEP):
            out.append((path, obj))
    return out


# --------------------------------------------------------------------------
# 静态判据（2026-10-07 审查 P1-2）
#
# 背景：`/api/tool/dsrna_panel`、`_mirna_lib_path`、`crispr_design` 的近缘库与
# BLAST 缓存曾用 `PLATFORM_ROOT/databases`。分发版由
# `config._detect_database_root()` 把数据库根重定向到同级
# `VirusPlatform-Database`，两者**不等价** —— dsrna 面板清单恒为空（且代码把
# "目录缺失"当正常，属静默失效）、miRNA 库直接 abort(500)。
#
# 上面那段响应式判据（"返回给前端的路径必须绝对"）**原理上抓不到**这类
# "后端从错误的根读盘"，故补一条 AST 判据。
#
# 判据走 AST 而非文本搜索：修复注释里必然会写"原来是 PLATFORM_ROOT/databases"，
# 文本搜索会把注释当代码（本项目已多次踩此坑）。
# --------------------------------------------------------------------------

# 白名单：这两个模块的职责**就是**操作平台内目录
_DB_ROOT_WHITELIST = {
    os.path.join('Virus_Platform_Core', 'config.py'),      # 定义 DIRS['databases'] 本身
    os.path.join('Virus_Platform_Core', 'db_migrate.py'),  # 迁移的**源**即平台内目录
}
_SCAN_ROOTS = ('Virus_Platform_Core',)
_SCAN_FILES = ('app.py', 'main.py')


def _is_platform_databases_join(node):
    """`os.path.join(PLATFORM_ROOT, 'databases', ...)` → True。

    ⚠️ `os.path.join` 的 AST 是 Attribute(attr='join',
    value=Attribute(attr='path', value=Name('os'))) —— 早期版本只判
    `func.value` 是 `Name('os')`，恒不匹配、守卫静默失灵（被本文件自测①抓住）。
    """
    if not isinstance(node, ast.Call):
        return False
    f = node.func
    if not (isinstance(f, ast.Attribute) and f.attr == 'join'
            and isinstance(f.value, ast.Attribute) and f.value.attr == 'path'
            and isinstance(f.value.value, ast.Name)
            and f.value.value.id == 'os'):
        return False
    a = node.args
    if len(a) < 2:
        return False
    if not (isinstance(a[0], ast.Name) and a[0].id == 'PLATFORM_ROOT'):
        return False
    return isinstance(a[1], ast.Constant) and a[1].value == 'databases'


def _find_in_source(src, fname):
    """在源码字符串里找出所有 PLATFORM_ROOT/databases 读盘点（AST，忽略注释）。"""
    out = []
    try:
        tree = ast.parse(src)
    except SyntaxError:
        return out
    for n in ast.walk(tree):
        if _is_platform_databases_join(n):
            out.append((fname, n.lineno))
    return out


def _scan_platform_databases():
    bad = []
    files = []
    for r in _SCAN_ROOTS:
        for dp, _d, fns in os.walk(os.path.join(ROOT, r)):
            if '__pycache__' in dp:
                continue
            files += [os.path.join(dp, f) for f in fns if f.endswith('.py')]
    files += [os.path.join(ROOT, f) for f in _SCAN_FILES]
    n_hit = 0
    for p in sorted(files):
        rel = os.path.relpath(p, ROOT)
        if rel in _DB_ROOT_WHITELIST:
            continue
        try:
            src = open(p, encoding='utf-8').read()
        except OSError:
            continue
        for fname, lineno in _find_in_source(src, rel):
            print(f'    △ {fname}:{lineno} 用 PLATFORM_ROOT/databases 读盘'
                  f'（分发版会读错目录）')
            bad.append(f'{fname}:{lineno} PLATFORM_ROOT/databases')
            n_hit += 1
    print(f'AST 扫描：{len(files)} 个文件，命中 {n_hit} 处'
          f'（白名单 {len(_DB_ROOT_WHITELIST)} 个模块）')
    return bad


def _selftest():
    """灵敏度自测：注入故意错误必须被抓住，注释/字符串里的同样字样不得误报。"""
    bad = []
    pos = ("import os\n"
           "from x import PLATFORM_ROOT\n"
           "p = os.path.join(PLATFORM_ROOT, 'databases', 'dsrna')\n")
    hit = _find_in_source(pos, '<selftest-pos>')
    print(f'  自测① 代码形式命中 {len(hit)} 处（应 1）')
    if len(hit) != 1:
        bad.append('自测①：注入的代码形式未被抓住 → 守卫失灵')

    neg = ("# 原来是 os.path.join(PLATFORM_ROOT, 'databases')，已改\n"
           'MSG = "os.path.join(PLATFORM_ROOT, \'databases\')"\n')
    hit2 = _find_in_source(neg, '<selftest-neg>')
    print(f'  自测② 注释/字符串形式命中 {len(hit2)} 处（应 0）')
    if len(hit2) != 0:
        bad.append('自测②：注释/字符串被误判为代码 → 假阳性')

    # ③ 判据在**真实代码**上确有牙齿：白名单里的 config.py 必须能被扫出命中
    cfg = os.path.join(ROOT, 'Virus_Platform_Core', 'config.py')
    try:
        real = _find_in_source(open(cfg, encoding='utf-8').read(),
                               'config.py')
    except OSError:
        real = []
    print(f'  自测③ 真实白名单文件 config.py 命中 {len(real)} 处（应 > 0）')
    if len(real) == 0:
        bad.append('自测③：真实代码里扫不出命中 → 判据写错（守卫无牙齿）')
    return bad


def main():
    c = app.app.test_client()
    bad = []
    for ep in ['/api/dbs', '/api/host_dbs', '/api/kv_index_list']:
        r = c.get(ep)
        if r.status_code != 200:
            print(f'{ep}: HTTP {r.status_code}（跳过）')
            continue
        ps = [x for x in walk_paths(r.get_json()) if not _is_label_field(x[0])]
        rel = [(k, v) for k, v in ps if not os.path.isabs(v)]
        print(f'{ep}: 路径字段 {len(ps)} 个，非绝对 {len(rel)} 个')
        for k, v in rel:
            print(f'    △ {k} = {v}')
            bad.append(f'{ep} {k}={v}')

    # ---- 静态判据：后端不得从 PLATFORM_ROOT/databases 读盘 ----
    print()
    bad += _scan_platform_databases()
    bad += _selftest()

    print()
    if bad:
        print(f'✘ 存在 {len(bad)} 个问题（分发版下会解析错库 / 守卫自身失灵）')
        return 1
    print('✔ 库/目录类 API 返回的路径全部为绝对路径；'
          '无 PLATFORM_ROOT/databases 读盘')
    return 0


if __name__ == '__main__':
    sys.exit(main())
