# -*- coding: utf-8 -*-
"""校验 spec 里 webapp 的 datas 展开逻辑：**既不漏文件，也不带工具残留**。

为什么单独测这个：`dev_tools/VirusPlatform.spec` 里把 `webapp/` 打进 `datas`
时，原先直接写 `(目录, 目标)`，PyInstaller 会整棵拷贝 —— 连
`webapp/templates/.mimosa`（编辑器 hook 缓存）一起进分发版。改成逐文件展开
（`_tree`）后，风险反转成"**漏掉某个前端文件 → 界面某个页面挂掉**"。

spec 是 exec 式的，没法直接 import（尾部会触发真正的 Analysis/EXE）。
所以这里用 AST 把 `_tree` **单独抽出来**编译执行 —— 测的是 spec 里那份真代码，
不是复制品。

判据：
  1. webapp 下所有非残留文件都在 datas 里（集合相等）→ 不漏
     （例外：*.bak / *.bak_* 是开发期手工备份，2026-09-18 起 spec 有意
      不打进分发版 —— 期望集合同样剔除，不算"漏"）
  2. 没有任何条目来自 `.mimosa` / `__pycache__` 等残留 → 不多
  3. 目标路径为 posix 风格、且都在 webapp/ 下
  4. 嵌套文件的目录映射正确
  5. 备份文件确实不在 datas 里（正面对照，防止排除逻辑被改坏）
"""
import ast
import fnmatch
import os
import sys
import tempfile

sys.stdout.reconfigure(encoding='utf-8', errors='replace')

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SPEC = os.path.join(ROOT, 'dev_tools', 'VirusPlatform.spec')
EXCLUDE = ('.mimosa', '__pycache__', '.git', '.github')
BAK_PATTERNS = ('*.bak', '*.bak_*')


def _is_bak(name):
    return any(fnmatch.fnmatch(name, p) for p in BAK_PATTERNS)

fails = []


def check(label, cond, detail=''):
    print(f'  {"OK  " if cond else "FAIL"} {label}' + (f'  {detail}' if detail else ''))
    if not cond:
        fails.append(label)


print('spec datas 展开校验：')

# --- 从 spec 里抽出 _tree 并编译执行 ---------------------------------
with open(SPEC, encoding='utf-8') as f:
    src = f.read()
tree = ast.parse(src)
fn = next((n for n in tree.body
           if isinstance(n, ast.FunctionDef) and n.name == '_tree'), None)
check('spec 中能找到 _tree', fn is not None)
if fn is None:
    print('\n✘ 无法继续')
    raise SystemExit(1)

seg = ast.get_source_segment(src, fn)
ns = {'os': os}
exec(compile(seg, SPEC, 'exec'), ns)
_tree = ns['_tree']

webapp = os.path.join(ROOT, 'webapp')
check('webapp 目录存在', os.path.isdir(webapp))

# --- 实际展开 --------------------------------------------------------
datas = _tree(webapp, 'webapp')
got_src = {os.path.relpath(s, ROOT).replace('\\', '/') for s, _d in datas}

# --- 期望集合：webapp 下所有非残留、非备份文件 --------------------------
want_src = set()
bak_src = set()
for cur, dns, files in os.walk(webapp):
    dns[:] = [d for d in dns if d not in EXCLUDE]
    for fnm in files:
        if fnm in EXCLUDE:
            continue
        rel = os.path.relpath(os.path.join(cur, fnm), ROOT).replace('\\', '/')
        if _is_bak(fnm):
            bak_src.add(rel)
            continue
        want_src.add(rel)

missing = sorted(want_src - got_src)
extra = sorted(got_src - want_src)
check(f'不漏文件（{len(want_src)} 个）', not missing,
      f'缺 {missing[:5]}' if missing else '')
check('不带残留', not extra, f'多 {extra[:5]}' if extra else '')

# --- 备份文件确实被排除（正面对照，防止 _tree 的排除逻辑被改坏）---------
bak_leak = sorted(got_src & bak_src)
check('备份文件（*.bak*）不在 datas 里', not bak_leak,
      f'混入 {bak_leak[:3]}' if bak_leak else '')
# 2026-09-27 起源码树不再有 *.bak_*（71 个手工备份移 archive/_bak_20260927，
# 规则见 .gitignore §8）——反向断言：若再出现，说明有人又复制了手工备份，
# 应归档而不是提交；此时「不在 datas 里」的排除逻辑重新变得可证伪。
#
# 2026-10-07 审查 P2 修正：这条断言原先只统计 **webapp/** 下的 bak（`bak_src`
# 来自上面的 `os.walk(webapp)`），标签却写"源码树" —— 覆盖面小于字面印象，
# 于是源码树里 17 个 .bak*（tests/ 12 个、Virus_Platform_Core/ 2 个、
# docs/ 1 个…）长期存在而本条始终绿灯。现改为扫**整棵源码树**，并加
# 「范围自证」+「灵敏度自测」，防止再次变成假绿。
TREE_EXCLUDE = {'archive', '3rd', 'dist', 'run', '.git', '__pycache__',
                'build', '.mimosa', '.pytest_cache', 'node_modules',
                'uploads', 'fastq', 'host-db', 'databases'}
TREE_SRC_DIRS = {'Virus_Platform_Core', 'tests', 'docs', 'webapp', 'dev_tools'}


def _scan_bak(root):
    """扫 root 下全部 *.bak*，返回 (相对 root 的命中列表, 文件数, 顶层目录集)。"""
    hits, n, tops = [], 0, set()
    for cur, dns, files in os.walk(root):
        dns[:] = [d for d in dns if d not in TREE_EXCLUDE]
        rel = os.path.relpath(cur, root)
        tops.add(rel.split(os.sep)[0] if rel != os.curdir else os.curdir)
        for fnm in files:
            n += 1
            if _is_bak(fnm):
                hits.append(os.path.relpath(os.path.join(cur, fnm), root)
                            .replace('\\', '/'))
    return sorted(hits), n, tops


tree_bak, tree_n, tree_tops = _scan_bak(ROOT)
check(f'源码树已无 .bak 手工备份（若再现请归档 archive/_bak_*）'
      f'（扫 {tree_n} 个文件）', not tree_bak,
      f'出现 {len(tree_bak)} 个：{tree_bak[:3]}' if tree_bak else '')
# 范围自证：必须真的走到了 webapp/ 之外的源码目录，否则上面那条是假绿
check('范围自证：扫描覆盖 webapp/ 之外的源码目录',
      TREE_SRC_DIRS <= tree_tops,
      f'未覆盖 {sorted(TREE_SRC_DIRS - tree_tops)}（已覆盖 {sorted(tree_tops)}）')
# 灵敏度自测：临时目录里放一个 .bak，必须被扫到（否则判据本身失灵）
with tempfile.TemporaryDirectory(prefix='vp_bak_probe_') as _td:
    open(os.path.join(_td, 'probe.py.bak_19700101_selftest'),
         'w', encoding='utf-8').close()
    _h, _n, _t = _scan_bak(_td)
check('灵敏度：临时目录里的 .bak 能被扫到', len(_h) == 1,
      f'扫到 {len(_h)} 个（应 1）')

# --- 残留确实被排掉（正面对照，避免"两边都错"）-------------------------
residue = [s for s, _d in datas if '.mimosa' in s or '__pycache__' in s]
check('残留被显式排除', not residue, f'{residue[:3]}' if residue else '')
check('源码树里确实存在 .mimosa（否则本条无意义）',
      any('.mimosa' in p for p in want_src) or
      os.path.exists(os.path.join(webapp, 'templates', '.mimosa')))

# --- 目标路径 --------------------------------------------------------
bad_dest = [d for _s, d in datas if '\\' in d or not d.startswith('webapp')]
check('目标路径为 posix 且都在 webapp/ 下', not bad_dest, f'{bad_dest[:3]}')

# 嵌套文件的目录映射：对**任意**一层嵌套做通用校验（不绑定具体文件名）
nested = [(s, d) for s, d in datas
          if os.path.relpath(s, webapp).replace('\\', '/').count('/') >= 1]
check('存在嵌套文件可供校验', bool(nested))
if nested:
    s0, d0 = nested[0]
    rel_dir = os.path.dirname(os.path.relpath(s0, webapp)).replace('\\', '/')
    check('嵌套文件目录映射正确',
          d0 == 'webapp/' + rel_dir,
          f'{os.path.basename(s0)} -> {d0!r}（期望 webapp/{rel_dir!r}）')

# --- 空目录边界：_tree 逐文件展开，不会创建空目录 -----------------------
# 原来写 (目录, 目标) 时 PyInstaller 会连空目录一起建；改成逐文件后不会。
# webapp 下若出现空目录（例如预留给上传的 static/uploads），产物里就会**静默消失**。
empty_dirs = []
for cur, dns, files in os.walk(webapp):
    dns[:] = [d for d in dns if d not in EXCLUDE]
    if not files and not dns and cur != webapp:
        empty_dirs.append(os.path.relpath(cur, webapp).replace('\\', '/'))
check('webapp 下无空目录（_tree 不创建空目录）', not empty_dirs,
      f'需处理: {empty_dirs}')

print()
if fails:
    print(f'✘ {len(fails)} 项未通过: {fails}')
    raise SystemExit(1)
print(f'✔ 全部通过（datas 覆盖 {len(datas)} 个前端文件，零残留）')
