# -*- coding: utf-8 -*-
"""依赖审计：扫描全部源码的 import，找出 requirements.txt 缺失的第三方依赖。

跑法：  C:\\Python312\\python.exe dev_tools\\audit_deps.py
输出：  已安装 / 缺失 / 未被 requirements 覆盖 三张表，退出码 1 表示有缺失。

打包前必跑：PyInstaller 不会自动带依赖，缺一个就是运行期 ModuleNotFoundError。
"""
import ast
import os
import sys
import importlib.util
# 控制台编码兜底：Windows 默认代码页是 GBK，本脚本的 ✔/✘/⚠ 等字符会让
# print 抛 UnicodeEncodeError（2026-09-11 实测多处踩过）。只改错误处理为
# replace（编码不动，中文照常可读），编不出的字符降级为 '?'。
for _stream in (sys.stdout, sys.stderr):
    try:
        _stream.reconfigure(errors='replace')
    except (AttributeError, OSError):
        pass


ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# 标准库（sys.stdlib_module_names 在 3.10+ 可用，这里补一份兜底）
try:
    STDLIB = set(sys.stdlib_module_names)
except AttributeError:                                   # pragma: no cover
    STDLIB = set()
STDLIB |= {'os', 'sys', 're', 'json', 'time', 'shutil', 'glob', 'csv',
           'math', 'random', 'uuid', 'threading', 'subprocess', 'collections',
           'argparse', 'dataclasses', 'gzip', 'io', 'itertools', 'functools',
           'pathlib', 'types', 'typing', 'warnings', 'unicodedata', 'hashlib',
           'tempfile', 'traceback', 'string', 'struct', 'copy', 'enum',
           'multiprocessing', 'socket', 'webbrowser', 'email', 'urllib',
           'http', 'xml', 'sqlite3', 'zipfile', 'tarfile', 'base64', 'errno',
           'platform', 'stat', 'textwrap', 'difflib', 'bisect', 'heapq',
           'operator', 'contextlib', 'inspect', 'abc', 'numbers', 'decimal',
           'datetime', 'calendar', 'locale', 'gettext', 'logging'}

# 本地模块（平台自己的包）
LOCAL = {'Virus_Platform_Core', 'app', 'main', 'dev_tools', 'tests'}


def _local_submodules():
    """扫描平台本地包目录下的模块名（含全部一级子包）。

    引擎已并入 Virus_Platform_Core/known_virus_suite/（2026-09-10 由顶层
    engines/ 迁入），改用包内相对导入后，此处的顶层名兜底主要用于扫描
    仍在用裸导入的遗留脚本。2026-09-20：mirna_target/consensus.py 的
    「独立脚本方式运行」回退分支用裸名导入 alignment/psrnatarget/rna22/
    tapir/rnahybrid/psrobot —— 只扫 known_virus_suite 时这些名字会被
    误判成第三方依赖，selfcheck ⑥ 因此恒红。改为遍历 Virus_Platform_Core
    下**所有一级子包**（带 __init__.py 的目录），新增子包不必再回来登记。
    """
    names = set()
    core = os.path.join(ROOT, 'Virus_Platform_Core')
    pkgs = [core]
    try:
        pkgs += [os.path.join(core, d) for d in os.listdir(core)
                 if os.path.isdir(os.path.join(core, d))
                 and os.path.isfile(os.path.join(core, d, '__init__.py'))]
    except OSError:
        pass
    # dev_tools/ 与 tests/ 下的本地脚本模块（test_mcp_e2e 等）也按本地模块
    # 处理：它们是仓库内的文件，不是 pip 包，判成"必需依赖未装"是假阳性。
    for extra in ('dev_tools', 'tests'):
        p = os.path.join(ROOT, extra)
        if os.path.isdir(p):
            pkgs.append(p)
    for pkg in pkgs:
        if not os.path.isdir(pkg):
            continue
        for fn in os.listdir(pkg):
            if fn.endswith('.py'):
                names.add(fn[:-3])
    return names


LOCAL_SUBMODULES = None  # 延迟到 ROOT 定义之后初始化

# 导入名 → pip 包名（不一致的才需要列）
PIP_NAME = {
    'flask': 'Flask', 'werkzeug': 'Werkzeug', 'jinja2': 'Jinja2',
    'plotly': 'plotly', 'matplotlib': 'matplotlib', 'numpy': 'numpy',
    'pandas': 'pandas', 'Bio': 'biopython', 'dna_features_viewer':
        'dna_features_viewer', 'pyrodigal': 'pyrodigal', 'pyhmmer': 'pyhmmer',
    'primer3': 'primer3-py', 'selenium': 'selenium', 'requests': 'requests',
    'openpyxl': 'openpyxl', 'PIL': 'Pillow', 'lxml': 'lxml',
    'sklearn': 'scikit-learn', 'scipy': 'SciPy', 'yaml': 'PyYAML',
    'psutil': 'psutil', 'tqdm': 'tqdm', 'xgboost': 'xgboost',
    'openai': 'openai', 'docx': 'python-docx', 'pptx': 'python-pptx',
    'pypdf': 'pypdf', 'fitz': 'PyMuPDF', 'cv2': 'opencv-python',
    'streamlit': 'streamlit', 'pycirclize': 'pycirclize',
    # 模块名 webview ↔ pip 包名 pywebview（app.py 的桌面窗口壳，import 时兜底回退浏览器）
    'webview': 'pywebview',
    # ViennaRNA 的 Python 绑定顶层模块名是 RNA（requirements.txt 已声明 viennarna）
    'RNA': 'viennarna',
}

# 已知可选依赖：缺失时功能降级，不阻断运行。
# 判定依据是「代码里有 try/except ImportError 兜底或功能可降级」，且
# requirements.txt 里通常只以注释给出安装命令。本表同时决定
# 「requirements.txt 未声明」的告警口径：表内的包不再计入。
# tkinterdnd2 只被 scripts/drag_proxy.py（独立的拖拽代理，非平台本体）使用，
# 且「启动拖拽代理.bat」检测到缺失时会自动 pip 安装 —— 它没装不该让
# 环境自检报 ✘（此前它不在本表，导致 selfcheck 恒为「存在未就绪项」）。
OPTIONAL = {'selenium', 'gbdraw', 'dna_features_viewer', 'pyhmmer',
            'primer3', 'openai', 'psutil', 'pycirclize',
            'distinctipy', 'taxburst', 'tkinterdnd2', 'adjustText',
            # treetime：phylodyn_trees.treetime_available() 探测，缺时走
            # 自研回归兜底（RTT 真引擎降级为本地实现），功能不中断
            'treetime',
            # seaborn：known_virus_suite/kv_fig.py 的 try/except ImportError
            # 导入（HAS_SEABORN），缺时回退 matplotlib 原生的箱线/色板实现
            # ——2026-10-03 该文件并入后未登记，导致自检 ⑥ 恒报
            # "requirements.txt 未声明 1 个: seaborn"（假警报）。
            'seaborn'}


def iter_py_files():
    # docs/bioaider_ref_src 是 BioAider 的参考源码（他项目 GUI，仅供对照），
    # 不是平台运行依赖；tests/ 的造数脚本同理不参与打包；
    # git-repo/ vendor/ _archive/ 是第三方镜像与归档，也不算平台依赖
    # （曾因此把 git-repo/MultiVirusConsensus 的 pysam 误报为必需依赖）。
    # _bioaider_re/ 是 BioAider.exe 逆向工程的工作目录（解包/解密/反编译），
    # 同样不是平台源码 —— 2026-09-11 它出现在仓库根目录后，把 BioAider 自己的
    # PySide2 / custom_libraries / sdt_identity 等误报成了"必需依赖未装"。
    # 3rd/ 是外部依赖（含 2026-09-11 起内置的绿色版 Python 3rd/python，其
    # site-packages 有几万个第三方 .py），里面的 import 不代表平台依赖。
    # 注意：**任何新增的仓库根级临时/参考目录都要在这里登记**，否则会污染审计。
    skip = {'dist', 'tools', 'open-virome', 'node_modules', '.git',
            '__pycache__', 'build', 'host-db', 'virus-db', 'databases',
            'results', 'downloads', 'logs', 'tasks', 'tool_runs',
            'meta_search', 'logan', 'submissions', 'fastq', 'uploads',
            'bin', 'logan', 'docs', 'tests',
            'git-repo', 'vendor', '_archive', '.pytest_cache',
            '_bioaider_re', '3rd',
            # 仓库根级临时/基准/备份目录（非平台运行依赖；2026-09-18
            # _bench_kvs/ 里的本地模块 bench_engines 曾被误报为必需依赖）
            '_bench_kvs', '_bench_kvs2', '_bench_cons',
            '_audit_20260916', '_backup_before_raxmlng_20260916',
            # analysis/ 是运维脚本与一次性实验工作目录（_gcva_phylogeo.py
            # 用 paramiko 连远端跑任务，平台代码不引用它）；2026-10-20 加入，
            # 此前它让 paramiko 被误报成"必需依赖未装"。
            'analysis',
            # dev_tools/archive_202610/ 是 2026-10 的归档探针脚本，平台不
            # 引用；里面 fill_species_gaps.py 还用着已移除的 polars，扫到会
            # 让自检误报"必需依赖未装 polars"（2026-10-20）。
            'archive_202610',
            'archive'}
    for cur, dirs, files in os.walk(ROOT):
        dirs[:] = [d for d in dirs if d not in skip and not d.startswith('.')]
        for fn in files:
            if fn.endswith('.py'):
                yield os.path.join(cur, fn)


def collect_imports():
    """返回 {顶层模块名: [引用它的文件相对路径]}，并跳过 try/except ImportError 块。"""
    found = {}
    for path in iter_py_files():
        try:
            with open(path, encoding='utf-8', errors='replace') as f:
                tree = ast.parse(f.read(), filename=path)
        except (SyntaxError, UnicodeDecodeError, ValueError):
            continue
        rel = os.path.relpath(path, ROOT).replace(os.sep, '/')
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                for a in node.names:
                    top = a.name.split('.')[0]
                    found.setdefault(top, []).append(rel)
            elif isinstance(node, ast.ImportFrom):
                if node.level:            # 相对导入（from .xxx）→ 本地
                    continue
                if node.module:
                    top = node.module.split('.')[0]
                    found.setdefault(top, []).append(rel)
    return found


def _read_req(path):
    """读取 requirements 文件里的包名（小写、`-` 归一为 `_`）。

    忽略注释行、空行与 `-r other.txt` 这类指令行。
    """
    out = set()
    if not os.path.isfile(path):
        return out
    with open(path, encoding='utf-8') as f:
        for line in f:
            line = line.split('#', 1)[0].strip()
            if not line or line.startswith('-'):
                continue
            name = (line.split('==')[0].split('>=')[0]
                        .split('<=')[0].split('~=')[0].split('[')[0])
            if name.strip():
                out.add(name.strip().lower().replace('-', '_'))
    return out


def analyze():
    """完整分析：返回 {third, missing, uncovered, installed, hard_missing}。"""
    found = collect_imports()
    local_sub = _local_submodules()   # 本地包内的子模块（kv_*.py 等）
    third = {m: fs for m, fs in found.items()
             if m not in STDLIB and m not in LOCAL
             and m not in local_sub and not m.startswith('_')}

    req_path = os.path.join(ROOT, 'requirements.txt')
    req = _read_req(req_path)
    # 开发/测试依赖另有一份清单（requirements-dev.txt）。dev_tools/ 与
    # tests/ 里的 import 也会被 collect_imports 扫到（dev_tools 未列入 skip），
    # 所以 PyInstaller / playwright / psutil 这类只属于开发链的包要按
    # dev 清单判定，否则会被误报成「requirements.txt 未声明」。
    dev_req = _read_req(os.path.join(ROOT, 'requirements-dev.txt'))

    installed, missing = {}, []
    for m in sorted(third):
        ok = importlib.util.find_spec(m) is not None
        installed[m] = ok
        if not ok:
            missing.append(m)

    uncovered = []
    for m in sorted(third):
        pip = PIP_NAME.get(m, m)
        key = pip.lower().replace('-', '_')
        # 可选依赖（OPTIONAL）本就不要求写进 requirements.txt —— 它们的功能
        # 缺失时自动降级，清单里只用注释给出安装命令。
        if m in OPTIONAL or key in dev_req:
            continue
        if key not in req:
            uncovered.append((m, pip))

    # hard_missing = 缺失且**运行时**需要。开发链依赖（requirements-dev.txt
    # 里声明的 pyinstaller/playwright/psutil/pytest 等）缺失不该算"必需
    # 未装"——源码分发模式下用户不需要它们，报红只会误导。
    # 2026-10-20 修：此前只在 uncovered 里排除 dev 依赖，hard_missing
    # 没排除，精简掉 PyInstaller 后自检就误报"必需依赖未装"。
    hard_missing = [m for m in missing
                    if m not in OPTIONAL
                    and PIP_NAME.get(m, m).lower().replace('-', '_') not in dev_req]
    return {'third': third, 'missing': missing, 'uncovered': uncovered,
            'installed': installed, 'hard_missing': hard_missing,
            'req': req,
            'all_ready': not hard_missing and not missing}


def missing_required():
    """返回缺失的必需依赖名列表（空 = 全部就绪）。供 selfcheck / CI 用。"""
    return analyze()['hard_missing']


def main():
    res = analyze()
    third, missing = res['third'], res['missing']
    uncovered, installed = res['uncovered'], res['installed']
    hard_missing = res['hard_missing']
    req = res['req']

    print('扫描到第三方依赖 %d 个（引用文件已去重统计）' % len(third))
    print('\n── 未安装 ──')
    if not missing:
        print('  （无）')
    for m in missing:
        tag = '可选' if m in OPTIONAL else '必需'
        print('  ✘ %-24s [%s]  pip install %s   ← %s'
              % (m, tag, PIP_NAME.get(m, m), ', '.join(sorted(set(third[m]))[:3])))

    print('\n── requirements.txt 未覆盖 ──')
    if not uncovered:
        print('  （无）')
    for m, pip in uncovered:
        tag = '可选' if m in OPTIONAL else '必需'
        print('  △ %-24s [%s]  建议加入: %s' % (m, tag, pip))

    print('\n── 已安装且已覆盖 ──')
    for m in sorted(third):
        pip = PIP_NAME.get(m, m)
        if installed[m] and pip.lower().replace('-', '_') in req:
            print('  ✔ %s' % m)

    hard_missing = [m for m in missing if m not in OPTIONAL]
    print('\n' + '=' * 46)
    if hard_missing:
        print('缺失必需依赖 %d 个: %s' % (len(hard_missing), ', '.join(hard_missing)))
        print('安装: pip install ' + ' '.join(PIP_NAME.get(m, m)
                                              for m in hard_missing))
        return 1
    if missing:
        print('仅缺可选依赖（对应功能会降级）: %s' % ', '.join(missing))
    else:
        print('全部依赖已安装')
    return 0


if __name__ == '__main__':
    sys.exit(main())
