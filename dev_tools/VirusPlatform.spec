# -*- mode: python ; coding: utf-8 -*-
# spec 位于 dev_tools/，平台根目录 = SPECPATH 的上级（由 dev_tools/package.py 调用）
import os

ROOT = os.path.abspath(os.path.join(SPECPATH, '..'))

from PyInstaller.utils.hooks import collect_data_files, collect_submodules

def _tree(src, dest, exclude=('.mimosa', '__pycache__', '.git', '.github')):
    """把目录展开成 datas 条目（逐文件），顺带跳过编辑器/AI 工具残留。

    为什么不能直接写 `(目录, 目标)`：那样 PyInstaller 会**原样整棵拷贝**，
    把 `webapp/templates/.mimosa`（编辑器 hook 缓存）这类残留一起打进分发版
    （实测确实在包里）。`datas` 不支持 ignore 参数，只能自己走 os.walk。

    目标路径统一用 `/` 分隔（PyInstaller 的 datas 目标用 posix 风格）。
    """
    import fnmatch
    # *.bak_YYYYMMDD_* / *.bak 是开发期的手工备份（webapp 下有 20+ 个），
    # 属于源码树但不该进分发版 —— 打进包只会白白增大体积、暴露内部迭代痕迹。
    bak_patterns = ('*.bak', '*.bak_*')
    out = []
    for cur, dns, files in os.walk(src):
        dns[:] = [d for d in dns if d not in exclude]
        rel = os.path.relpath(cur, src)
        tgt = dest if rel == '.' else dest + '/' + rel.replace(os.sep, '/')
        for fn in files:
            if fn in exclude:
                continue
            if any(fnmatch.fnmatch(fn, pat) for pat in bak_patterns):
                continue
            out.append((os.path.join(cur, fn), tgt))
    return out


datas = _tree(os.path.join(ROOT, 'webapp'), 'webapp')
datas += [(os.path.join(ROOT, 'README.md'), '.')]
datas += collect_data_files('plotly')
datas += collect_data_files('pycirclize')
datas += collect_data_files('gbdraw')
# pyrodigal_rv 的病毒 metagenomic 模型库是纯数据（meta.json），不带会被
# 运行时的 FileNotFoundError 挡住 —— 同样会被 orf.py 当"未安装"静默跳过。
datas += collect_data_files('pyrodigal_rv', includes=['meta.json'])

# 包内**非 .py 数据文件**不会被 PyInstaller 自动收集（PYZ 里只装代码，
# 包目录本身不会出现在 dist/ 下）。凡是模块用 `__file__` 相对路径定位的
# 资源，都必须显式登记，且**目标路径要与包结构一致**：
#   冻结后模块的 `__file__` = <_MEIPASS>/Virus_Platform_Core/<子包>/<模块>.pyc，
#   于是 `Path(__file__).parent / 'templates'` 会落到
#   <_MEIPASS>/Virus_Platform_Core/ncbi_submit/templates —— 只要把资源
#   按同名相对路径投放到那里即可命中，**不需要改业务代码**。
# 已实测会出问题的两处（均属「提交准备」页）：
#   - ncbi_submit/templates/submission_report.html
#       report_html.py 用 jinja2 FileSystemLoader 加载 → 缺失抛
#       jinja2.exceptions.TemplateNotFound（生成提交报告时）
#   - ncbi_submit/samples/sample_{public,selfseq}.csv
#       store._sample_dir() 定位 → 缺失时 list_samples() 静默返回空列表、
#       用样例新建提交项目则 FileNotFoundError
# 注：kv_config.yaml / _test_logs/*.log 等无人读取，故意不收。
datas += [
    (os.path.join(ROOT, 'Virus_Platform_Core', 'ncbi_submit', 'templates'),
     'Virus_Platform_Core/ncbi_submit/templates'),
    (os.path.join(ROOT, 'Virus_Platform_Core', 'ncbi_submit', 'samples'),
     'Virus_Platform_Core/ncbi_submit/samples'),
]

# pyrodigal 的 CPU 内核（avx2/avx512/sse2/generic/swar64）在
# `pyrodigal/impl/` 下，那是个**没有 __init__.py 的命名空间包**：
# PyInstaller 的静态分析看不到，`collect_submodules('pyrodigal')` 也不返回它们
# （实测只返回 pyrodigal / .cli / .lib / .tests.*）。漏掉则运行时 `import pyrodigal`
# 抛 `ModuleNotFoundError: No module named 'pyrodigal.impl'`，而 orf.py 把它当作
# "未安装"静默跳过 → 阶段④ ORF 预测整段失效，且日志只写"未安装"，极难排查。
# 按目录下实际存在的扩展模块动态登记，换 pyrodigal 版本/平台都不用改这里。
_pyro_impl = []
try:
    import glob as _glob
    import pyrodigal as _pd
    _impl_dir = os.path.join(os.path.dirname(_pd.__file__), 'impl')
    for _f in (_glob.glob(os.path.join(_impl_dir, '*.pyd'))
               + _glob.glob(os.path.join(_impl_dir, '*.so'))):
        _name = os.path.basename(_f).split('.')[0]      # avx2.cp312-win_amd64.pyd → avx2
        _pyro_impl.append('pyrodigal.impl.' + _name)
    _pyro_impl = sorted(set(_pyro_impl))
except Exception:
    pass

# pywebview 原生窗口壳（app.py 窗口模式；未安装 pywebview 的构建环境自动跳过，
# 运行时回退浏览器）。WebView2/WinForms 后端是 pywebview 运行时按平台动态
# import 的（importlib），静态分析看不到，必须显式登记，否则 exe 一启动就
# WebViewError 回退浏览器，窗口壳形同虚设。WebView2Loader.dll 等 lib 数据
# 文件不用管 —— hooks-contrib 的 hook-webview 会自动收集。
_webview_hidden = []
try:
    import webview as _wv  # noqa: F401
    _webview_hidden = ['webview.platforms.winforms',
                       'webview.platforms.edgechromium']
except Exception:
    pass

# MCP 控制台（2026-10-24）：mcp_server.py 不被 app.py 静态 import（冻结模式
# 经 `exe --mcp-stdio` 进程内 import），静态分析永远扫不到，必须显式登记；
# mcp_server 的 import 链会把 mcp.server.fastmcp / mcp.types / anyio /
# pydantic 全部静态可达地带入。**不能用 collect_all('mcp')**：它会枚举导入
# 全部子模块，mcp/cli/cli.py 在缺 typer 的构建环境 import 即 sys.exit(1)
# （2026-10-24 实测打包直接失败）；CLI 子包我们用不到，数据文件单独收。
_mcp_datas = collect_data_files('mcp')
datas += _mcp_datas

# ViennaRNA 的 _RNA.pyd 用 MSVC C++ 运行库（导入表依赖 MSVCP140.dll），
# 但平台其余扩展都只用 VCRUNTIME140 —— PyInstaller 的依赖分析因此没有把
# MSVCP140.dll 收进 _internal（2026-09-18 实测：冻结版 import RNA 报
# "Failed to load ViennaRNA RNAlib Python wrapper"，dsRNA 全链路 / miRNA
# 靶标预测均失效）。源码模式不受影响（System32 里有）。这里从 System32
# 显式随包携带（VC 运行库可随应用再分发，PyInstaller 对其它应用也是这么做的）。
import os as _os
import sys as _sys
_binaries = []
if _sys.platform == 'win32':
    _windir = _os.environ.get('WINDIR', r'C:\Windows')
    for _dll in ('msvcp140.dll', 'vcruntime140_1.dll'):
        _p = _os.path.join(_windir, 'System32', _dll)
        if _os.path.isfile(_p):
            _binaries.append((_p, '.'))

# mcp SDK 的二进制附件改由静态分析随 mcp_server 的 import 链带入
# （collect_all('mcp') 会触发 mcp.cli 的 import-time sys.exit(1)，已弃用）

# primer3-py（T7 引物设计）**故意不进分发版**（2026-09-19 实测）：
# 本机构建时装在 C:\ProgramData\VirusPlatformPrimer3（sys.path 注入），
# 尝试整体收进冻结包后，import primer3.bindings 直接 0xC0000005 段错误
# （p3helpers/thermoanalysis 两个 pyd 在冻结环境初始化即崩，原因未明）。
# 处理：primer_design.py 已改为惰性导入 + 优雅降级（不可用时 T7/引物
# 热力学功能提示 primer3 不可用，其余功能不受影响），分发版不携带
# primer3 —— 待查清 pyd 冻结崩溃原因后再恢复随包。


a = Analysis(
    [os.path.join(ROOT, 'app.py')],
    pathex=[],
    binaries=_binaries,
    datas=datas,
    hiddenimports=['Virus_Platform_Core.public_meta.search_engine',
                   'Virus_Platform_Core.public_meta.info_engine',
                   # 病毒样本预筛引擎（Web 侧以 --run-engine 子进程调用，
                   # 见 engine_entry.py 的 _ENGINES）
                   'Virus_Platform_Core.public_meta.prescreen_engine',
                   # 病毒谱：把预筛产物聚合成「病毒视角」（Web 侧直接 import，
                   # 不是子进程引擎；模块级 import 了 taxonomy，须显式登记）
                   'Virus_Platform_Core.public_meta.virome',
                   'Virus_Platform_Core.logan_submit', 'Virus_Platform_Core.lovis4u_run',
                   'Virus_Platform_Core.ncbi_submit.unified_metadata', 'Virus_Platform_Core.engine_entry',
                   'Virus_Platform_Core.cdd_search',   # CDD 低内存分片回退（函数内惰性导入）
                   # CRISPR 检测试剂设计（tool_jobs 内惰性导入，静态分析
                   # 不到；2026-10-07 生产审查补登。冻结版无 primer3 时
                   # 自动降级枚举路，无 ViennaRNA 时跳过结构过滤）
                   'Virus_Platform_Core.crispr_design',
                   # AI 助手（对话驱动智能代理，2026-10-08 补登）：web/ai_assist.py
                   # 顶层 import 了 chat/digest/broadcast/provider（静态分析可达），
                   # 但 **executors 只被函数内惰性 import 引用**
                   # （chat.py 的 propose_action 分派、tools.py 的 dispatch、
                   #   ai_assist.py 的两条 /api/ai/actions 路由）——
                   # 与本项目其它惰性导入引擎（crispr_design / dsrna_* / kv_*）
                   # 同因，须显式登记；否则冻结版里「AI 执行单 / 自动重试 /
                   # 自动调参」会 ImportError，而对话与播报看起来一切正常
                   # （静默失效，最难发现）。其余 5 个一并登记以防静态分析漏。
                   'Virus_Platform_Core.ai_agent.executors',
                   'Virus_Platform_Core.ai_agent.broadcast',
                   'Virus_Platform_Core.ai_agent.audit',
                   'Virus_Platform_Core.ai_agent.chat',
                   'Virus_Platform_Core.ai_agent.digest',
                   'Virus_Platform_Core.ai_agent.provider',
                   'Virus_Platform_Core.ai_agent.tools',
                   # 轨道 B（BEAST 交接，2026-10-08 接线）：beast_handoff 只被
                   # web/tool_jobs.py 的 _tool_job_beastexport 函数内 import 引用，
                   # beast_io 只被 phylogeo.import_beast 函数内 import 引用 ——
                   # 均为惰性导入，随本项目惯例显式登记，避免冻结版里
                   # 「打包 BEAST 作业包 / BEAST 产物导入」静默 ImportError。
                   'Virus_Platform_Core.beast_handoff',
                   'Virus_Platform_Core.beast_io',
                   # 已知病毒引擎（2026-09-10 由顶层 engines/ 并入包内子包；
                   # 用相对导入后 PyInstaller 静态分析不到，须显式登记）
                   'Virus_Platform_Core.known_virus_suite.kv_common',
                   'Virus_Platform_Core.known_virus_suite.kv_engines',
                   'Virus_Platform_Core.known_virus_suite.kv_identify',
                   'Virus_Platform_Core.known_virus_suite.kv_filter',
                   'Virus_Platform_Core.known_virus_suite.kv_consensus',
                   'Virus_Platform_Core.known_virus_suite.kv_plot',
                   'Virus_Platform_Core.known_virus_suite.kv_variant',
                   'Virus_Platform_Core.known_virus_suite.kv_variant_plot',
                   'Virus_Platform_Core.known_virus_suite.kv_variant_evo',
                   'Virus_Platform_Core.known_virus_suite.known_virus_suite',
                   'Virus_Platform_Core.kv_stage',
                   # dsRNA 设计全链路（2026-09 接线，2026-10-03 恢复）：三模块均为
                   # 函数内惰性 import，PyInstaller 静态分析不到，须显式登记；
                   # RNA 为 viennarna 的编译扩展（dsRIP 效价打分的折叠能项，
                   # 亦供 miRNA 靶标预测的 RNAhybrid/psRNATarget/TAPIR 复刻引擎）
                   'Virus_Platform_Core.dsrna_scoring',
                   'Virus_Platform_Core.dsrna_pipeline',
                   'Virus_Platform_Core.dsrna_offtarget',
                   'RNA',
                   'main',            # exe --cli <子命令> 用（打包后自检）
                   'mcp_server',      # exe --mcp-stdio 用（MCP 控制台进程内拉起）
                   'mcp', 'mcp.types', 'mcp.server.fastmcp',   # mcp SDK 主链显式兜底
                   'gbdraw', 'Virus_Platform_Core.public_meta.landscape_plot', 'seaborn']
                   + collect_submodules('gbdraw')
                   + _pyro_impl
                   + _webview_hidden,

    hookspath=[],
    hooksconfig={},
    # 冻结分发：引导期预载 ViennaRNA（时机相关的 DllMain 冲突，见钩子内注释）
    runtime_hooks=[os.path.join(ROOT, 'dev_tools', 'rna_early_import_hook.py')],
    # 体积治理：以下模块平台自身从不 import（2026-09-11 已逐项 grep 复核
    # Virus_Platform_Core/ + app.py + main.py，均 0 命中），只是被
    # pandas/plotly/gbdraw 的依赖链顺带拖进来。实测首版打包 2.05GB，
    # 其中 _internal 834MB（llvmlite 115MB / pyarrow 79MB / bokeh 21MB /
    # botocore 18MB …）。注意不要排除 selenium（LOGAN 批量提交必需）。
    #
    # ⚠️ 不要把 polars / _polars_runtime_32 加回来（2026-09-11 实测踩坑）：
    #   known_virus_suite/kv_filter.py 与 kv_identify.py 在**模块级**
    #   `import polars as pl`，缺失时 `raise SystemExit("需要 polars")`；
    #   而 known_virus_suite.py 又在模块级 import 这两个模块 —— 一旦排除，
    #   「已知病毒识别与定量」整段引擎在打包版里**完全失效且完全静默**。
    #   （当时还被 engine_entry.run_engine 吞成 return 0，双重静默；
    #    engine_entry 已修，但依赖本身仍必须打进包。）
    #   polars 本体仅 8.1MB，但 polars/_plr.py 硬依赖 _polars_runtime_32
    #   （174MB 的 .pyd，就是它的核心运行时，无回退），两者必须同进同出，
    #   合计约 +183MB —— 这是该功能的必要成本。
    #
    # ⚠️ kaleido 已在 2026-09-16 移出 excludes：迁移 GIF 导出
    #   （Virus_Platform_Core/phylogeo_gif.py）会**惰性 import** 它做逐帧出图，
    #   排掉后打包版点「导出 GIF 动画」只会得到 ImportError。kaleido 本体几
    #   百 KB，可忽略。它自身不随包带浏览器：kaleido 1.x 走 plotly 的
    #   choreographer，优先用系统已装的 Chrome/Edge；都没有的那台机器上，
    #   **首次**渲染会联网下载一份 Chromium（本机已下载到
    #   %LOCALAPPDATA%\plotly\choreographer\deps\chrome-win64\）。离线机器
    #   首次渲染失败时由 build_migration_gif 首帧给出人话提示。
    excludes=['pytest', '_pytest', 'IPython', 'jedi', 'parso', 'tkinter',
              '_tkinter', 'numba', 'llvmlite',
              'pyarrow', 'bokeh', 'botocore', 'boto3', 's3fs', 'sqlalchemy',
              'sympy', 'duckdb', 'reportlab', 'docutils', 'sphinx',
              'matplotlib.tests', 'numpy.testing',
              'pandas.tests', 'scipy.io.tests',
              # Qt 绑定（2026-10-24）：构建环境同时装了 PyQt5 与 PySide6
              #（传递安装），PyInstaller 拒绝多 Qt 绑定共存直接中止构建；
              # 平台是 Flask 网页 UI + pywebview(edgechromium/winforms)，
              # 零处 import Qt，显式排除两个绑定。
              'PyQt5', 'PyQt5_sip', 'PySide6', 'shiboken6', 'qtpy',
              # 开发/测试依赖（requirements-dev.txt）：真浏览器探针用，
              # 平台本体从不 import。显式排除，防止本机装了之后被动入包。
              'playwright',
              # primer3 本体不入包（2026-09-20）：site-packages 里的副本位于
              # 非 ASCII 路径，pyd 在冻结环境下构造 ThermoAnalysis 会段错误；
              # app.py 启动时会把 primer3_pkg\primer3 镜像到
              # %ALLUSERSPROFILE%\VirusPlatformPrimer3App 并插到 sys.path
              # 最前（见 package.py 的 primer3_pkg 拷贝步骤）。
              'primer3', 'primer3.bindings', 'primer3.argdefaults',
              'primer3.p3helpers', 'primer3.thermoanalysis'],
    noarchive=False,
    optimize=0,
)
pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name='VirusPlatform',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    console=True,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
)
coll = COLLECT(
    exe,
    a.binaries,
    a.datas,
    strip=False,
    upx=True,
    upx_exclude=[],
    name='VirusPlatform',
)
