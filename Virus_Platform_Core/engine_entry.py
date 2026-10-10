# -*- coding: utf-8 -*-
"""冻结分发（PyInstaller exe）的内部引擎入口。

exe 以 `VirusPlatform.exe --run-engine <脚本名> [args...]` 方式调用时，
app.py 顶部把控制权转交本模块，用 runpy 在进程内以 __main__ 语义执行
对应引擎（与源码模式下 `python 某引擎.py args...` 等价）。
引擎脚本必须在 VirusPlatform.spec 的 hiddenimports 中登记才会被打包。
"""
import runpy
import sys

_ENGINES = {
    'search_engine.py': 'Virus_Platform_Core.public_meta.search_engine',
    'info_engine.py': 'Virus_Platform_Core.public_meta.info_engine',
    'prescreen_engine.py': 'Virus_Platform_Core.public_meta.prescreen_engine',
    # 病毒谱只有排查用 CLI（Web 侧是直接 import，不走子进程）；
    # 登记它是为了冻结分发里也能 `--run-engine virome.py <项目>/prescreen` 复查
    'virome.py': 'Virus_Platform_Core.public_meta.virome',
    'logan_submit.py': 'Virus_Platform_Core.logan_submit',
    'unified_metadata.py': 'Virus_Platform_Core.ncbi_submit.unified_metadata',
    'known_virus_suite.py':
        'Virus_Platform_Core.known_virus_suite.known_virus_suite',
}


def run_engine(script_name):
    mod = _ENGINES.get(script_name)
    if not mod:
        print(f'未知引擎脚本: {script_name}', file=sys.stderr)
        return 2
    argv = [mod] + sys.argv[3:]
    try:
        sys.argv = argv
        runpy.run_module(mod, run_name='__main__', alter_sys=True)
        return 0
    except SystemExit as e:                      # 引擎主动退出码透传
        if e.code is None:
            return 0                             # sys.exit() = 正常结束
        if isinstance(e.code, int):
            return e.code
        # `sys.exit("消息")` / `raise SystemExit("消息")`：必须把消息打出来，
        # 并且必须返回非 0。**不能一律 return 0** —— 原先那种写法会把引擎
        # 启动期的致命错误变成"成功"：kv_filter / kv_identify 在模块级
        # `raise SystemExit("需要 polars: pip install polars")`，异常被这里
        # 捕获后既不打消息、又返回 0，而调用方 kv_stage.py 只在
        # returncode != 0 时报错 —— 结果「已知病毒识别与定量」整段静默失效
        # （零输出、退出码 0、界面显示成功、无任何产物）。实测于 2026-09-11。
        print(e.code, file=sys.stderr)
        return 1
    except ImportError as e:                     # 可选依赖未打包
        print(f'引擎 {script_name} 不可用（依赖未打包）: {e}', file=sys.stderr)
        return 3
