# -*- coding: utf-8 -*-
"""植物病毒分析平台 - 本地 Web GUI（组合根）。

启动（两种模式，可显式选）:
    python app.py            或 双击 启动平台-桌面窗口.bat → 默认独立桌面窗口
    python app.py --gui      或 双击 启动平台-桌面窗口.bat → 强制桌面窗口
    python app.py --web      或 双击 启动平台-网页.bat     → 强制系统浏览器
  桌面窗口由 pywebview(WebView2) 壳住本地页面，无地址栏，观感即本地应用；
  未安装 pywebview / 缺 WebView2 / 窗口启动失败时**自动回退**浏览器。
  兼容旧写法：环境变量 VP_GUI=browser|web 等价于 --web，VP_GUI=gui 等价于 --gui；
  命令行参数优先于环境变量。打包版同理：VirusPlatform.exe --web 即浏览器模式。
  仅监听本机回环地址；长任务在后台线程执行，页面实时轮询进度。
  所有文件写入经由 utils.safe_open（路径校验+限平台内），os.startfile 前一律 check_path。

本文件只做「组装」——创建 Flask 应用、注册 blueprint、启动服务。
业务逻辑分两层：
  - Virus_Platform_Core/        计算与流程（pipeline / kunpeng / phylo / orf_annot …）
  - Virus_Platform_Core/web/    HTTP 边界与任务工厂（各 blueprint）
**请勿在 app.py 新增路由**；新路由写进 Virus_Platform_Core/web/ 对应模块，在此登记 blueprint。
"""
import os
import sys
import threading
import time
import webbrowser
# Windows 注册表把 .svg 关联成 MIME 类型 "image/svg"，而标准类型是
# "image/svg+xml"；Python 的 mimetypes 会照搬注册表，Flask 便以 image/svg
# 发出去，Chromium 直接拒绝解码（<img> 拿到 naturalWidth=0 的裂图）。
# 这里在造 Flask app 之前显式注册，全站 SVG 才能当图片显示
# （首页封面主视觉、⑨ 基因组图、示例结果都用 <img> 引 SVG）。
import mimetypes
mimetypes.add_type('image/svg+xml', '.svg')

from flask import Flask, abort, request
from flask.json.provider import DefaultJSONProvider

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

if sys.platform == 'win32':
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')

# 冻结分发下 multiprocessing（SDT 精确矩阵的 ProcessPoolExecutor 等）子进程
# 会重新启动 exe，必须在任何业务代码前调用 freeze_support 拦截子进程引导；
# 非冻结环境是空操作
if __name__ == '__main__':
    import multiprocessing
    multiprocessing.freeze_support()

# 冻结分发：primer3-py 的 pyd 在非 ASCII 安装路径下构造 ThermoAnalysis 会
# 直接段错误（同 sitecustomize 的镜像机制口径：源码模式把包镜像到
# %ALLUSERSPROFILE%\VirusPlatformPrimer3）。程序包自带 primer3_pkg/primer3
# 副本，启动时镜像到纯 ASCII 的用户公共目录并插入 sys.path——必须在任何
# 业务模块可能触发 primer3 导入之前完成。镜像失败仅影响引物热力学功能。
# ⚠️ 位置必须在下方 --run-engine / --cli 派发之前：--run-engine 跑 dsRNA T7
# 引物设计、--cli importprobe 子进程探测 primer3，两者都依赖镜像路径已入 sys.path。
if getattr(sys, 'frozen', False):
    try:
        _p3_src = os.path.join(os.path.dirname(sys.executable),
                               'primer3_pkg', 'primer3')
        if os.path.isdir(_p3_src):
            import shutil as _sh
            _p3_base = os.path.join(os.environ.get('ALLUSERSPROFILE',
                                                   r'C:\ProgramData'),
                                    'VirusPlatformPrimer3App')
            _p3_dst = os.path.join(_p3_base, 'primer3')
            _ver_src = os.path.join(os.path.dirname(_p3_src), 'VERSION.txt')
            _ver_dst = os.path.join(_p3_base, 'VERSION.txt')
            _need = not os.path.isdir(_p3_dst)
            # 旧版程序包漏带 primer3_config（热力学参数表）时镜像会完整但
            # 段错误 —— 检测到配置目录缺失也强制重镜像（自愈升级路径）
            if not _need and not os.path.isdir(
                    os.path.join(_p3_dst, 'src', 'libprimer3',
                                 'primer3_config')):
                _need = True
            try:
                if os.path.isfile(_ver_src) and os.path.isfile(_ver_dst) and \
                   open(_ver_src, 'rb').read() != open(_ver_dst, 'rb').read():
                    _need = True
            except OSError:
                pass
            if _need:
                _sh.rmtree(_p3_dst, ignore_errors=True)
                os.makedirs(_p3_base, exist_ok=True)
                _sh.copytree(_p3_src, _p3_dst)
                if os.path.isfile(_ver_src):
                    _sh.copyfile(_ver_src, _ver_dst)
            if _p3_base not in sys.path:
                sys.path.insert(0, _p3_base)
    except Exception:
        pass  # 镜像失败仅引物热力学功能降级，平台其余功能不受影响

# 冻结分发：`VirusPlatform.exe --run-engine <脚本名> [args...]` → 进程内执行引擎
if getattr(sys, 'frozen', False) and len(sys.argv) > 2 and sys.argv[1] == '--run-engine':
    from Virus_Platform_Core.engine_entry import run_engine
    raise SystemExit(run_engine(sys.argv[2]))

# 冻结分发：`VirusPlatform.exe --cli <子命令> [args...]` → 进程内跑 CLI
# （打包后没有 python，无法 `python main.py selfcheck`；这条给打包验证和
#  用户自检用，例如 `VirusPlatform.exe --cli selfcheck`）
if getattr(sys, 'frozen', False) and len(sys.argv) > 2 and sys.argv[1] == '--cli':
    import main as _cli
    sys.argv = ['main.py'] + sys.argv[2:]
    _cli.main()
    raise SystemExit(0)

# MCP 控制台的服务器入口：`VirusPlatform.exe --mcp-stdio`（或源码模式
# `python app.py --mcp-stdio`）→ 进程内跑 mcp_server 的 stdio 服务。
# 冻结分发没有独立 python 可执行脚本，mcp_console._ensure 在 frozen 下
# 就是用本入口拉起（mcp_server 模块与 mcp SDK 已随 PyInstaller 入包）。
# 必须在 Flask / 网页启动之前接管进程：stdout 是 MCP 协议通道，本分支
# 之后不允许再有任何非协议输出。
if len(sys.argv) > 1 and sys.argv[1] == '--mcp-stdio':
    import mcp_server as _mcp_srv
    _mcp_srv.mcp.run()
    raise SystemExit(0)

# webapp 目录 / 配置单例 / tool_runs 根统一由 Virus_Platform_Core/web/state.py 提供，
# 各 blueprint 也从那里取，避免 app ↔ blueprint 循环导入。
from Virus_Platform_Core.web.state import _WWW  # noqa: E402

def _pause_before_exit():
    """退出前暂停（有控制台时等回车）；无 stdin 的场景（pythonw/计划任务/
    测试子进程）input 会抛 EOFError 直接崩栈——吞掉，安静退出。"""
    try:
        input('按回车键退出...')
    except (EOFError, OSError):
        pass


app = Flask(__name__, template_folder=os.path.join(_WWW, 'templates'),
            static_folder=os.path.join(_WWW, 'static'))
app.config['JSON_AS_ASCII'] = False
# 上传请求体硬上限（16GB）：此前未设，submit/meta 等上传点在请求线程整表
# 读进内存，误传超大文件会把本机内存/磁盘打满（2026-10-04 审查 E4）。
# 常规 FASTQ/FASTA（宿主基因组 ~2GB）覆盖足够；更大文件走「挂载数据目录」。
app.config['MAX_CONTENT_LENGTH'] = 16 * 1024 ** 3
# 浏览器固定请求站点根的 /favicon.ico；不挂路由的话每次页面加载都有一条
# 404 噪音（真身在 static/brand/favicon.ico）
@app.route('/favicon.ico')
def _favicon():
    from flask import send_file
    return send_file(os.path.join(_WWW, 'static', 'brand', 'favicon.ico'),
                     mimetype='image/x-icon')
# 模板热重载：debug=False 下 Jinja2 默认缓存模板，改 .html 后需重启才生效；
# 打开后每次请求按 mtime 检测，本地平台开销可忽略，改模板即刷新可见。
app.config['TEMPLATES_AUTO_RELOAD'] = True

# NaN/±Inf 清洗的唯一实现（原随病毒浏览器归档，现收编进 utils.py）；
# 非 JSON 出口（/tool_runs 原样发文件）的落盘清洗见 web/tool_jobs._dump_json。
from Virus_Platform_Core.utils import nan_to_null as _nan_to_null  # noqa: E402


class _CleanJSONProvider(DefaultJSONProvider):
    """全局 JSON 出口：清洗 NaN/±Inf 后再序列化（平台 JSON 铁律）。

    Flask 默认放行裸 `NaN`（`allow_nan=True`）——那是**非法 JSON**：
    Python 端 json.loads 能容忍（本地自测假通过），浏览器 JSON.parse 直接
    报错、整份响应被丢弃。来源多为空集合求均值/占比（pandas 聚合）。
    这里统一把 NaN/±Inf 换成 null，并以 `allow_nan=False` 兜底：若还有
    漏网之鱼会带类型报错，而不是静默产出非法 JSON。
    """

    ensure_ascii = False

    def default(self, o):
        try:
            return super().default(o)
        except TypeError:
            if hasattr(o, 'tolist'):        # numpy 数组
                return o.tolist()
            if hasattr(o, 'item'):          # numpy 标量
                return o.item()
            raise

    def dumps(self, obj, **kwargs):
        kwargs.setdefault('allow_nan', False)
        return super().dumps(_nan_to_null(obj), **kwargs)


app.json = _CleanJSONProvider(app)


# ------------------------------------------------------------------
# 本地 CSRF / DNS-rebinding 防护
# 服务只监听 127.0.0.1，但这挡不住两类远程攻击：
#   1) CSRF：用户浏览器里的任意恶意网页可静默 POST http://127.0.0.1:<port>
#      （删样品、全局重置、改设置…），响应虽不可读，破坏已经生效；
#   2) DNS rebinding：攻击者域名解析到 127.0.0.1 后，同源策略即被绕过，
#      恶意页面可读取响应（任意文件浏览/序列查看 → 数据外泄）。
# 两道校验（对每个请求，纯头部判断，开销可忽略）：
#   ① Host 必须是回环主机名 —— rebinding 时 Host 是攻击者域名，直接 403；
#   ② 带 Origin/Referer 的请求必须同源 —— 浏览器跨站 POST 必带 Origin，
#      恶意网页驱动不了任何端点。非浏览器客户端（curl/本机脚本）不带
#      Origin/Referer，正常放行。
# ------------------------------------------------------------------
_LOOPBACK_HOSTS = {'127.0.0.1', 'localhost', '::1'}

# 服务实际绑定的端口，由 main() 在 _pick_port() 之后、_start_server() 之前写入。
# None = 尚未起服务（test_client / 嵌入调用 / 打包自检），此时 Origin 校验
# 退化为主机名校验，保证既有测试与本地脚本行为不变（见下）。
_BOUND_PORT = None


def _url_hostname(v):
    try:
        from urllib.parse import urlparse
        return (urlparse(v).hostname or '').lower()
    except ValueError:
        return ''


def _origin_same_platform(value):
    """Origin/Referer 是否来自本平台源（http://<回环主机>:<绑定端口>）。

    只比 hostname 不够：本机**其它端口**上的任意页面同样是 127.0.0.1，
    浏览器视角与平台并不同源，却照样能驱动平台的破坏性端点（本地 CSRF）。
    端口在启动时已知，这里按「完整源」比对（scheme + host + port 三项）；
    未绑定端口（_BOUND_PORT is None）时退回主机名校验，避免 test_client
    与本地冒烟脚本被误伤。
    """
    try:
        from urllib.parse import urlparse
        u = urlparse(value)
    except ValueError:
        return False
    if (u.hostname or '').lower() not in _LOOPBACK_HOSTS:
        return False
    if _BOUND_PORT is None:
        return True                 # 未起服务：仅回环主机名校验
    if (u.scheme or '').lower() != 'http':
        return False                # 平台只跑 http，其它协议一律拒绝
    try:
        return u.port == _BOUND_PORT
    except ValueError:              # 端口非法（:abc / 超界）
        return False


@app.before_request
def _guard_local_origin():
    host_hdr = (request.host or '') if request else ''
    if host_hdr:
        # 尾部圆点（"localhost.:8765"）是合法 DNS 写法，剥点后再比对
        hostname = host_hdr.rsplit(':', 1)[0].strip('[]').lower().rstrip('.')
        if hostname and hostname not in _LOOPBACK_HOSTS:
            abort(403, '仅允许通过本机回环地址访问平台')
    origin = request.headers.get('Origin') or ''
    if origin:
        # 浏览器对同源 POST 也带 Origin（必须是本平台**完整源**才放行）
        if not _origin_same_platform(origin):
            abort(403, '已拦截跨站请求（Origin 校验失败）')
    else:
        referer = request.headers.get('Referer') or ''
        # 无 Origin 且非浏览器客户端（curl / 本机脚本 / 拖拽代理，都不发
        # Origin）保持既有放行逻辑——拖拽代理依赖它。
        # 威胁模型核对（2026-10-04 审查遗留项复核）：① 浏览器跨站请求
        # （含 CSRF）fetch/form POST 必带 Origin，已在上方完整源校验拦截；
        # ② 无 Origin 的只可能是本机进程，其权限等价于登录用户，回环
        # 绑定已是边界（进一步收紧对 CSRF 无增益，只会打断 test_client
        # 与 CLI）；③ 破坏性路由 14 个全为 POST、无 GET 型（否则跨站
        # 导航不带 Origin 可直达）——三项成立，本放行是安全决策而非疏漏。
        if referer and not _origin_same_platform(referer):
            abort(403, '已拦截跨站请求（Referer 校验失败）')

# ------------------------------------------------------------------
# Blueprint 注册
# 路由清单守卫：python tests/_route_inventory.py
#（锁定 229 条 rule + methods，拆分/新增后必须零差异）
# ------------------------------------------------------------------
from Virus_Platform_Core.web import (  # noqa: E402
    ai_assist as _bp_ai,
    build as _bp_build,
    cds_export as _bp_cds_export,
    download as _bp_download,
    examples_api as _bp_examples,
    io_api as _bp_io,
    logan as _bp_logan,
    mcp_console as _bp_mcp,
    meta as _bp_meta,
    pages as _bp_pages,
    prescreen as _bp_prescreen,
    refs as _bp_refs,
    samples as _bp_samples,
    settings as _bp_settings,
    submit as _bp_submit,
    tasks_api as _bp_tasks,
    tool_results as _bp_tool_results,
    tools_api as _bp_tools,
    virome as _bp_virome,
)

# 顺序与拆分前的 app.py 区块顺序一致（便于对照回溯）
#
# 2026-09-17：病毒浏览器（explorer）蓝图停止挂载 —— 平台不再提供 /explorer
# 页面与 /api/explorer/* 接口。
# 2026-09-27：**彻底归档** —— 引擎（Virus_Platform_Core/explorer/）、蓝图
# （web/explorer.py）、页面（explorer.html / app-explorer.js）、专属测试与
# 运行数据（run/explorer_exports、run/explorer_variation）全部移入
# archive/_retire_20260927/explorer/，代码与数据在本仓不再保留活动副本。
# 注意：app.py 的 JSON 清洗器原从 web/explorer.py 借 _nan_to_null（每个
# JSON 响应都走它）——现已内联为本文件顶部的 _nan_to_null，无外部依赖。
for _bp in (_bp_pages, _bp_ai, _bp_tool_results, _bp_refs, _bp_logan, _bp_io,
            _bp_mcp,
            _bp_samples, _bp_download, _bp_submit, _bp_meta, _bp_prescreen,
            _bp_settings, _bp_virome, _bp_tools, _bp_tasks, _bp_build,
            _bp_examples, _bp_cds_export):
    app.register_blueprint(_bp.bp)
del _bp


# 把不带 TaskManager 执行的子系统接进任务中心（统一列表/进度/停止/前往）：
# 下载批次（DownloadManager 自建线程）。执行仍归各自子系统，这里只是
# 注册 provider 回调做状态镜像；详见 web/download.py:bind_downloads。
from Virus_Platform_Core.web.download import bind_downloads as _bind_downloads  # noqa: E402
_bind_downloads()


# ------------------------------------------------------------------
# 兼容别名：app.py 拆分到 Virus_Platform_Core/web/ 之前，这些名字挂在 app 模块上，
# 仓库内测试（tests/_it_*.py、test_logan_trace.py）与外部脚本仍按老写法
# 引用（app.tm / app._KEEP_TOTAL / app._resolve_virus_db）。
# 这里只做转发，真实定义仍在各自模块，避免复制状态。
# ------------------------------------------------------------------
from Virus_Platform_Core.web.tasks import tm as tm  # noqa: E402,F401
from Virus_Platform_Core.web.tasks import _KEEP_TOTAL as _KEEP_TOTAL  # noqa: E402,F401
from Virus_Platform_Core.web.tools_api import _resolve_virus_db as _resolve_virus_db  # noqa: E402,F401


# ------------------------------------------------------------------
# 启动
# ------------------------------------------------------------------
def _pick_port():
    """选择可绑定的端口。

    Windows 上 Hyper-V/WSL 会动态保留端口段（netsh 可见），
    固定端口可能被系统排除导致 bind 被拒（访问权限不允许），
    因此逐个探测候选端口，全部失败再随机尝试。
    """
    import socket
    candidates = [8765, 8900, 8989, 9600, 8888, 5050, 5000]
    for p in candidates:
        try:
            with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
                s.bind(('127.0.0.1', p))
                return p
        except OSError:
            continue
    import random
    for _ in range(30):
        p = random.randint(10000, 19999)
        try:
            with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
                s.bind(('127.0.0.1', p))
                return p
        except OSError:
            continue
    raise RuntimeError('未找到可用的本地端口，请关闭占用端口的程序后重试')


def _open_browser(url):
    time.sleep(1.2)
    try:
        webbrowser.open(url)
    except Exception:
        pass


def _start_server(port):
    """在后台线程起 Flask，返回已绑定端口的 server 句柄。

    为什么不用 app.run：它会阻塞主线程，而 pywebview 的窗口必须开在主线程；
    换到后台线程后主线程才能建窗口，且拿到 server 句柄，窗口关闭时可停机。
    make_server 同步完成 bind，返回即意味着地址可访问，窗口首开不会白屏。
    """
    from werkzeug.serving import make_server
    srv = make_server('127.0.0.1', port, app, threaded=True)
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    return srv


def _gui_mode():
    """启动模式：'gui'（独立桌面窗口）| 'web'（系统浏览器）。

    优先级：命令行 `--web`/`--browser` 或 `--gui` > 环境变量 VP_GUI > 默认 gui。
    为什么要有命令行开关：原先只能用 `set VP_GUI=browser` 切模式，双击 .bat
    的人没有入口；打包后的 exe 也一样（`VirusPlatform.exe --web` 即可）。
    gui 模式下若 pywebview/WebView2 不可用，仍会自动回退浏览器（见 _open_window）。
    """
    for raw in sys.argv[1:]:
        a = raw.strip().lower()
        if a in ('--web', '--browser', '--no-gui'):
            return 'web'
        if a in ('--gui', '--window', '--webview'):
            return 'gui'
    v = os.environ.get('VP_GUI', '').strip().lower()
    if v in ('browser', 'web'):
        return 'web'
    if v in ('gui', 'window', 'webview'):
        return 'gui'
    return 'gui'


def _open_window(url, mode='gui'):
    """pywebview 独立桌面窗口（壳住本地页面，无浏览器地址栏）。

    返回 True=窗口已开过并关闭（进程应退出）；False=pywebview 不可用或
    指定了网页模式，调用方回退浏览器。Windows 固定 edgechromium 后端：
    系统缺 WebView2 运行时就走回退，绝不静默掉到老旧的 mshtml
    （IE 内核渲染会烂）。
    """
    if mode == 'web':
        return False
    try:
        import webview
    except Exception:
        print('  [提示] 未安装 pywebview，回退浏览器模式'
              '（pip install pywebview 后即为独立窗口）')
        return False
    try:
        webview.create_window('植物病毒分析平台', url,
                              width=1440, height=900, min_size=(1080, 700))
        webview.start(gui='edgechromium' if sys.platform == 'win32' else None)
        return True
    except Exception as e:
        print(f'  [提示] 桌面窗口启动失败（{type(e).__name__}: {e}），'
              f'回退浏览器模式')
        return False


def _pid_alive(pid):
    """pid 是否是一个仍在运行的进程（实例锁残留判定用）。"""
    try:
        pid = int(pid)
    except (TypeError, ValueError):
        return False
    if pid <= 0:
        return False
    if sys.platform == 'win32':
        import ctypes
        k32 = ctypes.windll.kernel32
        # 只 OpenProcess 成功不够（pid 被回收/僵尸也算存在）；
        # 再取退出码，等于 STILL_ACTIVE(259) 才算真活着。
        h = k32.OpenProcess(0x1000, False, pid)   # QUERY_LIMITED_INFORMATION
        if not h:
            return False
        try:
            code = ctypes.c_ulong()
            ok = k32.GetExitCodeProcess(h, ctypes.byref(code))
            return bool(ok) and code.value == 259
        finally:
            k32.CloseHandle(h)
    try:
        os.kill(pid, 0)
        return True
    except OSError:
        return False


def _acquire_instance_lock():
    """单实例锁：run/.platform.lock（O_CREAT|O_EXCL）。返回 (status, path)。

    为什么需要：多实例共享同一份 run/ 数据（queue.json 队列、tasks/*.json
    状态、results/ 结果目录），而 web/samples.py 的并发保护是**进程内**的
    set —— 跨进程完全失效，两个实例同时跑同一批样品会互相覆盖状态与结果。

    status: 'ok'=抢到锁（path 为锁文件）| 'busy'=已有实例在跑（调用方应退出）
            | 'skip'=锁不可用（磁盘/权限异常），放弃保护但不拦启动。
    锁内容为「pid 时间戳」：pid 仍存活 → 报错退出；pid 已消失（上次异常
    退出的残留锁）→ 覆盖接管。
    """
    from Virus_Platform_Core.config import PLATFORM_ROOT
    path = os.path.join(PLATFORM_ROOT, 'run', '.platform.lock')
    try:
        os.makedirs(os.path.dirname(path), exist_ok=True)
    except OSError as e:
        print(f'  [提示] 实例锁目录不可用（{e}），跳过并发保护。')
        return ('skip', None)
    for _ in range(2):
        try:
            fd = os.open(path, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
        except FileExistsError:
            old_pid = 0
            try:
                with open(path, encoding='utf-8') as f:
                    old_pid = int((f.read() or '0').split()[0])
            except (OSError, ValueError, IndexError):
                old_pid = 0
            if old_pid and _pid_alive(old_pid):
                print(f'  [错误] 检测到平台已在运行（pid={old_pid}）。')
                print('  同一份 run/ 数据不允许两个实例同时读写，请先关闭'
                      '已运行平台的窗口/控制台。')
                print(f'  确认该进程已不存在后，可删除锁文件再启动: {path}')
                _pause_before_exit()
                return ('busy', None)
            # 残留锁（上次异常退出，pid 已不存在）：覆盖接管
            print(f'  [提示] 发现残留锁（pid={old_pid or "未知"} 已不存在），'
                  '已接管。')
            try:
                os.remove(path)
            except OSError:
                pass
            continue
        except OSError as e:
            print(f'  [提示] 实例锁创建失败（{e}），跳过并发保护。')
            return ('skip', None)
        with os.fdopen(fd, 'w', encoding='utf-8') as f:
            f.write(f'{os.getpid()} {time.strftime("%Y-%m-%d %H:%M:%S")}\n')
        return ('ok', path)
    return ('skip', None)


def _release_instance_lock(path):
    """退出时删除自己写的锁（pid 不匹配则不动，避免误删他人锁）。"""
    if not path:
        return
    try:
        with open(path, encoding='utf-8') as f:
            if int((f.read() or '0').split()[0]) != os.getpid():
                return
        os.remove(path)
    except (OSError, ValueError):
        pass


def main():
    # 接管 SIGINT/SIGTERM/SIGBREAK：退出前把 AI 审计队列排空。Windows 下
    # SIGTERM **不触发 atexit**（实测），只靠 atexit 会丢最后几秒的记录。
    # 必须在主线程且尽早装（signal.signal 有线程限制，且装晚了会漏掉信号）。
    try:
        from Virus_Platform_Core.ai_agent import audit as _ai_audit
        _ai_audit.install_signals()
    except Exception:
        pass
    # 最早的启动动作：抢单实例锁（端口探测/起服务之前），防止第二个实例
    # 共享 run/ 目录把状态与结果互相覆盖。
    lock_status, lock_path = _acquire_instance_lock()
    if lock_status == 'busy':
        return
    mode = _gui_mode()
    try:
        port = _pick_port()
    except RuntimeError as e:
        print(f'  [错误] {e}')
        _release_instance_lock(lock_path if lock_status == 'ok' else None)
        _pause_before_exit()
        return
    # CSRF 端口校验的基准：起服务前写入，_guard_local_origin 按「完整源」
    # 比对 Origin/Referer（见 _origin_same_platform）。
    global _BOUND_PORT
    _BOUND_PORT = port
    # 实际端口落盘供拖拽代理等外部组件读取（此前代理只扫 8760-8789，
    # 平台备选端口 8900/9600/随机段全在扫描区间外——2026-10-04 审查 E4）
    try:
        import pathlib as _pl
        _port_file = _pl.Path(__file__).resolve().parent / 'run' / '.port'
        _port_file.parent.mkdir(parents=True, exist_ok=True)
        _port_file.write_text(str(port), encoding='ascii')
    except OSError:
        pass
    url = f'http://127.0.0.1:{port}'
    print('=' * 56)
    print('  植物病毒分析平台 启动中...'
          f'（模式: {"独立桌面窗口" if mode == "gui" else "系统浏览器"}）')
    print(f'  本机地址: {url}')
    if port != 8765:
        print(f'  （默认端口 8765 被系统占用/保留，已自动改用 {port}）')
    print('  关闭窗口/控制台即退出平台'
          if mode == 'gui' else
          '  平台在本控制台窗口内运行，关闭它即退出')
    print('=' * 56)
    try:
        srv = _start_server(port)
    except OSError as e:
        print(f'\n  [错误] 端口 {port} 监听失败: {e}')
        print('  可在管理员 PowerShell 运行以下命令查看被系统保留的端口段:')
        print('     netsh interface ipv4 show excludedportrange protocol=tcp')
        _release_instance_lock(lock_path if lock_status == 'ok' else None)
        _pause_before_exit()
        return
    try:
        if _open_window(url, mode):
            print('  桌面窗口已关闭，平台退出。')
            return
        # 浏览器模式：显式 --web / VP_GUI=browser，或 gui 模式但窗口不可用时回退
        if mode == 'gui':
            print('  （已回退浏览器模式；如需固定用网页模式，'
                  '下次用 启动平台-网页.bat 或 python app.py --web）')
        threading.Thread(target=_open_browser, args=(url,), daemon=True).start()
        while True:
            time.sleep(3600)    # 关闭本控制台窗口即退出（守护线程随之终止）
    except KeyboardInterrupt:
        print('\n  收到退出信号，平台关闭。')
    finally:
        try:
            srv.shutdown()
        except Exception:
            pass
        _release_instance_lock(lock_path if lock_status == 'ok' else None)


if __name__ == '__main__':
    main()
