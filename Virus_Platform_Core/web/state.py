# -*- coding: utf-8 -*-
"""Web 层共享状态与路径解析——app.py 拆分的「单一真源」。

为什么需要它：app.py 里 `cfg` / `_WWW` / `_tool_runs_root` 等模块级对象
被上百个路由共用。一旦把路由搬进独立模块，若各模块反过来 `import app`
取这些对象，就会形成 app → blueprint → app 的循环导入。因此统一收到
本模块，由 app.py 与各 blueprint 双向导入。

内容：
  - cfg         配置单例（Virus_Platform_Core.config.get_config() 本身即单例）
  - _WWW        webapp 目录（模板/静态资源根）
  - webapp_dir() 目录探测（源码运行 / PyInstaller 冻结分发两种布局）
  - tool_runs_root()  tool_runs 根（跟随自定义输出根）
"""
import json
import os
import sys

from Virus_Platform_Core.config import DIRS, PLATFORM_ROOT, get_config


def webapp_dir():
    """定位 webapp 目录。

    源码运行：平台根/webapp（app.py 与 Virus_Platform_Core/ 同级）；
    PyInstaller 冻结分发：优先内置资源 `sys._MEIPASS/webapp`。
    """
    # Virus_Platform_Core/web/state.py -> Virus_Platform_Core/web -> 平台根
    root = os.path.dirname(os.path.dirname(os.path.dirname(
        os.path.abspath(__file__))))
    src = os.path.join(root, 'webapp')
    if os.path.isdir(src):
        return src
    if getattr(sys, 'frozen', False):
        base = getattr(sys, '_MEIPASS',
                       os.path.dirname(os.path.abspath(sys.executable)))
        bundled = os.path.join(base, 'webapp')
        if os.path.isdir(bundled):
            return bundled
    return src


_WWW = webapp_dir()

# 配置单例：Virus_Platform_Core.config.get_config() 内部缓存，重复调用返回同一对象，
# 因此 app.py 与各 blueprint 拿到的是同一个 cfg（改语言/参数全局生效）。
cfg = get_config()


def tool_runs_root():
    """tool_runs 根目录（设置自定义输出根后由 DIRS 决定）。"""
    return DIRS.get('tool_runs') or os.path.join(PLATFORM_ROOT, 'run', 'tool_runs')


# 自定义输出位置的指针文件名：tool_runs/<run>/out_link.json 记录产物实际
# 所在目录（run_dir 外移后列表/下载/打开仍可寻址）。
RUN_OUT_LINK = 'out_link.json'


def run_real_dir(run):
    """运行的真实产物目录：tool_runs/<run>/out_link.json 存在且登记了
    out_dir 时返回外移目录，否则 tool_runs/<run> 本身。

    为什么需要指针：工具产物下载/历史列表/打开目录全都从 tool_runs 根拼
    路径，run_dir 外移到用户目录后没有这份登记就全部 404/失联。
    """
    from Virus_Platform_Core.utils import check_path
    base = os.path.join(tool_runs_root(), run)
    try:
        with open(os.path.join(base, RUN_OUT_LINK), encoding='utf-8') as f:
            out = (json.load(f) or {}).get('out_dir')
        if out:
            return check_path(out, must_exist=True, in_platform=True)
    except (OSError, ValueError):
        pass
    return base


# ------------------------------------------------------------------
# 批处理队列单例（SampleQueue 在 app.py 中实例化，2D 批次迁出）
# 用「注册 / 取用」而非 `import app`，避免 blueprint ↔ app 循环导入。
# ------------------------------------------------------------------
_sample_queue = None


def set_sample_queue(queue):
    """由 app.py 在创建 sample_queue 后调用。"""
    global _sample_queue
    _sample_queue = queue


def get_sample_queue():
    """供 blueprint 取用（download 的「转入分析流程」用）。"""
    return _sample_queue
