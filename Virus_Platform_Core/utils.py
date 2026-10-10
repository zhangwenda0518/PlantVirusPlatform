# -*- coding: utf-8 -*-
"""
通用工具：路径安全校验、命令执行、日志、FASTA/FASTQ 流式处理、步骤断点标记。
所有 subprocess 调用一律使用参数列表（shell=False），不拼接命令字符串。
所有文件访问统一经由 safe_open()/safe_remove()：先 check_path 规范化校验
（拒绝 .. 段、写模式限定平台根内、读模式要求存在），再执行 I/O。
"""
import math
import os
import re
import sys
import gzip
import time
import threading
import subprocess

if sys.platform == 'win32':
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
    sys.stderr.reconfigure(encoding='utf-8', errors='replace')


def nan_to_null(o):
    """把 NaN / ±Inf 递归换成 None（JSON 落盘 / 出口统一清洗）。

    Python 的 `json.dumps` 默认 `allow_nan=True`，会把 NaN 原样吐成裸字面量
    `NaN` —— 非法 JSON：Python 端 json.loads 能容忍（自测假通过），浏览器
    `JSON.parse` 直接报 `Unexpected token 'N'`、整份响应被丢弃。来源是
    pandas 聚合的均值/占比（空集合求均值 → NaN）。原实现随病毒浏览器
    归档（archive/_retire_20260927），现收编为本模块的公共工具。
    """
    if isinstance(o, float):          # np.float64 也是 float 子类，一并覆盖
        return o if math.isfinite(o) else None
    if isinstance(o, dict):
        return {k: nan_to_null(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [nan_to_null(v) for v in o]
    # plotly 的 fig.data 里 x/y/z 都是 numpy.ndarray；先 tolist 再递归清洗里面的
    # NaN。float 分支必须在前，否则 np.float64 会先走 tolist 丢掉标量语义。
    if hasattr(o, 'tolist'):
        return nan_to_null(o.tolist())
    if hasattr(o, 'item'):            # numpy 标量（np.int64 等不是 int 子类）
        return nan_to_null(o.item())
    return o



def decode_output(raw):
    """子进程输出解码：优先 UTF-8，失败退 GBK，最后才 replace。

    Windows 上不少工具（bcftools/minimap2/samtools 等）按系统 ANSI 代码页
    (936/GBK) 打印命令行，其中含中文路径。若直接用 encoding='utf-8',
    errors='replace' 解码，每字节变 U+FFFD 且不可逆。先按字节读入再逐
    编码试解，可无损保留中文。
    """
    if isinstance(raw, str):
        return raw
    for enc in ('utf-8', 'gbk'):
        try:
            return raw.decode(enc)
        except UnicodeDecodeError:
            continue
    return raw.decode('utf-8', errors='replace')


def to_utf8_file(path):
    """文本产物就地规范为 UTF-8（外部工具自写日志/报告专用）。

    文件整体已是合法 UTF-8（含纯 ASCII）→ 原样保留、返回 False；
    否则**按行**试解 utf-8 → GBK（decode_output 口径）后整体重写为
    UTF-8，返回 True。按行而不是整 buffer：混编码文件（如 cmd 重定向
    里 UTF-8 框线字符与 GBK 中文报错共存）里合法的 UTF-8 行逐字节保留、
    只有 GBK 行被转写，两头都不损。

    典型场景：RAxML-NG / LSD2 这类 C++ 工具把含中文路径的行按系统 ANSI
    代码页（中文 Windows = GBK）写进自产日志，与平台其余产物的编码
    不一致；跑完立即调用本函数即完成统一。IO/解码失败静默跳过。
    """
    try:
        with safe_open(path, 'rb') as f:
            raw = f.read()
    except OSError:
        return False
    try:
        raw.decode('utf-8')
        return False
    except UnicodeDecodeError:
        pass
    out, changed = [], False
    for line in raw.splitlines(keepends=True):
        try:
            line.decode('utf-8')
            out.append(line)
        except UnicodeDecodeError:
            changed = True
            out.append(decode_output(line).encode('utf-8', errors='replace'))
    if not changed:
        return False
    try:
        with safe_open(path, 'wb') as f:
            f.write(b''.join(out))
        return True
    except OSError:
        return False


# ------------------------------------------------------------------
# 路径安全
# ------------------------------------------------------------------
def _redirect_legacy_top(raw):
    """目录整理（2026-09-10）兼容：把旧的顶层目录重定向到新位置。

    运行期数据已搬进 run/、外部二进制与第三方搬进 3rd/、示例搬到
    <平台根>/examples/，但历史数据、前端拼接的路径、文档示例里仍可能写
    tool_runs/xxx、results/xxx、databases/examples/xxx… 这类旧路径。

    两类写法都要覆盖（很多 API 是**先 os.path.join(PLATFORM_ROOT, path)
    再 check_path**，只处理相对路径会漏掉它们）：
      1) 相对路径：tool_runs/x  → run/tool_runs/x
      2) 平台根下的绝对路径：<根>/tool_runs/x → <根>/run/tool_runs/x
    平台根之外的绝对路径原样返回。

    仅当旧位置不存在、新位置存在时才改写，避免误伤同名目录。
    """
    s = str(raw).replace('\\', '/')
    if not s:
        return raw
    try:
        from .config import LEGACY_TOP_DIRS, PLATFORM_ROOT
    except Exception:                                # 循环导入兜底
        return raw
    root_s = str(PLATFORM_ROOT).replace('\\', '/').rstrip('/')
    low_s, low_root = s.lower(), root_s.lower()
    if low_s.startswith(low_root + '/'):
        tail = s[len(root_s) + 1:]                   # 平台根下的绝对路径
        absolute = True
    elif os.path.isabs(str(raw)) or ':' in s.split('/')[0]:
        return raw                                   # 平台外的绝对路径 → 不动
    else:
        tail = s
        absolute = False

    def _finish(new_rel):
        """新相对路径 → 按原式（绝对/相对）重建。"""
        return (root_s + '/' + new_rel) if absolute else new_rel.replace('/', os.sep)

    def _usable(old_rel, new_rel):
        old_abs = os.path.join(PLATFORM_ROOT, *old_rel.split('/'))
        new_abs = os.path.join(PLATFORM_ROOT, *new_rel.split('/'))
        return (not os.path.exists(old_abs)) and os.path.exists(new_abs)

    # 多段前缀的旧→新映射（表见 Virus_Platform_Core/config.py 的 LEGACY_MULTI_DIRS）。
    # 例：示例数据 2026-09-10 由 databases/examples/ 迁到 <平台根>/examples/；
    # 用户数据 databases/misc/gb/、databases/ncbi_refs/ 迁入 run/。
    try:
        from .config import LEGACY_MULTI_DIRS
    except Exception:                                # 循环导入兜底
        LEGACY_MULTI_DIRS = {}
    for old_pre, new_pre in LEGACY_MULTI_DIRS.items():
        if not (tail == old_pre or tail.startswith(old_pre + '/')):
            continue
        rest = tail[len(old_pre):].lstrip('/')
        # 按**完整路径**判定（而非只看顶层目录是否还在）：根目录可能残留
        # 旧顶层空壳（tool_runs/、results/ 里只剩 _archive 等杂项），顶层级
        # 判定会让重映射被永久抑制、所有旧相对路径全部失效（2026-09-12 实测，
        # /api/align/file、submit /fasta 与前端 tool_runs/... 引用全部中招）。
        old_full = old_pre + ('/' + rest if rest else '')
        new_full = new_pre + ('/' + rest if rest else '')
        if not _usable(old_full, new_full):
            return raw
        return _finish(new_full)

    head, sep, tail_rest = tail.partition('/')
    if not sep:
        return raw
    new = LEGACY_TOP_DIRS.get(head)
    if not new or not _usable(tail, new + '/' + tail_rest):
        return raw                                   # 旧位置还在 / 新位置不存在 → 不动
    return _finish(new + '/' + tail_rest)


def check_path(path, must_exist=False, in_platform=False):
    """规范化并校验路径。

    - 拒绝路径中含 '..' 段（防穿越）
    - in_platform=True 时强制路径位于平台根目录内（用于输出/中间文件）
    - must_exist=True 时要求文件/目录已存在
    - 旧顶层目录（tool_runs/ results/ …）自动重定向到 run/ 3rd/（见
      _redirect_legacy_top）
    返回规范化后的绝对路径。
    """
    raw = _redirect_legacy_top(path)
    for seg in str(raw).replace('/', os.sep).split(os.sep):
        if seg == '..':
            raise ValueError(f"非法路径（含 ..）: {raw}")
    p = os.path.normpath(os.path.abspath(raw))
    if in_platform:
        # 平台根 + 自定义输出根（若有）均视为合法写入区
        from .config import write_roots
        for root in write_roots():
            root = os.path.normpath(root)
            if p == root or p.startswith(root + os.sep):
                break
        else:
            raise ValueError(f"路径越界（仅允许平台或输出目录内）: {raw}")
    if must_exist and not os.path.exists(p):
        raise FileNotFoundError(f"路径不存在: {p}")
    return p


_STEP_RE = re.compile(r'^[A-Za-z0-9_\-]+$')


def check_step(step):
    """步骤名白名单校验：仅字母/数字/下划线/连字符，防止拼入文件名做路径注入。"""
    s = str(step)
    if not _STEP_RE.match(s):
        raise ValueError(f"非法步骤名: {s!r}")
    return s


# ---- 统一文件 I/O 入口（open/os.remove 仅出现在以下两个函数内） ----
def safe_open(path, mode='rt', compress_level=6):
    """安全打开文件。

    读模式（mode 含 r）：路径必须已存在；
    写模式（mode 含 w/a/x）：路径必须位于平台根目录内，父目录自动创建。
    自动识别 .gz（读与写均支持；写 gz 需显式以 .gz 结尾）。
    .gz 文件在 crabz.exe 可用时改走 crabz 多线程管道（比 Python gzip
    快数倍，写默认级别 6 而非 gzip 模块的 9），失败自动回退。
    compress_level 仅写 .gz 时生效（crabz -l；回退 gzip.compresslevel）。
    """
    is_read = 'r' in mode
    is_write = any(c in mode for c in ('w', 'a', 'x'))
    is_bin = 'b' in mode
    if is_read and not is_write:
        p = check_path(path, must_exist=True)
        if p.endswith('.gz'):
            cz = crabz_path()
            if cz:
                try:
                    return _GzReadPipe(p, mode)
                except OSError:
                    pass     # crabz 启动失败 → 回退 gzip
            if is_bin:
                return gzip.open(p, mode)
            return gzip.open(p, mode, encoding='utf-8', errors='replace',
                             newline=None)
        if is_bin:
            return open(p, mode)
        return open(p, mode, encoding='utf-8', errors='replace', newline=None)
    # 写模式：限定平台内（追加模式 gz 不支持 crabz 管道，回退 gzip）
    p = check_path(path, in_platform=True)
    os.makedirs(os.path.dirname(p), exist_ok=True)
    if p.endswith('.gz'):
        cz = crabz_path()
        if cz and 'a' not in mode:
            try:
                return _GzWritePipe(p, mode, level=compress_level)
            except OSError:
                pass
        if is_bin:
            return gzip.open(p, mode, compresslevel=compress_level)
        # fastq/fasta 一律 LF 换行（CRLF 会导致 SPAdes 等解析失败）
        return gzip.open(p, mode, encoding='utf-8', errors='replace', newline='\n',
                         compresslevel=compress_level)
    if is_bin:
        return open(p, mode)
    return open(p, mode, encoding='utf-8', errors='replace', newline='\n')


def safe_remove(path, in_platform=True):
    """安全删除文件（默认限定平台内）。"""
    p = check_path(path, must_exist=True, in_platform=in_platform)
    os.remove(p)


# 兼容别名（历史调用点）
def open_maybe_gzip(path, mode='rt'):
    return safe_open(path, mode)


def open_write(path):
    return safe_open(path, 'wt')


def open_write_bin(path):
    return safe_open(path, 'wb')


# ------------------------------------------------------------------
# 日志
# ------------------------------------------------------------------
class TaskLogger:
    """线程安全日志器：写 UTF-8 文件 + 推送到回调（GUI 实时显示）。"""

    def __init__(self, log_file=None, callback=None, echo=False):
        self.log_file = log_file
        self.callback = callback          # callable(str) or None
        self.echo = echo
        self._lock = threading.Lock()
        self._fh = None
        if log_file:
            self._fh = safe_open(log_file, 'at')

    def log(self, msg, level="INFO"):
        line = f"[{time.strftime('%H:%M:%S')}] [{level}] {msg}"
        with self._lock:
            if self._fh:
                self._fh.write(line + "\n")
                self._fh.flush()
            if self.echo:
                try:
                    print(line)
                except Exception:
                    pass
        if self.callback:
            try:
                self.callback(line)
            except Exception:
                pass

    def close(self):
        with self._lock:
            if self._fh:
                self._fh.close()
                self._fh = None

    # 上下文管理器 + 析构兜底（2026-10-04 审查 E5）：CLI 命令中途抛异常时
    # 旧模式 `logger=make_logger(); …; logger.close()` 的 close 走不到，
    # 日志尾部丢失、句柄延迟释放。close 幂等，重复调用安全。
    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        self.close()
        return False

    def __del__(self):
        try:
            self.close()
        except Exception:
            pass


# ------------------------------------------------------------------
# crabz 加速的 gzip 流（可选；未装 crabz 时回退 Python gzip）
# ------------------------------------------------------------------
_CRABZ_PATH = None
_CRABZ_CHECKED = False


def crabz_path():
    """crabz.exe 路径（未安装返回 None）。Rust 多线程 gzip，替代 Python gzip。"""
    global _CRABZ_PATH, _CRABZ_CHECKED
    if not _CRABZ_CHECKED:
        _CRABZ_CHECKED = True
        try:
            from .config import get_config
            try:
                _CRABZ_PATH = get_config().tool('crabz')
            except (FileNotFoundError, RuntimeError):
                _CRABZ_PATH = None
        except Exception:
            _CRABZ_PATH = None
    return _CRABZ_PATH


def _crabz_threads():
    try:
        from .config import get_config
        return max(2, min(4, get_config().threads))
    except Exception:
        return 4


def _popen_hidden(cmd):
    flags = 0
    if sys.platform == 'win32':
        flags = getattr(subprocess, 'CREATE_NO_WINDOW', 0)
    return subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
                            stdin=subprocess.PIPE, creationflags=flags)


class _GzReadPipe:
    """crabz -d 解压流。用法与 gzip.open(path,'rb'/'rt') 等价的上下文对象。

    safe_open 的调用方直接使用返回对象（不经 with），故代理内部文件
    对象的读写方法（readline/read/迭代/...），close() 附加进程回收。
    """

    def __init__(self, path, mode):
        proc = _popen_hidden([crabz_path(), '-d', '-Q', '-p', '2', str(path)])
        self._proc = proc
        if 't' in mode:
            import io
            self._fh = io.TextIOWrapper(
                proc.stdout, encoding='utf-8', errors='replace', newline=None)
        else:
            self._fh = proc.stdout
        self._closed = False

    def __enter__(self):
        return self._fh

    def __exit__(self, *exc):
        self.close()

    def __iter__(self):
        return iter(self._fh)

    def __getattr__(self, name):
        return getattr(self._fh, name)

    def close(self):
        if self._closed:
            return
        self._closed = True
        try:
            self._fh.close()
        except OSError:
            pass
        try:
            rc = self._proc.wait(timeout=60)
        except subprocess.TimeoutExpired:
            # 网络盘/杀毒扫描可能让 crabz 收尾卡死：杀掉并给出可诊断错误，
            # 而不是向上抛裸 TimeoutExpired 让调用方莫名其妙崩溃
            self._proc.kill()
            self._proc.wait()
            raise RuntimeError('crabz 解压进程关闭超时（60s），已强制终止')
        if rc not in (0, None):
            raise RuntimeError(f"crabz 解压失败 (退出码 {rc})")


class _GzWritePipe:
    """crabz 压缩流（stdin 进，文件出）。等价 gzip.open(path,'wb'/'wt')。"""

    def __init__(self, path, mode, level=6):
        p = check_path(path, in_platform=True)
        os.makedirs(os.path.dirname(p), exist_ok=True)
        proc = _popen_hidden([crabz_path(), '-Q', '-l', str(level),
                              '-p', str(_crabz_threads()), '-o', p])
        self._proc = proc
        self._path = p
        if 't' in mode:
            import io
            self._fh = io.TextIOWrapper(
                proc.stdin, encoding='utf-8', errors='replace', newline='\n')
        else:
            self._fh = proc.stdin
        self._closed = False

    def __enter__(self):
        return self._fh

    def __exit__(self, *exc):
        self.close()

    def __iter__(self):
        return iter(self._fh)

    def __getattr__(self, name):
        return getattr(self._fh, name)

    def close(self):
        if self._closed:
            return
        self._closed = True
        try:
            self._fh.close()
        except (OSError, ValueError):
            pass
        try:
            rc = self._proc.wait(timeout=600)
        except subprocess.TimeoutExpired:
            self._proc.kill()
            self._proc.wait()
            try:
                os.remove(self._path)
            except OSError:
                pass
            raise RuntimeError('crabz 压缩进程收尾超时（600s），已强制终止')
        if rc not in (0, None):
            try:
                os.remove(self._path)
            except OSError:
                pass
            raise RuntimeError(f"crabz 压缩失败 (退出码 {rc})")


def dir_size(path):
    """目录总字节数（元数据遍历，用于实测 chunk 磁盘占用）。"""
    total = 0
    for dp, _dn, fns in os.walk(path):
        for fn in fns:
            try:
                total += os.path.getsize(os.path.join(dp, fn))
            except OSError:
                pass
    return total


def db_trace(paths, version_files=(), n_max=40):
    """数据库可追溯信息（写进 db_info.json / 阶段 summary.json 用）。

    paths: 库文件或库目录序列——目录取第一层文件（n_max 上限，防 mmseqs
           分片类大目录把 summary 撑爆）；不存在的项静默跳过。
    version_files: 版本锚点文件（如 bfvd.version），读首行内容。
    返回 (db_files, db_versions)：
      db_files    {文件名: [size_bytes, mtime 'YYYY-MM-DD HH:MM:SS']}
      db_versions {文件名: 首行内容}
    """
    import time as _time
    files = {}
    for p in paths:
        cands = []
        if os.path.isfile(str(p)):
            cands = [str(p)]
        elif os.path.isdir(str(p)):
            try:
                cands = [os.path.join(str(p), n)
                         for n in sorted(os.listdir(str(p)))[:n_max]]
            except OSError:
                continue
        for fp in cands:
            if not os.path.isfile(fp):
                continue
            try:
                st = os.stat(fp)
            except OSError:
                continue
            files[os.path.basename(fp)] = [
                st.st_size,
                _time.strftime('%Y-%m-%d %H:%M:%S',
                               _time.localtime(st.st_mtime))]
    versions = {}
    for p in version_files:
        if not os.path.isfile(str(p)):
            continue
        try:
            with safe_open(str(p)) as f:
                versions[os.path.basename(str(p))] = \
                    (f.readline() or '').strip()
        except (OSError, ValueError):
            pass
    return files, versions


def est_decompressed(paths):
    """估算一组（可能 .gz）输入解压后的总字节数。
    FASTQ 压缩比经验值 ~5.5:1，FASTA ~4:1（无质量行），普通文件 1:1。"""
    total = 0
    for p in paths:
        try:
            sz = os.path.getsize(str(p))
        except OSError:
            continue
        low = str(p).lower()
        if low.endswith(('.fastq.gz', '.fq.gz')):
            total += sz * 5.5
        elif low.endswith('.gz'):
            total += sz * 4.0
        else:
            total += sz
    return total


def log_res_plan(logger, title, threads=None, mem_gb=None, disk_gb=None,
                 note=''):
    """阶段资源预估日志（统一格式，便于核对各阶段内存/线程/磁盘）。"""
    parts = [f"线程 {threads if threads is not None else '自动'}"]
    if mem_gb is not None:
        parts.append(f"内存 ~{mem_gb:.1f} GB")
    if disk_gb is not None:
        parts.append(f"磁盘 ~{disk_gb:.1f} GB")
    if note:
        parts.append(note)
    if logger:
        logger.log(f"📊 资源预估 [{title}] " + "；".join(parts), "PLAN")


def fmt_eta(seconds):
    """秒 → 人读时长（2小时3分 / 5分10秒 / 40秒）。"""
    if seconds is None or seconds < 0:
        return ''
    s = int(seconds)
    if s >= 3600:
        return f"{s // 3600}小时{(s % 3600) // 60:02d}分"
    if s >= 60:
        return f"{s // 60}分{s % 60:02d}秒"
    return f"{s}秒"


# ------------------------------------------------------------------
# 命令执行
# ------------------------------------------------------------------
# 外部命令默认超时（秒）：任何 run_cmd/run_cmd_redirect 未显式传 timeout 时
# 生效。背景：SPAdes 死锁/杀毒软件锁文件曾让任务线程永久挂起，heavy 任务槽
# （默认 2）被占满后全平台分析排队瘫痪，只能重启。默认 6 小时——本平台输入
# 经子采样，正常阶段远低于此值；确有超长任务时用环境变量调整：
#   set VP_CMD_TIMEOUT_SEC=43200   （12 小时；0 = 关闭默认超时）
DEFAULT_CMD_TIMEOUT_SEC = int(os.environ.get('VP_CMD_TIMEOUT_SEC', '21600') or 0)


def total_ram():
    """物理内存总量（字节）；失败返回 0。"""
    import ctypes
    import sys as _sys
    if not _sys.platform.startswith('win'):
        try:
            with open('/proc/meminfo') as f:
                for line in f:
                    if line.startswith('MemTotal:'):
                        return int(line.split()[1]) * 1024
        except OSError:
            return 0
        return 0

    class MEMORYSTATUSEX(ctypes.Structure):
        _fields_ = [('dwLength', ctypes.c_ulong),
                    ('dwMemoryLoad', ctypes.c_ulong),
                    ('ullTotalPhys', ctypes.c_ulonglong),
                    ('ullAvailPhys', ctypes.c_ulonglong),
                    ('ullTotalPageFile', ctypes.c_ulonglong),
                    ('ullAvailPageFile', ctypes.c_ulonglong),
                    ('ullTotalVirtual', ctypes.c_ulonglong),
                    ('ullAvailVirtual', ctypes.c_ulonglong),
                    ('ullAvailExtendedVirtual', ctypes.c_ulonglong)]

    stat = MEMORYSTATUSEX()
    stat.dwLength = ctypes.sizeof(MEMORYSTATUSEX)
    if ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(stat)):
        return int(stat.ullTotalPhys)
    return 0


# ------------------------------------------------------------------
# 任务级取消上下文：TaskManager 在每个任务线程绑定 cancel 事件与子进程
# 注册表；run_cmd / run_cmd_redirect 自动注册所启动的进程并在取消时
# 由 TaskManager 杀进程树（硬停止）。非任务线程（CLI 等）绑定为 None，
# 行为与从前完全一致。
# ------------------------------------------------------------------
_TASK_STATE = threading.local()


def task_bind(cancel=None, procs=None):
    """绑定当前线程（任务）的取消事件与子进程注册表；不传即解绑。"""
    _TASK_STATE.cancel = cancel
    _TASK_STATE.procs = procs


def task_check_cancel():
    """任务取消事件置位时抛 RuntimeError（长循环内可随时调用）。"""
    ev = getattr(_TASK_STATE, 'cancel', None)
    if ev is not None and ev.is_set():
        raise RuntimeError('任务已停止（用户取消）')


def task_cancel_event():
    """返回当前线程（任务）绑定的取消事件（未绑定 → None）。

    _TASK_STATE 是 threading.local()：**子线程读不到调用线程的事件**。
    子线程看门狗（如 suvtk / LOGAN 子进程静默挂起时的轮询）必须先在任务
    线程里取出该事件引用再传进线程，否则取消永远不被感知
    （2026-10-03 修复 P3-18 踩坑）。
    """
    return getattr(_TASK_STATE, 'cancel', None)


def _kill_proc_tree(proc):
    """杀整棵进程树（Windows: taskkill /F /T；其它平台退化为 proc.kill()）。

    2026-10-04 审查 E3：run_cmd 的超时/异常路径此前只 proc.kill() 直接
    子进程——spades.bat→spades.py、kv 引擎→salmon 这类孙进程会成孤儿
    继续吃 CPU/内存，而任务闸门名额已释放。与 web/tasks._kill_tree
    （取消路径）统一为进程树口径。"""
    if proc.poll() is not None:
        return
    try:
        if sys.platform == 'win32':
            subprocess.run(['taskkill', '/F', '/T', '/PID', str(proc.pid)],
                           capture_output=True, timeout=15)
        else:
            proc.kill()
    except Exception:
        try:
            proc.kill()
        except Exception:
            pass


def _task_register(proc):
    procs = getattr(_TASK_STATE, 'procs', None)
    if procs is not None:
        procs.append(proc)
    # 竞态泄漏：Popen 到本函数之间有一个微秒级窗口，若 cancel 恰好在此期间
    # 置位，TaskManager.cancel 的 procs 快照看不到它，该进程会成为孤儿
    # （任务已结束但进程仍在跑，需手工 taskkill）。注册后立即复查一次。
    ev = getattr(_TASK_STATE, 'cancel', None)
    if ev is not None and ev.is_set():
        _kill_proc_tree(proc)


def task_register_proc(proc):
    """公开入口：run_cmd 之外自行 Popen 的代码也把子进程登记进任务，
    供 TaskManager.cancel 硬停止。非任务线程内调用是空操作。"""
    _task_register(proc)


def run_cmd(cmd, logger=None, timeout=None, cwd=None, check=True,
            env=None, silent=False, monitor_fn=None, monitor_interval=5.0):
    """执行外部命令（参数列表形式）。stdout/stderr 实时写入日志。

    monitor_fn: 可选回调，命令运行期间每 monitor_interval 秒调用一次
    （用于长命令的进度看门狗，如 chunk 目录增长、spades 日志解析）。
    返回 returncode；check=True 时非零退出抛 RuntimeError。
    timeout: 秒。缺省用 DEFAULT_CMD_TIMEOUT_SEC（外部工具挂死时不再永久
    占用 heavy 任务槽）；显式传 None 之外的值可覆盖，传 0 关闭超时。
    """
    cmd = [str(c) for c in cmd]
    cwd = check_path(cwd, must_exist=True) if cwd else None
    if logger:
        logger.log("$ " + subprocess.list2cmdline(cmd), "CMD")
    task_check_cancel()
    creationflags = 0
    if sys.platform == 'win32':
        creationflags = getattr(subprocess, 'CREATE_NO_WINDOW', 0)
    if timeout is None:
        timeout = DEFAULT_CMD_TIMEOUT_SEC or None
    elif timeout <= 0:
        timeout = None

    proc = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        text=False,
        cwd=cwd, env=env, creationflags=creationflags)
    _task_register(proc)

    stop_mon = threading.Event()
    mon_thread = None
    if monitor_fn:
        def _watch():
            while not stop_mon.wait(monitor_interval):
                try:
                    monitor_fn()
                except Exception:
                    pass
        mon_thread = threading.Thread(target=_watch, daemon=True)
        mon_thread.start()

    # timeout 看门狗：阻塞读 stdout 时 wait(timeout) 永远轮不到，
    # 由看门狗线程在超时后强杀进程，读循环随即拿到 EOF。
    # timed_out 标志让报错说清「是超时被杀」而非误导性的「退出码 <0」。
    timed_out = threading.Event()
    if timeout:
        def _kill_watchdog():
            if not stop_mon.wait(timeout):
                timed_out.set()
                _kill_proc_tree(proc)
        _wd = threading.Thread(target=_kill_watchdog, daemon=True)
        _wd.start()

    t0 = time.time()
    try:
        for raw_line in proc.stdout:
            # 按字节读入后 lossless 解码：工具若按 GBK 打印中文路径可完整保留
            line = decode_output(raw_line).rstrip()
            if line and not silent:
                if logger:
                    logger.log(line, "OUT")
        proc.wait(timeout=timeout)
    except subprocess.TimeoutExpired:
        _kill_proc_tree(proc)
        proc.wait()
        raise RuntimeError(f"命令超时被终止: {cmd[0]}")
    except KeyboardInterrupt:
        _kill_proc_tree(proc)
        proc.wait()
        raise
    finally:
        if mon_thread:
            stop_mon.set()
            mon_thread.join(timeout=2)

    if timed_out.is_set():
        raise RuntimeError(
            f"命令超时（运行 {timeout}s 后被强制终止）: {cmd[0]}"
            f"（可设环境变量 VP_CMD_TIMEOUT_SEC 调大，0=不限制）")
    rc = proc.returncode
    dt = time.time() - t0
    if logger:
        logger.log(f"[退出码 {rc}] 耗时 {dt:.1f}s", "OUT")
    task_check_cancel()
    if check and rc != 0:
        raise RuntimeError(f"命令失败 (退出码 {rc}): {subprocess.list2cmdline(cmd)}")
    return rc


def run_cmd_redirect(cmd, out_path, logger=None, timeout=None, cwd=None,
                     append=False):
    """执行外部命令并将 stdout 重定向到平台内文件（mafft/fasttree 等输出流）。

    append=True 时追加写入（kvsuite 引擎多个子命令共用一份日志）。
    返回 returncode；非零抛 RuntimeError。

    子进程环境强制 PYTHONUTF8=1 / PYTHONIOENCODING=utf-8：Python 子进程
    （如 kvsuite 引擎）把中文日志 print 到重定向文件时默认走系统 ANSI
    代码页（中文 Windows = GBK），会让同一份日志里 UTF-8 与 GBK 混存；
    声明 UTF-8 后落盘编码与平台自身日志一致。非 Python 工具忽略该变量。
    """
    cmd = [str(c) for c in cmd]
    out = check_path(out_path, must_exist=False, in_platform=True)
    os.makedirs(os.path.dirname(out), exist_ok=True)
    cwd = check_path(cwd, must_exist=True) if cwd else None
    if logger:
        logger.log(f"$ {subprocess.list2cmdline(cmd)} > {os.path.basename(out)}", "CMD")
    task_check_cancel()
    creationflags = 0
    if sys.platform == 'win32':
        creationflags = getattr(subprocess, 'CREATE_NO_WINDOW', 0)
    if timeout is None:
        timeout = DEFAULT_CMD_TIMEOUT_SEC or None
    elif timeout <= 0:
        timeout = None
    t0 = time.time()
    child_env = os.environ.copy()
    child_env['PYTHONUTF8'] = '1'
    child_env['PYTHONIOENCODING'] = 'utf-8'
    with safe_open(out, 'ab' if append else 'wb') as fo:
        proc = subprocess.Popen(
            cmd, stdout=fo, stderr=subprocess.PIPE,
            cwd=cwd, creationflags=creationflags, env=child_env)
        _task_register(proc)
        err_lines = []
        timed_out = threading.Event()
        stop_wd = threading.Event()
        if timeout:
            def _kill_watchdog():
                if not stop_wd.wait(timeout):
                    timed_out.set()
                    try:
                        proc.kill()
                    except Exception:
                        pass
            threading.Thread(target=_kill_watchdog, daemon=True).start()
        try:
            for line in iter(proc.stderr.readline, b''):
                # 与 run_cmd 同一解码口径：GBK stderr 中的中文路径不丢
                txt = decode_output(line).rstrip()
                if txt:
                    err_lines.append(txt)
                    if logger:
                        logger.log(txt, "OUT")
            proc.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            proc.kill()
            proc.wait()
            raise RuntimeError(f"命令超时被终止: {cmd[0]}")
        except KeyboardInterrupt:
            # 不杀会留孤儿（mafft/fasttree 挂后台继续吃 CPU）
            proc.kill()
            proc.wait()
            raise
        finally:
            stop_wd.set()
        if timed_out.is_set():
            raise RuntimeError(
                f"命令超时（运行 {timeout}s 后被强制终止）: {cmd[0]}"
                f"（可设环境变量 VP_CMD_TIMEOUT_SEC 调大，0=不限制）")
    task_check_cancel()
    rc = proc.returncode
    if logger:
        logger.log(f"[退出码 {rc}] 耗时 {time.time() - t0:.1f}s", "OUT")
    if rc != 0:
        raise RuntimeError(
            f"命令失败 (退出码 {rc}): {subprocess.list2cmdline(cmd)}\n"
            + "\n".join(err_lines[-10:]))
    return rc


# ------------------------------------------------------------------
# 外部工具版本探测（复现依据：同一样品不同机器/不同轮的版本留痕）
# ------------------------------------------------------------------
_VERSION_CMD_VARIANTS = (('--version',), ('-version',))
_TOOL_VERSION_CACHE = {}          # exe 路径 → 版本串（进程内缓存）


def _parse_tool_version(raw):
    """从 --version 输出里取一行有意义的版本文本；取不到返回 'unknown'。"""
    for line in decode_output(raw or b'').splitlines():
        s = line.strip()
        if s:
            return s[:200]
    return 'unknown'


def probe_tool_versions(tools):
    """探测外部工具版本。tools: {tool_key: exe 路径} → {tool_key: 版本串}。

    对每个工具依次试 `<exe> --version` / `<exe> -version`（10s 超时），
    把 stdout/stderr 首行作为版本串；两条命令都失败或超时记 'unknown'。
    进程内按 exe 路径缓存——同一次/同一天多次分析不重复 fork 子进程。
    只读探测，绝不因某个工具异常而抛出（版本记录是附加信息，不阻断流程）。
    """
    out = {}
    for key, exe in (tools or {}).items():
        exe = str(exe or '')
        if not exe:
            out[key] = 'unknown'
            continue
        if exe in _TOOL_VERSION_CACHE:
            out[key] = _TOOL_VERSION_CACHE[exe]
            continue
        ver = 'unknown'
        for args in _VERSION_CMD_VARIANTS:
            try:
                r = subprocess.run(
                    [exe] + list(args), capture_output=True, timeout=10,
                    creationflags=(getattr(subprocess, 'CREATE_NO_WINDOW', 0)
                                   if sys.platform == 'win32' else 0))
                raw = (r.stdout or b'') + b'\n' + (r.stderr or b'')
                cand = _parse_tool_version(raw)
                if cand != 'unknown':
                    ver = cand
                    break
            except Exception:
                # 超时/权限/非可执行（如 .bat 需 cmd 外壳）→ 换下一种写法；
                # 全失败记 unknown，调用方照常运行
                continue
        _TOOL_VERSION_CACHE[exe] = ver
        out[key] = ver
    return out


# ------------------------------------------------------------------
# 步骤断点标记（work_dir 均位于平台 results/ 内）
# ------------------------------------------------------------------
def step_done_file(work_dir, step):
    step = check_step(step)
    wd = check_path(work_dir, in_platform=True, must_exist=True)
    return os.path.join(wd, f".{step}.done")


def is_step_done(work_dir, step):
    try:
        return os.path.isfile(step_done_file(work_dir, step))
    except FileNotFoundError:
        # 目录尚不存在（如全新样品首次跑到该步）= 未完成，而非异常——
        # 调用方（如管道⑧的 conserved 守卫）据此走正常首次运行分支
        return False


def mark_step_done(work_dir, step):
    with safe_open(step_done_file(work_dir, step), 'wt') as f:
        f.write(time.strftime('%Y-%m-%d %H:%M:%S'))


def clear_step_done(work_dir, step=None):
    wd = check_path(work_dir, in_platform=True, must_exist=True)
    targets = []
    if step is None:
        for fn in os.listdir(wd):
            if fn.startswith('.') and fn.endswith('.done'):
                targets.append(check_path(os.path.join(wd, fn), in_platform=True))
    else:
        step = check_step(step)
        p = os.path.join(wd, f".{step}.done")
        if os.path.isfile(p):
            targets.append(check_path(p, in_platform=True))
    for p in targets:
        try:
            os.remove(p)
        except OSError:
            pass


# ------------------------------------------------------------------
# ------------------------------------------------------------------
# SVG 产物合法性自检
# ------------------------------------------------------------------
def verify_svg(path, logger=None):
    """SVG 产物必须是良构 XML（否则矢量交付物损坏、嵌入报告时浏览器只能
    容错渲染残图）。2026-10-04 审查实测：examples/results/synteny/
    compare.svg 曾有 18 处 </path 闭合缺 '>'，全部残在同一张图里。

    校验失败时自动尝试幂等修复（补齐缺 '>' 的闭合标签），仍失败则告警
    返回 False；成功/修复成功返回 True。只在产物落盘后调用。
    """
    import re as _re
    import xml.etree.ElementTree as _ET
    try:
        _ET.parse(path)
        return True
    except _ET.ParseError:
        pass
    try:
        with open(path, encoding='utf-8') as f:
            svg = f.read()
        fixed, n_fix = _re.subn(r'</([a-zA-Z][\w-]*)(?![\w>])', r'</\1>', svg)
        _ET.fromstring(fixed)
        with open(path, 'w', encoding='utf-8', newline='\n') as f:
            f.write(fixed)
        if logger:
            logger.log(f"SVG 非良构，已自动修复 {n_fix} 处闭合标签: {path}",
                       "WARN")
        return True
    except Exception as e:
        if logger:
            logger.log(f"SVG 非良构且无法自动修复（矢量图可能渲染残缺）: "
                       f"{path} · {type(e).__name__}: {e}", "WARN")
        return False


# FASTA / FASTA 流式处理
# ------------------------------------------------------------------
def iter_fasta(path):
    """迭代 FASTA：yield (header_no_arrow, seq)"""
    header, chunks = None, []
    with safe_open(path) as f:
        for line in f:
            line = line.rstrip('\n\r')
            if line.startswith('>'):
                if header is not None:
                    yield header, ''.join(chunks)
                header, chunks = line[1:], []
            elif line:
                chunks.append(line)
    if header is not None:
        yield header, ''.join(chunks)


def write_fasta_record(fh, name, seq, width=60):
    fh.write(f">{name}\n")
    if width and len(seq) > width:
        for i in range(0, len(seq), width):
            fh.write(seq[i:i + width] + "\n")
    else:
        fh.write(seq + "\n")


def iter_fastq(path):
    """迭代 FASTQ：yield (name_no_at, seq, plus, qual)。"""
    with safe_open(path) as f:
        while True:
            l1 = f.readline()
            if not l1:
                break
            l2 = f.readline()
            l3 = f.readline()
            l4 = f.readline()
            if not l4:
                raise ValueError(f"FASTQ 记录数不是 4 的倍数（文件截断?）: {path}")
            yield l1.rstrip('\n\r'), l2.rstrip('\n\r'), l3.rstrip('\n\r'), l4.rstrip('\n\r')


def iter_fastq_records(path):
    """迭代 FASTQ：yield 4 行原始文本（原样写出用）。

    path 可为文件路径或已打开的文本文件对象（调用方负责关闭）。
    """
    if hasattr(path, 'readline'):
        f = path
        outer = None
    else:
        f = outer = safe_open(path)
    try:
        while True:
            lines = [f.readline() for _ in range(4)]
            if not lines[0]:
                break
            if not lines[3]:
                raise ValueError(f"FASTQ 记录数不是 4 的倍数（文件截断?）: {path}")
            yield lines
    finally:
        if outer is not None:
            outer.close()


def count_fasta_seqs(path):
    n = 0
    with safe_open(path) as f:
        for line in f:
            if line.startswith('>'):
                n += 1
    return n


def load_alignment(path):
    """读比对 FASTA → dict {name: 大写序列}（保序，允许 - 与 ?）。

    name 取 FASTA 头**首个空白分隔 token** —— 与 phylo/viz 全平台的叶名口径
    一致（Newick 叶名也是按首个 token 净化的），保证树上叶名能回映射到序列。
    原为 rdp_detect.load_alignment；自研三序列法移除后上移到通用层，
    供系统地理 / 重组结果摘要共用。
    """
    out = {}
    for header, seq in iter_fasta(path):
        tokens = header.split()
        out[tokens[0] if tokens else header.strip()] = seq.upper()
    return out


def variable_cols(seqs):
    """比对的可变位点下标：该列存在 ≥2 种非 gap 字符。

    seqs 为等长序列列表。gap 含 '-' 与 '?'，任一序列为 gap 的列不参与判定。
    原为 rdp_detect.variable_cols（三序列法的特征位池）。
    """
    if not seqs:
        return []
    out = []
    for i in range(len(seqs[0])):
        seen = set()
        for s in seqs:
            c = s[i]
            if c in '-?':
                continue
            seen.add(c)
            if len(seen) >= 2:
                out.append(i)
                break
    return out


# ------------------------------------------------------------------
# FASTA 头 taxid 注入
# ------------------------------------------------------------------
def _clean_seq_id(header, fallback_no):
    """从 FASTA 头提取规范 seq_id：取首空白前 token、幂等去除历史
    kraken:taxid 标签；空 ID 自动编号 seq_1/seq_2/...（kunpeng 要求非空）。
    返回 (seq_id, 其余描述部分)。"""
    tokens = header.split(None, 1)
    seq_id = tokens[0] if tokens else ''
    desc = tokens[1] if len(tokens) > 1 else ''
    if '|kraken:taxid|' in seq_id:
        seq_id = seq_id.split('|kraken:taxid|')[0]
    if not seq_id:
        seq_id = f'seq_{fallback_no}'
    return seq_id, (' ' + desc if desc else '')


def inject_taxid_to_fasta(in_fasta, out_fasta, taxid, logger=None):
    """为 FASTA 每条序列头注入 taxid 标签（kunpeng add-library 要求）。

    kunpeng 解析要求标签位于第一个空白分隔 token 内: >ID|kraken:taxid|N
    """
    tag = f"|kraken:taxid|{taxid}"
    n = 0
    with safe_open(in_fasta) as fin, safe_open(out_fasta, 'wt') as fout:
        for line in fin:
            if line.startswith('>'):
                header = line[1:].rstrip('\n\r')
                n += 1
                if 'kraken:taxid|' in header:
                    fout.write('>' + header + '\n')   # 已有标签，保留
                    continue
                seq_id, rest = _clean_seq_id(header, n)
                fout.write(f'>{seq_id}{tag}{rest}\n'.rstrip() + '\n')
            else:
                fout.write(line if line.endswith('\n') else line + '\n')
    if logger:
        logger.log(f"taxid 注入完成: {n} 条序列 -> {os.path.basename(out_fasta)} (taxid={taxid})")
    return n


def inject_taxid_chunked(in_fasta, out_fasta, taxid, chunk_bp=1000000,
                         overlap_bp=34, logger=None):
    """注入 taxid 并把每条序列切成 ≤chunk_bp bp 的片段（kunpeng add-library 用）。

    kunpeng convert 阶段按 60 条/批整批载入序列且无字节上限：多条大染色体
    （如玉米 chr01 176MB）同批时内存随批内序列总长暴涨，宿主库 1.77GB 基因组
    实测触发 >20GB 分配失败。切成小片段后批内存有界，且 kmer 分类只依赖
    k-mer→taxid 映射，片段边界仅丢失每切口 k-1 个 k-mer，对去宿主无影响。

    overlap_bp: 相邻片段的重叠碱基数，默认 34 = k-1（kunpeng 默认 k=35），
    使跨切口 k-mer 全部保留、kmer 集合与不切割时完全一致（零损失）；0 关闭。

    片段标题: >{orig_id}.p{序号}|kraken:taxid|N {原描述}
    返回 (原始序列数, 输出片段数)。
    """
    tag = f"|kraken:taxid|{taxid}"
    n_seq = n_frag = 0
    seq_id = desc = ''
    part = 0
    acc = []          # 当前未写出序列的纯碱基块（无换行）
    acc_len = 0
    seen_ids = set()  # 记录已用 ID，重复 ID 自动加序号去重

    def flush():
        nonlocal n_frag, acc, acc_len
        if not acc_len:
            return
        # 尾部残余 ≤ overlap_bp 时内容已完全含于上一片段，无需重复写出
        if part > 0 and acc_len <= overlap_bp:
            acc, acc_len = [], 0
            return
        _write_piece(fout, seq_id, desc, tag, part, ''.join(acc))
        n_frag += 1
        acc = []
        acc_len = 0

    with safe_open(in_fasta) as fin, safe_open(out_fasta, 'wt') as fout:
        for line in fin:
            if line.startswith('>'):
                flush()
                header = line[1:].rstrip('\n\r')
                n_seq += 1
                seq_id, desc = _clean_seq_id(header, n_seq)
                if seq_id in seen_ids:          # 重复 ID 去重（如 draft 同名 scaffold）
                    seq_id = f'{seq_id}_r{n_seq}'
                seen_ids.add(seq_id)
                part = 0
                continue
            s = line.strip()
            if not s:
                continue
            acc.append(s)
            acc_len += len(s)
            if acc_len >= chunk_bp:
                # 依次切出整片段；下一片段从切点前 overlap_bp 处开始，
                # 跨切口的 k-mer 因此在相邻两片中各存在一份
                seq = ''.join(acc)
                while len(seq) >= chunk_bp:
                    _write_piece(fout, seq_id, desc, tag, part,
                                 seq[:chunk_bp])
                    n_frag += 1
                    part += 1
                    keep = chunk_bp - overlap_bp if overlap_bp else chunk_bp
                    seq = seq[keep:]
                acc = [seq]
                acc_len = len(seq)
        flush()
    if logger:
        logger.log(f"taxid 注入+分块完成: {n_seq} 条序列 -> {n_frag} 个片段"
                   f"（≤{chunk_bp // 1000}kb/片段, 重叠{overlap_bp}bp, "
                   f"taxid={taxid}）")
    return n_seq, n_frag


def _write_piece(fout, seq_id, desc, tag, part, seq):
    head = f"{seq_id}.p{part:05d}{tag}"
    fout.write(f">{head}{desc}\n".rstrip() + "\n")
    for i in range(0, len(seq), 60):
        fout.write(seq[i:i + 60] + "\n")


def inject_taxid_map(in_fasta, out_fasta, acc2taxid, logger=None, missing_fatal=True):
    """按 accession→taxid 映射逐条注入 kraken:taxid|N。

    acc2taxid: {accession: taxid(int)}；FASTA 头第一个空格前的 ID 用于查询。
    返回 (注入数, 未匹配列表)。
    """
    n, missing = 0, []
    with safe_open(in_fasta) as fin, safe_open(out_fasta, 'wt') as fout:
        for line in fin:
            if line.startswith('>'):
                header = line[1:].rstrip('\n\r')
                seq_id = header.split()[0] if header.split() else ''
                rest = header[len(seq_id):]
                taxid = acc2taxid.get(seq_id)
                if taxid is None:
                    missing.append(seq_id)
                    if missing_fatal:
                        raise KeyError(
                            f"序列 {seq_id} 在 info 表中未找到 Taxid（全部缺失项已记录）")
                    fout.write('>' + header + '\n')
                else:
                    fout.write(f'>{seq_id}|kraken:taxid|{taxid}{rest}\n'.rstrip() + '\n')
                n += 1
            else:
                fout.write(line if line.endswith('\n') else line + '\n')
    if logger:
        logger.log(f"taxid 映射注入: {n} 条, 未匹配 {len(missing)} 条")
        if missing:
            logger.log("未匹配示例: " + ", ".join(missing[:10]), "WARN")
    if missing and missing_fatal:
        raise KeyError(f"共 {len(missing)} 条序列缺少 Taxid 映射")
    return n, missing


# ------------------------------------------------------------------
# 杂项
# ------------------------------------------------------------------
def fmt_size(nbytes):
    for unit in ('B', 'KB', 'MB', 'GB', 'TB'):
        if nbytes < 1024:
            return f"{nbytes:.1f}{unit}"
        nbytes /= 1024
    return f"{nbytes:.1f}PB"


def sample_name_from_fastq(path):
    """NX-5_S2_L001_R1_001.fastq.gz -> NX-5"""
    base = os.path.basename(str(path))
    for suffix in ('.fastq.gz', '.fq.gz', '.fastq', '.fq', '.fasta.gz', '.fa.gz', '.fasta', '.fa', '.fna'):
        if base.lower().endswith(suffix):
            base = base[:-len(suffix)]
            break
    # 去 Illumina 命名尾部
    for pat in ('_R1_001', '_R2_001', '_R1', '_R2', '_1', '_2'):
        if base.endswith(pat):
            base = base[:-len(pat)]
            break
    base = base.split('_L00')[0]
    return base or 'sample'


# ------------------------------------------------------------------
# 样品名归一化 —— 全平台唯一真源
#
# 历史问题：样品名到目录名的转换曾经有 4 份各自为政的实现
#   pipeline._safe_sample_name   折叠非法字符 + 去首尾 + 截断 50（建目录用）
#   web/common._safe_sample      逐字符替换，不去首尾、不截断（读路径用）
#   logan_trace.safe_name        同 web/common
#   utils.sample_name_from_fastq 仅由文件名推导（用于生成初始候选名）
# 前两者的差异在中文/连续特殊字符/超长名下会显形：样品名「样品A」建成目录 `A`，
# 读路径却去找 `__A` —— 用户读不到自己刚建的样品，甚至撞名丢数据。
# 现在收敛成下面两个函数：safe_sample_name（纯转换）+ resolve_sample_name（查真实目录）。
# ------------------------------------------------------------------
SAFE_SAMPLE_MAX = 50


def safe_sample_name(name):
    """样品名 → 目录名（纯函数，唯一权威实现）。

    规则（与历史 pipeline._safe_sample_name 逐字节一致，保证已有目录名不漂移）：
      1. 非 [A-Za-z0-9_\\-.] 的**连续**字符折叠成一个 '_'
      2. 去掉首尾的 '.' / '_' / '-'
      3. 截断到 50 字符（Windows 全路径长度安全）；空则回退 'sample'

    为什么不做「保留中文」：样品名直接当目录名进 SPAdes/BLAST(LMDB)/minibwa/
    snpEff 等外部工具的路径（见 DEVELOPMENT_NOTES 环境与路径 1），
    非 ASCII 路径这些工具普遍不支持。折叠成 ASCII 是刻意的安全设计，
    代价只是"用户填的名字与实际目录名可能不同"——这一点必须显式告知用户
    （见 web/samples.api_pipeline_create 的 note 返回），不能静默。
    """
    s = re.sub(r'[^A-Za-z0-9_\-.]+', '_', str(name)).strip('._-')
    return s[:SAFE_SAMPLE_MAX] or 'sample'


def resolve_sample_name(name, base=None):
    """把请求里的样品名解析成 base 目录下**真实存在**的子目录名。

    为什么需要：光靠字符替换无法命中"已经躺在磁盘上的目录"。
    手工用资源管理器建的目录、历史遗留命名、旧版本写下的非规范名
    （含中文、空格、超长名）都会被字符替换改成另一个名字 → 读不到。
    这里改为以磁盘为准：
      1. 含路径分隔符 / 上跳 / 空的名字一律不原样接受，直接用规范名
      2. 原名命中 base 的子目录 → 用原名（人工/历史目录也能正常读写）
      3. 规范名命中 → 用规范名
      4. 都不命中 → 返回规范名，由调用方决定是新建还是报 404
    第 2 步返回的名字直接取自 os.listdir，天然不含分隔符，比字符替换更安全；
    路径越界仍由 check_path(in_platform=True) 兜底。
    """
    raw = str(name or '')
    canon = safe_sample_name(raw)
    if (not raw or raw in ('.', '..') or '/' in raw or '\\' in raw
            or (os.altsep and os.altsep in raw)):
        return canon
    if base is None:
        from .config import DIRS
        base = DIRS['results']
    try:
        names = set(os.listdir(base))
    except OSError:
        return canon
    return raw if raw in names else canon


def open_in_explorer(path):
    """在 Windows 资源管理器中打开目录/文件，并尽力把新窗口置顶。

    为什么不直接 os.startfile：目录已有窗口或 Shell 复用窗口时，窗口常
    开在当前前台窗口**后面**，用户点了「📂 打开」以为没反应（2026-09-20
    用户实测反馈「点击文件夹，展示在页面最前端」）。这里 startfile 之后
    用 EnumWindows 找新出现的 CabinetWClass 窗口 SetForegroundWindow；
    没等到新窗口（复用旧窗口）就把现有资源管理器窗口顶上来。
    非兜底不可时的 SetForegroundWindow 前台锁用「先敲一下 ALT」解锁。
    """
    path = str(path)
    if sys.platform != 'win32':
        try:
            os.startfile(path)          # type: ignore[attr-defined]
        except AttributeError:
            subprocess.Popen(['xdg-open', path])
        return

    import ctypes
    from ctypes import wintypes

    user32 = ctypes.windll.user32
    EnumProc = ctypes.WINFUNCTYPE(ctypes.c_bool, wintypes.HWND,
                                  wintypes.LPARAM)

    def _cabinet_windows():
        out = []
        buf = ctypes.create_unicode_buffer(64)

        def cb(hwnd, _lparam):
            user32.GetClassNameW(hwnd, buf, 64)
            if buf.value == 'CabinetWClass' and user32.IsWindowVisible(hwnd):
                out.append(hwnd)
            return True
        user32.EnumWindows(EnumProc(cb), 0)
        return out

    def _raise(hwnd):
        user32.keybd_event(0x12, 0, 0, 0)       # ALT down：解前台锁
        user32.SetForegroundWindow(hwnd)
        user32.keybd_event(0x12, 0, 2, 0)       # ALT up

    before = set(_cabinet_windows())
    try:
        os.startfile(path)
    except OSError:
        subprocess.Popen(['explorer', path])

    # 置顶等待放后台线程：HTTP 请求立即返回，不被最多 2s 的窗口探测拖住
    def _focus():
        deadline = time.time() + 2.0
        while time.time() < deadline:
            new = set(_cabinet_windows()) - before
            if new:
                _raise(next(iter(new)))
                return
            time.sleep(0.15)
        hwnds = _cabinet_windows()
        if hwnds:
            _raise(hwnds[0])
    threading.Thread(target=_focus, daemon=True,
                     name='explorer-focus').start()
