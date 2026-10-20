# -*- coding: utf-8 -*-
"""插件系统：外部工具 / 数据库包的安装、卸载、扫描与解析。

参考 PhyloSuite 的插件机制（zip 包 + 设置页 Install/Uninstall），但元数据
跟着包走：每个插件 zip 顶层带 plugin.json 清单，主程序无需重新编译即可
识别新插件——这是与 PhyloSuite（清单编译进 exe）的本质区别。

目录约定（与 3rd/ 平级，不进仓库）：
    <平台根>/plugins/<插件名>/plugin.json
    <平台根>/plugins/<插件名>/...            # 工具本体（相对路径由清单给出）

清单格式：
    {
      "name": "mafft",                 # 必填，安装后的目录名
      "version": "7.526",              # 必填
      "kind": "tool",                  # 必填：tool | database
      "description": "多序列比对",
      "provides": {"mafft": "mafft-win64/mafft.bat"}        # kind=tool
      "databases": {"virus": "db/plant"}                    # kind=database
    }

解析优先级（见 config.Config）：platform.json 显式覆盖 > 插件 > 自动探测。
插件只「填缺」：不动用户手工指定与 3rd/ 内已探测到的工具，卸载插件后
原探测结果自然回接，配置零残留。

zip 安全：解压前逐成员校验路径（拒绝绝对路径 / .. 穿越 / 盘符），
解压到临时目录、清单校验通过后再原子换入——坏包不会污染 plugins/。
"""
import hashlib
import json
import os
import re
import shutil
import stat
import tempfile
import time
import zipfile

if __name__ == '__main__':  # 直接脚本运行时保证可 import 平台包
    import sys
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from Virus_Platform_Core.config import PLATFORM_ROOT

PLUGIN_DIR = os.path.join(PLATFORM_ROOT, 'plugins')
MANIFEST_NAME = 'plugin.json'
_NAME_RE = re.compile(r'^[\w][\w.-]{0,63}$')     # 插件名即目录名，收紧字符集
_TMP_PREFIX = '.tmp-'
_OLD_PREFIX = '.old-'
KINDS = ('tool', 'database')

# 插件源索引（dev_tools/build_plugin_repo.py 生成，随仓库提交，raw URL 稳定；
# 可被 platform.json 的 plugin_index_url 覆盖）。
# 默认按顺序尝试：国内镜像站（自建，快且稳）→ GitHub raw（海外兜底）。
# 索引本身很小（~12KB），两处内容一致；镜像站不可达时自动切 GitHub。
DEFAULT_INDEX_URLS = [
    'http://39.106.101.94/plugins/index.json',
    'https://raw.githubusercontent.com/zhangwenda0518/'
    'plant-virus-analysis-platform/main/plugins-repo/index.json',
]
DEFAULT_INDEX_URL = DEFAULT_INDEX_URLS[0]


class PluginError(Exception):
    """安装/卸载失败（清单非法、zip 损坏、路径穿越等），消息直接展示给用户。"""


# ------------------------------------------------------------------
# 清单
# ------------------------------------------------------------------
def read_manifest(plugin_dir):
    """读插件目录的清单；无清单/清单非法返回 (None, 原因)。"""
    mp = os.path.join(plugin_dir, MANIFEST_NAME)
    if not os.path.isfile(mp):
        return None, f'缺少 {MANIFEST_NAME}'
    try:
        with open(mp, 'r', encoding='utf-8') as f:
            mf = json.load(f)
    except (OSError, ValueError) as e:
        return None, f'{MANIFEST_NAME} 解析失败: {e}'
    err = _validate_manifest(mf, plugin_dir)
    if err:
        return None, err
    return mf, None


def _validate_manifest(mf, base_dir):
    """清单最小校验；返回错误消息或 None。base_dir 用于核对 provides 路径存在。"""
    if not isinstance(mf, dict):
        return f'{MANIFEST_NAME} 必须是 JSON 对象'
    name = mf.get('name')
    if not isinstance(name, str) or not _NAME_RE.match(name):
        return f'name 缺失或非法（需匹配 {_NAME_RE.pattern}）: {name!r}'
    if not str(mf.get('version') or '').strip():
        return 'version 缺失'
    kind = mf.get('kind') or 'tool'
    if kind not in KINDS:
        return f'kind 必须是 {" / ".join(KINDS)}，当前: {kind!r}'
    mf['kind'] = kind
    key = 'provides' if kind == 'tool' else 'databases'
    mapping = mf.get(key)
    if not isinstance(mapping, dict) or not mapping:
        return f'kind={kind} 需要非空的 {key} 映射'
    for k, rel in mapping.items():
        if not isinstance(k, str) or not k.strip():
            return f'{key} 键非法: {k!r}'
        if not isinstance(rel, str) or not rel.strip():
            return f'{key}[{k}] 缺少相对路径'
        p = os.path.join(base_dir, rel)
        if not os.path.exists(p):
            return f'{key}[{k}] 指向的文件/目录不存在: {rel}'
    return None


# ------------------------------------------------------------------
# 扫描与解析
# ------------------------------------------------------------------
def installed_plugins(root=None):
    """扫描 plugins/ 下已安装的插件，返回清单列表（含 bad 标记的坏插件）。"""
    root = root or PLUGIN_DIR
    out = []
    if not os.path.isdir(root):
        return out
    for entry in sorted(os.listdir(root)):
        pdir = os.path.join(root, entry)
        if entry.startswith('.') or not os.path.isdir(pdir):
            continue
        mf, err = read_manifest(pdir)
        if mf is None:
            out.append({'name': entry, 'kind': 'unknown', 'version': '',
                        'description': '', 'dir': pdir,
                        'provides': {}, 'ok': False, 'error': err})
            continue
        out.append(_describe(mf, pdir))
    return out


def _describe(mf, pdir):
    """清单 → 对外描述（provides 键附解析后的绝对路径与在位状态）。"""
    kind = mf.get('kind') or 'tool'
    key = 'provides' if kind == 'tool' else 'databases'
    resolved = {}
    for k, rel in (mf.get(key) or {}).items():
        ap = os.path.normpath(os.path.join(pdir, rel))
        resolved[k] = ap
    ok = bool(resolved) and all(os.path.exists(v) for v in resolved.values())
    return {'name': mf['name'], 'kind': kind,
            'version': str(mf.get('version') or ''),
            'description': str(mf.get('description') or ''),
            'dir': pdir, 'map': resolved, 'ok': ok, 'error': None}


def resolve_plugins(cfg=None, root=None):
    """扫描插件 → (tool_map, db_map)：键 = Cfg.tools / Cfg.databases 的键，
    值 = 绝对路径。只收文件/目录真实在位的条目（坏插件静默失效，不拦启动）。"""
    tools, dbs = {}, {}
    for it in installed_plugins(root):
        if not it['ok']:
            continue
        if it['kind'] == 'tool':
            for k, ap in it['map'].items():
                if os.path.isfile(ap):
                    tools[k] = ap
        else:
            for k, ap in it['map'].items():
                if os.path.exists(ap):
                    dbs[k] = ap
    return tools, dbs


# ------------------------------------------------------------------
# 安装 / 卸载
# ------------------------------------------------------------------
def install_zip(zip_path, root=None, expected_sha256=None):
    """从 zip 安装插件；返回清单 dict。已有同名插件 = 升级覆盖。

    流程：解压到临时目录 → 定位并校验清单 → 原子换入 plugins/<name>/。
    清单位置：zip 顶层，或唯一顶层目录内（Windows 右键「压缩为 zip」
    常见布局都兼容）。expected_sha256 给出时先验哈希（插件源索引带
    sha256，防止下载损坏/被换包）。"""
    root = root or PLUGIN_DIR
    if not os.path.isfile(zip_path):
        raise PluginError(f'zip 不存在: {zip_path}')
    if expected_sha256:
        got = _sha256_file(zip_path)
        if got.lower() != str(expected_sha256).lower():
            raise PluginError(
                f'zip 校验失败：sha256 不匹配（期望 {expected_sha256[:12]}…，'
                f'实际 {got[:12]}…）——下载可能损坏或被替换，请重试')
    os.makedirs(root, exist_ok=True)
    tmp = tempfile.mkdtemp(prefix=_TMP_PREFIX, dir=root)
    try:
        _safe_extract(zip_path, tmp)
        mf, sub, err = _locate_manifest(tmp)
        if mf is None:
            raise PluginError(err or f'zip 内未找到 {MANIFEST_NAME}（顶层或唯一子目录内）')
        err = _validate_manifest(mf, sub)
        if err:
            raise PluginError(err)
        dest = os.path.join(root, mf['name'])
        _swap_into_place(sub, dest)
        return mf
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def _locate_manifest(extract_dir):
    """在解压结果里定位清单：顶层优先，其次唯一顶层目录。
    返回 (manifest, 插件根, 错误)——「有清单但非法」要报具体原因，
    不能混进「没清单」里（否则用户无从修包）。"""
    top = os.path.join(extract_dir, MANIFEST_NAME)
    if os.path.isfile(top):
        try:
            with open(top, 'r', encoding='utf-8') as f:
                return json.load(f), extract_dir, None
        except ValueError as e:
            return None, None, f'{MANIFEST_NAME} 解析失败: {e}'
    entries = [e for e in os.listdir(extract_dir) if not e.startswith('.')]
    if len(entries) == 1:
        sub = os.path.join(extract_dir, entries[0])
        if os.path.isdir(sub):
            mf, err = read_manifest(sub)
            if mf is not None:
                return mf, sub, None
            if err and '缺少' not in err:
                return None, None, err
    return None, None, None


def _safe_extract(zip_path, dest):
    """防御性解压：逐成员校验路径后写盘（zipfile.extract 的 Windows 补丁
    历史上放过绝对路径/盘符，这里自己控制，不依赖版本行为）。"""
    try:
        with zipfile.ZipFile(zip_path) as z:
            members = z.namelist()
            for m in members:
                if m.endswith('/'):
                    continue
                _check_member_path(m)
                z.extract(m, dest)
    except zipfile.BadZipFile as e:
        raise PluginError(f'zip 损坏或不是 zip: {e}') from e


def _check_member_path(member):
    """zip 成员路径校验：拒绝绝对路径、盘符、.. 穿越、NUL。"""
    if '\x00' in member:
        raise PluginError(f'zip 内含非法路径（NUL）: {member!r}')
    norm = member.replace('\\', '/')
    if re.match(r'^[A-Za-z]:', norm):
        raise PluginError(f'zip 内含绝对路径（盘符）: {member!r}')
    if norm.startswith('/'):
        raise PluginError(f'zip 内含绝对路径: {member!r}')
    parts = [p for p in norm.split('/') if p not in ('', '.')]
    if any(p == '..' for p in parts):
        raise PluginError(f'zip 内含 .. 路径穿越: {member!r}')


def _swap_into_place(src, dest):
    """src → dest 的原子换入：先挪走旧目录，成功后再删；中途失败可回滚。"""
    old = None
    if os.path.isdir(dest):
        old = os.path.join(os.path.dirname(dest),
                           _OLD_PREFIX + os.path.basename(dest)
                           + time.strftime('.%Y%m%d%H%M%S'))
        os.rename(dest, old)
    try:
        os.rename(src, dest)
    except OSError:
        if old and os.path.isdir(old):
            os.rename(old, dest)       # 回滚，旧插件原样保留
        raise
    if old:
        _rmtree_force(old)


def _rmtree_force(path):
    """shutil.rmtree 的只读位兜底版（Windows 上 zip 解出的文件常带只读位）。"""

    def _onerr(func, p, _exc):
        try:
            os.chmod(p, stat.S_IWRITE)
            func(p)
        except OSError:
            pass

    shutil.rmtree(path, onerror=_onerr)


def uninstall(name, root=None):
    """卸载插件：删除目录。运行中的任务若正引用该插件工具会失败——
    与 PhyloSuite 同口径（直接删），安装器不替代任务中心的并发管理。"""
    root = root or PLUGIN_DIR
    if not _NAME_RE.match(name or ''):
        raise PluginError(f'插件名非法: {name!r}')
    pdir = os.path.join(root, name)
    if not os.path.isdir(pdir):
        raise PluginError(f'插件未安装: {name}')
    _rmtree_force(pdir)


# ------------------------------------------------------------------
# 下载安装（对齐 PhyloSuite 的镜像/URL 直装）
# ------------------------------------------------------------------
def install_url(url, root=None, progress=None, expected_sha256=None):
    """下载 zip 到临时文件后走 install_zip。progress(done, total) 可选回调。"""
    return install_urls(url if isinstance(url, (list, tuple)) else [url],
                        root=root, progress=progress,
                        expected_sha256=expected_sha256)


def install_urls(urls, root=None, progress=None, expected_sha256=None):
    """按顺序尝试多个下载地址（国内镜像优先、GitHub 兜底），首个成功即返回。

    每个地址失败后重试一次（实测镜像站偶发瞬时 502，重试即可恢复），
    仍失败才切下一个地址。全部失败时抛 PluginError，消息里带上每个地址的
    失败原因——用户能看出是网络问题还是地址失效。"""
    if isinstance(urls, str):
        urls = [urls]
    urls = [u for u in (urls or []) if u]
    if not urls:
        raise PluginError('未提供下载地址')
    errors = []
    for u in urls:
        for attempt in (1, 2):
            try:
                return _install_one_url(u, root=root, progress=progress,
                                        expected_sha256=expected_sha256)
            except Exception as e:             # 先原地重试，再换下一个镜像
                if attempt == 2:
                    errors.append(f'{u} → {type(e).__name__}: {e}')
                else:
                    time.sleep(1.5)
    raise PluginError('所有下载地址均失败：\n' + '\n'.join(errors))


def _install_one_url(url, root=None, progress=None, expected_sha256=None):
    if not (isinstance(url, str) and url.lower().startswith(('http://', 'https://'))):
        raise PluginError(f'URL 非法（仅支持 http/https）: {url!r}')
    import urllib.request
    root = root or PLUGIN_DIR
    os.makedirs(root, exist_ok=True)
    fd, tmp_zip = tempfile.mkstemp(prefix='.dl-', suffix='.zip', dir=root)
    os.close(fd)
    try:
        req = urllib.request.Request(url, headers={'User-Agent': 'PVAP-Plugin/1.0'})
        with urllib.request.urlopen(req, timeout=60) as resp, \
                open(tmp_zip, 'wb') as f:
            total = int(resp.headers.get('Content-Length') or 0)
            done = 0
            while True:
                chunk = resp.read(1 << 20)
                if not chunk:
                    break
                f.write(chunk)
                done += len(chunk)
                if progress:
                    progress(done, total)
        return install_zip(tmp_zip, root, expected_sha256=expected_sha256)
    finally:
        try:
            os.remove(tmp_zip)
        except OSError:
            pass


# ------------------------------------------------------------------
# 插件源索引（自动检测缺失 → 一键下载解压）
# ------------------------------------------------------------------
_index_cache = {'ts': 0, 'data': None, 'url': ''}
_INDEX_TTL = 600            # 索引进程内缓存 10 分钟


def _sha256_file(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        while True:
            b = f.read(1 << 20)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def fetch_index(url=None, force=False, timeout=30):
    """拉取插件源索引（进程内缓存 10 分钟）。返回 (plugins 列表, 错误或 None)。

    url 为空时按 DEFAULT_INDEX_URLS 顺序回退（国内镜像 → GitHub raw）。
    url 显式给出时只用该地址（platform.json 覆盖/调试用）。"""
    urls = [url] if url else list(DEFAULT_INDEX_URLS)
    now = time.time()
    if (not force and _index_cache['data'] is not None
            and _index_cache['url'] == urls[0]
            and now - _index_cache['ts'] < _INDEX_TTL):
        return _index_cache['data'], None
    import urllib.request
    errors = []
    for u in urls:
        try:
            req = urllib.request.Request(u, headers={'User-Agent': 'PVAP-Plugin/1.0'})
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                data = json.loads(resp.read().decode('utf-8'))
            plugins = data.get('plugins')
            if not isinstance(plugins, list):
                raise ValueError("索引缺少 'plugins' 数组")
            clean = [p for p in plugins
                     if isinstance(p, dict) and p.get('name') and p.get('url')]
            _index_cache.update(ts=now, data=clean, url=u)
            return clean, None
        except Exception as e:
            errors.append(f'{u} → {type(e).__name__}: {e}')
    return None, '；'.join(errors)


def available_from_index(cfg, index_url=None, force=False):
    """索引 × 本机现状 → 可安装列表，标注自动检测结果。

    每条附：
      missing_keys   该插件提供、但本机当前缺失（探测不到也未由插件提供）的工具键
      installed      同名插件目录是否已存在
      needed         missing_keys 非空（缺失检测的核心标记）
    """
    plugins, err = fetch_index(index_url, force=force)
    if err:
        return None, err
    have = {k for k, v in cfg.tool_status().items() if v}
    installed = {it['name'] for it in installed_plugins()}
    out = []
    for p in plugins:
        keys = [k for k in (p.get('provides') or {})]
        missing = [k for k in keys if k not in have]
        out.append({**p,
                    'missing_keys': missing,
                    'installed': p.get('name') in installed,
                    'needed': bool(missing)})
    out.sort(key=lambda x: (not x['needed'], x['name']))
    return out, None


# ------------------------------------------------------------------
# Cfg 集成
# ------------------------------------------------------------------
def refresh_cfg_plugins(cfg):
    """重扫插件并填入 cfg.plugin_tools / cfg.plugin_databases（见 config.py）。"""
    cfg.plugin_tools, cfg.plugin_databases = resolve_plugins()
    return cfg.plugin_tools, cfg.plugin_databases


def source_of(cfg, key):
    """工具键的当前来源标签：platform / plugin / detect（诊断与 UI 用）。"""
    v = cfg.tools.get(key)
    if v and os.path.isfile(v):
        return 'platform' if key in getattr(cfg, '_platform_tool_keys', ()) \
            else 'detect'
    if cfg.plugin_tools.get(key):
        return 'plugin'
    return 'missing'


if __name__ == '__main__':
    for it in installed_plugins():
        tag = '✔' if it['ok'] else f"✘ {it.get('error')}"
        print(f"{tag} {it['name']} {it['version']} ({it['kind']})")
