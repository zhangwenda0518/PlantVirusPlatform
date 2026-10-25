# -*- coding: utf-8 -*-
"""设置中心 + 存储/磁盘水位（自 app.py 拆出）。

路由：
    GET  /settings              设置页
    GET  /api/settings          读取（语言/线程/邮箱/默认参数/目录/库状态）
    POST /api/settings          保存（逐项校验后持久化到 platform.json）
    GET  /api/storage           磁盘容量 + 目录占用（占用走后台扫描 + 缓存）
    POST /api/storage/rescan    强制重新统计目录占用
    GET  /api/selfcheck         环境自检（模块/工具/数据库/磁盘；?force=1 强制重跑）
    GET  /api/plugin/list       插件清单（plugins/ 目录扫描）
    GET  /api/plugin/available  插件源索引 × 本机现状（自动检测缺失工具）
    GET  /api/plugin/modules    模块视图（每模块依赖哪些工具、还缺哪些）
    POST /api/plugin/install_modules 按模块安装（工具去重、跳过已有）
    POST /api/plugin/install    安装插件 zip（multipart file= 或 JSON path=）
    POST /api/plugin/install_url 从 URL 下载安装（可带 sha256 校验）
    POST /api/plugin/uninstall  卸载插件（JSON name=）

设计要点：磁盘剩余用 shutil.disk_usage 毫秒级实时算；目录占用要遍历
6 万+ 文件（databases 12GB+），实测数十秒，因此走后台线程 + 10 分钟缓存，
绝不放在请求路径上同步扫描。
"""
import os
import tempfile
import threading
import time

from flask import Blueprint, abort, jsonify, render_template, request

from Virus_Platform_Core.config import DIRS, PLATFORM_ROOT
from Virus_Platform_Core.utils import fmt_size
from Virus_Platform_Core.web.state import cfg, tool_runs_root

bp = Blueprint('settings', __name__)

# ------------------------------------------------------------------
# 存储与磁盘水位
# ------------------------------------------------------------------
STORAGE_WARN_GB = 20          # 剩余低于此值触发水位横幅
_STORAGE_TTL = 600            # 目录占用缓存有效期（秒）
_storage = {'ts': 0, 'scanning': False, 'items': [], 'error': None,
            'scanned_at': None}
_storage_lock = threading.Lock()


def _dir_size(path):
    """累加目录内文件大小与文件数（os.walk，容忍无权限/已删除的文件）。"""
    total = 0
    count = 0
    for cur, _dirs, files in os.walk(path):
        for fn in files:
            try:
                total += os.path.getsize(os.path.join(cur, fn))
            except OSError:
                continue
            count += 1
    return total, count


def _scan_targets():
    """待统计的目录清单：(显示名, 路径, 类别)。"""
    items = [
        ('databases', DIRS['databases'], 'database'),
        ('host-db', DIRS['host_src'], 'database'),
        ('databases/virusref_db', DIRS['virus_src'], 'database'),
        ('results', DIRS['results'], 'output'),
        ('downloads', DIRS['downloads'], 'output'),
        ('tool_runs', tool_runs_root(), 'output'),
        ('meta_search', DIRS['meta_search'], 'output'),
        ('logan', DIRS['logan'], 'output'),
        ('logs', DIRS['logs'], 'output'),
        ('tasks', DIRS['tasks'], 'output'),
        ('fastq', DIRS.get('fastq', ''), 'input'),
        ('dist', os.path.join(PLATFORM_ROOT, 'dist'), 'build'),
    ]
    return [(n, p, k) for n, p, k in items if p and os.path.isdir(p)]


def _do_storage_scan():
    """后台扫描：平台各目录 + databases 下逐库细分。"""
    out = []
    try:
        for name, path, kind in _scan_targets():
            size, count = _dir_size(path)
            entry = {'name': name, 'path': path, 'kind': kind,
                     'size': size, 'size_h': fmt_size(size),
                     'files': count, 'children': []}
            # 数据库根太大，按子库拆开才有参考价值
            if name == 'databases':
                try:
                    subs = sorted(d for d in os.listdir(path)
                                  if os.path.isdir(os.path.join(path, d)))
                except OSError:
                    subs = []
                for sub in subs:
                    sp = os.path.join(path, sub)
                    ssz, scnt = _dir_size(sp)
                    entry['children'].append({
                        'name': sub, 'path': sp, 'kind': 'database',
                        'size': ssz, 'size_h': fmt_size(ssz), 'files': scnt})
                entry['children'].sort(key=lambda x: -x['size'])
            out.append(entry)
        out.sort(key=lambda x: -x['size'])
        with _storage_lock:
            _storage['items'] = out
            _storage['ts'] = time.time()
            _storage['scanned_at'] = time.strftime('%H:%M:%S')
            _storage['error'] = None
    except Exception as e:                      # 扫描失败不应影响平台
        with _storage_lock:
            _storage['error'] = str(e)
    finally:
        with _storage_lock:
            _storage['scanning'] = False


def _disk_snapshot():
    """盘容量快照（实时，含数据库根所在盘——可能与平台根不同盘）。"""
    import shutil
    from Virus_Platform_Core.config import get_database_root
    roots = [('平台', PLATFORM_ROOT)]
    db_root = get_database_root()
    if db_root:
        roots.append(('数据库', db_root))
    out = []
    for label, path in roots:
        try:
            u = shutil.disk_usage(path)
            out.append({'label': label, 'path': path,
                        'total': u.total, 'free': u.free, 'used': u.used,
                        'total_h': fmt_size(u.total),
                        'free_h': fmt_size(u.free),
                        'used_h': fmt_size(u.used),
                        'used_pct': round(u.used / u.total * 100, 1) if u.total else 0,
                        'warn': u.free < STORAGE_WARN_GB * 1e9})
        except OSError as e:
            out.append({'label': label, 'path': path, 'error': str(e),
                        'warn': False})
    return out


# ------------------------------------------------------------------
# 路由
# ------------------------------------------------------------------
@bp.route('/settings')
def page_settings():
    return render_template('settings.html')


@bp.route('/api/settings', methods=['GET'])
def api_settings_get():
    """设置中心读取：语言/线程/邮箱/默认参数/目录/数据库状态。"""
    from Virus_Platform_Core.config import Config
    fields = []
    for f in Config.DEFAULT_FIELDS:
        key, typ, dflt, lz, le = f[0], f[1], f[2], f[3], f[4]
        fields.append({'key': key, 'type': typ, 'def': dflt,
                       'label_zh': lz, 'label_en': le,
                       'opts': (f[5] if len(f) > 5 else [])})
    return jsonify({
        'language': cfg.language,
        'email': cfg.email,
        'threads': cfg.threads,
        'cpu_count': os.cpu_count() or 4,
        'defaults': cfg.defaults,
        'default_fields': fields,
        'databases': cfg.databases,
        'output_root': cfg.output_root,
        'input_root': cfg.input_root,
        'database_root': cfg.database_root,
        'examples_root': getattr(cfg, 'examples_root', ''),
        'examples_effective': DIRS.get('examples', ''),
        'plugin_index_url': getattr(cfg, 'plugin_index_url', ''),
        'plugin_source_pref': getattr(cfg, 'plugin_source_pref', 'auto'),
        'dirs': {'results': DIRS['results'],
                 'tool_runs': tool_runs_root(),
                 'databases': DIRS['databases'],
                 'fastq': DIRS.get('fastq')
                          or os.path.join(PLATFORM_ROOT, 'run', 'fastq'),
                 'uploads': DIRS.get('uploads')
                            or os.path.join(PLATFORM_ROOT, 'run', 'uploads'),
                 'host_src': DIRS.get('host_src')
                             or os.path.join(PLATFORM_ROOT, 'host-db'),
                 'virus_src': DIRS.get('virus_src')
                              or os.path.join(DIRS['databases'],
                                              'virusref_db'),
                 'logs': DIRS['logs']},
    })


@bp.route('/api/selfcheck', methods=['GET'])
def api_selfcheck():
    """环境自检（结构化）：模块导入 / 外部工具 / 数据库 / 磁盘空间。

    ?force=1 强制重跑（设置页「运行自检」按钮）；默认走进程内缓存
    （5 分钟）——首次收集要递归 import 全部子模块，实测数秒，总览页
    状态条与设置页共用同一份缓存，不重复跑。
    """
    from Virus_Platform_Core.selfcheck import collect
    force = request.args.get('force') in ('1', 'true')
    try:
        return jsonify(collect(force=force))
    except Exception as e:             # 自检自身出错也要给前端一个明确信号
        return jsonify({'ok': False, 'error': f'{type(e).__name__}: {e}',
                        'sections': []}), 500


@bp.route('/api/settings', methods=['POST'])
def api_settings_set():
    """设置中心保存：language/email/threads/defaults（逐项校验后持久化）。"""
    body = request.get_json(force=True) or {}
    from Virus_Platform_Core.config import Config, _ensure_not_system_dir
    known = {k: t for k, t, *_r in Config.DEFAULT_FIELDS}

    # ── ① 先整体校验，一个字节都不改 cfg ──
    # 原实现逐项就地写 cfg：中途某一项校验失败时，前面已生效的项既留在内存
    # 又可能已落盘（目录类 setter 内部各自 save 一次），配置处于"半套用"
    # 状态。现在先把全部入参校验/规整到 pending，全部通过才统一改 + 统一 save。
    pending = {}
    if 'language' in body:
        if body['language'] not in ('zh', 'en'):
            abort(400, 'language 仅支持 zh / en')
        pending['language'] = body['language']
    if 'email' in body:
        pending['email'] = str(body['email'] or '').strip()
    if 'threads' in body:
        try:
            pending['threads'] = max(1, int(body['threads']))
        except (TypeError, ValueError):
            abort(400, 'threads 必须是整数')
    if 'plugin_source_pref' in body:
        from Virus_Platform_Core.plugin_manager import SOURCE_PREFS
        v = str(body['plugin_source_pref'] or 'auto').strip().lower()
        if v not in SOURCE_PREFS:
            abort(400, 'plugin_source_pref 仅支持 ' + ' / '.join(SOURCE_PREFS))
        pending['plugin_source_pref'] = v
    if 'defaults' in body and isinstance(body['defaults'], dict):
        merged = dict(cfg.defaults)          # 基于现值合并，逐键覆盖
        for k, v in body['defaults'].items():
            if k not in known:
                continue
            kind = known[k]
            try:
                if kind == 'int':
                    merged[k] = int(v)
                elif kind == 'float':
                    merged[k] = float(v)
                elif kind == 'bool':
                    merged[k] = bool(v)
                else:
                    merged[k] = str(v)
            except (TypeError, ValueError):
                abort(400, f'参数 {k} 的值不合法: {v!r}')
        pending['defaults'] = merged

    # 目录类四项：set_*_root 会建目录/重写 DIRS（有副作用的写），先按同一套
    # 规则干校验（绝对路径 / 非系统关键目录 / 示例目录须存在 / 可创建），
    # 全部通过后才统一应用。
    root_setters = (
        ('output_root', cfg.set_output_root, '输出目录'),
        ('input_root', cfg.set_input_root, '输入目录'),
        ('database_root', cfg.set_database_root, '数据库目录'),
        ('examples_root', cfg.set_examples_root, '示例目录'),
    )
    roots_pending = []
    for key, _setter, label in root_setters:
        if key not in body:
            continue
        raw = str(body.get(key) or '').strip()
        if raw:
            if not os.path.isabs(raw):
                abort(400, f'{label}需为绝对路径')
            try:
                _ensure_not_system_dir(raw)
            except ValueError as e:
                abort(400, f'{label}不合法: {e}')
            if key == 'examples_root':
                if not os.path.isdir(raw):
                    abort(400, f'{label}不存在: {raw}')
            else:
                try:
                    os.makedirs(raw, exist_ok=True)
                except OSError as e:
                    abort(400, f'{label}无法创建: {e}')
        roots_pending.append((key, raw))

    # ── ② 统一应用（校验已全过；万一仍出错则整体回滚，不留半套用）──
    _snap = {'language': cfg.language, 'email': getattr(cfg, 'email', ''),
             'threads': cfg.threads, 'defaults': dict(cfg.defaults),
             'roots': {k: getattr(cfg, k, '') for k, _s, _l in root_setters}}

    def _restore():
        cfg.language = _snap['language']
        cfg.email = _snap['email']
        cfg.threads = _snap['threads']
        cfg.defaults = _snap['defaults']
        for k, setter, _l in root_setters:
            try:                     # 回滚失败（磁盘/权限）也只能尽力
                setter(_snap['roots'][k])
            except Exception:        # noqa: BLE001
                pass
        try:
            cfg.save()
        except OSError:
            pass

    _applying_roots = False
    try:
        for k, v in pending.items():
            setattr(cfg, k, v)
        _setter_map = {k: s for k, s, _l in root_setters}
        _applying_roots = True
        for key, raw in roots_pending:
            _setter_map[key](raw)        # setter 内部会 save 一次
        _applying_roots = False
        cfg.save()
    except (ValueError, OSError) as e:
        _restore()
        # 目录类失败多是用户输入不可用（不可创建/权限/系统目录）→ 400；
        # 落盘失败 → 500。两者都已整体回滚。
        abort(400 if _applying_roots else 500,
              f'配置保存失败（已回滚本次修改）: {e}')
    # 返回数据库就绪状态（前端切换 database_root 后可即时提示缺库）
    from Virus_Platform_Core.kunpeng import db_ready
    from Virus_Platform_Core.taxonomy import taxonomy_ready
    return jsonify({'ok': True, 'language': cfg.language,
                    'defaults': cfg.defaults,
                    'dbs': {'host': db_ready(cfg.databases['host']),
                            'virus': db_ready(cfg.databases['virus']),
                            'taxonomy': taxonomy_ready(),
                            'root': getattr(cfg, 'database_root', '')}})


def _bind_storage_scan():
    """把「存储占用扫描」挂进任务中心（外部任务）。

    该扫描是后台线程，走整棵目录树（含 databases/），可能跑几十秒；
    原先只在设置页显示一个 scanning 布尔，任务中心里完全看不见。
    扫描结束后 provider 返回 None，卡片自动消失。
    """
    from Virus_Platform_Core.web.tasks import tm

    def prov():
        with _storage_lock:
            if not _storage['scanning']:
                return None
            return {'status': 'running', 'pct': None, 'stage': '统计目录占用',
                    'msg': '正在统计各目录占用（结果会缓存，结束后本条自动消失）',
                    'log': []}

    tm.attach('存储占用扫描', prov, key='settings:storage-scan',
              link='/settings', weight='light')


@bp.route('/api/storage')
def api_storage():
    """磁盘与存储占用。

    ?detail=1 时返回目录占用明细；无缓存或缓存过期则后台起一次扫描，
    未就绪时 ready=false（前端可稍后再问），请求本身不阻塞。
    """
    detail = request.args.get('detail') in ('1', 'true', 'yes')
    disks = _disk_snapshot()
    resp = {'disks': disks,
            'warn': any(d.get('warn') for d in disks),
            'threshold_gb': STORAGE_WARN_GB,
            'scan': {'ready': False, 'scanning': False,
                     'scanned_at': None, 'items': [], 'error': None}}
    if not detail:
        return jsonify(resp)
    started = False
    with _storage_lock:
        fresh = (time.time() - _storage['ts']) < _STORAGE_TTL
        if not fresh and not _storage['scanning']:
            _storage['scanning'] = True
            threading.Thread(target=_do_storage_scan, daemon=True,
                             name='storage-scan').start()
            started = True
        resp['scan'] = {'ready': bool(_storage['items']) and fresh,
                        'scanning': _storage['scanning'],
                        'scanned_at': _storage['scanned_at'],
                        'items': _storage['items'] if fresh else [],
                        'error': _storage['error']}
    # 必须**在锁外**挂任务：_bind_storage_scan 会去拿 tm 的锁，而 _storage_lock
    # 是不可重入的普通 Lock —— 在里面挂就和 provider（快照时会拿 _storage_lock）
    # 形成了 storage→tm 的锁序，迟早咬上。
    if started:
        _bind_storage_scan()
    return jsonify(resp)


@bp.route('/api/storage/rescan', methods=['POST'])
def api_storage_rescan():
    """强制重新统计目录占用（忽略缓存）。"""
    with _storage_lock:
        if _storage['scanning']:
            return jsonify({'ok': False, 'msg': '统计进行中'})
        _storage['ts'] = 0
        _storage['items'] = []
        _storage['scanning'] = True
    threading.Thread(target=_do_storage_scan, daemon=True,
                     name='storage-scan').start()
    _bind_storage_scan()
    return jsonify({'ok': True})


# ------------------------------------------------------------------
# 插件管理（plugins/ 目录，机制见 plugin_manager.py）
# ------------------------------------------------------------------
def _plugin_describe_payload():
    """扫描结果 → 前端载荷（附当前生效路径与来源标记）。"""
    from Virus_Platform_Core import plugin_manager as pm
    items = []
    for it in pm.installed_plugins():
        items.append({
            'name': it['name'], 'kind': it['kind'],
            'version': it['version'], 'description': it['description'],
            'ok': it['ok'], 'error': it.get('error'),
            'map': {k: v for k, v in (it.get('map') or {}).items()},
            'active': {k: (cfg.plugin_tools.get(k) == v
                           or cfg.plugin_databases.get(k) == v)
                       for k, v in (it.get('map') or {}).items()},
        })
    return {'items': items, 'plugin_dir': pm.PLUGIN_DIR,
            'plugin_tools': cfg.plugin_tools,
            'plugin_databases': cfg.plugin_databases}


@bp.route('/api/plugin/list', methods=['GET'])
def api_plugin_list():
    return jsonify(_plugin_describe_payload())


@bp.route('/api/plugin/install', methods=['POST'])
def api_plugin_install():
    """安装插件 zip：multipart 上传（file=...）或 JSON {"path": 本地 zip 路径}。"""
    from Virus_Platform_Core import plugin_manager as pm
    tmp_zip = None
    try:
        f = request.files.get('file')
        if f is not None:
            os.makedirs(pm.PLUGIN_DIR, exist_ok=True)
            fd, tmp_zip = tempfile.mkstemp(prefix='.up-', suffix='.zip',
                                           dir=pm.PLUGIN_DIR)
            os.close(fd)
            f.save(tmp_zip)
            src = tmp_zip
        else:
            src = (request.json or {}).get('path') or ''
        if not src:
            return jsonify({'ok': False,
                            'error': '缺少插件包：请上传 zip（file=）或给出本地路径（path=）'}), 400
        mf = pm.install_zip(src)
        cfg.refresh_plugins()
        resp = _plugin_describe_payload()
        resp['ok'] = True
        resp['installed'] = mf.get('name')
        return jsonify(resp)
    except pm.PluginError as e:
        return jsonify({'ok': False, 'error': str(e)}), 400
    except Exception as e:
        return jsonify({'ok': False,
                        'error': f'{type(e).__name__}: {e}'}), 500
    finally:
        if tmp_zip:
            try:
                os.remove(tmp_zip)
            except OSError:
                pass


@bp.route('/api/plugin/available', methods=['GET'])
def api_plugin_available():
    """插件源索引 × 本机现状：自动检测缺失工具，标注 needed/missing_keys。

    ?force=1 跳过索引缓存。索引 URL 取 cfg.plugin_index_url（空 = 内置默认，
    即本仓库 plugins-repo/index.json 的 raw 地址）。"""
    from Virus_Platform_Core import plugin_manager as pm
    force = request.args.get('force') in ('1', 'true')
    url = (cfg.plugin_index_url or '').strip() or None
    items, err = pm.available_from_index(cfg, index_url=url, force=force)
    if err:
        return jsonify({'ok': False, 'error': err,
                        'index_url': url or pm.DEFAULT_INDEX_URL}), 502
    return jsonify({'ok': True, 'items': items,
                    'index_url': url or pm.DEFAULT_INDEX_URL})


@bp.route('/api/plugin/install_url', methods=['POST'])
def api_plugin_install_url():
    """从 URL 下载安装（对齐 PhyloSuite 的镜像直装）。下载失败时提示用户
    手动下载后走 /api/plugin/install 上传或指定本地路径。
    JSON 可带 sha256（索引条目自带）→ 下载后先验哈希再解压。"""
    from Virus_Platform_Core import plugin_manager as pm
    body = request.json or {}
    # urls = 按顺序回退的地址列表（国内镜像在前）；url = 单个地址（兼容旧调用）
    urls = body.get('urls') or ([body['url']] if body.get('url') else [])
    sha = body.get('sha256') or None
    try:
        mf = pm.install_urls(urls, expected_sha256=sha)
        cfg.refresh_plugins()
        resp = _plugin_describe_payload()
        resp['ok'] = True
        resp['installed'] = mf.get('name')
        return jsonify(resp)
    except pm.PluginError as e:
        return jsonify({'ok': False, 'error': str(e)}), 400
    except Exception as e:
        return jsonify({'ok': False, 'error': f'{type(e).__name__}: {e}',
                        'hint': '下载失败可手动下载 zip 后在「本地安装」里指定'}), 502


@bp.route('/api/plugin/modules', methods=['GET'])
def api_plugin_modules():
    """模块视图：每个功能模块依赖哪些工具包、本机还缺哪些。

    ?ids=a,b 只算指定模块；?force=1 跳过索引缓存。返回里 plan 已去重
    （同一工具被多模块共用只算一次）并跳过本机已有/已装。"""
    from Virus_Platform_Core import plugin_manager as pm
    force = request.args.get('force') in ('1', 'true')
    ids = [x for x in (request.args.get('ids') or '').split(',') if x.strip()]
    url = (cfg.plugin_index_url or '').strip() or None
    try:
        items, err = pm.available_from_index(cfg, index_url=url, force=force)
        if err:
            return jsonify({'ok': False, 'error': err}), 502
        plan = pm.plan_module_install(cfg, ids, items=items,
                                      index_url=url, force=force)
        # 每个模块单独标注「还需装几个 / 共几个」
        by_name = {p['name']: p for p in items}
        for m in plan['modules']:
            names = m.get('plugins') or []
            m['total'] = len(names)
            m['missing'] = [n for n in names
                            if (by_name.get(n) or {}).get('needed')
                            and not (by_name.get(n) or {}).get('installed')]
            m['ready'] = not m['missing']
        return jsonify({'ok': True, 'index_url': url or pm.DEFAULT_INDEX_URL,
                        **plan})
    except Exception as e:
        return jsonify({'ok': False, 'error': f'{type(e).__name__}: {e}'}), 500


@bp.route('/api/plugin/install_modules', methods=['POST'])
def api_plugin_install_modules():
    """按模块安装（自动去重 + 跳过已有）。JSON: {"modules": ["tree", ...]}"""
    from Virus_Platform_Core import plugin_manager as pm
    ids = (request.json or {}).get('modules') or []
    if not isinstance(ids, list) or not ids:
        return jsonify({'ok': False, 'error': '未指定模块'}), 400
    url = (cfg.plugin_index_url or '').strip() or None
    try:
        r = pm.install_modules(cfg, ids, index_url=url)
        resp = _plugin_describe_payload()
        resp.update({'ok': True, **r})
        return jsonify(resp)
    except pm.PluginError as e:
        return jsonify({'ok': False, 'error': str(e)}), 400
    except Exception as e:
        return jsonify({'ok': False,
                        'error': f'{type(e).__name__}: {e}'}), 500


@bp.route('/api/plugin/uninstall', methods=['POST'])
def api_plugin_uninstall():
    from Virus_Platform_Core import plugin_manager as pm
    name = (request.json or {}).get('name') or ''
    try:
        pm.uninstall(name)
        cfg.refresh_plugins()
        resp = _plugin_describe_payload()
        resp['ok'] = True
        resp['uninstalled'] = name
        return jsonify(resp)
    except pm.PluginError as e:
        return jsonify({'ok': False, 'error': str(e)}), 400
    except Exception as e:
        return jsonify({'ok': False,
                        'error': f'{type(e).__name__}: {e}'}), 500
