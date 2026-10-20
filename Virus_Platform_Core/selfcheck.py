# -*- coding: utf-8 -*-
"""环境自检（结构化版）：模块导入 / 外部工具 / 数据库 / 磁盘空间。

collect() 返回四段结构化结果，供两处共用：
  - web API /api/selfcheck（设置页「环境自检」卡片 + 总览页环境状态条）
  - main.py selfcheck（CLI 控制台版，另带 ⑤样品 / ⑥依赖 两段 CLI 专属检查）

结果按进程缓存（CACHE_TTL 秒内直接复用）：首次收集要递归 import 全部 80+
子模块，实测数秒；此后模块都已进 sys.modules，force=True 重跑也很快。
整体 ok 的口径与 CLI 一致：模块导入与数据库失败才算失败；工具缺失只算
降级（warn），磁盘不足只提示。
"""
import importlib
import shutil
import time

from Virus_Platform_Core.config import DIRS, PLATFORM_ROOT, get_config

CACHE_TTL = 300.0          # 缓存有效期（秒）；设置页按钮走 force=1 强制重跑
DISK_WARN_BYTES = 10e9     # 磁盘剩余低于此值提示（与 CLI 自检同一口径）
_CACHE = {'ts': 0.0, 'data': None}


def invalidate():
    """清掉进程内自检缓存，下次 collect() 立即重算。

    数据库状态在进程内发生变化（建库完成 / 加载宿主库 / 切换当前库 /
    登记病毒库）后应由改动方调一下：否则总览页环境状态条会把加载前的
    "数据库 未就绪"一直展示到缓存过期，看起来像"加载了但没同步"
    （2026-09-30 手动加载宿主库后即此现象）。模块探测较贵（子进程逐个
    import），所以不做短 TTL 轮询，只在状态变化点精确失效。"""
    _CACHE['ts'] = 0.0
    _CACHE['data'] = None


def iter_platform_modules():
    """枚举 Virus_Platform_Core 下的全部子模块（递归含子包）。

    返回 (模块名列表, 遍历期导入失败的子包名列表)。

    为什么遍历而不是硬编码清单：原先自检只 import 16 个写死的模块，而平台有
    80+ 个。**走子进程调用的引擎不在其中** —— kv 引擎由 kv_stage 以子进程启动
    （`--run-engine`），kv_stage 本身只 import config/utils，于是 kv_filter /
    kv_identify 在模块级的 `import polars` 缺失时，自检完全看不见。
    2026-09-11 就是这个盲区让「已知病毒识别与定量」整段静默失效。
    （2026-10-20 起该模块已改用 pandas + numpy，不再依赖 polars；
     这条历史教训仍成立——子进程调用的引擎必须在自检里覆盖到。）

    冻结分发下靠 PyInstaller 的 `pyi_rth_pkgutil` 运行时钩子让 pkgutil 可用，
    打包后同样能遍历。遍历本身若失败（返回空列表），调用方回退到手工清单。
    """
    import pkgutil
    mods, errs = [], []
    try:
        import Virus_Platform_Core as _pkg
        # 逐子包显式导入（带追踪）：pkgutil.walk_packages 遍历时会导入各
        # 子包 __init__，那里崩了（冻结版 C 扩展 0xC0000005）连一行痕迹
        # 都没有 —— 先自己按序导入并打印，肇事包一目了然。
        import os as _os
        import sys as _sys
        _trace = _os.environ.get('VP_IMPORT_TRACE') == '1'
        _subpkgs = []
        for _e in _pkg.__path__:
            for _n in sorted(_os.listdir(_e)):
                _p = _os.path.join(_e, _n)
                if _n.startswith('_'):
                    continue
                if _os.path.isfile(_os.path.join(_p, '__init__.py')):
                    _subpkgs.append(_pkg.__name__ + '.' + _n)
                elif _n.endswith('.py'):
                    _subpkgs.append(_pkg.__name__ + '.' + _n[:-3])
        for _full in _subpkgs:
            if _trace:
                print('[pkgimport] ' + _full, file=_sys.stderr, flush=True)
            try:
                importlib.import_module(_full)
            except Exception as e:
                errs.append(_full)
                if _trace:
                    print('[pkgimport FAILED] ' + repr(e), file=_sys.stderr,
                          flush=True)
        for info in pkgutil.walk_packages(
                _pkg.__path__, _pkg.__name__ + '.',
                onerror=lambda n: errs.append(n)):
            if info.ispkg:
                continue
            # 跳过**私有脚本**（`_smoke_test` 这类）。它们是"跑一次"的脚本，
            # 不是平台模块，且模块级代码带副作用 —— explorer/_smoke_test.py
            # 一被 import 就会完整跑一遍 explorer 冒烟测试（加载 19.9 万条
            # 元数据 + 457MB FASTA，冷启动分钟级）。自检本意只是"导入看能不能
            # 进来"，却因此把整个冒烟测试跑了一遍：既拖慢自检，又在 explorer
            # 数据移出 databases/ 后让这一项必然报错（2026-09-17 查实，调用栈
            # 为 api_selfcheck → _check_modules → importlib.import_module）。
            # `__init__.py` 不在 walk_packages 的返回里；这里仍显式放行，
            # 以免将来遍历实现变化时误伤包初始化。
            _leaf = info.name.rsplit('.', 1)[-1]
            if _leaf.startswith('_') and _leaf != '__init__':
                continue
            mods.append(info.name)
    except Exception as e:
        errs.append(f'<遍历失败: {e!r}>')
    return sorted(set(mods)), errs


# 遍历不可用时的兜底清单（原 CLI 内联清单，保持一致）
FALLBACK_MODULES = [
    'Virus_Platform_Core.config', 'Virus_Platform_Core.utils',
    'Virus_Platform_Core.taxonomy', 'Virus_Platform_Core.kunpeng',
    'Virus_Platform_Core.preprocess',
    'Virus_Platform_Core.host_removal',
    'Virus_Platform_Core.kv_stage', 'Virus_Platform_Core.assembly',
    'Virus_Platform_Core.host_analysis', 'Virus_Platform_Core.orf',
    'Virus_Platform_Core.orf_annot', 'Virus_Platform_Core.phylo',
    'Virus_Platform_Core.primer', 'Virus_Platform_Core.viz',
    'Virus_Platform_Core.pipeline', 'Virus_Platform_Core.msa_view',
]


def _check_modules():
    """逐模块导入探测。**必须走子进程**（importprobe 子命令）：

    冻结分发里某个 C 扩展的 import 可能直接段错误（0xC0000005，进程即死、
    无任何 Python 异常）。在 web 服务进程内做导入扫描时，一个坏模块就能把
    整个平台拖死（2026-09-19 实测）；放子进程里跑，段错误最多杀死探测进程，
    父进程靠协议行数定位肇事模块，平台自身毫发无损。
    """
    mods, walk_errs = iter_platform_modules()
    if not mods:
        mods = FALLBACK_MODULES
    total = len(mods) + len(walk_errs)
    t0 = time.time()
    import os
    import subprocess as _sub
    import sys as _sys
    if getattr(_sys, 'frozen', False):
        argv = [_sys.executable, '--cli', 'importprobe']
    else:
        argv = [_sys.executable, os.path.join(PLATFORM_ROOT, 'main.py'),
                'importprobe']
    mods_csv = ','.join(mods)
    fails = []
    crashed_at = None
    seen_p = 0
    _env = os.environ.copy()
    _env['VP_IMPORT_TRACE'] = '1'      # 子进程 stderr 留下逐包导入痕迹，崩溃可溯源
    try:
        # 显式 encoding='utf-8'：text=True 缺省走系统 ANSI 代码页（中文
        # Windows = GBK），探测子进程输出含非 GBK 字节时会解码报错/乱码
        proc = _sub.run(argv + [mods_csv], capture_output=True,
                        text=True, encoding='utf-8', errors='replace',
                        timeout=900, env=_env)
    except Exception as e:
        fails.append(('<importprobe>', f'探测子进程启动失败: {e!r}'))
        proc = None

    if proc is not None:
        done = False
        pending = None
        for ln in (proc.stdout or '').splitlines():
            if ln.startswith('P '):
                pending = ln[2:]
                seen_p += 1
            elif ln.startswith('F '):
                _m, _, _e = ln[2:].partition(' ')
                fails.append((_m, _e))
                pending = None
            elif ln.startswith('E '):
                fails.append((ln[2:], '子包导入失败，其下模块未能枚举'))
                pending = None
            elif ln.startswith('DONE'):
                done = True
                pending = None
        if pending:
            crashed_at = pending
        elif not done and seen_p == 0:
            crashed_at = '<枚举阶段>'
        if crashed_at:
            # 探测子进程被原生段错误杀死：肇事模块 + 其后未探测的都标失败
            _err_tail = ((proc.stderr or '').strip().splitlines() or [''])[-1]
            fails.append((crashed_at,
                          f'导入时进程崩溃（原生段错误）· 崩溃前最后痕迹: '
                          f'{_err_tail[:120]}'))
            idx = mods.index(crashed_at) if crashed_at in mods else 0
            for m in mods[idx + 1:]:
                fails.append((m, '未探测（探测进程在先前的模块崩溃后中止）'))

    return {'id': 'modules', 'title': '模块导入',
            'ok': not fails, 'warn': False,
            'summary': f'{total - len(fails)}/{total}',
            'seconds': round(time.time() - t0, 1),
            # 只列失败明细：86 个模块名全量回传没有信息量
            'items': [{'name': m, 'ok': False, 'detail': e}
                      for m, e in fails]}


def _check_tools():
    cfg = get_config()
    st = cfg.tool_status()
    missing = sorted(k for k, v in st.items() if not v)

    # 插件源能不能补：把缺失键和插件索引对一遍，能补的标出来。
    # 用户看到「缺 X」时最想知道的是「点一下能不能好」，而不是自己去
    # 翻插件列表找对应关系。索引不可达时降级为不标注（不报错）。
    fixable, unfixable = {}, []
    try:
        from Virus_Platform_Core import plugin_manager as pm
        items, _err = pm.available_from_index(cfg)
        if items:
            for it in items:
                prov = it.get('provides') or {}
                hit = [k for k in missing if k in prov]
                if hit:
                    fixable[it['name']] = {
                        'provides': hit,
                        'size': it.get('size', 0),
                        'installed': it.get('installed', False)}
            covered = {k for v in fixable.values() for k in v['provides']}
            unfixable = [k for k in missing if k not in covered]
        else:
            unfixable = list(missing)
    except Exception:
        unfixable = list(missing)

    # 本次会话自动装过什么（module_deps 落盘）——插件被悄悄装上时用户
    # 必须看得见，否则会以为"平台自己改了我的环境"。
    autofixed = {}
    try:
        from Virus_Platform_Core.module_deps import read_autofix_log
        autofixed = read_autofix_log()
    except Exception:
        pass

    return {'id': 'tools', 'title': '外部工具',
            'ok': not missing, 'warn': bool(missing),
            'summary': f'{len(st) - len(missing)}/{len(st)}',
            'missing': missing,
            # 插件源可一键补齐的（前端据此显示「可自动安装」按钮）
            'fixable': fixable,
            'unfixable': unfixable,
            'autofixed': autofixed}
    # 工具逐项明细不重复给：设置页上方「工具探测」卡片已有同源列表


def _check_databases():
    from Virus_Platform_Core.config import host_db_info
    from Virus_Platform_Core.kunpeng import db_ready
    from Virus_Platform_Core.taxonomy import taxonomy_ready
    cfg = get_config()
    items = [{'name': 'NCBI Taxonomy', 'ok': bool(taxonomy_ready())}]
    for key, label in (('host', '宿主库 host_db'),
                       ('virus', '病毒库 virus_db')):
        ok = bool(db_ready(cfg.databases[key]))
        item = {'name': label, 'ok': ok, 'detail': cfg.databases[key]}
        if key == 'host':
            # 宿主库"是谁"必须可见：库本体不记录物种，历史上出过
            # "目录名是枸杞、库却按 4081 番茄建"的事故，用错宿主是静默的。
            h = host_db_info(cfg.databases[key])
            if h:
                who = (f"{h.get('species')} (taxid={h.get('taxid')})"
                       if h.get('taxid') else '身份未记录')
                item['detail'] = f"{h.get('name')} | 物种 = {who}"
                if h.get('conflicted'):
                    item['detail'] += (' | ⚠ 元数据冲突 '
                                       f"{h.get('taxids_in_map')}")
        items.append(item)
    return {'id': 'databases', 'title': '数据库',
            'ok': all(i['ok'] for i in items), 'warn': False,
            'summary': f"{sum(i['ok'] for i in items)}/{len(items)}",
            'items': items}


def _check_disk():
    items = []
    for key, d in (('平台根', PLATFORM_ROOT), ('结果', DIRS['results'])):
        try:
            free = shutil.disk_usage(d).free
            items.append({'name': key, 'ok': free > DISK_WARN_BYTES,
                          'detail': f'剩余 {free / 1e9:.1f} GB'
                                    + ('' if free > DISK_WARN_BYTES
                                       else '（建议保留 >10GB）')})
        except OSError as e:
            items.append({'name': key, 'ok': False, 'detail': str(e)})
    return {'id': 'disk', 'title': '磁盘空间',
            'ok': all(i['ok'] for i in items),
            'warn': any(not i['ok'] for i in items),
            'summary': f"{sum(i['ok'] for i in items)}/{len(items)}",
            'items': items}


def collect(force=False):
    """跑一遍四段自检并返回结构化 dict；进程内缓存 CACHE_TTL 秒。

    返回形如：{'ok': bool, 'seconds': 3.2, 'cached': False,
               'sections': [ {id,title,ok,warn,summary,items?}, ... ]}
    """
    now = time.time()
    if (not force and _CACHE['data'] is not None
            and now - _CACHE['ts'] < CACHE_TTL):
        data = dict(_CACHE['data'])
        data['cached'] = True
        return data
    t0 = time.time()
    sections = [_check_modules(), _check_tools(),
                _check_databases(), _check_disk()]
    # 整体口径与 CLI 一致：模块导入与数据库失败才算失败（工具缺失=降级、
    # 磁盘不足=提示，都不算失败）
    modules, databases = sections[0], sections[2]
    data = {'ok': bool(modules['ok'] and databases['ok']),
            'seconds': round(time.time() - t0, 1),
            'ts': now, 'cached': False, 'sections': sections}
    _CACHE['ts'], _CACHE['data'] = now, data
    return data
