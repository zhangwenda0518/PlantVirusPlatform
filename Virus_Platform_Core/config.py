# -*- coding: utf-8 -*-
"""
平台配置：工具路径自动探测 / 手动覆盖、目录布局、线程数。
配置持久化到 <平台根>/platform.json，GUI 修改后立即生效。
"""
import os
import sys
import json
import glob
import shutil
import threading

if sys.platform == 'win32':
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
    sys.stderr.reconfigure(encoding='utf-8', errors='replace')


def _detect_platform_root():
    """平台根目录：PyInstaller 打包后 = exe 所在目录；源码运行 = Virus_Platform_Core/ 上级目录。"""
    if getattr(sys, 'frozen', False):          # PyInstaller 环境
        return os.path.dirname(os.path.abspath(sys.executable))
    return os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


PLATFORM_ROOT = _detect_platform_root()


def engine_cmd(script_abspath, *args):
    """构造「以子进程运行内部引擎脚本」的命令行。

    源码运行: [python, script.py, args...]
    冻结分发: [VirusPlatform.exe, --run-engine, script.py, args...]
    （打包环境无独立 python、无源码脚本路径；app.py 入口识别 --run-engine
    后经 Virus_Platform_Core/engine_entry.py 进程内执行对应引擎。）
    """
    if getattr(sys, 'frozen', False):
        return [sys.executable, '--run-engine',
                os.path.basename(script_abspath)] + list(args)
    # -u：子进程 stdout 无缓冲（管道默认块缓冲会让任务日志长时间空白）
    return [sys.executable, '-u', script_abspath] + list(args)

# ------------------------------------------------------------------
# 目录布局（全部限定在平台根目录内）
# ------------------------------------------------------------------
# 病毒参考库根 = databases/virusref_db（2026-09-11 从顶层 virus-db/ 迁入，
# 并按 <name>_db 惯例定名）。里面是运行期要用的病毒参考与派生索引：
#   <库名>/  —— 自包含鉴定库：reference.fasta + reference.ref_info.tsv +
#               salmon_k31/（定量唯一引擎 salmon；默认库目录名 plant_ref，
#               2026-09-29 前为 kv_index）
#   blast/   —— 本地检索用的预置 BLAST 库
# 2026-09-12 重构：根目录散置的 final.cluster.ref.* 已归档进 kv_index/
# （config.ref_annotation_path / kv_stage.default_reference 兼容回退旧布局）。
# 注：命名避开 virus_db —— 那是 kunpeng 植物病毒分类库的历史目录名。
# 目录已定局，不再保留旧布局兜底（virus_src / virus-db 均已废弃）。


DIRS = {
    'databases': os.path.join(PLATFORM_ROOT, 'databases'),   # kunpeng 库根目录
    'taxonomy':  os.path.join(PLATFORM_ROOT, 'databases', 'taxonomy'),
    'host_src':  os.path.join(PLATFORM_ROOT, 'host-db'),     # 宿主源基因组
    'virus_src': os.path.join(PLATFORM_ROOT, 'databases',
                              'virusref_db'),               # 病毒参考 + 派生索引
    # ---- 运行期数据统一收在 run/ 下（2026-09-10 目录整理）----
    'results':   os.path.join(PLATFORM_ROOT, 'run', 'results'),   # 样品结果
    'tool_runs': os.path.join(PLATFORM_ROOT, 'run', 'tool_runs'),  # 工具运行目录
    'downloads': os.path.join(PLATFORM_ROOT, 'run', 'downloads'),  # 公共数据下载
    'logan':     os.path.join(PLATFORM_ROOT, 'run', 'logan'),     # LOGAN 溯源任务
    'submissions': os.path.join(PLATFORM_ROOT, 'run', 'submissions'),  # NCBI 提交准备
    'meta_search': os.path.join(PLATFORM_ROOT, 'run', 'meta_search'),  # 公共数据检索
    'tasks':     os.path.join(PLATFORM_ROOT, 'run', 'tasks'),     # GUI 任务状态
    'logs':      os.path.join(PLATFORM_ROOT, 'run', 'logs'),      # 全局日志
    'ai_audit':  os.path.join(PLATFORM_ROOT, 'run', 'ai_audit'),  # AI 代理审计留痕
    'uploads':   os.path.join(PLATFORM_ROOT, 'run', 'uploads'),   # 粘贴/上传中转
    'fastq':     os.path.join(PLATFORM_ROOT, 'run', 'fastq'),     # 流程 FASTQ 中转
    # ---- 用户在平台内产生的「集合」类数据（不是分发的库）----
    'gb_collections': os.path.join(PLATFORM_ROOT, 'run', 'gb_collections'),
    'ncbi_refs':      os.path.join(PLATFORM_ROOT, 'run', 'ncbi_refs'),
    'webapp_static': os.path.join(PLATFORM_ROOT, 'webapp', 'static'),
    # 示例数据根（内置示例 FASTA/GenBank/树）：默认在程序目录内，
    # 但允许外置（打包时「程序 / 数据库 / 示例」三分离）。
    'examples':  os.path.join(PLATFORM_ROOT, 'examples'),
}

# 目录整理（2026-09-10）后的新旧顶层目录映射：
#   运行期数据 → run/，外部二进制与第三方 → 3rd/。
# 历史数据/前端/文档里仍可能写旧相对路径（tool_runs/xxx、results/xxx…），
# utils.check_path 会按这张表重定向，保证旧调用与旧记录不失效。
LEGACY_TOP_DIRS = {
    'results': 'run/results',
    'tool_runs': 'run/tool_runs',
    'logs': 'run/logs',
    'tasks': 'run/tasks',
    'uploads': 'run/uploads',
    'submissions': 'run/submissions',
    'meta_search': 'run/meta_search',
    'logan': 'run/logan',
    'fastq': 'run/fastq',
    'downloads': 'run/downloads',
    'tools': '3rd/tools',
    'bin': '3rd/bin',
    'vendor': '3rd/vendor',
    'open-virome': '3rd/open-virome',
}

# 多段前缀的旧→新映射（单段映射见上一张表）。
# 这些前缀随目录整理整体搬迁，历史记录 / 旧前端常量里仍可能写旧写法：
#   databases/examples/   2026-09-10 迁到 <平台根>/examples/
#   databases/misc/gb/    同上，与 ncbi_refs 一起迁入 run/（用户数据，
#   databases/ncbi_refs/  不是平台分发的库，不该占 databases/）
LEGACY_MULTI_DIRS = {
    'databases/examples':  'examples',
    'databases/misc/gb':   'run/gb_collections',
    'databases/ncbi_refs': 'run/ncbi_refs',
}

# ------------------------------------------------------------------
# 数据库分类路径注册表
# ------------------------------------------------------------------
# 每个条目 = 相对**数据库根**的子路径；只有宿主分类库用
# (子路径, base_key) 二元组 —— 它的根是 host-db/（按物种而异、单库约
# 1.5GB，不进 databases/ 分发包），其余都在 databases/ 下。
# 2026-09-11 目录整理定局，不再保留旧布局候选与 old_name 兜底：
#   annot_db / kunpeng_db / misc_db / tax_db / tree_db / virusref_db
# 用户数据（gb_collections / ncbi_refs）**不在此表** —— 它们不是平台
# 分发的库，已归入 run/，见上方 DIRS。
# 当前宿主分类库（绝对路径；'' = 按 active_host_db / 自动探测 / 旧槽位解析）。
# 提前声明：db_path() 里要读它，而本文件顶部有模块级 db_path() 调用。
_ACTIVE_HOST_DB = ''
DB_LAYOUT = {
    'virus': {
        'plant': 'kunpeng_db/plant',
        'ref':   'kunpeng_db/ref',
        # RVDB（C-RVDB 通用参考库）槽位：universal_ref 仍会按需建/取它
        'rvdb':  'kunpeng_db/rvdb',
    },
    'host': {
        # 宿主**分类库**（kunpeng hash）：host-db/host/classify
        'classify': ('host/classify', 'host_src'),
    },
    'annot': {
        # cdd 槽位已移除：本地 CDD 库 prefilter 索引 ~5GB 反复 OOM，
        # 2026-09-29 随功能一并下线
        'hmm':  'annot_db/hmm',
        'prot': 'annot_db/prot',
    },
    'tax': {
        # taxonomy 三件套（nodes/names/merged.dmp）扁平放在 databases/tax_db/
        'core': 'tax_db',
        # VMR 参考库与 plant 建树库并列住 tree_db/
        'ictv': 'tree_db/ictv_tree.db',
    },
    # 建树参考库（两套口径，同住 tree_db/，由库本身决定宿主范围）
    #   plant_tree.db : 植物病毒口径，**预下载**库
    #   ictv_tree.db  : 全病毒界口径，**按需下载**（VMR 谱系 + gb_cache）
    #                   （见 Virus_Platform_Core/acvirus.py 文件头「历史脉络」）
    'tree': {
        'plant': 'tree_db/plant_tree.db',
        'ictv':  'tree_db/ictv_tree.db',
    },
    'misc': {
        'viroids': 'misc_db/viroids',
        'suvtk':   'misc_db/suvtk',
        'prob':    'misc_db/prob',
    },
}

def _resolve_host_db():
    """当前宿主库，不依赖初始化顺序。

    为什么不能直接读 _ACTIVE_HOST_DB：它只有 Config.__init__ 跑过才填充，而
    db_path('host','classify') 既可能被模块级/脚本提前调用，也在 Config.__init__
    内部（apply_active_host_db 之前）被调用 —— 那时全局还是空串，会静默回落到
    旧槽位 host/classify（实测：删掉旧库后 db_path 仍打印 'classify'）。

    为什么不用 get_active_host_db()（它带惰性构造配置）：get_config() 是
    `_config = Config()`，构造完成才赋值，重入会**再构造一个 Config**。这里改用
    纯文件系统探测（扫描 host-db/*_host_db），无构造、无递归。
    """
    return _ACTIVE_HOST_DB or detect_active_host_db()


def db_path(category, name):
    """按分类取数据库路径（kunpeng_db/plant、annot_db/hmm 等）。

    单一路径解析：注册表已是当前布局，不做旧目录兜底（目录不存在时也返回
    该路径，供建库落盘用）。

    例外：宿主分类库按物种分目录（host-db/<taxid>_<slug>_host_db/），由
    platform.json 的 active_host_db 决定当前用哪个；未配置时自动探测，
    都没有才回落到注册表里的旧槽位 host/classify。
    """
    if category == 'host' and name == 'classify':
        resolved = _resolve_host_db()
        if resolved:
            return resolved
    entry = DB_LAYOUT[category][name]
    sub, base_key = entry if isinstance(entry, tuple) else (entry, 'databases')
    return os.path.join(DIRS[base_key], sub)


# taxonomy 目录也走注册表（databases/tax_db，三件套扁平放置）
DIRS['taxonomy'] = db_path('tax', 'core')


def ref_annotation_path():
    """参考病毒注释 TSV 的当前路径（accession → 物种/taxid/科/属）。

    注释随参考归档进默认鉴定库目录（2026-09-29 起为
    virusref_db/plant_ref/reference.ref_info.tsv，此前为 kv_index/）；
    旧布局在 virusref_db 根下的 final.cluster.ref_info.tsv 作为兼容回退。
    找不到时返回新布局期望路径，由调用方判存在。
    """
    base = DIRS['virus_src']
    for name in (os.path.join('plant_ref', 'reference.ref_info.tsv'),
                 os.path.join('kv_index', 'reference.ref_info.tsv'),
                 'final.cluster.ref_info.tsv'):
        p = os.path.join(base, name)
        if os.path.isfile(p):
            return p
    return os.path.join(base, 'plant_ref', 'reference.ref_info.tsv')


def host_source_dirs():
    """host-db 下所有宿主源目录（<taxid>_<slug>/ 且含 genome.fa）。"""
    base = DIRS['host_src']
    if not os.path.isdir(base):
        return []
    out = []
    for n in sorted(os.listdir(base)):
        p = os.path.join(base, n)
        if os.path.isdir(p) and os.path.isfile(os.path.join(p, 'genome.fa')):
            out.append(n)
    return out


def current_host_genome():
    """当前宿主源基因组。

    优先级：当前宿主库清单里记录的源基因组（host_db.json source_genome）
    → host-db 下第一个含 genome.fa 的宿主目录 → 旧默认 host-db/genome.fa。

    为什么要先看清单：宿主源目录可以有多个（多宿主共存），而"当前在用的
    宿主库"是显式选定的；跟随当前库的源基因组才不会出现"库是枸杞、回填的
    基因组却是另一个物种"。
    """
    m = host_db_manifest(get_active_host_db())
    src = str(m.get('source_genome') or '').strip()
    if src and os.path.isfile(src):
        return src
    dirs = host_source_dirs()
    if dirs:
        return os.path.join(DIRS['host_src'], dirs[0], 'genome.fa')
    old = os.path.join(DIRS['host_src'], 'genome.fa')
    return old if os.path.isfile(old) else None


# ------------------------------------------------------------------
# 宿主分类库（kunpeng hash）：host-db/<taxid>_<slug>_host_db/
# ------------------------------------------------------------------
# 历史上是单一固定槽位 host-db/host/classify，带来三个问题：
#   ① 只能存在一个宿主库 —— 建第二个宿主会覆盖/混进同一目录；
#   ② 库本身不记录物种，用错宿主去宿主去除是**静默**的；
#   ③ kunpeng 的 add-library 往 seqid2taxid.map 追加，而 build_db 的替换式
#      清理漏了该文件与 taxo.k2d → 陈旧 seqid→taxid 长期累积。
#      实测该库的 map 里同时有 4081（番茄 Solanum lycopersicum）与
#      112863（枸杞 Lycium barbarum）两个 taxid，而日志显示 2026-09-03 那次
#      建库把 689 条序列（= 枸杞基因组的序列数）注成了 taxid=4081。
# 现在改为：目录按物种命名 + 库内写 host_db.json 清单 + platform.json 的
# active_host_db 显式指定"当前用哪个"。求值优先级：
#   显式 active_host_db → 自动探测（唯一 / 最近构建的 *_host_db）
#   → 旧槽位 host/classify（兼容未迁移的部署）
_HOST_DB_SUFFIX = '_host_db'
# _ACTIVE_HOST_DB 在文件上方 db_path() 之前声明（那里要读它）
_PATH_BAD = set('<>:"/\\|?*') | {chr(_i) for _i in range(32)}


def host_db_manifest(db_dir):
    """读库清单 host_db.json（不存在/损坏一律返回 {}）。"""
    if not db_dir:
        return {}
    p = os.path.join(db_dir, 'host_db.json')
    try:
        with open(p, encoding='utf-8') as f:
            d = json.load(f)
        return d if isinstance(d, dict) else {}
    except (OSError, ValueError):
        return {}


def host_db_taxids(db_dir):
    """seqid2taxid.map 里出现的全部 taxid。

    正常应当**只有一个**（宿主物种本身）。出现多个 = 库的元数据被历史建库
    污染（见上方背景），调用方应据此告警而不是照跑。
    """
    out = set()
    if not db_dir:
        return out
    p = os.path.join(db_dir, 'seqid2taxid.map')
    if not os.path.isfile(p):
        return out
    try:
        with open(p, encoding='utf-8', errors='replace') as f:
            for line in f:
                parts = line.rstrip('\n').split('\t')
                if len(parts) >= 2 and parts[1].strip().isdigit():
                    out.add(int(parts[1]))
    except OSError:
        pass
    return out


def host_db_species(db_dir):
    """库对应的 (taxid, 物种拉丁名)。优先清单，其次 map 里唯一 taxid 反查。

    返回 (None, '') 表示无法判定（多 taxid 且无清单 = 身份不明）。
    """
    m = host_db_manifest(db_dir)
    if m.get('taxid'):
        try:
            return int(m['taxid']), str(m.get('species') or '')
        except (TypeError, ValueError):
            pass
    taxids = host_db_taxids(db_dir)
    if len(taxids) == 1:
        t = next(iter(taxids))
        name = ''
        try:
            from .taxonomy import TaxonomyGraph
            name = TaxonomyGraph().name(t)
        except Exception:
            name = ''
        return t, name
    return None, ''


def host_db_name(taxid, genome_fasta=''):
    """宿主库目录名：<源目录名或 taxid>_host_db。

    host-db/112863_Lycium_barbarum/genome.fa → 112863_Lycium_barbarum_host_db
    基因组不在 <taxid>_<slug>/ 里时退化为 <taxid>_host_db。
    """
    slug = ''
    if genome_fasta:
        parent = os.path.basename(os.path.dirname(os.path.abspath(str(genome_fasta))))
        if parent and parent.lower() not in ('host-db', 'host_db', 'fastq', ''):
            slug = parent
    if slug:
        base = slug if slug.startswith(f'{taxid}_') else f'{taxid}_{slug}'
    else:
        base = str(taxid)
    base = ''.join('_' if c in _PATH_BAD else c for c in base).strip('_')
    return f'{base}{_HOST_DB_SUFFIX}'


def legacy_host_db_dir():
    """旧布局的固定槽位（迁移前部署仍在用）。"""
    return os.path.join(DIRS['host_src'], 'host', 'classify')


def host_db_dirs():
    """host-db 与已登记外部库根下所有已建成的宿主库目录（含 hash_*.k2d）。

    外部库根（extra_db_roots，"加载已构建宿主库"登记平台外库时自动登记
    其父目录）也参与扫描：否则库在平台外时本函数返回空，建库页会一边
    显示"当前宿主库 ✓"（那走 databases.host 的 db_ready）一边又提示
    "还没有宿主库"（2026-09-30 实测）。外部根可能是病毒库的父目录，
    只认 *_host_db 命名的目录，免得 kunpeng 病毒库混进宿主列表。"""
    bases = [(DIRS['host_src'], False)]
    bases += [(r, True) for r in _EXTRA_DB_ROOTS]
    out, seen = [], set()
    for base, external in bases:
        if not base or not os.path.isdir(base):
            continue
        key = os.path.normcase(os.path.abspath(base))
        if key in seen:                       # 外部根与 host-db 同目录时去重
            continue
        seen.add(key)
        for n in sorted(os.listdir(base)):
            p = os.path.join(base, n)
            if n == 'host':            # 旧槽位单独处理，不混进列表
                continue
            if external and not n.endswith(_HOST_DB_SUFFIX):
                continue
            if os.path.isdir(p) and glob.glob(os.path.join(p, 'hash_*.k2d')):
                out.append(p)
    return sorted(out)


def detect_active_host_db():
    """自动挑一个：唯一候选 → 它；多个 → 最近构建；都没有 → 旧槽位（存在则用）。"""
    dirs = host_db_dirs()
    if len(dirs) == 1:
        return dirs[0]
    if dirs:
        def _rank(p):
            m = host_db_manifest(p)
            return str(m.get('built_at') or ''), os.path.getmtime(p)
        return max(dirs, key=_rank)
    legacy = legacy_host_db_dir()
    return legacy if os.path.isdir(legacy) else ''


def _host_db_selector(db_dir):
    """库的"选择器"：相对 host-db 的路径（写进 GUI/API，用于切库）。

    为什么不用 basename：旧槽位是 host-db/host/classify，basename 'classify'
    拼回去并不存在；相对路径 'host/classify' 才是无歧义的。
    """
    try:
        rel = os.path.relpath(db_dir, DIRS['host_src'])
    except ValueError:
        return db_dir
    return db_dir if rel.startswith('..') else rel.replace('\\', '/')


def apply_active_host_db(value):
    """设置当前宿主库；返回解析后的绝对路径（无可用库时返回 ''）。

    value: '' = 自动探测；相对路径（如 112863_Lycium_barbarum_host_db 或
    host/classify）= host-db 下同名路径；绝对路径 = 原样采用（须位于平台目录、
    已配置的数据库根或已登记的外部库根内——可用性由
    check_path(in_platform=True) 判定）。相对路径不存在时按**目录名**兜底
    匹配已知库（方便只记得库名的情况）。

    同时记录"是不是显式指定"（_ACTIVE_HOST_DB_EXPLICIT）：显式指定要压过
    platform.json 里旧的 databases.host，而自动探测不应压过它（否则用户显式
    配置的外部宿主库会被 host-db 下任何一个 *_host_db 抢走）。
    """
    global _ACTIVE_HOST_DB, _ACTIVE_HOST_DB_EXPLICIT
    v = (value or '').strip()
    if not v:
        _ACTIVE_HOST_DB = detect_active_host_db()
        _ACTIVE_HOST_DB_EXPLICIT = False
        return _ACTIVE_HOST_DB
    p = os.path.normpath(os.path.abspath(
        v if os.path.isabs(v) else os.path.join(DIRS['host_src'], v)))
    if not os.path.isdir(p):
        cands = [q for q in host_db_dirs() + [legacy_host_db_dir()]
                 if os.path.isdir(q) and os.path.basename(q) == v]
        if len(cands) == 1:
            p = cands[0]
    _ACTIVE_HOST_DB = p
    _ACTIVE_HOST_DB_EXPLICIT = True
    return p


def active_host_db_is_explicit():
    """当前宿主库是不是显式指定的（而非自动探测出来的）。"""
    return _ACTIVE_HOST_DB_EXPLICIT


_ACTIVE_LAZY = False       # 防重入：见 get_active_host_db()
_ACTIVE_HOST_DB_EXPLICIT = False   # 当前库是显式指定(True)还是自动探测(False)


def get_active_host_db():
    """当前宿主库绝对路径（'' = 无）。

    惰性初始化：_ACTIVE_HOST_DB 只有 Config.__init__ 跑过才填充，而 CLI/脚本
    片段常常直接调本函数（不先 get_config()），此时会拿到空串——实测踩过
    （清理中间产物时把空路径当成了库目录）。这里补一次配置加载。
    """
    global _ACTIVE_LAZY
    if not _ACTIVE_HOST_DB and not _ACTIVE_LAZY:
        _ACTIVE_LAZY = True
        try:
            get_config()
        except Exception:
            pass
        finally:
            _ACTIVE_LAZY = False
    return _ACTIVE_HOST_DB


def host_db_info(db_dir):
    """宿主库对外描述：GUI 列表、自检、阶段日志共用同一口径。"""
    if not db_dir:
        return {}
    db_dir = os.path.normpath(os.path.abspath(db_dir))
    m = host_db_manifest(db_dir)
    taxids = host_db_taxids(db_dir)
    taxid, species = host_db_species(db_dir)
    ready = bool(glob.glob(os.path.join(db_dir, 'hash_*.k2d'))) and all(
        os.path.isfile(os.path.join(db_dir, f))
        for f in ('opts.k2d', 'taxo.k2d', 'hash_config.k2d'))
    return {
        'name': os.path.basename(db_dir),
        'path': db_dir,
        'selector': _host_db_selector(db_dir),
        'ready': ready,
        'active': os.path.normcase(db_dir) == os.path.normcase(
            get_active_host_db() or ''),
        'legacy': os.path.normcase(db_dir) == os.path.normcase(legacy_host_db_dir()),
        'taxid': taxid,
        'species': species,
        'taxids_in_map': sorted(taxids),
        'conflicted': len(taxids) > 1,       # 元数据被污染：map 里不止一个 taxid
        'has_manifest': bool(m),
        'built_at': str(m.get('built_at') or ''),
        'source_genome': str(m.get('source_genome') or ''),
        'n_seq': m.get('n_seq'),
        'n_frag': m.get('n_frag'),
    }


CONFIG_FILE = os.path.join(PLATFORM_ROOT, 'platform.json')
# 配置文件写入互斥（Flask threaded=True，设置页可并发保存）
_SAVE_LOCK = threading.RLock()


def _rel_if_under(path):
    """写盘前把平台根之下的路径转成相对路径，保持 platform.json 跨机/跨盘可移植。

    平台外的路径（外置库等）与特殊值（如 'bundled:gbdraw'）原样保留；
    相对值在 _load 侧按 PLATFORM_ROOT 解析回来。
    """
    path = str(path)
    if not os.path.isabs(path):
        return path
    try:
        rel = os.path.relpath(path, PLATFORM_ROOT)
    except ValueError:               # Windows 跨盘符无法 relpath
        return path
    if rel.startswith('..'):
        return path
    return rel


# ------------------------------------------------------------------
# 自定义输出根目录：设置后所有产物（样品结果 / 工具运行 / 数据下载 /
# 公共数据检索 / LOGAN / 提交准备 / 日志）按推荐结构写入该目录下的
# 同名子目录（自动创建）；留空 = 平台目录内。
# ------------------------------------------------------------------
OUTPUT_SUBDIRS = ('results', 'tool_runs', 'downloads', 'meta_search',
                  'logan', 'submissions', 'logs', 'gb_collections',
                  'ncbi_refs', 'ai_audit')

_OUTPUT_ROOT = ''          # 模块级当前自定义输出根（'' = 默认平台目录）


def apply_output_root(root):
    """设置/切换输出根：原地重写 DIRS 输出类条目并自动创建推荐子目录。

    所有功能在运行时读取 DIRS，因此切换立即对全部模块生效，无需重启。
    """
    global _OUTPUT_ROOT
    root = (root or '').strip()
    if root:
        if not os.path.isabs(root):
            raise ValueError('输出目录需为绝对路径')
        root = os.path.normpath(os.path.abspath(root))
        os.makedirs(root, exist_ok=True)
    _OUTPUT_ROOT = root
    for name in OUTPUT_SUBDIRS:
        # 默认（root 为空）落 PLATFORM_ROOT/run/<name> —— 与 DIRS 初始化口径
        # 一致。原先写 os.path.join(PLATFORM_ROOT, name)，一旦用户把
        # 「输出目录」清空（set_output_root('')），就会把 results/tool_runs…
        # 悄悄退回平台根下的旧位置，绕开 run/ 整理。
        target = os.path.join(root, name) if root \
            else os.path.join(PLATFORM_ROOT, 'run', name)
        DIRS[name] = target
        os.makedirs(target, exist_ok=True)
    # 同步预建 tool_runs 固定子目录（_archive/_tmp/_scripts/_reports）；
    # 懒导入避免与 tool_runs_admin 循环依赖
    try:
        from .tool_runs_admin import ensure_fixed_dirs
        ensure_fixed_dirs(DIRS['tool_runs'])
    except Exception:
        pass
    return _OUTPUT_ROOT


def get_output_root():
    return _OUTPUT_ROOT


# ------------------------------------------------------------------
# 示例数据根（程序 / 数据库 / 示例 三分离打包时，示例可放在程序目录之外）
# ------------------------------------------------------------------
_EXAMPLES_ROOT = ''


def _default_examples_root():
    return os.path.join(PLATFORM_ROOT, 'examples')


def _detect_examples_root():
    """自动探测外置示例目录（无需用户手工配置）。

    探测顺序：
      1) <平台根>/examples（2026-09-10 起的内置布局，与数据库分离）
      2) <平台根>/databases/examples（旧内置布局，兼容未迁移的部署）
      3) <平台根>/../VirusPlatform-Examples/examples
      4) <平台根>/../VirusPlatform-Examples/databases/examples
    找到含 example_viral_contigs.fasta 的目录即采用。
    """
    cands = [
        _default_examples_root(),
        os.path.join(PLATFORM_ROOT, 'databases', 'examples'),
        os.path.join(os.path.dirname(PLATFORM_ROOT), 'VirusPlatform-Examples',
                     'examples'),
        os.path.join(os.path.dirname(PLATFORM_ROOT), 'VirusPlatform-Examples',
                     'databases', 'examples'),
    ]
    for c in cands:
        if os.path.isfile(os.path.join(c, 'example_viral_contigs.fasta')):
            return c
    return ''


def _detect_database_root():
    """自动探测与程序同级的 VirusPlatform-Database 数据库包（三分离分发）。

    程序/数据库分离打包后，两个目录同级放置即零配置可用。只在平台内
    databases/tax_db 尚无 nodes.dmp（新部署/未对接）时才探测，且探测结果
    不写回 platform.json——配置文件保持零硬编码路径、跨机可移植。"""
    cand = os.path.join(os.path.dirname(PLATFORM_ROOT),
                        'VirusPlatform-Database')
    if os.path.isfile(os.path.join(cand, 'databases', 'tax_db', 'nodes.dmp')):
        return cand
    return ''


def apply_examples_root(root):
    """设置/切换示例数据根（'' = 自动探测，找不到则用平台内默认路径）。"""
    global _EXAMPLES_ROOT
    root = (root or '').strip()
    if root:
        if not os.path.isabs(root):
            raise ValueError('示例目录需为绝对路径')
        root = os.path.normpath(os.path.abspath(root))
    else:
        root = _detect_examples_root()
    _EXAMPLES_ROOT = root
    DIRS['examples'] = root or _default_examples_root()
    return _EXAMPLES_ROOT


def get_examples_root():
    return _EXAMPLES_ROOT or DIRS['examples']


def write_roots():
    """允许写入的根目录集合：平台根 + 自定义输出/输入/数据库根 + 已登记的
    外部库根（若有）。check_path(in_platform=True) 据此放行写入类路径。"""
    roots = [PLATFORM_ROOT]
    for r in (*_EXTRA_DB_ROOTS, _OUTPUT_ROOT, _INPUT_ROOT, _DATABASE_ROOT):
        if r:
            roots.append(r)
    return roots


# 输入根：测序数据 fastq/ 与拖拽上传 uploads/ 的存放位置
INPUT_SUBDIRS = ('fastq', 'uploads')

# 数据库根：kunpeng 库/taxonomy/palmdb 等大数据目录 + 建库源数据
_DATABASE_KEYS = {                       # DIRS 键 → 数据库根下的子路径
    'databases': 'databases',
    # taxonomy 不在此表内：它必须走 DB_LAYOUT。写死 'databases/taxonomy' 会让
    # 外部数据库根下的 taxonomy 指向空目录（见 apply_database_root）。
    'host_src':  'host-db',
    'virus_src': os.path.join('databases', 'virusref_db'),
}

_INPUT_ROOT = ''           # 模块级当前自定义输入根（'' = 默认平台目录）
_DATABASE_ROOT = ''        # 模块级当前自定义数据库根（'' = 默认平台目录）
_EXTRA_DB_ROOTS = []       # 已登记的外部库根（自备宿主库等，可位于平台外）


def apply_extra_db_roots(roots):
    """用给定列表重置模块级外部库根（规范化，只保留真实存在的目录）。"""
    global _EXTRA_DB_ROOTS
    out = []
    for r in (roots or []):
        try:
            p = os.path.normpath(os.path.abspath(str(r).strip()))
        except (TypeError, ValueError, OSError):
            continue
        if p not in out and os.path.isdir(p):
            out.append(p)
    _EXTRA_DB_ROOTS = out
    return list(_EXTRA_DB_ROOTS)


def register_extra_db_root(path):
    """登记一个外部库根目录（含整棵子树）进入路径白名单并持久化。

    为什么需要：kunpeng 库的可用性检查（db_ready）内部走
    check_path(in_platform=True)，其放行范围就是 write_roots()。自备宿主库
    放在平台/数据库根之外时，若不先登记根目录，加载会在后续所有就绪检查
    里被静默拒绝。与「设置 → 数据库目录」同一信任级别：用户在 UI 里显式
    选定的目录即可信。重复登记幂等；目录消失时 _load 会自动剔除。"""
    p = os.path.normpath(os.path.abspath(str(path).strip()))
    if not os.path.isdir(p):
        raise ValueError(f'目录不存在: {path}')
    cfg = get_config()
    if p not in cfg.extra_db_roots:
        cfg.extra_db_roots.append(p)
        cfg.save()
    # 以配置实例为准重建模块级列表：get_config() 首次构造 Config 时会按
    # platform.json 的旧内容 apply_extra_db_roots，把先改的模块状态抹掉
    # （实测踩过：先 append 再 get_config，登记在进程内不生效）。
    apply_extra_db_roots(cfg.extra_db_roots)
    return p



def apply_input_root(root):
    """设置/切换输入根：重写 DIRS['fastq'/'uploads'] 并自动创建。"""
    global _INPUT_ROOT
    root = (root or '').strip()
    if root:
        if not os.path.isabs(root):
            raise ValueError('输入目录需为绝对路径')
        root = os.path.normpath(os.path.abspath(root))
        os.makedirs(root, exist_ok=True)
    _INPUT_ROOT = root
    for name in INPUT_SUBDIRS:
        DIRS[name] = os.path.join(root, name) if root \
            else os.path.join(PLATFORM_ROOT, name)
        os.makedirs(DIRS[name], exist_ok=True)
    return _INPUT_ROOT


def apply_database_root(root):
    """设置/切换数据库根：重写 DIRS 数据库类条目并自动创建。"""
    global _DATABASE_ROOT
    root = (root or '').strip()
    if root:
        if not os.path.isabs(root):
            raise ValueError('数据库目录需为绝对路径')
        root = os.path.normpath(os.path.abspath(root))
        os.makedirs(root, exist_ok=True)
    _DATABASE_ROOT = root
    for key, sub in _DATABASE_KEYS.items():
        DIRS[key] = os.path.join(root, sub) if root \
            else os.path.join(PLATFORM_ROOT, sub)
        os.makedirs(DIRS[key], exist_ok=True)
    # taxonomy 走注册表解析：新扁平 <root>/databases/tax_db 优先，
    # 旧布局 <root>/databases/tax/core 兜底（与平台内 DIRS 初始化口径一致）。
    DIRS['taxonomy'] = db_path('tax', 'core')
    os.makedirs(DIRS['taxonomy'], exist_ok=True)
    # virus_src 由上面的 _DATABASE_KEYS 循环统一落位（<root>/databases/virusref_db）
    _refresh_module_paths()
    return _DATABASE_ROOT


def _refresh_module_paths():
    """通知「import 时快照了 DIRS」的模块重算路径常量。

    否则切换数据库目录后，ictv_db / virus_ref / universal_ref /
    local_search / verify / hmm_annot / orf_annot / suvtk_submit
    仍指向旧库（长驻的 Web 进程尤其明显）。这类失配不会报错，
    只会让功能静默失效（库找不到 → 该层结果为空）。
    用延迟导入避免 config ↔ 这些模块的循环依赖。
    """
    for name in ('ictv_db', 'virus_ref', 'universal_ref', 'local_search',
                 'verify', 'hmm_annot', 'orf_annot', 'suvtk_submit',
                 'bfvd_search'):
        try:
            mod = __import__(f'Virus_Platform_Core.{name}', fromlist=['refresh_paths'])
            fn = getattr(mod, 'refresh_paths', None)
            if callable(fn):
                fn()
        except Exception as e:
            # 失配的后果是「该层结果为空」且无报错（见 docstring），至少留痕
            print(f'[config] 模块路径刷新失败 {name}: {e}', file=sys.stderr)


def get_input_root():
    return _INPUT_ROOT


def get_database_root():
    return _DATABASE_ROOT


def _ensure_not_system_dir(path):
    """拒绝把系统关键目录设为数据根。

    write_roots() 会把自定义输出/输入/数据库根视为合法写入区
    （check_path in_platform=True 放行其下的写入与删除），一旦被设成
    C:\\Windows 之类，删除类端点就能波及系统文件。"""
    drv = os.environ.get('SYSTEMDRIVE', 'C:').rstrip('\\') or 'C:'
    norm = os.path.normcase(os.path.normpath(os.path.abspath(path)))
    for seg in ('Windows', 'Program Files', 'Program Files (x86)',
                'ProgramData'):
        # 不能用 os.path.join(drv, seg)：drv='C:' 时 join 产出 'C:Windows'
        # （盘符相对路径），normpath 后与 'c:\\windows' 永不相等 → 整个
        # 守卫恒放行（实测过）。必须显式补上盘符后的分隔符。
        base = os.path.normcase(os.path.normpath(drv + os.sep + seg))
        if norm == base or norm.startswith(base + os.sep):
            raise ValueError(f'不允许把系统目录设为数据根: {path}')


def _safe_child(path, parent=PLATFORM_ROOT):
    """规范化路径并校验必须位于 parent 目录内，禁止 ../ 逃逸。"""
    norm = os.path.normpath(os.path.abspath(path))
    parent = os.path.normpath(os.path.abspath(parent))
    if norm == parent or norm.startswith(parent + os.sep):
        return norm
    raise ValueError(f"路径越界（仅允许平台目录内）: {path}")


def _first_exist(paths):
    for p in paths:
        if p and os.path.isfile(p):
            return os.path.abspath(p)
    return None


def _glob_first(pattern):
    hits = sorted(glob.glob(pattern))
    return os.path.abspath(hits[0]) if hits else None


def detect_tools():
    """自动探测外部工具路径（仅 PATH 与 3rd/bin、3rd/tools 目录，不访问其他用户目录）。
    3rd/bin/   = 单文件可执行（kunpeng、seqkit、crabz 等）
    3rd/tools/ = 带目录结构的工具套件（Blast、mafft-win、iQtree 等）
    下面各条目里沿用 PLATFORM_ROOT 的候选是「旧版平铺分发」回退（exe 与 app.py
    同级），仍然保留以兼容已分发的 exe 平台。"""
    R = PLATFORM_ROOT
    BIN = os.path.join(R, '3rd', 'bin')
    T = os.path.join(R, '3rd', 'tools')
    which = shutil.which
    cands = {
        'kunpeng': [
            _glob_first(os.path.join(BIN, 'kun_peng-v*-x86_64-pc-windows-msvc.exe')),
            _glob_first(os.path.join(BIN, 'kun_peng*.exe')),
            _glob_first(os.path.join(R, 'kun_peng-v*-x86_64-pc-windows-msvc.exe')),
            _glob_first(os.path.join(R, 'kun_peng*.exe')),
            which('kun_peng'), which('kun_peng.exe'),
        ],
        # 捆绑副本优先于 PATH：分发版始终用随包 SPAdes，不受对方机器已装版本影响
        'spades': [_glob_first(os.path.join(T, 'SPAdes*', 'bin', 'spades.bat')),
                   _glob_first(os.path.join(R, 'SPAdes*', 'bin', 'spades.bat')),
                   which('spades.bat'), which('spades.py')],
        'blastn':      [_glob_first(os.path.join(T, 'Blast', 'bin', 'blastn.exe')),
                        _glob_first(os.path.join(R, 'Blast', 'bin', 'blastn.exe')),      which('blastn')],
        'makeblastdb': [_glob_first(os.path.join(T, 'Blast', 'bin', 'makeblastdb.exe')),
                        _glob_first(os.path.join(R, 'Blast', 'bin', 'makeblastdb.exe')), which('makeblastdb')],
        'blastp':      [_glob_first(os.path.join(T, 'Blast', 'bin', 'blastp.exe')),
                        _glob_first(os.path.join(R, 'Blast', 'bin', 'blastp.exe')),      which('blastp')],
        'blastx':      [_glob_first(os.path.join(T, 'Blast', 'bin', 'blastx.exe')),
                        _glob_first(os.path.join(R, 'Blast', 'bin', 'blastx.exe')),      which('blastx')],
        # mafft/fasttree/trimal/seqkit 统一「捆绑优先、PATH 兜底」——与
        # spades/blast 同款口径：随包副本存在就用它，版本不随对方机器已装
        # 的版本漂移；_first_exist 逐项做 os.path.isfile 判断，无 3rd 的
        # 源码环境自然回退 which()。
        'mafft':   [os.path.join(T, 'mafft-win', 'mafft.bat'),
                    _glob_first(os.path.join(T, 'mafft-win', 'mafft*.bat')),
                    os.path.join(R, 'mafft-win', 'mafft.bat'),
                    _glob_first(os.path.join(R, 'mafft-win', 'mafft*.bat')),
                    which('mafft')],
        'fasttree': [_glob_first(os.path.join(T, 'FastTree', 'FastTree*.exe')),
                     _glob_first(os.path.join(R, 'FastTree', 'FastTree*.exe')),
                     which('fasttree'), which('FastTree')],
        # RAxML-NG：exe 必须与 3 个 msys DLL 同目录（捆绑副本按此布局打包）
        'raxml-ng': [which('raxml-ng'), which('raxml-ng.exe'),
                     _glob_first(os.path.join(T, 'raxml-ng*', 'raxml-ng*.exe')),
                     _glob_first(os.path.join(R, 'raxml-ng*', 'raxml-ng*.exe'))],
        # IQ-TREE 只认 v3（v2 已退役，不再探测，避免自检长期误报"缺 iqtree2"）
        # —— 建树引擎已换 RAxML-NG（实测 ML 搜索 143 s vs IQ-TREE 1049 s，快 ~7 倍；
        #    含 100 次 FBP 自举总耗时 730.6 s vs IQ-TREE 1000 次 UFBoot 1049 s，
        #    端到端快 ~1.4 倍），此处仅保留工具在位探测
        'iqtree3': [which('iqtree3'), which('iqtree3.exe'),
                    _glob_first(os.path.join(T, 'iQtree', '*', 'bin', 'iqtree3.exe')),
                    _glob_first(os.path.join(T, 'iQtree', '*', 'iqtree3*.exe')),
                    _glob_first(os.path.join(R, 'iQtree', '*', 'bin', 'iqtree3.exe')),
                    _glob_first(os.path.join(R, 'iQtree', '*', 'iqtree3*.exe'))],
        'trimal':  [os.path.join(T, 'trimAl', 'trimal.exe'),
                    _glob_first(os.path.join(T, 'trimAl*', 'bin', 'trimal*.exe')),
                    os.path.join(R, 'trimAl', 'trimal.exe'),
                    _glob_first(os.path.join(R, 'trimAl*', 'bin', 'trimal*.exe')),
                    which('trimal'), which('trimal.exe')],
        # clustalw / gblocks / muscle 不再探测（2026-09-19）：全平台无任何模块
        # 调用它们（比对走 mafft、修剪走 trimAl；RDP5 用的是 3rd/tools/rdp5
        # 自带的 clustalw2，不经过本表）。列在这里只会让自检长期红着
        # 「外部工具 31/34」误导用户，与上方退役 iqtree2 同理。
        'seqkit':   [_glob_first(os.path.join(BIN, 'seqkit*.exe')),
                     _glob_first(os.path.join(R, 'seqkit*.exe')),
                     which('seqkit')],
        'crabz':    [which('crabz'), which('crabz.exe'),
                     _glob_first(os.path.join(BIN, 'crabz*.exe')),
                     _glob_first(os.path.join(R, 'crabz*.exe'))],
        'gbdraw':   [which('gbdraw'), which('gbdraw.exe'),
                     os.path.join(os.path.dirname(sys.executable), 'Scripts',
                                  'gbdraw.exe'),
                     _glob_first(os.path.join(BIN, 'gbdraw*', 'gbdraw*.exe')),
                     _glob_first(os.path.join(R, 'gbdraw*', 'gbdraw*.exe'))],
        'orfipy':   [which('orfipy'),
                     # 打包运行时 Scripts 不一定在 PATH：按当前 Python 推断
                     os.path.join(os.path.dirname(sys.executable), 'Scripts', 'orfipy.exe')],
        'diamond':  [os.path.join(T, 'diamond', 'diamond.exe'),
                     os.path.join(R, 'diamond', 'diamond.exe'), which('diamond'),
                     which('diamond.exe')],
        'mmseqs':   [os.path.join(T, 'mmseqs', 'bin', 'mmseqs.exe'),
                     os.path.join(R, 'mmseqs', 'bin', 'mmseqs.exe'),
                     which('mmseqs'), which('mmseqs.exe')],
        'minimap2': [os.path.join(T, 'minimap2', 'minimap2.exe'),
                     os.path.join(T, 'minimap2', 'minimap2'),
                     os.path.join(R, 'minimap2', 'minimap2.exe'),
                     os.path.join(R, 'minimap2', 'minimap2'),
                     which('minimap2'), which('minimap2.exe')],
        # minibwa：minimap2 作者的下一代比对器，与 minimap2 同属单文件工具，
        # 故探测 bin/；随包分发不依赖 PATH。
        'minibwa': [_glob_first(os.path.join(BIN, 'minibwa*.exe')),
                     os.path.join(BIN, 'minibwa.exe'),
                     _glob_first(os.path.join(R, 'minibwa*.exe')),
                     which('minibwa'), which('minibwa.exe')],
        'fastp':    [which('fastp'), which('fastp.exe'),
                     _glob_first(os.path.join(T, 'fastp*', 'fastp.exe')),
                     _glob_first(os.path.join(R, 'fastp*', 'fastp.exe')),
                     _glob_first(os.path.join(R, 'biosoft', 'fastp*', 'fastp.exe')),
                     _glob_first(os.path.join(R, 'fastp*.exe'))],
        'aria2c':   [os.path.join(BIN, 'aria2c.exe'), os.path.join(R, 'aria2c.exe'),
                     which('aria2c'), which('aria2c.exe')],
        'sracha':   [os.path.join(BIN, 'sracha.exe'), os.path.join(R, 'sracha.exe'),
                     which('sracha'), which('sracha.exe')],
        # ---- 以下工具原先不在探测表里，只靠 platform.json 手填绝对路径；
        # 打包版用的是干净配置（tools: {}），于是"已随包"的它们反而探测不到
        # （实测发布版少了 salmon/samtools/table2asn 三个）。这里补齐候选，
        # 让源码模式与发布版行为一致、且不依赖手改配置。
        'samtools': [which('samtools'), which('samtools.exe'),
                     os.path.join(T, 'samtools', 'bin', 'samtools.exe'),
                     os.path.join(BIN, 'samtools.exe'),
                     os.path.join(R, 'samtools.exe')],
        'bcftools': [which('bcftools'), which('bcftools.exe'),
                     os.path.join(T, 'bcftools', 'bin', 'bcftools.exe'),
                     os.path.join(R, 'bcftools.exe')],
        'salmon':   [which('salmon'), which('salmon.exe'),
                     os.path.join(T, 'salmon2', 'salmon.exe'),
                     os.path.join(T, 'salmon', 'bin', 'salmon.exe'),
                     os.path.join(R, 'salmon.exe')],
        'table2asn': [which('table2asn'), which('table2asn.exe'),
                      os.path.join(T, 'table2asn', 'table2asn.exe'),
                      os.path.join(R, 'table2asn.exe')],
        'pandepth': [which('pandepth'), which('pandepth.exe'),
                     _glob_first(os.path.join(T, 'pandepth*', 'pandepth.exe')),
                     os.path.join(R, 'pandepth.exe')],
        'viral_consensus': [
            which('viral_consensus'), which('viral_consensus.exe'),
            _glob_first(os.path.join(T, 'viral_consensus*',
                                     'viral_consensus.exe')),
            os.path.join(R, 'viral_consensus.exe')],
        'lsd2':     [which('lsd2'), which('lsd2.exe'),
                     os.path.join(T, 'lsd2', 'lsd2.exe'),
                     os.path.join(BIN, 'lsd2.exe'),
                     os.path.join(R, 'lsd2.exe')],
        # RDP5CL（重组检测命令行）与 dsRNAmax（dsRNA 判定）原先只靠代码内的
        # 目录扫描兜底：设置页/自检看不到、platform.json 路径覆盖无效。
        # （dsRNAmax 条目 2026-09-29 随 dsRNA 迁出删除，2026-10-03 恢复。）
        'RDP5CL':   [os.path.join(T, 'rdp5', 'RDP5CL.exe'),
                     os.path.join(R, 'RDP5CL.exe'), which('RDP5CL')],
        'dsRNAmax': [os.path.join(T, 'dsrnamax', 'dsRNAmax_det.exe'),
                     os.path.join(T, 'dsrnamax', 'dsRNAmax.exe'),
                     os.path.join(BIN, 'dsrnamax.exe'),
                     os.path.join(R, 'dsRNAmax.exe'), which('dsRNAmax')],
        # 注：snpEff（jar）与 SNPGenie（perl 脚本）不是可直接执行的单文件工具，
        # 由 known_virus_suite 自己按 tools/snpeff/snpEff/snpEff.jar、
        # tools/snpgenie/snpgenie.pl 定位，不放进本探测表（避免"jar 当 exe"）。
    }
    out = {name: _first_exist(lst) for name, lst in cands.items()}
    if getattr(sys, 'frozen', False) and not out.get('gbdraw'):
        out['gbdraw'] = 'bundled:gbdraw'     # 随包 Python 引擎，进程内调用
    return out


def _ai_local_path():
    """AI API Key 的落盘位置：<平台根>/run/ai_local.json（固定，不随输出根迁移）。

    为什么不放 platform.json：那是要进 git 的配置文件，密钥进去就是事故。
    run/ 整体被 .gitignore 覆盖（.gitignore 第 25 行）。
    为什么不随自定义输出根迁移：Key 是机器本地的密钥，不是分析数据——
    配置的 _load 发生在 apply_output_root 重映射 DIRS **之前**，若跟随
    输出根会造成"存 A 读 B"；固定在平台根下则存读永远一致（platform.json
    就写在平台根，可写性有同样的保证）。
    """
    return os.path.join(PLATFORM_ROOT, 'run', 'ai_local.json')


class Config:
    """全局配置单例。tools 可被 platform.json 覆盖；未探测到的工具调用时报错。"""

    # 分析默认参数的定义（键 → 类型/默认/中英说明）；设置页据此渲染
    DEFAULT_FIELDS = [
        # key, 类型, 默认, 中文说明, 英文说明[, 选项]
        ('confidence',    'float',  0.0,        '分类置信度（0~1，越大越严格）',
         'Classification confidence (0-1)'),
        ('subsample',     'int',    0,          '默认子采样 reads 对数（0=不截取）',
         'Default subsample pairs (0 = all)'),
        ('fastp_dedup',   'bool',   False,      'Fastp 质控默认去重 (--dedup)',
         'Dedup by default in fastp QC'),
        ('fastp_min_len', 'int',    15,         'Fastp 最短保留读长 bp（sRNA 测序请保持 ≤18，原 50 会丢光小 RNA）',
         'Fastp min read length to keep (keep <=18 for sRNA)'),
        ('do_fq2fa',      'bool',   True,       '默认 FASTQ→FASTA 预转换（分类提速）',
         'FASTQ to FASTA pre-conversion by default'),
        ('assembly_mode', 'select', 'rnaviral', 'SPAdes 组装模式',
         'SPAdes assembly mode',
         ['rnaviral', 'metaviral', 'rna', 'meta', 'isolate', 'srna']),
        ('assembly_input','select', 'virus',    '组装输入 reads',
         'Assembly input reads',
         ['virus', 'kept', 'raw']),
        ('memory',        'int',    64,         'SPAdes 内存上限 (GB)',
         'SPAdes memory cap (GB)'),
        ('min_contig_len','int',    200,        '最小 contig 长度 (bp)',
         'Min contig length (bp)'),
        ('min_orf_aa',    'int',    100,        '最小 ORF 长度 (aa)',
         'Min ORF length (aa)'),
        ('top_n_refs',    'int',    10,         '进化分析近缘参考数',
         'Closest references for phylogenetics'),
        ('tree_tool',     'select', 'fasttree', '建树工具',
         'Tree building tool', ['fasttree', 'nj', 'raxml-ng']),
        ('tree_sampling', 'select', 'blast',    '建树参考抽样策略',
         'Tree reference sampling', ['blast', 'macro', 'genus', 'lineage']),
        ('primer_mode',   'select', 'conserved', '引物设计模式',
         'Primer design mode', ['conserved', 'plain']),
        ('specificity',   'bool',   False,      '引物宿主特异性检查（首次较慢）',
         'Primer host-specificity check'),
        ('gbdraw_max',    'int',    12,         '基因组图出图上限（条）',
         'Genome plot limit (contigs)'),
        ('plot_engine',   'select', 'auto',     '基因组图引擎',
         'Genome plot engine', ['auto', 'gbdraw', 'dfv']),
        ('force',         'bool',   False,      '默认强制重跑（忽略断点）',
         'Force re-run by default'),
        ('do_trim',       'bool',   True,       '⑦ 建树前用 trimAl 清剪比对',
         'Trim alignment with trimAl before tree building'),
    ]

    def __init__(self):
        self.threads = max(1, (os.cpu_count() or 4) - 1)
        self.tools = detect_tools()
        # 插件层（plugins/ 目录，机制见 plugin_manager.py）：独立于 tools/
        # databases 存放，tool()/tool_status() 只把它当「填缺」备份，save()
        # 不落盘——卸载插件后配置零残留。_platform_tool_keys 记录 platform.json
        # 显式指定的工具键（插件不得压过用户手工配置）。
        self.plugin_tools = {}
        self.plugin_databases = {}
        self._platform_tool_keys = set()
        # 插件源索引 URL（'' = 用 plugin_manager.DEFAULT_INDEX_URL，即本仓库
        # plugins-repo/index.json 的 raw 地址）；国内网络可换成镜像地址。
        self.plugin_index_url = ''
        self.email = ''            # NCBI EUtils 联系邮箱（可选，礼貌性参数）
        self.language = 'zh'       # 界面语言: zh | en
        # AI 助手（对话驱动智能代理）：OpenAI 兼容端点配置。键结构与预设
        # 清单见 ai_agent/provider.py（PROVIDER_PRESETS）；api_key 只存本地
        # run/ai_local.json（密钥分离，见 save()），platform.json 永不落 Key。
        # explain/broadcast/gate = 主动驱动三层开关（讲解/解读/决策闸门），
        # 详见 ai_agent/broadcast.py；gate 默认关——它会改变流水线行为
        # （关键阶段收尾后阻塞等确认），必须用户显式开启。
        self.ai = {'enabled': False, 'provider': 'deepseek', 'base_url': '',
                   'api_key': '', 'model': '', 'temperature': 0.3,
                   'max_tokens': 2048,
                   'explain': True, 'broadcast': True, 'gate': False,
                   'gate_timeout_min': 15,
                   'broadcast_stages': ['kvsuite', 'assembly', 'report'],
                   'gate_stages': ['kvsuite', 'assembly'],
                   # 执行权分级（见 ai_agent/executors.py）：
                   # exec_auto     = 第二级：阶段失败后原参数自动重试（限1次）
                   # exec_auto_tune= 第三级：AI 提参自动重跑闭环（限1次，
                   #                 调参失败回落第二级）。两者默认关。
                   'exec_auto': False, 'exec_auto_tune': False,
                   # 供应商档案（ZCode 式多端点并存）：每套 name/base_url/model，
                   # Key 永不进本文件（分离在 run/ai_local.json provider_keys）
                   'providers': [], 'active_profile': ''}
        # 分析默认参数（预填硬编码缺省，platform.json 的已存值在 _load_extras 覆盖）
        self.defaults = {k: v for k, _t, v, *_r in self.DEFAULT_FIELDS}
        self.databases = {
            'host':  db_path('host', 'classify'),
            'virus': db_path('virus', 'plant'),
        }
        self.output_root = ''      # 自定义输出根（'' = 平台目录内）
        self.input_root = ''       # 自定义输入根（fastq / uploads）
        self.database_root = ''    # 自定义数据库根（databases 含 virus_src/ + host-db）
        self.examples_root = ''    # 自定义示例数据根（'' = 自动探测）
        # 当前宿主分类库：留空 = 自动探测（唯一/最近构建的 *_host_db，再退旧槽位）
        self.active_host_db = ''
        self._db_override = {}
        # 使用者登记的外部已构建 kunpeng 病毒库（绝对路径，平台内/数据库根优先）
        self.extra_virus_libs = []
        # 使用者登记的外部库根（自备宿主库等可位于平台外，进入写白名单）
        self.extra_db_roots = []
        # 使用者登记的外部病毒鉴定库（salmon 比对索引目录，可为平台外）
        self.extra_kv_indexes = []
        self._load()
        for attr, apply_fn in (('input_root', apply_input_root),
                               ('database_root', apply_database_root),
                               ('output_root', apply_output_root)):
            path = getattr(self, attr)
            if path:
                try:
                    setattr(self, attr, apply_fn(path))
                except (OSError, ValueError):
                    setattr(self, attr, '')
        # 示例数据根：显式配置优先，否则自动探测（程序/数据库/示例三分离
        # 打包时示例目录位于程序目录之外，靠这里接上）。
        apply_examples_root(self.examples_root)
        # 外部库根先于宿主库解析登记：active_host_db 的可用性校验
        # （db_ready → check_path）依赖 write_roots() 已包含这些根。
        apply_extra_db_roots(self.extra_db_roots)
        # 与程序同级的数据库包自动接上（零配置分发）：仅当用户未显式配置
        # 数据库根、且平台内 databases/tax_db 还没有 nodes.dmp（新部署）
        # 时才探测；探测结果不落盘，用户随时可在设置页显式覆盖。
        if not self.database_root and not os.path.isfile(
                os.path.join(DIRS['databases'], 'tax_db', 'nodes.dmp')):
            _detected_root = _detect_database_root()
            if _detected_root:
                apply_database_root(_detected_root)
        # 数据库实际路径按（可能已迁移的）数据库根重算；
        # platform.json 的 databases 显式覆盖：仅当该路径真实存在时生效
        # （迁移后旧路径不存在则回退注册表解析，自动切到新分类路径）。
        # 宿主库要先应用 active_host_db，db_path('host','classify') 才算得对。
        apply_active_host_db(self.active_host_db)
        # 配置指向的库可能已被删除/换盘：不要死守一个坏路径，回落到自动探测
        # （否则 db_path 返回不存在的目录，selfcheck 只报 ✘ 看不出原因）。
        if self.active_host_db and not os.path.isdir(get_active_host_db()):
            apply_active_host_db('')
        self.databases = {
            'host':  db_path('host', 'classify'),
            'virus': db_path('virus', 'plant'),
        }
        for k, v in (getattr(self, '_db_override', {}) or {}).items():
            if not v:
                continue
            pv = os.path.normpath(os.path.abspath(v))
            if not os.path.isdir(pv):
                continue
            if not any(pv == os.path.normpath(r)
                       or pv.startswith(os.path.normpath(r) + os.sep)
                       for r in write_roots()):
                continue
            if k == 'host' and active_host_db_is_explicit():
                # 宿主库以**显式**的 active_host_db 为准：platform.json 的
                # databases.host 会被旧版本的常驻实例按自己的 schema 写回
                # （它不认识 active_host_db）。实测该键被写回已删除的
                # host-db/host/classify，若照它走，平台会静默换用另一个
                # （可能被污染的）宿主库。
                # 注意只压过"显式"的 active_host_db：自动探测出来的库不应
                # 抢走用户在 databases.host 里显式配置的宿主库。
                continue
            self.databases[k] = pv
        # 插件最后挂载：platform.json 覆盖与注册表解析已就位，插件只填缺
        # （tools 的填缺在 tool()/tool_status() 里做，见 plugin_manager.py）。
        try:
            self.refresh_plugins()
        except Exception as e:            # 坏插件/坏清单不拦启动
            sys.stderr.write(f'[config] 插件扫描失败（已跳过）: {e}\n')
        for d in ('databases', 'taxonomy', 'results', 'downloads', 'logan',
                  'submissions', 'meta_search', 'tasks', 'logs', 'ai_audit'):
            os.makedirs(DIRS[d], exist_ok=True)

    # ---- 持久化 -------------------------------------------------
    def _load(self):
        if not os.path.isfile(CONFIG_FILE):
            return
        try:
            with open(CONFIG_FILE, 'r', encoding='utf-8') as f:
                data = json.load(f)
            self._load_extras(data)
            self.output_root = str(data.get('output_root') or '').strip()
            self.input_root = str(data.get('input_root') or '').strip()
            self.database_root = str(data.get('database_root') or '').strip()
            self.examples_root = str(data.get('examples_root') or '').strip()
            self.active_host_db = str(data.get('active_host_db') or '').strip()
            # threads: 0 / 缺省 = 自动探测（cpu 核数 - 1，见 __init__）；
            # 显式 >0 才采用——干净发布配置里写 0 不应钳成单线程
            try:
                _t = int(data.get('threads') or 0)
            except (TypeError, ValueError):
                _t = 0
            if _t > 0:
                self.threads = _t
            if data.get('email'):
                self.email = str(data['email'])
            if data.get('plugin_index_url') is not None:
                self.plugin_index_url = str(data['plugin_index_url']).strip()
            for k, v in data.get('tools', {}).items():
                if not v:
                    continue
                # 相对路径按平台根解析（platform.json 可存 rel 路径，
                # 避免搬动目录后绝对路径失效）
                p = v if os.path.isabs(str(v)) else os.path.join(PLATFORM_ROOT,
                                                                 str(v))
                if os.path.isfile(p):
                    self.tools[k] = os.path.abspath(p)
                    self._platform_tool_keys.add(k)
            # databases 覆盖延后到根目录应用之后再校验（见 __init__）
            self._db_override = {}
            for k, v in data.get('databases', {}).items():
                v = str(v).strip()
                if not v:
                    continue
                # 相对路径按平台根解析（与 tools 一致，不依赖启动工作目录）
                self._db_override[k] = (v if os.path.isabs(v)
                                        else os.path.join(PLATFORM_ROOT, v))
            libs = []
            for v in (data.get('extra_virus_libs') or []):
                v = str(v).strip()
                if v and not os.path.isabs(v):
                    v = os.path.join(PLATFORM_ROOT, v)
                v = os.path.normpath(os.path.abspath(v))
                if v and os.path.isdir(v):
                    libs.append(v)
            self.extra_virus_libs = libs
            roots = []
            for v in (data.get('extra_db_roots') or []):
                v = str(v).strip()
                if v and not os.path.isabs(v):
                    v = os.path.join(PLATFORM_ROOT, v)
                v = os.path.normpath(os.path.abspath(v))
                if v and os.path.isdir(v):
                    roots.append(v)
            self.extra_db_roots = roots
            kvi = []
            for v in (data.get('extra_kv_indexes') or []):
                v = str(v).strip()
                if v and not os.path.isabs(v):
                    v = os.path.join(PLATFORM_ROOT, v)
                v = os.path.normpath(os.path.abspath(v))
                if v and os.path.isdir(v):
                    kvi.append(v)
            self.extra_kv_indexes = kvi
        except Exception as e:
            # 不要静默吞掉：配置损坏时用户"改了设置不生效"且无从发现。
            sys.stderr.write(f'[Virus_Platform_Core.config] 读取 platform.json 失败，'
                             f'已回退默认配置: {type(e).__name__}: {e}\n')

    def save(self):
        # AI API Key 分离落盘：真身写 run/ai_local.json（run/ 整体 gitignore，
        # 见 .gitignore），platform.json 里**永远写空串**——本文件进 git，
        # 密钥进了就是事故。密钥文件写失败时只报 stderr、平台内存态仍可用，
        # 绝不退回"写进 platform.json"的旧路。
        ai_out = dict(getattr(self, 'ai', {}) or {})
        try:
            lp = _ai_local_path()
            os.makedirs(os.path.dirname(lp), exist_ok=True)
            # 合并写：保留 ai_local.json 里的其他本地字段（如供应商档案
            # 的 provider_keys）——整文件覆盖会把它们抹掉。
            prev = {}
            if os.path.isfile(lp):
                try:
                    with open(lp, encoding='utf-8') as f:
                        prev = json.load(f)
                    if not isinstance(prev, dict):
                        prev = {}
                except (OSError, ValueError):
                    prev = {}
            prev['api_key'] = ai_out.get('api_key', '')
            with open(lp, 'w', encoding='utf-8') as f:
                json.dump(prev, f, ensure_ascii=False, indent=1)
        except OSError as e:
            if ai_out.get('api_key'):
                sys.stderr.write(f'[config] 警告：密钥文件写入失败（{e}），'
                                 f'API Key 本次不落盘，重启后需重新填写\n')
        ai_out['api_key'] = ''
        # 平台根之下的路径一律转相对再落盘：platform.json 搬目录/换电脑/换
        # 盘符仍可用（_load 侧按 PLATFORM_ROOT 解析回来）；平台外的路径
        # （外置库等）保持绝对，属用户显式登记的外部资源。
        data = {'threads': self.threads,
                'plugin_index_url': self.plugin_index_url,
                'tools': {k: _rel_if_under(v) for k, v in self.tools.items()},
                # 插件提供的库不落盘（值与 plugin_databases 相同即为插件来源）：
                # 卸载插件后 platform.json 里不能残留死路径
                'databases': {k: _rel_if_under(v) for k, v in self.databases.items()
                              if self.plugin_databases.get(k) != v},
                'defaults': self.defaults, 'language': self.language,
                'ai': ai_out,
                'output_root': self.output_root,
                'input_root': self.input_root,
                'database_root': self.database_root,
                'examples_root': self.examples_root,
                'active_host_db': self.active_host_db,
                'extra_virus_libs': [_rel_if_under(v)
                                     for v in self.extra_virus_libs],
                'extra_kv_indexes': [_rel_if_under(v)
                                     for v in self.extra_kv_indexes],
                'extra_db_roots': [_rel_if_under(v)
                                   for v in self.extra_db_roots]}
        # 保留本版本不认识的键：平台上"长驻实例"很常见（分析要跑很久），
        # 磁盘上的代码更新后旧进程仍在跑；它一旦保存设置就会按自己的 schema
        # 重写配置，把新版本新增的键**抹掉**。实测踩过：旧实例（不认
        # active_host_db）一次保存就把 active_host_db 抹了，databases.host 也
        # 被写回旧路径。这里先读回磁盘上的未知键再覆盖已知键——旧进程不会
        # 再吃掉新版本的配置项。
        try:
            if os.path.isfile(CONFIG_FILE):
                with open(CONFIG_FILE, 'r', encoding='utf-8') as f:
                    old = json.load(f)
                if isinstance(old, dict):
                    for k, v in old.items():
                        data.setdefault(k, v)
        except (OSError, ValueError):
            pass                     # 旧文件损坏/读不了：按全新写入即可
        if self.email:
            data['email'] = self.email
        # 原子写：先写临时文件再替换，断电/崩溃不会损坏现有配置。
        # 加锁 + 临时名带 pid/线程 id：Flask 以 threaded=True 运行，设置页
        # 并发保存会共用同一个 platform.json.tmp，交错写入会写出损坏 JSON。
        with _SAVE_LOCK:
            tmp = f'{CONFIG_FILE}.{os.getpid()}.{threading.get_ident()}.tmp'
            try:
                with open(tmp, 'w', encoding='utf-8') as f:
                    json.dump(data, f, ensure_ascii=False, indent=2)
                os.replace(tmp, CONFIG_FILE)
            finally:
                if os.path.isfile(tmp):
                    try:
                        os.remove(tmp)
                    except OSError:
                        pass

    # ---- 语言 ---------------------------------------------------
    @property
    def lang(self):
        return 'en' if self.language == 'en' else 'zh'

    def set_input_root(self, path):
        """设置自定义输入根（'' = 恢复默认平台目录），持久化并立即生效。"""
        path = (path or '').strip()
        if path and not os.path.isabs(path):
            raise ValueError('输入目录需为绝对路径')
        if path:
            _ensure_not_system_dir(path)
        self.input_root = apply_input_root(path)
        self.save()

    def set_database_root(self, path):
        """设置自定义数据库根（'' = 恢复默认平台目录），持久化并立即生效；
        数据库实际路径（host_db/virus_db）按新根重算。"""
        path = (path or '').strip()
        if path and not os.path.isabs(path):
            raise ValueError('数据库目录需为绝对路径')
        if path:
            _ensure_not_system_dir(path)
        self.database_root = apply_database_root(path)
        # host_src 随数据库根改变 → 当前宿主库要按新根重解析
        apply_active_host_db(self.active_host_db)
        self.databases = {
            'host':  db_path('host', 'classify'),
            'virus': db_path('virus', 'plant'),
        }
        self.save()

    def set_active_host_db(self, name):
        """切换当前宿主分类库并持久化；'' = 恢复自动探测。

        name 可以是 host-db 下的目录名（推荐，可移植）、绝对路径，或 ''。
        只写"用户显式指定"的值：自动探测出来的本机绝对路径不写回配置，
        否则 platform.json 没法在别的机器/盘符上用。

        校验在**落盘之前**完成并不留痕：早先的实现先应用再让调用方校验，
        结果是"接口返回 400、配置却已经被改成无效路径"（实测踩过——
        name='host' 会解析到 host-db/host 这个父目录而不是库目录）。
        """
        name = (name or '').strip()
        if name:
            resolved = apply_active_host_db(name)
            from .kunpeng import db_ready          # 延迟导入：kunpeng 依赖 config
            problems = []
            if not os.path.isdir(resolved):
                problems.append('目录不存在')
            elif not db_ready(resolved):
                problems.append('不是可用的宿主库（缺 hash_*.k2d / opts.k2d / '
                                'taxo.k2d / hash_config.k2d）')
            if problems:
                apply_active_host_db(self.active_host_db)     # 回滚内存状态
                avail = [os.path.basename(p) for p in host_db_dirs()]
                tip = ('；已有的宿主库: ' + '、'.join(avail)) if avail else \
                      '；还没有任何宿主库，请先构建'
                raise ValueError(f'{name}: {"，".join(problems)}{tip}')
            # host-db 下的库只记目录名（配置可移植，换盘/换机仍能解析）；
            # 放在其它位置（外置盘）才记绝对路径。
            under_base = (os.path.normcase(os.path.dirname(resolved))
                          == os.path.normcase(DIRS['host_src']))
            self.active_host_db = (os.path.basename(resolved) if under_base
                                   else resolved)
        else:
            self.active_host_db = ''
            apply_active_host_db('')
        self.databases['host'] = db_path('host', 'classify')
        self.save()
        return self.databases['host']

    def set_output_root(self, path):
        """设置自定义输出根（'' = 恢复默认平台目录），持久化并立即生效。"""
        path = (path or '').strip()
        if path and not os.path.isabs(path):
            raise ValueError('输出目录需为绝对路径（如 D:\\我的分析结果）')
        if path:
            _ensure_not_system_dir(path)
        self.output_root = apply_output_root(path)
        self.save()

    def set_examples_root(self, path):
        """设置示例数据根（'' = 自动探测），持久化并立即生效。

        程序 / 数据库 / 示例 三分离打包时，示例目录在程序目录之外，
        这里指定后各工具「✨ 示例」/「👁 示例结果」按钮立即指向新位置。
        """
        path = (path or '').strip()
        if path and not os.path.isabs(path):
            raise ValueError('示例目录需为绝对路径')
        if path and not os.path.isdir(path):
            raise ValueError(f'示例目录不存在: {path}')
        # 配置里只记用户显式指定的路径；留空 = 自动探测（保持 platform.json
        # 可移植，不要把探测出来的本机绝对路径写回去）
        apply_examples_root(path)
        self.examples_root = path
        self.save()

    def tr(self, zh, en):
        """按当前语言取字符串（后端动态文案用）。"""
        return en if self.lang == 'en' else zh

    def _load_extras(self, data):
        """非路径类配置：界面语言/AI 助手与分析默认参数（值均为标量，逐项校验类型）。"""
        try:
            if data.get('language') in ('zh', 'en'):
                self.language = data['language']
            # AI 助手段：缺省值已在 __init__ 预填，这里只合并合法项。
            ai = data.get('ai')
            if isinstance(ai, dict):
                merged = dict(self.ai)
                merged['enabled'] = bool(ai.get('enabled'))
                for k in ('provider', 'base_url', 'api_key', 'model'):
                    if k in ai:
                        merged[k] = str(ai.get(k) or '').strip()
                try:
                    merged['temperature'] = min(2.0, max(
                        0.0, float(ai.get('temperature', merged['temperature']))))
                except (TypeError, ValueError):
                    pass
                try:
                    merged['max_tokens'] = max(64, min(
                        32768, int(ai.get('max_tokens', merged['max_tokens']))))
                except (TypeError, ValueError):
                    pass
                # 供应商档案：白名单清洗（Key 永不从 platform.json 读）
                clean, _seen = [], set()
                for _pr in (ai.get('providers') or []):
                    if not isinstance(_pr, dict):
                        continue
                    _nm = str(_pr.get('name') or '').strip()
                    if not _nm or _nm in _seen:
                        continue
                    _seen.add(_nm)
                    try:
                        _tp = min(2.0, max(0.0, float(
                            _pr.get('temperature', 0.3))))
                    except (TypeError, ValueError):
                        _tp = 0.3
                    try:
                        _mt = max(64, min(32768, int(
                            _pr.get('max_tokens', 2048))))
                    except (TypeError, ValueError):
                        _mt = 2048
                    clean.append({'name': _nm,
                                  'base_url': str(_pr.get('base_url') or '').strip(),
                                  'model': str(_pr.get('model') or '').strip(),
                                  'temperature': _tp, 'max_tokens': _mt})
                merged['providers'] = clean
                merged['active_profile'] = str(
                    ai.get('active_profile') or '').strip()
                # 主动驱动三层开关与闸门超时
                for _bk in ('explain', 'broadcast', 'gate',
                            'exec_auto', 'exec_auto_tune'):
                    if _bk in ai:
                        merged[_bk] = bool(ai.get(_bk))
                if 'gate_timeout_min' in ai:
                    try:
                        merged['gate_timeout_min'] = max(1, min(
                            120, int(ai.get('gate_timeout_min'))))
                    except (TypeError, ValueError):
                        pass
                # 阶段清单：只留字符串、限长；合法性（是否真有该阶段）由
                # broadcast 在使用时按 STAGE_ORDER 过滤——config 不能反向
                # import pipeline（会成环），这里不做阶段名校验。
                for _lk in ('broadcast_stages', 'gate_stages'):
                    if _lk in ai and isinstance(ai.get(_lk), list):
                        vals = [str(x).strip() for x in ai[_lk]
                                if str(x).strip()][:15]
                        merged[_lk] = vals
                # API Key 与代码库分离存储：platform.json 进 git，密钥**绝不**
                # 落进它——真身在 run/ai_local.json（run/ 整体 gitignore），
                # 这里读回合并。platform.json 里即便残留历史 api_key（旧版
                # 直接落盘的），也以密钥文件为准，避免"两处真身"漂移。
                try:
                    lp = _ai_local_path()
                    if os.path.isfile(lp):
                        with open(lp, encoding='utf-8') as _f:
                            _local = json.load(_f)
                        if isinstance(_local, dict):
                            merged['api_key'] = str(_local.get('api_key')
                                                    or '').strip()
                except (OSError, ValueError):
                    pass
                self.ai = merged
            known = {kk: t for kk, t, *_r in self.DEFAULT_FIELDS}
            for k, v in (data.get('defaults') or {}).items():
                if k not in known or v is None:
                    continue
                kind = known[k]
                try:
                    if kind == 'int':
                        self.defaults[k] = int(v)
                    elif kind == 'float':
                        self.defaults[k] = float(v)
                    elif kind == 'bool':
                        self.defaults[k] = bool(v)
                    else:
                        self.defaults[k] = str(v)
                except (TypeError, ValueError):
                    # 用户在设置页填了类型不符的值：留痕，否则「改了不生效」
                    # 无从排查
                    print(f'[config] 跳过类型不合法的分析默认参数 {k}={v!r}',
                          file=sys.stderr)
        except Exception as e:
            # 外层兜底同样不许静默：整个 defaults 段被丢弃时用户会以为
            # 修改生效了
            print(f'[config] 分析默认参数段加载失败，已忽略: {e}',
                  file=sys.stderr)

    # ---- 工具 ---------------------------------------------------
    def tool(self, name):
        path = self.tools.get(name)
        if path and path.startswith('bundled:'):
            return path          # 随包 Python 引擎哨兵，由调用方进程内处理
        if not path or not os.path.isfile(path):
            # 插件只填缺：3rd/ 与 platform.json 都没有时才落到插件提供的副本
            path = self.plugin_tools.get(name)
        if not path or not os.path.isfile(path):
            raise FileNotFoundError(
                f"未找到工具 [{name}]。请在平台配置(platform.json 或 GUI 设置页)"
                "中手动指定路径，或在设置页「插件管理」安装对应插件。")
        return path

    def tool_status(self):
        """GUI 用：返回各工具探测状态（插件提供的键一并在列，来源见
        /api/plugin/list 的 plugin_tools 标记）。"""
        out = dict(self.tools)
        for k, v in self.plugin_tools.items():
            if not out.get(k):
                out[k] = v
        return out

    def refresh_plugins(self):
        """重扫 plugins/ 并挂载插件提供的工具/数据库（填缺，不覆盖既有）。

        安装/卸载插件后调用；启动时在 __init__ 末尾调用一次。"""
        from Virus_Platform_Core import plugin_manager as _pm
        self.plugin_tools, self.plugin_databases = _pm.resolve_plugins()
        # 插件数据库同样只填缺：platform.json 覆盖与注册表解析优先
        for k, v in self.plugin_databases.items():
            cur = self.databases.get(k)
            if not cur or not os.path.isdir(cur):
                self.databases[k] = v
        return self.plugin_tools, self.plugin_databases

    def ensure_dir(self, key):
        d = DIRS[key]
        os.makedirs(d, exist_ok=True)
        return d


_config = None
_config_lock = threading.Lock()


def get_config():
    """配置单例（双检锁）。

    Flask threaded=True 下两个请求可能并发首次触发构造，而 Config.__init__
    会原地重写全局 DIRS（apply_*_root），无锁时可能交错成半更新状态
    （output_root 的子目录只改了一半）。这类失配不报错、只会让功能静默
    失效，必须在入口处串行化。
    """
    global _config
    if _config is None:
        with _config_lock:
            if _config is None:
                _config = Config()
    return _config
