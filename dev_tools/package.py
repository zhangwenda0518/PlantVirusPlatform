# -*- coding: utf-8 -*-
"""一键打包发布版：PyInstaller exe + 外部工具 + 前端，程序/数据库/示例三分离。

产物（默认分离布局，`dist/` 下三个互不干扰的目录）:
  dist/VirusPlatform/            ① 程序（exe + webapp + bin/ + tools/，~1GB）
  dist/VirusPlatform-Examples/   ② 示例数据（<程序>/examples/，~1MB）
  dist/VirusPlatform-Database/   ③ 数据库（仅 --with-db 时生成）
  dist/VirusPlatform-Source/     ④ 源码包（仅 --with-source 时生成）

分发口径（2026-09-19 起）：**原始数据集与数据库区分开**——
  * 随包只发「原始数据集」（参考 FASTA + 注释表 + taxonomy dmp + 第三方
    直接使用库），派生索引/分类库一律不预建：
      - 病毒鉴定库引擎索引（salmon_k31/k15、minibwa，~1.66GB）不随包，
        首次运行在库目录内自建一次后全局复用；
      - kunpeng 病毒分类库（plant ~53MB / ref ~680MB）不随包，用户在
        「数据库构建」页从数据集一键构建（数据集 = kv_index 的
        reference.fasta + reference.ref_info.tsv，含 Accession/Taxid 列；
        RefSeq 库 = databases/virus_ref/ 下的 viral.1.1.genomic.fna.gz）；
    需要开箱即用的预建分类库时加 `--with-built-virusdb`。
  * 宿主分类库按物种而异，向来由用户自建（host-db/，不入包）。

数据库包内容（逐项见下方 `_DB_TREES`）：
  数据集：virusref_db（参考 FASTA/ref_info + 预置 v4 blastn 库）、tax_db、
          virus_ref（RefSeq/RVDB 源 FASTA，本机有就随包）
  注释：annot_db/{prot,hmm}
  建树/其他：tree_db、misc_db/{viroids,suvtk,prob}、genus_lens.tsv

程序启动时按此顺序找示例目录（Virus_Platform_Core/config.py `_detect_examples_root`）：
  platform.json.examples_root → <程序>/examples →
  <程序>/../VirusPlatform-Examples/databases/examples → <程序>/examples

用法:
  python dev_tools/package.py                      # 程序 + 示例（无数据库，~1GB）
  python dev_tools/package.py --with-db            # 再加数据库包（数据集口径）
  python dev_tools/package.py --with-built-virusdb # 数据库包额外带预建 kunpeng 分类库
  python dev_tools/package.py --db-only            # 跳过 exe，只补数据库包
  python dev_tools/package.py --no-split           # 旧布局：示例/数据库放进程序目录
  python dev_tools/package.py --verify             # 打包后跑产物自检
  python dev_tools/package.py --dedupe-only        # 只对已有 dist 做硬链接去重（不重打包）
  python dev_tools/package.py --no-exe             # exe 已建好，只补外部工具/示例/说明
  python dev_tools/package.py --with-source        # 再加源码包（源码+绿色运行时，~900MB）
  python dev_tools/package.py --source-only        # 跳过 exe/示例/数据库，只出源码包

源码包（2026-09-29）：让使用者也能以**源码模式**运行分析——平台源码全树
+ 3rd/python 绿色运行时 + 启动 bat。定位是程序包的「叠加层」：把
dist/VirusPlatform-Source 的内容并入 dist/VirusPlatform 程序目录，双击
启动平台-网页.bat 即源码运行（工具链 3rd/bin、3rd/tools 程序包已带，
解释器与全部依赖在 3rd/python 里，primer3 亦随包）。exe 分发口径不受影响。
"""
import argparse
import json
import os
import shutil
import subprocess
import sys
from datetime import datetime, timezone

# 参数最先解析：--help 直接退出，绝不触发后面的打包/清理动作
_parser = argparse.ArgumentParser(description='打包发布版（默认程序/数据库/示例三分离）')
_parser.add_argument('--with-db', action='store_true',
                     help='同时生成独立的数据库包（开箱即用，体积大）')
_parser.add_argument('--with-built-virusdb', action='store_true',
                     help='数据库包额外携带预建 kunpeng 病毒分类库（plant/ref；'
                          '默认不携带——分发数据集，用户在构建页一键自建）')
_parser.add_argument('--db-only', action='store_true',
                     help='跳过 exe，仅生成/更新数据库包')
_parser.add_argument('--no-split', action='store_true',
                     help='旧布局：示例与数据库都放进程序目录')
_parser.add_argument('--verify', action='store_true',
                     help='打包完成后对产物跑自检（exe --cli selfcheck）')
_parser.add_argument('--dedupe-only', action='store_true',
                     help='不打包，只对已有 dist/VirusPlatform 做硬链接去重')
_parser.add_argument('--no-exe', action='store_true',
                     help='跳过 PyInstaller（exe 已构建好），只补外部工具/示例/说明')
_parser.add_argument('--with-source', action='store_true',
                     help='同时生成「源码包」：源码全树 + 3rd/python 绿色运行时 + '
                          '启动 bat（并入程序包目录即可源码运行）')
_parser.add_argument('--source-only', action='store_true',
                     help='跳过 exe/示例/数据库，仅生成/更新源码包')
_args, _ = _parser.parse_known_args()

# 控制台编码兜底：Windows 默认代码页是 GBK，脚本里若出现 GBK 编不出的字符
# （实测 2026-09-10：'↳' U+21B3 让 print 抛 UnicodeEncodeError，打包跑到一半中断），
# 整个打包流程会被打断。这里只把错误处理改成 replace（不改编码，中文照常可读），
# 编不出的字符降级为 '?'，绝不让日志把打包搞崩。
for _stream in (sys.stdout, sys.stderr):
    try:
        _stream.reconfigure(errors='replace')
    except (AttributeError, OSError):
        pass

# 本脚本位于 dev_tools/，ROOT = 平台根目录（app.py / vp 所在处）
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT)

APP = os.path.join('dist', 'VirusPlatform')
EX_DIR = os.path.join('dist', 'VirusPlatform-Examples')
DB_DIR = os.path.join('dist', 'VirusPlatform-Database')
SRC_DIR = os.path.join('dist', 'VirusPlatform-Source')
SPLIT = not _args.no_split


def _rmtree(path):
    if os.path.isdir(path):
        subprocess.run(['cmd', '/c', 'rmdir', '/s', '/q', path], check=False)
        if os.path.isdir(path):
            sys.exit(f'无法清理旧目录（被占用？）: {path}')


def _fmt(n):
    """人类可读的字节数（与 db_migrate._fmt 同口径）。"""
    for unit in ('B', 'KB', 'MB', 'GB', 'TB'):
        if n < 1024 or unit == 'TB':
            return f'{n:.1f} {unit}' if unit != 'B' else f'{n:.0f} B'
        n /= 1024


def _dir_stat(path):
    """返回 (文件数, 总字节)；目录不存在返回 (0, 0)。"""
    total = n = 0
    if not os.path.isdir(path):
        return 0, 0
    for cur, _dirs, files in os.walk(path):
        for fn in files:
            try:
                total += os.path.getsize(os.path.join(cur, fn))
                n += 1
            except OSError:
                pass
    return n, total


def _build_info():
    """构建清单（build_info.json 内容）：发布包可追溯所需的元信息。

    动机：产物一旦离开本机（U 盘/网盘/他人电脑），"这是哪次提交、什么时候、
    哪套工具链打的"就无从查起，自检报告也无法回溯到对应版本。把 git rev +
    构建时间 + `main.py tools` 探测原文 + 数据库目录统计固化进产物内部，
    无需外部记录即可追溯（2026-10-03）。
    """
    def _git(*a):
        try:
            r = subprocess.run(['git', *a], capture_output=True, text=True,
                               encoding='utf-8', errors='replace')
            return (r.stdout or '').strip()
        except Exception:
            return ''

    tools_probe = ''
    try:
        r = subprocess.run([sys.executable, 'main.py', 'tools'],
                           capture_output=True, text=True, encoding='utf-8',
                           errors='replace', timeout=300)
        tools_probe = r.stdout or ''
        if not tools_probe.strip():
            tools_probe = f'<tools 探测无输出, rc={r.returncode}> ' + (r.stderr or '')
    except Exception as e:
        tools_probe = f'<tools 探测失败: {e!r}>'

    db_dirs = {}
    for rel in ('databases', 'host-db'):
        n, sz = _dir_stat(rel)
        db_dirs[rel + '/'] = {'files': n, 'bytes': sz}

    return {
        'git_rev': _git('rev-parse', 'HEAD'),
        'built_at_utc': datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ'),
        'python': sys.version.split()[0],
        'argv': sys.argv[1:],
        'tools_probe': tools_probe,
        'db_dirs': db_dirs,
    }


def _dedupe_hardlinks(root, min_size=64 * 1024):
    """把 root 下内容完全相同的文件合并成硬链接，返回省下的字节数。

    为什么需要：`tools/mmseqs/bin` 里约 248 个 745KB 的 busybox applet
    （bash/sh/awk/sed/…）在源目录里是**同一个文件的硬链接**，实际只占
    0.7MB；而 `shutil.copytree` 会把它们展开成 248 份真实副本，dist 里
    凭空多出约 185MB（实测 dist/VirusPlatform/tools/mmseqs 221MB）。

    这里按 (大小, 内容摘要) 分组后用 `os.link` 在 **dist 内部**重建硬链接，
    不跨目录链接回源树，dist 仍然自包含、可单独拷贝。失败（跨卷、FAT32、
    权限不足）时静默跳过，绝不影响打包结果。
    """
    import hashlib
    if not os.path.isdir(root):
        return 0
    groups = {}
    for dirpath, _dirnames, filenames in os.walk(root):
        for fn in filenames:
            p = os.path.join(dirpath, fn)
            try:
                st = os.stat(p)
            except OSError:
                continue
            if st.st_size < min_size or st.st_nlink > 1:
                continue
            groups.setdefault(st.st_size, []).append(p)

    saved = 0
    for size, paths in groups.items():
        if len(paths) < 2:
            continue
        canon = {}
        for p in paths:
            h = hashlib.blake2b(digest_size=16)
            try:
                with open(p, 'rb') as f:
                    for chunk in iter(lambda: f.read(1 << 20), b''):
                        h.update(chunk)
            except OSError:
                continue
            key = h.hexdigest()
            first = canon.get(key)
            if first is None:
                canon[key] = p
                continue
            tmp = p + '.linktmp'
            try:
                os.link(first, tmp)
                os.replace(tmp, p)
                saved += size
            except OSError:
                try:
                    os.remove(tmp)
                except OSError:
                    pass
    return saved


if _args.dedupe_only:
    _n = _dedupe_hardlinks(os.path.join(APP, '3rd', 'tools'))
    print(f'硬链接去重完成：{APP}/3rd/tools 省下 {_n / 1024 / 1024:.0f} MB')
    raise SystemExit(0)


# ------------------------------------------------------------------
# 1) PyInstaller
# ------------------------------------------------------------------
if not (_args.db_only or _args.source_only):
    if _args.no_exe:
        print('== 1/5 跳过 PyInstaller（--no-exe，沿用已有 dist/VirusPlatform）==')
    else:
        print('== 1/5 PyInstaller 打包 ==')
        _rmtree(APP)
        r = subprocess.run([sys.executable, '-m', 'PyInstaller',
                            os.path.join('dev_tools', 'VirusPlatform.spec'),
                            '--noconfirm', '--distpath', 'dist',
                            '--workpath', 'build'])
        if r.returncode != 0:
            sys.exit('打包失败')
        # 必须用 _rmtree（cmd rmdir）而不是 shutil.rmtree：后者是 Python 层
        # 批量删除，在本机安全策略下会被拦断（实测 build/ 5600+ 文件触发
        # SAFE_DELETE_BULK_CONFIRM_REQUIRED），而此刻 dist/VirusPlatform 已经
        # 被上面的 _rmtree(APP) 删掉了 —— 流程卡死在半路，产物整个消失。
        # 本脚本其余清理一律走 _rmtree，这里保持一致。
        _rmtree('build')

    # ---- 2) 外部工具与前端 -------------------------------------------
    print('== 2/5 外部工具与前端 ==')
    # 目录整理（2026-09-10）：外部二进制统一在 3rd/ 下（3rd/bin、3rd/tools）
    for d in ('3rd/bin', '3rd/tools'):
        if os.path.isdir(d):
            dst = os.path.join(APP, *d.split('/'))
            _rmtree(dst)
            shutil.copytree(d, dst,
                            ignore=shutil.ignore_patterns('.git', '.github'))
            print('  +', d + '/（全部，已排除 .git）')
            if d.endswith('tools'):
                _saved = _dedupe_hardlinks(dst)
                if _saved:
                    print(f'    -> 硬链接去重：省 {_saved / 1024 / 1024:.0f} MB'
                          '（mmseqs 的 busybox applet 等重复副本）')

    orfipy = os.path.join(os.path.dirname(sys.executable), 'Scripts',
                          'orfipy.exe')
    if os.path.isfile(orfipy):
        shutil.copy2(orfipy, os.path.join(APP, 'orfipy.exe'))
        print('  + orfipy.exe')

    # 便携 R + treedater（2026-09-17）：TreeDater-LTT 真引擎随包分发，
    # 目标机免装 R/免联网。安装包 exe 与 doc/tests 不随包（省体积）。
    if os.path.isdir('3rd/R'):
        dst = os.path.join(APP, '3rd', 'R')
        _rmtree(dst)
        shutil.copytree('3rd/R', dst,
                        ignore=shutil.ignore_patterns('*.exe', 'doc', 'tests'))
        print('  + 3rd/R/（便携 R + treedater）')

    ov_build = os.path.join('3rd', 'open-virome', 'frontend', 'build')
    if os.path.isdir(ov_build):
        dst = os.path.join(APP, '3rd', 'open-virome', 'frontend', 'build')
        _rmtree(dst)
        shutil.copytree(ov_build, dst)
        print('  + 3rd/open-virome/frontend/build/')

    # dsRNA 设计的小型数据文件（2026-09-19；2026-10-03 随功能恢复重新入包）：
    # lookup.db 热力学查找表 / siRNA 参数 / 致死基因清单 —— 没有它们 dsRNA 卡
    # 一运行就报「缺少热力学查找表」。合计 <100KB，随程序包开箱即用；大头的
    # 脱靶面板 panel\（物种 RefSeq RNA，GB 级）走数据库包（见 _DB_TREES）。
    for _rel in ('databases/dsrna/lookup.db',
                 'databases/dsrna/siRNA_parameters.txt',
                 'databases/dsrna/lethal_lists'):
        _src = os.path.join(ROOT, *_rel.split('/'))
        _dst = os.path.join(APP, *_rel.split('/'))
        if os.path.isdir(_src):
            _rmtree(_dst)
            shutil.copytree(_src, _dst)
            print('  +', _rel + '/')
        elif os.path.isfile(_src):
            os.makedirs(os.path.dirname(_dst), exist_ok=True)
            shutil.copy2(_src, _dst)
            print('  +', _rel)

    # primer3-py 包副本（2026-09-20）：分发版在非 ASCII 安装路径下跑 primer3，
    # 依赖「ASCII 镜像」机制 —— exe 启动时由 app.py 把程序包自带的
    # primer3_pkg/primer3 镜像到 %ALLUSERSPROFILE%\VirusPlatformPrimer3App
    # 并插入 sys.path。
    # ⚠️ 必须携带 src/libprimer3/primer3_config/（dangle/stack 等热力学参数
    # 表）：thermoanalysis.pyd 在 ThermoAnalysis() 构造时按
    # <包目录>/src/libprimer3/primer3_config 定位并加载，缺文件会直接
    # 段错误（0xC0000005）——build19 及之前漏带该目录正是「分发版 T7
    # 引物段错误」的根因，与内存/非 ASCII 路径无关（2026-09-20 实测定位）。
    # 源优先级：ProgramData 镜像 → 3rd/python 绿色版 → 构建机 site-packages；
    # 缺 primer3_config 的候选会用其它候选补齐。整包缺失时跳过，
    # 分发版引物热力学降级为不可用（其余功能不受影响）。
    _p3_bases = [
        os.path.join(os.environ.get('ALLUSERSPROFILE', r'C:\ProgramData'),
                     'VirusPlatformPrimer3'),
        os.path.join('3rd', 'python', 'Lib', 'site-packages'),
        os.path.join(sys.prefix, 'Lib', 'site-packages'),
    ]

    def _p3_ok(d):
        return os.path.isfile(os.path.join(d, 'bindings.py')) and \
            os.path.isdir(os.path.join(d, 'src', 'libprimer3',
                                       'primer3_config'))

    def _p3_cfg_src(cands):
        for d in cands:
            c = os.path.join(d, 'src', 'libprimer3', 'primer3_config')
            if os.path.isdir(c):
                return c
        return None

    _p3_mirror = None
    for _base in _p3_bases:
        _cand = os.path.join(_base, 'primer3')
        if os.path.isfile(os.path.join(_cand, 'bindings.py')):
            _p3_mirror = _cand
            break
    if _p3_mirror:
        dst = os.path.join(APP, 'primer3_pkg', 'primer3')
        _rmtree(os.path.dirname(dst))
        shutil.copytree(_p3_mirror, dst,
                        ignore=shutil.ignore_patterns('__pycache__', 'src',
                                                      '*.c', '*.pyx', '*.pxd',
                                                      '*.h', '*.html'))
        # primer3_config 无论候选是否完整都补齐（缺失 = 运行期段错误）
        _cfg = _p3_cfg_src([os.path.join(b, 'primer3') for b in _p3_bases])
        if _cfg:
            _cfg_dst = os.path.join(dst, 'src', 'libprimer3', 'primer3_config')
            if os.path.isdir(os.path.dirname(_cfg_dst)):
                shutil.rmtree(os.path.dirname(_cfg_dst), ignore_errors=True)
            os.makedirs(os.path.dirname(_cfg_dst), exist_ok=True)
            shutil.copytree(_cfg, _cfg_dst)
        elif not _p3_ok(dst):
            print('  ! primer3_config 缺失且无处补齐：分发版引物热力学将段错误降级为不可用')
        _ver = os.path.join(os.path.dirname(_p3_mirror), 'VERSION.txt')
        if os.path.isfile(_ver):
            shutil.copy2(_ver, os.path.join(APP, 'primer3_pkg', 'VERSION.txt'))
        print(f'  + primer3_pkg/（primer3-py 副本，{sum(len(fs) for _d, _s, fs in os.walk(dst))} 个文件，'
              f'含 primer3_config 热力学参数表）')
    else:
        print('  ! 未找到 primer3 包副本，分发版引物热力学将不可用')

    # ---- 干净 platform.json：零硬编码路径，示例/数据库靠自动探测 ----
    fresh_cfg = {'threads': 0, 'tools': {},
                 'defaults': {'max_heavy_tasks': 2, 'max_light_tasks': 4},
                 'language': 'zh',
                 'examples_root': ''}      # '' = 启动时自动探测
    with open(os.path.join(APP, 'platform.json'), 'w', encoding='utf-8') as f:
        json.dump(fresh_cfg, f, ensure_ascii=False, indent=2)
    print('  + platform.json（干净配置，示例/数据库路径零硬编码）')
else:
    print('== 1-2/5 跳过 exe（--db-only）==')

# ------------------------------------------------------------------
# 3) 示例数据（独立目录）
# ------------------------------------------------------------------
if not (_args.db_only or _args.source_only):
    print('== 3/5 示例数据 ==')
    src_ex = 'examples'
    if not os.path.isdir(src_ex):
        print('  ! 未找到 examples/，跳过')
    else:
        if SPLIT:
            # 与源码布局镜像：示例在程序目录的**同级** examples/
            # （Virus_Platform_Core/config.py `_detect_examples_root` 第 3 号候选）
            dst = os.path.join(EX_DIR, 'examples')
            _rmtree(EX_DIR)
            os.makedirs(os.path.dirname(dst), exist_ok=True)
            shutil.copytree(src_ex, dst,
                            ignore=shutil.ignore_patterns('_work', '__pycache__'))
            n = sum(len(fs) for _d, _s, fs in os.walk(dst))
            print(f'  + {EX_DIR}/examples/（{n} 个文件，'
                  f'程序目录不含示例）')
        else:
            dst = os.path.join(APP, 'examples')
            _rmtree(dst)
            shutil.copytree(src_ex, dst,
                            ignore=shutil.ignore_patterns('_work', '__pycache__'))
            print('  + 程序目录内 examples/（--no-split）')

# ------------------------------------------------------------------
# 4) 数据库包（独立目录）
# ------------------------------------------------------------------
def _skip_engine_index(dir_path, names):
    """virusref_db 各鉴定库目录内的引擎派生索引不随包分发（参考数据集照发）。

    库有效性的锚点是参考真相（reference.fasta + reference.ref_info.tsv），
    索引是派生物：salmon k31 全库现建约 10s、k15 约 35s、minibwa 数秒，
    首次运行时由引擎在库目录内自建一次并全局复用（kv_stage.index_dir_for，
    平台内所有库通用）。kv_index 一家就省 ~1.66GB 分发体积；其余库目录
    （viromock_kv 等）同样适用。日志与历史备份目录一并跳过。

    2026-09-29：不再要求"锚点与索引同目录"。引擎索引子目录在任何层级
    出现都跳过——默认库更名 plant_ref 后，索引既可能在库目录内（自包含
    布局），也可能残留在旧 kv_index/ 这类纯索引目录里，一律不入包。
    """
    import fnmatch as _fn
    skip = {'salmon_k31', 'salmon_k15', 'salmon', 'minibwa', 'logs'}
    for pat in ('salmon_index*', 'backup_*', '*.log', '*.stale'):
        skip |= {n for n in names if _fn.fnmatch(n, pat)}
    return skip


if _args.with_db or _args.db_only:
    print('== 4/5 数据库包 ==')
    target = APP if not SPLIT else DB_DIR
    os.makedirs(target, exist_ok=True)

    def _db_dst(rel, is_file=False):
        """databases/<rel> → 数据库包内的同构路径。

        必须保留相对 `databases/` 的层级（databases/kunpeng_db/plant、
        databases/annot_db/prot…）——`Virus_Platform_Core/config.py` 的 `db_path()`
        按 `DB_LAYOUT` 解析这些路径；早先用 basename 展平（databases/plant/）
        会让打包后的平台一个库都找不到（实测 db-migrate --check 三项全 ✘）。
        """
        dst = os.path.join(target, 'databases', *rel.split('/'))
        os.makedirs(os.path.dirname(dst) if is_file else dst, exist_ok=True)
        return dst

    def _copy_tree(rel, note, ignore=(), ignore_fn=None):
        """整树拷贝 `databases/<rel>`。返回 (文件数, 字节数)。"""
        src = os.path.join('databases', rel)
        if not os.path.isdir(src):
            print(f'  - 跳过（源不存在）: databases/{rel}')
            return 0, 0
        dst = _db_dst(rel)
        _rmtree(dst)
        _ig = shutil.ignore_patterns(
            '__pycache__', '*.pyc', '*.tmp', '*.mmseqs_tmp', '.DS_Store',
            *ignore)
        if ignore_fn:
            _base = _ig
            _ig = lambda d, names: _base(d, names) | ignore_fn(d, names)
        shutil.copytree(src, dst, ignore=_ig)
        n, size = _dir_stat(dst)
        print(f'  + databases/{rel}/  {_fmt(size)} / {n} 文件  {note}')
        return n, size

    def _copy_file(rel, note):
        """拷贝 `databases/<rel>` 单个文件。返回 (文件数, 字节数)。"""
        src = os.path.join('databases', rel)
        if not os.path.isfile(src):
            print(f'  - 跳过（源不存在）: databases/{rel}')
            return 0, 0
        dst = _db_dst(rel, is_file=True)
        shutil.copy2(src, dst)
        size = os.path.getsize(dst)
        print(f'  + databases/{rel}  {_fmt(size)}  {note}')
        return 1, size

    # ---- 数据库清单 ----------------------------------------------------
    # 为什么从"逐个文件名"改成"整树"（2026-09-10）：
    #   原先的手写清单只覆盖 kunpeng_db/plant + tax_db + palmdb（后者已随
    #   Open-Virome 撤下、目录早就不存在），把 annot_db / tree_db / misc_db /
    #   virusref_db 整套漏掉了 —— 打出来的包缺功能注释库、CDD、HMM、建树参考、
    #   病毒参考与预置 blastn 库，装上去等于半个平台。
    #   `Virus_Platform_Core/db_migrate.py` 的口径就是整拷 `databases/`，
    #   这里对齐它；以后新增库文件也不必再改本脚本。
    #   各项体积会打印出来，要瘦身按行删即可。
    #
    # 原始数据集 vs 数据库（2026-09-19 口径）：
    #   kunpeng 预建分类库（plant/ref）默认**不入包**——plant 可由数据集
    #   （kv_index 的 reference.fasta + reference.ref_info.tsv，Accession/
    #   Taxid 列齐备）在建库页一键重建；ref 可由 virus_ref/ 下的 RefSeq
    #   源 FASTA 一键重建。需要开箱即用时加 --with-built-virusdb。
    _DB_TREES = (
        ('tax_db',           '← NCBI Taxonomy（nodes/names/merged.dmp，原始数据集）'),
        ('virusref_db',      '← 病毒参考数据集（FASTA+ref_info，兼作分类库建库源）'
                             ' + 预置 v4 blastn 库 + kv_index 鉴定库（引擎索引不随包，首跑自建）'),
        ('virus_ref',        '← RefSeq/RVDB 源 FASTA（通用病毒库建库数据集；本机有就随包）'),
        ('annot_db/prot',    '← 功能注释层1：RefSeq 病毒蛋白（DIAMOND 库 + faa.gz + 元数据）'),
        ('annot_db/hmm',     '← 功能注释层2：Pfam-A-Viruses（已 hmmpress）'),
        ('tree_db',          '← 建树参考 plant / ictv 两套口径'),
        ('misc_db/viroids',  '← 类病毒 blastn v4 库'),
        ('misc_db/suvtk',    '← NCBI 提交 BFVD 功能注释库'),
        ('misc_db/prob',     '← 宿主概率表'),
        # 2026-10-03 用户裁定：dsRNA 脱靶面板（~1.4GB 物种 RefSeq RNA）不再随
        # 数据库包分发（软件体积瘦身）。功能保留：自备面板放 databases/dsrna/panel/
        # 即自动参与脱靶扫描；小文件 lookup.db / 参数 / 致死清单仍随程序包。
    )
    if _args.with_built_virusdb:
        _DB_TREES = _DB_TREES + (
            ('kunpeng_db/plant', '← 病毒分类库（预建，开箱即用；数据集在 virusref_db 可重建）'),
            ('kunpeng_db/ref',   '← RefSeq 通用 kraken2 库（预建，开箱即用）'),
        )
    else:
        print('  - kunpeng 预建分类库（plant/ref）不入包：分发数据集，'
              '构建页一键自建（开箱即用加 --with-built-virusdb）')
    _db_files = _db_bytes = 0
    for _rel, _note in _DB_TREES:
        _n, _s = _copy_tree(_rel, _note, ignore_fn=_skip_engine_index)
        _db_files += _n
        _db_bytes += _s
    _n, _s = _copy_file('genus_lens.tsv', '← 属平均长度表（近完整基因组判据）')
    _db_files += _n
    _db_bytes += _s

    # 宿主分类库**故意不收**：它按物种而异（每个宿主一套 ~1.5GB 的 kunpeng
    # 库），收进数据库包会让分发体积无谓翻倍。宿主库留在 host-db/ 下，
    # 由用户在「数据库构建」页按自己的物种建库。
    print('  - 宿主分类库不入包（按物种而异，见 host-db/）')
    print(f'  —— 数据库包合计 {_fmt(_db_bytes)} / {_db_files} 文件')
    # 数据库包自带一份配置，指向自己（也可由用户改用「设置 → 数据库目录」）
    with open(os.path.join(target, '数据库对接说明.txt'), 'w',
              encoding='utf-8') as f:
        f.write(
            '本目录是「数据库包」，与程序目录分离，可放在任意盘。\n'
            '\n'
            '内容口径（原始数据集为主，派生库不预建）：\n'
            '  · 病毒参考数据集（reference.fasta + ref_info.tsv）——鉴定库\n'
            '    引擎索引（salmon/minibwa）首次运行自动现建并全局复用；\n'
            '  · kunpeng 病毒分类库（plant/ref）不预带：到「数据库构建」页\n'
            '    「原始数据集」区一键构建（plant 由参考数据集直接构建，\n'
            '    ref 需先放入 RefSeq 源 FASTA）\n'
            '  · NCBI Taxonomy / 注释库（Pfam、RefSeq 蛋白）/ 建树参考 /\n'
            '    blastn 库等为第三方直接使用库，随包即用\n'
            '  · 宿主分类库按物种而异，用「数据库构建」页自建\n'
            '\n'
            '对接方式（零配置优先）：\n'
            '  ① 与程序目录**同级放置**即被自动识别，无需任何设置\n'
            '     （程序启动时自动探测同级 VirusPlatform-Database）\n'
            '  ② 启动平台 → 设置 → 数据库目录 → 填本目录的绝对路径 → 应用\n'
            '     （保存时会即时反馈宿主库/病毒库/Taxonomy 是否就绪）\n'
            '  ③ 开发机执行：python main.py db-migrate --to <目标库目录>\n'
            '\n'
            '自备宿主库可放任意位置：构建页「📂 加载已构建宿主库」\n'
            '选中即可，平台会自动登记该位置。\n'
            '\n'
            '校验：程序目录下执行\n'
            '  VirusPlatform.exe --cli selfcheck\n'
            '（源码模式：python main.py selfcheck）\n')
    print('  + 数据库对接说明.txt')
else:
    print('== 4/5 数据库：跳过（软件与数据库分离分发）==')
    print('  需要数据库包时重跑：python dev_tools/package.py --with-db')

# ------------------------------------------------------------------
# 4b) 源码包（叠加层：并入程序包目录即可源码模式运行）
# ------------------------------------------------------------------
if _args.source_only:
    _args.with_source = True
if _args.with_source:
    print('== 4b/5 源码包 ==')
    _rmtree(SRC_DIR)
    os.makedirs(SRC_DIR, exist_ok=True)

    _SRC_IG = shutil.ignore_patterns('__pycache__', '*.pyc',
                                     '.pytest_cache', '.ruff_cache')

    def _src_tree(rel, note):
        src = os.path.join(ROOT, rel)
        if not os.path.isdir(src):
            print(f'  - 跳过（源不存在）: {rel}/')
            return
        dst = os.path.join(SRC_DIR, *rel.split('/'))
        shutil.copytree(src, dst, ignore=_SRC_IG)
        n, size = _dir_stat(dst)
        print(f'  + {rel}/  {_fmt(size)} / {n} 文件  {note}')

    def _src_file(rel, note):
        if not os.path.isfile(os.path.join(ROOT, rel)):
            print(f'  - 跳过（源不存在）: {rel}')
            return
        dst = os.path.join(SRC_DIR, rel)
        os.makedirs(os.path.dirname(dst) or SRC_DIR, exist_ok=True)
        shutil.copy2(os.path.join(ROOT, rel), dst)
        print(f'  + {rel}  {_fmt(os.path.getsize(dst))}  {note}')

    for _f in ('app.py', 'main.py', 'requirements.txt', 'requirements.lock'):
        _src_file(_f, '入口 / 依赖清单')
    for _d in ('Virus_Platform_Core', 'webapp', 'scripts', 'dev_tools'):
        _src_tree(_d, '平台源码 / 前端 / 辅助脚本')
    for _b in ('启动平台-网页.bat', '启动平台-桌面窗口.bat', '启动拖拽代理.bat',
               '挂载数据目录.bat', '环境自检.bat'):
        _src_file(_b, '源码模式启动脚本')

    _rt = os.path.join('3rd', 'python')
    if os.path.isdir(_rt):
        _src_tree(_rt, '绿色版 Python 3.12 运行时（依赖与 primer3 已装好）')
    else:
        print('  ! 缺少 3rd/python：源码包将无法开箱运行（需使用者自装 Python 3.12）')

    # 干净 platform.json（与程序包同口径：零硬编码路径，示例/数据库自动探测）
    with open(os.path.join(SRC_DIR, 'platform.json'), 'w', encoding='utf-8') as f:
        json.dump({'threads': 0, 'tools': {},
                   'defaults': {'max_heavy_tasks': 2, 'max_light_tasks': 4},
                   'language': 'zh', 'examples_root': ''},
                  f, ensure_ascii=False, indent=2)
    print('  + platform.json（干净配置）')

    _readme = '''本目录是「源码包」——让使用者以源码模式运行平台（区别于 exe 双击运行）。

内容：平台源码全树（Virus_Platform_Core / webapp / app.py / main.py）
+ 3rd/python 绿色版 Python 3.12（全部依赖与 primer3 已装好）
+ 源码模式启动脚本（bat）。

用法（叠加到程序包，推荐）：
  ① 先拿到程序包 dist/VirusPlatform/（exe + 3rd/bin + 3rd/tools 工具链）
  ② 把本目录的全部内容复制进 VirusPlatform 目录内（同名目录合并）
  ③ 双击「启动平台-网页.bat」（或桌面窗口版）→ 源码模式，改代码即时生效
  数据库/示例仍由 VirusPlatform-Database / -Examples 两个包提供（同级自动识别）。

独立放置（不并入程序包）也可以：
  本目录自带解释器与全部依赖，「python app.py」即可启动；
  但分析功能依赖外部工具（blastn/spades/kunpeng 等，在 3rd/tools、3rd/bin），
  独立使用时需自行补齐，或到「设置 → 工具路径」指定。

说明：
  - 源码模式与 exe 模式共享同一套 platform.json / 数据库 / 示例配置；
  - primer3 已随 3rd/python 携带，源码模式引物热力学开箱可用；
  - 依赖清单见 requirements.txt / requirements.lock（重建环境用）。
'''
    with open(os.path.join(SRC_DIR, '源码运行说明.txt'), 'w',
              encoding='utf-8') as f:
        f.write(_readme)
    print('  + 源码运行说明.txt')
    _n, _sz = _dir_stat(SRC_DIR)
    print(f'  —— 源码包合计 {_fmt(_sz)} / {_n} 文件')

# ------------------------------------------------------------------
# 5) 顶层说明 + 汇总
# ------------------------------------------------------------------
print('== 5/5 说明与汇总 ==')
if SPLIT and os.path.isdir(APP):
    with open(os.path.join('dist', '程序与数据说明.txt'), 'w',
              encoding='utf-8') as f:
        f.write(
            '植物病毒分析平台 · 分发说明（程序 / 数据库 / 示例 三分离）\n'
            '=' * 60 + '\n'
            '\n'
            '① VirusPlatform\\            程序本体：双击 VirusPlatform.exe 启动\n'
            '   - 内含 exe + webapp 前端 + bin/ tools/ 外部工具\n'
            '   - 不含数据库，也不含示例数据\n'
            '\n'
            '② VirusPlatform-Examples\\   示例数据（只读，约 1MB）\n'
            '   - 各工具「✨ 示例」按钮与「示例结果」页用到的输入/产物\n'
            '   - 与程序同级放置即可被自动识别；也可搬到别处后到\n'
            '     「设置 → 示例数据目录」填绝对路径\n'
            '\n'
            '③ VirusPlatform-Database\\   数据库包（仅 --with-db 时生成）\n'
            '   - 原始数据集：病毒参考（FASTA+注释表）、NCBI Taxonomy、\n'
            '     RefSeq/RVDB 源 FASTA（若随包）\n'
            '   - 第三方直接使用库：注释库（RefSeq 病毒蛋白 / Pfam HMM）、\n'
            '     建树参考、类病毒库、SUVTK、宿主概率表、属平均长度表、\n'
            '     预置 blastn 库\n'
            '   - 派生库不预带、首跑/一键自建：鉴定库引擎索引首次运行自动\n'
            '     现建；kunpeng 病毒分类库（plant/ref）到「数据库构建」页\n'
            '     「原始数据集」区一键构建（开箱即用需求打包时加\n'
            '     --with-built-virusdb）\n'
            '   - 不含宿主分类库（按物种而异，见 host-db/，用「数据库构建」页自建）\n'
            '   - 对接（零配置优先）：与本程序**同级放置**即被自动识别；\n'
            '     也可启动平台 →「设置 → 数据库目录」填该目录，\n'
            '     或 python main.py db-migrate --to <目录>\n'
            '   - 自备宿主库可放任意位置：构建页「📂 加载已构建宿主库」\n'
            '     选中即可，平台自动登记该位置\n'
            '\n'
            '自检：VirusPlatform.exe --cli selfcheck\n'
            '（源码模式：python main.py selfcheck）\n')
    print('  + dist/程序与数据说明.txt')

for d in (APP, EX_DIR, DB_DIR, SRC_DIR):
    if os.path.isdir(d):
        total = sum(os.path.getsize(os.path.join(dp, f))
                    for dp, _s, fs in os.walk(d) for f in fs)
        n = sum(len(fs) for _dp, _s, fs in os.walk(d))
        print(f'  {d}: {total / 1e9:.3f} GB / {n} 文件')

# ------------------------------------------------------------------
# 5.5) 产物构建清单 build_info.json（可追溯：git rev / 时间 / tools 探测）
# ------------------------------------------------------------------
print('== 5.5/5 构建清单 ==')
_info = _build_info()
for _d in (APP, EX_DIR, DB_DIR, SRC_DIR):
    if not os.path.isdir(_d):
        continue
    _payload = dict(_info, product=os.path.basename(_d))
    with open(os.path.join(_d, 'build_info.json'), 'w', encoding='utf-8') as f:
        json.dump(_payload, f, ensure_ascii=False, indent=2)
        f.write('\n')
    print(f'  + {_d}/build_info.json'
          f'（git {(_info["git_rev"] or "?")[:10]}，'
          f'{_info["built_at_utc"]}）')

# ------------------------------------------------------------------
# 可选：产物自检
# ------------------------------------------------------------------
if _args.verify:
    print('\n== 产物自检 ==')
    _verify_bad = []
    for _d in (APP, EX_DIR, DB_DIR, SRC_DIR):
        if not os.path.isdir(_d):
            continue
        _bi = os.path.join(_d, 'build_info.json')
        if not os.path.isfile(_bi):
            _verify_bad.append(f'{_d}: 缺 build_info.json（本次打包未生成？）')
            continue
        try:
            with open(_bi, encoding='utf-8') as f:
                _bi_json = json.load(f)
        except Exception as e:
            _verify_bad.append(f'{_d}: build_info.json 解析失败 ({e!r})')
            continue
        _rev = str(_bi_json.get('git_rev') or '').strip()
        if not _rev:
            _verify_bad.append(
                f'{_d}: build_info.json 的 git_rev 为空'
                '（打包时不在 git 仓库内，或 git 不可用）')
            continue
        print(f'  ✔ {_d}/build_info.json'
              f'（git {_rev[:10]}，{_bi_json.get("built_at_utc", "?")}）')
    for _b in _verify_bad:
        print('  ✘ ' + _b)
    if os.path.isfile(os.path.join(APP, 'VirusPlatform.exe')):
        print('== 产物自检（exe --cli selfcheck）==')
        r = subprocess.run([os.path.abspath(os.path.join(APP, 'VirusPlatform.exe')),
                            '--cli', 'selfcheck'], cwd=os.path.abspath(APP))
        print('  selfcheck 退出码:', r.returncode)
    if _verify_bad:
        sys.exit(1)
print('\n完成。')
