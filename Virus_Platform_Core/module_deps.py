# -*- coding: utf-8 -*-
"""运行前依赖自愈：模块要用的外部工具缺失时，自动从插件源安装。

设计要点
--------
1. **每个模块声明自己要用哪些工具键**（MODULE_TOOLS），这是唯一真源：
   工具页与管道阶段都从这张表取需求，不再各自散落判断。
2. **先检测、再自愈、后运行**：`ensure_tools()` 在任务真正启动前调用；
   缺失项若能由插件源提供，就当场下载安装（带 sha256 校验），装完重扫
   `Cfg.tools` 再复核一遍。装不上才抛错——错误信息里给出可操作指引，
   而不是让任务跑到一半才 FileNotFoundError。
3. **全局可见**：自愈结果写入 `run/plugin_autofix.json`，`selfcheck` 的
   「外部工具」段会读它，把「本次会话自动装过什么」显式列给用户，
   避免"插件被悄悄装上、用户不知道平台改了什么"。

为什么不放在每个 job 里各写一遍
------------------------------
30 个 job × 各自的引擎，散着判断必然漏。集中在一张表 + 一个入口，
新增工具页时只加一行映射；漏声明由 tests/test_module_preflight.py 兜住。
"""
import json
import os
import threading
import time

from Virus_Platform_Core import plugin_manager as pm
from Virus_Platform_Core.config import DIRS, PLATFORM_ROOT

# ----------------------------------------------------------------------
# 模块 → 工具键（唯一真源）
# 键与 Cfg.tools / 插件 provides 一致；工具页键与 pipeline 阶段键都在这。
# ----------------------------------------------------------------------
MODULE_TOOLS = {
    # ── 工具页（TOOL_REGISTRY 的键）──────────────────────────────
    'convert': ['seqkit'],
    'fastp': ['fastp'],
    'hostremoval': ['seqkit', 'minibwa', 'kunpeng'],
    'orf': ['orfipy'],
    'orfa': ['blastn', 'makeblastdb'],
    'genoplot': ['gbdraw'],
    'primer': ['blastn', 'makeblastdb'],
    'identify': ['kunpeng', 'seqkit'],
    'assemble': ['spades', 'makeblastdb', 'blastn'],
    'contigs': ['kunpeng', 'seqkit'],
    'sdt': ['mafft'],
    'identity': ['mafft'],
    'structcmp': ['mafft', 'trimal', 'fasttree'],
    'verify': ['blastn', 'diamond', 'mmseqs'],
    'consensus': ['samtools', 'minimap2'],
    'kvsuite': ['minimap2', 'samtools', 'seqkit'],
    'align': ['mafft', 'trimal'],
    'quicktree': ['mafft', 'trimal', 'fasttree', 'raxml-ng'],
    'rdp': ['RDP5CL'],
    'phylogeo': ['lsd2'],
    'dsrna': ['dsRNAmax'],
    'crispr': ['blastn', 'makeblastdb'],
    'mirna': ['srna-target'],
    'phylodyn': ['lsd2'],          # rscript 是 TreeDater 可选分支，见 OPTIONAL_TOOLS
    'beastexport': [],
    'conan': [],
    'rtt': [],
    'pdprep': [],
    # pdrename / pdgroup 已于 2026-09-21 归档（注册表已摘），不再登记
    # ── 一键链式（web/chains.py）：各自是上面若干工具页的组合 ─────
    # virchain = ②提取序列组装 + ③组装后鉴定 + ④候选验证
    'virchain': ['kunpeng', 'seqkit', 'spades', 'makeblastdb',
                 'blastn', 'diamond', 'mmseqs'],
    # kvchain = kvsuite 全流程（鉴定→过滤→共识→绘图→变异）
    'kvchain': ['minimap2', 'samtools', 'seqkit'],
    # prepchain = 转换 → 质控 → 转 FASTA → 宿主去除
    'prepchain': ['seqkit', 'fastp', 'minibwa', 'kunpeng'],
    # annochain = ORF 预测 → 功能注释 → 基因组图谱
    'annochain': ['orfipy', 'blastn', 'makeblastdb', 'gbdraw'],
    # ── 管道阶段（pipeline.STAGE_NAMES 的键）─────────────────────
    'subsample': ['seqkit'],
    'fq2fa': ['seqkit'],
    'host': ['seqkit', 'minibwa', 'kunpeng'],
    'hostana': ['kunpeng'],
    # assembly 阶段内部还会用 makeblastdb 建骨架库（见 assembly.py）
    'assembly': ['spades', 'makeblastdb', 'blastn'],
    'phylo': ['mafft', 'trimal', 'fasttree', 'raxml-ng'],
    'gbdraw': ['gbdraw'],
    'report': [],
}

# 可选工具：缺了只降级对应子功能，不该阻断整个模块。
# 与 MODULE_TOOLS 分开：那里是「没有就跑不了」，这里是「没有就少一个分支」。
# 格式 {模块: {工具键: 缺了之后的影响}}——影响文案直接进日志/自检，让用户
# 知道少了什么能力，而不是莫名其妙看不到某个结果。
OPTIONAL_TOOLS = {
    'phylodyn': {
        'rscript': 'TreeDater 定年分支不可用（需 R + treedater 包）；'
                   'LSD2 定年与其余分析不受影响',
    },
    'genoplot': {
        'gbdraw': 'gbdraw 引擎不可用，自动回退纯 Python 引擎 dna_features_viewer'
                  '（出图样式略简，功能不受影响）',
    },
    'kvsuite': {
        'salmon': 'salmon 定量引擎不可用，需改用 minibwa 比对引擎',
    },
}

# 插件源里没有、只能人工装或本机自带的键（自愈时直接给出指引）
_NO_PLUGIN_HINT = {
    'rscript': 'R 脚本（LSD2 定年辅助）：请安装 R 并加入 PATH，'
               '或在 platform.json 的 tools.rscript 指定路径',
    'gbdraw': 'gbdraw：pip install git+https://github.com/satoshikawato/gbdraw.git'
              '（未装时基因组图会自动回退纯 Python 引擎，不影响运行）',
    'orfipy': 'orfipy：pip install orfipy（随 requirements.txt 安装）',
}

_AUTOFIX_LOG = os.path.join(DIRS.get('run') or
                            os.path.join(PLATFORM_ROOT, 'run'),
                            'plugin_autofix.json')
_lock = threading.Lock()


def _log_path():
    d = os.path.dirname(_AUTOFIX_LOG)
    if d:
        os.makedirs(d, exist_ok=True)
    return _AUTOFIX_LOG


def record_autofix(module, installed, failed):
    """把一次自愈结果落盘（供 selfcheck / 设置页展示）。"""
    if not installed and not failed:
        return
    with _lock:
        p = _log_path()
        try:
            data = json.load(open(p, encoding='utf-8')) if os.path.isfile(p) else {}
        except (OSError, ValueError):
            data = {}
        rec = data.setdefault('modules', {})
        rec[module] = {'installed': installed, 'failed': failed,
                       'at': time.strftime('%Y-%m-%d %H:%M:%S')}
        data['last'] = time.strftime('%Y-%m-%d %H:%M:%S')
        try:
            tmp = p + '.tmp'
            with open(tmp, 'w', encoding='utf-8') as f:
                json.dump(data, f, ensure_ascii=False, indent=1)
            os.replace(tmp, p)
        except OSError:
            pass


def read_autofix_log():
    """读自愈记录（selfcheck / 设置页用）。返回 {module: {...}} 或 {}。"""
    try:
        with open(_log_path(), encoding='utf-8') as f:
            return (json.load(f) or {}).get('modules') or {}
    except (OSError, ValueError):
        return {}


def missing_tools(cfg, module):
    """模块当前缺哪些**必需**工具键（可选工具不计入，见 OPTIONAL_TOOLS）。"""
    keys = MODULE_TOOLS.get(module, [])
    st = cfg.tool_status()
    return [k for k in keys if not st.get(k)]


def missing_optional(cfg, module):
    """缺哪些**可选**工具键，返回 {键: 影响说明}（只提示，不阻断）。"""
    opts = OPTIONAL_TOOLS.get(module, {})
    st = cfg.tool_status()
    return {k: why for k, why in opts.items() if not st.get(k)}


def ensure_tools(cfg, module, auto_install=True, log=None, timeout_s=600):
    """运行前自愈：缺的工具尽量用插件源补齐。

    返回 {'missing': 原始缺失键, 'installed': 装成功的插件,
          'still_missing': 装完仍缺的键, 'hints': {键: 指引}}。

    auto_install=False 时只检测不安装（供「仅检查」入口与单测用）。
    安装失败不抛异常——把结果如实返回，由调用方决定是否中止。
    """
    miss = missing_tools(cfg, module)
    result = {'missing': miss, 'installed': [], 'still_missing': [],
              'hints': {}, 'optional_missing': missing_optional(cfg, module)}
    if not miss:
        return result

    # 有指引的键（插件源里没有）先记下，避免白跑一趟下载
    hints = {k: _NO_PLUGIN_HINT[k] for k in miss if k in _NO_PLUGIN_HINT}

    if auto_install:
        try:
            items, err = pm.available_from_index(cfg, force=True)
            if err:
                if log:
                    log(f'插件源不可用（{err}），跳过自动安装')
            else:
                # 键 → 提供它的插件（一个插件可能提供多个键，如 blast 4 个）
                need_plugins = []
                for it in items:
                    prov = it.get('provides') or {}
                    if any(k in miss for k in prov):
                        need_plugins.append(it)
                for p in need_plugins:
                    if p.get('installed') and not p.get('needed'):
                        continue
                    try:
                        if log:
                            log(f'自动安装插件 {p["name"]}'
                                f'（{p.get("size", 0)/1048576:.1f}MB）…')
                        pm.install_urls(p.get('urls') or [p.get('url')],
                                        expected_sha256=p.get('sha256'))
                        result['installed'].append(p['name'])
                    except Exception as e:      # 单包失败不拦其他包
                        if log:
                            log(f'插件 {p["name"]} 安装失败: {e}')
        except Exception as e:
            if log:
                log(f'自动安装流程异常: {e}')

    # 装完重扫（refresh_plugins 会更新 Cfg.tools 的填缺层）
    if result['installed']:
        try:
            cfg.refresh_plugins()
        except Exception:
            pass

    still = missing_tools(cfg, module)
    result['still_missing'] = still
    result['hints'] = {k: hints.get(k) or
                       '插件源未提供该工具，请手动安装或到设置页指定路径'
                       for k in still}
    record_autofix(module, result['installed'], still)
    return result


def format_error(module, result):
    """自愈失败时的可操作错误文案（直接给用户看）。"""
    lines = [f'模块「{module}」缺少必要外部工具：'
             + '、'.join(result['still_missing'])]
    for k, hint in result['hints'].items():
        lines.append(f'  · {k}：{hint}')
    if result['installed']:
        lines.append('  已自动安装：' + '、'.join(result['installed']))
    lines.append('  可在「设置 → 插件管理 → 插件源」查看与手动安装。')
    return '\n'.join(lines)
