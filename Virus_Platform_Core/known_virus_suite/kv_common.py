#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
kv_common.py — known_virus_suite 公共层
=====================================
参考加载、日志、I/O、外部工具探测。

设计来源: virome_analysis_pipeline/batch_virus_depth.py 的
_load_reference_lengths / _load_ref_info_smart / setup_logging 等
（复制后改造，原文件不动）
"""

import gzip
import logging
import os
import shutil
import sys
import time
from pathlib import Path

# ── 平台适配 ──────────────────────────────────────────────
IS_WINDOWS = os.name == 'nt'


def smart_open(path, mode='rt'):
    """透明处理 .gz"""
    path = str(path)
    return gzip.open(path, mode, encoding=None if 'b' in mode else 'utf-8') if path.endswith('.gz') else open(path, mode, encoding=None if 'b' in mode else 'utf-8')


def which_tool(name, extra_dirs=None):
    """
    工具探测，兼顾 Windows 上的 .exe 和平台自带 tools 目录。
    返回绝对路径，找不到返回 None。
    """
    p = shutil.which(name)
    if p:
        return p
    # Windows: shutil.which 对没有 .exe 后缀的名字也能找到，但为稳妥再试一次
    if IS_WINDOWS:
        for suffix in ('.exe', '.bat', '.cmd'):
            p = shutil.which(name + suffix)
            if p:
                return p
    for d in (extra_dirs or []):
        d = Path(d)
        for cand in (d / name, d / f"{name}.exe"):
            if cand.is_file():
                return str(cand)
    return None


def platform_root() -> Path:
    """平台根目录（本引擎位于 <PLATFORM_ROOT>/Virus_Platform_Core/known_virus_suite/）。

    冻结分发（PyInstaller exe）：平台根 = **exe 所在目录**，与
    `config._detect_platform_root()` 同口径。这一支必须先判 —— 打包后本模块的
    `__file__` 落在 `<exe目录>/_internal/Virus_Platform_Core/...` 里，上溯既
    找不到 `Virus_Platform_Core/` 也找不到 `app.py`，最终兜底返回 `_internal`，
    于是 `ToolRegistry.probe()` 去找 `_internal/3rd/tools`（真实位置是
    `<exe目录>/3rd/tools`）→ 8 个外部工具全部探测失败、`tools.require()` 直接
    拒绝运行。注意：开发机上 dist/ 嵌在源码树内时，上溯会"碰巧"命中源码根而
    掩盖此问题，只有在用户机器（dist 单独拷走）才暴露 —— 实测于 2026-09-11。

    源码运行：从本文件所在目录逐级上溯，取同时含 Virus_Platform_Core/ 与 app.py
    的那一层。这样引擎换目录（包内子目录、或退回旧的顶层布局）都无需再改这里。
    """
    if getattr(sys, 'frozen', False):
        return Path(sys.executable).resolve().parent
    p = Path(__file__).resolve().parent
    for _ in range(6):
        if (p / 'Virus_Platform_Core').is_dir() and (p / 'app.py').is_file():
            return p
        if p.parent == p:
            break
        p = p.parent
    return Path(__file__).resolve().parent.parent.parent


class ToolRegistry:
    """
    记录本模块需要的所有外部工具，支持缺件降级判断。

    用法:
        reg = ToolRegistry(logger)
        reg.probe('minimap2', extra_dirs=[...])
        if reg.has('minimap2'): ...
    """

    NEEDED = {
        'salmon': 'salmon 定量引擎（identify；--writeBam 出映射位点）',
        'minibwa': 'minibwa 真比对（共识段内部固定使用）',
        'mafft': 'mafft 多序列比对（共识迭代）',
        'samtools': 'samtools（SAM/BAM 处理与覆盖度统计）',
        'pandepth': 'pandepth（覆盖度/深度统计，口径金标准）',
        'viral_consensus': 'viral_consensus（共识序列构建，对齐原管线）',
        'bcftools': 'bcftools（变异检出 caller，替代 freebayes/lofreq/ivar）',
    }

    # 工具可执行文件名与所在子目录（Windows 本机 tools\ 布局）
    TOOL_HINTS = {
        'salmon':   ('salmon',   r'salmon2'),
        'minibwa':  ('minibwa',  r'bin'),
        'mafft':    ('mafft',    r'mafft-win\ms\bin'),
        'samtools': ('samtools', r'samtools\bin'),
        'pandepth': ('pandepth', r'pandepth'),
        'viral_consensus': ('viral_consensus', r'viral_consensus'),
        'bcftools': ('bcftools', r'bcftools\bin'),
    }

    def __init__(self, logger=None):
        self.logger = logger or logging.getLogger('kv')
        self.paths = {}
        self.env = None

    def probe(self, extra_dirs=None):
        # 外部工具搜索目录：新布局 3rd/tools、3rd/bin 优先，旧根目录布局回退
        # （2026-09-10 目录整理把 tools/ bin/ 收进 3rd/）
        root = platform_root()
        bases = []
        for rel in ('3rd/tools', '3rd/bin', 'tools', 'bin'):
            d = root / rel
            if d.is_dir():
                bases.append(d)
        auto = []
        for base in bases:
            auto.append(str(base))
            for _exe, sub in self.TOOL_HINTS.values():
                d = base / sub
                if d.exists():
                    auto.append(str(d))
        dirs = list(dict.fromkeys(auto + [str(Path(d)) for d in (extra_dirs or [])]))
        for name in self.NEEDED:
            exe = self.TOOL_HINTS.get(name, (name, ''))[0]
            self.paths[name] = which_tool(exe, dirs)
        self.env = self._build_env(root)
        return dict(self.paths)

    def _build_env(self, root):
        """构造子进程环境：bcftools 插件目录 + 工具目录入 PATH。

        bcftools 的 +插件 与某些子命令需要 BCFTOOLS_PLUGINS 才能找到
        libexec/bcftools 下的插件；Windows 下不设会直接报插件未找到。
        PATH 里追加工具目录，让 samtools/bcftools 能互相调用。
        """
        env = dict(os.environ)
        for rel in ('3rd/tools', 'tools'):
            plugins = Path(root) / rel / 'bcftools' / 'libexec' / 'bcftools'
            if plugins.is_dir():
                env['BCFTOOLS_PLUGINS'] = str(plugins)
                break
        extra = [str(Path(p).parent) for p in self.paths.values() if p]
        if extra:
            env['PATH'] = os.pathsep.join(
                list(dict.fromkeys(extra)) + [env.get('PATH', '')])
        return env

    def has(self, name):
        return bool(self.paths.get(name))

    def require(self, *names):
        """缺件直接报错，附安装提示"""
        missing = [n for n in names if not self.has(n)]
        if missing:
            lines = [f"缺少必要工具: {', '.join(missing)}"]
            for m in missing:
                lines.append(f"  - {m}  ({self.NEEDED.get(m, '')})")
            raise RuntimeError('\n'.join(lines))
        return True

    def report(self):
        rows = []
        for name, desc in self.NEEDED.items():
            p = self.paths.get(name)
            rows.append(f"  [{'OK' if p else 'MISS':>4}] {name:<12} {desc}")
            if p:
                rows[-1] += f"\n         -> {p}"
        return '\n'.join(rows)


# ── 日志 ──────────────────────────────────────────────────
def setup_logger(out_dir, verbose=False, tag='kv'):
    """文件 + 控制台双通道日志"""
    out_dir = Path(out_dir)
    logs_dir = out_dir / 'logs'
    logs_dir.mkdir(parents=True, exist_ok=True)

    logger = logging.getLogger(tag)
    logger.setLevel(logging.DEBUG if verbose else logging.INFO)
    logger.handlers = []
    logger.propagate = False

    fh = logging.FileHandler(logs_dir / f"{tag}_{time.strftime('%Y%m%d_%H%M%S')}.log",
                             mode='a', encoding='utf-8')
    fh.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
    logger.addHandler(fh)

    ch = logging.StreamHandler(sys.stdout)
    ch.setFormatter(logging.Formatter('[%(asctime)s] %(levelname)s - %(message)s', '%H:%M:%S'))
    logger.addHandler(ch)

    return logger, logs_dir


def fmt_time(seconds):
    m, s = divmod(seconds, 60)
    h, m = divmod(m, 60)
    if h > 0:
        return f"{int(h)}h {int(m)}m {int(s)}s"
    if m > 0:
        return f"{int(m)}m {int(s)}s"
    return f"{s:.1f}s"


# ── 参考库加载 ────────────────────────────────────────────
def load_ref_lengths(ref_fasta, logger=None):
    """读参考 FASTA，返回 {accession: length}"""
    lengths = {}
    cur_id, cur_len = None, 0
    with smart_open(ref_fasta, 'rt') as f:
        for line in f:
            if line.startswith('>'):
                if cur_id:
                    lengths[cur_id] = cur_len
                cur_id = line.strip().split()[0][1:]
                cur_len = 0
            else:
                cur_len += len(line.strip())
        if cur_id:
            lengths[cur_id] = cur_len
    if logger:
        logger.info(f"参考库载入: {len(lengths)} 条序列, 合计 {sum(lengths.values())/1e6:.2f} Mb")
    return lengths


def load_ref_sequences(ref_fasta, only=None):
    """
    读参考 FASTA 序列本体。only 为 accession 集合时只保留这些。
    返回 {accession: sequence}
    """
    seqs = {}
    cur_id, buf = None, []
    with smart_open(ref_fasta, 'rt') as f:
        for line in f:
            if line.startswith('>'):
                if cur_id and (only is None or cur_id in only):
                    seqs[cur_id] = ''.join(buf).upper()
                cur_id = line.strip().split()[0][1:]
                buf = []
            else:
                buf.append(line.strip())
        if cur_id and (only is None or cur_id in only):
            seqs[cur_id] = ''.join(buf).upper()
    return seqs


# ref_info 列名同义词表（对齐原管线 _load_ref_info_smart）
ACC_SYNONYMS = ['Accession', 'accession', 'Virus GENBANK accession', 'ID']
TAX_SYNONYMS = ['Taxid', 'taxonomy_id', 'taxid']
SP_SYNONYMS = ['Species_NCBI', 'Species_ICTV', 'taxonomy', 'Species', 'description', 'Virus name(s)']
SEG_SYNONYMS = ['Segment', 'segment']


def read_text_lossless(path):
    """读文本文件，优先 UTF-8，其次 GBK，都失败才用 replace 兜底。

    同 kv_variant.read_text_lossless 口径（Windows 中文环境下外部工具/人工
    整理的表可能是系统 ANSI 代码页 936/GBK；只按 UTF-8 读 + replace 会把
    GBK 字节不可逆地替换成 U+FFFD，物种名「柑橘」变一串替换符）。
    """
    raw = Path(path).read_bytes()
    for enc in ('utf-8-sig', 'utf-8', 'gbk'):
        try:
            return raw.decode(enc)
        except UnicodeDecodeError:
            continue
    return raw.decode('utf-8', errors='replace')


def load_ref_info(ref_info_tsv, logger=None):
    """
    读参考注释表，返回
      {accession: {'taxid':..,'species':..,'segment':..,'molecule':..}}
    同时建立去版本号的别名索引。
    """
    info = {}
    if not ref_info_tsv or not os.path.exists(ref_info_tsv):
        if logger:
            logger.warning("未提供 ref_info，物种名将回退为 accession")
        return info

    idx = {}
    header = None
    # 编码兜底：GBK 表按 UTF-8+replace 读会把中文列变 U+FFFD（不可逆），
    # 这里用无损解码（utf-8 → gbk 回退，同 kv_variant.read_text_lossless 口径）
    for line in read_text_lossless(ref_info_tsv).splitlines():
        if line.startswith('####') or not line.strip():
            continue
        parts = line.rstrip('\n').split('\t')
        if header is None:
            header = [h.strip() for h in parts]

            def find(syns):
                return next((header.index(s) for s in syns if s in header), -1)

            idx = {
                'acc': find(ACC_SYNONYMS), 'tax': find(TAX_SYNONYMS),
                'sp': find(SP_SYNONYMS), 'seg': find(SEG_SYNONYMS),
                'mol': find(['Molecule_type2', 'Molecule_Type2', 'Molecule_type']),
                'host': find(['Host', 'host']),
                'geo': find(['Geo_Location', 'Geo_Location/Year', 'geo']),
                'gbtitle': find(['GenBank_Title', 'GenBank Title']),
                'mol1': find(['Molecule_type', 'Molecule_Type']),
            }
            if idx['acc'] == -1:
                if logger:
                    logger.error(f"ref_info 找不到 Accession 列，实际列: {header[:8]}")
                return {}
            continue

        if len(parts) <= idx['acc']:
            continue
        acc = parts[idx['acc']].strip()
        if not acc:
            continue

        def get(key):
            i = idx.get(key, -1)
            return parts[i].strip() if i != -1 and len(parts) > i else ''

        rec = {
            'taxid': get('tax') or 'Unannotated',
            'species': get('sp') or acc,
            'segment': get('seg'),
            'molecule': get('mol'),
            # 宿主/地理/GenBank 标题：识别表一并提供（此前只导出 4 列，
            # 用户反馈「提取信息不全」——库里 26 列注释大多被丢掉）
            'host': get('host'),
            'geo': get('geo'),
            'genbank_title': get('gbtitle'),
            'molecule_type': get('mol1'),
        }
        info[acc] = rec
        # 去版本号别名：同一 accession 多版本时后写入者胜出，
        # 会让物种名静默错配。检测到冲突就报警。
        bare = acc.split('.')[0]
        if bare in info and info[bare].get('species') != rec.get('species'):
            if logger:
                logger.warning(f"ref_info 多版本冲突: {bare} "
                               f"({info[bare].get('species')} vs "
                               f"{rec.get('species')})，使用后者")
        info[bare] = rec  # 去版本号别名

    if logger:
        logger.info(f"参考注释载入: {len(info)//2} 条（含去版本号别名）")
    return info


def ref_lookup(ref_info, accession):
    """按 accession 查注释，容错版本号差异"""
    if not ref_info:
        return {'taxid': 'Unannotated', 'species': accession, 'segment': '', 'molecule': ''}
    rec = ref_info.get(accession) or ref_info.get(accession.split('.')[0])
    return rec or {'taxid': 'Unannotated', 'species': accession, 'segment': '', 'molecule': ''}
