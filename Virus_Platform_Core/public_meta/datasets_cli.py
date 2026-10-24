#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
📥 NCBI datasets CLI 统一封装 v1.1
==================================
官方 datasets 命令行（download / summary 域）的平台化包装，四种数据类型全接：

  genome    按 taxon / accession 下基因组包（--chromosomes 染色体过滤）
  gene      按 gene-id / symbol+taxon / protein accession 下基因包
  taxonomy  按 taxon 下分类报告包
  virus     按 taxon(+host 宿主过滤) 下病毒基因组 / 蛋白包（accession 单条精确）

先查后下：query_records() 走 `datasets summary …` 预览匹配记录
（归一化成前端表格行），确认后再用 download 系列拉包——基因组查询下载
页的「🔍 查询」按钮即走这条链路。

平台化差异（与裸用 datasets 的差异）：
  - 工具定位：PATH → 3rd/tools/datasets/datasets.exe；缺失时
    `setup()` 可从 NCBI 官方（ftp.ncbi.nlm.nih.gov 白名单内）拉取 win64 版
  - 旧语法翻译：v14+ 已移除 --exclude-gff3 等旗标，统一翻译为 --include 组合
  - 输出落位：zip + 解包收进 <out>/，FASTA 合并去重复用 host_genome 的管线
  - virus 模式可直接产出平台鉴定库数据集布局（reference.fasta +
    reference.ref_info.tsv），建库页「原始数据集」卡自动识别，一键建索引
  - data_report.jsonl → reference.ref_info.tsv（Accession/Taxid/Species/
    Host/Geo_Location…列名对齐 kv_common.load_ref_info 同义词表）

用法（CLI 见 main.py datasets 子命令；本模块也可编程调用）:
  python main.py datasets genome --taxon "bos taurus"
  python main.py datasets genome --accession GCF_000001405.40 --chromosomes X,Y
  python main.py datasets gene --symbol brca1 --taxon "mus musculus"
  python main.py datasets taxonomy --taxon human,sars-cov-2
  python main.py datasets virus --taxon "tomato brown rugose fruit virus" --name tobrfv
"""

import os
import sys
import re
import json
import shutil
import zipfile
import argparse
import subprocess
from datetime import datetime

if sys.platform == 'win32':
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')

if __package__ in (None, ''):
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(
        os.path.abspath(__file__)))))

from Virus_Platform_Core.config import DIRS, PLATFORM_ROOT
from Virus_Platform_Core.utils import safe_open

# 官方分发（NCBI 域名白名单内；win64 直发 exe，无需解压）
DATASETS_EXE_URL = ('https://ftp.ncbi.nlm.nih.gov/pub/datasets/'
                    'command-line/LATEST/win64/datasets.exe')
_DATASETS_SUBDIR = os.path.join('3rd', 'tools', 'datasets')

# 引用注释表列名对齐 known_virus_suite.kv_common.load_ref_info 同义词表
_REF_INFO_HEADER = ('Accession', 'Taxid', 'Species', 'Host', 'Host_Taxid',
                    'Segment', 'Molecule_type', 'Geo_Location',
                    'Collection_Date', 'Isolate', 'Length', 'Completeness',
                    'GenBank_Title')


def _emitters(log=None):
    """统一输出出口（与 host_genome._emitters 同约定）。"""
    if log is None:
        class UI:
            CYAN = '\033[96m'; GREEN = '\033[92m'; YELLOW = '\033[93m'
            RED = '\033[91m'; PURPLE = '\033[95m'
            BOLD = '\033[1m'; RESET = '\033[0m'

            @staticmethod
            def ok(msg):    print(f"  {UI.GREEN}✓{UI.RESET} {msg}")
            @staticmethod
            def warn(msg):  print(f"  {UI.YELLOW}⚠{UI.RESET} {msg}")
            @staticmethod
            def err(msg):   print(f"  {UI.RED}✗{UI.RESET} {msg}")
            @staticmethod
            def info(msg):  print(f"  {UI.CYAN}→{UI.RESET} {msg}")

            @staticmethod
            def header(msg):
                print(f"\n{UI.PURPLE}{UI.BOLD}{'='*55}{UI.RESET}")
                print(f"{UI.PURPLE}{UI.BOLD} {msg}{UI.RESET}")
                print(f"{UI.PURPLE}{UI.BOLD}{'='*55}{UI.RESET}")
        return UI

    class _L:
        @staticmethod
        def ok(msg):    log('✓ ' + str(msg))
        @staticmethod
        def warn(msg):  log('⚠ ' + str(msg))
        @staticmethod
        def err(msg):   log('✗ ' + str(msg))
        @staticmethod
        def info(msg):  log('→ ' + str(msg))

        @staticmethod
        def header(msg):
            log('')
            log('════════ ' + str(msg) + ' ════════')

    return _L()


# ==========================================
# CLI 定位与安装
# ==========================================
def datasets_exe() -> str:
    """定位 datasets 可执行文件：PATH → 平台内置 3rd/tools/datasets/。"""
    import shutil as _sh
    found = _sh.which('datasets')
    if found:
        return found
    cand = os.path.join(PLATFORM_ROOT, _DATASETS_SUBDIR, 'datasets.exe')
    if os.path.isfile(cand):
        return cand
    return ''


def datasets_available() -> bool:
    return bool(datasets_exe())


def setup(dest_dir: str = '', log=None) -> str:
    """下载官方 datasets.exe（NCBI 官方域名，白名单内）到平台工具目录。

    返回可执行文件路径。已存在时直接返回（先删后装可强制更新）。"""
    from .host_genome import _http  # 复用 NCBI 域名白名单 HTTP 通道
    E = _emitters(log)
    if not dest_dir:
        dest_dir = os.path.join(PLATFORM_ROOT, _DATASETS_SUBDIR)
    os.makedirs(dest_dir, exist_ok=True)
    exe = os.path.join(dest_dir, 'datasets.exe')
    if os.path.isfile(exe) and os.path.getsize(exe) > 1_000_000:
        E.ok(f"datasets CLI 已就绪: {exe}")
        return exe
    E.info(f"下载官方 datasets CLI（win64）: {DATASETS_EXE_URL}")
    http = _http()                          # 默认白名单 = NCBI 域名
    http.download(DATASETS_EXE_URL, exe)
    E.ok(f"完成: {exe} ({os.path.getsize(exe) / (1024*1024):.1f} MB)")
    r = subprocess.run([exe, '--version'], capture_output=True, text=True,
                       timeout=30)
    ver = (r.stdout or r.stderr).strip().splitlines()[-1] if (
        r.stdout or r.stderr).strip() else '?'
    E.ok(f"自检: {ver}")
    return exe


# ==========================================
# 底层执行
# ==========================================
def _run_and_stream(cmd, log=None, timeout=7200):
    """执行子进程，stdout/stderr 逐行转发到 log（CLI 模式则直接透传）。"""
    if log is None:
        return subprocess.run(cmd, timeout=timeout).returncode
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, text=True,
                            encoding='utf-8', errors='replace')
    try:
        for line in proc.stdout or []:
            line = line.rstrip()
            if line:
                log('  ' + line)
        return proc.wait(timeout=timeout)
    except Exception:
        proc.kill()
        raise


def _download_package(args: list, out_dir: str, tag: str,
                      log=None) -> str:
    """执行 datasets download … --filename <out_dir>/<tag>.zip 并解包。

    返回解包后的 ncbi_dataset/data 目录。"""
    E = _emitters(log)
    exe = datasets_exe()
    if not exe:
        raise RuntimeError(
            '未找到 NCBI datasets CLI：先执行 python main.py datasets-setup '
            '（自动从 NCBI 官方下载到 3rd/tools/datasets/），'
            '或自行安装后确保 datasets 在 PATH')
    os.makedirs(out_dir, exist_ok=True)
    zip_path = os.path.join(out_dir, f'{tag}.zip')
    cmd = [exe] + args + ['--filename', zip_path]
    E.info(f"执行: {' '.join(args) or '(空)'}")
    rc = _run_and_stream(cmd, log=log)
    if rc != 0 or not os.path.isfile(zip_path):
        raise RuntimeError(f'datasets 下载失败（exit={rc}）: '
                           f"{' '.join(cmd)}")
    E.ok(f"数据包: {os.path.basename(zip_path)} "
         f"({os.path.getsize(zip_path) / (1024*1024):.1f} MB)")
    data_dir = os.path.join(out_dir, 'ncbi_dataset', 'data')
    if os.path.isdir(data_dir) and os.listdir(data_dir):
        E.ok("已解包，跳过")
    else:
        with zipfile.ZipFile(zip_path) as zf:
            zf.extractall(out_dir)
        E.ok(f"解包 → {os.path.relpath(data_dir, out_dir)}")
    return data_dir


def count_fasta_seqs(fasta_paths) -> int:
    """数 FASTA 总序列条数（'>' 头行计数）。"""
    n = 0
    for p in fasta_paths:
        with open(p, encoding='utf-8', errors='replace') as f:
            for line in f:
                if line.startswith('>'):
                    n += 1
    return n


def collect_reports(data_dir: str):
    """收集数据包内全部 data report（jsonl，每行一条记录）。"""
    reports = []
    if not os.path.isdir(data_dir):
        return reports
    for fn in sorted(os.listdir(data_dir)):
        if fn.endswith('.jsonl') and 'report' in fn:
            with open(os.path.join(data_dir, fn), encoding='utf-8') as f:
                for line in f:
                    line = line.strip()
                    if line:
                        try:
                            reports.append(json.loads(line))
                        except ValueError:
                            continue
    return reports


def collect_fasta(data_dir: str):
    """收集数据包内全部 FASTA（.fna/.fa/.fasta/.faa，含子目录）。"""
    out = []
    for root, _, files in os.walk(data_dir):
        for fn in files:
            if fn.lower().endswith(('.fna', '.fa', '.fasta', '.faa')):
                out.append(os.path.join(root, fn))
    return sorted(out)


def report_taxid(reports):
    """从 report 记录猜 taxid（assembly/virus/gene 报告字段名不一）。"""
    for r in reports:
        for path in (('organism', 'taxId'), ('organism', 'id'),
                     ('taxId',), ('tax_id',),
                     ('virus', 'taxId')):
            cur = r
            for k in path:
                if isinstance(cur, dict):
                    cur = cur.get(k)
                else:
                    cur = None
                    break
            if isinstance(cur, int):
                return cur
    return None


# ==========================================
# 旧语法翻译（v14+ 移除了 --exclude-*）
# ==========================================
def build_genome_args(taxon=None, accession=None, chromosomes=None,
                      include=None, exclude=()):
    """构造 `download genome …` 参数；exclude 旗标翻译为 include 组合。

    exclude 元素形如 'gff3'/'rna'/'protein'；include 与 exclude 同给时
    include 优先（先按 include 全集再剔除 exclude）。"""
    if not taxon and not accession:
        raise ValueError('genome 模式需要 --taxon 或 --accession')
    args = ['download', 'genome',
            'taxon' if taxon else 'accession', taxon or accession]
    if chromosomes:
        args += ['--chromosomes', chromosomes]
    if include:
        inc = [x.strip() for x in re.split(r'[,\s]+', include) if x.strip()]
    else:
        inc = ['genome', 'gff3', 'rna', 'protein']
    inc = [x for x in inc if x not in set(exclude)]
    if not inc:
        raise ValueError('--include 剔除 exclude 后为空')
    args += ['--include', ','.join(inc)]
    return args


def build_gene_args(gene_id=None, symbol=None, taxon=None, accession=None):
    """构造 `download gene …` 参数（三种定位方式三选一）。"""
    if gene_id:
        return ['download', 'gene', 'gene-id', str(gene_id)]
    if accession:
        return ['download', 'gene', 'accession', accession]
    if symbol and taxon:
        return ['download', 'gene', 'symbol', symbol, '--taxon', taxon]
    raise ValueError('gene 模式需要 --gene-id，或 --symbol + --taxon，'
                     '或 --accession（蛋白质 accession）')


def build_taxonomy_args(taxons):
    if not taxons:
        raise ValueError('taxonomy 模式需要 --taxon（逗号分隔多个）')
    return ['download', 'taxonomy', 'taxon',
            ','.join(t.strip() for t in taxons.split(',') if t.strip())]


def build_virus_args(taxon=None, protein=None, host=None, include=('genome',),
                     accession=None):
    """构造 `download virus …` 参数。

    taxon + include(genome/protein/annotation/biosample/all/none) →
      `virus genome taxon <t> [--host h] --include …`；
    accession（如 NC_045512.2，查询结果里挑单条时用）→
      `virus genome accession <a>`；
    protein（如 S / spike）→ `virus protein <p> [--host h]`。"""
    inc = [x.strip() for x in (include or []) if x.strip()] or ['genome']
    if accession and not taxon and not protein:
        args = ['download', 'virus', 'genome', 'accession', accession]
        if host:
            args += ['--host', host]
        args += ['--include', ','.join(inc)]
        return args
    if protein and not taxon:
        args = ['download', 'virus', 'protein', protein]
        if host:
            args += ['--host', host]
        return args
    if not taxon:
        raise ValueError('virus 模式需要 --taxon，或 --protein（+可选 --host），'
                         '或 --accession')
    args = ['download', 'virus', 'genome', 'taxon', taxon]
    if host:
        args += ['--host', host]
    args += ['--include', ','.join(inc)]
    return args


# ==========================================
# data_report → reference.ref_info.tsv
# ==========================================
def report_to_ref_info(reports, out_tsv: str, log=None) -> int:
    """病毒数据报告 → 平台鉴定库引用注释表（列名对齐 load_ref_info）。

    返回写入行数。"""
    E = _emitters(log)
    n = 0
    with safe_open(out_tsv, 'wt') as f:
        f.write('\t'.join(_REF_INFO_HEADER) + '\n')
        for r in reports:
            acc = str(r.get('accession') or '').strip()
            if not acc:
                continue
            virus = r.get('virus') or {}
            host = r.get('host') or {}
            isolate = r.get('isolate') or {}
            loc = r.get('location') or {}
            row = [
                acc,
                str(virus.get('taxId') or ''),
                str(virus.get('organismName') or ''),
                str(host.get('organismName') or ''),
                str(host.get('taxId') or ''),
                '',                                   # Segment（报告未给）
                '',                                   # Molecule_type
                str(loc.get('geographicLocation') or ''),
                str(isolate.get('collectionDate') or ''),
                str(isolate.get('name') or ''),
                str(r.get('length') or ''),
                str(r.get('completeness') or ''),
                str(r.get('sourceDatabase') or ''),
            ]
            f.write('\t'.join(row) + '\n')
            n += 1
    E.ok(f"引用注释表: {out_tsv}（{n} 条）")
    return n


# ==========================================
# 高层入口
# ==========================================
def _slug(s: str) -> str:
    return re.sub(r'[^\w\-]+', '_', str(s)).strip('_') or 'ds'


def _write_manifest(out_dir, mode, args, reports, fasta_files):
    meta = {
        'source': 'ncbi-datasets', 'mode': mode,
        'args': args, 'downloaded_at': datetime.now().isoformat(),
        'n_reports': len(reports), 'fasta_files': fasta_files,
    }
    with safe_open(os.path.join(out_dir, 'datasets_meta.json'), 'wt') as f:
        json.dump(meta, f, ensure_ascii=False, indent=2)


def download_genome(taxon=None, accession=None, chromosomes=None,
                    include=None, exclude=(), out_dir=None, log=None) -> dict:
    """基因组包下载（taxon / accession 双入口）。

    返回 dict：out_dir / data_dir / fasta / reports / taxid / n_seqs。"""
    E = _emitters(log)
    E.header(f"datasets genome: {accession or taxon}")
    if not out_dir:
        out_dir = os.path.join('datasets_dl', 'genome_' + _slug(
            accession or taxon))
    os.makedirs(out_dir, exist_ok=True)
    args = build_genome_args(taxon=taxon, accession=accession,
                             chromosomes=chromosomes, include=include,
                             exclude=exclude)
    data_dir = _download_package(args, out_dir, 'genome', log=log)
    reports = collect_reports(data_dir)
    fasta = collect_fasta(data_dir)
    n_seqs = count_fasta_seqs(fasta)
    for p in fasta:
        E.ok(f"FASTA: {os.path.relpath(p, out_dir)}")
    E.ok(f"记录 {len(reports)} 条 / FASTA {len(fasta)} 个 / 序列 {n_seqs} 条")
    taxid = report_taxid(reports)
    if taxid:
        E.ok(f"TaxID: {taxid}")
    _write_manifest(out_dir, 'genome', args, reports, fasta)
    return {'out_dir': out_dir, 'data_dir': data_dir, 'fasta': fasta,
            'reports': reports, 'taxid': taxid, 'n_seqs': n_seqs,
            'mode': 'genome'}


def download_gene(gene_id=None, symbol=None, taxon=None, accession=None,
                  out_dir=None, log=None) -> dict:
    """基因包下载（gene-id / symbol+taxon / protein accession）。

    默认包内含 rna.fna + protein.faa + data_report.jsonl。"""
    E = _emitters(log)
    label = (f"gene-id {gene_id}" if gene_id else
             f"accession {accession}" if accession else
             f"symbol {symbol} / {taxon}")
    E.header(f"datasets gene: {label}")
    if not out_dir:
        out_dir = os.path.join('datasets_dl', 'gene_' + _slug(
            gene_id or accession or symbol))
    os.makedirs(out_dir, exist_ok=True)
    args = build_gene_args(gene_id=gene_id, symbol=symbol, taxon=taxon,
                           accession=accession)
    data_dir = _download_package(args, out_dir, 'gene', log=log)
    fasta = collect_fasta(data_dir)
    reports = collect_reports(data_dir)
    for p in fasta:
        E.ok(f"FASTA: {os.path.relpath(p, out_dir)}")
    _write_manifest(out_dir, 'gene', args, reports, fasta)
    return {'out_dir': out_dir, 'data_dir': data_dir, 'fasta': fasta,
            'reports': reports, 'mode': 'gene'}


def download_taxonomy(taxons, out_dir=None, log=None) -> dict:
    """分类报告包下载（taxonomy_report.jsonl + taxonomy_summary.tsv）。"""
    E = _emitters(log)
    E.header(f"datasets taxonomy: {taxons}")
    if not out_dir:
        out_dir = os.path.join('datasets_dl', 'taxonomy_' + _slug(
            str(taxons).split(',')[0]))
    os.makedirs(out_dir, exist_ok=True)
    args = build_taxonomy_args(taxons)
    data_dir = _download_package(args, out_dir, 'taxonomy', log=log)
    reports = collect_reports(data_dir)
    E.ok(f"分类记录 {len(reports)} 条")
    _write_manifest(out_dir, 'taxonomy', args, reports, [])
    return {'out_dir': out_dir, 'data_dir': data_dir, 'reports': reports,
            'mode': 'taxonomy'}


def download_virus(taxon=None, protein=None, host=None,
                   include=('genome',), name=None, out_dir=None,
                   as_dataset=True, accession=None, log=None) -> dict:
    """病毒包下载（genome / protein；--host 宿主过滤；accession 单条精确）。

    as_dataset=True（默认）时按平台鉴定库数据集布局落位：
      databases/virusref_db/<name>/reference.fasta + reference.ref_info.tsv
    —— 建库页「原始数据集」卡自动识别，可直接一键建 kv 索引/kunpeng 病毒库。
    """
    E = _emitters(log)
    E.header(f"datasets virus: {taxon or protein or accession}"
             + (f" (host={host})" if host else ''))
    if as_dataset and not out_dir:
        if not name:
            name = _slug(taxon or protein or accession)
        out_dir = os.path.join(DIRS['virus_src'], name)
    if not out_dir:
        out_dir = os.path.join('datasets_dl',
                               'virus_' + _slug(taxon or protein or accession))
    os.makedirs(out_dir, exist_ok=True)
    args = build_virus_args(taxon=taxon, protein=protein, host=host,
                            include=include, accession=accession)
    data_dir = _download_package(args, out_dir, 'virus', log=log)
    reports = collect_reports(data_dir)
    fasta = collect_fasta(data_dir)
    n_seqs = 0
    for p in fasta:
        c = sum(1 for line in open(p, encoding='utf-8', errors='replace')
                if line.startswith('>'))
        n_seqs += c
        E.ok(f"FASTA: {os.path.relpath(p, out_dir)}（{c} 条序列）")
    if not n_seqs:
        raise RuntimeError('数据包内没有序列（换个 taxon/host/include 试试，'
                           '或加 --include protein 只取蛋白）')

    ref_info_path = ''
    reference_path = ''
    if reports and (taxon or accession):     # protein-only 包没有病毒报告
        ref_info_path = os.path.join(out_dir, 'reference.ref_info.tsv')
        report_to_ref_info(reports, ref_info_path, log=log)
    if n_seqs:
        # 平台数据集锚点名：reference.fasta（自包含库目录识别依据）
        merged = os.path.join(out_dir, 'reference.fasta')
        if len(fasta) == 1:
            shutil.copyfile(fasta[0], merged)
        else:
            with safe_open(merged, 'wb') as fo:
                for p in fasta:
                    with open(p, 'rb') as fi:
                        shutil.copyfileobj(fi, fo, 1 << 20)
        reference_path = merged
        E.ok(f"reference.fasta 就绪（可直接建库）: {reference_path}")

    _write_manifest(out_dir, 'virus', args, reports,
                    fasta + ([reference_path] if reference_path else []))
    return {'out_dir': out_dir, 'data_dir': data_dir, 'fasta': fasta,
            'reference': reference_path, 'ref_info': ref_info_path,
            'reports': reports, 'n_seqs': n_seqs, 'mode': 'virus',
            'name': name or os.path.basename(out_dir)}


def download(mode: str, log=None, **kw) -> dict:
    """统一分发：mode ∈ genome / gene / taxonomy / virus。"""
    table = {'genome': download_genome, 'gene': download_gene,
             'taxonomy': download_taxonomy, 'virus': download_virus}
    if mode not in table:
        raise ValueError(f'未知 datasets 模式: {mode}')
    return table[mode](log=log, **kw)


# ==========================================
# summary 查询（先查后下：预览匹配记录，再决定下什么）
# ==========================================
# 各 kind 默认预览条数（真实匹配总数用 limit+1 探测：返回了 limit+1 条
# 即说明被截断；datasets 的 total_count 会跟着 --limit 截断，探不出真总数）
_SUMMARY_DEFAULT_LIMIT = {'genome': 200, 'virus': 500, 'gene': 200,
                          'taxonomy': 50}
_SUMMARY_TIMEOUT = 300


def build_summary_args(kind, taxon=None, accession=None, host=None,
                       gene_id=None, symbol=None, protein_acc=None,
                       taxons=None, limit=None):
    """构造 `datasets summary …` 参数（与 build_*_args 同源的查询入口）。

    limit 传 None 时按 kind 用默认预览条数；taxonomy 记录数固定，不传 limit。"""
    args = ['summary']
    if kind == 'genome':
        if not taxon and not accession:
            raise ValueError('genome 查询需要 taxon 或 accession')
        args += ['genome', 'taxon' if taxon else 'accession',
                 taxon or accession]
    elif kind == 'virus':
        if not taxon and not accession:
            raise ValueError('virus 查询需要 taxon 或 accession')
        args += ['virus', 'genome', 'taxon' if taxon else 'accession',
                 taxon or accession]
        if host:
            args += ['--host', host]
    elif kind == 'gene':
        if gene_id:
            args += ['gene', 'gene-id', str(gene_id)]
        elif protein_acc:
            args += ['gene', 'accession', protein_acc]
        elif symbol and taxon:
            args += ['gene', 'symbol', symbol, '--taxon', taxon]
        else:
            raise ValueError('gene 查询需要 gene-id，或 symbol + taxon，'
                             '或蛋白 accession')
    elif kind == 'taxonomy':
        if not taxons:
            raise ValueError('taxonomy 查询需要 taxon（逗号分隔多个）')
        args += ['taxonomy', 'taxon', taxons]
    else:
        raise ValueError(f'未知查询类型: {kind}')
    if limit and kind in ('genome', 'virus'):
        args += ['--limit', str(int(limit))]
    return args


def _dig(d, *paths):
    """按多个候选点路径取第一个非空值（summary 输出 snake_case，
    download 的 data_report 是 camelCase，双风格都接）。"""
    for path in paths:
        cur = d
        for k in path.split('.'):
            if isinstance(cur, dict):
                cur = cur.get(k)
            else:
                cur = None
                break
        if cur not in (None, ''):
            return str(cur)
    return ''


def _norm_summary_row(kind, r):
    """单条 summary 记录 → 前端表格直渲染的扁平行（字段随 kind 而异）。"""
    if kind == 'genome':
        asm = r.get('assembly_info') or {}
        org = r.get('organism') or {}
        return {
            'accession': _dig(r, 'accession'),
            'name': _dig(asm, 'assembly_name'),
            'organism': _dig(org, 'organism_name', 'organismName'),
            'common': _dig(org, 'common_name', 'commonName'),
            'taxid': _dig(org, 'tax_id', 'taxId'),
            'length': _dig(r, 'assembly_stats.total_sequence_length',
                           'assemblyStats.totalSequenceLength'),
            'level': _dig(asm, 'assembly_level', 'assemblyLevel'),
            'source': _dig(r, 'source_database', 'sourceDatabase')
                      .replace('SOURCE_DATABASE_', ''),
            'released': _dig(asm, 'release_date', 'releaseDate')[:10],
            'submitter': _dig(asm, 'submitter'),
        }
    if kind == 'virus':
        return {
            'accession': _dig(r, 'accession'),
            'length': _dig(r, 'length'),
            'virus': _dig(r, 'virus.organism_name', 'virus.organismName'),
            'virus_taxid': _dig(r, 'virus.tax_id', 'virus.taxId'),
            'host': _dig(r, 'host.organism_name', 'host.organismName'),
            'host_taxid': _dig(r, 'host.tax_id', 'host.taxId'),
            'geo': _dig(r, 'location.geographic_location',
                        'location.geographicLocation'),
            'isolate': _dig(r, 'isolate.name'),
            'collected': _dig(r, 'isolate.collection_date',
                              'isolate.collectionDate')[:10],
            'completeness': _dig(r, 'completeness'),
            'source': _dig(r, 'source_database', 'sourceDatabase'),
            'released': _dig(r, 'release_date', 'releaseDate')[:10],
        }
    if kind == 'gene':
        g = r.get('gene') or {}
        chroms = g.get('chromosomes') or []
        return {
            'accession': _dig(g, 'gene_id'),
            'name': _dig(g, 'symbol'),
            'organism': _dig(g, 'taxname'),
            'common': _dig(g, 'common_name'),
            'taxid': _dig(g, 'tax_id'),
            'desc': _dig(g, 'description'),
            'chroms': ','.join(map(str, chroms[:4]))
                      + ('…' if len(chroms) > 4 else ''),
            'type': _dig(g, 'type'),
            'n_proteins': _dig(g, 'protein_count'),
            'n_transcripts': _dig(g, 'transcript_count'),
        }
    # taxonomy：嵌套在 taxonomy.* 下
    t = r.get('taxonomy') or r
    cls = t.get('classification') or {}
    rank = ''
    for rk in ('species', 'genus', 'family', 'order', 'class', 'phylum',
               'kingdom', 'realm', 'domain'):
        if isinstance(cls.get(rk), dict):
            rank = f"{rk}: {cls[rk].get('name')}"
            break
    return {
        'accession': _dig(t, 'tax_id', 'taxid'),
        'name': _dig(t, 'current_scientific_name.name',
                     'currentScientificName.name',
                     'scientific_name'),
        'rank': _dig(t, 'group_name'),
        'lineage': rank,
        'genetic_code': _dig(t, 'genetic_code.transl_table',
                             'geneticCode.translTable'),
    }


def query_records(kind, taxon=None, accession=None, host=None,
                  gene_id=None, symbol=None, protein_acc=None, taxons=None,
                  limit=None, log=None) -> dict:
    """datasets summary 预览查询（先查后下）。

    返回 {'kind', 'total_count', 'returned', 'truncated', 'rows', 'query'}。
    rows 已按 kind 归一化（_norm_summary_row），前端表格直接渲染；
    truncated=True 表示匹配数超出预览上限、只回传了前 limit 条。"""
    E = _emitters(log)
    exe = datasets_exe()
    if not exe:
        raise RuntimeError(
            '未找到 NCBI datasets CLI：先执行 python main.py datasets-setup '
            '（自动从 NCBI 官方下载到 3rd/tools/datasets/），'
            '或自行安装后确保 datasets 在 PATH')
    cap = int(limit) if limit else _SUMMARY_DEFAULT_LIMIT.get(kind, 200)
    args = build_summary_args(kind, taxon=taxon, accession=accession,
                              host=host, gene_id=gene_id, symbol=symbol,
                              protein_acc=protein_acc, taxons=taxons,
                              limit=cap + 1)
    E.info(f"执行: {' '.join(args)}")
    try:
        r = subprocess.run([exe] + args, capture_output=True, text=True,
                           encoding='utf-8', errors='replace',
                           timeout=_SUMMARY_TIMEOUT)
    except subprocess.TimeoutExpired:
        raise RuntimeError(f'datasets summary 查询超时（{_SUMMARY_TIMEOUT}s）：'
                           '换个更具体的 taxon 再试')
    if r.returncode != 0:
        tail = ((r.stderr or r.stdout or '').strip().splitlines() or [''])[-1]
        raise RuntimeError(f'datasets summary 查询失败（exit={r.returncode}）: '
                           f'{tail}')
    try:
        data = json.loads(r.stdout or '{}')
    except ValueError:
        raise RuntimeError('datasets summary 输出不是有效 JSON'
                           '（CLI 版本过旧？重跑 python main.py '
                           'datasets-setup 更新）')
    reports = data.get('reports') or []
    rows = [_norm_summary_row(kind, x) for x in reports[:cap]]
    # --limit 时 total_count 跟着截断：取了 cap+1 条、截到 cap 条，
    # total>len(rows) 即「还有更多没显示」
    total = int(data.get('total_count') or len(reports))
    return {'kind': kind, 'total_count': total, 'returned': len(rows),
            'truncated': total > len(rows), 'rows': rows,
            'query': ' '.join(args[1:])}


# ==========================================
# CLI
# ==========================================
def main():
    from .host_genome import UI
    parser = argparse.ArgumentParser(
        description='📥 NCBI datasets CLI 平台封装（genome/gene/taxonomy/virus）',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例（镜像官方 datasets download 语法）:
  python host_genome_datasets.py genome --taxon "bos taurus"
  python host_genome_datasets.py genome --accession GCF_000001405.40 ^
      --chromosomes X,Y --include genome
  python host_genome_datasets.py gene --gene-id 672
  python host_genome_datasets.py gene --symbol brca1 --taxon "mus musculus"
  python host_genome_datasets.py gene --accession NP_000483.3
  python host_genome_datasets.py taxonomy --taxon human,sars-cov-2
  python host_genome_datasets.py virus --taxon sars-cov-2 --host dog
  python host_genome_datasets.py virus --protein S --host dog --out SARS2-spike-dog
  python host_genome_datasets.py virus --taxon "tomato brown rugose fruit virus" --name tobrfv
        """)
    sub = parser.add_subparsers(dest='mode', required=True)

    def _out(p):
        p.add_argument('--out', default=None, help='输出目录（默认按模式命名）')

    g = sub.add_parser('genome', help='基因组包（taxon / accession）')
    g.add_argument('--taxon', default=None, help='物种名（科学名）')
    g.add_argument('--accession', default=None, help='组装 accession（GCF_/GCA_）')
    g.add_argument('--chromosomes', default=None,
                   help='染色体过滤（如 X,Y 或 NC_00xxxx.1）')
    g.add_argument('--include', default=None,
                   help='包含内容（genome,gff3,rna,protein；默认全量）')
    g.add_argument('--exclude', default='',
                   help='剔除内容（旧版 --exclude-* 翻译，如 gff3,rna）')
    _out(g)
    g.set_defaults(func=lambda a: download_genome(
        taxon=a.taxon, accession=a.accession, chromosomes=a.chromosomes,
        include=a.include, exclude=tuple(x for x in a.exclude.split(',')
                                         if x.strip()), out_dir=a.out))

    ge = sub.add_parser('gene', help='基因包（gene-id / symbol / accession）')
    ge.add_argument('--gene-id', default=None)
    ge.add_argument('--symbol', default=None)
    ge.add_argument('--taxon', default=None, help='symbol 定位时必给')
    ge.add_argument('--accession', default=None,
                    help='蛋白质 accession（如 NP_000483.3）')
    _out(ge)
    ge.set_defaults(func=lambda a: download_gene(
        gene_id=a.gene_id, symbol=a.symbol, taxon=a.taxon,
        accession=a.accession, out_dir=a.out))

    t = sub.add_parser('taxonomy', help='分类报告包')
    t.add_argument('--taxon', required=True, help='逗号分隔多个，如 human,sars-cov-2')
    _out(t)
    t.set_defaults(func=lambda a: download_taxonomy(a.taxon, out_dir=a.out))

    v = sub.add_parser('virus', help='病毒包（genome/protein + host 过滤）')
    v.add_argument('--taxon', default=None, help='病毒 taxon（如 sars-cov-2）')
    v.add_argument('--accession', default=None,
                   help='病毒核酸 accession（如 NC_045512.2，单条精确下载）')
    v.add_argument('--protein', default=None, help='蛋白名（如 S / spike）')
    v.add_argument('--host', default=None, help='宿主过滤（如 dog / 番茄拉丁名）')
    v.add_argument('--include', default='genome',
                   help='genome,protein,annotation,biosample,all,none')
    v.add_argument('--name', default=None,
                   help='数据集名（默认按 taxon 命名，落 virusref_db/<name>）')
    v.add_argument('--no-dataset', action='store_true',
                   help='不按平台数据集布局落位（存 datasets_dl/ 通用目录）')
    _out(v)
    v.set_defaults(func=lambda a: download_virus(
        taxon=a.taxon, protein=a.protein, host=a.host,
        include=tuple(x for x in a.include.split(',') if x.strip()),
        name=a.name, out_dir=a.out, as_dataset=not a.no_dataset,
        accession=a.accession))

    args = parser.parse_args()
    res = args.func(args)
    UI.header('下载完成')
    return res


if __name__ == '__main__':
    main()
