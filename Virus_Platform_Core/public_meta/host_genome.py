#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
📥 宿主参考基因组下载器 v1.1
============================
双源下载指定宿主物种的参考基因组（供建 kunpeng 宿主库）：
  - ncbi（默认）: NCBI datasets CLI / E-utilities + HTTPS 回退
  - gget: gget.ref → Ensembl / Ensembl Genomes（plants/fungi/metazoa…），
    覆盖 NCBI 之外的植物宿主（水稻/番茄/小麦等 Ensembl Plants 物种）

功能特性:
  - 自动解压与合并多文件
  - 序列名去重（处理多基因组合并时的重名问题）
  - 支持叶绿体/线粒体基因组的额外下载（仅 NCBI 源）
  - 可选基因组大小过滤
  - 输出统一 FASTA（all.genome.uniq.fasta + 平台锚点名 genome.fa）

用法:
  python host_genome.py --species "Lycium barbarum" --outdir ./host_genome
  python host_genome.py --species "Oryza sativa" --source gget --release 60
  python host_genome.py --species "Lycium barbarum" --include-organelles --ncbi-api xxx

（移植自 MMPV-RNA public_metadata_pipeline/download_host_genome.py；
 文件写入统一改走平台 Virus_Platform_Core.utils.safe_open；HTTP 统一走带域名校验的
 requests.Session；datasets CLI 不可用时自动改用 E-utilities + HTTPS 回退）
"""

import os
import sys
import re
import gzip
import json
import shutil
import argparse
import subprocess
from datetime import datetime
from typing import List, Optional, Set, Tuple

if sys.platform == 'win32':
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')

if __package__ in (None, ''):
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
    from Virus_Platform_Core.utils import safe_open
else:
    from ..utils import safe_open

ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
ESUMMARY_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
NCBI_FTP_HOST = "ftp://ftp.ncbi.nlm.nih.gov/"
NCBI_FTP_HTTPS = "https://ftp.ncbi.nlm.nih.gov/"

# gget/Ensembl 官方域名（gget.ref 返回的文件链接只落在这两台 FTP 主机上）
ENSEMBL_HOSTS = ('ftp.ensembl.org', 'ftp.ensemblgenomes.org')
GGET_INSTALL_HINT = 'pip install "gget==0.30.5"'


def _http(domains=('eutils.ncbi.nlm.nih.gov', 'ftp.ncbi.nlm.nih.gov')):
    """统一 HTTP 入口（仅允许指定官方域名，防 SSRF）。"""
    import requests

    allowed = tuple(domains)

    class _Scoped:
        """只放行白名单域名的极简 Session 包装。"""

        def __init__(self):
            self._s = requests.Session()

        @staticmethod
        def _check(url):
            host = url.split('//', 1)[-1].split('/', 1)[0]
            if host not in allowed:
                raise ValueError(f"仅允许 {'、'.join(allowed)} 域名，拒绝: {host}")
            return url

        def get(self, url, **kw):
            kw.setdefault('allow_redirects', False)
            return self._s.get(self._check(url), **kw)

        def download(self, url, dest, chunk=1 << 20, timeout=600,
                     exist_ok=True):
            """流式下载到 dest；已存在且非空则跳过；3xx 逐跳复核白名单后跟随。

            返回 True=实际下载，False=已存在跳过。"""
            if exist_ok and os.path.isfile(dest) and os.path.getsize(dest) > 1000:
                return False
            for _ in range(4):
                resp = self.get(url, stream=True, timeout=timeout)
                if resp.is_redirect or resp.is_permanent_redirect:
                    nxt = resp.headers.get('Location') or ''
                    resp.close()
                    if not nxt:
                        raise RuntimeError(f'重定向缺少 Location: {url}')
                    # Location 允许相对形式（'/path'、'../x'）：必须先 urljoin
                    # 再走白名单校验，否则 _check 把相对路径当主机名直接拒掉
                    url = urllib.parse.urljoin(url, nxt)
                    continue
                resp.raise_for_status()
                expected = resp.headers.get('Content-Length')
                got = 0
                with safe_open(dest, 'wb') as f:
                    for piece in resp.iter_content(chunk_size=chunk):
                        f.write(piece)
                        got += len(piece)
                # 截断校验（2026-10-06 实测事故：778MB 基因组 zip 下载中
                # 断，无 Content-Length 对照时静默落一个无 EOCD 的坏 zip，
                # 下游"未找到 FASTA"只剩一句 WARN，exit 0）
                if expected and got != int(expected):
                    try:
                        os.remove(dest)
                    except OSError:
                        pass
                    raise RuntimeError(
                        f'下载不完整（{got:,}/{int(expected):,} 字节），'
                        f'已丢弃残件，请重试: {url}')
                return True
            raise RuntimeError(f'重定向次数过多: {url}')

    return _Scoped()


def _emitters(log=None):
    """统一输出出口：log=None 打印控制台（CLI），否则推给回调（Web 任务面板）。"""
    if log is None:
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
# 终端UI
# ==========================================
class UI:
    CYAN = '\033[96m'; GREEN = '\033[92m'; YELLOW = '\033[93m'
    RED = '\033[91m'; PURPLE = '\033[95m'; GRAY = '\033[90m'
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


def check_datasets_cli() -> bool:
    """检查 NCBI datasets CLI 是否可用（PATH 或平台内置 3rd/tools/datasets/）"""
    try:
        from .datasets_cli import datasets_exe
        exe = datasets_exe()
    except ImportError:
        exe = 'datasets'
    if exe:
        try:
            result = subprocess.run(
                [exe, '--version'],
                capture_output=True, text=True, timeout=10
            )
            if result.returncode == 0:
                version = result.stdout.strip() or result.stderr.strip()
                UI.ok(f"NCBI datasets CLI 可用: {version}")
                return True
        except FileNotFoundError:
            pass
        except Exception:
            pass

    UI.warn("未找到 NCBI datasets CLI（将使用 E-utilities + HTTPS 直接下载；"
            "可执行 python main.py datasets-setup 自动安装）")
    return False


def open_fasta(path: str):
    """Context-managed opener for FASTA files (supports .gz compression)。

    errors='replace'：单条序列里一个坏字节不应让整份文件读崩（旧实现会直接
    抛 UnicodeDecodeError 丢掉整个源文件）。
    """
    if path.endswith('.gz'):
        return gzip.open(path, 'rt', encoding='utf-8', errors='replace')
    return open(path, 'r', encoding='utf-8', errors='replace')


def collect_fasta_files(directory: str) -> List[str]:
    """Recursively collect all FASTA/FA/FNA files (plain or .gz) under a directory."""
    files = []
    for root, _, filenames in os.walk(directory):
        for f in filenames:
            ext = os.path.splitext(f)[1]
            if f.endswith('.gz'):
                # 检查 .fna.gz 等双重后缀
                base = f[:-3]
                if any(base.endswith(e) for e in ['.fna', '.fasta', '.fa']):
                    files.append(os.path.join(root, f))
            elif ext in {'.fna', '.fasta', '.fa'}:
                files.append(os.path.join(root, f))
    return sorted(files)


def count_sequences(fasta_path: str) -> int:
    """Count sequences in a FASTA file by counting '>' header lines."""
    count = 0
    try:
        with open_fasta(fasta_path) as f:
            for line in f:
                if line.startswith('>'):
                    count += 1
    except Exception:
        pass
    return count


def merge_and_deduplicate(
    fasta_files: List[str],
    output_path: str,
    min_length: int = 0,
    prefix: str = '',
    log=None
) -> Tuple[int, int]:
    """
    合并多个 FASTA 文件并去重序列名

    返回值: (写入序列数, 重复序列名数)
    """
    E = _emitters(log)
    E.header("合并与去重")

    seen_ids: Set[str] = set()
    total_written = 0
    duplicate_count = 0
    total_bp = 0

    # FASTA 规范不允许 # 注释行（bowtie2-build 直接拒绝；kunpeng 注入 taxid
    # 时非 '>' 行一律当序列累积，# 行会变成垃圾碱基混进宿主库）——
    # 元数据写旁车 .meta.txt（移植自 MMPV-RNA 原版的同款修复）
    with safe_open(output_path + '.meta.txt', 'wt') as meta_f:
        meta_f.write("Merged host genome\n")
        meta_f.write(f"Created: {datetime.now().isoformat()}\n")
        meta_f.write(f"Source files: {len(fasta_files)}\n")

    with safe_open(output_path, 'wt') as out_f:
        for fa_path in fasta_files:
            basename = os.path.basename(fa_path)
            E.info(f"处理: {basename}")
            local_count = 0
            # 每源"输入序列数"与"写出序列数"分别计数：用于收尾核对
            # 输出总数 == 各源计数和（防空结果静默通过）
            n_in = 0

            try:
                with open_fasta(fa_path) as in_f:
                    current_seq = []
                    current_id = ''
                    current_full_header = ''

                    for line in in_f:
                        if line.startswith('>'):
                            # 写入上一条序列
                            if current_id and current_seq:
                                seq = ''.join(current_seq)
                                if len(seq) >= min_length:
                                    out_f.write(f">{current_full_header}\n")
                                    for i in range(0, len(seq), 60):
                                        out_f.write(seq[i:i+60] + '\n')
                                    total_written += 1
                                    total_bp += len(seq)

                            # 解析新序列头
                            full_header = line[1:].strip()
                            seq_id = full_header.split()[0]

                            # 序列名去重
                            if seq_id in seen_ids:
                                suffix = 1
                                while f"{seq_id}_dup{suffix}" in seen_ids:
                                    suffix += 1
                                new_id = f"{seq_id}_dup{suffix}"
                                seen_ids.add(new_id)
                                current_full_header = full_header.replace(seq_id, new_id, 1)
                                duplicate_count += 1
                            else:
                                seen_ids.add(seq_id)
                                current_full_header = full_header

                            current_id = seq_id
                            current_seq = []
                            local_count += 1
                            n_in += 1
                        else:
                            current_seq.append(line.strip())

                    # 处理最后一条序列
                    if current_id and current_seq:
                        seq = ''.join(current_seq)
                        if len(seq) >= min_length:
                            out_f.write(f">{current_full_header}\n")
                            for i in range(0, len(seq), 60):
                                out_f.write(seq[i:i+60] + '\n')
                            total_written += 1
                            total_bp += len(seq)

                # 核对口径：min_length 过滤会正常减少条数，这里只校验
                # "写出数 ≤ 读入数"，写出数多于读入数说明合并逻辑异常
                if local_count > n_in:
                    E.warn(f"  {basename}: 写出 {local_count} > 读入 {n_in}（异常）")
                E.ok(f"  {local_count} 条序列")

            except Exception as e:
                E.err(f"  处理 {basename} 时出错: {e}")

    file_size_mb = os.path.getsize(output_path) / (1024 * 1024)
    E.ok(f"合并完成: {total_written} 条序列, {total_bp:,} bp, {file_size_mb:.1f} MB")
    if duplicate_count:
        E.warn(f"重命名了 {duplicate_count} 个重复序列名")

    # 收尾核对：输出序列数必须等于各源计数和（去重只改名不丢序列）。
    # 不符说明某个源被静默跳过/截断（旧实现会产出空合并文件而不报错）。
    counted = sum(count_sequences(p) for p in fasta_files)
    out_counted = count_sequences(output_path)
    if not fasta_files:
        E.warn("合并核对：源文件列表为空，输出可能没有意义")
    elif out_counted != counted:
        msg = (f"合并核对不符：输出 {out_counted} 条 != 各源计数和 {counted} 条"
               f"（min_length={min_length} 过滤与去重改名会造成正常差异，"
               f"若差值明显请检查源文件）")
        E.warn(msg)
        try:                     # 不符时写进报告目录，便于事后追查
            with safe_open(os.path.join(os.path.dirname(output_path) or '.',
                                        'merge_check.txt'), 'wt') as f:
                f.write(msg + '\n')
                f.write(f"fasta_files={len(fasta_files)}\n")
                for p in fasta_files:
                    f.write(f"  {p}: {count_sequences(p)} 条\n")
                f.write(f"merged={output_path}: {out_counted} 条\n")
                f.write(f"total_written={total_written}\n")
        except (OSError, ValueError):
            pass
    else:
        E.ok(f"合并核对通过：输出 {out_counted} 条 == 各源计数和")

    return total_written, duplicate_count


def generate_genome_report(output_dir: str, fasta_files: List[str],
                           merged_fasta: str, log=None, source: str = 'ncbi',
                           meta: Optional[dict] = None):
    """生成基因组下载汇总报告"""
    E = _emitters(log)
    report_path = os.path.join(output_dir, 'genome_report.txt')

    total_seqs = count_sequences(merged_fasta)
    total_size = os.path.getsize(merged_fasta) if os.path.isfile(merged_fasta) else 0

    with safe_open(report_path, 'wt') as f:
        f.write("=" * 60 + "\n")
        f.write("  宿主参考基因组下载报告\n")
        f.write("=" * 60 + "\n")
        f.write(f"  生成时间: {datetime.now().isoformat()}\n")
        f.write(f"  下载来源: {source}\n")
        for k, v in sorted((meta or {}).items()):
            f.write(f"  {k}: {v}\n")
        f.write(f"  合并文件: {merged_fasta}\n")
        f.write(f"  平台锚点: {os.path.join(output_dir, 'genome.fa')}\n")
        f.write(f"  序列总数: {total_seqs}\n")
        f.write(f"  文件大小: {total_size / (1024*1024):.1f} MB\n")
        f.write(f"  源文件数: {len(fasta_files)}\n")
        f.write("\n  源文件列表:\n")
        for fa in fasta_files:
            f.write(f"    - {fa}\n")
        f.write("=" * 60 + "\n")

    E.ok(f"报告已保存: {report_path}")


def gunzip_full(path: str, log=None, max_layers: int = 3) -> str:
    """逐层解压 gzip 直到明文，返回明文文件路径；非压缩文件原样返回。

    为什么不假设只有一层：Ensembl Genomes 部分发布文件实测是**双层 .gz**
    （2026-09 拟南芥 TAIR10 release 63：外层 36MB → 内层 gzip → 121MB
    FASTA）。只解一层会留下一个仍带 gzip 魔数的文件，下游合并/建库全部
    静默产出空结果——按魔数逐层解到明文为止。"""
    E = _emitters(log)
    cur = path
    layers = 0
    while layers < max_layers:
        with open(cur, 'rb') as fh:
            if fh.read(2) != b'\x1f\x8b':
                break
        out = cur[:-3] if cur.lower().endswith('.gz') else cur + '.unz'
        if out.endswith('.'):
            out = out[:-1]
        with gzip.open(cur, 'rb') as fi, safe_open(out, 'wb') as fo:
            shutil.copyfileobj(fi, fo, 1 << 20)
        if layers:
            try:
                os.remove(cur)                    # 清掉中间层
            except OSError:
                pass
        cur = out
        layers += 1
    if layers:
        note = f"（检测到 {layers} 层嵌套 gzip）" if layers > 1 else ''
        E.ok(f"解压 → {os.path.basename(cur)} {note}")
    return cur


def finalize_host_dir(out_dir: str, merged_fasta: str, log=None) -> str:
    """把合并产物按平台约定落 genome.fa（建库锚点名），返回其路径。

    config.host_source_dirs / current_host_genome 与建库页都以
    host-db/<目录>/genome.fa 识别宿主源目录；旧流程只产
    all.genome.uniq.fasta，上下游在这里接通。"""
    E = _emitters(log)
    genome_fa = os.path.join(out_dir, 'genome.fa')
    shutil.copyfile(merged_fasta, genome_fa)
    E.ok(f"genome.fa 就绪（可直接构建宿主库）: {genome_fa}")
    return genome_fa


# ==========================================
# gget / Ensembl 源
# ==========================================
def normalize_species(species: str) -> str:
    """"Oryza sativa" → "oryza_sativa"（gget/Ensembl 命名）。"""
    return re.sub(r'\s+', '_', str(species).strip()).lower()


def _ensembl_https(url: str) -> str:
    """Ensembl/gget 链接强制 HTTPS。

    gget 返回的 genome_dna 链接实测有 http://（明文）与 ftp:// 形式；
    明文 http 会被中间人替换内容，而 _http() 只校验域名不校验协议，
    因此这里统一升级为 https（Ensembl 官方两台主机都支持 HTTPS）。
    """
    u = str(url or '').strip()
    if u.startswith('http://'):
        return 'https://' + u[len('http://'):]
    if u.startswith('ftp://'):
        return 'https://' + u[len('ftp://'):]
    return u


def parse_gget_ref(ref_dict: dict, species: str) -> dict:
    """解析 gget.ref() 返回值，取基因组 FASTA（dna）与注释（gtf）链接。"""
    info = (ref_dict or {}).get(species)
    if not info:
        raise RuntimeError(
            f'gget ref 返回中不含 {species}（Ensembl 无此物种或拼写不符？'
            '可用 gget.ref(None, list_species=True) 查看全部物种）')
    dna = info.get('genome_dna') or {}
    url = _ensembl_https(dna.get('ftp') or '')
    if not url:
        raise RuntimeError('Ensembl 未提供该物种的基因组 FASTA（dna）链接')
    return {
        'species': species,
        'dna_url': url,
        'gtf_url': _ensembl_https((info.get('annotation_gtf') or {}).get('ftp') or ''),
        'ensembl_release': dna.get('ensembl_release'),
        'release_date': dna.get('release_date', ''),
        'bytes': dna.get('bytes', ''),
    }


def gget_ref_lookup(species: str, release=None) -> dict:
    """gget.ref 封装：物种拉丁名（空格/下划线均可）→ 下载链接与版本元数据。

    gget 覆盖 Ensembl（脊椎动物）与 Ensembl Genomes（plants / fungi /
    metazoa / protists / bacteria），正好补 NCBI 覆盖不全的植物宿主。
    gget 自身的目录浏览请求固定访问 Ensembl 官方域名；文件下载在本模块
    的域名白名单内完成。"""
    try:
        import gget
    except ImportError as e:
        raise RuntimeError(
            f'gget 未安装（{GGET_INSTALL_HINT}）；gget 走 Ensembl / '
            'Ensembl Plants，适合 NCBI 覆盖不到的宿主物种') from e
    sp = normalize_species(species)
    try:
        ref = gget.ref(sp, which=['dna', 'gtf'], release=release, verbose=False)
    except ValueError as e:
        raise RuntimeError(
            f'Ensembl 未找到物种 {sp}（用拉丁学名，如 Oryza sativa / '
            f'oryza_sativa；可用 --source ncbi 走 NCBI）：{e}') from e
    return parse_gget_ref(ref, sp)


def download_genome_gget(species: str, out_dir: str, release=None,
                         taxid=None, log=None, min_length: int = 0) -> dict:
    """gget/Ensembl 宿主基因组下载（与 NCBI 流程同目录约定）。

    产物：extracted/（原始文件）、all.genome.uniq.fasta（合并去重）、
    genome.fa（平台建库锚点名）、genome_report.txt、gget_ref.json
    （Ensembl 版本与链接元数据，便于溯源/复现）。
    返回 dict：out_dir / genome_fa / taxid / ensembl_release / species。
    """
    E = _emitters(log)
    E.header(f"宿主参考基因组下载 (gget/Ensembl): {species}")
    os.makedirs(out_dir, exist_ok=True)

    meta = gget_ref_lookup(species, release)
    rel = meta['ensembl_release']
    E.ok(f"Ensembl release {rel} ({meta['release_date']})  {meta['bytes']}")
    E.info(f"基因组: {os.path.basename(meta['dna_url'])}")

    http = _http(ENSEMBL_HOSTS)
    extracted_dir = os.path.join(out_dir, 'extracted')
    os.makedirs(extracted_dir, exist_ok=True)
    fname = os.path.basename(meta['dna_url'].split('?')[0])
    dest = os.path.join(extracted_dir, fname)
    E.info(f"下载: {meta['dna_url']}")
    # 下载前记录 gget 元数据给出的字节数：下载后比对，不符 WARN（不阻断）
    expected_bytes = 0
    try:
        expected_bytes = int(str(meta.get('bytes') or '').replace(',', ''))
    except (TypeError, ValueError):
        expected_bytes = 0
    http.download(meta['dna_url'], dest)
    actual_bytes = os.path.getsize(dest)
    if expected_bytes and actual_bytes != expected_bytes:
        E.warn(f'字节数不符: 元数据 {expected_bytes} vs 实际 {actual_bytes}'
               f'（Ensembl 文件可能已更新或下载不完整，请核对）')
        meta['bytes_mismatch'] = {'expected': expected_bytes, 'actual': actual_bytes}
    else:
        E.ok(f"完成: {fname} ({actual_bytes / (1024*1024):.1f} MB)")

    plain = gunzip_full(dest, log=log)          # 逐层解到明文（含双层 .gz）
    if plain.endswith('.gz'):
        raise RuntimeError(f'解压后仍是 gzip: {plain}')
    fasta_files = [plain]                       # 只并明文文件，避免 .gz 重复计入

    merged_fasta = os.path.join(out_dir, 'all.genome.uniq.fasta')
    merge_and_deduplicate(fasta_files, merged_fasta,
                          min_length=min_length, log=log)
    genome_fa = finalize_host_dir(out_dir, merged_fasta, log=log)

    # Ensembl 不带 NCBI taxid：默认 NCBI taxonomy 在线查（建宿主库必需）
    if not taxid:
        E.info("查询 NCBI TaxID ...")
        taxid = _lookup_taxid(species)
        if taxid:
            E.ok(f"TaxID: {taxid}")
        else:
            E.warn('未能在线解析 TaxID（NCBI 不可达？），构建宿主库前请手工填写')

    meta.update({'source': 'gget', 'taxid': taxid,
                 'downloaded_at': datetime.now().isoformat(),
                 'genome_fa': genome_fa})
    with safe_open(os.path.join(out_dir, 'gget_ref.json'), 'wt') as f:
        json.dump(meta, f, ensure_ascii=False, indent=2)
    generate_genome_report(out_dir, fasta_files, merged_fasta, log=log,
                           source='gget',
                           meta={'Ensembl release': rel,
                                 'Release date': meta['release_date'],
                                 'URL': meta['dna_url'],
                                 'TaxID': taxid or '(未解析)'})

    E.header("下载完成")
    n_seqs = count_sequences(merged_fasta)
    size_mb = os.path.getsize(merged_fasta) / (1024 * 1024)
    E.ok(f"最终输出: {genome_fa}")
    E.ok(f"  {n_seqs} 条序列, {size_mb:.1f} MB")
    E.ok("下一步: 「数据库构建」页构建宿主库，或 "
         "python main.py build-host-db --genome ...")
    return {'out_dir': out_dir, 'genome_fa': genome_fa, 'taxid': taxid,
            'ensembl_release': rel, 'species': meta['species'],
            'source': 'gget'}


def download_genome_accession(accession: str, out_dir: str,
                              chromosomes: str = '', taxid=None,
                              log=None, min_length: int = 0) -> dict:
    """datasets/accession 模式宿主基因组下载（精确到组装号 + 染色体过滤）。

    例：--accession GCF_000001405.40 --chromosomes X,Y。目录约定与
    taxon/gget 模式一致；taxid 从 assembly_data_report.jsonl 自动解析，
    自动命名的 genome_<slug> 目录会改名为 <taxid>_<slug>（与 taxon 模式
    同约定）。"""
    E = _emitters(log)
    E.header(f"宿主参考基因组下载 (datasets/accession): {accession}")
    try:
        from . import datasets_cli as dsc
    except ImportError:
        import datasets_cli as dsc
    if not dsc.datasets_available():
        raise RuntimeError(
            'accession 模式需要 NCBI datasets CLI：执行 '
            'python main.py datasets-setup 自动安装（NCBI 官方 exe，'
            '落到 3rd/tools/datasets/）')
    os.makedirs(out_dir, exist_ok=True)

    res = dsc.download_genome(accession=accession, chromosomes=chromosomes,
                              include='genome,gff3', out_dir=out_dir, log=log)
    fasta_files = res['fasta'] or collect_fasta_files(
        os.path.join(out_dir, 'ncbi_dataset'))
    if not fasta_files:
        raise RuntimeError(
            f'数据包内没有基因组 FASTA（--chromosomes 过滤太严？当前过滤: '
            f'{chromosomes or "无"}）')

    merged_fasta = os.path.join(out_dir, 'all.genome.uniq.fasta')
    merge_and_deduplicate(fasta_files, merged_fasta,
                          min_length=min_length, log=log)
    genome_fa = finalize_host_dir(out_dir, merged_fasta, log=log)

    taxid = taxid or res.get('taxid')
    if taxid:
        E.ok(f"TaxID: {taxid}")
        # 自动命名的 genome_<slug> → <taxid>_<slug>（建库页按 taxid 识别）
        base = os.path.basename(os.path.normpath(out_dir))
        if base.startswith('genome_'):
            new_dir = os.path.join(os.path.dirname(os.path.normpath(out_dir)),
                                   f'{taxid}_{base[len("genome_"):]}')
            if os.path.normcase(new_dir) != os.path.normcase(out_dir) \
                    and not os.path.exists(new_dir):
                os.rename(out_dir, new_dir)
                out_dir = new_dir
                genome_fa = os.path.join(out_dir, 'genome.fa')
                E.ok(f"目录按 TaxID 重命名: {base} → {os.path.basename(new_dir)}")

    with safe_open(os.path.join(out_dir, 'datasets_host_meta.json'),
                   'wt') as f:
        json.dump({'source': 'ncbi-datasets', 'accession': accession,
                   'chromosomes': chromosomes, 'taxid': taxid,
                   'downloaded_at': datetime.now().isoformat(),
                   'genome_fa': genome_fa}, f, ensure_ascii=False, indent=2)
    generate_genome_report(out_dir, fasta_files, merged_fasta, log=log,
                           source='ncbi-datasets',
                           meta={'Accession': accession,
                                 'Chromosomes': chromosomes or '(全部)',
                                 'TaxID': taxid or '(未解析)'})

    E.header("下载完成")
    n_seqs = count_sequences(merged_fasta)
    E.ok(f"最终输出: {genome_fa}")
    E.ok(f"  {n_seqs} 条序列, "
         f"{os.path.getsize(merged_fasta) / (1024*1024):.1f} MB")
    E.ok("下一步: 「数据库构建」页构建宿主库，或 "
         "python main.py build-host-db --genome ...")
    return {'out_dir': out_dir, 'genome_fa': genome_fa, 'taxid': taxid,
            'accession': accession, 'source': 'ncbi-datasets'}


def _entrez_params(api_key=None, **extra):
    """E-utilities 参数（URL 由 requests params 序列化，不手工拼接）。"""
    p = dict(extra)
    if api_key:
        p['api_key'] = api_key
    return p


def _lookup_taxid(species: str, api_key=None):
    """用 NCBI esearch（taxonomy 库）查物种 taxid；失败返回 None。"""
    http = _http()
    term = f'"{species}"[Organism]'
    try:
        r = http.get(ESEARCH_URL,
                     params=_entrez_params(api_key, db='taxonomy', term=term,
                                           retmax=1, retmode='json'),
                     timeout=40)
        r.raise_for_status()
        ids = r.json().get('esearchresult', {}).get('idlist', [])
        return int(ids[0]) if ids else None
    except Exception:
        return None


def download_organelle_genome(species: str, organelle: str, out_dir: str,
                               ncbi_api: str = '') -> Optional[str]:
    """
    下载细胞器基因组（叶绿体/线粒体）

    参数:
        species: 物种名
        organelle: 'chloroplast' 或 'mitochondrion'
        out_dir: 输出目录
    """
    UI.info(f"正在下载{organelle}基因组...")

    org_dir = os.path.join(out_dir, organelle)
    os.makedirs(org_dir, exist_ok=True)

    if organelle == 'chloroplast':
        filter_word = 'chloroplast'
    else:
        filter_word = 'mitochondrion'
    search_term = f'"{species}"[Organism] AND {filter_word}[filter]'

    try:
        http = _http()
        resp = http.get(ESEARCH_URL,
                        params=_entrez_params(ncbi_api, db='nucleotide',
                                              term=search_term, retmode='json'),
                        timeout=30)
        data = resp.json()
        id_list = data.get('esearchresult', {}).get('idlist', [])

        if not id_list:
            UI.warn(f"  未找到 {species} 的 {organelle} 基因组")
            return None

        UI.ok(f"  找到 {len(id_list)} 条 {organelle} 记录")

        # 下载 FASTA
        ids_str = ','.join(id_list[:100])  # 限制100条
        resp = http.get(EFETCH_URL,
                        params=_entrez_params(ncbi_api, db='nucleotide',
                                              id=ids_str, rettype='fasta',
                                              retmode='text'),
                        timeout=120)
        fasta_content = resp.text

        out_name = f'{organelle}.fasta'
        out_path = os.path.join(org_dir, out_name)
        with safe_open(out_path, 'wt') as f:
            f.write(fasta_content)

        seq_count = fasta_content.count('>')
        UI.ok(f"  下载完成: {out_path} ({seq_count} 条序列)")
        return out_path

    except Exception as e:
        UI.err(f"  {organelle} 基因组下载失败: {e}")
        return None


def download_genome_https(species: str, out_dir: str, api_key=None,
                          max_assemblies: int = 5, log=None):
    """datasets CLI 不可用时的纯 Python 回退：E-utilities 搜 assembly →
    esummary 取 FTP 路径 → HTTPS 下载 RefSeq *_genomic.fna.gz 到 extracted/。

    filter 短语语法已实测核实（2026-10-03，db=assembly）：
      · `refseq[filter]`（有/无引号）在 esearch 里恒为 0 命中 → 旧查询条件
        实际上永远搜不到组装；
      · `"reference genome"[filter]` 是有效形式（全库 48826 条，
        N. benthamiana 命中 1 条）——加它即把"前 5 个任意组装"收紧为参考组装。
    下载到的组装 accession 清单会写入 out_dir/refseq_assemblies.json 供报告引用。
    """
    E = _emitters(log)
    http = _http()

    term = f'"{species}"[Organism] AND "reference genome"[filter]'
    r = http.get(ESEARCH_URL,
                 params=_entrez_params(api_key, db='assembly', term=term,
                                       retmax=max_assemblies, retmode='json'),
                 timeout=40)
    r.raise_for_status()
    ids = r.json().get("esearchresult", {}).get("idlist", [])
    if not ids:
        E.warn("E-utilities 未找到参考组装")
        return None
    E.ok(f"找到 {len(ids)} 个组装: {', '.join(ids)}")

    r = http.get(ESUMMARY_URL,
                 params=_entrez_params(api_key, db='assembly',
                                       id=",".join(ids), retmode='json'),
                 timeout=40)
    r.raise_for_status()
    docs = r.json().get("result", {})

    extracted_dir = os.path.join(out_dir, "extracted")
    os.makedirs(extracted_dir, exist_ok=True)
    n_downloaded = 0
    assembly_meta = []            # 合并清单：本次实际取用的组装 accession
    for uid in ids:
        d = docs.get(uid) or {}
        ftp = d.get("ftppath_refseq") or d.get("ftppath_genbank") or ""
        if not ftp:
            continue
        # 仅接受 NCBI 官方 FTP 主机（换 HTTPS），_http() 内再做域名校验
        if not ftp.startswith(NCBI_FTP_HOST):
            E.warn(f"跳过非 NCBI 路径: {ftp[:80]}")
            continue
        rel = ftp[len(NCBI_FTP_HOST):]
        dir_url = NCBI_FTP_HTTPS + rel
        try:
            lr = http.get(dir_url, timeout=60)
            lr.raise_for_status()
            fnames = re.findall(r'href="([^"]+_genomic\.fna\.gz)"', lr.text)
            if not fnames:
                continue
            fname = fnames[0]
            furl = dir_url + fname
            dest = os.path.join(extracted_dir, fname)
            E.info(f"下载: {furl}")
            if http.download(furl, dest):
                E.ok(f"完成: {fname}")
            else:
                E.ok(f"已存在: {fname}")
            n_downloaded += 1
            assembly_meta.append({
                'uid': uid,
                'assembly_accession': d.get('assemblyaccession') or '',
                'assembly_name': d.get('assemblyname') or '',
                'organism': d.get('organism') or species,
                'ftp_path': ftp,
                'file': fname,
            })
        except Exception as e:
            E.err(f"{dir_url} 下载失败: {e}")
    if not n_downloaded:
        E.warn("未下载到任何基因组 FASTA")
        return None
    # 合并清单（assembly accession 列表）写入输出 meta，报告/宿主库溯源都靠它
    if assembly_meta:
        try:
            with safe_open(os.path.join(out_dir, 'refseq_assemblies.json'),
                           'wt') as f:
                json.dump({'species': species, 'term': term,
                           'assemblies': assembly_meta}, f,
                          ensure_ascii=False, indent=1)
            E.ok('组装清单: ' + ', '.join(
                a['assembly_accession'] or a['uid'] for a in assembly_meta))
        except (OSError, ValueError):
            pass
    for fn in os.listdir(extracted_dir):
        if fn.endswith(".fna.gz"):
            src = os.path.join(extracted_dir, fn)
            dst_name = fn.replace('.fna.gz', '.fna')
            dst = os.path.join(extracted_dir, dst_name)
            if os.path.isfile(dst):
                continue
            with gzip.open(src, "rb") as fi, safe_open(dst, 'wb') as fo:
                fo.write(fi.read())
    E.ok(f"HTTPS 回退完成: {n_downloaded} 个组装 → {extracted_dir}")
    return extracted_dir


def main():
    parser = argparse.ArgumentParser(
        description="📥 宿主参考基因组下载器 — NCBI datasets 封装",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
  # 基础下载（NCBI）
  python host_genome.py --species "Lycium barbarum" --outdir ./host_genome

  # accession 精确下载（datasets，支持染色体过滤）
  python host_genome.py --accession GCF_000001405.40 --chromosomes X,Y

  # gget/Ensembl 源（含 Ensembl Plants，补 NCBI 覆盖不到的物种）
  python host_genome.py --species "Oryza sativa" --source gget
  python host_genome.py --species "Oryza sativa" --source gget --release 60 --taxid 39947

  # 包含细胞器基因组 + API key
  python host_genome.py --species "Lycium barbarum" \\
      --outdir ./host_genome --include-organelles --ncbi-api xxx

  # 仅验证/检查已有下载
  python host_genome.py --species "Lycium barbarum" \\
      --outdir ./host_genome --verify-only
        """
    )

    parser.add_argument('--species', default=None,
                        help='物种拉丁学名（taxon 模式必填）')
    parser.add_argument('--accession', default=None,
                        help='组装 accession（GCF_/GCA_，datasets 精确下载；'
                             '给出后忽略 --species）')
    parser.add_argument('--chromosomes', default=None,
                        help='染色体过滤（如 X,Y；仅 --accession 模式）')
    parser.add_argument('--source', choices=('ncbi', 'gget'), default='ncbi',
                        help='下载来源：ncbi=NCBI RefSeq（默认）；'
                             'gget=Ensembl/Ensembl Plants（补 NCBI 覆盖不到的物种）')
    parser.add_argument('--outdir', default=None,
                        help='输出目录（默认 host-db/<taxid>_<物种>/）')
    parser.add_argument('--release', type=int, default=None,
                        help='gget: Ensembl release（默认最新）')
    parser.add_argument('--taxid', type=int, default=None,
                        help='指定 NCBI TaxID（省去在线查询）')
    parser.add_argument('--ncbi-api', help='NCBI API Key (提升速率)')
    parser.add_argument('--include-organelles', action='store_true',
                        help='同时下载叶绿体和线粒体基因组')
    parser.add_argument('--min-length', type=int, default=0,
                        help='最小序列长度过滤 (bp)')
    parser.add_argument('--verify-only', action='store_true',
                        help='仅验证现有下载')
    parser.add_argument('--skip-datasets', action='store_true',
                        help='跳过 NCBI datasets 步骤 (使用已有文件)')

    args = parser.parse_args()

    if not args.species and not args.accession:
        parser.error('--species 与 --accession 至少给一个')
    if args.chromosomes and not args.accession:
        parser.error('--chromosomes 仅在 --accession 模式下可用')
    if not args.species:
        args.species = args.accession          # 兜底用于 slug/报告展示

    if not args.outdir:
        if args.accession:
            slug = re.sub(r'[^\w\-]+', '_', args.accession).strip('_')
            subdir = f'genome_{slug}'
        else:
            slug = re.sub(r'[^\w\-]+', '_', args.species).strip('_') or 'host'
            taxid = args.taxid or _lookup_taxid(args.species, args.ncbi_api)
            subdir = f'{taxid}_{slug}' if taxid else f'genome_{slug}'
        try:
            from Virus_Platform_Core.config import DIRS
            args.outdir = os.path.join(DIRS['host_src'], subdir)
        except Exception:
            args.outdir = subdir

    UI.header(f"宿主参考基因组下载: {args.species}")

    os.makedirs(args.outdir, exist_ok=True)

    # gget/Ensembl 源：独立入口（目录约定与 NCBI 源一致，产物直接供建库）
    if args.source == 'gget':
        if args.accession:
            UI.warn("gget/Ensembl 源不支持 accession（那是 NCBI 组装号），"
                    "已按 taxon 处理")
        if args.include_organelles:
            UI.warn("细胞器基因组仅 NCBI 源支持（Ensembl 不收录），已忽略")
        download_genome_gget(args.species, args.outdir, release=args.release,
                             taxid=args.taxid, min_length=args.min_length)
        return

    # accession 精确下载（datasets genome accession，支持染色体过滤）
    if args.accession:
        if args.include_organelles:
            UI.warn("细胞器基因组下载在 accession 模式下不适用，已忽略")
        download_genome_accession(args.accession, args.outdir,
                                  chromosomes=args.chromosomes,
                                  taxid=args.taxid,
                                  min_length=args.min_length)
        return

    # 验证模式
    if args.verify_only:
        UI.info("验证模式: 检查已有基因组文件...")
        merged_fasta = os.path.join(args.outdir, 'all.genome.uniq.fasta')
        genome_fa = os.path.join(args.outdir, 'genome.fa')
        if os.path.isfile(genome_fa):
            n_seqs = count_sequences(genome_fa)
            size_mb = os.path.getsize(genome_fa) / (1024 * 1024)
            UI.ok(f"已存在平台锚点 genome.fa: {n_seqs} 条序列, {size_mb:.1f} MB")
        elif os.path.isfile(merged_fasta):
            n_seqs = count_sequences(merged_fasta)
            size_mb = os.path.getsize(merged_fasta) / (1024 * 1024)
            UI.warn(f"已有合并文件但缺 genome.fa（旧版下载产物）: "
                    f"{n_seqs} 条序列, {size_mb:.1f} MB")
            UI.info("重跑一次下载即可自动补 genome.fa")
        else:
            extracted_dir = os.path.join(args.outdir, 'extracted')
            fasta_files = collect_fasta_files(extracted_dir)
            if fasta_files:
                UI.ok(f"找到 {len(fasta_files)} 个源 FASTA 文件，但尚未合并")
                for f in fasta_files:
                    UI.info(f"  {f}")
            else:
                UI.warn("未找到任何基因组文件")
        return

    # 步骤 1: 使用 NCBI datasets 下载
    if not args.skip_datasets:
        if not check_datasets_cli():
            UI.info("改用 E-utilities + HTTPS 直接下载 RefSeq 基因组 ...")
            try:
                download_genome_https(args.species, args.outdir, args.ncbi_api)
            except Exception as e:
                UI.err(f"HTTPS 回退下载失败: {e}")
        else:
            genome_zip = os.path.join(args.outdir, 'genome_down.zip')

            if os.path.isfile(genome_zip) and os.path.getsize(genome_zip) > 1000:
                UI.ok(f"基因组压缩包已存在: {genome_zip}")
            else:
                UI.header("步骤 1/3: NCBI datasets 下载")
                try:
                    from .datasets_cli import datasets_exe
                    ds_exe = datasets_exe() or 'datasets'
                except ImportError:
                    ds_exe = 'datasets'
                if args.accession:
                    cmd = [ds_exe, 'download', 'genome', 'accession',
                           args.accession]
                    if args.chromosomes:
                        cmd += ['--chromosomes', args.chromosomes]
                else:
                    cmd = [ds_exe, 'download', 'genome', 'taxon',
                           args.species]
                cmd += ['--filename', genome_zip,
                        '--include', 'genome,gff3,seq-report']
                if args.ncbi_api:
                    cmd.extend(['--api-key', args.ncbi_api])

                UI.info(f"执行: {' '.join(cmd)}")
                result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)

                if result.returncode != 0:
                    UI.err(f"datasets 下载失败: {result.stderr[:500]}")
                    UI.info("改用 E-utilities + HTTPS 直接下载 ...")
                    try:
                        download_genome_https(args.species, args.outdir, args.ncbi_api)
                    except Exception as e:
                        UI.err(f"HTTPS 回退下载失败: {e}")
                else:
                    UI.ok("基因组下载成功")

    # 步骤 2: 解压
    UI.header("步骤 2/3: 解压基因组")
    genome_zip = os.path.join(args.outdir, 'genome_down.zip')
    extracted_dir = os.path.join(args.outdir, 'extracted')

    if os.path.isfile(genome_zip):
        if not os.path.isdir(extracted_dir) or not os.listdir(extracted_dir):
            os.makedirs(extracted_dir, exist_ok=True)
            UI.info("正在解压...")
            # 先走 zipfile（标准库、无外部依赖）：旧实现直接调外部 unzip，
            # Windows 上 FileNotFoundError 没被捕获会直接崩，zipfile 兜底
            # 永不执行。zipfile 失败再退回 unzip 二进制（有则用，同样包 try）。
            import zipfile
            ok_unzip = False
            try:
                with zipfile.ZipFile(genome_zip, 'r') as zf:
                    zf.extractall(extracted_dir)
                UI.ok("解压完成 (Python zipfile)")
                ok_unzip = True
            except (OSError, zipfile.BadZipFile) as e:
                UI.warn(f"zipfile 解压失败({e})，尝试外部 unzip ...")
            if not ok_unzip:
                try:
                    result = subprocess.run(
                        ['unzip', '-o', genome_zip, '-d', extracted_dir],
                        capture_output=True, text=True
                    )
                    if result.returncode == 0:
                        UI.ok("解压完成 (unzip)")
                    else:
                        UI.err(f"解压失败: {(result.stderr or '')[:200]}")
                except (OSError, subprocess.SubprocessError) as e:
                    UI.err(f"外部 unzip 不可用: {e}")
        else:
            UI.ok("已解压，跳过")
    else:
        UI.warn(f"基因组压缩包不存在: {genome_zip}")

    # 步骤 3: 合并与去重
    UI.header("步骤 3/3: 合并与去重")

    fasta_files = collect_fasta_files(extracted_dir)
    if not fasta_files:
        UI.warn("未在解压目录中找到 FASTA 文件，尝试全目录搜索...")
        fasta_files = collect_fasta_files(args.outdir)

    merged_fasta = os.path.join(args.outdir, 'all.genome.uniq.fasta')
    if fasta_files:
        UI.ok(f"找到 {len(fasta_files)} 个 FASTA 文件:")
        for f in fasta_files:
            n_seqs = count_sequences(f)
            UI.info(f"  {os.path.basename(f)} — {n_seqs} 条序列")

        n_written, n_dups = merge_and_deduplicate(
            fasta_files, merged_fasta,
            min_length=args.min_length
        )
        generate_genome_report(args.outdir, fasta_files, merged_fasta)
    else:
        UI.warn("未找到任何 FASTA 文件")

    # 可选: 细胞器基因组（下载后重跑合并，把细胞器序列并入主文件）
    if args.include_organelles:
        UI.header("额外: 细胞器基因组下载")
        for org in ['chloroplast', 'mitochondrion']:
            org_fasta = download_organelle_genome(
                args.species, org, args.outdir, args.ncbi_api
            )
            if org_fasta and os.path.isfile(merged_fasta):
                UI.info(f"合并 {org} 基因组到主文件...")
                src_files = collect_fasta_files(extracted_dir)
                if org_fasta not in src_files:
                    src_files.append(org_fasta)
                merge_and_deduplicate(src_files, merged_fasta,
                                      min_length=args.min_length)
                UI.ok(f"{org} 基因组已合并")

    # 平台约定锚点：建库页 / config.host_source_dirs 以 genome.fa 识别宿主源目录
    if os.path.isfile(merged_fasta):
        finalize_host_dir(args.outdir, merged_fasta)

    UI.header("下载完成")
    if os.path.isfile(merged_fasta):
        n_seqs = count_sequences(merged_fasta)
        size_mb = os.path.getsize(merged_fasta) / (1024 * 1024)
        UI.ok(f"最终输出: {os.path.join(args.outdir, 'genome.fa')}")
        UI.ok(f"  {n_seqs} 条序列, {size_mb:.1f} MB")
        UI.ok("下一步: 数据库页构建宿主库（python main.py build-host-db --genome ...）")
    else:
        UI.warn("最终合并文件未生成，请检查错误日志")


if __name__ == '__main__':
    main()
