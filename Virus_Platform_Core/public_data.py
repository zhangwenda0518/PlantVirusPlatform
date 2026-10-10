# -*- coding: utf-8 -*-
"""
公共数据下载：SRA/ENA (SRR/ERR/DRR) 与 GSA/NGDC (CRR) 的链接解析与
批量下载管理（aria2c 优先、内置 HTTP 回退）。

安全约束（强制）：
- 仅允许 http/https；主机白名单（ENA 镜像 / NGDC 域）；所有出站请求前
  逐次校验：拒绝环回/私有/保留/链路本地地址，重定向每一跳重新校验
  （防 SSRF / DNS rebinding）。
- 批次目录与输出文件在创建时统一经 Virus_Platform_Core.utils.check_path 校验
  （拒绝 .. 段、限定平台 downloads/ 内），文件名严格白名单清洗。
- 完整性：ENA 提供的字节数比对 + gzip 魔数抽样。
- aria2c 低并发分段（-x4 -s4），NGDC 场景防 IP 封控（经验吸收自
  MMPV-RNA public_metadata_pipeline/gsa_sra.down.py）。
"""
import os
import re
import base64
import ipaddress
import json
import logging
import time
import gzip
import shutil
import socket
import random
import ssl
import subprocess
import threading
import urllib.request
import urllib.error

from .config import DIRS, PLATFORM_ROOT
from .utils import check_path, safe_open

# ------------------------------------------------------------------
# 常量与安全
# ------------------------------------------------------------------
ALLOWED_HOSTS = {
    'ftp.sra.ebi.ac.uk',            # ENA 读文件 HTTPS 镜像
    'www.ebi.ac.uk',                # ENA API（解析用）
    'ngdc.cncb.ac.cn',              # NGDC 搜索页（解析用）
    'download.cncb.ac.cn',          # GSA 下载
    'download.big.ac.cn',           # GSA 下载（big 镜像）
    'ddbj.nig.ac.jp',               # DDBJ（DRR 编号的原始 archive）
    'eutils.ncbi.nlm.nih.gov',      # NCBI eutils（DRR→DRA/DRX 元数据解析）
    'trace.ncbi.nlm.nih.gov',       # NCBI SRA STAT 预计算分类（病毒预筛，
                                    # 见 public_meta/prescreen_engine.py）
}
_IP_RE = re.compile(r'^\d{1,3}(\.\d{1,3}){3}$')
_UA = ('Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 '
       '(KHTML, like Gecko) Chrome/124.0 Safari/537.36')

DATA_FILE_RE = re.compile(
    r'\.(fastq\.gz|fq\.gz|fastq|fq|sra|fasta\.gz|fa\.gz|fasta|fa|fastq\.bz2)$', re.I)


def _is_private_ip(ip):
    """私有/保留地址判定（IPv4 + IPv6，同 ncbi_download.py 口径）。

    旧实现只按 IPv4 四段解析，IPv6 字符串必然 len(parts)!=4 → 一律判私有，
    导致 AAAA 记录（如 2600:1901:0:c831::）整条被拒。改用 ipaddress 统一判定：
    公网 IPv6 放行，内网/回环/链路本地/组播/保留/未指定一律拒绝；
    无法解析的字符串按不安全处理（True）。
    """
    try:
        a = ipaddress.ip_address(str(ip).strip())
    except ValueError:
        return True
    return bool(a.is_private or a.is_loopback or a.is_link_local
                or a.is_multicast or a.is_reserved or a.is_unspecified)


def check_url_host(url):
    """出站 URL 安全校验：仅 http(s)、主机在白名单、且不得解析到
    私有/保留地址。返回 (ok, msg)。所有网络请求前必须调用。"""
    m = re.match(r'^(https?)://([^/:?#]+)(?::(\d+))?([/?#].*)?$',
                 str(url or '').strip(), re.I)
    if not m:
        return False, f'仅支持 http/https URL: {url!r}'
    scheme, host = m.group(1).lower(), m.group(2).lower().rstrip('.')
    if _IP_RE.match(host) and _is_private_ip(host):
        return False, f'拒绝私有/保留 IP: {host}'
    if host == 'localhost' or host.endswith(('.local', '.internal')):
        return False, f'拒绝本机/内部主机名: {host}'
    if host not in ALLOWED_HOSTS:
        return False, (f'主机不在白名单内: {host}'
                       f'（允许: {", ".join(sorted(ALLOWED_HOSTS))}）')
    try:
        infos = socket.getaddrinfo(host, 443 if scheme == 'https' else 80,
                                   proto=socket.IPPROTO_TCP)
    except OSError as e:
        return False, f'主机解析失败: {host} ({e})'
    for info in infos:
        if _is_private_ip(info[4][0]):
            return False, f'{host} 解析到私有/保留地址: {info[4][0]}'
    return True, ''


class _SafeRedirectHandler(urllib.request.HTTPRedirectHandler):
    """重定向每一跳重新过白名单与私网校验。"""

    def redirect_request(self, req, fp, code, msg, headers, newurl):
        ok, why = check_url_host(newurl)
        if not ok:
            raise urllib.error.URLError(f'重定向目标被拒绝: {why}')
        return super().redirect_request(req, fp, code, msg, headers, newurl)


_OPENER = urllib.request.build_opener(_SafeRedirectHandler)

# NGDC 反爬验证页（2026-09-17 上线）：对无 Referer 的脚本请求返回一张
# 静态"正在进行安全检查"页（实测 nginx 按 Referer 放行，auth_token
# Cookie 是其页内 JS 的自检产物）。GSA 的 CRR 详情页正因此解析不出
# 下载链接——搜索页能过、browse 详情页被拦，症状恰是"未找到公开文件"。
_NGDC_CHALLENGE_MARK = '正在进行安全检查'
_NGDC_SALT_RE = re.compile(r'var\s+secret\s*=\s*"([^"]+)"\s*\+\s*dateStr')


def _ngdc_auth_token(salt):
    """按验证页 JS 的算法生成 auth_token：base64(salt + UTC 日期)，去
    padding。JS 用 toISOString() 取的是 UTC 日期——用本地日期算出的
    token 过不了当天 UTC 早于本地的时段。"""
    today = time.strftime('%Y-%m-%d', time.gmtime())
    return base64.b64encode((salt + today).encode()).decode().rstrip('=')


def _http_get(url, timeout=40, retries=3, _challenge_html=None):
    """受限 GET：请求前逐次校验（白名单 + 解析 IP + 重定向逐跳校验）。

    - NGDC 域一律带站内 Referer（无 Referer 会被反爬验证页拦截）；
      万一仍命中验证页，按页内 JS 的 salt 算 auth_token Cookie 重放一次。
    - SSL 瞬断 / 超时 / 连接重置按指数退避重试：批量解析时一次网络抖动
      不应把整条 accession 打成"解析失败"。
    """
    ok, why = check_url_host(url)
    if not ok:
        raise ValueError(f'出站请求被拒绝: {why}')
    headers = {'User-Agent': _UA, 'Accept': '*/*'}
    m = re.match(r'^https?://([^/:?#]+)', url, re.I)
    if m and m.group(1).lower().endswith(('.cncb.ac.cn', '.big.ac.cn')):
        headers['Referer'] = 'https://ngdc.cncb.ac.cn/gsa/search'
    if _challenge_html is not None:
        salt = _NGDC_SALT_RE.search(_challenge_html)
        token = _ngdc_auth_token(salt.group(1) if salt else 'my_nginx_salt_')
        headers['Cookie'] = 'auth_token=' + token
    last_err = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers=headers)
            with _OPENER.open(req, timeout=timeout) as r:
                html = r.read().decode('utf-8', 'replace')
            if _NGDC_CHALLENGE_MARK in html and _challenge_html is None:
                return _http_get(url, timeout=timeout, retries=retries,
                                 _challenge_html=html)
            return html
        except urllib.error.HTTPError:
            raise
        except (urllib.error.URLError, ssl.SSLError, TimeoutError,
                ConnectionError) as e:
            last_err = e
            if attempt < retries - 1:
                time.sleep(1.5 * (2 ** attempt) + random.random())
    raise last_err


def aria2c_path():
    """内置 aria2c.exe（bin/，兼容旧版根目录布局）；找不到再试 PATH。"""
    for p in (os.path.join(PLATFORM_ROOT, '3rd', 'bin', 'aria2c.exe'),
              os.path.join(PLATFORM_ROOT, 'aria2c.exe')):
        if os.path.isfile(p):
            return p
    return shutil.which('aria2c') or shutil.which('aria2c.exe')


def sracha_path():
    """内置 sracha.exe（bin/，官方 0.7.0；兼容旧版根目录布局）；
    找不到再试 PATH。"""
    for p in (os.path.join(PLATFORM_ROOT, '3rd', 'bin', 'sracha.exe'),
              os.path.join(PLATFORM_ROOT, 'sracha.exe')):
        if os.path.isfile(p):
            return p
    return shutil.which('sracha') or shutil.which('sracha.exe')


def fasterq_path():
    for p in (os.path.join(PLATFORM_ROOT, '3rd', 'bin', 'fasterq-dump.exe'),
              os.path.join(PLATFORM_ROOT, 'fasterq-dump.exe')):
        if os.path.isfile(p):
            return p
    return shutil.which('fasterq-dump') or shutil.which('fasterq-dump.exe')


def sra_convert_engine():
    """选定 .sra→FASTQ 转换引擎：优先 sracha（纯 Rust 单文件，比
    fasterq-dump 快 5-13 倍，直接输出 gzip），回退 fasterq-dump。
    返回 (引擎名, 可执行路径) 或 (None, None)。"""
    p = sracha_path()
    if p and os.path.isfile(p):
        return 'sracha', p
    w = shutil.which('sracha') or shutil.which('sracha.exe')
    if w:
        return 'sracha', w
    p = fasterq_path()
    if p:
        return 'fasterq-dump', p
    return None, None


def _safe_filename(url):
    """从 URL 取文件名并强制白名单清洗（杜绝路径段与保留名）。"""
    fn = os.path.basename(str(url).split('?')[0].split('#')[0])
    fn = re.sub(r'[^A-Za-z0-9._\-]+', '_', fn).strip('._')
    if not fn or fn in {'.', '..'} or not DATA_FILE_RE.search(fn):
        raise ValueError(f'URL 文件名不合规（需数据文件后缀）: {url!r}')
    return fn[:120]


# ------------------------------------------------------------------
# Accession → 文件 URL 解析
# ------------------------------------------------------------------
_ENA_API = ('https://www.ebi.ac.uk/ena/portal/api/filereport'
            '?accession={acc}&result=read_run'
            '&fields=run_accession,fastq_ftp,fastq_md5,fastq_bytes,'
            'scientific_name,sample_accession&format=tsv')


def _dedupe_urls(urls):
    """URL 去重 + NGDC 页面常见的 '//SRR' 双斜杠修正。

    旧实现是 re.sub(r'(?<=//SRR)', '/SRR', u)——在 '//SRR' 之后**插入** '/SRR'，
    '//SRR1000' 反而变成 '//SRR/SRR1000'（越改越错）。这里用负向后顾只吃掉
    多余的那一条斜杠：'//SRR1000' → '/SRR1000'，且不动 'https://' 的双斜杠。
    """
    out, seen = [], set()
    for u in urls:
        u = u.strip()
        if u.startswith('//'):
            u = 'https:' + u
        # 修正路径里的重复斜杠（//SRR / //ERR / //DRR / //CRR / //RR），保留协议 // 
        u = re.sub(r'(?<!:)//(?=S?RR|ERR|DRR|CRR|RR)', '/', u)
        if u not in seen:
            seen.add(u)
            out.append(u)
    return out


def _resolve_ena(acc):
    """SRR/ERR/DRR → ENA filereport → [{url, bytes}]。acc 已过严格正则。"""
    txt = _http_get(_ENA_API.format(acc=acc))
    lines = [l for l in txt.splitlines() if l.strip()]
    if len(lines) < 2:
        raise ValueError(f'ENA 未返回 {acc} 的数据行（run 可能不存在或未公开）')
    header = lines[0].split('\t')
    row = lines[1].split('\t')
    rec = dict(zip(header, row))
    ftps = [x for x in (rec.get('fastq_ftp') or '').split(';') if x]
    byts = [x for x in (rec.get('fastq_bytes') or '').split(';') if x]
    md5s = [x for x in (rec.get('fastq_md5') or '').split(';') if x]
    if not ftps:
        # ENA 未生成 FASTQ.GZ（常见于刚公布 / 发布者仅提供 SRA 原始数据）。
        # 只要 NCBI SRA 有数据（INSDC 三库元数据互通），就用 sracha
        # 从 NCBI .sra 拉取（--prefer-ena 会回落 NCBI）。
        # 不再直接判“受控数据”判死，避免“原来能下”的 run 因 ENA 无 fastq
        # 而整批失败。
        if re.fullmatch(r'[SED]RR\d+', acc) and sracha_path():
            return {'db': 'ENA→SRA', 'organism': rec.get('scientific_name', ''),
                    'bio_sample': rec.get('sample_accession', ''),
                    'files': [{'url': '', 'sracha': True}]}
        raise ValueError(f'ENA 无 {acc} 的 fastq 文件（可能为受控数据）')
    files = []
    for i, f in enumerate(ftps):
        url = 'https://' + f
        ok, msg = check_url_host(url)
        if not ok:
            raise ValueError(f'ENA 返回了不允许的地址: {msg}')
        try:
            size = int(byts[i]) if i < len(byts) else 0
        except ValueError:
            size = 0
        # ENA 已请求 fastq_md5，此前被丢掉（条目里 md5 恒为空 → 下载后无从校验）
        md5 = str(md5s[i]).strip() if i < len(md5s) else ''
        if not re.fullmatch(r'[0-9a-fA-F]{32}', md5 or ''):
            md5 = ''
        files.append({'url': url, 'bytes': size, 'md5': md5})
    # 单端库有时同时列出 .fastq.gz 与 _1/_2 —— _1/_2 齐全时只保留配对
    paired = [x for x in files if re.search(r'_1\.', x['url'])
              or re.search(r'_2\.', x['url'])]
    if len(paired) >= 2:
        files = paired
    return {'db': 'ENA', 'organism': rec.get('scientific_name', ''),
            'bio_sample': rec.get('sample_accession', ''), 'files': files}


def _resolve_ddbj(acc):
    """DRR → DDBJ 公共 FASTQ（.fastq.bz2）。

    DDBJ 与 ENA/NCBI 只共享元数据，原始文件留在提交方 archive，
    因此 DRR 编号必须走 DDBJ。布局：
      dra/fastq/<DRA前6>/<DRA>/<DRX>/<acc>.fastq.bz2        （单端）
      dra/fastq/<DRA前6>/<DRA>/<DRX>/<acc>_1/_2.fastq.bz2   （双端）
    DRA/DRX 号经 NCBI eutils 元数据解析（INSDC 三库元数据互通）。
    """
    import urllib.parse as _up
    term = _up.quote(acc)
    txt = _http_get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
                    f'esearch.fcgi?db=sra&term={term}')
    m = re.search(r'<Id>(\d+)</Id>', txt)
    if not m:
        raise ValueError(f'NCBI 元数据中未找到 {acc}（可能未公开）')
    xml = _http_get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
                    f'efetch.fcgi?db=sra&id={m.group(1)}&rettype=xml')
    dra = re.search(r'(DRA\d{6})', xml)
    drx = re.search(r'(DRX\d{6})', xml)
    if not (dra and drx):
        raise ValueError(f'NCBI 元数据缺少 {acc} 的 DRA/DRX 号')
    base = ('https://ddbj.nig.ac.jp/public/ddbj_database/dra/fastq/'
            f'{dra.group(1)[:6]}/{dra.group(1)}/{drx.group(1)}/')
    try:
        idx = _http_get(base)
    except Exception as e:
        raise ValueError(f'DDBJ 目录不可达: {e}') from e
    names = sorted(set(re.findall(
        r'href="(' + re.escape(acc) + r'[^"/]*\.fastq\.bz2)"', idx)))
    if not names:
        raise ValueError(f'DDBJ 暂无 {acc} 的 FASTQ（新提交的 bz2 生成有延迟，'
                         f'可稍后重试）')
    return {'db': 'DDBJ',
            'organism': '',
            'files': [{'url': base + n, 'bytes': 0, 'md5': ''}
                      for n in names]}


def _resolve_ngdc(acc):
    """CRR → NGDC 搜索页 → browse 详情页 → 下载 URL（吸收 gsa_sra.down.py 流程）。

    详情页抓取由 _http_get 自动过反爬验证页（Referer + auth_token）。"""
    html = _http_get(f'https://ngdc.cncb.ac.cn/gsa/search?searchTerm={acc}')
    urls = [u for u in _dedupe_urls(re.findall(
        r'(?:(?:ftp|https?):)?//download[0-9]*\.(?:cncb|big)\.ac\.cn/[^"\'<>\s]+',
        html)) if acc in u]
    browse_re = (rf'href="((?:https?://ngdc\.cncb\.ac\.cn)?'
                 rf'(?:/gsa/)?browse/[^"]*{acc}[^"]*)"')
    m = re.search(browse_re, html)
    if not m:
        html2 = _http_get(f'https://ngdc.cncb.ac.cn/gsa/search?db=GSA&term={acc}')
        m2 = re.search(browse_re, html2)
        path = m2.group(1) if m2 else None
    else:
        path = m.group(1)
    if path and not urls:
        if path.startswith('http'):
            detail = path
        elif path.startswith('/gsa/'):
            detail = 'https://ngdc.cncb.ac.cn' + path
        else:
            detail = 'https://ngdc.cncb.ac.cn/gsa/' + path.lstrip('/')
        html3 = _http_get(detail)
        urls = [u for u in _dedupe_urls(re.findall(
            r'(?:(?:ftp|https?):)?//download[0-9]*\.(?:cncb|big)\.ac\.cn/[^"\'<>\s]+',
            html3)) if acc in u]
    files = []
    rejected = []
    for u in urls:
        u = re.sub(r'^ftp://', 'https://', u)      # 统一 HTTPS（协议白名单）
        if not DATA_FILE_RE.search(u):
            continue
        ok, msg = check_url_host(u)
        if not ok:
            # 单个外站/非白名单 URL 不应整批判死：记录后继续看其余候选，
            # 只有全部候选都被拒才报错（页面里常混有 CDN/统计域名）
            rejected.append(f'{u}（{msg}）')
            continue
        if all(f['url'] != u for f in files):
            files.append({'url': u, 'bytes': 0})
    if rejected:
        logging.getLogger(__name__).warning(
            'NGDC %s：忽略 %d 个不在白名单/校验不通过的候选 URL: %s',
            acc, len(rejected), '; '.join(rejected[:3]))
    # cncb/big 双镜像去重：同文件名优先 cncb 主站
    by_name = {}
    for f in files:
        bn = os.path.basename(f['url'])
        cur = by_name.get(bn)
        if cur is None or ('cncb.ac.cn' in f['url'] and 'cncb.ac.cn' not in cur['url']):
            by_name[bn] = f
    files = list(by_name.values())
    if not files:
        if rejected:
            # 全部候选都被拒才报错（此时把拒绝原因带出来便于排查）
            raise ValueError(f'NGDC 返回的 {acc} 候选地址全部未通过校验：'
                             + '; '.join(rejected[:3]))
        raise ValueError(f'NGDC 未找到 {acc} 的公开下载文件（可能受控/未公开）')
    return {'db': 'GSA', 'organism': '', 'bio_sample': '', 'files': files}


_ACC_RE_ENA = re.compile(r'[SED]RR\d+')
_ACC_RE_CRR = re.compile(r'CRR\d+')


def resolve_accession(acc):
    """单条 accession → {db, organism, files:[{url,bytes}]}；
    不支持的编号抛 ValueError（CRA 项目请先取得 CRR 运行号列表）。"""
    acc = str(acc).strip().upper()
    if _ACC_RE_ENA.fullmatch(acc):
        if acc.startswith('DRR'):
            # DRR 归 DDBJ 管理，ENA 通常无数据；先走 DDBJ，失败回退 ENA
            try:
                return _resolve_ddbj(acc)
            except Exception:
                return _resolve_ena(acc)
        return _resolve_ena(acc)
    if _ACC_RE_CRR.fullmatch(acc):
        return _resolve_ngdc(acc)
    raise ValueError(f'不支持的编号格式: {acc}（支持 SRR/ERR/DRR/CRR）')


# ------------------------------------------------------------------
# 下载管理器（写入统一走 safe_open：内部再做一次路径与平台边界校验）
# ------------------------------------------------------------------
def _batch_root():
    return DIRS.get('downloads') or os.path.join(PLATFORM_ROOT, 'downloads')


def _sanitize(name):
    s = re.sub(r'[^A-Za-z0-9_\-.]+', '_', str(name)).strip('._-')
    return s[:40] or 'batch'


def _gzip_ok(path):
    """gzip 魔数 + 抽样解压：快速完整性检查。"""
    try:
        with open(path, 'rb') as f:
            if f.read(2) != b'\x1f\x8b':
                return False
        with gzip.open(path, 'rb') as f:
            f.read(1 << 16)
        return True
    except (OSError, EOFError):
        return False


def _bz2_ok(path):
    """bz2 流式解压校验（读全量）：返回 (ok, msg)。

    DDBJ 的 .fastq.bz2 此前没有任何完整性校验（既无 ENA 的 md5，也无数字节
    可比），半截文件会一路流到转换链才炸。这里按"能否完整解压"判定：
    读全量算通过，ReadError/EOFError/OSError 即 broken。
    """
    import bz2
    try:
        with bz2.open(path, 'rb') as f:
            while f.read(1 << 20):
                pass
        return True, ''
    except (OSError, EOFError, ValueError) as e:
        return False, f'{type(e).__name__}: {e}'


def _md5_ok(path, expect, max_bytes=2 * 1024 ** 3):
    """本地算 MD5 与期望值比对 → (ok, note)。

    ok: True=一致 / False=不一致 / None=未校验（超过 max_bytes 跳过或读失败，
    此时 note 给出原因，调用方写入条目标注）。
    """
    import hashlib
    try:
        if os.path.getsize(path) > max_bytes:
            return None, f'skipped_size>{max_bytes // 1024 ** 3}GB'
        h = hashlib.md5()
        with open(path, 'rb') as f:
            for chunk in iter(lambda: f.read(1 << 20), b''):
                h.update(chunk)
        got = h.hexdigest()
    except OSError as e:
        return None, f'error_{type(e).__name__}'
    if got.lower() == str(expect).strip().lower():
        return True, ''
    return False, f'md5 不符（期望 {expect}, 实际 {got}）'


def _sracha_owned_name(n, acc):
    """sracha 产物是否属于该 acc（边界明确的精确形式）。

    旧实现用 `n.startswith(acc)` 归属产物：SRR1000 会吃掉 SRR10001 的文件
    （'SRR10001_1.fastq.gz'.startswith('SRR1000') 为真），把别人的数据登记到
    本条目下。这里改成白名单式精确匹配：
      <acc>.sra / <acc>.fastq.gz / <acc>.fq.gz / <acc>_1|_2.fastq.gz|fq.gz
    """
    if not n or not acc:
        return False
    if n == f'{acc}.sra':
        return True
    for suf in ('.fastq.gz', '.fq.gz'):
        if n == acc + suf:
            return True
        if n in (f'{acc}_1{suf}', f'{acc}_2{suf}'):
            return True
    return False


def _known_total(fe):
    """条目总量（字节）：优先元数据精确值，未知时退回下载器实时近似值。

    为什么必须分开两个字段：`size` 早期被 aria2c 的整数 MiB 缩写值覆盖
    （'8.3MiB' = 8703180，真值 8772679），前端「已下载 / 总量」因此出现
    倒挂（8.8 MB / 8.7 MB），批次总量还会中途掉头。`expect_bytes` 来自
    ENA/DDBJ 元数据，全程不被覆盖，是唯一可信分母。
    """
    return int(fe.get('expect_bytes') or 0) or int(fe.get('size') or 0)


def _fe_weight(fe):
    """进度加权：有精确/近似总量的条目用字节数加权。

    下载批次里各文件体积可能差几个数量级，按条目数取平均会让大文件的
    进度被小文件稀释。没有总量的条目返回 0，由调用方用均值兜底。
    """
    t = _known_total(fe)
    return t if t > 0 else 0


def _sra_fastq_outputs(bdir, base, exts=('.fastq.gz', '.fq.gz')):
    """列出 .sra 解码产物文件名（相对批次目录），单端/双端都认。

    单端：`<base>.fastq.gz`（sracha）或 `<base>.fastq`（fasterq-dump）
    双端：`<base>_1/_2.fastq.gz` / `<base>_1/_2.fastq`

    只匹配 `<base>_` 是历史 bug：sracha 解码单端 run 写的是无下划线的
    `<acc>.fastq.gz`，于是转换链报「无输出」并把已经产出的 FASTQ 丢掉。
    """
    out = []
    for n in sorted(os.listdir(bdir)):
        stem = None
        for e in exts:
            if n.endswith(e):
                stem = n[:-len(e)]
                break
        if stem is None:
            continue
        if stem == base or stem.startswith(base + '_'):
            out.append(n)
    return out


def _retryable(fe):
    """该条目能否重跑下载。

    必须覆盖两条路径，缺一不可：
      · 有 URL 的普通条目（ENA/NGDC/DDBJ/自定义 URL）；
      · **无 URL 的 sracha 条目**（acc.sra，out 已定好）—— 这正是
        ENA 无预生成 FASTQ 时的主路径，只看 `url` 会让它永远无法重试。
    解析期就失败的条目 `out` 为空、也不是 INSDC 编号，不可重试。
    """
    if not (fe.get('out') and fe.get('acc')):
        return False
    if fe.get('url') or fe.get('sracha'):
        return True
    return bool(re.fullmatch(r'[SED]RR\d+', str(fe.get('acc') or '')))


# 批次终态：一旦落到这些状态就不再被 cancel 改写。
# 动因：阶段二把下载批次挂进任务中心后，「全局重置」会遍历运行中任务并
# 逐条转发停止 —— 若某张卡片的运行中判断有偏差，一次误转发就会把**已经
# 下载完成**的批次改写成「已取消」，用户会以为数据白下了。已完成/已取消/
# 失败/已中断（可续传）这四种状态都携带用户需要的信息，一律保护起来。
_TERMINAL_BATCH_STATES = ('completed', 'cancelled', 'failed', 'interrupted')


class DownloadManager:
    """批量下载：文件队列 + N 并发 aria2c（内置 HTTP 回退）+ 进度状态。

    batch = {id, name, status, created, dir(已校验), n_files, done, failed,
             overall, files: [{acc, url, out, size, done_bytes, progress,
                               status, md5, verified, error, db, organism}]}
    状态: resolving → downloading → (converting) → completed / failed /
          cancelled / interrupted
    """

    def __init__(self):
        self.lock = threading.RLock()
        self.batches = {}
        self.procs = {}
        self._spd = {}          # out -> [last_bytes, last_ts, ema_bytes_per_s]
        self._listeners = []    # 批次新增时的回调（web 层用来挂进任务中心）
        self._load_all()

    def add_listener(self, fn):
        """登记「新批次」回调，并对**已存在**的批次各回调一次。

        供 web 层把下载批次挂进任务中心（外部任务）。回调只传批次 id，
        web 层自己来取快照——core 不反向 import web，依赖保持单向。
        """
        with self.lock:
            if fn not in self._listeners:
                self._listeners.append(fn)
            ids = list(self.batches.keys())
        for bid in ids:
            try:
                fn(bid)
            except Exception:
                pass

    def _notify(self, bid):
        with self.lock:
            fns = list(self._listeners)
        for fn in fns:
            try:
                fn(bid)
            except Exception:
                pass

    # ---- 下载速率（字节/秒，EMA 平滑）----
    # 设计要点：不依赖 aria2c 的输出格式。任何一路下载只要周期性更新
    # done_bytes，就能算出速率 —— 内置 HTTP 回退路径同样生效。
    _SPD_ALPHA = 0.3

    @staticmethod
    def _aria_bytes(s):
        """aria2c 容量串（'12MiB' / '1.2GiB' / '0B'）→ 字节数。"""
        m = re.match(r'^\s*(\d+(?:\.\d+)?)\s*([KMGT]?)i?B\s*$', s or '', re.I)
        if not m:
            return 0
        mult = {'': 1, 'K': 1024, 'M': 1024 ** 2,
                'G': 1024 ** 3, 'T': 1024 ** 4}
        return int(float(m.group(1)) * mult[m.group(2).upper()])

    @staticmethod
    def _aria_eta(s):
        """aria2c ETA 串（'9m53s' / '25s' / '1h2m'）→ 秒；无法解析返回 None。"""
        t = (s or '').strip().lower()
        if not t:
            return None
        m = re.fullmatch(r'(?:(\d+)h)?(?:(\d+)m)?(?:(\d+)s)?', t)
        if not m or not any(m.groups()):
            return None
        h, mi, se = (int(x) if x else 0 for x in m.groups())
        return h * 3600 + mi * 60 + se

    def _tick_speed(self, out, done, alpha=None):
        """记录一次采样，返回平滑后的速率（B/s）。

        采样间隔 <0.5s 时不更新（避免抖动），但也不清零 —— 保留上一次
        平滑值，前端显示不会闪成 0。

        alpha: EMA 权重，缺省 _SPD_ALPHA。sracha 路径以 8MiB 分块完成
        为单位推进，瞬时差分天然是脉冲式的，需要更长的平滑窗口。
        """
        now = time.time()
        a = self._SPD_ALPHA if alpha is None else alpha
        prev = self._spd.get(out)
        if prev is None:
            self._spd[out] = [done, now, 0.0]
            return 0.0
        last_bytes, last_ts, ema = prev
        dt = now - last_ts
        if dt >= 0.5:
            if done >= last_bytes:
                inst = (done - last_bytes) / dt
                ema = inst if ema <= 0 else (a * inst + (1 - a) * ema)
            prev[0], prev[1], prev[2] = done, now, ema
        return prev[2]

    def _prime_speed(self, out, done):
        """同步 EMA 基准（采用外部速率时调用，避免回退时瞬间跳变）。"""
        prev = self._spd.get(out)
        now = time.time()
        if prev is None:
            self._spd[out] = [done, now, 0.0]
        else:
            prev[0], prev[1] = done, now

    def _drop_speed(self, out):
        self._spd.pop(out, None)

    # ---- 持久化 ----
    def _save(self, b):
        with self.lock:
            try:
                meta = check_path(b['dir'] + os.sep + 'batch.json',
                                  must_exist=False, in_platform=True)
                os.makedirs(b['dir'], exist_ok=True)
                # 原子落盘：前端/任务中心会并发读 batch.json，直接覆写会让读端
                # 拿到半截 JSON（下载页整页 500）。先写 .tmp 再 os.replace。
                tmp = meta + '.tmp'
                with safe_open(tmp, 'wt') as f:
                    json.dump(b, f, ensure_ascii=False, indent=1)
                os.replace(tmp, meta)
            except (OSError, ValueError):
                pass

    def _blog(self, b, msg):
        """追加一条批次过程日志到 <批次目录>/batch.log。

        任何写入失败都静默吞掉，绝不影响下载流程本身。"""
        try:
            d = b.get('dir') or ''
            if not d:
                return
            os.makedirs(d, exist_ok=True)
            with open(os.path.join(d, 'batch.log'), 'a',
                      encoding='utf-8') as f:
                f.write(f'[{time.strftime("%H:%M:%S")}] {msg}\n')
        except OSError:
            pass

    def _load_all(self):
        root = _batch_root()
        try:
            os.makedirs(root, exist_ok=True)
        except OSError:
            return
        for name in sorted(os.listdir(root)):
            try:
                meta = check_path(os.path.join(root, name, 'batch.json'),
                                  must_exist=False, in_platform=True)
            except ValueError:
                continue
            if not os.path.isfile(meta):
                continue
            try:
                with open(meta, encoding='utf-8') as f:
                    b = json.load(f)
                # batch.json 是可见、可手改、也可能被中断写坏的文件。
                # 下面每一步都必须容错：单个坏文件只跳过它自己，
                # **绝不能**让 _load_all 抛异常 —— 它跑在 DownloadManager()
                # 构造里、由 get_manager() 单例触发，一抛就是全部
                # /api/dl/* 端点 500、下载页整页不可用（实测：batch.json
                # 里缺一个 'id' 键即 KeyError；内容为 [] / "x" / 123 即
                # AttributeError，而原 except 只兜 OSError/ValueError）。
                if not isinstance(b, dict):
                    continue
                # 身份以**目录名**为准：缺 id、或 id 与目录不符（复制粘贴
                # 的产物）都不能造成批次错乱或整体丢失。
                b['id'] = name
                b.setdefault('name', name)
                b.setdefault('status', 'failed')
                b.setdefault('created', '')
                b.setdefault('dir', os.path.join(root, name))
                b.setdefault('error', '')
                if not isinstance(b.get('files'), list):
                    b['files'] = []
                if not isinstance(b.get('runs'), dict):
                    b['runs'] = {}
                if b.get('status') in ('downloading', 'resolving', 'converting'):
                    b['status'] = 'interrupted'
                    for fe in b['files']:
                        if isinstance(fe, dict) and fe.get('status') in (
                                'downloading', 'pending', 'verifying'):
                            fe['status'] = 'pending'
                b.setdefault('started', None)
                # 中断即已停机：让「耗时」停下，不要跨着重启一直涨
                b.setdefault('finished', time.time())
                self.batches[b['id']] = b
            except Exception:
                # 兜底：任何解析/结构异常都只影响这一个批次，不拖垮整个模块
                continue

    # ---- 创建 ----
    def create(self, name, accessions, concurrency=2, convert_sra=True,
               out_dir=None):
        """创建批次。out_dir：可选自定义输出位置（目录）；不填/留空则用
        默认 downloads/<批次id>/。指定目录即批次目录本身（不再嵌一层
        批次 id，用户指定哪里数据就在哪里）；经 check_path 校验——平台内
        或已挂载的输出根（write_roots）内，越界抛 ValueError。"""
        bid = f'{_sanitize(name)}_{time.strftime("%Y%m%d_%H%M%S")}'
        if str(out_dir or '').strip():
            bdir = check_path(str(out_dir).strip(), must_exist=False,
                              in_platform=True)
            os.makedirs(bdir, exist_ok=True)
            custom_dir = True
        else:
            bdir = check_path(os.path.join(_batch_root(), bid),
                              must_exist=False, in_platform=True)
            custom_dir = False
        b = {'id': bid, 'name': _sanitize(name), 'status': 'resolving',
             'created': time.strftime('%Y-%m-%d %H:%M'),
             'started': time.time(), 'finished': None,
             'concurrency': max(1, min(int(concurrency or 2), 4)),
             'convert_sra': bool(convert_sra), 'dir': bdir,
             'files': [], 'runs': {}, 'error': '',
             'converting': False, 'convert_msg': ''}
        with self.lock:
            self.batches[bid] = b
        self._save(b)
        self._blog(b, f'创建批次 {b["name"]}：共 {len(accessions)} 个条目，'
                      f'输出目录 {bdir}')
        if custom_dir and os.listdir(bdir):
            self._blog(b, '⚠ 指定目录非空：与本次数据同名的文件将被覆盖/续传')
        self._notify(bid)          # 让任务中心立刻出现这张卡

        def resolver():
            ok_n = 0
            n_effective = sum(1 for x in accessions
                              if str(x).strip() and not str(x).startswith('#'))
            n_done = 0
            for raw in accessions:
                raw = str(raw).strip()
                if not raw or raw.startswith('#'):
                    continue
                try:
                    if re.match(r'^https?://', raw, re.I):
                        # 完整数据文件 URL：主机白名单 + 文件名白名单
                        ok, why = check_url_host(raw)
                        if not ok:
                            raise ValueError(why)
                        fn = _safe_filename(raw)
                        out = check_path(bdir + os.sep + fn,
                                         must_exist=False, in_platform=True)
                        if not out.startswith(bdir + os.sep):
                            raise ValueError('输出文件名越界（已拦截）')
                        with self.lock:
                            acc = fn.split('.')[0] or raw[:24]
                            b['runs'][acc] = {'db': 'URL', 'organism': ''}
                            b['files'].append({
                                'acc': acc, 'url': raw, 'out': out,
                                'size': 0, 'done_bytes': 0, 'progress': 0.0,
                                'status': 'pending', 'md5': '',
                                'verified': False, 'error': '', 'db': 'URL',
                                'organism': ''})
                        self._blog(b, f'URL 条目: {raw[:80]}')
                        ok_n += 1
                        continue
                    acc = raw.upper()
                    info = resolve_accession(acc)
                    with self.lock:
                        b['runs'][acc] = {'db': info['db'],
                                          'organism': info.get('organism', '')}
                        for fe in info['files']:
                            # sracha 条目：无 ENA url，输出名用 acc.sra（sracha 落地）
                            if fe.get('sracha'):
                                fn = acc + '.sra'
                            else:
                                fn = _safe_filename(fe['url'])
                            out = check_path(bdir + os.sep + fn,
                                             must_exist=False,
                                             in_platform=True)
                            if not out.startswith(bdir + os.sep):
                                raise ValueError('输出文件名越界（已拦截）')
                            # size        = 展示/进度用总量，下载中会被
                            #               aria2c 实时值覆盖（可能按整数
                            #               MiB 缩写，仅供进度条，不可校验）
                            # expect_bytes = 校验用精确总量，只来自 ENA/
                            #               DDBJ 元数据，全程不被覆盖
                            b['files'].append({
                                'acc': acc, 'url': fe['url'], 'out': out,
                                'size': int(fe.get('bytes') or 0),
                                'expect_bytes': int(fe.get('bytes') or 0),
                                'done_bytes': 0, 'progress': 0.0,
                                'status': 'pending',
                                'md5': str(fe.get('md5', ''))[:64],
                                'verified': False, 'error': '',
                                'db': info['db'],
                                'organism': info.get('organism', '')})
                    self._blog(b, f'解析 {acc} -> {info["db"]} '
                                  f'{len(info["files"])} 个文件')
                    ok_n += 1
                except Exception as e:
                    with self.lock:
                        b['files'].append({
                            'acc': raw.upper()[:40], 'url': '', 'out': '',
                            'size': 0, 'done_bytes': 0, 'progress': 0.0,
                            'status': 'failed', 'md5': '', 'verified': False,
                            'error': f'解析失败: {e}', 'db': '?',
                            'organism': ''})
                    self._blog(b, f'解析失败 {raw}: {e}')
                n_done += 1
                # 落盘节流：原来每个条目都 self._save(b)，而 _save 要把整个
                # 批次（含全部已解析条目）重新 json.dump —— 500 条批次实测
                # 502 次落盘、约 12.5 万次条目写入、约 58MB 写盘（O(n²)）。
                # 每 20 条存一次 + 收尾必存，写入量回到 O(n)，进度不丢；
                # 内存态不受影响（前端看到的是内存快照，不是磁盘）。
                if n_done % 20 == 0 or n_done >= n_effective:
                    self._save(b)
            with self.lock:
                # cancel 可能落在解析期间：仅在仍是 resolving 时才转入下载
                if b['status'] == 'resolving':
                    b['status'] = 'downloading' if ok_n else 'failed'
                    if not ok_n:
                        b['error'] = '没有可下载的文件（全部解析失败）'
                        b['finished'] = time.time()
            if ok_n:
                self._blog(b, f'解析完成 {ok_n} 个条目，开始下载')
            else:
                self._blog(b, '解析全部失败')
            self._save(b)
            if ok_n and b['status'] == 'downloading':
                self._start_workers(bid)

        threading.Thread(target=resolver, daemon=True,
                         name=f'dl-resolve-{bid}').start()
        return bid

    # ---- 下载执行 ----
    def _start_workers(self, bid):
        threading.Thread(target=self._worker_loop, args=(bid,),
                         daemon=True, name=f'dl-{bid}').start()

    def _worker_loop(self, bid):
        with self.lock:
            b = self.batches.get(bid)
            conc = max(1, min(int((b or {}).get('concurrency', 2)), 4))
        threads = [threading.Thread(target=self._worker, args=(bid,),
                                    daemon=True) for _ in range(conc)]
        for t in threads:
            t.start()
        for t in threads:
            t.join()
        self._finalize(bid)

    def _next_pending(self, b):
        with self.lock:
            if b['status'] != 'downloading':
                return None
            for fe in b['files']:
                if fe['status'] == 'pending':
                    fe['status'] = 'downloading'
                    return fe
        return None

    def _worker(self, bid):
        while True:
            with self.lock:
                b = self.batches.get(bid)
                if not b or b['status'] != 'downloading':
                    return
            fe = self._next_pending(b)
            if not fe:
                return
            self._blog(b, f'开始下载 {fe["acc"]} ({fe["url"]})')
            try:
                self._download_file(b, fe)
            except Exception as e:                  # 单文件失败不拖垮批次
                with self.lock:
                    fe['status'] = 'failed'
                    fe['error'] = str(e)
                self._blog(b, f'失败 {fe["acc"]}: {e}')

    def _download_file(self, b, fe):
        try:
            out = check_path(fe['out'], must_exist=False, in_platform=True)
        except ValueError as e:
            with self.lock:
                fe['status'] = 'failed'
                fe['error'] = f'路径校验失败: {e}'
            self._blog(b, f'失败 {fe["acc"]}: {fe["error"]}')
            return
        # sracha 条目：无 ENA url，直接从 NCBI .sra 拉取（--prefer-ena 自动回落）
        if fe.get('sracha') or (not fe.get('url') and re.fullmatch(r'[SED]RR\d+', fe.get('acc') or '') and sracha_path()):
            os.makedirs(b['dir'], exist_ok=True)
            ok_dl = self._dl_sracha(b, fe)
            if ok_dl:
                return
            with self.lock:
                fe['status'] = 'failed'
                fe['error'] = fe.get('error') or 'sracha 拉取失败'
            self._blog(b, f'失败 {fe["acc"]}: {fe.get("error")}')
            return
        ok, why = check_url_host(fe['url'])
        if not ok:
            with self.lock:
                fe['status'] = 'failed'
                fe['error'] = f'URL 校验失败: {why}'
            self._blog(b, f'失败 {fe["acc"]}: {fe["error"]}')
            return
        if not out.startswith(b['dir'] + os.sep):   # 写盘前包含性检查
            with self.lock:
                fe['status'] = 'failed'
                fe['error'] = '输出路径越界（已拦截）'
            self._blog(b, f'失败 {fe["acc"]}: {fe["error"]}')
            return
        os.makedirs(b['dir'], exist_ok=True)
        a2 = aria2c_path()
        ok_dl = self._dl_aria2(a2, b, fe) if a2 else self._dl_python(fe)
        if not ok_dl:
            # 三级回退：INSDC run 编号（SRR/ERR/DRR）可用 sracha fetch
            # （--prefer-ena 优先 ENA FASTQ.GZ，失败自动回落 NCBI .sra，
            # 下载即断点续传 + MD5 校验；.sra 由收尾转换链自动转 FASTQ）。
            # GSA CRR / 任意 URL 不适用（NCBI 无此数据），维持内置 HTTP。
            if re.fullmatch(r'[SED]RR\d+', fe.get('acc') or '') \
                    and sracha_path():
                with self.lock:
                    fe['status'] = 'downloading'
                    fe['error'] = ''
                ok_dl = self._dl_sracha(b, fe)
        if not ok_dl:
            with self.lock:
                if fe['status'] != 'cancelled':
                    fe['status'] = 'failed'
                    fe['error'] = fe.get('error') or '下载失败'
            if fe['status'] == 'cancelled':
                self._blog(b, f'取消 {fe["acc"]}')
            else:
                self._blog(b, f'失败 {fe["acc"]}: {fe.get("error")}')
            return
        # 完整性：字节数比对（仅当元数据给出精确值时才做强校验）
        # ⚠️ 不能用下载中 aria2c 实时报的 size 做精确比对 —— aria2c 对
        #    ≥10MiB 的数按整数 MiB 缩写（36.84MiB 显示成 '36MiB'），
        #    误差可达 0.5MiB，用它卡字节数会把正常文件判成损坏。
        #    实测：DDBJ DRR003724.fastq.bz2 真实 38638145，aria2c 报
        #    37748736（=36.00MiB），差值 889409 → 误判 failed。
        #    expect_bytes 为空时退化为 ±1MiB 容差的粗校验。
        exp = int(fe.get('expect_bytes') or 0)
        try:
            actual = os.path.getsize(out)
        except OSError as e:
            with self.lock:
                fe['error'] = f'完整性检查出错: {e}'
            actual = 0
        if exp and actual != exp:
            with self.lock:
                fe['status'] = 'failed'
                fe['error'] = (f'字节数不一致（期望 {exp}, '
                               f'实际 {actual}）')
            self._blog(b, f'失败 {fe["acc"]}: {fe["error"]}')
            return
        if not exp and fe.get('size') and actual:
            if abs(actual - int(fe['size'])) > 1048576:
                with self.lock:
                    fe['status'] = 'failed'
                    fe['error'] = (f'字节数偏差过大（约 {fe["size"]}, '
                                   f'实际 {actual}）')
                self._blog(b, f'失败 {fe["acc"]}: {fe["error"]}')
                return
        verify_kind = 'bytes' if exp else 'size_approx'
        if out.lower().endswith('.gz') and not _gzip_ok(out):
            with self.lock:
                fe['status'] = 'failed'
                fe['error'] = 'gzip 文件损坏（魔数/抽样解压失败）'
            self._blog(b, f'失败 {fe["acc"]}: {fe["error"]}')
            return
        if out.lower().endswith('.bz2'):
            # DDBJ .fastq.bz2：无 md5/精确字节数可依，用 bz2 流式解压校验（读全量）
            ok_bz, why_bz = _bz2_ok(out)
            if not ok_bz:
                with self.lock:
                    fe['status'] = 'failed'
                    fe['error'] = f'bz2 文件损坏（流式解压失败: {why_bz}）'
                self._blog(b, f'失败 {fe["acc"]}: {fe["error"]}')
                return
            verify_kind = 'bz2_stream'
            self._blog(b, f'校验 {fe["acc"]}: bz2 流式解压通过')
        exp_md5 = str(fe.get('md5') or '')
        # 本地 md5 兜底（sracha 路径自己已做 MD5；.sra 的 md5 是 FASTQ 的，不适用）
        if (re.fullmatch(r'[0-9a-fA-F]{32}', exp_md5) and fe.get('url')
                and not out.lower().endswith('.sra')):
            ok_md5, note_md5 = _md5_ok(out, exp_md5)
            if ok_md5 is False:
                with self.lock:
                    fe['status'] = 'failed'
                    fe['error'] = f'MD5 校验失败: {note_md5}'
                self._blog(b, f'失败 {fe["acc"]}: {fe["error"]}')
                return
            if ok_md5 is None:
                # 大文件跳过/读取异常：如实标注，不阻断（WARN）
                with self.lock:
                    fe['md5_note'] = note_md5
                self._blog(b, f'警告 {fe["acc"]}: MD5 未校验（{note_md5}）')
            else:
                verify_kind = 'md5'
        with self.lock:
            fe['status'] = 'done'
            fe['verified'] = True
            fe['verify'] = verify_kind
            if not exp:
                # 元数据未给字节数：只做了 ±1MiB 粗校验，来源如实标注
                fe['size_source'] = 'unknown'
            fe['progress'] = 1.0
        self._blog(b, f'完成 {fe["acc"]}')
        self._save(b)

    def _dl_aria2(self, a2, b, fe):
        out = fe['out']
        extra = []
        if re.fullmatch(r'[0-9a-fA-F]{32}', str(fe.get('md5') or '')):
            # ENA 给了 md5 时让 aria2c 一并做完整性校验（有 .aria2 控制文件
            # 时按已记录的分块校验；本地收尾还会再算一遍 md5 兜底）
            extra = ['--check-integrity=true']
        cmd = [a2, '-x', '4', '-s', '4', '-k', '1M', '-c',
               '--file-allocation=none', '--console-log-level=warn',
               '--summary-interval=2', '--retry-wait=3', '--max-tries=5',
               *extra,
               '-d', b['dir'], '-o', os.path.basename(out), fe['url']]
        try:
            proc = subprocess.Popen(
                cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                text=True, encoding='utf-8', errors='replace')
        except OSError:
            return self._dl_python(fe)
        with self.lock:
            self.procs[fe['out']] = proc
        # aria2c 实测进度行（--summary-interval=2 每 2 秒一帧）：
        #   总量已知  [#289963 4.7MiB/0.9GiB(0%) CN:4 DL:1.5MiB ETA:9m53s]
        #   总量未知  [#289963 0B/0B CN:1 DL:0B]        ← 无百分比，size 保持 0
        # 解析要点（全部来自实测，见 handle 末尾注释）：
        #   · done_bytes 必须取 aria2c 自报的 DONE，**不能**用
        #     os.path.getsize(out)：-x4 -s4 -k1M 是分段落盘，aria2c 直接按
        #     偏移写，文件表观长度远早于真实进度就到满，实测 aria2c 自报 92%
        #     时 getsize 已经等于完整 8772679 → 进度条提前冲到 100% 再干等。
        #   · size 只在没有 ENA/DDBJ 精确元数据时才用 aria2c 的近似总量，
        #     绝不覆盖 expect_bytes；否则「已下载 / 总量」会倒挂。
        #   · progress 在下载过程中封顶 0.999：aria2c 的 DONE/TOTAL 按
        #     0.1MiB 四舍五入，最后一帧可能凑满而进程尚未收尾。
        #   · ETA 优先用 aria2c 的 ETA: 字段（落盘差分在分段下载下不可信）。
        pct_re = re.compile(r'\((\d+)%\)')
        sz_re = re.compile(r'\[#\w+\s+([\d.]+[KMGT]?i?B)/([\d.]+[KMGT]?i?B)')
        dl_re = re.compile(r'\bDL:\s*([\d.]+[KMGT]?i?B)')
        eta_re = re.compile(r'\bETA:\s*([0-9]+[hms][0-9hms]*)')
        try:
            for line in proc.stdout:
                m = pct_re.search(line)
                ms = sz_re.search(line)
                md = dl_re.search(line)
                me = eta_re.search(line)
                if not (m or ms or md or me):
                    continue
                with self.lock:
                    if ms:
                        aria_done = self._aria_bytes(ms.group(1))
                        aria_total = self._aria_bytes(ms.group(2))
                        fe['aria_done'] = aria_done
                        fe['aria_total'] = aria_total
                        if not fe.get('expect_bytes'):
                            fe['size'] = max(int(fe.get('size') or 0),
                                             aria_total, aria_done)
                        # 单调递增：每帧取历史最大值，进度条不回退
                        done = max(int(fe.get('done_bytes') or 0), aria_done)
                        fe['done_bytes'] = done
                        total = _known_total(fe)
                        if total > 0:
                            fe['progress'] = min(done / total, 0.999)
                        elif m:
                            fe['progress'] = min(int(m.group(1)), 100) / 100.0
                    elif m and not fe.get('done_bytes'):
                        fe['progress'] = min(int(m.group(1)), 100) / 100.0
                    if md:
                        fe['speed'] = self._aria_bytes(md.group(1))
                        self._prime_speed(out, fe.get('done_bytes') or 0)
                    elif fe.get('done_bytes'):
                        sp = self._tick_speed(out, fe['done_bytes'])
                        if sp:
                            fe['speed'] = sp
                    if me:
                        sec = self._aria_eta(me.group(1))
                        if sec is not None:
                            fe['eta'] = sec
        finally:
            rc = proc.wait()
            with self.lock:
                self.procs.pop(fe['out'], None)
                self._drop_speed(out)
                fe.pop('speed', None)
                fe.pop('eta', None)
                fe.pop('aria_done', None)
                fe.pop('aria_total', None)
        if rc != 0:
            with self.lock:
                fe['error'] = f'aria2c 退出码 {rc}'
            return False
        if not os.path.isfile(out) or os.path.getsize(out) == 0:
            with self.lock:
                fe['error'] = 'aria2c 完成但文件缺失/为空'
            return False
        with self.lock:
            # 收尾：进程已退出，此时 getsize 才是可信的真实字节数
            try:
                actual = os.path.getsize(out)
                fe['done_bytes'] = actual
                if not fe.get('expect_bytes'):
                    fe['size'] = max(int(fe.get('size') or 0), actual)
            except OSError:
                pass
        return True

    # sracha 断点续传文件：.<输出文件名>.sracha-progress，内容为 JSON
    #   {"url":…, "expected_size":N, "chunk_size":N, "total_chunks":N,
    #    "completed_chunks":[…], "expected_md5":"…"}
    # 这是 sracha 唯一可用的进度来源：目标文件在起手就被预分配到完整长度
    # （实测 175329159 字节从第 5 秒就是满值），且 indicatif 进度条在管道
    # （非 TTY）下不输出，所以 getsize 与 stdout 都拿不到真实进度。
    _SRACHA_PROG_SUFFIX = '.sracha-progress'

    def _sracha_progress(self, bdir, acc):
        """读 sracha 进度文件 → (done_bytes, total_bytes, md5) 或 None。

        归属判定同样必须精确：旧实现取前缀 '.'+acc，'.SRR1000' 会把
        SRR10001 的进度文件也算进来（startswith 前缀碰撞）。这里把
        '. <产物名> .sracha-progress' 中的产物名剥出来后按 _sracha_owned_name
        校验（产物名形如 SRR123.sra / SRR123.fastq.gz）。
        """
        cands = []
        try:
            for n in os.listdir(bdir):
                if not n.endswith(self._SRACHA_PROG_SUFFIX) or not n.startswith('.'):
                    continue
                stem = n[1:-len(self._SRACHA_PROG_SUFFIX)]
                # stem 是产物名（SRR123.sra / SRR123.fastq.gz）；也接受裸 acc
                # （startsWith 旧行为里 '.SRR123.sracha-progress' 这类精确形式）
                if stem == acc or _sracha_owned_name(stem, acc):
                    cands.append(n)
        except OSError:
            return None
        for n in cands:
            try:
                with open(os.path.join(bdir, n), encoding='utf-8') as f:
                    d = json.load(f)
            except (OSError, ValueError):
                continue
            try:
                exp = int(d.get('expected_size') or 0)
                csz = int(d.get('chunk_size') or 0)
                nch = len(d.get('completed_chunks') or [])
            except (TypeError, ValueError):
                continue
            if exp <= 0:
                continue
            done = min(csz * nch, exp) if csz else 0
            return done, exp, str(d.get('expected_md5') or '')
        return None

    def _dl_sracha(self, b, fe):
        """第三级回退：sracha fetch（仅 SRR/ERR/DRR）。

        --prefer-ena 优先下载 ENA 预生成的 FASTQ.GZ（与主路径同源同内容），
        ENA 无 FASTQ 时自动回落 NCBI 的 .sra（后续收尾链自动转 FASTQ）。
        自带断点续传与 MD5 校验。

        进度：sracha 在管道下不吐进度条，改用它的 .sracha-progress 续传
        文件轮询真实字节数（见 _sracha_progress），并同步刷新 size/speed。
        stderr 逐行转发到批次日志（-v）——重试与回退原因原先被整段吞掉。
        """
        exe = sracha_path()
        acc = fe['acc']
        bdir = b['dir']
        cmd = [exe, 'fetch', acc, '-O', bdir, '-f', '-v', '--no-progress',
               '--prefer-ena']
        try:
            proc = subprocess.Popen(
                cmd, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE,
                text=True, encoding='utf-8', errors='replace', bufsize=1)
        except OSError as e:
            with self.lock:
                fe['error'] = f'sracha 回退失败: {e}'
            return False
        self._blog(b, f'sracha: {" ".join(cmd[1:])}')
        tail = []
        stop = threading.Event()

        def _drain():
            """转发 sracha stderr（含重试告警），保留末 8 行用于失败报错。"""
            try:
                for line in proc.stderr:
                    line = re.sub(r'\x1b\[[0-9;]*m', '', line).rstrip()
                    if not line:
                        continue
                    tail.append(line)
                    if len(tail) > 8:
                        tail.pop(0)
                    self._blog(b, f'sracha[{acc}] {line}')
            except (OSError, ValueError):
                pass

        th = threading.Thread(target=_drain, daemon=True,
                              name=f'sracha-log-{acc}')
        th.start()
        deadline = time.time() + 24 * 3600
        sp_base = None            # (首次采样时刻, 首次采样已下载字节)
        while proc.poll() is None:
            if time.time() > deadline:
                proc.kill()
                with self.lock:
                    fe['error'] = 'sracha 超时（24 小时）'
                break
            pg = self._sracha_progress(bdir, acc)
            if pg:
                done, tot, md5 = pg
                now = time.time()
                with self.lock:
                    max_done = max(int(fe.get('done_bytes') or 0), done)
                    fe['size'] = max(int(fe.get('size') or 0), tot)
                    fe['done_bytes'] = max_done
                    if tot > 0:
                        fe['progress'] = min(done / tot, 0.999)
                    # 速率用「下载窗口内的累计均值」，不用逐块差分：
                    # 8MiB 分块完成是脉冲式的（实测同一秒内可能连落数块），
                    # 差分 EMA 会报出 14 MB/s 这种远高于真实吞吐的尖峰，
                    # 累计均值稳定且让 ETA 可信。
                    if max_done > 0:
                        if sp_base is None:
                            sp_base = (now, max_done)
                            fe.pop('speed', None)
                        elif now - sp_base[0] >= 0.5:
                            sp = (max_done - sp_base[1]) / (now - sp_base[0])
                            if sp > 0:
                                fe['speed'] = sp
                                rem = max(tot - max_done, 0)
                                if sp > 1024:
                                    fe['eta'] = int(rem / sp)
                    if md5 and not fe.get('md5'):
                        fe['md5'] = md5
            stop.wait(1.2)
        rc = proc.wait()
        stop.set()
        th.join(timeout=2)
        with self.lock:
            self._drop_speed(fe['out'])
            fe.pop('speed', None)
        if rc != 0:
            with self.lock:
                fe['error'] = (f'sracha 退出码 {rc}'
                               + ('；' + tail[-1] if tail else ''))
            return False
        # 识别产物：acc.sra（NCBI 路径）或 acc[_N].fastq.gz（ENA 路径）
        # 归属必须用 _sracha_owned_name 的精确形式，不用 startswith(acc)——
        # 否则 SRR1000 会把 SRR10001 的产物登记到自己名下
        names = set(os.listdir(bdir))
        sra_p = bdir + os.sep + acc + '.sra'
        fastqs = sorted(
            n for n in names
            if _sracha_owned_name(n, acc) and n.endswith(('.fastq.gz', '.fq.gz')))
        with self.lock:
            # ENA 路径已经产出 FASTQ（或 sracha 把 .sra 与 .sracha-progress
            # 之外的产物直接落好）时优先登记 FASTQ：原先先判 .sra，会把
            # 已经下好的 FASTQ.GZ 丢掉，只留一个还得再转换的 .sra。
            if fastqs:
                produced_paths = []
                for n in fastqs:
                    produced_paths.append(
                        check_path(bdir + os.sep + n, must_exist=True,
                                   in_platform=True))
                new_entries = []
                for e in b['files']:
                    if e['out'] == fe['out']:
                        for pp in produced_paths:
                            new_entries.append({
                                'acc': acc, 'url': '', 'out': pp,
                                'size': os.path.getsize(pp),
                                'expect_bytes': os.path.getsize(pp),
                                'done_bytes': os.path.getsize(pp),
                                'progress': 1.0, 'status': 'done',
                                'md5': fe.get('md5', ''), 'verified': True,
                                'error': '', 'db': fe.get('db', 'ENA'),
                                'organism': fe.get('organism', '')})
                    else:
                        new_entries.append(e)
                b['files'] = new_entries
                return True
            if os.path.isfile(sra_p):
                # 交给收尾转换链（条目 out 指向 .sra，_finalize 会触发转换）
                fe['out'] = sra_p
                fe['size'] = os.path.getsize(sra_p)
                fe['expect_bytes'] = fe['size']
                fe['done_bytes'] = fe['size']
                fe['status'] = 'done'
                fe['verified'] = True                # sracha 已做 MD5 校验
                fe['progress'] = 1.0
                return True
        with self.lock:
            fe['error'] = 'sracha 完成但未找到产物文件'
        return False

    def _dl_python(self, fe):
        """内置 HTTP 下载回退（Range 断点续传；URL 已在外层校验）。

        先写 .aria_part 临时文件（无 .gz 后缀 → safe_open 原样二进制写，
        不会对已压缩数据二次 gzip），完成后改名到位。"""
        out, url = fe['out'], fe['url']
        part = check_path(out + '.aria_part', must_exist=False,
                          in_platform=True)
        pos = os.path.getsize(part) if os.path.isfile(part) else 0
        headers = {'User-Agent': _UA}
        if pos:
            headers['Range'] = f'bytes={pos}-'
        try:
            req = urllib.request.Request(url, headers=headers)
            with _OPENER.open(req, timeout=60) as r:
                total = int(r.headers.get('content-length') or 0) + \
                    (pos if r.status == 206 else 0)
                if r.status != 206:
                    pos = 0
                with self.lock:
                    # 元数据已有精确总量时不覆盖（同 _dl_aria2：size 是
                    # 前端展示分母，被近似值覆盖会让「已下载/总量」倒挂）
                    if not fe.get('expect_bytes'):
                        fe['size'] = total or fe.get('size') or 0
                    denom = _known_total(fe) or total
                # 断点续传必须**追加**：'wb' 会把已下载的前缀截断，于是续传
                # 产出一个"只有尾巴"的坏文件，而 done_bytes 仍按累计字节数
                # 报成完整长度（进度条 100%、文件却不完整）。expect_bytes
                # 有值时会被完整性检查兜住判 failed；无值时就是静默交付坏文件。
                mode = 'ab' if (pos and r.status == 206) else 'wb'
                with safe_open(part, mode) as f:
                    while True:
                        chunk = r.read(1 << 20)
                        if not chunk:
                            break
                        f.write(chunk)
                        pos += len(chunk)
                        with self.lock:
                            # pos 是真实读到的字节数（非文件表观长度），可信
                            fe['done_bytes'] = pos
                            fe['progress'] = (min(pos / denom, 0.999)
                                              if denom else 0.0)
                            sp = self._tick_speed(out, pos)
                            if sp:
                                fe['speed'] = sp
            shutil.move(part, out)
            return True
        except urllib.error.HTTPError as e:
            # 416 = Range 起点超出文件尾 → 本地 .aria_part 其实已经完整，
            # 直接改名到位即可（这是"进程在改名之前挂掉"的续传收尾路径）。
            # 必须写在 except 里：urllib 的 opener 对非 2xx/3xx 是**先抛
            # HTTPError**、根本不会把响应对象交出来，所以原先写在 with 体内
            # 的 `if r.status == 416` 是**死代码**，416 会一路掉到下面那个
            # URLError 分支被当成下载失败 —— 完整文件永远收不了尾。
            if e.code == 416 and os.path.isfile(part):
                shutil.move(part, out)
                return True
            with self.lock:
                fe['error'] = f'HTTP 下载失败: {e}'
            return False
        except (urllib.error.URLError, OSError) as e:
            with self.lock:
                fe['error'] = f'HTTP 下载失败: {e}'
            return False
        finally:
            with self.lock:
                self._drop_speed(out)
                fe.pop('speed', None)

    # ---- 收尾：.sra 转换 ----
    def _finalize(self, bid):
        with self.lock:
            b = self.batches.get(bid)
            if not b:
                return
            done = sum(1 for f in b['files'] if f['status'] == 'done')
            failed = sum(1 for f in b['files'] if f['status'] == 'failed')
            b['done'], b['failed'] = done, failed
            sra_left = [f for f in b['files']
                        if f['status'] == 'done'
                        and f['out'].lower().endswith('.sra')]
            bz2_left = [f for f in b['files']
                        if f['status'] == 'done'
                        and f['out'].lower().endswith('.fastq.bz2')]
            engine, _path = sra_convert_engine()
            b['status'] = ('converting'
                           if ((sra_left and b.get('convert_sra') and engine)
                               or bz2_left)
                           else ('completed' if done else 'failed'))
            if b['status'] != 'converting':
                b['finished'] = b.get('finished') or time.time()
        if b['status'] == 'converting':
            self._blog(b, f'批次转换中（{len(sra_left)} 个 .sra 待转换）')
        elif b['status'] == 'completed':
            self._blog(b, f'批次完成 {b.get("done", 0)}/'
                          f'{len(b.get("files", []))} 成功')
        self._save(b)
        if b['status'] == 'converting':
            self._convert_sras(b)

    def _run_decode_cmd(self, cmd, b, timeout):
        """跑 .sra 解码子进程，并登记进 self.procs 以便 cancel() 终止它。

        原来用 subprocess.run：进程没人管 → 用户在收尾阶段点「停止」时
        既杀不掉进程（继续吃 CPU 解码 GB 级 .sra），批次状态还会在
        _convert_sras 末尾被覆写成 completed（取消被静默丢弃）。
        键放在**批次目录内**，cancel() 的前缀匹配才认得出它。
        """
        key = os.path.join(b['dir'], '.decode-proc')
        proc = subprocess.Popen(
            cmd, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE,
            text=True, encoding='utf-8', errors='replace')
        with self.lock:
            self.procs[key] = proc
        try:
            _out, err = proc.communicate(timeout=timeout)
        except subprocess.TimeoutExpired:
            proc.kill()
            try:
                proc.communicate(timeout=10)
            except Exception:
                pass
            raise RuntimeError('解码超时')
        finally:
            with self.lock:
                self.procs.pop(key, None)
        if proc.returncode != 0:
            tail = (err or '').strip().splitlines()
            raise RuntimeError(f'解码退出码 {proc.returncode}'
                               + (f'：{tail[-1][:200]}' if tail else ''))
        return proc.returncode

    def _convert_sras(self, b):
        engine, exe = sra_convert_engine()
        with self.lock:
            if not b or not engine:
                if b:
                    # 同理不能覆写取消态
                    if b.get('status') != 'cancelled':
                        b['status'] = 'completed'
                    b['convert_msg'] = ('无 .sra 需转换' if engine is None and b
                                        else b.get('convert_msg', ''))
                    self._save(b)
                    self._blog(b, f'批次完成 {b.get("done", 0)}/'
                                  f'{len(b.get("files", []))} 成功')
                return
            b['converting'] = True
            sras = [f['out'] for f in b['files']
                    if f['status'] == 'done'
                    and f['out'].lower().endswith('.sra')]
            threads = max(2, min((os.cpu_count() or 4) // 2, 8))
        msgs = []
        failed_conv = 0
        for sra in sras:
            # 收尾阶段也必须能取消：用户在大 .sra 解码中点了停止，
            # 剩余文件不该继续解码，状态也不该在末尾被覆写成 completed。
            if b.get('status') == 'cancelled':
                self._blog(b, '已取消：停止后续 .sra 解码')
                break
            try:
                if engine == 'sracha':
                    # 单文件静态引擎：双端输出 _1/_2（split-3 命名，与
                    # fasterq-dump/ENA 一致），单端输出 <acc>.fastq.gz。
                    # ⚠ 只匹配 `base + '_'` 会把单端产物判成「无输出」：
                    #   sracha 对单端 run 明确写道 "run is single-end …
                    #   output is a single file"，实测写出
                    #   SRR39909446.fastq.gz，而转换链只 glob 下划线名 →
                    #   抛错并在替换条目前中断，批次假报「已完成」而
                    #   产物被丢掉（ready_files 为空，无法转入流程）。
                    self._run_decode_cmd(
                        [exe, 'fastq', sra, '-O', b['dir'],
                         '-t', str(threads), '-f', '-q'], b, 8 * 3600)
                    produced = _sra_fastq_outputs(b['dir'],
                                                  os.path.basename(sra)[:-4])
                    if not produced:
                        raise RuntimeError('sracha 无输出')
                    os.remove(sra)
                    msgs.append(f'✔ {os.path.basename(sra)} → '
                                f'{len(produced)} FASTQ.GZ（sracha）')
                else:
                    self._run_decode_cmd(
                        [exe, '--split-files', '-e', str(threads), '-O',
                         b['dir'], sra], b, 4 * 3600)
                    produced = _sra_fastq_outputs(b['dir'],
                                                  os.path.basename(sra)[:-4],
                                                  exts=('.fastq', '.fq'))
                    if not produced:
                        raise RuntimeError('fasterq-dump 无输出')
                    for p in produced:
                        src = check_path(b['dir'] + os.sep + p,
                                         must_exist=True, in_platform=True)
                        with open(src, 'rb') as fin, \
                                safe_open(src + '.gz', 'wb', compress_level=6) \
                                as fout:
                            shutil.copyfileobj(fin, fout, 1 << 20)
                        os.remove(src)
                    produced = [p + '.gz' for p in produced]
                    os.remove(sra)
                    msgs.append(f'✔ {os.path.basename(sra)} → '
                                f'{len(produced)} FASTQ.GZ（fasterq-dump）')
                # 把 .sra 条目替换为产物 FASTQ 条目（供 ready_files 配对进流程）
                produced_paths = [check_path(b['dir'] + os.sep + p,
                                             must_exist=True,
                                             in_platform=True)
                                  for p in produced]
                with self.lock:
                    new_entries = []
                    for fe in b['files']:
                        if fe['out'] == sra:
                            acc = fe['acc']
                            for pp in produced_paths:
                                new_entries.append({
                                    'acc': acc, 'url': '', 'out': pp,
                                    'size': os.path.getsize(pp),
                                    'expect_bytes': os.path.getsize(pp),
                                    'done_bytes': os.path.getsize(pp),
                                    'progress': 1.0, 'status': 'done',
                                    'md5': '', 'verified': True,
                                    'error': '', 'db': fe.get('db', 'URL'),
                                    'organism': fe.get('organism', '')})
                        else:
                            new_entries.append(fe)
                    b['files'] = new_entries
            except Exception as e:
                if b.get('status') == 'cancelled':
                    # 是被 cancel() 杀掉的那一个，不算转换失败
                    self._blog(b, f'已取消 .sra 转换 {os.path.basename(sra)}')
                    break
                failed_conv += 1
                msgs.append(f'✘ {os.path.basename(sra)}: {e}')
                self._blog(b, f'✘ .sra 转换失败 {os.path.basename(sra)}: {e}')
        # DDBJ .fastq.bz2 条目（与 .sra 转换互不影响）
        bz2_entries = [f for f in b['files']
                       if f['status'] == 'done'
                       and f['out'].lower().endswith('.fastq.bz2')]
        if bz2_entries:
            self._blog(b, f'批次转换中（{len(bz2_entries)} 个 .fastq.bz2 待转换）')
            self._convert_bz2(b, bz2_entries)
        with self.lock:
            b['converting'] = False
            done = sum(1 for f in b['files'] if f['status'] == 'done')
            failed = sum(1 for f in b['files'] if f['status'] == 'failed')
            b['done'], b['failed'] = done, failed
            b['convert_failed'] = failed_conv
            if failed_conv:
                # 转换失败时批次仍是 completed（数据确实下完了），但必须
                # 显式带上失败数：只有 .sra 而没有 FASTQ 的批次是无法
                # 转入分析流程的，不能再显示成一片绿。
                msgs.insert(0, f'⚠ {failed_conv} 个 .sra 未转成 FASTQ'
                               f'（可点「重试转换」）')
            b['convert_msg'] = '；'.join(msgs) or '无 .sra'
            has_fastq = any(
                re.search(r'\.f(ast)?q\.gz$', f['out'] or '', re.I)
                for f in b['files'] if f['status'] == 'done')
            if failed_conv and not has_fastq:
                b['error'] = (b.get('error') or '') or (
                    f'{failed_conv} 个 .sra 转换失败，本批次暂无可用 FASTQ')
            if b.get('status') == 'cancelled':
                # 用户在收尾阶段点了停止：**不能**把它覆写成 completed，
                # 否则界面显示"已完成"而实际只解了一部分。
                b['convert_msg'] = ('；'.join(msgs + ['已取消（未完成的 .sra '
                                                     '可点「重试转换」）']))
            else:
                b['status'] = ('completed' if done or not b['files']
                               else 'failed')
            b['finished'] = b.get('finished') or time.time()
        self._blog(b, f'批次完成 {b.get("done", 0)}/'
                      f'{len(b.get("files", []))} 成功'
                      + (f'（{failed_conv} 个 .sra 转换失败）'
                         if failed_conv else ''))
        self._save(b)

    def retry_convert(self, bid):
        """重跑 .sra → FASTQ 转换（不动已下载的数据，供转换失败后补救）。

        为什么需要：`.sra` 下载成功后转换是独立一步，失败时批次已落
        `completed`，此前没有重试入口，用户只能整批重下。
        """
        with self.lock:
            b = self.batches.get(bid)
            if not b or b.get('status') in self._ACTIVE:
                return 0
            sra_left = [f for f in b['files']
                        if f['status'] == 'done'
                        and (f.get('out') or '').lower().endswith('.sra')]
            if not sra_left:
                return 0
            b['status'] = 'converting'
            b['converting'] = True
            b['convert_msg'] = ''
            b['error'] = ''
            b['finished'] = None        # 重新在跑：耗时继续累计
        self._blog(b, f'重试 .sra 转换（{len(sra_left)} 个）')
        self._save(b)
        threading.Thread(target=self._convert_sras, args=(b,), daemon=True,
                         name=f'dl-convert-{bid}').start()
        return len(sra_left)

    def _convert_bz2(self, b, entries):
        """DDBJ .fastq.bz2 → .fastq.gz（流式转压缩，gzip level 1）。
        单端单文件 → <acc>.fastq.gz；双端 _1/_2 → 同名 _1/_2.fastq.gz。
        完成后条目替换为 FASTQ（供 ready_files 配对进流程）。"""
        import bz2 as _bz2
        for fe in entries:
            if b.get('status') == 'cancelled':      # 同 .sra：收尾也要能停
                self._blog(b, '已取消：停止后续 .fastq.bz2 转换')
                break
            src = fe['out']
            try:
                base = os.path.basename(src)[:-len('.fastq.bz2')]
                outs = []
                with _bz2.open(src, 'rb') as fin:
                    if base.endswith('_1') or base.endswith('_2'):
                        dst = src[:-len('.fastq.bz2')] + '.fastq.gz'
                        with safe_open(dst, 'wb', compress_level=1) as fout:
                            shutil.copyfileobj(fin, fout, 1 << 20)
                        outs.append(dst)
                    else:
                        # 单端：直接转 gz，命名与 ENA/sracha 产物一致
                        dst = os.path.join(os.path.dirname(src),
                                           base + '.fastq.gz')
                        with safe_open(dst, 'wb', compress_level=1) as fout:
                            shutil.copyfileobj(fin, fout, 1 << 20)
                        outs.append(dst)
                os.remove(src)
                with self.lock:
                    new_entries = []
                    for e in b['files']:
                        if e is fe:
                            for pp in outs:
                                new_entries.append({
                                    'acc': e['acc'], 'url': '', 'out': pp,
                                    'size': os.path.getsize(pp),
                                    'done_bytes': os.path.getsize(pp),
                                    'progress': 1.0, 'status': 'done',
                                    'md5': '', 'verified': True,
                                    'error': '', 'db': e.get('db', 'DDBJ'),
                                    'organism': e.get('organism', '')})
                        else:
                            new_entries.append(e)
                    b['files'] = new_entries
                self._blog(b, f'✔ {os.path.basename(src)} → FASTQ.GZ（DDBJ）')
            except Exception as e:
                with self.lock:
                    fe['status'] = 'failed'
                    fe['error'] = f'bz2 转换失败: {e}'
                self._blog(b, f'✘ {os.path.basename(src)}: {e}')

    # ---- 状态 / 控制 ----
    # 逐文件明细里对前端可见的字段
    _SNAP_KEYS = ('acc', 'url', 'out', 'size', 'expect_bytes',
                  'done_bytes', 'progress', 'status', 'md5', 'verified',
                  'error', 'db', 'organism', 'speed', 'eta')

    def snapshot(self, bid, include_files=True):
        """批次快照。

        include_files=False 时聚合值照算，但不物化逐文件明细 ——
        list_snapshots() 要对最近 40 个批次各取一次快照，而它紧接着就把
        files 置空；40×500 条会把 2 万个 dict 造出来又丢掉
        （实测 32 ms/次轮询，下载页每 2 秒轮询一次）。
        """
        with self.lock:
            b = self.batches.get(bid)
            if not b:
                return None
            # 一律走 .get + 兜底：snapshot() 在请求热路径上（/api/dl/batches
            # 与 /api/dl/batch/<bid>），缺字段只应降级显示，不该 500。
            raw = b.get('files') or []   # 聚合直接用原始条目，不先物化
            total = len(raw)
            done = sum(1 for f in raw if f.get('status') == 'done')
            failed = sum(1 for f in raw if f.get('status') == 'failed')
            # 批次进度按字节加权：各文件体积可能差几个数量级，按条目数取
            # 平均会让大文件的进度被小文件稀释（100 个文件里 1 个大文件
            # 占 90% 数据量时，条目均值几乎不反映真实完成度）。
            # 没有总量的条目（如 sracha 尚未报出 expected_size）用已知
            # 总量的均值兜底；全部未知时退回条目均值。
            weights = [_fe_weight(f) for f in raw]
            known = [w for w in weights if w > 0]
            if total and known:
                default = sum(known) / len(known)
                denom = sum(w if w > 0 else default for w in weights)
                overall = (sum((f.get('progress') or 0)
                               * (weights[i] if weights[i] > 0 else default)
                               for i, f in enumerate(raw)) / denom
                           if denom else 0.0)
            else:
                overall = (sum((f.get('progress') or 0) for f in raw) / total
                           if total else 0.0)
            # 批次速率 = 正在下载条目的速率之和
            active = ('pending', 'downloading', 'verifying')
            speed = sum(float(f.get('speed') or 0) for f in raw
                        if f.get('status') == 'downloading')
            done_bytes = sum(int(f.get('done_bytes') or 0) for f in raw)
            # size 是前端展示分母：精确元数据优先，缺失时退回下载器近似值。
            # 仍取 max(size, done) 兜底，杜绝「已下载 > 总量」的倒挂。
            total_bytes = sum(max(_known_total(f),
                                  int(f.get('done_bytes') or 0))
                              for f in raw)
            remain = sum(max(_known_total(f)
                             - int(f.get('done_bytes') or 0), 0)
                         for f in raw if f.get('status') in active)
            # ETA：按剩余字节 / 当前总速率推算；速率不可用时退回各条目
            # 自报的 ETA 最大值（aria2c 的 ETA: 字段）。两者都没有才为 0，
            # 前端此时显示「计算中」——不会像以前那样因为 remain 被
            # 表观文件长度提前清零而永远显示「计算中」。
            eta = int(remain / speed) if speed > 1024 else 0
            if not eta:
                fe_etas = [int(f.get('eta') or 0) for f in raw
                           if f.get('status') == 'downloading']
                eta = max(fe_etas) if fe_etas else 0
            runs = b.get('runs') or {}
            files = None
            if include_files:
                files = [{k: f.get(k) for k in self._SNAP_KEYS} for f in raw
                         if isinstance(f, dict)]
            started = b.get('started')
            finished = b.get('finished')
            elapsed = None
            if started:
                try:
                    elapsed = max(0.0, (finished or time.time())
                                  - float(started))
                except (TypeError, ValueError):
                    elapsed = None
            return {'id': b.get('id') or bid, 'name': b.get('name') or bid,
                    'status': b.get('status') or 'failed',
                    'dir': b.get('dir') or '',
                    'created': b.get('created') or '',
                    'error': b.get('error') or '',
                    'started': started, 'finished': finished,
                    'elapsed': (round(elapsed, 1)
                                if elapsed is not None else None),
                    'concurrency': b.get('concurrency', 2),
                    'converting': b.get('converting', False),
                    'convert_msg': b.get('convert_msg', ''),
                    'convert_failed': int(b.get('convert_failed') or 0),
                    'n_files': total, 'done': done, 'failed': failed,
                    'n_runs': len(runs),
                    'speed': round(speed, 1), 'eta': eta,
                    'done_bytes': done_bytes, 'total_bytes': total_bytes,
                    'overall': round(min(overall, 1.0), 4), 'files': files}

    # 日志尾缓存：{(path, size, mtime_ns, lines): (行...)}
    # 与 TaskManager._tail_log_file 同一套做法（那边注释里记着同样的教训）。
    _logcache = {}

    @classmethod
    def _tail_log_file(cls, path, lines):
        """从 batch.log 尾部读 lines 行：seek 到尾部 + 结果缓存。

        为什么不能 `f.read().splitlines()[-n:]` 整文件读：
        list_snapshots() 会对最近 **40** 个批次各取一次日志尾，而下载页在
        任意批次进行中时**每 2 秒**轮询 /api/dl/batches。实测
        40 × 0.2MB = 47~59 ms、40 × 2MB = 383~442 ms、
        40 × 8MB = 1213~1236 ms（每次都重读数百 MB 磁盘，冷热无差别）——
        轮询本身成了负担。这里只读尾部必要字节，并以
        (size, mtime_ns, lines) 作键缓存，文件未变时零磁盘 IO。
        """
        try:
            st = os.stat(path)
        except OSError:
            return None
        key = (path, st.st_size, st.st_mtime_ns, lines)
        hit = cls._logcache.get(key)
        if hit is not None:
            return list(hit)
        try:
            size = st.st_size
            # 最小 256KB：单行日志可能很长（sracha -v 的报错/重试行会整行转发），
            # 块太小会把最后一行截断丢掉。256KB 对磁盘是噪声级开销。
            chunk = min(size, max(1 << 18, lines * 512))
            with safe_open(path, 'rb') as f:
                if size > chunk:
                    f.seek(size - chunk)
                data = f.read()
        except OSError:
            return None
        rows = data.decode('utf-8', errors='replace').splitlines()
        if size > chunk and rows:
            rows = rows[1:]              # 首行可能被字节截断，丢弃
        tail = tuple(rows[-lines:])
        if len(cls._logcache) > 64:      # 40 个批次 + 余量，防无限增长
            cls._logcache.pop(next(iter(cls._logcache)))
        cls._logcache[key] = tail
        return list(tail)

    def _read_log_tail(self, b, n=40):
        """读取批次过程日志 batch.log 的最后 n 行（读失败返回空列表）。"""
        d = b.get('dir') or ''
        if not d:
            return []
        return self._tail_log_file(os.path.join(d, 'batch.log'), n) or []

    def list_snapshots(self):
        with self.lock:
            ids = list(self.batches.keys())
        outs = []
        for bid in sorted(ids, reverse=True)[:40]:
            # include_files=False：本接口不下发逐文件明细（下载页点开 📁
            # 才单独拉 /api/dl/batch/<bid>），没必要物化 40×N 个 dict。
            s = self.snapshot(bid, include_files=False)
            if s:
                s['log_tail'] = self._read_log_tail(self.batches.get(bid) or {})
                outs.append(s)
        return outs

    def cancel(self, bid):
        """停止批次。已到终态的批次只清理残留子进程，**不改状态**。

        返回 True 表示「确实把在跑的批次停下了」；False 表示批次不存在或
        已经是终态（没什么可停的）。任务中心的「停止」与全局重置都靠这个
        返回值判断是否真的受理了。
        """
        with self.lock:
            b = self.batches.get(bid)
            if not b:
                return False
            bdir = b['dir']
            for out, proc in list(self.procs.items()):
                if os.path.normpath(out).startswith(
                        os.path.normpath(bdir) + os.sep):
                    try:
                        proc.terminate()
                    except OSError:
                        pass
            if b.get('status') in _TERMINAL_BATCH_STATES:
                return False
            for fe in b['files']:
                if fe['status'] in ('pending', 'downloading', 'verifying'):
                    fe['status'] = 'cancelled'
            b['status'] = 'cancelled'
            b['finished'] = time.time()
        self._blog(b, '批次已取消')
        self._save(b)
        return True

    def retry_failed(self, bid):
        """重试失败条目；批次处于 interrupted 时兼作「续传」。

        为什么要把 interrupted 纳进来：服务重启会把 downloading/resolving/
        converting 的批次标成 interrupted，其条目退回 pending（failed 为 0）。
        前端原本只在 failed>0 时给重试按钮，于是重启后的批次在界面上
        **没有任何出路**，只能删了重下。
        """
        n = 0
        with self.lock:
            b = self.batches.get(bid)
            if not b:
                return 0
            resume = b.get('status') == 'interrupted'
            for fe in b['files']:
                st = fe.get('status')
                if not (st == 'failed' or (resume and st == 'pending')):
                    continue
                if not _retryable(fe):
                    continue
                fe['status'] = 'pending'
                fe['error'] = ''
                fe['progress'] = 0.0
                # done_bytes 必须一起清零：进度分子现在就是它，
                # 只清 progress 会让重试瞬间显示上一次的高进度
                # （aria2c 续传后第一帧会报出真实 DONE，随即校正）。
                fe['done_bytes'] = 0
                fe.pop('speed', None)
                fe.pop('eta', None)
                n += 1
            if n:
                b['status'] = 'downloading'
                b['finished'] = None    # 重新在跑：耗时继续累计
        self._save(b)
        if n:
            self._blog(b, f'{"续传" if resume else "重试"} {n} 个文件')
            self._start_workers(bid)
        return n

    _RECORD_FILES = ('batch.json', 'batch.log')
    _ACTIVE = ('resolving', 'downloading', 'converting')

    def delete(self, bid):
        """全删：批次记录 + 日志 + 已下载文件（目录在创建时已校验）。"""
        self.cancel(bid)
        with self.lock:
            b = self.batches.pop(bid, None)
        if not b:
            return False
        if b.get('dir'):
            shutil.rmtree(b['dir'], ignore_errors=True)
        return True

    def delete_files(self, bid):
        """只删文件：清掉批次目录内的数据文件，保留记录(batch.json)与
        日志(batch.log)；文件条目标记 deleted（列表仍可见、可追溯）。"""
        with self.lock:
            b = self.batches.get(bid)
            if not b:
                return False
            if b.get('status') in self._ACTIVE:
                active = True
            else:
                active = False
        if active:
            self.cancel(bid)          # 运行中先停（cancel 会改状态并保存）
        bdir = b.get('dir') or ''
        n = 0
        if bdir and os.path.isdir(bdir):
            for name in os.listdir(bdir):
                if name in self._RECORD_FILES:
                    continue
                p_ = os.path.join(bdir, name)
                try:
                    if os.path.isdir(p_) and not os.path.islink(p_):
                        shutil.rmtree(p_, ignore_errors=True)
                    else:
                        os.remove(p_)
                    n += 1
                except OSError:
                    pass
        with self.lock:
            b = self.batches.get(bid)
            if b:
                for fe in b['files']:
                    if fe.get('status') in ('done', 'failed', 'cancelled',
                                            'pending'):
                        fe['status'] = 'deleted'
                        fe['progress'] = 0.0
                        # 文件已从磁盘删掉，done_bytes 不能再留在快照里
                        # （批次「已下载 X / 总量」会继续显示这些字节）
                        fe['done_bytes'] = 0
                        fe['size'] = 0
                        fe['expect_bytes'] = 0
                        fe.pop('speed', None)
                        fe.pop('eta', None)
                b['files_deleted'] = True
                self._save(b)
        self._blog(b, f'已删除 {n} 个数据文件（保留批次记录与日志）')
        return True

    def delete_record(self, bid):
        """只删记录：移除 batch.json/batch.log（批次从列表消失），已下载
        数据文件保留在 downloads/<批次>/ 磁盘原处。"""
        with self.lock:
            b = self.batches.get(bid)
            if not b:
                return False
            active = b.get('status') in self._ACTIVE
        if active:
            self.cancel(bid)
        with self.lock:
            b = self.batches.pop(bid, None)
        if not b:
            return False
        if b.get('dir'):
            for name in self._RECORD_FILES:
                try:
                    os.remove(os.path.join(b['dir'], name))
                except OSError:
                    pass
        return True

    # ---- 已完成文件的清单（进分析流程用） ----
    def ready_files(self, bid):
        snap = self.snapshot(bid)
        if not snap:
            return {}
        by_acc = {}
        for f in snap['files']:
            if f['status'] != 'done' or not f['out']:
                continue
            fn = os.path.basename(f['out'])
            if re.search(r'_1\.f(ast)?q\.gz$', fn, re.I):
                by_acc.setdefault(f['acc'], {})['r1'] = f['out']
            elif re.search(r'_2\.f(ast)?q\.gz$', fn, re.I):
                by_acc.setdefault(f['acc'], {})['r2'] = f['out']
            elif re.search(r'\.f(ast)?q\.gz$', fn, re.I):
                by_acc.setdefault(f['acc'], {}).setdefault('single', f['out'])
        return by_acc


_manager = None


def get_manager():
    global _manager
    if _manager is None:
        _manager = DownloadManager()
    return _manager
