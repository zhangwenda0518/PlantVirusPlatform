# -*- coding: utf-8 -*-
"""宿主基因组下载（host_genome）纯函数单测：gget/Ensembl 源接入层。

只测不联网的部分：物种名规范化、gget.ref 返回值解析、HTTP 域名白名单、
日志回调适配、genome.fa 锚点落盘。真实下载链路由 /api/host_genome_download
端点巡检（空参数 400）与人工联调覆盖。
"""
import os
import sys

import pytest

from Virus_Platform_Core.public_meta import host_genome as hg


CANNED_REF = {
    'oryza_sativa': {
        'genome_dna': {
            'ftp': 'http://ftp.ensemblgenomes.org/pub/release-60/plants/'
                   'fasta/oryza_sativa/dna/'
                   'Oryza_sativa.IRGSP-1.0.dna.toplevel.fa.gz',
            'ensembl_release': 60,
            'release_date': '2025-06-02',
            'bytes': '400.2 MB',
        },
        'annotation_gtf': {
            'ftp': 'http://ftp.ensemblgenomes.org/pub/release-60/plants/'
                   'gtf/oryza_sativa/Oryza_sativa.IRGSP-1.0.60.gtf.gz',
        },
    },
}


# ── normalize_species：空格/大小写 → Ensembl 命名 ─────────────
def test_normalize_species_space_and_case():
    assert hg.normalize_species('Oryza Sativa') == 'oryza_sativa'


def test_normalize_species_already_underscored():
    assert hg.normalize_species('oryza_sativa') == 'oryza_sativa'


# ── parse_gget_ref：链接/版本元数据提取 ───────────────────────
def test_parse_gget_ref_extracts_dna_and_gtf():
    meta = hg.parse_gget_ref(CANNED_REF, 'oryza_sativa')
    assert meta['ensembl_release'] == 60
    assert meta['dna_url'].endswith('.fa.gz')
    assert meta['gtf_url'].endswith('.gtf.gz')
    assert meta['bytes'] == '400.2 MB'


def test_parse_gget_ref_missing_species():
    with pytest.raises(RuntimeError):
        hg.parse_gget_ref(CANNED_REF, 'dryosaurus')


def test_parse_gget_ref_missing_dna_link():
    canned = {'x_y': {'genome_dna': {'ftp': ''}}}
    with pytest.raises(RuntimeError):
        hg.parse_gget_ref(canned, 'x_y')


def test_gget_ref_lookup_without_gget_installed(monkeypatch):
    """gget 未安装时给带安装提示的 RuntimeError，而不是裸 ImportError。"""
    import builtins
    real_import = builtins.__import__

    def fake_import(name, *a, **kw):
        if name == 'gget':
            raise ImportError('No module named gget')
        return real_import(name, *a, **kw)

    monkeypatch.setattr(builtins, '__import__', fake_import)
    with pytest.raises(RuntimeError, match='pip install'):
        hg.gget_ref_lookup('Oryza sativa')


# ── _http：域名白名单（SSRF 防护）────────────────────────────
def test_http_whitelist_allows_ensembl():
    s = hg._http(hg.ENSEMBL_HOSTS)
    assert s._check('http://ftp.ensembl.org/pub/release-60/') .startswith('http')


def test_http_whitelist_rejects_foreign():
    s = hg._http(hg.ENSEMBL_HOSTS)
    with pytest.raises(ValueError):
        s._check('http://evil.example.com/x.fa.gz')


def test_http_default_whitelist_is_ncbi():
    s = hg._http()
    assert s._check('https://ftp.ncbi.nlm.nih.gov/genomes/x.fna.gz')
    with pytest.raises(ValueError):
        s._check('http://ftp.ensembl.org/pub/')


# ── _emitters：log 回调适配（Web 任务面板链路）────────────────
def test_emitters_routes_to_log_callback():
    msgs = []
    E = hg._emitters(msgs.append)
    E.ok('下载完成')
    E.warn('注意')
    E.err('失败')
    E.header('标题')
    assert msgs[0] == '✓ 下载完成'
    assert msgs[1] == '⚠ 注意'
    assert msgs[2] == '✗ 失败'
    assert any('标题' in m for m in msgs)


def test_emitters_none_returns_ui():
    assert hg._emitters(None) is hg.UI


# ── finalize_host_dir：genome.fa 锚点落盘 ────────────────────
def _tiny_fasta(path, n=2):
    with open(path, 'wt', encoding='utf-8') as f:
        for i in range(n):
            f.write(f'>chr{i +1} test seq\nACGT{"ACGT" * 3}\n')
    return path


def test_finalize_host_dir_copies_genome_fa(tmp_path):
    merged = tmp_path / 'all.genome.uniq.fasta'
    _tiny_fasta(str(merged))
    msgs = []
    genome_fa = hg.finalize_host_dir(str(tmp_path), str(merged),
                                     log=msgs.append)
    assert os.path.isfile(genome_fa)
    assert os.path.basename(genome_fa) == 'genome.fa'
    with open(genome_fa, encoding='utf-8') as f:
        assert f.readline().startswith('>chr1')
    assert any('genome.fa' in m for m in msgs)


def test_merge_then_finalize_end_to_end(tmp_path):
    """合并去重 + genome.fa 落锚点的最小闭环（离线）。

    merge 走 safe_open（写模式限定平台内），产物放在平台 run/ 下临时目录，
    结束即清理。"""
    import shutil

    from Virus_Platform_Core.config import PLATFORM_ROOT
    work = os.path.join(PLATFORM_ROOT, 'run', '_test_host_genome_gget')
    os.makedirs(work, exist_ok=True)
    try:
        a = _tiny_fasta(os.path.join(work, 'a.fa'), n=2)
        b = _tiny_fasta(os.path.join(work, 'b.fa'), n=2)  # 同名序列 → 去重
        merged = os.path.join(work, 'all.genome.uniq.fasta')
        n_written, n_dups = hg.merge_and_deduplicate([a, b], merged,
                                                     log=lambda m: None)
        assert n_written == 4 and n_dups == 2
        genome_fa = hg.finalize_host_dir(work, merged, log=lambda m: None)
        assert hg.count_sequences(genome_fa) == 4
    finally:
        shutil.rmtree(work, ignore_errors=True)


# ── gunzip_full：逐层解压（Ensembl 双层 .gz 实测回归）─────────
def _gzip_file(src, dst):
    import shutil

    import gzip
    with open(src, 'rb') as fi, gzip.GzipFile(dst, 'wb', mtime=0) as fo:
        shutil.copyfileobj(fi, fo)
    return dst


def _workdir(name):
    import shutil

    from Virus_Platform_Core.config import PLATFORM_ROOT
    work = os.path.join(PLATFORM_ROOT, 'run', name)
    os.makedirs(work, exist_ok=True)
    return work, lambda: shutil.rmtree(work, ignore_errors=True)


def test_gunzip_full_single_layer(tmp_path):
    work, cleanup = _workdir('_test_gunzip_single')
    try:
        raw = _tiny_fasta(os.path.join(work, 'x.fa'), n=1)
        gz = _gzip_file(raw, raw + '.gz')
        msgs = []
        out = hg.gunzip_full(gz, log=msgs.append)
        assert out == os.path.join(work, 'x.fa')
        assert hg.count_sequences(out) == 1
        assert any('解压' in m for m in msgs)
    finally:
        cleanup()


def test_gunzip_full_peels_double_gzip():
    """Ensembl Genomes release 63 实测：部分文件是双层 gzip，
    只解一层会留下仍带 gzip 魔数的文件、下游静默产出空结果。"""
    work, cleanup = _workdir('_test_gunzip_double')
    try:
        raw = _tiny_fasta(os.path.join(work, 'y.fa'), n=2)
        gz1 = _gzip_file(raw, raw + '.gz')
        gz2 = _gzip_file(gz1, gz1 + '.gz')          # 双层嵌套
        out = hg.gunzip_full(gz2, log=lambda m: None)
        assert not out.endswith('.gz')
        assert hg.count_sequences(out) == 2
        with open(out, 'rb') as f:
            assert f.read(2) != b'\x1f\x8b'         # 已是明文，无 gzip 魔数
    finally:
        cleanup()


def test_gunzip_full_passthrough_plain(tmp_path):
    work, cleanup = _workdir('_test_gunzip_plain')
    try:
        raw = _tiny_fasta(os.path.join(work, 'z.fa'), n=1)
        msgs = []
        out = hg.gunzip_full(raw, log=msgs.append)
        assert out == raw                            # 非压缩文件原样返回
        assert not msgs
    finally:
        cleanup()
