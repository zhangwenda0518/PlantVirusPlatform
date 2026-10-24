# -*- coding: utf-8 -*-
"""跨样本 SNP 分型图（snipit 封装）单测：VCF 基因型解析 / 伪比对构建 /
编排入口的容错降级 / 端到端出图。

snipit 未安装时端到端用例自动跳过（kv_variant_plot 对缺包场景静默降级，
这里只保证「装了就必须产出」）；伪比对构建本身零 snipit 依赖，始终验证。
"""
import importlib.util
import logging
from pathlib import Path

import pytest

from Virus_Platform_Core.known_virus_suite import kv_variant_plot as kvp

REF = 'ACGT' * 15          # 60 bp 合成参考
SAMPLES = ['S1', 'S2', 'S3']


def _base(p: int) -> str:
    return REF[p - 1]


def _write_multisample_vcf(path: Path) -> None:
    """3 样本合成 VCF：SNV / 共享 SNV / 缺失 / 插入 / 无基因型各一。"""
    lines = [
        '##fileformat=VCFv4.2',
        f'##contig=<ID=REF1,length={len(REF)}>',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'
        + '\t'.join(SAMPLES),
        # S1 SNV
        f'REF1\t11\t.\t{_base(11)}\tT\t100\tPASS\t.\tGT\t1/1\t0/0\t0/0',
        # S1/S2 共享 SNV（S2 纯合）
        f'REF1\t21\t.\t{_base(21)}\tG\t100\tPASS\t.\tGT\t0/1\t1/1\t0/0',
        # S1 缺失 TT（REF 跨度 3bp → ALT 1bp）
        f'REF1\t31\t.\t{_base(31)}{_base(32)}{_base(33)}\t{_base(31)}'
        f'\t100\tPASS\t.\tGT\t1/1\t0/0\t0/0',
        # 插入行：伪比对表达不了，必须整行跳过
        f'REF1\t41\t.\t{_base(41)}\t{_base(41)}A\t100\tPASS\t.\tGT'
        f'\t1/1\t0/0\t0/0',
        # S3 SNV；S1 无基因型 → 记参考碱基
        f'REF1\t51\t.\t{_base(51)}\tC\t100\tPASS\t.\tGT\t./.\t0/0\t1/1',
    ]
    path.write_text('\n'.join(lines) + '\n', encoding='utf-8')


# ── _read_vcf_genotypes ─────────────────────────────────────────
def test_read_vcf_genotypes(tmp_path):
    vcf = tmp_path / 't.vcf'
    _write_multisample_vcf(vcf)
    samples, rows = kvp._read_vcf_genotypes(vcf)
    assert samples == SAMPLES
    assert len(rows) == 5
    assert rows[0]['GTS'] == ['1/1', '0/0', '0/0']
    assert rows[3]['ALTS'] == [_base(41) + 'A']      # 插入行原样保留


def test_read_vcf_genotypes_sites_only(tmp_path):
    vcf = tmp_path / 't.vcf'
    vcf.write_text(
        '##fileformat=VCFv4.2\n'
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
        f'REF1\t11\t.\t{_base(11)}\tT\t100\tPASS\t.\n', encoding='utf-8')
    samples, rows = kvp._read_vcf_genotypes(vcf)
    assert samples == [] and len(rows) == 1


# ── build_snipit_alignment ──────────────────────────────────────
def test_build_snipit_alignment(tmp_path):
    vcf = tmp_path / 'REF1.vcf'
    _write_multisample_vcf(vcf)
    aln = tmp_path / 'aln.fasta'
    st = kvp.build_snipit_alignment(vcf, 'REF1', REF, aln, chrom='REF1')

    assert st['n_samples'] == 3
    assert st['n_sites'] == 4                        # 11/21/31/51；41 被跳过
    assert st['n_rows_skipped_insertion'] == 1

    seqs = {}
    name = None
    for line in aln.read_text(encoding='utf-8').splitlines():
        if line.startswith('>'):
            name = line[1:]
            seqs[name] = ''
        else:
            seqs[name] += line.strip()
    assert len(seqs) == 4 and set(seqs) == {'REF1', *SAMPLES}
    for s in seqs.values():
        assert len(s) == len(REF)                    # 定长伪比对

    exp_s1 = list(REF)
    exp_s1[10] = 'T'                                 # 11 G>T
    exp_s1[20] = 'G'                                 # 21 共享 SNV
    exp_s1[30:33] = [_base(31), '-', '-']            # 31 缺失 TT → 补 gap
    assert seqs['S1'] == ''.join(exp_s1)
    exp_s2 = list(REF)
    exp_s2[20] = 'G'
    assert seqs['S2'] == ''.join(exp_s2)
    exp_s3 = list(REF)
    exp_s3[50] = 'C'                                 # 51 S3 SNV
    assert seqs['S3'] == ''.join(exp_s3)
    assert seqs['REF1'] == REF                       # 参考行原样


# ── run_snipit_plots：容错降级 ───────────────────────────────────
def test_run_snipit_plots_without_snipit(tmp_path, monkeypatch):
    """snipit 未安装 → available=False 静默返回，不抛异常。"""
    monkeypatch.setattr(kvp, '_find_snipit_runner', lambda: None)
    res = kvp.run_snipit_plots(tmp_path)
    assert res['available'] is False and res['n_accessions'] == 0


def test_run_snipit_plots_sites_only_skipped(tmp_path):
    """sites-only VCF → 按 accession 记 skipped，不产图。"""
    (tmp_path / 'vcf').mkdir()
    (tmp_path / 'vcf' / 'REF1.vcf').write_text(
        '##fileformat=VCFv4.2\n'
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
        f'REF1\t11\t.\t{_base(11)}\tT\t100\tPASS\t.\n', encoding='utf-8')
    res = kvp.run_snipit_plots(tmp_path)
    assert res['n_accessions'] == 0
    assert res['skipped'][0]['reason'] == 'VCF 无样本列（sites-only）'


def test_run_snipit_plots_skips_caller_intermediates(tmp_path):
    """all.variants.* / demo.variants.* 中间产物不进 skipped 噪音。"""
    (tmp_path / 'vcf').mkdir()
    body = ('##fileformat=VCFv4.2\n'
            '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n'
            f'REF1\t11\t.\t{_base(11)}\tT\t100\tPASS\t.\tGT\t1/1\t0/0\n')
    for name in ('all.variants.vcf', 'demo.variants.vcf'):
        (tmp_path / 'vcf' / name).write_text(body, encoding='utf-8')
    res = kvp.run_snipit_plots(tmp_path)
    assert res['n_accessions'] == 0 and res['skipped'] == []


# ── 端到端（需 snipit 已安装）────────────────────────────────────
needs_snipit = pytest.mark.skipif(
    importlib.util.find_spec('snipit') is None,
    reason='snipit 未安装')


@needs_snipit
def test_run_snipit_plots_end_to_end(tmp_path):
    (tmp_path / 'vcf').mkdir()
    _write_multisample_vcf(tmp_path / 'vcf' / 'REF1.vcf')
    vdir = tmp_path / 'virus-fasta' / 'ref_REF1'
    vdir.mkdir(parents=True)
    (vdir / 'ref_REF1.ref.fasta').write_text(
        '>REF1\n' + '\n'.join(REF[i:i + 60] for i in range(0, len(REF), 60))
        + '\n', encoding='utf-8')

    res = kvp.run_snipit_plots(tmp_path, logging.getLogger('test'))
    assert res['available'] is True and res['n_accessions'] == 1
    assert res['accessions']['REF1']['n_samples'] == 3
    assert res['accessions']['REF1']['n_snps'] == 4

    plot_dir = tmp_path / 'variant_plots'
    for fn in ('REF1_snipit.png', 'REF1_snipit.pdf', 'REF1_snipit.snps.csv'):
        assert (plot_dir / fn).is_file(), f'缺产物 {fn}'
    # snps.csv 首列是记录名（暂存目录改名归位后仍是本 accession 的样本）
    import csv
    with open(plot_dir / 'REF1_snipit.snps.csv', encoding='utf-8') as fh:
        rows = list(csv.DictReader(fh))
    assert {r['record'] for r in rows} <= set(SAMPLES)
