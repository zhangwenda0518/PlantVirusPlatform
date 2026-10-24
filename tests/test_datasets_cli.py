# -*- coding: utf-8 -*-
"""datasets_cli 纯函数单测：NCBI datasets CLI 包装层。

只测不联网的部分：命令构造（含 v14+ 旧语法翻译）、数据报告 → 引用注释表
转换、taxid 提取、FASTA 计数。真实下载链路由 CLI（main.py datasets …）
与 /api/datasets_virus 端点巡检覆盖。
"""
import os
import json

import pytest

from Virus_Platform_Core.public_meta import datasets_cli as dc


# ── 命令构造：镜像官方 datasets download 语法 ─────────────────
def test_build_genome_taxon():
    assert dc.build_genome_args(taxon='bos taurus') == [
        'download', 'genome', 'taxon', 'bos taurus',
        '--include', 'genome,gff3,rna,protein']


def test_build_genome_accession_chromosomes():
    args = dc.build_genome_args(accession='GCF_000001405.40',
                                chromosomes='X,Y')
    assert args == ['download', 'genome', 'accession', 'GCF_000001405.40',
                    '--chromosomes', 'X,Y',
                    '--include', 'genome,gff3,rna,protein']


def test_build_genome_exclude_translated_to_include():
    """v14+ 移除了 --exclude-gff3 等旗标：旧语法须翻译成 include 组合。"""
    args = dc.build_genome_args(accession='GCF_000001405.40',
                                exclude=('gff3', 'rna'))
    assert '--include' in args
    assert args[args.index('--include') + 1] == 'genome,protein'
    assert not any(a.startswith('--exclude') for a in args)


def test_build_genome_include_exclude_combined():
    args = dc.build_genome_args(taxon='x', include='genome,gff3',
                                exclude=('gff3',))
    assert args[-1] == 'genome'


def test_build_genome_requires_target():
    with pytest.raises(ValueError):
        dc.build_genome_args()


def test_build_gene_three_ways():
    assert dc.build_gene_args(gene_id=672) == [
        'download', 'gene', 'gene-id', '672']
    assert dc.build_gene_args(symbol='brca1', taxon='mus musculus') == [
        'download', 'gene', 'symbol', 'brca1', '--taxon', 'mus musculus']
    assert dc.build_gene_args(accession='NP_000483.3') == [
        'download', 'gene', 'accession', 'NP_000483.3']
    with pytest.raises(ValueError):
        dc.build_gene_args()
    with pytest.raises(ValueError):
        dc.build_gene_args(symbol='brca1')          # symbol 缺 taxon


def test_build_taxonomy_multi():
    assert dc.build_taxonomy_args('human,sars-cov-2') == [
        'download', 'taxonomy', 'taxon', 'human,sars-cov-2']
    with pytest.raises(ValueError):
        dc.build_taxonomy_args('')


def test_build_virus_genome_and_protein():
    assert dc.build_virus_args(taxon='sars-cov-2', host='dog') == [
        'download', 'virus', 'genome', 'taxon', 'sars-cov-2',
        '--host', 'dog', '--include', 'genome']
    assert dc.build_virus_args(protein='S', host='dog') == [
        'download', 'virus', 'protein', 'S', '--host', 'dog']
    with pytest.raises(ValueError):
        dc.build_virus_args()


def test_build_virus_accession_single():
    """查询表挑单条：accession 精确下载（taxon/protein 均缺省时生效）。"""
    assert dc.build_virus_args(accession='NC_045512.2') == [
        'download', 'virus', 'genome', 'accession', 'NC_045512.2',
        '--include', 'genome']
    assert dc.build_virus_args(accession='NC_045512.2', host='dog') == [
        'download', 'virus', 'genome', 'accession', 'NC_045512.2',
        '--host', 'dog', '--include', 'genome']
    # taxon 优先级高于 accession（旧行为不变）
    assert dc.build_virus_args(taxon='sars-cov-2', accession='NC_045512.2')[
        3] == 'taxon'


# ── summary 查询（先查后下）────────────────────────────────
def test_build_summary_args_all_kinds():
    assert dc.build_summary_args('genome', taxon='Oryza sativa',
                                 limit=200) == [
        'summary', 'genome', 'taxon', 'Oryza sativa', '--limit', '200']
    # limit 未给时不带 --limit（query_records 才会按 kind 默认值 +1 探测截断）
    assert dc.build_summary_args('genome', accession='GCF_000001405.40') == [
        'summary', 'genome', 'accession', 'GCF_000001405.40']
    assert dc.build_summary_args('virus', taxon='sars-cov-2', host='dog',
                                 limit=5) == [
        'summary', 'virus', 'genome', 'taxon', 'sars-cov-2',
        '--host', 'dog', '--limit', '5']
    assert dc.build_summary_args('gene', gene_id=672) == [
        'summary', 'gene', 'gene-id', '672']
    assert dc.build_summary_args('gene', symbol='brca1', taxon='mus') == [
        'summary', 'gene', 'symbol', 'brca1', '--taxon', 'mus']
    assert dc.build_summary_args('gene', protein_acc='NP_000483.3') == [
        'summary', 'gene', 'accession', 'NP_000483.3']
    # taxonomy 记录数固定，不带 --limit
    assert dc.build_summary_args('taxonomy', taxons='human,x') == [
        'summary', 'taxonomy', 'taxon', 'human,x']


def test_build_summary_args_validation():
    with pytest.raises(ValueError):
        dc.build_summary_args('genome')
    with pytest.raises(ValueError):
        dc.build_summary_args('virus')
    with pytest.raises(ValueError):
        dc.build_summary_args('gene', symbol='brca1')   # 缺 taxon
    with pytest.raises(ValueError):
        dc.build_summary_args('taxonomy')
    with pytest.raises(ValueError):
        dc.build_summary_args('bogus', taxon='x')


def test_dig_dual_style_paths():
    """summary 输出 snake_case、download 报告 camelCase，双风格都要接。"""
    assert dc._dig({'a': {'b_c': 'v'}}, 'a.b_c', 'a.bC') == 'v'
    assert dc._dig({'a': {'bC': 'v'}}, 'a.b_c', 'a.bC') == 'v'
    assert dc._dig({'a': {'b_c': ''}}, 'a.b_c') == ''
    assert dc._dig({}, 'x.y') == ''


def test_norm_summary_row_genome_snake():
    row = dc._norm_summary_row('genome', {
        'accession': 'GCA_1.1',
        'assembly_info': {'assembly_name': 'ASM1', 'assembly_level':
                          'Scaffold', 'release_date': '2020-01-02T00:00:00',
                          'submitter': 'NCBI'},
        'assembly_stats': {'total_sequence_length': '307224746'},
        'organism': {'organism_name': 'Oryza sativa', 'tax_id': 4530},
        'source_database': 'SOURCE_DATABASE_GENBANK'})
    assert row['accession'] == 'GCA_1.1'
    assert row['name'] == 'ASM1'
    assert row['organism'] == 'Oryza sativa'
    assert row['taxid'] == '4530'
    assert row['length'] == '307224746'
    assert row['source'] == 'GENBANK'                # 前缀被剥掉
    assert row['released'] == '2020-01-02'           # 截到日期


def test_norm_summary_row_virus_camel_fallback():
    """download 风格的 camelCase 报告同样能归一化（容错双风格）。"""
    row = dc._norm_summary_row('virus', {
        'accession': 'KT383474.1', 'length': 6393,
        'completeness': 'COMPLETE', 'sourceDatabase': 'GenBank',
        'virus': {'organismName': 'Tomato brown rugose fruit virus',
                  'taxId': 1761477},
        'host': {'organismName': 'Solanum lycopersicum', 'taxId': 4081},
        'isolate': {'name': 'Tom1-Jo', 'collectionDate': '2015-04-17'},
        'location': {'geographicLocation': 'Jordan'}})
    assert row['virus'] == 'Tomato brown rugose fruit virus'
    assert row['host'] == 'Solanum lycopersicum'
    assert row['geo'] == 'Jordan'
    assert row['collected'] == '2015-04-17'
    assert row['completeness'] == 'COMPLETE'


def test_norm_summary_row_gene_and_taxonomy():
    grow = dc._norm_summary_row('gene', {'gene': {
        'gene_id': '672', 'symbol': 'BRCA1', 'taxname': 'Homo sapiens',
        'description': 'BRCA1 DNA repair associated',
        'chromosomes': ['17'], 'type': 'protein-coding',
        'protein_count': 1, 'transcript_count': 5}})
    assert grow['accession'] == '672'
    assert grow['name'] == 'BRCA1'
    assert grow['organism'] == 'Homo sapiens'
    assert grow['chroms'] == '17'
    trow = dc._norm_summary_row('taxonomy', {'taxonomy': {
        'tax_id': 2697049,
        'current_scientific_name': {'name':
            'Severe acute respiratory syndrome coronavirus 2'},
        'group_name': 'viruses',
        'classification': {'family': {'id': 11118, 'name':
                                      'Coronaviridae'}}}})
    assert trow['accession'] == '2697049'
    assert trow['name'] == 'Severe acute respiratory syndrome coronavirus 2'
    assert trow['lineage'] == 'family: Coronaviridae'


# ── 数据报告 → 引用注释表（列名对齐 kv_common.load_ref_info）──
def _report_line(**over):
    rec = {
        'accession': 'KT383474.1',
        'completeness': 'COMPLETE',
        'length': 6393,
        'sourceDatabase': 'GenBank',
        'virus': {'organismName': 'Tomato brown rugose fruit virus',
                  'taxId': 1761477},
        'host': {'organismName': 'Solanum lycopersicum', 'taxId': 4081},
        'isolate': {'name': 'Tom1-Jo', 'collectionDate': '2015-04-17'},
        'location': {'geographicLocation': 'Jordan'},
    }
    rec.update(over)
    return rec


def test_report_to_ref_info_columns(tmp_path, monkeypatch):
    """safe_open 限定平台内写 → 产物放平台 run/ 临时目录。"""
    import shutil

    from Virus_Platform_Core.config import PLATFORM_ROOT
    work = os.path.join(PLATFORM_ROOT, 'run', '_test_ds_refinfo')
    os.makedirs(work, exist_ok=True)
    try:
        out = os.path.join(work, 'reference.ref_info.tsv')
        n = dc.report_to_ref_info([_report_line(), {'accession': ''}],
                                  out, log=lambda m: None)
        assert n == 1                            # 空 accession 行被剔除
        header = open(out, encoding='utf-8').readline().rstrip('\n')
        cols = header.split('\t')
        assert 'Accession' in cols and 'Taxid' in cols and 'Host' in cols
        row = open(out, encoding='utf-8').readlines()[1].rstrip('\n').split('\t')
        d = dict(zip(cols, row))
        assert d['Accession'] == 'KT383474.1'
        assert d['Taxid'] == '1761477'
        assert d['Species'] == 'Tomato brown rugose fruit virus'
        assert d['Host_Taxid'] == '4081'
        assert d['Geo_Location'] == 'Jordan'
    finally:
        shutil.rmtree(work, ignore_errors=True)


def test_report_taxid_paths():
    assert dc.report_taxid([{'organism': {'taxId': 123}}]) == 123
    assert dc.report_taxid([{'virus': {'taxId': 456}}]) == 456
    assert dc.report_taxid([{'accession': 'x'}]) is None


def test_count_fasta_seqs(tmp_path):
    p = tmp_path / 'a.fna'
    p.write_text('>a\nACGT\n>b\nACGT\n', encoding='utf-8')
    assert dc.count_fasta_seqs([str(p)]) == 2
    assert dc.count_fasta_seqs([]) == 0
