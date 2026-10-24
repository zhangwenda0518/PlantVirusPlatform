# -*- coding: utf-8 -*-
"""为 examples/results/variant/kvsuite 合成演示变异并跑真注释链。

背景：examples 里 variant 运行的 VCF 是 0 条数据线的空壳（demo 无真 BAM），
变异棒棒糖图（kv_variant_plot.run_variant_plots）因此没有
annotated/*.ann.tsv 输入、无法复跑与图形体检。

本脚本补齐整条演示数据：
1. 从 kvsuite 示例的 virus-fasta 拼合三参考 FASTA；
2. 以固定种子合成演示 VCF（SNV 为主 + 少量 indel，位置取自真实参考序列）；
   VCF 为**多样本**（S1/S2/S3 基因型列，固定携带先验），跨样本 SNP 分型图
   （snipit，kv_variant_plot.run_snipit_plots）据此演示「一行一样本」视图；
3. 走平台真链路 VariantStage(input_vcf) → snpEff 建库/注释
   → annotated/{gid}.ann.tsv + virus-annotations/ + virus-fasta/，
   随后 run_variant_plots 出全套图件（含 snipit 分型图）。

产出是"合成变异 × 真实注释/绘图代码"，演示与体检两用。重跑幂等（覆盖）。
"""
from __future__ import annotations

import argparse
import logging
import random
import sys
from pathlib import Path

_ROOT = Path(__file__).resolve().parents[1]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from Virus_Platform_Core.known_virus_suite.kv_variant import VariantStage  # noqa: E402

KVSUITE = _ROOT / 'examples' / 'results' / 'kvsuite' / 'kvsuite'
VARIANT_RUN = _ROOT / 'examples' / 'results' / 'variant' / 'kvsuite'
SEED = 42
SNV_PER_ACC = 40        # 每参考的 SNV 数（演示量级）
INDEL_PER_ACC = 3       # 每参考的短 indel 数
DEMO_SAMPLES = ('S1', 'S2', 'S3')   # 演示样本名（与 BAM 命名口径一致）
CARRY_P = {'S1': 0.60, 'S2': 0.50, 'S3': 0.40}   # 各样本携带变异的先验


def collect_reference() -> tuple[dict[str, str], Path]:
    """拼合三参考 FASTA → ({seq_id: seq}, 合并文件路径)。"""
    refs: dict[str, str] = {}
    for sub in sorted(KVSUITE.glob('virus-fasta/ref_*')):
        fa = next(sub.glob('*.fasta'), None)
        if fa is None:
            continue
        seq_id, chunks = None, []
        for line in fa.read_text(encoding='utf-8').splitlines():
            if line.startswith('>'):
                seq_id = line[1:].split()[0]
            elif seq_id:
                chunks.append(line.strip())
        if seq_id and chunks:
            refs[seq_id] = ''.join(chunks)
    if not refs:
        raise SystemExit(f'未在 {KVSUITE}/virus-fasta 找到参考 FASTA')
    combined = VARIANT_RUN / 'reference.fasta'
    combined.parent.mkdir(parents=True, exist_ok=True)
    with combined.open('w', encoding='utf-8', newline='\n') as f:
        for sid, seq in refs.items():
            f.write(f'>{sid} demo combined reference\n')
            for i in range(0, len(seq), 60):
                f.write(seq[i:i + 60] + '\n')
    return refs, combined


def synth_vcf(refs: dict[str, str], out_vcf: Path) -> int:
    """合成演示 VCF：多样本（S1/S2/S3 基因型列），SNV 为主 + 少量 1-3 bp
    indel，固定种子可复现。每行至少一个样本携带（全 0/0 的行分型图没得画）。"""
    rng = random.Random(SEED)
    alt_of = {'A': 'CGT', 'C': 'AGT', 'G': 'ACT', 'T': 'ACG'}
    lines = ['##fileformat=VCFv4.2', '##source=demo-synthetic']
    for sid, seq in refs.items():
        lines.append(f'##contig=<ID={sid},length={len(seq)}>')
    lines.append('##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">')
    lines.append('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    cols = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO',
            'FORMAT', *DEMO_SAMPLES]
    lines.append('\t'.join(cols))
    n_total = 0
    for sid, seq in refs.items():
        positions = sorted(rng.sample(range(1, len(seq) - 4), SNV_PER_ACC
                                      + INDEL_PER_ACC))
        snv_pos = set(positions[:SNV_PER_ACC])
        for p in positions:
            base = seq[p - 1].upper()
            if base not in alt_of:
                continue
            if p in snv_pos:
                ref, alt = base, rng.choice(alt_of[base])
            elif rng.random() < 0.5:                    # 缺失 1-2 bp
                ref = base + ''.join(rng.choice('ACGT')
                                     for _ in range(rng.randint(1, 2)))
                alt = base
            else:                                       # 插入 1-2 bp
                ref = base
                alt = base + ''.join(rng.choice('ACGT')
                                     for _ in range(rng.randint(1, 2)))
            # 基因型：按先验抽携带样本，强制至少一个携带者
            carriers = [s for s in DEMO_SAMPLES if rng.random() < CARRY_P[s]] \
                or [rng.choice(DEMO_SAMPLES)]
            gts = '\t'.join('1/1' if s in carriers else '0/0'
                            for s in DEMO_SAMPLES)
            dp = rng.randint(20, 800)
            lines.append('\t'.join((sid, str(p), '.', ref, alt,
                                    f'{rng.uniform(80, 900):.1f}', 'PASS',
                                    f'DP={dp}', 'GT', gts)))
            n_total += 1
    out_vcf.write_text('\n'.join(lines) + '\n', encoding='utf-8')
    return n_total


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--out', default=str(VARIANT_RUN),
                    help='variant 运行目录（缺省 examples demo run）')
    a = ap.parse_args()
    out_dir = Path(a.out)

    logging.basicConfig(level=logging.INFO,
                        format='%(levelname)s %(message)s')
    log = logging.getLogger('make_variant_demo')

    refs, combined = collect_reference()
    log.info('参考拼合：%d 条（%s）', len(refs),
             '、'.join(f'{s}({len(q)})' for s, q in refs.items()))

    vcf_dir = out_dir / 'vcf'
    vcf_dir.mkdir(parents=True, exist_ok=True)
    demo_vcf = vcf_dir / 'demo.variants.vcf'
    n = synth_vcf(refs, demo_vcf)
    log.info('合成演示变异 %d 条 -> %s', n, demo_vcf)

    class _Args:
        reference = str(combined)
        ncbi_email = None
        ncbi_api_key = None

    stage = VariantStage(_Args(), tools={}, logger=log, out_dir=out_dir,
                         cfg={'input_vcf': str(demo_vcf),
                              'gbk_dir': str(KVSUITE / 'gbk_files'),
                              'threads': 2})
    results = stage.run()
    ok = [r for r in results if r.get('n_coding_ann')]
    log.info('注释完成：%d / %d 条参考产出编码区注释', len(ok), len(results))
    for r in results:
        log.info('  %s: 变异 %s，编码注释 %s', r.get('genome'),
                 r.get('n_variants'), r.get('n_coding_ann'))

    # 出全套图件（与 suite 的变异可视化同一路径，含 snipit 跨样本分型图）
    from Virus_Platform_Core.known_virus_suite.kv_variant_plot import (
        run_variant_plots,
    )
    vpres = run_variant_plots(out_dir, log,
                              gbk_dir=str(KVSUITE / 'gbk_files'))
    sn = vpres.get('snipit') or {}
    log.info('变异绘图：%d 条参考 → %d 张图；跨样本分型图（snipit）'
             '%d 条参考 / %d 个产物',
             vpres.get('n_accessions', 0), vpres.get('n_plots', 0),
             sn.get('n_accessions', 0), len(sn.get('plots', [])))
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
