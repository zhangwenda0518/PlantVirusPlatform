#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
known_virus_suite.py — 已知病毒五段整合模块 主入口
=================================================
鉴定 → 过滤 → 共识 → 绘图 → 变异注释

用法:
    # 一键全跑
    python known_virus_suite.py all --reference ref.fasta --ref-info ref_info.tsv \\
        -1 sample_R1.fq.gz -2 sample_R2.fq.gz --engine salmon --out out_dir

    # 建索引（可复用）
    python known_virus_suite.py index --reference ref.fasta --out out_dir

    # 只跑鉴定
    python known_virus_suite.py identify --reference ref.fasta -1 R1.fq.gz -2 R2.fq.gz \\
        --out out_dir

    # 只跑过滤（读已有鉴定结果）
    python known_virus_suite.py filter --out out_dir --min-cov 10

    # 默认两段（平台「已知病毒识别与定量」卡的默认）：
    # 识别/定量（identify）+ 过滤（filter）；共识/绘图/变异由下游模块各自完成
    python known_virus_suite.py identify+filter --reference ref.fasta \\
        -1 R1.fq.gz -2 R2.fq.gz --out out_dir

    # 只跑共识
    python known_virus_suite.py consensus --reference ref.fasta -1 R1.fq.gz --out out_dir

    # 只跑变异注释（读共识段产出的 BAM + GenBank 缓存，caller 为 bcftools mpileup+call）
    python known_virus_suite.py variant --reference ref.fasta --out out_dir \
        --bam-dir out_dir/bam --variant-gbk-dir out_dir/gbk_files

引擎分工（定量引擎 --engine 二选一，默认 salmon；共识段固定 minibwa）:
    salmon  —— 定量（identify）：quant.sf 出 EM 定量，--writeBam 出映射位点
               （排序 BAM + pandepth 覆盖度/深度，命名 <样本>_pandepth.*）
    minibwa —— 定量（identify）真比对计数（idxstats），同一覆盖度后端与
               命名；2026-09-13 恢复（原管线 traditional 路径）
    minibwa —— 共识（consensus）内部固定使用的真比对引擎，产共识/变异用 BAM
               （自建单参考索引重新比对，与定量引擎互不混用）
"""

import argparse
import os
import sys
import time
from pathlib import Path

from .kv_common import ToolRegistry, fmt_time, load_ref_info, setup_logger
from .kv_engines import make_engine
from .kv_filter import FilterStage
from .kv_identify import IdentifyStage

# 平台「已知病毒识别与定量」卡的默认阶段：只跑 识别/定量（identify）+ 过滤（filter）。
# 共识序列、深度绘图、变异注释不属于「识别与定量」，由下游模块
# （t-consensus 共识序列分析 / t-variant 病毒变异分析 / t-kvchain 一键全流程）各自完成
# —— 这与管道阶段 Virus_Platform_Core/kv_stage.py 的口径一致（那里也只跑这两段）。
DEFAULT_STAGE = 'identify+filter'


# ── 工具目录自动发现（Windows 本机 tools/）────────────────
def _extra_tool_dirs(base=None):
    """外部工具目录的补充候选（base 非空时视为平台根覆盖）。

    `ToolRegistry.probe()` 内部已用 `kv_common.platform_root()` 扫
    `3rd/tools`、`3rd/bin` 及其各工具子目录，本函数只补它不覆盖的位置。

    2026-09-11 修正：原先按 `here/'..'/'tools'` 拼路径，那是 2026-09-10
    目录整理**之前**的顶层布局。整理后 here =
    <平台根>/Virus_Platform_Core/known_virus_suite，`here/'..'/'tools'`
    指向并不存在的 Virus_Platform_Core/tools，整个列表在 Windows 上恒为
    空 —— 属静默降级（probe() 自己找得到，故一直没暴露）。
    """
    if base:
        root = Path(base)
    else:
        try:
            from .kv_common import platform_root
            root = platform_root()
        except Exception:
            root = Path(__file__).resolve().parent.parent.parent
    cands = [
        root / '3rd' / 'tools',                     # 植物病毒分析平台/3rd/tools
        root / '3rd' / 'bin',                       # minibwa 等单文件 exe
        root / '3rd' / 'tools' / 'minimap2',
        root / '3rd' / 'tools' / 'salmon2',
        # mafft-win/ms/bin（MSYS 运行时里的 mafft）：目录已随 2026-10 瘦身移除，
        # 候选保留仅为向后兼容老包/自装 MSYS 的机器 —— 下方 exists() 会静默过滤，
        # 不存在时自动跳过（不影响实际探测结果）。
        root / '3rd' / 'tools' / 'mafft-win' / 'ms' / 'bin',
        Path('/usr/local/bin'), Path('/usr/bin'),
    ]
    return [str(c.resolve()) for c in cands if c.exists()]


def collect_samples(args):
    """
    支持两种样本指定方式：
      1. 单样本: -1 R1 -2 R2 (--sample-name 指定名字)
      2. 批量:   --sample-sheet tsv（列: name, r1, r2）
    """
    if getattr(args, 'sample_sheet', None):
        sheet = Path(args.sample_sheet)
        if not sheet.exists():
            sys.exit(f"样本表不存在: {sheet}")
        samples = []
        with open(sheet, 'r', encoding='utf-8') as f:
            head = f.readline().rstrip('\n').split('\t')
            idx = {h.strip().lower(): i for i, h in enumerate(head)}
            for line in f:
                p = line.rstrip('\n').split('\t')
                if len(p) < 2 or not p[0].strip():
                    continue
                name = p[idx.get('name', 0)].strip()
                r1 = p[idx.get('r1', 1)].strip()
                r2 = p[idx['r2']].strip() if 'r2' in idx and len(p) > idx['r2'] else None
                samples.append({'name': name, 'r1': r1, 'r2': r2 or None})
        return samples

    if not args.r1:
        sys.exit("必须提供 -1/--r1（或 --sample-sheet）")
    name = args.sample_name or Path(args.r1).name.split('.')[0]
    return [{'name': name, 'r1': args.r1, 'r2': args.r2}]


def _index_dir(args, out_dir):
    """索引目录：--index-dir 优先，否则 <out>/index。

    指向已有索引目录可跨 run 复用（build_index 内部已有存在即复用逻辑），
    避免每次新建 run 都重建 8 千条参考的 60 MB 索引。
    """
    d = getattr(args, 'index_dir', None)
    return Path(d) if d else (out_dir / 'index')


def build_parser():
    p = argparse.ArgumentParser(
        description='已知病毒五段整合模块：鉴定 → 过滤 → 共识 → 绘图 → 变异注释',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__)
    sub = p.add_subparsers(dest='stage', required=True)

    def common(sp):
        sp.add_argument('--out', required=True, help='输出目录')
        sp.add_argument('--reference', help='参考序列 FASTA')
        sp.add_argument('--ref-info', dest='ref_info', help='参考注释 TSV（Accession/Taxid/Species/Segment）')
        sp.add_argument('--engine', default='salmon', choices=['salmon', 'minibwa'],
                        help='定量引擎（库与引擎正交：salmon=EM 定量默认；'
                             'minibwa=真比对计数；共识段内部固定 minibwa 不受此参影响）')
        sp.add_argument('--threads', type=int, default=max(1, (os.cpu_count() or 4) // 2),
                        help='样本级并行/批次大小')
        sp.add_argument('--align-threads', dest='align_threads', type=int, default=8,
                        help='引擎内部线程数')
        sp.add_argument('--batch-size', dest='batch_size', type=int, default=1,
                        help='每批处理样本数（断点续跑粒度）')
        sp.add_argument('--min-mapq-identify', dest='min_mapq_identify', type=int, default=10,
                        help='（已弃用/保留兼容）鉴定段比对 MAPQ 下限：pandepth 后端'
                             '照原管线 MAPQ 0 起算、本参数无消费者，仅内部保留传参')
        sp.add_argument('--coverage-tool', dest='coverage_tool', default='pandepth',
                        choices=['pandepth', 'samtools', 'builtin'],
                        help='覆盖度/深度统计后端。默认 pandepth：严格照 '
                             'virome_analysis_pipeline 的方式（samtools view -F 0x04 '
                             '建 BAM + pandepth -a，不传 -q）。缺件时自动降级。'
                             'builtin 已弃用（自研解析器已删除，收到即告警回退 pandepth）')
        sp.add_argument('--index-dir', dest='index_dir', default=None,
                         help='比对索引目录（默认 <out>/index）；指向已有索引可跨运行复用，'
                              '避免每次重建 8 千条参考的索引')
        sp.add_argument('--verbose', action='store_true')

    # identify+filter —— 平台的默认两段（见 DEFAULT_STAGE）。
    for name in ('all', DEFAULT_STAGE, 'index', 'identify', 'filter',
                 'consensus', 'plot', 'variant'):
        help_txt = ('只跑 识别/定量 + 过滤 两段（平台默认；共识/绘图/变异交给下游模块）'
                    if name == DEFAULT_STAGE else f'只跑 {name} 段')
        sp = sub.add_parser(name, help=help_txt)
        common(sp)

    # 输入
    # variant 段用不到 fastq 输入，但平台 _tool_job_kvsuite 不分 stage 统一拼参，
    # 因此这里也注册（接受即丢），避免 argparse 报 unrecognized arguments。
    for name in ('all', DEFAULT_STAGE, 'identify', 'consensus', 'variant'):
        sp = sub.choices[name]
        sp.add_argument('-1', '--r1', help='R1 fastq(.gz)')
        sp.add_argument('-2', '--r2', help='R2 fastq(.gz)')
        sp.add_argument('--sample-name', dest='sample_name', help='样本名（默认取文件名前缀）')
        sp.add_argument('--sample-sheet', dest='sample_sheet', help='批量样本表 TSV(name,r1,r2)')

    # 鉴定阈值
    # 同理，variant 段不用这些阈值（它直接吃 BAM），注册仅为容错平台统一拼参。
    for name in ('all', DEFAULT_STAGE, 'identify', 'filter', 'variant'):
        sp = sub.choices[name]
        g = sp.add_argument_group('过滤/鉴定阈值')
        g.add_argument('--min-cov', dest='min_cov', type=float, default=10.0, help='覆盖率下限 %% (默认 10)')
        g.add_argument('--min-depth', dest='min_depth', type=float, default=0.5, help='平均深度下限 (默认 0.5)')
        g.add_argument('--min-reads', dest='min_reads', type=int, default=10, help='最少 reads (默认 10)')
        g.add_argument('--min-tpm', dest='min_tpm', type=float, default=1.0, help='最低 TPM (默认 1.0)')
        g.add_argument('--min-poisson', dest='min_poisson', type=float, default=0.3, help='泊松比值下限 (默认 0.3)')
        g.add_argument('--ratio', type=float, default=0.3, help='双轨 A 轨泊松比 (默认 0.3)')
        g.add_argument('--min-ani', dest='min_ani', type=float, default=0.0,
                       help='（已弃用/保留兼容）最低 ANI：当前无消费者'
                            '（ANI 分种用 --ani-thresh），收到非 0 值即告警并忽略')
        g.add_argument('--ani-thresh', dest='ani_thresh', type=float, default=95.0, help='ANI 分种阈值 (默认 95)')
        g.add_argument('--genes-cov', dest='genes_cov', help='基因覆盖表 TSV，启用 B 轨')
        g.add_argument('--min-gene-total-cov', dest='min_gene_total_cov', type=float, default=80.0)
        g.add_argument('--min-gene-avr-cov', dest='min_gene_avr_cov', type=float, default=5.0)

    # 共识参数
    for name in ('all', 'consensus'):
        sp = sub.choices[name]
        g = sp.add_argument_group('共识参数（照搬 batch_virus_variants.py 默认值）')
        g.add_argument('--vc-qual', dest='vc_qual', type=int, default=20,
                       help='最低碱基质量 Phred 值 (默认 20，对齐原管线 -q)')
        g.add_argument('--vc-depth', dest='vc_depth', type=int, default=5,
                       help='最低覆盖深度 (默认 5，对齐原管线 -d；传 0 则按 MeanDepth/2 动态算)')
        g.add_argument('--vc-freq', dest='vc_freq', type=float, default=0.5,
                       help='变异最低频率阈值 (默认 0.5，对齐原管线 -f)')
        g.add_argument('--vc-ambig', dest='vc_ambig', default='N',
                       help='低于阈值的碱基用此字符代替 (默认 N，对齐原管线 -a)')
        g.add_argument('--jobs', '-j', type=int, default=4,
                       help='共识段并行任务数 (默认 4，对齐原管线 --jobs)')
        g.add_argument('--resume', action='store_true', default=False,
                       help='断点续传：已有共识产物则跳过 (对齐原管线 --resume)')

    # 绘图参数
    for name in ('all', 'plot'):
        sp = sub.choices[name]
        g = sp.add_argument_group('绘图参数（照搬 batch_plot_virus_depth.py --mode depth）')
        g.add_argument('--window', type=int, default=10,
                       help='深度曲线滑动平均窗口 (默认 10)')
        g.add_argument('--plot-fontsize', dest='plot_fontsize', type=int, default=9,
                       help='统计信息框字号 (默认 9)')
        g.add_argument('-g', '--add-genes', dest='add_genes', action='store_true',
                       default=True,
                       help='叠加 GenBank 基因轨道 (默认开)')
        g.add_argument('--no-genes', dest='add_genes', action='store_false',
                       help='关闭基因轨道（不联网下载注释）')
        g.add_argument('--gbk-dir', dest='gbk_dir', default=None,
                       help='GenBank 缓存目录 (默认 <out>/gbk_files)')
        g.add_argument('--depth-dir', dest='depth_dir', default=None,
                       help='per-position 深度表目录 (默认 <out>/align)')
        g.add_argument('--summary', default=None,
                       help='鉴定/过滤结果表 (默认 <out>/filter/filtered.tsv)')
        g.add_argument('--ncbi-email', dest='ncbi_email', default=None,
                       help='NCBI Entrez 邮箱 (或设环境变量 NCBI_EMAIL)')
        g.add_argument('--ncbi-api-key', dest='ncbi_api_key',
                       default=os.environ.get('NCBI_API_KEY'),
                       help='NCBI API key (或设环境变量 NCBI_API_KEY；'
                            '推荐环境变量，key 不进命令行/日志)')

    # 变异检出与注释参数
    for name in ('all', 'variant'):
        sp = sub.choices[name]
        g = sp.add_argument_group('变异检出与注释（bcftools caller + SnpEff 5.4c）')
        g.add_argument('--bam-dir', dest='bam_dir', default=None,
                       help='BAM 目录 (默认 <out>/bam，共识段的比对产物)')
        g.add_argument('--input-vcf', dest='input_vcf', default=None,
                       help='直接指定 VCF 文件作为变异输入（跳过 bcftools caller）；'
                            '与 --bam-dir 二选一，同时给出时 VCF 优先')
        g.add_argument('--variant-gbk-dir', dest='variant_gbk_dir', default=None,
                       help='GenBank 缓存目录 (默认 <out>/gbk_files；'
                            '与绘图段的 --gbk-dir 独立)')
        g.add_argument('--variant-ncbi-email', dest='variant_ncbi_email', default=None,
                       help='NCBI Entrez 邮箱')
        g.add_argument('--variant-ncbi-api-key', dest='variant_ncbi_api_key',
                       default=os.environ.get('NCBI_API_KEY'),
                       help='NCBI API key (或设环境变量 NCBI_API_KEY)')
        g.add_argument('--variant-qual', dest='variant_qual', type=float,
                       default=3.5,
                       help='QUAL 下限 (默认 3.5，bcftools 实测定稿；'
                            '原管线的 20 是 freebayes 标度，不可移植)')
        g.add_argument('--min-freq', dest='min_freq', type=float,
                       default=0.05, help='最小等位频率 (默认 0.05)')
        # 变异图氨基酸标签（参照 246 plot_snpeff.py 口径）
        g.add_argument('--aa-label-cutoff', dest='aa_label_cutoff', type=float,
                       default=0.50,
                       help='变异图上标注氨基酸变化的 AF 下限 (默认 0.50；'
                            '设为 1.01 可关闭标签，设为 0 则全部标注)')
        g.add_argument('--max-aa-labels', dest='max_aa_labels', type=int,
                       default=40,
                       help='单张图上氨基酸标签数量上限 (默认 40，防拥挤)')
        g.add_argument('--all-variants', dest='nonsyn_only',
                       action='store_false', default=True,
                       help='变异图画全部变异（默认只画非同义）')
        g.add_argument('--bcftools', dest='bcftools', default=None,
                       help='bcftools 可执行文件路径 (默认平台 tools/bcftools)')
        g.add_argument('--samtools', dest='samtools', default=None,
                       help='samtools 可执行文件路径 (默认平台 tools/samtools)')

        # 扩展变异分析（三模块）
        g2 = sp.add_argument_group('扩展变异分析（分子谱 / 变异密度 / 群体遗传）')
        g2.add_argument('--no-variant-evo', dest='run_variant_evo',
                        action='store_false', default=True,
                        help='关闭扩展变异分析（分子谱/密度/群体遗传成图）')
        g2.add_argument('--no-snpgenie', dest='run_snpgenie',
                        action='store_false', default=True,
                        help='关闭 SNPGenie（不做 πN/πS 与 dN/dS）')
        g2.add_argument('--evo-accession', dest='evo_accession', action='append',
                        default=None,
                        help='只对指定参考做扩展变异分析（可重复；默认全部）')

    # 平台兼容层：_tool_job_kvsuite 不区分 stage 统一拼全套参数，
    # 而各 stage 只注册自己用得到的。这里给每个 stage 补齐缺失的开关
    # （已注册的跳过），避免 argparse 报 unrecognized arguments。
    _register_platform_compat(sub)

    return p


def _register_platform_compat(sub):
    """让每个子命令都能接受平台统一拼出的全套参数。

    已存在的参数跳过；补齐的用 argparse.SUPPRESS 作 default，
    不覆盖任何显式语义（取值会被忽略）。
    """
    compat = [
        (('-1', '--r1'), {'dest': 'r1', 'help': argparse.SUPPRESS}),
        (('-2', '--r2'), {'dest': 'r2', 'help': argparse.SUPPRESS}),
        (('--sample-name',), {'dest': 'sample_name', 'help': argparse.SUPPRESS}),
        (('--sample-sheet',), {'dest': 'sample_sheet', 'help': argparse.SUPPRESS}),
        (('--min-cov',), {'dest': 'min_cov', 'type': float, 'help': argparse.SUPPRESS}),
        (('--min-depth',), {'dest': 'min_depth', 'type': float, 'help': argparse.SUPPRESS}),
        (('--min-reads',), {'dest': 'min_reads', 'type': int, 'help': argparse.SUPPRESS}),
        (('--min-tpm',), {'dest': 'min_tpm', 'type': float, 'help': argparse.SUPPRESS}),
        (('--min-poisson',), {'dest': 'min_poisson', 'type': float, 'help': argparse.SUPPRESS}),
        (('--ratio',), {'type': float, 'help': argparse.SUPPRESS}),
        (('--min-ani',), {'dest': 'min_ani', 'type': float, 'help': argparse.SUPPRESS}),
        (('--ani-thresh',), {'dest': 'ani_thresh', 'type': float, 'help': argparse.SUPPRESS}),
        (('--genes-cov',), {'dest': 'genes_cov', 'help': argparse.SUPPRESS}),
        (('--min-gene-total-cov',), {'dest': 'min_gene_total_cov', 'type': float,
                                     'help': argparse.SUPPRESS}),
        (('--min-gene-avr-cov',), {'dest': 'min_gene_avr_cov', 'type': float,
                                   'help': argparse.SUPPRESS}),
        (('--vc-qual',), {'dest': 'vc_qual', 'type': int, 'help': argparse.SUPPRESS}),
        (('--vc-depth',), {'dest': 'vc_depth', 'type': int, 'help': argparse.SUPPRESS}),
        (('--vc-freq',), {'dest': 'vc_freq', 'type': float, 'help': argparse.SUPPRESS}),
        (('--vc-ambig',), {'dest': 'vc_ambig', 'help': argparse.SUPPRESS}),
        (('--jobs', '-j'), {'dest': 'jobs', 'type': int, 'help': argparse.SUPPRESS}),
        (('--resume',), {'dest': 'resume', 'action': 'store_true', 'help': argparse.SUPPRESS}),
        (('--window',), {'type': int, 'help': argparse.SUPPRESS}),
        (('--plot-fontsize',), {'dest': 'plot_fontsize', 'type': int,
                                'help': argparse.SUPPRESS}),
        (('-g', '--add-genes'), {'dest': 'add_genes', 'action': 'store_true',
                                 'default': argparse.SUPPRESS, 'help': argparse.SUPPRESS}),
        (('--no-genes',), {'dest': 'add_genes', 'action': 'store_false',
                           'default': argparse.SUPPRESS, 'help': argparse.SUPPRESS}),
        (('--gbk-dir',), {'dest': 'gbk_dir', 'help': argparse.SUPPRESS}),
        (('--depth-dir',), {'dest': 'depth_dir', 'help': argparse.SUPPRESS}),
        (('--summary',), {'help': argparse.SUPPRESS}),
        (('--ncbi-email',), {'dest': 'ncbi_email', 'help': argparse.SUPPRESS}),
        (('--ncbi-api-key',), {'dest': 'ncbi_api_key', 'help': argparse.SUPPRESS}),
        (('--bam-dir',), {'dest': 'bam_dir', 'help': argparse.SUPPRESS}),
        (('--input-vcf',), {'dest': 'input_vcf', 'help': argparse.SUPPRESS}),
        (('--variant-gbk-dir',), {'dest': 'variant_gbk_dir', 'help': argparse.SUPPRESS}),
        (('--variant-ncbi-email',), {'dest': 'variant_ncbi_email', 'help': argparse.SUPPRESS}),
        (('--variant-ncbi-api-key',), {'dest': 'variant_ncbi_api_key',
                                       'help': argparse.SUPPRESS}),
        (('--variant-qual',), {'dest': 'variant_qual', 'type': float,
                               'help': argparse.SUPPRESS}),
        (('--min-freq',), {'dest': 'min_freq', 'type': float, 'help': argparse.SUPPRESS}),
        (('--aa-label-cutoff',), {'dest': 'aa_label_cutoff', 'type': float,
                                  'help': argparse.SUPPRESS}),
        (('--max-aa-labels',), {'dest': 'max_aa_labels', 'type': int,
                                'help': argparse.SUPPRESS}),
        (('--all-variants',), {'dest': 'nonsyn_only', 'action': 'store_false',
                               'default': argparse.SUPPRESS, 'help': argparse.SUPPRESS}),
        (('--bcftools',), {'dest': 'bcftools', 'help': argparse.SUPPRESS}),
        (('--samtools',), {'dest': 'samtools', 'help': argparse.SUPPRESS}),
        (('--no-variant-evo',), {'dest': 'run_variant_evo', 'action': 'store_false',
                                 'default': argparse.SUPPRESS, 'help': argparse.SUPPRESS}),
        (('--no-snpgenie',), {'dest': 'run_snpgenie', 'action': 'store_false',
                              'default': argparse.SUPPRESS, 'help': argparse.SUPPRESS}),
        (('--evo-accession',), {'dest': 'evo_accession', 'action': 'append',
                                'help': argparse.SUPPRESS}),
    ]
    for sp in sub.choices.values():
        existing = set()
        for a in sp._actions:  # noqa: SLF001 - argparse 无公开 API 枚举已注册选项
            existing.update(a.option_strings)
        for opt, kw in compat:
            if any(o in existing for o in opt):
                continue
            sp.add_argument(*opt, **kw)


def main():
    parser = build_parser()
    args = parser.parse_args()
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    logger, _ = setup_logger(out_dir, verbose=args.verbose, tag='known_virus_suite')

    # 弃用参数告警（audit 第 16 项）：参数保留兼容（web 层/旧脚本可能还在传），
    # 但收到即显式告警，不再让它们静默失效。
    if getattr(args, 'coverage_tool', '') == 'builtin':
        logger.warning('--coverage-tool builtin 已弃用（自研 SAM 解析器已删除），'
                       '本次回退 pandepth（缺件时自动降级 samtools）')
    if float(getattr(args, 'min_ani', 0.0) or 0.0) > 0:
        logger.warning('--min-ani 已弃用（无消费者；ANI 分种请用 --ani-thresh）：'
                       '传入值 %s 被忽略', args.min_ani)
    if int(getattr(args, 'min_mapq_identify', 10)) != 10:
        logger.warning('--min-mapq-identify 已弃用（pandepth 后端不消费，'
                       '照原管线 MAPQ 0 起算）：传入值 %s 被忽略',
                       args.min_mapq_identify)

    print("=" * 64)
    print(" 已知病毒五段整合模块  known_virus_suite")
    print(f" 阶段: {args.stage}   引擎: {args.engine}   输出: {out_dir}")
    print("=" * 64)

    # 工具探测
    tools = ToolRegistry(logger)
    tools.probe(extra_dirs=_extra_tool_dirs())
    print(tools.report())
    print()

    t_all = time.time()

    # ── index ──
    if args.stage == 'index':
        tools.require(args.engine)
        eng = make_engine(args.engine, tools, args.align_threads, logger)
        p = eng.build_index(args.reference, _index_dir(args, out_dir),
                            args.align_threads)
        logger.info(f"索引完成 -> {p}")
        print(f"\n索引路径: {p}")
        return

    # ── identify / all / identify+filter ──
    confirmed = None
    ref_info = {}
    if args.stage in ('all', DEFAULT_STAGE, 'identify'):
        tools.require(args.engine)
        if args.engine == 'minibwa':
            tools.require('samtools')   # 真比对定量产 BAM 依赖 samtools
        if not args.reference:
            sys.exit("鉴定段需要 --reference")
        samples = collect_samples(args)
        logger.info(f"样本数: {len(samples)}")

        eng = make_engine(args.engine, tools, args.align_threads, logger,
                          min_mapq=getattr(args, 'min_mapq_identify', 10),
                          coverage_tool=getattr(args, 'coverage_tool', 'pandepth'))
        stage = IdentifyStage(args, eng, tools, logger, out_dir)
        ref_info = stage.ref_info
        df, confirmed, novel = stage.run(samples)

    # ── filter / all / identify+filter ──
    if args.stage in ('all', DEFAULT_STAGE, 'filter'):
        if args.reference and not ref_info:
            ref_info = load_ref_info(args.ref_info, logger)
        fstage = FilterStage(args, logger, out_dir, ref_info=ref_info)
        if args.stage == 'filter':
            kept, dropped = fstage.run()
        else:
            # 走到这里必然跑过识别段（stage 同为 'all'/DEFAULT_STAGE）。
            # 此前退化成分支用空 DataFrame 兜底——真被触发时会写出一份
            # "全部剔除"的 filter 结果冒充过滤结论；宁可直接失败。
            if confirmed is None:
                raise RuntimeError(
                    'filter 段缺少识别段结果（confirmed 为空）：'
                    '请先跑 identify，或使用 --stage filter 从磁盘结果表过滤')
            kept, dropped = fstage.apply(confirmed)
        confirmed = kept

    # ── consensus / all ──
    if args.stage in ('all', 'consensus'):
        from .kv_consensus import ConsensusStage
        if not args.reference:
            sys.exit("共识段需要 --reference")
        # 共识段固定 minibwa 真比对（定量引擎 salmon 无真实 CIGAR/NM，
        # 撑不起共识位点）；samtools 负责管道与索引，viral_consensus 出共识。
        tools.require('minibwa', 'samtools', 'viral_consensus')
        samples = collect_samples(args)
        cstage = ConsensusStage(args, tools, logger, out_dir)
        cstage.run(samples, confirmed)

    # ── plot / all ──
    if args.stage in ('all', 'plot'):
        from .kv_plot import PlotStage
        pstage = PlotStage(args, tools, logger, out_dir, cfg={
            'window': getattr(args, 'window', 10),
            'fontsize': getattr(args, 'plot_fontsize', 9),
            'add_genes': getattr(args, 'add_genes', True),
            'gbk_dir': getattr(args, 'gbk_dir', None),
            'ncbi_email': getattr(args, 'ncbi_email', None),
            'ncbi_api_key': getattr(args, 'ncbi_api_key', None),
        })
        pstage.run()

    # ── variant / all ──
    if args.stage in ('all', 'variant'):
        from .kv_variant import VariantStage
        if not args.reference:
            sys.exit("变异段需要 --reference")
        vstage = VariantStage(args, tools, logger, out_dir, cfg={
            'bam_dir': getattr(args, 'bam_dir', None),
            'input_vcf': getattr(args, 'input_vcf', None),
            'gbk_dir': getattr(args, 'variant_gbk_dir', None) or getattr(args, 'gbk_dir', None),
            'ncbi_email': getattr(args, 'variant_ncbi_email', None) or getattr(args, 'ncbi_email', None),
            'ncbi_api_key': getattr(args, 'variant_ncbi_api_key', None) or getattr(args, 'ncbi_api_key', None),
            'qual': getattr(args, 'variant_qual', 3.5),
            'min_minor_freq': getattr(args, 'min_freq', 0.05),
            'threads': getattr(args, 'align_threads', 4),
            'bcftools': getattr(args, 'bcftools', None) or tools.paths.get('bcftools'),
            'samtools': getattr(args, 'samtools', None) or tools.paths.get('samtools'),
            'env': getattr(tools, 'env', None),
        })
        vstage.run()
        ann_stat = getattr(vstage, '_annotations', None) or {}
        if ann_stat:
            logger.info(
                "注释产物: virus-annotations/ %d 个 .gb, %d 个 .gtf"
                % (ann_stat.get('gb', 0), ann_stat.get('gtf', 0)))
        ref_fastas = getattr(vstage, '_ref_fastas', None) or {}
        if ref_fastas:
            logger.info("单参考 FASTA: virus-fasta/ %d 条（SNPGenie 输入）"
                        % len(ref_fastas))

        # 下游两段与本段共用同一个 GenBank 源：显式传下去，避免
        # --variant-gbk-dir 只对第一段生效（否则基因轨道/SNPGenie 静默丢失）。
        _gbk_src = (getattr(args, 'variant_gbk_dir', None)
                    or getattr(args, 'gbk_dir', None))

        # 变异可视化（紧跟变异段，产物落到 <out>/variant_plots/）
        try:
            from .kv_variant_plot import run_variant_plots
            vpres = run_variant_plots(
                out_dir, logger,
                nonsyn_only=getattr(args, 'nonsyn_only', True),
                label_af_cutoff=getattr(args, 'aa_label_cutoff', 0.50),
                max_aa_labels=getattr(args, 'max_aa_labels', 40),
                gbk_dir=_gbk_src)
            logger.info(f"变异绘图：{vpres.get('n_accessions', 0)} 条参考 → "
                        f"{vpres.get('n_plots', 0)} 张图")
        except Exception as e:
            logger.warning(f"变异绘图失败（不影响变异段产物）: {e}")

        # 扩展变异分析：分子突变谱 / 变异密度 / 群体遗传学（SNPGenie）
        if getattr(args, 'run_variant_evo', True):
            try:
                from .kv_variant_evo import run_variant_evo
                eres = run_variant_evo(
                    out_dir, logger,
                    run_snpgenie_stage=getattr(args, 'run_snpgenie', True),
                    accession_filter=getattr(args, 'evo_accession', None),
                    label_af_cutoff=getattr(args, 'aa_label_cutoff', 0.50),
                    max_aa_labels=getattr(args, 'max_aa_labels', 40),
                    gbk_dir=_gbk_src)
                logger.info(
                    f"扩展变异分析：{eres.get('n_accessions', 0)} 条参考 → "
                    f"{eres.get('n_plots', 0)} 张图，"
                    f"{len(eres.get('snpgenie', []))} 个 SNPGenie 结果")
                _sg_skip = eres.get('snpgenie_skipped') or []
                if _sg_skip:
                    logger.warning(
                        "SNPGenie 未产出 %d 条（是工具未运行，非无多态）："
                        % len(_sg_skip))
                    for _s in _sg_skip:
                        logger.warning("  %s: %s"
                                       % (_s.get('accession'), _s.get('reason')))
            except Exception as e:
                logger.warning(f"扩展变异分析失败（不影响变异段产物）: {e}")

    logger.info(f"全流程耗时 {fmt_time(time.time()-t_all)}")
    print("\n" + "=" * 64)
    print(f" 完成。产物目录: {out_dir.resolve()}")
    print("=" * 64)


if __name__ == '__main__':
    main()
