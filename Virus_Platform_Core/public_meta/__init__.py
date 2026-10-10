# -*- coding: utf-8 -*-
"""Virus_Platform_Core.public_meta — 公共数据检索与元数据模块。

源自 MMPV-RNA public_metadata_pipeline 的功能子集：

  - search_engine.py   SRA + GSA 双引擎物种检索（合并去重 → 检索结果 CSV）
  - info_engine.py     Run 元数据深度提取 → Core14 / Full 统一元数据表
  - prescreen_engine.py 病毒样本预筛（NCBI STAT / CNCB GSA **预计算**分类，
                       不下载测序数据即可判断哪些 run 含病毒）
  - virome.py          病毒谱：把预筛产物从「样本视角」翻成「病毒视角」
                       （检出清单 / 分类归属 / 占比 / 跨样本发生率）。
                       纯读产物、不联网；Web 侧直接 import，不是子进程引擎。
  - landscape_plot.py  检索/元数据 SCI 级可视化
  - host_genome.py     宿主参考基因组下载（NCBI datasets/E-utilities + gget·Ensembl 双源）
  - datasets_cli.py    NCBI datasets CLI 全模式封装（genome/gene/taxonomy/virus+宿主过滤）

下载环节复用平台既有「数据下载」模块（ENA/NGDC 直下 + aria2c + fasterq-dump），
不引入原管线的 gsa_sra.down.py / sra2fastx.py。建宿主库请用平台
build-host-db（kunpeng），不引入原管线的 kraken2/bowtie2 建库脚本。

CLI：python main.py meta-… / host-genome；Web：导航「公共数据检索」页。
病毒预筛走 Web 导航「病毒样本预筛」页（引擎 CLI 见 prescreen_engine.py --help）；
病毒谱是同页结果区的「病毒谱」标签（排查用 CLI 见 virome.py --help）。
"""
