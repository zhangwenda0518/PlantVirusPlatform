# Public-Data · 公共数据检索与下载

> 返回 [Home](Home) ｜ 上一页：[Phylogeography](Phylogeography) ｜ 下一页：[Toolbox](Toolbox)

对应引擎：`public_meta/`（search_engine / info_engine / landscape_plot / host_genome）、`public_data.py`（下载管理）。

形成 **「检索 → 下载 → 分析」闭环**：按物种检索公共测序数据 → 提取统一元数据 → 勾选 Run 一键转入下载 → 自动配对建样品进分析流程。

---

## 公共数据检索（meta_search）

源自 MMPV-RNA public_metadata_pipeline：

- **SRA + GSA 双引擎**：NCBI SRA（esearch/esummary）+ CNCB NGDC GSA 按物种检索，结果合并；
- **详细模式**：解析 Tissue / Location / 发育阶段（可选 `--deepseek-api` 做 AI 元数据清洗，未配 Key 自动纯规则模式）；
- 产物（`run/meta_search/<物种>/`）：
  - `search/SRA_GSA_Merged_Final.csv`（检索表）+ `download_links.csv`（直链）；
  - `info/Global_Unified_Metadata_Core14.csv|tsv`（**14 列核心表**）与 `Full.csv|tsv`（**34 列全维表**）；
  - `plot/`（时间/机构/组织/地理出版级图组）。

## 统一元数据提取（meta-info）

Run 列表（SRR/ERR/DRR/CRR）→ Core14/Full 统一元数据；`--fill-date` 采集日期缺失时用发布日期兜底；local/api 双模式。

## SCI 级元数据可视化（meta-plot）

对检索表或 Core14/Full 出**出版级图组**：采样时间分布、机构排名、组织类型、地理分布等。

## 批量下载（数据下载页 / public_data）

| 能力 | 说明 |
|---|---|
| 输入 | SRR/ERR/DRR（ENA 直下 FASTQ.GZ）、CRR（NGDC/GSA）、任意数据文件 URL |
| 引擎 | **aria2c 多线程**批量下载（断点续传）；内置 HTTP 回退；NGDC 低并发防封控 |
| 完整性 | 字节数 + gzip 完整性校验 |
| GSA .sra | **内置 sracha**（纯 Rust 引擎，比 fasterq-dump 快 5-13 倍，Windows 版已随平台编译）自动转 FASTQ.GZ——无需安装 sra-tools |
| 转入分析 | 一键转入分析流程（自动配对 R1/R2 建样品） |
| 进度 | 独立进度窗口可弹出挂机观察 |

## 宿主基因组下载（host-genome，双源）

- **NCBI 源**（默认）：按物种下载宿主参考基因组 FASTA（datasets CLI 优先、E-utilities 回退；`--include-organelles` 同时下叶绿体/线粒体）；另支持 `--accession GCF_...` 按组装号精确下载（可 `--chromosomes X,Y` 染色体过滤，taxid 从数据报告自动解析、目录自动重命名）；
- **gget/Ensembl 源**（`--source gget`）：走 [gget](https://github.com/pachterlab/gget) → Ensembl / Ensembl Genomes（含 **Ensembl Plants**，覆盖水稻、番茄、小麦等 NCBI 覆盖不全的植物宿主）；`--release` 可钉 Ensembl 版本；需先 `pip install "gget==0.30.5"`（依赖较重，作为可选依赖维护在 requirements-dev.txt）；
- 两源产物同目录约定：`host-db/<taxid>_<物种>/` 下 `all.genome.uniq.fasta`（合并去重）+ **`genome.fa`（平台建库锚点名，建库页/配置按它识别宿主源目录）** + `genome_report.txt`（gget 源另有 `gget_ref.json` 记录 Ensembl 版本与链接，便于溯源）；
- 下载域名校验：NCBI 源仅允许 ncbi.nlm.nih.gov，gget 源仅允许 ftp.ensembl.org / ftp.ensemblgenomes.org；Ensembl 部分发布文件为**双层 gzip**，下载器按魔数逐层解到明文；
- Web 入口：数据资源组「基因组查询下载」独立卡（/refdata，2026-09-30 自建库页迁出）——下载完成显示产物路径与 TaxID，一键跳「数据库构建」页建库；
- `--verify-only` 只检查已有下载。

## NCBI datasets CLI 全模式下载（datasets / datasets-setup）

- 四种数据类型全接：`genome`（taxon / accession+染色体过滤）、`gene`（gene-id / symbol+taxon / 蛋白 accession）、`taxonomy`、`virus`（taxon + **宿主过滤** `--host`，或直接取蛋白 `--protein S`）；
- 工具定位：PATH → 平台内置 `3rd/tools/datasets/datasets.exe`；缺失时 `python main.py datasets-setup` 从 NCBI 官方（域名白名单内）自动下载；
- v14+ 已移除 `--exclude-*` 旗标，包装层自动翻译为 `--include` 组合（旧命令照写即可）；
- **virus 模式与建库上下游直通**：默认按平台数据集布局落位（`databases/virusref_db/<名>/reference.fasta + reference.ref_info.tsv`，列名对齐 `load_ref_info` 同义词表），数据库构建页「原始数据集」卡即时识别——下载 → 一键建鉴定库索引 / kunpeng 分类库，全链路不出平台；Web 入口即数据资源组「基因组查询下载」卡的「⬇ 下载病毒数据集」表单（同组还有 genome/gene/taxonomy 通用下载）。

## 界面与 CLI

- 导航「公共数据检索」「数据下载」两页；检索结果勾选 Run → 转入下载 → 分析。

```bat
python main.py meta-search --species "Lycium chinense"
python main.py meta-search --species "Lycium chinense" --db sra --source TRANSCRIPTOMIC --no-detailed
python main.py meta-info --runs meta_search\Lycium_chinense\search\sra.list --fill-date
python main.py meta-info --runs SRR24100141
python main.py meta-plot --input meta_search\...\SRA_GSA_Merged_Final.csv
python main.py host-genome --species "Lycium barbarum" --include-organelles
python main.py host-genome --species "Oryza sativa" --source gget
python main.py datasets virus --taxon "tomato brown rugose fruit virus" --name tobrfv
```

## 相关

- 下载后进流程 → [Pipeline-Overview](Pipeline-Overview)
- 宿主基因组建库 → [Databases](Databases)
