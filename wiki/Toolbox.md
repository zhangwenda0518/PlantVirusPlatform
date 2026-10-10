# Toolbox · 专项工具箱

> 返回 [Home](Home) ｜ 上一页：[Public-Data](Public-Data) ｜ 下一页：[NCBI-Submission](NCBI-Submission)

导航「专项分析」页：**不分样品流程的独立分析**，两区组织——🧬 单步工具 与 🧬 参考数据与比较分析。产物统一落 `run/tool_runs/`（规范化双层结构，见 [Settings-Ops](Settings-Ops)）。

---

## 单步工具（把主管道拆开单跑）

| 工具 | 说明 |
|---|---|
| 质控 | fastp 独立运行 |
| 病毒鉴定提取 | kvsuite 识别 + 病毒 reads 提取 |
| 组装 | SPAdes 单跑 |
| contig 分类与深度分析 | kunpeng 分类 / contig 深度注释（Taxonomy 8 级谱系统列）/ 离线检索（BLASTN·DIAMOND BLASTX，`local_search.py`）/ contig 功能画像（含 BFVD 病毒域） |

## 参考数据与比较分析

### NCBI 参考序列下载（ncbi_download）

- Entrez 检索式 → 会话式批量下载（esearch usehistory=y → WebEnv/QueryKey 翻页，借鉴 PhyloSuite）；
- accession 去重、断点续传、含宿主/分离物元数据 TSV；
- 集合名即「⑦ 进化分析 → NCBI 参考集合」，可追加进每组比对（[Phylo](Phylo)）。

```bat
python main.py ncbi-dl "Tobamovirus[ORGN] AND complete genome[TITL]" -n TobamoRefs [--db nucleotide|protein] [--max 100]
python main.py ncbi-list
```

### GenBank 参考集合（gb_collection）

与 ncbi_download 的关键差异：**保留完整 GenBank 记录与 CDS 注释**（供建树/共线性/特征表）：

- 四种来源：**ICTV 界→门→纲→目→科→属→种递进级联下拉**（选得越深集合越聚焦；计数=accession 条数；本地序列优先、缺的 NCBI 补齐，预览离线可用）/ accession 列表 / Entrez 检索式 / 本机 .gb 导入（多记录自动拆分）；
- 每属/每种上限防大属淹没整科集合；
- 一键 **🧬 提取 CDS/PEP**（`cds_export.py`）：genome / CDS / PEP 分类分目录 + 按基因拆分，直接送序列比对或基因建树；
- 集合巡检：记录数/CDS 数/警告。

```bat
python main.py gb-dl -n PotyvirusSet --term "Potyvirus[ORGN] AND complete genome[TITL]" --max 50
python main.py gb-import -n MySet --files a.gb,b.gb
python main.py gb-list
python main.py gb-check -n PotyvirusSet
```

### MSA 查看（msa_view）

SNP-only 变异热图（详见 [Phylo](Phylo)）；可直接输入任意 FASTA 比对。

## dsRNA 全链路设计（dsrna_pipeline / dsrna_scoring / dsrna_offtarget）

RNAi 防控用的 dsRNA 臂设计，四段链路一次跑完：

1. **选窗**（dsRNAmax 思路，maximin 多候选）：扫描病毒基因组选出最优 dsRNA 臂窗口；
2. **效价评分**（`dsrna_scoring.py`）：dsRIP 论文公式的平台统一实现（Cedden et al., BMC Biology 2025, 23:114 忠实移植）——臂/siRNA 效价与安全评分，合并原管线两份复制实现为唯一口径；
3. **宿主脱靶扫描**（`dsrna_offtarget.py`）：k-mer 排序 uint64 数组 + 二分，纯 numpy 实现的大规模扫描（对宿主转录组/基因组排查 siRNA 脱靶）；
4. **引物**：Primer3 设计扩增臂的引物。

> 脱靶面板（`databases/dsrna/panel/`，物种 RefSeq mRNA，GB 级）**不随软件/数据库包分发**（2026-10 瘦身裁定）：不放置面板时自动跳过脱靶扫描与安全评分（效价照常），自备面板放入该目录后即自动参与扫描。

## miRNA 靶标预测（mirna_target/）

5 引擎并行预测病毒序列 ↔ 植物 miRNA 相互作用：

| 引擎 | 说明 |
|---|---|
| psRNA Target | 植物专用 miRNA 靶标预测 |
| psRobot | 精细打分 |
| RNA22 | 模式驱动 |
| RNAhybrid | 最小自由能 |
| Tapir | 打分 + 位点 |

输出统一靶标表（位点、打分、抑制证据），供病毒-宿主互作分析。

## 其他工具

- **SnpEff**（`3rd/tools/jre-snpeff/`，jlink 精简 Java 21 运行时 47MB）：变异注释（③c 变异位点 → 基因/蛋白影响）；
- **环境自检 / 内存自检**：`selfcheck.py` / `mem_check.py`（见 [Settings-Ops](Settings-Ops)）。

## 相关

- 比较基因组四步（参考获取 → 比对 → 建树 → SDT）→ [Phylo](Phylo)
- 产物目录管理 → [Settings-Ops](Settings-Ops)
