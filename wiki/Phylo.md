# Phylo · ⑦ 进化树与 SDT

> 返回 [Home](Home) ｜ 上一页：[ORF](ORF) ｜ 下一页：[Primer](Primer)

对应引擎：`phylo.py`（⑦ 编排）、`sdt_exact.py`（SDT 精确引擎）、`msa_view.py`（MSA 查看）。

---

## 建树流程

1. **分组**：按病毒物种分组（BLAST top hit），每组 = 病毒 contigs + top-N 近缘参考（`--top-n-refs`，默认 10）；
2. **参考来源**：
   - 默认比对 hits（`--tree-sampling blast`）；
   - `macro` = 同科建树（目标属 + 同科各属 3 条背景）；`genus` = 属级树；`lineage` = 种级树（层级抽样策略）——后三者用 `databases/tree_db/plant_tree.db` 植物病毒参考库（48 科 / 452 属 / 6,167 条）；
   - **NCBI 参考集合**追加：`ncbi-dl` 下载的集合（`--ncbi-refs 集合名`）并入每组比对（见 [Toolbox](Toolbox)）；
   - ICTV 口径选参：`ictv-refs --genus Tobamovirus [--download]`（本地优先、缺的 NCBI 按需下载，见 [Databases](Databases)）；
3. **比对**：MAFFT 全长比对 → **trimAl（automated1）清剪**，过度修剪自动回退原比对；Gblocks 可用；
4. **建树**：

| 工具 | 口径 | 适用 |
|---|---|---|
| **FastTree**（默认） | 类 ML，秒级 | 快速查看、序列少 |
| **RAxML-NG** | GTR+G，ML 搜索 + 100 次 FBP 自举支持值；固定种子可复现、线程自动封顶 8 | 发表级 |
| **NJ** | 距离法 | 探索 |

- 序列 <4 条时 RAxML-NG 自动降级 NJ（RAxML-NG 自举要求 ≥4 条）；比对无信息（所有列全 N / 全 gap）导致 RAxML-NG 底层断言崩溃时降级 FastTree 兜底保树——两类降级原因都写入任务日志；最优模型与对数似然写入 05_phylo summary.json。
- 实测口径（209 条 × 966 nt）：RAxML-NG 含自举 730.6s；IQ-TREE 3（GTR+G + 1000 UFBoot）1049s——RAxML 端到端快约 1.4 倍。

## SDT 同一性分析

- 平台用 **SDT 精确引擎**（`sdt_exact.py`，SDTv1.3 的纯 Python 复刻）：MAFFT 比对 + 成对 gap 删除口径重算全长 identity 矩阵（与 SDT 算法一致）；
- 产物：`sdt_matrix.csv` + 热图；同时生成 `sdt_input.fas`，可在结果中心点「启动 SDT」后拖入 SDT v1.3 GUI 交互查看。

## MSA 查看（SNP-only 变异热图）

结果中心内置（`msa_view.py`）：选样品与 05_phylo 分组，只显示比对中的**变异位点**——ACGT 彩色字符热图（A 绿/C 蓝/G 橙/T 红）、共识行、每列变异度柱、悬浮位点计数、大比对自动分页；优先展示建树所用清剪后比对（aln.trim）。参考 PhyloSuite MSA Viewer 的 SNP-only 设计，平台自有前端实现。

## 产物（05_phylo/<组>/）

```
aln.fasta / aln.trim.fasta    比对与清剪后比对
tree.nwk (+tree.png)          Newick 与树图
sdt_matrix.csv                identity 矩阵
sdt_input.fas                 可拖入 SDT GUI
summary.json                  建树参数/模型/耗时
```

## 界面与 CLI

- 分析管道 ⑦ 卡片：参考数/建树工具/抽样策略 + `📊 查看`；结果中心有**交互树查看器**（Newick 渲染）与 MSA 查看。

```bat
python main.py analyze --r1 ... --sample S1 --stages phylo ^
       --tree-tool raxml-ng --tree-sampling lineage --top-n-refs 10 ^
       --ncbi-refs TobamoRefs
python main.py ictv-refs --genus Tobamovirus --download    :: ICTV 选参
python main.py ncbi-dl "Tobamovirus[ORGN] AND complete genome[TITL]" -n TobamoRefs
```

## 相关

- 更高层级的进化动力学（分子钟/地理）→ [Phylodynamics](Phylodynamics) / [Phylogeography](Phylogeography) / [BEAST](BEAST)
- 参考库细节 → [Databases](Databases)
