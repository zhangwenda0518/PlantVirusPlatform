# Assembly · ③ 组装与 contig 分类 · ③b 候选验证 · ③c 共识与变异

> 返回 [Home](Home) ｜ 上一页：[Kvsuite](Kvsuite) ｜ 下一页：[Host-Prediction](Host-Prediction)

对应引擎：`assembly.py`（③）、`verify.py`（③b）、`consensus.py`（③c）、`contig_annot.py`（深度注释）、`local_search.py`（离线比对）。

---

## ③ 组装与分类

**组装（SPAdes）**

- 模式：`rnaviral`（默认）/ `metaviral` / `meta` / `rna` / `isolate`（卡片下拉 + `--assembly-mode`）；
- 组装输入三选一（卡片下拉）：**病毒 reads**（默认，②b 比对上的病毒 reads——宿主/杂菌污染最少、最聚焦；不足 500 对自动回退去宿主 reads/原始 reads 并写日志）/ 去宿主 reads / 原始 reads；
- 自动降级：metaviral 低深度拼不出 contigs → 自动改 rna 模式重试；单端数据 metaviral/meta → 自动转 rna（宏基因组模式不用单端 reads）；
- 内存上限 `--memory`（默认 64GB）；已有组装结果重跑直接复用；
- **中文路径兼容**：SPAdes 经 `%TEMP%\vp_spades`（纯 ASCII）中转，完成后拷回。

**contig 分类（kunpeng + BLAST）**

- kunpeng 对 contigs 分类（kraken2 同构 C 行输出，LCA 判到上级节点时经 nodes/names.dmp 谱系还原完整分类）；
- 对病毒 contig 跑一次本地 BLASTN（vs `databases/virusref_db/` 参考核酸库）填最近参考（accession/identity/覆盖度/物种/科）。

产物（`03_assembly/`）：`spades/`、`contigs.filtered.fasta`（`--min-contig-len` 默认 500）、`viral_contigs.fasta`（病毒子集，④⑦⑧⑨/LOGAN 共用）、`virus_contigs.tsv`（11 列 = kunpeng 5 列 + blast 6 列）、`contig_blast.tsv`。

> **契约**：`virus_contigs.tsv` 是 ④ 宿主预测的输入——blast 六列供 kunpeng 未判到 taxid 的 contig 走「最近参考 → 物种/科 → 宿主概率」回退；手工表缺列时 ④ 自动补齐。

## ③b 候选序列验证

对齐 MMPV-RNA 02b_Filter / 09b_Analysis_Verify（`verify.py`）：

- 宿主污染筛选 → 长度分流 → **blastx**（蛋白级）+ **BFVD 病毒域**（mmseqs2，与 NCBI 提交链同库；原 CDD B 路因 ~5GB 索引反复 OOM 于 2026-09-29 换为 BFVD） + **类病毒 blastn**（环状 RNA 病原专项）；
- 产物（`03b_verify/`）：验证通过的候选集与证据明细。

## ③c 共识序列与变异

- reads 回贴参考 → 逐参考位置计数（ACGT/覆盖度）→ 生成**共识序列**与变异位点表；
- 产物（`03c_consensus/`）：共识 FASTA、变异表；供进化分析（⑦）与提交准备使用。

## contig 深度注释（工具箱配套）

`contig_annot.py`：NCBI Taxonomy 8 级谱系统列（realm→kingdom→…→species）+ 覆盖度/深度汇总，工具箱「contig 分类与深度分析」四件套之一；本地版由 `local_search.py` 提供同一套命中结构（BLASTN / DIAMOND BLASTX，前端无需区分数据来源）。

## CLI

```bat
python main.py analyze --r1 ... --sample S1 --stages assembly ^
       --assembly-mode rnaviral --memory 64 --min-contig-len 500
```

## 相关

- 病毒 contigs 的去向 → [Host-Prediction](Host-Prediction)（④）、[ORF](ORF)（⑥/⑥b）、[Phylo](Phylo)（⑦）、[Primer](Primer)（⑧）、[Genome-Plots](Genome-Plots)（⑨）、[Logan](Logan)（溯源）
