# Pipeline-Overview · 分析流程总览

> 返回 [Home](Home) ｜ 上一页：[Architecture](Architecture) ｜ 下一页：[Preprocess](Preprocess)

流程编排器 `Virus_Platform_Core/pipeline.py`：按阶段依赖链调度、断点续跑、加权总进度、预计剩余时间与逐阶段资源预估。GUI 与 CLI 共用。

---

## 14 步级联（按功能模块分组显示）

| 组 | 步骤 | 说明 | 缺依赖时 |
|---|---|---|---|
| 🧹 测序数据预处理 | ⓪ fastp 质控 | 自动去接头+低质量裁剪（双端/单端） | 灰显可跳过 |
| | ⓪b FASTQ→FASTA | seqkit fq2fa → conv_R1/R2.fa.gz，①分类直接复用提速 | 自动绕过 |
| | ① 宿主去除 | kunpeng 宿主库分类 → 剔除宿主 read 对 | 需先建宿主库 |
| 🦠 病毒鉴定 | ②b 识别与定量 | kvsuite：识别 + 二次过滤定量 + 比对深度 + 病毒 reads 提取 | — |
| 🧬 病毒组装 | ③ 组装·分类·病毒contigs | SPAdes（rnaviral 默认/metaviral/meta/rna/isolate）+ kunpeng+BLAST 分类 | — |
| | ③b 候选序列验证 | 宿主筛选 → 长度分流 → blastx/BFVD + 类病毒 blastn | — |
| | ③c 共识序列与变异 | reads 回贴参考 → 逐参考位置计数 → 共识 | — |
| 🧲 宿主预测 | ④ ICTV 宿主判定 | kunpeng 分类 → ICTV 宿主概率级联 + NCBI 元数据交叉 | — |
| 🔬 下游分析 | ⑥ ORF 预测 | pyrodigal（meta）+ pyrodigal-rv（RNA 病毒）+ orfipy（显式指定） | — |
| | ⑥b ORF 功能注释 | DIAMOND/MMseqs2/blastp × RefSeq 蛋白 + Pfam HMM + BFVD 病毒域 | — |
| | ⑦ 进化树与 SDT | MAFFT+trimAl + FastTree/RAxML-NG/NJ；SDT identity 矩阵 | — |
| | ⑧ 引物设计 | primer3 保守区/全长 + 热力学评分 + 特异性检查（可选） | — |
| | ⑨ 基因组图 | gbdraw 首选 / dna_features_viewer 回退 | 两引擎全缺才灰显 |
| 📄 报告 | ⑩ 可视化报告 | plotly/pycirclize/matplotlib 全图汇总 report.html | — |

- **每步输出自动作为下一步输入**；上游完成后下游自动解锁，任务结束页面自动刷新。
- **组装输入可选**：默认「病毒 reads」（②b 比对上的病毒 reads，宿主/杂菌污染最少）；病毒 reads 不足 500 对自动回退去宿主 reads/原始 reads 并在日志说明。
- **组装自动降级**：metaviral 低深度拼不出 contigs 自动改 rna 模式重试；单端数据 metaviral/meta 自动转 rna。
- **操作**：`▶ 运行此步`（单跑）/ `⏩ 依次运行到此` / `↻ 重跑此步`（配合强制重跑忽略断点）/ 顶部 `▶ 依次运行剩余步骤`。
- **卡片内直接暴露参数**：组装模式/输入/内存/最小 contig、最小 ORF、参考数/建树工具、引物模式/特异性检查、出图上限/自备 fasta+gff/gb；右上全局参数（线程数/分类置信度/chunk 目录/强制重跑）。设置页「分析默认参数」为初值，逐步可覆盖。

## 断点与重跑

- 每阶段完成落 `.done` 标记；`--force` 或卡片勾选「强制重跑」忽略断点。
- SPAdes 已有组装结果时重跑直接复用，不重复计算。

## 进度模型：全局进度 · 预计剩余时间 · 资源预估

- 任务卡显示**全局进度**——按阶段耗时**加权**（非步数均分）+ 已运行时长 + **预计剩余时间**。
- 剩余时间来自每阶段耗时**自学习模型**（`run/logs/stage_perf.json`，随使用越来越准；首次用缺省粗估）。
- 无原生进度的引擎也有近似上报：kunpeng 分类按 chunk 中间盘增长、SPAdes 按 spades.log 的 k-mer 阶段、过滤/提取/转换按记录数。
- 每阶段开始/结束输出 `📊 资源预估`（线程/内存/磁盘）与实际耗时；分类阶段结束额外输出 `📊 磁盘核对`（预估 vs 实测 chunk 占用）。

## 并发闸门

| 类别 | 同时上限 | 涵盖 |
|---|---|---|
| 重任务 | 2（可配 1-8） | kunpeng 分类 / SPAdes / DIAMOND 注释 / 建库 / 建树 / SDT |
| 轻任务 | 4（可配 1-16） | 下载 / 转换 / 绘图 / 检索 |

超出自动排队，任务卡显示「排队中 · 前面还有 N 个任务」。在 `platform.json` 的 `defaults` 设 `max_heavy_tasks` / `max_light_tasks` 调整。批处理队列（分析流程页「批量导入」TSV）顺序执行，同时只跑一个样品。

## crabz 加速（可选）

平台根目录有 `crabz.exe`（已内置）时，所有 .gz 读写自动走 crabz 多线程管道（kept/viral reads、子样本、fq2fa 产物等），实测写出较 Python gzip 级别 9 **快约 10 倍**（195MB 模拟 FASTQ：21.5s → 2.1s），压缩率相当；无 crabz 自动回退，无功能差异。

## 相关页面

- 各步骤细节 → [Preprocess](Preprocess) / [Kvsuite](Kvsuite) / [Assembly](Assembly) / [Host-Prediction](Host-Prediction) / [ORF](ORF) / [Phylo](Phylo) / [Primer](Primer) / [Genome-Plots](Genome-Plots) / [Report](Report)
- 运维（tool_runs 归档、自检）→ [Settings-Ops](Settings-Ops)
