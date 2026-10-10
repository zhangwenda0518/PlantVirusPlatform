# 平台病毒检出性能评测（Bari 小RNA + VIROMOCK 挑战集，全面版）

生成日期：2026-09-16（数据快照）；**2026-10-03 按 `results/results.json` /
`metrics_summary.csv` / `perf.csv` 逐项核对修订**（修订点见文末"修订记录"）。
共 **46 个实验**（`run_bench.py` 现行配置集），由 `run_bench.py` 统一执行：
psutil 全程监控进程树峰值内存、minibwa 单实验 20 分钟限时（驱动侧 monitor 超时
强杀并记 `timed_out`——**归档记录中无一条带该标记**，见 §4 关于 D5/D6 的说明）、
每个实验把关键产物（确诊表/定量表/日志）保留到 `results/artifacts/<exp>/` 后
清理其余中间文件。原始数据 `results/`，图 `figs/`（PNG+SVG，中文，共 10 张）。
工具/环境版本见 `ENVIRONMENT.md`；性能表字段口径见 `results/PERF_NOTES.md`。

> **记录口径**：`results/results.json` 现存 **50 条记录 = 46 个配置实验
> + 2 个早期对照实验**（ctrl_s9_k31 / ctrl_s9_bwa，已不在脚本配置里）
> **+ 2 条断点续跑残留的重复旧 error 记录**（chain_s3 / chain_s9 各 1 条；
> `run_bench.py` 已补 (exp,error) 去重逻辑）。本文数字均以每实验的最新有效
> 记录为准，重复的旧 error 记录不参与统计。

> **参考库版本（重要）**：本文全部基准数字绑定 **2026-09 的 `kv_index` 快照库**
> （物种级，EN 稿记 8,464 条参考；该库未随仓库归档、无 manifest 原件）。
> 2026-09-29 起平台默认库已由 **`virusref_db/plant_ref/`** 取代
> （manifest `built_at=2026-09-19`，12,159 条，含 TSWV L/M/S 注入记录），
> 旧 `kv_index/` 目录已移除。**重跑本文实验需先重建等价库**，或用
> `VP_BENCH_DB` / `--db-dir` 显式指定库目录（各脚本已参数化，缺失时报中文错误）。
> 参见 `ENVIRONMENT.md` §4（库哈希与清单）。

---

## 1. 评测设计

### 1.1 数据集与真值（18 个 VIROMOCK 数据集中的 18 个子集 + Bari 2 样品）

| 数据集 | 类型 | 规模 | 真值 |
|---|---|---|---|
| **Bari Sample_3 / Sample_9** | 小RNA（21–24nt，250K，马铃薯/葡萄） | 各 25 万 reads | Massart et al. 2019 (Phytopathology 109:488) 表 1（2 + 9 项病原，含逐病毒映射 reads 数） |
| **VIROMOCK D11–D18** | 完全人工（150bp PE、纯病毒 reads） | **4.7–5.2 万对/集**（按各集实际 reads 数 ÷2） | 各 Dataset*.md 组成表（accession 级） |
| **VIROMOCK D1–D10** | 半人工/真实（50–301bp，含真实宿主背景） | 6.5 万–4900 万对/集 | README 表 1 + 各 md（背景病毒 + 人工添加病毒，物种级） |

### 1.2 方法矩阵（46 组实验）

| 维度 | 取值 |
|---|---|
| 引擎 | salmon（伪比对 EM 定量）× minibwa（真比对计数）——全部常规 reads 数据集双引擎都测 |
| 小RNA路线 | salmon-k15 适配、srna 组装+BLAST（≥60bp）、全链路（fastp→鉴定+组装）、50K 深度稀释 |
| 对照 | salmon-k31 原始配置（小RNA 上 0 映射）、minibwa 小RNA 拒绝执行 |

检出判定用平台自身口径（确诊表；阈值以平台 argparse 为准，`kv_config.yaml:29-53`
为文档快照）：**A 轨** cov≥10% & mean depth≥0.5 & reads≥10 & TPM≥1 & 泊松比≥0.3；
**B 轨**（RNA 且提供基因覆盖）gene_total_cov≥80 & gene_avr_cov≥5；
**多节段病毒另加"全段检出"规则**（参考库含 Species/Segment 列时自动生效，
任一段未过阈则整组剔除）。
物种级 TP 按真值物种去重，**真值别名表纳入 ICTV 改名**
（CTV=*Closterovirus tristezae*、PVY=*Potyvirus yituberosi*、EMDV=*Alphanucleorhabdovirus melongenae*、
BPEV=*Bell pepper alphaendornavirus*、BYDV=*Luteovirus pavhordei* 等，见 §7.1）。

---

## 2. 检出性能

### 2.1 常规 reads（VIROMOCK）

| 数据集组 | 引擎 | 灵敏度 | 精确率 | 说明 |
|---|---|---|---|---|
| D11–D18 人工集（8 种病毒） | salmon k31 | **8/8 = 100%** | 62% | FP 为近缘物种溢出 |
| D11–D18 人工集 | minibwa | 7/8 = 88% | 58% | 漏检 1 集（D17 *Little cherry virus 1*） |
| D1–D10 半人工/真实集（28 项真值） | salmon | **20/28 = 71.4%** | 32.8% | FP 以同属/同科近缘溢出与疑似参考伪影为主（§7.4） |
| D1–D10 | minibwa | **19/27 = 70.4%** | 35.2% | D5 记录作废（见下），剔除后计；D6 正常完成 1/1 |

双引擎在半人工集灵敏度基本持平——minibwa 按 27 项有效真值计 70.4%，
salmon 按 28 项计 71.4%（**D5 记录作废剔除**；若把 D5 的异常记录也算进去，
minibwa 表面为 20/28 = 71.4%，但该结果依赖已作废的记录，不采用）。
差异点在资源（§4）：耗时中位数约 3.4×（逐集 1.8–7.6×）、峰值内存最大差约 14×
（随输入增大而放大，小集相当）。

**D5/D6 的 minibwa 核查结论（以 `results.json` 为准，2026-10-03 修订）**：
- **D6 正常完成**：274.0 s、31.3M×50bp、mapped 18,038,128 / total 31,327,327，
  物种级检出 1/1——**不存在"超时强杀"**；
- **D5 运行异常提前退出**：0.9 s、峰值 2 MB、**无 `mapped`/`total` 字段**
  （监测失真），该集结果作废、不计入任何指标；
- 归档结果中**没有任何 `timed_out` 标记**（20 分钟限时机制存在但本轮从未触发）。

低灵敏度的归因（均有实据）：**D7 TSWV 库内完全缺失**（双引擎 0%，真实数据库缺口——
TSWV 是最重要的植物病毒之一，建议优先补库）〔**时效标注**：该结论与建议针对
2026-09 评测快照；2026-09-19 起 `plant_ref` 库已注入 TSWV L/M/S
（NC_002052.1/002050.1/002051.1），**TSWV 现已入库**，重测预期可检出〕；
**D3** 是"故意稀释/不完整覆盖"挑战集
（GRSPaV/GLRaV-2/GRVFV 被刻意去除部分 reads， salmon 20%、minibwa 40%）；
**D9 PiVB** 库内无该物种，命中来自疑似误标的 *Pistacia emaravirus* 参考（100% 覆盖，
判 TP 但存疑——数据库准确性问题）。

### 2.2 Bari 小RNA（图 fig_b1 / fig_b2）

| 方法 | 灵敏度 | TP/FP/FN |
|---|---|---|
| salmon-k15（小RNA适配） | **6/11 = 54.5%** | **6/4/5**（合并 S3+S9，精确率 60%） |
| 全链路（fastp→鉴定+组装并集） | S3 2/2 + S9 4/9 | S9 精确率 **28.6%**（TP4/FP10/FN5，组装 contig 归属更模糊） |
| srna组装+BLAST(≥60bp) | 4/11 = 36.4% | 4/5/7 |
| salmon-k31（原始配置） | 0/11 | 小RNA 下完全不可用（见下对照说明） |
| minibwa | 0/11 | 平台入口明确拒绝（rejected=true） |

**对照说明（2026-10-03 补充，数据取自 `results.json` 的 ctrl_* 记录）**：
salmon-k31 对照跑在 **3 参考小库 `idx_k31`**（NC_011620.1 PVX + NC_043447.1 /
NC_043448.1 PVB，`num_refs=3`，2026-09-15 构建），**不是平台全库**——0 映射既是
k=31 > 21–24nt 读长的必然结果，也受小库构成限制；两样品（S3/S9）均测，
`num_mapped=0 / percent_mapped=0.0`。minibwa 对照（S3/S9）为平台入口直接拒绝
（`rejected=true`，message 见 `results.json`）。0/11 = 两样品合计 0/(2+9) 项真值。

**方法变更记录（深度稀释）**：50K 稀释的抽样实现自 2026-10-03 起为
**随机抽样（`random.Random(seed=42).sample`，seed 写入输出元数据）**；
本文 50K 数字基于旧的"顺序截取前 50K 条 reads"实现（无种子、受 fastq 排序影响）。
抽样实现不同、结论不受影响（见 §10）。

对照 Massart 2019：21 实验室 250K 平均 76%（宿主过滤 + 自定阈值 + 专家判读），
本评测按平台默认诊断阈值且未去宿主。S9 漏检的 5 项在 reads 层均有信号，
但**"放宽阈值可到 8/9"的说法已撤回**：该判断所需的阈值扫描产物未随本文档归档，
本次修订无法复核，待专门的阈值扫描实验补齐前不给出数字（阈值与检出率的关系
不属于本文已验证结论）。

### 2.3 深度效应（图 fig_b7）

S3：250K 100% → 50K **100%**（高丰度病原不受深度影响）；
S9：250K 44.4% → 50K **33.3%**（低丰度病原随深度下降）——
与论文"深度每降一个数量级灵敏度大降"的结论一致（46% @50K）。

---

## 3. 定量准确性（图 fig_b4）

预期 vs 观测（平台 EM reads）log-log：**Pearson r = 0.968**（n=17）。
PVB RNA1 8,109 vs 真值 8,333（**97%**）；RNA2 4,197 vs 4,363（**96%**）；
PVX 偏低源于参考选择（泛型 NC_011620 vs 样品自身两分歧分离物），非工具偏差。

**单位口径与系统偏移（2026-10-03 后验核算披露，不重跑实验）**：
- 真值按 **reads** 计（VIROMOCK 少量集为 150bp PE），平台 EM 计的是 **片段
  （fragments）**；PE 数据 **reads = 2×fragments**，故观测值系统性偏小；
- 后验核算（对归档结果复算，log10 空间）：斜率 1.02、**截距 −0.56**、
  观测/真值比值 **中位数 0.43**（范围 0.07–0.97）；审计独立复算为截距 −0.54、
  中位比值 0.40（点位选取微差），结论一致；
- **r = 0.968 只表征相关（方向与强度），不表征准确度**——系统偏移在相关系数里
  不可见；跨引擎/跨单位比较定量数字时必须先统一到同一单位口径。
- 图 fig_b4 已加同一注记（y=x 参考线之外）。

---

## 4. 运行开销（图 fig_b5 / fig_b6）

表内 wall_s / 峰值内存取自 `results/perf.csv`（2026-10-03 对齐 CSV 数字）。

> **字段口径警告**：perf.csv 中 salmon 行的 `total` 列语义有缺陷（取自已比对
> BAM 记录数而非输入片段数，23 行出现 mapped>total），该列及其派生比值不可用；
> 详见 `results/PERF_NOTES.md`。

| 数据 | salmon | minibwa |
|---|---|---|
| D1（43.4M reads） | 127.4 s / 233.6 MB | 384.9 s / **2,474 MB** |
| D3（49.1M reads） | 101.8 s / 563.5 MB | 453.2 s / **7,875 MB** |
| D6（31.3M×50bp） | 35.9 s / 619.4 MB | **274.0 s（正常完成）** / 8,783.7 MB |
| D8（13 万对） | 8.1 s / 277.9 MB | 14.4 s / 297.9 MB |
| D11–D18（人工集） | 8.7–10.9 s / 234 MB | 9.3–12.3 s / 236–372 MB |
| Bari S3/S9（k15 已缓存） | 9.5 / 8.7 s | 拒绝执行（小RNA） |
| Bari S3/S9（srna 组装+BLAST） | 76.7 / 71.7 s | — |
| Bari 全链路（fastp→鉴定+组装） | 69.2 / 74.1 s | — |
| Bari S3/S9（salmon-k31 对照小库 / minibwa 对照） | 0.2 / 0.2 s（0 映射） | 0.0 / 0.0 s（拒绝，rejected） |

D5 的 minibwa 行（0.9 s / 2.0 MB）为异常作废记录，不列表比较（见 §2.1）。

结论：**常规 reads 上 salmon 在耗时与内存上全面占优**——配对耗时比值
1.8–7.6×（**中位数 3.4×**；小集 D8 仅 1.8×，大 PE 集最高 7.6×）；
峰值内存差距随输入增大而放大：D1/D2/D3/D6 低 **10.6–14.2×**
（D1 233.6 MB vs 2,474 MB；D3 563.5 MB vs 7,875 MB；D6 619.4 MB vs 8,783.7 MB），
D4/D9/D10 低 3.2–5.6×，小集两引擎相当（D8 1.1×）；
灵敏度持平或更高。minibwa 的价值在真比对产物（CIGAR/NM/ANI、共识支持），
但其定量+汇总路径对大样本的内存（最高 8.8 GB）与耗时需要限额保护（本评测 20 分钟，
本轮未触发）。

## 5. 全链路 sRNA 实测（fastp→鉴定→组装）

| 样品 | fastp | 鉴定(k15) | 组装 | QC 存活 reads | 确诊（stage1 表行数） |
|---|---|---|---|---|---|
| S3 | 4.6 s | 8.5 s | 55.6 s | **249,802/250,000**（99.9%） | 3 行（含 1 行重复 PVB） |
| S9 | 4.5 s | 8.8 s | 60.1 s | 249,889/250,000 | 8 行 |

（确诊列 2026-10-03 按 `results/chain_results.json` 的 stage1_confirmed 行数改为
"3 行 / 8 行"——原稿 3 / 9 与归档数据不符，疑为早期一轮链运行的旧值。）

fastp `min_len` 已自适应小RNA（默认 15，`fastp_min_len` 可配），小RNA reads
几乎无损通过质控（仅 ~0.05% 丢弃）；全链路端到端约 70 秒/样品。

### 5.1 工具链全流程实测（图 fig_b10，run_chain_bench.py）

四段工具链在全部 20 个样本上端到端串测（①kvsuite 识别+双轨过滤+抽病毒 reads →
②SPAdes 组装 → ③kunpeng 再鉴定 + 谱系表 → ④blastx+CDD 验证）。
Bari sRNA 用 srna 模式 + 60bp；VIROMOCK D1–D18 用 rnaviral + 200bp（minibwa 双引擎
已在 §2 验证，此处 salmon 单引擎跑链）：

| 样本 | ① 识别提取 | ② 组装 | ③ 再鉴定 | ④ 验证 | 合计(s) | 链级检出 |
|---|---|---|---|---|---|---|
| S3 | 13.6s | 29.3s | 16.8s | 0.8s | 60.5 | **2/2 = 100%**（PVX+PVB；28/28 contig 判病毒） |
| S9 | 12.6s | 24.4s | 0.4s | 1.4s | 38.8 | 4/9 = 44.4%（GLRaV-1/GRSPaV/GRVFV/GYSVd1） |
| D1 Citrus | 250.2s | 6.1s | 20.2s | 360.7s | 637.2 | **5/5 = 100%** |
| D2 Citrus | 287.5s | 7.1s | 0.3s | 376.4s | 671.3 | **5/5 = 100%** |
| D3 Grapevine | 181.0s | 11.5s | 0.3s | 339.5s | 532.3 | 3/5 = 60% |
| D4 Grapevine | 48.6s | 5.6s | 0.3s | 321.4s | 375.9 | **4/5 = 80%** |
| D5 Potato | 62.1s | 6.2s | 0.3s | 319.2s | 387.8 | **1/1 = 100%** |
| D6 Potato | 59.0s | 6.0s | 0.3s | 313.4s | 378.7 | **1/1 = 100%** |
| D7 Tobacco | 36.0s | 2.8s | 0.3s | 0.7s | 39.8 | 0/1（TSWV 库缺） |
| D8 Chenopodium | 11.0s | 29.0s | 0.2s | 0.9s | 41.1 | 1/2（PFBV 检出；CqMV1 0.5% 低浓度） |
| D9 Pistachio | 45.0s | 6.3s | 0.3s | 313.3s | 364.9 | **1/1 = 100%** |
| D10 Prunus | 52.4s | 13.2s | 0.3s | 314.7s | 380.6 | **2/2 = 100%** |
| D11–D18 (8) | 9.6–13.7s | 4.5–8.0s | 0.2–0.7s | **326.9–431.1s** | 347–450 | **各 1/1 = 100%**（8 种病毒全检出） |
| **合计** | | | | | | D1–D10 链路 **23/28 = 82.1%**；含 D11–D18 的 VIROMOCK 全量 **31/36 = 86.1%**；全 20 样本无运行失败 |

要点：①链级检出（①确诊 ∪ ③再鉴定）与单路线一致——S3 100%、S9 受低丰度病原与
默认确诊阈值限制；②S3 全部 28 条 contig 均判为病毒（无宿主污染混入）；
③S9 抽取的病毒 reads 分散在 9 个病原上，≥60bp contig 仅 2 条——组装路线对
低丰度多样品场景覆盖有限（与 §6.1 门槛扫描结论一致：瓶颈在组装而非门槛）；
④整链峰值内存 2.2 GB。过程中新增两处平台加固：srna 组装对小子集输入的
k-mer 覆盖模型崩溃自动降 k 阶梯重试；samtools sort 显式内存限幅。
说明：D1–D10 的 ④验证段含联网/数据库检索型步骤（每样本 313–431 s），
D11–D18 的 ④验证同为 326.9–431.1 s——这是链时长的主要构成；D7/D8 的 ④为
秒级（无候选或跳过）。原稿"VIROMOCK 24/28 = 85.7%"与逐行数据不符
（D1–D10 实为 23/28），已按 `results/chain_results.json` 更正（2026-10-03）。

## 6. 共识序列分析比较（图 fig_b8，run_cons_bench.py）

**实验口径披露（2026-10-03）**：①"模拟自洽"路线为**无错模拟**（150bp PE、
平铺、**无测序错误、无污染**）——所得"零引入错误"是**自洽性测试**，只证明管线
自身不引入变异，不代表真实数据（含测序错误/低质量碱基）的处理能力；
②三路实验均**显式指定 vc_depth=5**（vc_qual=20 / vc_freq=0.5）；**评测当时
（2026-09）平台默认 vc_depth 实为 1**（帮助文本误写 5，2026-10 已修正默认为 5），
故本结果不代表当时平台默认口径；**平台默认 vc_depth=1（旧默认）下的表现未测**。

三路对比（Bari S3 / 平台共识段 = minibwa 真比对 + viral_consensus）：

| 路线 | 目标 | 覆盖率 | 一致率 | 说明 |
|---|---|---|---|---|
| **模拟自洽**（库内参考→模拟 150bp PE 无错误 reads 200×→平台共识） | PVB RNA1 (7,148bp) | **99.82%** | **100.0%**（0 错 0 缺口） | 管线自身零引入错误，覆盖损失仅端点 |
| 同上 | PVB RNA2 (4,527bp) | **99.69%** | **100.0%** | 同上 |
| 同上 | GLRaV-1 (18,659bp) | **99.93%** | **100.0%** | 同上 |
| **sRNA 实测**（S3 真实 21–24nt reads → PVB 两段） | PVB RNA1/RNA2 | **0%** | 产物为 **100% N** | 工具**未拒绝**输入：minibwa stderr 自报"mapped 5,613,132 bp in 250,000 sequences"（恰等于输入 reads 总长度），而产出的 BAM 实测 **250,000 条记录全为 unmapped**（`samtools idxstats`：两参考各 0 mapped）——stderr 与 BAM 互相矛盾；viral_consensus 消费 0 映射 BAM 输出参考长度全 N 序列且**不报错**（静默假产物） |
| **组装替代路线**（S3 srna 组装 contigs vs PVB 两段） | PVB RNA1 / RNA2 | **14.7% / 11.8%** | **100%** | siRNA 组装 contig 的高一致片段，可作 sRNA 共识的低覆盖替代 |

**sRNA 失败根因（2026-10-03 按证据改写）**：原稿写"minibwa 拒绝 21–24nt"，
与实测不符——共识段内部固定 minibwa 自建索引重比对，**工具并未报拒绝**：stderr
自报已映射（bp 数=输入全长），BAM 却 0 mapped，二者矛盾，随后产出全 N 共识。
这是"0 映射无守卫"的静默失败（与定量段入口的 `rejected=true` 明确拒绝是两种不同
表现）。**平台侧核查**：`Virus_Platform_Core/known_virus_suite/kv_consensus.py`
现无 0 映射守卫（2026-10-03 grep 核实）——该守卫属 kvsuite 修复范围
（`docs/audit_20261003_修复对照表/02-kvsuite.md` 经查**未含此项 → 待处理**），
本文档不跨子系统改动平台代码，仅在此登记"0 映射守卫为平台侧待办"。

结论：①共识管线在常规 reads 上**零引入错误**（200× 自洽 0 SNP/0 indel），可放心用于
提交级共识；②**小RNA 数据的共识段会静默产出全 N 序列**——与"定量全零"同类的
静默失败，应在共识段加 0 映射守卫（拒绝或告警）；③sRNA 场景的共识替代以 siRNA
组装 contigs 一致率 100%（PVB 覆盖 12–15%）。

### 6.1 min contig 门槛扫描（图 fig_b9，run_contig_sweep.py）

同一套 S3 组装 contigs（48 条）按 30–200bp 八档阈值过滤后对 PVB 两段评估：

| 发现 | 数据 |
|---|---|
| 放宽门槛的**覆盖收益有限** | RNA1 全档恒为 14.7%；RNA2 仅从 7.0%（200bp）→ 11.8%（≤40bp），即 +4.8 个百分点 |
| **无一致率代价** | 全部门槛下命中一致率均为 100%，无错误片段混入 |
| 门槛敏感区在 150–200bp | 低于 150bp 后 contig 数稳定（36→46），覆盖不再增长 |

修正此前的推断：**共识覆盖的瓶颈在组装产出的 contigs 本身（未去宿主、株系分歧、
siRNA 片段天然短），而不是过滤门槛**。min contig 下调的主要收益在**物种检出**
（S9 案例：200bp 门槛丢失 HSVd/GYSVd1，60bp 可保），且无一致率代价——
建议 sRNA 场景默认 60bp（与 Massart 建议一致），但不要期待它提高共识覆盖率；
共识覆盖要提升需从组装端入手（去宿主后组装、k 值/模式调优）。

## 7. 数据库与判读发现

1. **ICTV 改名是自动判读的头号陷阱**：本轮评测 6 个物种需按新名匹配才算正确
   （清单见 §1.2）。平台确诊表带 NCBI/ICTV 双物种列，自动比对需双列 + 别名表。
2. **真实数据库缺口**：TSWV（D7）库内 0 条参考，双引擎必然漏检——**建议优先补录
   Orthotospovirus 属**〔时效标注：该建议针对 2026-09 快照；2026-09-19 起
   `plant_ref` 已注入 TSWV L/M/S，**TSWV 现已入库**〕；PiVB（D9）缺失且命中的
   *Pistacia emaravirus* 参考 100% 覆盖疑似误标，建议核对来源。
3. **株级 accession 缺失是设计使然**（物种 curated 库）：VIROMOCK 38 个株级
   accession 不在库内，物种级检出不受影响。
4. **假阳性（FP）构成：近缘溢出只是其中一类（2026-10-03 按数据改写）**。
   D1–D10 双引擎合计 81 个 FP 实例、22 个物种（`results.json` 重新统计）：
   - **同属/同科近缘溢出**：*Physalis floridana potyvirus*、*Spiraea potyvirus*
     （PVY 集，各 4 次）等——原稿"FP 全为同属/同科近缘"**不成立**；
   - **跨背景复发的待查实体**：*Citrus exocortis Yucatan viroid* 在 **21 个实验**
     的确诊表里复发（含 D5/D6 马铃薯、D7 烟草等非柑橘背景——柑橘类病毒在生物学上
     不可能出现在马铃薯/烟草中），*Dapple peach fruit disease viroid* /
     *Dapple plum and peach fruit disease viroid*（合计 14 次）、
     *Citrus gummy bark viroid*（7 次）同样跨背景出现。此类物种（约 300nt 类病毒，
     种间同源性高），**跨背景复发原因待查**（索引/污染/参考交叉映射三选一，需专项
     排查）。注意：CEVd-Yucatan 并非纯误报——官方 Dataset1.md 的实验室记录中
     USDA-APHIS 亦报告检出该实体（FJ751933.1）并注明"与 CEVd 不同"，说明挑战赛
     真实材料中确有记录；但平台在马铃薯/烟草等背景的命中无法据此解释，仍为待查；
   - **不相关物种命中**：*Prunus virus F*（6 次）、*Apple ourmia-like virus 1*
     （4 次）、*Elephantopus scaber closterovirus* 等与数据集真值无近缘关系，
     属门槛（cov/depth/poisson）与 ANI 分流松紧问题。
   结论：FP 由**近缘溢出 + 跨背景复发实体（待查）+ 无关命中**三类构成，
   门槛与 ANI 分流的松紧应随诊断场景调整；跨背景复发需专项排查（**待查项**）。
5. **官方"Real but unconfirmed"条目的从紧口径**（2026-10-03 补引官方文档）：
   官方 `Dataset1.md` 表 1 将 *Citrus endogenous pararetrovirus* (CitPRV)、
   *Citrus dwarfing viroid* (CDVd)、*Apple hammerhead viroid-like RNA* (AHVd RNA)
   3 个条目列为 "Real but unconfirmed"（脚注¹："There are a few reads of this
   virus/viroid present, but these could not be confirmed by qPCR and are most
   probably the result of cross-contamination"）。本评测**按从紧口径**：未确诊
   条目不列入真值，平台若检出即计 FP——D1/D2 的 FP 表里 CitPRV（2 次）、CDVd
   （4 次）正来自该口径；若改用松口径（计为真背景），灵敏度分母与 FP 口径都会变化。
   该口径差异应随诊断场景（监管从紧 vs 研究从松）显式选择。

## 8. 评测过程发现并已修复的平台缺陷

| 缺陷 | 影响 | 修复 |
|---|---|---|
| `kv_identify.summarize` 用 polars 默认 schema 推断（前 100 行） | **minibwa 引擎**：Avg_Read_ANI 浮点未出现在前 100 行时汇总崩溃，确诊表静默为空（D11–D18 全军覆没即此因） | `infer_schema_length=None` 全量扫描（kv_identify.py） |
| fastp JSON 报告含未转义 Windows 路径 | 评测 chain 实验解析崩溃 | 改数 fastq 记录数（评测侧） |
| 评测框架：Windows 孙进程继承管道致挂起；页面文件耗尽打死监控 | 评测自身 | 文件重定向 + os._exit + 20 分钟限时 + OSError 容错 + 结果提取后自动清理中间文件（2026-10-03 改为先保留关键产物到 `results/artifacts/`，不再销毁原始证据） |

## 9. 局限

- 半人工集 FP 判定含官方 "Real but unconfirmed" 背景条目的从紧口径（见 §7.5）；
  FP 的"近缘溢出 vs 跨背景复发（待查）vs 无关命中"三分需专家判读与专项排查
  （CEVd-Yucatan 跨背景复发为待查项）；
- Bari 参试实验室基线（76%）含宿主过滤与自定阈值，与平台默认口径不同，仅作量级对照；
- 定量对比存在**单位口径差异**：PE 数据集真值为 reads（=2×fragments），平台 EM
  计数为片段——log-log 拟合整体下移（后验核算：截距约 −0.54～−0.56、比值中位数
  0.40～0.43），**r=0.968 只表征相关、不表征准确度**（详见 §3 与 fig_b4 注记）；
- **D5 集 minibwa 记录作废**（运行异常提前退出，无 mapped/total），其灵敏度不计入
  合计；**D6 minibwa 正常完成**（274 s，1/1）；归档记录中无任何 `timed_out` 标记；
- 真值本身的已知缺口（不虚构补值）：`truth_viromock.json` 的 D18 组成合计 85.27%
  （约 14.7% reads 未见于源表）、D14 合计 105.5%（源表问题）；BARI_TRUTH 的 S9
  缺 *Grapevine virus A* 的 reads 数（定量点不含该物种，检出判定不受影响）；
  D1–D10 分组真值见 `truth_viromock_d110_groups.json`（与 `truth_viromock_d110.json`
  物种集合一致，别名更全——含 ICTV 改名，评估采用前者）；
- 基准数字全部绑定 2026-09 的 `kv_index` 快照库（现已被 plant_ref 取代、未随仓库
  归档）——**未重建等价库前，本文不可原位复跑**。

## 10. 复现

```bash
cd docs/paper/benchmark
python run_bench.py            # 断点续跑（results.json 已有实验跳过，含 (exp,error) 去重）
python compute_and_plot.py     # 指标 CSV + figs/fig_b1..b7（输入默认取仓库内 results/）
python run_cons_bench.py       # 共识比较实验 → results/cons_results.json
python make_fig_cons.py        # figs/fig_b8_consensus
python run_contig_sweep.py     # min contig 门槛扫描 → results/contig_sweep.json
python make_fig_sweep.py       # figs/fig_b9_sweep（缺 sweep.log 时跳过并提示）
python run_chain_bench.py      # 工具链全流程 → results/chain_results.json
python make_fig_chain.py       # figs/fig_b10_chain
```

**路径参数化（2026-10-03）**：鉴定库目录优先级 `spec['db_dir']` > 环境变量
`VP_BENCH_DB` > 旧默认 `databases/virusref_db/kv_index/`（已不存在，缺库时脚本报
中文修复指引）；BLAST 库 `VP_BENCH_BLAST`；Bari 原始数据 `VP_BENCH_BARI`；
VIROMOCK/Bari 数据集根 `VP_BENCH_VM`；对照输入根（in1.fastq / idx_k31，
评测机外部产物，未随仓库归档）`VP_BENCH_CTRL`；`compute_and_plot.py` 输入目录
`VP_BENCH_RUNS`（默认仓库内 `results/`）；真值 `VP_BENCH_TRUTH` /
`VP_BENCH_D110_GROUPS`；扫描日志 `VP_BENCH_SWEEP_LOG`；输出目录
`VP_BENCH_CONS_OUT` / `VP_BENCH_SWEEP_OUT`。缺输入时的行为：库缺→清晰中文报错；
扫描日志缺→跳过并打印说明。

数据依赖：Bari 样品 + VIROMOCK D1–D18 fastq.gz（评测机外部目录
`E:/谷歌下载/测试数据验证/`，可用 `VP_BENCH_BARI` 覆盖）；鉴定库为 2026-09 的
`kv_index` 快照（**需重建等价库**，见头部"参考库版本"与 `ENVIRONMENT.md` §4）。
共识评测数据：`results/cons_results.json`；性能字段口径：`results/PERF_NOTES.md`。
真值三件（2026-10-03 起角色明确）：
- `truth_viromock.json` —— D11–D18（accession 级，含 `_notes` 缺口披露），
  脚本默认入口（compute_and_plot 可用 `VP_BENCH_TRUTH` 覆盖）；
- `truth_viromock_d110.json` —— D1–D10 **原始解析**（别名较简）；
- `truth_viromock_d110_groups.json` —— D1–D10 **评估用分组别名**（含 ICTV 改名，
  与上者物种集合一致；`compute_and_plot.py` 已改为读取本文件，脚本内常量已删除）。
方法变更：50K 深度稀释自 2026-10-03 起为 seed=42 随机抽样（旧数字基于顺序截取）。

## 修订记录（2026-10-03）

本次修订逐项核对 `results/results.json` / `metrics_summary.csv` / `perf.csv` /
`chain_results.json` 与平台代码，主要更正：①D5/D6 minibwa 结论（原文"超时强杀
记 0 分"与数据不符：D6 正常完成、D5 记录作废）；②Bari 表（6/1/5→6/4/5、
S9 全链路 40%→28.6%）；③实验数 50→46（并说明重复记录来源）；④定量单位口径
与系统偏移披露（截距/中位比值/r 的解读边界）；⑤FP 构成如实改写（近缘溢出并非
全部，CEVd-Yucatan 待查）；⑥"放宽阈值可到 8/9"撤回（依据未归档）；
⑦参考库版本与失效披露（kv_index→plant_ref）；⑧D7 TSWV 时效标注；
⑨链级合计 24/28→23/28（D1–D10），并补 D11–D18 验证段秒数；
⑩对照说明（3 参考小库 idx_k31 / minibwa 拒绝）与 ctrl_s9 结果入表；⑪性能表数字对齐 perf.csv 并补
PERF_NOTES；⑫共识实验口径（无错模拟、vc_depth=5 非当时默认）披露。
脚本侧：图 b1/b2/b4/b6 已重生成并目检通过。
