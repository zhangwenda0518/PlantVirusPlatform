# CRISPR-Design（CRISPR/Cas12a 检测试剂设计）

> 位置：**检测与防治 → CRISPR 检测试剂设计**（工具卡 `t-crispr`）
> 引擎：`Virus_Platform_Core/crispr_design.py` · CLI：`main.py crispr-design`
> 版本：1.0.0（2026-10-06）

一次产出 **Cas12a crRNA 候选 × RPA（或 PCR）引物对** 的协同设计、复合评分排序、
订购表与 L1–L5 实验验证清单。全程本地离线，复用平台已纳管的 MAFFT /
primer3 / ViennaRNA / BLAST，不引入新外部工具。

## 输入

| 输入 | 必选 | 说明 |
|---|---|---|
| 目标 FASTA | ✅ | 多分离物（建议 ≥10 株，覆盖度评估才可靠）。header 含 `RNA1`/`S2`/`segment 3` 等标签时自动按片段分组，多分体基因组各片段独立设计；无标签视为同一分子 |
| 近缘病毒 FASTA | 可选 | 同属/同科近缘种 + 常见复合侵染毒源。也可固定放入 `databases/crispr/near_neighbors/`（放即生效） |
| 宿主基因组筛查 | 默认开 | 用 host-db 当前宿主源基因组（无配置则跳过并如实注明）。FASTA ≤1.5GB 走 k-mer 索引路（首次建索引较慢，之后有缓存）；更大走 BLAST 路 |

目标序列可用「比较基因组分析 → 参考序列获取」按 ICTV 科/属一键拉取。

## 流程（设计文档见 `docs/CRISPR检测设计模块规划_20261006.md`）

1. **分组与比对**：按片段标签分组 → 组内多序列 MAFFT → 逐列保守性 +
   多数共识（平票记 N 不猜）；
2. **crRNA 漏斗**：双链 TTTV PAM 扫描（严格档；放松档前 3 位容 1 错配，
   需实验验证）→ 过滤：spacer 平均保守分（严格 0.95/放松 0.85）、GC 30–70%、
   同聚物 ≥5 剔除、poly-T≥4 标记 T7 风险、RNAfold MFE ≥ −5 kcal/mol；
   PAM 优先级 TTTA/TTTC/TTTG > TTTT > 容错 PAM（LbCas12a 活性口径）；
3. **引物设计**：每个存活位点取侧翼 ±250bp 窗口，primer3 设计 RPA
   （30–35nt / Tm 60–70 / 产物 100–200bp）或 PCR 引物对；primer3 出对
   不足时自动枚举兜底。**扩增子必须完整包含 PAM+spacer 且两侧各留 ≥25bp**；
4. **双库特异性**：crRNA spacer / 引物对宿主库 + 近缘库逐条筛查
   （20-mer，mm0/mm1/mm2 分层）——spacer 精确命中一票否决；≤1 错配标记
   「近缘≤1mm」；引物精确命中整对剔除、3′ 端 ≤1 错配重罚；
5. **复合评分 × 群体覆盖率**：21 项特征效用加权（脱靶/二聚体/3′ 端/Tm 平衡/
   GC/种子区/结构，初始权重来自 CASPER 拟合 100 组实验对的公开配置），
   总分再乘覆盖率因子（mm0 覆盖 100% 不折减、0% 折半——覆盖不了目标群体
   的序列特征再好也不是可用试剂，PathoGD prevalence 口径）；
6. **输出**（落 `run/tool_runs/<run>/crispr/`）：
   `ranked_sets.tsv`（top 组合全特征）、`crrna_candidates.tsv`（漏斗全程
   明细含被否候选）、`order_sheet.xlsx`（订购表：spacer / DR-spacer /
   T7 体外转录模板三种 crRNA 格式 + 引物）、`validation_checklist.md`
   （L1–L5 植物病毒语境验证清单）、`report.json/html`、`crispr_design.json`
   （参数快照，可复现）。

## 效应酶（NucleaseProfile，2026-10-07 二期）

引擎的酶参数已抽象为档案对象（借鉴 crisprDesign `CrisprNuclease`
的轻量净室版）：PAM/PFS、spacer 长度区间、DR 序列、seed 长度、poly-T
口径全部随档案切换，新增酶 = 新增档案条目。

| 档案 | 适用 | 关键差异 |
|---|---|---|
| **LbCas12a**（默认） | DNA 靶标（RPA/PCR 扩增后检测） | TTTV PAM；spacer 20–24；poly-T≥4 仅标记 |
| **LwaCas13a** | RNA 病毒直接检测（免逆转录路线） | 无 PAM（保守窗口即位点）；spacer 27–30；**poly-T≥4 硬否决**（crRNA 一律 T7 体外转录，36 nt DR：Kellner et al. 2017 SHERLOCK 口径） |

Cas13a 模式的订购表只出 T7 体外转录模板（T7 启动子 + DR + spacer）；
对 RNA 植物病毒可走「粗提液 → RT-free 检测」，需要信号放大时仍配
RT-RPA 引物（设计逻辑同 Cas12a，位点=spacer 本身）。

## 二期能力（2026-10-07 实装）

- **EXO/NFO 探针**（RPA 模式可选）：在扩增子内部（引物结合区之间，
  留 5bp 边距）扫描 30–35nt 窗口，纯 ACGT、GC 30–70、同聚物 <4、
  达保守档阈值、Tm≥60℃，尽量居中。EXO = 5′-FAM—THF(中位)—BHQ1-3′
  （ExoIII 实时荧光）；NFO = 5′-FAM—THF—C3-3′（EndoIV 侧流，
  **需配套一端 5′ 生物素引物**，订购表已注明）。空间不足时该组合标记
  跳过原因。订购表与组合表自动带修饰标记序列。
- **intentional mismatch 建议**：近缘库 ≤1–2 错配命中的候选，自动定位
  天然差异位，并建议在 PAM 远端（避开 seed）引入人工错配实现差分去
  稳定（靶标 1 个人工错配 vs 近缘 ≥2）。序列与位置依赖，**须实验验证
  区分度增益与活性损失**（综述 §guide engineering 口径，建议随报告展示）。
- **跨片段 panel 汇总**：`panel.tsv` + 报告/订购表 panel 页——每个片段
  取评分最高的一套组合，多分体病毒多联检视图。
- **IUPAC 兼并码共识**：平票列不再记 N 而是记兼并码（W/S/M/K/R/Y），
  引物模板区不再被 N 打洞；crRNA 与引物交付序列仍有纯 ACGT 硬守卫。

## 对照标准操作清单（2026-10-07 补齐后全覆盖）

| 操作清单要求 | 平台实现 |
|---|---|
| MSA 保守区（CP/MP/RdRp，多分离株） | MAFFT 逐列保守性 + 多数共识（参考序列可直接从「参考序列获取」按 ICTV 科/属拉取） |
| Cas12a TTTV PAM + 下游 20–24nt | 双链 PAM 扫描（严格/容错档），spacer 20–24 |
| Cas13a 无 PAM，28–30nt 直接选 | NucleaseProfile cas13a 档案，保守窗口即位点，27–30 |
| **Cas13a DNA 体系前向引物加 T7 启动子** | cas13a 模式自动输出 forward+T7（订购表/详情/报告） |
| spacer GC 40–60% | 漏斗硬界 30–70%（宽），评分钟形最优 40–60 |
| 避免与 DR 骨架形成稳定结构/自二聚 | spacer 单独 MFE≥−5 **+ DR+spacer 联合折叠 ΔMFE 闯门**（全长比 DR 单折低 >3 kcal 判互作剔除） |
| 种子区（PAM 近端 8–10nt）绝对保守 | seed 区逐列保守硬约束（严格档 1.0 / 放松 0.95），漏斗独立统计 |
| RPA 产物 100–250bp、位点距引物 ≥10–20bp | 严格档 100–200（放松 90–300），位点边距默认 25bp |
| 引物 30–35nt、GC 30–70% | primer3 主路 + 枚举兜底同口径 |
| 3′ 端避 GC 富集/回文 | 末三碱基全 G/C 一票否决（primer3 后置 + 枚举路 + 探针），3′ 同碱基 run 与交叉二聚体进评分 |
| 多 crRNA × 多引物对候选 | 每位点 ≤2 组合、top_n 输出、21 特征评分排序 + 覆盖率因子 |
| 体外荧光报告正交筛选 | L1 清单（aLoD/time-to-positive/多 crRNA 并排） |
| Cas12a/13a crRNA 结构（DR+spacer） | 订购表三格式：spacer / DR+spacer / T7 转录模板（DR 取自档案：Cas12a UAAUUUCUACUAAGUGUAGAU；Cas13a GAUUUAGACUACCCCAAAAACGAAGGGGACUAAAAC） |

## 生产验证记录（2026-10-07 真实数据冒烟）

15 株真实 PVY 完整基因组（9.7kb，NCBI）+ 6 株 PVA 近缘 + 番茄真实
宿主基因组（810MB）全配置跑通，Cas12a 14s / Cas13a 14s。冒烟抓出并
修复三个只在真实数据下暴露的问题：

1. **primer3 不认 IUPAC 兼并码**：真实 15 株在若干列平票出 M/K/W/Y，
   主路整批报错降级——模板送 primer3 前兼并位转 N（PRIMER_MAX_NS=0
   自动阻止引物落位，兼并码保模板区长度的本意不变）；
2. **大宿主 k-mer 索引内存爆**：810MB 基因组 np.unique 峰值 ~40GB，
   本机 21GB 可用必然失败且曾静默跳过——k-mer 上限降到 250MB，失败/
   超限一律降级 BLAST 路（blastn-short 批查），绝不静默放弃筛查；
3. **BLAST+ LMDB 中文路径崩溃**：平台根含中文时 mdb_env_open 失败
   ——BLAST 库改建在 %TEMP% ASCII 中转目录（_ascii_work_base 同
   MAFFT 方案），host.fa 按大小比对增量重拷。

真实数据下各门控的工作量：Cas12a 309 位点（保守-182/种子-107/DR互作
-1）→ 17 存活；Cas13a 15,686 窗口（保守-9,168/种子-5,631/polyT-74/
DR互作-233）→ 384 存活。宿主 BLAST 筛查在 17 个候选中标记 8 个
≤1mm 背景命中。番茄 BLAST 库 810MB 实建 6.4s。

## 规模化生产验证（2026-10-07 第二轮）

**分发冒烟**：dist/VirusPlatform-Source（源码分发形态）——① 分发自带
解释器（3rd/python，含 primer3/RNA/openpyxl）跑主树真实数据 11.8s ✓；
② 叠加用户形态（dist 独立目录 + platform.json 指外部工具路径）10.6s ✓；
③ MAFFT 缺失时清晰报错"未找到工具 [mafft]，请在 platform.json 或 GUI
设置页指定路径" ✓；④ 最新代码已同步进 dist 快照（164 文件）。

**四病毒 × 双酶真实数据矩阵**（全部番茄真实宿主 + BLAST/k-mer 双路 +
真实近缘属内种 + EXO/NFO 探针）：

| 病毒 | 基因组 | 耗时 | 关键验证点 |
|---|---|---|---|
| PVY（马铃薯 Y 病毒） | 单体 9.7kb ×15 株 | 14s | 全配置基线；IUPAC→primer3 修复后 route 恢复 |
| ToMV（番茄花叶病毒） | 单体 6.4kb ×8 株 | 37s | **近缘 TMV（同属 ~80% 同源）veto 9 候选**——交叉反应高危区实证 |
| CMV（黄瓜花叶病毒） | **三分体** RNA1/2/3 ×7 株 | 45s / 164s（Cas13a） | 多分体分组、panel 按段取优；RNA2 亚组分裂覆盖 57% 被正确降权；Cas13a×多分体组合路径 |
| TSWV（番茄斑萎病毒） | **字母式 L/M/S 分段** | 121s | 字母分段命名支持；INSV（同属）veto 7 候选 |

第二轮修掉的真实数据 bug：**多分体分组失效**（真实 NCBI 标题的
`segment RNA1` 在标题中段而非首词——记录名保留完整 header）；**字母式
L/M/S 分段不识别**（TSWV 命名）；同时修掉字母分支对 isolate 编号
（A-25）的误伤。UX：任务 stage 全中文化 + UI 契约测试
（tests/test_crispr_ui.py 锁定卡片控件/i18n/渲染分支）。

真实数据集与产物保留在 `run/tests_crispr/real/`（PVY/CMV/TSWV/ToMV
四组 fasta + out_* 七个运行目录），可直接打开 report.html 复核。

## 参数档

| 档 | 适用 | 关键差异 |
|---|---|---|
| **严格**（默认） | 低变异病毒 / 分型精确 | 保守 0.95、仅 TTTV、产物 100–200bp |
| **放松** | 高变异 RNA 病毒 | 保守 0.85、容错 PAM、产物 90–300bp（PathoGD 黄病毒实测先例） |

## ⚠ 证据边界（ACS Sensors 2026 综述口径）

计算设计**仅为序列级候选优先级排序**，不构成 crRNA 活性、错配区分度、
扩增兼容性或基质性能的证据。实验结论必须按清单逐级建立：

- **L1** 纯化靶标：合成靶标测 aLoD / time-to-positive；
- **L2** 加标基质：健康植株粗提液（纸片/碱裂解）加标 + **近缘病毒交叉反应**；
- **L3** 回顾性田间样本（平台 02_kvsuite 阳性可召回，RT-qPCR 对照）；
- **L4** 前瞻性多地点田间（非开发者操作）；
- **L5** 实际使用（储存稳定性 / 批次一致性 / 单样成本）。

## 结果可视化

三层可视化（2026-10-06 实装，经 visual-jude 五轮目检交付）：

1. **工具卡结果区**（t-crispr）：位点轨道图（每片段一条 SVG 泳道——保守性
   折线 + 漏斗存活位点三角标记（红=入选/灰=未入选/空心=已否决）+ top 组合
   扩增子条带）；点击表格行展开**评分分解面板**（7 类效用条形 + 三种 crRNA
   订购序列与 T7 提示）。
2. **report.html**（独立交互报告，plotly）：① 各片段保守性剖面与位点
   散点；② Top 组合评分构成堆叠条形（条形总长 = 基础分/100，悬浮含各类
   权重）；③ 评分 vs 群体覆盖率散点（x 轴按数据自适应）；④ 明细表。
   通过平台服务打开自动加载 plotly；本地双击打开的副本表格仍可读。
3. **order_sheet.xlsx / ranked_sets.tsv / crrna_candidates.tsv**：订购与
   归档级明细（含被否决候选与否决原因）。

## CLI

```bash
python main.py crispr-design -i targets.fa --near nearvirus.fa \
    --amp RPA --profile relaxed --top-n 20
python main.py crispr-design -i rna_virus.fa --nuclease cas13a \
    --profile strict --top-n 20
python main.py crispr-design -i targets.fa --no-host -o run/my_crispr
python main.py crispr-design -i targets.fa --probe exo --top-n 10
```

## 引用

设计输出报告自动附带：Wang H. et al. *ACS Sensors* 2026（框架与 L1–L5）；
Low S.J. et al. *Commun Biol* 2025 PathoGD（群体覆盖口径）；CASPER
（初始评分权重）；AutoCORDSv2 / RPA-Primer-Design / CRISPRPrimerDesigner
（规则与参数参考，全部净室实现，未复制任何 GPL/无许可代码）。
