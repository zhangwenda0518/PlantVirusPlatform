# VirusPlatform 脚本审查报告：审稿人视角

**审查日期：** 2026-10-10  
**对象：** 当前工作树中的样品主流水线、病毒鉴定/组装/验证、进化与群体遗传分析、公共数据检索、基准脚本及论文 Methods。排除 `3rd/` 第三方源码、生成产物本身与归档快照；只在需要时看调用方和测试入口。  
**结论标记：** ✅ 已实测（本轮可复跑）／推断（源码与调用链证据，未运行真实样品或外部工具）。

> **工作树状态说明：** 审查开始时 `git status` 已显示大量已修改、已删除、已暂存和未跟踪文件，包含核心脚本、论文、benchmark、样例结果与工具文件。本报告只针对当时可见的工作树，不把这些改动归因于本轮审查；本轮没有修改业务代码或样例数据。除下文单项隔离探针外，未运行会重写 `examples/`、`run/` 的测试套件。

## 0. 总体评价

**审稿建议：大修后再接受任何“鉴定能力/新病毒发现/诊断性能”层面的结论。** 作为工程平台论文，宜继续限定为实现与可复现性展示；当前稿件已在 Methods/Discussion 中承认尚无独立生物学验证，这一点是正确的。

平台的工程基础并不薄弱：主流程有明确阶段登记；多个统计步骤保留阈值、警告和结构化产物；kvsuite 会记录失败样本而非一律当阴性；TreeTime/DRT 有固定 seed、经验 p 的 +1 校正和低功效提示。主要问题集中在四个审稿人会卡住的边界：

1. **旧产物可能越过“本次运行”边界，造成假阳性；**
2. **不少 checkpoint 只凭 `.done` 或不完整指纹复用结果；**
3. **候选序列在进入验证前已受分类器过滤，限制新病毒发现的可解释范围；**
4. **基准真值、失败判定与论文参数口径还不够严格。**

## 1. P0 — 可能直接改变结论或静默失效

### P0-1. 新输入零命中时，旧阳性鉴定表可作为本次阳性结果返回 ✅ 已实测（隔离合成探针）

- `Virus_Platform_Core/known_virus_suite/kv_identify.py:237-240`：`all_rows` 为空时只返回空 DataFrame，不重写或清除此前的 `all_viruses.summary.tsv` / `all_viruses.best.summary.tsv`。
- `Virus_Platform_Core/known_virus_suite/kv_filter.py:148-167`：只有非空表才写 `filtered.tsv` / `discarded.tsv`；磁盘过滤时又优先读旧 `best.summary.tsv`，不存在才读全量表。
- `Virus_Platform_Core/kv_stage.py:592-597,672-680`：外层每次依次运行 identify、filter，再从磁盘汇总结果；即使本轮 FASTQ 已变化，最终 `summary.json` 仍会写入本轮输入指纹。

**触发：** 同一样品输出目录先有阳性结果，随后换成无命中输入或本轮没有可确认命中。识别脚本虽重新运行，空结果没有覆盖旧文件；过滤/汇总仍可能读到旧表，之后还把它与新输入指纹绑定。**影响：** 页面与摘要可能把上次阳性误当成本次检出，属于可以改变生物学结论的错误。

**本轮隔离探针：** 在系统临时目录放入一行旧阳性 TSV，调用当前 `IdentifyStage.summarize([])` 后再调用 `_summarize()`，结果为 `old_best_file_still_exists=True`、`species_detected=['Pistachio virus B']`、`passed_rows=1`。探针未触碰仓库样例或 `run/`。

**建议：** 每次重跑先将该阶段产物写入新的运行子目录，成功后原子替换；或确保所有“空结果”也写带表头的空 TSV 并清除过期通过表。汇总时校验每个产物的运行 ID/输入指纹，缺失或不一致即报错，不能降级读取旧表。

## 2. P1 — 可靠性、筛选偏倚与推断边界

### P1-1. 多处断点复用未验证全部输入和有效参数【源码静态推断】

- 子采样：`Virus_Platform_Core/pipeline.py:339-347` 仅比较旧 `n_pairs` 与当前采样数，并检查 `sub_R1/R2` 存在；不比较原始 FASTQ 指纹。换源 FASTQ 但 `n_pairs` 不变时可复用旧子样本。
- Fastp：`Virus_Platform_Core/preprocess.py:184-193,207-220,234-237` 只比较/保存源文件指纹；`dedup`、`min_len` 等处理参数未进入指纹。参数变了但输入路径未变，仍可能跳过重做。
- 宿主去除：`Virus_Platform_Core/host_removal.py:359-363` 仅凭 `.done` 即返回旧 stats；没有在该分支核验输入 reads、宿主数据库或 `confidence/min_hit_groups`。
- 组装：`Virus_Platform_Core/assembly.py:702-707` 在检查 `.done` 后立即返回，绕过后续 `715-734` 已实现的输入/组装参数指纹验证。
- kvsuite 外层：`Virus_Platform_Core/kv_stage.py:537-542` 的快速复用指纹仅含 R1/R2/参考 FASTA，不含 `ref_info`、engine、过滤参数、覆盖后端/工具版本及 `need_viral_reads`；内层批次虽然有更细指纹，但外层命中时不会进入内层。

**影响范围：** 同名样品重跑、改参数、切换数据库注释、更新输入或算法后，如果没有手动强制重跑，可能悄悄混用前一轮结果。这里不是“没有指纹”，而是 `.done` 快速返回绕过已有指纹，或指纹未覆盖实际影响结果的参数。

**建议：** 用统一的 stage fingerprint（输入内容哈希或稳健元数据 + 数据库/注释版本 + 全部生效参数 + 工具版本）作为 `.done` 的判定条件；摘要记录 fingerprint 与 run ID。至少增加“只改变一个输入/参数/数据库时必须重算”的守卫。

### P1-2. 子采样是按文件顺序取前 N 条，不具代表性保证（仅在启用子采样时）【偏倚风险推断】

`Virus_Platform_Core/pipeline.py:174-225` 的 `subsample_fastq()` 明确逐条读入并在 `n >= n_pairs` 时停止；`main.py:944-945` 也把它定义为“仅取前 N 对 reads（快速测试用）”。

**风险：** FASTQ 顺序若与 lane、测序批次、拼接文件或上机时间相关，前缀子样本可能集中来自一个子群；低丰度病毒和 reads 组成估计会受影响。平台把它限定为快速测试是合理的，但不能把前缀子样本产物当成代表性定量或最终阴性证据。

**建议：** 保留“前缀模式”作为快速冒烟选项；另提供固定 seed 的随机配对抽样（R1/R2 同步），记录 seed、抽取数、输入哈希，最终分析默认使用全量 reads 或随机抽样，并在摘要醒目标记抽样状态。

### P1-3. “新病毒候选”输入在组装阶段已被分类器门控【代码路径已核实；漏检影响推断】

- `Virus_Platform_Core/assembly.py:808-814`：`viral_contigs.fasta` 只收 `flag == 'C' and taxid > 0` 的 contig。
- `Virus_Platform_Core/pipeline.py:521-545`：③b 验证阶段只把该 `viral_contigs.fasta` 传给 blastx/BFVD 验证。

**含义：** kunpeng 标为 `U` 或没有正 taxid 的组装 contig，在进入 blastx/BFVD 之前就已排除。于是“③b `unclassified` 保留为候选池”只表示**分类器已判定候选**在后续蛋白/域证据上的状态，不等同于“所有未知 contig 都进入新病毒筛查”。这会系统性漏掉对参考分类库差异较大的候选序列。

**建议：** 将“分类已知病毒”与“未知/未分类 contig 候选”分成独立路径；至少对长度合格的未分类 contig 做蛋白域/核酸近邻筛查或采样复核，并分别报告进入、排除数量与原因。若不打算做这一步，论文应把“新病毒发现”限定为分类器先行筛选后的候选验证。

### P1-4. Salmon 缺少覆盖度后端时，双端样本的深度估算可能少算约一端【公式静态推断】

- `Virus_Platform_Core/kv_stage.py:627-634` 明确 Salmon `NumReads` 以 fragment/片段计数，PE 一对算 1。
- `Virus_Platform_Core/known_virus_suite/kv_engines.py:365-374`：samtools/pandepth 缺失时不产比对 BAM，退回 reads×读长估算。
- `Virus_Platform_Core/known_virus_suite/kv_identify.py:169-188`：只用 R1 算 `avg_len`，fallback 为 `qreads * avg_len / L`；R2 的长度未纳入。

对双端片段，估算覆盖的测序碱基应包含两端；当前公式在两端长度相近时可能低估约一半。`min_depth/min_cov` 阈值附近的参考可能因此被误筛掉。触发范围限于 coverage 后端缺失的 Salmon 路径。

**建议：** PE fallback 使用 `mean_len_R1 + mean_len_R2`，并将来源标为估算；靠近阈值时不给出“通过实测深度”的表述。补一组长度相同/不同的模拟 PE 真值对照。

### P1-5. 公共数据“检索完成”不能保证两数据源成功【源码静态推断】

- `Virus_Platform_Core/public_meta/search_engine.py:147-149`：SRA 请求异常被捕获后返回空 DataFrame。
- `search_engine.py:700-712,775-793`：合并结果为空时直接返回；主程序仍继续执行且没有将源级失败转换为非零退出码。
- `Virus_Platform_Core/web/meta.py:249-282`：任务包装层只依据子进程退出码判断成功；退出码 0 即返回 `done_msg`。

若 SRA/GSA 请求都失败，任务可显示完成但没有结果文件；若只有一个来源失败，表格可能是部分检索结果，却没有清晰的失败来源/覆盖范围字段。用户会把“暂时取数失败”误读成“无公共样本”。

**建议：** 每个来源输出 `success / no_hits / failed` 与错误详情；全部失败时返回非零退出码；部分成功则任务状态为 partial，CSV/query 快照中记录来源状态和结果完整性。

### P1-6. 默认 Salmon 路径没有 ANI，却仍把筛选通过的条目放进“confirmed”【代码与稿件对账】

- `Virus_Platform_Core/pipeline.py:442-447` 调用 `run_kvsuite_stage()` 时未传 engine；`Virus_Platform_Core/kv_stage.py:506-508` 默认 engine 为 `salmon`。
- `Virus_Platform_Core/known_virus_suite/kv_identify.py:277-287,352-363`：Salmon 无真实比对位点，`Avg_Read_ANI` 为空；代码会告警，但筛选逻辑把 `ANI is null` 与 ANI 达阈值一起放入 `confirmed`，并且 novel 列只收有 ANI 且低于阈值的行。
- 稿件 `docs/paper/manuscript_body.md:83` 表述为“probabilistic engine 下缺失 ANI 报为 undetermined”。

**影响：** 当前表格中的“confirmed”可表示“计数/覆盖度/Poisson 门槛通过”，而非该 reference 的 ANI 已核验；默认 Salmon 下疑似新变异不会进入 ANI 分流。虽然日志有告警，但“confirmed”字段与论文的“undetermined”语义不一致，容易被下游报告读成确定鉴定。

**建议：** 状态拆成 `screen_positive_ani_unavailable / ANI_confirmed / ANI_below_threshold`；未经 ANI 核验不得统一称“confirmed”。如论文需要 ANI 结论，明确运行 minibwa 或单独的真实比对复核流程，并把未测比例作为结果报告。

### P1-7. 基准脚本可能把运行失败当作有效实验，且续跑按实验名而非输入版本【源码静态推断】

- `docs/paper/benchmark/run_one.py:103-140`：保存 `returncode`，但无论退出码是否为 0，仍尝试读取 `confirmed` 表并生成 `RESULT`。
- `run_bench.py:188-200,223-238`：续跑时只用实验名 `exp` 判定 done；拿到 `RESULT` 就保存，未基于 `returncode` 设 `error`。输入、参考库、参数或脚本变更而实验名不变时，旧结果仍跳过。
- `compute_and_plot.py:157-163` 仅对 minibwa 缺 mapped/total 的特定结果作废；Salmon 非零退出不在此作废逻辑内。

**后果：** 部分产物或退出失败的运行可能被记成“无检出/低灵敏度”；如果中断后目录残留旧 TSV，还可能把上轮结果读为本轮。benchmark 的“已完成”状态也不能证明与当前脚本/输入匹配。

**建议：** `returncode != 0` 一律作废；每次执行前使用干净 scratch；结果记录含 input/index/spec/代码版本哈希，`done` 按哈希而非名称判断；另外以外部真值独立复核产物，而不只解析平台自己的输出。

### P1-8. D9 真值别名把“存疑参考命中”直接计成 TP【文档与评分代码已核实】

- `docs/paper/benchmark/truth_viromock_d110_groups.json:180-187` 把 `pistacia emaravirus` 作为 PiVB 别名，并将目标标注为“库存疑参考命中”。
- `docs/paper/benchmark/BENCHMARK_zh.md:85-88` 已写明命中来自疑似误标参考、100% 覆盖，“判 TP 但存疑”。
- `compute_and_plot.py:111-135,251-257` 对别名做双向子串匹配，命中任一别名即计真阳性。

**影响：** D9 的 TP/FN 受不确定数据库注释驱动；这会抬高该数据集灵敏度，并使汇总数字依赖一个尚未核实的别名。文档披露疑点是优点，但当前评分仍把它作为确定 TP。

**建议：** D9 单列为 unresolved/ambiguous，不进入主灵敏度分子分母，或用可靠 accession/物种判定来源核实后再锁定别名；再报告“计入/剔除 D9”两种敏感性结果。别名匹配应优先精确规范化名称或 taxid/accession，避免一般性子串误匹配。

### P1-9. 单样本 Tajima’s D 把测序深度当作有效样本量【统计有效性审稿意见；需外部验证】

`Virus_Platform_Core/known_virus_suite/kv_variant_evo.py:1156-1193` 明确在单样本 VCF（平台常态）下用窗口内变异位点平均 DP 作为 `n_eff`，再代入 Tajima’s D 样本量方差公式；调用路径可达：`known_virus_suite.py:541-544` → `kv_variant_evo.run_variant_evo()` → `plot_popgen()` (`kv_variant_evo.py:1448-1452,1788-1791`)。

**审稿关注：** read depth 是技术抽样数，不等同于独立病毒基因组/生物学个体数；PCR duplicates、扩增偏好、映射偏差与测序错误会使读段不独立。这个量可作为池化测序的探索性近似，但不能未经误差模型/独立样本验证就支持种群史或选择压力结论。

**建议：** 明确称为探索性 Pool-Seq 指标；若要推断，使用经验证的 Pool-Seq Tajima’s D 估计、纳入有效群体大小/测序误差与覆盖度不均校正，并在模拟数据与生物学重复上估计偏差/区间。

### P1-10. BEAST 树样本可脱离日志诊断独立导入【源码静态推断】

- `Virus_Platform_Core/phylogeo.py:2149-2176`：只有提供 `beast_log` 才计算 ESS/低 ESS 诊断；`2197-2208` 可以仅凭 `beast_trees` 计算 posterior corridor。
- `Virus_Platform_Core/beast_io.py:554-565`：burn-in 后若树数超过上限，直接取前 `max_trees` 棵；默认上限 2000。

**风险：** 用户只给 `.trees` 或链尚未收敛时，后验树自相关/有效样本量无法从树数本身判断；按顺序取前 2000 棵也没有抽稀/ESS 证据。结果仍会形成迁移走廊频率，容易被误读为稳定后验支持。

**建议：** 只有 `.trees` 时明确标记“未做收敛诊断、仅探索性”；正式频率报告要求配套 log 和 ESS 诊断，记录树采样间隔、burn-in、有效树数；超上限用可复现等距抽样/按 ESS 控制而不是无说明地截取前段。

## 3. P2 — 需要补强的验证、统计表达与文档

1. **DRT 默认置换数只 20。** `rtt_drt.py:54-56,297-345` 默认 seed=42、20 次；p 值分辨率下限 `1/21=0.0476`。模块已诚实说明小置换数限制，因此不是静默 bug；但“显著/不显著”几乎只能由最极端的一次结果决定。发表级结论建议至少 999 次（或给出更高分辨率要求），并固定记录有效置换数与 seed。
2. **没有独立诊断真值验证。** `docs/paper/manuscript_body.md:199-207,225` 已明确性能数据属于工程 benchmark，尚无 RT-PCR 等独立 ground truth、敏感度/特异度/LoD 验证。建议保持这一定位；将来测试需预先锁定阈值、盲法分析、阳/阴性对照、多批次与独立样本，并报告置信区间。不要从当前小型挑战集外推临床/检疫性能。
3. **Benchmark 指标目前以点估计为主。** `docs/paper/benchmark/compute_and_plot.py:138-144` 计算 PRF 点估计；现有数据集数量和真值项较少，未见分层置信区间/重复抽样不确定性。若用于比较引擎/阈值，应按样本/病毒组为重采样单位给区间，并处理多阈值、多引擎比较；不能把同一批 reads 的多个参数组合当独立生物学重复。
4. **Methods 参数与代码有可定位的口径漂移。**

| 论文/说明承诺 | 当前代码 | 审稿影响 |
|---|---|---|
| `manuscript_body.md:89` 写最小 contig 长度默认 500 bp | `main.py:946`、`pipeline.py:862`、`config.py:1008`、`platform.json:49` 当前均为 200 bp | 200–499 bp 序列的保留率及后续候选集不同；需要确认正文是否应更新或代码参数是否恢复。 |
| `manuscript_body.md:91` 写 mmseqs2 对 CDD 域模型 | `verify.py:265-285` 实际搜索 BFVD 库；`pipeline.py:536` 也记录 BFVD | 数据库来源、领域覆盖与可重复性说明不一致，Methods 必须改成实际库或代码改回承诺的方法。 |
| `manuscript_body.md:83` 把缺 ANI 描述为 undetermined | 当前 Salmon null ANI 会进入 `confirmed`，日志另提示 ANI 未筛 | 应统一状态语义，避免表格标签与稿件说法不一致。 |
| README/Pipeline-Overview 写“14 步”；论文方法写 15 阶段 | `pipeline.py:30-32` 的 `STAGE_ORDER` 含 15 项（包括子采样） | 可能只是是否计入 prep 的口径差异，但应明确口径，避免阶段数和流程图对不上。 |

## 4. 修复优先级建议

**第一档：先阻断错误结果进入报告**
1. 修复 P0-1：零命中/异常结果必须以“本轮新产物”覆盖旧文件；汇总端强制校验产物指纹。
2. 修复 P1-6：把筛选通过与 ANI 核验分开命名；Salmon 无 ANI 不得暗含“已确认”。
3. 修复 P1-7/8：benchmark 失败退出作废、真值不确定项隔离、续跑按数据/参数哈希判定。

**第二档：保证运行与分析口径可复现**
4. 用统一 fingerprint 覆盖 subsample、fastp、host removal、assembly、kvsuite 的全部输入/数据库/参数。
5. 重新设计未知 contig 的验证入口，或严格缩小“新病毒候选”表述。
6. 修正 PE coverage fallback；明确子采样前缀仅用于快速测试。
7. 公共数据检索记录每个 source 的成功/失败/零命中状态。

**第三档：统计与论文提升**
8. 将读段深度版 Tajima’s D 和无 ESS 的 BEAST 树频率限制为探索性输出，补充独立验证或收敛门槛。
9. 对 benchmark 增加重复/置信区间，Methods 按当前实际默认值、数据库和输出标签逐项对账。

## 5. 值得保留的设计

- 主流程有单一阶段注册表和清晰的样品级调用链；不是只看函数存在与否。
- Fastp 对零 reads/不可解析报告会失败，kvsuite 有 `failed_samples` 并提醒“失败不能当阴性”；这类语义闸门应继续扩展到公共数据检索和所有阶段产物。
- kvsuite 批次指纹包含 engine、版本、参考、ref_info 与阈值；问题主要是外层快速缓存和其他阶段没有统一使用它，值得收敛而非推倒重来。
- DRT 使用固定 seed、经验 p 的 +1 修正、有效置换数/低功效提示，并暴露根搜索的反保守风险；这些是严肃的统计防护。
- 稿件主动承认“尚无独立生物验证”，将实现级证据与诊断级证据区分开来，是正确的审稿姿态。
- 自动化测试有真实的旧行为对照臂和敏感性守卫；下一步应把相同纪律推广到失败路径、cache provenance、候选筛选偏倚与基准真值。

---

**本轮验证范围：** 静态检查调用链、代码与 Methods 对账；对 P0-1 做了一次临时目录隔离的合成回归。未运行全量测试、外部工具、完整样品流水线或 benchmark；未修改业务代码和样例数据。