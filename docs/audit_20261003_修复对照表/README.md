# 修复对照表总索引（2026-10-03 全库体检 → 全量修复）

上一轮只读审计（11 个子系统、150+ 条发现）的全部实质问题已完成核验与修复。
本目录 11 份分表为**逐条对照**（问题 → 修复方式 → 改动位置 → 验证结果），供人工审核。
所有改动**未提交 git**，可整体 `git diff` 复查或逐文件回退。

## 修复统计

| 分表 | 范围 | 清单项 | 已修 | 核实后不改/证伪 | 结果将变化 |
|---|---|---|---|---|---|
| [01-主管道](01-主管道.md) | pipeline/consensus/verify/assembly 等 | 25 | 25 | （2 处审计前提修正） | ✅ 共识 min_depth 1→4（第二批裁定"深度>3"） |
| [02-kvsuite](02-kvsuite.md) | 已知病毒识别/定量/变异 | 20 | 20 | 0 | ✅ TPM/DNA0.5/双链/-d5 四道守门恢复 |
| [03-系统发育I](03-系统发育I.md) | 树/分子钟/SDT/BEAST/RDP5 | 21 | 19 | 2（beast with 块证伪等） | — |
| [04-系统发育II](04-系统发育II.md) | phylogeo/phylodyn | 23 | 22 | 1（文案已由 03 落地） | — |
| [05-注释与可视化](05-注释与可视化.md) | ORF/注释/viz/参考库 | 20 | 19 | 1（orf.py min_nt 证伪） | — |
| [06-miRNA与引物](06-miRNA与引物.md) | mirna_target/primer/LOGAN | 24 | 24 | 0 | ✅ 比对方向修复后命中数 2→4（24 对） |
| [07-提交与公共数据](07-提交与公共数据.md) | NCBI 提交/检索/元数据 | 37 | 36 | 1（33 号已达标） | — |
| [08-Web后端](08-Web后端.md) | web/* + app.py | 21 | 21 | 0 | — |
| [09-Web前端](09-Web前端.md) | webapp/* | 9 | 8 | 1（loadNcbiCollections 证伪） | — |
| [10-入口与文档](10-入口与文档.md) | main.py/.bat/scripts/dev_tools/tests/README/wiki | 22 | 22 | 0 | — |
| [11-论文基准](11-论文基准.md) | docs/paper/benchmark | 23 | 23 | 0（另 4 处附加修正） | ✅ 基准数字改为与归档数据一致 |
| [12-dsRNA恢复](12-dsRNA恢复.md) | dsRNA 设计（迁出后恢复 + 独立审计 + 修复优化） | 16 | 13 | 3（保留原语义/环境前置） | ✅ 见该文第四节：引物限定区间、0 错配档语义变化 |

另加 **4 项跨代理交接缝隙修复**（主会话直接完成，见文末）与 **2 项清单外附带修复**（10 号分表）。
**合计处理 250+ 项；P0 全部 7 项、P1 全部 22 项闭环。**

> 2026-10-03 追加：dsRNA 设计子系统（2026-09-29 迁出，同日由主会话从 HEAD 恢复核心
> 五件套）已由 **12 号分表** 完成接线恢复 + 1846 行核心代码独立审计 + 修复优化 +
> 全量验证。**10 号分表第 17 项 a)/b)**（删 README 的 dsRNA 短语、删工具表 `dsrnamax`）
> 已被 12 号分表的恢复动作反转，以 12 号分表为准（10 号分表作为上一轮历史记录未改）。

## P0 修复对照（速览）

| # | 问题 | 修复 | 验证 |
|---|---|---|---|
| P0-1 | salmon 索引折叠 931 条重复参考（192 对异名物种）静默丢定量 | `--keepDuplicates` + 旧索引 WARN | salmon 2.7.0 实测支持；plant_ref 旧索引命中 WARN |
| P0-2 | identify 批次缓存无指纹，旧 FASTQ 冒充新结果 | 批次文件带 input_fingerprint（r1/r2+库+引擎+阈值） | 4 轮实验：未变跳过/改 mtime 重算/旧格式重算 |
| P0-3 | 共识 min_qual 对无质量 BAM 空转 + min_depth=1 单 read 定 SNV | QUAL='*' 检测+`min_qual_effective` 字段；MIN_DEPTH 默认提高（后续第二批裁定为 4，"深度>3"） | FASTA/FASTQ 双输入端到端 |
| P0-4 | miRNA 比对方向反（镜像位点），坐标/seed/identity 全错 | DP 改反转帧（miRanda 书写约定）+ 三处消费端适配 | 正控 identity=100/duplexfold -36.50 一致；反控不判中；新增 17 例回归测试 |
| P0-5 | featuretable.tbl 丢链方向（反向 CDS 全写成正向） | st=='-' 时大坐标在前 + 部分标记 | 3 CDS 用例：`700 500` 正确 |
| P0-6 | 提交 FASTA 把 U 换成 N（应 U→T） | 先 upper+U→T 再白名单清洗；>1% 替换告警区块 | ACGU×8 → u_to_t=8、0 个 N |
| P0-7 | 基准：D5/D6 叙事与数据矛盾、fig_b6 丢 8/10 柱、指标写死、库已消失 | 以 results.json 重写、过滤器修、真值外置、库路径参数化+版本披露 | figs b1-b7/b9 重生成并目检 |

## 回归验证记录（修复后全绿）

- `py_compile` 285 个自有 py（含恢复的 dsRNA 四件套）：**0 失败**
- `pytest tests/ -q`：**75 passed**（修复前 58；新增 miRNA 比对回归 17 例）
- `_check_js_syntax` ✔ / `_check_i18n` ✔（新增弹窗 CJK 扫描，裸中文字面量=0）
- `_route_inventory` ✔（228 条路由基线已含 dsRNA 恢复新增，一致）`/ _audit_routes` 230 请求 0 异常
- `main.py selfcheck`：依赖全部就绪（③ 宿主库未建为本机环境态，与修复无关）；dsRNAmax 探测条目已恢复
- `_it_dsrna` DSRNA_IT_OK（2.3s，合成面板自足化）/ `_check_phylodyn_kit` / `_check_phylodyn_e2e`（真管线+真浏览器 32 项）/ `_check_pdprep_chain` / `_check_phylogeo_fitch/geo` / `_it_platform` PLATFORM CHECKS PASSED
- `test_web_guard` 5 passed（CSRF 端口级校验）；dsRNA 端到端冒烟 SMOKE_OK（产物与历史示例逐字节一致）

## 第二批调整（用户裁定，2026-10-03 下午）

| # | 裁定 | 落地 | 验证 |
|---|---|---|---|
| 调① | **类病毒支 >200bp；组装输出下限默认 200** | `verify.py` LEN_VIROID_MIN 250→200（[200,1000) 走类病毒支，<200 不处理）；`config.py` DEFAULT_FIELDS / `platform.json` / `main.py --min-contig-len` 三处默认统一 **500→200** | 分流实测：1200bp→病毒支、250bp→类病毒支、150bp→跳过 |
| 调② | **共识/变异：深度 >3** | `consensus.py` MIN_DEPTH 5→**4**（含注释口径"深度>3"）；`web/tool_jobs.py` ③c 面板默认值同步 5→4 | `MIN_DEPTH==4` 断言；py_compile 通过 |
| 调③ | **dsmax（dsRNA 设计套件）恢复 + 核查优化** | 见 [12-dsRNA恢复](12-dsRNA恢复.md)：核心五件套从 HEAD 恢复 + 12 处接线 + 1846 行独立审计（16 项，修 13）+ 端到端冒烟 | `_it_dsrna` OK、SMOKE_OK、路由基线一致 |
| 调④ | **脱靶面板不随软件分发**（1.4GB 瘦身） | `package.py` `_DB_TREES` 移除 `dsrna/panel` 条目（数据库包不再含）；i18n zh/en 面板标签 + wiki/Toolbox.md 注明"不随软件分发、自备放置"；功能保留——不放置时三态逻辑自动跳过脱靶扫描，效价链路不受影响 | py_compile + `node --check i18n.js` 通过 |

> 说明：kvsuite ②b 共识段的 `--vc-depth`（原管线/基准口径=5）不在本批"深度>3"调整范围内，
> 保持 5；如需同步为 4，改 `known_virus_suite.py --vc-depth` 默认值一行即可。

## ⚠️ 结果将变化（有意为之，审核重点）

以下修复恢复了文档承诺的守门，**重新分析旧样品会得到更严的数字**：

1. **共识/变异**：深度 ≤3 的位点不再定共识碱基、不再产出 SNV（P0-3 + 第二批调②，MIN_DEPTH=4 即"深度>3"）。
2. **kvsuite**：TPM≥1 门槛真实生效；DNA 病毒泊松比改用 ≥0.5；变异需正负链**各自**达标；共识 `-d` 默认 1→5（02 号分表第 3/4/5/6 项）。
3. **组装输出下限 200**（第二批调①）：200–1000bp 的 contig（kunpeng 未判毒）进类病毒支体检；**已判病毒的 contig 无论长度都走病毒支**（01 号分表第 2 项），不再因长度丢失身份。
4. **miRNA 命中数上升**（真实位点不再被判镜像失败），坐标/种子口径与物理互补一致。
5. **论文基准头条数字**改为与归档数据一致（minibwa 19/27=70.4%、S9 全链路 28.6% 等，详见 11 号分表）。
6. **dsRNA 设计**（12 号分表 D1/D2）：① 引物/T7 模板改按**臂（sense）帧**最优窗切片（此前取到
   guide 帧镜像区间，离中窗时错位）——交付扩增子落回真正的最优窗；② 界面"0 错配"档（`max_mm=0`）
   不再被静默放宽成 ≤1 错配——该档下的脱靶计数/红线判定会变严。`report.json` 新增
   `best_window.sense_nt_*` 与 `frame` 字段（纯新增，旧字段语义与旧运行回退不变）。

## 未修复 / 遗留清单（如实披露，均已在分表注明原因）

1. `drt_permutation` 置换检验开关已实现（默认关），调用方未接线——不影响现有行为，启用属后续增强。
2. pandas 3.x 下 `standardize_metadata_df` 的 object-dtype 门控不触发（既有问题；修改会变更数据输出口径，需专项决定）。
3. `main.py:487` cmd_host_genome 的 `--ncbi-api` 为进程内传参（不进 OS 进程表，风险显著低于子进程场景，未改）。
4. kvsuite B 轨仍无 genes_cov 生产者（已加显式 WARN"B 轨不可达"；补生产者或下线该轨属产品决策）。
5. 旧 salmon 索引不自动重建（复用处 WARN 提示需重建后 `--keepDuplicates` 才生效）。
6. `package.py` 新增 build_info.json 为函数级实测，全量打包（GB 级）未演练。
7. 论文基准如实披露项（11 号分表末）：`_bench_kvs2/_bench_cons` 缺失待重跑、kv_index 未重建、perf.csv 23 行缺陷数字只披露未重算、8/9 阈值扫描未补做、fig_b8/b10 未重生成。
8. `3rd/` 的 122 项删除（mafft ms/ 等）保留工作区现状，待作者确认取舍。
9. miRNA seed 窗保持本模块 2–8 口径（未按审计建议改 2–12；改口径会连带变过滤结果，推导见 06 号分表）。
10. `.bat` 在 cmd.exe 下的 UTF-8 注释碎片怪癖（布局性，退出码与功能不受影响，已记录）。
11. 3 处需浏览器/真实结果环境的端到端展示（LOGAN per-contig 报告、UNKNOWN 兜底、plain 窗口 WARN）以数据契约+单测替代验证。
12. 审计误报存档：beast_handoff with 块（本就规范）、orf.py min_nt（orfipy 官方语义核实）、loadNcbiCollections（调用有效）——3 条"问题"经核实不成立，未改。

## 跨代理交接缝隙（主会话收尾，均已验证）

| 位置 | 问题 | 修复与验证 |
|---|---|---|
| `config.py:_ensure_not_system_dir` | `os.path.join('C:',seg)`='C:Windows' 盘符相对 → 系统目录守卫恒放行 | 改 `drv+os.sep+seg`；C:\Windows\System32\Program Files\ProgramData 全拦截、正常目录放行 |
| `kv_consensus.py` | sRNA 0 映射静默产全 N 共识（stderr 自报成功 vs BAM 0 mapped） | 新增 `_bam_mapped_count`+样品级守卫：0 映射 ERROR 并跳过该样品；真实 mini BAM 判别 1/0 正确 |
| `phylodyn_local.py:_slim_clock` | 白名单丢 argv 与 clock_filter/time_marginal/reroot | 白名单补 4 键；断言通过 |
| `phylodyn_trees.py:_parse_translate` | `[,;]` 强制结尾吃掉块尾 `;`，最后一条叶名映射永丢 | 分隔符改 `[,;]?`；'1 A,\n2 B;' → 两项全取 |
