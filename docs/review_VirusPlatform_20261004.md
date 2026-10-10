# VirusPlatform 全面审查报告（2026-10-04）

> 审查方式：环境自检 + 测试套件（--quick）+ 实际运行（report 重生成）+ 三路子代理代码审查 + 11 张产物图视觉验收（visual-judge ×3）+ reviewer-two / reproducibility-reviewer 双章程审稿。
> 所有问题均有实测证据，标注文件:行号。
>
> **➡ 2026-10-06 修复记录见文末附录（全部 P0/P1/P2 已修 + 守卫测试 + 全量回归）。**

## 结论

**建议决定：大修（Major Revision）**——作为软件工具类稿件（Bioinformatics Applications 级别）。
平台骨架优秀（防御性编程、种子可复现、示例链路设计、文档-实现对齐度高），但存在 4 个实测复现的核心功能缺陷，其中 **AA 同一性数值错误直接动摇 SDT 核心卖点的可信度**，必须修复后重审。

---

## 一、致命/主要问题（必须修）

### M1. AA 同一性出现生物学不可能值【P0·数据正确性】
`examples/results/identity/identity_table.csv`：
- ORSV vs ORSV_：NT 88.52% / **AA 33.53%**
- TBRFV_1 vs TMMV_va：NT 75.65% / **AA 27.95%**

核酸 76–89% 一致的 RNA 病毒，蛋白不可能只有 28–34% 一致。疑蛋白比对未做 MSA（逐位置硬比对）或翻译/ORF 提取错位。SDT 是平台 ⑦ 步核心输出，此错直接影响物种划界结论。**修法：AA 矩阵改用与 NT 同源的 MSA→翻译流程，加"AA≤NT-30 时告警"的自检；补一个已知病毒对的回归测试。**

### M2. 基因组圈图静默全灭【P1·功能可用性】
实测 `python main.py report --sample smoke_e2e`：3/3 contig 圈图失败，WARN 回退线性图。
根因：`Virus_Platform_Core/viz.py:496` 调 `sector.text(..., deg=90, ...)`，而 pycirclize 1.10.1 的 `Sector.text()` 已删除 `deg` 形参（现签名 `text(text, x=None, r=105, *, adjust_rotation, orientation, ...)`），`deg` 经 `**kwargs` 传入 matplotlib `Text` 抛 `AttributeError`。
`requirements.txt:23` 仅钉 `pycirclize>=1.0`，lock 锁 1.10.1 与代码不兼容，且无圈图回归测试。**修法：删 `deg=90` 改 `x=...` 居中；requirements 收紧为实测兼容区间；tests 增加圈图冒烟断言（产物存在且非回退标记）。**

### M3. 共线性 compare.svg 是非法 XML【P1·交付物损坏】
`examples/results/synteny/compare.svg` 第 150–151 行：`<title` 开标签后直接 `</g>`，title 元素未闭合，XML parse error。浏览器容错渲染出残图。矢量图是平台宣称交付物。**修法：排查生成端（matplotlib SVG backend 交互或文本注入），生成后加 `xml.etree.ElementTree.parse` 自检。**

### M4. NCBI 提交链 U→N 替换失效【P1·功能回归】
`tests/_it_submit.py:132` FAIL：期望替换 20 个 U，实际 0。提交准备是 ⑩ 之后的研究闭环终点，RNA 病毒序列带 U 提交 NCBI 会被拒。**修法：定位 submit-fill/replace 调用链，修复后该测试自然转绿。**

## 二、绘图质量（11 张产物验收）

剔除审查方截图裁切因素后，平台自身问题如下：

| 图组 | 问题 |
|---|---|
| 热图 ×3（identity/nt/sdt，平台直接产出 PNG） | ① 黑色网格线错位半格、穿过色块中部并外伸（渲染缺陷）；② 行列标签程序截断至 ~40 字符，`..._vir` 与 `..._vir_1` 同前缀样本不可区分；③ 色标 20–100/50–100 不从 0 起；④ YlOrRd 低端 50/55/60 档色盲不可辨 |
| sdt_distribution | 仅 13/15 组合入图（2 个低值被 x 轴起点裁掉），y 轴 Proportion 合计≈0.87 ≠ 1，误导 |
| tree.png（真实运行） | 支持率数值压在枝线上难辨；叶标签越面板边框被穿字/词中截断；横轴无 substitutions/site 单位 |
| kv_depth | 图例框+信息框叠放遮断深度峰；未平滑灰色曲线无图例；标题 "Sample: EXAMPLE" 占位名 |
| orf1（dna_features_viewer） | domain 命中条仅存在于 SVG `<title>` 悬浮提示，静态导出/打印不可读；ORF 名省略号截断 |
| e2e 圈图（历史产物） | 标题直接用 `NODE_1_length_6240_cov_29.664667`，成图标题不规范 |

正面：中文渲染全部正常无方框；数据层（GC 曲线、基因方向、对角线 100%、深度均值/阈值线）本身合理。
**修法顺序：热图错位网格线 → 标签不截断（或缩写+图注全名）→ 树图避让与单位 → 深度图图例避让 → 命中条可见图例。**

## 三、工程健壮性（代码审查摘要）

1. 【P1】`utils.py:724,740` 超时看门狗只 `proc.kill()` 直接子进程，与取消路径的 `taskkill /F /T`（web/tasks.py:570）语义不一致——spades.bat→spades.py、salmon 孙进程成孤儿。
2. 【P2】`assembly.py:677` `memory_gb=64` 硬编码且不参照 `utils.py:605` 的物理内存探测；非 shim 分支成功后不清理 K* 中间目录（可达数十 GB）。
3. 【P2】并发闸门只限任务数（heavy=2）不限线程（每任务默认 19）→ 38 线程超订 20 核。
4. 【P2】`app.py:306` 端口候选 `[8765,8900,8989,9600,8888,5050,5000]+随机` vs `scripts/drag_proxy.py:41` 只扫 8760–8789——8765 被占则代理必失联。
5. 【P2】无全局 `MAX_CONTENT_LENGTH`；`submit.py:95`、`meta.py:990` 上传无大小预检（`io_api.py:68` 已有 64GB 预检可复用）。
6. 【P2】`app.py:226-231` 无 Origin/Referer 即放行，本机任意进程可驱动删样品等破坏性端点。
7. 【P3】`kv_engines.py` 15 处裸 subprocess 无超时；`consensus.py:458-470` CIGAR M 段 SEQ 越界无保护；`main.py` 无顶层异常包装、logger 无 try/finally；`webapp/static/i18n.js.bak_20261001_prescreen` 等备份残留触发 `_check_spec_datas` FAIL。

正面：`shell=False` 全程、`check_path` 防穿越集中化、上传名清洗、AI 执行单三重白名单+审计、密钥分离（api_key 落 run/ai_local.json）、并发互斥锁、测试运行器超时分级设计——均为同类项目少见的工程规范。

## 四、示例数据治理

- ✅ TMV/PVY/PSTVd/CMV 三段、PZ 系列登录号与 NCBI RefSeq 逐碱基 100% 一致（抽样实测）。
- ❌ `example_acmv.fasta`：登录号 NC_001467.2 **不存在**（现行 ACMV DNA-A 为 .1/2779bp，该序列 2801bp 且仅 62% 一致）。
- ❌ `example_conserved_set.fasta`：5 条全挂同一真实登录号 PZ530881.1，isolate1–4 系虚构，其中 2 条完全相同——"保守区"示例退化。
- ❌ `example_pdprep.meta.csv`（name 列）与 `example_pdprep.fasta`（头部名）**零精确匹配**；实际配的是 phylogeo 对；`_check_examples_inputs.py` 恰好只校验 phylogeo 对，漏掉此错配；伴生裸 pdprep.fasta 与 phylogeo.fasta 逐条 100% 相同属冗余。
- ❌ `example_treedater.nwk` 为 476 叶**人流感 H3N2** 树，叶名与任何示例 FASTA 不匹配。
- ❌ `examples/results/` 34 个文件含本机绝对路径（混 `tool_runs\` 与 `run\tool_runs\` 两个根），不宜对外分发。
- ❌ `manifest.json` schema 混用（6 模块 files 是字符串数组，26 个是 {path,size}），45 个声明文件不存在、73 个实际文件未声明。

## 五、可复现性审查

**结论：源码级可复现性好；产物快照与二进制分发层有缺口。**

正面：
- 降采样随机种子固定并写进结果（`phylodyn_kit.py:489` `random.Random(seed)`，结果含 `'seed': seed`，注释明示对齐 GeoSubsampler 口径）。
- `tests/make_example_results.py` 用与用户点击完全相同的 `/api/tool/run` 链路固化示例——示例即真实链路，设计优秀。
- `requirements.lock` 106 包；`__version__="1.0.0"`；外部工具版本实测记录在 README 附录。
- 测试运行器：独立超时+互斥锁+陈旧锁回收，注释记录了历史事故根因。

缺口：
1. `pyrodigal-rv`、`viennarna` 出现在 requirements.txt 但未在 lock 中锁定（Flask/pyCirclize 为大小写差异非缺失）。
2. `3rd/` 2GB 预编译二进制**无 checksum 清单**——克隆者无法校验完整性，与"二进制审计"诉求矛盾（README 自述源码构建一节即为此留档，但无对照基准）。
3. `examples/results/` 快照与当前代码输出脱节（45 缺 73 多），`_check_examples_manifest.py` 明文"允许落后"——快照失去"预期产物"的对账价值。
4. 数据库版本（taxonomy/ICTV 下载日期戳）未在 manifest 层面记录。
5. `ncbi_submit/_server_collect.py` 含 `/home/zhangwenda/...` 服务器绝对路径（文档性质，但建议参数化默认值）。

## 六、文档与测试覆盖

- README 6 个 wiki 链接 slug 失效（Assembly-Verify-Consensus / ORF-Annotation / Phylo-SDT / Primer-Design / BEAST-Interface / Logan-Trace；实际文件为 Assembly.md/ORF.md/Phylo.md/Primer.md/BEAST.md/Logan.md）。
- Primer 模块 wiki 强调的热力学评分/推荐分级/宿主特异性**无专项测试**（仅 `_it_annotate.py` 断言 n_primers>0）。
- 测试 51/56 通过（--quick、跳过在线），失败 5 项根因：M4（U→N）、`_it_dl_to_pipeline`（外部任务状态读取，待查）、`_check_hostdb`（宿主库旧布局未迁移，环境项）、`_check_spec_datas`（.bak 残留）、`_check_pages_console`（已清理旧任务的数据表 404，环境残留）。

## 七、修复优先级路线

1. **立即（数据正确性）**：M1 AA 同一性 → 加回归测试。
2. **立即（功能）**：M2 圈图 deg（一行修复）、M4 U→N、M3 SVG 合法性自检。
3. **本周（绘图）**：热图错位网格线与标签截断 → 树图避让+单位 → 分布图 x 轴范围。
4. **短期（示例治理）**：ACMV 换真实序列、conserved_set 去虚构、pdprep 对键名、treedater 树换植物病毒、results 快照去本机路径并重生成。
5. **迭代（工程）**：进程树统一终止、闸门-线程联动、端口发现统一、MAX_CONTENT_LENGTH、lock 补全 + 3rd checksum 清单。

---
*审查产物：渲染对比图 11 张位于 `run/_audit_figs/`；测试明细 `tests/_run_all_report.txt`。*

---

# 附录：修复记录（2026-10-06）

> 审查发现的全部 P0/P1/P2 与大部分 P3 已修复，每项修复均有验证；新增 4 个守卫测试钉住回归。

## 核心缺陷（M1-M4）

| 项 | 修复 | 验证 |
|---|---|---|
| M1 AA 同一性伪值 | `sdt_exact.py` 新增**分簇锚定翻译**（`_anchored_aa_seqs`）：NT≥70% 并簇 → 簇内最长 ORF 定 (链,框) → 整框宽容翻译（终止→X 不切断）→ 多轮重簇锚定 → 末位回退独立 ORF；AA 生成移至 NT 矩阵后 | ORSV 对 AA 33.53→**74.86**，15/15 对生物学合理；新守卫 `_check_identity_aa.py` 通过 |
| M2 圈图静默回退 | `viz.py` 删 `deg=90`（pycirclize≥1.9 已移除该形参）；标题经 `_pretty_contig_name` 美化并移至几何圆心（r=0） | report 重生成 3/3 真圈图；守卫 `_check_genome_circos.py` 通过；视觉终检 pass |
| M3 synteny SVG 非法 XML | 18 处 `</path` 缺 `>` 已补齐；`utils.verify_svg()` 幂等自检挂到 dfv_plot 两处 savefig；守卫 `_check_svg_wellformed.py` 扫描 examples 全部 SVG | ET.parse 全部合法 |
| M4 U→N "失效" | 实现本为正确口径（U→T，NCBI DNA 等价形式，单独计 u_to_t）；**过时断言**修正：`_it_submit.py` 改断言 u_to_t==20、fsa 无残留 U、replaced==0；顺带删恒真断言（`or True`） | SUBMIT INTEGRATION PASSED |

## 绘图（D1-D6）

| 项 | 修复 | 验证 |
|---|---|---|
| D1 热图错位网格线 | `plot_heatmap`/`plot_heatmap_composite` 网格线移到格边界 ±0.5（三角模式按可见阶梯精确收线） | 像素级核验：线全部落格边界，4/4 pass |
| D2 标签截断/色档 | `_disp_labels`：n≤12 全名、大 n 截断保唯一；`_cell_values`：n≤14 格内标数值（白/深字按底色）；数据层名称截断 40→56 | `_vir` vs `_vir_1` 可区分；36 格数值与矩阵逐一吻合 |
| D3 分布图丢组合 | bins 从 `floor(nanmin)` 起覆盖全域（原 `max(30,min)` 裁掉低值） | 15/15 组合入图、Proportion 合计 1.0 |
| D4 树图 | 支持率加白色 halo 描边；去上/左/右边框（穿字根因）；横轴补 Substitutions per site；标签 48→64 字符 + 画布自适应加宽 | 两轮目检：halo/单位/标签 `_1`/`_2` 完整可辨 |
| D5 深度图 | 示例快照重生成（新版代码：单行主标题+灰副标题+despine）；reads 口径修正（MappedReads 优先，Uniq_Reads×2 折算）；coverage 两位小数；图例挪左上 framealpha 0.45 | 终检 pass（2,400 reads、99.98%、峰尖无遮挡） |
| D6 ORF 图例/圈图标题 | `genome_diag.py` 图底加结构域命中图例（含深浅分道说明）；圈图标题去双后缀+圆心居中 | 图例存在且不遮图形；圈图终检 pass |

## 示例数据（E1）

- `example_acmv.fasta` → 真实 **NC_001467.1**（efetch 下载，2779bp，登录号核验存在）。
- `example_conserved_set.fasta` → 1 真实 ORSV + 4 合成派生（SIMisolate1-4，突变率 0.5-1.5%，**seed=20261004 可复现**，defline 诚实标注 synthetic）。
- `example_pdprep.meta.csv` name 列与 phylogeo.fasta 头逐字一致（治理链路 15/15 通过）；UI 卡 `t-pdprep` 改配 phylogeo.fasta + pdprep.meta.csv（与 make_example_results 权威链路对齐）。
- `example_treedater.nwk` → 用平台 MAFFT+FastTree 对 example_phylogeo.fasta 现建的 15 叶植物病毒树（叶名与比对一致）。
- `example_mirna.fasta` miRN73 注释能量 -26.2 → **-22.8**（真实 ACMV 上 consensus 引擎实测）。

## 工程健壮性（E3-E5）

- **进程树终止**：`utils._kill_proc_tree`（taskkill /F /T）替换 run_cmd 三处裸 kill（超时看门狗/TimeoutExpired/KeyboardInterrupt），与取消路径口径统一。
- **资源钳制**：POST threads ≤ 物理核、memory ≤ 物理内存 60%（web/samples.py）；CLI 直调路径同样钳制（assembly.py，`-m 64` 不再硬通）。
- **组装中间目录**：非 shim 分支成功后清理 K* 目录（原仅 shim 丢弃，可达数十 GB）。
- **端口发现**：app.py 绑定后落盘 `run/.port`；drag_proxy 优先读端口文件 + 扫描全候选集（原只扫 8760-8789，错配必失联）。
- **上传上限**：`MAX_CONTENT_LENGTH=16GB`（原无限，防误传打满内存）。
- **consensus.py**：CIGAR M 段 SEQ 越界跳过（原只护 QUAL）。
- **kv_engines.py**：6 处阻塞 subprocess.run 补 `timeout=7200`（管道 Popen 留外层 6h+进程树兜底）。
- **main.py**：顶层异常包装（0 成功/1 selfcheck 未就绪/2 运行失败/130 中断；VP_CLI_TRACEBACK=1 可看回溯）；`TaskLogger` 加上下文管理器+析构兜底（close 幂等）。
- **.bak 残留**：4 个手工备份归档 `archive/_bak_20261004/`，`_check_spec_datas` 转绿。

## 可复现性与文档（E6 + E2）

- `examples/results/` 全量重生成（25 模块 + pdprep 修复重跑）：本机绝对路径 **0 残留**；manifest schema 归一为 {path,size}（0 混用）、对账 **185 声明 0 缺失**（自动剔除缺失声明并注明）。
- README.md 6 处失效 wiki slug 修正（README_EN 本就正确）。
- `requirements.lock` 补 ViennaRNA==2.7.2（pyrodigal_rv 为下划线/连字符同包，原判缺失系误报）。
- 新增 `3rd/CHECKSUMS.sha256`（1613 个 exe/dll/bat/pyd/jar 的 SHA-256，180KB）——2GB 预编译二进制可完整性校验。

## 遗留（非代码缺陷，供运维参考）—— 2026-10-06 第二轮已全部处置

1. **宿主库旧布局 → 已重建闭环**：下载番茄参考基因组（NCBI datasets CLI，RefSeq GCF_036512215.1，核 + 叶绿体 + 线粒体共 54 条/849MB），`build-host-db` 建库至新布局 `host-db/4081_Solanum_lycopersicum_host_db/`（1.15GB hash 表，自动写 host_db.json 并设为 active）；selfcheck 宿主库 ✔，`_check_hostdb` 全绿（其中"kreport 不得出现 Solanum"断言原为枸杞污染事故的写死口径，已普适化为"物种级节点不得偏离声明 taxid"）。
2. **`_it_dl_to_pipeline` "时序敏感" → 实为测试桩契约 bug，已修**：FakeQueue.add 返回字符串列表，而真实 SampleQueue.add 返回条目 dict；download.to_pipeline 的 bind_queue_entry(it) 拿字符串调 .get('id') 500。测试桩对齐真实契约后 ALL PASS。
3. **旧检索任务 404 → 已修**：/meta 页对未提取元数据的项目（n_core14=0）仍自动请求 core14 表，浏览器把 404 记为 console.error。meta.html 现按 n_core14 跳过请求并显示"先提取元数据"提示（i18n 双语）；`_check_pages_console` 22/22 页面通过。
4. **闸门-线程完整联动 → 已实现**：TaskManager 维护 heavy 名额实时占用计数（acquire/release 配对），pipeline 启动时按 物理核 ÷ 占用数 钳制本任务线程（CLI 直调不变）；单元验证 20 核 ÷ 2 任务 = 每任务 10 线程（原 2×19 超订）。
5. **Origin/Referer 口子 → 威胁模型复核后维持并文档化**：浏览器跨站请求（含 CSRF）fetch/form POST 必带 Origin，已被完整源校验拦截；无 Origin 的仅剩本机进程（权限等价登录用户，回环绑定即边界）；破坏性路由 14 个全为 POST、无 GET 型。收紧对 CSRF 无增益且会打断全部 test_client/CLI，故将三项依据写进 app.py 注释固化决策。
6. **新发现并修复：host-genome 下载无截断校验**——778MB 基因组 zip 中断下载后静默落盘坏文件、下游只剩 WARN 且 exit 0；`_Scoped.download` 现对照 Content-Length 校验，不完整即丢弃残件并报错。
7. **批量回归波动（`_it_async_ui` 偶发超时）→ 已加余量**：浏览器元素等待写死 40s，批量满载时偶发不够（单跑 67s 全过、批量 105s 翻车）；等待时长改走 `VP_TEST_WAIT` 环境变量（默认 40 不变），`_run_all.py` 批量调度时注入 90s。

**第二轮处置后回归：56/56 全绿**（第一轮修复后 53/56 → 遗留处置后全量通过）。

---

# 附录 2：P3 清零轮（2026-10-06 第三轮）

> 审查清单中剩余的 12 个 P3 小项全部处置；其中发现并修复 1 个真实功能退化（Primer 示例 0 引物）。3 项大重构类建议明确不做并附理由。

## 已修（12 项）

| 项 | 修复 |
|---|---|
| kv_engines 管道 Popen 未注册取消表 | minibwa→samtools 管道两进程注册 `task_register_proc`（硬停止可杀） |
| kv_engines 命令日志不可复制执行 | `' '.join` → `shlex.join`（含空格路径加引号） |
| host_removal 配对后缀探测静默失败 | except 补 WARN 日志（不再无痕退化为无后缀匹配） |
| consensus 进度分母量纲错 | 文本 SAM 统一字符口径；gz 场景改按行数报进度（不给假百分比） |
| app.py Host 头 `localhost.:8765` 被 403 | 剥尾部圆点后比对 |
| tasks.py 排队取消竞态 | acquire 成功瞬间复查 cancel，被取消则归还名额走取消分支 |
| tasks.py 进度持久化写盘风暴 | 进度帧 ≥1.2s 落盘一次（终态/状态变化仍立即写；节流键不落盘） |
| submit 上传后缀未校验 | 白名单 csv/txt/tsv/xlsx/xls，其余统一 .csv |
| virome 资产路由风格不一 | 补 check_path 与全站下载路由统一口径 |
| drag_proxy 文档 `--host` 参数未实现 | 注册参数并贯通 send/find_platform |
| add_segment_std/normalize_ref_info 二跑固化改写产物 | filecmp 对比：备份只保留首跑原稿；硬编码路径可 `VP_DBPIPE_DATA_DIR` 覆盖 |
| Primer 无专项测试 | 新增 `_check_primer.py`（走 /api/tool/run 同链路：n≥3、score 0-100、分级枚举、长度/Tm），并注册进 `_run_all.py` |

## 顺带发现并修复的真实退化

**Primer 示例 0 引物**：E1 轮的 conserved_set 用均匀随机点突变（0.5-1.5%），把最长连续保守段打散到 233bp（< min_len=400），保守区提取为 0 → primer3 无从设计，示例 summary 记 `n_primers: 0`。改为**块状突变**（模拟真实分离物成簇变异：每条 3-6 个 250bp 窗口内 5-8% 突变、seed=20261006），最长连续保守段恢复 489bp，示例重生成得 **9 对引物**（score 73.3-95.6、RECOMMENDED/USABLE 分级正常）。

## 明确不做（3 项，理由）

1. **main.py 40 命令拆分子模块**（1268 行单文件）：纯组织性重构，无行为收益；现有 argparse 全在同一 `main()` 内、参数命名已一致，拆分反而增加 CLI 回归面。建议随下一次功能性改版一并做。
2. **dev_tools/package.py 模块顶层执行改 main()**：仅维护者工具，import 期执行的唯一代价是不可被单测导入；为其重构打包脚本的收益不抵验证成本。
3. **CSP 去 unsafe-inline/eval**：前端模板体系（plotly 注入、内嵌配置块）全面依赖内联脚本，且已有 Origin 完整源校验 + 仅回环监听双层防护；改造需重排全部静态资源引用，属独立工程。

**第三轮后全量回归：60/60 全绿**（新增 4 守卫并入 _run_all）。
