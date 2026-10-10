# miRNA 靶标 / 引物 / LOGAN 提交溯源修复对照表（2026-10-03）

范围：上一轮只读审计清单第 1–24 项（miRNA 靶标工具组 + 引物工具组 + LOGAN 提交/溯源 + 平台级
联网管控的配套项）。**本组最敏感项是第 1 项（P0-4 比对方向反向）**，已按"正控 / 反控 /
ViennaRNA 交叉验证 / 官方 miranda.exe 实测输出"四重证据核验。

约定（与审计要求一致）：
- 每项先读代码核实再改；改动均为精准 Edit（未整文件重写，保持 UTF-8 与中文注释）。
- 改动文件全部 `python -m py_compile` 通过；`ruff check --select F` 仅剩**修复前即存在**的既有噪音（文末列名）。
- 验证实验一律在 `%TEMP%` 或内存中完成；只调用了 `3rd/` 里 miranda.exe / PITA perl 管线的**短调用**（单条序列），未跑全量批处理、未动 `3rd/`、未做 git commit。
- 行号以本次修复后的文件为准（工作区存在其它子系统并行改动，行号可能漂移，以函数名/锚点为准）。

## 对照表

| 编号 | 级别 | 位置(修复前) | 问题摘要 | 修复方式 | 改动位置(修复后) | 验证方式与结果 |
|---|---|---|---|---|---|---|
| 1 | **P0** | `mirna_target/alignment.py:129-159`（联动 `consensus.py:105-181,304-308`） | 核心比对器链向反：DP 拿 miRNA 5'→3' 原向做同索引互补匹配 = 找 `comp(miRNA)` 的**同向镜像**序列，把反平行物理位点判反。实测（修复前）：target 含 `comp(mir)` → identity **100.0**；含 `revcomp(mir)` → **47.6**；官方 miranda.exe 对同一 revcomp target 报 100.00% | ① DP 改用反转帧 `mir_rev = mirna[::-1]`（列 = miRNA 3'→5' vs target 5'→3'，与官方输出书写约定一致），traceback/打分同步；② `t_start=j` 左偏 1 → `t_start=j+1`；③ 位置图/种子/query 坐标按 `k ↔ L-k+1` 换算回原始 5'→3' 帧；④ 能量在**原始帧**做 cofold（`aln_a[::-1]`，两链各按 5'→3' 读入）；⑤ 比对串标签改 `miRNA 3' … target 5' …`；⑥ docstring 注明打分矩阵非 miRanda 原值（**数值未改**，避免连带改结果）；⑦ consensus 三处消费（`_validate_psrna` / `_validate_rnahybrid` / `compute_seed_quality`）改为自右向左遍历 + 解析新标签；⑧ tapir 调用点改用比对列（见第 2 项）；⑨ 前端比对渲染/Query2 判据同步适配（见"联动项 A2"） | alignment.py:11-19（种子常量+注释）、:68-83 / :107-124（帧内位置换算）、:127-233（DP/traceback/坐标/种子/能量/query）、:286-296（能量帧说明）；consensus.py:37-91（seed_quality）、:108-197（psrna 列序）、:199-211（rnahybrid 列序） | ① **正控**：target=`GGGGG+revcomp(mir)+AAAAA` → identity=**100.0**、seed_perfect=True、`target=6-26`（精确对准插入段）、`query=1-21`，tier_strict=True；② **反控**：target 含 `comp(mir)` → identity=**61.9**、seed=False、tier_strict=False（不再判中；审计旧值 100.0）；③ **ViennaRNA**：`RNA.duplexfold(mir, revcomp(mir))=-36.5` vs `duplexfold(mir, comp(mir))=-9.4`（审计 -37.3/-9.4 同结论），且比对能量 **-36.50 == duplexfold 位点值**（一致）；④ **官方 miranda.exe v3.3a 实测**：同一 revcomp target → `Score 195 Q:2 to 20 R:6 to 26 (100.00%)`，与我们的书写约定/坐标语义一致（若 target 只含 comp(mir)，miranda 亦只给 57.89%）；⑤ 回归测试 `tests/test_mirna_alignment.py` 17 passed；⑥ 示例端到端见"结果变化说明"第 1 条 |
| 2 | P1 | `mirna_target/tapir.py:77-130,254-279` | ① 列语义错：`validate_with_tapir` 传的是"无缺口 miRNA vs 原始目标片段"，反平行位点被当同向序列比；② target 短于 miRNA 时 `tline[i]` → IndexError；③ 种子区定位反向（在旧帧里自右向左数，2-12 位被当 10-20 位） | ① 调用点改传 **miRanda 比对列**（新帧上行=miRNA 3'→5'；兼容旧格式：标签为 5' 时先翻转上行）；② `score_alignment_columns` 长度不一致时按 `-` 补齐（不再 IndexError），返回 `score/…/seed_mismatch` 语义不变；③ 种子列定位改为"先翻回原始 5'→3' 帧（`mline[::-1]`）再自左向右数第 2..12 个非缺口字符"（等价于新帧自右向左数，但语义可自证） | tapir.py:36-38（`SEED_POS_START/END`）、:79-114（列语义+种子+补齐）、:189-270（`tapir_predict` 反转帧+种子窗修正）、:278-306（`validate_with_tapir` 用比对列+旧帧兼容） | ① 手工构造比对列（含 1 处真错配）：**第 3 位错配 → seed_mismatch=1、score=2.0**（判"种子错配"）；**第 16 位错配 → seed_mismatch=0、nonseed_mismatch=1、score=1.0**（判"种子保留"）；完美列 → score=0/无错配；② `target` 仅 8nt（<21）→ 返回 score=21.0，**无 IndexError**；③ 完美位点上 `validate_with_tapir` → tapir_score=0.0；④ pytest 用例（test_tapir_seed_band_positions / test_tapir_short_target_no_indexerror / test_tapir_validate_uses_alignment_columns）通过 |
| 3 | P1 | `mirna_target/psrnatarget.py:160-190` | `score_alignment` 缺 `'wc'` 分支：`_column_class` 返回 `'wc'|'gu'|'mm'`，else 把 **WC 完美配对也按错配 1.0** 计分（self 比对 penalty=25.0、完美命中被自己否掉） | 在 `gu` 分支前补 `elif _column_class(q,h)=='wc': curr=0.0`（并关闭缺口状态），与 MiRNATarget 口径一致 | psrnatarget.py:184-190 | `score_alignment(q,q)`（19 列）→ **penalty=0.0、seed_mm=0、total_mm=0、gaps=0、gus=0**（修复前 25.0/25）；pytest `test_psrna_wc_branch_perfect_pairing` 通过 |
| 4 | P1 | `consensus.py:279,285`（+ `web/tool_jobs.py:2531`） | 裸 `except: pass` 吞掉 TAPIR/psRobot 引擎失败；且 summary 里的 `algorithms` 记录的是**请求**引擎而非实际跑成，引擎崩溃会被伪装成"跑过但无命中" | ① 改 `except Exception as e` + `log.warning("TAPIR/psRobot 校验失败跳过: %r")`；② `ConsensusResult` 新增 `engines_actual`（实际投票引擎）/`engines_skipped`（失败引擎+原因）/`engine_notes`，`to_dict` 导出；③ tool_jobs 的 summary **保留** `algorithms`（字段兼容）并新增 `engines_actual`/`engines_skipped`/`engine_notes`（跨命中并集），有跳过时另打 WARN 日志 | consensus.py:245-250（dataclass 字段）、:276-315（try/except+记账）、:358-410（to_dict 导出）；web/tool_jobs.py:2623-2655（summary） | ① pytest `test_consensus_engine_failure_recorded`（打桩 `validate_with_tapir` 抛 ValueError）→ `engines_actual` 无 tapir、`engines_skipped=[{'engine':'tapir','reason':'ValueError: 模拟引擎失败'}]`，且走的是 `log.warning` 而非静默；② 正常路径 `engines_actual=['miranda','psrna']`、`engines_skipped=[]`；③ tool_jobs 汇总块**按源码原文抽取后用真实 rows 执行**（`exec` 去缩进片段）：`algorithms=['miranda','psrna','rnahybrid']`（请求）vs `engines_actual=['miranda','psrna']`（实际）、`engines_skipped=[]`、无 WARN——两字段语义已分离 |
| 5 | P1 | `rna22.py:59-70,148`；`psrobot.py:64-97`；`rnahybrid.py:98` | 口径不诚实三连：① RNA22 判定**不看 miRNA** 却把 `passes` 硬编码 True；② psRobot docstring 声称 `gap=2.0` 罚分（等长滑窗从不产生缺口，永不生效）；③ `rnahybrid_scan` 返回对象 `passes` 恒 False（死函数） | ① `RNA22Result` 增字段 `miRNA_independent: bool = True`（判定只分析靶标自身结构：GC 富集种子+hairpin），docstring/返回处标注，`to_dict` 导出 `rna22_miRNA_independent`；② psRobot docstring 删除 gap 声明并注明"num_gaps 恒 0"，`PSROBOT_SCORING['gap']` 加"未实现"标注；③ `rnahybrid_scan` docstring 标注 **deprecated（死函数）**并指向实际使用的 `consensus._validate_rnahybrid` | rna22.py:74-95（dataclass）、:101-113（docstring）、:157-159（返回标注）；psrobot.py:12-27（docstring）、:196-202（表注）；rnahybrid.py:44-49（deprecated 说明） | ① `validate_with_rna22`（构造 GC 富集种子+hairpin 窗口 `AAAAAAAACCCCGGGGAAAA`）→ `passes=True, pattern_score=7.0, miRNA_independent=True`；pytest `test_rna22_flagged_mirna_independent` 通过；② consensus 汇总处 `engine_notes` 含"RNA22：仅靶标结构特征…未使用 miRNA 序列"（pytest 断言 + 示例端到端实测），且 summary 透出；③ psRobot docstring/表注 grep 确认 gap 声明已删；④ rnahybrid 死函数保留但已标注（删除会破坏旧调用名兼容） |
| 6 | P2 | `consensus.py:96,184,216-221` | `Optional` 未导入（F821×7，`get_type_hints` 必炸）；`from __future__` 写在字符串**之前**导致 `__doc__` 为 None | ① 模块 docstring 提到最前、`from __future__ import annotations` 紧随其后，`from typing import Optional` 补上；② 顺带 `import logging` + 模块级 `log`（第 4 项用） | consensus.py:1-13 | 修复前实测：`__doc__ = None`、`get_type_hints(_validate_psrna)` → `NameError: name 'Optional' is not defined`；修复后：`__doc__` 非 None、`get_type_hints` 对 `ConsensusResult`（23 个键）与 `_validate_psrna` 均通过；ruff `--select F` 该文件 F821 从 7 → 0（pytest `test_consensus_module_typing_and_docstring`） |
| 7 | P2 | `logan_trace.py:800-802` | `cross_top` 把**全局 Top5 物种**塞给每个 contig（与片段序号 s 无关）——不同 contig 显示同一列，属冒充 per-contig 归属 | 改为按该 contig 关联片段**已导入的结果行**统计 Top5（多片段时按匹配样本数累加）；表头改"LOGAN Top 物种（本 contig）"、正文注明"无导入结果显示 —" | logan_trace.py:822-835（统计）、:886-894（表头/文案） | 代码核实：`agg['per_seg'][seg['index']]` 的每行确有 `organism` 字段 → per-contig 归属**数据支持**，故做真 per-contig（非改文案）；渲染路径未跑真报告（需真实 Logan 结果表），以源码锚点 + 数据契约核对（`_load_rows`→`organism`）为准 |
| 8 | P2 | `bin_engines.py:100-106` | miRanda `Energy:` 行在窗口外：官方 v3.3a 输出块含**空行**，Energy 实际在第 6 行（i+6），原窗口只到 i+5 → 能量恒 None；同因，比对三行（Query/杠线/Ref）也取错（取到空行/标签行） | 窗口放宽到 `i+8`；比对三行按实测格式取 `i+2/i+3/i+4`（含注释说明 7 行块结构） | bin_engines.py:145-165 | 真 miranda.exe 单条短序列：`energy=-39.0`（修复前 None）；合成块解析 `Energy=-39.0`、`query_aln` 为真 `Query:` 行（pytest `test_miranda_block_parsing`） |
| 9 | P2 | `bin_engines.py:89-90` | Q 区间解析取到 Score 数字：`_ints_after(body.split('R:')[0])` 在 `"Score: 195.000000 Q:2 to 20"` 上得 `[195,2]` → `q=(195,2)` | 先按 `'Q:'` 定位再 `_ints_after`：`body.split('Q:')[1].split('R:')[0]` | bin_engines.py:140-144 | `'Score: 144.000000  Q:2 to 20  R:6 to 26 …'` → **q=(2,20)**、r=(6,26)；真 miranda.exe 输出 → `q_start=2,q_end=20,t_start=6,t_end=26`（修复前 q=(195,2)） |
| 10 | P2 | `primer_thermo.py:115-118` | `calc_end_stability(fwd, revcomp(rev))` 参数方向错：把 R 摆成模板链方向，真实的 3' 端互配被算成**不结合**（正 ΔG） | 改为 `calc_end_stability(fwd, rev)`（primer3 语义：seq1 的 3' 端退火到 seq2，seq2 按自身 5'→3' 读入），docstring 记录实测对拍值 | primer_thermo.py:104-123 | 实测 F=`ACGTA` / R=`TACGT`（互为 3' 互补）：`calc_end_stability(F, R).dg = -2749 cal/mol`（结合）；旧写法 `calc_end_stability(F, revcomp(R)).dg = +1700 cal/mol`（不结合）；11-mer 用例 -3421 vs +116 cal/mol。pytest `test_end_stability_direction` 通过（断言新值<0、旧写法>0） |
| 11 | P2 | `bin_engines.py:179-198` | PITA 空格路径必挂：busybox `sh/cat` 以空白分词，实测 `cat: '.../Temp/vp pita space test/...'` 断在空格 → 命中 0 | ① 新增 `_ascii_space_free` / `_temp_free_base` / `_pita_root`：输入 fasta 目录含空格 → 迁到 `%TEMP%` 下 ASCII 无空格目录（跑完清理）；PITA_ROOT 自身含空格 → 最小文件集复制到临时目录（整目录 ~103MB，按路径哈希缓存复用）；② 产物仍在 pita_root，`PITA_PERL/PITA_SCRIPT` 改走运行时解析 | bin_engines.py:42-87（三个工具函数）、:231-300（`run_pita_batch`） | 真 PITA perl 管线 + **含空格 workdir**（`…\Temp\vp pita space test`）：修复前 0 hit（日志 `cat: 'C:/Users/17711/AppData/Local/Temp/vp': No such file or directory`）→ 修复后 **1 hit**（19.6s）；常规路径（无空格）行为不变 |
| 12 | P2 | `bin_engines.py:68-69,127-128` | miRanda/RNAhybrid 输入未做 U→T（PITA 侧已做 T→U），RNA 字母表输入被 DNA 引擎当非法字符 | 新增 `_to_dna()`（去缺口 + `U/u→T/t`，**保留原大小写**），miRanda/RNAhybrid 两个 `_write_fasta` 前统一调用；docstring 注明与 PITA 侧对称 | bin_engines.py:33-39、:125-126、:193-194 | `_to_dna('AUGCu-')='ATGCt'`、`_to_dna('UUU-ggg')='TTTggg'`（pytest `test_bin_engines_u_to_t`）；用 U 输入的 miRNA 真跑 miranda → 命中能量 -39.0 正常解析 |
| 13 | P2 | `suvtk_submit.py:57-77,362-366` | MIUVIG 元数据硬编码工具版本与项目名 `'Lycium virome project'`（项目私有字面量会被当真实项目名提交） | ① 新增 `default_miuvig(project_name, samp_taxon_id, log)`：`feat_pred` 从 **pyrodigal 的 importlib.metadata 版本**读取；`sim_search_meth` 优先读 suvtk 库目录 `db_info.json` 引擎记录、否则读 **mmseqs 二进制自身版本**；读不到保留字面量并 WARN；② `project_name` 默认值改中性占位 `plant virus metagenome project`，未显式覆盖时 WARN"必须人工核对/替换"；③ `build_sqn`/`/api/submit/sqn`/CLI 全链接收 `project_name`/`samp_taxon_id` | suvtk_submit.py:55-160（常量+`default_miuvig`）、:461-466（build_sqn 调用）、:678-684（CLI）；web/submit.py:410-421 | 实测：`feat_pred='pyrodigal;3.7.1;single mode'`（环境读取）、`sim_search_meth='MMseqs2;8cc5ce367b5638c4306c2d7cfc652dd099a4643f;default parameters'`（读 mmseqs 二进制）、未传 project_name → 占位 + `⚠ MIUVIG 方法学字段为默认值…请人工核对` WARN；传 `project_name='My virome 2026'` → 原样采用、不再报 project_name WARN。**注**：项目内唯一 `db_info.json` 是 `annot_db/prot`（④b ORF 注释的 diamond 记录，与 BFVD/mmseqs 检索不是同一条链），直接采信会把引擎名写错 → 故未用它填 `sim_search_meth`（见"未修复项"说明） |
| 14 | P2 | `primer.py:181,256`（+ `main.py:894` / `web/tool_results.py:1190` 交互路径） | ① 特异性检查默认关且产物**无"未检查"声明**（host_hit 列全 NA 易被读成"已检查无脱靶"）；② 检查失败也写 NO（冒充已检查）；③ `-max_target_seqs 1` 低估脱靶；④ 交互式 primer_design 路径完全没接特异性参数 | ① summary 恒写 `specificity_checked`（两条出口都写；无引物出口也写）；② 显式 false 时日志 WARN"未做脱靶特异性检查"；检查异常 → `host_hit='FAILED'` + summary `specificity_checked=false`（跟踪 `spec_ok`）；③ `-max_target_seqs` 1→10；④ `primer_design.run_design` 增 `specificity=False, logger=None` 参数：先判 host-db 是否可用（不可用→如实"未配置宿主基因组"），可用才生成扩增子 FASTA 调 `check_specificity_blast` 并逐对回填 `host_hit`，summary 恒写 `specificity_requested/checked/note`；`/api/tool/primer_design` 透传 `specificity` | primer.py:164-168（max_target_seqs）、:238-244（窗口 WARN）、:256-262（无引物出口声明）、:265-330（WARN/FAILED/spec_ok）、:328（summary）；primer_design.py:376-382（签名）、:397-407（RuntimeError 兜底）、:416-490（特异性块+summary）；web/tool_results.py:1153-1195 | ① 交互路径打桩（伪造宿主库 + 假 blast）：5 对引物 → 扩增子 FASTA 5 条（`>t1_P1…`）、`host_hit` = YES×1 + NO×4、summary `checked=True`（note：1/5 命中）；② 默认关：`requested=False, checked=False`、note 含"未做"、**逐对不标注 host_hit**；③ 无宿主库但勾选：note="未配置宿主基因组（host-db）或未构建 BLAST 库" + WARN（不冒充已检查）；④ 检查中间抛异常 → `host_hit='FAILED'`、note 含异常类型（实测触发 KeyError 路径）；pytest `test_interactive_primer_specificity_declared` 通过 |
| 15 | P3 | `consensus.py:296-301,303` | `best_energy` 把**未通过过滤**的引擎能量一起取 min（拉低成非命中值）；`passes_pamirdb_filter` 恒 True 且不导出 | ① `best_energy` 只取通过过滤的引擎（`result.miranda/psrna/rnahyb`）；通过数=0 时输出 None（不再伪造 0.0）；② `passes_pamirdb_filter` 如实计算并**导出到 `to_dict`** | consensus.py:336-347（能量）、:406-409（导出） | pytest 断言 `best_energy == min(miranda_energy, psrna_energy)`（只含通过引擎）且与修复前 min 口径一致；`passes_pamirdb_filter` 出现在 to_dict（新增导出）；下游 `(r.get('best_energy') or 0)` 对 None 安全（tool_jobs:2333 已核实） |
| 16 | P3 | `alignment.py:216`（F841）、`:191-193` | `passes` 赋值后从不使用；`query_s` 死分支（`if query_s is None` 永假） | 删除 `passes=…` 死变量（tier_* 才是真实过滤）；query 坐标段重写为"k↔L-k+1 换算"（死分支随之消失） | alignment.py:225-233 | ruff `--select F`：该文件 F841 从 1 → 0（修复前 `git show HEAD` 抓取复核）；query 坐标在正控上 = (1,21) |
| 17 | P3 | `logan_submit.py:389,421,226-240,311-317` | ① `--email-rotate` 帮助声称"每 N 条顺序轮换"，实际是**每条随机轮换**；② Submit 后 600s 无 session 记 UNKNOWN 后不补漏、不汇总；③ curl 缺失要到下载时才炸；④ zip 落点含非 ASCII 路径 | ① 帮助文本 + 启动日志改如实（"每次提交随机轮换"，`--email-rotate` 标"[已弃用/不生效]"）；② UNKNOWN 记入 `unconfirmed`，主循环后兜底扫 localStorage 列出"未归属 session 候选下载链接"并写入条目 message，结束前 WARN 汇总清单；③ `probe_curl()` 启动探测（`curl.exe`→`curl`），缺失 `exit(3)`；④ 输出目录含非 ASCII → 先下到 `%TEMP%` ASCII 临时 zip 再拷回（`_ret()` 在所有返回路径清理） | logan_submit.py:64-82（probe/ascii 判定）、:322-337（中转+`_ret`）、:350-372（HEAD/下载走 `_CURL`）、:428-436（帮助）、:464-470（启动探测）、:513-514 / :584-614（UNKNOWN 汇总）、:650-656（结束 WARN） | ① `probe_curl()` → `D:\Program Files\Git\mingw64\bin\curl.exe`（本机 PATH）；`_ascii_space_free('D:/桌面/x')=False` → 本机（中文路径）必走中转分支；② 4 处 `return` 路径均替换为 `_ret(...)`（源码 grep 确认 9 处）；③ UNKNOWN/兜底扫描为新增分支，需真 selenium 会话才能端到端跑（本环境不跑浏览器），以源码 + `ast.parse`/py_compile 复核 |
| 18 | P3 | `logan_trace.py:595-621` + `suvtk_submit.py:173-181` | 取消/超时只在**读到新输出行**时检查：子进程静默挂起（等服务器/等 mmseqs）时 cancel 与 timeout 双双失效 | 两处都加独立看门狗线程（1s 轮询 → 越限即 `kill()`，阻塞读循环随即 EOF，循环后按原因报错）：logan_trace 轮询传入的 `cancel` 事件；suvtk 轮询"取消事件 + timeout"。**关键坑**：`utils._TASK_STATE` 是 `threading.local()`，子线程读不到调用线程的事件 → 新增公开助手 `utils.task_cancel_event()`，在任务线程先取出引用再传进看门狗 | logan_trace.py:596-611（看门狗）、:837-841（`_wd.set()` 清理）；suvtk_submit.py:247-281（看门狗）；utils.py:642-651（`task_cancel_event`） | 真子进程验证（假"静默挂起"命令）：① suvtk timeout=4s → **4.0s** 抛 `suvtk 子命令超时`（旧实现要等子进程自己结束）；② suvtk 取消（经 `utils.task_bind(cancel=ev)` 真实事件，2s 置位）→ **2.0s** 抛 `任务已停止（用户取消）`；③ logan_trace `batch_submit`（真实 Popen 挂起子进程 + 2s 取消）→ **2.0s** 抛 `批量提交已取消`（旧实现须等下一行输出，实测会等满 120s）。**修正记录**：首版看门狗直接调 `task_check_cancel()`（子线程读 thread-local 恒为无事件）→ 取消通路实测失效，已按上述改法修复并复测通过 |
| 19 | P3 | `primer_thermo.py:255-259,132-149` | Tm 盐条件与 primer3 不同尺度混用无声明；`gc_fraction` 对含 N 引物只算非 N 子集无说明 | 取值与字段名**保持不变**，`gc_tm_check` 输出新增两个旁注字段：`tm_note`（Biopython Tm_NN 默认盐条件 50mM Na/1.5mM Mg/0.6mM dNTP/50nM DNA，与 primer3 内部条件不同、仅供排序）与 `gc_note`（`ambiguous='remove'`：含 N 只统计非 N 子集）；docstring 同步 | primer_thermo.py:125-142（docstring+两字段） | 源码/字段核对：`evaluate_pair` 的 `out.update(gctm)` 会把两个字符串字段带进 `thermo`（前端 JSON 可显示，不参与数值计算——`compute_score` 只读 gc_fw/gc_rev/tm_*）；实测 `gc_fraction(Seq('NNACGT'))=0.5`（remove 口径）/ `ambiguous='ignore'`=0.333…，与旁注一致 |
| 20 | P3 | `primer.py:232-235` | plain 模式末端 <600bp 窗口被 `break` 静默丢弃（用户以为全长已覆盖） | 丢弃时 WARN 一次，含 contig 名、坐标与丢弃长度 | primer.py:238-244 | 源码确认 WARN 文案含 `{cid}` / 窗口范围 / `len(sub)`；此分支需 plain+短 contig 真跑（本环境未跑全量管线），以 py_compile + 分支审读为准 |
| 21 | P3 | `primer_design.py:28-33,395` | `_p3()` 不可用时抛 `RuntimeError`，而 `run_design` 只接 `OSError` → 交互页整体 500 而非逐序列软跳过 | `except (OSError, RuntimeError) as e:`，提示文案改"primer3 不可用或参数不合法" | primer_design.py:397-407 | py_compile 通过；类型语义核对（`_p3()` 抛 RuntimeError → 现被捕获转 `explain`）；pytest 中 primer3 可用路径回归通过 |
| 22 | P3 | `tool_runs_admin.py:139-159` | organize 无单条容错：一条移动失败（占用/权限）整批抛异常中断，且无成功/失败统计 | 每条移动包 `try/except OSError`（含 `makedirs`），失败条目记入 `failed` + WARN，末尾输出"成功 N 条 / 跳过 N 条 / 失败 N 条"，返回值新增 `failed` | tool_runs_admin.py:128-175 | py_compile + 逻辑审读（`failed` 为新增键，`main.py:650` 只读 `moved/skipped` 不破）；`dry_run` 路径不变（不建目录不移动、仍计 moved） |
| 23 | P3 | `tapir.py:180,187,227-240`；`rna22.py:66-69` | tapir：`revmir` 死变量、种子窗 `12..23` 死算式、`try` 外引用 `mfe2/mfe1` 的 NameError、`import RNA` 未用；rna22：嵌套重复 if | tapir：`tapir_predict` 改用反转帧 + 正确种子窗（`SEED_POS_START/END`），`mfe1=mfe2=0.0` 预置（try 外引用安全），删 `import RNA` 与死变量；rna22：`if score_seed(rc)>4: if find_hairpin(...)` 合并为单条 `if … and …`（语义不变） | tapir.py:189-270；rna22.py:59-66 | py_compile；`rna22.RNA22()` 行为回归（同输入同输出，纯合并）；tapir 用 pytest 的三条用例覆盖（seed 窗/短 target/比对列） |
| 24 | P3 | `psrnatarget.py:305,320-321` | ① revcomp 前未 `_dna()`：含 U 输入经 `_revcomp_dna` 的未知碱基映射整条变 N；② 择优选中更差比对（实测完美命中被按 25.0 计分而落败）；③ 坐标区间倒置（`t_start=22, t_end=1`，start>end） | ① `_revcomp_dna(_dna(target_seq))`；② 维持"通过优先、同状态 expectation（越小越好）者胜"的比较方向（方向本身正确，选错的根因是第 3 项的计分 bug），并在注释中记录该结论；③ `_sw_align` 额外返回比对覆盖区间 `(t_lo,t_hi)`，`psrna_score_alignment` 新增 `target_range` 参数据此写坐标（未给时退化为旧口径但 `start≤end`） | psrnatarget.py:80-117（返回区间）、:239-268（target_range）、:311-350（evaluate） | ① `_revcomp_dna(_dna('ACGU'))='ACGT'`；U 版 target 与 T 版结果一致（坐标 6-26、passes=True）；② 修复后 forward 完美命中 penalty=0 被正确选中（修复前选中 reverse 链 penalty=10.5 的次优比对）；③ 坐标 **(6,26)** 精确对准插入段（修复前 (22,1)）；pytest `test_psrna_evaluate_coords_and_u_input` 通过 |

### 联动项（第 1 项直接引起的跨文件改动，超出清单但必须同步，否则前端/引擎会解析错帧）

| 编号 | 级别 | 位置 | 问题摘要 | 修复方式 | 改动位置(修复后) | 验证方式与结果 |
|---|---|---|---|---|---|---|
| A2 | 联动(P0 配套) | `webapp/static/app-mirna.js:283-309,504-521` | 前端按旧帧写死解析比对串：`^miRNA\s+5'` 正则匹配不到新标签 → 整行当序列渲染；Query2 过滤（10-11 位 WC + A/U）在新帧下从 3' 端数位次 | `mirAlnHTML`/`mirIsQuery2` 改为**按上行标签定帧**（3' 起 → 种子位/位次自右向左换算），并兼容旧格式；Query2 的配对判定改用碱基互补（原只看 `'|'`，会把 G:U wobble 误判为 WC，与后端切割口径不一致） | app-mirna.js:282-310（渲染）、:503-527（Query2） | `node --check webapp/static/app-mirna.js` 通过；`python tests/_check_js_syntax.py` → "前端 JS 语法检查通过 ✔" |
| A3 | 联动(P0 配套) | `mirna_target/alignment.py:198` | 靶起点整体左偏 1（`t_start=j` 而非 `j+1`）：即便方向修正，区间仍会指错 1 bp | 改 `t_start=j+1`（回溯止步列即"最后一个未参与比对的列"） | alignment.py:196-198 | 正控 `target=6-26` 精确等于插入段（len('GGGGG')+1 .. +21）；示例真命中位点能量与 duplexfold 逐位一致（见"结果变化说明"第 2 条） |
| A4 | 联动(P0 配套) | `consensus.py:37-91` | `compute_seed_quality` 依赖比对串且只按 `'|'` 分类：本模块 `'|'` 同时覆盖 WC 与 G:U wobble → 字段名 `seed_wc`（承诺 WC 计数）会把 wobble 计成 WC | 改双行解析（miRNA+target）、自右向左数位次、逐列用 `_is_wc/_is_wobble` 分类；窗口保持本模块既有 **2-8** 口径（不在本次改窗口，见"未修复项"） | consensus.py:37-91 | 完美位点 → `seed_wc=7, seed_wobble=0, seed_perfect=True`；人工构造种子第 3 位 (G:U) wobble → `seed_wobble=1, seed_wc=6, is_seed_perfect=False`（修复前 wobble 会被算成 WC=7）；pytest `test_seed_quality_new_frame` |
| A5 | 联动(第 18 项配套) | `utils.py:642-651` | `_TASK_STATE` 是 `threading.local()`，子线程看门狗读不到任务线程的取消事件（取消通路形同虚设） | 新增公开助手 `task_cancel_event()`（返回当前线程绑定的取消事件），调用方在任务线程取出后传进看门狗 | utils.py:642-651；suvtk_submit.py:252 | 见第 18 项验证（取消通路 2.0s 生效；修复前后对比实测） |
| A6 | 联动(第 4 项配套) | `consensus.py:358-410`、`web/tool_jobs.py:2623-2655` | 新增的 `engines_actual/engines_skipped/engine_notes` 需要端到端可见（per-row 与 summary 两级） | `to_dict` 导出三个字段（+ `passes_pamirdb_filter`、`rna22_miRNA_independent`）；tool_jobs summary 聚合（并集）+ 跳过时 WARN | consensus.py:400-409；web/tool_jobs.py:2623-2655 | 示例端到端跑通并落 summary（`engines_actual=['miranda','psrna','rnahybrid']` 等）；pytest 两条用例断言 |

## 结果变化说明（刻意的行为变化）

1. **第 1 项（P0-4 方向修正）——命中集合会变，这是修复的目的**
   - 示例端到端（`examples/example_mirna.fasta` × `examples/example_cmv.fasta`，8×3=24 对，快 3 引擎）：
     **修复前 2 条命中 → 修复后 4 条**。修复前 2 条中 1 条（Theobroma_cacao_miRNA_00667→NC_002035）
     在新方向下不再成立（旧命中来自镜像链判读）；另新增 3 条。
   - 新增命中已用 ViennaRNA 交叉验证：`Oryza_sativa-miRN73-5p→NC_002034`（t=822-843，能量 -26.7 ==
     `duplexfold` 位点值 -26.7；±15nt 扫描最优点在 823）、`Rhamphospermum_arvense_miRNA_00033→NC_002035`
     （-21.5 vs -21.6）、`Brachypodium_distachyon_miRNA_00736→NC_001440`（-23.1 vs -23.5）——位点区间与
     能量自洽（最优点落在区间 ±1nt 内）。
   - `examples/example_viral_target_sites.fasta` × 同一 miRNA 集（16 对）修复前后均为 **0 命中**：实测该
     演示文件的位点既非 `revcomp(miR159a/miR408-3p)` 也非 `comp(...)` 的精确子串（是旧平台的演示构造），
     对本次方向修正不敏感——如需"命中型"演示，建议另建含 revcomp 位点的示例文件（未擅自改动示例资产）。
2. **第 3/24 项（psRNATarget 计分与坐标）**：`score_alignment` 的 penalty 数值会整体下降（完美比对从
   列数→0），`psrna_score/evaluate` 的通过了率上升、靶坐标从倒置区间变为正常回文区间。下游
   `consensus._validate_psrna`（走比对列自打分，不经 `score_alignment`）不受影响。
3. **第 5 项（RNA22/psRobot/rNAhybrid 口径）**：RNA22 的 `miRNA_independent` 为**新增字段**（恒 True），
   不代表算法变化；psRobot 的 gap 声明删除（数值从未生效，无行为变化）；rnahybrid 死函数仅加标注。
4. **第 13 项（MIUVIG）**：`sim_search_meth` 由字面量 `MMseqs2;14.7e284;…` 变为**本机二进制实读值**
   （本机为开发构建 `MMseqs2;8cc5ce36…`）、`feat_pred` 由字面量变为实读 pyrodigal 版本；
   `project_name` 默认值由项目私有字面量变为中性占位 + WARN（提交前必须人工覆盖）。
5. **第 17 项（logan_submit）**：`--email-rotate` 从"声称生效"变为"标注已弃用"（行为本来就是随机轮换，
   未改代码逻辑）；UNKNOWN 条目会多出"候选 session 链接 + 汇总 WARN"；含非 ASCII 输出目录时 zip 会经
   `%TEMP%` 中转（产物路径与文件名不变）。
6. **第 15 项（best_energy）**：只在通过过滤的引擎中取 min（此前会把被过滤引擎的 mfe 拉进来）→ 该字段
   可能变"小"（更负）或从 0.0 变为真实值；示例命中 best_energy 由 -28.5/-25.0 变为 -26.7/-29.5/-22.9/-23.0
   （命中集合本身也变了，见第 1 条）。
7. **第 14 项（引物特异性）**：新增 `host_hit='FAILED'` 取值（检查异常时）与 summary 三字段；默认关路径
   的产物/日志多出"未检查"声明（数值不变）。

## 未修复项 + 原因

| 项 | 说明 | 原因 |
|---|---|---|
| 第 1 项：种子窗 **2-8 → 2-12** | 任务描述写"seed 判定：miRNA 5' 第 2-12 位"。核实：本模块（`alignment.py`）既有口径是 **2-8**（与 `consensus.compute_seed_quality` 同口径、与 miRanda 经典种子一致；2-12 是**tapir** 的 `get_upper_lower` 口径、2-13 是 psRNATarget 口径）。**未改窗口**，改为命名常量 `SEED_START=2/SEED_END=8` 并在 docstring 完整写清两套口径与反转帧换算（2-8 ↔ L-7..L-1；若改 2-12 即 ↔ L-11..L-1，已推导核对） | 改窗口会连带改变 `seed_ok → tier_strict → passes_filters`（miRanda 命中是否进入共识）从而改变结果，与"本次不改数值（会变结果）"的克制原则冲突，且与模块自身/frontend 的 2-8 口径不一致；如需 2-12 请作为独立科学决策项（改常量一处即可） |
| 第 13 项：`db_info.json` 引擎记录来源 | 审计建议"读 `db_info.json` 的引擎记录"填 `sim_search_meth`。核实：项目内唯一 `db_info.json` = `databases/annot_db/prot/db_info.json`，记录的是 **diamond（④b ORF 注释）**，与 MIUVIG 的 BFVD/mmseqs 相似性检索不是同一条链 | 直接采信会把"实际用 mmseqs 检索"写成 diamond（引入新的口径错误）；改为"suvtk 库目录 `db_info.json`（若存在）→ mmseqs 二进制实读版本 → 字面量 + WARN"三级回退，语义与来源一致 |
| 第 4 项：`engines_actual` 的"故障可见性"到前端 | 后端 summary/results.json 已如实落三个字段；前端表格（`app-mirna.js` 引擎状态表）仍按 `consensus_detail` 展示，未新增"引擎失败"提示行 | 清单只要求 tool_jobs summary 对接；前端工具页归"Web 前端"子系统（另一代理负责），跨文件边界未越界 |
| 第 7 项 / 第 17 项 / 第 20 项：真跑端到端 | ① `cross_top` per-contig 需真实 Logan 导入结果表 + 平台样品目录；② UNKNOWN 兜底需真 selenium 浏览器会话；③ plain 模式短 contig WARN 需真跑引物管线 | 本环境不跑浏览器自动化、不跑全量管线（审计约定：不跑外部引擎全量批处理）。三处均以"数据契约/源码锚点 + 相关单测"验证，py_compile 通过；建议在联网+浏览器环境的联调轮补跑 |
| 既有 lint 噪音（非清单项，修复前后一致） | `consensus.py:5`(`math` 未用 F401)、`psrnatarget.py:28,30`(`field`/`math`)、`rna22.py:19`(`field`)、`tapir.py:25`(`field`)、`bin_engines.py:308`(F541 无占位 f-string) | 经 `git show HEAD` 抓取修复前版本复跑 ruff 确认**修复前即存在**（同一批 6 条）。属既有噪音/死导入，按"纯风格问题不动"处理（第 16/23 项点名要删的 F841/未用 `import RNA` 已删） |

## 验证记录（可直接复跑）

```bash
# 1) 回归测试（P0-4 正控/反控 + 全部联动项）
python -m pytest tests/test_mirna_alignment.py -q          # 17 passed
python -m pytest tests/test_mirna_alignment.py tests/test_utils_pure.py -q   # 36 passed

# 2) 编译与 lint（对比基线：git show HEAD:<file> | ruff check --select F）
python -m py_compile Virus_Platform_Core/mirna_target/*.py Virus_Platform_Core/{primer,primer_design,primer_thermo,logan_trace,logan_submit,suvtk_submit,tool_runs_admin,utils}.py Virus_Platform_Core/web/{tool_jobs,tool_results,submit}.py

# 3) 示例端到端（8 miRNA × 3 病毒 = 24 对）：修复前 2 命中 → 修复后 4 命中
node --check webapp/static/app-mirna.js && python tests/_check_js_syntax.py
```

关键实测数值（正文引用出处）：

| 场景 | 修复前 | 修复后 |
|---|---|---|
| target 含 `revcomp(miR159a)` → identity / seed / tier_strict | 47.6 / False / False | **100.0 / True / True**（t=6-26 对准插入段） |
| target 含 `comp(miR159a)` → identity | 100.0（判中） | **61.9**（seed=False，不判中） |
| `RNA.duplexfold(mir, revcomp)` vs `duplexfold(mir, comp)` | — | **-36.5 vs -9.4**；比对能量 -36.50 == duplexfold |
| 官方 `miranda.exe` 对含 revcomp 的 target | — | `Score 195 Q:2 to 20 R:6 to 26 (100.00%)`（书写约定/坐标语义一致） |
| `psrnatarget.score_alignment(q,q)` | penalty 25.0 / mm 25 | **0.0 / 0** |
| `psrnatarget.evaluate` 靶坐标 | (22, 1)（倒置） | **(6, 26)** |
| tapir 第 3 位 / 第 16 位错配 | seed 判定错位 | **seed_mismatch=1** / **nonseed=1, seed=0** |
| `primer3.calc_end_stability(F=ACGTA, R=TACGT)` | +1700 cal/mol（判不结合） | **-2749 cal/mol** |
| PITA 含空格 workdir | 0 hit（`cat` 断在空格） | **1 hit**（19.6s） |
| `consensus.__doc__` / `get_type_hints` | `None` / `NameError: Optional` | 非 None / 通过 |
| suvtk 子进程静默挂起（timeout=4s / cancel@2s） | 等子进程自己结束 | **4.0s / 2.0s** 抛错 |
| logan_trace 静默挂起 + cancel@2s | 等下一行输出（≥120s） | **2.0s** 抛 `批量提交已取消` |
