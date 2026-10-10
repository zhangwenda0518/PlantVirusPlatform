# dsRNA 设计恢复 + 独立审计 + 修复优化 对照表（2026-10-03）

范围：dsRNA 设计子系统（引擎三模块 1846 行 + web 接线面 + 测试/文档/示例）。
背景：2026-09-29 dsRNA 设计整体迁出（工作区级改动，未提交 git），核心五件套
先由主会话从 HEAD 恢复到工作区；本轮把**接线面按当前代码结构重新接入**（不整块
拷贝旧版），并对恢复的 1846 行核心代码做**独立审计**（体检层 + 审稿人层双线），
再修复/优化、全量验证。

约定：所有改动未提交 git（`git diff` 可整体复查）；核心三模块自 HEAD 恢复后
本轮**按审计结论修改**，故相对 HEAD 显示为 modified（内容见下表，可逐项回退）。
行号以本轮修复后的文件为准。

---

## 一、恢复清单（恢复了什么 / 从哪 / 怎么接线）

| # | 对象 | 来源 | 接线方式 |
|---|---|---|---|
| 1 | `Virus_Platform_Core/dsrna_pipeline.py`(829行) `dsrna_scoring.py`(430) `dsrna_offtarget.py`(587) | `git show HEAD:<path>`（HEAD = 53bc2863 收编版，其间无提交改动） | 主会话已恢复；本轮按审计结论修改（见第三节） |
| 2 | `3rd/tools/dsrnamax/dsRNAmax_det.exe`（v1.1.15，带 `-seed`） | 同上 | 直接可用；`config.detect_tools()` 恢复 `dsRNAmax` 探测条目 |
| 3 | 工具卡 `t-dsrna`（HTML + 内嵌 JS） | `git show HEAD:webapp/templates/tools.html` 的卡块与 dsrna JS 块 | 卡插在 `t-mirna` 前（原位置），去掉旧的 phylodyn 预填充 IIFE（现行版本已无该结构）；表格"最优窗"列改用**臂帧**坐标（旧运行回退 nt_*） |
| 4 | i18n 键 `tk.ds*` ×20（zh+en） + `hm.entry.detect.d` | `git show HEAD:webapp/static/i18n.js` | 插在 `tk.mirH2` 前（与 miRNA 同组相邻）；`alert` 文案改用既有键 `ui.tk.needFasta`（i18n 扫描器要求） |
| 5 | 工具注册 + 面板接口 | `git show HEAD:Virus_Platform_Core/web/{tools_api,tool_jobs}.py` | `_tool_job_dsrna` 重新入 tool_jobs（按现行规范重写：显式参数解析 / safe_open / 取消透传）；`TOOL_REGISTRY['dsrna']` + `/api/tool/dsrna_panel` 重新入 tools_api |
| 6 | `config.py` 的 `dsRNAmax` 工具条目 | HEAD | 探测表恢复（设置页/自检/platform.json 覆盖重新可见） |
| 7 | 示例：`examples/results/dsrna/`（10 文件）+ `tests/fixtures/dsrna/`（3）+ `tests/_it_dsrna.py` | 主会话已从 HEAD 恢复 | `examples.js` 的 `MAP['t-dsrna']`（`ds_input:CMV`）与 `RESULT_MAP['t-dsrna']='dsrna'` 恢复；`manifest.json` 的 dsrna 条目按 HEAD 原样回填；`make_example_results.py` 恢复 dsrna 生成条目（`panel:''`，约 3s） |
| 8 | 打包接线 | HEAD | `VirusPlatform.spec` hiddenimports 恢复 `dsrna_{scoring,pipeline,offtarget}`；`package.py` 恢复 `databases/dsrna/{lookup.db,siRNA_parameters.txt,lethal_lists}` 入包 + `dsrna/panel` 数据库树 |
| 9 | 文档回正 | HEAD | `README.md`/`README_EN.md` 工具箱行恢复 dsRNA 短语、④仅二进制表恢复 `dsrnamax`；`home.html` 入口描述与注释改回"已恢复"；`app.py` 注释恢复"dsRNA T7 引物"；`pages.py` 迁出注释改恢复记录 |

**未恢复/无需恢复**：`_nav.html`（注释中提及 dsRNA 属历史口径，仍正确）；
wiki `Home.md/Toolbox.md/Primer.md`（dsRNA 文档自 771daa67 起一直在库，本轮
逐条核对与实现口径一致，未改）；`kv_identify.py`（`dsRNA` 仅出现在分子类型
正则，与本子系统无关）。

## 二、接线清单（文件:行，修复后行号）

| 文件 | 位置 | 内容 |
|---|---|---|
| `Virus_Platform_Core/web/pages.py` | 147-155 | detect 组恢复 `t-dsrna` 卡 + 恢复注释 |
| `Virus_Platform_Core/web/tools_api.py` | 44 / 113-114 / 617-640 | import `_tool_job_dsrna` / `TOOL_REGISTRY['dsrna']` / `GET /api/tool/dsrna_panel` |
| `Virus_Platform_Core/web/tool_jobs.py` | 1881-1960 | `_tool_job_dsrna`（ctx.req 目标 FASTA；0-safe 参数钳制；cancel 透传；safe_open 写 summary.json） |
| `Virus_Platform_Core/config.py` | 956-959 | `dsRNAmax` 工具探测条目 |
| `webapp/templates/tools.html` | 2426-2494 / 4174-4338 | 卡 HTML（目标/排除/面板/臂长·kmer·迭代·种子·候选·臂数·错配档/范围·产物区间/分离株）与 JS（`loadDsrnaPanel` / `runDsrna` / `loadDsrnaResult` + 载入时拉面板清单） |
| `webapp/static/i18n.js` | 202 / 2362 / 1580-1599 / 3704-3723 | `hm.entry.detect.d`（中/英）+ `tk.ds*` 20 键 ×2 |
| `webapp/static/examples.js` | 57 / 346 | `MAP` 与 `RESULT_MAP` 恢复 |
| `webapp/templates/home.html` | 99-106 | 入口注释 + `data-i18n` 回退文案 |
| `README.md` / `README_EN.md` | 82 / 172（两版对称） | 工具箱行 + ④仅二进制表 |
| `app.py` | 54-55 | 注释 |
| `dev_tools/VirusPlatform.spec` | 154-161 | hiddenimports（dsRNA 三模块 + RNA） |
| `dev_tools/package.py` | 299-317 / 514-515 | 小数据文件入包 + `dsrna/panel` 数据库树 |
| `examples/results/manifest.json` | `dsrna` 条目 | 10 文件快照（HEAD 原文） |
| `tests/_it_platform.py` | 231-250 | 断言由"已迁出"反转为"已恢复"（卡/接口/组内三项） |
| `tests/_run_all.py` | 105-111 | 批量入编 `_it_dsrna.py`（300s 预算，实测 2.3s） |
| `tests/make_example_results.py` | 116-131 + MANIFEST_NOTES | dsrna 示例生成条目 + note |
| `tests/_route_baseline.json` | `+"/api/tool/dsrna_panel"` | `python tests/_route_inventory.py --save` 更新（228 路由，与基线一致） |

## 三、审计发现与修复对照表

审稿人层依据：与上游 dsRIP 参考实现（`抗病毒dsRNA分析/dsRIP/dsRNA_efficiency.py`）
逐函数对拍效价公式；与 C-host 侧 `safety_scoring` 口径（经模块 docstring 记录
的 2026-09-18 数值对拍）保持 `nt_*` 坐标系不动；扫描器以其自带 selftest 的
五条红线为验收。

| 编号 | 级别 | 位置 | 问题 | 修复 | 验证 |
|---|---|---|---|---|---|
| D1 | **P1** | `dsrna_pipeline.design_primers` / `window_overlap`（原 487-491、735-737） | **坐标系错位**：`pick_best_window` 的 `nt_start/nt_end` 是 guide（反义链）帧（dsRIP 口径），但引物按 `arm_seq[nt_start-1:nt_end]` 切**sense 臂** → 取到**镜像区间**（臂长 > 窗口+20 时错位）；QC 的"窗内占比"用同一组镜像坐标自比，恒 ≈1，掩盖错位 | `pick_best_window` 增 `sense_nt_start/sense_nt_end`（+`frame` 说明、`efficacy_best_window` 同步）；`design_primers`、`window_overlap` 改用臂帧（旧报告无 sense_nt_* 时按 `L-nt_end+1 .. L-nt_start+1` 就地换算）；前端"最优窗"列改臂帧（旧运行回退） | 离中窗实例 TEF-17（300nt）：guide 窗 4..280 ↔ 臂窗 21..297；修复后引物 scope=`best_window(臂 nt 21-297)`、切片 off0=20、引物落臂坐标 24..242（窗内）；**交付区间 guide 均分 89.046 vs 修复前镜像区间 88.649**（窗口自身 89.061）；CMV 示例 QC 窗内占比仍 1.0（现在是真测量） |
| D2 | **P1** | `dsrna_pipeline.run_chain`（原 2 处 `int(params.get('max_mm') or 1)`）、`safety_counts`（原 `scanres.get('max_mm') or 1`）、`_tool_job_dsrna` 原版 | **`max_mm=0` 被 `or 1` 吞成 1**：界面明确提供"0 错配"档（前端注释也写明不能用 `|| 1`），后端却把 0 当缺省 → 静默放宽为 ≤1 错配，1mm 命中也计入脱靶加权/红线 | 三处全部显式解析（0 为合法值、越界钳 0..2）；job 层用 `_int_param`（None/空串才取默认） | `_it_dsrna.py` 步骤 9（max_mm=0：11 窗=None、869 窗=0、mm1 计数 0）通过；`safety_counts` 阈值改显式 `mm_cap` |
| D3 | P2 | `dsrna_pipeline.run_dsrnamax_once`（原 `subprocess.run`） | 子进程**不可取消**（单次最长 3600s，多候选×多窗更久，用户点"停止"无效）；未设 `CREATE_NO_WINDOW`（GUI 模式弹黑框）；未登记任务进程（`TaskManager.cancel` 杀不到） | 新增 `_run_exe`：`Popen` + `task_register_proc`（平台口径）+ 读线程 + cancel/timeout 看门狗轮询 + 硬杀 + `decode_output`；超时可 `VP_DSRNA_TIMEOUT_SEC` 覆盖 | 集成测试全过（含取消=None 路径）；`-h` 探测同通道；取消/超时路径与 `run_cmd` 同构 |
| D4 | P2 | `probe_exe`（原 `text=True, encoding='utf-8', errors='replace'`） | 与平台解码口径不一致：dsRNAmax 报错/回显含中文路径时固定 utf-8+replace 会打成 U+FFFD，异常不可读 | 统一走 `_run_exe`（按字节 + `decode_output`：utf-8→GBK） | `probe_exe()` 正常返回 1.1.15/支持 -seed |
| D5 | P2 | `_mem_available_gb`（原只读 `/proc/meminfo`） | **Windows 恒 None** → 内存守卫"按可用内存 70% 收紧预算"在主平台是**死代码**，只剩固定 32GB 上限 | Windows 走 ctypes `GlobalMemoryStatusEx`（无新依赖），Linux 保持 /proc | 单测日志：`预算 24.1 GB`（=0.7×本机 ~34.4GB 可用），不再是固定 32GB |
| D6 | P2 | `run_chain` 面板三态（原 757-769） | 用户勾了面板但**一个都没匹配上**时静默降级为"未启用"（只有任务日志一行），report 里 `offtarget_panel=None` 与"本就跳过"不可区分 | 显式勾选而未扫成 → `logs` 追加 `⚠ …`（前端告警列表会显示），日志同步 | 代码路径 + 文案；面板接口/清单联动 |
| D7 | P2 | `dsrna_offtarget.build_species_index` / `build_lethal_index`（原无 cancel） | 真实面板 GB 级、建索引可能数分钟且**不可中断**；`scan_panel` 只在文件之间查 cancel | 两函数新增 `cancel` 参数并逐记录检查；`scan_panel` 传入 | 集成测试（cancel=None）+ 审阅；与 `scan_arms` 的取消语义一致 |
| D8 | P2 | 引擎全部写盘点（`report.json`/`primers.tsv`/`offtarget_lethal.tsv`/`arms.fa`/`amplicons.fa`/`targets_wN.fa`；job 的 `summary.json`） | 裸 `open`/`Path.write_text`，无平台写入区校验、编码/行尾口径与现行 `safe_open`（UTF-8 + LF + 平台根内）不一致 | 全部改 `safe_open`；`write_fasta` 同步 | 集成测试 + 端到端冒烟产物清单/内容与恢复示例一致（safe_open 只影响写入通道：编码/行尾/根内校验） |
| D9 | P2 | `tests/_it_dsrna.py` | ①硬依赖仓库外 GB 级面板（缺失即失败，且**从不在 `_run_all` 批量里** → 迁出期间无人报警）；②工作目录在 `%TEMP%`（真实运行不存在的配置，与 D8 后的写入口径冲突） | ①真实面板优先，缺失时在临时目录**合成小面板**并临时重指 `PANEL_DIR`（扫描代码路径照跑，命中数不参与断言）；②工作目录改平台 `run/` 下（gitignored）并清理；③数据面前置检查（缺 lookup.db 给出可执行说明）；④入编 `_run_all` | `python tests/_it_dsrna.py` → `DSRNA_IT_OK`（2.3s）；`python tests/_run_all.py _it_dsrna` → PASS |
| D10 | P2 | `run_chain` 的 `report['params']` | 复现口径不全：缺 `safety_essential_coeff` / `primer_num_return` / `isolates` | 补齐（+job 层参数钳制：`primer_product_min>max` 自动交换、kmer 12..32、候选 1..10 等） | 冒烟：report 结构对照仅"新增"，无缺失 |
| D11 | P3 | 三模块 ruff（`E4,E7,E9,F,B,PLE`）7 条 | B905（zip 无 strict）、B904×2（except 内 raise 无 from）、F541×2（无占位符 f-string）、B007（未用循环变量）、F401（未用 `import json`） | 逐条修正 | `ruff check` 三模块 **0 告警** |
| D12 | P3 | `ThermoLookup.__init__` | `sqlite3` 连接在 `dict(...)` 抛错时不经 close（异常路径泄漏句柄） | `try/finally` 关闭 | 审阅 + 单测 |
| D13 | P3 | `_run_exe` 管道 | `proc.stdout` 未显式 close | join 后 close（异常安全） | 审阅 |

**核实后不改（如实记录）**：`candidate_key` 对表格行解析失败静默 skip、`load_params`
忽略非法值（均为文档化设计，改会变语义）；`scan_panel` 的"双向前缀"宽松面板匹配
（改会变更选面板行为，属产品决策，见第六节）。

## 四、优化项（审计之外主动做）

1. **取消/超时链路贯通**（D3+D7）：exe 单次调用、面板索引构建、扫描循环三处均可中断；
   超时可用 `VP_DSRNA_TIMEOUT_SEC` 配置（默认 3600s/次）。
2. **job 层参数健壮性**：`_int_param` 统一显式解析（0 合法）+ 上下界钳制；产物区间反向自动交换。
3. **前端坐标口径**（D1 配套）：结果表"最优窗(nt·臂坐标)"；旧运行（无 sense_nt_*）自动回退。
4. **示例可再生产**：`make_example_results.py` 恢复 dsrna 条目（`panel:''`，约 3s）+
   `MANIFEST_NOTES`，示例结果与 manifest 条目可一键重生成。
5. **临时目录 ASCII 化 —— 实测不需要**（不做）：`dsRNAmax_det.exe` 对中文根、含空格
   与中文的文件名读写均正常（实测 `run/_probe_cn/dir with space/目标 文件.fa` → rc=0、
   CSV 正常落盘），且用户勾选面板与示例均不经过额外中转 → 不引入 `%TEMP%` 中转层
   （少一层清理面/少一处路径翻译）。
6. **断点续跑指纹**：本工具为单遍链（无断点语义），不做续跑；已有的面板索引缓存
   带 `mtime_ns + size + k + n` 指纹（任一不匹配即重建），符合审计要求。

## 五、验证记录（全部实跑）

| 项 | 命令 | 结果 |
|---|---|---|
| 语法 | `py_compile` 全部改动 py | 通过 |
| 语法（全库） | 编译 471 个自有 py（跳过 3rd/.git/node_modules） | **0 失败** |
| 静态检查 | `ruff check --select E4,E7,E9,F,B,PLE` 三模块 | **All checks passed** |
| 扫描器自测 | `python -m Virus_Platform_Core.dsrna_offtarget` | `SELFTEST_OK`（双链/分层/N 断窗/编码序/e2e 全过） |
| 集成测试 | `python tests/_it_dsrna.py` | `DSRNA_IT_OK`（端到端 9 段，2.3s；含确定性、内存守卫、max_mm=0、isolates、QC） |
| 批量入编 | `python tests/_run_all.py _it_dsrna` | PASS（2.3s） |
| 单元/回归 | `python -m pytest tests/ -q` | **75 passed**（无回归） |
| 平台集成 | `python tests/_it_platform.py` | `PLATFORM CHECKS PASSED`（含反转后的 dsRNA 恢复断言） |
| i18n | `python tests/_check_i18n.py` | ✔（模板锚点 1316 个全有定义、中英成对；裸 CJK 弹窗=0） |
| 前端语法 | `python tests/_check_js_syntax.py` | ✔（56/56 块语法通过；tools.html 无重复 id） |
| 路由守卫 | `python tests/_route_inventory.py` | ✔ 228 路由与基线一致（--save 已更新） |
| 示例清单 | `python tests/_check_examples_manifest.py` | ✔ 无硬性差异（dsrna 条目 + 10 文件对齐） |
| 端到端冒烟 | Flask test_client → `POST /api/tool/run{tool:dsrna}`（CMV 示例，同示例参数） | `SMOKE_OK`：产物 10 项与 `examples/results/dsrna/` **逐项一致**；臂序列与恢复示例**逐字节相同**（同 seed 同 exe）；report.json 结构差异**仅新增**（sense 帧字段 + params 补全）；summary.json 结构零差异；QC 窗内占比 1.0（真测量）；运行目录已清理 |
| 面板接口 | `GET /api/tool/dsrna_panel` | 200，`files: []`（本机未装真实面板）+ panel_dir 正确 |
| 中文/空格路径 | exe 直跑 `run/_probe_cn/dir with space/目标 文件.fa` | rc=0，CSV 正常 |
| 坐标系 | 离中窗 TEF-17 全链（见 D1） | 引物落臂窗内；交付区间 guide 均分 89.046 > 镜像 88.649 |
| 内存守卫 | 集成测试日志 | 预算 24.1 GB（Windows API 生效，非固定 32） |

## 六、未修复 / 遗留（如实披露）

1. **真实脱靶面板不在本机/仓库**（8 物种 RefSeq RNA，约 1.4G）：本轮以合成小面板验证
   panel 扫描**代码路径与产物结构**；真实面板上的命中数/性能未复现（数据面环境前置）。
   **2026-10-03 用户裁定：面板不再随软件/数据库包分发**（`package.py` `_DB_TREES` 的
   `dsrna/panel` 条目已移除，软件瘦身；UI 中英标签与 wiki/Toolbox.md 已注明"不随软件
   分发、自备放置"）。自备面板放入 `databases/dsrna/panel/` 即自动启用（含致死表命名约定）；
   不放置时三态逻辑自动跳过脱靶扫描与安全评分，效价链路不受影响。
2. **C-host 侧引物实现不在本机**（仓库已归档）：D1 的坐标系问题只能与本仓单测钉住的
   `nt_*`（C-host safety_scoring 同口径）对拍，无法与 C-host 的 primer 段端到端对拍；
   修复保持 `nt_*` 不变 + 新增臂帧字段，向后兼容。
3. **打包链路未真跑**（PyInstaller 全量 GB 级）：`VirusPlatform.spec` hiddenimports 与
   `package.py` 数据项按 HEAD 恢复，语法编译通过，未实际构包演练。
4. `scan_panel` 的宽松面板匹配（双向前缀）保留原语义：多物种面板且名称前缀相近时可能
   误选；改动会变更"选哪些面板"的行为，属产品决策，未动（仅把"完全没匹配上"升级为显式告警）。
5. `python tests/_check_spec_datas.py` 有 **1 项既有红灯**："源码树已无 .bak 手工备份"
   —— `Virus_Platform_Core/web/pages.py.bak_20261001_prescreen`（2026-10-01 前一轮残留，
   **非本轮产物**）；按该检查的提示应归档到 `archive/_bak_*`，但可能属并行工作留存，
   本轮未动（不在本任务范围）。
6. `10-入口与文档.md` 第 17 项 a)/b)（"删 dsRNA 短语、删 dsrnamax"）**已被本轮恢复所反转**：
   该文作为上一轮历史记录保留未改；以本文为准。

## 七、与既有审计表的关系

- 本轮**不动** `01`–`11` 分表与其中已闭环的结论；dsRNA 恢复带来的文档口径变化见
  `10-入口与文档.md` 第 17/64 行（已在本文第六节第 6 条声明反转）。
- 新增回归项已入编 `tests/_run_all.py`（`_it_dsrna.py`），后续迁出/恢复类改动会直接红灯。
