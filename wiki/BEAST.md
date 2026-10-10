# BEAST · BEAST 接口（交接包 · MJRM · MCC 解析）

> 返回 [Home](Home) ｜ 上一页：[Phylodynamics](Phylodynamics) ｜ 下一页：[Phylogeography](Phylogeography)

**本地平台不做 MCMC**——贝叶斯链交给服务器 BEAST 跑；平台负责**前后两头**：把输入打包送出去（beast_handoff）、把结果读回来核验复用（beast_io），外加两个专用解析器（mjrm / mot）。

涉及引擎：`beast_handoff.py`、`beast_io.py`、`mjrm.py`、`mot.py`、`bsp.py`（ESS 口径复用）。

> ✅ **现状（2026-10-08 接线完成，源码实读 + 端到端实测）**：本仓**已补齐作业接线**，导出与回读两条链都能跑：
> - **导出（上半程）**：新增工具 `beastexport`（`web/tool_jobs.py::_tool_job_beastexport` → `beast_handoff.export_job`），已进 `TOOL_REGISTRY`；「系统地理」卡新增「📦 打包 BEAST 作业包」按钮。产物落 `run/tool_runs/<任务>/beastexport/`（`aln.fasta` / `metadata.csv` / `README.txt` / `job.json`）。
> - **回读（下半程）**：`phylogeo.analyze` 的 `beast_log` / `beast_trees` 已由 `_tool_job_phylogeo` 透传（卡内两栏 + burn-in / 后验阈值 / 最多树数三个参数），`import_beast()` 真正可达；`summary.json` 的 `beast` 段也已透传 —— 前端 BEAST 面板与 AI digest 都读它。
> - 前端新增 6 个控件（`#pg_beast_log` / `#pg_beast_trees` / `#pg_beast_burnin` / `#pg_beast_min_prob` / `#pg_beast_max_trees` / `#pg_export_tag`）与 `#pgBeast` / `#pgExport` 两个结果位；`i18n.js` 的 `tk.pdTrack` 悬空文案已改为指向真实存在的两个控件。
>
> 接线前这里是「有模块、无接线」（`export_job` 无调用方、参数没透传、前端无控件）；本仓实现的是与姊妹平台「植物病毒进化分析平台」同一条接缝（对照见 `docs/进化动力学实现详解与对比_20260917.md` 第二部分，其 `_tool_job_pgjob` + `dev_tools/export_beast_job.py`）。

---

## 工作面分工

```
本地平台                          服务器 BEAST 1.x
──────────                        ────────────────
beast_handoff  作业包生成  ──→     跑 MCMC（链长 1×10⁸ 等）
beast_io       回读+核验   ←──     .log / .trees / MCC 树
mjrm           Markov jump 计数解析
mot            MCC 节点状态解析
（可选）beast_io.run_beast  本机装有 BEAST CLI 时可直接代跑一份 XML
```

## beast_handoff · 作业包生成

与 virome_phylo_pipeline 的交接通道：从平台产物（比对、树、元数据）生成 BEAST 输入作业包（XML 模板 + 数据 + 说明），交服务器或本机 BEAST 执行。

## beast_io · 回读与核验

| 能力 | 说明 |
|---|---|
| `probe_beast` | 探测本机 BEAST 可执行（beast/beast2/打包版 beast2.exe），读版本 |
| `read_table` | 读 .log 表格日志（优先 state 列头口径；mjrm 严格口径失败退通用口径） |
| `log_summary` | 逐列后验统计（mean/median/CI95/**ESS**）+ 自动分组（posterior/likelihood/prior/rate/indicator/jump/dwell）+ 低 ESS 告警 |
| `jump_corridors` | 跳转计数列（`c_A-to-B[1]`）→ **走廊后验表**（均值/中位/95% CI）+ 停留时间 |
| `tree_corridors` | `.trees`/MCC 树 → 走廊后验（sample/mcc 两模式）+ 树高直方 + 动画事件采样 |
| `xml_outputs` | 从 XML 读 `fileName=` 声明，定位全部产物路径 |
| `run_beast` | 可选代跑：stdout/stderr 逐行回吐、可取消、可超时 |

> ⚠️ **口径提示（平台原注）**：Markov jump 计数是**采样历史下的期望跳转次数**（不是树上边数、也不是 BSSVS 的速率指示），与 Fitch 观测事件数是两条口径——面板上分开摆。

## mjrm · Markov Jump Randomization Model（BEAST 1.x）

**只做三件本地能做的事**（跑链需 BEAST，其余全本地）：

1. **生成**：把 N×(N−1) 个 `A-to-B` 指示矩阵 + `<rewards>` 块插入已有 BEAST XML 模板，并**补上让计数真的出现在日志里的 `<log>` 引用**（实测 VirPhyKit 原版插法是"哑"的——矩阵声明了但没有 log 引用，跑几小时也拿不到计数；本模块顺手修复，BEAST 1.10.4 实测列名 `c_Region.count[1]` / `c_AS-to-OC[1]` / `c_AS_reward[1]`）；
2. **审计**：`audit_mjrm_xml` 跑 BEAST **之前**检查 XML 能否产出跳转计数——别白等几小时；
3. **解析**：`parse_beast_log`（裁烧弃）→ `jump_matrix_from_log`（迁移矩阵 + 停留时间，逐路线均值/中位/95% CI）→ `write_jump_csvs`（routes/matrix/dwell 三张表）。

另含 `make_dating_config`：从含离散性状的 XML 生成"只靠序列似然"的定年版配置；`build_pair` 一次生成定年/迁移两份配置。

## mot · MCC 树节点状态

MCC 注释树的节点状态/高度解析（供走廊动画、树高估计；与 beast_io 的 `tree_corridors` mcc 模式共用口径）。

## 界面与产物

- **界面**：不在独立卡片，而在「专项分析 → 进化动力学 → 系统地理分析」卡内的一段「🅑 BEAST 产物导入」区：上半程「📦 打包 BEAST 作业包」按钮（+ 作业名），下半程两栏（`.log` / `.trees`）+ 三个导入参数。结果面板的 BEAST 段给出后验样本数 / burn-in 丢弃 / ESS 诊断 / 树模式与用树数 / 走廊数 / 分带口径。
- **产物**：
  - 导出：`run/tool_runs/<任务>/beastexport/{aln.fasta, metadata.csv, README.txt, job.json}`（+ 可选 `iqtree.log` / `coords.tsv`）；`README.txt` 里是粘贴即用的服务器命令。
  - 回读：`beast_io` 的全量产物落 `run/tool_runs/<任务>/phylogeo/beast_log.json` / `beast_trees.json`；`summary.json` 的 `beast` 段是精简版（逐列诊断截断，逐树事件不进 summary）。
- **闭环**：本卡导出 → 拷到服务器跑（README 三条命令）→ 把 `merged/{mcc.tree, merged.log, merged.trees}` 拷回 → 填进本卡两栏 → 运行本卡，BEAST 面板与弧线分带即用真后验口径。

## 相关

- 本地 ML 轨道 → [Phylodynamics](Phylodynamics)
- 走廊 → 地图动画 → [Phylogeography](Phylogeography)
