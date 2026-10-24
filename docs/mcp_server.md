# MCP 服务器（virusplatform_mcp）

平台根的 [`mcp_server.py`](../mcp_server.py) 把平台暴露为 **MCP（Model Context Protocol）服务器**，
Claude Desktop、ZCode、Codex 等 MCP 客户端可直接驱动；
**网页端总览页也有「🛰 MCP 控制台」卡片**——点「连接」后从下拉选任意
`vp_*` 工具、按 schema 填参、就地运行看结果（后端 `/api/mcp/*` 以最小
stdio 客户端拉起本服务器，走同一 MCP 协议；首次连接才拉起进程）：

- **工具箱全量可执行**：与网页端同源的 `TOOL_REGISTRY`（32 张卡：质控 / 宿主去除 /
  鉴定定量 / 组装 / 建树 / SDT / 引物 / 基因组图 / dsRNA / CRISPR / CONAN-SWIFT（Cas3）/
  四条一键链……）
- **流水线全流程驱动**：`pipeline.run_analysis` 编排，⓪~⑩ 15 个阶段任选
- **任务与产物查询**：TaskManager 快照、运行清单、产物文本读取

任务引擎与网页端共用同一套代码路径（heavy/light 并发闸门、run.log、结果可视化补全都保留）。

## 快速开始

```bat
:: 自检（确认平台模块 / 双库就绪 / 工具卡数量）
3rd\python\python.exe mcp_server.py --check

:: 全链路自测（stdio 握手 → 清单 → 错误路径 → convert 真实运行）
3rd\python\python.exe dev_tools\test_mcp_e2e.py --pipeline

:: streamable-HTTP 模式自测（断线任务存活 + 重连轮询）
3rd\python\python.exe dev_tools\test_mcp_e2e.py --http 8899
```

首次运行缺 `mcp` 包会自动安装（锁 `mcp<2`，清华镜像，与平台启动脚本同一套兜底）。

## 两种传输模式

| | stdio（默认） | `--http [端口] [--token X]`（默认 8765） |
| --- | --- | --- |
| 任务生命周期 | 随客户端进程：客户端断开 → 服务退出 → 任务终止 | 服务端独立进程：**断开照跑、重连可查** |
| 接入方式 | 每客户端各自配置命令 | 局域网多客户端共享 `http://127.0.0.1:8765/mcp` |
| 鉴权 | 无需（进程级隔离） | 可选 `--token X`：所有请求须带 `Authorization: Bearer X`，否则 401 |
| 适用 | Claude Desktop / ZCode 单机挂载 | 长跑任务、多客户端、远程机器挂载 |

HTTP 模式参考 awesome-mcp-fastapi 的形态落地（官方 SDK `streamable-http`
传输），默认只绑 `127.0.0.1`；暴露到局域网请务必加 `--token` 并自备反代 TLS。

## 客户端挂载

**Claude Desktop**（`claude_desktop_config.json`，stdio）：

```json
{
  "mcpServers": {
    "virusplatform": {
      "command": "D:\\桌面\\植物病毒分析平台\\3rd\\python\\python.exe",
      "args": ["D:\\桌面\\植物病毒分析平台\\mcp_server.py"]
    }
  }
}
```

**ZCode / Codex 系（TOML，stdio）**：

```toml
[mcp_servers.virusplatform]
command = "D:\\桌面\\植物病毒分析平台\\3rd\\python\\python.exe"
args = ["D:\\桌面\\植物病毒分析平台\\mcp_server.py"]
```

**HTTP 模式**：先起服务 `mcp_server.py --http --token s3cret`，客户端按
URL 型 MCP 挂 `http://127.0.0.1:8765/mcp` 并配置请求头
`Authorization: Bearer s3cret`（支持 `"url": ...` + `"headers": {...}` 的
客户端直接填）。任务提交后即可断开客户端，回来重连 `vp_get_task` 续查。

## 工具一览（22 个）

| 工具 | 说明 |
| --- | --- |
| `vp_platform_info` | 环境总览：目录、双库就绪、就绪体检、**版本 + 启用工具清单 + 能力指纹**（首接先调） |
| `vp_list_tools` | 工具箱 32 张卡清单 + **逐参数 schema**（键/类型/必填/默认/枚举；`gated=true` 的卡已被下线，勿选） |
| `vp_run_tool` | 执行任一工具卡（提交即返回 task_id） |
| `vp_pipeline_stages` | 流水线 15 阶段键 / 中文名 / 默认全流程 |
| `vp_pipeline_start` | 启动样品全流程（run_analysis 全参数：kv_lib/组装输入与回退门槛/子采样/建树采样/清剪/验证组合/gbdraw 覆盖…） |
| `vp_get_task` | 任务快照（稳定投影）：状态 / 阶段 / 百分比 / 日志尾 / 终态结果（可截断） |
| `vp_list_tasks` | 近期任务列表（分页；`include_archived=true` 并入磁盘归档历史） |
| `vp_cancel_task` | 取消任务（当前外部步骤硬停止进程树） |
| `vp_get_task_log` | 任务全量日志尾（含归档任务三级回源） |
| `vp_list_samples` | 已接入样品清单（分页） |
| `vp_list_tool_runs` | 近期工具/链式运行清单（分页，带文件清单） |
| `vp_read_run_file` | 读运行目录内文本产物（tsv/json/log/fasta/nwk/svg） |
| `vp_delete_tool_run` | 删除运行目录（**必须 `confirm: true`**；运行中拒删） |
| `vp_ov_search_runs` | **OpenVirome（Serratus 云端）**：按物种/组织/地理等查 palmprint 阳性 SRA Run（公开 API 直连，无需凭据） |
| `vp_ov_counts` | 对一批 Run 做元数据分组统计（宿主物种/组织/地理/病毒科分布） |
| `vp_ov_mwas` | 对一批 Run 做 MWAS 关联分析（哪些因素与病毒丰度显著相关） |
| `vp_cli_help` | **37 个 CLI 命令**清单：用途 / sync·task 模式 / 逐参数 schema |
| `vp_cli_run` | 执行 CLI 快查（同步）：ictv-status/refs、ncbi-list、gb-list/check、submit-list/init/validate/fill/sbt、logan-jobs、selfcheck、ref-status |
| `vp_cli_submit` | 提交 CLI 长任务：**数据库构建**（init-taxonomy/build-host-db/build-virus-db）、**公共数据**（ncbi-dl/gb-dl/gb-import/meta-search/meta-info/meta-plot/host-genome/datasets）、**NCBI 提交**（submit-export/adopt）、**LOGAN**（create/import/batch）、**ICTV 更新**、report/host-analysis/orfa 重跑、db-migrate/tool-runs/samples-backfill |

## CLI 能力面（vp_cli_*）

`vp_cli_help` 返回全部命令的参数 schema（与 main.py argparse 同源；参数键
下划线自动转 `--中划线`，bool 传 true 即加旗标）。两类执行：

- **sync → `vp_cli_run`**：秒~分钟级快查，同步返回 `output`（尾部 ≤64KB），
  超时可调（10~3600s）
- **task → `vp_cli_submit`**：建库/下载/检索/批量提交等长任务，进 TaskManager，
  `vp_get_task` 轮询、日志实时可读；敏感键（`ncbi_api`/`deepseek_api`）
  在任务日志中**自动打码**

破坏面守卫：`db-migrate`（落盘切换平台库根）与 `tool-runs clean`（删超期运行）
必须 `confirm: true`（预检 `check`/`dry_run` 免）。被 gate 下线的命令随工具
整体摘除。

## 输出信封（2026-10-16）

查询/执行类工具返回统一信封：**`summary`（一句话）/ `result`（结构化数据）/
`interpretation`（规则化解读，由真实字段生成、不经 LLM）/ `next`（下一步
指引含工具名）**，原有关键字段全部保留。失败任务（vp_get_task）的
interpretation 会按错误关键词给出可操作建议（缺库/缺参/路径/超时…）。

## 与平台内置 AI（对话坞）的关系

内置对话坞的 LLM 现在也**能驱动 MCP**：工具循环新增 `mcp_tool_catalog`
（发现 19 个工具与参数口径）与 `mcp_call`（执行任意 vp_* 工具，结果自动
截断到 4KB 防撑上下文；`vp_delete_tool_run` 破坏面禁用）。即：外部客户端
经 MCP 协议驱动、网页对话坞经工具循环驱动，**两套 AI 共用同一个执行层**。

## 调用模型（重要）

长任务一律 **提交 → 轮询**，绝不同步阻塞：

```
vp_run_tool / vp_pipeline_start  →  { task_id, run, out_dir }
vp_get_task(task_id)             →  status: running → done / failed / cancelled
vp_read_run_file(run, path)      →  产物文本（限运行目录内，二进制拒收并提示路径）
```

`vp_get_task` 返回**固定键集合**的稳定投影（id/name/status/stage/pct/msg/
error/eta/out_dir/link/started_h/finished_h/log/result），键恒在、顺序恒定；
`result_max_kb`（默认 256）超限时 result 截断为 `{_truncated, _note, _head}`
并指回 `run/tasks/<tid>.json`。

## 工具下线（gate，借鉴 wisp-bio）

不删代码即可摘除指定 `vp_*` 工具，两个来源取并集：

- 环境变量 `VP_MCP_DISABLED_TOOLS=vp_a,vp_b`
- 平台根 `deferred.json`：`{"disabled_tools": ["vp_a"]}`

未知的 `vp_` 名忽略并告警；**deferred.json 存在但解析失败时拒绝启动**
（fail-closed——写了闸门文件就必须生效）。下线清单与能力指纹会体现在
`vp_platform_info` 的 `tools_gated` / `fingerprint` 字段里。

- 提交类工具缺必填参数 / 库未就绪会**当场**报中文原因（与网页端同文案）；
  `vp_list_tools` 的 params 字段即每张卡完整参数口径（31 卡约 150+ 参数），
  填 params 前先查它，不必猜键名
- 终态 `result` 注入了 `run` 与 `out_dir`；工具运行产物同时在
  `run/tool_runs/<run>/`（可被网页端历史区直接看到）
- 输入文件路径：绝对路径（电脑任意位置只读）或平台根相对路径；产物一律落平台内

## 注意事项

- stdio 模式任务在客户端进程内执行：断开即终止；HTTP 模式任务在独立服务
  进程内，与客户端解耦。两种模式与网页端同时运行时并发闸门各自计数
- `vp_read_run_file` 只读文本产物（≤4 MB，超出截断标注）；二进制
  （fastq.gz/bam/png…）请按返回的磁盘路径直接取用
- 未暴露删除类操作（删运行目录 / 清任务记录仍走网页端），`vp_cancel_task`
  只停任务不删产物
- 流水线取消语义：外部工具步骤（kunpeng/SPAdes/MAFFT…）被硬停止后，
  编排层收到失败即停；纯 Python 段之间的取消粒度与 CLI 一致
