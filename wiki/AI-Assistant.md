# AI-Assistant · AI 助手（对话驱动智能代理）

> 返回 [Home](Home) ｜ 上一页：[FAQ](FAQ) ｜ 下一页：无

对应引擎：`Virus_Platform_Core/ai_agent/`（`provider` / `digest` / `tools` / `chat` / `broadcast` / `executors` / `audit`，共 7 模块）。
Web 层：`Virus_Platform_Core/web/ai_assist.py`（9 条路由）；前端：`webapp/static/app-ai.js`（随 `_nav.html` 注入每个页面）；设置入口：设置中心「AI 助手」卡片（`/settings#aiCard`）。

回答一个问题：**「这一堆结果说明什么？下一步该做什么？」**

AI 助手是平台的**可选顾问层**：所有分析计算仍由既有引擎完成，助手只读结果、给出解读与建议，默认不改变任何判定。未配置模型服务时，平台功能一个不少。

---

## 三条设计原则

| 原则 | 含义 |
|---|---|
| **顾问而非替代** | 不参与任何分析判定；工具集**全部只读**（列样品 / 读摘要 / 读任务 / 读结果表），没有"重跑阶段 / 改配置 / 删数据"的通道 |
| **数字来自代码** | 检出物种、读段数、一致性百分比等**一切进入模型上下文的数值**，都由确定性代码（`digest.py`）从落盘产物算出后注入；模型只解读、不生成数字——科研场景防幻觉的底线 |
| **可插拔** | 未配 Key 平台功能不减；配 Key 后走 OpenAI 兼容协议，DeepSeek / 智谱 GLM / 通义 Qwen / Kimi / OpenAI / 本地 Ollama 通吃 |

## 界面入口

- 每个页面右下角的 **✨ 悬浮按钮** → 展开对话坞（`_nav.html` 全局注入，同一页只加载一次）；
- 对话坞内含：样品下拉（选定后该样品摘要注入上下文）、清空对话、跳设置、未读红点；
- 事件流 `/api/ai/events`（SSE 常驻），提问 `/api/ai/chat`（POST + fetch 流式）。

## 主动驱动三层

| 层 | 触发 | 消耗 | 说明 |
|---|---|---|---|
| 第一层 · 讲解 | 阶段**开始**时 | **零**（静态知识库，中英双语） | 推一张"这一步在做什么"卡片 |
| 第二层 · 播报 | 关键阶段**收尾**（可勾选）/ 阶段失败（恒播） | 少量 token | 后台线程取该样品摘要调模型生成"结果速读—注意—下一步"，发布为 `ai_read` 事件；AI 不可用时静默跳过，**绝不阻塞流水线** |
| 第三层 · 决策闸门 | 配置的阶段收尾后 | 少量 token + **阻塞** | 给出"继续 / 先检查"评估并暂停，等你在对话坞点继续或终止；**超时自动继续**（默认 15 分钟，可调 1~120）。⚠ 会改变流水线行为，默认关 |

事件类型：`stage_start` / `stage_done` / `stage_failed` / `ai_read` / `run_done` / `ai_gate`。

## 对话式追问

助手经 **function calling** 循环检索上下文后作答，工具集见 `ai_agent/tools.py`：

| 工具 | 作用 | 根白名单 |
|---|---|---|
| `list_samples` | 列出全部样品与各阶段完成度（done/total） | — |
| `get_sample_digest` | 读指定样品的结果摘要（各阶段状态、检出病毒、组装/注释/宿主/进化/引物关键数字与表头） | — |
| `get_task_overview` | 读任务中心概览 | — |
| `read_result_table` | 读结果表 | `run/results/<样品>/` |
| `read_tool_run_file` | 读工具运行产物 | `run/tool_runs/<任务>/` |
| `read_prescreen_file` | 读病毒预筛产物（`stat_run_summary.tsv` 等） | `run/meta_search/<项目>/prescreen/` |
| `propose_action` | **起草**动作（见下） | — |

三个读表工具共用同一份路径校验（`_safe_read_table_under`）：根与目标都过 `check_path`（平台内）+ 规范化前缀比对（防 `..` 与符号链接逃逸）+ 体积闸，行/列/单元格三重截断。预筛那条的根**只到 `prescreen/`**，不放开整个 `meta_search/`（那里还有原始下载件与检索缓存）。

## 执行权三级

| 级别 | 默认 | 行为 | 约束 |
|---|---|---|---|
| **第一级 · 执行单** | **常开** | 助手只能**起草**动作 → 待批卡片（含参数 diff 与理由）→ 你批准后经**与界面按钮完全相同的代码路径**执行；驳回即作废 | 参数走白名单 + 类型钳制，未知键直接拒绝 |
| **第二级 · 自动重试** | 关 | 阶段失败后按**原参数**自动重跑（即断点续跑机制本身） | 每样品每阶段限 1 次 |
| **第三级 · 自动调参闭环** | 关（需手动开启） | 失败后由助手在**限定键 + 数值边界**内提参并自动重跑 | 每样品限 1 次；失败回落第二级 |

硬编码安全闸（不受配置影响）：同一样品同时只允许一个运行；助手发起的运行**全平台并发上限 1**。

**永不进动作注册表的动作**：删除样品/数据、全局重置、修改 AI 配置自身、修改设置中心、平台目录之外的任何写操作。

## 审计留痕

`ai_agent/audit.py` 把全部动作与事件按天轮转落盘：

- 产物：`run/ai_audit/actions-YYYY-MM-DD.jsonl`、`run/ai_audit/events-YYYY-MM-DD.jsonl`（**保留 90 天**，超期自动清理）；
- 脱敏：`key` / `token` / `secret` / `password` 等字段递归替换为 `••••`，文本 / 列表 / 深度三重截断；
- **同步落盘**：执行单走同步写（进程被强制终止也不丢）；播报事件走后台异步写（不阻塞流水线）；
- **重启重建**：`list_actions` 从盘上重建历史，进行中的动作标 `interrupted`，动作编号推高防撞号；
- 读取：界面回看，或 `GET /api/ai/audit?kind=actions|events&limit=&day=`。

> 事件总线的**内存环形缓冲**（最近 100 条）重启即失；审计盘才是跨重启可追溯的来源。

## 设置项（设置中心 · AI 助手卡片）

| 字段 | 默认 | 说明 |
|---|---|---|
| 启用 AI 助手 | 关 | 总开关 |
| 提供方预设 | DeepSeek | DeepSeek / 智谱 GLM / 通义 Qwen / Kimi / OpenAI / Ollama（本地离线）/ 自定义 |
| API Base URL · Key · 模型 · 温度 · 最大生成长度 | — · — · — · 0.3 · 2048 | Key **分离落盘** `<平台根>/run/ai_local.json`（`platform.json` 永不落 Key） |
| 第一层 主动讲解 | 开 | 零 token |
| 第二层 结果播报（+ 触发阶段勾选） | 开 | 阶段失败恒播 |
| 第三层 决策闸门（+ 超时分钟 + 触发阶段勾选） | 关 | 1~120 分钟 |
| 第二级 自动重试 | 关 | — |
| 第三级 自动调参闭环 | 关 | 最高权限 |

## API

| 路由 | 方法 | 说明 |
|---|---|---|
| `/api/ai/config` | GET / POST | 读 / 写配置（Key 回显打码） |
| `/api/ai/test` | POST | 测试连接 |
| `/api/ai/context` | GET | 取当前上下文（样品清单 / 摘要） |
| `/api/ai/gate` | POST | 决策闸门：继续 / 终止 |
| `/api/ai/actions` | GET | 执行单列表 |
| `/api/ai/actions/resolve` | POST | 批准 / 驳回执行单 |
| `/api/ai/events` | GET | 事件流（SSE 常驻 / JSON 轮询快照） |
| `/api/ai/audit` | GET | 读审计盘（`kind=actions\|events`，`limit` 1~1000，`day=YYYY-MM-DD`） |
| `/api/ai/chat` | POST | 对话（流式） |

## 隐私与边界

- 出本机的数据**只有统计摘要与表格行**，不含任何序列本体；
- 助手与 [Public-Data](Public-Data) 的元数据清洗是平台**仅有的两个联网模块**，均需用户显式配置方才生效；分析核心（管道 ⓪→⑩）永不联网；
- 选**本地 Ollama** 可完全离线；
- 所有 AI 触碰的输出可经审计盘追溯；任何 AI 故障都静默降级，不影响分析任务本身。

## 相关

- 配置入口 → [Settings-Ops](Settings-Ops)
- 元数据清洗（另一个联网模块）→ [Public-Data](Public-Data)
- 后台任务与断点续跑（"重试"即复用此机制）→ [Pipeline-Overview](Pipeline-Overview)
