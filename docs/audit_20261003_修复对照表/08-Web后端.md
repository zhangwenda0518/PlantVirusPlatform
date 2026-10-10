# 08 — Web 后端修复对照表（2026-10-03）

- 审计日期：2026-10-03；修复日期：2026-10-03
- 范围：`Virus_Platform_Core/web/*`（HTTP 边界与任务工厂）+ `app.py`（组合根/守卫/启动锁）
- 约定：位置行号指**修复前**（工作区存在其它子系统并行改动，行号可能小幅漂移，以函数名/锚点为
  准）；每项先读代码核实再改；全部为精准 Edit（未整文件重写，保持 UTF-8 与中文注释）；
  所有改动文件 `python -m py_compile` 通过；路由/校验类修复一律用 Flask `test_client` 冒烟
  （`python -B`，不起真实服务）。
- 回归基线：`python tests/_route_inventory.py` → 「✔ 路由清单与基线一致」；
  `python tests/_audit_routes.py` → 230 请求 / 500 错误 0 / 未捕获异常 0；
  `pytest tests/test_web_guard.py tests/test_utils_pure.py` → 24 passed；
  `python tests/_it_platform.py` → PLATFORM CHECKS PASSED。
- 本次未新增/删除任何路由（第 7 项用查询参数扩展），故 `tests/_route_baseline.json` **无需更新**
  （已跑基线脚本确认零差异）。
- 验证实验一律在平台 `run/` 下临时目录或系统 temp 内完成并清理（`run/uploads`、`run/results`、
  `run/tool_runs`、`run/meta_search` 已核对无 `_smoke*` 残留）；未做 git commit。

## 对照表

| 编号 | 级别 | 位置(修复前) | 问题摘要 | 修复方式 | 改动位置(修复后) | 验证方式与结果 |
|---|---|---|---|---|---|---|
| 1 | P1 | `web/tool_jobs.py:1088`（kvsuite `job` 内 Popen）；配合 `web/tasks.py:600`、`utils.py:611` | kvsuite 是全平台唯一裸 `subprocess.Popen` 的任务：未走 `task_register_proc`，取消无法硬停；且无超时看门狗，引擎挂死即永久占满 heavy 名额 | ① Popen 后立即 `task_register_proc(proc)`（`utils.task_register_proc` 顶层导入）；② 加守护线程看门狗：`DEFAULT_CMD_TIMEOUT_SEC`（默认 6h，`VP_CMD_TIMEOUT_SEC` 可调）到点置 `_timed_out` 并 `taskkill /F /T` 杀整棵树，读循环随即 EOF，收尾抛「超时…末几行」文案；`wd_stop` 在 finally 释放，取消路径不泄漏线程 | `tool_jobs.py:20-24`（import）、`:1152`（登记）、`:1155-1175`（看门狗）、`:1185-1199`（读循环 try/finally + 超时判定）；`tasks.py:600` 未改（`utils._task_register` 已有 cancel 竞态复查） | py_compile 通过。功能实验（`engine_cmd` 打桩为「写 pid 文件后 sleep 120s / 持续打印」的 python 子进程）：① 看门狗 4s 生效 → `RuntimeError` 含「超时」、耗时 4.6s（非 120s）；② 取消置位 → `RuntimeError('用户取消')`；两种情形子进程均被 `task_register_proc` 登记且 `tasklist` 查证已被杀（`被登记? True`／`已被杀? True`） |
| 2 | P1 | `web/samples.py:909`、`web/pages.py:398,406,415`、`web/tools_api.py:151`、`web/logan.py:39` | 5 条文件服务路由的 `<path:...>` 传绝对路径时 `os.path.join` 整段替换 base，`check_path(in_platform=True)` 只保证「平台内」→ 可读平台内任意文件（含 `run/ai_local.json` 的 API Key） | 照抄 `meta.py:1080` 样板逐条加校验：`base = os.path.abspath(check_path(base,…))`；`p = check_path(os.path.join(base, rel), must_exist=True, in_platform=True)`；`if not p.startswith(base + os.sep): abort(400,'非法路径')`。`pages.py:398/415` 本身无 `<path:>` 参数（原不可利用），按清单一并加包含校验（防御纵深） | `samples.py:944-953`（api_sample_file）；`pages.py:410-441`（page_report / page_report_file / page_archive_report）；`tools_api.py:151-170`（page_tool_run_file）；`logan.py:39-50`（page_logan_report_file）；`pages.py:14`（补 `abort` 导入） | py_compile 通过。test_client（造真实样品目录 / logan job 目录 / tool run 目录 + 平台内敏感文件 `run/_smoke_secret.txt`）：5 条路由传绝对路径（`.../run/_smoke_secret.txt`、`C:/Windows/win.ini`）全部 **400「非法路径」**；正常相对路径全部 **200**（内容正确） |
| 3 | P2 | `web/tasks.py:222-227` | 任务收尾在全局锁内做目录遍历（`build_result_preview` 全树 os.walk），大运行目录会把 `snapshot()`/`_progress()` 等所有轮询卡住 | preview 组装移到锁外先算好，再进锁做终态赋值（`status`/`result`/`result_preview` 仍一次性原子写入） | `tasks.py:220-231` | py_compile 通过；改动为临界区收窄，终态写入仍在同一把锁内（原子性语义不变） |
| 4 | P2 | `app.py:175-198` | CSRF 守卫只比 hostname 不比端口：本机其它端口上的页面（另一个本地服务/演示页）可驱动平台破坏性端点 | 新增模块级 `_BOUND_PORT`（`main()` 在 `_pick_port()` 后、起服务前写入）与 `_origin_same_platform()`：Origin/Referer 需 scheme=http 且 host∈回环 且 **port==绑定端口**；`_BOUND_PORT is None`（test_client/嵌入态）退化为原主机名校验，保证既有测试与本机脚本不受影响；无 Origin/Referer 的放行逻辑保持（拖拽代理 urllib 不发 Origin） | `app.py:172-219`（`_BOUND_PORT`/`_origin_same_platform`/守卫分支）、`app.py:506-511`（`global _BOUND_PORT` 赋值） | py_compile 通过；`pytest tests/test_web_guard.py` 5 passed（原有 5 例全过）。冒烟：`_BOUND_PORT=8765` 时 `Origin: http://127.0.0.1:9999` → **403**、`:80`/https → **403**、`Referer http://127.0.0.1:9999/…` → **403**；`Origin: http://127.0.0.1:8765` / `http://localhost:8765` → 非 403；无 Origin/Referer → 放行 |
| 5 | P2 | `web/build.py:829-1002`（build_kv_index）等 5 个建库类接口 | ① 无并发去重：同一目标连点两次起两个任务同写一个库目录（一个删旧索引另一个在建）；② `api_build_kv_index` 在**排队前**就删旧库：任务排队期间用户已失去可用旧库，任务失败则「旧的没了、新的没建成」 | ① 新增 `_busy_build_task(target,name)`/`_reject_busy_build(...)`：扫 `tm.tasks` 中 running/queued 且 `out_dir` 归一后相等或任务名相同 → 409（附带 busy 任务 id/名，便于去任务中心取消）；5 个接口提交时调用，并用 `out_dir=` 把目标目录登记进任务卡；② `build_kv_index`：提交时只保留「已存在且未勾 overwrite → 409」判定，删除动作移进 job 体（任务开始后才删，删失败 `RuntimeError`） | `build.py:31-60`（两个 helper）、`:534`（host_genome_download，按任务名）、`:586`（datasets_virus，有 name 时按目录否则按名）、`:701`+`:718`（convert_kraken2，含 falsy 口径）、`:951-965`（kv_index 覆盖清理入 job）、`:1046`（build_kv_index 去重）、`:1245-1247`（build_virus_db 去重） | py_compile 通过。test_client（`tm.start` 打桩，不真跑）：注入同目标 running/queued 记录后 5 个接口均 **409**（文案含 busy id）；host 下载用不同物种 → 200（不过度拦截）；build_kv_index 正常提交 200 且 `tm.start` 收到 `out_dir`、旧库目录**未被提前删** |
| 6 | P2 | `web/meta.py:556-567,141-157,705-768,887-928` | 元数据表编辑整表读改写（并发编辑互相整表覆盖）；AI 清洗任务持开始时的旧 DataFrame 收尾整表覆盖，静默抹掉用户编辑 | ① 模块级 per-name RLock（`_meta_edit_locks` + `_meta_edit_lock()`）与编辑事务上下文 `_meta_edit_ctx()`（校验请求→持锁→读表，锁覆盖「读-改-写」全程）；cell / row / rows_delete / fill_cache 改用它，import 整表替换持同一把锁；② AI 清洗：提交时记表 csv 的 `mtime` 为基线、在 `df.copy()` 快照上工作；收尾在锁内重读 mtime，与基线不同 → **不覆盖**，把结果另存 `<表名>.ai_result.csv/.tsv`（`_meta_store_df(stem_suffix=...)`，且不污染内存缓存）并返回 `conflict/ai_result/message` + 任务日志 `[WARN]`；一致才正常落盘 | `meta.py:49-87`（锁与编辑事务）、`:141-160`（`_meta_store_df` 支持后缀副本、返回 csv 路径）、`:611-700`（cell/row/rows_delete 用 `with _meta_edit_ctx()`）、`:735-800`（fill_cache 持锁；ai_fill 基线/快照）、`:832-855`（ai_fill 收尾冲突分支）、`:1000-1010`（import 持锁） | py_compile 通过。test_client（`_meta_ai_chat` 打桩，不联网）：① 编辑/整行/删行/连续 5 次编辑均 200，无死锁；② AI 清洗期间外部改表 → 任务 `done` 且结果含 `conflict=true`、写出 `SRA_GSA_Merged_Final.ai_result.csv/.tsv`、**表保留用户的编辑行且未被 AI 覆盖**（日志 `[WARN] 表在 AI 清洗期间被编辑过…`）；③ 无冲突 → 正常落盘（表出现 AI 结果） |
| 7 | P2 | `web/tasks_api.py:57-86` | 每任务一条 SSE 流，批量任务易触顶浏览器同源 6 连接上限（把轮询也堵死） | `/api/task/<tid>/stream` 支持多路复用参数 `?also=<tid2>,<tid3>`（逗号分隔、单连接上限 8）：一次连接轮询并推送多个任务快照，**多任务时**每条 data 事件带 `tid` 字段、结束发 `event: gone` 带 tid；所有任务到终态/消失才关流；docstring 写明用法。不传 `also` 时事件负载与旧版逐字节一致（前端零改动）；无新路由（查询参数扩展），基线无需更新 | `tasks_api.py:57-116` | py_compile 通过。test_client（真起 2 个 `tm.start` 轻任务）：单任务首帧**不含** `tid` 字段（向后兼容），`?also=<t2>` 时帧含 `tid` 且两任务都被推送，`also=` 20 项被截到 8；mimetype `text/event-stream` 正常 |
| 8 | P2 | `web/samples.py:307-322` + `app.py` | 多实例并发保护失效：`_RUNNING_SAMPLES` 是**进程内** set，多实例共享同一份 `run/`（queue.json、tasks/*.json）会互相覆盖状态与结果 | `app.py` 新增 `_acquire_instance_lock()`：`run/.platform.lock` 以 `O_CREAT|O_EXCL` 写入「pid 时间戳」；已被占且 pid 存活（Windows: OpenProcess+GetExitCodeProcess==259；其它平台 os.kill(pid,0)）→ 控制台打印「检测到平台已在运行（pid=…)」后退出；pid 已消失（残留锁）→ 提示接管并覆盖；锁不可创建时降级放行并提示。`main()` 最早处调用（端口探测前），退出/端口失败路径 `_release_instance_lock()`（pid 不匹配不动，防误删他人锁） | `app.py:396-488`（`_pid_alive`/`_acquire_instance_lock`/`_release_instance_lock`）、`app.py:491-496`（main 最早处调用）、`:502`、`:526`、`:547`（三处释放） | py_compile 通过。冒烟（把 `config.PLATFORM_ROOT` 指向 temp 目录，不落项目目录）：首次抢锁 → `('ok', path)` 且锁文件存在；锁内 pid=1（不存在）→ 残留锁接管 `ok`；锁内为本进程 pid（存活）→ `busy`（并打印提示）；`_release_instance_lock` 后锁文件消失 |
| 9 | P3 | `web/samples.py:231-234` | 队列把 cancelled 记成 failed 且 error 常为空 —— 卡片显示「失败」却无任何原因 | 三态显式映射：done → `status='done'`、`error=''`；cancelled → `status='cancelled'`、`error=快照 error 或「用户取消（<样品> 的任务已停止）」`；其余（failed / snap is None）→ `failed`，error 依次回退「快照 error → 任务记录已消失（可能被中断或清理）」 | `samples.py:236-256` | py_compile 通过。桩测试（`tm.start`/`tm.snapshot` 打桩 + 真实样品目录 + 免 sleep）：快照 done → `done/''`；cancelled(空 error) → `cancelled/用户取消（…）`；cancelled(带 error) → 原文；failed → 原文；snap=None → `failed/任务记录已消失…`。前端 `static/app-batch.js` 已有 `cancelled: '已取消'` 标签，无需改动 |
| 10 | P3 | `web/tools_api.py:390-399` | `out_dir` 只校验平台内，可建到 `results/` 等下变成幽灵样品（列表里删不掉、点开 404） | 提交时对 `out_dir` 加结构性目录黑名单（前缀命中即 400，文案含目录名与路径）：`DIRS['results']`、`DIRS['databases']`、`DIRS['tasks']`、`PLATFORM_ROOT/webapp`、`PLATFORM_ROOT/3rd`（两侧 `normcase(normpath(abspath))` 归一，避免大小写/分隔符绕过） | `tools_api.py:414-435` | py_compile 通过。test_client：5 个黑名单目录请求全部 **400「输出目录不能位于…内」**；`run/uploads/<新目录>` 通过黑名单（后续 400 来自 convert 工具自身的输入格式校验，非黑名单） |
| 11 | P3 | `web/tools_api.py:242-268` | 删除运行不检查是否有任务在跑：删掉正在写入的结果目录，任务随后写盘失败/产物残缺 | 删除前按 `out_dir` 匹配 `tm.tasks` 中 running/queued 任务（路径 `normcase(abspath)` 归一比较），命中 → 409（带 busy 任务 id/名） | `tools_api.py:255-268` | py_compile 通过。test_client：注入 `out_dir==运行目录` 的 running 任务 → **409**；置为 done 后 → 200 且目录被真删（`目录已删? True`） |
| 12 | P3 | `web/tool_results.py:1013-1021`（写）+ 1077-1078（缓存分支） | 结果 JSON 非原子写（并发读可拿到截断内容）；缓存命中不校验（截断/坏 JSON 会被永久返回且不重算） | ① 写入改先写 `out_path + '.tmp'` 再 `os.replace`（同卷原子）；② 缓存命中先用 `json.loads` 探针，`OSError/ValueError` → 删除坏缓存并继续走重算 | `tool_results.py:1080-1092`（缓存探针）、`:1136-1143`（原子写） | py_compile 通过。test_client：① 手工放入非法 JSON 的 `analysis/ctg1_local_primer.json` → `POST /api/tool/analyze` 返回 `{"cached":false,"task":…}`（不再假命中）且坏缓存文件**已被删除**；② 原子写端到端（`primer.design_primers_for_seq` 打桩）：任务 `done`、`analysis/` 仅 `ctg1_local_primer.json`、JSON 合法、**无 `.tmp` 残留**，二次调用返回 `{"cached":true,"task":null}`（合法缓存正常命中） |
| 13 | P3 | `web/tool_jobs.py:2744,2381-2388,795-833`（映射到修复前实际行号：`_phylodyn_resume`、mirna job、consensus job） | job 体内 `abort(400)` 在任务线程执行（无请求上下文），异常仅剩「BadRequest: 400 Bad Request」，提示文案全丢 | 用 AST 扫描全文件 job 子树，把 job 体内全部 `abort(400, …)` 改 `raise ValueError(…)/RuntimeError(…)`（文案逐字保留）：consensus 5 处、contigs 1 处、gb 1 处、mirna 3 处、`_phylodyn_resume` 1 处；另新增 `_job_req()`（与 `ctx.req` 同语义但不 abort），替换 3 处**可达**的 job 体内 `ctx.req`（identify 的 R2、pdprep 的 fasta/meta）。工厂期 abort 保留（仍在请求线程，400 语义正确）。扫描后 job 体内 abort **0 处** | `tool_jobs.py:54-69`（`_job_req`）、`:414`（identify）、`:650-654`（contigs）、`:812-830/869-891`（consensus）、`:2088-2089`（pdprep）、`:2124`（gb）、`:2479-2486`（mirna）、`:2867`（phylodyn_resume） | py_compile 通过；AST 复扫：`tool_jobs.py: job 体内 abort 0 处`；`app` 可正常导入、路由 231 条 |
| 14 | P3 | `web/tool_jobs.py:579`（当前 601-607） | 同一性矩阵用 `zip(cols[i], cols[j])`，比对列长不等时按短序列静默截断，identity 偏大且无提示 | 计算前先校验全部序列等长：`len({len(c) for c in cols}) > 1` → `RuntimeError`，报出两条序列名与各自长度并提示检查 MAFFT | `tool_jobs.py:600-610` | py_compile 通过；逻辑为「先校验再算」，`aln_cols` 仍取 `len(cols[0])`（校验通过后各列必然等长） |
| 15 | P3 | `web/tool_jobs.py:1815`（当前 1906 起） | `_tool_job_rtt` 丢弃 `analyze()` 返回的 `warnings`（未解析 trait、同名净化、无时间信号等），前端只看到一个斜率数字 | `summary` 增加 `'warnings': warns`，并在任务日志逐条 `[WARN] …` | `tool_jobs.py:1950-1962` | py_compile 通过；`phylogeo.analyze` 返回 dict 确含 `warnings` 键（`phylogeo.py:3352`），字段名对齐 |
| 16 | P3 | `web/settings.py:198-262` | 逐项就地改 cfg 后统一 save：中途某项失败时前面已生效（目录类 setter 内部还各自 save），配置「半套用」既在内存又在磁盘 | 改为两阶段：① 先整体校验（language/email/threads/defaults 只写本地 `pending`；四个目录根按同一套规则干校验：绝对路径 / `_ensure_not_system_dir` / 示例目录须存在 / 其余尝试 `makedirs` 验可创建）；② 统一应用 + 统一 save，任一异常 → `_restore()` 回滚全部已改项（含用旧值调回 setter）并回滚式 save；目录类失败 400、落盘失败 500，文案均注明「已回滚本次修改」 | `settings.py:198-300` | py_compile 通过。test_client（`cfg.save` 打桩计数）：① `language=en + threads=99999 + output_root=相对路径` → **400**，且 language/threads **未变**、save **未被调用**（无半套用）；② defaults 非法值 → 400 且 threads 未变；③ `output_root=C:\Windows` → 400 回滚（见文末备注）；④ 全合法 → 200，language/threads/defaults 全部生效 |
| 17 | P3 | `web/pages.py:398-423` | 样品 report 路由无 CSP 头（与 logan 溯源报告不一致）：报告 HTML 含样本衍生字符串与内联脚本，属半可信输入 | 新增 `_REPORT_CSP`（与 `logan.py` 逐字同口径）与 `_send_report()`（send_file + CSP 头），`/report/<sample>/`、`/report/<sample>/<path:filename>`、`/archive_report/<sample>/` 三条报告路由统一使用（归档报告同为 report.html，一并覆盖） | `pages.py:397-441` | py_compile 通过。test_client：三条路由均 **200 且 `Content-Security-Policy` 存在**（含静态文件路由） |
| 18 | P3 | `web/build.py:662` | `overwrite` 判定用 `is None`，前端传 `false`/`""` 时被当成「要覆盖」放行，与 `:1144` 的 falsy 口径不一致 | 改 `if not body.get('overwrite') and _db_dir_ready(db_dir):`（并加注释说明口径对齐） | `build.py:698-702` | py_compile 通过。test_client（`_db_dir_ready` 打桩为 True）：`overwrite=false` / `""` / 缺省三种均 **409**（修复前 false/"" 会放行） |
| 19 | P3 | `web/download.py:243-245` | 「转入分析流程」入队后未 `bind_queue_entry`：条目只在分析流程页可见，任务中心整批「消失」 | `_queue().add(...)` 的返回值逐条 `bind_queue_entry(it)`（延迟导入 `web.samples`，对齐 `/api/queue/add`） | `download.py:242-252` | py_compile 通过。test_client（`get_dl_manager`/`save_sample_input` 打桩，`enqueue=true`）：响应 `{"created":["SRR_SMOKE"],"queued":1}`，任务中心出现外部卡 **「排队·SRR_SMOKE」**（测试后已清卡与样品目录） |
| 20 | P3 | `web/chains.py:83-87` | 链内步骤 run 名为 `{step_key}_{秒级时间戳}`，同秒两次链会复用同一 run 目录互相覆盖产物 | run 名追加 4 位随机后缀 `uuid.uuid4().hex[:4]`（口径同 `tools_api.py` 的 `/api/tool/run`）；顶部补 `import uuid` | `chains.py:11`（import）、`:83-89`（run_name） | py_compile 通过。功能实验（临时 chain 目录 + 桩 job_factory）：同一秒内两次调用得 `identify_20261003_131907_2821` 与 `identify_20261003_131907_db74`（后缀 4 位、**不撞名**、两个目录均存在） |
| 21 | P3 | `web/io_api.py:67-69` | 上传只预检 Content-Length：`Transfer-Encoding: chunked` 且无长度时可无上限流式落盘（写满磁盘） | 选择「直接拒绝」分支（注释写明理由）：无 `Content-Length` 且头含 chunked → **411 Length Required**「上传请求需带 Content-Length（不支持 chunked 无长度上传）」；页面 fetch/XHR+FormData 与拖拽代理均带 Content-Length，chunked 只可能来自异常/恶意客户端，拒绝比读满 64GB 再断更省资源 | `io_api.py:65-76` | py_compile 通过。test_client：`EnvironBuilder`（无 Content-Length + `Transfer-Encoding: chunked`）→ **411**；正常 multipart 上传 → 200 且落盘（测试文件已删） |

## 未修复项 + 原因

1. **`Virus_Platform_Core/config.py:778-790` `_ensure_not_system_dir()` 实际是空操作**（第 16 项验证时发现，
   非本清单项）：`os.environ.get('SYSTEMDRIVE','C:')` 的值是 `C:`（无尾反斜杠），
   `os.path.join('C:', 'Windows')` 得到**盘符相对**路径 `C:Windows`，`normpath` 后仍无分隔符，
   于是 `c:\windows`.startswith(`c:windows\`) 恒 False —— 该守卫对
   `C:\Windows`/`C:\Program Files`/`C:\ProgramData` 全部放行（实测返回未抛错）。属 Core/config.py
   范围（本代理只负责 Web 后端；该文件极可能正被其它代理修改，避免冲突未动）。影响：设置页把
   `C:\ProgramData` 这类**用户可写**的系统目录设为数据根时不会被拦（只靠后续 makedirs 是否失败兜底），
   而 config 注释明确该守卫是防止 `write_roots()` 把系统目录纳入可写/可删区。建议负责人修复为
   `drv = (os.environ.get('SYSTEMDRIVE') or 'C:').rstrip('\\') + os.sep` 或直接 `os.path.splitdrive(os.path.abspath(path))[0] + os.sep`。
2. **第 5 项的「临时目录构建 + 原子替换」未做**（只做了清单允许的最小方案「把删除动作挪进 job」+ 并发去重）：
   `build_kv_index` 支持**分次补建**引擎索引（先 salmon 后 minibwa，读的是同一库目录里已有的
   另一引擎索引）。若改成「新目录构建后整目录替换」，补建第二批引擎时会丢掉第一批的索引，
   需先整库复制再构建 —— 代价与风险明显高于收益。现方案已消除「排队前删库」与「同目标并发」两个实际危害。
3. **第 7 项前端未改**（按分工由 Web 前端代理处理）：后端已支持 `?also=a,b,c` 多路复用，前端
   `webapp/templates/logan.html` 仍是一任务一条 `EventSource`（`subscribeBatch`）。**建议前端代理**：
   把同一批次的多个 tid 合并为 `/api/task/<首个tid>/stream?also=<其余tid逗号拼接>` 一条连接，
   收到事件时用 `snap.tid` 路由到对应卡片（`updateBatchLive`），全部任务结束后关流。
4. **第 21 项选择了「拒绝 chunked」而非「流式读取并硬上限」**（清单允许二选一，注释已写明）：
   本平台上传方（页面 fetch/XHR + FormData、拖拽代理 urllib）都会带 Content-Length，
   拒绝策略更省资源、错误信息更明确；若未来出现必须支持 chunked 的上传端，再补「边读边计数、
   超 64GB 即断」的实现（`request.stream` 分块写 + 长度累加）。
5. **第 6 项 AI 冲突副本的人工合并是「文件级」的**：冲突时前端绿色标记所用的 `cells`（行号+列名）
   仍指向 AI 快照的行号；若用户在任务期间删过行，行号与当前表不再一一对应，
   `.ai_result.csv` 需按 `Run`/`Sample` 等键人工比对（任务结果与日志已明确提示）。
6. **`run/.platform.lock` 与「自定义输出根」的关系**：锁固定位于 `PLATFORM_ROOT/run/`（清单要求）。
   两个实例若配置了不同的自定义输出根，仍会共享 `run/tasks`、`run/logs`、`run/ai_local.json`，
   故依旧视为同一实例互斥（这是有意的保守选择）；如需支持「多输出根并行」，需把锁改为按输出根分片
   （清单未要求，未做）。

## 备注（透明说明）

- `Virus_Platform_Core/web/tool_jobs.py` 与 `app.py` 在 git HEAD 中的 blob 是 CRLF，而**本次改动前后
  工作区文件均为 LF**（`app.py` 的 2026-10-01 备份、`dist/VirusPlatform-Source/` 的 2026-09-29
  打包副本均为 LF；`tool_jobs.py` 由 9-29 副本 162353B/LF → 9-30 修改后 163318B/LF，字节差 965
  与行数差 +16 一致，排除「该次修改把 LF 翻成 CRLF」）。因此 `git diff` 相对 HEAD 显示整文件差异
  是**既有状态**（`git diff --ignore-cr-at-eol` 可只看真实改动），不是本次修复引入的噪声；
  本次全部改动保持原文件换行符与 UTF-8 未变。
- 行号一栏「位置(修复前)」沿用审计清单编号（含其实际指向的函数）；第 13 项清单给的三处行号
  （2744 / 2381-2388 / 795-833）按函数锚点核实为 `_phylodyn_resume`、mirna job、consensus job，
  已全部处理。
