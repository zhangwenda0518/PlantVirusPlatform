# 开发笔记 · 踩坑记录与目录规范

> 本文档面向平台维护者：汇总开发过程中踩过的坑（及修复方式）、
> 目录规范约定、以及一轮全面回顾后的验证基线。使用手册见根目录 README.md。

## 一、踩坑记录（重要，改代码前先看）

### 环境与路径
1. **中文路径**：SPAdes 与 BLAST(LMDB) 不支持含中文/非 ASCII 的路径。
   平台的对策：SPAdes 经 `%TEMP%\vp_spades`（纯 ASCII）中转后拷回；
   BLAST 库建在 `%TEMP%\vp_blast`。新增外部工具调用时必须考虑这一点。
2. **Windows 端口保留段**：Hyper-V/WSL 会动态保留端口段，固定端口可能
   bind 失败（"访问权限不允许"）。`app.py::_pick_port` 逐个探测候选端口。
3. **控制台编码**：所有入口（app.py/main.py/Virus_Platform_Core.config）在 win32 下
   `sys.stdout.reconfigure(encoding='utf-8')`，否则日志中文乱码/报错。
4. **CRLF**：FASTA/FASTQ 写出一律 LF（`safe_open` 已统一），
   CRLF 会导致 SPAdes 等解析失败。
5. **orfipy 日志目录**：orfipy 默认把日志写到输入文件旁边的
   `<输入名>_out/` 目录（曾在平台根目录留下 `orfipy_viral_contigs.fasta_out`）。
   修复：调用时显式传 `--outdir` 到样品目录内（见 Virus_Platform_Core/orf.py）。
6. **kunpeng 大基因组建库 OOM**：convert 阶段 60 条/批整批载入，
   多条大染色体同批会触发 ~20.8GB 确定性巨型分配。对策：注入 taxid 时
   把每条序列切 ≤1MB 片段（相邻 34bp 重叠保 k-mer），见 Virus_Platform_Core/kunpeng.py。
7. **线程数与建库内存**：kunpeng 建库内存峰值 ≈ 1GB/线程，GUI 默认限 8 线程。

### 任务与并发
8. **幽灵任务**：服务重启后 running 状态任务永远不结束。对策：
   `TaskManager._recover_interrupted` 启动时把遗留 running 标记为 failed。
   **2026-09-11 补第二层坑**：`run/tasks/` 里除了任务记录还住着
   `queue.json`（SampleQueue 的批处理队列）。恢复逻辑原先只按 `.json`
   后缀筛文件，于是 ① 修剪到最近 100 个时会把最老的 `queue.json` 当旧任务
   删掉（整条队列静默丢失）② `json.load` 得到 list，`d.get()` 抛
   AttributeError 被外层 except 吞掉，**整个恢复循环中断、幽灵任务照样留着**。
   对策：文件列表加 `self._tid_ok(n[:-5])`（12 位十六进制才算任务）。
   `cancel_all/purge_all` 同样用 `_tid_ok` 挡 `queue.json`。
9. **同样品并发**：同一样品并发分类会写爆磁盘。对策：`/api/pipeline/<s>/run`
   拒绝同名样品的第二个运行中任务（任务名以 `<样品> ·` 前缀识别）。
10. **任务结果丢失**：早期 `snapshot()` 不返回 `result` 字段，任务完成后
    前端拿不到任何结果摘要（用户第 6 条痛点的根源）。已修复：
    snapshot/持久化均含 `result` 与 `result_preview`。

### 数据与安全
11. **路径穿越**：所有文件读写必须经 `utils.check_path` / `safe_open`
    （拒绝 `..`、输出限定平台根内）。新增 API 时禁止绕过直接 open()。
12. **NCBI 下载**：域名白名单（eutils.ncbi.nlm.nih.gov）+ 限速重试 +
    会话式翻页（借鉴 PhyloSuite），改 Virus_Platform_Core/ncbi_download.py 时保持这些约束。
13. **GUI 仅监听 127.0.0.1**，不要改成 0.0.0.0。

### 前端
14. **服务失联诊断**：连续 3 次轮询失败显示红色横幅（app.js `setConnBanner`），
    恢复后自动消失——修改轮询逻辑时保留该机制。
15. **LOGAN 批量面板**用 SSE（`/api/task/<id>/stream`）实时刷新，
    不整页重渲染（防止打断正在输入的邮箱框）。新增面板交互时同样
    避免全量 innerHTML 覆盖正在编辑的区域。
16. **进化树查看器用 Archaeopteryx.js**（`webapp/static/vendor/archaeopteryx/`，
    LGPL-3.0，npm 包 `archaeopteryx`）。渲染入口是 app.js `renderTreeTo()`：
    同一容器重新 `launchArchaeopteryx` 即换树，换树前先 `viewer.destroy()`；
    一页只能有一个 viewer。FastTree/RAxML-NG（以及 IQ-TREE 旧产物）把支持值
    写成内部节点名，
    解析需 `nhConfidenceValuesAsInternalNames: true`（此时必须同时
    `nhConfidenceValuesInBrackets: false`，两者互斥，见 app.js
    `_treeInternalLabelsAllNumeric()`）。面板自带布局/支持值/搜索/导出，
    PNG 导出依赖 `window.Canvg`（ESM 模块桥接），PDF 依赖 `window.jspdf` + svg2pdf。
17. **切模块就失忆**：一级导航全是裸 `<a href>`（无 SPA 路由），每次切模块都是
    整页刷新，页面级 JS 变量必然归零。表现最重的是 `/pipeline` 的当前样品：
    只有 URL 带 `?sample=` 才恢复，而导航栏的 `/pipeline` 不带参数；曾经写进
    `localStorage['vp_last_sample']` 的"上次样品"**全仓没有任何读取点**，
    于是"切回来忘了在分析哪个样品"。对策：`app.js` 新增全局上下文
    （`VP_CTX`，键 `vp_ctx_sample` / `vp_ctx_project`）+ 顶部导航常驻胶囊，
    变更经 `document` 上的 `vp-ctx` 事件广播给各页重挂自己。
    **改这一带时注意**：项目下拉的变更必须走模板 `onchange="gxProjChanged(this)"`，
    不能只在 `DOMContentLoaded` 里 `addEventListener('change')`——属性处理器在
    解析期就注册，页面自身的 onchange（会重刷下拉）先跑，会把用户刚选的值刷回去。
18. **`applyI18n` 会覆盖动态文案**：`i18n.js` 在 `DOMContentLoaded` 里对所有
    `[data-i18n]` 元素重写 `innerHTML`。想在静态标签后拼动态内容（例如
    「只看当前样品（NX-1）」）时，动态部分要放进**不带 `data-i18n` 的独立元素**，
    否则会被整段覆盖掉。
19. **顶部导航必须能"跟着内容长高"，且不能让它把品牌挤走**（2026-09-11 用户截图报障）。
    `.topnav` 原是 `height: 54px` 固定高度，而 `.navlinks` 是 `flex-wrap: wrap`；
    13 个一级链接约占 1150px，窗口一窄（实测 ≲1500px，尤其加了右侧全局上下文条
    之后）链接内部就换行，可绿色横条仍是 54px —— **第二行直接溢出到横条外的白底上**。
    `height: auto` 的兜底只写在 ≤720px，远低于实际窗口宽度，所以平时必然踩到。
    修的时候踩到两个反直觉点，都靠真浏览器量出来（`tests/_it_nav_responsive.py`）：
    - 横条自己**不能**用 `flex-wrap: wrap`：宽度不够时它会把整个 `.navlinks`
      挤到第二行（品牌独占一行、高度直接 72~80px），比原来更丑。正确做法是
      `flex-wrap: nowrap` + `.navlinks{flex:1 1 auto; min-width:0}`，
      让品牌与导航同处一行、由 `.navlinks` 内部消化换行；
      到 ≤1180px 再显式切成"品牌一行 + 导航整行"。
    - **末尾控制簇要绑成一组**：13 个链接刚好差几十像素时，flex 换行只会把
      最后的 📂 一个图标挤到第二行，孤零零挂在右边像坏了（实测 1884px 正是
      这个症状）。故 `app.js` 把「上下文条 + #langBtn + 📂」收进 `.navtail`
      （`flex-wrap: nowrap` 的单个 flex 项），要么整簇留第一行、要么整簇换行。
    - 量"行数"不能用 `round(top)` 去重：`align-items: center` 下同一行里不同
      高度的元素 `top` 本来就不同（把单行数成 5 行）。要按"中心 Y 分簇"。
    - 上下文胶囊宽度有限会截断长样品名，`title` 里必须给出**完整**名字。

### 打包
20. **PyInstaller**：webapp 模板/静态资源经 `_webapp_dir()` 多路径探测
    （源码同级 → 上级 → `sys._MEIPASS`）。新增静态目录时同步改
    package.py / VirusPlatform.spec。
21. `dist/VirusPlatform/` 是打包产物（含运行时 results/logs/tasks），
    不要把开发期临时文件混进去；`build/` 为 PyInstaller 中间目录，可随时删除。

### 破坏性操作与巡检工具
22. **破坏性接口必须要求显式确认，且巡检工具必须跳过它**（2026-09-11 实测事故）。
    `POST /api/global/reset` 最初只用请求体里的 `include_archive` 做可选项，
    空 POST 就会执行「取消全部任务 + 清队列 + 清任务记录」。而
    `tests/_audit_routes.py` 的设计是"对每条路由发一次空 POST"来找 500，
    于是它把 `run/tasks/*.json` 里约 100 条历史任务记录**全部删掉了**。
    该脚本当时的 docstring 写着"全程只读，不删除任何文件"——**这句话没有实现**：
    它只比对 tasks/ 目录的"新增"，完全不检查"删除"。
    两条修复，缺一不可：
    - 接口侧：加 `confirm: 'reset'` 必填令牌（**防误触，不是鉴权**），
      空 body 一律 400；
    - 工具侧：`_audit_routes.DESTRUCTIVE_POST` 显式列出跳过的破坏性路由
      （`/api/global/reset`、`/api/tasks/clear_finished`、
      `/api/queue/clear_finished` —— 后两条早就是空 POST 即删的状态），
      并把快照比对从"只查新增"改成"新增与删除都比对"。
    **新增任何破坏性路由时，必须同时登记到 `DESTRUCTIVE_POST`。**
    损失范围：只丢了任务中心的记录索引（名称/状态/日志尾部/时间），
    任务**日志文件本体**（`run/logs/tasks/` 698 个）与全部科研结果
    （results/ tool_runs/）都在——删记录不删日志是 `tasks.py` 的既有设计。

## 二、目录规范

**2026-09-10 目录整理**：顶层按「代码 / 数据库 / 运行期数据 / 外部依赖」四层
收敛，运行期数据进 `run/`、外部二进制与第三方进 `3rd/`，数据库保持原位。

```
<平台根>/
├─ app.py / main.py        GUI 与 CLI 入口
├─ Virus_Platform_Core/    核心 pipeline 包（GUI 与 CLI 共用，勿在 app.py 写业务逻辑）；
│  │                       2026-09-10 由 vp/ 更名而来（原 vp = Virus Platform 缩写）
│  ├─ known_virus_suite/   已知病毒识别与定量引擎（五段整合，自带 CLI）；
│  │                       同日由顶层 engines/known_virus_suite/ 迁入包内，
│  │                       模块间改用相对导入（from .kv_xxx）
│  ├─ web/                 HTTP 边界与任务工厂（Flask blueprint）
│  ├─ public_meta/         公共数据检索引擎
│  └─ ncbi_submit/         NCBI 提交准备
├─ webapp/
│  ├─ templates/           页面模板（Jinja2）
│  └─ static/              app.css / app.js / i18n.js / examples.js / vendor/
├─ databases/              kunpeng 库、taxonomy、tree_db、annot、misc（大数据）
├─ examples/               内置各工具示例数据 + results/（示例运行结果，与真实结果隔离）；
│                          2026-09-10 由 databases/examples/ 迁到根目录，整体纳管 git
├─ host-db/                建库源数据（宿主源基因组）
├─ databases/virusref_db/  病毒参考 + 派生索引（kv_index/ 鉴定库：salmon_k31/ + minibwa/）
│                          2026-09-11 由顶层 virus-db/ 迁入 databases/（命名随 _db 惯例）
├─ run/                    运行期数据（可清理重建，git 忽略）
│  ├─ results/<样品>/       每样品产物（00_prep … 09_genome_plots、logs）
│  ├─ tool_runs/           工具箱独立运行产物
│  ├─ logan/<查询名>/       LOGAN 溯源任务
│  ├─ tasks/ · logs/       GUI 任务状态与运行日志
│  ├─ uploads/ · submissions/ · meta_search/ · fastq/ · downloads/
├─ tests/                  自检/造数脚本与烟测数据
├─ docs/                   开发文档（本文档）
├─ scripts/                辅助脚本（build_virus_db.py / package.py / VirusPlatform.spec）
├─ 3rd/                    外部依赖（git 忽略）
│  ├─ bin/                 单文件外部可执行（kunpeng / seqkit / crabz / clustalw2 /
│  │                       muscle / aria2c / sracha；*.orig.exe 备份同放）
│  ├─ tools/               带目录结构的工具套件（Blast / mafft-win / FastTree /
│  │                       iQtree / trimAl / Gblocks / diamond / mmseqs /
│  │                       fastp / SPAdes / salmon2 等）
│  └─ open-virome/         Open-Virome 前端构建源（/virome 页 `_FRONTEND_BUILD`
│                          引用，勿改名/搬动）
├─ git-repo/               第三方仓库镜像（ViralConsensus 等，仅供参考）
└─ _archive/               归档区（历史产物、外平台工具、迁移备份）
```

**旧路径兼容（重要）**：`Virus_Platform_Core/utils.py::check_path` 会把旧的平台相对顶层目录
自动重定向到新位置（`tool_runs/`、`results/`、`uploads/`… → `run/`；
`tools/`、`bin/`、`vendor/`、`open-virome/` → `3rd/`），映射表见
`Virus_Platform_Core/config.py::LEGACY_TOP_DIRS`。因此历史数据、前端拼接的路径、旧文档示例
都无需手工改；**新增代码请直接用 `DIRS[...]`，不要再写字面量路径**。

工具路径探测集中在 `Virus_Platform_Core/config.py::detect_tools()`：先 `3rd/bin/`、
`3rd/tools/`，根目录平铺布局保留为回退（兼容已分发的 exe 平台）；
`platform.json` 的 `tools` 覆盖优先级最高。

**`platform.json` 的 tools 路径约定（2026-09-09 起）**：平台根内的工具写
**相对路径**（如 `bin/minibwa.exe`），由 `Virus_Platform_Core/config.py::Config._load()`
按 `PLATFORM_ROOT` 解析；平台根外的工具（Python Scripts、外部 SPAdes）
仍写绝对路径。这样整个平台目录搬到别的盘或改名，工具探测不会连带失效。

约定：
- **版本管理改用 git（2026-09-09 起）**：本仓库已 `git init` 并推送
  GitHub（`zhangwenda0518/plant-virus-analysis-platform`）。
  改代码前不要复制 `*.bak_<改动名>_<日期>`——历史用
  `git log -p -- <path>` 看演进、`git show <sha>:<path>` 取回旧版；
  `.gitignore` 已禁止提交 `*.bak_*` / `*.bak` / `*.orig`。
  2026-09-09 一次性清理了 74 个历史手工备份（含 `tools.html` 的 8 个版本），
  内容全部在 git 历史里（如 `git show 38243d2:webapp/templates/tools.html.bak_nodeid_20260909`）。
  回退：`git revert <sha>` / `git checkout <sha> -- <path>`。
- 临时/中间文件一律写进对应任务目录（results/<样品>/… 或 tool_runs/<run>/），
  禁止落在平台根目录（历史上曾出现 `C:` 空目录、`orfipy_*_out` 等遗留物）。
- `*.orig.exe` 等二进制备份保留在原处需注明用途，长期不用应移入 docs/ 归档说明。
  **当前唯一备份**：`bin/kun_peng-*.orig.exe` 是 kunpeng 打补丁前的原版
  （MD5 与在用 exe 不同，勿当重复文件删）。
  `bin/sracha.exe` 与 `bin/sracha-0.6.0-local-patch.exe` 也不是同一文件，
  前者为在用版本。
- minibwa 单文件 exe 统一放 `bin/minibwa.exe`（2026-09-09 前曾重复放在
  `minibwa_win_build/`，已合并）；构建说明在 `docs/minibwa-build/`。
- **`tools/strawberry-perl/` 不要清理**：它是 SNPGenie 的 Perl 解释器
  （便携版 5.42.3，`perl/bin/perl.exe`），被
  `Virus_Platform_Core/known_virus_suite/kv_variant_evo.py::_snpgenie_exe()` 引用。
  2026-09-09 曾因目录瘦身把它归档到 `_archive/`，导致 SNPGenie 链路断掉，
  现已恢复。使用要点见 `tools/snpgenie/DEPLOY_NOTE.md`（工作目录必须纯
  英文路径，且要显式传 `--workdir` 避开 MSYS 的 `pwd` 污染）。
- 日志只增不改；排查问题先看 `logs/` 与 `results/<样品>/logs/`。

## 三、验证基线（2026-09-05 全面回归）

- `python main.py tools`：全部内置工具探测通过（单文件工具在 bin/）。
- `python -c "import app"` / 各 vp 模块 import：无错误。
- Flask test client 冒烟：`/` `/pipeline` `/build` `/results` `/tools` `/logan`
  `/virome` `/help` `/settings` `/hostremoval` 均 200；核心 API
  `/api/tools` `/api/dbs` `/api/samples` `/api/tasks` `/api/settings` 正常。
- 独立模块烟测：`Virus_Platform_Core.msa_view`、`Virus_Platform_Core.phylo._safe_name`、`Virus_Platform_Core.pipeline` 视图函数、
  设置中心读写（platform.json `defaults`/`language`）。
- 上下游依赖链：fq2fa→host→virus→assembly→(hostana/orf)→orfa→
  phylo→primer→gbdraw→report 的依赖表见 `Virus_Platform_Core/pipeline.py::STAGE_DEPS`，
  管道页卡片按依赖解锁；工具缺失自动降级（fastp/seqkit/gbdraw/orfa）。

## 四、更新日志（第一轮 · 凌晨）

- 清理：根目录 `C:` 空目录、`orfipy_viral_contigs.fasta_out`、`__pycache__`、
  PyInstaller `build/` 中间目录、散落的 `.mimosa` 插件状态目录。
- 新增：设置中心（/settings，语言/线程/邮箱/默认参数/路径）、
  任务结果预览与完成通知、前端 i18n（zh/en）、orfipy --outdir 修复、
  `python main.py selfcheck` 环境自检命令。
- 修复（浏览器实机验证发现）：
  - settings.html 内联脚本重复声明 `const $`（与 app.js 冲突）导致整页脚本
    失效、默认参数网格空白——内联脚本一律复用全局 `$`，禁止重复声明；
  - LOGAN 页 `.step b` 选择器误伤正文内联 `<b>`（渲染成巨型绿色梭形），
    改为 `.step > b` 仅命中步骤编号，并删除模板内重复的本地样式块；
  - 静态资源缓存：app.css/app.js/i18n.js 引用统一带 `?v={{ asset_v }}`
    （app.js 的 mtime），升级后浏览器不会再用旧缓存（本次 LOGAN 修复
    曾因缓存迟迟不生效，即为该问题的实证）。
- 验证基线补充：Flask test client 全 9 页 + 9 API 通过；设置读写往返通过；
  真实任务（contigs 分类）完成后任务卡结果预览面板渲染通过（统计 chips +
  7 个产物下载链接）；`python main.py selfcheck` 全绿。

## 五、第二轮深度检查（2026-09-05 下午）

### 关键发现：坏修复被实测拦截
- **orfipy `--outdir` 修复本身是坏的**：相对路径输出目录 + 绝对路径输入时
  orfipy 直接退出码 1（首轮"验证"只测了绝对路径场景，漏掉了
  `predict_orfs → run_orfipy(viral_fa, out_dir)` 传相对路径的真实调用形态）。
  修复：`run_orfipy` 内部 `os.path.abspath(out_dir)` 后再拼参。
  **教训**：改外部工具调用参数后，必须用生产调用形态（相对路径）实测。
- **端到端全流程**：新建样品 E2E-VERIFY 用 CLI 跑完 13 阶段
  （fastp→fq2fa→host→virus→assembly→hostana→virome→orf→orfa→phylo→
  primer→gbdraw→report）全部 done、日志 0 ERROR，保留作为回归基线。

### 新增能力
- **阶段结果预览**：`STAGE_VIEW` 从仅 report 扩展到 9 个阶段
  （fastp 报告/host stats/virus 汇总/组装 contigs/宿主桑基/RdRP 汇总/
  ORF 注释表/引物表/主报告）；管道卡片「查看」按钮打开预览弹窗
  （TSV→分页表格、JSON→格式化、HTML→新窗口），带下载按钮。
- **样品管理**：结果页每行加「🗑 删除」；后端 `/api/samples/<s>/delete`
  三重防护（样品名规范化 + 必须 results/ 一级子目录 + 有运行中任务时拒绝）。
- **参数开放**：⑦ 卡片暴露 `tree_sampling`（blast/macro/genus/lineage）、
  ⑨ 卡片暴露 `plot_engine`（auto/gbdraw/dfv）——此前仅 CLI 可用。
- **任务状态文件自动修剪**：服务启动时 tasks/*.json 只保留最近 100 个。
- **i18n 补盲**：结果页 MSA/树视图说明与下拉/复选框标签接词典（16→30 处）。

### 目录清整
- 删除 9 个昨日开发期 tool_runs（保留最新 1 个作预览面板示例）与
  results/NX-5（0 文件的空骨架）。
- **results/NX-6（61GB）是真实样品数据，永不清整**；01_host_removal 52GB
  为其宿主去除中间产物，属用户数据。
- `.zcode/ui_shots/` 为验收截图存档，保留供对照。

## 六、布局重构（2026-09-05 晚）

问题：结果入口分散在 4 个页面（样品报告在"结果"、共线性图在"工具箱"、
溯源报告在"LOGAN"）；工具箱是杂物抽屉（单步工具与进阶比较混排）；
导航顺序不符合工作流。

调整（URL 全部不变，仅重组导航与页面内容）：
- 导航按工作流重排：总览 → 分析流程 → 结果中心 → 专项分析 → 数据与库
  → 公共病毒组 → LOGAN 溯源 → 设置 → 手册（更名三处）。
- 总览"分析工作流"入口按 建库→分析→结果→专项→溯源 排序。
- 结果中心新增"专项结果入口"卡：聚合共线性比较图（gb_collections）与
  LOGAN 溯源报告（logan jobs），一页可达全部结果面。
- 专项分析分两区：单步工具 / 参考数据与比较分析。
- 细节：.btn 补 text-decoration:none（聚合入口链接按钮带下划线）。

## 七、公共数据闭环 + 项目批处理（2026-09-05 下午）

新增两大能力（对应路线图第 7、8 项）：

### 公共数据下载（/download，导航"数据下载"）
- **解析双通道**：SRR/ERR/DRR → ENA filereport API（直接 FASTQ.GZ + 字节数）；
  CRR → NGDC 搜索页 → browse 详情页提取下载链接（cncb/big 双镜像按主站去重，
  ftp 统一改写 https）。CRA 项目级编号不支持，需先取得 CRR 列表。
- **引擎**：平台根目录内置 aria2c 1.37.0（-x4 -s4 -c 断点续传，解析进度百分比）；
  缺失时回退内置 HTTP（Range 续传，经 .aria_part 临时文件 + 完成后改名，
  避免 safe_open 对 .gz 后缀二次压缩）。并发默认 2（NGDC 防封控，吸收自
  MMPV-RNA public_metadata_pipeline/gsa_sra.down.py）。
- **安全**（Mimosa 多轮拦截后的最终形态）：出站 URL 逐次校验（仅 http/https、
  主机白名单、解析 IP 非私有/保留、重定向逐跳校验）；落盘路径全部经
  check_path；**写入一律 safe_open**（代码库强制约定：raw open 变量路径
  写入会被扫描器拦截；读取可原生 open）。
- **完整性**：ENA 字节数比对 + gzip 魔数/抽样解压（md5 串仅存档展示，
  不做弱哈希运算——Mimosa 禁 hashlib.md5）。
- **.sra→FASTQ**：fasterq-dump（平台根目录或 PATH）自动转换；未装则保留
  .sra 并在批次信息中提示。
- **独立进度窗口**：`/download?pop=1`（无导航精简视图，window.open 弹出）。
- **闭环**：批次完成后「🧪 进分析流程」→ 按 R1/R2 配对建样品（项目标签）→
  可选直接入批处理队列。已实测：ERR7586041（SARS-CoV-2，2×8.8MB）下载
  60s 完成、校验通过、队列自动跑完 13 阶段全 done。

### 项目管理与批处理队列
- **项目标签**：input.json 增 project 字段；/api/samples 返回；
  分析流程页样品列表显示 🏷 项目 + 下拉筛选。
- **批量导入**：分析流程页左栏"批量导入"卡，TSV（样品名/R1/R2）逐行建样，
  可选建后自动入队。
- **批处理队列**：SampleQueue 顺序执行（同一时间仅跑一个样品，防资源打满），
  状态持久化 tasks/queue.json；页面队列面板 4s 轮询；无 stages 时按
  "依次运行剩余步骤"语义自动补全；运行中条目不可删（先取消任务）。
- **避坑记录**：本模块开发中被 Mimosa 反复拦截 6 次——最终结论是三条硬规则：
  ① 出站请求必须"白名单 + 解析 IP + 重定向逐跳"三层校验；② 写文件一律
  safe_open（raw open 变量路径=拦截）；③ hashlib.md5 不可出现（完整性
  校验用"字节数 + gzip 抽样"替代）。后续写网络/文件类新模块请直接按此写。

## 八、.sra 转换引擎选型：sracha（2026-09-05 晚）

对比结论（针对 GSA .sra → FASTQ 的本地转换需求）：
- **xsra（ArcInstitute）不选**：CI 只测 ubuntu/macOS，从未支持 Windows；
  依赖 ncbi-vdb-sys（C 库捆绑编译），README 自述"构建系统相关、不可移植"。
- **sracha（rnabioco）选定**：纯 Rust、零 C 依赖、零 unsafe，单文件静态
  二进制。**2026-09-05 起改用官方 0.7.0 Windows 版**（bin/sracha.exe，
  CLI 与旧版兼容：fetch -O/-f/-q/--no-progress/--prefer-ena、
  fastq -O/-t/-f/-q、split-3 数字后缀、默认 gzip level 1 均保留；
  SRR28588231 实测 fetch 24s / fastq 0.2s，产物命名 `_1/_2.fastq.gz` 一致）。
- **历史：本地 0.6.0 Windows 补丁版**（已退役，备份于
  bin/sracha-0.6.0-local-patch.exe）：官方早期 release 只发 linux/macOS，
  曾自行在 Windows 编译（cargo 1.98 + MSVC）；唯一补丁是
  `crates/sracha-core/src/download/mod.rs` 的 Unix `write_all_at` →
  Windows `seek_write`（OVERLAPPED 定位写）。补丁源码曾存于
  thirdparty/git-repo/sracha-rs（该目录 2026-09-06 已清理删除，
  官方 0.7.0 后无需本地补丁）。
- **性能实测**（SRR28588231，23MiB，66220 spots/132440 reads）：
  sracha fastq → 0.26s 输出 `_1/_2.fastq.gz`（gzip level 1），
  命名与 fasterq-dump/ENA 惯例一致，直接对接 ready_files 配对。
- **平台接入**：`sra_convert_engine()` 优先 sracha.exe（bin/，探测已兼容旧根目录布局）→
  回退 fasterq-dump；下载批次收尾自动转换，并把 .sra 条目替换为产物
  FASTQ 条目（修复过 ready_files 看不到产物的 bug）。
- **附带修复**：cancel 落在解析阶段时 resolver 会把 cancelled 覆盖回
  downloading 并并发触发转换——加状态守卫（仅 resolving 才转 downloading）。
- sracha 还自带 `get/fetch`（直连 NCBI 下载 .sra，URL 模式
  sra-pub-run-odp.s3.amazonaws.com/sra/{acc}/{acc}），未来若要
  "NCBI 直下 .sra" 可用它；当前 ENA 直下 FASTQ 路径更快更省，暂不启用。

## 九、下载引擎三级回退（aria2c → 内置 HTTP → sracha）

sracha 本身是"下载+转换"一体工具：`fetch` 走 NCBI S3
（sra-pub-run-odp.s3.amazonaws.com/sra/{acc}/{acc}）+ 可 `--prefer-ena`
优先拿 ENA 预生成 FASTQ.GZ（无则自动回落 NCBI .sra，由收尾链转 FASTQ）。
自带断点续传与 MD5 校验。

- 适用范围：仅 INSDC run 编号（SRR/ERR/DRR）；GSA CRR 与任意 URL 不适用
  （NCBI 无此数据），这两类维持内置 HTTP 回退。
- 接入点：`_download_file` 中 aria2c/内置 HTTP 均失败后，若条目 acc 匹配
  `[SED]RR\d+` 且 sracha 存在 → `_dl_sracha`。产物识别：acc.sra →
  条目 out 指向 .sra（收尾转换链接管）；acc*.fastq.gz → 条目替换（同
  转换链逻辑）。verified=True（sracha 已做 MD5）。
- 实测：禁用 aria2c + 主 URL 404（白名单主机）→ 自动 sracha 回退，
  ENA 三个文件全部取回、配对正确。
- 引擎链全貌：下载 aria2c → 内置 HTTP → sracha fetch（仅 INSDC）；
  转换 sracha fastq → fasterq-dump。

## 十、BioAider 展示结构借鉴落地（2026-09-05 深夜，6 项全做）

数据流导向范式保持不变，落地 6 项展示增强（详见 BIOAIDER_REFERENCE.md 第五节）：

1. **专项分析目录侧栏**：`.tools-layout` 双栏 + 粘性 `.toc`（MENU 锚点跳转，
   scroll-margin-top 已设），组内条目覆盖 ①-⑥/序列查看器/最近运行。
2. **常用工具收藏**：工具卡标题 ⭐（localStorage `vp_fav_tools`）；
   总览页"⭐ 常用工具"条（favSection，无收藏时自动隐藏）。
   工具注册表 TOOL_REGISTRY 在 app.js（图标/中英名/href，跨页共享）。
3. **拖拽文件上传**：所有 `.filerow` 输入框接 drop 事件 → `/api/upload`
   流式落盘到 `uploads/`（重名加序号不覆盖）→ 自动填路径。
   注意：浏览器出于安全不给完整路径，所以方案是"拖拽=上传到平台再引用"。
4. **工具卡内嵌运行状态**：`injectStageLogs` 扩展 TOOL_LABELS 匹配
   （工具·/Tool·/NCBI/GenBank/同属 前缀）→ 各工具卡内 `#toolrun-<key>`
   容器（状态徽标+进度条+日志尾+取消+结果预览）。
5. **深色主题**：CSS 全量令牌化（surface/input 令牌补齐）后，
   `html[data-theme=dark]` 覆盖调色板；设置页 浅色/深色 分段切换
   （localStorage `vp_theme`）；12 个模板 head 加防闪脚本
   （注意：Flask 非 debug 模板有缓存，改模板要重启服务才可见）。
6. **轻量序列查看器**：`/api/seqview` 流式统计（FASTA/.gz，序列数/总碱基/
   GC%/简并%，分页 300bp 预览，2 万条封顶）+ 专项分析"📄 序列查看器"卡
   （svPath + 分页器）。

实测截图：tools_sidebar_fav / home_favstrip / settings_dark / tools_dark /
seqview / toolrun_inline（.zcode/ui_shots/）。
最终回归：12 页面 + 6 API 全 200；selfcheck 全绿。

## 十一、RdRP 模块移除 + trimAl 开关（2026-09-05 深夜，用户确认）

- **⑤ Open-Virome（RdRP/palmdb）彻底删除**：STAGE_ORDER/名称表/分组/权重/
  估时/PIPELINE_STAGES/STAGE_DEPS/STAGE_OUTPUTS/STAGE_VIEW/汇总函数（中英）
  全部移除；viz 报告章节、main.py CLI（cmd_virome + 子命令 + selfcheck 清单）、
  app.py 工具④的 rdrp 佐证块、tools.html 复选框、Virus_Platform_Core/virome.py 模块文件均删。
  README 同步（功能表/CLI/目录树/整章）。
  注意：导航"公共病毒组"（Open-Virome 雷达网页）是独立数据挖掘入口，保留；
  历史样品的 09_virome 目录已于 2026-09-08 用户确认后彻底删除
  （ERR7586041 / _archive:E2E-VERIFY / REGRESS 三处，各 5 文件）。
- **trimAl 开关**：⑦卡片新增 do_trim 复选（默认开）。关闭后直接用
  aln.fasta 建树（跳过清剪，日志明确提示）；开启时保持原行为（automated1
  + 过度修剪回退）。全链打通：卡片 → collectParams → _analysis_kwargs
  （_chk 默认 True）→ run_analysis(do_trim) → build_phylo(do_trim)。
  设置页"分析默认参数"新增 do_trim（共 18 项），可改全局默认。

## 五、模块化架构（2026-09-05 深夜改造）

三层接入点全部声明式，新增功能/调整顺序不再改核心代码：

1. **导航栏组件化**：所有页面 `{% include "_nav.html" %}`（唯一来源），
   当前页高亮由 `request.path` 自动判定。改导航 = 改 `_nav.html` 一处
   + i18n.js 两行。
2. **工具箱注册表（app.py `TOOL_REGISTRY`）**：
   `{'tool': {'title', 'title_en', 'need_db', 'job'}}`，
   job 构造器 `_tool_job_<name>(ctx)` 返回 `job(log, prog, cancel)`；
   ctx 提供 `p / run_dir / threads / db_virus / req / opt`。
   新增工具 = 实现一个函数 + 登记一行 + 前端卡片。
3. **管道阶段注册表（Virus_Platform_Core/pipeline.py `STAGE_REGISTRY`）**：
   `{key: {'deps': [...], 'fn': _stage_xxx}}`，阶段函数签名 `(C)`，
   通过 `C`（_Ctx）读写跨阶段产物（`C.cur_r1 / C.host_stats / C.vs / C.asm`）。
   执行器 `_stage_execution_order()` 对已选阶段做**稳定拓扑排序**——
   改变阶段顺序/插入新阶段只改注册表；未选中的依赖不补跑，
   沿用断点衔接（阶段函数读上一轮产物摘要）。
   `STAGE_DEPS` 兼容导出给管道页卡片解锁。

约定：
- 阶段函数只经 `C` 传数据，禁止模块级可变状态；输入解析（如
  `_assembly_inputs` 的 病毒reads→去宿主→原始 回退链）放在消费阶段内。
- 工具 job 内新增产物请写入 `ctx.run_dir` 并放进返回 dict（前端
  「文件」标签自动列出）。

## 十二、宿主去除独立模块化（2026-09-06，用户要求）

- **导航不再跳管道**：侧栏「宿主去除与序列提取（kunpeng 宿主库分类）」原
  `href: '/pipeline'`，点击落在自动化分析流程页，用户要求改为独立模块。
  现在 NAV_GROUPS 指向新页 `/hostremoval`（`_PATH_TO_GROUP` 归入 sample 组，
  二级侧栏高亮沿用 request.path 精确匹配）。
- **独立页** `webapp/templates/host_removal.html`：输入类型（双端/单端）、
  R1/R2 选择、置信度、线程、可选自备宿主库目录（留空=平台库），
  库状态徽章复用 `loadDbs()` 的 `#hostStat`。
- **后端**：TOOL_REGISTRY 新增 `hostremoval`（`_tool_job_hostremoval`），
  复用 `Virus_Platform_Core.host_removal.remove_host`（run_dir 充当 sample_dir，
  产物在 `<run>/01_host_removal/`：kept_R1/R2.fastq.gz + stats.json +
  分类报告）。启动前校验宿主库就绪（`Virus_Platform_Core.kunpeng.db_ready`），未就绪 400。
  样品管道的 ① 宿主去除阶段不受影响，两处共存。
- **前端匹配**：app.js `TOOL_LABELS` 增加 `hostremoval: ['宿主去除与序列提取',
  'Host removal']` → 模块卡 `#toolrun-hostremoval` 内嵌状态/结果；
  注意任务名含「工具·」前缀才参与 TOOL_LABELS 匹配，与管道阶段
  `① 宿主去除`（STAGE_LABELS，进 `#log-host`）不冲突。
- **回归**：py_compile 通过；test client `/hostremoval` 200 且侧栏正确指向；
  合成 50 对 reads 端到端跑通（任务 done、kept=50/50、stats.json 正确，
  随机序列不命中宿主库属预期），测试产物已清理。

### 十二（续）：其余管道入口全部独立模块化（2026-09-06，用户确认）

NAV_GROUPS 中剩余 4 个 `href: '/pipeline'` 条目全部改为独立页，
`_PATH_TO_GROUP` 同步归组，侧栏不再有任何跳管道的入口（管道页仅从
顶部导航 / 总览进入）：

| 模块 | 页面 | tool key | 复用后端 |
|------|------|----------|----------|
| 样品创建 / 批量导入 | `/samples` | —（纯前端，复用 `/api/pipeline/create` `/api/samples` `/api/queue/add`） | — |
| 宿主预测（ICTV 级联） | `/hostpredict` | `hostpredict` | `Virus_Platform_Core.host_analysis.predict_hosts` |
| ORF 预测 / 功能注释 | `/orf` | `orf` | `Virus_Platform_Core.orf.predict_orfs` + `Virus_Platform_Core.orf_annot.run_orf_annotation` |
| 基因组图谱 | `/genome` | `genoplot` | `gbdraw_plot.run_genome_plots` / `dfv_plot.run_dfv_plots`（fasta_in/ann_in 原生 standalone 模式） |
| 引物设计 | `/primer` | `primer` | `Virus_Platform_Core.primer.design_primers` |

实现约定（与 hostremoval 相同）：
- **伪样品目录**：stage 函数都是 sample_dir 口径，tool job 在 run_dir 下
  摆出最小布局再调原函数——orf/hostpredict 造 `03_assembly/`
  （summary.json+viral_contigs.fasta+contigs.filtered.fasta）；
  引物 conserved 模式造 `05_phylo/G1/aln.fasta`+summary.json，plain 模式
  造 `03_assembly/viral_contigs.fasta`。
- **分类表归一**：hostpredict 输入自动识别列名；工具④ metabuli 风格
  （taxid/taxon）自动转管道③口径（kunpeng_flag/kunpeng_taxid）。
- **引擎回退**：genoplot auto=gbdraw 优先、缺则 DFV；显式指定缺引擎时报
  可读错误。
- **TOOL_LABELS 新增** hostpredict/orf/genoplot/primer（任务名带「工具·」
  前缀才参与匹配，与管道阶段 `⑥ ORF 预测`/`⑧ 引物设计` 等无冲突）。
- **坑**：`extract_viral_contigs` 除 viral_contigs.fasta 外还需要
  `contigs.filtered.fasta` 兜底（首次端到端即被拦下），伪目录两件都要写。
- **回归**：16 页全 200；四工具合成数据端到端全 done（ORF 2 contigs/14
  ORFs、引物 2 对、宿主预测③④两种格式、图谱 gbdraw 出图）；内联 JS
  node --check 通过；测试产物已清理。

## 十三、检索页修复 + 全站按钮防重复 + 任务硬停止（2026-09-06 凌晨，用户反馈）

- **「开始检索」无反应根因**：meta.html 的 doSearch 调用全局 `val()`，
  但该函数不存在（app.js/i18n.js 均未定义，页面也没本地副本）→
  ReferenceError，fetch 根本没发出。download.html / submit.html 同病。
  修复：app.js 增加全局 `val(id)`。
- **任务硬停止（真停）**：原先 TaskManager.cancel 只置 Event，无任何任务
  检查，外部工具子进程照跑。现在 Virus_Platform_Core/utils 新增线程级 `task_bind /
  task_check_cancel` 上下文，run_cmd / run_cmd_redirect 启动的子进程自动
  注册到任务 rec['procs']；TaskManager.cancel 置位事件 + taskkill /F /T
  杀进程树；_run 结束统一解绑。取消后的任务 status=cancelled、
  error=「任务已停止（用户取消）」。单元验证：run_cmd 起长 ping →
  cancel → 进程树被杀、无残留。
- **按钮防重复点击**：app.js 新增 lockBtn / withBtn / watchTaskBtn
  （锁定按钮直到后台任务终态）。全站接线：meta（开始检索/提取/转入下载，
  任务卡新增「⏹ 停止」）、tools 工作台（runTool 五卡 + ncbi/gb/比较，
  _watchTask 扩展 btn 参数）、build 三按钮、管道页创建/批量、五个独立
  模块页全部运行按钮。app.css 增加 .btn:disabled/.btn-busy 禁用态。
- **回归**：py_compile 通过；16 页 200 + 接线抽查通过；全部模板内联 JS
  node --check 通过（help/tools 的 {{ tojson }} Jinja 表达式误报除外）；
  硬停止单元测试通过。启动服务后生效。

## 十四、检索页任务卡不显示 + 日志不流动修复（2026-09-06 凌晨，浏览器实测）

- **任务卡从未显示的根因（历史遗留）**：meta.html refreshTasks 渲染用
  `mine.map(t => …)`，参数 `t`（任务对象）**遮蔽全局 i18n 函数 `t()`**，
  回调内 `t('c.st.running')` 抛 TypeError 被 `catch(e){}` 静默吞掉。
  浏览器实测复现 → 参数改名 `task`。该页任务卡自上线以来从未渲染过。
- **运行日志空白根因**：engine_cmd 以 `python search_engine.py` 起子进程，
  stdout 接管道默认**块缓冲**，print 全攒在缓冲区。修复：源码模式加
  `-u`；_meta_stream_task 子进程 env 注入 PYTHONUNBUFFERED=1（冻结分发
  也覆盖）；Popen 登记进任务进程表（task_register_proc，取消即硬停止）。
- **停止按钮生效验证**（真实浏览器）：检索 7323 runs 任务 → 卡内
  「⏹ 停止」→ 立刻终止，日志尾 `[ERROR] 任务已停止（用户取消）`；
  状态映射补 `cancelled → 已停止`（原先误显示为失败）。
- 回归：py_compile 通过；live 服务重启后 /meta 200；任务卡 + 实时日志 +
  停止 全链路浏览器验证通过。

## 十五、全站任务日志审计：上屏 + 落盘（2026-09-06 凌晨，双代理并行审计）

**落盘（用户要求：日志保留到输出目录，重启/出错可查）**：
- TaskManager.start 新增 `log_file` 参数：任务每行日志（带 [HH:MM:SS] 前缀）
  同步追加写盘；缺省落 `logs/tasks/<任务名>_<时间>_<tid>.log`（文件名支持中文）。
- 显式路径：工具任务 → `tool_runs/<run>/run.log`；管道运行/样品队列 →
  `results/<样品>/logs/run_<ts>.log` / `queue_<ts>.log`；样品分析 →
  `results/<样品>/logs/analyze_<ts>.log`；公共检索 →
  `meta_search/<项目>/search/run.log`；元数据提取 → `.../info/run.log`。
  tasks/<id>.json 增记 log_file 路径。至此全部 15 类任务日志均落盘
  （此前只有 LOGAN 批量自带 run.log，其余重启即丢）。

**上屏（审计发现并修复）**：
- /api/tasks 日志尾 5 → 30 行（此前阶段卡/工具卡只能看到 5 行）。
- app.js：日志框运行中自动跟随尾部（autoscrollLogs，阶段卡/工具卡/任务卡
  统一生效）；refreshTasks 渲染体包 try/catch（单次渲染异常不再打断轮询）；
  断连横幅自愈探针（无 #tasks 的页面横幅也能自动消失）；
  STAGE_LABELS 补 `subsample`（子采样阶段卡此前永远无实时日志）。
- logan.html：任务名匹配正则补英文前缀（英文界面下批量任务完全不可见）。
- download.html：refresh 轮询链一次瞬时失败即永久停止 → 重排移入 finally。
- build.html：`</main>` 误入文件选择弹窗内部的标签嵌套修正。

**已知遗留（低优先，未改）**：多阶段管道任务名会同时命中链上多张阶段卡
（日志串卡）；下载批次无过程日志（aria2c 百分比被丢弃，仅 batch.json 状态）；
若干良性变量遮蔽（`const t =`）与硬编码中文待清理；settings 页"打开目录"
按钮全部只开平台根目录。

回归：py_compile + node --check + 16 页 200 + ORF 端到端（任务卡 14 行日志、
run.log 带时间戳落盘）+ 默认路径单测全部通过。

## 十六、三项遗留修复（2026-09-06 凌晨，用户点名 + 并行代理实施）

1. **多阶段管道日志串卡**：injectStageLogs 阶段卡注入从"任务名标签匹配"
   （多阶段任务名含链上每个阶段标签 → 整条链每张卡都塞同一任务）改为
   **按后端 prog 上报的当前阶段（task.stage）精确归属**——日志只进正在跑
   的那张卡，跑完自动挪到下一阶段卡；stage 未上报的启动瞬间退化为名称
   匹配且仅单阶段命中才注入；另按 curSample 过滤，其它样品的管道任务不再
   串到当前样品的卡上。附带修复：队列·样品 / 样品分析 任务此前不进任何
   阶段卡，现在同样按 stage 正确归属。
2. **下载批次过程日志**（并行代理实施，仅改 Virus_Platform_Core/public_data.py +
   download.html）：DownloadManager 新增 `_blog`（[HH:MM:SS] 前缀追加写
   `downloads/<bid>/batch.log`，OSError 静默）；创建/解析成功失败/开始下载/
   回退链失败原因（区分取消）/完成/转换/取消/重试 全程落盘；
   `list_snapshots()` 每批带 `log_tail`（40 行）；下载页每批卡片新增
   「📋 日志」展开按钮（展开状态跨轮询保持）。
3. **变量遮蔽清理**：app.js（watchTaskBtn snap、toolName/renderFavStrip
   tool）、tools.html（showRows ty、三处 watcher snap）、logan.html
   （bt/snap/task/prev 共 6 处）、results.html（tr）全部改名；app.js 中
   `let/const t =` 已清零。results.html 原有"良性"遮蔽（回调内 t 属性访问
   + 回调外 t() 混用）消除。

回归：py_compile（app/public_data/config/utils）+ 全部模板内联 JS
node --check + 16 页 200。服务重启后生效。

## 十七、全站结果表格分页（2026-09-06 凌晨，用户要求：每页 20 行 + 翻页）

- **共享分页器**：app.js 新增 `PER_PAGE=20` + `pagerHtml(page, pages, fnExpr)`
  （fnExpr 支持 `{p}` 占位，适配按批/按样品的翻页函数；省略号折叠长页码）。
- **meta 检索/元数据表（服务端分页）**：/api/meta/table 接受 page/per_page
  （默认 20，上限 200），返回 page/pages/total；废除原"仅显示前 500 行"
  截断（15851 行实测 → 793 页正确翻页）；Run 勾选跨页保留（SEL 按表
  身份重置）。
- **客户端分页 20 行**：工具④病毒 contigs 表（vcGo）、下载批次文件列表
  （dlFileGo，页码跨 2s 轮询保持）、LOGAN contigs 表（ctGo + ctSel Set
  勾选跨页保留，全选=全部页）、阶段文件预览 TSV/CSV（pvGo，废除 500 行
  截断）。
- **缓存版本号补齐**：meta/download/submit 三页 script 标签补 ?v={{asset_v}}
  （此前浏览器会拿旧缓存 app.js，本次及以往前端修复可能"不生效"即此因）。
- 序列查看器/引物/BLAST/CDD 等小表已自带分页或行数封顶，不重复处理。
- 回归：py_compile + 全模板 JS node --check + 16 页 200 + meta 分页端到端。

## 十八、导航调整 + 默认样品命名 + 跨模块文件交接（2026-09-06 凌晨，用户要求，双代理并行）

- **导航**：「公共数据下载」从数据资源组移入样本处理组（首位）——
  工作流顺序：拿数据 → 建样品 → 质控/去宿主 → 分析。_PATH_TO_GROUP 同步。
- **默认样品命名**：新增 `_default_sample_name()`，不填样品名时按文件名
  推导：去 .fastq(.gz) 等扩展名 + 循环剥技术后缀（_R1/_R2/_1/_2/_001/
  _L001，仅下划线数字避免误删 -01 真名）。/api/pipeline/create 与
  /api/analyze 两处旧 `split('_')[0]`（Lycium_barbarum 会变成 Lycium）
  统一替换。samples.html 补自动命名规则提示。
- **跨模块文件交接（并行代理实施）**：app.js 新增 sendTo/takePrefill
  （sessionStorage 暂存 → 目标页预填后销毁）；下载页每个完成文件行加
  「🧾 样品 / 🚫 去宿主 / 🧹 质控 / 🦠 鉴定」四按钮，按 acc 自动配对
  R1/R2（支持 _1/_2 与 _r1/_r2 形态）；hostremoval/samples/tools 三页
  支持预填。反斜杠路径统一转正斜杠避免 JS 转义破坏。
- **results 深链**：compare 组两条 /results 改 /results#msa、#sdt，
  结果页卡片补 id 锚点。
- **审计遗留（低优先待办，见报告）**：产物 files 页签加"发送到模块"
  按钮（toolrun files / 模块页最近运行）；t-contigs 独立成页；
  /api/analyze 孤儿端点去留；_safe_sample 双实现合并；管道阶段产物
  补下载链接。
- 回归：py_compile + 全模板 JS + 16 页 200 + 命名 7 用例 + 建样品端到端。

## 十九、文件浏览「此电脑」盘符修复（2026-09-06 上午，用户报告）

- 现象：文件浏览对话框停在「此电脑」盘符列表，点 C:\ / D:\ / E:\ 一律
  alert「无法打开目录」。
- 根因：`_joinBrowse(cwd, name)` 在 cwd='此电脑' 时仍执行 `cwd + name`，
  盘符被拼成 `此电脑D:\` → 服务端 400。盘符本身就是绝对路径，应直接返回。
- 修复：`_joinBrowse` 对 '此电脑' 层级直接返回子项名；盘根（尾反斜杠）
  直接拼接；其余 / 连接。node 单测 5 用例通过。

## 二十、contig 分类结果报告（对标 32-server）+ 下载失败排查（2026-09-06 凌晨）

- **文件下不下来**：实测 /tool_runs/<run>/<file> 路由 200 且内容完整——
  用户遇到的「无法从网站上提取文件」发生在服务重启窗口期（多次重启期间
  点了下载），重试即可，非代码问题。
- **contig 分类结果报告（对标 32-server 的 Metabuli 结果页）**：
  Virus_Platform_Core/viz.py 新增 `build_contig_report(run_dir)`：解析运行目录 kreport →
  分类谱系桑基图（root→界→…→种，按片段数）+ 分类旭日图 + 分类表
  （级别/分类单元/TaxID/%/片段数）+ contig 明细表（近完整绿色）+ 产物
  下载按钮；离线 plotly.min.js，固定浅色背景。写入 <run>/report.html。
- 端点 GET /api/tool/report?run=<name>：生成（每次重生成）并 302 到
  report.html；tools.html 工具④卡片新增「📊 分类报告」按钮。
- 浏览器实测：桑基/旭日 plotly 渲染正常、分类表 59 行、下载按钮齐全。

## 二十一、全站结果表格排序（2026-09-06 上午，用户要求）

点击表头切换升/降序（▲/▼ 指示，再次点击同一列反向）：
- **meta 检索/元数据表（服务端排序）**：/api/meta/table 接受 sort/dir；
  数值列全为数字时按数值排（pd.to_numeric），否则字典序；排序后再分页，
  翻页保持排序。前端表头 class=sortable，列下标定位（避免列名转义问题）。
- **客户端排序**：工具④病毒 contigs 表（contig/taxon/length/genus_avg_len/
  ratio/score）、LOGAN contigs 表（contig/length/species/family，勾选
  跨页跨排序保留）、阶段文件预览 TSV/CSV（任意列）。数字优先 parseFloat、
  否则 localeCompare('zh')。
- 样式：th.sortable（手型 + hover 提亮）。
- 未排（小表/编辑器）：提交网格（可编辑）、引物/CDD/BLAST 结果（≤20 行）、
  队列表。回归：py_compile + 全模板 JS + 16 页 200 + meta 排序端到端
  （asc/desc 首行不同且 asc 与 sorted 一致）。

## 二十二、工具页布局重排 + 结果中心整合（2026-09-06 凌晨，用户反馈，双代理并行）

- **工具卡布局重排**：toolRunHtml 从「结果/日志/文件」三标签切换改为纵向
  平铺：状态+进度条 → 📜 运行日志（常显、自动滚尾）→ 📦 结果 → 🗂 产物
  文件。refreshTasks 容错：无 #tasks 的页面跳过面板渲染但通知/
  injectStageLogs 照常。
- **工具页去乱**：删除底部「🗂 最近专项运行」与「⏳ 任务」区块；
  空 #moduleTree 左栏移入 main 顶部改横向切换条（消除空白列）。
- **独立模块页（5 页）**：删除右列「任务」+「本模块最近运行」卡片，
  合并为单列全宽。
- **结果中心新增「🧪 专项分析运行」区块**：/api/tool/runs 的 30 条记录
  迁移展示，文件行可展开下载（/tool_runs/...?dl=1），contigs_ 运行带
  「📊 分类报告」按钮。
- **内嵌分类报告**：工具④卡片新增 #contigReport iframe 内嵌
  report.html（运行完成自动内嵌/切换运行自动内嵌/📊 按钮手动内嵌，
  去重防轮询闪烁）。**坑**：/tool_runs 全局 Content-Disposition:attachment
  会让 iframe 空白——改为仅 ?dl=1 时 attachment，默认 inline。
- **下载失败说明**：/tool_runs 路由实测 200 完整；用户遇到的
  「无法从网站上提取文件」集中在服务重启窗口期（下载进行中服务被重启），
  重试即可；现已加 attachment 头使浏览器下载行为更明确。
- 回归：py_compile + 全模板 JS + 16 页 200 + 浏览器实测（工具页单列
  无空白栏、卡内日志平铺、内嵌报告 6 个 plotly SVG + 347 行表格渲染）。

## 二十三、桑基/旭日对标 32-server 实现（2026-09-06 上午，用户提供本地管线源码）

对照 `C-host_classify/plant_virus_db_pipeline/9.metabuli`（metabuli_api.py +
metabuli_page.html，即 <DEMO-IP> 服务器同源代码）重写报告图表：
- **数据流一致**：kreport 行（name/rank/count/depth）内嵌为 KROWS JSON，
  页面 JS 渲染（32-server 同款：SANKEY_DATA + limitSpeciesPerGenus）。
- **过滤规则一致**：剔除 no-rank 层级与 unclassified 大灰流（可勾选显示）；
  「每属仅显示前 5 的种」开关（默认开，丢弃种的深层后代）。
- **渲染细节一致**：plotly 桑键横向 pad20/thickness25、重复边合并累加、
  固定配色表（root/Viruses/Riboviria/…）+ 默认蓝、链接淡蓝 0.15 透明度、
  字号 11；旭日图同数据源 depth 栈链。
- 旭日图用 plotly 替代服务器端 KronaTools（Windows 无 kronagraph），
  交互能力等价（点击下钻/缩放）。
- 坑位记录：占位符 __KROWS__ 在 head 末尾，需对 (head+tail) 整体替换。
- 回归：report 90KB、KROWS 内嵌、浏览器实测桑基 188 链接/旭日正常、
  无 unclassified 灰流。

## 二十四、报告树形分类表 + 宿主预测整合（2026-09-06 上午，用户要求）

- **分类表树形展示**：kreport 本身就是树（kunpeng 输出）——分类表按
  缩进深度渲染层级，有子节点的行显示 ▾/▸ 可点击折叠/展开全部后代。
- **宿主预测整合进报告**：若运行目录存在
  08_host_analysis/host_prediction.tsv，contig 明细表自动在 taxon 后
  插入「宿主(预测)」列（final_host + 置信度），并在桑基图前绘制
  「宿主预测统计」柱状图（按 final_host 类别计数）；无宿主数据时显示
  引导提示。
- **原地宿主预测端点**：POST /api/tool/hostpredict_run {run}——对既有
  contigs 运行目录原地做 ICTV 宿主预测（含工具④ metabuli 风格 TSV →
  ③ 口径归一），日志落 <run>/hostpredict.log。tools.html 工具④卡片新增
  「🧲 宿主预测」按钮（运行完成自动刷新表格与报告）。
- **viral_contigs API 合并 host 列**：/api/tool/viral_contigs 若存在
  host_prediction.tsv 自动并入 host（final_host (置信度)），
  工具④页面表格新增「宿主(预测)」可排序列。
- 实测：148 条 contigs 全部判定（如 TRINITY_DN71 → Plant (Medium)），
  报告含树形分类表 + 宿主列 + 宿主统计柱状图。

## 二十五、分类报告 1:1 对齐 <DEMO-IP>/metabuli（2026-09-06 上午，用户明确要求"别自己发挥"）

- 已核对：线上页面与本地 `C-host_classify/plant_virus_db_pipeline/9.metabuli/metabuli_page.html`
  完全一致（空白差异除外），全部照本地源码逐行移植。
- build_contig_report v4 全面重写，报告结构/交互与 32-server 结果页一致：
  1. Taxonomy Sankey（"Only top 5 species per genus" 开关；过滤 no-rank
     保留 root/unclassified；重复边合并累加；固定配色 CM + 淡蓝链接）
  2. Classification Table（Rank 彩色 no-rank 灰/std 深蓝、Taxon 按深度
     缩进 8px 且 std 加粗、TaxID 灰、% / Reads；440px 滚动盒）
  3. Virus Sequence Classification（"Total N contigs, M classified as
     Virus"；列 Contig + Realm→Species + Length + Score + E-value +
     Genus Avg Len + Analyze；近完整(0.8≤len/avg≤1.2)绿底/长度绿粗；
     genus/species 加粗；sortVirusTable 点击表头排序默认降序；
     Analyze 按钮 CDD/BLASTN/BLASTX/Primer 经 parent.runAnalyze 调
     工具④四件套）
  4. Krona Taxonomy Plot（"Show unclassified" 开关；Windows 无
     KronaTools，用 plotly 旭日图等价替代，位置/交互一致）
  5. Downloads 底部按钮组（各产物 ?dl=1 强制下载）
- 保留平台自有增强：宿主预测 host 列 + 宿主统计柱状图（有
  08_host_analysis 时）；无则提示引导。
- 工具②（identify_*，测序 reads）报告：Classified Sequences 表
  （Seq ID/TaxID/Score）+ 相同桑基/旭日/分类表；tools.html 工具②卡片
  新增「📊 分类报告」按钮 + identRun 输入 + 内嵌 iframe + 完成自动打开。
- 实测：contigs 运行（148 contigs/199 谱系节点）与 identify 运行（真实
  病毒参考序列）报告均正常。

补充（同日）：分类报告已从 iframe 嵌套改为**页面内联渲染**——
/api/tool/report_data 返回 JSON，tools.html 用自带 plotly.min.js 页内绘制
桑基/分类表/宿主统计；工具④的「🔬 Contig 深度分析四件套」表格与 contig
明细合并（每行 CDD/BLASTN/BLASTX/Primer 按钮 + 宿主列 + 表头排序 +
每页 20 行分页）；iframe 全部移除；宿主预测完成后表格与报告自动刷新。
工具②（测序 reads）同样支持（Classified Sequences 表）。

## 二十六、工具④报告区合并重构（2026-09-06 上午，用户反馈）

- 移除重复的「🔬 Contig 深度分析四件套」独立标题区；工具④卡片统一为
  一条流：运行分类 → 输出面板 → 📊 分类报告（桑基 / 分类表 / contig
  明细+四件套按钮 / 宿主统计 / 旭日图 / 下载按钮组）。
- 宿主预测统计不再绘图，改为**表格**（宿主类别 / Contigs 数）。
- 补上内联报告缺失的**旭日图**（plotly sunburst，点击下钻）与
  **Downloads 下载按钮组**（各产物 ?dl=1）。
- report_data 接口增加 downloads 清单（按运行类型生成、校验文件存在）。

## 二十七、Krona 旭日图改用 taxburst（2026-09-06 上午，用户指定 github.com/taxburst/taxburst）

- `pip install taxburst`（离线 d3 内嵌，无 CDN 依赖）。
- Virus_Platform_Core/viz.py 新增 `_taxburst_nodes(krows)`（kreport 累计行 → name/rank/count/
  children 节点树）与 `build_taxburst(run_dir, krows)`：生成 taxburst.html
  （仅已分类）与 taxburst_unc.html（含未分类）两个离线交互 HTML。
- report.html 的 Krona 区块改为 taxburst iframe（Show unclassified 开关
  切换两个变体）；taxburst 不可用时回退 plotly 旭日图（renderSun）。
- 工具页内联（renderSunC）同样优先 taxburst iframe（report_data 附带
  taxburst 标记并负责生成文件），无文件时回退 plotly。
- 实测：浏览器内 taxburst 正常渲染（Viruses / Count: 148 / 交互控件齐全）。

## 二十八、组导航改顶部横向切换条（2026-09-06 上午，用户选择）

- 全部 11 个带组导航的页面：左侧栏移除，`_group_nav.html` 重写为
  顶部横向切换条（🏠 组名 + 各二级功能 pill，当前页高亮，悬停显示说明）。
- app.css：.group-layout 改 block、新增 .group-topbar 样式；
  .group-toc 旧规则弃用（保留不删）。
- 条目本身不变（用户确认）。
- 回归：16 页 200、全部页 topbar 在/旧侧栏无。

## 二十九、全量自查报告（2026-09-06 上午，用户休息前最终检查）

13 项全量检查全部通过：编译（6 文件）/ JS 语法（全部模板）/ 16 页 200 /
侧栏（顶部切换条全覆盖，旧侧栏归零，无管道跳转）/ 无嵌套 iframe /
report_data 双模式 / meta 排序+分页（793 页）/ 日志落盘（run.log +
hostpredict.log）/ 下载响应头（?dl=1 attachment）/ 工具②报告按钮/
宿主预测按钮/树形分类表/旭日图/下载容器/sendTo 跨模块/report.html
（taxburst iframe + 客户端渲染 + Analyze 按钮）。

## 三十、UI 精修（2026-09-06 上午，用户要求）

- **Toast 改屏幕中央**：CSS `#toastBox` 从 `right/bottom` 改为
  `left:50%; top:40%; transform:translate(-50%,-50%)`，动画从平移
  改缩放，更醒目。
- **线程数入各工具卡参数区**：移除工具页顶部全局线程独立区；每张
  工具卡（fastp/identify/assemble/contigs）新增各自的「线程数」输入框。
- **min_hit_groups 暴露**：工具②鉴定卡新增「最小命中组数」参数
  （kunpeng -g），不再硬编码 2。
- **E-value 列**：viral_contigs API 和 report_data 均从 contig_blast.tsv
  合并最优命中 E-value。
- 回归：py_compile + 全模板 JS + 16 页 200 + 全部标记。

## 三十、粘贴序列输入（2026-09-06 上午，用户要求）

- 后端 POST /api/paste_input：粘贴文本 → 写 uploads/paste_<ts>.<ext>，
  支持 .fasta/.fa/.fna/.fas/.fastq/.fq/.tsv/.txt；FASTA 校验 > 开头、
  FASTQ 校验 @ 开头；空内容 400。
- app.js 新增 pasteSeq(inputId, ext)：弹出粘贴对话框（modal），确认后
  POST 写文件并回填 input 路径。全站各序列输入框旁加 📋 按钮。
- 覆盖：工具④contigs FASTA、工具②鉴定输入、工具①质控 R1/R2、
  工具③组装 R1/R2、管道页样品创建 R1/R2、宿主去除 R1/R2、宿主预测
  TSV+FASTA、ORF FASTA、图谱引物 FASTA/GenBank。

## 三十一、通用病毒参考库（ref-virus / RVDB）（2026-09-06 下午）

> **2026-09-10 更新**：RVDB 库（databases/virus/rvdb/，1.6GB Kraken2 库）与
> 构建入口已移除——病毒验证/蛋白注释两条链路只依赖 nnot/prot + nnot/cdd
> + nnot/hmm/pfam，与 RVDB 无关；保留 RefSeq 通用库（databases/virus/ref/）。
> 以下为该功能的实现记录，留作历史参考。

用户在 databases/virus_ref/ 补充：NCBI RefSeq Viral
（viral.1.1.genomic.fna.gz，1.9 万条）、RVDB C-RVDBv32.1（132 万条，
头格式 acc|GENBANK|ACC|desc）、RVDB_Taxon_Current.tab.gz（1100 万行
accession→taxid+谱系，覆盖 RefSeq/RVDB accession）。

新增 Virus_Platform_Core/universal_ref.py：
- build_universal_meta(source)：流式扫 RVDB_Taxon 表 + FASTA 内 accession
  交集 → virus_ref/universal/<source>_acc2taxid.tsv（复用缓存）。
- build_universal_db(source)：注入 |kraken:taxid|N → kunpeng build_db →
  独立库 databases/refvirus_db（~1 分钟建完）或 rvdb_db（大库）。
- 端点 POST /api/build_universal_db {source}；/api/dbs 增加 refvirus/rvdb
  就绪徽章；数据库构建页新增「🌐 通用病毒参考库」卡片（两按钮 + 状态）。
- 用法：分析时把工具②/管道的病毒库路径指向 databases/refvirus_db 或
  rvdb_db 即切换到通用参考；与植物病毒库并存可按样品选择。
- 实测：RefSeq 库端到端建库成功（12,105/19,625 有 taxid 映射，其余为
  MAP/WT 未有序列条目属正常），库文件齐全、/api/dbs ready=True。
  RVDB 库按钮已就位（构建约 132 万条，用户按需启动）。

## 三十二、kunpeng 三种建库方式全支持 + Kraken2 库转换（2026-09-06 下午）

用户下载 Kraken2 官方预构建病毒库（databases/k2_viral_20260626.tar.gz，
含 hash.k2d/opts.k2d/taxo.k2d/nodes.dmp/names.dmp）。kunpeng 三种建库
方式现状与本平台支持：

- **A. 从基因组下载建库**（build --download-dir）：本平台等同能力为
  数据库构建页 Taxonomy + 宿主库/病毒库/通用库（走 add-library+build-db）。
- **B. add-library 自备 FASTA**：已有（build_virus_db / build_host_db /
  build_universal_db）。
- **C. Kraken2 库转换（hashshard）**：新增 Virus_Platform_Core/kunpeng.convert_kraken2()：
  支持 .tar.gz 包或已解包目录；解包（自动定位 hash.k2d）→
  `kunpeng hashshard --db <k2目录>`（**就地转换**，输出在 k2 目录内）→
  移动 hash_*.k2d/hash_config.k2d/opts.k2d/taxo.k2d 到目标库目录 +
  复制 nodes/names.dmp。端点 POST /api/convert_kraken2 {tar,name,
  hash_capacity}；数据库构建页「🌐 通用病毒参考库」卡片新增
  「🔄 从 Kraken2 库转换」按钮 + k2Stat 徽章 + /api/dbs k2viral 条目。
- **坑**：hashshard 不能重复跑（hash_config 已存在即 panic）——重试前
  需清理 k2 目录内的 hash_*.k2d/hash_config.k2d；输出位置是 k2 目录
  而非目标目录，需自行搬运。
- 实测：k2_viral tar.gz → k2viral_db 转换成功（/api/dbs ready=True）；
  用库内 RefSeq 真实序列分类冒烟：C t1 265522（Ichnoviriform fugitivi），
  kreport 层级完整（root→Viruses→…→Species）。

### 三十一（补）：refvirus 映射覆盖率 61.7% → 100%（2026-09-06 下午）

- **用户质疑 12,105/19,625 不合理 —— 属实**。根因：RVDB_Taxon_Current 表
  只含 RVDB 自身收录集，RefSeq 大量官方基因组（噬菌体/卫星等，7,520 条）
  不在其中。
- **修复（三源优先级）**：
  0. **NCBI 官方 nucl_gb.accession2taxid**（权威，覆盖全部 RefSeq/GenBank；
     aria2c 8 线程下载 2.7GB gz 到 virus_ref/universal/（80 分钟→约 12 分钟）；
     流式解析 4 亿行按 want 过滤）
  1. RVDB_Taxon（12,105）2. Kraken2 库包 seqid2taxid（18,567）兜底。
- 服务器 <DEMO-IP> 无 /home/USER（已核实；其 taxonomy 目录仅
  dmp，无 accession2taxid）。
- 实测：NCBI 官方源命中 19,625/19,625 = **100%**；库重建（24s，
  hash 0.71GB）；分类冒烟：原无映射的 NC_010393.1 Phage Gifsy-2 正确
  分类到种级 Essonnevirus Gifsy2（C t1 129862，kreport 全层级）。
- 注意：nucl_gb.accession2taxid.gz 为 NCBI 周更文件，DATA_VERSION 里
  可记录下载日期；过期只影响新增 accession，已有映射不受影响。

## 三十三、246 服务器三库构建 + 拉回本地（2026-09-06 下午，用户指定 246）

- **246**（zhangwenda@<HPC-IP>，256 核/629GB/12TB）：kunpeng 0.7.11
  已装（~/.cargo/bin/kun_peng）；源数据由本地上传（RVDB 919MB + RefSeq
  170MB + Taxon 表 102MB + Plant_Virus complete_ref + NCBI 映射用服务器
  现成的 /home/USER/database/taxonomy/nucl_gb.accession2taxid.gz，
  2026-03 版）。
- scripts/build_virus_dbs.py（纯标准库独立脚本）上传 246 构建三库：
  - plant_db：5,772 条，覆盖率 100%（NCBI 源），44s
  - refvirus_db：19,145/19,625 = 97.6%（3 月版 NCBI 映射缺 480 条新序列；
    上传本地 9 月版映射后复用重建 → 100%）
  - rvdb_db：1,321,456/1,321,608 = 100%，991s，hash 1.55GB
  - 坑：kunpeng build-db 需要库目录内 taxonomy/nodes.dmp+names.dmp
    （246 脚本从 /home/USER/database/taxonomy/ 复制）。
- 三库 tar 包拉回本地解压至 databases/{plant_db,refvirus_db,rvdb_db}，
  db_ready 全 True；分类冒烟三库各自参考序列全部种级命中
  （Janusivirus portis / Essonnevirus Gifsy2 / Salisharnavirus britensis）。
- 246 端三库保留在 /home/USER/kunpeng_dbs/（含 tar 包），可独立使用。
- 本地库全家福：virus_db（植物 55MB）+ plant_db（44MB）+ refvirus_db
  （713MB，100% 映射）+ rvdb_db（1.65GB，100% 映射）全部 ready。

## 三十四、外部审查修复（2026-09-06 下午，4 路并行审查 + 用户要求修复）

**已修复（高危 4 + 中危 6 + 低危 2）**：
1. 🔴 app.py /api/tool/viral_contigs NameError（run_dir 未定义，接口 500，
   contig 分类表全挂）——系 E-value 合并补丁引入，改回 _tool_runs_root 路径。
2. 🔴 Virus_Platform_Core/config.py tool() 放行 bundled: 哨兵（冻结分发 gbdraw 引擎恢复）。
3. 🔴 /api/seqview 允许平台外绝对路径（只读统计；与文件浏览对话框对齐）。
4. 🔴 app.js watchTaskBtn 轮询加 10 次不可达上限（服务重启后按钮不再永久
   卡「运行中」）。
5. 🟡 platform.json / queue.json 原子写（tmp + os.replace），断电不再丢配置。
6. 🟡 样品删除守卫改用 _task_matches_sample（前缀匹配对不上任务名的漏洞）。
7. 🟡 run_cmd timeout 看门狗线程（阻塞读 stdout 时也能超时杀进程）。
8. 🟡 工具运行 threads 钳制（1 ~ cpu*4）。
9. 🟡 esc() 补单引号转义；plotly.min.js 双加载去重（带 ?v=）；
   universal_ref NCBI 下载改 urlopen timeout=60 + 分块写。
10. 🟢 results.html checkbox data-i18n 误用改 span；暗色模式
    badge/sbadge/tstat 三组配色固定深底亮字；宿主去除 kept=0 早停；
    ORF 0 contigs 早停；selfcheck 失败退出码 1。

**核实无问题**：pyrodigal 3.7.1（v3 坐标，审查担心的 v2 移码不存在）；
STAGE_DEFS 重复已不存在；路径穿越防护/subprocess 参数/safe_open 均完好。

**遗留（择机，记录在案）**：orf_annot 无命中 ORF 被丢弃（n_orfs 语义，
需产品决策）；MAFFT/trimAl 中文路径 ASCII 中转（需实机 MAFFT 测试）；
API key 走命令行参数；assembly BLAST 库版本戳；i18n 硬编码中文若干；
kaleido 未使用。

## 三十五、审查遗留四项修复（2026-09-06 下午续）

1. **MAFFT/trimAl 中文路径**：Virus_Platform_Core/phylo.py 新增 _ascii_stage/_ascii_fetch
   （复用 assembly._ascii_work_base）；_run_mafft 输入非 ASCII 时复制到
   %TEMP%/vp_mafft 中转后运行；_run_trimal in/out 均非 ASCII 时中转并取回。
   实机验证：中文目录下 MAFFT 比对成功（mafft.bat 退出码 0）。
2. **API key**：/api/meta/info 的 DeepSeek key 优先级改为
   请求体 > 环境变量 DEEPSEEK_API_KEY > local 模式，密钥不再必须走命令行。
3. **BLAST 库版本戳**：ensure_virus_blast_db 记录参考 FASTA mtime+大小到
   virus.stamp，源更新后自动删除重建（此前一旦建好永不更新）。
4. **i18n 硬编码清理**：app.js 运行日志/结果/产物文件/位点/共识/页码/
   比全长/无变异/暂无数据等 10 处改 t()；i18n.js 补 c.runlog 等 9 键
   （中英各一份）。曾把英文键误插中文段——已修复（node --check 通过）。

## 三十六、注释组模块归位 + 保守域/同源复核独立成卡（2026-09-06，用户反馈）

用户指出：①「保守域 / 同源搜索复核」不应只是工具④的导航别名，应像
ORF/图谱/引物设计一样独立成模块；②「结构比较」放病毒注释组不合适，
应归比较基因组组。

- **NAV_GROUPS 调整**（app.py）：
  - annotate 组：`t-contigs-annot`（ref 别名 → 工具④）删除，新增真实
    条目 `t-annot`（独立卡片）；`t-contigs-struct`（结构比较）移出。
  - compare 组：`t-contigs-struct` 插到 `t-sdt` 之前（同为 identity
    矩阵/热图口径，SDT 卡紧随其后）。
- **tools.html**：新增 `<section id="t-annot">` 独立工作区——选 contigs_
  运行 → 病毒 contig 表（近完整绿底）→ 逐条 CDD/BLASTN/BLASTX/Primer；
  结构比较卡片原样移到 SDT 卡前。四件套 JS 改造：`runAnalyze/showCached/
  showAnaResult` 增加 scope 参数（`ANA_SCOPES.ana`=工具④内嵌表格 /
  `.ca`=独立卡片，共用 /api/tool/analyze 与结果缓存），结果 HTML 组装
  抽成 `anaResultHtml()` 两处共用；新增 `loadCaRuns/loadCaContigs/
  renderCaTable`（30s 轮询与工具④一致）。
- **参数钳制补齐**（后端兜底，此前只有前端 min/max）：orf min_aa
  1~5000、genoplot max_plots 1~200、primer num_return 1~20。
- **文档修正**：上文模块表「引物设计」页面误写 `/genome`，实为 `/primer`。
- **回归**：py_compile 通过；test client 全页 200；tools 页内联 JS
  node --check 通过；旧别名 `#t-contigs-annot` 链接失效属预期（落地卡片
  兜底），工具④卡内四件套按钮不受影响（scope 默认 ana）。

## 三十七、246 三库清理重建版拉回本地 + ②卡复选框归位（2026-09-07 凌晨）

- 246 端三库清理重建（零残留确认）后重新打包拉回（plant_db 107MB /
  refvirus_db 752MB / rvdb_db 3.0GB，构建时间戳 09-07 00:27-30）。
- 本地替换：`databases/{plant_db,refvirus_db,rvdb_db}` 三库 db_ready 全
  True（hash: 42MB / 676MB / 1.5GB）；旧 rvdb_db 已删除。
- 分类冒烟（Virus_Platform_Core.kunpeng.classify，NC_116488.1 前 20kb 切 150bp reads×15）：
  三库全部 C=15 种级命中——plant_db→Janusivirus portis、refvirus_db、
  rvdb_db 符合 §33 判据。注意 plant_db 是「植物病毒参考库」（4,286
  taxid 中 3,806 为植物病毒，源 plant_tagged.fa=Plant_Virus complete_ref），
  不是宿主库，勿用于宿主去除。
- tar 包移回 /tmp/dbpull/（246 端 /home/USER/kunpeng_dbs/ 另有全套）。
- 工具②卡内嵌报告复选框归位：`skTop5I` 移到「Taxonomy Sankey」标题下、
  `skUncI` 移到「Classification Table」标题下（此前堆在报告区顶部）。
- 服务端已验证：curl /tools?g=virus 返回新版布局；/api/dbs 六项全 ready。

## 三十八、公共数据检索页三问题修复（2026-09-07 凌晨，用户反馈）

**现象**：Core14/Full 元数据"都显示 1 个"；检索结果点击没反应；绘图没有
分开的图。**根因与修复**：
1. **陈旧元数据误导**：Lycium_chinense 的 Core14/Full 是 09-05 旧检索
   （1 条 Run）的产物，09-07 01:02 重新检索得 69 条后未作废 → 界面照常
   显示"1 个"。修复：`/api/meta/collections` 返回 `n_core14/n_full/n_plots/
   stale`（快算 CSV 数据行数，stale=元数据条数≠检索条数），卡片显示条数
   徽章 + 「⚠ 请重新提取」黄色警示。
2. **点击无反应**：showTable 无 catch，结果文件被新一轮检索重建期间
   404 → 静默失败。修复：所有请求失败在卡片内显示红色错误条 + 「刷新
   项目列表」按钮；空项目（如 Nicotiana_tabacum，search 目录为空）显示
   「无产物」提示。
3. **绘图拆分**：landscape_plot.py 面板重构为 `_pnl_a~f` 六个独立绘制
   函数（合并图与单图共用），`plot_sci_landscape` 现输出 3×2 合并总览图
   + A~F 六张独立图（各 PNG 600dpi + PDF 矢量，Panel_A_temporal /
   B_database / C_organizations / D_tissues / E_regions / F_stages）。
   Lycium_barbarum 已实跑验证（14 个文件）。
4. **布局重构**：meta.html ② 检索项目由"扁平行 + 全局结果卡"改为每项目
   一张独立卡片——头部统计徽章、操作按钮、内嵌结果表（排序/分页/勾选/
   下载，✕ 收起）、内嵌图库（fig-grid 网格缩略 + PDF 下载行）；排序状态
   按 `name:which` 独立存储；全局 plotPanel/tblCard 移除。
- 回归：py_compile / node --check 通过；/api/meta/collections 新字段、
  /meta 新布局、/api/meta/table 均实测正常；服务已重启。

## 三十九、注释组四卡拆分 + 全模块内置示例（2026-09-07 凌晨，用户反馈）

- **拆分**：「ORF 预测 / 功能注释」→「ORF 预测」(/orf) +「功能注释」
  (/annotation 新页)；「保守域 / 同源搜索复核」(t-annot) → 「保守域搜索
  (CDD)」(t-cdd) + 「同源搜索复核 (BLAST)」(t-hom)。注释组现为六条独立
  入口：ORF 预测 / 功能注释 / CDD / BLAST / 基因组图谱 / 引物设计。
- **后端**：TOOL_REGISTRY 新增 `orfa`（_tool_job_orfa）——参数 run=已有
  orf_ 运行（校验 `orf_` 前缀 + 目录存在）直接 run_orf_annotation；
  参数 fasta=输入 FASTA 时复用 _tool_job_orf（annotate=True）一步预测+注释。
  'orf' 标题改为「ORF 预测」；新增 /annotation 路由；
  NAV_GROUPS + _PATH_TO_GROUP 同步拆分。
- **annotation.html**：方式 A=选 orf_ 运行注释（列表按 files 含
  03_assembly 过滤）；方式 B=FASTA 预测+注释；均带 ✨示例。
- **t-cdd / t-hom**：ANA_SCOPES 扩展 cdd/hom 两作用域（cdd 卡只出 CDD
  按钮、hom 卡只出 BLASTN/BLASTX 按钮），loadCaRuns/loadCaContigs/
  renderCaTable 全部 scope 参数化（表格 id 由 run select id 推导），
  caExample(scope) 自动选最新 contigs_ 运行载入表。
- **内置示例**：databases/examples/example_viral_contigs.fasta（NC_077216.1
  6096bp + NC_027131.1 6988bp 两条植物病毒完整基因组）；app.js 新增
  EXAMPLE_FASTA + fillExample(id)；orf / annotation / genome / primer
  四页输入行均有「✨ 示例」按钮（t-cdd/t-hom 为运行选择型示例）。
- 回归：py_compile、tools/annotation/orf 三页内联 JS node --check 通过；
  /orf、/annotation、/tools?g=annotate 实测新布局；orfa 空参与非法 run
  名均 400 带明确错误；服务已重启。

## 四十、SDT exe 移除 → 纯 Python 精确复刻（2026-09-07 凌晨，移植 MMPV-RNA）

- **来源**：D:\桌面\延伸基因组\MMPV-RNA（sdt_genus_matrix.py +
  virus_auto_pipeline.py 的 build_mat_sdt_exact / plot_heatmap /
  plot_distribution / get_safe_leaf_order），SDT v1.3（Brejnev & Muhire
  2014）口径完全复刻：**每对序列 MAFFT 独立全局比对（--localpair）→
  Get_Similarity 公式**（分母只计两序列均有碱基的列），而非单次 MSA 近似。
- **Virus_Platform_Core/sdt_exact.py（新）**：_sdt_pair_worker（临时 FASTA → mafft.bat
  --quiet --localpair → 解析 → 原公式）经 ProcessPoolExecutor 并行
  （initializer 传序列全局，避免逐对 pickle）；.resume_cache 原子刷盘
  断点续传（缓存键 = n + 全序列 sha1 前 12 位）；21-mer 自动定向
  （与最长 3 条参考比 k-mer 共享，rev> fwd×1.1 翻正，revcomp 手写不依赖
  Biopython）；cluster_order（scipy average linkage，NaN 按 100 距离）；
  plot_heatmap（SDT 经典红黄 10 档 ListedColormap + BoundaryNorm 硬分界、
  自适应下限 95/90/85/80/十位、三角模式、手绘网格、cividis 色盲安全可选，
  PNG 300dpi + PDF fonttype42）；plot_distribution（SDT 原版红色分布曲线）；
  run_sdt_exact 主入口（≤max_n 截断、矩阵 CSV/JSON、六产物）。
  注意 matplotlib set_xticklabels 不接受 pad（用 tick_params(pad=)）。
- **接线**：TOOL_REGISTRY 新增 'sdt'（_tool_job_sdt，参数 seqs/max_n/
  orient/palette，threads=并行进程数）；t-sdt 卡片重写——「▶ 运行 SDT
  精确分析」跑新工具，结果区嵌入热图/分布图 PNG + 矩阵表 + 6 个下载，
  保留「查看已有矩阵」plotly 视图。
- **SDT exe 移除**：config.py 'sdt' 探测项、/api/open_sdt 端点、
  openSDT() JS、results.html「启动 SDT v1.3」按钮、i18n rs.sdtBtn 全部
  删除；results/tools 相关提示文本改为在线查看口径。
  identityMatrixTable 兼容 null（SDT JSON 对角与缺失为 null）。
- **回归**：6 条病毒基因组长序列实测（15 对，近缘 90.85% / 远缘 57.84%，
  聚类排序正常，热图/分布图质量目检通过）；续传缓存命中生效；API 端到端
  跑通 example FASTA（sdt_20260907_024615，产物齐全）；py_compile +
  app.js/tools 页 JS node --check 通过；服务已重启。

## 四十一、SDT 已比对直算 + 多色阶 + NT+AA 同一性表（2026-09-07，用户需求）

- **已比对输入**：SDT 卡新增「输入已比对 MSA」勾选——跳过 MAFFT 与自动
  定向，每对直接从 MSA 两行按 Get_Similarity 公式计算（_msa_pair_worker，
  同等并行/续传）；长度不一致时报错提示。API 参数 aligned。
- **热图色阶可选**：_cmap_colors(name) 把任意 matplotlib colormap 采样成
  10 档硬分界色阶；支持 sdt / cividis / viridis / RdYlBu / Spectral /
  YlGnBu / coolwarm / magma，SDT 与同一性表两卡共用下拉。
- **NT+AA 同一性表**（t-identity 新卡，tool 'identity'，BioAider
  Sequence Identity Matrix 口径）：run_identity_table——
  - 输入：核苷酸 FASTA（必选）+ 氨基酸 FASTA（可选，按名称匹配，缺失条
    自动 6-frame 最长 ORF 翻译兜底，Virus_Platform_Core/contig_annot.longest_orf_protein）；
  - NT / AA 两套矩阵（同一 aligned 参数口径，AA 长度不齐自动退回逐对
    MAFFT）；产物：identity_table.csv（逐对长表，NT 降序）、nt/aa_matrix.csv、
    identity_composite.png/pdf（**复合热图：上三角 NT / 下三角 AA**）、
    nt_heatmap.png、identity.json（前端逐对表直读）；
  - 前端：复合热图预览 + 逐对同一性表（# / 序列 A / B / NT% / AA%）+ 下载区。
- **缓存污染修复**：build_mat_sdt_exact——比对器缺失/不可达时快速失败；
  全部新算对均 NaN 时抛错且不写续传缓存（避免一次 mafft 故障把 NaN 固化，
  之后永远复用）。
- **回归**：4 条基因组建 MSA 后 aligned 模式 5.5s 出 NT+AA 全表（6/6 对
  有 AA 值，近缘 90.92/97.81）；未比对模式 10.1s；复合热图 viridis 目检
  通过；API 端到端 identity 任务完成；服务已重启。
- **合并（用户反馈）**：NT+AA 本就是 SDT 同源结果，独立 t-identity 卡撤销
  并回 t-sdt——「模式」下拉切换 NT（SDT 矩阵/热图/分布图）与 NT+AA
  （同一性表/复合热图），AA FASTA 输入行仅在 NT+AA 模式使用；前端
  sdtRunMatrix 按模式派发 tool 'sdt' / 'identity'，结果区共用一张卡
  （#sdExact / #idtResult 二选一显示）。NAV_GROUPS 恢复单项；工具
  'identity' 后端保留不变；API 端到端复测通过。
- **双轨补齐（用户指出 MMPV 原生支持核酸+蛋白）**：detect_seqtype——
  统计蛋白特征残基 E/F/I/L/P/Q 占比（核酸歧义码中不存在，蛋白序列通常
  >20%），>3% 判蛋白；run_sdt_exact 增加 seqtype='auto'|'nt'|'aa'——
  AA 模式自动跳过 21-mer 定向，标题/JSON 标注 AA Identity；SDT 卡模式
  下拉扩为三档 NT / AA / NT+AA（前端 seqtype 透传，结果 meta 显示
  序列类型与直算标记）；run_identity_table 增加蛋白输入守卫（NT+AA
  模式①误投蛋白序列时报错并提示切 AA 模式）。AA 实跑：2 条 ORF 蛋白
  3.2s 完成（detect=aa）；API 端到端 NT 自动判别正常。

## 四十二、全模块示例补齐 + 两轮集成实测揪出五个真 bug（2026-09-07，用户通宵验证）

### 修复的高危 bug

1. 🔴 **pyrodigal 坐标口径错位（Virus_Platform_Core/orf.py）**——`g.begin/g.end` 是 1-based
   closed（v2/v3 文档一致，本次实测 3.7.1 实证），代码却按 0-based 半开
   区间切片 `seq[begin:end]`：正链基因整体移码，faa 全是含几十个内部 `*`
   的垃圾序列，DIAMOND 命中为 0 → ⑥b 功能注释空表。负链基因因
   revcomp 方向巧合不受影响，历史 E2E（全是负链）因此没暴露——
   第三十四节外部审查"核实无问题"结论系误判，已在实机推翻。
   修复：`begin0 = begin-1` 切片；GFF 写回 `begin..end`（原 `begin+1` 也
   偏 1）。实测 2 条示例 contig：修复前注释 0 行，修复后 7/7 ORF 全部
   命中产物+类别。
2. 🔴 **gbdraw 纯 FASTA 模式坏死 + 假成功**——gbdraw CLI 要求 `--fasta`
   必须配 `--gff`（空 GFF 又报 No valid records），裸骨架路径逐 contig
   失败后任务仍报"完成 0 张图"。修复：`_write_skeleton_gff` 合成跨全长
   `region` 特征的骨架 GFF（实测 gbdraw 接受）；`_tool_job_genoplot`
   auto 引擎 0 张图时自动回退 DFV，仍为 0 则报错。
3. 🔴 **比对定向缺失（Virus_Platform_Core/phylo.py `_run_mafft`）**——de novo 组装常出
   反向互补 contig，样品流程比对/树/SDT 矩阵不翻正：IT-E2E 的自身
   contig 对自己的参考只有 36.4% identity、⑧保守区引物 0 对。
   平台 MAFFT 是 v6.864b（2011）无 `--adjustdirection`；新增
   `_orient_normalize`（比对前 21-mer 定向翻正，复用 sdt_exact）。
   修复后 identity 100.0%、保守区引物 3 对。
4. 🔴 **`auto_orient` 自锚失效（Virus_Platform_Core/sdt_exact.py）**——旧实现把全部序列
   k-mer 并集当锚：反向序列自身贡献的 k-mer 恒在锚集内（fwd≈rev），
   永远判不出反向——SDT 卡的「21-mer 自动定向」实际从未生效。
   改为固定全局锚（ref_seqs 最长者，否则 recs 内最长者）。
5. 🟡 **`/api/tool/runs` 500**——os.walk 遍历遇 MMseqs2 临时库特殊文件
   `*.mmseqs_tmp/latest`（WinError 1920 不可 stat）直接炸接口；
   getsize 加 OSError 兜底跳过。顺手 `_run_search`（mmseqs 分支）补
   `shutil.rmtree(tmp)` 清残留临时目录（存量 10 处已清理）。
   另：该接口原按目录名倒序取 30 条，字母序挤掉真正新运行
   （structcmp > genoplot），改按 mtime 排序。
6. 🟡 **冻结分发 multiprocessing**——SDT 精确矩阵/同一性表用
   ProcessPoolExecutor；打包 exe 子进程会重启 exe，app.py 顶部补
   `multiprocessing.freeze_support()`（源码模式空操作）。

### 测试与示例

- **tests/_it_annotate.py**（新增）：注释组全链路集成——orf → orfa
  （run 模式写回 orf_ 目录 + FASTA 一步法）→ genoplot（FASTA 骨架图 +
  GenBank 输入）→ primer（plain + conserved）→ contigs 分类 →
  viral_contigs/report_data/analyze(primer, blastn 在线) → runs 列表。
  实测 BLASTN 在线命中 20 条。修复前该脚本抓出 bug 1/2/5。
- **tests/_it_compare.py**（新增）：比较组全链路集成——structcmp（近缘
  98.5% / 跨参考 34.4% 校验）→ msa_data → SDT 精确矩阵（15 对 +
  对角线 100 + 近缘>跨参考）→ NT+AA 同一性表 → quicktree（NJ/FastTree）
  → 树文件 API + 路径穿越拒绝 → Entrez 在线检索（1030 条）→ 示例
  GenBank 集合导入 → 共线性比较（html/svg/png/clusters/similarity/
  faa/m8）→ 集合建树 → gb: 伪样品 MSA/树/SDT API。**测试体必须收进
  main()+__name__ 保护**——sdt/identity 的 ProcessPoolExecutor 在
  Windows spawn 下会重导入 __main__，无保护递归崩溃（早期运行即踩）。
- **tests/_it_platform.py**（新增）：21 个页面全 200 + 渲染含各示例按钮
  + 示例文件齐全 + 15 个关键 API + 样品报告可访问 + 静态资源。
- **tests/_it_submit.py**：修复过时 API（create_table(demo=True) →
  sample='demo'；demo 数据现为 6 行），改受控 3 行表续跑。
- **databases/examples/** 新增（tests/make_examples.py 确定性生成）：
  example_virus_set.fasta（3 条真实 tobamovirus 参考 + 3 条 2%/8%/15%
  受控突变衍生，MSA/结构比较/SDT/快速建树示例）、
  example_conserved_set.fasta（同种近缘 0.1%~1.2%，保守区引物示例——
  保守区设计要求全表 ≥90% 一致连续列 ≥400，跨物种集达不到）、
  example_tree.nwk（上集 FastTree 树）、example_genome.gb（pyrodigal
  坐标 + 占位 product，基因组图谱 GenBank 模式）、
  example_synteny_A/B/C.gb（3 条同属小基因组各 6 CDS，共线性离线示例）。
- **示例按钮补位**（tools.html / genome.html / primer.html + app.js）：
  t-msa、t-tree（快速建树 + 本机树文件）、t-contigs-struct、t-sdt 指向
  近缘集、t-ncbi（ncbiExample 预填检索式）、t-synteny（syntenyExample
  预填内置 3 条 .gb 离线导入）、genome 页「✨ 示例 .gb」、primer 页按
  模式切换示例文件（app.js 新增 EXAMPLE_SET_FASTA /
  EXAMPLE_CONSERVED_FASTA / EXAMPLE_TREE_NWK / EXAMPLE_GENBANK_GB /
  EXAMPLE_SYNTENY_GBS 常量 + ncbiExample/syntenyExample）。

### 全平台验证结论（模拟 + 真实数据）

- **主流程 E2E**：tests/make_synthetic.py 重生成 PMMoV 8000 对高 Q reads
  （旧数据 Q 全 ≥60 触发 SPAdes "Failed to determine offset"，fake_qual
  下限改 55；fastp/SPAdes 对合成数据的边界又踩清一次）→ IT-E2E 全 12
  阶段（fastp→fq2fa→宿主→筛查→组装→宿主预测→ORF→注释→树→引物→图→
  报告）一次通过；8000/8000 reads 分类命中，宿主预测 Tombusviridae/
  Plant 0.998 High；报告 200KB 内嵌 SVG/plotly。
- **公共检索**：meta-search 真查 5 条 Run（SRA_GSA_Merged_Final.csv）。
- **下载**：/api/dl/create→cancel→delete 生命周期实测通过（真实大批量
  下载此前已有 downloads/ 证据）。
- **LOGAN**：logan-create（2 片段切分）→ import_result（ fabricated
  Logan-Search 表）→ trace_report.html 28KB 生成；真实外网提交不代测。
- **回归**：_it_msa / _it_synteny / _it_phylo / _it_submit / _smoke_phylo
  / _it_annotate / _it_compare / _it_platform 全部 PASSED。

## 四十三、比较基因组组整改：删「结构比较」、拆「GenBank 集合管理」（2026-09-07，用户反馈）

- **删除「结构比较」卡**（t-contigs-struct）：与 ⑨ SDT 卡功能重叠（同为
  两两 identity 矩阵/热图，仅"一次共享 MSA"与"逐对独立比对"之差），按用户
  要求整卡移除——NAV_GROUPS 条目、tools.html `<section>`、专用 JS
  （loadStructRuns/scUseRun/loadStructData/renderStruct/SC_DATA）、runTool
  structcmp 分支与初始化调用全部清除。**后端 TOOL_REGISTRY 的
  'structcmp' 工具与 /api/tool/structcmp_data、/api/tool/msa_data 保留**：
  ⑦ MSA 卡「运行比对并查看」仍以 tool=structcmp 出比对（msaRunAlign），
  集成测试也经 API 复用。
- **「GenBank 集合下载 / 导入 / 巡检」拆为独立卡片**（t-synteny-gb）：
  原 t-synteny-gb 是 `ref: t-synteny` 别名菜单（点击落在共线性卡，用户
  明确不要）——现改为落地卡片：集合名/检索式/accession/条数 + ⬇ 下载 +
  ✨ 示例（预填内置 3 条 .gb）+ 📥 本机 .gb 导入 + 输出区 +
  **#gbCollectionsMgmt 集合列表（🔍 巡检）** + #gbWarnings 警告明细。
  loadGbCollections 改为双渲染：同一份 /api/gb/collections 数据，
  共线性卡 #gbCollections 出 ▶ 比较 / 🌳 建树，管理卡 #gbCollectionsMgmt
  出 🔍 巡检。共线性卡（t-synteny）只留比较参数（一致性/覆盖度/着色/
  LoVis4u）+ 建树工具 + 集合列表；失去入口的 runCompare（读 s_name 的
  手动按钮）删除——集合名输入移去管理卡，比较统一走列表「▶ 比较」。
- **回归**：app.js node --check 通过；tools 页渲染后两段内联 JS
  node --check 通过；py_compile 通过；导航 JSON 确认 compare 组 =
  t-ncbi / t-synteny-gb / t-msa / t-tree / t-sdt / t-synteny 且页面
  无 t-contigs-struct 残留；_it_platform 断言同步更新（结构比较移除 +
  双卡存在 + 双列表渲染标记）后 PASSED。

## 四十四、比较基因组组按分析目的重组 + ICTV 科属选参 + 旧项目归档（2026-09-07，用户反馈）

用户定调四模块（按分析目的而非工具堆叠）：
**参考序列获取 → 进化树构建（科/属级）→ SDT 同一性（属级）→ 同属共线性（LoVis4u）**。

- **t-seqprep「参考序列获取」**（合并原 t-ncbi + t-synteny-gb，都是序列准备）：
  - 新增 <b>ICTV 科/属选参下载</b>：GET /api/ictv/taxa（taxa.txt 按宿主含
    plant 过滤，科/属下拉含 accession 数，实测 236 属）；POST /api/ictv/preview
    （select_refs 本地优先预览，不下载）；POST /api/ictv/download（后台任务：
    select_refs → download_gb_collection(accessions) → 同步 FASTA 参考集到
    databases/ncbi_refs/<coll>/refs.fa）。实测 Dianthovirus 属：6 条 GenBank +
    6 条 FASTA，直接喂「进化树构建」建树成功。
  - accession 列表 / Entrez 检索式下载（/api/gb/download 也补了 FASTA 同步）；
    本机 .gb 导入（✨ 示例 .gb 预填 3 条内置文件 + 集合名）；
    集合列表（🔍 巡检）+ FASTA 参考集列表。
  - 原 t-ncbi 的"仅 FASTA 集合"下载入口撤销——现在所有来源统一产
    GenBank 集合 + FASTA 双产物（ncbi_refs 集合名即样品流程可用的参考集合名）。
- **t-treebuild「进化树构建（科/属级）」**（合并原 t-msa + t-tree）：
  ① 集合建树（tbColl 下拉 + NJ/FastTree/RAxML-NG，/api/gb/phylo）；
  ② FASTA 直接建树（quicktree）；③ 树查看（本机 Newick / 样品与集合树）；
  ④ MSA 查看（SNP-only）。所有元素 id 原样保留（qt_fa/tv_*/ms_*），
  msaToolInit/treeToolInit 等 JS 零改动。
- **t-sdt**：标题改「SDT 同一性分析（属级）」。
- **t-synteny「同属共线性比较（LoVis4u）」**：/api/compare/run 的 lovis4u 默认
  True（前端不再有开关/着色风格选择）；内置引擎降级为 LoVis4u 不可用时的回退
  （clusters/similarity 表照出）。实测 Dianthovirus 集合 → lovis4u.pdf 51KB。
  集合列表按钮 = ▶ 比较 + 🌳 建树（tbBuildFor 跳建树卡并选中集合）。
- **旧项目归档**：E2E-VERIFY/REFACTOR-E2E/REGRESS/SYN-smoke/IT-E2E/dfv_demo/
  _it_* 等 9 个测试样品 → results/_archive/（保留 NX-6、ERR7586041、
  SRR39909438 三个真实样品）；gb_collections 的 _it_synteny/_example_synteny/
  tobamo_test → databases/gb_collections/_archive/。列表口径统一过滤下划线
  前缀（api_samples / api_msa_samples / api_gb_collections /
  api_ncbi_collections / page_results）；新增 /api/samples/archived +
  /archive_report/<sample>/，结果中心样品表下新增折叠「📁 归档项目」区
  （可打开报告，提示如何恢复）。meta_search 测试项目 Tobacco_mosaic_virus
  删除。
- **坑**：tm.start 的 prog 是 (stage, pct, msg) 三参，与
  download_gb_collection 内部 prog 约定一致——中间包一层 2 参 _prog 必炸；
  collection_records 返回 (acc, organism, seq) 三元组而非 SeqRecord。
- **回归**：_it_platform（含新模块断言 + 归档断言）/ _it_msa / _it_synteny /
  _it_submit / _smoke_phylo / _it_compare / _it_annotate 全部 PASSED；
  渲染内联 JS node --check 通过；ICTV 预览/下载/建树/LoVis4u 比较真实链路实测。

## 四十五、ICTV 级联选参 + 基因级建树 + 独立序列比对模块（2026-09-07，用户反馈）

用户三点需求：① ICTV 科属选择改为界门纲目科属种**递进级联**（并核实计数）；
② GenBank 集合支持提取 CDS/PEP 做比对建树；③ 比对+trimAl+查看编辑独立成模块
（参照本机 PhyloSuite——其为 PyInstaller+PyQt，内置 Lg_msa_view 彩色表格
查看器支持编辑保存，无外部 Jalview 依赖）。

- **计数核实**：taxa.txt 每行=1 个 accession（多 accession 已拆行），各级
  计数即 accession 条数：Geminiviridae 986（植物口径）/988（全部宿主）、
  Begomovirus 846、属 236、种 2478——与 select_refs 实取数一致。旧版预览
  按植物过滤而下载不过滤的不一致已统一为植物口径（plant_only）。
- **级联（Virus_Platform_Core/ictv_db）**：RANK_COLS 九级（含 Kingdom/Subfamily 列，主级
  Realm/Phylum/Class/Order/Family/Genus/Species）；rank_options(levels) 给
  已选过滤下每级的 [{name,n}]；select_refs 增加 ranks={列:值} 交集过滤 +
  plant_only。API：/api/ictv/cascade（GET，已选级作过滤返回全部级选项）、
  preview/download 改吃任意级组合（_ictv_levels 取最深层定位）。前端
  spCascade 七个下拉逐级联动（变更即清更深级并重取），示例预填
  Riboviria→Tombusviridae→Dianthovirus。
- **基因级建树（Virus_Platform_Core/gb_collection）**：_gene_feature_seqs 按 product/gene/
  note 关键词匹配 CDS（每记录取首个命中），PEP 优先 translation 限定符、
  缺则按 transl_table 翻译；build_collection_phylo 增加 molecule=
  genome|cds|pep + gene，基因级产物落 gene_trees/<slug>_<molecule>/；
  PEP 的 FastTree 用蛋白模式（-gtr 不适用）。/api/gb/phylo 与建树卡同步
  （分子下拉 + 关键词输入，genome 时禁用）。实测示例集合 coat protein：
  PEP/CDS 各 3 条提取建树通过；无关键词/无命中均有清晰报错。
- **独立序列比对模块（t-align，导航第 2 位）**：TOOL_REGISTRY 新增
  'align'——_run_mafft 增加 strategy=auto|linsi|fast（L-INS-i=--localpair
  --maxiterate 1000；fast=FFT-NS-1），trimAl automated1/gappyout/strict/
  none，产物 input/aln/aln.trim。查看器 API：/api/align/file（平台内或
  绝对路径只读，返回 names/seqs/aligned/type）与 /api/align/save（FASTA
  校验后写 <名>.edit.fasta 副本；平台外源落 tool_runs/align_edit_*/）。
  前端：彩色等宽查看器（核酸/蛋白双配色 + 分页 + 图例 + 长度不一致警示）、
  ✏️ 编辑模式（textarea FASTA 化）→ 💾 保存副本 → 自动重载，结果可
  一键「送进化树构建 / 送 SDT」。t-treebuild 的 MSA 查看保留（建树上下文）。
- **测试**：tests/_it_align.py（级联/比对/查看/编辑/基因树全链路，PASSED）；
  _it_platform 断言更新（五模块导航 + 级联/比对/分子标记）；_it_compare/
  _it_msa/_it_synteny 回归通过；渲染内联 JS node --check 通过。

## 四十六、参考序列获取卡：CDS/PEP 提取 + ICTV 抽样上限（2026-09-07，用户反馈）

- **CDS/PEP 提取（PhyloSuite 布局）**：集合管理列表新增「🧬 提取」按钮 →
  后台任务 extract_collection_features →
  gb_collections/<名>/extract/{genome.fa, CDS.fa, PEP.fa, genes.tsv,
  CDS/<基因>.fa, PEP/<基因>.fa}——分类分目录 + 按基因拆分（同源基因每
  基因组一份，序列 id=<acc>|<基因>；PEP 优先 translation 限定符）。
  产物面板（集合下拉 + 清单表）每项带 ⬇下载 / →送比对 / →送建树；
  /api/gb/collections 增加 has_extract 标记；/api/gb/extract(_files)
  端点。实测示例集合：3 基因组 × 6 CDS/PEP、6 个按基因文件、
  genes.tsv 7 行。
- **ICTV 每属/每种抽样上限**：ictv_db.cap_per_rank(rows, per_genus,
  per_species)（保持本地优先序封顶）；preview/download 接受
  per_genus(0-50)/per_species(0-20)，先封顶再取 limit（select_refs 内层
  limit 放大到 max(limit*20, 2000) 保证封顶前有足够候选）。UI 两个下拉
  （每属 不限/3/5/10，每种 不限/1/2/3），预览标注抽样口径与丢弃数。
  实测：Tombusviridae 101→64（每属≤5）→62（+每种≤1）；Potyviridae
  290→42（max/属=5）。
- **测试**：tests/_it_extract.py（抽样效果 + 提取布局 + 清单 API + 404）
  PASSED；_it_platform / _it_align / _it_compare 回归通过。

## 四十七、下载批次删除拆分 + 数据库构建页去任务面板（2026-09-07，用户反馈）

- **公共数据下载 · 三模式删除**：原单个 🗑（记录+文件一并删）拆为
  📂 删文件（清数据文件，保留 batch.json 记录与 batch.log 日志，文件
  条目标记 deleted 仍可在列表追溯）/ 🧾 删记录（删 batch.json+batch.log，
  数据文件保留在 downloads/<批次>/ 磁盘原处）/ 🗑 全删（原行为）。
  后端 DownloadManager.delete_files / delete_record（运行中批次先自动
  cancel；已完结批次不再被误标 cancelled）；delete/delete_record 对
  不存在批次改返回 False（API 400）。确认弹窗按模式区分文案。
  测试：tests/_it_dl_delete.py（三模式目录断言 + 日志落笔 + API 链路 +
  400 守卫）PASSED。
- **数据库构建页**：底部「任务」汇总卡移除——各构建卡下方 buildrun 区
  已实时输出日志，不再重复展示；refreshTasks 的空面板守卫保证连接
  横幅/toast/卡内日志注入不受影响。

## 四十八、公共数据检索 ↔ 下载联动补齐（2026-09-07，用户反馈）

- **Info 元数据表（Core14/Full，ib 视图）**：勾选栏原只有计数，补
  「⬇ 勾选 Run 转入数据下载」——intoDownload 泛化 which 参数
  （缺省 'search' 兼容原 S-浏览调用），runsOfRows 本就支持任意 which。
- **行内单 Run 快捷下载**：通用表格渲染（renderView，sb/ib/local 三视图
  共用）在 Run 列追加 ⬇ 锚点（正则校验 SRR/ERR/DRR/CRR 形态才显示），
  dlOneRun 单 Run 直接 POST /api/dl/create（批次名 <项目>_runs）并跳转
  /download。实测单 Run 建批 → 取消 → 全删链路 200。
- 回归：/meta 渲染含新按钮与 dlOneRun；两段内联 JS node --check 通过。

## 四十九、工具卡历史口径统一说明（2026-09-07，用户问询）

用户问：下载/样品/格式转换/质控/宿主去除页为何有的有旧记录有的没有。
根源与口径：
- **公共数据下载**：批次落盘 downloads/<id>/batch.json，服务重启也在 → 一直有记录（该页本身就是批次管理器）。
- **样品创建/批量导入**：列表 = results/ 目录；旧测试样品已按需求归档进
  结果中心「归档项目」，故只剩真实样品（NX-6 / ERR7586041 / SRR39909438）。
- **格式转换 / 序列质控 / 宿主去除（及识别 / 组装）**：卡片 .toolrun 输出区
  只渲染 TaskManager 内存中的本会话任务（injectStageLogs），服务重启即清空；
  历史一直都在 **结果中心 → 专项运行**（/api/tool/runs 按目录 mtime 取
  最近 30 条，含全部产物文件下载与打开目录）。
统一处理：五处输出区下追加固定提示行「🕘 本输出区仅显示本次会话任务；
历史运行与产物 → 结果中心 · 专项运行（重启后到此查看）」（tools.html
四卡 + host_removal.html 独立页），口径对齐"模块页干净、历史集中结果中心"。

## 五十、全模块历史运行区：持久保留 + 上下游衔接 + 折叠/删除（2026-09-07，用户需求）

用户定调：不只下载数据，**所有模块都保留运行记录**，且支持**上下游衔接
复用**；为防过多只要求**可折叠、可删除**。方案（取代第四十九节的"结果中心
指引"静态提示）：

- **数据源**：/api/tool/runs（tool_runs/ 目录持久，重启不丢），按每卡
  prefix 过滤（convert_/fastp_/hostremoval_/identify_/assemble_/contigs_/
  align_/structcmp_+quicktree_/sdt_+identity_/orf_/orfa_/genoplot_/primer_/
  hostpredict→contigs_）。
- **组件**（app.js 通用 RH_DEFS/rhRefresh/rhToggle/rhFiles/rhDel/rhGo）：
  卡输出区下方插入 🕘 历史运行（N）——默认折叠、计数徽标；条目 = 运行名 +
  文件数 + **衔接按钮组** + 📁 文件（展开逐文件下载）+ 🗑 删除
  （POST /api/tool/runs/<run>/delete，名称白名单 + 平台内限定 + 不存在 400）。
- **上下游衔接**：每工具定义 chains(files)，按产物文件名（fastp_R1/R2、
  kept_R1/R2、viral_sequences、contigs.filtered、viral_contigs、
  pyrodigal.faa/ffn、orf_annotation.gff3、aln(.trim)、tree/nj.nwk、
  virus_classification.tsv）生成"→ 下游模块"按钮；点击经 sessionStorage
  vp_runfill 暂存 {输入框id: 路径} → 跳目标页 → rhApplyStashed 回填
  （支持 input/select/checkbox，如 i_type=fasta、sd_aligned 勾选、
  oa_run 选运行）；同页目标立即填充不跳转。
- **宿主预测页**无自有运行目录，历史 = 有分类结果的 contigs_ 运行，
  一键「载入分类表」（hp_tsv + hp_fa）。
- 任务刚结束时（refreshTasks 钩子 __rhOnTasks）自动刷新各历史区。
- **回归**：_it_platform 新增容器数量与删除 API 断言 PASSED；全部页面
  内联 JS node --check 通过；删除 API 实测（建目录→删→400 守卫）。

## 五十一、④分类卡重排 + 宿主自动化 + ID 复制 + 四件套跳转与 metabuli 风格结果（2026-09-07，用户需求）

- **布局**：Krona Taxonomy Plot（旭日图）移到 Virus Sequence Classification
  表之前（Sankey → 分类表 → 宿主统计 → Krona → 分类明细表）。
- **宿主预测自动化**：新增「宿主预测自动运行」开关（默认开）——分类表加载
  后若无宿主列（autoHostIfMissing，按运行去重）→ 自动 POST
  /api/tool/hostpredict_run → 完成后重载表格，宿主列并入（实测
  'Plant (High)'）；宿主库未就绪时静默跳过（手动 🧲 按钮会提示）。
  表头新增**宿主筛选下拉**（vcHostF，从行数据去重生成，按 "类别 (置信)"
  前缀匹配），计数行显示筛选前后条数。
- **ID 复制序列**：分类表 Contig ID 单元格可点（⧉ 前缀 + 虚线下划线）→
  copyContigSeq 从运行目录 viral_contigs.fasta（回退 contigs.filtered）
  取该记录，clipboard.writeText（execCommand 兜底）写 FASTA，toast 报长度。
- **四件套跳转病毒注释分析**：行内按钮改为 🧩CDD / N / X / 🧪引物 →
  anaJump：cdd→#t-cdd，blastn/blastx/primer→#t-hom；sessionStorage
  vp_anajump 暂存（60s 有效），processAnaJump 在 tools 页加载/跳转后
  自动选运行（loadCaRuns→set→loadCaContigs）→ runAnalyze(action,
  contig, scope) → 滚动到卡片；t-hom 卡升级为
  「同源搜索复核与引物（BLASTN·BLASTX·Primer）」+ metabuli 式 tab 徽标，
  CA_ACTIONS.hom 增加 primer。
- **结果 metabuli 风格**（anaResultHtml 全作用域共用重写，参考
  <DEMO-IP>/metabuli）：BLAST 命中表 #/Accession(链接 NCBI nuccore|
  protein)/Description/Organism(斜体)/Identity(≥90/60/40 分档着色)/
  Align Len/E-value/Bit Score + 顶部 Best E-value/Best identity；
  CDD 命中前加**域条可视化**（按 from/to 定位到最长 ORF，10 色段，
  段名内嵌）+ 表 #/Accession(cddsrv 链接)/Short Name/E/Bit/From/To/
  Superfamily；Primer 表 #/引物对/F与R序列(等宽+⧉复制按钮+Tm/GC)/
  产物/罚分。表头 sticky、行 hover。
- **验证**：伪造 contigs 运行实测宿主自动补列与筛选候选、primer 字段
  （name/F_seq/F_tm/F_gc/R_*/product/penalty 全齐）；_it_platform 新增
  布局顺序/控件/跳转/metabuli 标记断言 PASSED；各页内联 JS
  node --check 通过。

## 五十二、Metabuli 示例病毒引入（2026-09-07，用户："看 metabuli 的示例病毒，我们的用"）

- 从 <DEMO-IP>/metabuli /examples 端点固化五个示例到
  databases/examples/：TMV（NC_001367.1，6395nt）、PVY（NC_001616.1，
  9704nt）、CMV（NC_002034/2035/1440 三分体 RNA 共 8623nt）、
  PSTVd（NC_002030.1 类病毒 359nt）、Mix All（四毒 6 条混合 25081nt）。
- ④「病毒 contig 深度分析」卡输入框下新增 metabuli 同款小按钮行
  [TMV][PVY][CMV][PSTVd][Mix All]（vEx，app.js EXAMPLE_METABULI 映射；
  PSTVd 自动把最短 contig 调到 200——359nt 低于默认 500 过滤线）。
- 端到端实测（TMV）：kunpeng+BLAST 分类 → Virgaviridae / Tobamovirus /
  Tobamovirus tabaci（E=0）→ 宿主自动补列 Plant (High) → 四件套 primer
  5 对 → 运行清理。_it_platform 补示例按钮/文件断言 PASSED。

## 五十三、分类表去除「选择运行」下拉（2026-09-07，用户："很不明白，去除"）

- Virus Sequence Classification 的手动「选择运行」下拉移除——表格自动
  跟随**最新一次**含 virus_classification.tsv 的 contigs_ 运行
  （loadAnaRuns 重写：anaCurRun + #anaRunLabel 标签显示当前运行名；
  无新运行不重复加载；30s 轮询/任务完成回调触发）。
- 所有 val('anaRun') 改读 anaCurRun；runAnalyze 的 ana 作用域动态取值
  特判（sc.run==='anaRun' → anaCurRun，cdd/hom 不受影响）；宿主自动/
  ID 复制/四件套跳转全部随最新运行。历史运行区（rh-contigs）仍可
  复用/删除旧运行。实测：TMV 示例跑完 → 最新运行即表格显示对象。

## 五十四、输出区 UX 三修：日志滚动/折叠、Downloads 胶囊、历史运行沉底（2026-09-07，用户反馈）

- **日志滚动不再跳顶**：injectStageLogs 每 2.5s 重绘是根因——重绘前
  snapshotLogScroll 记录每个 pre.logbox 的 scrollTop 与"是否贴底"，
  重绘后 restoreLogScroll：运行中且原本贴底 → 跟随新日志滚底；
  其余（已结束任务、用户上翻阅读中）→ 精确恢复原位置（WeakMap 记忆，
  取代旧的全局 data-autoscroll 一刀切）。
- **日志可折叠**：日志头新增 ▾/▸ 日志 按钮（logToggle，
  sessionStorage vp_logclosed 记忆状态；i18n c.logBtn 中英）。
- **Downloads 胶囊化**：taskResultHtml 的文件链接加 btn small 胶囊样式
  （tr-links 本就 flex wrap）。
- **历史运行沉底**：tools.html 八卡的 rh 容器从 tool-out 内移到各
  </section> 前（卡片最底部）；orf/annotation/host_removal/host_predict
  独立页移到产物 hint 之后；genome/primer 的 tool-out 本就是卡片末块。
- 回归：_it_platform PASSED（rh 数量断言不受移动影响）；app.js/i18n.js/
  各页内联 JS node --check 通过。

## 五十五、④分类报告区顺序调整（2026-09-07，用户："宿主预测统计放到病毒序列分类(contig明细)下面"）

④卡分类报告折叠区最终顺序：桑基流向图 → 分类表 → **病毒序列分类
（contig 明细，含宿主筛选/四件套跳转）→ 宿主预测统计 → 分类旭日图**。
中途曾把宿主统计误插进 ② identify 卡（sections 全在同一渲染页，需按
源码 section 边界核对），已纠正为唯一副本。_it_platform 顺序断言
（明细 < 宿主 < 旭日）+ 唯一性断言 PASSED。

## 五十六、并发闸门 + 存储面板 + 数据库分离（2026-09-07 晚，用户评估后要求实施）

### 并发闸门（P0-1，app.py TaskManager）
- heavy/light 两档 Semaphore（缺省 2/4，platform.json defaults 可调
  max_heavy_tasks / max_light_tasks，钳制 1-8/1-16）。
- start(name, fn, log_file=None, weight='heavy')：默认 heavy 保守；
  _run() 先 0.5s 轮询 acquire 排队，期间响应 cancel。
- **排队状态对外映射 running**（snapshot()），否则前端 status==='running'
  判断会让取消按钮消失。排队提示走 msg（"排队中 · 前面还有 N 个任务"）。
- LIGHT_TOOLS={'convert','genoplot','primer'} + 10 处下载类标 light，其余默认 heavy。
- cancel() 支持 queued；_recover_interrupted 把遗留 queued 一并标记失败。
- _trim_memory()：_run finally 调用；>60 个终态任务清 log/result/thread，
  >300 移出内存（磁盘 tasks/*.json 仍可回溯）。running/queued 跳过。
- list_all(limit=120)：running 30 行日志，终态 5 行。

### 轮询瘦身（app.js）
- setInterval(2.5s) → setTimeout 递归自适应：有任务 2.5s、空闲 8s、
  断连 2.5s 快速重试（connDown 标志）。
- **坑**：递归调度丢 setInterval 的异常免疫，schedulePoll 必须 try/catch
  兜 refreshTasks，否则一次异常永久断链。

### 存储与磁盘（P0-4）
- /api/storage（盘容量实时 + 目录占用后台扫描 + 10 分钟缓存）；
  /api/storage/rescan 强制重扫。目录统计走后台线程，绝不阻塞请求路径
  （databases 47GB/6 万文件 os.walk 实测数十秒）。
- /settings 新增「存储与磁盘」卡（storagePanel）；盘剩余 <20GB 全站顶部
  黄横幅（checkStorageBanner，插在 conn-banner 之后）。

### 数据库分离（任务：软件与数据库分开分发）
- 基础设施原本就有：config.apply_database_root + platform.json database_root +
  设置页「数据库目录」卡。本轮补：
- **Virus_Platform_Core/db_migrate.py + main.py db-migrate**：robocopy(/E /MT /Z 断点) →
  文件数+总字节校验 → 通过才 set_database_root。--mode copy|move、--dry-run、
  --check。目标须在平台根外、盘剩余 ≥ 源+10GB。写失败配置不动源目录原样。
- /api/settings 返回 dbs 就绪状态，前端切换 database_root 后即时提示缺哪个库。
- **打包改造**：scripts/package.py 默认不预置数据库（~1.9GB 软件包），
  --with-db 恢复开箱即用，--db-only 对已有发布目录补库。
  **坑**：argparse 曾写在执行路径中部，--help 会真的触发清理+打包（误删
  dist/）。参数解析必须放文件最前、任何副作用之前。
- **config threads=0 bug**：_load 曾 `max(1, int(threads))`，干净配置写 0 →
  新机被钳成 1 线程。已改 0/缺省=自动探测，显式 >0 才采用。

### 依赖审计（scripts/audit_deps.py + selfcheck ⑥）
- ast 扫描全部 import，对比 requirements.txt 与 installed，输出
  未装/未覆盖/可选三表，退出码非 0 表示缺必需依赖。打包前必跑。
- requirements.txt 补齐直接依赖声明：numpy/matplotlib/scipy/Jinja2/
  Werkzeug/urllib3/PyYAML/defusedxml（此前靠传递依赖）；可选注释加
  distinctipy/taxburst。
- 排除 docs/（bioaider_ref_src 是他项目参考源码）与 tests/。

### 回归
- tests/_it_concurrency.py：heavy 峰值≤2、light≤4、排队取消、名额归还
  无死锁、日志 30/5 行 —— ALL PASS。
- 28 项页面/API 冒烟 200；全量 py_compile + node --check 通过。
- dist/VirusPlatform 重新打包 1.9GB（无库），含 数据库对接说明.txt。

## 五十七、tool_runs 输出目录规范化（2026-09-07）

### 背景
tool_runs/ 积累 100+ 条目，混有补丁脚本、冒烟测试目录、工具日志残留
（orfipy_..._out），GUI 列表噪音大，无归档/清理机制。

### 设计
- 固定目录（`_` 前缀）：`_archive`（按工具分层归档）、`_tmp`（冒烟/
  一次性实验）、`_scripts`（运维脚本）、`_reports`（跨运行汇总）。
  各含 README.txt 说明用途与规则；GUI `api_tool_runs()` 跳过 `_` 前缀。
- 动态活动区：扁平 `<工具>_<YYYYMMDD_HHMMSS>`，正则
  `^([a-zA-Z0-9]+)_(\d{8})_(\d{6})$` 判定；不匹配者 organize 按规则归位
  （`.py`→_scripts，`_xxx`/`*_out`/含 test|smoke→_tmp，其余保留人工判断）。

### 落地
- `Virus_Platform_Core/tool_runs_admin.py`：ensure_fixed_dirs / status / organize /
  archive / clean；clean 的删除不可逆（shutil.rmtree），CLI 强制提示
  先 --dry-run。
- `main.py tool-runs <status|organize|archive|clean>` 四个动作。
- `Virus_Platform_Core/config.py`：默认 DIRS 补 `'tool_runs'` 键（此前只有设置自定义
  输出根时 apply_output_root 才写入该键，app.py 靠 `DIRS.get(...) or`
  兜底——隐性缺口）；apply_output_root 内懒导入 ensure_fixed_dirs，
  切换输出根时在新根自动预建固定目录。
- 存量归位：3 个补丁脚本→_scripts/，_k2_smoke2/_sracha_test/
  orfipy_viral_contigs.fasta_out→_tmp/。

### 验证
- `tool-runs status`：活动区 101 个/434.1MB，固定区四类统计正常。
- organize/archive/clean --dry-run 均 0 误报（存量已规范）；
  clean 缺期限参数正确报错退出码 1；--help 正常。

## 五十八、CDD / BLAST 语义归位（第一步：文案）+ 验证模块规划（2026-09-08，用户提问）

**起因**：用户问「保守域搜索（CDD）/ 同源搜索复核（BLAST）本质是病毒序列的
验证，应整合到『病毒识别和分类分析』还是保留在『病毒注释分析』」。

**结论：不二选一，按「点击目的」分流。** CDD / BLAST 天然横跨「验证」与
「注释」两阶段，按工具归类必丢一半语义。识别组做批量裁决，注释组做单条
深挖，共用同一套 `/api/tool/analyze` 与 `analysis/` 缓存，仅入口语义与结果
视图不同。规划文档：`docs/候选序列验证模块规划.md`。

**佐证（代码自陈）**：`tools.html` `anaJump()` 内注释明写「t-cdd / t-hom 属于
annotate 组，而 contig 明细表在 virus 组……必须把 g 参数也改成 annotate」——
用户看分类表想复核可疑 contig 却被弹到另一组，靠改 URL 参数才能跳过去，
正是归属错位的证据。

### 第一步落地（已完成）

- **app.py NAV_GROUPS**（annotate 组两条）：
  - `t-cdd` 「保守域搜索（CDD）」→「保守域与功能元件注释（CDD）」
  - `t-hom` 「同源搜索复核（BLAST）」→「同源比对与引物（BLASTN · BLASTX · Primer）」
  - desc 同步改为注释语义。「复核」是识别阶段的词，出现在注释组制造歧义。
- **tools.html**：`t-cdd` / `t-hom` 两卡 h2、卡片说明、两处 JS 注释同步改名；
  交叉引用文案一并更新；hom 卡补注「identity 分档仅反映最近同源水平，正式
  物种界定用 compare 组 SDT」。
- **未动（留第二步）**：`anaJump()` 跳转目标、分类表批量验证按钮——依赖
  验证模块存在，现阶段改动会成半成品。
- **回归**：`py_compile app.py` 通过；tools.html 内联 JS（13 个 block）
  `node --check` 通过。

### HMM 现状核实（规划输入）

- 三库齐备：VOG 4.7GB/49,116 + RVDB 2.7GB/13,679 + vFam 679MB/5,585，
  VOG 元数据与索引就位；本地 CDD 库 419MB + 精选病毒域表 1,580 条。
- pyhmmer **0.12.1 已装**（`C:/Python312/python.exe`，即 `启动平台.bat` 所用
  环境；managed 3.13.12 未装，不影响平台运行）。mmseqs2 已配置
  `tools/mmseqs/bin/mmseqs.exe`。
- `Virus_Platform_Core/hmm_annot.py` 完整（26KB），但**仅被 `orf_annot.py:705` 调用做层 2 功能
  兜底，无独立路由**。
- 关键提示：`annotate_orfs_hmm` 吃 **ORF 蛋白 faa**，不是 contig 核酸，
  验证模块需先补一步 ORF 预测（口径对齐现有「6-frame 最长 ORF ≥50aa」）。

### 第二步方案要点（待用户确认后实施）

**v2 修正（用户指认参照系）**：本模块等价于 MMPV-RNA
`virome_discovery_pipeline/filter_virus.py` 的 **02b_Filter**
（流程位：`02a_Identification → 02b_Filter → 03a_COBRA`，见
`virome_pipeline.py:11`、`:2630`）。据此重写 `docs/候选序列验证模块规划.md`。

**最重要修正：这是「过滤器」不是「裁决器」。** v1 的"判定四档含排除档"思路
作废——02b 哲学是**不确定时倾向保留**，保留规则四条：①已知病毒 taxid 命中
②关键词抢救（phage/virus/capsid/terminase/integrase/prophage…）③比对不可信
抢救（length<50 或 pident<30 或 evalue>1e-5，比对差不足以据此排除）④**全局
no-hit 也保留**。第④条是关键：no-hit contig 恰恰是新病毒来源，是最高价值
候选而非"排除"。改为 02b 口径 filter/strict 双模式。

- 证据链：DIAMOND vs viral_prot（本地）→ HMM 三库 → mmseqs2 CDD → NCBI 按需。
- 宿主归属走平台自有路线，比 02b 更直接：DIAMOND 命中 → `organism_tax.tsv`
  → family → `ictv_family_host.tsv` → plants/fungi/bacteria/invertebrates。
  02b 需 `cdd_classified_taxid.tsv`（从 CDD 域间接推断），平台无此表且
  `cddid_all.tbl` 无分类字段无法自生成 → **决定不补，用 ICTV 路线替代**。
- call 八类：viral_plant / viral_fungi / viral_bacteria / viral_invertebrate /
  viral / ambiguous / non_viral / unclassified；filter 模式保留除 non_viral
  外全部，strict 仅保留 viral*。unclassified 界面单独高亮（新病毒候选池）。
- 三步实施：后端 `Virus_Platform_Core/verify.py` → 前端 `t-verify` 卡 → 打通跳转（消除跨组 hack）。

### 参照系资产现状（核实）

- 02b 依赖 vs 平台：UniRef90 → `viral_prot.dmnd` 239MB（已有，更聚焦）；
  virus.taxid.txt → `taxonomy/nodes.dmp` 可生成；CDD 库+白名单 1,580 条已有；
  `cdd_classified_taxid.tsv` 缺失（走 ICTV 路线替代）。
- diamond `tools/diamond/diamond.exe`、mmseqs2 `tools/mmseqs/bin/mmseqs.exe`
  均已配置；pyhmmer 0.12.1 已装于 `C:/Python312`。
- 平台特有增强：02b 无 HMM 层，平台三库齐备，可捕获远源病毒。

### v3 补入 09b 严格层（用户："09b 还有 cdd 的过滤更严格些"）

MMPV-RNA 有**两个**过滤点，本模块需同时覆盖：
- **L1 = 02b_Filter**（`filter_virus.py`）：早期粗筛，CDD tier + 宿主归属，
  倾向保留（含 ambiguous / unclassified）。
- **L2 = 09b_Analysis_Verify**（`validate_rescue_cdd.py`，见
  `virome_pipeline.py:2244` Part V，由 9a 迁入）：后期精验，**严格在多了
  「属/科核心域」交叉验证**——不只问"像不像病毒"，还问"像不像它声称的那个属"。

L2 六档 EVIDENCE（`validate_rescue_cdd.py:259-330`）：PASS_CORE（属核心命中）/
PASS_CORE_FAM（仅科核心）/ PASS_VIRAL（有 T1 无核心，可能片段）/ AMBIGUOUS /
REVIEW / NON_VIRAL（疑似假阳性宿主）。核心域表含 core100 / core90 / aux 三集合。

**实现必带的坑**：CDD 多标签效应使同一家族多模型重复计数（**实测膨胀 ~6×**），
且跨界同源超家族（RT / 解旋酶 / pol / HSP70）有 viral 与宿主双版本模型，制造
虚假外源域信号。09b v2 解法：`convertalis` 取 qstart/qend 核酸坐标 → 聚类成
位点 → 按成员构成分 viral_only / dual / proximal_dual / confirmed_capture /
pure_clean 五类。

**资产现成，无需自建**（用户指认"应该有的"后核实——此前找漏了，
v3 初稿误判为"缺失需自建"）。全部在 `D:\桌面\延伸基因组\Cenote-Taker3\`：
- `docs/genera_core_aux_species_v3.tsv` 199 属 / 46 科，188 属有 core100
- `docs/family_core_aux_species_v1.tsv`（科级兜底）、`genera_core90_domains.tsv`、
  `all_211_genera_domains.tsv`、`suprafam_sites_final.tsv`（2.6MB，超家族位点）
- `results/cdd_classified_taxid.tsv` 43.9MB / 98,863 域（列 accession/category/
  taxname/taxid/title/div/smp_lineage）← **此前误判为平台缺失**
- `results/cdd_virus_final_v4.txt` 31,636 条白名单（比平台现有 1,580 条全 20×）

**编号兼容性已验证**：核心域表用 `pfam02407`/`cd23247`，而平台 `cddid_all.tbl`
确有 `425413  pfam00001  ...` 记录（NCBI CDD 整合 Pfam，accession 保留 pfam
前缀）→ 与 mmseqs2 CDD 搜索结果直接对得上，无需编号映射。

覆盖实测：Potexvirus(51) / Tobamovirus(34) / Potyvirus(199) / Tombusvirus(15) /
Cucumovirus(4) 均有 core100；**Carlavirus(71) core100 为空 → 必须回退 core90**。
注意 n_species 中位数仅 3（最大 467），≥5 的属仅 72 个 → 低样本属 core100 统计
意义有限，按 09b 原逻辑取 core100|core90 并集，界面标注 n_species 提示参考级。
落地：复制到 `databases/cdd/` 与 `databases/viral_prot/`，约 45MB。
实施步骤：第 0 步资产就位 → 第 1 步 L1 → 第 2 步 L2 → 第 3 步前端 → 第 4 步跳转。

## 五十九、类病毒库下载 + 验证模块按宿主/长度分流（2026-09-08，用户需求）

### 类病毒库就位

从 246 下载 11 个文件到 `databases/viroids-db/`（约 11MB），源
`/home/USER/database/virus-db/viroids-db/`：
`viroids.fasta`（9,353 条，header 带描述）/ `viroids.acc.fasta`（仅 accession）/
`viroids.taxonomy_info.tsv`（2.1MB，accession → realm..species 完整分类）/
`viroids.fasta.blast.db.*` 8 文件；另写 `db_info.json`。
打包 tar 传输（远程 753KB），非逐个 scp。

### 踩坑：BLAST v5 库在中文路径下不可用（平台级，全局适用）

- **现象**：`blastn -db databases/viroids-db/viroids.fasta.blast.db` →
  `LMDB runtime error: mdb_env_open: 系统找不到指定的路径`
- **排查链**：① 复制到 `C:\temp\` 纯 ASCII 路径 → **可用**（100% self-hit）
  ② 移除 `.njs` → 仍失败（不是 .njs 单独的问题）
  ③ 用平台 `makeblastdb 2.15 -blastdb_version 4` 本地重建 → **中文路径下可用**
- **根因：BLAST v5 库的 LMDB 后端无法在含中文的路径下打开。**
- **结论**：平台根目录 `D:\桌面\植物病毒分析平台\` 含中文，**所有 v5 格式的
  blast 库都不能直接用**，必须 `makeblastdb -blastdb_version 4` 重建。
  后续接入任何外部 blast 库一律照此处理。
- 已建 `viroids_v4`（9,353 条，耗时 0.43s），验证：self-hit
  `NC_039241.1 → ref|NC_039241.1| 100.000/333`、近缘 `gb|KT901877.1| 99.700/333`。
- 附：`makeblastdb 2.15` 即使指定 `-blastdb_version 4` 仍会生成 `.njs`，
  本例带不带均可用，最终保留为 `.njs.lmdb-disabled`（保险）。
  另 `blastdbcmd.exe` 存在于 `tools/Blast/bin/` 但未登记进 config。

### 验证模块设计变更（规划 v4）

用户需求：① 宿主分类选择（默认全部）② 按长度分流。

| 长度 | 归属 | 方法 |
|---|---|---|
| ≥1000 bp | 病毒 | DIAMOND **blastp** + mmseqs2 **CDD**（单选/双跑，并集或交集） |
| 200–1000 bp | 类病毒 | **blastn** vs `viroids_v4` |
| <200 bp | 不处理 | — |

分流理由：类病毒全长仅 246–401 nt、**无蛋白编码能力**，blastp/CDD 蛋白层方法
对其完全无效；反之 blastn 查类病毒库对 ≥1000bp 病毒 contig 无意义。

前端：宿主下拉（默认全部）→ 方法勾选 → 并集/交集 → 运行；结果分病毒 /
类病毒两个 Tab。后端 `Virus_Platform_Core/verify.py` 的 `filter_contigs(run_dir, host='all')`
先筛选后分流，两支分别 `run_virus_arm()` / `run_viroid_arm()`。
类病毒判据建议：identity ≥80% 且覆盖 ≥60%（200–300bp 放宽到 ≥40%）。

## 六十、候选序列验证模块 L1 实施（2026-09-08，用户需求）

验证模块（对齐 02b/09b）的 L1 后端 + 前端已落地，L2（09b 核心域精验）待实施。

### 关键设计决策（用户拍板）

- **不预测 ORF**：DIAMOND blastx + mmseqs 均自带翻译，直接吃核酸 contig
  （与 02b 同构，不落 ORF 中间文件）。
- **宿主归属唯一来源 `virus_classification.tsv`**（工具④产物，8 级谱系）；
  **不引入 organism_tax.tsv**——用户："完全不需要，不兜底没意义"。
  `ictv_family_host.tsv` 做科→宿主映射。筛选语义为「含」匹配
  （host.source 是复合值如 `algae+fungi+plants`）。
- 长度分流：≥1000bp 病毒支（blastx+CDD，并集/交集）；250–1000bp 类病毒支
  （blastn vs viroids_v4）；<250bp 跳过。下限 250 由用户从 200 提上来。

### 实现

- **Virus_Platform_Core/verify.py**（新增，约 450 行）：`verify(run_dir, fasta, host, methods,
  combine)` 主编排；`split_by_length`（重复 header 追加 `_2` 去重）/
  `run_blastx` / `run_cdd` / `run_viroid_blastn` / `viroid_call` / `virus_call` /
  `combine_pass` / `_load_cls_family` / `_filter_by_host`。
  过滤器语义：no-hit 保留为 `unclassified`（新病毒/新类病毒候选池）。
- **app.py**：`_tool_job_verify`（run 引用优先 viral_contigs.fasta，否则独立
  fasta；引用 run 时复制其 virus_classification.tsv 供宿主归属）；
  `TOOL_REGISTRY['verify']`；`NAV_GROUPS.virus` 第 4 条；新增
  `/api/tool/verify_result` 路由（calls.tsv + summary.json）。
- **tools.html**：`<section id="t-verify">` 卡（宿主下拉/方法/并集交集 + 结果
  分色表）+ JS（loadVerifyRuns / runVerify / loadVerifyResult / renderVerifyTable）。

### 坑与验证

- blastn 的 `-outfmt '6 qseqid sseqid ...'` 必须作为**单个参数**传（否则
  "Too many positional arguments"）。
- DIAMOND blastx stitle 的 organism 名在方括号 `[Tobacco mosaic virus]` 里——
  本想据此推 family，后因 organism_tax 弃用而移除（见上）。
- `ictv_family_host.tsv` 的 host.source 是复合值，精确 `==` 会漏判，改子串包含。
- **端到端通过**：example_mix.fasta（5 病毒 + 1 类病毒 PSTVd 359bp）正确分流，
  病毒支 blastx/CDD 双跑、宿主筛选 plants、类病毒支 blastn 均验证。
- 回归：py_compile app.py + Virus_Platform_Core/verify.py 通过；tools.html 内联 JS（13 block）
  node --check 通过；test client `/tools?g=virus` 200 + t-verify 卡 + 导航条目。

### 待办（L2 + 跳转）

- L2：`validate_cdd`（CDD tier + 属/科核心域 → 6 级 EVIDENCE）+ 超家族归并。
- 第 4 步：`anaJump()` 跳本组 `#t-verify`；分类表加「批量验证」按钮。
- 第 0 步收尾：复制 Cenote-Taker3 七个文件（L2 需要）。

## 六十一、验证模块收尾：跳转 + 粘贴 + 桑基/旭日重绘（2026-09-08，用户需求）

用户决策：**L2 不做**（09b 核心域精验 + 超家族归并取消）；第 4 步跳转 + 分类表
批量验证按钮 + 独立 fasta/粘贴输入 + 结果重绘桑基图和旭日图。

### 实现

- **批量验证按钮**：contig 分类报告区顶部「🛡 候选序列验证」→ `jumpVerify()`，
  跳本组（virus）`#t-verify`，无跨组 hack；`processVerifyJump` 自动选中该运行。
  （四件套的 `anaJump` 保持跳注释组不动——那是单条深挖，语义不同。）
- **输入三方式**：① 选 contigs 运行（viral_contigs.fasta）② 独立 FASTA ③ 粘贴
  （`pasteSeq` → `/api/paste_input` 已有，卡片直接复用）。
- **结果重绘**：`#t-verify` 结果区加两张 Plotly 图——桑基图（宿主→判定，
  含 unclassified 独立流）+ 旭日图（分支→判定→科）。

### 全流程验证（真实 contigs 运行）

用真实运行 `contigs_20260906_040347`（148 条 viral contig）端到端验证：
- 分流正确：121 条（27 条 <250bp 跳过），110 类病毒支 + 11 病毒支。
- **重要观察（非 bug）**：该批数据 94% 判 unclassified——因为 kunpeng 分类
  score 仅 0.001-0.003 的低置信短片段，blastx（viral_prot 库对远源覆盖弱，仅
  2 条命中 identity 45%/33%）与 CDD（白名单模式，PRK/COG/KOG 细菌/真核域不在
  病毒白名单内，正确被过滤）都拿不到强病毒域证据。只有 2 条命中 Potexvirus
  核心域（cd23246/47、pfam01443）判 viral。符合「过滤器不裁决」——unclassified
  保留给人看，不强判。
- 宿主筛选：host=plants 时 host_dropped=0 正确——该批数据所有科的 host.source
  都含 plants 或复合值（Rhabdoviridae=fungi+invertebrates+plants+soil+vertebrates），
  无纯真菌/细菌科可剔除。「含」匹配语义符合"找植物病毒"。
- 类病毒链路：完整 PSTVd 359bp → 正确判 viroid/Pospiviroidae；180bp 短序列
  正确被 <250bp 跳过。
- 粘贴链路：paste_input 200 → verify 全通。
- 回归：py_compile + node --check + test client 全绿。

### 遗留认知（重要，后续可参考）

- `viral_prot.dmnd` 对远源病毒灵敏度有限（blastx 只命中 identity 33-45% 的
  极近缘）；真正远源的要用 HMM 层（本模块默认关，未接）。若后续需要，把
  `annotate_orfs_hmm` 作为第三路接入即可。
- ICTv `host.source` 是"该科可感染宿主范围"而非"实际宿主"，复合值居多，
  宿主筛选是粗筛而非精确。
- 低置信 contig（kunpeng score<0.01）大量存在时，验证模块会如实返回
  unclassified——这是对的，但用户需理解 unclassified ≠ 非病毒。

---

## 2026-09-08 · 阶段 ③b verify 接入管道（闭环闭合）

### 背景与决策

此前 kunpeng contig 分类（③ 组装）与 verify 双证据验证是两个独立工具，
未串联。本轮把 verify 整合为管道阶段，置于 `03_assembly` 之后：

- **产物目录 `03b_verify`**（大王拍板，紧邻上游便于对照）
- **orf 输入收窄**：只吃 verify 判定 `known`/`novel` 的 contig（大王拍板）
- `domain_only` / `unclassified` 不下传注释，但 `unclassified` 保留为
  新病毒候选池，不删除
- 管道 ③ 补产 `virus_classification.tsv`，使 verify 不再依赖工具④
- verify 缺工具（DIAMOND / MMseqs2）时优雅剔除：WARN 跳过，不阻断下游

### 本轮修复的 3 个真 bug

**① `STAGE_VIEW` / `STAGE_OUTPUTS` 路径多写一层阶段目录前缀**

原写 `'verify': '03b_verify/calls.tsv'`，但 `_list_outputs()` 与 `view_ok`
已按 `os.path.join(sample_dir, dirname, view)` 拼路径（`dirname` 本身即
`03b_verify`），实际变成 `03b_verify/03b_verify/calls.tsv` → 文件不
存在 → 卡片 `outputs=[]`、`view=null`。
**修复**：改为 `'verify': 'calls.tsv'` / `['calls.tsv', 'summary.json']`。
教训：`STAGE_OUTPUTS`/`STAGE_VIEW` 的所有路径一律**相对阶段目录**。

**② `n_calls` / `out_dir` 在写盘之后才赋值**

`Virus_Platform_Core/verify.py` 里 `summary['out_dir']`、`summary['n_calls']` 原在
`json.dump` 之后 → 落盘 `summary.json` 缺这两键 → `_sum_verify` 显示
「共 0 条：known 1」。内存对象（C.verify）有值，落盘文件没有，这个
差异非常隐蔽。
**修复**：把两个赋值移到 `json.dump` 之前。

**③ `_load_cls(run_dir)` 只读根目录，管道布局下谱系全空**

分类表在管道模式下落于 `<run_dir>/03_assembly/virus_classification.tsv`，
而 `_load_cls` 只找 `<run_dir>/virus_classification.tsv` → 所有谱系字段
（taxid / taxon / realm..species）与 `host_source` 全空。
这个 bug 不报错、不阻断，只让 verify 的核心信息静默损失。
**修复（追溯版）**：`verify()` 与 `_load_cls()` 均加显式 `assembly_dir`
参数，管道调用点传 `a_dir`（`<sample_dir>/03_assembly`）。
与平台既有约定一致——`orf.py:L217` / `phylo.py:L492` / `primer.py:L182`
全用 `assembly_dir or os.path.join(sample_dir, '03_assembly')` 模式，
verify 是唯一的例外，本轮归队。`assembly_dir` 缺省时仍探测
`run_dir/` → `run_dir/03_assembly/`，兼容未传参的旧调用。

**两种入口的语义差异（关键背景）**

| 入口 | 分类表位置 | 谁负责放置 |
|---|---|---|
| 管道 ③b | `<sample_dir>/03_assembly/` | ③ 组装阶段产出 |
| 工具④ | `<run_dir>/` 根 | `app.py:2028` 主动复制 |

因此该表位置无法靠猜统一，必须由调用者声明。

**同类读取点全景（4 处，另 3 处本就正确、不动）**

| 位置 | 拼法 | 布局 | 判定 |
|---|---|---|---|
| `verify.py:248` | 显式参数 + 缺省探测 | 双 | 本轮修 |
| `phylo.py:511` | `a_dir/` | 管道 | 正确 |
| `viz.py:782` | `run_dir/` | 工具 | 正确 |
| `ncbi_submit/store.py:823` | `root/run/` | 工具 | 正确 |

### 实测证据

- `_stage_verify` 端到端实跑（真实 contig，118.9s / 165.4s 两次）：
  blastx 命中 NP_619711.1、CDD 命中 `pfam00998;pfam00729` tier=1 →
  判定 `known`。修复③后 calls.tsv 谱系完整：
  `taxid=12268 / taxon=Carnation ringspot virus / family=Tombusviridae /
  genus=Dianthovirus / species=Dianthovirus dianthi / host_source=plants`
- `_load_cls` **五场景单测**：①显式 assembly_dir（管道）②缺省探测
  （工具根目录）③缺省探测（03_assembly 兜底）④显式时优先于根、
  不误读根目录的表 ⑤两处皆无 → 返回空 dict
- **纯管道布局端到端**（77.3s）：分类表**仅在** `03_assembly/`、根目录
  无副本时，谱系 / host_source / n_calls 全部正确
- **工具模式端到端**（76.9s）：分类表在 `run_dir` 根时未受影响，谱系完整
- `verified_contigs()` + `extract_viral_contigs(only=...)` 收窄四场景：
  4 条 contig（known/novel/domain_only/unclassified）→ 收窄 2 条，
  fasta 实际内容仅 `NODE_1`/`NODE_2`
- HTTP API（中英双语）：14 阶段、verify 紧跟 assembly、
  `summary=共 1 条：known 1` / `1 total: known 1`、`outputs` 2 项、
  `view` 实测 HTTP 200 + `text/tab-separated-values` + 35 列
- 三个真实样品回归（ERR7586041 / NX-6 / SRR39909438）中英双语全绿
- `py_compile` + `node --check` + 登记项一致性检查全绿

### 已知未改（记录）

- 约 50 个任务归档 json（9/6-9/8 批次）丢失，不在回收站
- `results/ERR7586041/03_assembly/contig_classify` 是空目录，`spades/` 含
  10 个真实产物（spades.log / contigs.fasta / GFA 等）；两者均为 ③ 组装
  阶段的工作目录（`assembly.py:203/217` 写 spades、`:362` 写
  contig_classify），**不是残留**，未删。
- 管道页 `STAGE_PARAMS.subsample` 卡片是死控件：`collectParams()` 不读
  `#subsample` 元素，填值不生效（后端永远用设置页默认值，默认 0）。

---

## 2026-09-08 · 完整注释档去除子采样

### 需求

大王指令：「完整注释 去除子采样」。

### 改动

`webapp/static/app.js:1795` `PIPE_TEMPLATES` 的 `full` 档由
`stages: null`（全部可用阶段）改为保留 `stages: null` +
新增 `exclude: ['subsample']` 字段；`runTemplate()`（L1800）在展开前
按 `exclude` 过滤 `usable`。

用排除法而非显式全列表，好处是将来新增阶段仍自动进 full 档，
只有明确该排除的才写进 `exclude`。

```js
{ id: 'full', label: '🔬 完整注释',
  tip: '全部可用阶段（含宿主预测 / 功能注释 / 引物设计 / 基因组图），不含子采样',
  stages: null, exclude: ['subsample'] },
```

### 为什么前端改就够

`subsample` 阶段的执行体 `_stage_subsample`（`Virus_Platform_Core/pipeline.py:269`）是
参数驱动的：`if not (C.subsample and C.subsample > 0): 跳过`。
执行循环是 `for key in _stage_execution_order(stages)`，只跑**已选阶段**；
full 档排除后 `subsample` 不在 `stages` 里，函数根本不被调用。

`C.subsample` 的全部读取点（`pipeline.py:272/284/290/292/294/298`）
都在 `_stage_subsample` 内部，无其他阶段消费该参数，因此无需后端配合。

### 实测证据

- 展开验证（复刻 `runTemplate` 逻辑）：
  - `fast` 4 阶段、`std` 9 阶段、`full` **14 → 13 阶段**
  - full 展开为 `fastp → fq2fa → host → virus → assembly → verify →
    hostana → orf → orfa → phylo → primer → gbdraw → report`
  - 三档均不含 `subsample`；full 的 `hostana`/`orfa`/`primer`/`gbdraw` 保留
  - 模拟 `hostana` 不可用：full 降为 12 阶段，仍不含 subsample
- 后端 `_stage_execution_order` 三档实算：阶段数与传入一致，
  无 `subsample`，顺序为稳定拓扑序
- `node --check webapp/static/app.js` 通过

### 顺带查清（未改）

- 管道页 `STAGE_PARAMS.subsample` 卡片是死控件：`collectParams()`
  （`app.js:2048`）不含 `subsample` 键，页面填的数值不会传到后端；
  后端 `app.py:6129` 用 `_val('subsample', 0)` 取设置页默认值（默认 0）。
  三档现已全不含 subsample，该卡片只在手动单阶段运行时出现。
- 快速筛查档 `report` 声明 deps `['virus','assembly']`，档内无 assembly。
  后端不补齐依赖（deps 只约束已选集合内先后），report 照跑且
  `build_report` 各处有 `os.path.isfile` 保护可降级；但跑完后卡片
  `deps_ok` 判定会显示 blocked（`pipeline.py:1193`）。待大王定夺是否
  改 deps 或调整档位定义。

---

## 2026-09-08 · 修复 subsample 死控件（collectParams 未接线）

### 问题

管道页 `subsample` 阶段卡片里的数字输入框是死控件：`collectParams()`
（`app.js:2048`）不读 `#subsample` 元素，页面填的数值永远传不到后端，
后端一律用设置页默认值。

### 差集审计（精确口径）

抠出 `STAGE_PARAMS`（`app.js:1926`）声明的全部 `id` 与 `collectParams()`
实际读取的键做双向差集：

- 声明 24 个 → 收集 27 个
- **死控件：仅 `subsample` 一个**
- 反向 4 个（`threads` / `confidence` / `chunk_dir` / `force`）是全局参数，
  不在阶段卡片里，正常

### 改动

`webapp/static/app.js:2052` 在 `collectParams()` 补一行：

```js
subsample: +($('subsample')?.value) || undefined,
```

沿用平台既有惯例（`+($('id')?.value) || undefined`）：元素不存在或值为空
时归 `undefined`，`runStages()`（L2083）的 `filter(([_k, v]) => v !== undefined)`
会剔除该键，后端 `_val('subsample', 0)` 回落到设置页默认值。

### 实测证据

**行为矩阵**（构造受控 DOM 执行 `collectParams`）：

| 场景 | 输出 | 发送 |
|---|---|---|
| 卡片存在，值 100000 | `100000` | 发送 100000 |
| 卡片存在，值 50000 | `50000` | 发送 50000 |
| 卡片存在，值清空 | `undefined` | 不发，用后端默认 |
| 卡片不存在（阶段未渲染） | `undefined` | 不发，用后端默认 |

**端到端四场景**（真实 `fastp_R1.fastq.gz` 4.42MB，直调 `_stage_subsample`）：

- `n=0` → 正确跳过，零产物
- `n=20000` → 实际截取 4.42MB → 0.70MB，0.48s，`subsample.json` 记 `n_pairs=20000`
- 同参数复跑 → 断点复用，0.00s，日志「子采样已完成（20,000 对），复用」
- 改参数 `n=5000` → 触发重跑，0.17MB，0.12s

**API 层解析**（`_analysis_kwargs`，设置页默认 `subsample=0`）：

| 请求体 | `run_analysis(subsample=)` |
|---|---|
| `{'subsample': 100000}` | 100000 |
| `{'subsample': 5000}` | 5000 |
| `{'subsample': 0}` | 0 |
| `{}` / `null` / `''` | 0（回落默认） |

`node --check webapp/static/app.js` + `py_compile` 全绿。

### 关联

三档模板现已全部不含 `subsample` 阶段（详见上一节），该卡片只在手动
单阶段运行时出现；接线后手动运行时填值才真正生效。

---

## 2026-09-08 · 修复 report 阶段依赖声明过宽

### 问题

`STAGE_REGISTRY` 里 `report` 的 deps 声明为 `['virus', 'assembly']`，
但快速筛查档（`fastp → host → virus → report`）不含 assembly。

后果：跑完快筛后，若 report 尚未执行，`pipeline_overview`
（`Virus_Platform_Core/pipeline.py:1259`）判定 `status = 'ready' if deps_ok else 'blocked'`，
而 `deps_ok = all(ran_map.get(d) for d in STAGE_DEPS[key])`
（L1232，`ran_map` 按 PIPELINE_STAGES 顺序累积），assembly 未跑
→ `deps_ok=False` → report 卡片显示 **blocked**，用户误以为不可运行。

### 审计：全表 deps 逐个复核

| 阶段 | deps | 判定 |
|---|---|---|
| verify / consensus / hostana / orf / phylo / primer / gbdraw | `['assembly']` | 合理（输入都是 contig） |
| orfa | `['orf']` | 合理 |
| **report** | ~~`['virus','assembly']`~~ → `['virus']` | **过宽** |

报告是汇总展示层，核心输入是病毒筛查结果；assembly 产物只是让其更丰富
（多几个卡片）。`build_report`（`Virus_Platform_Core/viz.py:446`）对每个数据源都有独立
`os.path.isfile` 保护，缺 assembly 时逐项降级，不会失败。

### 改动

`Virus_Platform_Core/pipeline.py:635` 一行：

```python
'report':    {'deps': ['virus'],              'fn': _stage_report},
```

### 实测证据

**状态判定**（`pipeline_overview` 完整路径调用，NX-6：virus 已 done、
assembly 未跑）：

| 阶段 | 改动后 |
|---|---|
| virus | `done`（检出物种 2064，提取病毒 reads 14,826,661 对） |
| assembly | `ready` |
| **report** | **`ready`**（改动前 blocked） |
| verify/hostana/orf/phylo/primer/gbdraw | 仍正确 `blocked` |

**HTTP API**（服务重启后，中英双语全绿）：

| 样品 | virus | assembly | report |
|---|---|---|---|
| NX-6 | done | ready | **ready** |
| ERR7586041 | done | done | done |
| SRR39909438 | ready | ready | blocked（依赖确未满足，正确） |

**降级实测**：构造只含 `02_virus_screen`（无 assembly、无 host_removal）
的测试样品，`build_report` 0.45s 生成 50KB `report.html`，含病毒卡片与
分类树、HTML 结构完整；桑基图因缺 `01_host_removal/stats.json` 自动跳过。

**闭环实测**：同一测试样品跑 report 前 `ready` → 跑后 `done`
（摘要「报告已生成」、产物 `07_report/report.html`），全程 assembly 未跑。

### 附带清理

- 停掉两个测试遗留实例（PID 284256→8989、286292→8900）
- 重启主服务（8765）以加载新代码：Flask `debug=False` 不热重载，
  源码改动必须重启才生效（本轮一度因未重启导致 API 仍返回旧值）

### 排查教训

`pipeline_overview(sample_dir, ...)` 形参是**完整路径**，传样品名
（如 `'NX-6'`）会静默拼成相对路径并全部返回 `ready`，不报错。
调试时须传绝对路径，否则得到假阴性。

---

## 2026-09-08 · 修复 FASTQ read ID 后缀不匹配导致的宿主去除/病毒提取静默失效

### 问题

用 MGI 数据（`fastq/GQMIX.R1.fq.gz` + `GQMIX.R2.fq.gz`）建测试样品
跑全流程，出现一组自相矛盾的症状：

| 阶段 | 症状 |
|---|---|
| ① 宿主去除 | `kept_R1` 507,368 reads = 输入 reads（**一条都没剔除**），但 `stats.json` 声称 `dropped_pairs: 452400`、`host_ratio: 0.471364` |
| ② 病毒筛查 | `viral_R1/R2.fastq.gz` 均 **0 字节**，但 `summary.json` 声称 `viral_pairs: 11104` |
| ③ 组装 | SPAdes 报 `== Error == file is empty: ...\in_R1.fastq.gz`，退出码 64 |

三个阶段都「声称成功、产物为空」，属于典型的静默失效。

### 根因：read ID 后缀在两条比对路径上不一致

FASTQ 头有三种格式，平台此前只遇过第一种：

| 来源 | FASTQ 头 | 首 token 有无 `/1` `/2` |
|---|---|---|
| 纯 Illumina | `@A00151:1303:H527CDSXF:1:1101:4901:1597032` | 无（NX-6 即此类，故从未暴露） |
| **MGI/华大** | `@LH00278:577:23VC2NLT4:8:1101:1205:26992/1` | **有**（后缀与 ID 同 token） |
| SRA 下载 | `@SRR39909438.9 A00459:...:1000/1` | 无（后缀在第二 token，天然被丢弃） |

kunpeng/kraken 输出的 read 名一律**无后缀**
（`LH00278:577:23VC2NLT4:8:1368:3171:24708`），而 MGI 原始 FASTQ 头带
`/1` `/2`。seqkit `-f` 精确匹配是整行/整 token 匹配，两边对不上：
**实测前 100,001 个 read 名与 keep_ids 交集为 0**。

两个 bug 叠加掩盖了彼此：

1. `_seqkit_filter` 用 `-v` 反选，匹配失败 = 全部保留，`kept` 数看着正常；
   而 `dropped` 直接拿 `len(ids)` 充数，掩盖了「一条都没剔除」。
2. `_extract_reads_seqkit` 返回 `len(keep_ids)`（输入 ID 数）而非实际
   提取条数；seqkit 无匹配**不报错**，静默产出空文件 + 虚假计数，
   直到 SPAdes 才炸出来。

### 审计：匹配方案实测选型

修复方向有两个候选，用真实数据（452,400 个宿主 ID）实测对比：

| 方案 | 5,000 个模式 | 452,400 个模式 |
|---|---|---|
| **精确匹配 + 写进后缀** | 2.4s | **2.4s** |
| `seqkit grep -r` 锚定正则 | 197.6s | 卡死（实测 15 分钟未结束被 kill） |

精确匹配是 hash 查表 O(1)，耗时与模式数无关；`-r` 是逐 read 逐模式的正则
扫描 O(n×m)，45 万条必然灾难。**结论：放弃 `-r`，改用「探测输入后缀 +
把后缀写进模式文件 + 精确匹配」**。

顺带实测否掉的两个选项：

- 子串模式（不带 `-r` 的裸 ID）会**误匹**：`READ:1:100` 匹到 `READ:1:1000`
- 原样精确（不带后缀）：**0/5 命中**，输出 30 字节空 gz

### 改动

**`Virus_Platform_Core/host_removal.py`** — 新增 ID 规范化与后缀探测：

| 函数 | 作用 |
|---|---|
| `_strip_pe_suffix(s)` | 剥掉尾部的 `/1` `/2` |
| `_norm_read_id(name_line)` | 取首 token 并剥后缀，**两条比对路径共用** |
| `_id_patterns(ids, pe_suffix='')` | 生成 seqkit 模式行（逐行一个，精确匹配用） |
| `detect_pe_suffix(fastq_path)` | 读首条头探测实际后缀，返回 `''`/`'/1'`/`'/2'` |

`_seqkit_filter` 重写为「探测后缀 → 写模式 → 精确匹配（去掉 `-r`）」，
`dropped` 改为 `n_in1 - kept1`（seqkit `stats` 实测差值），不再用
`len(ids)` 充数；`filter_paired_fastq` / `filter_single_fastq` 的内置
Python 回退路径统一走 `_norm_read_id`。

**`Virus_Platform_Core/virus_screen.py`** — `_extract_reads_seqkit` 同样改为后缀探测 +
精确匹配（去掉 `-r`），并新增 `_count_reads_seqkit` 用 `seqkit stats`
实测产物条数，无匹配/配对不等时抛异常，不再静默虚报。内置回退路径
改用 `_norm_read_id`。

### 实测证据

**单元级**（真实数据，452,400 宿主 ID / 2,448 病毒 ID）：

| 环节 | 输入 | 结果 | 耗时 |
|---|---|---|---|
| 宿主过滤 | 507,368 reads | kept 54,968，**实测剔除 452,400**（与 stats 一致） | 13.2s |
| 病毒提取 | 507,368 reads | 返回 2,448 对 = 实测 R1/R2 条数 = 目标 ID 数 | 5.3s |
| 后缀探测 | — | MGI R1→`/1`、R2→`/2`、Illumina/SRA→`(无)` | — |

**全流程**（`/api/pipeline/GQMIX/run`，host → virus → assembly）：

| 指标 | 修复前 | 修复后 |
|---|---|---|
| 宿主剔除 | 0 条（声称 452,400） | **452,400 条**（实测一致） |
| `host_ratio` | 0.471364（虚报） | **0.89166**（真实） |
| viral fastq | 0 字节（声称 11,104） | **2,448 对**（实测一致，0.083MB） |
| 组装 | SPAdes 退出码 64 | **11 条 contig**，12,611 bp |

宿主占比从虚报 47.1% 修正到 89.2%——枸杞数据宿主 reads 占绝大多数，
后者才合理。三个阶段状态全部 `done`。

### 关联

- 修复前所有样品均为 Illumina 数据（无 `/1` 后缀），故此 bug 长期未暴露
- 新建 GQMIX 测试样品而非覆盖 NX-6：NX-6 的 virus 产物是真实跑出的
  （占 2.6GB），且 GQMIX 是不同数据
- NX-6 `input.json` 相对路径指向的原始 fastq 已删，无法实跑；
  `/api/pipeline/NX-6/run` 报 HTTP 400「路径不存在」

### 排查教训

- **同一份数据在两条比对路径上必须共用同一个规范化函数**。此前
  `_read_id`（不剥后缀）散落在多处，新增 `_norm_read_id` 后已全量替换
  （`grep _read_id\(` 复核确认无遗漏）。
- **统计值必须来自实测**。`dropped = len(ids)`、`return len(keep_ids)`
  这类「拿输入数当输出数」的写法，在工具静默失败时会产出完美自洽的
  假报告，比直接报错危险得多。
- **seqkit 无匹配不报错**，任何 `grep` 提取后都应实测产物条数。
- **性能方案要看复杂度而不只是能不能跑通**。`-r` 在小样本上只是「慢」，
  在真实规模上是「卡死」；5,000 vs 452,400 的耗时对比（197.6s vs 2.4s）
  才是有决策力的数据。
- **`force` 不会清 checkpoint 标志**。重跑前须手动删除 `.host_removal.done`
  / `.virus_screen.done` 与旧产物，否则阶段会被 skip 而显示 `done`。

---

## 全局上下文（当前样品 / 项目）+ 全局清理（2026-09-11，用户要求）

### 需求与背景

用户提出两点：①「样品创建 / 批量导入」里的**样品名 / 项目**如何与下游关联，
怎么更好地做成全局使用；②为避免切换模块导致各模块运行进度与关联丢失，
要一个**全局清理**按钮，一次清掉所有模块的相互关联与运行。

### 现状勘察结论（改造前）

- **样品名 = 唯一主键**：`_safe_sample_name()` 归一化后即目录名
  `run/results/<样品名>/`，全平台靠"扫目录 + 拼路径"找它，**没有注册表/数据库**。
  真正的消费者只有：管道页、结果中心、批处理队列、任务结果预览、
  **kvsuite 工具卡选样**（`tool_jobs.py` 拼 `sample_sheet.tsv`）、
  MSA/建树下拉（`refs.py`）、LOGAN 溯源（`logan_trace.py`）。
  其余 21 个专项工具只认显式文件路径，与样品名无关。
- **项目名 = 纯展示标签**：只写进 `00_prep/input.json` 与 `project.json`，
  仅在样品卡片 🏷、项目筛选下拉、kvsuite 选样弹窗里显示；
  **不参与任何路径拼接/计算/模块关联**，且改造前没有任何改名接口。
- **切模块丢的不是任务**：任务在服务端线程里跑，`startPolling()` 会在每页重新接上；
  真正丢的是页面级 UI 上下文（`curSample` / `curStages` / `curBranch` 等），
  因为一级导航是整页刷新。`vp_last_sample` 写了却从不读，是关键缺口。
- **清理入口全是局部的**：`/api/queue/clear_finished`（只清终态）、
  `/api/tasks/clear_finished`（不碰磁盘归档）、逐样品 `/clear` `/delete`、
  逐运行 `/api/tool/runs/<run>/delete`。全仓搜 `cancel_all|purge|stop_all` → 0 命中。

### 实现（用户选定口径）

用户明确选择：**非破坏性**清理（永不删结果文件）+ 样品与项目**都做全局** +
按钮放**顶部导航常驻**。

**后端**

| 文件 | 改动 |
|---|---|
| `web/tasks.py` | 新增 `active_ids()` / `cancel_all()` / `purge_all(include_archive=)`；并修 `_recover_interrupted` 误删 `queue.json` 的坑（见坑记录 8） |
| `web/samples.py` | `SampleQueue.reset_all()`：连 queued/running 一起清空并落盘 |
| `web/tasks_api.py` | 新增 `POST /api/global/reset`：cancel_all → queue.reset_all → 等待收敛 ≤2s → purge_all；返回各项计数 + `active_left`，并显式回 `kept_results: true` |
| `web/pages.py` | `/results` 每行补 `project`（新增 `_sample_project()`，口径与 `/api/samples` 一致：`input.json` 优先、清单兜底） |

`purge_all` 的边界写死：**只动任务状态**（内存记录 + `tasks/<tid>.json`），
`results/` `tool_runs/` `downloads/` `submissions/` 一律不碰；
`include_archive=False` 时复用 `delete()` 是错的（它会顺手 unlink json），
故改为直接 `pop` 内存 + 单独遍历文件。

**前端**

- `app.js` 新增 `VP_CTX`（真源 = `localStorage['vp_ctx_sample'/'vp_ctx_project']`）、
  `gxModal()`、`gxSyncProjSelect()`、`gxPickSample()`、`gxPickProject()`、
  `gxGlobalReset()`、`gxProjChanged()`、`gxApplyResultsFilter()`、`gxApplyPage()`；
  `document` 上的 `vp-ctx` 事件做跨页广播。
- 顶部导航由 JS 注入（与 `injectLangSwitch` 同法，零模板改动）：样品胶囊 /
  项目胶囊 / 🧹 全局清理。清理确认后会清浏览器关联键 + 整页刷新复位 DOM 态，
  结果经 `sessionStorage['vp_reset_result']` 带到刷新后以 toast 回显。
- 接入点：`/pipeline`（无 `?sample=` 时用全局样品自动接上，这是用户主诉）、
  `/samples` 与 `/pipeline` 的项目下拉、`/results`（新增项目列 + 筛选）、
  `/tasks`（新增「只看当前样品（<样品名>）」）、kvsuite 选样弹窗
  （默认勾选全局样品、按全局项目筛选、全选只作用于可见行）。
- `/results` 样品表新增项目列，`filerow-tr` 的 `colspan` 4 → 5。

### 验证

- `tests/_it_global_reset.py`：临时目录造数据（沙箱要求路径必须在平台内，
  故用 `run/_it_global_reset_tmp/`，跑完整体删除），覆盖 `_recover_interrupted`
  防幽灵与 `queue.json` 存活、`cancel_all`、`reset_all`、`purge_all` 边界
  （含 `include_archive=False` 真正保留归档，以及 `results/`/`tool_runs/` 分毫不动）、
  `/api/global/reset` 端到端、`/results` 项目字段。
- `tests/_it_global_ctx_ui.py`：起临时实例 + selenium(Edge headless) 做 39 项
  浏览器断言（导航注入、样品/项目弹窗、管道页自动接上并刷新后仍恢复、
  任务中心筛选、结果中心筛选、清理按钮 fetch 打桩并校验请求体 + 关联键清除 +
  刷新后 toast 回显、kvsuite 选样弹窗接入）。
  **坑**：服务日志不能挂 `subprocess.PIPE`——werkzeug 每请求一行访问日志，
  64KB 管道写满后服务阻塞在 write 上，页面加载随之超时（排查花了很久）。
- `python tests/_route_inventory.py`：新增 1 条路由，已 `--save` 更新基线（192 条）。

### 遗留 / 未做

- **项目名仍不可改**：没有 `POST /api/samples/<s>/project`，创建后只能手改
  `00_prep/input.json`。全局项目筛选因此只能筛选"创建时填对"的样品。
- `/api/tool/runs` 硬截断 30 条，第 31 条以后无接口可清（全局清理刻意不含它）。
- 另有独立发现（本次未改）：`consensus` 缺在 `PIPELINE_STAGES` /
  `STAGE_SUMMARIES` / `STAGE_OUTPUTS` / `STAGE_VIEW` 里，虽在 `STAGE_ORDER`
  且会写 `03c_consensus/`，但永不计入 `stages_done`、不显示进度卡；
  `_safe_sample_name`（建目录）与 `web/common._safe_sample`（读路由）两套
  归一化口径不一致，中文样品名会被压成 `A`/`1` 之类甚至撞名；
  `backfill_manifest_from_products` 无任何调用点（死代码）。
- **运行中的平台需重启一次**：前端静态资源会热更新，但新路由
  `/api/global/reset` 要重启服务才加载；未重启时点清理会收到 404，
  前端已给出「请重启平台」的定向提示。

---

## 阶段表治理 · 样品名归一化 · 清单回填（2026-09-11，承接上一节）

用户要求「扩大改动面，做好检查，好好修复」上一节勘察出的三处独立问题。
实际修的过程又牵出**四个同族 bug**——都是"改了一处、漏了另一处"。

### 1. ③c 共识阶段"半接入"（真 bug）

`consensus` 早已在 `STAGE_ORDER` / `STAGE_REGISTRY` / `STAGE_GROUPS` /
`DEFAULT_ANALYZE_STAGES` / `DEFAULT_STAGE_EST` 里，前端 `PIPE_FANOUT`
与 `STAGE_LABELS` 也早等着它——**唯独漏了 `PIPELINE_STAGES`**（以及配套的
`PIPELINE_STAGE_NAMES_EN` / `STAGE_SUMMARIES(_EN)` / `STAGE_OUTPUTS` /
`STAGE_VIEW` / `_MISSING_TOOL_MSG`）。后果：`pipeline_overview` 从不产出该卡、
`/api/samples` 的 `done/total` 漏算、`stages_done` 永不记录、跑完看不见结果。
磁盘上至今没有任何样品有 `03c_consensus/`，所以一直没暴露。

修法：补全上述 6 张表；`STAGE_VIEW` 指向 `coverage.tsv`（逐参考覆盖 + 是否
判定存在，最可读）；`_MISSING_TOOL_MSG` 加 minibwa 探测
（`consensus.consensus_engine()`），缺 minibwa 时显示 `unavailable` 并自动
排除出"依次运行剩余步骤"。
**刻意不给 consensus 加 skipped 规则**：`present` 为空是"真实阴性"结论，
不是"跑失败"，标成 skipped 会把它踢出 `stages_done`。

### 2. 前端仍用退役阶段键 `virus`（同族，两个真实故障）

② 病毒筛查（kraken2）2026-09-10 退役，阶段键由 `virus` 换成 `kvsuite`，
但 `app.js` 三处漏改：

| 位置 | 后果 |
|---|---|
| `renderPipe` 的 `known` 集合 + `byKey.virus` | ②b 卡片既不匹配显式分支、又被判为 extras → 渲染到 **⑩报告之后**，顺序错乱 |
| `STAGE_PARAMS.virus` 持有 `db_virus` 输入 | 没有任何阶段叫 `virus` → **病毒参考库目录输入框从不渲染**，`collectParams` 取到 null，用户无法在管道页指定病毒库 |

修法：`known` 用 `'kvsuite'`、改 `byKey.kvsuite`、`STAGE_PARAMS.db_virus` 挂到
`kvsuite`。浏览器测试断言顶层卡顺序为
`预处理 → kvsuite → assembly → report` 且 `kvsuite` 卡内有 `#db_virus`。

### 3. 样品名归一化：4 套实现收敛成 1 套 + 按磁盘真实目录名解析（真 bug）

历史上有 4 份各自为政的转换：

| 实现 | 规则 | 用途 |
|---|---|---|
| `pipeline._safe_sample_name` | 连续非法字符折叠成一个 `_`、去首尾 `._-`、截断 50 | 建目录 |
| `web/common._safe_sample` | 逐字符替换，不去首尾、不截断 | **读路径** |
| `logan_trace.safe_name` | 同 `web/common` | LOGAN 读路径 |
| `utils.sample_name_from_fastq` | 由文件名推导 | 生成候选名 |

前两者在中文/连续特殊字符/超长名下不一致：「样品A」建成目录 `A`，
读路径却去找 `__A`；带空格的手工目录、`_` 开头的目录同理读不到。
且规范化是**静默**的——前端 toast 显示用户输入的名字，磁盘上却是另一个。

修法：
- `utils.safe_sample_name()` = 唯一权威纯转换（**逐字节保持旧规则**，
  保证已有目录名不漂移）；`pipeline._safe_sample_name` 变成它的别名；
  `logan_trace.safe_name` 同样委托。
- 新增 `utils.resolve_sample_name(name, base=None)`：**以磁盘真实目录名为准**
  ——① 含分隔符/上跳的名字不原样接受 ② 原名命中 `os.listdir(base)` 就用原名
  ③ 否则退回规范名。读路径一律改用它，于是手工创建/历史遗留/中文/带空格的
  目录都能正常读写；第 ② 步的名字直接来自 listdir，天然不含分隔符，
  比字符替换更安全；越界仍由 `check_path(in_platform=True)` 兜底。
- 改动点：`web/common._safe_sample`（加 `base` 参数供 `results/_archive` 用）、
  `web/samples._sample_dir` / `api_queue_add` / `SampleQueue._run_entry` /
  `api_analyze`、`pipeline.run_analysis` / `run_report_only`、
  `web/tasks._sample_result_preview`、`web/refs._msa_group_dir`、
  `web/tool_jobs` 的 kvsuite 选样、`logan_trace` 的四处路径拼装与
  `meta.sample` 快照。
- **不保留中文**是刻意决定（`safe_sample_name` docstring 写明理由）：样品名
  要当目录名进 SPAdes/BLAST(LMDB)/minibwa 等不支持非 ASCII 路径的工具。
  代价是"用户填的名字 ≠ 实际目录名"，所以改为**显式告知**：
  `POST /api/pipeline/create` 回 `{sample, typed, note}`，撞名时报
  「样品名「样品A」会被规范化为「A」，而「A」已存在」；前端
  `createSample`/`batchCreate`/`scCreate`/`scBatch` 全部改用后端返回的
  **实际目录名**展示与入队（原先用用户输入的原文）。
- 刻意**不动** `gbdraw_plot` / `cds_export` / `gb_collection` / `phylo._safe_name`
  / `tool_results` 里同形的 `re.sub`：它们处理的是 FASTA 头 / contig ID /
  分类单元 / 基因名，不是样品目录，统一会改掉别人 glob 的文件名。

### 4. `backfill_manifest_from_products` 是死代码 → 接成显式 CLI（真 bug）

`project.json` 只由 `update_project_manifest` 在跑流程时写；老样品（或手工
拷进 `results/` 的目录）没有清单，项目名/输入/stages_done 全缺。反推函数早就
写好，grep 全仓**零调用点**，且注释里承诺的 `ensure_project_manifest()`
**根本不存在**（名字都是错的）。

修法：新增 `python main.py samples-backfill [--sample X] [--dry-run]`
→ `pipeline.backfill_all_manifests()`，只处理"清单缺失或来自回填"的样品，
已有真实清单的一律不动；`backfill_manifest_from_products` 加 `dry_run` 参数。
**刻意不在 `GET /api/samples` 里自动触发**：读接口不该写盘，且回填会用目录名
兜底出"项目名"，静默给所有老样品造出伪项目会污染项目筛选器。
边界如实写明：只填 `project` / `input` / `stages_done`，**不猜
`last_run` / `last_status`**（凭产物推断运行状态属臆测）。

### 5. 顺带清掉的测试残留（同族：假样品污染 UI）

`tests/_it_synteny.py` 把工作目录建在 `results/it_synteny_work`（不带 `_`
前缀），而 `/api/samples`、`/results`、项目筛选器、kvsuite 选样弹窗**只跳过
`_` 开头的目录** → 每次跑测试都在用户样品列表里多一个假样品
（本项目其它测试用的是 `_it_phylo` / `_smoke_phylo` 这种规范命名）。
改为 `_it_synteny_work` 并在结束时 `rmtree`；磁盘上遗留的那个（5 个合成
fixture、33KB、无任何产物）已删除。

### 6. `STAGE_WEIGHTS` 是没有读取点的死表（§4.4 更正）

它只有定义、没有任何调用（`_RunState.est` 走 `_est_stage_sec` →
`DEFAULT_STAGE_EST`），且只列了 14 项（漏 `subsample`），而
`docs/全链路复核报告_20260910_1545.md` §4.4 曾声称它"覆盖全部 15 阶段"。
已删除该表、把 `subsample` 补进 `DEFAULT_STAGE_EST`（原值靠 `('plain', 60)`
兜底，补上后取值不变），并更正该文档。

### 7. 新增的不变量守卫：`tests/_it_stage_tables.py`

阶段信息散落在后端 9 张表 + 前端 3 张表，靠人工对照必然漏。本测试用集合
相等/子集断言把这类漏项变成可执行检查：

- 后端 8 张表（含 `PIPELINE_STAGES`）与 `STAGE_ORDER` 双向相等；
  `STAGE_SUMMARIES(_EN)` 等价于 `STAGE_ORDER - {report}`（report 走特殊路径）
- `STAGE_GROUPS(_EN)` 覆盖且不重复、中英同构；`STAGE_VIEW` /
  `DEFAULT_ANALYZE_STAGES` / `STAGE_DEPS` 只引用已知阶段
- `PIPELINE_STAGES` 里 deps 必须先出现（拓扑顺序）
- 前端 `STAGE_LABELS` 必须覆盖全部阶段；`PIPE_FANOUT` 里的阶段都必须存在；
  `renderPipe` 的 `known` 集合必须覆盖全部阶段（否则有卡掉进 extras 乱序）
- **代码行里不得再出现 `virus` 阶段键 / `byKey.virus`**
  （注释与"数据库名 virus"不算）
- 实测 `pipeline_overview` 对含 `03c_consensus/` 的样品产出 done 卡
  （摘要 / 查看入口 / 产物齐全），且 `/api/samples` 的 `total` == 阶段数

### 8. 本轮新增/修改的测试

| 文件 | 覆盖 |
|---|---|
| `tests/_it_stage_tables.py` | 阶段表不变量 + ③c 可见性（新） |
| `tests/_it_sample_names.py` | 归一化规则/幂等/磁盘解析/防穿越 + 非规范命名端到端读写 + note 与撞名文案（新） |
| `tests/_it_manifest_backfill.py` | 回填 dry_run / 写盘 / 跳过已有清单 / `--sample` 过滤 / CLI smoke（新） |
| `tests/_it_global_ctx_ui.py` | 追加 [9] 管道卡顺序与 `#db_virus`、[10] 规范化提示（51 项） |
| `tests/_it_synteny.py` | 工作目录改 `_` 前缀 + 结束时清理 |
| `tests/_it_nav_responsive.py` | 顶部导航在 12 档宽度（1920→640）× 3 个页面下的自适应：横条不纵向溢出、无子元素跑到条外、页面无横向溢出、上下文/语言/📂 都可见、13 个一级链接都在、链接最多 2 行；宽屏必须回到单行 54px（新） |

回归结果：`_it_platform.py` PASSED、`_verify_fixes.py` 全通过、
`_verify_fixes_b.py` 通过、`_audit_routes.py` 193 条路由 0 个 500、
23 个页面全 200、5 个新增测试全通过、路由基线与 192 条一致（本轮未新增路由）。

### 遗留（本轮仍未动）

- 项目名依旧不可改（无 `POST /api/samples/<s>/project`）。
- `/api/tool/runs` 硬截断 30 条。
- 样品目录名改名后，目录内的**内嵌旧名副本**
  （`input.json.sample`、`project.json.project`、报告标题、
  `02b_kvsuite/plots/<样品>/`、`tool_runs/*/sample_sheet.tsv`、
  `logan/<job>/meta.json.sample`）不会跟着更新。
  `resolve_sample_name` 只让读路径能命中目录，消不掉这些陈旧副本。

### 补充：顶部导航宽度自适应（2026-09-11，用户截图报障）

用户反馈"页面上面的导航目录没有自动适应窗口"。根因与修法见坑记录 19。
实测结果（`tests/_it_nav_responsive.py`，真浏览器量 12 档宽度）：

| 窗口宽度 | 横条高度 | 排布 |
|---|---|---|
| ≥1900px | 54px | 品牌 + 13 链接 + 上下文簇，全在一行 |
| 1680–1180px | 80px | 链接一行，上下文簇（胶囊+EN+📂）整体换到第二行 |
| ≤1180px | 103–106px | 品牌独占一行，导航整行起第二行 |
| ≤720px | 139px | 链接 2 行 |

任何宽度下：横条不纵向溢出、无子元素跑到条外、页面无横向溢出、
上下文/语言/📂 入口都仍可见可点。

### 事故记录：本轮我自己造成的任务记录清空

跑 `tests/_audit_routes.py` 做回归时，它按设计对每条路由发了一次空 POST，
其中包含新加的 `POST /api/global/reset` → 200，
于是 `run/tasks/*.json` 里约 100 条历史任务记录被删除（详见坑记录 21）。
- 未受影响：`results/`（样品与产物）、`tool_runs/`、`run/logs/tasks/`
  下的 698 个任务日志文件、队列文件本身。
- 已受影响：任务中心的"历史归档"列表（记录索引），无法恢复。
- 已修复：接口加 `confirm` 必填令牌（空 body → 400）；巡检脚本加
  `DESTRUCTIVE_POST` 跳过表，并把快照比对补上"删除检测"。
  `tests/_it_global_reset.py` 新增「无 confirm → 400 且不清任何记录」断言。

## 一键分析（全流程）卡补齐输入 / 参数 / 输出（2026-09-11，用户截图报障）

用户截图报障：「一键分析（全流程）页面没有输入，也没有输出，也没有参数，啥都没」。

### 根因（三条叠加，缺一都不会是"整页空白"的观感）

1. **没有输入、没有参数**：`t-kvchain` 卡只有「⚡ 一键运行 / 📂 最近汇总」两个按钮。
   `runKvchain()` 读的是 **t-kvsuite 卡的 `kv_*` 字段**（`kv_samples` / `kv_ref` /
   `kv_min_cov` …）——而那张卡属于同一模块组的**另一个二级模块**，`showModule()`
   一次只显示一个卡片，所以站在本页既看不到那些字段、也无从修改；表单值全靠
   另一张卡里"上次填过什么"，本页表现为"没有任何参数"。
2. **没有输出**：本卡没有 `.toolrun` 容器（`#toolrun-kvchain`）。平台的口径是
   "运行日志/结果/产物平铺在工具卡内"（app.js `injectStageLogs` + `toolRunHtml`），
   而 `TOOL_LABELS` 里**从未登记 kvsuite / kvchain（virchain 同样漏了）**，
   于是即便有容器也匹配不到任务名，日志只会进全局任务抽屉，卡片输出区恒空。
3. **没有结果**：跑完只有 `#kvChainResult` 里的"流程各段"表 + 汇总目录链接
   （且 `_chain.json` 里 `steps[0].summary` 只有 stage/engine/n_confirmed 等标量），
   看不到任何鉴定/过滤/共识/变异表与图——相比另一张卡"什么都没有"。

### 修法

- **本卡自带一整套表单**（`kvc_*` 命名，17 个 input + 2 个 select）：样品（📋 多选，
  `kvPickSamples(targetId)` 支持指定回填框）、参考库、参考注释 TSV、比对引擎、
  鉴定库、覆盖率/深度/reads/泊松四阈值、变异 QUAL/AF/氨基酸标签、NCBI 邮箱、
  线程数、扩展变异 / SNPGenie / 不下载注释 / 画全部变异开关。
  默认值与 t-kvsuite 卡一致；`runKvchain()` 只读本卡 `kvc_*`。
- **参数互取**：`kvchainPullKvsuite()`（按钮「⤓ 取该卡参数」）按
  `KVCHAIN_FIELD_MAP` 把 t-kvsuite 卡当前表单值整份搬过来。两套表单不做静默同步，
  避免"改了这边没改那边"的隐性漂移；固定跑全五段（`stage=all`）由后端
  `_tool_job_kvchain` 保证，前端以文案明示，不设"跑哪些段"下拉。
- **输出区**：新增 `<div class="tool-out">` + `#toolrun-kvchain`；app.js
  `TOOL_LABELS` 补 `kvsuite` / `kvchain` / `virchain` 三个键，
  任务名（工具·已知病毒识别与定量 / 工具·病毒定量与共识·一键）即可匹配到卡内输出区，
  实时日志与产物下载就地可见（此前 virchain 卡同样是空的）。
- **内联结果**：`renderKvsuite` / `renderKvVariantPlots` / `renderKvEvo` 增加
  `ids` 参数（默认 `KV_IDS` 指向 t-kvsuite 卡的 `kv*` 容器）。一键流程跑完或点
  「📂 最近汇总」时，`kvchainLoadInline()` 从 `_chain.json` 的
  `steps[0].run`（= kvsuite 子运行）拉 `/api/tool/kvsuite_result`，
  用 `KVC_IDS` 渲染进本卡的 `kvc*` 容器：计量徽标 + 鉴定表 + 过滤表 + 共识 QC +
  变异注释 + 变异可视化 + 演化分析 + 深度图，并给「📊 在已知病毒识别与定量卡中打开」
  （`kvchainOpenKvsuite()` 走左侧模块条点击 → `loadKvsuiteResult(run)`）。
- `loadKvIndexOptions()` 同时填充两个鉴定库下拉（`kv_index_dir` / `kvc_index_dir`）。
- 示例体系：`examples.js` 的 `fillKvsuiteExample()` 同时填 `kv_samples` 与
  `kvc_samples`，本卡点「✨ 示例」不再毫无反应。
- i18n：新增/改写 `tk.kvChainHint`、`tk.kvChainStages`、`tk.kvChainPull(_Title)`、
  `tk.kvChainRunHint`、`tk.kvChainStepsSec`、`tk.kvChainSubRunLbl`、
  `tk.kvChainOpenKv`、`tk.kvChainNoSub`（zh/en 同步）。

### 回归

- `tests/_it_platform.py` 新增 **3c-2** 段：断言本卡含样品/参考/注释/鉴定库/引擎/
  阈值/开关/线程数 9 类输入、`#toolrun-kvchain` 输出区、两张结果表容器、
  运行与参数互取按钮；断言 `runKvchain` 只读 `kvc_*`（不再出现 `val('kv_samples')`）、
  内联走 `KVC_IDS`、`TOOL_LABELS` 登记 kvchain/kvsuite。
- `tests/_check_kvchain_ui.py`（真浏览器，自起临时实例，已挂进 `_run_all.py`，新）：
  本卡可见输入 ≥15 项、14 个容器齐全、「✨ 示例」落到 `#kvc_samples`、
  「⤓ 取该卡参数」搬运值/开关、未选样品点运行有提示；并用 `page.route` 桩住
  `chain_result` / `kvsuite_result`，验证「汇总清单 → 子运行 → 内联 kvc* 表格」
  这条链路渲染正确且**不污染 t-kvsuite 卡容器**，全程无 pageerror。
  （也可 `--port 8765` 直接对已在运行的平台跑。）
- `tests/_audit_contract.py`：i18n zh/en 1718 键对齐、使用键全部有定义，无硬性差异。
- 渲染回归：`/tools` 各组页面内联 JS 经 `node --check` 全部通过；
  `app.js` / `examples.js` / `i18n.js` 同步通过。
- 注：本机 `run/results/` 无样品（用户已归档），`_it_platform.py` 的
  "样品列表 > 0" 断言属既有数据依赖项，非本轮引入；补一个临时样品目录后
  全量 PASSED。

## 裁剪版系统地理移植 + 两个"永远报绿"的巡检陷阱（2026-09-16）

完整清单/证据见 `docs/主平台同步_裁剪移植_20260916.md`；这里只记**踩坑**。

### 1. 有组无入口 = 卡片到不了（组落地路由的隐含前提）

09-15 拆分把 `compare` / `phylodyn` 从 `NAV_GROUPS` 摘掉后，`#t-rdp` / `#t-rtt` /
`#t-phylogeo` 三张卡变成"路由在、卡片在、界面上到不了"：tools 页的
`currentGroupFromURL()` → `showLanding/showModule` 要求组在 `NAV_GROUPS` 里，
否则 `?g=phylodyn` 掉回 sample 组、`#t-phylogeo` 深链落不到卡片。
**只补组还不够** —— `_nav.html` 的一级导航项是硬编码列表，不补就是"只能手敲 URL"。

### 2. `tests/_check_pages_console.py` 的判据曾经恒真（巡检静默全绿）

两个缺陷叠加，缺一都不会"全绿"：

1. 单页模式末尾写的是 `main()` 而不是 `sys.exit(main())` → 子进程**永远 exit 0**，
   而父进程只看 `returncode` ⇒ 页面坏成什么样都打 ✔。
2. 父进程 `subprocess.run(..., text=True)` 按本地代码页（GBK）解码子进程 stdout，
   而平台自己的模块会把**子进程** stdout 重配成 UTF-8（`/settings` 加载全量数据时
   打印的自检报告就是 UTF-8）→ 读线程抛 `UnicodeDecodeError`、stdout 变空 ⇒
   连"结论行"都读不到，而缺陷 1 又把空输出当成功。

现在：`sys.exit(main())` + `encoding='utf-8', errors='replace'` + **以子进程自报的
末行结论（`OK ` / `BAD `）为主判据**、`returncode` 兜底，两者不满足一律 ✘。
**负对照**（把子进程结论强改 `BAD` 的临时副本）：修前"全部 ✔"，修后 `rc=1` 且逐页 ✘。

### 3. 脚本化整文件写入会静默改换行符（CRLF → LF）

用 Python 脚本做批量文本替换时，`io.open(..., newline='')` 读进来是 `\n`，
写回就把 CRLF 文件整篇变成 LF ⇒ `git diff` 出现上千行噪声（正是 `.gitattributes`
里 `* -text` 那段注释警告的事）。**对策**：改完用字节比对核一遍 CRLF/LF 计数，
必要时 `.replace(b'\n', b'\r\n')` 还原。本次涉及 `i18n.js` / `tools.html` / `phylogeo.py`。


## 2026-09-17 dsRNA 设计全链路（工具箱 t-dsrna）改造

原 `_tool_job_dsrna` 只是 dsRNAmax_det 裸调（产出 result.csv），本次升级为
完整链路，引擎侧三个新模块（均为函数内惰性 import，**打包 spec 的
hiddenimports 必须登记**，已加）：

| 模块 | 职责 |
|---|---|
| `dsrna_pipeline.py` | 编排：exe 探测/内存守卫/maximin 选窗/面板扫描/引物/成品 QC |
| `dsrna_scoring.py` | dsRIP 效价 + 安全合成评分（**平台唯一口径**，合并 C-host 管线的 dsrip_scoring/safety_scoring 两份复制实现） |
| `dsrna_offtarget.py` | uint64 排序数组 21-mer 脱靶扫描器（双链入索引、0→1→2mm 分层、磁盘缓存） |

数据资产（发行包需带）：
- `databases/dsrna/lookup.db` + `siRNA_parameters.txt`（dsRIP 权重/热力学表）
- `databases/dsrna/lethal_lists/*.tsv`（5 物种致死基因，按**裸种名**匹配，
  面板文件名里的 `_refseq_rna` 后缀会被剥掉再匹配）
- `databases/dsrna/panel/*.fna`（RefSeq mRNA 脱靶面板，8 物种约 750MB，
  首扫建索引并缓存到 `panel/_cache/`，之后秒级；**冻结分发可不带**，
  空面板时工具自动降级为纯效价选窗并在日志说明）

链路设计要点（背景：2026-09 对 C-host dsRNAmax 管线的审计）：
1. **打分对象=交付分子**：引物限定在效价×安全最优窗内设计（scope=
   best_window），引物定稿后对真实扩增子重跑效价 + 面板脱靶复扫（QC）。
   旧控制台/管线在 300nt 臂上打分、却交付 Primer3 随手挑的 215bp 产物，
   两者的脱靶面完全不同。
2. **median 遮蔽离群株**：多候选按（最差目标，中位，总和）maximin 重排；
   命中 0 的"盲目标"与 <25% 覆盖的目标显式告警（单臂只能覆盖共享 21-mer
   的同源目标，异源目标须各自出臂做鸡尾酒）。
3. **绝对红线**：致死基因 0 错配命中一票否决（`red_lines`），不随臂内
   相对归一化被稀释；安全分仍是 dsRIP 口径的臂内归一（最安全=100）。
4. exe 兼容四件套：`-otKmerLen` 恒显式传（kmerLen<21 不传直接 fatal）；
   `-seed` 经 `probe_exe()` 探测后才传（stock exe 会 exit 2）；`-csv` 写
   失败时 exe 退出码仍为 0 且不再输出臂序列，靠解析失败兜底；内存守卫按
   **160B × 去重 k-mer × iterations** 估算（iterations=100 实测 13-16KB/
   kmer 反推），超预算抛 `MemoryBudgetExceeded` 而不是等 OOM。iterations
   默认从 1000 回落到 100——每个迭代完整拷贝一份 k-mer 表且全部并发存活，
   1000 次 ≈ 10 倍内存（18 条高分化 15kb 靶标 ≈ 40GB）。
5. 同 seed 两次运行臂逐字节一致（det 构建，tag `v1.1.15-det.1`，源码在
   C-host dsRNAmax 仓库 commit 1e0c19e）。

测试：`python tests/_it_dsrna.py`（端到端/确定性/内存守卫/与 C-host 管线
v4 口径对拍 w933=7..296）；致死红线阳性用例见该次改造的会话记录
（Apis NM_001160052.1 → 280/280 lethal mm0 → red_line）。

已知边界：
- 植物特有生物学未建模：DCL2 22nt / RDR6 transitivity 会次级放大脱靶，
  21-mer 模型在植物里偏乐观（报告未声明，属已知科学局限）。
- `offtarget_lethal.tsv` 只含致死命中行（全量命中行体积大，扫描器返回值
  里 `sirna_tiers` 已含逐 siRNA 分层，需要时可再导出）。

### dsRNA 测试夹具：论文已发表序列（tests/fixtures/dsrna/）

`_it_dsrna.py` 第 5–7 节用 dsRNAmax 原论文（Fletcher et al. 2025, NAR Genom
Bioinform, lqaf064）补充材料做真实数据夹具：

- `tef17_tef21_arms.fa` = Suppl. Table 1 的两条已发表臂（300nt，有三种
  基础生测数据背书），经 PyMuPDF 从 PDF 文本层抽取，非手工抄写；
- `ce_ef1a_alnpartners.fa` = Suppl. Figure 3 的 C. elegans EF1a 比对段
  （图无文本层，来自图像 OCR）。**保真锚点是同源度**：TEF-17 前 168nt
  对 Ce-EF1a 实测 85.12%（图注 85.1%，几乎精确），TEF-21 实测 79.00%
  （图注 79.3%，差 1 个碱基，容差内）。

⚠ 论文内部口径差异（勿当 bug 修）：Figure 1/4 截图里的运行产物与
Table 1 最终发表序列**不是同一分子**——共享保守前缀，但 GC 不同
（图内 45.9%/45.7% vs Table 1 的 45.00%/42.67%）。测试断言只用
Table 1 长度 + Figure 3 同源度，不用图注 GC。

语义断言：发表设计 = 与 C. elegans 无连续 ≥17nt 精确匹配 → 我们的
扫描器在 21-mer 口径下 mm0 必为 0（TEF-17/21 实测均 0，mm1 9/18）。
第 7 节把 TEF-17 喂回全链路：v1.1.15 的 constructLen==目标长 边界
（off-by-one 修复的标志性场景）出臂为目标精确子串，引物/扩增子 QC 全通过。

#### 2026-09-17 代码审查修复（第二轮回归）

审查发现并修复 4 个 P1（隔离株覆盖日志 TypeError——`isolate_coverage` 返回
结构改为 `{overall, per_isolate}` 后 `min(m.values())` 会在 float/dict 间比较；
空面板选择语义——`panel=''` 此前被归一化成 None 而误扫全部 750MB 面板，
现三态化：不给=全部 / ''或'none'=跳过 / 列表=所选；JS `+val('ds_mm')||1`
把"0 错配"吞成 1；result.csv 固定拷 0 号候选的 CSV 而非 maximin 中选者——
候选现在携带 `orig` 原始编号）+ 2 个口径问题（报告里面板 bp/n_records 恒 0
——indexes 未透传；lethal_rows 的 sirna_idx 0 基与 TSV 契约不一致，改 1 基）
+ 卫生项（删 `extract_amplicon` 死代码与 `assign_safety_scores` 未用形参；
`build_lethal_index` 进程内缓存——同任务臂扫描+扩增子复扫免整文件重读；
`scan_arms` 支持 cancel 逐 siRNA 中断；revcomp 文档注明仅 DNA 字母表）。
`tests/_it_dsrna.py` 增至 9 节：新增 isolates 分支端到端（含 5% 突变衍生株
的逃逸信号断言）、`panel=''` 跳过、result.csv 中选臂校验、max_mm=0 分层
（11 跨突变窗 None + 869 精确窗 0）。

## 2026-09-17 病毒浏览器「输入包」导出：两个按钮、一份信封

需求（用户）：「停止，还是本地部署吧，修改 explore 导出 fasta 和 VirPhyKit 需要的
元数据信息文件，并且进行预过滤」→ 随后追问「我们之前的脚本呢，导出我们之前流程
支持的格式，是不是更好些」→ 定案**拆两个按钮**、**不信 GenBank 标记、按长度判完整**。

### 唯一的真冲突是 date 列编码

信封本来就已经一样了 —— 我们自己的 BEAST 交接包（进化平台
`beast_handoff.py:471`）写的表头也是 `name,date,location`，FASTA 头也与 `name`
列逐字相等。冲突只在 `date` 这一列的**写法**，而且**单向不可兼容**：

| `date` 写法 | VirPhyKit/TreeTime | 自家 `virome_phylo_pipeline` |
|---|---|---|
| `2018-06-29` | ✔ 2018.49178 | ✔ 2018.49444 |
| `2018.49589` | ✔ 2018.49589 | ✘ **2018.95145，偏 0.46 年且不报错** |

原因：`utils/decimal_year.py` 的 `to_decimal_year(strict=True)` 只认
`YYYY[-MM[-DD]]`，把小数年当成「年=2018.58 + 6 月 15 日」再折一次；越界年份还会
直接 raise，让服务器整个 BEAST 作业中止。TreeTime 的 `parse_dates()` 两种都收。
→ `YYYY-MM-DD` 是唯一两边都读对的写法，但 VirPhyKit 的 Example 用的是小数年
（`2007.569473`），所以**各出一个按钮，不混用**：

| 按钮 | 路由 | date 列 | 产物名 |
|---|---|---|---|
| 导出 VirPhyKit 输入包 | `POST /api/explorer/export/virphykit` | 小数年 `2007.56947` | `virphykit_*.zip` |
| 导出管线输入包 | `POST /api/explorer/export/pipeline` | `YYYY[-MM[-DD]]` | `pipeline_*.zip` |

两包**除 date 列外逐字相同**（name/location/FASTA 全等），巡检里有断言锁这条，
因为它是「两个按钮」这个设计的全部前提。

### 完整性：不信 `Nuc_Completeness` 标记，按长度判

标记 `complete` 直接放行（这信号只在「说完整」时可靠）；标记非 `complete` 时按长度判：
落在同物种 `complete` 记录长度**中位数**的 `[0.90, 1.20]` 倍内就留。
该物种一条 `complete` 都没有时（全库 6,508 个物种里有 1,985 个），退化为
**它自己最长序列**的长度。

数据依据：PVY 标 `complete` 的 685 条长度 9376–9800（很紧），而标 `partial`
的最长到 9762 —— 只认标记会白丢 269 条本该留的全长序列。全库同口径 66,552 →
81,137（+14,585）。跑在真实筛选上：PVY 588 → **805 条**（靠长度救回 217 条）。
档位取的是本平台既有判据：下界 0.90 = `align_qc.py:38` 的 `DEFAULT_MIN_LENGTH`，
上界 1.20 = `contig_annot.py:245` 的「近完整」。
**不能用属平均做基准**：`databases/genus_lens.tsv` 的 Potyvirus 平均只有 5292 bp
（真值 ~9700，被片段拉低），拿它判只命中 2 条。

### 顺手修掉一个互通缺陷：`metadata.csv` 里不能有引号

TreeTime 读元数据走 `treetime/utils.py:251` 的 `parse_dates()`，对非 `.tsv` 用
`sep=r'\s*,\s*'` + `engine='python'`，**不按引号包裹切分**。所以 `"Anshun, Guizhou"`
这种字段会被从逗号劈开 —— PVY 那个包里 43/588 条 `location` 就这样被截断（还带个
孤立引号），截断后的区域在 Mapping.txt 里永远匹配不上 → 全被标成 Unknown。
现在所有会破坏 CSV 的字符（逗号/引号/制表/换行）一律压成空格（`csv_safe`），
写出来**一个引号都没有**，TreeTime 的正则、pandas 默认、R `read.csv`、
最朴素的 `split(',')` 读出来完全一致。巡检里既有断言也有相应的灵敏度负控。

### 契约来源（外部目录，只按方法对齐，不复用代码）

VirPhyKit 是 GPL v3，本平台**只做方法对齐、不引入其代码**。三个约束从源码读出：
`Treetime/function_treetime.py:131-134` 读 `name/date/location` 并在 `:141-148`
**硬性中止**于「叶缺日期/缺地点」；`Treedater/function_treedater.py:25` 用 R
`read.csv`；`Subsample/fuction_subsample.py:37` 的
`re.search(r"(?:^|>|_)([A-Za-z0-9]+)(?:_|$)")` 从头里取**第一个** alnum token。
所以序列名口径定为 `{区域}_{登录号}_{小数年}`、区域去掉全部非 ASCII 字母数字字符。

### 巡检与验证

新增 `tests/_check_virphykit_export.py`（8 节，已登记进 `_run_all`）：小数年解析
（8 正 8 负）、区域/name token、长度判据（中位数 vs 均值 + 闭区间边界 + 基准缺失负控）、
预过滤夹具（含「标记 partial 但长度到位必须留下」的正控与 6 条负控）、
两种日期口径对拍真下家、契约灵敏度自检、两路由端到端真库、可选 GeoSubsampler 实跑。

**负控（证明会红）**：A 把 `statistics.median` 改成 `mean` → 1 项失败；
B 让 `_metadata_text` 忽略 `date_mode` → 4 项失败、退出码 1。改回后 md5 逐字节复原。
`tests/_check_explorer.py` 重跑全绿（证明 `filter_records` 抽取与路由重构没扰动
原有两个导出与 5000 行表格口径）。浏览器端真点击两个按钮，提示条数字与
`Content-Disposition` 文件名均正确（4,806 → 3,867；PVY 805），
再把两个 zip 从**运行中的实例**上取回来复检：契约 0 问题、两包除 date 列外全等。

### 环境提醒（本次踩到）

- `dev_tools/run_sanitized.py` **只存在于进化平台**（`D:\桌面\植物病毒进化分析平台`），
  所以 8900 是**进化平台**的实例；主平台在 **8765**（`python -u app.py --web`）。
  改主平台代码后要重启 8765 才生效（Flask 无 reload），别对着 8900 测主平台。
- 主平台的 `/api/explorer/query` 响应里 `rows_total` 是**全库行数**（`len(df)`），
  不是筛选后的条数 —— 别拿它当筛选结果计数（筛选后的量看 `stats.total` / `table`）。
- 前端 `filters()` 会把上一次的 `state.table` 一并回传，但服务端按筛选条件重跑，
  不吃这张表（导出尤其如此，见 `explorer.py` 模块头「已知差异 1」）。

## 2026-09-17 更正：比较基因组分析撤错了，已复原（导航层）

**谁改的**：2026-09-15 的「三块业务迁出」拆分（本次会话之外的操作，见
`docs/平台拆分记录_20260915.md`）。**为什么错**：拆分把「比较基因组分析」与
「进化动力学分析」当成同一件事捆着撤出导航，但两者是不同的业务对象 ——
比较基因组 = **不同病毒之间**的比较（跨科/属取参考序列 → MAFFT 比对 → 建树 / SDT
同一性）；进化动力学 = **同一个基因 / 序列集内部**的时间与地理信号
（重组 / 根到尾回归 / 系统地理）。09-16 恢复 `phylodyn` 时又用「覆盖」而不是
「并列」，把 `compare` 的一级导航位吞了，于是它落到「**有组、有卡、有路由，
但界面上到不了**」的状态。

**复原改动（三处，仍仅导航层）**：

| 文件 | 改动 |
|---|---|
| `Virus_Platform_Core/web/pages.py` | `NAV_GROUPS` 在 `annotate` 与 `phylodyn` **之间**插回 `compare` 组（5 张卡：`t-seqprep` / `/cds-export` / `t-align` / `t-treebuild` / `t-sdt`）；`_PATH_TO_GROUP` 恢复 `'/cds-export': 'compare'` |
| `webapp/templates/_nav.html` | `_nav_items` 在 `annotate` 之后、`phylodyn` 之前补回 `('/tools?g=compare', 'nav.g.compare', '比较基因组分析', 'compare')` |
| `webapp/templates/home.html` | 模块卡区补回「🧬 比较基因组分析」 |

`i18n.js` 的 `nav.g.compare` / `hm.entry.compare.d`（中英）**一直在**，无需补。
**未动 `phylodyn`** —— 09-16 的裁剪移植保持原样，两组现在并列。
`/explorer` **未恢复**进导航：病毒浏览器确实迁走了，与本次更正无关。

**关键澄清（推翻 09-16 文档里的一条说法）**：compare 的**后端在本站一直都好**，
从未迁走。`TOOL_REGISTRY`（`web/tools_api.py`）里的 `align` / `quicktree` / `sdt` /
`identity` / `structcmp` 与 `Virus_Platform_Core/{align_qc,sdt_exact,gb_collection,
cds_export,msa_view,phylo}.py` 都在本站；`tools.html` 里 4 张卡的 `section.card[id]`
（`:1322` / `:1430` / `:1546` / `:1649`）也都在。09-16 那条「比较基因组仍只在
进化平台」对本站不准确 —— 拆掉的只有导航入口。

**验证（全绿）**：
- 真浏览器点击链路：首页模块卡 → 一级导航（该项高亮）→ `/tools?g=compare` 组树 →
  4 个 hash 条目逐一点击，每次**只**该卡可见（`getComputedStyle(s).display`），
  `/cds-export` 可达；控制台仅 1 条 ArchaeopteryxJS 常规 log，无 pageerror。
- `tests/_it_compare.py`（SDT 精确矩阵 / NT+AA 同一性 / NJ+FastTree / GenBank 集合
  离线全流程）、`tests/_it_platform.py`（含「`NAV_GROUPS` 每个 id 必须是 `main > *`
  直接子元素」的结构断言 `:246`）、`tests/_it_align.py`、`tests/_it_msa.py`、
  `tests/_check_pages_console.py`（22 页）—— 全部通过。

**已定位但未修（与本次改动无关）**：`tests/_check_restore.py` 步 [3]「通用卡恢复」
会红 —— 它等 `#toolrun-contigs` 出现「上次运行」，而 `/api/tool/runs` 是
**按 mtime 倒序只取前 30 条**（`tools_api.py:190` 的 `names[:30]`）；最新一条
`contigs_*` 运行在 222 个运行目录里排第 **168** 位（2026-09-12），永远进不了这个窗口。
它是独立诊断脚本、**未登记进 `_run_all.py`**，不是门禁；修它属于改测试语义，
本次未动（步 [1][2] 走的正是 compare 组，均通过）。

**顺带发现并已定案（2026-09-17）**：主平台的 `/explorer` 也是「有页面、导航里没有」，
但**与 compare 性质不同** —— 它是**活的 7 页签浏览器**（trends / mutation / table /
primers / host / vector / profile），且这个「服务器版 7 页签」是**已提交**的工作
（`6c59773` 引入、`ef07b55` 修跨站深链 + 真点击验证）。用户明确决定**暂不恢复**
该入口、保持现状（只作深链访问）。所以：**compare 要进导航，explorer 不进** ——
别看形态一样就一起改。

## 2026-09-17 进化树 + 基因组叠加上线：补丁早就在，缺的是「接线」

**用户的观察是对的**：进化树查看器里确实做过「叶标签之外的蛋白结构与相似性连线」，
但 ③ 树查看的示例里没有。定位结论：能力在、**数据链路没有**。

**能力在哪**：`webapp/static/vendor/archaeopteryx/archaeopteryx.js` 里的自研补丁
「Comparative Synteny & Homology as tree overlays」（`~:5507` 起）—— 只画在矩形布局
（rectangular），三层叠加：属色带（`g.aptx-genus-bands`）、每叶基因轨道
（`g.aptx-genome`，叶标签右侧的箭头 gene box + 底部 `Genome Coordinates (bp)` 标尺）、
同源连线（`g.aptx-genome-links` 的半透明 ribbon，按蛋白同一性走 5 档色阶
30/55/75/90/100，右上角 `Protein Identity (%)` 色标）。数据入口是
`launchArchaeopteryx(box, name, nwk, {genomeOverlay: {...}})` 或句柄 `setGenomeOverlay()`；
面板里自动多出 `Genome Tracks` / `Genus Bands` / `Homology Links` / `Identity %` 开关。

**为什么示例没有**：`genomeOverlay` 在库和 `demo-synteny.html` **之外零引用** —— 主平台、
进化平台、`dist/` 全都没有生产方，两个渲染入口
（`webapp/static/app-compare.js: renderTreeTo`、`webapp/templates/results.html: treeRender`）
都只传 `{layout:'rectangular'}`。没有数据，补丁就永远不画。

**接线（本轮）**：

| 层 | 改动 |
|---|---|
| 数据约定 | `Virus_Platform_Core/web/common.py` 新增 `tree_overlay_for(nwk_path)`：叠加 JSON 放在**与树同目录同名**的位置（`tree.nwk` → `tree.overlay.json`）。缺失 / 损坏 / 非 dict / 空 tracks / 超 16MB **一律返回 None**，树照常渲染（叠加是可选产物，不能把树查看卡变成报错） |
| 后端 | 三个树入口都带上 `overlay` 字段：`web/tool_results.py` 的 `/api/tree/file`、`/api/tool/quicktree_data`（取运行目录里实际选中的那棵树）；`web/refs.py` 的 `/api/tree/data` |
| 前端 | `renderTreeTo(box, newick, opts)` 透传 `opts.genomeOverlay`；`treeToolLoadFile` / `treeToolView` / `quickTreeLoad` 与 `results.html treeLoad/treeRender` 都把 `d.overlay` 带上；`tvMeta` / `treeMeta` 多一行「基因组叠加：N 条轨道 / M 个基因 / K 条同源连线（按蛋白同一性着色）（仅矩形布局绘制）」 |
| 示例 | `tests/make_example_overlay.py` 产出 `examples/example_synteny_tree.nwk` + `examples/example_synteny_tree.overlay.json`（3 轨道 × 6 CDS，18 条蛋白同一性连线 68.6%~93.2%，NJ 树走平台自己的 `phylo._nj_distance_matrix` + `_neighbor_joining`）；③ 树查看新增按钮「✨ 示例·基因组叠加」（`data-i18n-title="tk.fillExampleOverlayTitle"`，两个按钮的 onclick 改成 `fillExample(...).then(treeToolLoadFile)`，否则示例根外置时相对路径会落空） |

**巡检**：`tests/_check_tree_overlay.py`（已登记 `_run_all.py`）—— 后端契约 +
**五类坏 overlay 的降级负控**（无同名 JSON / JSON 损坏 / JSON 是数组 / 空 tracks /
17MB 超限：都必须 overlay=null 且树仍 200）+ 真浏览器点示例数 SVG
（属色带 3 / 基因框 18 / 同源连线 18 / 标尺 1 / 色标 1，注意宽度 ≥7px 的 gene 画成
箭头 `polygon`、窄的才退化成 `rect`）+ **换纯树示例后四类元素必须清零**的负控 +
i18n 键 + 无 pageerror。已按纪律做**断电负控**：临时让 `tree_overlay_for` 直接
`return None` → 巡检 10 项变红，恢复后全绿。

**补丁留档（2026-09-17 已做）**：`webapp/static/vendor/` 在 `.gitignore:39` —— 这个补丁
**从未入库**，只存在于本机与该机 `dist/` 副本，而丢了之后平台只会「静默不画叠加」
（没有数据不算错，所以不会报错）。排查确认**本机已无上游原件**（源码、`dist/.../_internal/`、
进化平台三份 `archaeopteryx.js` md5 全同 `18e59d72…`），做不出可信 diff，
故按用户选择留**字节级副本**：`docs/vendor_patches/archaeopteryx/{archaeopteryx.js,
demo-synteny.html}` + `README.md`（md5、上游版本与 LGPL-3.0 说明、补丁三层与代码锚点、
恢复步骤与冒烟方法）。**`.gitignore` 未改** —— vendor 整体仍不纳管，只有这两个文件的
副本入库。

**本次不做（用户已定案）**：① `dev_tools/package.py` **暂不重打包** `dist/` —— 源码模式
8765 已带接线与新示例，打包版留待改动攒批再打；② 进化平台 8900 **暂不同步**这套接线
（那边有同一份补丁与 `demo-synteny.html`，接线时逐文件留 `.bak`，因为它无 git）。

**附带发现（2026-09-17 已修，见下一节）**：`examples/example_synteny_{A,B,C}.gb` 的 CDS 曾是
840bp @ 500bp 间隔（互相重叠），且第 6 个 CDS 结束在 3490 > `LOCUS` 声明的 3300。
生成器因此把标尺取 `max(len(seq), max(end))`，避免最右一个基因框被裁掉。

## 2026-09-17 「全部修复」：夹具自洽性 / 行尾 / 残留进程 / 一个被砸坏的表单标签

用户对上一轮报告的三个知情项 + 排查中新发现的缺陷做了统一修复。

**① 示例夹具 `example_synteny_{A,B,C}.gb` 改自洽（真实缺陷）**

| 缺陷 | 修法 |
|---|---|
| 6 个 CDS 按 500bp 间距排但每个 840bp → 相邻重叠 340bp（40%） | `tests/make_examples.py:make_synteny_trio()` 改为 840bp 编码 + 60bp 间隔、首尾各留 150bp |
| 第 6 个 CDS 收在 **3490** 而 `LOCUS`/序列只有 **3300bp** → 越界，按坐标切片被静默截断 | 基因组长度显式推导：`2*150 + 6*840 + 5*60 = ` **5640bp**（< 6.8kb 仍贴近 Potexvirus 口径） |
| 生成器用 `max(len(seq), max(end))` 把越界**兜住**了（缺陷被掩盖，没人再看得见） | `tests/make_example_overlay.py` 改为**断言**：CDS 越界直接报错并点名文件/基因 |
| 巡检不检查这条不变量 | `tests/_check_tree_overlay.py` 新增两条：每个 CDS 落在轨道长度内、同轨 CDS 互不重叠 |

重生成链（幂等）：`tests/make_examples.py:make_synteny_trio` → `tests/make_example_synteny_gb.py`
（ORIGIN 写真序列，5640bp / GC≈58% / EXAMPLE_SET 集合重建）→ `tests/make_example_overlay.py`
（树 + `overlay.json`）。**蛋白层面完全不变**：3 轨道 × 6 CDS、18 条连线 68.6%~93.2%（中位 75.0%）、
树拓扑 `(C,(A,B))` 照旧 —— 变的是坐标、轨道长度与（DNA 派生的）枝长。
`examples/results/synteny/summary.json` 里的 `length` 同步为 5640，并把指向已删工作目录
（`databases/examples/_work/synteny`）的绝对路径改为相对路径、补 `note` 说明该快照由已摘除的
LoVis4U 口径引擎产出（故不可原样重跑）；`ident` 数值来自蛋白、未受影响，保持原样。

**② 行尾归一化（9 个文件，只动行终止符）**
工作区里 9 个改动文件的行尾与 `HEAD` 约定不一致，`git diff` 被整文件噪音淹没（`app.js`
4058 行、`tools.html` 4489 行）。按各自 `HEAD` 约定重写：`app.js`/`assembly.py`/`db_migrate.py`/
`kv_engines.py`/`kv_stage.py`/`io_api.py` → LF；`kv_config.yaml`/`pages.py`/`tools.html` → CRLF
（后两个原本混入 45 / 2 行纯 LF）。改写前后以「按 LF 归一后内容逐字节一致」断言把关，原文件
备份在 `run/_eol_backup_20260917/`。归一后 `app.js` 的 diff 降到 110 行（只剩真实改动）。
**`.gitattributes` 的 `* -text` 未改** —— 纪律仍是「编辑时保持文件原有行尾」。

**③ 清理 6 个残留进程**（只留 8765 主平台 / 8900 进化平台）
`:8791` `:8792`（进化平台临时实例，父进程已退）、`:8818`（主平台临时实例）、
`:8931`（`http.server` 目录列表 —— 会暴露整个平台目录）、`:8766`（nohup 拉起的进化平台副本）、
PID 409424（`tests/_diag_pages.py` 跑了两天，占 `:18777`）。
**`:_pick_port` 的连带发现**：`app.py:_pick_port()` 按 `[8765, 8900, 8989, 9600, 8888, 5050, 5000]`
顺序取第一个可绑定端口，所以 8765/8900 被占用时，用户双击 `启动平台-网页.bat` 会静默落到
8989/9600 —— `:8989`、`:9600` 正是这样来的用户实例，**未动**（用户自己关）。

**④ 排查中新发现并修掉的真实 UI 缺陷：`tools.html` 定年卡的 `<label>` 被砸坏**
`tests/_it_platform.py` 的「导航项均有 `<main>` 直接子卡片」断言报 15 个 id 缺失
（`t-dsrna` / `t-pd*` / `t-phylogeo`）。根因不是行尾也非本次改动，而是 `tools.html:1795`：

```html
<label data-i18n="tk.rttMetaLbl"
<button ... onclick="pdPickDataset('rtt_input','rtt_meta','date')" ...>📦 从数据接入 run 选择</button> data-page-node-id="ew5WnP90Z0UaA8mckngklu">元数据 CSV / TSV（必选：seq_id + 日期列）</label>
```

按 `webapp/templates/tools.html.bak_20260916_geo移植前` 的原文可知是 09-16 那次
`data-page-node-id` 注入把「📦 从数据接入 run 选择」按钮**插进了 `<label>` 的属性区**：元数据标签
文本丢失、「选择」按钮退化成一串文字、页面上还漏出裸的 `data-page-node-id="…">`，并且未闭合的
`<label>` 让后面 15 张卡片在 DOM 里**嵌套错位**（点导航项即空白页 —— 正是该断言防的那类 bug）。
已按原文恢复 label + 按钮为兄弟节点。全模板扫描（「行内标签未闭合且下一行直接开新标签」）
**只此一处**；`_it_platform.py` 恢复 `PLATFORM CHECKS PASSED`，真浏览器复核：标签文本回来、
按钮存在可见可点、无游离文本、`main > #t-*` 抽查全 True、22 页无 JS 报错。

**⑤ `tools.html` 又被改回混行尾 → 再次归一（并更正归因）**
归一化完成后复核，`tools.html` 又变成 `worktree=MIXED`（CRLF 6318 + **恰好 4 行**纯 LF，
mtime 15:26:01）：`pdltt` 卡片里的 `<b id="pdlttTitle">` / `<span style="flex:1">` /
`<a id="pdlttDated">` / `<a id="pdlttPhyloPdf">`。已按 CRLF 再次归一（备份
`run/_eol_backup_20260917/tools.html.mixed_20260917b`），全仓 130 个改动文件复核 **0 处不一致**。

**归因更正（当日稍后查明）**：这 4 行不是"树外 IDE 插件"写的，而是**另一个会话在同一工作区
就地编辑**时新加的行 —— 当天有人在给 `t-pdltt` 卡加「真 TreeDater 引擎」，`pdlttTitle` /
`pdlttDated` / `pdlttPhyloPdf` 正是那批新元素（同一小时里 `run/` 下出现 7 个就地改写
`phylodyn_trees.py` 的补丁脚本，还有两个 `_check_phylodyn_example.py` 并行在跑）。
**教训**：看到 MIXED 先查时间线与"有没有别人的 python 在跑"，别抢在别人正在写的文件上做行尾
归一化（会互相覆盖）；09-16 那次 `data-page-node-id` 砸坏 `<label>` 的来源仍未证实。
真正的守卫仍是 `tests/_it_platform.py` 的「导航项均为 `<main>` 直接子元素」断言
—— 它上次就是这样抓出 `<label>` 被砸坏的。

**⑥ 示例结果清单里 26 处 `source` 指向已不存在的路径（改对）**
`examples/results/manifest.json` 的 `source` 全部写着 `databases/examples/…`，而该目录
2026-09-10 起已不存在（示例根搬到了 `examples/`）；生成器 `tests/make_example_results.py:179`
早就是 `'source': 'examples/'`，只有这份旧快照留在旧值。`source` 不经
`/api/examples` 暴露（只出 title / run / files / generated_at），故无用户可见影响，属
"死元数据指错路"。已按 26 处文本替换为 `examples/`（备份
`run/_eol_backup_20260917/manifest.json.before_source_fix`），JSON 仍可解析、条目数仍 26。
**未删任何示例结果目录** —— `identity` / `structcmp` / `synteny` 三套产物是"卡片已下架的
历史快照"（前者现为 SDT 卡模式二、中者降为 MSA 卡内部比对、后者是已摘除的 LoVis4U 引擎），
它们在「👁 示例结果总览」里仍能正常打开（后端按目录列产物，与有没有卡片无关），
所以**删不删是用户的决定**，本轮只报告不动手。

**验证**：`_check_tree_overlay.py`（含两条新断言，exit 0）、`_it_platform.py`（PASSED）、
`_check_pages_console.py`（22 页 ✔）、`_check_examples_manifest.py`（✔ 无硬性差异）、
`_audit_contract.py`（✔ 无硬性差异）。两条新守卫都做了**负控**：把 ORF6 改成 `4651..6490`
→ `make_example_overlay.py` 以「CDS 超出记录长度 5640 bp」退出 1；把 overlay 首轨基因改成
越界+重叠 → 巡检恰好多出 2 条 FAIL。两处还原后 md5 与改动前一致。**未改 `dist/`、未同步 8900、未提交 git。**

## 2026-09-18 引物设计 / dsRNA 设计从「病毒注释分析」拆出，新组「检测与防治」

**用户提出的判据**：「引物设计、dsRNA 设计是不是独立功能？从功能注释移除，归到检测与防治。」
核过一遍，**是独立的** —— 而且判据不止"看起来不像注释"，有三条硬理由：

1. **回答的问题不同**。注释组（`annotate`）回答的是「这条序列**是什么**」——
   ORF / 数据库比对 / 保守域 / 图谱，产物都是**对序列的描述**；
   引物设计回答「怎么把它**检出来**」（产物 = 引物对，服务 PCR / qPCR），
   dsRNA 设计回答「怎么把它**防住**」（产物 = RNAi 分子，服务防治）。
2. **产物的消费方不同**。注释产物喂给「结果中心 / 提交」这条链（读、看、投稿）；
   这两项的产物是**可下单的湿实验材料**，必须真做实验才能闭环，与注释不是一条链。
3. **依赖与参数集完全不相交**。注释靠病毒库 / CDD / nt·nr；引物靠 primer3 热力学
   与模板保守性；dsRNA 靠 dsRNAmax 选窗 + dsRIP 效价 + 脱靶面板 → T7 引物。
   三者没有共享输入，挂在一个组下只是历史偶然（旧页面把它们放在注释页侧栏）。

**改动（5 处，仅导航/文案层，**不动任何分析代码**）**：

| 文件 | 改动 |
|---|---|
| `Virus_Platform_Core/web/pages.py` | `NAV_GROUPS` 里 `annotate` 组从 7 项减到 5 项（移出 `/primer` 与 `t-dsrna`），紧随其后新建 `{'id':'detect','label':'检测与防治'}`（2 项）；`_PATH_TO_GROUP['/primer']` 由 `annotate` 改为 `detect` |
| `webapp/templates/_nav.html` | `_nav_items` 在 annotate 与 compare 之间插 `('/tools?g=detect','nav.g.detect','检测与防治','detect')` |
| `webapp/templates/home.html` | 「病毒注释分析」卡 desc 去掉「引物设计」；新增「🎯 检测与防治」卡 |
| `webapp/static/i18n.js` | 中英各加 `nav.g.detect`（检测与防治 / Detection & control）与 `hm.entry.detect.d`；`hm.entry.annotate.d` 中英同步删掉引物（**中英必须同改，`_audit_contract.py` 会卡 zh/en 键数一致**） |
| `tests/` | `_it_platform.py` 加 `/tools?g=detect` 到 PAGES + 4 条新断言；`_it_nav_responsive.py` 把写死的「13 个一级链接」改成按 `NAV_GROUPS` 派生 |

**`_PATH_TO_GROUP` 是这层最容易漏的一行**：`_group_nav.html` 用 `current_group`
拼 `/tools?g=<组id>#<卡id>`，漏改则 `/primer` 侧栏仍列注释组、dsRNA 入口指错组。
故新断言直接切 `/primer` 渲染出的 `class="group-topbar"` 区块来验，
而不是照抄 `NAV_GROUPS` 自证（那等于把测试绑在实现上）。

**验证**：`_it_platform.py`（PASSED，含 4 条新断言）、`_audit_contract.py`
（✔ 无硬性差异，i18n zh/en 各 2154 键一致）、`_check_pages_console.py`
（20 页真浏览器 ✔，含 `/primer`）、`_it_nav_responsive.py`（✔ 全部通过；**首次跑
`/pipeline 640px` 报过一条 linkRows=3，复跑即 2 行** —— 640px 属极窄窗口、
折行处在临界，与本次 +1 项无关：真浏览器里在同一页面内逐档摘掉 detect / compare /
phylodyn 三个组，15 / 14 / 13 / 12 个链接在 640px **一律折 2 行**，1920px 一律 1 行）。
线上实例已重启（8765）实测：首页导航 15 项、`/primer` 组顶栏「🏠 检测与防治」。
**未改 `dist/`、未同步 8900（该平台无 annotate 组，无需镜像）、未提交 git。**

**没做但可以做的**（留给用户定）：本组目前只有这 2 项。同属"检测"语义、现挂在别组的
还有 `t-variant`（病毒变异分析，现属「病毒定量与共识」，产物是 SNP/变异图 = 检测）
与 `/hostpredict`（现属「病毒识别和分类分析」）—— 要不要一并搬进来，等用户口径。


## 2026-09-30 导航：「公共数据下载」回归「数据资源」组 + GSA 反爬坑

### 导航调整（用户提出）

「公共数据下载」（/download）从「样本处理」组移回「数据资源」组，排在
「公共数据检索」正后。历史上它本来就在数据资源组，2026-09-06（本笔记
第十八节）应要求移入样本处理组首位，理由是工作流顺序"拿数据→建样品→
质控→分析"——但该衔接靠的是下载页文件行的跨模块按钮（🧾样品/🚫去宿主/
🧹质控/🦠鉴定），并不依赖导航分组。语义上样本处理组其余各卡都是"处理
本地已有文件"，只有下载是"从公共库获取"；且"公共数据检索/公共数据下载"
这对名字对仗的页面分居两组确实别扭。

改动两处（`Virus_Platform_Core/web/pages.py`）：NAV_GROUPS 的 download 卡
从 sample.items 移入 resource.items（检索之后），`_PATH_TO_GROUP` 的
`'/download': 'sample'` 改 `'resource'`（漏改则 /download 侧栏仍列样本组）。
`webapp/templates/download.html` 无需动（无写死组名）。

验证：test client 实测 /download 组顶栏=数据资源、resource 工作台 4 卡含
下载、sample 5 卡不含；`_it_platform.py` PASSED（含导航↔卡对应检查）。
`_it_nav_responsive.py` 报的 `/pipeline 640px` 折行与本次无关——stash 回退
旧版复跑同样失败（甚至多一项 1920px，属真浏览器测量波动），见 2026-09-18
节"640px 属极窄窗口、折行处在临界"的既有结论。

### GSA（CRR）下载解析失败根因：NGDC 反爬验证页（重要踩坑）

用户报告"CRR708330 解析失败（SSL EOF）"，排查出**两层**：

1. **SSL 瞬断**（04:41 那次直接原因）：`_http_get` 原先一次失败即判死。
   现已加指数退避重试 3 次（URLError/SSLError/TimeoutError/ConnectionError）。
2. **根本原因**：NGDC 于 **2026-09-17 上线 JS 反爬验证页**（nginx 静态页，
   ETag/Last-Modified 可证），对**无 Referer** 的脚本请求，
   `/gsa/browse/` 详情页一律返回"正在进行安全检查"页。GSA 解析流程
   （搜索页→详情页→提取下载链接）中搜索页能过、详情页被拦 → 症状是
   "NGDC 未找到公开下载文件"。实测放行条件：**带站内 Referer 即可**
   （`Referer: https://ngdc.cncb.ac.cn/gsa/search`）；页内 JS 算的
   `auth_token` Cookie（base64(salt+UTC日期)，salt 从页内 `var secret=`
   提取）作为兜底二次重放。**注意**：MMPV-RNA `gsa_sra.down.py` 的
   requests.Session 写法（无 Referer）2026-09-17 之后同样会被拦。
   修复在 `public_data.py::_http_get`：NGDC 域自动带 Referer，命中验证页
   按页内 salt 算 token 重放。下载域名 download.cncb.ac.cn **不需要**
   任何 Cookie，aria2c 直接下。

另：meta 页「④ 运行任务」卡片加了日志折叠块（懒加载
`/api/task/<id>/log`，展开状态跨 10s 轮询重绘保持）与「任务中心 →」
入口——此前已完成的历史任务只有状态没有回看入口；完整历史（含重启前
归档）一直在 /tasks 任务中心页。浏览器验收脚本 `tests/_check_meta_tasklog.py`。

## 2026-09-30 续：公共数据下载支持自定义输出位置

用户要求：下载任务可指定输出位置（不填 = 默认），下载后的批次卡能看到
数据在哪。

- **创建时**：`DownloadManager.create(..., out_dir=None)`（public_data.py）。
  out_dir 非空 → `check_path(in_platform=True)` 校验后**该目录即批次目录
  本身**（不再嵌 `<批次id>` 一层，用户指定哪里数据就在哪里）；目录自动
  创建，已存在且非空时批次日志提醒"同名文件将被覆盖/续传"。越界（平台外
  且未挂载输出根）抛 ValueError，`/api/dl/create` 捕获后 400 原因透出
  （想用平台外目录：设置页挂载输出根即可，check_path 的 write_roots 自动
  放行）。worker 落盘全走 `fe['out']`/`b['dir']`，自定义目录零适配。
- **可见性**：snapshot 原先**不带 dir**（前端根本拿不到输出目录），已补
  `'dir'` 字段；批次卡进度条下方显示 📁 完整路径 + 「📂 打开」按钮
  （openOutDir，app-jobs.js 已有）。
- **表单**：download.html 任务名/并发行下加「输出位置（可选）」输入 +
  📁 browseDir 浏览按钮；创建成功后清空（每批一目录心智，防下一批误混）。
  **download.html 原先没 include `_browse_dlg.html`** —— 页面一加
  browseDir 就会踩"死键"坑（见 2026-09-11 节），已补 include。
- **样品链路无需改**：文件行「🧾 样品」按钮预填的是文件实路径
  （sendTo→sessionStorage→samples 页 R1/R2），自定义位置自动适配；
  样品目录本身仍固定 `results/<样品名>/`（全平台样品发现/分析/结果的
  根基假设，不宜随单次下载漂移）。
- 验证：离线 create 4 行为（默认/自定义/越界 400/空白回落）+
  API test client 5 断言 + 浏览器 6 断言（tests/_check_dl_outdir.py）+
  `_it_dl_progress` / `_it_download_ui` / `_it_platform` 全 PASS。

## 2026-09-30 续②：工具卡输出区被旧任务回显覆盖（格式转换实测）

用户报告：格式转换连跑两次，第二次已完成（右上「任务完成」弹窗都是新
run），输出区的运行日志却还停在第一次（旧输入文件的 CMD）。

**根因**：`injectStageLogs`（app-jobs.js）对 `/api/tasks` 里**每个**名称
匹配「格式转换」的任务都执行 `toolrun-convert 容器.innerHTML = ...`，而
任务列表是**新→旧**序（TaskManager.order.insert(0)）——同卡旧任务数组
靠后、后写覆盖新任务。凡是能对同名任务跑多次的卡（转换/质控/鉴定…全部
TOOL_LABELS 卡、管道阶段日志 log-<key>、建库日志 buildrun-<key>）都有此
隐患，转换只是最先被撞见。

**修复**：injectStageLogs 里加 `claimed` Set，每个容器只吃数组里**第一个**
（=最新）匹配任务；三处注入（stage/toolrun/buildrun）统一。运行中优先于
已完成（顺序天然保证：新任务在前）。

**测试**：tests/_check_toolrun_latest.py —— 真浏览器连跑两次转换
（example_R1 → example_R2，seqkit 本地秒级），断言输出区最终是第二次的
run 名与日志；回退修复复跑 → 第二项 TimeoutException（=用户所见 bug 复现），
修复后 4 项全 PASS。测试基建备忘：
- 工作台卡片默认 display:none，需 `showModule('t-convert')` 激活；
- headless 全新 profile 必中 examples.js「首访自动弹示例结果」
  （.dlgmask z-index 320），JS 直调 `runTool('convert', btn)` 可绕开一切
  遮挡/点击拦截。

## 2026-09-30 续③：格式转换支持双端（R1/R2 分别转换输出）

用户需求：格式转换卡支持双端数据，双端则分别转换输出。

- **前端**（tools.html t-convert 卡）：`cv_input` 下新增可选 `cv_input2`
  （R2 文件）；`runTool('convert')` params 带 `input2`（trim，空则不传）。
- **后端**（tool_jobs.py `_tool_job_convert`）：
  - FASTQ→FASTA：`inputs = [inp, inp2] if inp2 else [inp]` 逐条独立跑
    seqkit fq2fa，产物按输入文件名各自命名（`GQMIX.R1.fq.gz →
    GQMIX.R1.fa.gz`），不合并；进度阶段 R1=0.3 / R2=0.7。
  - `.sra` 输入自含双端（sracha 解码天然出 _1/_2），带 input2 直接 400
    「.sra 文件自含双端信息，无需另选 R2」。
  - R2 后缀与 R1 不一致 → 400。`ctx.opt('input2')` 自带存在性校验。
- **验证**：API 层双端（R1/R2 两个 .fa.gz）/单端（行为不变）/
  sra+R2 400 / R2 格式不符 400 全过；浏览器 5 断言
  （tests/_check_convert_paired.py）含输出区双产物与两条命令日志；
  `_check_toolrun_latest`、`_it_platform` 回归通过。

## 2026-09-30 续④：样本处理四卡支持自定义输出目录（可选，不填默认）

用户需求：⚡一键前处理 / 格式转换 / 序列质控 / 宿主去除与序列提取 都加
输出目录设置（不填默认）。

**方案（工具运行框架层统一，四卡一次覆盖，全 registry 工具同样受益）**：
- 创建层（tools_api `/api/tool/run`）：body 顶层可选 `out_dir`，非空 →
  `check_path(in_platform=True)` 校验（平台内/挂载输出根内），产物写
  **<out_dir>/<run_name>/**（换父目录不换结构，多次运行不互混）；
  `tool_runs/<run_name>/` 仍创建并写 `out_link.json` 指针登记实际目录。
- 寻址层：`state.run_real_dir(run)` 统一解析（有指针→外移目录，无→原
  目录）；`/api/tool/runs` 列表、产物下载 `/tool_runs/<run>/<file>`、
  打开目录、删除运行、rdp/phylogeo/quicktree/contig 报告等 run 寻址全部
  换用。**没有指针登记，外移运行在历史区就是失联卡片**——这是本设计
  的关键不变量。
- 签名刷新：`/api/tool/runs` 的缓存签名基于 tool_runs/<run> 的 mtime，
  外移运行产物写在外部目录时列表会停在"只剩指针"——任务结束
  `_job_with_run` 里 `os.utime(base_dir)` touch 一下。
- 删除：清外移目录 + 指针目录两处；`_rmtree_retry`（Windows 下刚被
  send_file 下载过的文件句柄未释放时，gc + 0.3s 后重试一趟）。
- job 零改动：`ctx.run_dir` 已是外移后目录。
- 前端：四卡各加「输出目录（可选）」输入 + browseDir（cv_outdir /
  f_outdir / pc_outdir / hr_outdir）；runTool 的 convert/fastp/prepchain
  分支与 hostremoval 页提交带 `out_dir`。
- 验证：API 层 14 断言（tests/_check_tool_outdir.py：外移/指针/列表可见/
  下载可达/touch 刷新/删除双清 + 四工具各跑）；浏览器 7 断言
  （tests/_check_tool_outdir_ui.py：四卡字段渲染 + convert 端到端外部
  out_dir 显示 + 产物落盘）；`_check_toolrun_latest`、`_it_platform`
  回归通过。宿主库未建时 hostremoval 创建即 400（既有 db_ready 拦截，
  与本功能无关）。test client 教训：send_file 的 response 不 close 会
  占住文件句柄，让后续删目录在 Windows 上静默残留。

## 2026-10-21：主导航改版（顶部横条 → 左侧常驻栏）+ 首页分层架构总览

用户提供了一张智能辅助育种平台的架构图与工作台截图，要求借鉴其布局。
可迁移的三点里，落地了两点；第三点（右侧常驻文件面板）暂不做，理由见下。

**① 一级导航：顶部 14 项横条 → 左侧固定栏 + 顶部工具条**
- 旧方案的固有缺陷：14 个一级链接在 ≤1500px 就把横条撑成两行，
  `tests/_it_nav_responsive.py` 从 2026-09-11 起一直在守这个换行/溢出；
  用 CDP 实测 1920px 下内容宽 1674px 对盒宽 1675px，只差个位数像素，
  浏览器字体度量一抖动就掉行。竖排后该约束从根上消失。
- `_nav.html` 重写：`.sidenav`（品牌 / 一级链接 / 组内二级内联展开 /
  底部任务·结果·设置·手册）+ `.topbar`（☰ / 当前组标题 / 上下文胶囊 /
  语言 / 📂）。数据源仍是 `pages.py` 的 NAV_GROUPS，与 /tools 组落地、
  `_group_nav.html` 共用同一份定义。
- **结果中心组不在侧栏渲染**：底部已有直达入口（旧导航同口径），且该组
  没有 `nav.g.result` 的 i18n 键，渲染出来会显示成裸 key。
- 组内二级入口改为**内联展开**（`.sn-caret` 点击 + 桌面悬停），不再依赖
  悬停浮层菜单——触屏与窄窗口同样可用。展开态存 localStorage
  （`vp_nav_open`），跨页面保持。
- 窄屏 ≤1024px 收成抽屉：☰ 展开、遮罩/Esc/侧栏内 ✕ 收起、点链接自动关。
- 新增 `webapp/static/app-nav.js`（手风琴 / 抽屉 / 顶栏标题）；不依赖
  `$`/`t`/`esc`，纯原生 API，先于 app.js 执行也安全。
- `app.js injectGlobalCtx` 挂载点从 `.navlinks` 改为 `.navctx-host`；
  `i18n.js injectLangSwitch` 从 `.navlinks` 改为 `.navtail`。两处都保留
  旧选择器兜底，避免第三方模板漏改时胶囊/language 按钮凭空消失。
- **DOM 契约未变**：`.navlinks` / `.nav-drop > a[data-group]` /
  `#navTasksBadge` / `#navCtx*` / `#langBtn` / `.navtail` / `.iconbtn`
  全部保留，故 `_check_phylogeo_ui.py`、`_it_global_ctx_ui.py`、
  `_check_ai_dock_ui.py` 零改动通过。
- `tests/_it_nav_responsive.py` 按新语义重写：桌面档断言侧栏常驻/不遮挡
  主区/无横向溢出/链接数（NAV_GROUPS - result + 总览）/顶栏入口可见；
  窄屏档断言默认收起/主区全宽/☰ 与遮罩交互。旧口径的"横条单行"已无意义。

**② 首页新增「平台架构总览」（借鉴其分层架构图）**
- 位置：四张叙事卡之后、数据库模块之前。
- 结构：数据底座 → 分析引擎（五条业务线，底色加重作为视觉重心）→
  结果解读，右侧一列竖排「贯穿能力」（AI 助手 / 任务中心）。
- **每一块都是真实链接**，不是静态配图——点卡片直接进对应模块。
- 网格显式定位（`.arch-lane/.arch-arrow` 锁第 1 列，`.arch-side` 占第 2 列
  全部行）：不写死的话 auto-placement 会把第二条泳道塞进右侧列，整图错位
  （实测高度从 587px 涨到 1676px）。
- ≤1180px 时右侧贯穿列改为横排落到图下方，避免窄窗口"架构图看不全"。
- i18n 键 `hm.arch*` 中英成对（`_check_i18n.py` 1405 锚点全通过）。

**③ 未落地：右侧常驻文件/数据面板**
- 参考平台的右栏是"个人文件 / 公共文件"目录树，与其数据底座强绑定。
- 本平台现有 `/api/sample_files/<sample>` 已能列样品产物，但各模块的
  产物寻址是"运行目录 + tool_runs"而非统一文件树；做成全局右栏会与
  「结果中心」「任务中心」功能重叠，且挤压主区宽度。留待确有需求时再议。

**回归**：`_route_inventory` 246 条零差异 / `_check_js_syntax` 59 块通过 /
`_check_i18n` 1405 锚点全通过 / `_check_pages_console --isolated` 22 页
零 JS 错误 / `_it_nav_responsive` 全通过 / `_it_global_ctx_ui` 全通过 /
`_check_ai_dock_ui`、`_check_ai_global_ui` 全通过 / `_check_home_cover` 通过 /
`_check_phylogeo_ui` 通过。

### 2026-10-21 修复：分组点击后收不起来（悬停规则压过折叠态）

现象：点分组箭头折叠后，二级菜单仍然显示。

原因：CSS 里同时有 `.sn-group.open > .sn-sub{display:block}` 和
`.sn-group:hover > .sn-sub{display:block}`。鼠标必然停在分组上才点得到箭头，
此时 `:hover` 恒成立，把刚被 `.open` 移除的菜单又强制显示出来 —— 点击等于
没生效。:hover 那条是"桌面悬停也能看"的顺手便利，与点击开合语义直接冲突，
已删除；展开态现在只由 `.open` 类决定。

回归：`tests/_it_nav_responsive.py` 新增 `check_accordion()`，在 4 个页面 ×
1600px 下用**真实鼠标**先 `move_to_element`（让 :hover 激活）再点击，断言
1) 收起后 `.open` 为假且 `.sn-sub` 无渲染盒；2) 再点一次能重新展开。
必须用真实鼠标：只调 `element.click()` 不会让 :hover 成立，测不出这个缺陷。

### 2026-10-23 首页精简：移除静态封面叙事与重复架构图

按需求删除两块「只占版面、不产生操作」的内容：

- **四张叙事卡（Ⅰ 背景 / Ⅱ 方法 / Ⅲ 流程 / Ⅳ 意义）**：纯静态文案，
  与 hero 的科研定位句、下方模块入口信息重复。
- **平台架构总览（2026-10-21 新增）**：其节点与下方模块入口卡、
  左侧一级导航三处同源，属于第三份重复；删掉后模块入口仍由
  `hm.entry.*` 卡片承担。

同步清理：`hm.arch*` / `hm.cover.*` 中英各 22 / 8 条 i18n 键、
`.arch*` 与 `.cover-card/.cover-grid` CSS。

顺手修的四处首页缺陷：

- **深色主题白卡片**：`home.html` 内联样式用了 `var(--bg2)` / `var(--line)`，
  两个变量在 `app.css` 字典里根本不存在，深色主题下回落到 `#f8fafc` /
  `#cbd5e1` 浅色兜底值，快速开始与 MCP 控制台整块发白。现改用
  `--surface-2 / --line-200 / --input-bg / --blue-100` 等真实令牌。
- **回车不扫描**：输入框 placeholder 承诺「回车扫描」，但此前只有按钮绑定
  `QS.scan()`；补 `onkeydown` 后与文案一致。
- **首页漏掉进化动力学**：架构图删除后，模块入口卡是主区唯一入口，
  但卡片清单里没有 `phylodyn`（只有左栏有）。已按 `NAV_GROUPS` 顺序
  补卡与 `hm.entry.phylodyn.d` 中英键；现为 10 卡。
- **MCP 控制台移位**：原位置卡在「快速开始」与模块入口之间，把
  「粘贴路径 → 选模块」的主线截断；移到模块入口之后、最近样品之前。
  它是开发者调试面板，不该占据首屏主线位置。

`tests/_check_home_cover.py` 判据同步反转：断言 `.hero.cover` 存在、
`.cover-card` / `.arch` 不存在、封面文案非空、无横向溢出。
回归：`_check_i18n` 1376 锚点全通过 / `_check_js_syntax` 59 块通过 /
`_check_home_cover` 浅色+深色通过 / 1440·1280·1024·900·768 五档无横向滚动 /
回车触发 `scan_folder` 实测 1 次请求 / 首页控制台零 JS 错误。


### 2026-10-21 修复②：分组改为「整行标题」开关（箭头热区太小）

用户反馈：只有右侧小箭头能点，热区太小；应当点标题一次展开、再点收起。

改动：
- `_nav.html`：箭头从 `<button>` 降级为 `<span aria-hidden>`（纯装饰）；
  `aria-expanded` 移到 `a.sn-link` 上（整行才是真正的开关）。
- `app-nav.js`：点击整行标题 → `preventDefault()` + 切换 `.open`。
  href 仍保留组概览地址（JS 失效时退化为正常跳转）；有 JS 时不跳转，
  组概览走展开后的「🏠 组概览」项 —— 保持"点标题=开合"的单一语义。
- `app.css`：`.sn-caret{pointer-events:none}` 让点击穿透到整行。

**两个连带踩到的坑（都记在 CSS 注释里）**：
1. 箭头 `top:50%` 是相对 `.sn-group` 算的，而分组盒子 = 标题行 + 展开的
   二级菜单，高度随开合变化 —— 实测箭头偏出 104px、压到下一个分组上。
   改为相对分组顶部固定 `top:14px`（标题行 ≈38px 高，取半）。
2. 箭头**不能**挪进 `<a data-i18n>` 里：`applyI18n` 用 `el.innerHTML = v`
   覆盖，切语言时会把箭头整段抹掉。

回归：`tests/_it_nav_responsive.py` 的 `check_accordion()` 改为点 `a.sn-link`
（整行，不是 .sn-caret），并新增一条「点标题不跳转」断言 —— 若日后改回
只有箭头能点，或点标题误触发跳转，这里都会红。

### 2026-10-21 修复③：设置页侧栏点击跳转错位

用户报障：点左侧功能，跳转位置很乱、没有对齐。实测是**三个独立缺陷叠加**，
每一个单独看都像"跳转不准"：

1. **多列网格**：`.settings-grid` 是 `repeat(auto-fit, minmax(320px,1fr))`
   三列，而侧栏按文档顺序逐个列 —— 「输出/输入/数据库目录」三张卡在同一行，
   `top` 完全相同（实测都是 167），点哪个都滚到同一处，分不出目标。
2. **锚点被顶栏压住**：`.se-main section{scroll-margin-top:14px}`，而吸顶
   顶栏高 54px —— 滚动后卡片顶边贴在视口 0 处，标题被吃掉 40px
   （实测 `top=14` vs `topbarBottom=54`）。
3. **侧栏顺序 ≠ 文档顺序**：`sec-storage` / `sec-paths` 在 DOM 里排在
   「通用设置 / 数据库 / 分析默认参数」之后，侧栏却按 GROUPS 的书写顺序
   把它们归在「目录与存储」组前面 —— 点「平台目录」时按位置算出的当前
   区块是夹在中间的 `sec-defaults`，高亮跳到别人身上。

修复：
- 设置卡片改**单列流式**（`.se-main > .settings-grid{display:block}`），
  卡片内部的小网格 `#defaultsGrid` 不受影响。
- 新增 `--topbar-h: 54px`（app.css 根变量），`scroll-margin-top` 与
  sticky 侧栏 `top` 都引用它 —— 以后改顶栏高度只改一处。
- 侧栏条目改为**按 DOM 实际位置排序**（`docTop()` 累加 offsetParent），
  并把模板里 `sec-storage` / `sec-paths` 的物理位置前移到「通用设置」之前，
  让「目录与存储」组连续。侧栏顺序 = 视觉顺序 = 高亮判定顺序，三者同源。
- 高亮判定从 IntersectionObserver 改为**按滚动位置算**：观察器在多个区块
  同时进入观察带时后触发者覆盖前者（实测点「输入目录」高亮跳到
  「数据库目录」）。另加 700ms 点击锁 —— 平滑滚动途中 scroll 事件会连续
  触发，按位置算出的仍是上一屏的区块，会把刚点亮的目标覆盖掉。
- 点击时 `history.replaceState` 同步 `#hash`（不新增历史），带 hash 直接
  进入也按同一套偏移滚一次。

**过程中发现的两个坑（都写进了代码注释）**：
- `pickActive` 里"后面的都在下方就 break"是错的：从下往上跳时（sec-plugins
  → sec-paths）目标上方仍有多个区块越过顶栏，提前退出会停在上一个错项上。
  必须扫完全部再取最后一个。
- 首页 `home.html` 的 `VC` 对象被一行多余的 `};` 提前闭合（第 481 行），
  导致 `models: []` 变成非法语句、整个首页脚本 `Unexpected identifier`。
  该行是另一会话编辑首页时留下的残缺，已删。

回归：新增 `tests/_check_settings_nav.py`（已加入 `_run_all`）——13 个条目
逐个点击，断言落点偏移在 12~30px、高亮正确、各条目滚动位置互不相同、
侧栏顺序 == 文档顺序、`#hash` 深链也对齐、全程无 JS 报错。
`_check_pages_console` 23 页零报错、`_check_js_syntax`、`_check_i18n` 全通过。

### 2026-10-24 发布前 P0 批量修复（三档分发就绪）

按发布前体检结论逐项修复，全部有实测对照：

1. **分发说明文件损坏（打包脚本 bug）**：`package.py` 顶层说明的
   `'=' * 60 + '\n'` 写在隐式字符串拼接块里，Python 把它解析成「对前面
   整个字面量重复 60 遍」——产物 `程序与数据说明.txt` 首行重复 60 次。
   改为行列表 `'\n'.join()` 拼接，并顺手补 MCP / 插件两条说明。
2. **MCP 控制台在分发版静默失效**：`mcp_console._ensure` 固定
   `[python, mcp_server.py]`，冻结分发没有独立 python（exe 会把参数当网页
   参数）；`mcp_server.py` 不被 app.py 静态 import，PyInstaller 也扫不到，
   `_internal` 里连 mcp SDK 都没有——exe 版 MCP 控制台「连接」永远失败。
   修复三件套：app.py 顶部新增 `--mcp-stdio` 入口（在 Flask 启动前接管
   进程）；`mcp_console` 冻结模式改拉 `[exe, --mcp-stdio]`；spec 对 mcp
   `collect_all` + `hiddenimports += ['mcp_server']`。源码模式行为不变。
   `mcp_server.py` 同时随程序包/源码包各放一份。
3. **t-mirna 分发版必 500**：`databases/mirna_lib`（~13MB）从未入包，
   `tools_api._mirna_lib_rows` 缺库直接 abort。现随程序包拷贝，
   `dist_completeness.DEFAULT_PAIRS` 同步登记（+`__pycache__` 过滤）。
4. **数据库包混入 357MB 备份**：`_skip_engine_index` 的 `backup_*` 匹配
   不到带前导下划线的 `_backup_d17_20261019_*`（7 个目录），补
   `_backup_*` / `_bak_*` 两个模式。
5. **插件完整版打包入口**：`package.py --with-plugins` —— 32 个插件包
   （本地 `plugins-dist/` 优先，缺的按索引多地址下载并 sha256 校验）经
   `plugin_manager.install_zip(root=产物)` 预装进产物 `plugins/`。
   连带修 `install_zip` 的 `install_to` 基准 bug：原写死
   `PLATFORM_ROOT`，以 root=产物 调用时 6 个硬编码路径工具会写进**开发机
   仓库**而不是产物——改为 `dirname(root)`（默认调用两者等价）。
6. **测试两红 + 一挂死**：
   - `_it_platform`：detect 组期望集合补 `t-conan`（10-20 新增卡未同步断言）；
   - `_it_async_ui`：headless Edge 同会话连跑 4 个重页面在流程 4 原生崩溃
     （与 `_check_pages_console` 的 chromium 段错误同因）。改为**每流程一个
     全新会话**，首个会话失败仍 SKIP；实测批量满载下会话切换还会踩两个坑，
     一并修掉：流程 2 换会话后必须重新导航 /submit 并恢复页面级 CUR 状态
     （否则 no such element）；切换用 `_recycle()`——旧会话 quit 容错
     （批量下第 4 次切换 quit 会抛空消息 WebDriverException）、新会话 5 次
     递增退避重试，main() 的重试信息带异常类型（原来是空 Message，没法排）。
     VP_TEST_WAIT=90 批量口径连跑 3 次全绿；
   - `_check_pages_console`：父进程 `capture_output` 管道被 playwright 的
     chrome 孙进程继承，子进程超时被 kill 后管道写端仍被握着，
     `subprocess.run` 等 EOF **永久挂死**（实测卡 6 分钟无输出）。子进程
     输出改落文件 + `taskkill /F /T` 收树；18766 被占自动换临时端口；
     goto 错误截断 120→400 字符（原截断点恰好在 URL 端口前，误导排查）。
7. **依赖声明补齐**：requirements.txt 显式声明 mcp<2 / anyio / httpx /
   pydantic / uvicorn / cycler（此前仅靠传递依赖，干净克隆后 MCP 全断）；
   requirements.lock 以 3rd/python 绿色运行时 `pip freeze` 重新生成。
8. **RDP5 运行态文件出库**：`LastSave.rdp5` / `tempfile` 是用户会话快照，
   加入 .gitignore 并 `git rm --cached`（工具本体照常纳管）。
9. **第三方许可索引**：新增 `dev_tools/third_party_notices.py` 扫描
   3rd/bin+3rd/tools 随附许可文件与 3rd/python 全部 dist-info，生成
   `THIRD_PARTY_NOTICES.md`（33 工具 / 129 包）；README 分发形态章节改写
   为三档口径。

回归：见当日测试记录（目标测试 + 全量 quick + 产物自检）。

### 2026-10-24 补②：三档分层轴改定「插件当轴」+ 体积二轮瘦身

发布形态讨论定案（方案 B）：**"插件"从"完整版的标志"改为真正的分层轴**——

- **① 独立版** `--slim-tools`：3rd/tools 剔除 6 个大件进阶工具
  （srna-target / strawberry-perl / snpeff / jre-snpeff / snpgenie / datasets，
  ≈275MB）。它们的权威分发形态本就是 install_to 插件包，用户用到时
  「插件管理」一键补齐（模块运行前依赖自愈也会自动装）；
- **② +数据库版** = ① + 数据库包；
- **③ 完整版** = 全量 26 工具 + 数据库包，**不再预装 plugins/**——程序包
  自带全量 3rd 时插件层永远休眠（解析优先级 platform.json > 3rd 探测 > 插件），
  旧口径"完整版预装 32 包"实测是 1GB 纯冗余。`--with-plugins` 保留但降级为
  非默认的离线修复口径。
- 三档共用同一 exe/_internal，只有 3rd/tools 不同 → 档位切换用 `--no-exe`
  补拷，不必重跑 PyInstaller。

体积二轮瘦身（全档受益）：
1. **polars/_polars_runtime_32 剔除（-183MB）**：2026-09-11 它是 kv 引擎的
   必要成本（spec 旧注释）；f02736e 起 kv 全套迁 pandas（对拍
   tests/test_kv_pandas.py），全库 0 处 import——构建机残留的 polars 1.41
   连带把 pyd 拖进冻结包。spec excludes 落闸，注释同步改写；改动后须重跑
   test_kv_pandas + 已知病毒定量冒烟。
2. **OpenVirome SPA 双份收敛（-71MB）**：virome.py 优先级
   frontend/build(25MB 生产构建) > webapp/static/virome(71MB 开发副本)，
   冻结版此前两份都带、71MB 那份永远读不到。spec._tree 排除 static/virome，
   package.py 继续拷 frontend/build（两种运行模式的实际服务源都不变）。
