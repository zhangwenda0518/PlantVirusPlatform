# Settings-Ops · 设置 · 任务中心 · 运维

> 返回 [Home](Home) ｜ 上一页：[Databases](Databases) ｜ 下一页：[CLI-Reference](CLI-Reference)

对应引擎：`selfcheck.py`、`tool_runs_admin.py`、`db_migrate.py`、`web/state.py`（设置与任务状态）。

---

## 设置中心

| 区块 | 内容 |
|---|---|
| 界面语言 | 🌐 一键切换中/英（导航栏或设置页）；前端界面与后端阶段名/卡片摘要同步切换；偏好写入 platform.json 重开保持 |
| 全局线程数 | 影响所有阶段默认线程 |
| NCBI 邮箱 | Entrez 工具必填 |
| 分析默认参数 | **17 项**：confidence / subsample / 组装模式 / 组装输入 / 内存 / 最小 contig / 最小 ORF / 参考数 / 建树工具 / 抽样策略 / 引物模式 / 特异性检查 / 出图上限 / 绘图引擎 / 强制重跑 / max_heavy_tasks / max_light_tasks——管道页与 CLI 未显式传参时自动采用 |
| 存储与磁盘 | 盘剩余容量条 + 各目录占用明细（databases 按子库拆分）；后台任务 + 10 分钟缓存；**平台盘剩余 <20GB 全站顶部黄色水位横幅** |
| 数据库目录 | 外置数据库根对接（即时反馈宿主库/病毒库/Taxonomy 就绪状态；切换后相关路径常量同步刷新，无需重启） |
| 示例数据目录 | ✨示例按钮的数据根（自动探测 + 手动指定） |

## 任务中心

- 全平台统一任务列表（样品分析 / 工具 / 下载 / 比较 / LOGAN / 建库 / 提交），导航栏入口带**运行中任务数角标**；
- 任务卡：实时全局进度（加权）+ 已运行时长 + **预计剩余时间**（自学习模型 `run/logs/stage_perf.json`）+ 排队信息；
- **结果预览面板**：任务结束自动展开——样品任务显示各阶段关键数字 +「打开报告」；工具任务显示统计值 + 产物清单（可直接下载）；弹出自消失 toast；页面在后台时发系统通知；
- 服务失联自动显示红色诊断横幅，恢复自动消失。

## 环境自检

```bat
python main.py selfcheck
```

五段结构化结果（模块导入 / 外部工具 / 数据库 / 磁盘 / 依赖），同时供两处使用：Web API `/api/selfcheck`（设置页卡片 + 总览页环境状态条）与 CLI（存在未就绪项时以非零码退出，供脚本判断）。`mem_check.py` 单独验证进程可用内存上限（识别受限沙箱）。

## 输出目录管理（tool_runs）

「固定目录 + 动态运行区」双层结构：

```
run/tool_runs/
├─ _archive/   历史运行归档（_archive/<工具>/<运行>/），长期留存
├─ _tmp/       冒烟测试与一次性实验（随时可整体清空）
├─ _scripts/   运维/补丁脚本
├─ _reports/   跨运行汇总产物
└─ <工具>_<YYYYMMDD_HHMMSS>/   动态活动区（新运行都在这层）
```

`_` 前缀 = 固定目录：不出现在 GUI「运行目录」列表、永不被 archive/clean 当作运行处理；无时间戳命名的目录视为杂项。

```bat
python main.py tool-runs status                        :: 活动区/固定区统计
python main.py tool-runs organize [--dry-run]          :: 杂项自动归位
python main.py tool-runs archive --before 20260901     :: 归档旧运行
python main.py tool-runs clean --active-days 30 --archive-days 90 --dry-run
```

建议节奏：每月 `archive --older-days 30`；磁盘紧张时 `clean --active-days 30 --archive-days 90`。删除不可逆，**务必先 `--dry-run` 预演**。

## 数据库迁移（db_migrate）

软件与数据库分离部署：

```bat
python main.py db-migrate --to D:\库目录              :: 复制后切换（源保留）
python main.py db-migrate --to D:\库目录 --mode move  :: 复制校验后删源
python main.py db-migrate --to D:\库目录 --dry-run    :: 只预检
python main.py db-migrate --to D:\库目录 --check      :: 只校验目标已有库
```

- 迁移安全：robocopy（多线程+断点续传）→ 文件数与总字节校验 → **通过后才切配置**；任一失败配置不动、源原样，可断点重跑；
- 目标须在平台目录外，目标盘剩余 ≥ 源体积 + 10GB。

## 打包分发（已停用）

> **2026-10-11 定案：平台不再打包分发**（`VirusPlatform.exe` 安装包路线停用），以源码仓库为唯一发布形态——克隆、阅读、修改、定制都在源码上完成，也更契合智能体（AI agent）驱动的工作方式。
> `dev_tools/package.py` 仍保留在仓库中，供有自分发需求的用户自行使用；数据库外置迁移（上文）照常可用。

## 品牌资产

```bat
python tests/make_brand.py            :: 徽标 SVG → logo-512/192/96.png + favicon.ico
python tests/_check_home_cover.py     :: 真浏览器回归（需服务在跑）
```

首页封面为文字封面（平台名 + 科研定位 + 「背景/方法/流程/意义」四张叙事卡，i18n 键驱动双语）。

## 相关

- 进度与并发模型 → [Pipeline-Overview](Pipeline-Overview)
- 数据库内容 → [Databases](Databases)
