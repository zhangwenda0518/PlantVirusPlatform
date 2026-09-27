# Quick-Start · 快速开始

> 返回 [Home](Home) ｜ 上一页：无 ｜ 下一页：[Architecture](Architecture)

本文覆盖：环境要求、依赖安装、两种启动模式、首次数据库配置、跑通第一个样品。

---

## 环境要求

| 项目 | 要求 |
|---|---|
| 操作系统 | Windows 10 / 11（平台按 Windows 原生设计，含中文路径兼容） |
| Python | 3.12+（打包版内置运行时，无需安装） |
| 内存 | 建议 ≥16GB（宿主库构建、SPAdes、DIAMOND 为内存大户） |
| 磁盘 | 平台本体 ~1.5GB；数据库包 ~3.4GB；宿主库 ~1.5GB/物种；样品结果另计 |
| 网络 | 仅在下载/检索/按需选参时需要（全部走 NCBI/ENA/NGDC/ICTV 官方域名白名单） |

## 安装依赖

```bat
git clone https://github.com/zhangwenda0518/VirusPlatform.git
cd VirusPlatform
python -m pip install -r requirements.txt
```

**外部工具**已内置在平台目录 `3rd/`（bin/ 与 tools/，kunpeng、seqkit、crabz、Blast、mafft、FastTree、iQtree、trimAl、Gblocks、diamond、mmseqs、fastp、salmon、aria2c、sracha 等），**2026-09-28 起全量随 git 仓库分发，克隆即用**；仅两个超 100MB 的文件走 [Release（external-tools）](https://github.com/zhangwenda0518/VirusPlatform/releases/tag/external-tools)，下载后放回对应路径即可（`3rd/tools/rdp5/3seqTable`、`3rd/python/Lib/site-packages/_polars_runtime_32/_polars_runtime.pyd`，各目录内有同名 README 说明）。也可自行放置后用 `platform.json` 的 `tools` 段指定路径。

**需要单独安装的两项**：

- **SPAdes**：安装 `SPAdes-Windows-4.3.0-dev-Setup.exe` 并确保 `spades.bat` 在 PATH；
- **gbdraw**（基因组图首选引擎）：`pip install git+https://github.com/satoshikawato/gbdraw.git`。未安装时自动回退纯 Python 引擎 dna_features_viewer（随 requirements.txt 安装），功能不缺、风格不同。

**验证安装**：

```bat
python main.py selfcheck    :: 模块导入 / 外部工具 / 数据库 / 磁盘 / 依赖 五段结构化体检
python main.py tools        :: 只看工具探测状态
```

## 启动（两种模式）

| 模式 | 命令 | 说明 |
|---|---|---|
| 独立桌面窗口 | `python app.py --gui` 或双击 `启动平台-桌面窗口.bat` | pywebview/WebView2 壳住本地页面，无地址栏；缺 WebView2 自动回退浏览器 |
| 网页模式 | `python app.py --web` 或双击 `启动平台-网页.bat` | 默认浏览器打开；地址打印在控制台，可复制到其他浏览器 |
| 自动（默认） | `python app.py` | 优先桌面窗口，不可用回退浏览器 |

- 端口从 **8765** 起依次探测（8765→8900→8989→…），被占用自动顺延，**以控制台打印的地址为准**；
- 只监听 `127.0.0.1`，不对外网开放；
- 桌面窗口模式关闭窗口即退出；网页模式关闭控制台窗口即退出；
- 兼容旧写法：环境变量 `VP_GUI=browser|gui` 等价于 `--web|--gui`。

辅助脚本：`环境自检.bat`（selfcheck + 内存上限检查）、`挂载数据目录.bat`（数据库外置挂载）、`启动拖拽代理.bat`（浏览器拖拽文件支持）。

## 首次配置（仅一次）

进入 **数据库** 页按顺序完成：

1. **下载 NCBI Taxonomy**（~57MB，`python main.py init-taxonomy`）——nodes/names/merged.dmp，宿主判定与谱系解析的基础；
2. **构建宿主库**：选宿主基因组 FASTA（NCBI 下载、自行组装均可，`python main.py host-genome --species "..."` 可直接下载）+ 填宿主 **NCBI TaxID**（如 4081 番茄、112863 枸杞）→ 构建。TaxID 建库前自动校验（不存在/已合并直接提示）；序列名自动规范化；大基因组自动切片（见 [Databases](Databases)）；
   - 输出默认 `host-db/<TaxID>_<源目录名>_host_db/`，建完即设为「当前宿主库」；多宿主库可共存，列表里一键「设为当前」；
   - 病毒库**无需构建**：平台已预置 `databases/virusref_db`（鉴定+定量+BLAST）与 `databases/kunpeng_db`（分类库）。

## 跑第一个样品

**图形界面**：分析流程页 → 新建样品（名称 + R1/R2；单端数据 R2 留空）→ 点 **▶ 依次运行剩余步骤**。
卡片按功能模块分组（预处理→病毒鉴定→组装→宿主预测→下游分析→报告），每步显示摘要（宿主占比/病毒种数/contigs/ORF 数等）与产物文件。

> 💡 首次建议把「子采样」设为 **100000**（10 万对 reads），几分钟验证全流程。

**命令行等价**：

```bat
python main.py analyze --r1 R1.fastq.gz --r2 R2.fastq.gz --sample DEMO
python main.py analyze --r1 ... --sample T1 --stages virus,orf --subsample 200000   :: 只跑部分阶段
python main.py report --sample DEMO                                                 :: 重新生成报告
```

阶段可自由组合：`subsample,fastp,fq2fa,host,kvsuite,assembly,hostana,orf,orfa,phylo,primer,gbdraw,report`；每阶段有 `.done` 断点标记，中断后重跑自动跳过已完成阶段（`--force` 强制重跑）。

## 下一步

- 了解平台如何组织代码与数据 → [Architecture](Architecture)
- 了解 14 步级联与进度模型 → [Pipeline-Overview](Pipeline-Overview)
- 各阶段细节 → 左侧导航「主分析管道」各页
