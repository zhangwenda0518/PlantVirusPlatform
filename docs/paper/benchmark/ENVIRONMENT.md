# ENVIRONMENT —— 基准评测的版本与环境记录

生成日期：2026-10-03（本节所有版本均在当前仓库工作区**实测采集**，采集命令见文末）。
适用对象：`BENCHMARK_zh.md`（Bari 小RNA + VIROMOCK 全面版，2026-09-16 数据快照）
与 `docs/paper/viromock_benchmark_manuscript_EN.md`（VIROMOCK 稿件）。
原始实验跑于 2026-09-15/16，**工具版本以外置 `3rd/` 目录内二进制的现值为准**
（目录随仓库纳管，未升级）；若曾升级，性能数字需重比。

## 1. 外部工具（3rd/ 预编译，实测 `--version`）

| 工具 | 版本 | 路径 | 采集方式 |
|---|---|---|---|
| salmon | **2.7.0** | `3rd/tools/salmon2/salmon.exe` | `salmon --version` |
| minimap2 | **2.31-r1302** | `3rd/tools/minimap2/minimap2.exe` | `minimap2 --version` |
| minibwa | **0.7-r424-dirty** | `3rd/bin/minibwa.exe` | `minibwa version`（`--version` 不支持） |
| SPAdes | **v4.3.0-dev** | `3rd/tools/SPAdes-Windows/`（Python 前端，Windows 版） | `python spades.py --version` |
| kunpeng | **v0.7.12**（文件名标注） | `3rd/bin/kun_peng-v0.7.12-x86_64-pc-windows-msvc.exe` | 文件名；二进制无 `--version` |
| DIAMOND | **2.2.6** | `3rd/tools/diamond/diamond.exe` | `diamond --version` |
| blastn (BLAST+) | **2.15.0+**（build 2023-10-19） | `3rd/tools/Blast/bin/blastn.exe` | `blastn -version` |
| samtools | **1.24**（htslib 1.24） | `3rd/tools/samtools/bin/samtools.exe` | `samtools --version` |
| pandepth | 版本号未暴露（`--version` 不支持，`--help` 无版本行；文件日期 2026-09-09） | `3rd/tools/pandepth/pandepth.exe` | 见左；如需精确版本请查上游发布说明对照构建 |
| TreeTime | **0.12.1** | 站点包 `C:\Python312\Lib\site-packages\treetime`（`treetime/__init__.py` 内 `version = '0.12.1'`；无 dist-info，属手工/vendored 安装） | 读源码常量 |

## 2. 运行时与硬件

| 项 | 值 | 采集方式 |
|---|---|---|
| OS | Windows 11 专业工作站版，10.0.26200 | `Win32_OperatingSystem` |
| Python | 3.12.10（64-bit MSC v.1943） | `python -VV` |
| CPU | 12th Gen Intel Core i7-12700K（12 核 / 20 逻辑处理器） | `Win32_Processor` |
| 内存 | 68,475,338,752 B ≈ **63.8 GiB** | `Win32_ComputerSystem` |
| Python 依赖（关键） | matplotlib 3.11.0、polars 1.41.0、psutil 7.2.2 | `importlib.metadata` |

## 3. 线程配置与平台默认的偏离（重要）

评测脚本统一按 **4 线程**跑（`run_one.py` 传 `--align-threads 4`；SPAdes
`threads=4`；samtools 4）：

- **平台默认 `--align-threads 8`**（`known_virus_suite.py:153` argparse 默认；
  `kv_config.yaml` 同值）→ 基准的绝对耗时偏保守（线程数低于默认）；
  §4 表格与 fig_b5/fig_b6 的墙钟时间均应理解为 4 线程口径。
- 全链路评测（`run_chain_bench.py`）同用 4 线程；EN 稿件（VIROMOCK 篇）用
  19 线程口径（`_bench_kvs2/`，未随本仓库归档，见该稿数据可用性节）——
  **两套口径不可直接比较**。

## 4. 参考库版本快照（基准数字的绑定物）

| 库 | 用途 | 构建时间 | 记录数 | 哈希/清单 |
|---|---|---|---|---|
| `virusref_db/kv_index` | **基准默认库**（BENCHMARK_zh.md 全部数字 + EN 稿物种级对照） | 2026-09-19 快照（含 TSWV 注入记录，见下） | 未在库内留 manifest；EN 稿 §2.2 记 **8,464 参考** | **未随仓库归档、无原件哈希**（`databases/` 不入 git；2026-09-29 起被 plant_ref 取代，目录已移除） |
| `virusref_db/plant_ref` | 2026-09-29 起平台默认库（取代 kv_index） | manifest `built_at=2026-09-19`，12,159 条 | 12,159 | manifest sha256 `d180f726…5261`；reference.fasta sha256 `3cec4271…9c10`；ref_info sha256 `d9cb8949…262f` |
| `virusref_db/viromock_kv` | EN 稿任务定制库（株级，61 条） | manifest `built_at=2026-09-15 21:39:23` | 61 条（+表头） | manifest sha256 `f4ba53c2…66c2`；reference.fasta sha256 `9d1c5f35…5816`（516,344 B） |

关于 kv_index 的可追溯信息：`plant_ref/manifest.json` 的 notes 字段记载了
2026-09-19 hotfix（accession 版本号归一 + **注入 TSWV L/M/S**：
NC_002052.1 / NC_002050.1 / NC_002051.1，*Orthotospovirus tomatomaculae*，
taxid 3052585）——BENCHMARK_zh.md "参考库版本"节引用的即此记录。

**复现提示**：BENCHMARK_zh.md 的所有数字绑定 kv_index 快照（非 plant_ref）。
重跑需先重建等价库（物种级 8,464 参考、含 TSWV 注入），或显式指定
`VP_BENCH_DB`（详见 BENCHMARK_zh.md §参考库版本 与各脚本 `--db-dir`）。

## 5. 采集命令（可复核）

```bash
3rd/tools/salmon2/salmon.exe --version
3rd/bin/minibwa.exe version
3rd/tools/minimap2/minimap2.exe --version
python 3rd/tools/SPAdes-Windows/bin/spades.py --version
3rd/tools/diamond/diamond.exe --version
3rd/tools/Blast/bin/blastn.exe -version
3rd/tools/samtools/bin/samtools.exe --version
python -c "import treetime,pprint;print(open(treetime.__file__).read()[:200])"
powershell -NoProfile -Command "Get-CimInstance Win32_Processor | Select Name,NumberOfCores,NumberOfLogicalProcessors; Get-CimInstance Win32_ComputerSystem | Select TotalPhysicalMemory"
```
