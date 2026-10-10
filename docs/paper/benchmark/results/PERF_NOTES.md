# PERF_NOTES —— results/perf.csv 字段口径说明与已知缺陷

生成日期：2026-10-03（配合 docs/audit_20261003_修复对照表 的论文基准修复项）。
适用数据：本目录 `results.json` / `perf.csv` / `metrics_summary.csv`，均为
**2026-09 评测快照**（库 = 当时的 kv_index，见 `../BENCHMARK_zh.md` 参考库版本节）。

## 1. salmon 行 `total` 字段语义错误（已知缺陷，未重跑修正）

### 现象
`perf.csv` 中 23 行的 `mapped > total`（比对率 >100%），全部为 **salmon 行**
（外加 `vm_d12_bwa` 一行差 2 条，见 §3）：

```
bari_s3_50k,11464,10894   bari_s3_k15,29067,27257   bari_s9_50k,3272,2319
bari_s9_k15,16617,11794   chain_s3,32984,29988      chain_s9,19544,13948
vm_d1,116289,109316       vm_d2,158910,150408       vm_d3,1252442,1251398
vm_d4,1540508,1474283     vm_d5,2550972,934506      vm_d6,2190374,842341
vm_d7,7983,7918           vm_d8,56151,53821         vm_d9,117506,117422
vm_d10,1279325,1213362    vm_d11,69988,69641        vm_d12,270648,95902
vm_d13,80723,80571        vm_d14,175830,93578       vm_d15,77778,75741
vm_d18,92248,84920        vm_d12_bwa,96446,96444
```

### 根因（代码级，2026-10-03 核实）
评测当时（2026-09-15/16 跑的 run_one.py → 平台 kv_engines.py）salmon 引擎的
两列来自不同口径：

- `mapped` = 各参考 `Mapped_Reads` 求和（`count_mapped_by_ref`）——**多映射
  记录按参考重复计数**；
- `total` = `samtools view -F 0x4 -F 0x100 -c sorted_bam`
  （`BamCoverageBackend.count_records`，参见 kv_engines.py:777）——这是
  **"已比对记录"数**，不是输入片段数（名字叫 total，实际是 mapped-only 的
  分母）。

两个口径叠加导致 `mapped` 可能超过 `total`（多映射越多的集越明显，如
D12 95902→270648、D5 934506→2550972）。**`mapped/total` 比值不能当作比对率
使用**；本仓库文档未使用该比值（性能对比只用 wall_s / rss_mb）。

### 正确口径
- salmon：输入片段数取 `quant/aux_info/meta_info.json` 的 `num_processed`，
  比对数取 `num_mapped`（salmon 自己写的权威统计，不依赖 BAM flag）；
- minibwa：`total` 用输入 reads 数（或 `samtools view -c -F 0x900` 全主记录，
  含未比对），`mapped` 用去重复的多映射分摊计数（现平台为每参考计数求和，
  多映射仍会略偏大，见 §3）。

### 平台侧修复状态
2026-10 平台 `kv_engines.py` 已改为 `_meta_counts()`（优先读
`meta_info.json` 的 num_processed/num_mapped，BAM 兜底路径保留但记 warning；
代码注释明确写了"多映射记录……比对率虚高甚至 >100%"）。**本基准归档未重跑**，
`perf.csv` / `results.json` 中的旧字段保持原样，故本说明随数据发布。

### 影响与行动
- 影响：仅 `perf.csv` 的 `mapped`/`total` 两列（及其派生比值）。
  `wall_s`、`rss_mb`、`metrics_summary.csv`（基于确诊表，不经 mapped/total）
  不受影响。
- 重跑修正方式：重建 2026-09 等价库（见 `../BENCHMARK_zh.md` 参考库版本节）
  后重跑 `run_bench.py`，新记录将写入修复后的 meta_info 口径。

## 2. minibwa 行 `mapped` 的轻微重复计数

`vm_d12_bwa` 出现 `mapped=96446 > total=96444`（差 2 条）。minibwa 的
`total` 取 `-F 0x900` 主记录数（含未比对，本身口径正确），`mapped` 为逐参考
比对记录求和——同一 read 若同时比对到两条参考会被计两次。差值极小（该集 2 条），
但说明 minibwa 的 `mapped` 也**不是严格的 read 去重计数**，跨引擎比较映射率时
应使用去重后的 read 数。

## 3. 与本说明配套的文档标注

`BENCHMARK_zh.md` §4（运行开销）表下方已加脚注指向本文件；
`results/perf.csv` 仅增加"备注"列（对照行口径），未改动历史数字。
