# 本地时间与地理推断 · 运行报告

生成时间：2026-10-06 16:27:02

> 本报告由平台「本地时间与地理推断流水线」（轨道 A）自动生成。**全程本机完成，未跑 MCMC**；凡涉及后验显著性的结论仍需 BEAST（轨道 B）。

## 一、阶段状态

| 阶段 | 名称 | 状态 |
|---|---|---|
| A0 | 数据准备 | ✅ 完成 |
| A1 | 时间推断 | ✅ 完成 |
| A2 | 地理推断 | ✅ 完成 |
| A3 | 群体动态 | 未运行 |
| A4 | 祖先序列 | 未运行 |
| A5 | 汇总溯源 | ✅ 完成（本报告） |

## 二、定年结果（**每条都带方法名，不可互相比较**）

| 方法 | 钟模型 | 速率（每位点每年） | R²（时钟拟合, TreeTime 口径） | tMRCA | CoV | 备注 |
|---|---|---|---|---|---|---|
| TreeTime | 严格钟（默认） | 0.003241 | 0.44 | 1978.43 |  | ML 分子钟；r² 为根到尾回归（与 DRT 不是一回事） |

时间树接口产物来源：**TreeTime**（A2/A3 用的就是这一棵）。

时间信号检验（DRT，200 次有效置换）：R²（DRT 最大根位）=0.6171，经验 p=0.005，判据 **通过**。（此 R² 与上表"时钟拟合"列的 R² 不是同一个量，不可互相比较。）

## 三、地理迁移（两条口径**并列**）

- TreeTime mugration（**ML 点估计**）：迁移 **4** 次，走廊 3 条，事件 4 条
- 平台 Fitch 简约（同树）：迁移 **4** 次

> ⚠️ 两列**不是同一个统计量**，只能看趋势一致性：**不可相加、不可平均、不可排名**。实测 RSV 209 总数 70 vs 71（差 1.4%），但逐走廊方向差别很大（如 NCS→KOR 24 vs 2）—— 方向性结论务必谨慎。

> ⚠️ mugration 的 `confidence.csv` 是 **ML 边际概率，不是贝叶斯后验**（无先验、无拓扑不确定性），与 BEAST 后验概率**不可比大小**。

- 枝长口径对拍（换 `divergence_tree.nwk` 再跑一次 mugration）：迁移 4 次，相对主口径 0.0%

- 采样偏倚：最多区划占 0.4，最少占 0.2667；已按系数 0.992767 校正

## 六、告警（逐条看，别跳过）

- **[A1]** rooting="optimize"：零分布同样逐轮取最大根位 R²（与参照发布版管道 temporal_signal.py 一致）。真实数据"最大根位 R²" = 0.6171，而树自带根位只有 0.6132 —— 两者差距就是"根位搜索"贡献的量，**单看 R² 数值会高估信号**。
- **[A2]** 2/14 个内部节点的ML 最大边际概率 < 0.7（最低 0.43）—— 涉及这些节点的迁移方向不可靠，做方向性结论前先看置信分布

## 八、产物清单（含校验和）

| 阶段 | 产物 | 大小 | sha256(16) |
|---|---|---|---|
| A1 | `clock/clock_stats.tsv` | 253 | `8bd7217671b8c268` |
| A1 | `clock/divergence_tree.nwk` | 3577 | `a22ed856b58da5b7` |
| A1 | `clock/clock_qc/drt.json` | 8725 | `e5b8b05a866763cb` |
| A1 | `clock/clock_qc/rtt_scatter.tsv` | 2509 | `39b312fd3d180ca0` |
| A1 | `clock/timetree.nwk` | 3521 | `329576de507913e5` |
| A1 | `clock/treetime/ancestral_sequences.fasta` | 36573 | `111aca98e8a5168b` |
| A1 | `clock/treetime/auspice_tree.json` | 11846 | `7dbed31f71e07d23` |
| A1 | `clock/treetime/divergence_tree.nexus` | 4045 | `119f83ea5fd82e25` |
| A1 | `clock/treetime/molecular_clock.txt` | 60 | `8d6eaa26f112dc48` |
| A1 | `clock/treetime/root_to_tip_regression.pdf` | 16028 | `9be830d084e508ee` |
| A1 | `clock/treetime/timetree.nexus` | 3990 | `d6ff91a17ca20f5c` |
| A1 | `clock/treetime/timetree.nwk` | 3521 | `329576de507913e5` |
| A1 | `clock/treetime/timetree.pdf` | 13025 | `459e14f012e636e3` |
| A2 | `phylogeo/corridors.tsv` | 109 | `4f794af77a9a2817` |
| A2 | `phylogeo/events_year.tsv` | 488 | `92ba4ed3ab73d17d` |
| A2 | `phylogeo/ml/annotated_tree.nexus` | 1789 | `feb8c6f44dc86245` |
| A2 | `phylogeo/ml/confidence.csv` | 1797 | `8e05652db863e79a` |
| A2 | `phylogeo/ml/GTR.txt` | 409 | `a8c908dee5ea80dd` |
| A0 | `prep/dates.csv` | 668 | `22e1e19adb24b97d` |
| A0 | `prep/prep.fasta` | 18658 | `7518725a9cbfed3e` |
| A0 | `prep/prep.nwk` | 677 | `140139383068ae96` |
| A0 | `prep/states.csv` | 488 | `3ffe1f96ed9c3156` |

## 九、实际生效的参数

```json
{
 "trait": "region",
 "date_col": null,
 "date_trait": null,
 "methods": [
  "treetime"
 ],
 "relax": null,
 "reroot": "least-squares",
 "keep_root": false,
 "do_drt": true,
 "drt_perm": 200,
 "drt_seed": 42,
 "clock_rate": null,
 "prune": false,
 "prune_method": "fps",
 "prune_date_resolution": "month",
 "prune_max_reps": 3,
 "prune_min_snp_diff": 2,
 "prune_clade_cutoff": 5,
 "sampling_bias_correction": "auto",
 "pc": null,
 "cross_check": true,
 "rssp_bs": 0,
 "seed": 0,
 "n_skyline": 10,
 "gen_per_year": 50.0,
 "rng_seed": 0,
 "timeout": 3600,
 "divergence_tree": null,
 "marginal": false,
 "coord_table": null
}
```
