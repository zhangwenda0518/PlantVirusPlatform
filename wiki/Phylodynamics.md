# Phylodynamics · 进化动力学工具组

> 返回 [Home](Home) ｜ 上一页：[Report](Report) ｜ 下一页：[BEAST](BEAST)

进化动力学是平台在主管道（⑦ 建树）之上的扩展层：**分子钟定年、祖先状态与地理重构、种群动态、时空采样设计**。与 VirPhyKit（Yin et al., 2025, Ecol Evol）**方法学对齐、代码全部自研**（VirPhyKit 为 GPL-3.0，不复制其代码，只对齐输入/输出契约与算法口径）。

涉及引擎：`phylodyn_kit.py`（序列/元数据工具箱）、`phylodyn_trees.py`（树类引擎）、`phylodyn_local.py`（本地编排）、`phylodyn_alignqc.py`（比对 QC）、`phylodyn_govern.py`（元数据治理）、`treetime_ml.py`（TreeTime ML）、`treetime_rtt.py`、`rtt_drt.py`（RTT/DRT）、`bsp.py`（天际线）。

---

## 能力总览

### A. 序列/元数据工具（phylodyn_kit）

| 能力 | 说明 |
|---|---|
| SeqIDRenamer | 按 TSV 映射表批量重命名 FASTA 序列 ID |
| SeqGrouper | GenBank 记录表 + 分组映射表 → Group 列/分组计数 |
| **GeoSubsampler** | FASTA 时空降采样：均等采样 / 定区剔除 / 随机抽取（控制采样偏倚） |
| VirSpaceTime | Temporal.txt / Spatial.txt 解析 → 前端 plotly 绘图（时间轴/地点分布/地图点位） |
| 数据接入 | 三来源（explorer 导出包 / 手动输入 / 在线下载）+ 时间地点同检 + 标准化数据集 |

### B. 元数据治理与比对 QC（入口把关）

- **phylodyn_govern**：时间 + 地理元数据归一/标注/体检——**不猜、不静默、不改语义**（平台统一日期口径 YYYY / YYYY-MM / YYYY-MM-DD；小数年折算 ISO 入数据集）；
- **phylodyn_alignqc**：比对质量检查（长度 / gap% / 简并字符% / identity）→ 干净集 + 剔除集 + 报告（重组分析前预筛同源 `align_qc.py`）。

### C. 树类引擎（phylodyn_trees）

| 引擎 | 功能 |
|---|---|
| MCC 注释树解析 | BEAST MCC 树的节点状态/高度读取（`mot.py` 口径） |
| TempMig | 离散性状迁移时序 |
| **RRT**（root-to-tip 回归） | 分子钟检验：TreeTime ML 同口径重实现（`treetime_rtt.py` 走 `treetime_ml.run_clock_api`，速率/R² 直接读 TreeTime 收敛产物 `molecular_clock.txt`）。刻意**不**照抄上游的「正则抠日志 β/R² + scipy.linregress 自算 p」——实测那是数学恒等式（R²≡1、p≡0，见模块 docstring ①）；时间信号显著性改由 `rtt_drt.py` 的零分布检验给出，双口径修正 AUDIT 缺陷 |
| **DRT**（dating regression） | 定年回归诊断（与 RRT 配对使用） |
| RSPP | 序列抽样代表性 |
| **BSP**（Bayesian skyline plot） | 种群动态历史（`bsp.py`，含 ESS 有效样本量诊断） |
| LTT | lineage-through-time 曲线 |

### D. TreeTime ML 引擎（轨道 A，本地 ML）

`treetime_ml.py`（~1900 行）封装 **TreeTime 0.12.1**：

| 子能力 | 内容 |
|---|---|
| **clock 分子钟定年** | 时间树、进化速率、R²；可选**天际线**（coalescent 模型 + 置信带）；`audit_clock` 自动体检（叶数/速率/R² 阈值判定） |
| **mugration 地理重构** | 离散性状 ML 祖先态重构 + GTR 迁移模型 + 置信度；`audit_mugration` 体检；**采样偏倚校正系数**建议 `(1−Σpᵢ²)/(1−Σtᵢ²)`；**对拍交叉验证**（对拍树状态重算走廊，rel_diff > 0.2 判 divergent）；ML 与 Fitch 两条口径走廊并列成表 |
| **ancestral 祖先序列** | 祖先序列重构 + 逐枝替换 |
| **homoplasy 同塑性** | 替换饱和/同塑性检验 |
| 环境治理 | 逐级找 treetime 入口（env → platform.json → 3rd/tools → Scripts → PATH）、locale 剥离、不抛异常的超时封装、注释 Newick 解析 |

## 典型工作流

```
FASTA + 时间地点元数据
  → 元数据治理（govern：归一/体检）
  → 比对 QC（alignqc：剔除坏序列）
  → （可选）GeoSubsampler 时空降采样
  → 建树（复用平台 NJ/FastTree/RAxML-NG）
  → RRT/DRT 分子钟检验 → TreeTime clock 定年（+天际线 BSP）
  → TreeTime mugration 地理重构（迁移走廊 ML/Fitch 双口径 + 对拍验证）
  → 迁移弧线动画 GIF（见 Phylogeography）
```

## 界面与 CLI

- 专项分析页「进化动力学」工具组卡片（每能力一卡，任务中心统一调度）；
- 结果预览面板展示统计值 + 产物清单。

## 相关

- BEAST 轨道（贝叶斯 MCMC）→ [BEAST](BEAST)
- 地图动画 → [Phylogeography](Phylogeography)
- 建树基础 → [Phylo](Phylo)
