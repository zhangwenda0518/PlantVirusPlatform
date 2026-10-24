# MMPV-RNA 结果解析与绘图体系学习笔记

> 源仓库：`D:\桌面\延伸基因组\MMPV-RNA`（大规模植物病毒宏基因组挖掘平台 v3.1，七大管线）
> 整理日期：2026-10-12。四路并行精读：virome_analysis_pipeline（已知病毒深度分析）、virome_phylo_pipeline（系统发育/群体遗传）、public_metadata + virome_discovery 管线、论文图工程 + 结果分析 skills。
> 用途：为植物病毒分析平台的"结果解析 + 绘图 + 报告"模块提供借鉴蓝本。

---

## 0. 全局架构：结果解析怎么组织

MMPV 的核心理念是**三层分离**：

```
管线代码（跑工具、产 TSV/图） → 汇总表/结果文件体系（TSV/CSV 矩阵） → 视图层（静态 SCI 图 + 交互 HTML 报告 + 论文图脚本）
```

- **目录键名注册表**：所有管线不写死输出目录名，统一由 `mmpv_common/io_layout.py` 的 `layout_dirs()` 解析（standard/legacy 双布局，环境变量 `MMPV_IO_LAYOUT` 切换，键名如 `a_detect/a_variants/a_report`）。解析脚本开头 `sys.path.insert(0, _REPO_ROOT)` 后 import。
- **收集层与渲染层解耦**：报告只是汇总 TSV 的视图；每张图、每个数字都能指回源文件。
- **结果分析"技能化"**：`virome_analysis_pipeline/skills/` 下 8 个纯 Markdown SKILL.md，把"如何解读结果"写成可执行知识库（与流水线代码完全解耦）。

---

## 1. 绘图统一规范（最值得照抄的部分）

### 1.1 配色：单一真源 `sci_palette.py`

位置：`virome_discovery_pipeline/utils/sci_palette.py` 与 `virome_phylo_pipeline/utils/sci_palette.py`（同源，两侧改色必须同步）。

```python
OKABE_ITO = ["#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","#000000","#999999"]  # 色盲安全，顺序固定
OKABE_BRANCH = [BLUE,VERM,GREEN,ORANGE,PURPLE]  # 去掉最浅的 yellow/sky，保证 0.5pt 细线可辨
SEQ_CMAP = "viridis"    # 连续数据
DIV_CMAP = "RdBu_r"     # 发散数据
SEMANTIC = {...}        # 业务语义色表：verdict(KEEP绿/REVIEW橙/DROP朱红)、genome(DNA蓝/RNA朱红/Viroid绿)、novelty、quality(CheckV五档)、stage
```

铁律：**所有业务色一律从 OKABE_ITO 取，禁止新增游离 hex**；语义色表保证"同一含义跨图同色"（静态 matplotlib、plotly 交互、HTML chart.js 三个渲染载体共用）。

### 1.2 样式：`sci_rcparams()` 返回整份 rcParams 字典（强制共享）

- 字体 `Arial/Helvetica/DejaVu Sans`；`pdf.fonttype=42`（TrueType 嵌入，期刊要求可编辑）；`svg.fonttype='none'`；`axes.unicode_minus=False`；白底；`spines.top/right=False`；`legend.frameon=False`。
- **必须接管 `axes.prop_cycle=cycler(OKABE_ITO)`**（否则混入 matplotlib 默认 tab10）；seaborn violin/box 要 `saturation=1`。
- 手抄 rcParams 必漏键——实测漏 `spines.top|right` 会导致"同一份代码单独跑 150dpi/Type3、跟在别的模块后跑 300dpi/Type42"的诡异不一致。所以规定：任何模块只许调 `apply_sci_style()`，不许自己 update。

### 1.3 尺寸：画布宽 = 最终排版栏宽

```python
COL_W_2 = 4.58  # 单栏（英寸）
COL_W_1 = 9.32  # 通栏
figsize_for(ncol, ratio=0.62)
```

核心逻辑：**画布宽=目标印刷栏宽 → 印刷缩放≈1.000**。旧版 13/16.5 in 画布排进通栏缩放 0.717/0.565，8.6pt 标签印出来只剩 4.9pt——已全部返工。印刷字号 ≈ 源字号 ×(面板宽/源画布宽)，所以 ≥8.8 in 的源图一律占通栏。

验收手段：`pypdf` 读 PDF MediaBox 算 `9.32/width` 必须落在 0.9–1.1；检查字体嵌入（FontFile2>0、Type3=0）；`bbox_inches="tight"` 兼作"元素越界探测器"（画布异常变大=有东西画到框外）。

### 1.4 dpi 与格式

- 常规 `dpi=300`，PDF+PNG 双写；SDT 热图 600 dpi；Circos `png(3000,3000,res=300)`。
- 发表图 PDF 为主 + 补充 SVG（可编辑文字）；探索图 PNG；交互图单文件 HTML（plotly `write_html(strict=True)` 做零外链自检）。

### 1.5 画布封顶 + 标签分级（图随数据放大必须配刹车）

`batch_plot_virus_depth.py` 的教训：262 样本时产出 23550×1768 px 横条、6809×47127 px 竖条，无法阅读也无法嵌报告。解决：

```python
SAMPLE_BAR_MAX_W_IN / ABUND_HEAT_MAX_H_IN / STACKED_MAX_W_IN  # 三个封顶常量
DENSE_ELEM_THRESHOLD = 60   # 元素>60 自动关闭逐元素文字标注、缩小行标签字号
```

### 1.6 版式自动体检（把"图好不好"变成可回归测试）

`scripts/barbarum/_check_discovery_overlap.py`：劫持 `plt.Figure.savefig` → `fig.canvas.draw()` → 收集所有 Text 图元的显示 bbox 两两求交，判据 **0 overlapping pairs（TOL 0.5px）**，退出码即结论。26 figures、0 重叠通过。

---

## 2. 结果解析模块逐个拆解

### 2.1 virome_analysis_pipeline（已知病毒深度分析）

| 模块 | 输入 | 输出 | 要点 |
|---|---|---|---|
| `filter_summary.py` | `all_viruses.best.summary.tsv`（样本×病毒最佳命中） | high_conf/discarded/per_sample/per_virus 表 + 双面板图 | `KNOWN_COLS` 把 CLI 阈值映射为列名；统一 `thresholds` 字典累加 mask；列缺失只 WARN 跳过不崩。HQ 口径：coverage≥50%、depth≥5×、reads≥100、TPM≥1、Poisson≥0.3、ANI≥95% |
| `stage_funnel_summary.py` | 五阶段产物 | `stage_funnel_summary.tsv` + 漏斗图 + 气泡留存图 | 五阶段 raw→best→highconf→assembly→extracted；`build()` 返回 dict 可被报告直接 import；漏斗图手绘梯形（`Blues(np.linspace(0.88,0.42,n))` 五层色阶+Polygon 连接带）；气泡图有检出=圆泡（大小/颜色随值）、淘汰=灰底空格 |
| `virus_detection_summary.py` | 三处产物 | 一张 per-virus 综合表 | 以 `Rep_Accession` 为键拼 best/high_conf 计数、consensus 数（含 N 含量分箱）、组装成功数、提取数+失败原因回填 note——一张表讲清「检出→过滤→组装→提取」的衰减 |
| `virus_variants_analyzer.py` | iVar VCF | 9+1 图引擎 | **iVar 输出 8 列 sites-only 无 GT** → 必须用 allele_frequencies.tsv 重建 0/1 矩阵（修复前 Hamming 全零矩阵） |
| `snpeff_analysis.py` | `*.ann.vcf` | 六矩阵 + OncoPrint | `determine_oncoprint_state(af, impact)`：连续 AF 离散成 9 态（Minor<0.5 / Major 0.5-0.95 / Fixed≥0.95 × Truncating/Missense/Synonymous/Regulatory）；Quasispecies=Consensus(>80%)/iSNV(5-80%)；`min_occ = 1 if n_samples<2 else max(2, 3%×n)` 修"单样本矩阵必空→R oncoPrint 静默失败" |
| `panvirome_landscape.py` | 基因坐标+π+dN/dS | 三轨进化全景图 | 上轨基因箭头、中轨 SNP 密度+π 双轴、下轨 dN/dS 柱（>1 红 #E31A1C，<1 绿 #74C476） |
| `panvirome_novel_known.py` | blastn/mmseqs/CDD | 三层证据链判定 | nt（qcovs≥50% 才采信）→蛋白→CDD 结构域；nt≥85%=已知、50-85%=新种候选；判定阈值原文写进 HTML 报告的 Verdict Rules 卡片 |
| `virema_summary_report.R` | ViReMa BED | Circos 4-track + 断点偏好 | A/C/G/T 固定色 #E41A1C/#377EB8/#4DAF4A/#984EA3；`circos.genomicLink(arr.type="triangle", directional=1)` |
| `batch_draw_pymol.py` | AlphaFold CIF + DRHIP csv | 正选择位点 3D 映射 | 顶刊风格白底 orthoscopic + 四象限 offset 矩阵散开重叠标签；落盘 `*_ManualReplay_Settings.pml` 供 GUI 微调 |
| `generate_pipeline_report.py`（2580 行） | 全管线产物 | 自包含单文件 HTML | 不用 Jinja2，纯 f-string 拼 HTML+内联 CSS/JS；`img_to_base64`（>2000-3500KB 不内嵌）；侧栏导航+指标卡+漏斗章+逐 stage+每病毒小节+自动生成 Methods；**AI 按钮**（前端直连 OpenAI/DeepSeek/Kimi/Ollama，key 存 sessionStorage）+ 防幻觉系统提示词（只用所给数据、标明 Figure X、审慎措辞、数据不支持明说"未检出"） |

### 2.2 virome_phylo_pipeline（系统发育/群体遗传）

三大绘图模块：

- **`publication_figures.py`**：Fig1 Root-to-Tip 回归（`Bio.Phylo` 累加枝长——注意 `unit_branch_lengths=False`，旧参数 True 是"边数"导致 R² 从 0.966 错成 0.028；R²/slope/TMRCA 写进左上 monospace 统计框）、Fig2 MCC 时间树按地点着色、Table S1（ESS 用 Geyer initial positive sequence 与 Tracer 逐参数对账 + LaTeX 版输出）。
- **`popgen_figures.py`**：逐比对列 Shannon 熵面积图 + 基因区 axvspan。**核心正确性机制**：熵谱横轴是比对列、基因坐标是参考序列坐标，MAFFT 含 gap 时两者不等——必须重建"参考位置→比对列"映射，否则基因标注整体错位；命中 0 或 ≥2 条参考序列就**不画基因区**（宁可不画）。
- **`phylogeo_figures.py`**：五面板 GridSpec(3,2)（MCC 树/圆周布局迁移网络/skyline/时钟速率直方/迁移 BF 热图）。skyline 时间轴由 `groupSize1..k` 累加还原（旧版 `np.linspace(0,1,n)` 是编造刻度），并用 rootHeight 做量纲自检（差>5% 降级为等距）。**有效性闸门**：<2 面板有数据**拒不出图**（SystemExit(2)）、2-4 面板打红色 INCOMPLETE 横幅——防止五面板全灰占位被当成"内容完整"收录。

解析模块：

- `popgen_analysis.py`：不依赖 R，自实现 π（pairwise deletion 只算 ACGT，防 pypopart 只删 `-` 不删 N 高估）、Tajima's D（Nei 1987 完整方差项）、Hd、Fst 双法（单倍型频率法 + 999 次置换 seed=42）、AMOVA。
- `selection_analysis.py`：解析 HyPhy BUSTED/FEL/MEME/FUBAR/RELAX JSON + PAML codeml；FEL/MEME 过 BH-FDR；FUBAR 取 `Prob[alpha<beta]≥0.9`（旧版取反列结论全反了）。
- `recombination_analysis.py`：RDP5 csv（≥3 方法支持才 high-confidence）+ 自研 DnaSP 复刻算 Rm → 100nt 窗口分箱取 hotspot → 与 RDP5 断点取交集 → 注释到基因名。
- `geo_analysis.py`：RRT 距离置换（seed 固定保可复现）、父子节点地点变化=迁移事件、RRT-tree 用真实根状态 vs 随机树同名状态（空随机集→防假性 PASS 守卫）。

### 2.3 public_metadata_pipeline（公共数据检索）

- `gsa_sra.plot.py`：SCI 级 6 面板（年度折线/环形/机构/组织/地区/发育期 Top10 横条），统一 run 计数口径；缺列只跳该面板不整图崩；机构名 title-case+错拼修复；PDF+PNG(600dpi)。
- `stat_prescreen.py`（1496 行，刻意无图）：NCBI/GSA 预筛 k-mer 分类判定链；**关键口径补丁**——NCBI 分母是已鉴定 reads、GSA 是总 reads，不可比 → 额外输出归一化 `virus_pct_ident` 再出 decision（DOWNLOAD/REVIEW/SKIP）。
- `stat_virus_profile.py`：**`mark_leaf()` 父子节点去重**（同分支科/属/种同现会让父节点检出率虚高）→ 同时给两套口径 `n_samples/prevalence_pct` 与 leaf 版；4 图+1 组合面板。

### 2.4 virome_discovery_pipeline（病毒发现）

- `utils/taxonomic_sankey.py`（1203 行）：plotly `go.Sankey`，7 级层级（Realm→Genus），缺失填 `Unknown_Genus`，只有末级应用 min-flow 过滤；40+ 内嵌色板；交互 HTML 或 kaleido 静态导出。
- `utils/acvirus_tree_pro.py`：MAFFT→trimAl→FastTree/RAxML-NG/IQ-TREE；`root_tree` 三级兜底（taxonomy 引导→MAD→midpoint）；`normalize_support` 统一 0-100；输出 annotated.nwk + phyloxml + **auspice.json（Nextstrain v2，可拖入 auspice.us）**；绘图 = 树+共线性双面板 GridSpec(1,3)。
- `report_pipeline.py`（3265 行）：16 个 `_collect_*` 收集函数各自"找目录→解析产物→写 TSV→登记 KPI"；Chart.js/UpSet.js 全文 inline 进 HTML（离线可用）；Sankey/旭日图 base64→Blob URL→iframe（绕 data URI 大小限制）；表内 filter/排序/导出 CSV/分页/列显隐全客户端；漏斗 KPI 从 Key_Metric 字符串正则回捞，`Raw Reads → Clean → Host-free → Contigs → Viral Seqs → vOTU → HQ vOTU` 横条+log 轴。

---

## 3. "结果→结论"解读规则（skills 里的数字化阈值）

`virome_analysis_pipeline/skills/` 下 8 个 SKILL.md，分两类：

**结果分析型**（virome-results-analysis 等 3 个）骨架：When to use → Inputs you must locate（缺失记 `N/A (stage not run)`，绝不推断）→ Core statistics to compute → Standard audits（PASS/WARN/FAIL 判据）→ 双交付物（`RESULTS_*.md` 叙事 + `AUDIT_*.md` 质检表）→ Manuscript alignment。

**论文写作型**（mmpv-paper-* 5 个）：权威来源表 + **Lycium result anchors（把已定稿数字写死成锚点）** + 反造假规则（不得发明的量、`[N — TO CONFIRM]` 占位、冲突标 `[WARN] manuscript-anchor conflict` 而非改叙事）。

数字化阈值举例：

- 新颖度三级：★已知 = nt 一致性≥90% 且覆盖≥90%；★★潜在新种 = aa≥70% 但 nt<90%；★★★全新 = aa<70%。
- 证据链评分 = 0.30·BLASTN + 0.30·BLASTX + 0.40·domain → KEEP/REVIEW/DROP。
- RDP5 重组：7 方法联合，≥3 方法支持才算 high-confidence。
- dN/dS：dS≈0 会让比值爆炸（实测 49 个 dS=0 对把 y 轴拉到 29）→ 画 dN-vs-dS 双对数散点+中性线，看**落在线上方的比例**；比值只对 dS>0 子集算并记录被剔对数。
- 时间信号硬门：TreeTime RTT R²>0.3 + DRT（20 次重排，真实 R² 需>95 百分位）+ BETS 三道门全过才允许 BEAST 定年；未过就不报 TMRCA。
- Tajima's D 负值须同时给"扩张"与"纯化选择"两种解释，不单选。
- 审计判据：漏斗单调性（每级下游≤上游，违反即 FAIL）；CDD 通过率 sanity（很高→污染、近零→过滤没咬）；SDT 口径 >95% 同株 / 85-95% 种内 / <85% 异种。

---

## 4. 论文图工程化

- **架构示意类图**：`doc/figures/make_figure1_full.py` 用 `fig.add_axes([0,0,1,1])` 单轴全幅 + 毫米数据坐标手工摆放（190×228mm 对应期刊尺寸），自建 `rbox()/arrow()/chip()/panel()` 图元函数；面板字母 `ax.text` 手工放（与色带联动精确定位）。
- **数据多面板图**：`scripts/barbarum/make_discovery_sci_figs.py` 用 **GridSpec(5,3) + 声明式 LAYOUT 表** `[("A",0,0,1,2), ...]`；每个面板函数只画一张轴、异常被 try/except 捕获画占位不拖垮整图；面板字母挂 `ax.transAxes` 左上外。`make_sci_composite_figs.py` 多个手工 add_gridspec + 字母绝对锚定。
- **可重跑重画纪律**：全 argparse 参数化；硬断言守卫（文字高度换算到数据坐标，溢出 assert 报错不产图）；**symlink 防写穿**（savefig 前 unlink，避免写回源数据）；分母 FATAL 校验（分母<检出样本数直接退出）；每轮修订 `rev_X` 后缀保留历史；`reports/*.md` 记录 md5 清单+复现命令。
- **图注三层**：图内面板标题 / `fig_manifest.json`（机读溯源）/ 稿件 Figure legends 集中管理（每条附 `[Source: <脚本/产物目录>]`）。

---

## 5. 值得借鉴的设计模式（汇总清单）

1. **样式单一真源 + 顺序无关性**：`sci_rcparams()` 返回整份字典强制共享；每个 savefig 都显式传 dpi。
2. **退化守卫 / "声明即存在"**：文件非空才登记交付物；多面板<2 拒出；空图例不画；失败删陈旧图；空随机集防假性 PASS。核心信条："算不出来就不产出对应图元，绝不静默降级"。
3. **口径同源 + 历史坑注释体系**：小数年、序列名、ESS/HPD 统一走同一 utils；源码注释写"症状+根因+实测数字+修复日期"。
4. **两套口径并列输出**（父/叶节点、NCBI/GSA 分母）——不替用户选，写明口径。
5. **漏斗三件套**（表+双面板漏斗图+留存气泡图）+ 可 import 的 `build()` 接口，跨项目通用质控叙事。
6. **容错三件套** `_fnum/safe_read_csv/safe_plot`：`-`/NA/空文件到处容错；R 侧单点输入 `density()` 抛错降级跳过仍出主图。
7. **批次绘图封装**：`ProcessPoolExecutor` + 样本级任务元组 + 模式开关 + 失败清单汇总。
8. **统计口径写进产物**：队列分母+Wilson 95% CI+Fisher+BH-FDR+共检零模型（log2 obs/exp），把"原始计数误导"从代码层面消掉。
9. **离线自包含 HTML**：Chart.js inline、图 base64、TSV data URI、大图 Blob iframe——单文件交付/存档。
10. **AI 按钮模式**：前端直连多家 LLM（key 存 sessionStorage）+ 每 stage 组装"数据+图清单+写作纪律"提示词 + 防幻觉系统提示词。
11. **结果分析技能化**：输入文件表→核心统计→固定审计→双交付物（RESULTS/AUDIT）骨架，阈值全部数字化。
12. **配置驱动**：`datasets.yaml` 管数据集注册（坐标/配色/基因区），绘图 CLI 由它遍历；config 键必须同步登记映射否则静默失效（踩过坑）。
13. **交互与静态分层**：plotly/Leaflet 用于探索，PDF 用于发表，同一数据双路并存。
14. **`sys.executable` 替代裸 `python`** 调子进程，避免落到无 matplotlib 的解释器静默不出图。
15. **精确 Stage→章节映射表**替代目录名前缀猜测（旧实现串章 bug）。

---

## 6. R 与 Python 分工

- **Python/matplotlib 绝对主力**（30+ 模块）；R 只在四处：ComplexHeatmap oncoPrint、circlize Circos、ggplot2 棒棒糖/ VirSpaceTime、RDP5 Rgrapher。
- R 全部走 `subprocess` + `Rscript`（无 reticulate）；第三方 R 源码**一字节不改**，运行时字符串替换路径/配色/尺寸（`utils/rdp5_bdp_plot.py`：临时目录复制输入——R 不认 Unicode 路径，`encoding='latin-1'`）。
- R 侧脚本用 JSON 哨兵回传结构，依赖缺失只告警给安装清单。

---

## 7. 技术债与反面教材（同样有学习价值）

- `publication_figures.py`（NATURE_COLORS）与 `phylogeo_figures.py`（NATURE_PALETTE）尚未迁到 Okabe-Ito 真源——同库内可见的配色不一致。
- 旧版 skyline 用 `np.linspace(0,1,n)` 编造时间刻度（已废弃）；旧版 Fig1 root-to-tip 用错 `unit_branch_lengths=True`（R² 0.966→0.028）。
- 旧版 FUBAR 取错列（正选择↔纯化选择结论全反）；旧版 phylogeo significant 只认 'Yes' 丢掉全部箭头；旧版画布 13-30 in 未按栏宽返工。
- 未接线的占位脚本 `popgen_r_backend.R`（文件头明确写是假值）——占位要诚实标注。
