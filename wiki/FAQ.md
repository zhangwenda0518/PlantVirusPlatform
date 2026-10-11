# FAQ · 常见问题与已知事项

> 返回 [Home](Home) ｜ 上一页：[CLI-Reference](CLI-Reference) ｜ 下一页：[AI-Assistant](AI-Assistant)

---

## 常见问题

### 页面报 "浏览失败: TypeError: Failed to fetch" 或点击无反应？

浏览器连不上平台服务（服务端本身正常），即**服务进程已退出或无响应**：

1. 黑色控制台窗口被关闭了——关闭它就等于退出平台。重新启动平台（双击 `启动平台-桌面窗口.bat` 或 `python app.py`），浏览器会自动打开新页面；
2. 双击了多次平台，浏览器停在已关闭实例的旧页面——关掉旧标签页即可；
3. 大任务把内存/CPU 占满导致服务暂时无响应——等任务结束。

新版页面失联时自动显示红色诊断横幅，恢复后自动消失。

### 输入文件支持哪些格式？支持压缩吗？

- 测序数据：FASTQ / FASTQ.gz（`.fastq` `.fq` `.fastq.gz` `.fq.gz`），双端选 R1/R2，单端只选 R1；
- 基因组/参考：FASTA / FASTA.gz；GenBank `.gb/.gbk`（基因组图与参考集合）；
- `.zip` / `.rar` / `.tar` 请先解压出数据文件再选择；
- 文件浏览对话框只显示数据文件类型，浏览范围限平台目录内。

### 提示「宿主库不可用」？

到「数据库构建」页构建宿主库（需 TaxID），或取消勾选「① 宿主去除」阶段。

### SPAdes 报错 67 / non-ASCII？

旧版本问题，现已自动中转（`%TEMP%\vp_spades` 纯 ASCII）。若仍出现，检查 `%TEMP%` 路径是否含中文。

### 想用 ML 更严谨建树？

分析页把「建树工具」切到 **RAxML-NG**（随包 `3rd/tools/raxml-ng/`，GTR+G，ML 搜索 + 100 次 FBP 自举；固定种子可复现、线程自动封顶 8）。实测口径（209 条 × 966 nt）：RAxML-NG 含自举 730.6s，IQ-TREE 3 同数据 1049s——端到端快约 1.4 倍。只要拓扑不要支持值可直接 FastTree（秒级）；序列 <4 条或比对全 N 自动降级 FastTree 保树（写入任务日志）。

### 进化树想加入更多近缘参考？

工具箱「NCBI 参考序列下载」：Entrez 检索式（如 `Tobamovirus[ORGN] AND complete genome[TITL]`）→ 下载为参考集合 → 分析管道 ⑦ 的「NCBI 参考集合」填集合名。CLI：`python main.py ncbi-dl "<检索式>" -n <集合名>`，分析时 `--ncbi-refs <集合名>`。也可用 ICTV 级联下拉整科整属拿序列（见 [Toolbox](Toolbox)）。

### 想在网页上直接看比对差异（不用 SDT）？

结果中心 **MSA 查看（SNP-only 变异热图）**：只显示变异位点，ACGT 彩色热图 + 共识行 + 变异度柱，大比对自动分页；优先展示建树所用的清剪后比对。见 [Phylo](Phylo)。

### 引物对宿主特异性检查很慢？

首次需对 1.8GB 宿主基因组建 BLAST 库（10-30 分钟），此后复用。

### 想直接比较同属病毒的基因组结构？

比较基因组按四步组织：**参考序列获取（ICTV 级联/检索式/accession）→ 序列比对（MAFFT+trimAl，可编辑查看器）→ 进化树构建（科/属级，全基因组/CDS/PEP）→ SDT 同一性（属级）**。GenBank 集合一键提取 CDS/PEP 送基因建树。（同属共线性比较模块已移除；如需出版级共线性图建议用 clinker / LoVis4u 独立运行，GenBank 集合可直接作为其输入。）

## 已知事项与设计说明

### 中文路径兼容

SPAdes 与 BLAST(LMDB) 不支持含中文路径。平台自动处理：SPAdes 经 `%TEMP%\vp_spades`（纯 ASCII）中转拷回；BLAST 库建在 `%TEMP%\vp_blast`。其余工具（kunpeng/mafft/FastTree/RAxML-NG）原生支持。

### 宿主库大基因组内存问题（已修复）

kunpeng convert 阶段整批载入序列且无字节上限，多条大染色体同批会触发确定性 ~20.8GB 巨型分配（与机器内存无关）。平台在注入 taxid 时自动把每条序列切成 **≤1MB 片段**（`>{id}.p00001|kraken:taxid|N`），相邻片段重叠 34bp（k-1）——跨切口 k-mer 全部保留，全量实测 k-mer 数与整条建库**完全一致**（286,327,790）；2GB 内存即可完成 1.8GB 基因组建库（约 1.5 分钟，库 ~1.6GB）。建库为替换式：重试自动清理旧 library。

### 分类树来源

报告的科属种层级来自 ②b 鉴定表（逐参考 Uniq_Reads 按种/属/科累计）+ 参考库 ref_info 的 VMR_Family/VMR_Genus，**不加载全量 NCBI taxonomy**（避免内存压力；未收录 VMR 谱系的参考只画到种）。

### SDT 口径

平台用 MAFFT 比对 + 成对 gap 删除口径重算全长 identity 矩阵（与 SDT 算法一致），输出 CSV/热图；同时生成 `sdt_input.fas` 供 SDT v1.3 GUI 交互查看。

### 断点续跑

每阶段 `.done` 标记；SPAdes 已有组装结果时重跑直接复用；`--force` 强制重跑。

### 数据安全

GUI 仅监听 127.0.0.1；文件读写限制在平台目录内；下载仅限 NCBI/ENA/NGDC 官方域名白名单；subprocess 一律参数列表（shell=False）。

### SVG 的 MIME 类型（Windows 特有）

Windows 注册表把 `.svg` 关联成非标准的 `image/svg`，Flask 照搬后 Chromium 拒绝解码（`<img>` 裂图）。`app.py` 建应用前已显式注册 `image/svg+xml`——**该修复在 app.py 里，源码运行即时生效**。

### Markov jump 计数口径（BEAST）

跳转计数是**采样历史下的期望跳转次数**（不是树上边数、也不是 BSSVS 速率指示），与 Fitch 观测事件数是两条口径——面板上分开摆。详见 [BEAST](BEAST)。

### 参考库不做在线自动更新

需要更新时把新版 4 件套（两个 FASTA + 两个 Info.tsv + DATA_VERSION）整体替换进 `databases/virus_ref/`，平台按 DATA_VERSION 自动重建元数据缓存与 BLAST 库。ICTV VMR 例外：`ictv-update` 在线/本地 xlsx 均可。

## 相关

- 启动与安装 → [Quick-Start](Quick-Start)
- 平台安全边界 → [Architecture](Architecture)
