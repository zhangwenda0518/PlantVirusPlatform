# ORF · ⑥ ORF 预测 与 ⑥b 功能注释

> 返回 [Home](Home) ｜ 上一页：[Host-Prediction](Host-Prediction) ｜ 下一页：[Phylo](Phylo)

对应引擎：`orf.py`（⑥）、`orf_annot.py` + `hmm_annot.py` + `bfvd_search.py` + `genome_diag.py`（⑥b）。

---

## ⑥ ORF 预测

- **pyrodigal**（meta 模式）：基因预测主引擎，输出 faa/ffn/gff；
- **pyrodigal-rv**：RNA 病毒优化模型（核糖体滑移、重叠 ORF 等），产物 `pyrodigal_rv.*`；
- **orfipy**：仅在显式指定时运行（pep/nt/bed）；
- 参数：最小 ORF 长度 `--min-orf-aa`（默认 100 aa）；
- 产物（`04_orf/`）。

## ⑥b ORF 功能注释（三层）

回答「每个 ORF 是什么蛋白、属于哪个科属、基因组什么类型」。

### 层 1 · 序列同源（主策略）

- **DIAMOND（首选）→ MMseqs2 → BLAST+ blastp** 自动降级，对照 **NCBI RefSeq 病毒蛋白库**（release viral 全量，首次运行自动从 NCBI 官方域名下载 ~107MB 并建库，之后全缓存）；
- **E-value 加权类别投票**（依据 LanderDC/annotation_benchmark 基准：序列同源是病毒蛋白注释最准确的主策略——BLASTp 78.4% > 结构折叠 72.2% > 蛋白语言模型 63.6%）；有效类别权重须 ≥ 最佳「假想蛋白」类别的 50%；
- 命中物种经 NCBI taxonomy 谱系解析到属/科，再查 **ICTV 科级对照表**（265 科，源自 LazypipeX）得**宿主类别**与**基因组类型**（ssRNA(+)/dsDNA/…），与 ④ 宿主预测交叉印证。

### 层 2 · HMM / 结构域（远缘兜底）

- **pyhmmer 流式扫描 Pfam 病毒库**（`Pfam-A-Viruses.hmm`，virsorter2 病毒子集 1,074 profiles；内存恒定 ~50MB）；
- 双门槛过滤：域级 i-Evalue ≤ 1e-3 **且** 模型覆盖率 ≥ 0.5；HMM 内嵌 NAME/ACC/DESC 直接得可读产物名；
- 序列层已注释 ORF 保留层 1 结果、HMM 命中记入 `hmm_hits` 列；未命中 ORF 由 HMM 兜底（evidence=hmm）；
- **BFVD 病毒域层**：mmseqs2 搜索 BFVD（`bfvd_search.py`，-s 7.5，evalue<1e-3）；BFVD 为病毒专属域库、命中即病毒证据（免白名单），product/Organism 取库内 `bfvd_uniprot_names.tsv`——与 SUVTK 提交链 featuretable 的 product **同库同命名**；（原 CDD 结构域层因 ~5GB 索引反复 OOM，2026-09-29 换为 BFVD）；
- 库缺失时对应层整体跳过，不报错；可用 `--libs prot,pfam,bfvd` 指定启用层。

## 产物（04b_orf_annot/）

| 文件 | 内容 |
|---|---|
| `orf_annotation.tsv` | 逐 ORF：产物/物种/identity/覆盖度/E-value/科属/宿主/基因组类型/功能类别/证据/hmm_hits/bfvd_hits |
| `orf_function_summary.tsv` / `orf_family_summary.tsv` | 功能类别与科分布汇总 |
| `contig_function_profile.tsv` | 逐 contig 功能画像 |
| `orf_annotation.gff3` | pyrodigal GFF 追加 product/organism/category/evidence 属性 |
| `genome_diagrams/<contig>.svg` | littlegenomes 风格线性示意图：ORF 按功能类别着色方向箭头 + HMM/BFVD 结构域窄条 + 比例尺（零依赖 SVG） |

类别/科分布图与 Top 注释表自动嵌入 ⑩ 报告。

## 界面与 CLI

- 分析管道 ⑥/⑥b 卡片：`📊 查看`（注释表就地浏览）；结果中心可看功能示意图。

```bat
python main.py orfa --sample S1 [--threads 8] [--libs prot,pfam,bfvd] [--force]
python main.py analyze --r1 ... --sample S1 --stages orf,orfa --min-orf-aa 100
```

## 参考库维护

替换 `databases/annot/prot/viral_prot.faa` 后删除 `organism_tax.tsv`、`db_info.json` 重跑即自动重建索引与搜索库。（BFVD 库随包在 `databases/misc_db/suvtk/`，与 SUVTK 提交链共用；可用 platform.json `databases.bfvd` 指定其他前缀。原 CDD 库位说明已随 CDD 下线移除。）

## 相关

- 基因组类型/宿主类别 ↔ ④ → [Host-Prediction](Host-Prediction)
- ORF 注释如何被用 → [Phylo](Phylo)（PEP/CDS 建树）、[Genome-Plots](Genome-Plots)（着色注释）
