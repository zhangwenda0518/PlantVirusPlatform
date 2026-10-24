# 变异分析工具 Windows 可行性调研（实测版）

> 调研时间：2026-09-08
> 调研方式：全部本机实测，不接受文档推断

## 一、结论速览

| 工具 | Windows 可行性 | 本机状态 | 结论 |
|---|---|---|---|
| **SnpEff 5.4c** | ✅ 已跑通 | 已部署，端到端验证通过 | **采用** |
| iVar 1.4.4 | ❌ 不可得 | 无包、无源码构建链 | 放弃 |
| SNPGenie | ⚠️ 需 Perl | 本机无 perl | 暂缓 |

---

## 二、SnpEff（已落地）

### 2.1 部署事实

| 项 | 值 |
|---|---|
| 版本 | SnpEff 5.4c (build 2026-02-23) |
| 位置 | `tools/snpeff/snpEff/snpEff.jar` (27.3 MB) |
| Java | 平台自带 Temurin **JDK 21.0.12.1 LTS** (`tools/jdk21/jdk-21.0.12.1+1/`) |
| 下载源 | `https://snpeff-public.s3.amazonaws.com/versions/snpEff_latest_core.zip` (63.6 MB) |
| JDK 源 | `https://api.adoptium.net/v3/binary/latest/21/ga/windows/x64/jdk/hotspot/normal/eclipse` (195.6 MB) |
| 磁盘占用 | jdk21 523.5 MB + snpeff 229.4 MB = 753 MB |

### 2.2 为什么必须自带 JDK 21

本机系统 `java` 是 **1.8.0_341（class 52）**，SnpEff 5.4c 编译为 **class 65（Java 21）**。
实测报错：

```
java.lang.UnsupportedClassVersionError:
org/snpeff/SnpEff has been compiled by a more recent version of the Java Runtime
(class file version 65.0), this version of the Java Runtime only recognizes
class file versions up to 52.0
```

SnpEff 官方文档也已明确「requires Java 21 or later」。旧版 4.3t 只存在于 SourceForge，
实测直链返回 HTML 重定向页（0.1 MB）、`downloads.sourceforge.net` 返回 403，
GitHub releases 资产为空数组，取用不可靠。**自带 JDK 21 是唯一稳妥路线。**

### 2.3 端到端实测结果

测试对象：`OR489165.1`（Cytorhabdovirus sp. 'lycii'，枸杞 crinkle 病毒，14812 bp，6 个 CDS：N/P/P4/M/G/L）

```
GenBank (Bio.Entrez.efetch)
   ↓ SnpEff build -genbank
snpEffectPredictor.bin + sequence.bin
   ↓ poscounts.tsv → VCF
SnpEff ann → ANN 字段 → 结构化 TSV
```

建库输出：
```
Protein check: OR489165  OK: 6  Not found: 0  Errors: 0  Error percentage: 0.0%
```

注释结果（4 条编码区变异）：

| POS | 变异 | 基因 | 效果 | 影响 | HGVS.c | HGVS.p |
|---|---|---|---|---|---|---|
| 438 | A>G | N | start_retained_variant | LOW | c.3A>G | p.Met1? |
| 8000 | T>C | L | missense_variant | MODERATE | c.439T>C | p.Ser147Pro |
| 8002 | T>A | L | synonymous_variant | LOW | c.441T>A | p.Ser147Ser |
| 8005 | T>G | L | synonymous_variant | LOW | c.444T>G | p.Leu148Leu |

单库建库耗时 **0.27 秒**。

### 2.4 三个必须记住的坑

1. **中文路径 + Java 编码**
   不加 `-Dfile.encoding=UTF-8` 时路径显示为 `D:\����\ֲ�ﲡ������ƽ̨`，
   文件读取随即失败。所有调用必须带该参数。

2. **config 里 `data.dir` 必须用相对路径**
   `data.dir = D:\...\data` 会被 SnpEff 当成转义序列吃掉，报：
   ```
   Cannot read file '...\snpeff_work/D:����ֲ�ﲡ...data/OR489165/genes.gbk'
   ```
   正确写法 `data.dir = ./data`，且**工作目录必须切到 config 所在目录**执行。

3. **染色体名必须与库内一致（保留版本号）**
   SnpEff 库内染色体名是 `OR489165.1`，VCF 写 `OR489165` 会报
   `ERROR_CHROMOSOME_NOT_FOUND`，ANN 字段降级成空白注释。
   平台端必须做自动对齐（从 GenBank `VERSION` 行取 ID）。

### 2.5 附带收获：SnpEff 会校验参考碱基

人工构造的 VCF 中 REF 与库内序列不符时，输出
`WARNING_REF_DOES_NOT_MATCH_GENOME`。这说明它不是单纯按坐标套注释，
而是真的做了参考序列一致性检查，对结果可信度是正面信号。

---

## 三、iVar（不可行，已定性）

| 项 | 实测结果 |
|---|---|
| bioconda 平台 | 仅 `linux-64` / `linux-aarch64` / `osx-64` / `osx-arm64`（61 个文件），**无 win-64** |
| 官方二进制 | 不提供 Windows 版 |
| 源码构建 | 需 autotools（autoconf/automake/libtool）+ HTSlib + GCC≥5.0 |

本机实测：gcc 16.2.0 ✅、samtools 1.24 + htslib 1.24 ✅、make ✅，
但 **autoconf / automake / libtool / perl / cpan 全部缺失**，且无 bash / Git Bash / WSL。
按官方构建路径堵死。

**功能重叠判断**：平台现有 `viral_consensus/variants.tsv` 已产出
`contig/pos/ref_base/cons_base/depth/A/C/G/T/N/gap/major_base/major_freq/minor_base/minor_freq/type`，
iVar 的主要能力（iSNV 检出）已覆盖，边际价值低。

---

## 四、SNPGenie（需 Perl）

官方描述：「command-line application written in Perl, with no additional dependencies
beyond the standard Perl package」。

- **输入**：单序列 FASTA + GTF 注释 + VCF/CLC/Geneious SNP 报告
- **输出**：πN/πS、dN/dS、gene diversity
- **三种分析**：within-pool (`snpgenie.pl`)、within-group、between-group
- **VCF 必须指定 `--vcfformat=4`**

本机 perl 缺失。若要用，需额外部署 Perl（Strawberry Perl / MSYS2）。
与 SnpEff 的分工很清楚：SnpEff 做功能注释，SNPGenie 做群体遗传学指标（πN/πS）。

---

## 五、平台集成建议

### 5.1 数据链路

```
样本 poscounts.tsv (viral_consensus 产物，已有)
   ↓ poscounts_to_vcf()          ← kv_variant.py 已实现并验证
VCF (CHROM 自动对齐库内命名)
   ↓ SnpEff build -genbank       ← GenBank 复用绘图段缓存
SnpEff 库
   ↓ SnpEff ann
ANN 注释 VCF
   ↓ parse_ann()
结构化变异注释 TSV (基因/效果/影响/HGVS.c/HGVS.p)
```

### 5.2 已验证的数字

- 平台自产 poscounts → VCF：NODE_1 得 **43 条变异**，
  与 `summary.json` 的 `n_isnv: 43` **完全一致**
- 单库建库 0.27 秒，8464 条参考库理论约 38 分钟
- 建库必须有 GenBank 注释，`databases/virusref_db/final.cluster.ref.fasta` 只有序列，
  注释需走 NCBI（复用 `kv_plot.fetch_gb`）

### 5.3 待定问题

- 批量建库策略：全量 8464 条 vs 按样本检出病毒按需建库
  （倾向按需，避免 8464 次 NCBI 请求）

---

## 六、snipit 1.7（跨样本 SNP 分型图——已落地）

> 补充调研：2026-10-13。上文 SnpEff 链路回答「变异的功能效应」，
> snipit 回答「同一种病毒的 N 个样本谁携带哪些共享/私有 SNP」——
> 一行一样本的分型视图，是既有三张图（位点分布/效应谱/AF-DP）不覆盖的图型。

### 6.1 部署事实（实测）

| 项 | 值 |
|---|---|
| 版本 | snipit 1.7（PyPI 最新；Bioinformatics 2024, btae510 配套工具） |
| 许可证 | **GPL-3.0**。本仓库源码分发 + `3rd/` 本就随仓库带 bcftools 等 GPL 工具，作为独立进程调用无传染问题；README 致谢已登记 |
| 依赖 | 仅 `biopython>=1.78` + `matplotlib>=3.5`（纯 Python wheel，29KB）——平台既有版本即满足，**零新增传递依赖**，无需像 SnpEff 自带 JDK |
| 安装 | `pip install snipit`（console script 在用户 Scripts 目录，不一定上 PATH；平台用 `python -m snipit.command` 调用） |
| 中文路径 | ✅ 无坑（纯 Python 不走 JVM，输出到 `…\中文目录测试\` 实测成功） |
| 输出 | png/jpg/pdf/svg/tiff + `snps.csv` |

### 6.2 实测坑

1. **序列 ID 全串匹配**：平台 FASTA 头是 `|` 分隔长 ID，`-r` 必须传
   整串（`NC_001440.1|Cucumber_mosaic_virus_RNA3_RefSeq`），只传
   accession 报 `input reference not found in alignment`。
2. **snps.csv 是固定文件名**（`write_out_snps` 硬编码 `output_dir/snps.csv`），
   多条参考连续出图会互相覆盖——实现里先落暂存目录再改名归位
   （`kv_variant_plot._snipit_one`）。
3. **vcf/ 目录混着 caller 中间产物**（`all.variants.raw/soft/filtered.vcf`、
   `demo.variants.vcf`），glob `*.vcf` 会捞进来；靠 `_SNIPIT_VCF_SKIP`
   正则显式跳过 + 参考序列查找天然排除双保险。

### 6.3 接入设计（已实现，kv_variant_plot.py）

```
<out>/vcf/{gid}.vcf（BcftoolsCaller 全部样本 BAM 一起 calling → 多样本 VCF）
   ↓ build_snipit_alignment()：按基因型把 REF/ALT 套回参考序列
     （SNV/MNP 直接套换；纯缺失以 '-' 补齐跨度→snipit 读作 deletion；
      含插入的变异行整行跳过——定长伪比对表达不了插入；
      无基因型 ./. 记参考碱基，与共识序列口径一致）
伪比对（无需 MAFFT，参考坐标系天然对齐）
   ↓ snipit CLI（--sort-by-mutation-number --show-indels --solid-background）
<out>/variant_plots/{gid}_snipit.{png,pdf} + {gid}_snipit.snps.csv
```

- 不重比对的理由：VCF 本身在参考坐标系（同 collect_bams 多样本口径），
  套回即对齐；MAFFT 会引入与 VCF 无关的差异，反而和其余三张图不一致。
- 样本列 <2（含 sites-only）跳过：单样本场景 variant_map 已覆盖。
- snipit 未安装静默跳过（adjustText 同款容错），manifest 记 `snipit` 段。
- web 端 `tool_results.py` 平铺扫描 `variant_plots/*.png`，新图自动展示。

### 6.4 端到端实测

CMV RNA3（`examples/example_recomb_set.fasta`，10 条）真实分离物
法国/德国/斯洛文尼亚基因型写成多样本 VCF → `run_snipit_plots`：
3 样本 × 44 位点 → png/pdf/csv，图面与 MAFFT 直比对版完全吻合
（法国 C@104/106、德国 T@106/A@107、斯洛文尼亚 del@2 等模式一致）。

示例运行数据同步升级：`dev_tools/make_variant_demo.py` 的合成 VCF 从
sites-only 改为 S1/S2/S3 三样本基因型（固定种子），demo 目录三条参考
各出 snipit 分型图（41-43 位点/条，插入行按设计跳过并计数）；
无编码注释的参考 snipit 照常出图（分型图不依赖 SnpEff 注释）。

### 6.5 真实公共数据实测（2026-10-13）

数据：ToBRFV 爆发研究 PRJEB42296（约旦）两个真实 Illumina paired 样本
**ERR5032739 / ERR5032741**（NCBI SRA，sracha 经 s3-direct 下载，~1GB/对），
参考 ToBRFV RefSeq NC_028478.1（efetch FASTA+GBK，4 CDS）。
链路全程真实：步进抽样 ~30k 对 → minibwa 比对（病毒占比 94.2% / 87.9%）
→ bcftools 多样本 calling → SnpEff → run_snipit_plots。
结果：**22 条真实变异、2 样本 × 22 位点分型图**，共享位点（爆发克隆系
对 RefSeq 的系内差异）与样本私有位点清晰分列。

实测暴露并修复两个缺陷：

1. **老记录蛋白校验卡死建库**：NC_001367（TMV U1，1982）复制酶含通读
   终止（leaky stop），SnpEff 蛋白校验报 `'*' : Error` → 
   `Database check failed`（rc=4294967295），整库拒建。已在
   `SnpEffRunner.build` 加降级：失败且报该错时用 `-noCheckCds
   -noCheckProtein` 重建一次（旗标是单数 Protein，复数会报用法错误）。
2. **样本行名是 BAM 全路径**：bcftools 无 @RG/SM 时用 BAM 完整路径当
   样本列名，snipit 行名成了一长串 `C:\...\ERR5032739.bam`。已在
   `build_snipit_alignment` 写 FASTA 头时清洗（去目录、去 .bam 后缀、
   撞名加序号）。

选样教训（供后续复用）：ENA `tax_tree` 能按病毒 taxid 检索 runs 的
只有少数病毒（TMV 6 个、ToBRFV 217 个、PVY 300+）；NCBI SRA
`efetch rettype=runinfo` 实测返回空表不可用；宿主深 RNA-seq
（polyA 文库）里病毒 reads 趋近于零（TMV"感染组"实测 0% 比对），
要选病毒富集/爆发克隆系研究——克隆系内株间差异 <0.1%，任意两样本
reads 可互比，且保留真实样本间微多样性。sracha 支持 NCBI
accession 直下（8 连接 s3-direct，1GB 约 1.5 分钟）。

扩展回归（同日第二轮）：① 真编排入口
`known_virus_suite.py variant`（含 variant_evo/SNPGenie）与 snipit
共存无冲突，全流程 22.5s——SNPGenie 链路自带多样本 VCF 样本列清洗
（2 列 -> SAMPLE）；② 第三个阳性真实样本 ERR5032747（87.9%）拼入后
3 样本 × 22 位点图直接读出亚系结构：ERR5032747 与 5032741 变异谱
完全一致（同亚系），5032739 独有 6 处差异（另一亚系）；
③ 同研究其余候选（2763/2764/2767/2783/2784）病毒占比 0-1.5%
（"Nextstrain update" 样本实为低占比宏基因组），弃用。
