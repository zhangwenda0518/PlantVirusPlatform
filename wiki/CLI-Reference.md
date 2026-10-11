# CLI-Reference · 命令行完整参考

> 返回 [Home](Home) ｜ 上一页：[Settings-Ops](Settings-Ops) ｜ 下一页：[FAQ](FAQ)

入口：`python main.py <子命令> [参数]`。GUI 能做的 CLI 都能做（共用同一引擎）。

---

## 建库与初始化

### `init-taxonomy` — 下载/准备 NCBI Taxonomy（首次必做）

```bat
python main.py init-taxonomy [--force]
```

### `build-host-db` — kunpeng 宿主数据库构建

```bat
python main.py build-host-db --genome host-db\genome.fa --taxid 4081 ^
       [--hash-capacity 256M] [--threads 0] [--out-dir DIR] [--clean-mid] [--rebuild]
```

| 参数 | 说明 |
|---|---|
| `--genome` | 宿主基因组 FASTA（必填） |
| `--taxid` | 宿主 NCBI TaxID（必填，建库前自动校验） |
| `--threads 0` | 线程（0=自动 8；内存不足用 4） |
| `--out-dir` | 输出库目录（默认 `host-db/<taxid>_<源目录名>_host_db`） |
| `--clean-mid` | 建库后清理中间文件（省 ~3GB/库） |
| `--rebuild` | 重建（替换式，自动清理旧 library） |

### `build-virus-db` — kunpeng 病毒数据库构建

```bat
python main.py build-virus-db --fasta final.cluster.ref.fasta --info final.cluster.ref_info.tsv
       [--hash-capacity 64M] [--threads N] [--rebuild]
```

## 样品分析

### `analyze` — 样品全流程分析

```bat
python main.py analyze --r1 R1.fastq.gz [--r2 R2.fastq.gz] --sample NAME ^
       [--stages subsample,fastp,fq2fa,host,kvsuite,assembly,hostana,orf,orfa,phylo,primer,gbdraw,report] ^
       [--db-host DIR] [--db-virus DIR] [--threads N] [--confidence F] ^
       [--assembly-mode rnaviral|metaviral|meta|rna|isolate] ^
       [--subsample N] [--min-contig-len 500] [--min-orf-aa 100] [--memory 64] ^
       [--top-n-refs 10] [--tree-tool fasttree|raxml-ng|nj] ^
       [--tree-sampling blast|macro|genus|lineage] [--ncbi-refs NAME[,NAME2]] ^
       [--primer-mode conserved|plain] [--specificity] ^
       [--no-fq2fa] [--gbdraw-max 12] [--gbdraw-fasta x.fa --gbdraw-ann x.gff|x.gb] ^
       [--plot-engine auto|gbdraw|dfv] [--force]
```

- 阶段自由组合，`.done` 断点自动跳过，`--force` 强制重跑；
- `--subsample N` 只取前 N 对 reads（快速验证）。

### `report` — 重新生成可视化报告

```bat
python main.py report --sample NAME [--force]
```

### `host-analysis` — 对已有 ③ 结果做 ICTV 宿主预测

```bat
python main.py host-analysis --sample NAME [--force]
```

### `orfa` — ⑥b ORF 功能注释

```bat
python main.py orfa --sample NAME [--threads N] [--libs prot,pfam,bfvd] [--force]
```

## CRISPR 检测设计

### `crispr-design` — Cas12a/Cas13a 检测位点与引物设计（2026-10-06 新增）

```bat
python main.py crispr-design -i isolates.fasta [--near close_relatives.fasta] ^
       [--no-host] [--no-group] ^
       [--amp RPA|PCR] [--profile strict|relaxed] ^
       [--nuclease cas12a|cas13a] [--spacer-len N] [--probe none|exo|nfo] ^
       [--top-n 20] [--flank 250] [--margin 25] ^
       [--product-min N] [--product-max N] ^
       [-o run/tool_runs/crispr_cli_<时间戳>] [-t N]
```

- `-i` 目标 FASTA（多分离物；header 含 `RNA1` / `S2` 等标签时**自动按片段分组**，`--no-group` 关掉）；
- `--near` 近缘病毒 FASTA（特异性排除用）；宿主基因组筛查默认启用 `host-db`，`--no-host` 跳过；
- `--profile strict`（默认）为通用口径，`relaxed` 对高变异 RNA 病毒放松；
- `--spacer-len 0` 表示按酶档案默认（Cas12a 20 / Cas13a 30，允许区间 20–24 / 27–30）；
- `--product-min/max 0` 表示用档位默认产物长度区间；
- 产物含 crRNA 候选表、引物对、评分与订购信息，落盘在 `-o` 指定目录（默认 `run/tool_runs/crispr_cli_<时间戳>`）。

## 参考序列

### `ncbi-dl` / `ncbi-list` — NCBI Entrez 参考下载（进化树扩充）

```bat
python main.py ncbi-dl "Tobamovirus[ORGN] AND complete genome[TITL]" -n 集合名 ^
       [--db nucleotide|protein] [--max 100]
python main.py ncbi-list
```

### `gb-dl` / `gb-import` / `gb-list` / `gb-check` — GenBank 集合

```bat
python main.py gb-dl -n 集合名 [--term "Potyvirus[ORGN] AND complete genome[TITL]" | --acc acc1,acc2] [--max 50]
python main.py gb-import -n 集合名 --files a.gb,b.gb
python main.py gb-list
python main.py gb-check -n 集合名
```

### 建树参考库

```bat
python main.py ref-status                                    :: 参考库版本与统计
python main.py ictv-status                                   :: ICTV VMR 库状态
python main.py ictv-update [--xlsx VMR.xlsx]                 :: 更新 VMR
python main.py ictv-refs [--genus G] [--family F] [--species S] ^
       [--db plant|ictv] [--limit 20] [--any] [--download]   :: 选参（缺的按需下载）
```

## 公共数据

### `meta-search` — SRA+GSA 双引擎检索

```bat
python main.py meta-search --species "Lycium chinense" ^
       [--source TRANSCRIPTOMIC|GENOMIC|All] [--db sra|gsa|both] ^
       [--out DIR] [--no-detailed] [--ncbi-api KEY] [--deepseek-api KEY]
```

### `meta-info` — Run 列表 → 统一元数据（Core14/Full）

```bat
python main.py meta-info --runs <列表文件|单个Run> [--out DIR] [--mode local|api|both] ^
       [-t 4] [--fill-date] [--ncbi-api KEY] [--deepseek-api KEY]
```

### `meta-plot` — 元数据 SCI 可视化

```bat
python main.py meta-plot --input <检索表|Core14|Full.csv> [--out SCI_Figures_Output]
```

### `host-genome` — 宿主参考基因组下载（NCBI taxon/accession / gget·Ensembl）

```bat
python main.py host-genome --species "Lycium barbarum" [--source ncbi|gget] ^
       [--out DIR] [--release 60] [--taxid 112863] [--include-organelles] ^
       [--ncbi-api KEY] [--min-length N] [--verify-only] [--skip-datasets]

rem accession 精确下载（datasets，支持染色体过滤；taxid 从报告自动解析）
python main.py host-genome --accession GCF_000001405.40 --chromosomes X,Y

rem gget/Ensembl Plants 源（补 NCBI 覆盖不到的植物宿主；需先 pip install "gget==0.30.5"）
python main.py host-genome --species "Oryza sativa" --source gget
```

### `datasets` — NCBI datasets CLI 全模式下载（genome / gene / taxonomy / virus）

```bat
python main.py datasets genome --taxon "bos taurus"
python main.py datasets genome --accession GCF_000001405.40 --chromosomes X,Y --include genome
python main.py datasets gene --gene-id 672
python main.py datasets gene --symbol brca1 --taxon "mus musculus"
python main.py datasets gene --accession NP_000483.3
python main.py datasets taxonomy --taxon human,sars-cov-2
python main.py datasets virus --taxon sars-cov-2 --host dog
python main.py datasets virus --protein S --host dog --out SARS2-spike-dog
python main.py datasets virus --taxon "tomato brown rugose fruit virus" --name tobrfv

rem virus 模式默认落平台数据集布局（virusref_db/<名>/reference.fasta +
rem reference.ref_info.tsv）→ 数据库构建页「原始数据集」卡即时识别、一键建库；
rem 加 --no-dataset 存通用目录 datasets_dl/。

### `datasets-setup` — 安装/更新 NCBI datasets CLI

```bat
python main.py datasets-setup
rem 从 NCBI 官方（域名白名单内）下载 win64 exe 到 3rd/tools/datasets/，
rem PATH 与平台内置位置都会被自动探测；旧版 --exclude-* 旗标自动翻译为 --include
```

## LOGAN 溯源

```bat
python main.py logan-create --name 查询名 [--sample S1 | --fasta x.fa] ^
       [--contigs c1,c2] [--segments 1-4]
python main.py logan-batch --name 查询名 --email a@qq.com[,b@qq.com] ^
       [--group Fast_No_human|All|All_No_viral_human|Fast|Fast_No_RefSeq|Transcriptomic|Metatranscriptomic|Metagenomic|GenBank_RefSeq] ^
       [--show-browser] [--first-wait 300] [--max-wait 1800]
python main.py logan-import --name 查询名 --segment 1 --result 结果表.csv
python main.py logan-jobs
```

## NCBI 提交准备

```bat
python main.py submit-list
python main.py submit-init --name 项目名 [--taxonomy t.tsv] [--metadata m.csv] ^
       [--import-csv x.csv] [--demo] [--authors "Last, First"] [--title "..."] ^
       [--bioproject PRJNA...] [--host ...] [--lat-lon "38.47 N 106.27 E"] ^
       [--sequencer "Illumina NovaSeq 6000"] [--assembler "SPAdes;4.3.0;rnaviral"] [--coverage 42.5x]
python main.py submit-validate --name 项目名
python main.py submit-fill --name 项目名 --column 列 --value 新值 [--old 旧值]
python main.py submit-export --name 项目名 [--assembler ...] [--sequencer ...] [--enrichment "rRNA depletion"]
python main.py submit-adopt --name 项目名 (--out-dir DIR | --from-server DATASET) ^
       [--samples 样本表] [--dataset 名] [--dry-run] [--no-export] [--min-length N] ^
       [--geo-loc "China:Ningxia"] [--collection-date YYYY-MM-DD] [--host ...] [...]
python main.py submit-sbt --name 项目名 --last Zhang --first Wenda ^
       --affil 机构 --city 城市 --country 国家 --email me@x.com [--div] [--sub] [--title]
```

## 平台运维

```bat
python main.py tools                        :: 工具探测状态
python main.py selfcheck                    :: 环境自检（未就绪时非零退出）
python main.py db-migrate --to DIR [--mode copy|move] [--dry-run] [--check]
python main.py tool-runs status|organize|archive|clean ^
       [--before YYYYMMDD] [--older-days N] [--active-days N] [--archive-days N] [--dry-run]
python main.py samples-backfill [--sample 名]... [--dry-run]   :: 老样品回填清单
python main.py importprobe <模块清单>        :: 内部命令（selfcheck 子进程用）
```

## 全局约定

- 路径：相对路径以平台根为基准；产物一律限定平台目录内（`run/`）；建库/分析输入支持平台外绝对路径；
- 旧路径（`results/`、`tool_runs/`、`tools/`、`bin/`…）自动重定向，无需改写；
- `--force`：忽略阶段 `.done` 断点强制重跑（各子命令通用语义）。

## 相关

- GUI 操作对应 → [Pipeline-Overview](Pipeline-Overview) 与各模块页
- 环境自检细则 → [Settings-Ops](Settings-Ops)
