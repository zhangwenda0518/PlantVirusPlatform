# Databases · 数据库与参考库

> 返回 [Home](Home) ｜ 上一页：[Logan](Logan) ｜ 下一页：[Settings-Ops](Settings-Ops)

平台数据库总体 ~3.4GB，全部在 `databases/`（宿主数据独立在 `host-db/`，~1.5GB/物种）。数据库页负责 Taxonomy 下载与宿主库构建；其余库平台预置或按需自动建。

---

## 总览

| 目录 | 内容 | 体积 | 来源 |
|---|---|---|---|
| `tax_db/` | NCBI Taxonomy：nodes/names/merged.dmp | ~545MB | 数据库页下载（~57MB 压缩包） |
| `kunpeng_db/` | kunpeng 分类库：plant / ref / rvdb | ref 680MB | 预置 |
| `virusref_db/` | 已知病毒**鉴定库** kv_index/（自包含：reference.fasta + ref_info + salmon_k31 索引）+ blast/ 本地检索库 | ~850MB | 预置 |
| `annot_db/` | 注释库：prot/（RefSeq 病毒蛋白 ~550MB）· hmm/pfam/（Pfam 病毒子集 ~119MB）。（cdd/ CDD 结构域库已于 2026-09-29 移除：prefilter 索引 ~5GB 小内存机器反复 OOM；病毒域注释改用 misc_db/suvtk/ 内的 BFVD 库，与 SUVTK 提交链共用） | ~0.67GB + BFVD 212MB | prot 首次运行自动下载；pfam/BFVD 预置 |
| `tree_db/` | 建树参考库：**plant_tree.db**（植物口径预下载）· **ictv_tree.db**（ICTV VMR 元数据 + 按需序列缓存） | ~37MB | 预置 + ictv-update |
| `misc_db/` | host_prob（宿主概率表 5.9 万条）· suvtk · viroids | — | 预置 |
| `virus_ref/`（**可选**） | 非冗余植物病毒参考库：8,465 条 98% ANI 聚类代表 + 5,773 条完整基因组 + ~199K 全量 + DATA_VERSION | 自备 | 外部获取后整体放入即启用 |

## Taxonomy（tax_db）

- 数据库页一键下载 / `python main.py init-taxonomy`（new_taxdump.tar.gz → 解压 nodes/names/merged）；
- 用途：宿主判定谱系、④/⑥b 的宿主类别归类、contig 深度注释 8 级谱系统列、TaxID 校验（merged.dmp 处理旧编号并回新编号）。

## 病毒鉴定与分类库（virusref_db + kunpeng_db）

- **kv_index 自包含**：reference.fasta + reference.ref_info.tsv（VMR_Family/VMR_Genus 列供 ②b 分类树）+ salmon_k31 索引；②b 识别/定量、③ BLASTN 回退、⑦ 参考池兜底共用；
- 重建：`build-virus-db` 或 `dev_tools/build_virus_db.py`；
- 更新：新版 4 件套（两个 FASTA + 两个 Info.tsv + DATA_VERSION）整体替换，平台按 DATA_VERSION 自动重建元数据缓存与 BLAST 库。

## 建树参考库（tree_db）

### plant_tree.db（plant 口径，预下载）

植物病毒参考库 **48 科 / 452 属 / 6,167 条**（45 科来自 acvirus 谱系 6,150 条 + 补齐 Ourmiaviridae/Ambiguiviridae/Pestiviridae）；`plant_meta.tsv` 的 `Seq_Source` 列标明来源（ictv_db / ref_db）。⑦ 的 `--tree-sampling macro/genus/lineage` 从这里取参考；植物库未收录的由选参层联网补齐（单科缺口通常个位数）。

### ictv_tree.db（ictv 口径，按需）

**ICTV VMR 参考库**（`ictv_db.py`）：官方 [VMR 当前版 xlsx](https://ictv.global/vmr/current?fid=15873) 解析全病毒分类元数据（MSL41：22,785 accession / 4,068 属 / 393 科）+ 按需下载的参考序列缓存。分层架构：

1. **解析层**：VMR xlsx → `taxa.txt`（谱系列 + 病毒名/基因组完整性/Baltimore/宿主组扩展列，多 accession 拆分）+ DATA_VERSION；
2. **选参层**：按属/科/种过滤，本地已有（plant_tree.db 命中 > gb 缓存）优先、Complete genome 优先；
3. **下载层**：缺的 accession 经 NCBI efetch 落 `gb_cache/<acc>.gb`（幂等缓存，保留特征注释）并汇总 `gb_refs.fa`；
4. **兜底**：⑦ 只读接入——谱系并入分类索引、缓存序列并入参考池；**分析中不联网**，未下载自动回退 plant_tree.db。预览列 Source 显示命中库名。

```bat
python main.py ictv-status                  :: 库状态（MSL 版本/覆盖）
python main.py ictv-update                  :: 在线下载 VMR → 解析
python main.py ictv-update --xlsx VMR_MSL42.xlsx
python main.py ictv-refs --genus Tobamovirus            :: 选参预览
python main.py ictv-refs --genus Nepovirus --download   :: 缺的按需下载
python main.py ref-status                   :: 参考库版本与统计
```

> ictv 口径**不再预下载**全病毒界序列（原 559MB all_virus.fasta 已按设计移除）。

## virus_ref（可选自备数据源）

非冗余植物病毒参考库（源自 plant_virus_db_pipeline）。目录存在即启用（无需改配置）：

- **④** 用其 ICTV 谱系 + 已知宿主交叉验证（预测宿主 ∉ 已知宿主类别 → host_check=WARN）；
- **⑦** 参考池优先 RefSeq/完整基因组；`ref_meta.tsv` 为规范化元数据缓存（首次使用自动构建：Segment 归一 + 本地 taxonomy 谱系补全，ICTV 种覆盖 10%→100%；`segment_norm.py` 把 299 种原始 Segment 写法收敛到受控词表）；
- 不存在时平台照常运行，各调用点回退自带数据（`ref-status` 会提示「参考库未部署」）。

`universal_ref.py` 另支持通用病毒库（ref-virus / RVDB 口径）建库，按需重建 `universal/` 子目录。

## 宿主库（host-db，独立于 databases/）

- 按物种一库：`host-db/<TaxID>_<源目录名>_host_db/`，库内 `host_db.json` 清单记录目标物种；多宿主库共存，`active_host_db` 指定当前（页面一键切换）；旧布局固定槽位 host/classify/ 仍识别（legacy）；
- **建库要点**：TaxID 建库前校验；FASTA 头自动规范化注入 taxid；**大基因组自动切片**（≤1MB 片段、相邻重叠 34bp，k-mer 数与整条建库完全一致——2GB 内存完成 1.8GB 基因组建库）；替换式重建；
- 建库耗材（library/prep/taxonomy 副本）可清：勾「构建后清理中间文件」或 `--clean-mid`，省 3-4GB/库，重建时自动再生成；
- 自备宿主库放任意位置：构建页「📂 加载已构建宿主库」选中即登记路径白名单（extra_db_roots）。

## 宿主概率表（misc_db/host_prob）

`{species,genus,family,order}_host_probability.tsv` 共 **5.9 万条目**，源自 ICTV 分类 × 宿主记录交叉统计——④ 宿主预测的级联查表基础（[Host-Prediction](Host-Prediction)）。

## 注释库（annot_db）

- **prot**：NCBI RefSeq 病毒蛋白全量（首次自动下载 ~107MB 并建库）；换库：替换 `viral_prot.faa` 后删 `organism_tax.tsv`/`db_info.json` 自动重建；
- **hmm/pfam**：Pfam-A-Viruses.hmm（1,074 profiles，自描述 DESC）；缺失时 HMM 层整体跳过；
- **bfvd**：mmseqs 格式 BFVD 病毒域库（`misc_db/suvtk/bfvd`，212MB，与 SUVTK 提交链共用；⑥b 注释层2 与 ③b 验证 B 路共用，product 命名同 `bfvd_uniprot_names.tsv`）。原 cdd（NCBI Cdd，~5GB 索引反复 OOM）已于 2026-09-29 移除。

## 数据库迁移与外置

- `python main.py db-migrate --to D:\库目录 [--mode move] [--dry-run] [--check]`：robocopy（多线程+断点续传）→ 文件数与总字节校验 → **通过后才切配置**；任一失败配置不动、源原样，可断点重跑；目标须在平台目录外且剩余 ≥ 源体积 +10GB；
- 设置页「数据库目录」直接填路径（即时反馈宿主库/病毒库/Taxonomy 就绪状态，切换后 ictv_db/virus_ref/universal_ref/local_search/verify 路径常量同步刷新，无需重启）；
- 从其他机器迁移：整体拷贝 `databases/`（与 `host-db/`）到目标机同位置即可，探测结果不写回 platform.json 保持可移植。

## 相关

- 各库在分析里的用法 → [Kvsuite](Kvsuite) / [Host-Prediction](Host-Prediction) / [ORF](ORF) / [Phylo](Phylo)
- 数据库迁移与外置 → [Settings-Ops](Settings-Ops)
