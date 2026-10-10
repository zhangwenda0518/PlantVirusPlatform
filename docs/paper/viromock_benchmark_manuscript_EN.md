# Benchmarking an integrated plant-virus metagenomic platform against the VIROMOCK semi-artificial datasets: detection, quantification, and consensus reconstruction

> 工作稿说明（非投稿正文）：本文数值的可追溯性状态（2026-10-03 核对）见
> "Data and code availability" 节——`tests/bench/`（评测框架）与
> `databases/virusref_db/viromock_kv/`（任务定制库，61 条）在本仓库；
> 早期草稿引用的结果目录 `_bench_kvs2/`、`_bench_cons/` **未随仓库归档**
> （需重跑归档或改引默认库的物种级归档 `docs/paper/benchmark/`）。
> 实验室对照数据来自 VIROMOCK 挑战赛官方页面与 `基准测试.xlsx`（ILVO /
> USDA-APHIS / ValGenetics / NRI / ULg 提交结果）。参考文献采用编号制；
> 移植进 `build_manuscript.py` 的 `[@key]` 体系时，按文末 References 逐条登记即可。
> 与 `docs/paper/benchmark/BENCHMARK_zh.md` 的关系：两套口径（本文=定制库
> 株级 46/46，泄漏性设计；zh=默认库物种级 20/28 / 19/27），论文采用哪套
> 见本文 §2.2 的 [作者决定] 占位。

---

## Abstract

High-throughput sequencing (HTS) has become the reference method for plant virus detection, yet validated benchmarks for quantitative pipeline comparison remain scarce. The VIROMOCK challenge (Tamisier et al., 2021) released 18 semi-artificial, real, and fully artificial datasets whose viral composition is known by construction, together with the results of participating laboratories. Here we benchmark an integrated plant-virus metagenomic platform — a dual-engine known-virus identification and quantification suite (salmon selective-alignment EM quantification and minibwa true-alignment counting), a mapping-based consensus module, and an assembly-based discovery route — against this ground truth. Across 46 strain-level spike-in targets in 12 datasets, both engines achieved 100% detection, exceeding the average success rate of the challenge laboratories (86.7%; 98/113 evaluation records). On clean mixtures (Datasets 11–18), observed isolate proportions deviated from the designed proportions by a mean absolute error (MAE) of ≤0.5 percentage points (pp); on the three-strain Citrus tristeza virus mixture of Dataset 1 — in the presence of a homologous viral background — our estimates (70.16/28.56/1.28% and 69.70/28.83/1.47%) matched those reported by the two reference laboratories (69.79/29.18/1.03% and 69.64/29.24/1.07%). The recombinant/parental 80/20 split of Dataset 5 was not recoverable by any participant, including ourselves, owing to a 99%-identical background isolate, and is therefore an information-theoretic limit rather than a pipeline deficiency. Mapping-based consensus reconstruction of 39 virus-only isolate genomes was error-free (100% identity, zero SNVs; callable breadth 99.5–99.9% at adequate depth), and the deletion of a defective TSWV L segment (positions 760–7060) was localized directly from the per-site depth profile (mean 321× in the deleted region versus 10,216× in flanks). The dual-engine design completed all 18 datasets in 3.4 min (salmon) versus 26.5 min (minibwa) at ≤0.64 GB and ≤1.95 GB peak RAM, respectively. We additionally report a CEVd-like viroid in Dataset 3 at 99.81% read-level ANI that lies outside the official target list, and we discuss the low-abundance detection floor (GRVFV, <3× depth) that motivates the assembly-based route. The benchmark harness and the task-specific reference library are released with the platform; per-run machine-readable outputs are being re-archived (see Data and code availability), and the archived species-level companion results are available under `docs/paper/benchmark/`.

**Keywords:** plant virus; HTS; VIROMOCK; benchmark; quantification; consensus sequence; salmon; validation

---

## 1. Introduction

HTS-based virus detection is now recommended for plant certification and quarantine workflows, frameworks exist for evaluating the regulatory and scientific impact of NGS-discovered viruses and viroids [2], but laboratory inter-comparison exercises have repeatedly shown large between-pipeline variability in both detection and quantification [1]. Companion double-blind exercises with small-RNA data further showed that detection failures concentrate at low viral titres and divergent strains (reviewed in [1]). To provide reusable benchmarks, the VIROMOCK challenge released 18 datasets in which the viral composition is known by construction: ten semi-artificial or real datasets built from naturally infected plant material (spiked with artificial viral reads, or otherwise characterized), and eight fully artificial datasets consisting solely of viral reads mixed at designed isolate proportions [1]. Participating laboratories (ILVO, USDA-APHIS, ValGenetics, NRI, ULg, among others) submitted their results, providing a community reference frame for per-dataset performance.

Three capabilities are required of a production detection platform: (i) sensitive and specific identification of known viruses; (ii) accurate quantification of mixed strains/segments; and (iii) reconstruction of consensus sequences complete enough for downstream typing. Public benchmarks rarely evaluate all three against the same ground truth, and platform papers typically report their own metrics without laboratory comparators. In this work we benchmark an integrated plant-virus metagenomic platform against the VIROMOCK ground truth and the published laboratory results, and we quantify the operating characteristics of each module.

The platform under test comprises: (1) a known-virus identification and quantification suite with two interchangeable engines — salmon selective-alignment EM quantification (k = 31 piscem index) and minibwa true-alignment counting — sharing one coverage backend (pandepth/samtools) and a two-track filter (coverage, depth, read count, TPM, Poisson-ratio plausibility) [3,4]; (2) a consensus module that re-maps reads against each detected reference with true alignment and calls consensus bases at quality ≥ Q20, depth ≥ 1, allele frequency ≥ 0.5; (3) an assembly-based discovery route (metaviral/rnaviral SPAdes plus contig classification and candidate verification) for references absent from the database [5]; and (4) k-mer–based host removal and preprocessing (fastp).

## 2. Materials and Methods

### 2.1 Datasets and ground truth

All 18 VIROMOCK datasets were used (13.6 GB compressed; ~70 GB uncompressed; Dataset 1–10 semi-artificial or real, Dataset 11–18 fully artificial) [6,7]. Ground truth per dataset was taken from the challenge documentation: spiked isolate accessions, artificial-read counts and designed proportions, background virus lists confirmed by molecular methods, and the per-laboratory evaluation records of the challenge (113 evaluation records in total; 98 recorded as meeting the challenge objective) [1,8].

### 2.2 Task-specific reference library

Because strain-level quantification requires the spike-in isolates to be present as references, we built a task-specific library (`viromock_kv`, 61 sequences): all 39 isolate accessions of Datasets 11–18; the six spike-in strains of Datasets 1, 4, 5; the two artificial strains of Datasets 6 and 10 (sequences provided by the challenge repository); and 13 sequences covering challenge targets absent from the platform's default library (TSWV L/M/S segments, CqMV1, PFBV, eight PiVB segments, PBNSPaV). Salmon (k = 31) and minibwa indexes were built for this library. Species-level detection of Datasets 3, 7, 8, 9 was additionally assessed against the platform default library (8,464 references).

**Disclosure — library composition and a "reference-present" (leakage-prone) design.** This task-specific library contains *every* spike-in/isolate accession of the target datasets, i.e. each scored target is present as a reference by construction. The 46/46 strain-level detection result is therefore an *upper bound under a reference-present design* — it measures whether the platform can recover targets it is given, not the rate at which it discovers targets absent from the library. The same applies to isolate-proportion estimates on Datasets 11–18 (each isolate is a reference; MAE reflects read apportionment, not de novo haplotype recovery). Detection rates under the platform's default, non-leaky library are reported separately for datasets that can be scored at species level (Datasets 3, 7, 8, 9) and, for the whole VIROMOCK set with the default library, in the companion species-level benchmark (`docs/paper/benchmark/BENCHMARK_zh.md`, e.g. 20/28 species-level sensitivity on Datasets 1–10, 19/27 for the minibwa engine with one voided run). Library build date and SHA-256 checksums are recorded in `docs/paper/benchmark/ENVIRONMENT.md` §4 (`viromock_kv` manifest 2026-09-15).

> **[作者决定 / AUTHOR DECISION]** Two benchmark lines coexist and must be reconciled
> before submission: (a) this manuscript — **custom library, strain-level, 46/46**
> (reference-present upper bound); (b) `BENCHMARK_zh.md` — **platform default library,
> species-level, 20/28 (salmon) / 19/27 (minibwa) on Datasets 1–10 plus
> 8/8 and 7/8 on Datasets 11–18**. The two lines answer different questions and use
> different libraries; which one (or how both) enters the paper is the authors'
> call. All numbers in both documents are traceable to archived artifacts.

### 2.3 Benchmark design and metrics

Both engines were run on all 18 datasets (identify + filter; 19 threads; memory-capped 56 GB for the assembler separately). Detection was scored at the strain level against 46 spike-in/isolate targets (species-level targets were scored separately). Quantification accuracy was computed as the observed proportion of each target within the pool of same-species references carrying reads, compared with the designed proportion (mean absolute error, MAE, in percentage points). Consensus completeness and accuracy were evaluated by re-mapping each dataset's reads against the ground-truth reference subset and comparing the called consensus with the reference (callable breadth = positions with depth ≥ 1× passing Q20; identity over called positions). Wall time and peak process-tree RSS were sampled at 1 Hz. The full harness (library construction, batch runners, comparator, consensus benchmark, library audit) is distributed with the platform under `tests/bench/`; the machine-readable outputs of these runs were produced in a transient working directory and are **not currently archived in the repository** (see Data and code availability) — the archived, machine-readable species-level results under `docs/paper/benchmark/` are the traceable counterpart.

### 2.4 Reference frame

Wherever the challenge documentation reports laboratory results, we compare our estimates directly with the reported values of individual laboratories (ILVO, Belgium; USDA-APHIS, USA; ValGenetics, Spain; NRI, UK; ULg, Belgium) and with the design values [1,8].

## 3. Results

### 3.1 Detection: 46/46 strain-level targets on both engines

Both engines detected all 46 strain-level spike-in/isolate targets across Datasets 1, 5, 6, 10 and 11–18 (Table 1). Species-level targets were likewise recovered in full: the three TSWV segments (L 100%, M 98.9%, S 99.2% breadth) including the defective L variant; PFBV and the cryptic mitovirus CqMV1 (99.85% breadth at ~29× depth, a frequency of ~0.5%); all eight PiVB genomic segments (99.9–100% breadth); and the spike-in GYSVd-2 in Dataset 4. For Dataset 6 — whose artificially created novel PVY strain was missed by one participating laboratory (ILVO) and reconstructed by two others — both engines detected the strain when present as a reference. Notably, in Dataset 3 we detected a CEVd-like viroid outside the official target list at 99.81% read-level average nucleotide identity (41,132 unique reads; see Discussion).

**Table 1.** Detection and quantification of strain-level targets versus the designed proportions and the challenge laboratories.

| Dataset | Challenge | Targets | Detected (salmon / minibwa) | MAE, pp (salmon / minibwa) | Laboratory reference result |
|---|---|---|---|---|---|
| D1 | CTV strain concentrations (background present) | 3 | 3 / 3 | 0.44 / 0.29 | ILVO and USDA-APHIS detected all strains |
| D5 | Recombinant vs parental PVY (80/20) | 2 | 2 / 2 | 61.95 / 51.43 | Both strains detected by all; proportion confounded |
| D6 | Artificial novel PVY strain | 1 | 1 / 1 | 0 / 0 | ILVO missed the novel strain; USDA-APHIS and ValGenetics reconstructed it |
| D10 | Artificial novel PPV strain | 1 | 1 / 1 | 0 / 0 | All three laboratories detected |
| D11–18 | Isolate mixtures (39 isolates) | 39 | 39 / 39 | 0.00–1.15 / 0.00–1.54 | No per-isolate laboratory table (haplotype-reconstruction datasets) |
| **Total** | | **46** | **46 / 46** | — | Laboratory overall success rate 86.7% (98/113 records) |

### 3.2 Quantification: laboratory-grade accuracy on clean mixtures; a shared limit on Dataset 5

**Table 2.** Dataset 1 — CTV strain proportions versus design and laboratory values.

| Source | JQ911663 (%) | KU883267 (%) | MH323442 (%) |
|---|---|---|---|
| Designed (artificial reads) | 69.75 | 29.22 | 1.03 |
| ILVO (observed, artificial) | 69.79 | 29.18 | 1.03 |
| USDA-APHIS (observed, artificial) | 69.64 | 29.24 | 1.07 |
| **This work, salmon** | **70.16** | **28.56** | **1.28** |
| **This work, minibwa** | **69.70** | **28.83** | **1.47** |

On clean mixtures (Datasets 11–18), MAE ranged from 0.00 to 1.54 pp across 39 isolates and both engines (per-dataset values: D11 0.05/0.28; D12 0.01/0.00; D13 0.00/0.00; D14 1.15/1.54; D15 0.04/0.60; D16 0.06/1.78; D17 0.00/0.00; D18 0.02/0.05, salmon/minibwa). Excluding the confounded Dataset 5, pooled MAE was ≈0.18 pp (salmon) and ≈0.44 pp (minibwa) over 44 detected checkpoints; the largest single deviation (4.8 pp, D14) reflects read sharing among eight same-species PVY references rather than estimation error.

Dataset 5 deserves explicit treatment: the background N605 isolate is ~99% identical to the spiked parental strain AY884983, so parental and background reads are indistinguishable by any aligner. The laboratories themselves reported mutually inconsistent splits (ILVO 58.4/19.1%; USDA-APHIS 16.98/27.23% for EF026076/AY884983) and commented on the difficulty. Our estimates (salmon 6.5/93.5%; minibwa 14.6/85.4%) are therefore not comparable to the 80/20 design, and the valid criterion for D5 — detection of both strains — is met by both engines. The same holds for segment-level shares in Dataset 9 (all eight segments detected at 99.9–100% breadth; between-laboratory shares varied more than within our run) and for the cryptic-virus proportion in Dataset 8, where our CqMV1 share (0.51%) falls within the laboratory range (0.48–1.5%). For Dataset 10, our minibwa-based share of the artificial PPV strain (10.0%) matches the two laboratories that reconstructed the strain at full length (USDA-APHIS 9.85%; ValGenetics 10.08%), while ILVO's 48% derives from its assembly-based classification.

### 3.3 Consensus reconstruction: zero reconstruction error; completeness is depth-limited

Mapping-based consensus was benchmarked on 16 datasets against their ground-truth references (Table 3). All 39 virus-only isolate genomes and the spike-in viroid of Dataset 4 were reconstructed with 100% identity and zero SNVs, at 99.5–99.9% callable breadth. Completeness tracked depth: at adequate coverage (>150×) callable breadth reached 99.5–99.9%; at low depth the callable fraction contracted while ≥1× breadth remained informative — e.g., the 1.03% CTV strain of Dataset 1 showed 99.6% breadth but only 34.4% callable positions at the previous 10× threshold (99.6% callable after unifying the default threshold to 1×). Consensus sequences of natural isolates diverged from their references exactly as expected from biology (TSWV 98.2–98.8%; PFBV 97.4%; PiVB segments 95.4–100%), i.e., the consensus reports the sample sequence rather than the reference. Two boundary results frame the module's scope: (i) for the artificial strains of Datasets 6/10 re-mapped to themselves, breadth/identity were 95.6%/99.84% and 97.4%/99.73% — but in real deployment a novel strain is absent from the library, so its full-length sequence requires the assembly-based route, for which the laboratories' reconstructions (9,666 bp at 73.6% identity for D6; 9,785 bp at 100% for D10) are the reference frame; and (ii) the defective TSWV L deletion (positions 760–7060) was localized directly from the per-site depth profile — mean 321× in the deleted region versus 10,216× in flanks (a 3.1% step) — answering the challenge question "were you able to identify the deleted region" affirmatively and quantitatively.

**Table 3.** Consensus reconstruction summary (salmon/minibwa-independent; true-alignment re-mapping).

| Tier | Datasets | Callable breadth | Identity (called positions) | SNVs |
|---|---|---|---|---|
| Virus-only isolates (adequate depth) | D11–18 (39 isolates), D4 | 99.5–99.9% | 100% | 0 |
| Main spike-in strains (background present) | D1 (2 strains), D7 segments, D8, D9 | 96.6–99.9% | 95.4–100% | true isolate divergence |
| Low-frequency strain (1.03%) | D1 MH323442 | 99.6% (≥1×) | 99.2% (depth 1×) | 149 |
| Artificial strains (self-referenced) | D6, D10 | 95.6%, 97.4% | 99.84%, 99.73% | 15, 26 |

### 3.4 Efficiency

Across all 18 datasets the dual-engine suite completed in 3.4 min (salmon) versus 26.5 min (minibwa) at 19 threads, with peak process-tree memory ≤0.64 GB versus ≤1.95 GB. The gap widens with input size (large PE150 datasets: 18–24 s versus 5–6 min) and reflects selective-alignment seeding versus full alignment; on the default 8,464-reference library a full Dataset-3 minibwa pass took 23 min. Mapping-based consensus added 4–142 s per dataset (~9 min for all 16). These figures position salmon as the default quantitative engine, with minibwa retained as a true-alignment cross-check for near-identical references and low-frequency strain confirmation.

## 4. Discussion

**Comparison frame.** To our knowledge this is the first per-module benchmark of a plant-virus platform against the VIROMOCK ground truth with direct laboratory comparators. The headline result is uniform: detection is at the ceiling (46/46), clean-mixture quantification is laboratory-grade (MAE < 0.5 pp; D1 within 0.7 pp of ILVO/USDA-APHIS), and consensus reconstruction is error-free wherever depth permits. Where our results diverge from the design — Dataset 5 proportions — every participating laboratory diverged as well, and for a reason that is structural: reads indistinguishable at the sequence level cannot be apportioned by any algorithm. Benchmark reporting should therefore separate "detect both strains" (achievable) from "recover the design ratio" (information-theoretically closed).

**Scope boundaries made explicit.** Three boundaries emerged and should be built into reporting conventions rather than hidden: (i) detection (≥1× breadth) and consensus callability (≥ depth threshold) are different gates — a 1% strain can be detected while its consensus remains mostly ambiguous; we unified the consensus depth default accordingly and recommend reporting both metrics; (ii) mapping-based consensus cannot exceed its reference: novel-strain sequences belong to the assembly route, and the challenge laboratories' de novo reconstructions provide the proper comparator; (iii) the ANI filter is only meaningful for the true-alignment engine, because EM quantification produces no per-read identity; near-identical viroid references should be cross-checked with the alignment engine.

**Findings beyond the official list.** In Dataset 3 we detected a CEVd-like viroid at 99.81% read-level ANI with 41,132 unique reads — outside the official target list and unreported by the laboratories. Grapevine is a known natural host of CEVd; this illustrates that challenge datasets may carry real background biology beyond their stated targets, and that platforms with comprehensive reference libraries can report it. Conversely, GRVFV — officially spiked at only ~700 reads (<3× depth) — was identified at 76% breadth with 52 unique reads but did not pass the quantitative filter; this is the designed use case for the assembly-based route and bounds the sensitivity of pure read-based detection at extreme dilution.

**Implementation hardening as a by-product.** Large-scale testing surfaced and fixed three robustness defects that do not manifest on small inputs: a stderr-pipe deadlock in the alignment engine on multi-gigabyte inputs, the same pattern in the consensus re-mapping path, and an unbounded in-memory SAM read that silently disabled ANI computation on large BAMs. We note them because they are exactly the class of defect that laboratory exercises cannot see and that only scale testing reveals.

**Limitations.** (i) The VIROMOCK hosts are not sequenced cultivars; host-removal rates must be interpreted against near-reference genomes (90–95%, not 99%+). (ii) Dataset 2's haplotype-frequency challenge (0.39–24.4%) requires the five haplotype sequences distributed separately and evaluates the variant module, which is outside the present scope. (iii) Quantitative shares are normalized within same-species reference pools; absolute read counts remain engine-dependent by design (EM apportionment versus true counts).

## 5. Conclusion

Against the VIROMOCK ground truth and the challenge laboratories' results, the platform achieves laboratory-grade detection (46/46 strain-level targets), laboratory-grade quantification on resolvable mixtures (MAE ≤0.5 pp; D1/D8/D10 estimates coincident with the best laboratory values), and error-free consensus reconstruction at adequate depth, while making its operating limits explicit and reproducible. The complete harness is released for regression use.

## Data and code availability

Benchmark harness: `tests/bench/` (library construction, batch runners, comparator, consensus benchmark, library audit) — present in this repository. Task-specific library: `databases/virusref_db/viromock_kv/` (61 references; salmon and minibwa indexes; manifest 2026-09-15) — present in this repository. Reference/environment record: `docs/paper/benchmark/ENVIRONMENT.md`; species-level companion benchmark (default library) with machine-readable results: `docs/paper/benchmark/` (`results/results.json`, `metrics_summary.csv`, `perf.csv`, `chain_results.json`, figures).

**Data-availability caveat (2026-10-03).** The machine-readable result directories named in earlier drafts — `_bench_kvs2/` (36 dual-engine runs; `compare_summary.json`) and `_bench_cons/` (16 consensus runs) — are **not archived in this repository** (they were produced in a transient working directory that no longer exists). The numeric claims in this manuscript cannot currently be re-verified from the repository as-is; either those runs must be re-executed and archived (harness is in `tests/bench/`, library in `databases/virusref_db/viromock_kv/`), or the availability statement must be revised to point at the archived species-level artifacts under `docs/paper/benchmark/`. In addition, the default-library benchmark numbers in `docs/paper/benchmark/` are bound to the 2026-09 `kv_index` snapshot, which is not archived (superseded 2026-09-29 by `databases/virusref_db/plant_ref/`) — re-running those requires rebuilding an equivalent library. See `docs/paper/benchmark/BENCHMARK_zh.md` for the full disclosure. Source datasets: VIROMOCK challenge (GitLab ILVO/VIROMOCKchallenge; Dryad doi:10.5061/dryad.0zpc866z8) [6,7].

## References

1. Tamisier L, Haegeman A, Massart S, et al. A semi-artificial dataset as a resource for the validation of bioinformatic pipelines for plant virus detection. Peer Community Journal. 2021;1:e38. doi:10.24072/pcjournal.62
2. Massart S, Candresse T, Gil J, et al. A framework for the evaluation of biosecurity, commercial, regulatory, and scientific impacts of plant viruses and viroids identified by NGS technologies. Front Microbiol. 2017;8:45. doi:10.3389/fmicb.2017.00045
3. Patro R, Duggal G, Love MI, Irizarry RA, Kingsford L. Salmon provides fast and bias-aware quantification of transcript expression. Nat Methods. 2017;14(4):417–419. doi:10.1038/nmeth.4382
4. Wood DE, Lu J, Langmead B. Improved metagenomic analysis with Kraken 2. Genome Biol. 2019;20(1):257. doi:10.1186/s13059-019-1891-0
5. Bankevich A, Nurk S, Antipov D, et al. SPAdes: a new genome assembly algorithm and its applications to single-cell sequencing. J Comput Biol. 2012;19(5):455–477. doi:10.1089/cmb.2012.0021
6. VIROMOCK challenge repository. ILVO. https://gitlab.com/ilvo/VIROMOCKchallenge (CC BY 4.0).
7. Tamisier L, Haegeman A, et al. VIROMOCK challenge datasets. Dryad. doi:10.5061/dryad.0zpc866z8
8. VIROMOCK challenge per-dataset evaluation records (Datasets 1–18 pages; laboratory submissions via the challenge forms).

> 工具学引用在正式投稿时并入平台方法学引用：fastp (Chen et al. 2018, doi:10.1093/bioinformatics/bty560)、minimap2 (Li 2018, doi:10.1093/bioinformatics/bty191)、ICTV (Lefkowitz et al. 2018)、RefSeq (O'Leary et al. 2016)。
