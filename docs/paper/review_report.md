# Peer Review Report — VirusPlatform manuscript (draft v1)

> **审稿执行摘要（中文）**：五位评审一致认为论文选题有实际价值、工程事实扎实（代码库可核实），但存在三类必须解决的问题：① 缺少生物学验证（仅有工程性能基准，无灵敏度/特异性数据）；② "完全离线"与"公共数据库检索+AI接口需要联网"存在内部矛盾，必须澄清；③ 参考库（8,464条）的构建标准、ICTV版本、软件版本号与代码可用性声明缺失。编辑决定：**Major Revision（大修）**。修订路线图含 9 项 Major、7 项 Minor，已按优先级排序，见文末。

---

## Phase 0 — Field Analysis & Reviewer Panel Configuration

| Item | Assessment |
|---|---|
| Primary discipline | Bioinformatics (viroinformatics / software platform paper) |
| Secondary disciplines | Plant virology, plant quarantine diagnostics, phylodynamics, software engineering |
| Paper type | Tool/platform description (implementation paper) |
| Target journal tier | Mid-tier SCI: *Virology Journal*, *BMC Bioinformatics*, *Frontiers in Plant Science*, *Phytopathology Research*, *Plant Methods* |
| Maturity | Implementation complete and internally verified; biological evaluation absent |

**Panel configuration:**
- **EIC** — Editor, virology-methods journal; cares about fit, claims discipline, availability compliance.
- **R1 (Methodology)** — Bioinformatician/statistician: identification model validity, benchmark design, reproducibility.
- **R2 (Domain)** — Plant virologist: positioning vs. existing plant-virus pipelines, reference database curation, ICTV/VMR usage.
- **R3 (Perspective)** — Research software engineer / One-Health infra: deployment reality, security, sustainability, FAIR compliance.
- **DA (Devil's Advocate)** — Challenges the core contribution claim.

---

## Phase 1 — Independent Review Reports

### Report 1 — Editor-in-Chief

**Recommendation: Major Revision** | Overall: 6/10

**Strengths.** The manuscript addresses a real, well-motivated gap (air-gapped Windows plant-quarantine laboratories) and reports an unusually complete implementation: 15 pipeline stages, 25 tools, a dual-engine identification model with an explicit statistical safeguard, wrapped RDP5 recombination analysis, dependency-free phylodynamics, and a documented packaging pipeline. The negative result on the in-house recombination prototype (18 detections vs. 1 true event on synthetic chimeras) is commendable transparency that raises trust. Claims are, with one exception, disciplined and tied to verifiable implementation facts.

**Major concerns.**
1. **The "entirely offline" claim conflicts with the paper's own Section 2.9** (NCBI SRA/GSA retrieval, Europe PMC, and OpenAI-compatible AI endpoints). This is not pedantry: it is the platform's headline differentiator and it must be stated precisely (analysis offline; metadata retrieval is an optional networked module; AI cleaning optional and disabled by default?). Location: Abstract, Introduction contribution 5, Section 2.9.
2. **No biological evaluation.** For a diagnostic-adjacent platform, the absence of any sensitivity/specificity or concordance data is the single largest obstacle to acceptance. If full validation is out of scope, the paper must (a) reframe claims accordingly, and (b) present at least one worked end-to-end case on a public dataset with known content.
3. **Availability is incomplete** (code "to be completed", no version number anywhere, no license). Most journals will desk-reject a software paper without a repository DOI. Assign a version, deposit code, and cite it.

**Minor.** Figures are referenced (Figure 1, Figure 2) but not provided; references [18] and [40] appear in the reference list but are never cited; the abstract should state the platform's software version.

---

### Report 2 — Reviewer 1 (Methodology & Statistics)

**Recommendation: Major Revision** | Rigor: 5/10, Reproducibility: 6/10, Clarity: 7/10

**Strengths.** The two-track identification model is clearly specified (base gate; Track A coverage/Poisson/depth; Track B gene-region rescue), thresholds are enumerated, and outputs include machine-readable discard reasons — good auditability. The ANI formula (aln_len − NM)/aln_len is explicit, and coverage-backend equivalence testing (bit-identical outputs) is a meaningful dependency-hygiene check.

**Major issues.**

1. **The Poisson safeguard is asserted, not validated (Section 2.3).** The model λ = N·L_read/L_ref assumes uniform random read placement, but real viral genomes exhibit strong coverage non-uniformity (terminal under-coverage, secondary structure, degradation); consequently R = (Coverage/100)/(1−e^(−λ)) is not distribution-free. The paper must show the empirical distribution of R for accepted vs. rejected references on real or spiked samples, and justify (or at least sensitivity-analyze) the 0.3 cut-off. If threshold sensitivity is unknown, state it as a limitation and note that the threshold is user-configurable.
2. **The 95% ANI species/strain boundary is family-dependent.** ICTV species demarcation thresholds vary widely across virus families (e.g., potyviruses ~76–82% for species, geminivirid thresholds differ again). A single fixed 95% cut-off needs justification (what is it calibrated against — strain vs. isolate?) and the "suspected novel" label should be softened or made family-aware, or explicitly documented as a pragmatic screening heuristic.
3. **Benchmark reporting is insufficient.** The 5.1 s / 7.9 s figures (Section 2.3, 3.3) come from one sample (GQMIX.q06) on unspecified hardware with no replicate variance. Report: CPU/RAM/storage, read count and length, number of references actually receiving reads, and ≥3 replicates. Even better, add a small panel of public SRA samples.
4. **Segment-completeness enforcement could mask biology.** Requiring all reference segments to be present will reject genuinely incomplete (partially degraded) multipartite infections. State whether partial detections are reported anywhere (they appear to be dropped from "confirmed" but should surface as a distinct category or flag).
5. **RDP5 event parsing: significance union.** Defining event significance as "methods with p < 0.05" and best-p as their minimum is RDP5's own display convention; fine — but the paper should note that RDP5 itself recommends confirmation by phylogenetic evidence, and that the exported masked alignment supports exactly that. One sentence will do.

**Minor issues.** (a) Track B thresholds (80%/5%) lack provenance — empirical? inherited? (b) "bit-identical on 3 samples" — give sample IDs/read counts in a supplementary table. (c) Alignment QC's identity criterion is correctly gated on explicit reference designation — good; consider reporting the number of removed/kept sequences as standard output metrics. (d) Reproducibility: requirements.lock is mentioned nowhere in the manuscript — cite it under reproducibility.

---

### Report 3 — Reviewer 2 (Plant Virology Domain)

**Recommendation: Major Revision** | Literature: 5/10, Domain contribution: 7/10

**Strengths.** The manuscript is grounded in plant-virology practice (quarantine, certification, viroids, multipartite viruses, vector transmission views in the browser), and the CMV RNA3 example used to justify gating the alignment-QC identity criterion shows genuine hands-on familiarity. The segment-designation normalization across 299 public spellings addresses a real, chronically annoying problem in plant-virus metadata.

**Major issues.**

1. **Related work is too thin on plant-virus-specific tools.** Beyond VirusDetect/VirusSeeker/BioAider, the field has metagenomic classifiers and virus-identification tools that readers will expect to see discussed as alternatives: Kraken 2 (used by many plant-virome studies), Kaiju, geNomad/VirSorter2 (virus identification from contigs), and end-to-end virome pipelines applied to plants. A short paragraph contrasting design choices (k-mer classification vs. alignment vs. translated search) and explaining why the platform pivoted from Kraken2/kunpeng-based screening to reference-alignment quantification (the pivot is mentioned in passing at Section 2.2) would materially strengthen the paper.
2. **Reference database provenance is underspecified.** 8,464 sequences: selected how, from what source snapshot, filtered by which completeness/length criteria, deduplicated how, and against which ICTV/VMR release? This is the platform's core asset; it needs a construction paragraph (or a supplementary table) and a stated update policy.
3. **No figures.** A platform paper without its architecture and workflow figures is hard to evaluate; Figure 1 (pipeline) and Figure 2 (known-virus suite with Poisson model) are referenced but absent.
4. **Viroid handling deserves a sentence of method detail.** Viroids are stated to be handled (blastn route; SnpEff skipped for CDS-less genomes) — a reader in quarantine will want to know the detection path explicitly, since viroids are PRA-relevant organisms.
5. **Primer design module** is listed (conserved/full-length modes, thermodynamic scoring) but never described; either describe it in Methods (2–3 sentences) or remove it from the abstract-level claims.

**Minor.** (a) State ICTV Master Species List / VMR release dates used. (b) "8-level lineage classification" of contigs — say what the levels are. (c) The claim that platform phylogenetics uses "dual branch support" (UFBoot + SH-aLRT) is current best practice — good; cite SH-aLRT (Guindon et al. 2010) or the IQ-TREE 2 paper where the combination is benchmarked (Minh et al. 2020). (d) Consider citing the ICTV Report (current) alongside Lefkowitz 2018.

---

### Report 4 — Reviewer 3 (Research Software Engineering / Deployment Perspective)

**Recommendation: Major Revision** | Practical impact: 8/10, Sustainability: 5/10

**Strengths.** The deployment story is unusually credible: loopback-only binding with Host/Origin/Referer validation, path-traversal guards, shell-free subprocess execution with timeout watchdogs and cancellation, two-tier semaphore task manager with SSE progress and restart recovery, three-way distribution split with post-build self-check, and a route-inventory regression guard. The bilingual interface and bundled examples directly serve the target user group. This is more operational rigor than most academic platform papers demonstrate.

**Major issues.**

1. **Versioning and release discipline are absent from the manuscript.** No version number, no changelog reference, no release artifact identifiers. For software with regulatory-adjacent uses, versions of the *bundled tools* should also be enumerated (the platform pins IQ-TREE 3.0.1, bcftools 1.24, kunpeng 0.7.12 — say so in a supplementary table).
2. **The AI metadata layer needs a governance paragraph.** Which model (DeepSeek/Kimi) is documented, but not: default on/off, data sent off-premises (a hard problem for the stated air-gapped audience!), failure behavior when no API key is configured, and audit trails of AI-touched fields. The provenance suffix is good; the default-off and network-boundary statements are missing.
3. **Sustainability.** A 90,066-line core with 30 bundled third-party binaries is a maintenance liability; reviewers at software journals will ask about update policy (tool upgrades, database refresh cadence, OS-target policy, and what happens when a bundled binary is CVE-patched upstream). One paragraph in Discussion or Availability suffices.

**Minor.** (a) Engineering density in the main text (route counts, semaphore values) may exceed the interest of the virology readership — consider moving detail to "Software metadata" supplementary. (b) State explicitly that the application must never be exposed beyond loopback. (c) The six-hour command timeout default and its configurability are worth documenting for users processing large runs. (d) Cite the front-end test automation (85 scripts) as evidence but clarify which portions run in CI vs. manually.

---

### Report 5 — Devil's Advocate

**Strongest counter-argument.** *"This is a wrapper, not a contribution."* Every individual component — fastp, salmon/minibwa, SPAdes, DIAMOND, MMseqs2, MAFFT, IQ-TREE, RDP5, Fitch parsimony — is a published third-party tool; the manuscript's own negative result (18:1 false detections from the in-house recombination prototype) shows the integrator does not outperform the tools it wraps. On this view, VirusPlatform is competent systems integration with a GUI, and the scientific content reduces to one small statistical heuristic (the Poisson ratio) whose validity the paper never demonstrates. Under this reading, the paper belongs in a software venue, not a science venue, and its "lowers the barrier" claim is untested — there is no user study, no deployment record, no diagnostic concordance data.

**Issue list.**

| # | Severity | Dimension | Location | Issue |
|---|---|---|---|---|
| DA-1 | MAJOR | Validity of headline claim | Abstract, §5 | "Entirely offline" is contradicted by the online metadata/AI module — an internal-consistency defect in the paper's core positioning |
| DA-2 | MAJOR | Contribution substantiation | §2.3, §3 | The one candidate for genuine scientific novelty (Poisson safeguard) is never validated; if it is not validated, the paper must honestly reframe as engineering integration |
| DA-3 | MAJOR | Unfalsifiable benefit claim | §5 | "Lowers the barrier" is asserted without any user-experience or adoption evidence |
| DA-4 | MINOR | Selection effects | §3.3 | All quantitative evidence is self-reported from development benchmarks; no external benchmark or third-party reproduction |
| DA-5 | MINOR | Rhetorical risk | §4 | "To our knowledge, unusual" dual-engine claim — soften or support with a survey of related tools |

**Ignored alternatives.** (1) A containerized (Docker/Singularity) Linux distribution — cheaper to build, standard in bioinformatics; the paper should argue why Windows-native beats containers *for this user base* (it plausibly does: no WSL on locked-down government machines — but say it). (2) Contributing the statistical safeguard upstream into existing pipelines as a small library, rather than a full platform.

**Missing stakeholders.** Diagnostics regulators (validation dossiers), non-Windows users, educators (the example sets have pedagogical value — say so).

**"So what?" test.** Passes for the target community: quarantine labs on air-gapped Windows machines genuinely cannot use current cloud/Linux pipelines. The need is real; the solution is real; the *evidence* that it works as a diagnostic aid is what is missing.

**Observations (non-defects).** The transparent reporting of the failed in-house recombination detector is the single most trust-building element in the manuscript; keep it prominent. The explicit "failure ≠ negative" batch semantics is the kind of operational honesty regulators notice.

---

## Phase 2 — Editorial Decision & Revision Roadmap

### Editorial Decision: **MAJOR REVISION**

**Rationale.** Five reviews converge: the implementation is credible and unusually well engineered; the manuscript is honest about its negative results; but (i) the core evaluation evidence (biological validation) is absent, (ii) the offline/online positioning is internally inconsistent, and (iii) software-publishing compliance items (version, repository, license, figures) are incomplete. No CRITICAL fabrication-level defect was found, so major revision (rather than reject/restart) is appropriate.

**Consensus issues (multiple reviewers):** offline/online contradiction (EIC, R3, DA-1); missing validation (EIC, R1, DA-2); reference-database provenance (R1, R2); availability/versioning (EIC, R3); missing figures (EIC, R2).

### Revision Roadmap (prioritized)

| ID | Priority | Requirement | Addressed by |
|---|---|---|---|
| R-1 | Major | Resolve the offline/online contradiction: reposition as "offline-capable analysis; optional networked metadata/AI module (off by default, clearly network-boundaried)" across Abstract/Intro/§2.9/§5 | Revision |
| R-2 | Major | Add a validation subsection: worked end-to-end demonstration on public CMV/PVY data with known composition + report identification results; explicitly scope claims to "implementation evaluation" where validation is absent | Revision |
| R-3 | Major | Reference-database provenance: construction criteria, source snapshot, dedup, ICTV/VMR release versions, update policy (new §2.3 subsection) | Revision |
| R-4 | Major | Benchmark reporting: hardware, sample composition, replicates; qualify the single-sample caveat; move details to a table | Revision |
| R-5 | Major | Poisson model: state assumptions, sensitivity caveat, configurability; present R-distribution evidence if available or declare as future validation; ANI 95% justification + family-dependence caveat | Revision |
| R-6 | Major | Software compliance: assign version (e.g., v1.0), complete Availability (repository/license), enumerate pinned tool versions in supplementary table | Revision |
| R-7 | Major | Expand related work: Kraken 2/Kaiju/geNomad/VirSorter2 contrast; explain the Kraken2→alignment pivot; containers-vs-Windows-native justification | Revision |
| R-8 | Minor | Add Figures 1–3 (architecture; known-virus suite; evolutionary-dynamics module) with full captions | Revision |
| R-9 | Minor | Viroid detection path detail; primer-design module description (or drop from claims); ICTV MSL release dates; cite SH-aLRT/IQ-TREE 2; cite requirements.lock; soften "to our knowledge, unusual"; segment-incompleteness surfacing; move engineering minutiae to supplementary; sustainability paragraph; AI governance paragraph (default-off, data boundary); fix uncited refs [18]/[40]; mention user-study/adoption as future work | Revision |

### Scores summary

| Reviewer | Recommendation | Key score |
|---|---|---|
| EIC | Major Revision | Overall 6/10 |
| R1 Methodology | Major Revision | Rigor 5/10 |
| R2 Domain | Major Revision | Literature 5/10 |
| R3 RSE | Major Revision | Impact 8/10 |
| Devil's Advocate | — (5 issues: 4 MAJOR, 1 MINOR) | No CRITICAL |

*Review performed per academic-paper-reviewer skill v1.10.0 (full mode). Reports are read-only with respect to the manuscript; revision is executed separately in the pipeline's revision stage.*
