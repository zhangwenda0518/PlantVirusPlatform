<div align="center">

[English](README_EN.md) · [简体中文](README.md)

# 🌿 VirusPlatform — Plant Virus Analysis Platform

**A Windows-local plant-virus diagnostics & deep-analysis platform built on [kunpeng](https://github.com/afshinokh/kunpeng) (an ultra-low-memory metagenomic classifier)**

From raw FASTQ to NCBI submission files, **entirely point-and-click**:
Host removal → known-virus identification & quantification → assembly → ORF prediction → functional annotation → phylogeny / SDT / molecular clock / phylogeography → primer design → genome plots → interactive reports → public-data tracing → NCBI submission prep

![Platform](https://img.shields.io/badge/platform-Windows%2010%2B-0078D6?logo=windows11&logoColor=white)
![Python](https://img.shields.io/badge/python-3.12%2B-3776AB?logo=python&logoColor=white)
![GUI](https://img.shields.io/badge/GUI-Flask%20%2B%20pywebview-8B5CF6)
![CLI](https://img.shields.io/badge/CLI-%E2%88%9A%20main.py-16A34A)
![i18n](https://img.shields.io/badge/UI-Bilingual%20CN%2FEN-EAB308)
![Docs](https://img.shields.io/badge/docs-Wiki%2025%20pages-blue)

</div>

---

## ✨ At a Glance

- **One pipeline, end to end** — a 14-step cascaded analysis (⓪–⑩) where each step's output feeds the next, with resume-from-checkpoint, forced reruns, and automatic graceful degradation of optional steps.
- **Local-first, offline-friendly** — all external tools (kunpeng / SPAdes / MAFFT / RAxML-NG / DIAMOND / salmon / primer3 …) ship with the platform or are auto-detected; databases are built and searched locally; downloads hit only an NCBI/ENA/NGDC official-domain whitelist.
- **Dual entry: GUI + CLI** — a pywebview desktop window or browser UI (bilingual CN/EN one-click switch), backed by the same engine that powers the full `main.py` CLI for scripted batch runs.
- **Phylogenetics beyond tree building** — built-in TreeTime ML molecular clock, discrete-trait geographic reconstruction (mugration), skyline plots, RTT/DRT regression, BSP, BEAST hand-off with Markov-jump parsing, and phylogeographic migration animations (GIF) — plant-virus epidemiology in one place.
- **It gets smarter with use** — per-stage runtimes feed a self-learning model that drives real-time progress and ETA on every task card; plus concurrency gates for heavy jobs, a disk-watermark banner, and resource-estimate logging.
- **A closed research loop** — public-data search (NCBI SRA + NGDC GSA) → bulk download → analysis → LOGAN tracing across all of SRA → NCBI GenBank/BioSample submission prep, without ever leaving the platform.

## 🧭 Analysis Workflow

```mermaid
flowchart LR
    A["⓪ QC<br/>fastp · fq2fa"] --> B["① Host removal<br/>kunpeng host DB"]
    B --> C["②b ID & quant<br/>kvsuite · salmon"]
    C --> D["③ Assembly<br/>SPAdes rnaviral"]
    D --> D2["③b/③c Verify<br/>consensus & variants"]
    D2 --> E["④ Host prediction<br/>ICTV cascade"]
    E --> F["⑥ ORF calls<br/>pyrodigal"]
    F --> G["⑥b Annotation<br/>DIAMOND · Pfam · CDD"]
    G --> H["⑦ Tree/SDT<br/>MAFFT · FastTree/RAxML"]
    H --> I["⑧ Primers<br/>primer3"]
    I --> J["⑨ Genome plots<br/>gbdraw / DFV"]
    J --> K["⑩ Report<br/>plotly · pycirclize"]
    K --> L["Tracing / clock /<br/>phylogeo / submission"]
```

## 📦 Module Overview

> Each module's **detailed documentation (inputs/outputs · parameters · CLI · product folders · caveats) lives in the [Wiki](https://github.com/zhangwenda0518/VirusPlatform/wiki)** (pages in Chinese).

### Main analysis pipeline (per-sample, ⓪ → ⑩)

| Step | Wiki | Engine | Output |
|---|---|---|---|
| ⓪ fastp QC (optional) | [Preprocess](https://github.com/zhangwenda0518/VirusPlatform/wiki/Preprocess) | fastp | clean reads, QC report |
| ⓪b FASTQ→FASTA (optional) | same page | seqkit fq2fa | conv_R1/R2.fa.gz (speeds up classification) |
| ① Host removal | same page | kunpeng (host DB) | kept_R1/R2, host fraction |
| ②b Known-virus ID & quant | [Kvsuite](https://github.com/zhangwenda0518/VirusPlatform/wiki/Kvsuite) | salmon (EM quant) / minibwa | ID table, filtered quant table, BAM+depth, viral reads |
| ③ Assembly & classification | [Assembly](https://github.com/zhangwenda0518/VirusPlatform/wiki/Assembly) | SPAdes metaviral/rna + BLAST | contigs, viral-contig annotation table |
| ③b Candidate verification | same page | host filter → length split → blastx/CDD/viroid | verified candidate set |
| ③c Consensus & variants | same page | reads re-mapped to reference → per-position counts | consensus sequences, variant table |
| ④ Host prediction | [Host-Prediction](https://github.com/zhangwenda0518/VirusPlatform/wiki/Host-Prediction) | kunpeng + ICTV host-probability cascade | host_prediction.tsv, sankey/sunburst |
| ⑥ ORF calling | [ORF](https://github.com/zhangwenda0518/VirusPlatform/wiki/ORF) | pyrodigal + pyrodigal-rv | faa/ffn/gff |
| ⑥b ORF annotation | same page | DIAMOND/MMseqs2/blastp × RefSeq viral proteins + Pfam HMM + CDD | orf_annotation.tsv, GFF3, feature diagrams |
| ⑦ Phylogeny & SDT | [Phylo](https://github.com/zhangwenda0518/VirusPlatform/wiki/Phylo) | MAFFT + trimAl + FastTree/RAxML-NG/NJ; exact SDT engine | tree.nwk, identity matrix + heatmap |
| ⑧ Primer design | [Primer](https://github.com/zhangwenda0518/VirusPlatform/wiki/Primer) | primer3 (conserved/full-length) + thermodynamic scoring | primers.tsv |
| ⑨ Genome plots | [Genome-Plots](https://github.com/zhangwenda0518/VirusPlatform/wiki/Genome-Plots) | gbdraw / dna_features_viewer | circular + linear SVG |
| ⑩ Visual report | [Report](https://github.com/zhangwenda0518/VirusPlatform/wiki/Report) | plotly + pycirclize + matplotlib | report.html |

### Extended analyses

| Module | Wiki | Description |
|---|---|---|
| Phylodynamics toolkit | [Phylodynamics](https://github.com/zhangwenda0518/VirusPlatform/wiki/Phylodynamics) | TreeTime ML clock / mugration / skyline / ancestral sequences, RTT·DRT regression, BSP, MCC-tree parsing, spatiotemporal downsampling, metadata governance, alignment QC |
| BEAST interface | [BEAST](https://github.com/zhangwenda0518/VirusPlatform/wiki/BEAST) | server job package generation/read-back verification, Markov-jump (MJRM) XML injection & log parsing, MCC state parsing |
| Phylogeography | [Phylogeography](https://github.com/zhangwenda0518/VirusPlatform/wiki/Phylogeography) | metadata-driven migration reconstruction + migration-arc GIF animation; RDP5 recombination; pre-recombination sequence pre-filtering |
| Public data search & download | [Public-Data](https://github.com/zhangwenda0518/VirusPlatform/wiki/Public-Data) | SRA+GSA dual-engine search, unified metadata (Core14/Full), publication-grade metadata figures, aria2c bulk download + built-in sracha for GSA .sra, host-genome download |
| Toolbox | [Toolbox](https://github.com/zhangwenda0518/VirusPlatform/wiki/Toolbox) | end-to-end dsRNA design (windowing → efficacy → off-target → primers), miRNA target prediction (5 engines), NCBI reference download, GenBank collections & CDS/PEP export, MSA viewer, contig depth annotation, one-step tools |
| NCBI submission prep | [Submission](https://github.com/zhangwenda0518/VirusPlatform/wiki/NCBI-Submission) | one unified_metadata.csv table drives it all: online editing → validation → source.src / miuvig / assembly / BioSample / template.sbt |
| LOGAN tracing | [Logan](https://github.com/zhangwenda0518/VirusPlatform/wiki/Logan) | submit viral contigs to Logan-Search (~23.4M public SRA samples), batch submission + SSE live progress + tracing report |

### Platform infrastructure

| Module | Wiki | Description |
|---|---|---|
| Databases & references | [Databases](https://github.com/zhangwenda0518/VirusPlatform/wiki/Databases) | Taxonomy, kunpeng DBs, virus reference/ID DB, annotation DBs (RefSeq proteins / Pfam / CDD), tree DBs (plant_tree.db / ICTV VMR), host-probability tables, host-DB building |
| Pipeline orchestration | [Pipeline-Overview](https://github.com/zhangwenda0518/VirusPlatform/wiki/Pipeline-Overview) | dependency scheduling, checkpoint resume, weighted progress with self-learning ETA, resource estimates, concurrency gates |
| Settings · tasks · ops | [Settings-Ops](https://github.com/zhangwenda0518/VirusPlatform/wiki/Settings-Ops) | bilingual UI, 17 default parameters, disk watermark, environment self-check, tool_runs archiving/cleanup, database migration, packaging |
| Full CLI reference | [CLI-Reference](https://github.com/zhangwenda0518/VirusPlatform/wiki/CLI-Reference) | all 30+ subcommands with per-parameter quick reference |
| FAQ & known issues | [FAQ](https://github.com/zhangwenda0518/VirusPlatform/wiki/FAQ) | common questions, non-ASCII path handling, large-genome DB building, SDT semantics, ML tree building… |

## 🚀 Quick Start

### 1. Install

```bat
git clone https://github.com/zhangwenda0518/VirusPlatform.git
cd VirusPlatform
python -m pip install -r requirements.txt
```

- Requires **Windows 10+ and Python 3.12+**. **All external tools ship with the repository** (`3rd/`, ~2GB: single-file tools + tool suites + a bundled Python runtime) — clone and go.
- **SPAdes** must be installed separately (`SPAdes-Windows-4.3.0-dev-Setup.exe`; put `spades.bat` on PATH).
- gbdraw is the preferred genome-plot engine: `pip install git+https://github.com/satoshikawato/gbdraw.git` (falls back automatically to the pure-Python dna_features_viewer).
- Only two files exceed GitHub's 100MB per-file hard limit and are served from the [Release (external-tools)](https://github.com/zhangwenda0518/VirusPlatform/releases/tag/external-tools) instead — download and put them back in place:
  `3rd/tools/rdp5/3seqTable` (105MB) and `3rd/python/Lib/site-packages/_polars_runtime_32/_polars_runtime.pyd` (168MB).
- Don't want the ~2GB of prebuilt tools? Build only the subset you need — see **Building the External Tools from Source (Optional)** below.
- Databases (~3.4GB) are not in git: copy them from a database package / a packaging machine, or build with `dev_tools/package.py` (see [Packaging](#-packaging--distribution-software--examples--databases-in-three-separate-folders)).

### 2. Launch (two modes)

| Mode | Command | Notes |
|---|---|---|
| Desktop window | `python app.py --gui` | pywebview/WebView2 shell around the local UI; no address bar; auto-falls back to browser |
| Web mode | `python app.py --web` | opens in the default browser; the URL is printed to the console |
| Auto (default) | `python app.py` | prefers the desktop window, falls back to browser |

Port probing starts at **8765** and increments automatically — **use the URL printed in the console**. Binds to `127.0.0.1` only; never exposed to the network.

### 3. First-time setup (once)

1. **Databases page** → download NCBI Taxonomy (~57MB);
2. **Databases page** → pick a host genome FASTA + the host TaxID (e.g. 4081 for tomato) → build the host DB (the virus DB ships prebuilt — nothing to do);
3. **Settings page** → threads, e-mail, and 17 analysis defaults as needed (defaults work out of the box).

### 4. Run your first sample

**Pipeline page** → new sample (name + R1/R2) → ▶ Run remaining steps.
For a first pass, set subsampling to 100,000 read pairs to validate the whole flow in minutes.

CLI equivalent:

```bat
python main.py selfcheck                                   :: environment self-check
python main.py init-taxonomy                               :: download taxonomy
python main.py build-host-db --genome host-db\genome.fa --taxid 4081
python main.py analyze --r1 R1.fastq.gz --r2 R2.fastq.gz --sample NX-5
python main.py report --sample NX-5                        :: regenerate the report
```

## 🔧 Building the External Tools from Source (Optional)

**First, the common question: compiling from source does NOT reduce runtime memory.** Runtime RAM is determined by algorithms and data scale (kunpeng's hash capacity, SPAdes' assembly graph, DIAMOND's search index), not by where the binary came from — a locally built binary behaves identically to the official release of the same source.
What source builds actually buy you is **clone size** (2GB of prebuilt binaries → tens of MB of sources) and a **need-based subset** (build only the chain you use). The price is paid at build time: linking big C++ projects (SPAdes/BLAST) can peak at several GB of RAM, with another 5–10GB of disk for intermediate artifacts.
Also, `3rd/python/` (the 868MB bundled interpreter) is optional — a system Python 3.12 plus `python -m pip install -r requirements.txt` works just as well (the bundled interpreter only serves the packaged exe).

### Toolchain setup

```bat
winget install MSYS2.MSYS2 Rustlang.Rustup GoLang.Go Kit.CMake
:: inside the MSYS2/UCRT64 shell:
pacman -S mingw-w64-ucrt-x86_64-gcc make mafft samtools bcftools
:: on networks where proxy.golang.org is unreachable:
set GOPROXY=https://goproxy.cn,direct
```

### Tiered build guide

| Tier | Tools (put built binaries in) | Language | Notes |
|---|---|---|---|
| ① minutes | kunpeng / crabz / sracha → `3rd/bin/` | Rust | `cargo build --release` |
| ① | seqkit → `3rd/bin/` | Go | `go install github.com/shenwei356/seqkit/v2/seqkit@latest` (main package lives in the module's `seqkit/` subfolder) |
| ① | FastTree / trimAl / minimap2 / minibwa | C | straight `gcc`/`make`; the minibwa Windows port lives in `docs/minibwa-build/` (mmap-stub patch + parity report) |
| ② ready-made MSYS2 packages | mafft / samtools / bcftools | C | one `pacman -S` each |
| ② CMake | diamond / iQ-TREE / RAxML-NG / lsd2 / SPAdes | C++ | CMake + VS Build Tools; for SPAdes the official installer is easier |
| ③ not recommended locally | Blast+ / salmon / mmseqs | C++ | extremely heavy toolchains or degraded Windows support (salmon pushes WSL/Docker; mmseqs native build is experimental) — use the prebuilt `3rd/tools/` copies |
| ④ binary-only | RDP5 (VB6, no public source) / table2asn / Gblocks / strawberry-perl / dsrnamax / pandepth / snpgenie / viral_consensus | — | no upstream sources (or platform-internal); rely on the prebuilt `3rd/` copies |
| side | SnpEff | Java | source is public, jar can be self-built; the slim JRE in `3rd/tools/jre-snpeff/` only serves runtime |
| side | open-virome frontend | Node | source ships in-repo: `npm ci && npm run build` (rename `.git-upstream/` back to `.git` for upstream git work) |

### ✅ Build test log (2026-09-28, on this machine: Win11 · MinGW gcc 16.2 / cargo 1.98 / go 1.27 / cmake 4.4)

| Tool | Result | Notes |
|---|---|---|
| crabz | ✅ v0.10.1 | straight from crates.io, ~2 min |
| sracha | ✅ v0.7.0 | `rnabioco/sracha-rs`, ~3 min |
| seqkit | ✅ v2.14.0 | `go install .../v2/seqkit@latest`; use goproxy.cn where proxy.golang.org is unreachable |
| FastTree | ✅ v2.2.0 | author mirror `morgannprice/fasttree`, single-file gcc build |
| trimAl | ✅ | `inab/trimal`, `make -C source` |
| minimap2 | ✅ v2.31 | `lh3/minimap2`, `make` |
| minibwa | ✅ 0.7-r424-dirty | apply `docs/minibwa-build/windows-mmap-stub.patch` then **`make mimalloc=0`** (the vendored mimalloc is trimmed and has no Windows backend; disabling falls back to kalloc); index+map smoke test passed |
| kunpeng | ⚠️ blocked | crates.io only has 0.6.9–0.7.5; the bundled v0.7.12 carries **locally modified sources** — rebuilding requires that modified source tree, the vanilla crate is no substitute |
| RAxML-NG | ❌ MinGW route blocked | cleared 4 hurdles (coraxlib's `__declspec(dllexport)` on `__thread`, `M_PI`, `asprintf`, size_t width) then hit `sysutil.cpp`'s POSIX header `sys/resource.h` — a real port, not a patch. **The realistic tier-② route = official MSVC releases/installers** (SPAdes installer, RAxML-NG official release), exactly as the table advises |
| MSYS2 pacman route | ❌ not verified | the silent installer popped a GUI here and never finished; documented route kept |
| RDP5.93 (tier ④ install test) | ❌ not adopted | its RDP5CL is actually **5.84**: -f/-ofp honored, all nine methods run, exit 0 — but **zero output files written**; bundled **5.69** stays. Both 3seqTable copies are byte-identical (110,174,436 B). Install kept at `3rd/tools/rdp5.93/` (gitignored) |

- **Paths**: drop binaries into the matching `3rd/` location and they are auto-detected; anywhere else, register them in the `tools` section of `platform.json` (takes effect immediately).
- **Verify**: `python main.py tools` (detection status) → `python main.py selfcheck` (full check-up). A missing tool only affects its own step — optional steps (fastp / fq2fa / gbdraw) gray out and get skipped automatically; required steps clearly name what's missing.

## 🖥️ UI Tour

The top navigation follows the workflow (🌐 switches CN/EN; the task-center entry carries a running-tasks badge):

**Overview** (cover + DB status + workflow entries) · **Pipeline** (14 cascaded step cards) · **Results** (reports / extended results / tree·MSA·SDT viewers) · **Toolbox** (one-step tools + reference data) · **Downloads** (SRR/CRR/URL bulk) · **Public Data Search** · **Databases** (Taxonomy / host DB) · **Public Virome** · **LOGAN Tracing** · **Submission Prep** · **Task Center** · **Settings** · **Manual**

Highlights:
- Every finished task auto-expands a **result-preview panel** (key numbers + product list + open report); system notifications when the page is in the background;
- 9 pipeline step cards carry a `📊 view` button — tables open in-place, HTML reports in a new window;
- Storage page: free-space bar + per-subdatabase usage breakdown; a site-wide yellow banner appears below 20GB free;
- Concurrency gates: ≤2 heavy / ≤4 light jobs, with visible queue positions.

## 📁 Directory Layout

```
VirusPlatform/
├─ app.py                    Web GUI service (localhost only)
├─ main.py                   CLI entry (30+ subcommands)
├─ Virus_Platform_Core/      Core pipeline package (shared by GUI & CLI, 70+ modules)
│  ├─ web/                   Flask blueprints (pages/API/tasks)
│  ├─ known_virus_suite/     ②b known-virus ID & quant engine
│  ├─ ncbi_submit/           NCBI submission-prep engine
│  ├─ public_meta/           public-data search engine
│  └─ mirna_target/          miRNA target prediction (5 engines)
├─ webapp/                   page templates & static assets (i18n bilingual)
├─ tests/                    self-checks, data generators, integration tests
├─ dev_tools/                packaging & DB-build scripts (package.py etc.)
├─ docs/                     development docs & module inventories
├─ wiki/                     source Markdown of this repo's Wiki (kept in sync)
├─ examples/                 built-in example data (✨Example buttons, ~8MB)
├─ 3rd/                      external dependencies (bin/ tools/ python/ open-virome/,
│                            all tracked in git; only the two >100MB files come from
│                            the Release)
├─ databases/                virus classification/reference/annotation/tree DBs (~3.4GB, not in git)
├─ host-db/                  host genomes & host classification DBs (one DB per species)
├─ platform.json             tool paths / language / analysis defaults
└─ run/                      runtime data (results/logs/tasks/…, rebuildable)
```

Legacy path compatibility: platform-relative paths like `tool_runs/`, `results/` are auto-redirected into `run/`; `tools/`, `bin/` → `3rd/`. No old calls need rewriting.

## 📤 Packaging & Distribution (software / examples / databases in three separate folders)

```bat
python dev_tools/package.py            :: ① program + ② examples (~1.5GB)
python dev_tools/package.py --with-db  :: plus ③ database pack (~3.4GB)
python dev_tools/package.py --db-only  :: database pack only
python dev_tools/package.py --verify   :: run exe --cli selfcheck after packing
```

The program pack (`VirusPlatform.exe`) upgrades/distributes independently; the examples and database packs are auto-detected when placed next to the program — zero configuration. Databases can also live on any drive: `python main.py db-migrate --to D:\dbdir` (robocopy + byte verification before switching config), or set the path in Settings → Database directory.

## 📚 Documentation

- **[Wiki home](https://github.com/zhangwenda0518/VirusPlatform/wiki)** — 25 pages, one per module (quick start / pipeline stages / phylodynamics / BEAST / phylogeography / public data / toolbox / submission / LOGAN / databases / CLI / FAQ). Pages are written in Chinese.
- Wiki sources are kept in the [`wiki/`](wiki/) folder of this repository.
- Development notes: [`docs/DEVELOPMENT_NOTES.md`](docs/DEVELOPMENT_NOTES.md) (Chinese).

## ❓ FAQ (quick picks)

<details>
<summary><b>Page shows "Failed to fetch" / clicks do nothing?</b></summary>

The server process is gone — closing the black console window quits the platform; or the browser tab belongs to a dead instance; or a heavy job is saturating memory/CPU. The UI shows a red diagnostic banner while disconnected and clears it automatically on recovery.
</details>

<details>
<summary><b>Which input formats are supported?</b></summary>

FASTQ / FASTQ.gz (paired: R1+R2; single-end: R1 only); FASTA / FASTA.gz; GenBank .gb/.gbk (genome plots & reference collections). Unpack .zip/.rar/.tar archives first.
</details>

<details>
<summary><b>SPAdes / BLAST errors (non-ASCII paths)?</b></summary>

Handled automatically: SPAdes runs through an ASCII-safe `%TEMP%\vp_spades` staging dir; BLAST databases are built under `%TEMP%\vp_blast`.
</details>

<details>
<summary><b>Host-DB building fails on large genomes?</b></summary>

Fixed: sequences are auto-split into ≤1MB chunks (34bp overlap preserves cross-cut k-mers). A 1.8GB genome builds with 2GB of RAM, with k-mer counts identical to whole-sequence indexing.
</details>

More in the [Wiki FAQ](https://github.com/zhangwenda0518/VirusPlatform/wiki/FAQ) (Chinese).

## 🙏 Credits & Acknowledgements

The platform integrates many excellent open-source tools and data resources, including (but not limited to):

- **kunpeng** — ultra-low-memory metagenomic classifier (host removal / virus classification core)
- **SPAdes** (metaviral/rna), **salmon**, **minibwa**, **fastp**, **seqkit**, **crabz**
- **pyrodigal / pyrodigal-rv / orfipy** (ORF calling), **DIAMOND / MMseqs2 / BLAST+** (homology search)
- **pyhmmer × Pfam**, NCBI **CDD** (domain annotation), NCBI **RefSeq Viral** (protein DB)
- **MAFFT / trimAl / Gblocks / FastTree / RAxML-NG / IQ-TREE 3** (phylogenetics), **TreeTime** (molecular clock)
- **BEAST 1.x** (Markov-jump counting semantics; no local MCMC), **Logan-Search** (Chikhi et al. 2025, bioRxiv 10.1101/2024.07.30.605881)
- **ICTV VMR / MSL** (taxonomy metadata), NCBI **Taxonomy / SRA / GenBank / BioSample**, NGDC **GSA**
- **primer3** (primers), **gbdraw / dna_features_viewer** (genome plots), **plotly / pycirclize / matplotlib** (visualization)
- **RDP5** (recombination detection), **sracha** (GSA .sra conversion, pure Rust)
- Methodological alignment: VirPhyKit (Yin et al., 2025, Ecol Evol; code fully re-implemented, only input/output contracts aligned), MMPV-RNA pipelines (host-prediction C9 cascade / SDT semantics / submission prep)
