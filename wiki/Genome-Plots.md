# Genome-Plots · ⑨ 基因组图

> 返回 [Home](Home) ｜ 上一页：[Primer](Primer) ｜ 下一页：[Report](Report)

对应引擎：`gbdraw_plot.py`（gbdraw）、`dfv_plot.py`（dna_features_viewer 回退）。

---

## 双引擎

| 引擎 | 定位 | 注释来源 |
|---|---|---|
| **gbdraw**（默认首选） | 圈图更精美（SVG 圈图 + 线图） | 自动取 ⑥ ORF 的 pyrodigal GFF3 |
| **dna_features_viewer（DFV）** | gbdraw 不可用时自动顶上（纯 Python，无外部 CLI） | 优先 ⑥b 的 `orf_annotation.gff3`（带 product/category，按结构蛋白/聚合酶等功能类别着色），回退 pyrodigal GFF3；支持 GenBank 输入 |

- 引擎选择：卡片参数或 CLI `--plot-engine auto|gbdraw|dfv`（默认 auto）；
- 两引擎都未安装时 ⑨ 才灰显（DFV 随 requirements.txt 安装，因此实际很少触发）。

## 自备文件出图

卡片参数（或 CLI `--gbdraw-fasta/--gbdraw-ann`）改用**自备 FASTA + GFF3，或直接 GenBank `.gb/.gbk` 文件**——可以对任意病毒基因组出图，不限于本样品。

- 出图上限默认 12 条（按长度取最长，`--gbdraw-max`）；
- SVG 内嵌 ⑩ 报告；结果文件在 `09_genome_plots/<contig>.circular.svg / .linear.svg`。

## Windows 注意：SVG MIME

Windows 注册表把 `.svg` 关联成非标准的 `image/svg`，Python mimetypes 照搬后 Flask 发出的 MIME 会让 Chromium 拒绝解码（`<img>` 裂图 naturalWidth=0）。`app.py` 建应用前已 `mimetypes.add_type('image/svg+xml', '.svg')`——该修复在 app.py 里，源码运行即时生效。

## CLI

```bat
python main.py analyze --r1 ... --sample S1 --stages gbdraw --plot-engine auto --gbdraw-max 12
python main.py analyze ... --gbdraw-fasta x.fa --gbdraw-ann x.gff   :: 自备文件
python main.py analyze ... --gbdraw-ann x.gb                        :: GenBank 直接出图
```

## 相关

- 注释从哪来 → [ORF](ORF)（⑥ GFF3 / ⑥b 带功能类别的 GFF3）
- 图如何汇总 → [Report](Report)
