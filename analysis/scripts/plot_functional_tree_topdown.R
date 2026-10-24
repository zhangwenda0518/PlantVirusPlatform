# =====================================================================
# VirusPlatform 功能树（SCI 发表级）—— 顶向下布局变体
# 树朝下生长：根在顶部，50 项末端功能横排在下部（纵排标签），
# 平台真实输出示例图排在最底部、与各功能域的水平跨度对齐
# 用法: Rscript plot_functional_tree_topdown.R EN|CN
# 输出: figures/functional_tree_topdown_EN|CN .pdf/.png
# =====================================================================

suppressPackageStartupMessages({
  library(ggplot2); library(showtext); library(png)
})

lang <- toupper(commandArgs(trailingOnly = TRUE)[1])
DENS <- toupper(commandArgs(trailingOnly = TRUE)[2])
if (!lang %in% c("EN", "CN")) lang <- "EN"
# 第二参数：FULL 全面 / SIMPLE 简单
if (!DENS %in% c("FULL", "SIMPLE")) DENS <- "FULL"

# ---- 字体 ------------------------------------------------------------
if (lang == "EN") {
  font_add("mainfont",
           regular = "C:/Windows/Fonts/arial.ttf",
           bold    = "C:/Windows/Fonts/arialbd.ttf",
           italic  = "C:/Windows/Fonts/ariali.ttf")
  fam <- "mainfont"
} else {
  font_add("MSYaHei",
           regular = "C:/Windows/Fonts/msyh.ttc",
           italic  = "C:/Windows/Fonts/msyh.ttc",
           bold    = "C:/Windows/Fonts/msyhbd.ttc")
  fam <- "MSYaHei"
}
showtext_auto()
DPI <- 600
showtext_opts(dpi = DPI)

tip <- function(name) list(name = name, children = NULL)

# ---- 树结构定义（数据驱动：FULL 全面 / SIMPLE 简单两种密度）----------
DENS <- toupper(commandArgs(trailingOnly = TRUE)[2])
if (!DENS %in% c("FULL", "SIMPLE")) DENS <- "FULL"

tip <- function(name) list(name = name, children = NULL)

# 每行 "功能域|模块|末端功能"（模块留空 = 直接挂在功能域下）
SPEC_FULL_CN <- c(
"数据资源与获取|数据资源|公共数据检索（SRA + GSA）",
"数据资源与获取|数据资源|公共数据批量下载（aria2c + sracha）",
"数据资源与获取|数据资源|基因组/宿主基因组下载（NCBI datasets）",
"数据资源与获取|数据资源|Open-Virome 公共病毒组",
"数据资源与获取|数据资源|病毒样本预筛",
"数据资源与获取|数据资源|本地数据库构建",
"病毒检测与鉴定|样本处理|样品创建与批量导入",
"病毒检测与鉴定|样本处理|一键前处理",
"病毒检测与鉴定|样本处理|测序质控（fastp）",
"病毒检测与鉴定|样本处理|格式转换（fq2fa）",
"病毒检测与鉴定|样本处理|宿主去除（kunpeng）",
"病毒检测与鉴定|定量与共识|已知病毒识别与定量（salmon）",
"病毒检测与鉴定|定量与共识|覆盖深度与 BAM",
"病毒检测与鉴定|定量与共识|共识序列",
"病毒检测与鉴定|定量与共识|变异分析",
"病毒检测与鉴定|识别与组装|一键分析（②→④）",
"病毒检测与鉴定|识别与组装|宏病毒组组装（SPAdes）",
"病毒检测与鉴定|识别与组装|候选序列验证（BLASTx / BFVD / 类病毒）",
"病毒检测与鉴定|识别与组装|组装结果再鉴定",
"病毒检测与鉴定|识别与组装|宿主预测（ICTV 级联）",
"注释与比较基因组|基因组注释|ORF 预测（pyrodigal）",
"注释与比较基因组|基因组注释|功能注释（DIAMOND / Pfam / BFVD）",
"注释与比较基因组|基因组注释|同源比对（BLASTN / BLASTX）",
"注释与比较基因组|基因组注释|基因组图谱（gbdraw / DFV）",
"注释与比较基因组|比较基因组|参考获取与 CDS/PEP 提取",
"注释与比较基因组|比较基因组|多序列比对（MAFFT + trimAl）",
"注释与比较基因组|比较基因组|进化树构建（RAxML-NG / IQ-TREE / FastTree）",
"注释与比较基因组|比较基因组|SDT 同一性分析（精确引擎）",
"进化与动力学分析|进化动力学|元数据治理与比对 QC",
"进化与动力学分析|进化动力学|时间信号与定年",
"进化与动力学分析|进化动力学|ML 分子钟（TreeTime）",
"进化与动力学分析|进化动力学|天际线与贝叶斯天际线图",
"进化与动力学分析|进化动力学|RTT / DRT 回归",
"进化与动力学分析|进化动力学|系统地理（mugration + 迁移动画）",
"进化与动力学分析|进化动力学|重组分析（RDP5）",
"进化与动力学分析|进化动力学|BEAST 交接（Markov jump / MCC）",
"进化与动力学分析|进化动力学|一键定年→地理→溯源",
"检测与防治设计|检测与防治试剂设计|引物设计（primer3 + 热力学评分）",
"检测与防治设计|检测与防治试剂设计|CRISPR 检测试剂设计",
"检测与防治设计|检测与防治试剂设计|dsRNA 设计全链路（含脱靶评估）",
"检测与防治设计|检测与防治试剂设计|miRNA 靶标预测（5 引擎）",
"溯源与提交|溯源与提交|LOGAN 溯源（2340 万公开 SRA）",
"溯源与提交|溯源与提交|NCBI 提交准备（GenBank / BioSample / MIUVIG）",
"溯源与提交|溯源与提交|AI 助手（结果解读 · 检索执行）"
)

SPEC_FULL_EN <- c(
"Data acquisition|Data resources|Public data search (SRA + GSA)",
"Data acquisition|Data resources|Bulk download (aria2c + sracha)",
"Data acquisition|Data resources|Genome / host genome download (NCBI datasets)",
"Data acquisition|Data resources|Open-Virome public viromes",
"Data acquisition|Data resources|Sample virome prescreening",
"Data acquisition|Data resources|Local database construction",
"Virus Detection & Identification|Sample processing|Sample creation & batch import",
"Virus Detection & Identification|Sample processing|One-click preprocessing",
"Virus Detection & Identification|Sample processing|Read QC (fastp)",
"Virus Detection & Identification|Sample processing|Format conversion (fq2fa)",
"Virus Detection & Identification|Sample processing|Host read removal (kunpeng)",
"Virus Detection & Identification|Quantification & consensus|Known virus ID & quantification (salmon)",
"Virus Detection & Identification|Quantification & consensus|Coverage depth & BAM",
"Virus Detection & Identification|Quantification & consensus|Consensus sequence calling",
"Virus Detection & Identification|Quantification & consensus|Variant analysis",
"Virus Detection & Identification|Identification & assembly|One-click pipeline (stage 2-4)",
"Virus Detection & Identification|Identification & assembly|Metaviral assembly (SPAdes)",
"Virus Detection & Identification|Identification & assembly|Candidate validation (BLASTx / BFVD / viroid)",
"Virus Detection & Identification|Identification & assembly|Assembly re-identification",
"Virus Detection & Identification|Identification & assembly|Host prediction (ICTV cascade)",
"Annotation & genomics|Genome annotation|ORF prediction (pyrodigal)",
"Annotation & genomics|Genome annotation|Functional annotation (DIAMOND / Pfam / BFVD)",
"Annotation & genomics|Genome annotation|Homology search (BLASTN / BLASTX)",
"Annotation & genomics|Genome annotation|Genome maps (gbdraw / DFV)",
"Annotation & genomics|Comparative genomics|Reference acquisition & CDS/PEP export",
"Annotation & genomics|Comparative genomics|Multiple alignment (MAFFT + trimAl)",
"Annotation & genomics|Comparative genomics|Tree inference (RAxML-NG / IQ-TREE / FastTree)",
"Annotation & genomics|Comparative genomics|SDT pairwise identity (exact engine)",
"Evolution & Phylodynamics|Phylodynamics|Metadata governance & alignment QC",
"Evolution & Phylodynamics|Phylodynamics|Temporal signal & dating",
"Evolution & Phylodynamics|Phylodynamics|ML molecular clock (TreeTime)",
"Evolution & Phylodynamics|Phylodynamics|Skyline & Bayesian skyline plot",
"Evolution & Phylodynamics|Phylodynamics|RTT / DRT regression",
"Evolution & Phylodynamics|Phylodynamics|Phylogeography (mugration + migration GIF)",
"Evolution & Phylodynamics|Phylodynamics|Recombination (RDP5)",
"Evolution & Phylodynamics|Phylodynamics|BEAST interface (Markov jump / MCC)",
"Evolution & Phylodynamics|Phylodynamics|One-click dating - geography - tracing",
"Assay design|Assay & control design|Primer design (primer3 + thermodynamics)",
"Assay design|Assay & control design|CRISPR detection assay design",
"Assay design|Assay & control design|dsRNA design pipeline (off-target scoring)",
"Assay design|Assay & control design|miRNA target prediction (5 engines)",
"Tracing & submission|Tracing & submission|LOGAN tracing (23.4M public SRA runs)",
"Tracing & submission|Tracing & submission|NCBI submission (GenBank / BioSample / MIUVIG)",
"Tracing & submission|Tracing & submission|AI assistant (interpret / retrieval)"
)

SPEC_SIMPLE_CN <- c(
"数据资源与获取||公共数据检索（SRA + GSA）",
"数据资源与获取||参考与宿主基因组下载",
"数据资源与获取||病毒样本预筛",
"数据资源与获取||数据库构建",
"病毒检测与鉴定||宿主去除",
"病毒检测与鉴定||已知病毒识别与定量",
"病毒检测与鉴定||宏基因组组装",
"病毒检测与鉴定||宿主预测",
"注释与比较基因组||ORF 预测与功能注释",
"注释与比较基因组||基因组图谱",
"注释与比较基因组||进化树构建",
"注释与比较基因组||SDT 同一性分析",
"进化与动力学分析||分子钟定年",
"进化与动力学分析||系统地理分析",
"进化与动力学分析||重组分析（RDP5）",
"进化与动力学分析||天际线与 BSP",
"检测与防治设计||引物设计",
"检测与防治设计||CRISPR 检测试剂设计",
"检测与防治设计||dsRNA 设计",
"检测与防治设计||miRNA 靶标预测",
"溯源与提交||LOGAN 溯源",
"溯源与提交||NCBI 提交准备",
"溯源与提交||AI 助手"
)

SPEC_SIMPLE_EN <- c(
"Data acquisition||Public data search (SRA + GSA)",
"Data acquisition||Reference & host genome download",
"Data acquisition||Sample virome prescreening",
"Data acquisition||Database construction",
"Virus Detection & Identification||Host read removal",
"Virus Detection & Identification||Known virus ID & quantification",
"Virus Detection & Identification||Metaviral assembly",
"Virus Detection & Identification||Host prediction",
"Annotation & genomics||ORF prediction & annotation",
"Annotation & genomics||Genome maps",
"Annotation & genomics||Tree inference",
"Annotation & genomics||SDT identity analysis",
"Evolution & Phylodynamics||Molecular clock dating",
"Evolution & Phylodynamics||Phylogeography",
"Evolution & Phylodynamics||Recombination (RDP5)",
"Evolution & Phylodynamics||Skyline & BSP",
"Assay design||Primer design",
"Assay design||CRISPR assay design",
"Assay design||dsRNA design",
"Assay design||miRNA target prediction",
"Tracing & submission||LOGAN tracing",
"Tracing & submission||NCBI submission",
"Tracing & submission||AI assistant"
)
spec_to_tree <- function(rows, root_name) {
  sp <- data.frame(do.call(rbind, strsplit(rows, "|", fixed = TRUE)),
                   stringsAsFactors = FALSE)
  names(sp) <- c("domain", "module", "tip")
  sp$module[sp$module == ""] <- NA
  doms <- unique(sp$domain)
  doms_l <- lapply(doms, function(dm) {
    sub <- sp[sp$domain == dm, ]
    mods <- unique(sub$module[!is.na(sub$module)])
    ml <- lapply(mods, function(m) {
      list(name = m, children = lapply(sub$tip[sub$module == m], tip))
    })
    direct <- sub$tip[is.na(sub$module)]
    if (length(direct)) {
      ml <- c(ml, lapply(direct, function(t) {
        list(name = t, children = NULL)
      }))
    }
    list(name = dm, children = ml)
  })
  list(name = root_name, children = doms_l)
}

tree <- spec_to_tree(
  if (lang == "EN") (if (DENS == "FULL") SPEC_FULL_EN else SPEC_SIMPLE_EN)
  else              (if (DENS == "FULL") SPEC_FULL_CN else SPEC_SIMPLE_CN),
  if (lang == "EN") "VirusPlatform" else "植物病毒分析平台")

# ---- 嵌套列表 → 顶向下布局表 ------------------------------------------
# x = 末端列号（1..50，按定义顺序，内部节点取子树均值）；y = -深度 × 步长
DEPTH_STEP <- 2.0
layout_tree <- function(root) {
  nodes <- list(); nid <- 0; col <- 0
  recur <- function(nd, parent, depth) {
    nid <<- nid + 1; id <- nid
    if (is.null(nd$children)) {
      col <<- col + 1
      c0 <- c1 <- col
    } else {
      kids <- vapply(nd$children, recur, numeric(2), parent = id, depth = depth + 1)
      c0 <- min(kids[1, ]); c1 <- max(kids[2, ])
    }
    xc <- (c0 + c1) / 2
    nodes[[length(nodes) + 1]] <<- data.frame(
      node = id, parent = parent, depth = depth,
      label = nd$name, x = xc, c0 = c0, c1 = c1)
    c(c0, c1)
  }
  recur(root, NA_integer_, 0)
  do.call(rbind, nodes)
}
d <- layout_tree(tree)
maxdepth <- max(d$depth)
ntips <- sum(d$depth == maxdepth)
d$y <- -d$depth * DEPTH_STEP          # 0(根) → -6(末端)，向下生长
d$x <- d$x + 0.5                       # 列留出半列边距
cat("tips:", ntips, " nodes:", nrow(d), "\n")

# ---- 配色（Okabe-Ito 色盲安全 6 色 + 浅色底带）------------------------
domain_names <- d$label[d$depth == 1]
pal_main <- c("#0072B2", "#009E73", "#E69F00", "#CC79A7", "#D55E00", "#4D4D4D")
pal_tint <- c("#EAF2F9", "#E7F4EE", "#FBF2DF", "#F8EDF4", "#FAECE3", "#EEEEEE")
names(pal_main) <- domain_names
names(pal_tint) <- domain_names
dd <- d[d$depth == 1, ]

# ---- 画布几何 ----------------------------------------------------------
w_in <- 14.2
XMIN <- -0.2; XMAX <- ntips + 1.2      # 列坐标范围
unit_x <- (w_in - 0.223) / (XMAX - XMIN)

# 纵向分段（单位 = y 数据单位；unit_y 由总高换算）
Y_ROOT <- 0; Y_DOM <- -DEPTH_STEP; Y_MOD <- -2 * DEPTH_STEP; Y_TIP <- -3 * DEPTH_STEP
Y_TIP_LABEL <- Y_TIP - 0.55            # 末端纵排标签起点
LABEL_UNITS <- 6.4                     # 标签区高度（约 2.5in）
Y_IMG_TOP <- Y_TIP_LABEL - LABEL_UNITS - 0.5
IMG_H_UNITS <- 5.2                     # 示例图区高度（约 2.2in）
Y_CAP <- Y_IMG_TOP - IMG_H_UNITS - 0.55
Y_BOT <- Y_CAP - 0.45
h_in <- 8.9
unit_y <- (h_in - 0.223) / (1.1 - Y_BOT)
cat(sprintf("unit_x=%.3f in  unit_y=%.3f in\n", unit_x, unit_y))

# 模块标签自动换行（按子树宽度估算每行字符数）
has_cjk <- function(s) grepl("[\u4e00-\u9fff]", s)
wrap_label <- function(s, span_in) {
  ch_w <- if (has_cjk(s)) 0.105 else 0.064   # 单字符估计宽度（英寸）
  mx <- max(4, floor(span_in / ch_w))
  if (nchar(s) <= mx) return(s)
  if (has_cjk(s)) {
    chs <- strsplit(s, "")[[1]]
    lat <- grepl("[A-Za-z0-9]", chs)
    n <- length(chs)
    parts <- character(0); st <- 1
    while (st <= n) {
      en <- min(st + mx - 1, n)
      if (en < n && lat[en + 1]) {
        wb <- en
        while (wb > st && lat[wb]) wb <- wb - 1
        if (wb == st) {
          we <- en + 1
          while (we <= n && lat[we]) we <- we + 1
          en <- we - 1
        } else en <- wb
      }
      parts <- c(parts, paste(chs[st:en], collapse = ""))
      st <- en + 1
    }
    if (length(parts) >= 2 && nchar(parts[length(parts)]) <= 1) {
      prev <- parts[length(parts) - 1]
      parts[length(parts)] <- paste0(substring(prev, nchar(prev)),
                                     parts[length(parts)])
      parts[length(parts) - 1] <- substring(prev, 1, nchar(prev) - 1)
    }
    paste(parts, collapse = "\n")
  } else {
    words <- strsplit(s, " ")[[1]]
    lines <- character(0); cur <- ""
    for (wd in words) {
      cand <- if (nzchar(cur)) paste(cur, wd) else wd
      if (nchar(cand) <= mx || !nzchar(cur)) cur <- cand else { lines <- c(lines, cur); cur <- wd }
    }
    paste(c(lines, cur), collapse = "\n")
  }
}

# ---- 绘图 --------------------------------------------------------------
p <- ggplot()

# 浅色底带：每个功能域一根纵向立柱（从域行到末端行）
p <- p + geom_rect(data = data.frame(
                    xmin = dd$c0 + 0.5 - 0.45, xmax = dd$c1 + 0.5 + 0.45,
                    ymin = Y_DOM + 0.75, ymax = Y_TIP + 0.75,
                    fill = pal_tint[dd$label]),
                   aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                       fill = fill), color = NA, show.legend = FALSE) +
  scale_fill_identity()

# 直角树边：父节点竖降 + 子节点横移
edges <- d[!is.na(d$parent), ]
px <- d$x[match(edges$parent, d$node)]; py <- d$y[match(edges$parent, d$node)]
seg_v <- data.frame(x = px, y = py,     xend = px,     yend = edges$y)
seg_h <- data.frame(x = px, y = edges$y, xend = edges$x, yend = edges$y)
p <- p + geom_segment(data = seg_v, aes(x = x, y = y, xend = xend, yend = yend),
                      color = "grey35", linewidth = 0.35) +
  geom_segment(data = seg_h, aes(x = x, y = y, xend = xend, yend = yend),
               color = "grey35", linewidth = 0.35)

# 根主标题（顶部居中）
rx <- d$x[d$depth == 0]
p <- p + geom_text(data = data.frame(x = rx, y = Y_ROOT + 0.55, label = d$label[d$depth == 0]),
                   aes(x = x, y = y, label = label), size = 4.8,
                   fontface = "bold", family = fam, color = "black")

# 功能域标签（横排、域色、加粗，位于域节点上方）
p <- p + geom_label(data = data.frame(x = dd$x, y = Y_DOM + 0.55,
                                      label = dd$label, color = pal_main[dd$label]),
                    aes(x = x, y = y, label = label, color = color),
                    size = 2.9, fontface = "bold", family = fam,
                    fill = "white", linewidth = 0, alpha = 1,
                    show.legend = FALSE, na.rm = TRUE) +
  scale_color_identity()

# 模块标签（横排加粗、白底，按子树宽度自动换行）
mod <- d[d$depth == 2 & d$depth < maxdepth, ]
if (nrow(mod)) {                     # 单模块域不画模块标签（域标签即可代表）
  tb <- table(mod$parent)
  mod <- mod[mod$parent %in% names(tb)[tb >= 2], ]
}
if (nrow(mod)) {
mod$wrapped <- mapply(function(s, c0, c1) {
  span_in <- (c1 - c0 + 1) * unit_x
  wrap_label(s, span_in)
}, mod$label, mod$c0, mod$c1)
p <- p + geom_label(data = data.frame(x = mod$x, y = Y_MOD + 0.45, label = mod$wrapped),
                    aes(x = x, y = y, label = label),
                    size = 2.5, fontface = "bold", family = fam,
                    color = "#2B2B2B", fill = "white", linewidth = 0,
                    alpha = 0.9, label.padding = unit(0.12, "lines"),
                    lineheight = 0.9, na.rm = TRUE)
}

# 末端功能标签（纵排，首字在上向下读）
tps <- d[d$depth == maxdepth, ]
p <- p + geom_text(data = data.frame(x = tps$x, y = Y_TIP_LABEL, label = tps$label),
                   aes(x = x, y = y, label = label),
                   angle = 90, hjust = 1, vjust = 0.5, size = 2.3, na.rm = TRUE,
                   family = fam, color = "#1A1A1A")

# ---- 底部示例图行（与各功能域水平跨度对齐）-----------------------------
assets_def <- if (lang == "EN") list(
  list(domain = "Virus Detection & Identification", items = list(
    list(file = "figures/assets/kvsuite_depth.png", cap = "Depth plot (kvsuite)"))),
  list(domain = "Annotation & genomics", items = list(
    list(file = "figures/assets/gcva_genome_circular.png", cap = "GCVA genome map (gbdraw)"))),
  list(domain = "Evolution & Phylodynamics", items = list(
    list(file = "figures/assets/gcva_rdp_breakpoints.png", cap = "GCVA recombination (RDP5)"))),
  list(domain = "Results, Sharing & Platform", items = list(
    list(file = "figures/assets/gcva_sdt_heatmap.png", cap = "GCVA SDT identity heatmap"),
    list(file = "figures/assets/gcva_sdt_distribution.png", cap = "GCVA identity distribution")))
) else list(
  list(domain = "病毒检测与鉴定", items = list(
    list(file = "figures/assets/kvsuite_depth.png", cap = "覆盖深度图（kvsuite）"),
    list(file = "figures/assets/viromock_quant.png", cap = "检出与定量基准（VIROMOCK）"))),
  list(domain = "注释与比较基因组", items = list(
    list(file = "figures/assets/gcva_genome_circular.png", cap = "GCVA 基因组圈图（gbdraw）"))),
  list(domain = "进化与动力学分析", items = list(
    list(file = "figures/assets/gcva_rdp_breakpoints.png", cap = "GCVA 重组分析（RDP5）"))),
  list(domain = "结果、共享与平台", items = list(
    list(file = "figures/assets/gcva_sdt_heatmap.png", cap = "GCVA SDT 同一性热图"),
    list(file = "figures/assets/gcva_sdt_distribution.png", cap = "GCVA 同一性分布")))
)

for (a in assets_def) {
  r <- dd[dd$label == a$domain, ]
  items <- Filter(function(it) file.exists(it$file), a$items)
  if (nrow(r) == 0 || length(items) == 0) next
  n <- length(items)
  cell_x0 <- r$c0 + 0.5 - 0.45; cell_x1 <- r$c1 + 0.5 + 0.45
  cell_w_in <- (cell_x1 - cell_x0) * unit_x
  max_h_in <- IMG_H_UNITS * unit_y
  if (n == 1) cols <- list(c(cell_x0, cell_x1)) else {
    mid <- (cell_x0 + cell_x1) / 2
    cols <- list(c(cell_x0, mid - 0.15), c(mid + 0.15, cell_x1))
  }
  for (k in seq_len(n)) {
    it <- items[[k]]; xr <- cols[[min(k, length(cols))]]
    img <- readPNG(it$file)
    asp <- dim(img)[2] / dim(img)[1]
    sub_w_in <- (xr[2] - xr[1]) * unit_x
    h_in_fit <- min(sub_w_in / asp, max_h_in)
    w_units <- h_in_fit * asp / unit_x
    h_units <- h_in_fit / unit_y
    xc <- (xr[1] + xr[2]) / 2
    ytop <- Y_IMG_TOP
    p <- p + annotation_raster(img,
                               xmin = xc - w_units / 2, xmax = xc + w_units / 2,
                               ymin = ytop - h_units, ymax = ytop,
                               interpolate = TRUE)
    p <- p + geom_rect(data = data.frame(xmin = xc - w_units / 2, xmax = xc + w_units / 2,
                                         ymin = ytop - h_units, ymax = ytop),
                       aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                       fill = NA, color = "grey55", linewidth = 0.3)
    p <- p + geom_text(data = data.frame(x = xc, y = ytop - h_units - 0.55, label = it$cap),
                       aes(x = x, y = y, label = label), size = 2.05,
                       fontface = "italic", family = fam, color = "grey30")
  }
}

p <- p + coord_cartesian(xlim = c(XMIN, XMAX), ylim = c(Y_BOT, 1.1),
                         expand = FALSE, clip = "off") +
  theme_void(base_family = fam) +
  theme(plot.background = element_rect(fill = "white", color = NA),
        plot.margin = margin(8, 10, 8, 6))

# ---- 输出 --------------------------------------------------------------
outdir <- "figures"; dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
tag <- paste0(if (lang == "EN") "EN" else "CN", "_", tolower(DENS))

ggsave(file.path(outdir, paste0("functional_tree_topdown_", tag, ".pdf")), p,
       width = w_in, height = h_in, device = cairo_pdf)
ggsave(file.path(outdir, paste0("functional_tree_topdown_", tag, ".png")), p,
       width = w_in, height = h_in, dpi = DPI, bg = "white")
cat("saved:", file.path(outdir, paste0("functional_tree_topdown_", tag)),
    " size:", w_in, "x", h_in, "in  dpi:", DPI, "\n")
