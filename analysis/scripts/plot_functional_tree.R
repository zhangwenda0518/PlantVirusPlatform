# =====================================================================
# VirusPlatform 功能树（SCI 发表级树形思维导图）
# 用法: Rscript plot_functional_tree.R EN|CN
# 输出: figures/functional_tree_EN.pdf/.png  或  functional_tree_CN.pdf/.png
# 依赖: ggplot2 + showtext（已装）；布局由嵌套列表直接计算，不经 ape/read.tree，
#       保证「定义顺序 = 图中从上到下」
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

# ---- 嵌套列表 → 布局表（y 翻转后：定义顺序 = 图中从上到下）------------
layout_tree <- function(root) {
  nodes <- list(); nid <- 0; y <- 0
  recur <- function(nd, parent, depth) {
    nid <<- nid + 1; id <- nid
    if (is.null(nd$children)) {
      y <<- y + 1
      ymin <- ymax <- y
    } else {
      kys <- vapply(nd$children, recur, numeric(2), parent = id, depth = depth + 1)
      ymin <- min(kys[1, ]); ymax <- max(kys[2, ])
    }
    yc <- (ymin + ymax) / 2
    nodes[[length(nodes) + 1]] <<- data.frame(
      node = id, parent = parent, depth = depth,
      label = nd$name, y = yc, ymin = ymin, ymax = ymax)
    c(ymin, ymax)  # 返回子树末端范围，供父节点聚合
  }
  recur(root, NA_integer_, 0)
  do.call(rbind, nodes)
}
d <- layout_tree(tree)
maxdepth <- max(d$depth)
ntips <- sum(d$depth == maxdepth)
# 模块列与末端列间距加倍（x: 根0 → 域1 → 模块2 → 功能4），给标签留空间
d$x <- ifelse(d$depth == maxdepth, 4, d$depth)
# 翻转 y：让先定义的分支显示在图的顶部
Y <- ntips
oy <- d$y; oymin <- d$ymin; oymax <- d$ymax
d$y    <- Y + 1 - oy
d$ymin <- Y + 1 - oymax
d$ymax <- Y + 1 - oymin
cat("tips:", ntips, " nodes:", nrow(d), "\n")

# ---- 配色（Okabe-Ito 色盲安全 6 色 + 对应浅色底带）--------------------
domain_names <- d$label[d$depth == 1]
pal_main <- c("#0072B2", "#009E73", "#E69F00", "#CC79A7", "#D55E00", "#4D4D4D")
pal_tint <- c("#EAF2F9", "#E7F4EE", "#FBF2DF", "#F8EDF4", "#FAECE3", "#EEEEEE")
names(pal_main) <- domain_names
names(pal_tint) <- domain_names
dd <- d[d$depth == 1, ]   # 功能域行

# ---- 绘图 -------------------------------------------------------------
p <- ggplot()

# 浅色底带：每个功能域从其节点 x 延伸到末端
hi <- data.frame(xmin = dd$x, xmax = max(d$x),
                 ymin = dd$ymin - 0.42, ymax = dd$ymax + 0.42,
                 fill = pal_tint[dd$label])
p <- p + geom_rect(data = hi,
                   aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                       fill = fill), color = NA, show.legend = FALSE) +
  scale_fill_identity()

# 直角树边：父节点竖线 + 子节点横线
edges <- d[!is.na(d$parent), ]
px <- d$x[match(edges$parent, d$node)]; py <- d$y[match(edges$parent, d$node)]
seg_v <- data.frame(x = px, y = py,     xend = px,     yend = edges$y)
seg_h <- data.frame(x = px, y = edges$y, xend = edges$x, yend = edges$y)
p <- p + geom_segment(data = seg_v, aes(x = x, y = y, xend = xend, yend = yend),
                      color = "grey35", linewidth = 0.35) +
  geom_segment(data = seg_h, aes(x = x, y = y, xend = xend, yend = yend),
               color = "grey35", linewidth = 0.35)

# 根部纵排主标题（思维导图根）
rx <- d$x[d$depth == 0]
p <- p + geom_text(data = data.frame(x = rx - 0.45, y = (min(d$y) + max(d$y)) / 2,
                                     label = d$label[d$depth == 0]),
                   aes(x = x, y = y, label = label),
                   angle = 90, hjust = 0.5, vjust = 0.5, size = 4.6,
                   fontface = "bold", family = fam, color = "black")

# 功能域纵排标签（域色、加粗、白底遮住穿过的边线）
p <- p + geom_label(data = data.frame(x = dd$x - 0.18, y = dd$y,
                                      label = dd$label, color = pal_main[dd$label]),
                    aes(x = x, y = y, label = label, color = color),
                    angle = 90, hjust = 0.5, size = 2.8, fontface = "bold",
                    family = fam, fill = "white", linewidth = 0,
                    alpha = 1, show.legend = FALSE, na.rm = TRUE) +
  scale_color_identity()

# 模块标签（节点右侧，白底避免与边线冲突）
mod <- d[d$depth == 2 & d$depth < maxdepth, ]
if (nrow(mod)) {                     # 单模块域不画模块标签（域标签即可代表）
  tb <- table(mod$parent)
  mod <- mod[mod$parent %in% names(tb)[tb >= 2], ]
}
if (nrow(mod)) {
p <- p + geom_label(data = data.frame(x = mod$x + 0.06, y = mod$y, label = mod$label),
                    aes(x = x, y = y, label = label),
                    hjust = 0, vjust = 0.5, size = 2.6, fontface = "bold",
                    family = fam, color = "#2B2B2B", fill = "white",
                    linewidth = 0, alpha = 0.85, na.rm = TRUE)
}

# 末端功能标签
tps <- d[d$depth == maxdepth, ]
p <- p + geom_text(data = data.frame(x = tps$x + 0.08, y = tps$y, label = tps$label),
                   aes(x = x, y = y, label = label),
                   hjust = 0, vjust = 0.5, size = 2.45, na.rm = TRUE,
                   family = fam, color = "#1A1A1A")

# ---- 画布尺寸 ----------------------------------------------------------
h_in  <- 0.185 * ntips + 1.7
w_in  <- 11.6
XMAX  <- 9.6

# ---- 画布尺寸 ----------------------------------------------------------
h_in  <- 0.185 * ntips + 1.7
w_in  <- 13.6
XMAX  <- 12.4

# ---- 右侧示例图列（平台真实输出，按功能域对齐；每域 1-2 张并排）---------
AX0 <- 6.6; AX1 <- 12.1   # 缩略图列左右边界（数据坐标）
unit_w <- (w_in - 0.223) / (XMAX + 0.9)   # 每数据单元的英寸宽
unit_h <- (h_in - 0.223) / (ntips + 0.2)          # 每数据单元的英寸高

assets_def <- if (lang == "EN") list(
  list(domain = "Virus Detection & Identification", items = list(
    list(file = "figures/assets/kvsuite_depth.png", cap = "Depth plot (kvsuite)"))),
  list(domain = "Annotation & genomics", items = list(
    list(file = "figures/assets/gcva_genome_circular.png", cap = "GCVA genome map (gbdraw)"),
    list(file = "figures/assets/synteny_compare.png", cap = "Synteny comparison"))),
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
    list(file = "figures/assets/gcva_genome_circular.png", cap = "GCVA 基因组圈图（gbdraw）"),
    list(file = "figures/assets/synteny_compare.png", cap = "共线性比较"))),
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
  ytop <- r$ymax - 0.45; ybot0 <- r$ymin + 0.35   # 面板可用范围
  avail_h <- (ytop - ybot0 - 0.7) * unit_h        # 英寸（预留图注）
  # 子列划分：1 张占满整列，2 张各占一半（留 0.12u 列间距）
  if (n == 1) {
    cols <- list(c(AX0, AX1))
  } else {
    mid <- (AX0 + AX1) / 2
    cols <- list(c(AX0, mid - 0.06), c(mid + 0.06, AX1))
  }
  for (k in seq_len(n)) {
    it <- items[[k]]; xr <- cols[[min(k, length(cols))]]
    img <- readPNG(it$file)
    asp <- dim(img)[2] / dim(img)[1]              # 宽/高（dim = 行,列,通道）
    col_w_in <- (xr[2] - xr[1]) * unit_w
    h_in_fit <- min(col_w_in / asp, avail_h)
    h_units <- h_in_fit / unit_h
    w_units <- h_in_fit * asp / unit_w
    xc <- (xr[1] + xr[2]) / 2
    ymid <- (ytop + ybot0) / 2
    p <- p + annotation_raster(img,
                               xmin = xc - w_units / 2, xmax = xc + w_units / 2,
                               ymin = ymid - h_units / 2, ymax = ymid + h_units / 2,
                               interpolate = TRUE)
    # 边框：显式 data.frame 传值，避免 aes 懒求值吃到循环变量的最终值
    p <- p + geom_rect(data = data.frame(xmin = xc - w_units / 2,
                                         xmax = xc + w_units / 2,
                                         ymin = ymid - h_units / 2,
                                         ymax = ymid + h_units / 2),
                       aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                       fill = NA, color = "grey55", linewidth = 0.3)
    p <- p + geom_text(data = data.frame(x = xc, y = ymid - h_units / 2 - 0.45,
                                         label = it$cap),
                       aes(x = x, y = y, label = label), size = 2.05,
                       fontface = "italic", family = fam, color = "grey30")
  }
}

p <- p + coord_cartesian(xlim = c(-0.9, XMAX), ylim = c(0.4, ntips + 0.6),
                         expand = FALSE, clip = "off") +
  theme_void(base_family = fam) +
  theme(plot.background = element_rect(fill = "white", color = NA),
        plot.margin = margin(8, 10, 8, 6))

# ---- 输出 --------------------------------------------------------------
outdir <- "figures"; dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
tag <- paste0(if (lang == "EN") "EN" else "CN", "_", tolower(DENS))

ggsave(file.path(outdir, paste0("functional_tree_", tag, ".pdf")), p,
       width = w_in, height = h_in, device = cairo_pdf)
ggsave(file.path(outdir, paste0("functional_tree_", tag, ".png")), p,
       width = w_in, height = h_in, dpi = DPI, bg = "white")
cat("saved:", file.path(outdir, paste0("functional_tree_", tag)),
    " size:", w_in, "x", round(h_in, 2), "in  dpi:", DPI, "\n")
