# =====================================================================
# VirusPlatform 功能树（SCI 发表级树形思维导图）
# 用法: Rscript plot_functional_tree.R EN|CN
# 输出: figures/functional_tree_EN.pdf/.png  或  functional_tree_CN.pdf/.png
# 依赖: ggplot2 + showtext（已装）；布局由嵌套列表直接计算，不经 ape/read.tree，
#       保证「定义顺序 = 图中从上到下」
# =====================================================================

suppressPackageStartupMessages({
  library(ggplot2); library(showtext)
})

lang <- toupper(commandArgs(trailingOnly = TRUE)[1])
if (!lang %in% c("EN", "CN")) lang <- "EN"

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
DPI <- 150
showtext_opts(dpi = DPI)

tip <- function(name) list(name = name, children = NULL)

# ---- 树结构定义（根 → 6 功能域 → 模块 → 功能）-----------------------
if (lang == "EN") {
  tree <- list(name = "VirusPlatform", children = list(
    list(name = "Data Resources & Acquisition", children = list(
      list(name = "Data resources", children = list(
        tip("Public data search (SRA + GSA)"),
        tip("Bulk download (aria2c + sracha)"),
        tip("Genome / host genome download (NCBI datasets)"),
        tip("Open-Virome public viromes"),
        tip("Sample virome prescreening"),
        tip("Local database construction"))))),
    list(name = "Virus Detection & Identification", children = list(
      list(name = "Sample processing", children = list(
        tip("Sample creation & batch import"),
        tip("One-click preprocessing"),
        tip("Read QC (fastp)"),
        tip("Format conversion (fq2fa)"),
        tip("Host read removal (kunpeng)"))),
      list(name = "Quantification & consensus", children = list(
        tip("Known virus ID & quantification (salmon)"),
        tip("Coverage depth & BAM"),
        tip("Consensus sequence calling"),
        tip("Variant analysis"))),
      list(name = "Identification & assembly", children = list(
        tip("One-click pipeline (stage 2-4)"),
        tip("Metaviral assembly (SPAdes)"),
        tip("Candidate validation (BLASTx / BFVD / viroid)"),
        tip("Assembly re-identification"),
        tip("Host prediction (ICTV cascade)"))))),
    list(name = "Annotation & Comparative Genomics", children = list(
      list(name = "Genome annotation", children = list(
        tip("ORF prediction (pyrodigal)"),
        tip("Functional annotation (DIAMOND / Pfam / BFVD)"),
        tip("Homology search (BLASTN / BLASTX)"),
        tip("Genome maps (gbdraw / DFV)"))),
      list(name = "Comparative genomics", children = list(
        tip("Reference acquisition & CDS/PEP export"),
        tip("Multiple alignment (MAFFT + trimAl)"),
        tip("Tree inference (RAxML-NG / IQ-TREE / FastTree)"),
        tip("SDT pairwise identity (exact engine)"))))),
    list(name = "Evolution & Phylodynamics", children = list(
      list(name = "Phylodynamics", children = list(
        tip("Metadata governance & alignment QC"),
        tip("Temporal signal & dating"),
        tip("ML molecular clock (TreeTime)"),
        tip("Skyline & Bayesian skyline plot"),
        tip("RTT / DRT regression"),
        tip("Phylogeography (mugration + migration GIF)"),
        tip("Recombination (RDP5)"),
        tip("BEAST interface (Markov jump / MCC)"),
        tip("One-click dating - geography - tracing"))))),
    list(name = "Detection & Control", children = list(
      list(name = "Assay & control design", children = list(
        tip("Primer design (primer3 + thermodynamics)"),
        tip("CRISPR detection assay design"),
        tip("dsRNA design pipeline (off-target scoring)"),
        tip("miRNA target prediction (5 engines)"))))),
    list(name = "Results, Sharing & Platform", children = list(
      list(name = "Results & reporting", children = list(
        tip("Sample & module result center"),
        tip("Tree / SDT / MSA viewers"),
        tip("Interactive reports (plotly / pycirclize)"))),
      list(name = "Tracing & submission", children = list(
        tip("LOGAN tracing (23.4M public SRA runs)"),
        tip("NCBI submission (GenBank / BioSample / MIUVIG)"))),
      list(name = "Platform capabilities", children = list(
        tip("AI assistant (interpret / broadcast / gate)"),
        tip("Task center (self-learning ETA)"),
        tip("Full CLI (40 subcommands)"),
        tip("Bilingual interface (EN / CN)")))))))
} else {
  tree <- list(name = "植物病毒分析平台", children = list(
    list(name = "数据资源与获取", children = list(
      list(name = "数据资源", children = list(
        tip("公共数据检索（SRA + GSA）"),
        tip("公共数据批量下载（aria2c + sracha）"),
        tip("基因组/宿主基因组下载（NCBI datasets）"),
        tip("Open-Virome 公共病毒组"),
        tip("病毒样本预筛"),
        tip("本地数据库构建"))))),
    list(name = "病毒检测与鉴定", children = list(
      list(name = "样本处理", children = list(
        tip("样品创建与批量导入"),
        tip("一键前处理"),
        tip("测序质控（fastp）"),
        tip("格式转换（fq2fa）"),
        tip("宿主去除（kunpeng）"))),
      list(name = "定量与共识", children = list(
        tip("已知病毒识别与定量（salmon）"),
        tip("覆盖深度与 BAM"),
        tip("共识序列"),
        tip("变异分析"))),
      list(name = "识别与组装", children = list(
        tip("一键分析（②→④）"),
        tip("宏病毒组组装（SPAdes）"),
        tip("候选序列验证（BLASTx / BFVD / 类病毒）"),
        tip("组装结果再鉴定"),
        tip("宿主预测（ICTV 级联）"))))),
    list(name = "注释与比较基因组", children = list(
      list(name = "基因组注释", children = list(
        tip("ORF 预测（pyrodigal）"),
        tip("功能注释（DIAMOND / Pfam / BFVD）"),
        tip("同源比对（BLASTN / BLASTX）"),
        tip("基因组图谱（gbdraw / DFV）"))),
      list(name = "比较基因组", children = list(
        tip("参考获取与 CDS/PEP 提取"),
        tip("多序列比对（MAFFT + trimAl）"),
        tip("进化树构建（RAxML-NG / IQ-TREE / FastTree）"),
        tip("SDT 同一性分析（精确引擎）"))))),
    list(name = "进化与动力学分析", children = list(
      list(name = "进化动力学", children = list(
        tip("元数据治理与比对 QC"),
        tip("时间信号与定年"),
        tip("ML 分子钟（TreeTime）"),
        tip("天际线与贝叶斯天际线图"),
        tip("RTT / DRT 回归"),
        tip("系统地理（mugration + 迁移动画）"),
        tip("重组分析（RDP5）"),
        tip("BEAST 交接（Markov jump / MCC）"),
        tip("一键定年→地理→溯源"))))),
    list(name = "检测与防治设计", children = list(
      list(name = "检测与防治试剂设计", children = list(
        tip("引物设计（primer3 + 热力学评分）"),
        tip("CRISPR 检测试剂设计"),
        tip("dsRNA 设计全链路（含脱靶评估）"),
        tip("miRNA 靶标预测（5 引擎）"))))),
    list(name = "结果、共享与平台", children = list(
      list(name = "结果与报告", children = list(
        tip("样品与专项结果中心"),
        tip("树 / SDT / MSA 查看器"),
        tip("交互式报告（plotly / pycirclize）"))),
      list(name = "溯源与提交", children = list(
        tip("LOGAN 溯源（2340 万公开 SRA）"),
        tip("NCBI 提交准备（GenBank / BioSample / MIUVIG）"))),
      list(name = "平台能力", children = list(
        tip("AI 助手（解读/播报/闸门/检索执行）"),
        tip("任务中心（耗时自学习）"),
        tip("全量 CLI（40 子命令）"),
        tip("中英双语界面")))))))
}

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
                    angle = 90, hjust = 0.5, size = 2.9, fontface = "bold",
                    family = fam, alpha = 0, show.legend = FALSE, na.rm = TRUE) +
  scale_color_identity()

# 模块标签（节点右侧，白底避免与边线冲突）
mod <- d[d$depth == 2, ]
p <- p + geom_label(data = data.frame(x = mod$x + 0.06, y = mod$y, label = mod$label),
                    aes(x = x, y = y, label = label),
                    hjust = 0, vjust = 0.5, size = 2.75, fontface = "bold",
                    family = fam, color = "#2B2B2B", fill = "white",
                    linewidth = 0, alpha = 0.85, na.rm = TRUE)

# 末端功能标签
tps <- d[d$depth == maxdepth, ]
p <- p + geom_text(data = data.frame(x = tps$x + 0.08, y = tps$y, label = tps$label),
                   aes(x = x, y = y, label = label),
                   hjust = 0, vjust = 0.5, size = 2.45, na.rm = TRUE,
                   family = fam, color = "#1A1A1A")

p <- p + geom_point(data = dd, aes(x = x - 0.18, y = y), color = "red", size = 4)
p <- p + coord_cartesian(xlim = c(-0.9, 7.0), ylim = c(0.4, ntips + 0.6),
                         clip = "off") +
  theme_void(base_family = fam) +
  theme(plot.background = element_rect(fill = "white", color = NA),
        plot.margin = margin(8, 10, 8, 6))

# ---- 输出 --------------------------------------------------------------
h_in  <- 0.185 * ntips + 1.7
w_in  <- 8.6
outdir <- "figures"; dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
tag <- if (lang == "EN") "EN" else "CN"

ggsave(file.path(outdir, paste0("functional_tree_", tag, ".pdf")), p,
       width = w_in, height = h_in, device = cairo_pdf)
ggsave(file.path(outdir, paste0("functional_tree_", tag, ".png")), p,
       width = w_in, height = h_in, dpi = DPI, bg = "white")
cat("saved:", file.path(outdir, paste0("functional_tree_", tag)),
    " size:", w_in, "x", round(h_in, 2), "in  dpi:", DPI, "\n")
