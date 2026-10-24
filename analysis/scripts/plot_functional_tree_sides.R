# =====================================================================
# VirusPlatform 功能树 —— 对称双向思维导图 + 底部结果图网格
# 根居中在上，域1-2 向左、域3-6 向右展开（树冠等高对称），
# 底部 3 行 × 6 列 = 18 格平台结果图（按功能域顺序、图注域色）
# 用法: Rscript plot_functional_tree_sides.R <EN|CN> <FULL|SIMPLE>
# 输出: figures/functional_tree_sides_<lang>_<dens>.pdf/.png
# =====================================================================

suppressPackageStartupMessages({
  library(ggplot2); library(showtext); library(png)
})

lang <- toupper(commandArgs(trailingOnly = TRUE)[1])
if (!lang %in% c("EN", "CN")) lang <- "CN"
DENS <- toupper(commandArgs(trailingOnly = TRUE)[2])
if (!DENS %in% c("FULL", "SIMPLE")) DENS <- "FULL"

if (lang == "EN") {
  font_add("mainfont", regular = "C:/Windows/Fonts/arial.ttf",
           bold = "C:/Windows/Fonts/arialbd.ttf",
           italic = "C:/Windows/Fonts/ariali.ttf")
} else {
  font_add("mainfont", regular = "C:/Windows/Fonts/msyh.ttc",
           bold = "C:/Windows/Fonts/msyhbd.ttc",
           italic = "C:/Windows/Fonts/msyh.ttc")
}
fam <- "mainfont"
showtext_auto()
DPI <- 300
showtext_opts(dpi = DPI)

tip <- function(name) list(name = name, children = NULL)

# ---- 树结构（"功能域|模块|末端功能"，模块留空 = 直挂域下）-------------
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
"溯源与提交|溯源与提交|AI 助手（结果解读 · 检索执行）")

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
"Tracing & submission|Tracing & submission|AI assistant (interpret / retrieval)")

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
"溯源与提交||AI 助手")

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
"Tracing & submission||AI assistant")

spec_rows <- if (lang == "EN") {
  if (DENS == "FULL") SPEC_FULL_EN else SPEC_SIMPLE_EN
} else {
  if (DENS == "FULL") SPEC_FULL_CN else SPEC_SIMPLE_CN
}

sp <- data.frame(do.call(rbind, strsplit(spec_rows, "|", fixed = TRUE)),
                 stringsAsFactors = FALSE)
names(sp) <- c("domain", "module", "tip")

# C 区域设置下 Rscript 会把 UTF-8 源码按本地字节流解析：
# 字节本身没坏，把 Encoding 标志标回 UTF-8 即可（ASCII 串不受影响）
mark_utf8 <- function(v) { Encoding(v) <- "UTF-8"; v }
sp$domain <- mark_utf8(sp$domain)
sp$module <- mark_utf8(sp$module)
sp$tip    <- mark_utf8(sp$tip)

# ---- 对称双向布局：域1-3 走左、域4-6 走右，两半等高 --------------------
X_DOM <- 1.2; X_MOD <- 2.4; X_TIP <- 3.6; X_LAB <- 3.78
GAPY  <- 1.6

dom_levels <- unique(sp$domain)
side_of <- setNames(rep(c(-1, 1), each = 3)[seq_along(dom_levels)], dom_levels)
tpd <- table(sp$domain)[dom_levels]
SPAN <- if (DENS == "FULL") 30 else 14.5
side_tips  <- sapply(c("-1", "1"), function(s) sum(tpd[side_of[dom_levels] == as.integer(s)]))
side_ndom  <- sapply(c("-1", "1"), function(s) sum(side_of[dom_levels] == as.integer(s)))
pitch_side <- (SPAN - (side_ndom - 1) * GAPY) / side_tips

nid <- 0L; di <- 0L; d <- list()
for (s in c(-1, 1)) {
  pc <- pitch_side[[as.character(s)]]
  y_cur <- 0
  for (dm in dom_levels[side_of[dom_levels] == s]) {
    rows <- sp[sp$domain == dm, , drop = FALSE]
    n <- nrow(rows)
    ys <- y_cur + seq_len(n) * pc
    nid <- nid + 1L; dom_id <- nid
    d[[di <<- di + 1]] <- data.frame(node = dom_id, parent = NA_integer_,
      depth = 1L, label = dm, dom = dm, y = (ys[1] + ys[n]) / 2, side = s,
      x = s * X_DOM, c0 = ys[1], c1 = ys[n], stringsAsFactors = FALSE)
    parent_of_tip <- character(n)
    mod_ids <- list()
    for (m in unique(rows$module[nzchar(rows$module)])) {
      sel <- which(rows$module == m)
      nid <- nid + 1L; mod_ids[[m]] <- nid
      d[[di <<- di + 1]] <- data.frame(node = nid, parent = dom_id,
        depth = 2L, label = m, dom = dm, y = (min(ys[sel]) + max(ys[sel])) / 2,
        side = s, x = s * X_MOD, c0 = min(ys[sel]), c1 = max(ys[sel]),
        stringsAsFactors = FALSE)
      parent_of_tip[sel] <- m
    }
    for (i in seq_len(n)) {
      nid <- nid + 1L
      pid <- if (nzchar(parent_of_tip[i])) mod_ids[[parent_of_tip[i]]] else dom_id
      d[[di <<- di + 1]] <- data.frame(node = nid, parent = pid,
        depth = 3L, label = rows$tip[i], dom = dm, y = ys[i], side = s,
        x = s * X_TIP, c0 = NA_real_, c1 = NA_real_, stringsAsFactors = FALSE)
    }
    y_cur <- ys[n] + GAPY
  }
}
d <- do.call(rbind, d)
d$y  <- -d$y                     # 翻转：根在上（y>0），树向下延展
d$c0 <- -d$c0; d$c1 <- -d$c1     # 范围同步翻转（c0=上缘，c1=下缘）

mod_disp <- d[d$depth == 2, ]
if (nrow(mod_disp)) {
  tb <- table(mod_disp$parent)
  mod_disp <- mod_disp[mod_disp$parent %in% as.integer(names(tb))[tb >= 2], ]  # 单模块域不画模块标签
}

tree_name <- mark_utf8(if (lang == "CN") "VirusPlatform · 植物病毒分析平台" else "VirusPlatform · Plant Virus Analysis Platform")

# ---- 配色（Okabe-Ito 色盲安全 + 浅色底带）------------------------------
pal_main <- c("#0072B2", "#009E73", "#E69F00", "#CC79A7", "#D55E00", "#4D4D4D")
pal_tint <- c("#EAF2F9", "#E7F4EE", "#FBF2DF", "#F8EDF4", "#FAECE3", "#EEEEEE")
names(pal_main) <- dom_levels
names(pal_tint) <- dom_levels

# ---- 画布 --------------------------------------------------------------
unit_x <- 1.15
unit_y <- 0.29
XMIN <- -6.55; XMAX <- 6.55
w_in <- 0.223 + (XMAX - XMIN) * unit_x

# ---- 文本换行（中英混排，CJK 记 2 个半角宽）---------------------------
tok_w <- function(s) {
  cs <- strsplit(s, "", fixed = TRUE)[[1]]
  if (!length(cs)) return(0L)
  sum(ifelse(grepl("[\u2E80-\u9FFF\uFF00-\uFFEF]", cs, perl = TRUE), 2L, 1L))
}
wrap_label <- function(s, budget) {
  if (tok_w(s) <= budget) return(s)
  words <- strsplit(s, " ", fixed = TRUE)[[1]]
  pieces <- list()
  for (w in words) {
    while (tok_w(w) > budget) {
      acc <- 0L; cut <- 0L
      for (ci in seq_len(nchar(w))) {
        cw <- if (grepl("[\u2E80-\u9FFF\uFF00-\uFFEF]", substr(w, ci, ci),
                        perl = TRUE)) 2L else 1L
        if (acc + cw > budget) break
        acc <- acc + cw; cut <- ci
      }
      if (cut == 0L) cut <- 1L
      pieces[[length(pieces) + 1L]] <- substr(w, 1L, cut)
      w <- substr(w, cut + 1L, nchar(w))
    }
    pieces[[length(pieces) + 1L]] <- w
  }
  lines <- character(0); cur <- ""; curw <- 0L
  for (pc in pieces) {
    pw <- tok_w(pc)
    add <- if (nzchar(cur)) 1L + pw else pw
    if (curw + add > budget && nzchar(cur)) {
      lines <- c(lines, cur); cur <- pc; curw <- pw
    } else {
      cur <- if (nzchar(cur)) paste(cur, pc) else pc
      curw <- curw + add
    }
  }
  if (nzchar(cur)) lines <- c(lines, cur)
  if (length(lines) >= 2 && nchar(lines[length(lines)]) <= 2 &&
      tok_w(lines[length(lines) - 1]) > 4) {
    prev <- lines[length(lines) - 1]
    lines[length(lines)] <- paste0(substr(prev, nchar(prev), nchar(prev)),
                                   lines[length(lines)])
    lines[length(lines) - 1] <- substr(prev, 1L, nchar(prev) - 1L)
  }
  paste(lines, collapse = "\n")
}
bud_tip <- if (lang == "EN") 46 else 30
bud_mod <- if (lang == "EN") 24 else 14

# ---- 绘图 --------------------------------------------------------------
p <- ggplot()

# 功能域浅色底带（左右各立柱）
bands <- do.call(rbind, lapply(which(d$depth == 1), function(i) {
  r <- d[i, ]
  data.frame(xmin = ifelse(r$side == 1, X_DOM - 0.15, -X_TIP - 0.30),
             xmax = ifelse(r$side == 1, X_TIP + 0.30, -(X_DOM - 0.15)),
             ymin = r$c1 - 0.50, ymax = r$c0 + 0.50,
             fill = pal_tint[[r$label]], stringsAsFactors = FALSE)
}))
p <- p + geom_rect(data = bands,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
        color = NA, show.legend = FALSE) +
  scale_fill_identity()

# 直角树边（父节点竖线 + 子节点横线）
edges <- d[!is.na(d$parent), ]
px <- d$x[match(edges$parent, d$node)]; py <- d$y[match(edges$parent, d$node)]
p <- p +
  geom_segment(data = data.frame(x = px, y = py, xend = px, yend = edges$y),
               aes(x = x, y = y, xend = xend, yend = yend),
               color = "grey35", linewidth = 0.35) +
  geom_segment(data = data.frame(x = px, y = edges$y, xend = edges$x, yend = edges$y),
               aes(x = x, y = y, xend = xend, yend = yend),
               color = "grey35", linewidth = 0.35)

# 根（顶部居中）+ 中轴主干 + 域分支横线
y_dom_all <- d$y[d$depth == 1]
p <- p +
  geom_segment(data = data.frame(x = 0, y = 0.26, xend = 0, yend = min(y_dom_all)),
               aes(x = x, y = y, xend = xend, yend = yend),
               color = "grey35", linewidth = 0.5) +
  geom_segment(data = data.frame(x = 0, y = y_dom_all,
                                 xend = d$side[d$depth == 1] * X_DOM,
                                 yend = y_dom_all),
               aes(x = x, y = y, xend = xend, yend = yend),
               color = "grey35", linewidth = 0.5) +
  geom_text(data = data.frame(x = 0, y = 0.62, label = tree_name),
            aes(x = x, y = y, label = label), size = 4.6,
            fontface = "bold", family = fam, color = "black")

# 末端功能标签（水平向外，域色）
tips <- d[d$depth == 3, ]
tips$wrapped <- vapply(tips$label, wrap_label, character(1), budget = bud_tip)
p <- p + geom_text(data = tips,
        aes(x = side * X_LAB, y = y, label = wrapped),
        family = fam, size = 2.6, lineheight = 0.9,
        hjust = ifelse(tips$side == 1, 0, 1),
        color = pal_main[tips$dom])

# 模块标签（白色小牌压在总线上；单模块域不显示）
if (nrow(mod_disp)) {
  mod_disp$wrapped <- vapply(mod_disp$label, wrap_label, character(1),
                             budget = bud_mod)
  p <- p + geom_label(data = mod_disp,
          aes(x = x, y = y, label = wrapped),
          family = fam, size = 2.5, lineheight = 0.9,
          color = "grey15", fill = "white", fontface = "plain",
          label.size = 0.22, label.padding = unit(0.18, "lines"))
}

# 功能域标签（白色粗体大牌压在域总线中点）
doms <- d[d$depth == 1, ]
p <- p + geom_label(data = doms,
        aes(x = x, y = y, label = label),
        family = fam, size = 3.15, lineheight = 0.9,
        color = pal_main[doms$label], fill = "white",
        fontface = "bold", label.size = 0.30,
        label.padding = unit(0.28, "lines"))

# ---- 底部结果图网格：3 行 × 6 列 = 18 格（按功能域分组，图注域色）------
img_root <- "figures/assets/grid"
IMG_TXT <- if (lang == "CN") c(
  homepage                = "平台总览首页",
  meta_temporal           = "公共数据时间分布",
  meta_regions            = "公共数据地域分布",
  meta_organizations      = "公共数据机构分布",
  meta_stages             = "公共数据生育期分布",
  prescreen_hist          = "预筛·病毒占比分布",
  prescreen_scatter       = "预筛·测序量×病毒信号",
  prescreen_taxa          = "预筛·Top 病毒类群",
  prescreen_funnel        = "预筛漏斗与未决构成",
  stage_funnel            = "分析流程漏斗",
  kvsuite_depth           = "病毒覆盖深度",
  viromock_quant          = "病毒定量验证",
  variant_overview        = "变异全景总览",
  variant_map             = "变异基因组分布",
  variant_effect          = "变异效应预测",
  classification_overview = "候选验证判定构成",
  gcva_genome_circular    = "GCVA-32 基因组圈图",
  gcva_sdt_heatmap        = "SDT 同一性热图",
  gcva_identity_composite = "同一性矩阵组合",
  gcva_nt_heatmap         = "NT 同一性热图",
  gcva_tree               = "NJ 系统进化树",
  identity_matrix_heatmap = "结构比较·同一性矩阵",
  gcva_rdp_breakpoints    = "RDP5 重组断点",
  gcva_sdt_distribution   = "同一性分布",
  rdp_events_overview     = "RDP5 事件与断点区间",
  align_qc_overview       = "比对预筛 QC 分布",
  metadata_coverage_overview = "元数据治理覆盖",
  mot_route_year_heatmap  = "MOT 迁移路线×年",
  phylogeo_migration_overview = "系统地理迁移概览",
  migration_map           = "迁移地图（中国）",
  primer_design_overview  = "引物设计概览",
  crispr_design_overview  = "CRISPR 候选评分",
  crispr_top_sets         = "CRISPR Top 引导组合",
  mirna_targets_overview  = "miRNA 靶标概览",
  mirna_relations         = "miRNA–蛋白关系图",
  mirna_duplex            = "miRNA 双链配对结构",
  logan_hits              = "LOGAN 溯源命中分布（示例）"
) else c(
  homepage                = "Platform home",
  meta_temporal           = "Public data by year",
  meta_regions            = "Public data by region",
  meta_organizations      = "Public data by organization",
  meta_stages             = "Public data by growth stage",
  prescreen_hist          = "Prescreen: viral fraction",
  prescreen_scatter       = "Prescreen: depth vs signal",
  prescreen_taxa          = "Prescreen: top taxa",
  prescreen_funnel        = "Prescreen funnel + unresolved",
  stage_funnel            = "Pipeline funnel",
  kvsuite_depth           = "Viral coverage depth",
  viromock_quant          = "Viral quantification (mock)",
  variant_overview        = "Variant landscape",
  variant_map             = "Variants along genome",
  variant_effect          = "Variant effect prediction",
  classification_overview = "Verification verdicts",
  gcva_genome_circular    = "GCVA-32 genome map",
  gcva_sdt_heatmap        = "SDT identity heatmap",
  gcva_identity_composite = "Identity matrix composite",
  gcva_nt_heatmap         = "NT identity heatmap",
  gcva_tree               = "NJ phylogenetic tree",
  identity_matrix_heatmap = "Pairwise identity matrix",
  gcva_rdp_breakpoints    = "RDP5 breakpoints",
  gcva_sdt_distribution   = "Identity distribution",
  rdp_events_overview     = "RDP5 events & breakpoints",
  align_qc_overview       = "Alignment prescreen QC",
  metadata_coverage_overview = "Metadata coverage",
  mot_route_year_heatmap  = "MOT route × year",
  phylogeo_migration_overview = "Phylogeography overview",
  migration_map           = "Migration map (China)",
  primer_design_overview  = "Primer design overview",
  crispr_design_overview  = "CRISPR candidate scores",
  crispr_top_sets         = "CRISPR top guide sets",
  mirna_targets_overview  = "miRNA targets overview",
  mirna_relations         = "miRNA–protein relations",
  mirna_duplex            = "miRNA:target duplex",
  logan_hits              = "LOGAN tracing hits (demo)")
IMG_TXT <- mark_utf8(IMG_TXT)

img_dom <- c(rep(dom_levels[1], 9), rep(dom_levels[2], 7),
             rep(dom_levels[3], 6), rep(dom_levels[4], 8),
             rep(dom_levels[5], 6), rep(dom_levels[6], 1))
NCOL_ <- 6L
col_pitch <- 2.1; img_w <- 1.92
img_dy <- 4.9; cap_dy <- 1.05; row_gap <- 0.55
row_pitch <- img_dy + cap_dy + row_gap
Y_GRID_TOP <- min(d$y) - 1.3

imgdf <- data.frame(name = names(IMG_TXT), cap = unname(IMG_TXT),
                    dom = img_dom, stringsAsFactors = FALSE)
imgdf$path <- file.path(img_root, paste0(imgdf$name, ".png"))
stopifnot(all(file.exists(imgdf$path)))
for (k in seq_len(nrow(imgdf))) {
  r <- (k - 1L) %/% NCOL_; cc <- (k - 1L) %% NCOL_
  imgdf$x0[k] <- -6.3 + cc * col_pitch
  imgdf$x1[k] <- imgdf$x0[k] + img_w
  imgdf$y1[k] <- Y_GRID_TOP - r * row_pitch
  imgdf$y0[k] <- imgdf$y1[k] - img_dy
  imgdf$cy[k] <- imgdf$y0[k] - 0.62
}
Y_BOT <- min(imgdf$y0) - cap_dy

for (k in seq_len(nrow(imgdf))) {
  im <- readPNG(imgdf$path[k])
  a <- dim(im)[2] / dim(im)[1]                 # 像素宽高比
  w_inch <- img_w * unit_x; h_inch <- img_dy * unit_y
  w2 <- min(w_inch, h_inch * a); h2 <- w2 / a   # 单元内按比例适配
  cx <- (imgdf$x0[k] + imgdf$x1[k]) / 2
  cy <- (imgdf$y0[k] + imgdf$y1[k]) / 2
  p <- p + annotation_raster(im,
        cx - (w2 / unit_x) / 2, cx + (w2 / unit_x) / 2,
        cy - (h2 / unit_y) / 2, cy + (h2 / unit_y) / 2,
        interpolate = TRUE)
}
p <- p +
  geom_rect(data = imgdf,
            aes(xmin = x0, xmax = x1, ymin = y0, ymax = y1),
            fill = NA, color = "grey78", linewidth = 0.25) +
  geom_text(data = imgdf,
            aes(x = (x0 + x1) / 2, y = cy, label = cap),
            family = fam, size = 2.35, fontface = "bold",
            color = pal_main[imgdf$dom])

# ---- 画幅与输出 --------------------------------------------------------
h_in <- 0.223 + (1.15 - Y_BOT) * unit_y
p <- p + coord_cartesian(xlim = c(XMIN, XMAX),
                         ylim = c(Y_BOT - 0.35, 1.15), expand = FALSE) +
  theme_void() +
  theme(text = element_text(family = fam))
cat(sprintf("unit_x=%.2f unit_y=%.3f w=%.1f h=%.1f\n",
            unit_x, unit_y, w_in, h_in))

dir.create("figures", showWarnings = FALSE)
base <- file.path("figures",
                  sprintf("functional_tree_sides_%s_%s", lang, DENS))
ggsave(paste0(base, ".pdf"), p, width = w_in, height = h_in,
       device = cairo_pdf, limitsize = FALSE)
ggsave(paste0(base, ".png"), p, width = w_in, height = h_in,
       dpi = 300, limitsize = FALSE)
cat("OK:", base, "\n")
