# GCVA-32 进化树（平台 quicktree 产物 tree.nwk 的发表级渲染）
suppressPackageStartupMessages({
  library(ggtree); library(ape); library(ggplot2); library(showtext)
})
font_add("mainfont", regular = "C:/Windows/Fonts/arial.ttf",
         bold = "C:/Windows/Fonts/arialbd.ttf")
showtext_auto(); showtext_opts(dpi = 300)

nwk <- "run/tool_runs/quicktree_20261012_051924_9812f6/tree.nwk"
tr <- read.tree(nwk)
tr <- ladderize(tr)
# 去掉建树脚本加的后缀，保留 run 编号
tr$tip.label <- sub("^(\\w+)_OR489165\\.1.*$", "\\1", tr$tip.label)
tr$tip.label <- sub("_+", "_", tr$tip.label)

p <- ggtree(tr, layout = "rectangular", color = "grey30", size = 0.5) +
  geom_tiplab(size = 2.6, family = "mainfont", color = "#1A1A1A", offset = 0.2) +
  xlim(0, 0.17) +
  theme_tree2() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(family = "mainfont", face = "bold", size = 13),
        axis.text.x = element_text(size = 6)) +
  labs(title = "GCVA-32 · Cytorhabdovirus sp. lycii (OR489165.1) — 32 genomes, FastTree NJ",
       x = "Substitutions per site")

ggsave("figures/assets/gcva_tree.png", p, width = 7.5, height = 5.6, dpi = 300, bg = "white")
cat("saved gcva_tree.png\n")
