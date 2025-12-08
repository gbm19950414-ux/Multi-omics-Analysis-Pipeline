#!/usr/bin/env Rscript

# scripts/plot/figure_7_a.R
#
# 目的：生成“发表级” Figure 7A 流程图
#       概念化展示 EphB1 中介分析的整体 pipeline
#       （不包含具体脚本名 / 文件路径，只保留生物学与分析逻辑）
#
# 输出：
#   results/figs/figure_7_a.pdf
#
# 使用方式：
#   Rscript scripts/plot/figure_7_a.R

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

# -------- 1. 输出路径 --------
out_pdf <- "results/figs/figure_7_a.pdf"

# -------- 2. 定义节点（boxes） --------
# 坐标是人为指定的逻辑布局，可在需要时微调
nodes <- tibble::tibble(
  id = c(
    "genotype",
    "rna_seq", "lipidomics",
    "pathway_scores", "lipid_indicators",
    "integration",
    "mediation"
  ),
  label = c(
    # X：设计层面
    "Genotype\nWT vs Ephb1−/−",
    # 组学测量
    "RNA-seq\nperitoneal macrophages",
    "Lipidomics\nmitochondrial / membrane lipids",
    # 中间层：特征构建
    "EphB1 downstream signalling modules (M)\npathway-level activity scores",
    "CL-related lipid phenotypes (Y)\naxes & indicators (sum CL, CL fraction,\nmito lipid mass, PC/PE/CL)",
    # 样本整合
    "Sample-level integration\npaired RNA + lipid per mouse",
    # 统计模型
    "Mediation analysis\nX (genotype) → M (pathway) → Y (axis / lipid)\nbootstrap a×b, total/direct effects (c, c')"
  ),
  # 手动画一个“左右流动”的 pipeline 布局
  x = c(
    0,   # genotype
    2, 2, # RNA-seq, lipidomics
    4, 4, # pathway_scores, lipid_indicators
    6,    # integration
    8     # mediation
  ),
  y = c(
    2,   # genotype
    3, 1, # RNA-seq, lipidomics
    3, 1, # pathway_scores, lipid_indicators
    2,    # integration
    2     # mediation
  )
)

box_width  <- 1.8
box_height <- 0.9

nodes <- nodes %>%
  mutate(
    xmin = x - box_width / 2,
    xmax = x + box_width / 2,
    ymin = y - box_height / 2,
    ymax = y + box_height / 2
  )

# -------- 3. 定义边（arrows） --------
edges <- tibble::tibble(
  x    = c(
    # 设计 → 测量
    nodes$x[nodes$id == "genotype"],
    nodes$x[nodes$id == "genotype"],
    # RNA → 通路分数
    nodes$x[nodes$id == "rna_seq"],
    # Lipid → 脂质指标
    nodes$x[nodes$id == "lipidomics"],
    # 通路 + 指标 → 整合
    nodes$x[nodes$id == "pathway_scores"],
    nodes$x[nodes$id == "lipid_indicators"],
    # 整合 → 中介模型
    nodes$x[nodes$id == "integration"]
  ),
  y    = c(
    nodes$y[nodes$id == "genotype"],
    nodes$y[nodes$id == "genotype"],
    nodes$y[nodes$id == "rna_seq"],
    nodes$y[nodes$id == "lipidomics"],
    nodes$y[nodes$id == "pathway_scores"],
    nodes$y[nodes$id == "lipid_indicators"],
    nodes$y[nodes$id == "integration"]
  ),
  xend = c(
    nodes$x[nodes$id == "rna_seq"],
    nodes$x[nodes$id == "lipidomics"],
    nodes$x[nodes$id == "pathway_scores"],
    nodes$x[nodes$id == "lipid_indicators"],
    nodes$x[nodes$id == "integration"],
    nodes$x[nodes$id == "integration"],
    nodes$x[nodes$id == "mediation"]
  ),
  yend = c(
    nodes$y[nodes$id == "rna_seq"],
    nodes$y[nodes$id == "lipidomics"],
    nodes$y[nodes$id == "pathway_scores"],
    nodes$y[nodes$id == "lipid_indicators"],
    nodes$y[nodes$id == "integration"],
    nodes$y[nodes$id == "integration"],
    nodes$y[nodes$id == "mediation"]
  )
)

# 为了让箭头不要正好贴在 box 边上，做一点点缩短
shrink_edge <- function(x, y, xend, yend, shrink = 0.25) {
  dx <- xend - x
  dy <- yend - y
  len <- sqrt(dx^2 + dy^2)
  # 避免除零
  len[len == 0] <- 1
  x_new    <- x    + dx / len * shrink
  y_new    <- y    + dy / len * shrink
  xend_new <- xend - dx / len * shrink
  yend_new <- yend - dy / len * shrink
  list(x = x_new, y = y_new, xend = xend_new, yend = yend_new)
}

sh <- shrink_edge(edges$x, edges$y, edges$xend, edges$yend, shrink = 0.3)
edges$x    <- sh$x
edges$y    <- sh$y
edges$xend <- sh$xend
edges$yend <- sh$yend

# -------- 4. 绘图 --------
p <- ggplot() +
  # 箭头
  geom_segment(
    data = edges,
    aes(x = x, y = y, xend = xend, yend = yend),
    arrow = grid::arrow(length = grid::unit(0.18, "cm"), type = "closed"),
    linewidth = 0.4
  ) +
  # box
  geom_rect(
    data = nodes,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "white",
    linewidth = 0.4
  ) +
  # 文本
  geom_text(
    data = nodes,
    aes(x = x, y = y, label = label),
    size = 3,
    lineheight = 0.9
  ) +
  coord_equal(xlim = c(-0.5, 8.5), ylim = c(0, 4), expand = FALSE) +
  theme_void(base_size = 11) +
  labs(
    title = "Figure 7A – Conceptual pipeline for EphB1-dependent mediation analysis"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

message("[Figure 7A] 导出发表级流程图: ", out_pdf)
ggsave(out_pdf, p, width = 8, height = 4.5, device = cairo_pdf)

message("[Figure 7A] 完成。")