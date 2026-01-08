#!/usr/bin/env Rscript

## ============================================================
## 06_ephb1_downstream_plots.R
##
## 目的：
##   基于 05_rna_ephb1_downstream_scores.R 的结果表，
##   可视化 “EphB1 下游信号通路在各 batch 中的变化”：
##     1) 按跨 batch 综合效应 z_meta 排序的条形图
##     2) pathway × batch 的 z_pathway 热图
##
## 输入（默认，可通过 config 覆盖）：
##
##   1) ephb1_downstream_signaling_scores_per_batch.tsv
##        results/multiomics/tables/ephb1_downstream_signaling_scores_per_batch.tsv
##
##   2) 可选配置文件（命令行第 1 个参数）：
##        scripts/plot/figure_7_c.yaml
##
##      支持字段：
##        ephb1_downstream_scores_per_batch: "results/multiomics/tables/ephb1_downstream_signaling_scores_per_batch.tsv"
##        plots_dir: "results/multiomics/plots"
##
## 输出：
##
##   1) EphB1 下游通路整体条形图（跨 batch 综合 Z）：
##        results/multiomics/figs/figure_7_c_ephb1_downstream_bar_overall.png
##
##   2) EphB1 下游通路 pathway × batch 热图：
##        results/multiomics/figs/figure_7_c_ephb1_downstream_heatmap_per_batch.png
##
## 说明：
##   - z_pathway > 0：该 EphB1 下游通路在 HO 中整体激活；
##   - z_pathway < 0：该通路在 HO 中整体受抑；
##   - z_meta：对每条通路跨 batch 进行 Stouffer 合并的综合 Z，
##             绝对值越大，说明该通路在所有 batch 中变化越稳定、效应越强。
##
## ============================================================

suppressPackageStartupMessages({
  library(yaml)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(stringr)
  library(forcats)
  library(grDevices)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

## --------- style helpers (Nature-like) ----------------------
pt_to_mm <- function(pt) pt * 25.4 / 72

get_yaml_value <- function(x, path, default = NULL) {
  # path: character vector, e.g. c("typography","sizes_pt","axis_tick_default")
  cur <- x
  for (k in path) {
    if (is.null(cur) || is.null(cur[[k]])) return(default)
    cur <- cur[[k]]
  }
  cur
}

## --------- 0. 解析命令行参数 & 配置 -------------------------

args <- commandArgs(trailingOnly = TRUE)

config_path <- if (length(args) >= 1) {
  args[1]
} else {
  "scripts/plot/figure_7_c.yaml"
}

if (file.exists(config_path)) {
  cat("[INFO] 使用配置文件: ", config_path, "\n", sep = "")
  cfg <- yaml::read_yaml(config_path)
} else {
  cat("[WARN] 未找到配置文件: ", config_path,
      "，将使用脚本内置默认路径。\n", sep = "")
  cfg <- list()
}

scores_path_default <- "results/multiomics/tables/ephb1_downstream_signaling_scores_per_batch.tsv"
plots_dir_default   <- "results/figs"

style_yaml_default  <- "/Volumes/Samsung_SSD_990_PRO_2TB_Media/EphB1/02_protocols/figure_style_nature.yaml"
max_width_mm_default <- 84
combined_width_mm_default <- 178
combined_height_mm_default <- 118

scores_path <- cfg$ephb1_downstream_scores_per_batch %||% scores_path_default
plots_dir   <- cfg$plots_dir                          %||% plots_dir_default
pathway_rename <- cfg$pathway_rename %||% list()

style_yaml_path <- cfg$style_yaml %||% style_yaml_default
max_width_mm    <- cfg$max_width_mm %||% max_width_mm_default
combined_width_mm  <- cfg$combined_width_mm  %||% combined_width_mm_default
combined_height_mm <- cfg$combined_height_mm %||% combined_height_mm_default

if (file.exists(style_yaml_path)) {
  cat("[INFO] 使用 figure style YAML: ", style_yaml_path, "\n", sep = "")
  fig_style <- yaml::read_yaml(style_yaml_path)
} else {
  cat("[WARN] 未找到 figure style YAML: ", style_yaml_path,
      "\n       将使用脚本内置 ggplot2 默认样式。\n", sep = "")
  fig_style <- list()
}

# 强制单栏宽度限制（Nature 常见单栏 ~84 mm）
max_width_mm <- min(as.numeric(max_width_mm), 84)

cat("============================================================\n")
cat("[INFO] 06_ephb1_downstream_plots.R\n")
cat("  scores_path     : ", scores_path, "\n", sep = "")
cat("  plots_dir       : ", plots_dir,   "\n", sep = "")
cat("  style_yaml_path : ", style_yaml_path, "\n", sep = "")
cat("  max_width_mm    : ", max_width_mm, "\n", sep = "")
cat("  combined_width_mm  : ", combined_width_mm, "\n", sep = "")
cat("  combined_height_mm : ", combined_height_mm, "\n", sep = "")
cat("  pathway_rename : ", length(pathway_rename), " entries\n", sep = "")
cat("============================================================\n\n")

## --------- 1. 读取 EphB1 下游通路 per-batch Z 表 ------------

if (!file.exists(scores_path)) {
  stop("[ERROR] 找不到 ephb1_downstream_signaling_scores_per_batch.tsv: ",
       scores_path, "\n请先运行 05_rna_ephb1_downstream_scores.R")
}

cat("[STEP1] 读取 EphB1 下游通路 per-batch Z 表 ...\n")

scores <- readr::read_tsv(scores_path, show_col_types = FALSE)

cat("  [scores] nrow = ", nrow(scores),
    " ; ncol = ", ncol(scores), "\n", sep = "")

required_cols <- c("pathway", "axis", "batch",
                   "n_genes", "mean_LFC", "sd_LFC", "se_LFC", "z_pathway")

if (!all(required_cols %in% colnames(scores))) {
  stop("[ERROR] scores 表缺少必要列。缺少: ",
       paste(setdiff(required_cols, colnames(scores)), collapse = ", "),
       "\n当前列名: ", paste(colnames(scores), collapse = ", "))
}

## 过滤掉没有基因 / 没有 Z 的组合
scores_use <- scores %>%
  filter(n_genes > 0, !is.na(z_pathway))

# 可选：对通路名称做重命名（用于出图显示）
# 在 figure_7_c.yaml 里用 `pathway_rename:` 提供映射：旧名: 新名
if (length(pathway_rename) > 0) {
  rename_vec <- unlist(pathway_rename)
  idx <- match(scores_use$pathway, names(rename_vec))
  scores_use$pathway <- ifelse(!is.na(idx), as.character(rename_vec[idx]), scores_use$pathway)
}

cat("  [scores_use] nrow = ", nrow(scores_use),
    " ; pathway 数 = ", length(unique(scores_use$pathway)),
    " ; batch 数 = ", length(unique(scores_use$batch)), "\n", sep = "")

## --------- 1.1 统一 Nature-like 主题 ------------------------

style_font_family <- get_yaml_value(fig_style, c("typography","font_family_primary"), "Helvetica")

pt_axis_tick   <- as.numeric(get_yaml_value(fig_style, c("typography","sizes_pt","axis_tick_default"), 5.5))
pt_axis_label  <- as.numeric(get_yaml_value(fig_style, c("typography","sizes_pt","axis_label_default"), 6.5))
pt_legend_text <- as.numeric(get_yaml_value(fig_style, c("typography","sizes_pt","legend_text_default"), 6))
pt_legend_title<- as.numeric(get_yaml_value(fig_style, c("typography","sizes_pt","legend_title_default"), 6.5))
pt_title       <- as.numeric(get_yaml_value(fig_style, c("typography","sizes_pt","title_optional"), 7))

pt_axis_line   <- as.numeric(get_yaml_value(fig_style, c("lines","axis_line_default_pt"), 0.25))
pt_grid_major  <- as.numeric(get_yaml_value(fig_style, c("lines","major_grid_default_pt"), 0.35))
pt_grid_minor  <- as.numeric(get_yaml_value(fig_style, c("lines","minor_grid_default_pt"), 0.25))

col_grid_major <- get_yaml_value(fig_style, c("colors","grid_major"), "grey85")
col_grid_minor <- get_yaml_value(fig_style, c("colors","grid_minor"), "grey92")

theme_nature_like <- function() {
  theme_bw(base_size = pt_axis_label, base_family = style_font_family) +
    theme(
      axis.text  = element_text(size = pt_axis_tick),
      axis.title = element_text(size = pt_axis_label),
      legend.text  = element_text(size = pt_legend_text),
      legend.title = element_text(size = pt_legend_title),
      plot.title   = element_text(size = pt_title, hjust = 0.5, face = "bold"),
      axis.line    = element_line(linewidth = pt_to_mm(pt_axis_line), colour = "black"),
      panel.grid.major = element_line(linewidth = pt_to_mm(pt_grid_major), colour = col_grid_major),
      panel.grid.minor = element_line(linewidth = pt_to_mm(pt_grid_minor), colour = col_grid_minor)
    )
}

## --------- 2. 计算每条通路跨 batch 的综合 Z -----------------

cat("\n[STEP2] 计算每条 EphB1 下游通路跨 batch 的综合 Z (z_meta) ...\n")

overall <- scores_use %>%
  group_by(pathway, axis) %>%
  summarise(
    n_batch       = n_distinct(batch),
    n_batch_nonNA = sum(!is.na(z_pathway)),
    mean_z_pathway = ifelse(
      n_batch_nonNA > 0,
      mean(z_pathway, na.rm = TRUE),
      NA_real_
    ),
    ## Stouffer 风格综合 Z（只看 z_pathway 本身）
    z_meta = ifelse(
      n_batch_nonNA > 0,
      sum(z_pathway, na.rm = TRUE) / sqrt(n_batch_nonNA),
      NA_real_
    ),
    mean_LFC_overall = ifelse(
      n_batch_nonNA > 0,
      mean(mean_LFC, na.rm = TRUE),
      NA_real_
    ),
    .groups = "drop"
  ) %>%
  arrange(desc(z_meta))

cat("  [overall] nrow = ", nrow(overall),
    " ; pathway 数 = ", length(unique(overall$pathway)), "\n", sep = "")

## 为可视化准备因子顺序（按 z_meta 从小到大，这样条形图从下到上强度递增）
overall <- overall %>%
  mutate(
    pathway_factor = forcats::fct_reorder(pathway, z_meta, .desc = FALSE)
  )

## 用于条形图：统一对称 x 轴范围（让 0 作为视觉中心，便于比较强弱）
max_abs_z_meta <- max(abs(overall$z_meta), na.rm = TRUE)
if (!is.finite(max_abs_z_meta) || max_abs_z_meta <= 0) max_abs_z_meta <- 1
# 取到 0.5 的倍数，视觉更干净
max_abs_z_meta <- ceiling(max_abs_z_meta * 2) / 2

## --------- 3. 画条形图：EphB1 下游通路整体效应 -------------

cat("\n[STEP3] 绘制 EphB1 下游通路整体综合 Z 条形图 ...\n")

dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

p_bar <- ggplot(overall,
                aes(y = pathway_factor, x = z_meta, fill = z_meta)) +
  geom_col(width = 0.72) +
  geom_vline(xintercept = 0, linewidth = pt_to_mm(pt_axis_line)) +
  scale_x_continuous(
    limits = c(-max_abs_z_meta, max_abs_z_meta),
    breaks = scales::pretty_breaks(n = 4)
  ) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B",
    midpoint = 0,
    name = "z_meta"
  ) +
  labs(
    x = "Cross-batch meta Z",
    y = NULL
  ) +
  theme_nature_like() +
  theme(
    # 右侧作为热图的数值注释：隐藏通路名称与 y 轴
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    # 去掉图例：让右侧更像注释条
    legend.position = "none",
    # 减少网格干扰
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(t = 2, r = 2, b = 2, l = 0.5, unit = "mm")
  )

out_bar <- file.path(plots_dir, "figure_7_c_ephb1_downstream_bar_overall.png")

# 单栏宽度 84 mm；高度可按内容调整
height_bar_mm <- combined_height_mm

ggsave(out_bar, p_bar,
       width = max_width_mm, height = height_bar_mm,
       units = "mm", dpi = 300)

out_bar_pdf <- file.path(plots_dir, "figure_7_c_ephb1_downstream_bar_overall.pdf")

ggsave(out_bar_pdf, p_bar,
       width = max_width_mm, height = height_bar_mm,
       units = "mm", device = cairo_pdf)

cat("  [OK] 已保存条形图: ", out_bar, "\n", sep = "")
cat("  [OK] 已保存条形图 PDF: ", out_bar_pdf, "\n", sep = "")

## --------- 4. 画热图：pathway × batch 的 z_pathway -----------

cat("\n[STEP4] 绘制 EphB1 下游通路径 × batch z_pathway 热图 ...\n")

## 为热图统一 pathway 顺序（沿用 overall 的顺序）
scores_heat <- scores_use %>%
  mutate(
    pathway = factor(pathway, levels = levels(overall$pathway_factor))
  )

p_heat <- ggplot(scores_heat,
                 aes(x = batch, y = pathway, fill = z_pathway)) +
  geom_tile(color = "white", linewidth = pt_to_mm(pt_grid_minor)) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B",
    midpoint = 0,
    name = "z_pathway"
  ) +
  labs(
    x     = "Batch",
    y     = NULL
  ) +
  theme_nature_like() +
  theme(
    axis.text.x = element_text(size = pt_axis_tick, angle = 0, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = pt_axis_tick),
    legend.position = "bottom",
    panel.grid  = element_blank(),
    plot.margin = margin(t = 2, r = 2, b = 2, l = 2, unit = "mm")
  )

out_heat <- file.path(plots_dir, "figure_7_c_ephb1_downstream_heatmap_per_batch.png")

height_heat_mm <- combined_height_mm

ggsave(out_heat, p_heat,
       width = max_width_mm, height = height_heat_mm,
       units = "mm", dpi = 300)

out_heat_pdf <- file.path(plots_dir, "figure_7_c_ephb1_downstream_heatmap_per_batch.pdf")

ggsave(out_heat_pdf, p_heat,
       width = max_width_mm, height = height_heat_mm,
       units = "mm", device = cairo_pdf)

cat("  [OK] 已保存热图: ", out_heat, "\n", sep = "")
cat("  [OK] 已保存热图 PDF: ", out_heat_pdf, "\n", sep = "")

## --------- 5. 导出组合图：左热图 + 右条形注释（同高左右排列） ---------

cat("\n[STEP5] 导出组合图（左热图 + 右条形注释，同高左右排列） ...\n")

# 组合：左为热图（主图），右为条形注释（z_meta）
# widths 控制左右占比：热图更宽，条形注释更窄
p_combo <- p_heat + p_bar +
  plot_layout(widths = c(3.2, 1))

out_combo_png <- file.path(plots_dir, "figure_7_c_ephb1_downstream_heatmap_plus_meta_bar.png")
out_combo_pdf <- file.path(plots_dir, "figure_7_c_ephb1_downstream_heatmap_plus_meta_bar.pdf")

ggsave(out_combo_png, p_combo,
       width = combined_width_mm, height = combined_height_mm,
       units = "mm", dpi = 300)

ggsave(out_combo_pdf, p_combo,
       width = combined_width_mm, height = combined_height_mm,
       units = "mm", device = cairo_pdf)

cat("  [OK] 已保存组合图: ", out_combo_png, "\n", sep = "")
cat("  [OK] 已保存组合图 PDF: ", out_combo_pdf, "\n", sep = "")

cat("============================================================\n")
cat("[DONE] 06_ephb1_downstream_plots.R 完成。\n")
cat("  - 请查看:\n")
cat("      ", out_bar,  "\n", sep = "")
cat("      ", out_bar_pdf,  "\n", sep = "")
cat("      ", out_heat, "\n", sep = "")
cat("      ", out_heat_pdf, "\n", sep = "")
cat("      ", out_combo_png, "\n", sep = "")
cat("      ", out_combo_pdf, "\n", sep = "")
cat("============================================================\n")