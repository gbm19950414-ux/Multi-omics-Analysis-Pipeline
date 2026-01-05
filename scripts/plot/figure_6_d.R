#!/usr/bin/env Rscript

## ============================================================
## figure_6_d.R
##
## Panel 1D – 机制轴 × 批次 × 组学 heatmap
##
## 读入：
##   - results/multiomics/mechanistic_axis_effects_long.tsv
##     （由 07_mechanistic_axis_effects_matrix.R 生成）
##
## 画图：
##   - 行：机制轴 axis
##   - 列：omics_batch（如 RNA_batch1, RNA_batch2, Lipid_batch1, ...）
##   - 填色：effect_z（蓝负/红正）
##   - 覆盖显著性标记：
##       FDR < 0.05 或（FDR 为 NA 且 |effect_z| >= 1.96）
##
## 用法：
##   Rscript scripts/plot/figure_6_d.R \
##     scripts/plot/figure_6_d_config.yaml
##
## 若未提供 config，则使用默认路径。
## ============================================================

suppressPackageStartupMessages({
  library(yaml)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(purrr)
  library(grid)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

## --------- 0. 解析命令行参数 & 配置 -------------------------

args <- commandArgs(trailingOnly = TRUE)

config_path <- if (length(args) >= 1) {
  args[1]
} else {
  "/Volumes/Samsung_SSD_990_PRO_2TB_Media/EphB1/02_protocols/figure_style_nature.yaml"
}

if (file.exists(config_path)) {
  cat("[INFO] 使用配置文件: ", config_path, "\n", sep = "")
  cfg <- yaml::read_yaml(config_path)
} else {
  cat("[WARN] 未找到配置文件: ", config_path,
      " ，将使用脚本内置默认设置。\n", sep = "")
  cfg <- list()
}

# --------- Global figure style (Nature) ---------------------
style_path <- "/Volumes/Samsung_SSD_990_PRO_2TB_Media/EphB1/02_protocols/figure_style_nature.yaml"
if (!file.exists(style_path)) {
  stop("[ERROR] 找不到 figure style 文件: ", style_path)
}
style <- yaml::read_yaml(style_path)

## 兼容两种 YAML 写法：
# A) typography$sizes_pt$legend_text_default (推荐：本项目 figure_style_nature.yaml)
# B) legend_text_default$family / legend_text_default$size_pt (旧版)
font_family <- style$typography$font_family_primary %||%
  style$legend_text_default$family %||% "Helvetica"

font_size_pt <- style$typography$sizes_pt$legend_text_default %||%
  style$legend_text_default$size_pt %||% 10

if (is.null(font_size_pt) || is.na(font_size_pt)) {
  stop("[ERROR] figure_style_nature.yaml 中未找到 typography:sizes_pt:legend_text_default（或旧版 legend_text_default:size_pt）")
}

# 进一步从 figure_style_nature.yaml 读取外观参数（字号/线宽/边距）
axis_tick_pt  <- style$typography$sizes_pt$axis_tick_default  %||% font_size_pt
axis_label_pt <- style$typography$sizes_pt$axis_label_default %||% font_size_pt
legend_text_pt  <- style$typography$sizes_pt$legend_text_default  %||% font_size_pt
legend_title_pt <- style$typography$sizes_pt$legend_title_default %||% font_size_pt

line_width_pt <- style$lines$line_width_pt %||% 0.5
axis_line_pt  <- style$lines$axis_line_default_pt %||% 0.25

plot_margin_pt <- style$layout$plot_margin_pt %||% list(top = 0, right = 0, bottom = 0, left = 0)

# pt → mm（ggplot 的 linewidth 用 mm）
pt_to_mm <- function(pt) pt * 25.4 / 72
line_width_mm <- pt_to_mm(line_width_pt)
axis_line_mm  <- pt_to_mm(axis_line_pt)

effects_long_path <- cfg$effects_long_path %||%
  "results/multiomics/mechanistic_axis_effects_long.tsv"

output_dir  <- cfg$output_dir  %||% "results/figs"
output_file <- cfg$output_file %||% "figure_6_d_mechanistic_axis_heatmap.pdf"

z_thresh   <- cfg$z_threshold   %||% 1.96
fdr_thresh <- cfg$fdr_threshold %||% 0.05

## 默认轴顺序（可在 config 中通过 axis_order 覆盖）
default_axis_order <- c(
  "Synthesis",
  "Remodeling",
  "Oxidation",
  "Transport",
  "Supply",
  "Membrane context"
)
axis_order <- cfg$axis_order %||% default_axis_order

cat("============================================================\n")
cat("[INFO] figure_6_d.R\n")
cat("  input  : ", effects_long_path, "\n", sep = "")
cat("  output : ", file.path(output_dir, output_file), "\n", sep = "")
cat("============================================================\n\n")

## --------- 1. 读取 long 表 ----------------------------------

if (!file.exists(effects_long_path)) {
  stop("[ERROR] 找不到机制轴效应 long 表: ", effects_long_path)
}

cat("[STEP1] 读取机制轴效应 long 表 ...\n")

dat <- readr::read_tsv(effects_long_path, show_col_types = FALSE)

cat("  [dat] nrow = ", nrow(dat), " ; ncol = ", ncol(dat), "\n", sep = "")

needed_cols <- c("omic", "batch", "axis", "effect_z", "FDR")
miss_cols <- setdiff(needed_cols, colnames(dat))
if (length(miss_cols) > 0) {
  stop("[ERROR] 输入表缺少必要列: ", paste(miss_cols, collapse = ", "))
}

## --------- 2. 构建 omics_batch & 因子顺序 -------------------

cat("[STEP2] 构建 omics_batch 列，并设置行列顺序 ...\n")

dat2 <- dat %>%
  mutate(
    omics_batch = paste(omic, batch, sep = "_")
  )

## 列顺序：默认 RNA_* 在前，Lipid_* 在后，同组内按 batch 排序
omics_batch_levels <- dat2 %>%
  distinct(omic, batch, omics_batch) %>%
  mutate(
    omic = factor(omic, levels = c("RNA", "Lipid"), ordered = TRUE)
  ) %>%
  arrange(omic, batch) %>%
  pull(omics_batch)

## 若配置中指定 column_order，则覆盖
if (!is.null(cfg$column_order)) {
  omics_batch_levels <- cfg$column_order
}

## 行顺序：按 axis_order 否则按出现顺序
axis_levels <- axis_order
missing_axes <- setdiff(unique(dat2$axis), axis_levels)
if (length(missing_axes) > 0) {
  axis_levels <- c(axis_levels, missing_axes)
}

dat2 <- dat2 %>%
  mutate(
    axis        = factor(axis, levels = axis_levels),
    omics_batch = factor(omics_batch, levels = omics_batch_levels)
  )

## --------- 3. 计算显著性标记 -------------------------------

cat("[STEP3] 计算显著性标记 (FDR / |Z|) ...\n")

dat2 <- dat2 %>%
  mutate(
    is_sig_fdr = !is.na(FDR) & FDR < fdr_thresh,
    is_sig_z   = is.na(FDR) & !is.na(effect_z) & abs(effect_z) >= z_thresh,
    is_sig     = is_sig_fdr | is_sig_z,
    sig_label  = dplyr::if_else(is_sig, "*", "")
  )

## --------- 4. 绘制 heatmap ----------------------------------

cat("[STEP4] 绘制 Panel 1D heatmap ...\n")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

p <- ggplot(dat2, aes(x = omics_batch, y = axis, fill = effect_z)) +
  geom_tile(color = "grey90", linewidth = axis_line_mm) +
  ## 显著性用小星号标注
  geom_text(
    aes(label = sig_label),
    family = font_family,
    size = legend_text_pt / ggplot2::.pt,
    color = "black"
  ) +
  scale_fill_gradient2(
    name = "Effect Z",
    low = "#2166AC", mid = "white", high = "#B2182B",
    midpoint = 0, na.value = "grey90"
  ) +
  coord_fixed() +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = font_size_pt, base_family = font_family) +
  theme(
    text = element_text(family = font_family, size = legend_text_pt),
    axis.text.x  = element_text(angle = 45, hjust = 1, vjust = 1, family = font_family, size = axis_tick_pt),
    axis.text.y  = element_text(color = "black", family = font_family, size = axis_tick_pt),
    axis.title.x = element_text(family = font_family, size = axis_label_pt),
    axis.title.y = element_text(family = font_family, size = axis_label_pt),
    legend.text  = element_text(family = font_family, size = legend_text_pt),
    legend.title = element_text(family = font_family, size = legend_title_pt),
    panel.grid   = element_blank(),
    legend.position = "right",
    plot.margin  = margin(plot_margin_pt$top, plot_margin_pt$right, plot_margin_pt$bottom, plot_margin_pt$left, unit = "pt")
  )

out_path <- file.path(output_dir, output_file)

# 输出尺寸：宽度严格 84 mm；高度按 tile 数量比例自适应（并留出一定空间给坐标文字）
width_mm <- 84
n_x <- nlevels(dat2$omics_batch)
n_y <- nlevels(dat2$axis)
height_mm <- width_mm * (n_y / max(1, n_x)) + 10  # +10mm 给旋转 x 轴文字/边距

ggsave(out_path, p, width = width_mm, height = height_mm, units = "mm", device = grDevices::cairo_pdf)

cat("  [OK] 输出图像: ", out_path, "\n", sep = "")
cat("============================================================\n")
cat("[DONE] figure_6_d.R 完成，可用于 Panel 1D heatmap。\n")
cat("============================================================\n")