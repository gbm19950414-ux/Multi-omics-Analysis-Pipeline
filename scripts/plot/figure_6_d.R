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
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

## --------- 0. 解析命令行参数 & 配置 -------------------------

args <- commandArgs(trailingOnly = TRUE)

config_path <- if (length(args) >= 1) {
  args[1]
} else {
  "scripts/plot/figure_6_d_config.yaml"
}

if (file.exists(config_path)) {
  cat("[INFO] 使用配置文件: ", config_path, "\n", sep = "")
  cfg <- yaml::read_yaml(config_path)
} else {
  cat("[WARN] 未找到配置文件: ", config_path,
      " ，将使用脚本内置默认设置。\n", sep = "")
  cfg <- list()
}

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
  geom_tile(color = "grey90") +
  ## 显著性用小星号标注
  geom_text(
    aes(label = sig_label),
    size = 2.8,
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
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.text.y  = element_text(color = "black"),
    panel.grid   = element_blank(),
    legend.position = "right",
    plot.margin  = margin(5.5, 5.5, 5.5, 5.5)
  )

out_path <- file.path(output_dir, output_file)

ggsave(out_path, p, width = 6, height = 4, device = "pdf")

cat("  [OK] 输出图像: ", out_path, "\n", sep = "")
cat("============================================================\n")
cat("[DONE] figure_6_d.R 完成，可用于 Panel 1D heatmap。\n")
cat("============================================================\n")