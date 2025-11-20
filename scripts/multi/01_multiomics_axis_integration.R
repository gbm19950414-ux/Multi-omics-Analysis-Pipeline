#!/usr/bin/env Rscript

## ============================================================
## 01_multiomics_axis_integration.R
##
## 目的：
##   在机制轴（axis）× 批次（batch）层面整合
##   转录组与脂质组的机制 Z 分数，得到：
##     1）per-axis × batch 的 multi-omics Z
##     2）per-axis 的整体 multi-omics Z（跨 batch meta）
##
## 输入（默认路径，可在 config 中覆盖）：
##   - results/rna/mechanistic_axis/rna_mechanistic_axis_scores_per_batch.tsv
##       必须包含列：axis, batch, z_score
##   - results/lipid/tables/lipid_axis_effect_ho_vs_wt.tsv
##       （由 06a_lipid_axis_effect_ho_vs_wt.R 生成）
##       必须包含列：axis, batch, z_lipid_axis
##
## 可选配置文件（YAML）：
##   默认：scripts/multiomics/00_multiomics_axis_config.yaml
##
##   支持字段：
##     rna_axis_scores_path:   RNA 轴 Z 输入表路径
##     lipid_axis_effect_path: 脂质轴 Z 输入表路径
##     tables_dir:             综合结果输出目录（默认：results/multiomics/tables）
##     plots_dir:              绘图输出目录（默认：results/multiomics/plots）
##     axis_levels:            轴名顺序（字符向量），例如：
##                             ["Synthesis", "Remodeling", "Oxidation", "Transport", "Supply"]
##
## 输出：
##   1) results/multiomics/tables/multiomics_axis_Z_per_batch.tsv
##        axis, batch,
##        z_rna, p_rna,
##        z_lipid, p_lipid,
##        n_omics, Z_combined, p_combined
##
##   2) results/multiomics/tables/multiomics_axis_Z_overall.tsv
##        axis,
##        n_batch_rna, Z_rna_overall, p_rna_overall,
##        n_batch_lipid, Z_lipid_overall, p_lipid_overall,
##        n_batch_multiomics, Z_multiomics_overall, p_multiomics_overall
##
##   3) results/multiomics/plots/multiomics_axis_Z_overall_barplot.pdf
##        各轴整体 multi-omics Z 的柱状图
##
## 整合方法：
##   - axis × batch 层面：
##       若同时有 z_rna 和 z_lipid：
##         Z_combined = (z_rna + z_lipid) / sqrt(2)
##       若只存在一个组学：
##         Z_combined = 该组学的 Z
##
##   - axis 层面（跨 batch meta）：
##       对每个 axis，分别对各来源 Z 做 Stouffer 等权合并：
##         Z_overall = sum(Z_i) / sqrt(k)   (k = 非 NA 的 Z 数量)
##
## ============================================================

suppressPackageStartupMessages({
  library(yaml)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

## --------- 0. 解析命令行参数 & 读取配置 ----------------------

args <- commandArgs(trailingOnly = TRUE)

config_path <- if (length(args) >= 1) {
  args[1]
} else {
  "scripts/multiomics/00_multiomics_axis_config.yaml"
}

if (file.exists(config_path)) {
  cat("[INFO] 使用配置文件: ", config_path, "\n", sep = "")
  cfg <- yaml::read_yaml(config_path)
} else {
  cat("[WARN] 未找到配置文件: ", config_path,
      "，将使用脚本内置默认路径。\n", sep = "")
  cfg <- list()
}

tables_dir_default <- "results/multiomics/tables"
plots_dir_default  <- "results/multiomics/plots"

tables_dir <- cfg$tables_dir %||% tables_dir_default
plots_dir  <- cfg$plots_dir  %||% plots_dir_default

rna_axis_scores_path <- cfg$rna_axis_scores_path %||%
  "results/rna/mechanistic_axis/rna_mechanistic_axis_scores_per_batch.tsv"

lipid_axis_effect_path <- cfg$lipid_axis_effect_path %||%
  "results/lipid/tables/lipid_axis_effect_ho_vs_wt.tsv"

axis_levels <- cfg$axis_levels %||% NULL  # 可选：控制绘图顺序

cat("============================================================\n")
cat("[INFO] 01_multiomics_axis_integration.R\n")
cat("  RNA axis scores  : ", rna_axis_scores_path, "\n", sep = "")
cat("  Lipid axis effect: ", lipid_axis_effect_path, "\n", sep = "")
cat("  tables_dir       : ", tables_dir, "\n", sep = "")
cat("  plots_dir        : ", plots_dir, "\n", sep = "")
cat("============================================================\n\n")

## --------- 1. 读取 RNA & lipid 轴 Z 表 -----------------------

if (!file.exists(rna_axis_scores_path)) {
  stop("[ERROR] 找不到 RNA 轴分数表: ", rna_axis_scores_path)
}
if (!file.exists(lipid_axis_effect_path)) {
  stop("[ERROR] 找不到脂质轴效应表: ", lipid_axis_effect_path)
}

cat("[STEP1] 读取 RNA 轴 Z 表 ...\n")
rna_axis <- readr::read_tsv(rna_axis_scores_path, show_col_types = FALSE)
cat("  [rna_axis] nrow = ", nrow(rna_axis),
    " ; ncol = ", ncol(rna_axis), "\n", sep = "")

required_rna_cols <- c("axis", "batch", "z_score")
missing_rna <- setdiff(required_rna_cols, colnames(rna_axis))
if (length(missing_rna) > 0) {
  stop("[ERROR] RNA 轴表缺少必要列: ", paste(missing_rna, collapse = ", "))
}

## 只保留需要的列，并重命名 z_score -> z_rna
rna_axis_use <- rna_axis %>%
  transmute(
    axis,
    batch,
    z_rna = z_score
  )

cat("[STEP1] 读取脂质轴 Z 表 ...\n")
lipid_axis <- readr::read_tsv(lipid_axis_effect_path, show_col_types = FALSE)
cat("  [lipid_axis] nrow = ", nrow(lipid_axis),
    " ; ncol = ", ncol(lipid_axis), "\n", sep = "")

required_lipid_cols <- c("axis", "batch", "z_lipid_axis")
missing_lipid <- setdiff(required_lipid_cols, colnames(lipid_axis))
if (length(missing_lipid) > 0) {
  stop("[ERROR] 脂质轴表缺少必要列: ", paste(missing_lipid, collapse = ", "))
}

lipid_axis_use <- lipid_axis %>%
  transmute(
    axis,
    batch,
    z_lipid = z_lipid_axis
  )

## 如有 axis_levels 配置，则统一 axis 因子顺序
if (!is.null(axis_levels)) {
  rna_axis_use <- rna_axis_use %>%
    mutate(axis = factor(axis, levels = axis_levels))
  lipid_axis_use <- lipid_axis_use %>%
    mutate(axis = factor(axis, levels = axis_levels))
}

## --------- 2. axis × batch 层面的 Z 合并 ---------------------

cat("\n[STEP2] 在 axis × batch 层面合并 RNA 与脂质 Z ...\n")

axis_batch <- full_join(
  rna_axis_use,
  lipid_axis_use,
  by = c("axis", "batch")
)

cat("  [axis_batch] nrow (union of RNA & lipid) = ",
    nrow(axis_batch), "\n", sep = "")

axis_batch <- axis_batch %>%
  mutate(
    ## 为了方便查看，也给出各自的 p 值（双侧）
    p_rna   = if_else(!is.na(z_rna),
                      2 * pnorm(-abs(z_rna)), NA_real_),
    p_lipid = if_else(!is.na(z_lipid),
                      2 * pnorm(-abs(z_lipid)), NA_real_),

    n_omics = as.integer(!is.na(z_rna)) + as.integer(!is.na(z_lipid)),

    Z_combined = dplyr::case_when(
      n_omics == 2L ~ (z_rna + z_lipid) / sqrt(2),
      n_omics == 1L ~ dplyr::coalesce(z_rna, z_lipid),
      TRUE          ~ NA_real_
    ),

    p_combined = if_else(
      !is.na(Z_combined),
      2 * pnorm(-abs(Z_combined)),
      NA_real_
    )
  )

## 输出 axis × batch 层面的 multi-omics Z
dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)

per_batch_out_path <- file.path(tables_dir, "multiomics_axis_Z_per_batch.tsv")
readr::write_tsv(axis_batch, per_batch_out_path)

cat("  [OK] 写出 axis × batch multi-omics Z 表: ",
    per_batch_out_path, "\n", sep = "")

## --------- 3. axis 层面（跨 batch meta）Z 合并 ---------------

cat("\n[STEP3] 在 axis 层面做跨 batch meta（Stouffer 等权）...\n")

## Stouffer 等权合并：Z_overall = sum(Z_i) / sqrt(k)
stouffer_equal <- function(z_vec) {
  z_vec <- z_vec[is.finite(z_vec)]
  k <- length(z_vec)
  if (k == 0) return(NA_real_)
  sum(z_vec) / sqrt(k)
}

axis_overall <- axis_batch %>%
  group_by(axis) %>%
  summarise(
    ## RNA-only meta
    n_batch_rna    = sum(!is.na(z_rna)),
    Z_rna_overall  = stouffer_equal(z_rna),
    p_rna_overall  = if_else(
      !is.na(Z_rna_overall),
      2 * pnorm(-abs(Z_rna_overall)),
      NA_real_
    ),

    ## Lipid-only meta
    n_batch_lipid       = sum(!is.na(z_lipid)),
    Z_lipid_overall     = stouffer_equal(z_lipid),
    p_lipid_overall     = if_else(
      !is.na(Z_lipid_overall),
      2 * pnorm(-abs(Z_lipid_overall)),
      NA_real_
    ),

    ## Multi-omics meta（基于 Z_combined）
    n_batch_multiomics  = sum(!is.na(Z_combined)),
    Z_multiomics_overall = stouffer_equal(Z_combined),
    p_multiomics_overall = if_else(
      !is.na(Z_multiomics_overall),
      2 * pnorm(-abs(Z_multiomics_overall)),
      NA_real_
    ),
    .groups = "drop"
  )

overall_out_path <- file.path(tables_dir, "multiomics_axis_Z_overall.tsv")
readr::write_tsv(axis_overall, overall_out_path)

cat("  [OK] 写出 axis 层面 multi-omics meta Z 表: ",
    overall_out_path, "\n", sep = "")

## --------- 4. 绘制整体 multi-omics Z 柱状图 ------------------

cat("\n[STEP4] 绘制各机制轴整体 multi-omics Z 柱状图 ...\n")

plot_df <- axis_overall %>%
  mutate(
    axis = if (!is.null(axis_levels)) {
      factor(axis, levels = axis_levels)
    } else {
      axis
    }
  )

p_bar <- ggplot(plot_df, aes(x = axis, y = Z_multiomics_overall)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_col(width = 0.6) +
  geom_text(
    aes(
      label = sprintf("%.2f", Z_multiomics_overall),
      vjust = ifelse(Z_multiomics_overall >= 0, -0.3, 1.1)
    ),
    size = 3
  ) +
  labs(
    title = "Multi-omics integrated axis Z-score (HO vs WT)",
    x = "Mechanistic axis",
    y = "Z_multiomics_overall"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)
plot_path <- file.path(plots_dir, "multiomics_axis_Z_overall_barplot.pdf")
ggsave(plot_path, p_bar, width = 6, height = 4)

cat("  [OK] 写出柱状图: ", plot_path, "\n", sep = "")

cat("\n============================================================\n")
cat("[DONE] 01_multiomics_axis_integration.R 完成。\n")
cat("  说明：\n")
cat("    - Z_multiomics_overall > 0 表示在该轴上综合证据支持 HO 瓶颈增强；\n")
cat("    - Z_multiomics_overall < 0 表示综合证据支持 HO 瓶颈减轻或补偿；\n")
cat("    - 绝对值越大，效应越强，p_multiomics_overall 越小。\n")
cat("============================================================\n")