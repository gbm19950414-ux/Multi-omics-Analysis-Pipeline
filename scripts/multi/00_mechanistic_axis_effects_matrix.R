#!/usr/bin/env Rscript

## ============================================================
## 07_mechanistic_axis_effects_matrix.R
##
## 目的：
##   读取转录组 (RNA) 与脂质组 (Lipid) 的机制轴 **per-batch** 效应，
##   在每个 omic × batch × axis 上整理 KO vs WT 的效应 (effect_z)，
##   生成一个适合画 Panel 1D heatmap 的长表和宽表矩阵。
##
## 重要变化：
##   - RNA：仍然直接使用 rna_mechanistic_axis_scores_per_batch.tsv 中的 z_score。
##   - Lipid：不再从 per-sample 表重新做 t 检验，
##            而是直接读取 06a_lipid_axis_effect_ho_vs_wt.R 输出的
##            lipid_axis_effect_ho_vs_wt.tsv 中的 z_lipid_axis。
##   这样可以保证：
##     * 00 脚本中的脂质 per-batch Z 与 01 脚本中的完全一致；
##     * 去掉重复的 t 检验计算。
##
## 默认输入 (可通过 config 覆盖)：
##   - RNA per-batch 轴分数表：
##       results/rna/mechanistic_axis/rna_mechanistic_axis_scores_per_batch.tsv
##   - Lipid per-batch 轴效应表 (来自 06a)：
##       results/lipid/tables/lipid_axis_effect_ho_vs_wt.tsv
##
## 输出：
##   1) results/multiomics/mechanistic_axis_effects_long.tsv
##      每行对应一个 omic × batch × axis 的 KO vs WT 效应：
##        - effect   = mean_ko - mean_ref
##        - effect_z = 对应的 Z 统计量 (RNA: z_score, Lipid: z_lipid_axis)
##        - p_value, FDR 等
##   2) results/multiomics/mechanistic_axis_effects_matrix.tsv
##      行：axis
##      列：omics_batch（如 RNA_batch1, Lipid_batch1）
##      值：effect_z，作为 Panel 1D heatmap 输入。
##
## 用法：
##   Rscript scripts/multi/07_mechanistic_axis_effects_matrix.R \
##     scripts/multi/00_mechanistic_axis_effects_config.yaml
##
## 若未提供 config，则使用脚本内置默认路径。
##
## 作者：ChatGPT（根据用户需求定制，2025-11-25 版本：去除脂质重复 t 检验）
## ============================================================

suppressPackageStartupMessages({
  library(yaml)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
  library(purrr)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

## --------- 0. 解析命令行参数 & 配置 -------------------------

args <- commandArgs(trailingOnly = TRUE)

config_path <- if (length(args) >= 1) {
  args[1]
} else {
  "scripts/multi/00_mechanistic_axis_effects_config.yaml"
}

if (file.exists(config_path)) {
  cat("[INFO] 使用配置文件: ", config_path, "\n", sep = "")
  cfg <- yaml::read_yaml(config_path)
} else {
  cat("[WARN] 未找到配置文件: ", config_path,
      " ，将使用脚本内置默认路径。\n", sep = "")
  cfg <- list()
}

## RNA per-batch 轴分数表路径

RNA_tables_dir <- cfg$RNA_tables_dir %||% "results/rna/mechanistic_axis"
RNA_axis_file  <- cfg$RNA_axis_per_batch %||% "rna_mechanistic_axis_scores_per_batch.tsv"
RNA_axis_path  <- file.path(RNA_tables_dir, RNA_axis_file)

## Lipid per-batch 轴效应表路径（来自 06a_lipid_axis_effect_ho_vs_wt.R）

Lipid_tables_dir        <- cfg$Lipid_tables_dir %||% "results/lipid/tables"
Lipid_axis_effect_file  <- cfg$Lipid_axis_effect_file %||% "lipid_axis_effect_ho_vs_wt.tsv"
Lipid_axis_effect_path  <- file.path(Lipid_tables_dir, Lipid_axis_effect_file)

## 输出路径

output_dir      <- cfg$output_dir %||% "results/multiomics"
long_out_path   <- file.path(output_dir, "mechanistic_axis_effects_long.tsv")
matrix_out_path <- file.path(output_dir, "mechanistic_axis_effects_matrix.tsv")

cat("============================================================\n")
cat("[INFO] 07_mechanistic_axis_effects_matrix.R\n")
cat("  RNA per-batch axis:    ", RNA_axis_path,          "\n", sep = "")
cat("  Lipid axis effect (06a): ", Lipid_axis_effect_path, "\n", sep = "")
cat("  output long  : ", long_out_path,   "\n", sep = "")
cat("  output matrix: ", matrix_out_path, "\n", sep = "")
cat("============================================================\n\n")

## --------- 1. 读取 RNA / Lipid per-batch 轴表 -----------

if (!file.exists(RNA_axis_path)) {
  stop("[ERROR] 找不到 RNA per-batch 轴分数表: ", RNA_axis_path)
}
if (!file.exists(Lipid_axis_effect_path)) {
  stop("[ERROR] 找不到 Lipid per-batch 轴效应表: ", Lipid_axis_effect_path)
}

cat("[STEP1] 读取 RNA / Lipid per-batch 轴表 ...\n")

RNA_df          <- readr::read_tsv(RNA_axis_path,          show_col_types = FALSE)
lipid_effect_df <- readr::read_tsv(Lipid_axis_effect_path, show_col_types = FALSE)

cat("  [RNA_df]          nrow = ", nrow(RNA_df),          " ; ncol = ", ncol(RNA_df),          "\n", sep = "")
cat("  [lipid_effect_df] nrow = ", nrow(lipid_effect_df), " ; ncol = ", ncol(lipid_effect_df), "\n", sep = "")

## 简单检查必要列是否存在

required_RNA_cols <- c("batch", "axis", "mean_LFC", "z_score")
miss_RNA <- setdiff(required_RNA_cols, colnames(RNA_df))
if (length(miss_RNA) > 0) {
  stop("[ERROR] RNA per-batch 表缺少必要列: ", paste(miss_RNA, collapse = ", "))
}

required_Lipid_cols <- c(
  "batch", "axis",
  "ref_group", "ko_group",
  "n_ref", "n_ko",
  "mean_ref", "mean_ko",
  "delta_axis", "z_lipid_axis", "p_value"
)
miss_Lipid <- setdiff(required_Lipid_cols, colnames(lipid_effect_df))
if (length(miss_Lipid) > 0) {
  stop("[ERROR] Lipid per-batch 轴效应表缺少必要列: ", paste(miss_Lipid, collapse = ", "))
}

## --------- 2. RNA / Lipid per-batch 轴效应 → 统一格式 long 表 ---------

cat("[STEP2] 整理 RNA / Lipid 轴效应为统一 long 表 ...\n")

## RNA：直接使用 rna_mechanistic_axis_scores_per_batch.tsv 中的 mean_LFC / z_score

axis_effects_RNA <- RNA_df %>%
  transmute(
    omic      = "RNA",
    batch     = batch,
    axis      = axis,
    group_ko  = "HO",          # 约定：z_score 定义为 HO vs WT
    group_ref = "WT",
    n_ko      = NA_integer_,
    n_ref     = NA_integer_,
    mean_ko   = NA_real_,
    mean_ref  = NA_real_,
    effect    = mean_LFC,       # log2( HO / WT )
    t_stat    = z_score,        # RNA Z 统计量
    df_t      = NA_real_,
    p_value   = NA_real_,
    FDR       = NA_real_,
    effect_z  = z_score
  )

## Lipid：直接使用 06a 计算好的 delta_axis / z_lipid_axis / p_value

axis_effects_Lipid <- lipid_effect_df %>%
  transmute(
    omic      = "Lipid",
    batch     = batch,
    axis      = axis,
    group_ko  = ko_group,
    group_ref = ref_group,
    n_ko      = n_ko,
    n_ref     = n_ref,
    mean_ko   = mean_ko,
    mean_ref  = mean_ref,
    effect    = delta_axis,      # mean_ko - mean_ref
    t_stat    = z_lipid_axis,    # 06a 中定义的 Z 统计量
    df_t      = NA_real_,
    p_value   = p_value,
    FDR       = ifelse(is.na(p_value), NA_real_, p.adjust(p_value, method = "BH")),
    effect_z  = z_lipid_axis
  )

axis_effects <- bind_rows(axis_effects_RNA, axis_effects_Lipid)

cat("  [axis_effects] nrow = ", nrow(axis_effects),
    " ; ncol = ", ncol(axis_effects), "\n", sep = "")

## --------- 3. 输出 long 表 & 构建 Panel 1D 矩阵 --------------

cat("[STEP3] 写出 long 表，并构建 axis × (omics_batch) 矩阵 ...\n")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

readr::write_tsv(axis_effects, long_out_path)

## 构建矩阵：行 = axis，列 = omic_batch，值 = effect_z
axis_effects_wide <- axis_effects %>%
  mutate(omics_batch = paste(omic, batch, sep = "_")) %>%
  select(axis, omics_batch, effect_z) %>%
  tidyr::pivot_wider(
    names_from  = omics_batch,
    values_from = effect_z
  )

readr::write_tsv(axis_effects_wide, matrix_out_path)

cat("  [OK] 写出 long 表:   ", long_out_path,   "\n", sep = "")
cat("  [OK] 写出矩阵表: ", matrix_out_path, "\n", sep = "")

cat("============================================================\n")
cat("[DONE] 07_mechanistic_axis_effects_matrix.R 完成。\n")
cat("  下一步：可使用该矩阵作为 Panel 1D heatmap 的输入 (例如 figure_6_d.R)。\n")
cat("============================================================\n")