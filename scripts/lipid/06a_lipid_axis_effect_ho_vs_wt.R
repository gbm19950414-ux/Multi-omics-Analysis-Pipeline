#!/usr/bin/env Rscript

## ============================================================
## 06a_lipid_axis_effect_ho_vs_wt.R
##
## 目的：
##   从现有的 per-batch × group × axis 表
##   （lipid_mechanistic_axis_scores_per_batch_group.tsv）
##   推导出每个 batch、每条 axis 的：
##     - HO vs WT 轴级效应量 (delta_axis)
##     - 差值的标准误 (SE_delta)
##     - 轴级 Z 分数 (z_lipid_axis)
##     - 对应双侧 p 值 (p_value)
##
## 输入：
##   - 00_lipid_downstream_config.yaml （可选）
##       * tables_dir: 结果表路径（默认：results/lipid/tables）
##       * ref_group:  参考组名（默认："WT"）
##       * ko_group:   敲除组名（默认："HO"）
##       * lipid_axis_scores_per_batch_group:  (可选) 输入表路径覆盖
##       * lipid_axis_effect_ho_vs_wt:         (可选) 输出表路径覆盖
##
##   - results/lipid/tables/lipid_mechanistic_axis_scores_per_batch_group.tsv
##       （由 05a_lipid_mechanistic_indicators_per_batch.R 生成）
##       必须至少包含列：
##         batch, group, axis, n_sample, mean_score, sd_score
##
## 输出：
##   - results/lipid/tables/lipid_axis_effect_ho_vs_wt.tsv
##       列包括：
##         batch, axis,
##         ref_group, ko_group,
##         n_ref, n_ko,
##         mean_ref, mean_ko,
##         sd_ref, sd_ko,
##         delta_axis, se_delta, z_lipid_axis, p_value
##
## 解释：
##   delta_axis   = mean_score(ko_group) - mean_score(ref_group)
##   se_ref       = sd_ref / sqrt(n_ref)
##   se_ko        = sd_ko  / sqrt(n_ko)
##   se_delta     = sqrt(se_ref^2 + se_ko^2)
##   z_lipid_axis = delta_axis / se_delta
##   p_value      = 2 * pnorm(-abs(z_lipid_axis))
##
##   其中 mean_score 来自 05a 生成的 axis_score（样本瓶颈 z 的平均），
##   因此 z_lipid_axis 表示：在“样本间噪音”的单位下，
##   HO vs WT 在该机制轴上的瓶颈效应强度。
##
## ============================================================

suppressPackageStartupMessages({
  library(yaml)
  library(dplyr)
  library(readr)
  library(tidyr)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

## --------- 0. 解析命令行参数 & 配置 -------------------------

args <- commandArgs(trailingOnly = TRUE)

config_path <- if (length(args) >= 1) {
  args[1]
} else {
  "scripts/lipid/00_lipid_downstream_config.yaml"
}

if (file.exists(config_path)) {
  cat("[INFO] 使用配置文件: ", config_path, "\n", sep = "")
  cfg <- yaml::read_yaml(config_path)
} else {
  cat("[WARN] 未找到配置文件: ", config_path,
      "，将使用脚本内置默认路径。\n", sep = "")
  cfg <- list()
}

tables_dir_default <- "results/lipid/tables"
tables_dir <- cfg$tables_dir %||% tables_dir_default

## 输入输出路径
in_path <- cfg$lipid_axis_scores_per_batch_group %||%
  file.path(tables_dir, "lipid_mechanistic_axis_scores_per_batch_group.tsv")

out_path <- cfg$lipid_axis_effect_ho_vs_wt %||%
  file.path(tables_dir, "lipid_axis_effect_ho_vs_wt.tsv")

## 参考组 / 敲除组名称
ref_group <- cfg$ref_group %||% "WT"
ko_group  <- cfg$ko_group  %||% "HO"

cat("============================================================\n")
cat("[INFO] 06a_lipid_axis_effect_ho_vs_wt.R\n")
cat("  input  : ", in_path, "\n", sep = "")
cat("  output : ", out_path, "\n", sep = "")
cat("  ref_group = ", ref_group, " ; ko_group = ", ko_group, "\n", sep = "")
cat("============================================================\n\n")

## --------- 1. 读取输入表 & 检查 -----------------------------

if (!file.exists(in_path)) {
  stop("[ERROR] 找不到输入表 lipid_mechanistic_axis_scores_per_batch_group.tsv: ", in_path)
}

axis_grp <- readr::read_tsv(in_path, show_col_types = FALSE)

cat("[STEP1] 读取 per-batch × group × axis 表 ...\n")
cat("  nrow = ", nrow(axis_grp), " ; ncol = ", ncol(axis_grp), "\n", sep = "")

required_cols <- c("batch", "group", "axis", "n_sample", "mean_score", "sd_score")
missing_cols <- setdiff(required_cols, colnames(axis_grp))
if (length(missing_cols) > 0) {
  stop("[ERROR] 输入表缺少必要列: ", paste(missing_cols, collapse = ", "))
}

## 检查 ref_group / ko_group 是否存在
groups_present <- unique(axis_grp$group)
if (!ref_group %in% groups_present) {
  warning("[WARN] ref_group = ", ref_group, " 在数据中不存在。实际 group 值: ",
          paste(groups_present, collapse = ", "))
}
if (!ko_group %in% groups_present) {
  warning("[WARN] ko_group = ", ko_group, " 在数据中不存在。实际 group 值: ",
          paste(groups_present, collapse = ", "))
}

axis_grp_sub <- axis_grp %>%
  filter(group %in% c(ref_group, ko_group))

cat("  [INFO] 过滤只保留 ref_group / ko_group 后行数: ",
    nrow(axis_grp_sub), "\n", sep = "")

if (nrow(axis_grp_sub) == 0) {
  stop("[ERROR] 过滤 ref_group / ko_group 后没有任何数据，请检查 group 命名。")
}

## --------- 2. 宽表化：一行 = batch × axis，包含两组信息 ---------

cat("\n[STEP2] 宽表化 batch × axis × (ref, ko) ...\n")

axis_wide <- axis_grp_sub %>%
  select(batch, axis, group, n_sample, mean_score, sd_score) %>%
  tidyr::pivot_wider(
    names_from  = group,
    values_from = c(n_sample, mean_score, sd_score),
    names_sep   = "_"
  )

cat("  [axis_wide] nrow = ", nrow(axis_wide),
    " ; ncol = ", ncol(axis_wide), "\n", sep = "")

## 构造列名字符串，方便后面通过 .data[[col]] 访问
col_n_ref   <- paste0("n_sample_",   ref_group)
col_n_ko    <- paste0("n_sample_",   ko_group)
col_mean_ref<- paste0("mean_score_", ref_group)
col_mean_ko <- paste0("mean_score_", ko_group)
col_sd_ref  <- paste0("sd_score_",   ref_group)
col_sd_ko   <- paste0("sd_score_",   ko_group)

## --------- 3. 计算 delta / SE_delta / Z / p ------------------

cat("\n[STEP3] 计算 HO vs WT 的轴级效应量 / 标准误 / Z 分数 ...\n")

if (!all(c(col_n_ref, col_n_ko,
           col_mean_ref, col_mean_ko,
           col_sd_ref, col_sd_ko) %in% colnames(axis_wide))) {
  warning("[WARN] 某些 batch × axis 缺少 ref 或 ko 组信息，将在相应行给出 NA。")
}

axis_effect <- axis_wide %>%
  mutate(
    ## 各组样本数
    n_ref = .data[[col_n_ref]],
    n_ko  = .data[[col_n_ko]],

    ## 组内均值 / SD
    mean_ref = .data[[col_mean_ref]],
    mean_ko  = .data[[col_mean_ko]],
    sd_ref   = .data[[col_sd_ref]],
    sd_ko    = .data[[col_sd_ko]],

    ## 组内标准误
    se_ref = if_else(!is.na(sd_ref) & !is.na(n_ref) & n_ref > 0,
                     sd_ref / sqrt(n_ref),
                     NA_real_),
    se_ko  = if_else(!is.na(sd_ko) & !is.na(n_ko) & n_ko > 0,
                     sd_ko / sqrt(n_ko),
                     NA_real_),

    ## 轴级差异（ko - ref）
    delta_axis = if_else(!is.na(mean_ref) & !is.na(mean_ko),
                         mean_ko - mean_ref,
                         NA_real_),

    ## 差异的标准误（假设两组估计独立）
    se_delta = if_else(
      !is.na(se_ref) & !is.na(se_ko),
      sqrt(se_ref^2 + se_ko^2),
      NA_real_
    ),

    z_lipid_axis = if_else(
      is.finite(se_delta) & se_delta > 0,
      delta_axis / se_delta,
      NA_real_
    ),

    p_value = if_else(
      !is.na(z_lipid_axis),
      2 * pnorm(-abs(z_lipid_axis)),
      NA_real_
    ),

    ref_group = ref_group,
    ko_group  = ko_group
  ) %>%
  select(
    batch, axis,
    ref_group, ko_group,
    n_ref, n_ko,
    mean_ref, mean_ko,
    sd_ref, sd_ko,
    delta_axis, se_delta, z_lipid_axis, p_value
  )

cat("  [INFO] 计算完毕：共有 ",
    nrow(axis_effect), " 条 batch × axis 记录。\n", sep = "")

## --------- 4. 写出结果 ---------------------------------------

dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)
readr::write_tsv(axis_effect, out_path)

cat("\n[OK] 写出 HO vs WT 轴级效应表: ", out_path, "\n", sep = "")
cat("============================================================\n")
cat("[DONE] 06a_lipid_axis_effect_ho_vs_wt.R 完成。\n")
cat("  说明：\n")
cat("    - delta_axis > 0 表示在该轴上 HO 的 axis_score 高于 WT，\n")
cat("      若 axis_score 已按瓶颈方向校正，则表示 HO 瓶颈更严重。\n")
cat("    - z_lipid_axis 越大（绝对值越大），表示效应越强、越显著。\n")
cat("============================================================\n")