#!/usr/bin/env Rscript

## ============================================================
## 05a_mechanistic_indicators_per_batch.R
##
## 目的：
##   在「单批次 PQN 矫正」后的长表基础上，
##   为每个批次、每个生物学样本计算脂质学机制指标：
##     - 按 Class 汇总 ΣCL, ΣPG, ΣMLCL, ΣoxCL, ΣPA, ΣPC, ΣPE
##     - 计算机制指标（CL/PG、MLCL/CL、oxCL/CL、PG/PA、ΣCL 等）
##     - 根据方向性矩阵（sign matrix）做 z-score + 符号校正，
##       得到每个机制轴（合成、重塑、氧化、转运、供给）的瓶颈分数
##
## 关键特征：
##   - 完全不依赖 batch 间 scaling（不使用 intensity_pqn_bb）
##   - 每个 batch 单独标准化（z-score 在 batch 内计算）
##   - 输出：
##       1）per sample 指标 + 轴分数
##       2）per batch × group × axis 的汇总表
##
## 依赖：
##   - lipid_long_dedup_pqn.tsv
##   - sample_info_lipid.tsv
##   - lipid_indicator_sign_matrix.yaml
##
## 作者：ChatGPT（根据用户需求定制）
## ============================================================

suppressPackageStartupMessages({
  library(yaml)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
  library(purrr)
})

## --------- 0. 解析命令行参数 & 配置 -------------------------

args <- commandArgs(trailingOnly = TRUE)

## 可以可选传入一个 config（目前主要用于保持路径一致）
## 若未提供，则使用默认路径
config_path <- if (length(args) >= 1) {
  args[1]
} else {
  "scripts/lipid/00_lipid_downstream_config.yaml"
}

tables_dir_default <- "results/lipid/tables"
qc_dir_default     <- "results/lipid/qc"
plots_dir_default  <- "results/lipid/plots"

## 读取 00_config（若存在），否则使用默认路径
if (file.exists(config_path)) {
  cat("[INFO] 使用配置文件: ", config_path, "\n", sep = "")
  cfg <- yaml::read_yaml(config_path)
} else {
  cat("[WARN] 未找到配置文件: ", config_path,
      "，将使用脚本内置默认路径。\n", sep = "")
  cfg <- list()
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

tables_dir <- cfg$tables_dir %||% tables_dir_default
qc_dir     <- cfg$qc_dir     %||% qc_dir_default
plots_dir  <- cfg$plots_dir  %||% plots_dir_default

## 机制指标 & 符号矩阵配置路径（如未在 cfg 中定义，则使用默认）
indicator_cfg_path <- cfg$lipid_mechanistic_indicator_config %||%
  "scripts/lipid/lipid_mechanistic_indicators.yaml"

sign_matrix_path   <- cfg$lipid_indicator_sign_matrix %||%
  "scripts/lipid/lipid_indicator_sign_matrix.yaml"

## 输入文件路径（本脚本只用 PQN，不用 PQN+bb）
long_pqn_path <- file.path(tables_dir, "lipid_long_dedup_pqn.tsv")
sample_info_path <- file.path(tables_dir, "sample_info_lipid.tsv")

## 输出文件路径
out_per_sample_path <- file.path(
  tables_dir, "lipid_mechanistic_indicators_per_sample.tsv"
)
out_per_batch_group_path <- file.path(
  tables_dir, "lipid_mechanistic_axis_scores_per_batch_group.tsv"
)

cat("============================================================\n")
cat("[INFO] 10_lipid_mechanistic_indicators_per_batch.R\n")
cat("  input long : ", long_pqn_path, "\n", sep = "")
cat("  input sample_info: ", sample_info_path, "\n", sep = "")
cat("  sign matrix : ", sign_matrix_path, "\n", sep = "")
cat("  output per-sample  : ", out_per_sample_path, "\n", sep = "")
cat("  output per-batch×group×axis: ", out_per_batch_group_path, "\n", sep = "")
cat("============================================================\n\n")

## --------- 1. 读取数据 & 检查 ------------------------------

if (!file.exists(long_pqn_path)) {
  stop("[ERROR] 找不到 lipid_long_dedup_pqn.tsv: ", long_pqn_path)
}

if (!file.exists(sample_info_path)) {
  stop("[ERROR] 找不到 sample_info_lipid.tsv: ", sample_info_path)
}

if (!file.exists(sign_matrix_path)) {
  stop("[ERROR] 找不到 sign matrix 配置文件: ", sign_matrix_path)
}

cat("[STEP1] 读取 PQN 后长表 & sample_info ...\n")

long_pqn <- readr::read_tsv(long_pqn_path, guess_max = 1e6, show_col_types = FALSE)
sample_info <- readr::read_tsv(sample_info_path, show_col_types = FALSE)

cat("  [long_pqn] nrow = ", nrow(long_pqn),
    " ; ncol = ", ncol(long_pqn), "\n", sep = "")
cat("  [sample_info] nrow = ", nrow(sample_info),
    " ; ncol = ", ncol(sample_info), "\n", sep = "")

required_cols_long <- c("batch", "sample_id", "Class", "intensity_pqn")
missing_long <- setdiff(required_cols_long, colnames(long_pqn))
if (length(missing_long) > 0) {
  stop("[ERROR] lipid_long_dedup_pqn.tsv 缺少必要列: ",
       paste(missing_long, collapse = ", "))
}

required_cols_si <- c("batch", "sample_id", "group")
missing_si <- setdiff(required_cols_si, colnames(sample_info))
if (length(missing_si) > 0) {
  stop("[ERROR] sample_info_lipid.tsv 缺少必要列: ",
       paste(missing_si, collapse = ", "))
}

## 若存在 is_qc 列，则转换为逻辑，并尽量排除 QC 样本
to_logical <- function(x) {
  if (is.logical(x)) return(x)
  if (is.numeric(x)) return(x != 0)
  if (is.character(x)) {
    return(toupper(x) %in% c("TRUE", "T", "1", "YES", "Y"))
  }
  as.logical(x)
}

if ("is_qc" %in% colnames(long_pqn)) {
  long_pqn <- long_pqn %>% mutate(is_qc = to_logical(is_qc))
} else {
  long_pqn$is_qc <- FALSE
  cat("  [INFO] 长表中没有 is_qc 列，将默认所有记录为非 QC。\n")
}

## 将 group 信息合并到长表
long_anno <- long_pqn %>%
  left_join(
    sample_info %>% select(batch, sample_id, group),
    by = c("batch", "sample_id")
  )

## 只保留生物学样本（group 非 NA）
long_bio <- long_anno %>%
  filter(!is.na(group))

cat("  [INFO] 生物学样本记录数（过滤 group 非 NA 之后）: ",
    nrow(long_bio), "\n", sep = "")

## --------- 2. 读取 sign matrix（方向性信息） -----------------

cat("\n[STEP2] 读取 directionality sign matrix ...\n")

sign_cfg <- yaml::read_yaml(sign_matrix_path)

## sign_cfg 的结构：
##   Axis:
##     indicator_name:
##       bottleneck_sign: +1/-1
## 需要把它展开为 data.frame: axis, indicator, bottleneck_sign

sign_df <- purrr::imap_dfr(sign_cfg, function(axis_list, axis_name) {
  purrr::imap_dfr(axis_list, function(ind_cfg, ind_name) {
    tibble::tibble(
      axis = axis_name,
      indicator = ind_name,
      bottleneck_sign = as.numeric(ind_cfg$bottleneck_sign)
    )
  })
})

cat("  [sign_df] 行数 = ", nrow(sign_df),
    " ; 涉及 axis 数 = ", length(unique(sign_df$axis)),
    " ; 指标数 = ", length(unique(sign_df$indicator)), "\n", sep = "")

## --------- 3. 按 Class 统计 ΣCL/ΣPG/ΣMLCL/...（per sample） ----

cat("\n[STEP3] 按 Class / sample 统计总强度 ...\n")

## 这里假设 long_bio$Class 已经规范为:
##   "CL", "PG", "MLCL", "oxCL", "PA", "PC", "PE" 等
## 如果实际有别名（如 "OxCL"、"Cardiolipin" 等），可在此处做简单映射。

class_map <- list(
  CL   = c("CL"),
  PG   = c("PG"),
  MLCL = c("MLCL"),
  oxCL = c("oxCL", "OxCL"),
  PA   = c("PA"),
  PC   = c("PC"),
  PE   = c("PE")
)

## 把 Class 映射为简化的 class_key
long_bio_mapped <- long_bio %>%
  mutate(
    class_key = dplyr::case_when(
      Class %in% class_map$CL   ~ "CL",
      Class %in% class_map$PG   ~ "PG",
      Class %in% class_map$MLCL ~ "MLCL",
      Class %in% class_map$oxCL ~ "oxCL",
      Class %in% class_map$PA   ~ "PA",
      Class %in% class_map$PC   ~ "PC",
      Class %in% class_map$PE   ~ "PE",
      TRUE                      ~ NA_character_
    )
  )

## 只对我们关心的这些 class_key 做聚合
long_bio_use <- long_bio_mapped %>%
  filter(!is.na(class_key), !is.na(intensity_pqn), intensity_pqn > 0)

class_totals <- long_bio_use %>%
  group_by(batch, sample_id, group, class_key) %>%
  summarise(
    sum_intensity = sum(intensity_pqn, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  tidyr::pivot_wider(
    names_from = class_key,
    values_from = sum_intensity,
    values_fill = 0
  )

cat("  [class_totals] nrow = ", nrow(class_totals),
    " ; ncol = ", ncol(class_totals), "\n", sep = "")

## 确保 Σ 列存在（即使全 0）
safe_add_col <- function(df, col) {
  if (!col %in% colnames(df)) df[[col]] <- 0
  df
}

for (nm in c("CL", "PG", "MLCL", "oxCL", "PA", "PC", "PE")) {
  class_totals <- safe_add_col(class_totals, nm)
}

## --------- 4. 计算脂质机制指标（raw 指标） -------------------

cat("\n[STEP4] 计算脂质机制指标（raw 指标）...\n")

indicators_df <- class_totals %>%
  mutate(
    ## 合成轴
    CL_PG_ratio = dplyr::if_else(PG > 0, CL / PG, NA_real_),
    sum_CL      = CL,
    sum_PG      = PG,

    ## 重塑轴
    MLCL_CL_ratio = dplyr::if_else(CL > 0, MLCL / CL, NA_real_),
    sum_MLCL      = MLCL,

    ## 氧化轴
    oxCL_CL_ratio = dplyr::if_else(CL > 0, oxCL / CL, NA_real_),
    sum_oxCL      = oxCL,

    ## 转运 / 供给轴
    PG_PA_ratio = dplyr::if_else(PA > 0, PG / PA, NA_real_),

    mito_lipid_mass = PC + PE + CL,
    CL_fraction     = dplyr::if_else(
      (PC + PE + CL) > 0,
      CL / (PC + PE + CL),
      NA_real_
    )
  )

## 当前脚本实现的指标名（可以与 sign_matrix.yaml 对应）
available_indicators <- setdiff(
  colnames(indicators_df),
  c("batch", "sample_id", "group")
)

cat("  [INFO] 当前可计算的指标列:\n    ",
    paste(available_indicators, collapse = ", "), "\n", sep = "")

## --------- 5. 在 batch 内对指标做 z-score + 符号校正 -----------

cat("\n[STEP5] 在每个 batch 内进行 z-score & 瓶颈方向校正...\n")

## sign matrix 中出现的指标
sign_indicators <- unique(sign_df$indicator)

## 取交集：既在 sign matrix 中又在实际数据中存在的指标
use_indicators <- intersect(sign_indicators, available_indicators)

if (length(use_indicators) == 0) {
  stop("[ERROR] sign matrix 中的指标在当前数据中都不存在，请检查命名是否一致。")
}

cat("  [INFO] 将用于瓶颈评分的指标:\n    ",
    paste(use_indicators, collapse = ", "), "\n", sep = "")

## 把 wide 格式转为 long 格式：batch × sample × indicator
ind_long <- indicators_df %>%
  select(batch, sample_id, group, all_of(use_indicators)) %>%
  pivot_longer(
    cols = all_of(use_indicators),
    names_to = "indicator",
    values_to = "value"
  )

## 合并 sign 信息
ind_long_sign <- ind_long %>%
  left_join(sign_df, by = "indicator")

## 在 batch 内对每个 indicator 做 z-score
ind_long_z <- ind_long_sign %>%
  group_by(batch, indicator) %>%
  mutate(
    mean_val = mean(value, na.rm = TRUE),
    sd_val   = sd(value, na.rm = TRUE),
    z_raw    = dplyr::if_else(
      is.finite(sd_val) & sd_val > 0,
      (value - mean_val) / sd_val,
      NA_real_
    )
  ) %>%
  ungroup() %>%
  mutate(
    signed_z = z_raw * bottleneck_sign
  )

## --------- 6. 聚合到机制轴（axis-level bottleneck score） -------

cat("\n[STEP6] 计算每个机制轴（axis）的瓶颈分数 ...\n")

axis_scores <- ind_long_z %>%
  group_by(batch, sample_id, group, axis) %>%
  summarise(
    n_indicators = sum(!is.na(signed_z)),
    axis_score   = dplyr::if_else(
      n_indicators > 0,
      mean(signed_z, na.rm = TRUE),
      NA_real_
    ),
    .groups = "drop"
  )

## wide 化：每个 axis 一列，例如 Synthesis_score, Remodeling_score, ...
axis_scores_wide <- axis_scores %>%
  mutate(axis_col = paste0(axis, "_score")) %>%
  select(batch, sample_id, group, axis_col, axis_score) %>%
  pivot_wider(
    names_from = axis_col,
    values_from = axis_score
  )

## --------- 7. 合并 raw 指标 + 轴分数，并输出 per-sample 表 -----

cat("\n[STEP7] 合并 raw 指标与 axis 分数，输出 per-sample 指标表 ...\n")

per_sample_out <- indicators_df %>%
  left_join(axis_scores_wide,
            by = c("batch", "sample_id", "group"))

dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)

readr::write_tsv(per_sample_out, out_per_sample_path)

cat("  [OK] 写出 per-sample 指标表: ",
    out_per_sample_path, "\n", sep = "")

## --------- 8. 按 batch × group × axis 做汇总（均值/标准差） -----

cat("\n[STEP8] 按 batch × group × axis 聚合瓶颈分数 ...\n")

per_batch_group_axis <- axis_scores %>%
  group_by(batch, group, axis) %>%
  summarise(
    n_sample   = n(),
    mean_score = mean(axis_score, na.rm = TRUE),
    sd_score   = sd(axis_score, na.rm = TRUE),
    .groups = "drop"
  )

readr::write_tsv(per_batch_group_axis, out_per_batch_group_path)

cat("  [OK] 写出 per-batch × group × axis 汇总表: ",
    out_per_batch_group_path, "\n", sep = "")

cat("============================================================\n")
cat("[DONE] 10_lipid_mechanistic_indicators_per_batch.R 完成。\n")
cat("  下一步：\n")
cat("    - 可以在转录组机制 gene set 上构建类似的 axis score\n")
cat("    - 然后在多批次层面对 axis score 做 meta 分析 / 关联分析。\n")
cat("============================================================\n")