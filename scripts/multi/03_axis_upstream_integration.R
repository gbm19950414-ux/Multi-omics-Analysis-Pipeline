#!/usr/bin/env Rscript

## ============================================================
## 03_axis_upstream_integration.R
##
## 目的：
##   将机制轴 multi-omics Z（脂质+RNA）与 RNA 上游通路 Z 整合，
##   在 axis × batch 层面对齐，并在 axis × pathway 层面计算
##   “上游支持度”与相关性，输出可用于挑选关键上游通路的表。
##
## 输入（默认，可通过 config 覆盖）：
##
##   1) 机制轴 multi-omics Z（per batch）：
##        results/multiomics/tables/multiomics_axis_Z_per_batch.tsv
##      要求至少包含列：
##        axis, batch, 以及一个“综合 Z”列，列名优先匹配：
##          - Z_multiomics_combined
##          - Z_multiomics
##          - Z_axis_multiomics
##          - Z_combined
##          - Z_meta
##
##   2) 上游通路 Z（per batch）：
##        results/multiomics/tables/rna_upstream_pathway_scores_per_batch.tsv
##      要求至少包含列：
##        pathway, axis, batch, z_pathway
##
##   3) 可选配置文件（命令行第 1 个参数）：
##        scripts/multiomics/00_axis_upstream_integration_config.yaml
##
##      支持字段：
##        multiomics_axis_Z_per_batch: "results/multiomics/tables/multiomics_axis_Z_per_batch.tsv"
##        rna_upstream_pathway_Z_per_batch: "results/multiomics/tables/rna_upstream_pathway_scores_per_batch.tsv"
##        tables_dir: "results/multiomics/tables"
##
## 输出：
##
##   1) axis × pathway × batch 的合并表：
##        results/multiomics/tables/axis_upstream_Z_per_batch.tsv
##      列：
##        axis, batch, pathway,
##        Z_axis, z_pathway
##
##   2) axis × pathway 跨 batch 汇总表：
##        results/multiomics/tables/axis_upstream_support_overall.tsv
##      列（示例）：
##        axis, pathway,
##        n_batch, n_batch_nonNA,
##        mean_Z_axis, mean_z_pathway,
##        cor_Z_axis_pathway,
##        support_Z   （Stouffer 风格的“方向一致支持度”）
##
##   解释：
##     - Z_axis：多组学层面该机制轴在该 batch 的效应（正=机制增强，负=机制受抑）
##     - z_pathway：该上游通路在该 batch 的转录层面效应（正=通路激活，负=通路受抑）
##     - support_Z：
##         对每个 axis–pathway，跨 batch 计算：
##           signed_z_i = sign(Z_axis_i) * z_pathway_i
##         然后：
##           support_Z = sum(signed_z_i) / sqrt(n_batch_nonNA)
##         若：
##           - 轴是负向（被抑制），通路也负向（被抑制） -> support_Z > 0
##           - 轴是正向（增强），通路也正向（激活）   -> support_Z > 0
##           - 方向相反或不一致                   -> support_Z ~ 0 或 < 0
##
##     support_Z 越大、相关性越高的通路，
##     越可能是该机制轴变化的上游驱动候选。
##
## ============================================================

suppressPackageStartupMessages({
  library(yaml)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(purrr)
  library(stringr)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

## --------- 0. 解析命令行参数 & 配置 -------------------------

args <- commandArgs(trailingOnly = TRUE)

config_path <- if (length(args) >= 1) {
  args[1]
} else {
  "scripts/multiomics/00_axis_upstream_integration_config.yaml"
}

if (file.exists(config_path)) {
  cat("[INFO] 使用配置文件: ", config_path, "\n", sep = "")
  cfg <- yaml::read_yaml(config_path)
} else {
  cat("[WARN] 未找到配置文件: ", config_path,
      "，将使用脚本内置默认路径。\n", sep = "")
  cfg <- list()
}

axis_Z_per_batch_default   <- "results/multiomics/tables/multiomics_axis_Z_per_batch.tsv"
upstream_Z_per_batch_default <- "results/multiomics/tables/rna_upstream_pathway_scores_per_batch.tsv"
tables_dir_default          <- "results/multiomics/tables"

axis_Z_path   <- cfg$multiomics_axis_Z_per_batch        %||% axis_Z_per_batch_default
upstream_path <- cfg$rna_upstream_pathway_Z_per_batch   %||% upstream_Z_per_batch_default
tables_dir    <- cfg$tables_dir                         %||% tables_dir_default

cat("============================================================\n")
cat("[INFO] 03_axis_upstream_integration.R\n")
cat("  axis Z per-batch      : ", axis_Z_path,   "\n", sep = "")
cat("  upstream Z per-batch  : ", upstream_path, "\n", sep = "")
cat("  tables_dir            : ", tables_dir,    "\n", sep = "")
cat("============================================================\n\n")

## --------- 1. 读取 multi-omics axis Z（per batch） ----------

if (!file.exists(axis_Z_path)) {
  stop("[ERROR] 找不到 multi-omics 轴 Z 表: ", axis_Z_path)
}

cat("[STEP1] 读取 multi-omics 轴 Z 表 (per batch) ...\n")
axis_Z <- readr::read_tsv(axis_Z_path, show_col_types = FALSE)

cat("  [axis_Z] nrow = ", nrow(axis_Z),
    " ; ncol = ", ncol(axis_Z), "\n", sep = "")

if (!all(c("axis", "batch") %in% colnames(axis_Z))) {
  stop("[ERROR] axis_Z 表中缺少 axis 或 batch 列。实际列: ",
       paste(colnames(axis_Z), collapse = ", "))
}

## 选一个“综合 Z”列
candidate_Z_cols <- c(
  "Z_multiomics_combined",
  "Z_multiomics",
  "Z_axis_multiomics",
  "Z_combined",
  "Z_meta"
)

Z_col <- intersect(candidate_Z_cols, colnames(axis_Z))

if (length(Z_col) == 0) {
  stop("[ERROR] 在 multi-omics 轴 Z 表中找不到综合 Z 列。\n",
       "请确认是否存在以下列名之一: ",
       paste(candidate_Z_cols, collapse = ", "), "\n",
       "当前列名: ", paste(colnames(axis_Z), collapse = ", "))
} else if (length(Z_col) > 1) {
  cat("  [WARN] 检测到多个候选综合 Z 列，将使用第一个: ",
      Z_col[1], "\n", sep = "")
}

Z_col <- Z_col[1]

axis_Z_use <- axis_Z %>%
  dplyr::select(axis, batch, Z_axis = !!sym(Z_col))

cat("  [INFO] 使用综合 Z 列: ", Z_col, " 作为 Z_axis\n", sep = "")

## --------- 2. 读取 upstream pathway Z（per batch） ----------

if (!file.exists(upstream_path)) {
  stop("[ERROR] 找不到 RNA 上游通路 Z 表: ", upstream_path)
}

cat("\n[STEP2] 读取 RNA 上游通路 Z 表 (per batch) ...\n")
upstream_Z <- readr::read_tsv(upstream_path, show_col_types = FALSE)

cat("  [upstream_Z] nrow = ", nrow(upstream_Z),
    " ; ncol = ", ncol(upstream_Z), "\n", sep = "")

required_cols_upstream <- c("pathway", "axis", "batch", "z_pathway")

if (!all(required_cols_upstream %in% colnames(upstream_Z))) {
  stop("[ERROR] upstream Z 表中缺少以下必要列: ",
       paste(setdiff(required_cols_upstream, colnames(upstream_Z)),
             collapse = ", "),
       "\n当前列名: ", paste(colnames(upstream_Z), collapse = ", "))
}

## --------- 3. 在 axis × batch 层面合并 ----------------------

cat("\n[STEP3] 合并 axis Z 与 upstream pathway Z (axis × batch) ...\n")

## 只保留需要的列，防止重复列名干扰
upstream_for_join <- upstream_Z %>%
  dplyr::select(axis, batch, pathway, z_pathway)

axis_upstream_per_batch <- upstream_for_join %>%
  dplyr::inner_join(axis_Z_use, by = c("axis", "batch")) %>%
  dplyr::relocate(axis, batch, pathway, Z_axis, z_pathway)

cat("  [axis_upstream_per_batch] nrow = ", nrow(axis_upstream_per_batch),
    " ; axis 数 = ", length(unique(axis_upstream_per_batch$axis)),
    " ; pathway 数 = ", length(unique(axis_upstream_per_batch$pathway)),
    " ; batch 数 = ", length(unique(axis_upstream_per_batch$batch)), "\n", sep = "")

## --------- 4. 在 axis × pathway 层面汇总 --------------------

cat("\n[STEP4] 在 axis × pathway 层面计算跨 batch 总结指标 ...\n")

axis_upstream_overall <- axis_upstream_per_batch %>%
  group_by(axis, pathway) %>%
  summarise(
    n_batch        = n_distinct(batch),
    n_batch_nonNA  = sum(!is.na(Z_axis) & !is.na(z_pathway)),
    mean_Z_axis    = ifelse(n_batch_nonNA > 0,
                            mean(Z_axis, na.rm = TRUE),
                            NA_real_),
    mean_z_pathway = ifelse(n_batch_nonNA > 0,
                            mean(z_pathway, na.rm = TRUE),
                            NA_real_),
    cor_Z_axis_pathway = ifelse(
      n_batch_nonNA >= 2,
      suppressWarnings(stats::cor(Z_axis, z_pathway,
                                  use = "complete.obs")),
      NA_real_
    ),
    ## Stouffer 风格的“方向一致支持度”
    support_Z = ifelse(
      n_batch_nonNA > 0,
      {
        signed_z <- sign(Z_axis) * z_pathway
        sum(signed_z, na.rm = TRUE) / sqrt(n_batch_nonNA)
      },
      NA_real_
    ),
    .groups = "drop"
  ) %>%
  arrange(axis, desc(support_Z))

cat("  [axis_upstream_overall] nrow = ", nrow(axis_upstream_overall),
    " ; axis 数 = ", length(unique(axis_upstream_overall$axis)),
    " ; pathway 数 = ", length(unique(axis_upstream_overall$pathway)), "\n", sep = "")

## --------- 5. 写出结果 --------------------------------------

cat("\n[STEP5] 写出结果表 ...\n")

dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)

out_per_batch  <- file.path(tables_dir, "axis_upstream_Z_per_batch.tsv")
out_overall    <- file.path(tables_dir, "axis_upstream_support_overall.tsv")

readr::write_tsv(axis_upstream_per_batch, out_per_batch)
readr::write_tsv(axis_upstream_overall,   out_overall)

cat("  [OK] axis × pathway × batch 合并表  : ", out_per_batch, "\n", sep = "")
cat("  [OK] axis × pathway 跨 batch 汇总表 : ", out_overall,   "\n", sep = "")

cat("============================================================\n")
cat("[DONE] 03_axis_upstream_integration.R 完成。\n")
cat("  使用说明：\n")
cat("    - 对每个 axis，按 support_Z 从大到小查看最可能的上游通路；\n")
cat("    - support_Z > 0 且 |support_Z| 较大，说明通路变化方向与轴效应方向一致，\n")
cat("      且跨 batch 较为稳定，可作为“上游驱动候选”。\n")
cat("    - cor_Z_axis_pathway 可辅助判断 Z 序列是否同步（仅在 batch 数 ≥ 2 时有意义）。\n")
cat("============================================================\n")