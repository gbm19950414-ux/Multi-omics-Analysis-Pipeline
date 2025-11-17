#!/usr/bin/env Rscript

## ============================================================
## 07_lipid_DE_allfat_and_CL.R
##
## 用 limma 对 allfat / CL 两个矩阵做差异分析（默认 HO vs WT）
## - 使用 06_prepare_lipid_matrices.R 生成的：
##     * sample_info_lipid.tsv（含 global_sample_id / batch / group）
##     * lipid_matrix_allfat_log2_filtered.tsv
##     * lipid_matrix_CL_log2_filtered.tsv
## - 表达量已经是 log2(intensity + 1)
##
## 输出：
##   results/lipid/de/de_allfat_HO_vs_WT.tsv
##   results/lipid/de/de_CL_HO_vs_WT.tsv
##
## 逻辑：
##   1. 从矩阵中解析样本列（global_sample_id）
##   2. 用 sample_info_lipid.tsv 解析 group（WT / HO）
##   3. 仅保留 group 为 WT / HO 的生物学样本，去掉 QC / NA group
##   4. 使用 limma: log2Expr ~ group（WT 为基准，系数 groupHO）
##   5. 输出完整 DE 表（logFC, pval, FDR 等）并在终端汇总 up/down 数量
## ============================================================

suppressPackageStartupMessages({
  library(yaml)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(limma)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

log_msg <- function(...){
  cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), ..., "\n")
}

## ------------------------------------------------------------
## 0. 读配置
## ------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1) args[1] else "scripts/lipid/00_lipid_downstream_config.yaml"

if (!file.exists(config_path)) {
  stop("Config file not found: ", config_path,
       "\nUsage: Rscript 07_lipid_DE_allfat_and_CL.R path/to/00_lipid_downstream_config.yaml")
}

log_msg("[INFO] 使用配置文件: ", config_path)
cfg   <- yaml::read_yaml(config_path)
paths <- cfg$paths

de_dir <- paths$de_dir %||% "results/lipid/de"
dir.create(de_dir, recursive = TRUE, showWarnings = FALSE)

## ------------------------------------------------------------
## 1. 小工具：做一次 HO vs WT 的 DE 分析
## ------------------------------------------------------------
run_de_for_matrix <- function(matrix_path,
                              sample_info_path,
                              universe_label = "allfat",
                              out_prefix = "de_allfat_HO_vs_WT") {

  log_msg("------------------------------------------------------------")
  log_msg("[INFO] DE 分析开始: ", universe_label)
  log_msg("[INFO] 读取矩阵: ", matrix_path)

  mat <- readr::read_tsv(matrix_path, show_col_types = FALSE)

  ## 1.1 解析 feature 注释和样本列（global_sample_id）
  annot_cols <- c("feature_id", "lipidName", "LipidIon", "Class", "CalcMz")
  sample_cols <- setdiff(colnames(mat), annot_cols)

  if (length(sample_cols) == 0) {
    stop("在矩阵中没有检测到样本列（除注释列外）。")
  }

  log_msg("[INFO] 检测到样本列数: ", length(sample_cols))

  ## 1.2 读取 sample_info（含 global_sample_id）
  log_msg("[INFO] 读取 sample_info: ", sample_info_path)
  sample_info <- readr::read_tsv(sample_info_path, show_col_types = FALSE)

  if (!all(c("global_sample_id", "group") %in% colnames(sample_info))) {
    stop("sample_info_lipid.tsv 中必须包含列: global_sample_id, group")
  }

  ## 1.3 只保留 WT / HO 的生物学样本（排除 QC / 其他）
  meta <- sample_info %>%
    filter(global_sample_id %in% sample_cols) %>%
    mutate(group = as.character(group))

  log_msg("[INFO] 样本总数（矩阵中）: ", length(sample_cols))
  log_msg("[INFO] 样本在 sample_info 中匹配到的数量: ", nrow(meta))

  meta_de <- meta %>%
    filter(!is.na(group), group %in% c("WT", "HO"))

  log_msg("[INFO] 用于 DE 的生物学样本数 (WT/HO): ", nrow(meta_de))
  table_group <- table(meta_de$group)
  log_msg("[INFO] 组别分布: ",
          paste(names(table_group), table_group, collapse = " | "))

  if (length(unique(meta_de$group)) < 2) {
    stop("用于 DE 的样本中，group 少于 2 个水平（需要 WT 和 HO）。")
  }

  ## 1.4 构建表达矩阵（行 = feature, 列 = global_sample_id）
  expr_mat <- as.matrix(mat[, sample_cols])
  rownames(expr_mat) <- mat$feature_id

  ## 按 meta_de$global_sample_id 重新排序列
  expr_de <- expr_mat[, meta_de$global_sample_id, drop = FALSE]

  ## 1.5 设计矩阵（WT 为基准）
  group_factor <- factor(meta_de$group)
  if (!("WT" %in% levels(group_factor))) {
    stop("group 中没有 WT，无法设置 WT 为对照。")
  }
  group_factor <- relevel(group_factor, ref = "WT")

  design <- model.matrix(~ group_factor)
  colnames(design) <- c("Intercept", "groupHO")

  log_msg("[INFO] 设计矩阵列: ", paste(colnames(design), collapse = ", "))

  ## 1.6 limma 拟合
  log_msg("[INFO] 使用 limma 拟合线性模型 ...")
  fit <- limma::lmFit(expr_de, design)
  fit <- limma::eBayes(fit)

  ## 1.7 导出 DE 结果（对系数 groupHO）
  log_msg("[INFO] 提取 DE 结果（系数: groupHO, HO 相对 WT 的 logFC）...")
  tt <- limma::topTable(
    fit,
    coef = "groupHO",
    number = Inf,
    sort.by = "none"
  )

  ## 确保行顺序和 feature 对上
  tt$feature_id <- rownames(tt)

  annot <- mat %>%
    select(all_of(annot_cols))

  de_res <- annot %>%
    right_join(tt, by = "feature_id") %>%
    relocate(feature_id, lipidName, LipidIon, Class, CalcMz)

  ## 添加方便的 FDR 和一些阈值判断列
  de_res <- de_res %>%
    rename(
      logFC   = logFC,
      AveExpr = AveExpr,
      t       = t,
      P.Value = P.Value,
      adj.P.Val = adj.P.Val
    ) %>%
    mutate(
      sig_FDR_0.05   = adj.P.Val < 0.05,
      sig_FDR_0.05_logFC1 = adj.P.Val < 0.05 & abs(logFC) >= 1
    )

  ## 1.8 简单汇总
  n_sig <- sum(de_res$sig_FDR_0.05, na.rm = TRUE)
  n_sig_lfc1 <- sum(de_res$sig_FDR_0.05_logFC1, na.rm = TRUE)

  log_msg("[INFO] universe = ", universe_label,
          " ; FDR < 0.05 的 feature 数: ", n_sig)
  log_msg("[INFO] universe = ", universe_label,
          " ; FDR < 0.05 且 |logFC| >= 1 的 feature 数: ", n_sig_lfc1)

  ## 上下调各多少
  de_sig <- de_res %>%
    filter(sig_FDR_0.05)

  if (nrow(de_sig) > 0) {
    up   <- sum(de_sig$logFC > 0, na.rm = TRUE)
    down <- sum(de_sig$logFC < 0, na.rm = TRUE)
    log_msg("[INFO] universe = ", universe_label,
            " ; FDR < 0.05: 上调(HO>WT) = ", up,
            " ; 下调(HO<WT) = ", down)
  }

  ## 1.9 写出结果
  out_tsv <- file.path(de_dir, paste0(out_prefix, ".tsv"))
  readr::write_tsv(de_res, out_tsv)
  log_msg("[OK] 写出 DE 结果: ", out_tsv)

  invisible(de_res)
}

## ------------------------------------------------------------
## 2. 对 allfat / CL 分别做 HO vs WT 差异分析
## ------------------------------------------------------------
log_msg("============================================================")
log_msg("[INFO] 07_lipid_DE_allfat_and_CL.R 开始执行")
log_msg("============================================================")

## 路径从 config 中读取（与 06 脚本一致）
matrix_allfat <- paths$matrix_allfat_filtered %||% paths$matrix_allfat
matrix_CL     <- paths$matrix_CL_filtered     %||% paths$matrix_CL
sample_info   <- paths$sample_info_output     %||% "results/lipid/tables/sample_info_lipid.tsv"

## allfat
if (!is.null(matrix_allfat) && file.exists(matrix_allfat)) {
  run_de_for_matrix(
    matrix_path       = matrix_allfat,
    sample_info_path  = sample_info,
    universe_label    = "allfat",
    out_prefix        = "de_allfat_HO_vs_WT"
  )
} else {
  log_msg("[WARN] 找不到 allfat 矩阵文件: ", matrix_allfat,
          " ，跳过 allfat DE。")
}

## CL
if (!is.null(matrix_CL) && file.exists(matrix_CL)) {
  run_de_for_matrix(
    matrix_path       = matrix_CL,
    sample_info_path  = sample_info,
    universe_label    = "CL",
    out_prefix        = "de_CL_HO_vs_WT"
  )
} else {
  log_msg("[WARN] 找不到 CL 矩阵文件: ", matrix_CL,
          " ，跳过 CL DE。")
}

log_msg("============================================================")
log_msg("[DONE] 07_lipid_DE_allfat_and_CL.R 完成。")
log_msg("  下一步：可以基于 de_*.tsv 做火山图 / 热图 / 富集分析等。")
log_msg("============================================================")