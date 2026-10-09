#!/usr/bin/env Rscript

## ============================================================
## 07_feature_filter_QA.R
##
## 功能：
##   1) 对 06_prepare_lipid_matrices.R 产生的 feature_filter_stats_*.tsv
##      做 QA（检测率 / 方差 / 保留原因）。
##   2) 合并 lipid_feature_qc_summary.tsv 中的 QC 指标，生成扩展 QA 表。
##   3) 画 feature-level 的 missing rate / variance 直方图。
##
## 约定：
##   - 配置文件：scripts/lipid/00_lipid_downstream_config.yaml
##   - 输入：
##       * feature_filter_stats_allfat.tsv
##       * feature_filter_stats_CL.tsv
##       * lipid_feature_qc_summary.tsv
##   - 输出：
##       * results/lipid/qc/qa_feature_filter_allfat.tsv
##       * results/lipid/qc/qa_feature_filter_CL.tsv
##       * results/lipid/qc/qa_feature_filter_allfat_with_qc.tsv
##       * results/lipid/qc/qa_feature_filter_CL_with_qc.tsv
##       * results/lipid/qc/qc_feature_missingrate_hist_allfat.png
##       * results/lipid/qc/qc_feature_missingrate_hist_CL.png
##       * results/lipid/qc/qc_feature_variance_hist_allfat.png
##       * results/lipid/qc/qc_feature_variance_hist_CL.png
## ============================================================

suppressPackageStartupMessages({
  library(yaml)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(ggplot2)
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
       "\nUsage: Rscript 07_feature_filter_QA.R path/to/00_lipid_downstream_config.yaml")
}

log_msg("[INFO] 使用配置文件: ", config_path)
cfg <- yaml::read_yaml(config_path)

paths <- cfg$paths
qc_dir <- paths$qc_dir %||% "results/lipid/qc"

dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

## 自动推断 feature_filter_stats 路径
get_filter_stats_path <- function(label){
  # 优先从 config: paths$feature_filter_stats_allfat / CL
  key <- paste0("feature_filter_stats_", label)
  p_cfg <- paths[[key]]
  if (!is.null(p_cfg)) return(p_cfg)

  # 否则根据矩阵路径推断
  if (label == "allfat") {
    mat_path <- paths$matrix_allfat_filtered %||% paths$matrix_allfat
  } else if (label == "CL") {
    mat_path <- paths$matrix_CL_filtered %||% paths$matrix_CL
  } else {
    stop("Unknown label: ", label)
  }
  file.path(dirname(mat_path), paste0("feature_filter_stats_", label, ".tsv"))
}

stats_path_allfat <- get_filter_stats_path("allfat")
stats_path_CL     <- get_filter_stats_path("CL")

feature_qc_path <- paths$feature_qc_summary %||% "results/lipid/qc/lipid_feature_qc_summary.tsv"

log_msg("[INFO] feature_filter_stats_allfat: ", stats_path_allfat)
log_msg("[INFO] feature_filter_stats_CL    : ", stats_path_CL)
log_msg("[INFO] lipid_feature_qc_summary   : ", feature_qc_path)

if (!file.exists(feature_qc_path)) {
  stop("找不到 feature_qc_summary 文件: ", feature_qc_path)
}

## ------------------------------------------------------------
## 1. 读取 QC summary，并汇总到 feature 级别
## ------------------------------------------------------------
log_msg("[STEP1] 读取 lipid_feature_qc_summary.tsv 并汇总到 feature_id 级别 ...")

feature_qc_raw <- readr::read_tsv(feature_qc_path, show_col_types = FALSE)

# 构造与 06 中相同的 feature_id：lipidName | LipidIon | CalcMz(4位小数)
feature_qc_raw <- feature_qc_raw %>%
  mutate(
    feature_id = paste0(
      lipidName, " | ", LipidIon, " | ",
      sprintf("%.4f", CalcMz)
    )
  )

# 汇总为 per-feature 的 QC 概览（跨 batch）
feature_qc_agg <- feature_qc_raw %>%
  group_by(feature_id) %>%
  summarise(
    n_batches_qc          = n_distinct(batch[n_qc_non_na > 0]),
    any_keep_for_corr     = any(keep_for_correction, na.rm = TRUE),
    qc_mean_median        = median(qc_mean, na.rm = TRUE),
    qc_cv_median          = median(qc_cv,   na.rm = TRUE),
    qc_sn_median_median   = median(qc_sn_median, na.rm = TRUE),
    qc_rt_sd_median       = median(qc_rt_sd, na.rm = TRUE),
    .groups = "drop"
  )

log_msg("[INFO] feature_qc_agg 维度: ", nrow(feature_qc_agg),
        " features × ", ncol(feature_qc_agg), " 列")

## ------------------------------------------------------------
## 2. 一个通用函数：对某个 universe 做 QA
## ------------------------------------------------------------
qa_one_universe <- function(label, stats_path){

  if (!file.exists(stats_path)) {
    warning("统计文件不存在，跳过: ", stats_path)
    return(invisible(NULL))
  }

  log_msg("------------------------------------------------------------")
  log_msg("[INFO] QA for universe: ", label)
  log_msg("[INFO] 读取 feature_filter_stats_", label, " ...")

  stats <- readr::read_tsv(stats_path, show_col_types = FALSE)

  log_msg("[INFO] 行数: ", nrow(stats), " ; 列数: ", ncol(stats))

  if (!"keep_feature" %in% colnames(stats)) {
    warning("stats 中缺少 keep_feature 列，无法 QA: ", stats_path)
    return(invisible(NULL))
  }

  ## 2.1 概要计数：保留/剔除
  n_total <- nrow(stats)
  n_keep  <- sum(stats$keep_feature, na.rm = TRUE)
  n_drop  <- n_total - n_keep

  log_msg("[INFO] 总 feature 数: ", n_total,
          " ; 保留: ", n_keep,
          " ; 剔除: ", n_drop)

  ## 2.2 剔除原因拆解（基于 pass_global / pass_group / pass_var）
  if (!all(c("pass_global", "pass_group", "pass_var") %in% colnames(stats))) {
    log_msg("[WARN] 缺少 pass_global / pass_group / pass_var，无法拆解剔除原因，仅输出 keep_feature 情况。")
    stats_reason <- stats %>%
      mutate(reason = if_else(keep_feature, "kept", "dropped")) %>%
      count(reason)
  } else {
    stats_reason <- stats %>%
      mutate(
        reason = case_when(
          keep_feature ~ "kept",
          !pass_global ~ "fail_global_detect_rate",
          pass_global & !pass_group ~ "fail_group_detect_rate",
          pass_global & pass_group & !pass_var ~ "fail_low_variance",
          TRUE ~ "dropped_other"
        )
      ) %>%
      count(reason)
  }

  log_msg("[INFO] 剔除原因分布：")
  print(stats_reason)

  ## 2.3 与 QC summary 合并，写出扩展 QA 表
  stats_qc <- stats %>%
    left_join(feature_qc_agg, by = "feature_id")

  qa_out_path <- file.path(qc_dir, paste0("qa_feature_filter_", label, ".tsv"))
  qa_out_qc   <- file.path(qc_dir, paste0("qa_feature_filter_", label, "_with_qc.tsv"))

  # 不带 QC 的 QA（其实 stats 本身已经含大部分信息，这里只保证路径统一）
  readr::write_tsv(stats, qa_out_path)
  log_msg("[OK] 写出 QA 表（不含 QC 汇总）: ", qa_out_path)

  # 含 QC 汇总的 QA
  readr::write_tsv(stats_qc, qa_out_qc)
  log_msg("[OK] 写出 QA 表（含 QC 汇总）: ", qa_out_qc)

  ## 2.4 画 histogram：detect_rate / variance
  # 为了避免 NA，对指标先做过滤
  # 1) detect_rate_global
  if ("detect_rate_global" %in% colnames(stats)) {
    p1 <- ggplot(stats, aes(x = detect_rate_global, fill = keep_feature)) +
      geom_histogram(position = "identity", alpha = 0.6, bins = 40) +
      scale_fill_manual(values = c("TRUE" = "#1b9e77", "FALSE" = "#d95f02")) +
      labs(
        title = paste0("Feature global detect rate (", label, ")"),
        x = "detect_rate_global",
        y = "count",
        fill = "keep_feature"
      ) +
      theme_bw()

    out_png1 <- file.path(qc_dir, paste0("qc_feature_missingrate_hist_", label, ".png"))
    ggsave(out_png1, p1, width = 7, height = 5, dpi = 300)
    log_msg("[OK] 写出 detect_rate_global 直方图: ", out_png1)
  } else {
    log_msg("[WARN] stats 中没有 detect_rate_global，跳过 missing rate 图: ", label)
  }

  # 2) var_log（log2 方差）
  if ("var_log" %in% colnames(stats)) {
    # 为了视角更清晰，可以画原始 var_log，也可以 log10(var_log+1e-8)
    p2 <- ggplot(stats, aes(x = var_log, fill = keep_feature)) +
      geom_histogram(position = "identity", alpha = 0.6, bins = 40) +
      scale_fill_manual(values = c("TRUE" = "#1b9e77", "FALSE" = "#d95f02")) +
      labs(
        title = paste0("Feature variance (log2 scale) (", label, ")"),
        x = "var_log",
        y = "count",
        fill = "keep_feature"
      ) +
      theme_bw()

    out_png2 <- file.path(qc_dir, paste0("qc_feature_variance_hist_", label, ".png"))
    ggsave(out_png2, p2, width = 7, height = 5, dpi = 300)
    log_msg("[OK] 写出 variance 直方图: ", out_png2)
  } else {
    log_msg("[WARN] stats 中没有 var_log，跳过 variance 图: ", label)
  }

  invisible(list(
    stats = stats,
    stats_qc = stats_qc,
    reason = stats_reason
  ))
}

## ------------------------------------------------------------
## 3. 分别对 allfat / CL 做 QA
## ------------------------------------------------------------
log_msg("============================================================")
log_msg("[INFO] 07_feature_filter_QA.R 开始 QA ...")
log_msg("============================================================")

res_allfat <- qa_one_universe("allfat", stats_path_allfat)
res_CL     <- qa_one_universe("CL",     stats_path_CL)

log_msg("============================================================")
log_msg("[DONE] 07_feature_filter_QA.R QA 流程完成。")
log_msg("  请查看 ", qc_dir, " 下的 qc_feature_*_hist_*.png 和 qa_feature_filter_*.tsv")
log_msg("============================================================")