#!/usr/bin/env Rscript

# 06_explore_feature_qc_thresholds.R
# 用途：
#   - 读取 lipid_feature_qc_summary.tsv
#   - 对 qc_mean / qc_cv / qc_sn_median / qc_rt_sd 画直方图或密度图
#   - 在图上标出当前阈值位置，方便你判断阈值是否合理
# 注意：
#   - batch1 的 QC 你可以只当“监测”，不要用来驱动矫正；
#     本脚本仍会画出 batch1，方便你比较。
# cd /Volumes/Samsung_SSD_990_PRO_2TB_Media/multiomics_mech

# Rscript scripts/lipid/06_explore_feature_qc_thresholds.R \
  # results/lipid/qc/lipid_feature_qc_summary.tsv \
  # results/lipid/qc
  
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(stringr)
})

##--------------------------------------------------
## 1. 解析命令行参数
##--------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("用法: Rscript 06_explore_feature_qc_thresholds.R lipid_feature_qc_summary.tsv [out_dir]",
       call. = FALSE)
}

infile  <- args[1]
out_dir <- ifelse(length(args) >= 2, args[2], "results/lipid/qc")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("============================================================\n")
cat("[INFO] Explore QC feature thresholds\n")
cat("  infile : ", infile,  "\n")
cat("  outdir : ", out_dir, "\n")
cat("============================================================\n\n")

##--------------------------------------------------
## 2. 读取 feature 级 QC summary
##--------------------------------------------------

qc <- readr::read_tsv(infile, show_col_types = FALSE)
cat("[INFO] feature_qc 维度:", nrow(qc), "行 ×", ncol(qc), "列\n\n")

required_cols <- c("batch",
                   "qc_mean", "qc_cv",
                   "qc_sn_median", "qc_rt_sd")

missing <- setdiff(required_cols, names(qc))
if (length(missing) > 0) {
  stop("缺少必要列: ", paste(missing, collapse = ", "), call. = FALSE)
}

##--------------------------------------------------
## 3. 当前使用中的默认阈值（和 04_dedup_and_qc_filter.R 一致）
##    你可以在这里手动修改阈值，然后重跑画图。
##--------------------------------------------------

min_qc_mean      <- 1e5   # QC 平均强度 < 1e5 视为太弱
max_qc_cv        <- 0.3   # QC CV > 0.3 视为不稳定
min_qc_sn_median <- 10    # QC 中位 S/N < 10 视为噪声
max_qc_rt_sd     <- 0.05  # QC RT SD > 0.05 视为漂移大

cat("[INFO] 当前阈值设置:\n")
cat("  min_qc_mean      =", format(min_qc_mean, scientific = TRUE), "\n")
cat("  max_qc_cv        =", max_qc_cv, "\n")
cat("  min_qc_sn_median =", min_qc_sn_median, "\n")
cat("  max_qc_rt_sd     =", max_qc_rt_sd, "\n\n")

##--------------------------------------------------
## 4. 通用画图函数
##--------------------------------------------------

plot_metric_by_batch <- function(df, metric, threshold = NULL,
                                 direction = c(">=", "<="),
                                 log10_x = FALSE,
                                 filename_prefix = NULL,
                                 xlab = NULL,
                                 note = NULL) {
  direction <- match.arg(direction)

  if (!metric %in% names(df)) {
    warning("列不存在: ", metric, ", 跳过绘图。")
    return(invisible(NULL))
  }

  df_plot <- df %>%
    select(batch, !!sym(metric)) %>%
    rename(value = !!sym(metric)) %>%
    filter(!is.na(value) & is.finite(value))

  if (nrow(df_plot) == 0) {
    warning("列 ", metric, " 全部 NA 或非有限值，无法画图。")
    return(invisible(NULL))
  }

  if (log10_x) {
    df_plot <- df_plot %>%
      mutate(value_log10 = log10(value))
    x_var <- "value_log10"
    xlab_default <- paste0("log10(", metric, ")")
    thr_plot <- if (!is.null(threshold) && threshold > 0) log10(threshold) else NA_real_
  } else {
    x_var <- "value"
    xlab_default <- metric
    thr_plot <- threshold
  }

  if (is.null(xlab)) xlab <- xlab_default
  if (is.null(filename_prefix)) filename_prefix <- metric

  # 文本中说明阈值含义
  thr_text <- if (is.null(threshold)) {
    "no explicit threshold"
  } else if (direction == ">=") {
    paste0("threshold: ", metric, " ≥ ", threshold)
  } else {
    paste0("threshold: ", metric, " ≤ ", threshold)
  }

  if (!is.null(note)) {
    thr_text <- paste(thr_text, " | ", note)
  }

  p <- ggplot(df_plot, aes(x = .data[[x_var]])) +
    geom_histogram(bins = 50, alpha = 0.7) +
    facet_wrap(~ batch, scales = "free_y") +
    labs(
      title = paste0("Distribution of ", metric, " by batch"),
      x = xlab,
      y = "Count",
      caption = thr_text
    ) +
    theme_bw()

  if (!is.null(threshold) && is.finite(thr_plot)) {
    p <- p + geom_vline(xintercept = thr_plot,
                        linetype = "dashed",
                        color = "red")
  }

  out_file_png <- file.path(out_dir,
                            paste0("qc_", filename_prefix, "_by_batch.png"))
  out_file_pdf <- file.path(out_dir,
                            paste0("qc_", filename_prefix, "_by_batch.pdf"))

  ggsave(out_file_png, p, width = 8, height = 5, dpi = 300)
  ggsave(out_file_pdf, p, width = 8, height = 5)

  cat("[INFO] 已输出图: ", out_file_png, "\n")
  invisible(list(p = p,
                 file_png = out_file_png,
                 file_pdf = out_file_pdf))
}

##--------------------------------------------------
## 5. 对各个 QC 指标画图
##--------------------------------------------------

cat("==== 画 qc_mean 分布（log10 标尺）====\n")
plot_metric_by_batch(
  qc,
  metric          = "qc_mean",
  threshold       = min_qc_mean,
  direction       = ">=",
  log10_x         = TRUE,
  filename_prefix = "qc_mean",
  xlab            = "log10(qc_mean)",
  note            = "features with mean below threshold will be filtered"
)
cat("\n")

cat("==== 画 qc_cv 分布 ====\n")
plot_metric_by_batch(
  qc,
  metric          = "qc_cv",
  threshold       = max_qc_cv,
  direction       = "<=",
  log10_x         = FALSE,
  filename_prefix = "qc_cv",
  xlab            = "qc_cv",
  note            = "features with CV above threshold will be filtered"
)
cat("\n")

cat("==== 画 qc_sn_median 分布 ====\n")
# 先看看是否存在非 NA 的 qc_sn_median
if (any(!is.na(qc$qc_sn_median))) {
  plot_metric_by_batch(
    qc,
    metric          = "qc_sn_median",
    threshold       = min_qc_sn_median,
    direction       = ">=",
    log10_x         = FALSE,  # 如果分布跨度太大，可以手动改 TRUE
    filename_prefix = "qc_sn_median",
    xlab            = "qc_sn_median",
    note            = "features with median S/N below threshold will be filtered"
  )
} else {
  cat("[WARN] qc_sn_median 全部 NA，跳过图。\n")
}
cat("\n")

cat("==== 画 qc_rt_sd 分布 ====\n")
if (any(!is.na(qc$qc_rt_sd))) {
  plot_metric_by_batch(
    qc,
    metric          = "qc_rt_sd",
    threshold       = max_qc_rt_sd,
    direction       = "<=",
    log10_x         = FALSE,
    filename_prefix = "qc_rt_sd",
    xlab            = "qc_rt_sd",
    note            = "features with RT SD above threshold will be filtered"
  )
} else {
  cat("[WARN] qc_rt_sd 全部 NA，跳过图。\n")
}
cat("\n")

cat("============================================================\n")
cat("[DONE] QC metric exploration finished.\n")
cat("  请查看 ", out_dir, " 中的 qc_*.png / pdf 图，\n")
cat("  结合各 batch 的分布，决定是否调整阈值。\n")
cat("============================================================\n")