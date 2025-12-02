#!/usr/bin/env Rscript

## ============================================================
## figure_7_d.R
##
## 目的：
##   可视化中介分析筛选后的“top hits”结果：
##   - M1: Ephb1 下游信号通路 (ephb1_downstream_signaling_sets)
##   - M2: CL 轴上游 RNA 调控通路 (rna_upstream_pathway_sets)
##
## 输入（固定路径）：
##   - results/multiomics/tables/mediation/mediation_top_hits_M1.tsv
##   - results/multiomics/tables/mediation/mediation_top_hits_M2.tsv
##
## 输出：
##   - results/figs/figure_7_d_mediation_top_hits_M1.pdf
##   - results/figs/figure_7_d_mediation_top_hits_M2.pdf
##
## 图形：
##   - 每个点 = 一个 (mediator_col, outcome_col) 组合
##   - x 轴：间接效应 (优先 ab_boot_mean，其次 ab)
##   - 水平误差条：95% CI (优先 ab_boot_ci_low/high，其次 ab ± 1.96*se_ab_sobel)
##   - y 轴：paste0(axis, " : ", mediator_pretty)
##   - 颜色：CL 机制轴（由 outcome_col 的 *_score 推断）
## ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(ggplot2)
})

## ---------------- 辅助：推断轴名 / mediator 简名 ----------------

infer_axis_from_outcome <- function(outcome_col) {
  # outcome 列如 "Synthesis_score" -> "Synthesis"
  axis <- sub("_score$", "", outcome_col)
  axis
}

pretty_mediator_name <- function(mediator_col) {
  # mediator 列如 "Supply_mTOR_SREBP_lipogenesis"
  # 先尝试去掉前缀 "Supply_", 如果没有 "_" 就原样返回
  if (grepl("_", mediator_col)) {
    parts <- strsplit(mediator_col, "_", fixed = TRUE)[[1]]
    # 第一个部分大概率是轴名，如 "Supply" / "Synthesis"
    paste(parts[-1], collapse = "_")
  } else {
    mediator_col
  }
}

axis_level_order <- c(
  "Synthesis",
  "Remodeling",
  "Oxidation",
  "Transport",
  "Supply",
  "Membrane context"
)

## ---------------- 辅助：准备数据 ----------------

prepare_mediation_hits <- function(tsv_path) {

  if (!file.exists(tsv_path)) {
    stop("[ERROR] 找不到 top hits 表：", tsv_path)
  }

  dat <- readr::read_tsv(tsv_path, show_col_types = FALSE)

  required_cols <- c("mediator_col", "outcome_col")
  missing_cols  <- setdiff(required_cols, colnames(dat))
  if (length(missing_cols) > 0) {
    stop("[ERROR] top hits 表缺少必要列: ",
         paste(missing_cols, collapse = ", "),
         "\n实际列名: ", paste(colnames(dat), collapse = ", "))
  }

  ## 选择用于效应与 CI 的列
  effect_col <- if ("ab_boot_mean" %in% names(dat)) {
    "ab_boot_mean"
  } else if ("ab" %in% names(dat)) {
    "ab"
  } else {
    stop("[ERROR] top hits 表中既没有 ab_boot_mean 也没有 ab，无法作图。")
  }

  has_boot_ci <- all(c("ab_boot_ci_low", "ab_boot_ci_high") %in% names(dat))
  has_se_sobel <- "se_ab_sobel" %in% names(dat)

  dat2 <- dat %>%
    mutate(
      axis = infer_axis_from_outcome(outcome_col),
      axis = factor(axis,
                    levels = axis_level_order,
                    ordered = TRUE),
      mediator_pretty = vapply(mediator_col, pretty_mediator_name,
                               FUN.VALUE = character(1)),
      effect = .data[[effect_col]]
    )

  if (has_boot_ci) {
    dat2 <- dat2 %>%
      mutate(
        ci_low  = ab_boot_ci_low,
        ci_high = ab_boot_ci_high
      )
  } else if (has_se_sobel) {
    dat2 <- dat2 %>%
      mutate(
        ci_low  = effect - 1.96 * se_ab_sobel,
        ci_high = effect + 1.96 * se_ab_sobel
      )
  } else {
    # 没有 CI 信息，给 NA，错误条不会画出来
    dat2 <- dat2 %>%
      mutate(
        ci_low  = NA_real_,
        ci_high = NA_real_
      )
  }

  ## y 轴标签：轴名 : mediator 简名
  dat2 <- dat2 %>%
    mutate(
      pair_label = paste0(axis, " : ", mediator_pretty)
    )

  ## 按效应大小排序（y 轴从下到上）
  dat2 <- dat2 %>%
    mutate(
      pair_label = factor(pair_label,
                          levels = dat2$pair_label[order(dat2$effect)])
    )

  dat2
}

## ---------------- 辅助：画图 ----------------

plot_mediation_hits <- function(dat_hits, title_text) {

  p <- ggplot(dat_hits,
              aes(x = effect,
                  y = pair_label,
                  color = axis)) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_errorbarh(aes(xmin = ci_low, xmax = ci_high),
                   height = 0.2, alpha = 0.7, na.rm = TRUE) +
    geom_point(size = 2) +
    scale_color_brewer(palette = "Set1", na.value = "grey50") +
    labs(
      title = title_text,
      x = "Indirect effect a*b (Ephb1_KO → M → Y)",
      y = NULL,
      color = "Axis (Y)"
    ) +
    theme_bw(base_size = 10) +
    theme(
      plot.title = element_text(hjust = 0, face = "bold"),
      legend.position = "right"
    )

  p
}

## ---------------- 主流程：读入 M1/M2 并输出 PDF ----------------

cat("============================================================\n")
cat("[INFO] figure_7_d.R: 可视化中介分析 top hits\n")
cat("============================================================\n")

in_M1  <- "results/multiomics/tables/mediation/mediation_top_hits_M1.tsv"
in_M2  <- "results/multiomics/tables/mediation/mediation_top_hits_M2.tsv"
out_M1 <- "results/figs/figure_7_d_mediation_top_hits_M1.pdf"
out_M2 <- "results/figs/figure_7_d_mediation_top_hits_M2.pdf"

dir.create("results/figs", showWarnings = FALSE, recursive = TRUE)

cat("[STEP1] 读取并整理 M1 top hits: ", in_M1, "\n", sep = "")
hits_M1 <- prepare_mediation_hits(in_M1)

cat("[STEP2] 读取并整理 M2 top hits: ", in_M2, "\n", sep = "")
hits_M2 <- prepare_mediation_hits(in_M2)

cat("[STEP3] 生成并保存图形...\n")

p_M1 <- plot_mediation_hits(hits_M1,
                            title_text = "M1: Ephb1 downstream signaling as mediator")
p_M2 <- plot_mediation_hits(hits_M2,
                            title_text = "M2: Upstream CL-axis RNA pathways as mediator")

ggsave(out_M1, p_M1,
       width = 7, height = max(3, 0.3 * nrow(hits_M1)), units = "in")
ggsave(out_M2, p_M2,
       width = 7, height = max(3, 0.3 * nrow(hits_M2)), units = "in")

cat("  [OK] 已输出: ", out_M1, "\n", sep = "")
cat("  [OK] 已输出: ", out_M2, "\n", sep = "")
cat("============================================================\n")
cat("[DONE] figure_7_d 完成。\n")
cat("============================================================\n")