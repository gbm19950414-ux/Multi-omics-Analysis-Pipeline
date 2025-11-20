#!/usr/bin/env Rscript

## ============================================================
## 06_ephb1_downstream_plots.R
##
## 目的：
##   基于 05_rna_ephb1_downstream_scores.R 的结果表，
##   可视化 “EphB1 下游信号通路在各 batch 中的变化”：
##     1) 按跨 batch 综合效应 z_meta 排序的条形图
##     2) pathway × batch 的 z_pathway 热图
##
## 输入（默认，可通过 config 覆盖）：
##
##   1) ephb1_downstream_signaling_scores_per_batch.tsv
##        results/multiomics/tables/ephb1_downstream_signaling_scores_per_batch.tsv
##
##   2) 可选配置文件（命令行第 1 个参数）：
##        scripts/multiomics/00_ephb1_downstream_plots_config.yaml
##
##      支持字段：
##        ephb1_downstream_scores_per_batch: "results/multiomics/tables/ephb1_downstream_signaling_scores_per_batch.tsv"
##        plots_dir: "results/multiomics/plots"
##
## 输出：
##
##   1) EphB1 下游通路整体条形图（跨 batch 综合 Z）：
##        results/multiomics/plots/ephb1_downstream_bar_overall.png
##
##   2) EphB1 下游通路 pathway × batch 热图：
##        results/multiomics/plots/ephb1_downstream_heatmap_per_batch.png
##
## 说明：
##   - z_pathway > 0：该 EphB1 下游通路在 HO 中整体激活；
##   - z_pathway < 0：该通路在 HO 中整体受抑；
##   - z_meta：对每条通路跨 batch 进行 Stouffer 合并的综合 Z，
##             绝对值越大，说明该通路在所有 batch 中变化越稳定、效应越强。
##
## ============================================================

suppressPackageStartupMessages({
  library(yaml)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
  library(forcats)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

## --------- 0. 解析命令行参数 & 配置 -------------------------

args <- commandArgs(trailingOnly = TRUE)

config_path <- if (length(args) >= 1) {
  args[1]
} else {
  "scripts/multiomics/00_ephb1_downstream_plots_config.yaml"
}

if (file.exists(config_path)) {
  cat("[INFO] 使用配置文件: ", config_path, "\n", sep = "")
  cfg <- yaml::read_yaml(config_path)
} else {
  cat("[WARN] 未找到配置文件: ", config_path,
      "，将使用脚本内置默认路径。\n", sep = "")
  cfg <- list()
}

scores_path_default <- "results/multiomics/tables/ephb1_downstream_signaling_scores_per_batch.tsv"
plots_dir_default   <- "results/multiomics/plots"

scores_path <- cfg$ephb1_downstream_scores_per_batch %||% scores_path_default
plots_dir   <- cfg$plots_dir                          %||% plots_dir_default

cat("============================================================\n")
cat("[INFO] 06_ephb1_downstream_plots.R\n")
cat("  scores_path : ", scores_path, "\n", sep = "")
cat("  plots_dir   : ", plots_dir,   "\n", sep = "")
cat("============================================================\n\n")

## --------- 1. 读取 EphB1 下游通路 per-batch Z 表 ------------

if (!file.exists(scores_path)) {
  stop("[ERROR] 找不到 ephb1_downstream_signaling_scores_per_batch.tsv: ",
       scores_path, "\n请先运行 05_rna_ephb1_downstream_scores.R")
}

cat("[STEP1] 读取 EphB1 下游通路 per-batch Z 表 ...\n")

scores <- readr::read_tsv(scores_path, show_col_types = FALSE)

cat("  [scores] nrow = ", nrow(scores),
    " ; ncol = ", ncol(scores), "\n", sep = "")

required_cols <- c("pathway", "axis", "batch",
                   "n_genes", "mean_LFC", "sd_LFC", "se_LFC", "z_pathway")

if (!all(required_cols %in% colnames(scores))) {
  stop("[ERROR] scores 表缺少必要列。缺少: ",
       paste(setdiff(required_cols, colnames(scores)), collapse = ", "),
       "\n当前列名: ", paste(colnames(scores), collapse = ", "))
}

## 过滤掉没有基因 / 没有 Z 的组合
scores_use <- scores %>%
  filter(n_genes > 0, !is.na(z_pathway))

cat("  [scores_use] nrow = ", nrow(scores_use),
    " ; pathway 数 = ", length(unique(scores_use$pathway)),
    " ; batch 数 = ", length(unique(scores_use$batch)), "\n", sep = "")

## --------- 2. 计算每条通路跨 batch 的综合 Z -----------------

cat("\n[STEP2] 计算每条 EphB1 下游通路跨 batch 的综合 Z (z_meta) ...\n")

overall <- scores_use %>%
  group_by(pathway, axis) %>%
  summarise(
    n_batch       = n_distinct(batch),
    n_batch_nonNA = sum(!is.na(z_pathway)),
    mean_z_pathway = ifelse(
      n_batch_nonNA > 0,
      mean(z_pathway, na.rm = TRUE),
      NA_real_
    ),
    ## Stouffer 风格综合 Z（只看 z_pathway 本身）
    z_meta = ifelse(
      n_batch_nonNA > 0,
      sum(z_pathway, na.rm = TRUE) / sqrt(n_batch_nonNA),
      NA_real_
    ),
    mean_LFC_overall = ifelse(
      n_batch_nonNA > 0,
      mean(mean_LFC, na.rm = TRUE),
      NA_real_
    ),
    .groups = "drop"
  ) %>%
  arrange(desc(z_meta))

cat("  [overall] nrow = ", nrow(overall),
    " ; pathway 数 = ", length(unique(overall$pathway)), "\n", sep = "")

## 为可视化准备因子顺序（按 z_meta 从小到大，这样条形图从下到上强度递增）
overall <- overall %>%
  mutate(
    pathway_factor = forcats::fct_reorder(pathway, z_meta, .desc = FALSE)
  )

## --------- 3. 画条形图：EphB1 下游通路整体效应 -------------

cat("\n[STEP3] 绘制 EphB1 下游通路整体综合 Z 条形图 ...\n")

dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

p_bar <- ggplot(overall,
                aes(x = pathway_factor, y = z_meta, fill = z_meta)) +
  geom_col() +
  coord_flip() +
  labs(
    title = "EphB1 downstream signaling pathways (cross-batch meta Z)",
    x     = "Pathway",
    y     = "z_meta (Stouffer 合并 z_pathway)",
    fill  = "z_meta"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(size = 9),
    axis.text.y = element_text(size = 9),
    plot.title  = element_text(hjust = 0.5, face = "bold")
  )

out_bar <- file.path(plots_dir, "ephb1_downstream_bar_overall.png")
ggsave(out_bar, p_bar, width = 8, height = 5, dpi = 300)

cat("  [OK] 已保存条形图: ", out_bar, "\n", sep = "")

## --------- 4. 画热图：pathway × batch 的 z_pathway -----------

cat("\n[STEP4] 绘制 EphB1 下游通路径 × batch z_pathway 热图 ...\n")

## 为热图统一 pathway 顺序（沿用 overall 的顺序）
scores_heat <- scores_use %>%
  mutate(
    pathway = factor(pathway, levels = levels(overall$pathway_factor))
  )

p_heat <- ggplot(scores_heat,
                 aes(x = batch, y = pathway, fill = z_pathway)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B",
    midpoint = 0,
    name = "z_pathway"
  ) +
  labs(
    title = "EphB1 downstream signaling pathways per batch (z_pathway)",
    x     = "Batch",
    y     = "Pathway"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 9),
    plot.title  = element_text(hjust = 0.5, face = "bold"),
    panel.grid  = element_blank()
  )

out_heat <- file.path(plots_dir, "ephb1_downstream_heatmap_per_batch.png")
ggsave(out_heat, p_heat, width = 7, height = 5, dpi = 300)

cat("  [OK] 已保存热图: ", out_heat, "\n", sep = "")

cat("============================================================\n")
cat("[DONE] 06_ephb1_downstream_plots.R 完成。\n")
cat("  - 请查看:\n")
cat("      ", out_bar,  "\n", sep = "")
cat("      ", out_heat, "\n", sep = "")
cat("============================================================\n")