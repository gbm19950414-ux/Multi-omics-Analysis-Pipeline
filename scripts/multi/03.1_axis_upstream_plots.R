#!/usr/bin/env Rscript

## ============================================================
## 04_axis_upstream_plots.R
##
## 目的：
##   基于 03_axis_upstream_integration.R 的结果表，
##   可视化 “axis ↔ upstream pathway” 的整体关系：
##     1) 每个 axis 的 topN 上游通路支持度条形图
##     2) axis × pathway 的支持度热图（可限制为每个 axis 的 topN）
##
## 输入（默认，可通过 config 覆盖）：
##
##   1) axis_upstream_support_overall.tsv
##        results/multiomics/tables/axis_upstream_support_overall.tsv
##
##   2) axis_upstream_Z_per_batch.tsv（目前仅用于备份，可选）：
##        results/multiomics/tables/axis_upstream_Z_per_batch.tsv
##
##   3) 可选配置文件（命令行第 1 个参数）：
##        scripts/multiomics/00_axis_upstream_plots_config.yaml
##
##      支持字段：
##        axis_upstream_overall: "results/multiomics/tables/axis_upstream_support_overall.tsv"
##        axis_upstream_per_batch: "results/multiomics/tables/axis_upstream_Z_per_batch.tsv"
##        plots_dir: "results/multiomics/plots"
##        topN_per_axis: 6
##
## 输出：
##
##   1) 每个 axis 的 topN 上游通路条形图：
##        results/multiomics/plots/axis_upstream_bar_topN.png
##
##   2) axis × pathway 支持度热图（topN × axis）：
##        results/multiomics/plots/axis_upstream_heatmap_topN.png
##
## 说明：
##   - support_Z > 0 且绝对值较大：通路变化方向与 axis 效应方向一致，
##                                   是“上游驱动”的候选。
##   - 条形图：x = pathway（按 support_Z 排序），y = support_Z，按 axis 分 facet。
##   - 热图：y = axis，x = pathway，填色 = support_Z。
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
  "scripts/multiomics/00_axis_upstream_plots_config.yaml"
}

if (file.exists(config_path)) {
  cat("[INFO] 使用配置文件: ", config_path, "\n", sep = "")
  cfg <- yaml::read_yaml(config_path)
} else {
  cat("[WARN] 未找到配置文件: ", config_path,
      "，将使用脚本内置默认路径。\n", sep = "")
  cfg <- list()
}

overall_path_default   <- "results/multiomics/tables/axis_upstream_support_overall.tsv"
per_batch_path_default <- "results/multiomics/tables/axis_upstream_Z_per_batch.tsv"
plots_dir_default      <- "results/multiomics/plots"
topN_default           <- 6L

overall_path   <- cfg$axis_upstream_overall   %||% overall_path_default
per_batch_path <- cfg$axis_upstream_per_batch %||% per_batch_path_default
plots_dir      <- cfg$plots_dir               %||% plots_dir_default
topN_per_axis  <- cfg$topN_per_axis           %||% topN_default

cat("============================================================\n")
cat("[INFO] 04_axis_upstream_plots.R\n")
cat("  axis_upstream_overall : ", overall_path,   "\n", sep = "")
cat("  axis_upstream_per_batch: ", per_batch_path, "\n", sep = "")
cat("  plots_dir             : ", plots_dir,      "\n", sep = "")
cat("  topN_per_axis         : ", topN_per_axis,  "\n", sep = "")
cat("============================================================\n\n")

## --------- 1. 读取 axis_upstream_support_overall.tsv ---------

if (!file.exists(overall_path)) {
  stop("[ERROR] 找不到 axis_upstream_support_overall.tsv: ", overall_path,
       "\n请先运行 03_axis_upstream_integration.R")
}

cat("[STEP1] 读取 axis_upstream_support_overall.tsv ...\n")

overall <- readr::read_tsv(overall_path, show_col_types = FALSE)

cat("  [overall] nrow = ", nrow(overall),
    " ; ncol = ", ncol(overall), "\n", sep = "")

required_cols <- c("axis", "pathway", "n_batch_nonNA",
                   "mean_Z_axis", "mean_z_pathway",
                   "cor_Z_axis_pathway", "support_Z")

if (!all(required_cols %in% colnames(overall))) {
  stop("[ERROR] overall 表缺少必要列。缺少: ",
       paste(setdiff(required_cols, colnames(overall)), collapse = ", "),
       "\n当前列名: ", paste(colnames(overall), collapse = ", "))
}

## 过滤掉没有有效 batch 的组合
overall_use <- overall %>%
  filter(n_batch_nonNA > 0, !is.na(support_Z))

cat("  [overall_use] nrow = ", nrow(overall_use),
    " ; axis 数 = ", length(unique(overall_use$axis)),
    " ; pathway 数 = ", length(unique(overall_use$pathway)), "\n", sep = "")

## --------- 2. 为每个 axis 选出 topN pathway -----------------

cat("\n[STEP2] 为每个 axis 按 support_Z 选取 topN pathway ...\n")

topN_per_axis <- as.integer(topN_per_axis)

topN_tbl <- overall_use %>%
  group_by(axis) %>%
  arrange(desc(support_Z), .by_group = TRUE) %>%
  slice_head(n = topN_per_axis) %>%
  ungroup()

cat("  [topN_tbl] nrow = ", nrow(topN_tbl),
    " ; axis 数 = ", length(unique(topN_tbl$axis)),
    " ; pathway 数 = ", length(unique(topN_tbl$pathway)), "\n", sep = "")

## 为可视化整理因子顺序：
##   - 在每个 axis 内，pathway 按 support_Z 从小到大排列，这样条形图从左到右递增。
topN_tbl <- topN_tbl %>%
  group_by(axis) %>%
  mutate(
    pathway_factor = forcats::fct_reorder(pathway, support_Z, .desc = FALSE)
  ) %>%
  ungroup()

## 对 axis 也按整体最大 support_Z 排序（可选）
axis_order <- topN_tbl %>%
  group_by(axis) %>%
  summarise(max_support = max(support_Z, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(max_support)) %>%
  pull(axis)

topN_tbl <- topN_tbl %>%
  mutate(
    axis = factor(axis, levels = axis_order)
  )

## --------- 3. 画条形图：每个 axis 的 topN pathway -----------

cat("\n[STEP3] 绘制每个 axis 的 topN pathway 支持度条形图 ...\n")

dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

p_bar <- ggplot(topN_tbl,
                aes(x = pathway_factor, y = support_Z, fill = support_Z)) +
  geom_col() +
  facet_wrap(~ axis, scales = "free_x") +
  coord_flip() +
  labs(
    title    = "Top upstream pathways per axis (by support_Z)",
    x        = "Upstream pathway",
    y        = "support_Z (方向一致支持度)",
    fill     = "support_Z"
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.background = element_rect(color = NA, fill = "grey90"),
    axis.text.x      = element_text(size = 9),
    axis.text.y      = element_text(size = 9),
    plot.title       = element_text(hjust = 0.5, face = "bold")
  )

out_bar <- file.path(plots_dir, "axis_upstream_bar_topN.png")
ggsave(out_bar, p_bar, width = 10, height = 6, dpi = 300)

cat("  [OK] 已保存条形图: ", out_bar, "\n", sep = "")

## --------- 4. 画热图：axis × pathway（topN 内） -------------

cat("\n[STEP4] 绘制 axis × pathway 支持度热图 (topN per axis) ...\n")

## 为了在热图上横轴顺眼一些，把 pathway 做一个全局顺序：
## 可以用 support_Z 的全局均值或中位数排序
heat_tbl <- topN_tbl %>%
  mutate(
    pathway_global = pathway  # 可以保留原名
  )

## 全局 pathway 排序：按 support_Z 的绝对值或均值，这里用均值
pathway_order <- heat_tbl %>%
  group_by(pathway_global) %>%
  summarise(
    mean_support = mean(support_Z, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(mean_support) %>%
  pull(pathway_global)

heat_tbl <- heat_tbl %>%
  mutate(
    pathway_global = factor(pathway_global, levels = pathway_order)
  )

p_heat <- ggplot(heat_tbl,
                 aes(x = pathway_global, y = axis, fill = support_Z)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B",
    midpoint = 0,
    name = "support_Z"
  ) +
  labs(
    title = "Axis–upstream pathway support (topN per axis)",
    x     = "Upstream pathway",
    y     = "Axis"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, vjust = 1, size = 8),
    axis.text.y  = element_text(size = 10),
    plot.title   = element_text(hjust = 0.5, face = "bold"),
    panel.grid   = element_blank()
  )

out_heat <- file.path(plots_dir, "axis_upstream_heatmap_topN.png")
ggsave(out_heat, p_heat, width = 10, height = 5, dpi = 300)

cat("  [OK] 已保存热图: ", out_heat, "\n", sep = "")

cat("============================================================\n")
cat("[DONE] 04_axis_upstream_plots.R 完成。\n")
cat("  - 请查看:\n")
cat("      ", out_bar,  "\n", sep = "")
cat("      ", out_heat, "\n", sep = "")
cat("============================================================\n")