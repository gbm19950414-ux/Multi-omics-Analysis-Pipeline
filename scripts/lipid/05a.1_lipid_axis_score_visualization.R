#!/usr/bin/env Rscript

## ============================================================
## 11_lipid_axis_score_visualization.R
##
## 目的：
##   使用 10_lipid_mechanistic_indicators_per_batch.R 的输出：
##     - lipid_mechanistic_indicators_per_sample.tsv
##     - lipid_mechanistic_axis_scores_per_batch_group.tsv
##   对机制轴分数做可视化：
##     1) per-sample 箱线图（axis × batch × group）
##     2) per-batch × group × axis 的热图
##     3) 按 batch 的“雷达风格”极坐标折线图（axis 做角度，score 做半径）
##
## 使用方式：
##   Rscript scripts/lipid/11_lipid_axis_score_visualization.R \
##     scripts/lipid/00_lipid_downstream_config.yaml
##
## 若不传 config，将使用脚本内默认路径：
##   tables_dir = results/lipid/tables
##   plots_dir  = results/lipid/plots
## ============================================================

suppressPackageStartupMessages({
  library(yaml)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(pheatmap)
  library(tibble)
})

## --------- 0. 解析命令行参数 & 配置 -------------------------

args <- commandArgs(trailingOnly = TRUE)

config_path <- if (length(args) >= 1) {
  args[1]
} else {
  "scripts/lipid/00_lipid_downstream_config.yaml"
}

tables_dir_default <- "results/lipid/tables"
plots_dir_default  <- "results/lipid/plots"

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
plots_dir  <- cfg$plots_dir  %||% plots_dir_default

## 输入文件
per_sample_path <- file.path(
  tables_dir, "lipid_mechanistic_indicators_per_sample.tsv"
)
per_batch_group_axis_path <- file.path(
  tables_dir, "lipid_mechanistic_axis_scores_per_batch_group.tsv"
)

## 输出图路径
plot_box_path <- file.path(
  plots_dir, "lipid_axis_score_boxplot_by_batch_group.png"
)
plot_heatmap_path <- file.path(
  plots_dir, "lipid_axis_score_heatmap_batch_group.png"
)
plot_radar_path <- file.path(
  plots_dir, "lipid_axis_score_radar_by_batch.png"
)

cat("============================================================\n")
cat("[INFO] 11_lipid_axis_score_visualization.R\n")
cat("  input per-sample      : ", per_sample_path, "\n", sep = "")
cat("  input per-batch-group : ", per_batch_group_axis_path, "\n", sep = "")
cat("  output boxplot  : ", plot_box_path, "\n", sep = "")
cat("  output heatmap  : ", plot_heatmap_path, "\n", sep = "")
cat("  output radar(polar) : ", plot_radar_path, "\n", sep = "")
cat("============================================================\n\n")

## --------- 1. 读取数据 --------------------------------------

if (!file.exists(per_sample_path)) {
  stop("[ERROR] 找不到 per-sample 指标表: ", per_sample_path)
}

if (!file.exists(per_batch_group_axis_path)) {
  stop("[ERROR] 找不到 per-batch × group × axis 汇总表: ",
       per_batch_group_axis_path)
}

cat("[STEP1] 读取 per-sample & per-batch-group 表...\n")

per_sample <- readr::read_tsv(per_sample_path, show_col_types = FALSE)
per_batch_group_axis <- readr::read_tsv(per_batch_group_axis_path,
                                        show_col_types = FALSE)

cat("  [per_sample] nrow = ", nrow(per_sample),
    " ; ncol = ", ncol(per_sample), "\n", sep = "")
cat("  [per_batch_group_axis] nrow = ", nrow(per_batch_group_axis),
    " ; ncol = ", ncol(per_batch_group_axis), "\n", sep = "")

required_cols_ps <- c("batch", "sample_id", "group")
missing_ps <- setdiff(required_cols_ps, colnames(per_sample))
if (length(missing_ps) > 0) {
  stop("[ERROR] per-sample 表缺少必要列: ",
       paste(missing_ps, collapse = ", "))
}

required_cols_bg <- c("batch", "group", "axis", "mean_score")
missing_bg <- setdiff(required_cols_bg, colnames(per_batch_group_axis))
if (length(missing_bg) > 0) {
  stop("[ERROR] per-batch-group 表缺少必要列: ",
       paste(missing_bg, collapse = ", "))
}

## 识别 axis_score 列（形如 Synthesis_score, Remodeling_score ...）
axis_cols <- grep("_score$", colnames(per_sample), value = TRUE)

if (length(axis_cols) == 0) {
  stop("[ERROR] 在 per_sample 表中未找到 *_score 列，请检查 10 脚本输出。")
}

cat("  [INFO] 识别到 axis_score 列:\n    ",
    paste(axis_cols, collapse = ", "), "\n", sep = "")

## 只保留 group 非 NA 的生物学样本
per_sample_bio <- per_sample %>%
  filter(!is.na(group))

## --------- 2. 箱线图：axis × batch × group ------------------

cat("\n[STEP2] 生成 axis score 的 per-batch × group 箱线图...\n")

## long 格式：batch × sample × axis × score
axis_long_ps <- per_sample_bio %>%
  select(batch, sample_id, group, all_of(axis_cols)) %>%
  pivot_longer(
    cols = all_of(axis_cols),
    names_to = "axis_col",
    values_to = "axis_score"
  ) %>%
  mutate(
    axis = stringr::str_remove(axis_col, "_score$")
  )

## 排序 axis（固定顺序：Synthesis, Remodeling, Oxidation, Transport, Supply）
axis_levels <- c("Synthesis", "Remodeling", "Oxidation", "Transport", "Supply")
axis_long_ps <- axis_long_ps %>%
  mutate(axis = factor(axis, levels = axis_levels))

## 简单过滤掉全 NA 的 axis/batch/ group 组合
axis_long_ps_non_na <- axis_long_ps %>%
  filter(!is.na(axis_score))

## 绘图
p_box <- ggplot(axis_long_ps_non_na,
                aes(x = group, y = axis_score, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 1) +
  facet_grid(axis ~ batch, scales = "free_y") +
  theme_bw(base_size = 12) +
  labs(
    title = "Mechanistic axis scores by batch and group",
    x = "Group",
    y = "Axis score (higher = bottleneck heavier)"
  ) +
  theme(
    strip.background = element_rect(fill = "grey92"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)
ggsave(plot_box_path, p_box, width = 9, height = 6, dpi = 300)

cat("  [OK] 写出箱线图: ", plot_box_path, "\n", sep = "")

## --------- 3. 热图：axis × (batch_group) ---------------------

cat("\n[STEP3] 生成 axis mean score 的热图 (axis × batch_group)...\n")

## 构造列名 batch_group：如 batch1_WT, batch1_HO ...
heat_df <- per_batch_group_axis %>%
  filter(!is.na(group)) %>%
  mutate(
    axis = factor(axis, levels = axis_levels),
    batch_group = paste0(batch, "_", group)
  ) %>%
  select(axis, batch_group, mean_score)

## wide：行 = axis，列 = batch_group
heat_mat_df <- heat_df %>%
  pivot_wider(
    names_from = batch_group,
    values_from = mean_score
  ) %>%
  arrange(axis)

## rownames = axis
heat_mat <- heat_mat_df %>%
  column_to_rownames(var = "axis") %>%
  as.matrix()

## --- 关键修正：处理 heat_mat 中的 NA，避免 hclust/dist 出错 ---

n_na_heat <- sum(is.na(heat_mat))

if (n_na_heat > 0) {
  cat("  [WARN] Heatmap 矩阵中检测到 ", n_na_heat,
      " 个 NA 单元格。仅用于可视化，",
      "这些 NA 将被填为 0（轴的“中性”水平）。原始 TSV 保留 NA。\n",
      sep = "")
  heat_mat[is.na(heat_mat)] <- 0
}

## pheatmap
pheatmap::pheatmap(
  heat_mat,
  scale = "none",
  clustering_method = "complete",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  main = "Axis mean scores (batch × group)",
  fontsize = 10,
  angle_col = 45
)

ggsave(plot_heatmap_path, width = 7, height = 5, dpi = 300)

cat("  [OK] 写出热图: ", plot_heatmap_path, "\n", sep = "")

## --------- 4. 雷达风格极坐标折线图（per-batch） --------------

cat("\n[STEP4] 生成按 batch 的雷达风格极坐标折线图...\n")

## 从 per_batch_group_axis 里提取每个 batch × group × axis 的 mean_score
radar_df <- per_batch_group_axis %>%
  filter(!is.na(group)) %>%
  mutate(
    axis = factor(axis, levels = axis_levels)
  )

## 为 polar plot 准备：每个 batch 一块 facet，x=axis, y=mean_score, group=group
## 为了让线闭合，可选地把首个 axis 再复制一遍到末尾
radar_df_closed <- radar_df %>%
  group_by(batch, group) %>%
  arrange(axis, .by_group = TRUE) %>%
  mutate(axis_index = as.numeric(axis)) %>%
  bind_rows(
    radar_df %>%
      filter(!is.na(group)) %>%
      mutate(
        axis = factor(axis, levels = axis_levels),
        axis_index = as.numeric(axis)
      ) %>%
      group_by(batch, group) %>%
      arrange(axis_index, .by_group = TRUE) %>%
      slice(1L) %>%
      mutate(axis_index = axis_index + length(axis_levels))
  ) %>%
  ungroup()

## 画极坐标折线图
p_radar <- ggplot(radar_df_closed,
                  aes(x = axis_index, y = mean_score,
                      group = group, color = group)) +
  geom_line(size = 0.8) +
  geom_point(size = 1.8) +
  coord_polar() +
  facet_wrap(~ batch) +
  theme_bw(base_size = 11) +
  labs(
    title = "Mechanistic axis scores (radar-like) per batch",
    x = NULL,
    y = "Axis mean score"
  ) +
  theme(
    axis.text.x = element_blank(), ## 用下面自定义标签
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

## 给极坐标添加 axis 名字标签：使用 annotation_custom 会比较麻烦，
## 这里直接在外圈加文本（在原 axis_index 位置）
## 简单做法：生成一个辅助数据框，单独 geom_text

axis_label_df <- data.frame(
  axis = factor(axis_levels, levels = axis_levels),
  axis_index = seq_along(axis_levels)
)

## 取所有 batch 的 y 最大值，作为标签半径
max_score <- max(radar_df$mean_score, na.rm = TRUE)
label_radius <- max_score * 1.2

p_radar <- p_radar +
  geom_text(
    data = axis_label_df,
    aes(x = axis_index, y = label_radius, label = axis),
    inherit.aes = FALSE,
    size = 3
  )

ggsave(plot_radar_path, p_radar, width = 8, height = 6, dpi = 300)

cat("  [OK] 写出雷达风格极坐标图: ", plot_radar_path, "\n", sep = "")

cat("============================================================\n")
cat("[DONE] 11_lipid_axis_score_visualization.R 完成。\n")
cat("  产出：\n")
cat("    - 箱线图：  batch × axis × group 层面的 score 分布\n")
cat("    - 热图  ：  axis × (batch_group) 的 mean score\n")
cat("    - 极坐标：  每个 batch 的 HO vs WT 机制谱对比\n")
cat("============================================================\n")