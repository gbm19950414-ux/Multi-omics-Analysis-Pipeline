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
##   ★ 额外集成：
##     4) 每个指标在 HO vs WT 的差异检验（per-batch）
##     5) 每个指标的箱线图（indicator × batch × group）
##     6) 指标 Δmean( HO − WT ) 的热图 + 显著性星号
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

## --- NEW: 指标级输出路径 ---
indicator_de_path <- file.path(
  tables_dir, "lipid_mechanistic_indicator_DE_per_batch.tsv"
)
plot_indicator_box_path <- file.path(
  plots_dir, "lipid_indicator_boxplot_by_batch_group.png"
)
plot_indicator_heatmap_path <- file.path(
  plots_dir, "lipid_indicator_heatmap_logFC_per_batch.png"
)

cat("============================================================\n")
cat("[INFO] 11_lipid_axis_score_visualization.R\n")
cat("  input per-sample      : ", per_sample_path, "\n", sep = "")
cat("  input per-batch-group : ", per_batch_group_axis_path, "\n", sep = "")
cat("  output boxplot  : ", plot_box_path, "\n", sep = "")
cat("  output heatmap  : ", plot_heatmap_path, "\n", sep = "")
cat("  output radar(polar) : ", plot_radar_path, "\n", sep = "")
cat("  output indicator DE table : ", indicator_de_path, "\n", sep = "")
cat("  output indicator boxplot  : ", plot_indicator_box_path, "\n", sep = "")
cat("  output indicator heatmap  : ", plot_indicator_heatmap_path, "\n", sep = "")
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

## 方便后续：axis 固定顺序
axis_levels <- c("Synthesis", "Remodeling", "Oxidation", "Transport", "Supply")

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

## 处理 NA，避免 pheatmap 内部 hclust 出错
n_na_heat <- sum(is.na(heat_mat))
if (n_na_heat > 0) {
  cat("  [WARN] Heatmap 矩阵中检测到 ", n_na_heat,
      " 个 NA 单元格。仅用于可视化，这些 NA 将被填为 0。\n", sep = "")
  heat_mat[is.na(heat_mat)] <- 0
}

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

radar_df <- per_batch_group_axis %>%
  filter(!is.na(group)) %>%
  mutate(
    axis = factor(axis, levels = axis_levels)
  )

## 复制首个 axis 到末尾，闭合曲线
radar_df_closed <- radar_df %>%
  group_by(batch, group) %>%
  arrange(axis, .by_group = TRUE) %>%
  mutate(axis_index = as.numeric(axis)) %>%
  ungroup() %>%
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
      mutate(axis_index = axis_index + length(axis_levels)) %>%
      ungroup()
  )

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
    axis.text.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

axis_label_df <- data.frame(
  axis = factor(axis_levels, levels = axis_levels),
  axis_index = seq_along(axis_levels)
)

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

## --------- 5. NEW: 指标级 HO vs WT 差异检验 ----------------

cat("\n[STEP5] 指标级 HO vs WT 差异检验 (per-batch)...\n")

## 识别“指标列”：所有 numeric 列里，排除 *_score（axis 分数）和显然不是指标的列
numeric_cols <- names(per_sample_bio)[sapply(per_sample_bio, is.numeric)]
indicator_cols <- setdiff(numeric_cols, axis_cols)

## 如果你有一些“明显非指标的 numeric 列”（比如某些 ID），可在这里排除：
indicator_cols <- setdiff(
  indicator_cols,
  c()  ## 如需排除的列名写在这里
)

cat("  [INFO] 识别到指标列 (", length(indicator_cols), "): ",
    paste(indicator_cols, collapse = ", "), "\n", sep = "")

## long 格式：batch × sample × group × indicator × value
indicator_long <- per_sample_bio %>%
  select(batch, sample_id, group, all_of(indicator_cols)) %>%
  pivot_longer(
    cols = all_of(indicator_cols),
    names_to = "indicator",
    values_to = "value"
  )

## 仅保留 HO / WT 两组，用于差异检验
indicator_long_hw <- indicator_long %>%
  filter(group %in% c("HO", "WT"))

## 小工具：对某个 batch × indicator 做检验
do_test_one <- function(df) {
  vals_HO <- df$value[df$group == "HO"]
  vals_WT <- df$value[df$group == "WT"]
  
  mean_HO <- if (length(vals_HO) > 0) mean(vals_HO, na.rm = TRUE) else NA_real_
  mean_WT <- if (length(vals_WT) > 0) mean(vals_WT, na.rm = TRUE) else NA_real_
  
  if (length(vals_HO) >= 2 && length(vals_WT) >= 2) {
    ## 这里用 t.test；如果你更喜欢 Wilcoxon，可改成 wilcox.test
    tt <- try(stats::t.test(vals_HO, vals_WT), silent = TRUE)
    pval <- if (inherits(tt, "try-error")) NA_real_ else tt$p.value
  } else {
    pval <- NA_real_
  }
  
  tibble(
    mean_HO = mean_HO,
    mean_WT = mean_WT,
    diff_HO_WT = mean_HO - mean_WT,
    p_value = pval
  )
}

indicator_de <- indicator_long_hw %>%
  group_by(batch, indicator) %>%
  group_modify(~ do_test_one(.x)) %>%
  ungroup() %>%
  mutate(
    FDR = p.adjust(p_value, method = "BH")
  )

dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)
readr::write_tsv(indicator_de, indicator_de_path)

cat("  [OK] 指标级差异检验结果写出: ", indicator_de_path, "\n", sep = "")

## --------- 6. NEW: 指标级箱线图 (indicator × batch × group) ----

cat("\n[STEP6] 生成指标级箱线图 (indicator × batch × group)...\n")

## 为了可读性，对 indicator 做一个固定顺序（按字母排序）
indicator_levels <- sort(unique(indicator_long_hw$indicator))

indicator_long_hw <- indicator_long_hw %>%
  mutate(
    indicator = factor(indicator, levels = indicator_levels)
  )

p_ind_box <- ggplot(indicator_long_hw,
                    aes(x = group, y = value, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 0.7) +
  facet_grid(indicator ~ batch, scales = "free_y") +
  theme_bw(base_size = 10) +
  labs(
    title = "Lipid mechanistic indicators by batch and group",
    x = "Group",
    y = "Indicator value"
  ) +
  theme(
    strip.background = element_rect(fill = "grey92"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

ggsave(plot_indicator_box_path, p_ind_box, width = 10, height = 10, dpi = 300)

cat("  [OK] 写出指标箱线图: ", plot_indicator_box_path, "\n", sep = "")

## --------- 7. NEW: 指标 Δmean 热图 + 显著性星号 ---------------

cat("\n[STEP7] 生成指标 Δmean( HO − WT ) 热图 + 显著性星号...\n")

## 准备绘图表：indicator × batch
indicator_de_plot <- indicator_de %>%
  mutate(
    batch = factor(batch, levels = sort(unique(batch))),
    indicator = factor(indicator, levels = indicator_levels),
    sig_label = case_when(
      is.na(p_value) ~ "",
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE ~ ""
    )
  )

p_ind_heat <- ggplot(indicator_de_plot,
                     aes(x = batch, y = indicator, fill = diff_HO_WT)) +
  geom_tile(color = "grey80") +
  geom_text(aes(label = sig_label), size = 2) +
  scale_fill_gradient2(
    name = "Δmean (HO − WT)",
    low = "blue", mid = "white", high = "red",
    midpoint = 0
  ) +
  theme_bw(base_size = 10) +
  labs(
    title = "Indicator-level HO − WT differences per batch",
    x = "Batch",
    y = "Indicator"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

ggsave(plot_indicator_heatmap_path, p_ind_heat, width = 8, height = 10, dpi = 300)

cat("  [OK] 写出指标 Δmean 热图: ", plot_indicator_heatmap_path, "\n", sep = "")

## --------- DONE ------------------------------------------------

cat("============================================================\n")
cat("[DONE] 11_lipid_axis_score_visualization.R 完成。\n")
cat("  产出：\n")
cat("    - Axis 级：\n")
cat("        · 箱线图：  batch × axis × group 层面的 score 分布\n")
cat("        · 热图  ：  axis × (batch_group) 的 mean score\n")
cat("        · 极坐标：  每个 batch 的 HO vs WT 机制谱对比\n")
cat("    - Indicator 级：\n")
cat("        · DE 表  ：  每个指标在 HO vs WT 的 per-batch 检验 (mean_HO, mean_WT, diff, p, FDR)\n")
cat("        · 箱线图：  indicator × batch × group 的原始值分布\n")
cat("        · 热图  ：  indicator × batch 的 Δmean(HO − WT) + 显著性星号\n")
cat("============================================================\n")