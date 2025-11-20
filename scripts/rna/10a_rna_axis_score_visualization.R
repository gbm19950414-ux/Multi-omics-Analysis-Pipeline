#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(ggplot2)
  library(pheatmap)
})

cat("============================================================\n")
cat("[INFO] 09b_rna_axis_score_visualization.R\n")
cat("  input axis_scores : results/rna/mechanistic_axis/rna_mechanistic_axis_scores_per_batch.tsv\n")
cat("  input gene_LFC    : results/rna/mechanistic_axis/rna_mechanistic_axis_gene_LFC_long.tsv\n")
cat("  outdir            : results/rna/mechanistic_axis\n")
cat("============================================================\n\n")

## ------------------------------------------------------------
## 基本路径（如需将来可传参，可以在这里加 commandArgs 逻辑）
## ------------------------------------------------------------
axis_score_path <- "results/rna/mechanistic_axis/rna_mechanistic_axis_scores_per_batch.tsv"
gene_lfc_path   <- "results/rna/mechanistic_axis/rna_mechanistic_axis_gene_LFC_long.tsv"
outdir          <- "results/rna/mechanistic_axis"

if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}

## ------------------------------------------------------------
## 读取表
## ------------------------------------------------------------
axis_scores <- read_tsv(axis_score_path, show_col_types = FALSE)
gene_lfc    <- read_tsv(gene_lfc_path,   show_col_types = FALSE)

cat("[STEP1] 读取轴级分数与基因 LFC 表...\n")
cat("  [axis_scores] nrow =", nrow(axis_scores), "; ncol =", ncol(axis_scores), "\n")
cat("  [gene_lfc]    nrow =", nrow(gene_lfc),    "; ncol =", ncol(gene_lfc),    "\n")

## 期望列检查（温和检查，如果缺失就报错退出）
required_axis_cols <- c("axis", "batch", "mean_LFC")
missing_axis_cols  <- setdiff(required_axis_cols, colnames(axis_scores))
if (length(missing_axis_cols) > 0) {
  stop("[FATAL] axis_scores 缺少以下列: ", paste(missing_axis_cols, collapse = ", "))
}

required_gene_cols <- c("axis", "batch", "GeneID", "log2FC")
missing_gene_cols  <- setdiff(required_gene_cols, colnames(gene_lfc))
if (length(missing_gene_cols) > 0) {
  stop("[FATAL] gene_lfc_long 缺少以下列: ", paste(missing_gene_cols, collapse = ", "))
}

## 统一 axis 排序
axis_levels <- c("Synthesis", "Remodeling", "Oxidation", "Transport", "Supply")

axis_scores <- axis_scores %>%
  mutate(
    axis  = factor(axis, levels = axis_levels),
    batch = factor(batch, levels = sort(unique(batch)))
  ) %>%
  arrange(axis, batch)

gene_lfc <- gene_lfc %>%
  mutate(
    axis  = factor(axis, levels = axis_levels),
    batch = factor(batch, levels = sort(unique(batch)))
  )

cat("[INFO] axis_levels: ", paste(axis_levels, collapse = ", "), "\n")
cat("[INFO] batches    : ", paste(levels(axis_scores$batch), collapse = ", "), "\n\n")

## ------------------------------------------------------------
## STEP2: 轴级均值热图 (axis × batch, value = mean_LFC)
## ------------------------------------------------------------
cat("[STEP2] 绘制轴级均值热图 (axis × batch, mean_LFC)...\n")

heat_df <- axis_scores %>%
  select(axis, batch, mean_LFC) %>%
  filter(!is.na(axis), !is.na(batch))

# wide -> matrix
heat_mat <- heat_df %>%
  mutate(
    axis  = as.character(axis),
    batch = as.character(batch)
  ) %>%
  pivot_wider(
    names_from  = batch,
    values_from = mean_LFC
  ) %>%
  column_to_rownames("axis") %>%
  as.matrix()

# 处理 NA：如果存在 NA，用 0 填充并给出提示
if (any(is.na(heat_mat))) {
  n_na <- sum(is.na(heat_mat))
  cat("[WARN] heat_mat 中存在 NA，数量 =", n_na, "，将以 0 替代以保证聚类可运行。\n")
  heat_mat[is.na(heat_mat)] <- 0
}

heatmap_out <- file.path(outdir, "rna_axis_scores_heatmap_batch_mean_LFC.png")

pheatmap(
  heat_mat,
  scale            = "none",
  clustering_method = "complete",
  main             = "Transcriptomic mechanistic axis scores (mean log2FC)\n(HO vs WT per batch)",
  filename         = heatmap_out,
  width            = 6,
  height           = 4
)

cat("  [OK] 热图写出: ", heatmap_out, "\n\n")

## ------------------------------------------------------------
## STEP3: 轴级柱状图（axis × mean_LFC，按 batch 分色）
## ------------------------------------------------------------
cat("[STEP3] 绘制轴级柱状图 (mean_LFC by axis × batch)...\n")

bar_out <- file.path(outdir, "rna_axis_scores_barplot_by_batch.png")

p_bar <- ggplot(axis_scores, aes(x = axis, y = mean_LFC, fill = batch)) +
  geom_col(position = position_dodge(width = 0.7)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    x = "Mechanistic axis",
    y = "Mean log2FC (HO vs WT)",
    title = "Transcriptomic mechanistic axis scores per batch"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1),
    plot.title   = element_text(hjust = 0.5)
  )

ggsave(bar_out, p_bar, width = 7, height = 4.5, dpi = 300)
cat("  [OK] 柱状图写出: ", bar_out, "\n\n")

## ------------------------------------------------------------
## STEP4: 轴内基因 LFC 分布箱线图 (每个 axis 在不同 batch 中的分布)
## ------------------------------------------------------------
cat("[STEP4] 绘制轴内基因 LFC 分布箱线图...\n")

boxplot_out <- file.path(outdir, "rna_axis_gene_LFC_boxplot_by_axis_batch.png")

# 只保留 axis 非 NA、LFC 非 NA
gene_lfc_clean <- gene_lfc %>%
  filter(!is.na(axis), !is.na(log2FC))

# 若行数太少给个友好提示
if (nrow(gene_lfc_clean) == 0) {
  cat("[WARN] gene_lfc_long 中没有有效的 axis × LFC 记录，跳过箱线图绘制。\n")
} else {
  p_box <- ggplot(gene_lfc_clean, aes(x = axis, y = log2FC, fill = batch)) +
    geom_boxplot(outlier.size = 0.8, position = position_dodge(width = 0.7)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(
      x = "Mechanistic axis",
      y = "Per-gene log2FC (HO vs WT)",
      title = "Distribution of gene-level log2FC within each axis and batch"
    ) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title  = element_text(hjust = 0.5)
    )

  ggsave(boxplot_out, p_box, width = 8, height = 5, dpi = 300)
  cat("  [OK] 箱线图写出: ", boxplot_out, "\n\n")
}
## ------------------------------------------------------------
## STEP5: 轴级分数雷达图 (per-batch radar / polar plot)
## ------------------------------------------------------------
cat("[STEP5] 绘制轴级分数雷达图 (per-batch radar plot)...\n")

radar_df <- axis_scores %>%
  # 只用有效值
  filter(!is.na(axis), !is.na(batch), !is.na(mean_LFC)) %>%
  mutate(
    axis  = factor(axis, levels = axis_levels),
    batch = factor(batch, levels = levels(axis_scores$batch))
  ) %>%
  # 在每个 batch 内，对各轴的 mean_LFC 做 z-score
  group_by(batch) %>%
  mutate(
    mean_LFC_z = as.numeric(scale(mean_LFC))
  ) %>%
  ungroup()

if (nrow(radar_df) == 0) {
  cat("  [WARN] 没有可用于雷达图的轴级分数记录，跳过雷达图绘制。\n\n")
} else {
  range_z <- range(radar_df$mean_LFC_z, na.rm = TRUE)
  # 给 y 轴留一点 buffer，避免线贴边
  y_min <- range_z[1] - 0.2
  y_max <- range_z[2] + 0.2

  radar_out <- file.path(outdir, "rna_axis_scores_radar_by_batch.png")

  p_radar <- ggplot(radar_df,
                    aes(x = axis,
                        y = mean_LFC_z,
                        group = batch,
                        color = batch)) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    coord_polar() +
    ylim(y_min, y_max) +
    labs(
      x = NULL,
      y = "Z-scored mean log2FC (per batch)",
      title = "Transcriptomic mechanistic axis radar plot (HO vs WT)"
    ) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x     = element_text(size = 10),
      axis.title.y    = element_text(vjust = 1.5),
      plot.title      = element_text(hjust = 0.5),
      panel.grid.major = element_line(linetype = "dotted", linewidth = 0.3),
      panel.grid.minor = element_blank()
    )

  ggsave(radar_out, p_radar, width = 6, height = 6, dpi = 300)
  cat("  [OK] 雷达图写出: ", radar_out, "\n\n")
}
cat("============================================================\n")
cat("[DONE] 09b_rna_axis_score_visualization.R 完成。\n")
cat("  输出图像位于: ", outdir, "\n")
cat("============================================================\n")