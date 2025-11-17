#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(ggplot2)
})

cat("============================================================\n")
cat("[INFO] 05.1_between_batch_panel_QA.R (panel-based QA)\n")
cat("============================================================\n\n")

## ==== 路径设定 ====
in_long_pqn      <- "results/lipid/tables/lipid_long_dedup_pqn.tsv"
in_long_pqn_bb   <- "results/lipid/tables/lipid_long_dedup_pqn_bb.tsv"
in_feature_qc    <- "results/lipid/qc/lipid_feature_qc_summary.tsv"
in_factors_panel <- "results/lipid/qc/lipid_between_batch_panel_scaling_factors.tsv"

out_dir_qc <- "results/lipid/qc"
dir.create(out_dir_qc, showWarnings = FALSE, recursive = TRUE)

## ==== 1. 读取 PQN 前/后长表 + feature_qc + panel factors ====
cat("[STEP1] 读取 PQN 前/后长表 & feature_qc & panel factors...\n")
long_pqn <- read_tsv(in_long_pqn, show_col_types = FALSE)
long_pqn_bb <- read_tsv(in_long_pqn_bb, show_col_types = FALSE)
feature_qc <- read_tsv(in_feature_qc, show_col_types = FALSE)
factors_panel <- read_tsv(in_factors_panel, show_col_types = FALSE)

cat(sprintf("  [long_pqn]    nrow = %d  ncol = %d\n", nrow(long_pqn), ncol(long_pqn)))
cat(sprintf("  [long_pqn_bb] nrow = %d  ncol = %d\n", nrow(long_pqn_bb), ncol(long_pqn_bb)))
cat(sprintf("  [feature_qc]  nrow = %d  ncol = %d\n", nrow(feature_qc), ncol(feature_qc)))
cat(sprintf("  [factors_panel] nrow = %d  ncol = %d\n\n", nrow(factors_panel), ncol(factors_panel)))

## ==== 2. 从长表构建 QC 层面的 pqn 中心值（per batch, per lipid）====
cat("[STEP2] 构建 QC 层面的 pqn 中心值（per batch, per lipid）...\n")

# 统一 lipid id，避免 m/z 小数误差
make_lipid_id <- function(df) {
  df %>%
    mutate(
      lipid_id = paste(
        lipidName,
        LipidIon,
        sprintf("%.4f", round(CalcMz, 4)),
        sep = "|"
      )
    )
}

qc_long <- long_pqn %>%
  filter(is_qc == TRUE, keep_for_correction == TRUE) %>%
  make_lipid_id() %>%
  group_by(batch, lipid_id) %>%
  summarise(
    qc_center_pqn = median(intensity, na.rm = TRUE),
    n_qc_pqn = sum(!is.na(intensity)),
    .groups = "drop"
  )

cat(sprintf("  [qc_long] 有 QC 记录的 (batch, lipid) 数量: %d\n\n", nrow(qc_long)))

batches <- sort(unique(long_pqn$batch))
cat("  批次列表:", paste(batches, collapse = ", "), "\n\n")

## ==== 3. Panel A: 重建与 05_between_batch_panel_scaling.R 一致的 panel 逻辑 ====
cat("[STEP3] Panel A QA（batch1–3 共同 lipid 的 ratio 分布）...\n")

# 这里必须和 05_between_batch_panel_scaling.R 一致：
THR_MEAN_B1 <- 0     # batch1 qc_mean 不作筛选
THR_SN_B1   <- 10    # batch1 qc_sn_median ≥ 10

# 为 batch1 构建 keep_for_panel_b1
b1_panel <- feature_qc %>%
  filter(batch == "batch1") %>%
  mutate(
    keep_for_panel_b1 = (
      n_qc_non_na >= 1 &
        !is.na(qc_mean) & qc_mean >= THR_MEAN_B1 &
        !is.na(qc_sn_median) & qc_sn_median >= THR_SN_B1
    )
  )

# 合并回 feature_qc，得到 keep_for_panel（Panel A 规则）
feature_qc_panel <- feature_qc %>%
  left_join(
    b1_panel %>%
      select(batch, lipidName, LipidIon, CalcMz, keep_for_panel_b1),
    by = c("batch", "lipidName", "LipidIon", "CalcMz")
  ) %>%
  mutate(
    keep_for_panel = dplyr::case_when(
      batch == "batch1"               ~ keep_for_panel_b1,
      batch %in% c("batch2", "batch3") ~ keep_for_correction,
      TRUE                            ~ FALSE
    )
  ) %>%
  make_lipid_id()

# Panel A 候选：batch1–3 内 keep_for_panel == TRUE 的 lipid
panelA_candidates <- feature_qc_panel %>%
  filter(batch %in% c("batch1", "batch2", "batch3"),
         keep_for_panel == TRUE) %>%
  distinct(batch, lipid_id)

# 真正 Panel A：在 batch1–3 都出现的 lipid_id
panelA_common_ids <- panelA_candidates %>%
  count(lipid_id, name = "n_batch") %>%
  filter(n_batch == 3) %>%
  pull(lipid_id)

cat(sprintf("  Panel A 公共 lipid 数: %d\n", length(panelA_common_ids)))

if (length(panelA_common_ids) > 0) {
  qc_panelA <- qc_long %>%
    filter(batch %in% c("batch1", "batch2", "batch3"),
           lipid_id %in% panelA_common_ids)
  
  # 为每个 lipid 构造参考中心值：三批次 QC center 的 median
  panelA_ref <- qc_panelA %>%
    group_by(lipid_id) %>%
    summarise(
      ref_center = median(qc_center_pqn, na.rm = TRUE),
      .groups = "drop"
    )
  
  panelA_ratio <- qc_panelA %>%
    inner_join(panelA_ref, by = "lipid_id") %>%
    mutate(
      ratio = ref_center / qc_center_pqn,
      log10_ratio = log10(ratio)
    )
  
  cat(sprintf("  Panel A ratio 记录数: %d\n", nrow(panelA_ratio)))
  
  # 画 boxplot：log10(ratio) 按 batch
  pA <- ggplot(panelA_ratio, aes(x = batch, y = log10_ratio, fill = batch)) +
    geom_boxplot(outlier.size = 0.8, alpha = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(
      title = "Panel A (batch1–3) ratio 分布：log10(ref / qc_center_pqn)",
      x = "batch",
      y = "log10(ratio)"
    ) +
    theme_bw() +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5)
    )
  
  out_png_A <- file.path(out_dir_qc, "qc_panelA_ratio_boxplot.png")
  ggsave(out_png_A, pA, width = 6, height = 4, dpi = 300)
  cat("[OK] Panel A ratio boxplot 写出:", out_png_A, "\n\n")
} else {
  cat("  [INFO] Panel A 公共 lipid 数为 0，跳过 Panel A QA。\n\n")
}

## ==== 4. Panel B QA（batch2–4 CL 公共 lipid 的 ratio 分布）====
cat("[STEP4] Panel B QA（batch2–4 公共 CL 的 ratio 分布）...\n")

# 用长表里的 Class 信息定义 CL
long_qc_panel <- long_pqn %>%
  filter(is_qc == TRUE, keep_for_correction == TRUE) %>%
  make_lipid_id() %>%
  select(batch, lipid_id, Class, intensity)

# 限定 CL + batch2–4
long_qc_CL <- long_qc_panel %>%
  filter(batch %in% c("batch2", "batch3", "batch4"),
         Class == "CL")

if (nrow(long_qc_CL) > 0) {
  qc_CL_center <- long_qc_CL %>%
    group_by(batch, lipid_id) %>%
    summarise(
      qc_center_pqn_CL = median(intensity, na.rm = TRUE),
      .groups = "drop"
    )
  
  # 找到在 batch2–4 都出现的 CL lipid
  panelB_ids <- qc_CL_center %>%
    count(lipid_id, name = "n_batch") %>%
    filter(n_batch == 3) %>%
    pull(lipid_id)
  
  cat(sprintf("  Panel B 公共 CL lipid 数: %d\n", length(panelB_ids)))
  
  if (length(panelB_ids) > 0) {
    qc_panelB <- qc_CL_center %>%
      filter(lipid_id %in% panelB_ids)
    
    panelB_ref <- qc_panelB %>%
      group_by(lipid_id) %>%
      summarise(
        ref_center = median(qc_center_pqn_CL, na.rm = TRUE),
        .groups = "drop"
      )
    
    panelB_ratio <- qc_panelB %>%
      inner_join(panelB_ref, by = "lipid_id") %>%
      mutate(
        ratio = ref_center / qc_center_pqn_CL,
        log10_ratio = log10(ratio)
      )
    
    cat(sprintf("  Panel B ratio 记录数: %d\n", nrow(panelB_ratio)))
    
    pB <- ggplot(panelB_ratio, aes(x = batch, y = log10_ratio, fill = batch)) +
      geom_boxplot(outlier.size = 0.8, alpha = 0.7) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      labs(
        title = "Panel B (batch2–4 CL) ratio 分布：log10(ref / qc_center_pqn)",
        x = "batch",
        y = "log10(ratio)"
      ) +
      theme_bw() +
      theme(
        legend.position = "none",
        plot.title = element_text(hjust = 0.5)
      )
    
    out_png_B <- file.path(out_dir_qc, "qc_panelB_ratio_boxplot.png")
    ggsave(out_png_B, pB, width = 6, height = 4, dpi = 300)
    cat("[OK] Panel B ratio boxplot 写出:", out_png_B, "\n\n")
  } else {
    cat("  [INFO] Panel B 公共 CL lipid 数为 0，跳过 Panel B QA。\n\n")
  }
} else {
  cat("  [INFO] 在长表中未检测到 batch2–4 的 CL QC 记录，跳过 Panel B QA。\n\n")
}

## ==== 5. Panel bridge QA：factor_A vs factor_B ====
cat("[STEP5] Panel bridge QA（factor_A vs factor_B）...\n")

bridge_df <- factors_panel %>%
  filter(!is.na(factor_between_panelA),
         !is.na(factor_between_panelB),
         batch %in% c("batch2", "batch3"))

if (nrow(bridge_df) > 0) {
  k_bridge <- median(bridge_df$factor_between_panelA / bridge_df$factor_between_panelB,
                     na.rm = TRUE)
  cat(sprintf("  [bridge] nrow = %d  k_bridge(median A/B) = %.3f\n",
              nrow(bridge_df), k_bridge))
  
  p_bridge <- ggplot(bridge_df,
                     aes(x = factor_between_panelB,
                         y = factor_between_panelA,
                         label = batch)) +
    geom_point(size = 3) +
    geom_text(vjust = -1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    labs(
      title = "Panel bridge: factor_A vs factor_B (batch2/3)",
      x = "factor_between_panelB",
      y = "factor_between_panelA"
    ) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  
  out_png_bridge <- file.path(out_dir_qc, "qc_panel_bridge_A_vs_B.png")
  ggsave(out_png_bridge, p_bridge, width = 5, height = 5, dpi = 300)
  cat("[OK] Panel bridge 图写出:", out_png_bridge, "\n\n")
} else {
  cat("  [INFO] 无法找到同时拥有 PanelA/PanelB 因子的 batch2/3 记录，跳过 bridge QA。\n\n")
}

## ==== 6. QC-only PCA（PQN vs PQN+panel）====
cat("[STEP6] QC-only PCA（PQN vs PQN+panel）...\n")

do_qc_pca <- function(long_df, value_col, out_png, title_prefix) {
  df_qc <- long_df %>%
    filter(is_qc == TRUE, keep_for_correction == TRUE) %>%
    make_lipid_id() %>%
    select(batch, sample_id, lipid_id, !!sym(value_col)) %>%
    mutate(val = log10(!!sym(value_col) + 1)) %>%
    select(batch, sample_id, lipid_id, val)
  
  mat <- df_qc %>%
    pivot_wider(names_from = lipid_id, values_from = val) %>%
    as.data.frame()
  
  if (nrow(mat) < 3 || ncol(mat) < 4) {
    cat("  [WARN] QC-only PCA: 行或列太少，跳过 PCA。\n")
    return(invisible(NULL))
  }
  
  row_info <- mat[, c("batch", "sample_id")]
  X <- mat[, !(names(mat) %in% c("batch", "sample_id")), drop = FALSE]
  X <- X[, colSums(!is.na(X)) > 0, drop = FALSE]
  X_scaled <- scale(X, center = TRUE, scale = TRUE)
  X_scaled[is.na(X_scaled)] <- 0
  
  pca <- prcomp(X_scaled, center = FALSE, scale. = FALSE)
  pca_df <- as.data.frame(pca$x[, 1:2])
  pca_df$batch <- row_info$batch
  pca_df$sample_id <- row_info$sample_id
  
  p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = batch, label = sample_id)) +
    geom_point(size = 3) +
    geom_text(vjust = -0.8, size = 3) +
    labs(
      title = paste0(title_prefix, " - QC-only PCA"),
      x = "PC1",
      y = "PC2"
    ) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  
  ggsave(out_png, p, width = 6, height = 5, dpi = 300)
  cat("  [OK] QC-only PCA 写出:", out_png, "\n")
}

out_pca_pqn   <- file.path(out_dir_qc, "qc_pca_qc_only_pqn.png")
out_pca_pqnbb <- file.path(out_dir_qc, "qc_pca_qc_only_pqn_bb.png")

do_qc_pca(long_pqn,    "intensity",        out_pca_pqn,   "PQN only")
do_qc_pca(long_pqn_bb, "intensity_pqn_bb", out_pca_pqnbb, "PQN + panel scaling")
cat("\n")

## ==== 7. 样本层面 median intensity QA（PQN vs PQN+panel）====
cat("[STEP7] 样本层面 median intensity QA...\n")

compute_sample_median <- function(long_df, value_col) {
  long_df %>%
    filter(keep_for_correction == TRUE,
           !is.na(!!sym(value_col)),
           !!sym(value_col) > 0) %>%
    group_by(batch, sample_id) %>%
    summarise(
      median_intensity = median(!!sym(value_col), na.rm = TRUE),
      .groups = "drop"
    )
}

sample_med_pqn <- compute_sample_median(long_pqn, "intensity") %>%
  mutate(type = "PQN")

sample_med_pqnbb <- compute_sample_median(long_pqn_bb, "intensity_pqn_bb") %>%
  mutate(type = "PQN+panel")

sample_med_all <- bind_rows(sample_med_pqn, sample_med_pqnbb) %>%
  mutate(log10_median = log10(median_intensity + 1))

p_med <- ggplot(sample_med_all,
                aes(x = batch, y = log10_median, fill = type)) +
  geom_boxplot(position = position_dodge(width = 0.7)) +
  labs(
    title = "样本 median intensity：PQN vs PQN+panel",
    x = "batch",
    y = "log10(median_intensity + 1)"
  ) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

out_png_med <- file.path(out_dir_qc, "qc_sample_median_before_after_boxplot.png")
ggsave(out_png_med, p_med, width = 7, height = 4, dpi = 300)
cat("[OK] 样本 median intensity QA 图写出:", out_png_med, "\n\n")

cat("============================================================\n")
cat("[DONE] 05.1_between_batch_panel_QA.R QA 流程完成。\n")
cat("  请查看 results/lipid/qc/ 下的 qc_*.png 图，\n")
cat("  重点关注：Panel A/B 分布、bridge 关系、QC-only PCA、\n")
cat("  以及样本 median intensity 在 PQN 与 PQN+panel 下的变化。\n")
cat("============================================================\n")