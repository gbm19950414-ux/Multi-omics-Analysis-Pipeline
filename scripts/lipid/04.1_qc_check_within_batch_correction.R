#!/usr/bin/env Rscript

## ============================================================
## 04.1_qc_check_within_batch_correction.R
##
## 目的：对“批次内 PQN 矫正”的质量与完整性进行 QA
##
## 依赖输入（相对项目根目录）：
##   results/lipid/tables/lipid_long_dedup.tsv
##   results/lipid/tables/lipid_long_dedup_pqn.tsv
##   results/lipid/qc/lipid_pqn_factors_per_sample.tsv
##
## 主要检查内容：
##   1) PQN factor 合理性（按 batch 的分布、极端值、QC factor 是否 ≈1）
##   2) 每个样本用于 PQN 的 feature 数量 n_features_used
##   3) 矫正前后 sample-level 强度统计 (total / median / p90) 对比
##   4) QC 样本在矫正前后的变化是否 ~1（不被错误矫正）
##
## 输出：
##   results/lipid/qc/qa_pqn_factor_summary.tsv
##   results/lipid/qc/qa_pqn_extreme_factors.tsv
##   results/lipid/qc/qa_pqn_sample_stats.tsv
##   results/lipid/qc/pqn_factor_by_batch.png
##   results/lipid/qc/pqn_n_features_used_by_batch.png
##   results/lipid/qc/pqn_ratio_total_intensity_by_batch.png
##   results/lipid/qc/pqn_ratio_total_intensity_qc_vs_nonqc.png
## ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(ggplot2)
})

## --------- 小工具函数 --------------------------------------------------------

stop_if_missing <- function(path) {
  if (!file.exists(path)) {
    stop("[ERROR] 找不到文件: ", path,
         "\n       请确认 04_dedup_and_qc_filter.R 与 04_within_batch_pqn.R 已成功运行。",
         call. = FALSE)
  }
}

## --------- 1. 读取输入 -------------------------------------------------------

message("============================================================")
message("[INFO] QA for within-batch PQN correction")
message("============================================================")

base_dir <- "results/lipid"
tables_dir <- file.path(base_dir, "tables")
qc_dir <- file.path(base_dir, "qc")

dir.create(qc_dir, showWarnings = FALSE, recursive = TRUE)

raw_path    <- file.path(tables_dir, "lipid_long_dedup.tsv")
pqn_path    <- file.path(tables_dir, "lipid_long_dedup_pqn.tsv")
factor_path <- file.path(qc_dir, "lipid_pqn_factors_per_sample.tsv")

stop_if_missing(raw_path)
stop_if_missing(pqn_path)
stop_if_missing(factor_path)

message("[STEP1] 读取长表（矫正前/后）和 PQN factor ...")

long_raw <- readr::read_tsv(raw_path, guess_max = 1e6, show_col_types = FALSE)
long_pqn <- readr::read_tsv(pqn_path, guess_max = 1e6, show_col_types = FALSE)
pqn_factors <- readr::read_tsv(factor_path, guess_max = 1e6, show_col_types = FALSE)

message("  [raw] nrow = ", nrow(long_raw), " ; ncol = ", ncol(long_raw))
message("  [pqn] nrow = ", nrow(long_pqn), " ; ncol = ", ncol(long_pqn))
message("  [factors] nrow = ", nrow(pqn_factors), " ; ncol = ", ncol(pqn_factors))

## 基础列检查
required_cols_long <- c("batch", "sample_id", "lipidName", "intensity")
missing_raw <- setdiff(required_cols_long, colnames(long_raw))
missing_pqn <- setdiff(c(required_cols_long, "intensity_pqn"), colnames(long_pqn))

if (length(missing_raw) > 0) {
  stop("[ERROR] lipid_long_dedup.tsv 缺少列: ", paste(missing_raw, collapse = ", "), call. = FALSE)
}
if (length(missing_pqn) > 0) {
  stop("[ERROR] lipid_long_dedup_pqn.tsv 缺少列: ", paste(missing_pqn, collapse = ", "), call. = FALSE)
}

required_cols_factor <- c("batch", "sample_id", "pqn_factor", "pqn_factor_raw", "n_features_used", "is_qc")
missing_factor <- setdiff(required_cols_factor, colnames(pqn_factors))
if (length(missing_factor) > 0) {
  stop("[ERROR] lipid_pqn_factors_per_sample.tsv 缺少列: ",
       paste(missing_factor, collapse = ", "), call. = FALSE)
}

## --------- 2. PQN factor 合理性检查 -----------------------------------------

message("")
message("==== 2. PQN factor 合理性检查 ====")

## 防止 is_qc 不是逻辑型
to_logical <- function(x) {
  if (is.logical(x)) return(x)
  if (is.numeric(x)) return(x != 0)
  if (is.character(x)) return(toupper(x) %in% c("TRUE", "T", "1", "YES", "Y"))
  as.logical(x)
}

pqn_factors <- pqn_factors %>%
  mutate(is_qc = to_logical(is_qc))

## 按 batch 汇总 factor 分布
factor_summary <- pqn_factors %>%
  group_by(batch) %>%
  summarise(
    n_samples       = n(),
    n_qc_samples    = sum(is_qc, na.rm = TRUE),
    median_factor_qc = median(pqn_factor_raw[is_qc], na.rm = TRUE),
    factor_min      = min(pqn_factor, na.rm = TRUE),
    factor_q1       = quantile(pqn_factor, 0.25, na.rm = TRUE),
    factor_median   = median(pqn_factor, na.rm = TRUE),
    factor_q3       = quantile(pqn_factor, 0.75, na.rm = TRUE),
    factor_max      = max(pqn_factor, na.rm = TRUE),
    n_features_used_median = median(n_features_used, na.rm = TRUE),
    .groups = "drop"
  )

readr::write_tsv(factor_summary,
                 file.path(qc_dir, "qa_pqn_factor_summary.tsv"))

message("[INFO] PQN factor summary by batch:")
print(factor_summary)

## 标记“极端 factor”的样本，用于后续人工查看
extreme_samples <- pqn_factors %>%
  filter(is.finite(pqn_factor)) %>%
  mutate(
    is_extreme = (pqn_factor < 0.2 | pqn_factor > 5)
  ) %>%
  filter(is_extreme) %>%
  arrange(batch, pqn_factor)

if (nrow(extreme_samples) > 0) {
  readr::write_tsv(
    extreme_samples,
    file.path(qc_dir, "qa_pqn_extreme_factors.tsv")
  )
  message("[WARN] 发现 PQN factor 极端样本（pqn_factor < 0.2 或 > 5），已写出: ",
          file.path(qc_dir, "qa_pqn_extreme_factors.tsv"))
} else {
  message("[INFO] 未发现 pqn_factor < 0.2 或 > 5 的极端样本。")
}

## 画出 pqn_factor 的分布（按 batch）
p_factor <- ggplot(pqn_factors, aes(x = pqn_factor, fill = is_qc)) +
  geom_histogram(alpha = 0.5, position = "identity", bins = 30) +
  facet_wrap(~ batch, scales = "free_y") +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_fill_manual(values = c("grey70", "tomato"),
                    name = "is_qc") +
  theme_bw() +
  labs(
    title = "PQN factor distribution by batch",
    x = "pqn_factor",
    y = "Count"
  )

ggsave(
  filename = file.path(qc_dir, "pqn_factor_by_batch.png"),
  plot = p_factor,
  width = 9, height = 5, dpi = 300
)

## 画出 n_features_used 的分布
p_nfeat <- ggplot(pqn_factors, aes(x = n_features_used, fill = is_qc)) +
  geom_histogram(alpha = 0.5, position = "identity", bins = 30) +
  facet_wrap(~ batch, scales = "free_y") +
  scale_fill_manual(values = c("grey70", "steelblue"),
                    name = "is_qc") +
  theme_bw() +
  labs(
    title = "Number of features used for PQN per sample",
    x = "n_features_used (keep_for_correction & valid ratio)",
    y = "Count"
  )

ggsave(
  filename = file.path(qc_dir, "pqn_n_features_used_by_batch.png"),
  plot = p_nfeat,
  width = 9, height = 5, dpi = 300
)

message("[INFO] PQN factor 与 n_features_used 分布图已输出。")

## --------- 3. 矫正前后 sample-level 强度统计对比 ---------------------------

message("")
message("==== 3. 矫正前后 sample-level intensity 统计对比 ====")

## 保证 raw / pqn 两张表的 (batch, sample_id, lipidName) 对齐
common_keys <- intersect(
  long_raw$batch %>% paste(long_raw$sample_id, long_raw$lipidName, sep = "||"),
  long_pqn$batch %>% paste(long_pqn$sample_id, long_pqn$lipidName, sep = "||")
)

if (length(common_keys) == 0) {
  stop("[ERROR] raw 与 pqn 表在 (batch, sample_id, lipidName) 上没有交集，请检查。", call. = FALSE)
}

## 为了安全起见，用 inner_join 对齐后再统计
long_joined <- long_raw %>%
  select(batch, sample_id, lipidName, intensity, is_qc) %>%
  inner_join(
    long_pqn %>%
      select(batch, sample_id, lipidName, intensity_pqn),
    by = c("batch", "sample_id", "lipidName")
  )

## 计算 sample-level 统计量（before & after）
sample_stats <- long_joined %>%
  group_by(batch, sample_id, is_qc) %>%
  summarise(
    n_feat          = sum(!is.na(intensity)),
    total_raw       = sum(intensity, na.rm = TRUE),
    median_raw      = median(intensity, na.rm = TRUE),
    p90_raw         = quantile(intensity, 0.9, na.rm = TRUE),
    total_pqn       = sum(intensity_pqn, na.rm = TRUE),
    median_pqn      = median(intensity_pqn, na.rm = TRUE),
    p90_pqn         = quantile(intensity_pqn, 0.9, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    ratio_total  = total_pqn / total_raw,
    ratio_median = median_pqn / median_raw,
    ratio_p90    = p90_pqn / p90_raw
  )

readr::write_tsv(
  sample_stats,
  file.path(qc_dir, "qa_pqn_sample_stats.tsv")
)

message("[INFO] sample-level before/after intensity 统计已写出: ",
        file.path(qc_dir, "qa_pqn_sample_stats.tsv"))

## --------- 4. 画矫正前后强度比值分布（重点看 QC vs 非 QC）------------------

## 只考虑 total_intensity 的比值（作为代表），按 batch & is_qc 看分布
## 去掉 total_raw 非正的样本，防止比值问题
sample_stats_valid <- sample_stats %>%
  filter(is.finite(ratio_total), total_raw > 0, total_pqn > 0)

p_ratio_batch <- ggplot(sample_stats_valid,
                        aes(x = ratio_total, fill = is_qc)) +
  geom_histogram(alpha = 0.5, position = "identity", bins = 30) +
  facet_wrap(~ batch, scales = "free_y") +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_fill_manual(values = c("grey70", "tomato"),
                    name = "is_qc") +
  theme_bw() +
  labs(
    title = "Ratio of total intensity (PQN / raw) by batch",
    x = "total_intensity_pqn / total_intensity_raw",
    y = "Sample count"
  )

ggsave(
  filename = file.path(qc_dir, "pqn_ratio_total_intensity_by_batch.png"),
  plot = p_ratio_batch,
  width = 9, height = 5, dpi = 300
)

## 拆开 QC vs 非 QC，在全体 batch 上看 boxplot
p_ratio_qc_box <- ggplot(sample_stats_valid,
                         aes(x = is_qc, y = ratio_total, fill = is_qc)) +
  geom_boxplot(outlier.alpha = 0.4) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_fill_manual(values = c("grey70", "tomato"),
                    name = "is_qc") +
  theme_bw() +
  labs(
    title = "Total intensity ratio (PQN / raw): QC vs non-QC",
    x = "is_qc",
    y = "total_intensity_pqn / total_intensity_raw"
  )

ggsave(
  filename = file.path(qc_dir, "pqn_ratio_total_intensity_qc_vs_nonqc.png"),
  plot = p_ratio_qc_box,
  width = 6, height = 5, dpi = 300
)

## --------- 5. 对 batch1 的特殊检查（如存在）-------------------------------

if ("batch1" %in% unique(pqn_factors$batch)) {
  message("")
  message("==== 5. batch1 特殊检查（预期“只监测不矫正”） ====")

  b1_factors <- pqn_factors %>% filter(batch == "batch1")
  b1_stats   <- sample_stats %>% filter(batch == "batch1")

  if (nrow(b1_factors) == 0) {
    message("[INFO] pqn_factors 中没有 batch1 记录（可能已在主脚本中排除），跳过该检查。")
  } else {
    message("[INFO] batch1 pqn_factor 摘要：")
    print(summary(b1_factors$pqn_factor))

    ## 检查 batch1 的 ratio_total 是否全部接近 1
    if (nrow(b1_stats) > 0) {
      message("[INFO] batch1 total_intensity ratio (PQN / raw) 摘要：")
      print(summary(b1_stats$ratio_total))
    }
  }
}

message("============================================================")
message("[DONE] 04.1_qc_check_within_batch_correction.R QA 流程完成。")
message("  请结合 qc 目录下的 qa_*.tsv 和 png 图，人工评估 PQN 批次内矫正是否合理。")
message("============================================================")