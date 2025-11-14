#!/usr/bin/env Rscript

## ============================================================
## 批次内矫正：PQN（以 QC median profile 为参考）
##
## 输入（相对当前工作目录）：
##   results/lipid/tables/lipid_long_dedup_with_qcflag.tsv
##   results/lipid/tables/lipid_long_for_batch_correction.tsv
##
## 输出：
##   results/lipid/qc/lipid_pqn_factors_per_sample.tsv
##   results/lipid/tables/lipid_long_dedup_pqn.tsv
##
## 约定的关键列名：
##   batch, sample_id, lipidName, intensity,
##   is_qc, keep_for_correction
##
## 特别约定：
##   - 现在所有批次 (batch1–batch4) 的 QC 都参与 PQN 因子估计
## ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
  library(purrr)
})

## -------- 0. 一些路径和设置 ---------------------------------------------------

tables_dir <- "results/lipid/tables"
qc_dir     <- "results/lipid/qc"

dedup_qcflag_path   <- file.path(tables_dir, "lipid_long_dedup_with_qcflag.tsv")
for_correction_path <- file.path(tables_dir, "lipid_long_for_batch_correction.tsv")

## 哪些批次的 QC 用于驱动 PQN 矫正
##   - batch1: 仅监测，不驱动 => 不在这里列出
batches_use_qc_for_pqn <- c("batch1", "batch2", "batch3", "batch4")

## -------- 1. 读取数据 -------------------------------------------------------

message("[INFO] Reading input tables ...")

dedup_qcflag <- readr::read_tsv(dedup_qcflag_path, guess_max = 1e6)
for_corr     <- readr::read_tsv(for_correction_path, guess_max = 1e6)

## ---- 基本列检查（尽量早暴露问题）-------------------------------------------
required_cols_for_corr <- c(
  "batch", "sample_id", "lipidName", "intensity",
  "is_qc", "keep_for_correction"
)

missing_cols <- setdiff(required_cols_for_corr, colnames(for_corr))
if (length(missing_cols) > 0) {
  stop(
    "[ERROR] '", for_correction_path, "' 缺少必要列： ",
    paste(missing_cols, collapse = ", ")
  )
}

required_cols_dedup <- c("batch", "sample_id", "lipidName", "intensity")
missing_cols2 <- setdiff(required_cols_dedup, colnames(dedup_qcflag))
if (length(missing_cols2) > 0) {
  stop(
    "[ERROR] '", dedup_qcflag_path, "' 缺少必要列： ",
    paste(missing_cols2, collapse = ", ")
  )
}

## 确保 is_qc / keep_for_correction 为逻辑型，方便后面筛选
to_logical <- function(x) {
  if (is.logical(x)) return(x)
  if (is.numeric(x)) return(x != 0)
  if (is.character(x)) {
    return(toupper(x) %in% c("TRUE", "T", "1", "YES", "Y"))
  }
  as.logical(x)
}

for_corr <- for_corr %>%
  mutate(
    is_qc = to_logical(is_qc),
    keep_for_correction = to_logical(keep_for_correction)
  )

## 构造用于 PQN 的数据：
## - batch1：使用 dedup_qcflag 中的全部样本行（QC + 非 QC），并在本步骤中强制 keep_for_correction = TRUE
## - 其他批次：使用 for_corr（已经是 keep_for_correction == TRUE 的子集）
for_pqn_batch1 <- dedup_qcflag %>%
  filter(batch == "batch1") %>%
  mutate(keep_for_correction = TRUE)

for_pqn <- bind_rows(
  for_pqn_batch1 %>% select(colnames(for_corr)),
  for_corr
) %>%
  distinct()

## -------- 2. 构建 batch 内的 QC reference profile ---------------------------
##    只在:
##      - keep_for_correction == TRUE
##      - is_qc == TRUE
##      - batch ∈ batches_use_qc_for_pqn
##    的数据上构建 reference
##    => batch1 的 QC 不参与 ref 估计，只监测

message("[INFO] Building QC reference profiles for PQN ...")

ref_profile <- for_pqn %>%
  filter(
    keep_for_correction,
    is_qc,
    batch %in% batches_use_qc_for_pqn
  ) %>%
  group_by(batch, lipidName) %>%
  summarise(
    ref_intensity = median(intensity, na.rm = TRUE),
    n_qc = sum(!is.na(intensity)),
    .groups = "drop"
  ) %>%
  ## 去掉 ref 为 NA 或 <=0 的 feature
  filter(!is.na(ref_intensity), ref_intensity > 0)

if (nrow(ref_profile) == 0) {
  stop("[ERROR] 在 batches_use_qc_for_pqn 中没有有效的 QC reference feature（检查 keep_for_correction / is_qc / intensity）")
}

## -------- 3. 计算每个样本的 PQN factor -------------------------------------

message("[INFO] Computing PQN factors per sample (within batch) ...")

## 样本层信息（是否 QC）
sample_meta <- for_pqn %>%
  distinct(batch, sample_id, is_qc)

## 在 keep_for_correction 的 feature 上，与 reference 做比值
## 仅对 batches_use_qc_for_pqn 进行
ratios <- for_pqn %>%
  filter(
    keep_for_correction,
    batch %in% batches_use_qc_for_pqn
  ) %>%
  inner_join(ref_profile, by = c("batch", "lipidName")) %>%
  mutate(
    ratio = intensity / ref_intensity
  )

## 对每个 (batch, sample_id) 求 ratio 的中位数 = PQN 原始因子
pqn_factors_raw <- ratios %>%
  group_by(batch, sample_id) %>%
  summarise(
    pqn_factor_raw = median(
      ratio[is.finite(ratio) & ratio > 0],
      na.rm = TRUE
    ),
    n_features_used = sum(is.finite(ratio) & ratio > 0),
    .groups = "drop"
  )

## 合并 is_qc 信息
pqn_factors <- pqn_factors_raw %>%
  left_join(sample_meta, by = c("batch", "sample_id"))

## 对每个 batch，用 QC 样本的 median factor 作为基准，将 factor 重新居中到 1
##   - 仅对 batches_use_qc_for_pqn 生效
##   - batch1 没有出现在 pqn_factors_raw 中，后面会统一设为 factor = 1

pqn_factors <- pqn_factors %>%
  group_by(batch) %>%
  mutate(
    median_factor_qc = median(pqn_factor_raw[is_qc], na.rm = TRUE),
    pqn_factor = dplyr::if_else(
      batch %in% batches_use_qc_for_pqn &
        is.finite(median_factor_qc) & median_factor_qc > 0,
      pqn_factor_raw / median_factor_qc,
      pqn_factor_raw
    )
  ) %>%
  ungroup()

## 所有 batch（包括 batch1）的 factor 已包含在 pqn_factors 中
pqn_factors_all <- pqn_factors %>%
  arrange(batch, sample_id, is_qc)

## 简单提示一下每个 batch 的统计情况
batch_stats <- pqn_factors_all %>%
  group_by(batch) %>%
  summarise(
    n_samples       = n(),
    n_qc_samples    = sum(is_qc, na.rm = TRUE),
    n_with_factor   = sum(!is.na(pqn_factor_raw)),
    median_factor_qc = unique(median_factor_qc),
    .groups = "drop"
  )

message("[INFO] PQN factor summary by batch:")
print(batch_stats)

## -------- 4. 将 PQN 因子写入 QC 目录，便于检查 -----------------------------

dir.create(qc_dir, showWarnings = FALSE, recursive = TRUE)

pqn_factor_out <- pqn_factors_all %>%
  arrange(batch, is_qc, sample_id)

readr::write_tsv(
  pqn_factor_out,
  file.path(qc_dir, "lipid_pqn_factors_per_sample.tsv")
)

message("[INFO] Written PQN factors: ", file.path(qc_dir, "lipid_pqn_factors_per_sample.tsv"))

## -------- 5. 将样本级 PQN 因子回写到完整长表，并计算 intensity_pqn -------

message("[INFO] Applying PQN factors to full deduplicated table ...")

dedup_with_pqn <- dedup_qcflag %>%
  left_join(
    pqn_factors_all %>%
      select(batch, sample_id, pqn_factor, n_features_used),
    by = c("batch", "sample_id")
  ) %>%
  mutate(
    ## 对于没有算出 factor 的样本（极端异常或缺 feature），默认 factor = 1
    pqn_factor = ifelse(is.na(pqn_factor) | !is.finite(pqn_factor) | pqn_factor <= 0,
                        1, pqn_factor),
    intensity_pqn = intensity / pqn_factor
  )

## 简单检查一下是否有样本 n_features_used 特别小（不含 batch1，因为它是 NA）
problem_samples <- pqn_factors_all %>%
  filter(batch %in% batches_use_qc_for_pqn,
         n_features_used < 20)

if (nrow(problem_samples) > 0) {
  message("[WARN] 以下样本用于 PQN 的 feature 数 < 20，建议人工检查：")
  print(problem_samples %>% arrange(n_features_used))
}

## -------- 6. 输出带 PQN 矫正结果的长表 -------------------------------

dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)

dedup_out_path <- file.path(tables_dir, "lipid_long_dedup_pqn.tsv")

readr::write_tsv(
  dedup_with_pqn,
  dedup_out_path
)

message("[OK] Written PQN-normalized table: ", dedup_out_path)
message("[DONE] 批次内 PQN 矫正完成。")