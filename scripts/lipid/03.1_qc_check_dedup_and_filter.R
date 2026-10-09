#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(tibble)
})

## ------------------------------------------------------------
## 0. 解析参数 & 路径
## ------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
out_dir <- ifelse(length(args) >= 1, args[1], "results/lipid")

tables_dir <- file.path(out_dir, "tables")
qc_dir     <- file.path(out_dir, "qc")

file_all      <- file.path(tables_dir, "lipid_long_all_batches.tsv")
file_dedup    <- file.path(tables_dir, "lipid_long_dedup.tsv")
file_dedup_flag <- file.path(tables_dir, "lipid_long_dedup_with_qcflag.tsv")
file_for_corr <- file.path(tables_dir, "lipid_long_for_batch_correction.tsv")
file_feat_qc  <- file.path(qc_dir, "lipid_feature_qc_summary.tsv")
file_conflict <- file.path(qc_dir, "dedup_conflicts.tsv")

cat("============================================================\n")
cat("[INFO] QA for Dedup + QC filtering\n")
cat("  out_dir:", out_dir, "\n")
cat("============================================================\n\n")

## 简单检查文件存在性
files_to_check <- c(
  all_batches   = file_all,
  dedup         = file_dedup,
  dedup_withflag = file_dedup_flag,
  for_corr      = file_for_corr,
  feature_qc    = file_feat_qc
)

missing_files <- names(files_to_check)[!file.exists(files_to_check)]
if (length(missing_files) > 0) {
  cat("[FATAL] 以下关键文件不存在:\n")
  for (nm in missing_files) {
    cat("  -", nm, ":", files_to_check[[nm]], "\n")
  }
  stop("缺少关键输入文件，无法进行 QA 检查。", call. = FALSE)
}

## ------------------------------------------------------------
## 1. 读取主要表 & 基础结构检查
## ------------------------------------------------------------

cat("==== 1. 读取长表 & 基础结构检查 ====\n")

df_all   <- readr::read_tsv(file_all, show_col_types = FALSE)
df_dedup <- readr::read_tsv(file_dedup, show_col_types = FALSE)

cat("[all_batches] 行数:", nrow(df_all), "列数:", ncol(df_all), "\n")
cat("[dedup]       行数:", nrow(df_dedup), "列数:", ncol(df_dedup), "\n\n")

required_cols <- c("batch", "sample_id", "intensity",
                   "lipidName", "LipidIon", "CalcMz", "is_qc")

missing_all  <- setdiff(required_cols, names(df_all))
missing_dedup <- setdiff(required_cols, names(df_dedup))

if (length(missing_all) > 0) {
  stop("[all_batches] 缺少必要列: ", paste(missing_all, collapse = ", "), call. = FALSE)
}
if (length(missing_dedup) > 0) {
  stop("[dedup] 缺少必要列: ", paste(missing_dedup, collapse = ", "), call. = FALSE)
}

# 再保险：确认没有 "Group" 行
n_group_all   <- df_all  %>% filter(!is.na(lipidName) & lipidName == "Group") %>% nrow()
n_group_dedup <- df_dedup %>% filter(!is.na(lipidName) & lipidName == "Group") %>% nrow()

cat("  [CHECK] all_batches 中 lipidName == 'Group' 的行数:", n_group_all, "\n")
cat("  [CHECK] dedup 中 lipidName == 'Group' 的行数:", n_group_dedup, "\n\n")

## ------------------------------------------------------------
## 2. Dedup 质量检查
## ------------------------------------------------------------

cat("==== 2. Dedup 质量检查 ====\n")

key_cols <- c("batch", "sample_id", "lipidName", "LipidIon", "CalcMz")

## 2.1 在原始长表上重新统计重复组（理论值）
dup_groups_all <- df_all %>%
  group_by(across(all_of(key_cols))) %>%
  summarise(
    n_rows      = n(),
    intensity_n = n_distinct(intensity, na.rm = FALSE),
    qc_grade_n  = if ("qc_grade"  %in% names(.)) n_distinct(qc_grade,  na.rm = FALSE) else NA_integer_,
    qc_sn_n     = if ("qc_sn"     %in% names(.)) n_distinct(qc_sn,     na.rm = FALSE) else NA_integer_,
    qc_rt_n     = if ("qc_rt"     %in% names(.)) n_distinct(qc_rt,     na.rm = FALSE) else NA_integer_,
    qc_height_n = if ("qc_height" %in% names(.)) n_distinct(qc_height, na.rm = FALSE) else NA_integer_,
    .groups = "drop"
  ) %>%
  mutate(
    is_dup = n_rows > 1,
    all_equal = is_dup &
      intensity_n <= 1 &
      (is.na(qc_grade_n)  | qc_grade_n  <= 1) &
      (is.na(qc_sn_n)     | qc_sn_n     <= 1) &
      (is.na(qc_rt_n)     | qc_rt_n     <= 1) &
      (is.na(qc_height_n) | qc_height_n <= 1)
  )

n_dup_total_all    <- dup_groups_all %>% filter(is_dup) %>% nrow()
n_dup_equal_all    <- dup_groups_all %>% filter(all_equal) %>% nrow()
n_dup_conflict_all <- dup_groups_all %>% filter(is_dup & !all_equal) %>% nrow()

cat("  [原始表 Dedup 统计]\n")
cat("    - 有重复 key 组数:", n_dup_total_all, "\n")
cat("      · 其中完全相同(可安全合并)组数:", n_dup_equal_all, "\n")
cat("      · 其中存在数值冲突组数:", n_dup_conflict_all, "\n")

# 理论上减少的行数
rows_in_dup_groups <- dup_groups_all %>%
  filter(is_dup) %>%
  summarise(total_rows_dup = sum(n_rows)) %>%
  pull(total_rows_dup)

expected_dedup_nrow <- nrow(df_all) - (rows_in_dup_groups - n_dup_total_all)

cat("  [行数比对]\n")
cat("    - 原始总行数         :", nrow(df_all), "\n")
cat("    - 理论 Dedup 后行数  :", expected_dedup_nrow, "\n")
cat("    - 实际 Dedup 后行数  :", nrow(df_dedup), "\n")
cat("    - 实际 vs 理论是否一致:", ifelse(nrow(df_dedup) == expected_dedup_nrow, "YES", "NO"), "\n\n")

## 2.2 检查 Dedup 后是否仍有重复键
dup_after <- df_dedup %>%
  group_by(across(all_of(key_cols))) %>%
  summarise(n_rows = n(), .groups = "drop") %>%
  filter(n_rows > 1)

n_dup_after <- nrow(dup_after)
cat("  [Dedup 后剩余重复键检查]\n")
cat("    - Dedup 后仍有重复键组数:", n_dup_after, "\n")
if (n_dup_after > 0) {
  cat("    [WARN] 仍存在重复键，前几条如下:\n")
  print(dup_after %>% head(10))
}
cat("\n")

## 2.3 冲突组选峰策略 QA（如果有 dedup_conflicts 文件）
if (file.exists(file_conflict)) {
  cat("  [冲突组选峰 QA]\n")
  df_conf <- readr::read_tsv(file_conflict, show_col_types = FALSE)
  if (nrow(df_conf) == 0) {
    cat("    - dedup_conflicts.tsv 为空，无需检查\n\n")
  } else {
    has_qc_sn <- "qc_sn" %in% names(df_conf)

    # 为每个冲突 key 计算理论最佳行的 qc_sn / intensity
    best_theoretical <- df_conf %>%
      mutate(
        .score_sn  = if (has_qc_sn) qc_sn else NA_real_,
        .score_int = intensity
      ) %>%
      group_by(across(all_of(key_cols))) %>%
      summarise(
        qc_sn_max  = if (has_qc_sn) max(.score_sn, na.rm = TRUE) else NA_real_,
        int_max    = max(.score_int, na.rm = TRUE),
        .groups = "drop"
      )

    # 从 dedup 结果中取出对应 key 的行，看是否匹配理论最佳
    dedup_for_conf <- df_dedup %>%
      inner_join(best_theoretical,
                 by = key_cols) %>%
      mutate(
        match_sn  = if (has_qc_sn) (qc_sn == qc_sn_max) else TRUE,
        match_int = (intensity == int_max)
      )

    n_conf_key <- nrow(best_theoretical)
    n_match_sn <- if (has_qc_sn) sum(dedup_for_conf$match_sn, na.rm = TRUE) else NA_integer_
    n_match_int <- sum(dedup_for_conf$match_int, na.rm = TRUE)

    cat("    - 冲突 key 组总数:", n_conf_key, "\n")
    if (has_qc_sn) {
      cat("    - 选中 qc_sn 最大行的 key 数:", n_match_sn, "\n")
    }
    cat("    - 选中 intensity 最大行的 key 数:", n_match_int, "\n\n")
  }
} else {
  cat("  [冲突组选峰 QA]\n")
  cat("    - 未找到 dedup_conflicts.tsv，跳过该步\n\n")
}

## ------------------------------------------------------------
## 3. Feature 级 QC 统计 QA
## ------------------------------------------------------------

cat("==== 3. Feature 级 QC 统计 QA ====\n")

feature_qc <- readr::read_tsv(file_feat_qc, show_col_types = FALSE)
dedup_flag <- readr::read_tsv(file_dedup_flag, show_col_types = FALSE)
for_corr  <- readr::read_tsv(file_for_corr, show_col_types = FALSE)

cat("[feature_qc] 行数:", nrow(feature_qc), "列数:", ncol(feature_qc), "\n")
cat("[dedup_withflag] 行数:", nrow(dedup_flag), "列数:", ncol(dedup_flag), "\n")
cat("[for_batch_correction] 行数:", nrow(for_corr), "列数:", ncol(for_corr), "\n\n")

## 3.1 有 QC 记录的 feature 是否都出现在 feature_qc 里？

# 在 dedup 表中：有 is_qc == TRUE 的 feature
feat_with_qc_in_dedup <- dedup_flag %>%
  filter(is_qc) %>%
  distinct(batch, lipidName, LipidIon, CalcMz)

n_feat_with_qc_in_dedup <- nrow(feat_with_qc_in_dedup)

# 在 feature_qc 中：真正有 QC 行的 feature（n_qc_total > 0 或 !is.na(n_qc_total)）
if ("n_qc_total" %in% names(feature_qc)) {
  feat_with_qc_in_summary <- feature_qc %>%
    filter(!is.na(n_qc_total) & n_qc_total > 0) %>%
    distinct(batch, lipidName, LipidIon, CalcMz)
} else {
  # 保险兜底：用 n_qc_non_na
  feat_with_qc_in_summary <- feature_qc %>%
    filter(!is.na(n_qc_non_na) & n_qc_non_na > 0) %>%
    distinct(batch, lipidName, LipidIon, CalcMz)
}

n_feat_with_qc_in_summary <- nrow(feat_with_qc_in_summary)

cat("  [有 QC 记录的 feature 数]\n")
cat("    - 在 dedup 表中 (is_qc == TRUE):", n_feat_with_qc_in_dedup, "\n")
cat("    - 在 feature_qc 中 (n_qc_total > 0):", n_feat_with_qc_in_summary, "\n")

# 检查二者集合是否一致
join_qc_sets <- feat_with_qc_in_dedup %>%
  full_join(feat_with_qc_in_summary,
            by = c("batch", "lipidName", "LipidIon", "CalcMz"),
            suffix = c(".dedup", ".summary")) %>%
  mutate(
    in_dedup   = !is.na(batch),
    in_summary = !is.na(batch)
  )

n_only_in_dedup   <- join_qc_sets %>% filter(in_dedup & !in_summary) %>% nrow()
n_only_in_summary <- join_qc_sets %>% filter(!in_dedup & in_summary) %>% nrow()

cat("    - 仅在 dedup 中出现但不在 summary 中的 feature 数:", n_only_in_dedup, "\n")
cat("    - 仅在 summary 中出现但不在 dedup 中的 feature 数:", n_only_in_summary, "\n\n")

## 3.2 各 QC 条件 & keep_for_correction 的统计（按 batch）

qc_cond_cols <- c("pass_n_qc", "pass_mean", "pass_cv", "pass_sn", "pass_rt", "keep_for_correction")
missing_qc_cond <- setdiff(qc_cond_cols, names(feature_qc))
if (length(missing_qc_cond) > 0) {
  cat("[WARN] feature_qc 中缺少以下 QC 条件列，将跳过部分统计:\n  ",
      paste(missing_qc_cond, collapse = ", "), "\n\n")
}

# 按 batch 做一个 summary 表
feat_qc_by_batch <- feature_qc %>%
  group_by(batch) %>%
  summarise(
    n_feature_total  = n(),
    n_feature_with_qc = sum(!is.na(n_qc_non_na) & n_qc_non_na > 0, na.rm = TRUE),
    n_pass_n_qc      = if ("pass_n_qc" %in% names(.)) sum(pass_n_qc, na.rm = TRUE) else NA_integer_,
    n_pass_mean      = if ("pass_mean" %in% names(.)) sum(pass_mean, na.rm = TRUE) else NA_integer_,
    n_pass_cv        = if ("pass_cv" %in% names(.)) sum(pass_cv, na.rm = TRUE) else NA_integer_,
    n_pass_sn        = if ("pass_sn" %in% names(.)) sum(pass_sn, na.rm = TRUE) else NA_integer_,
    n_pass_rt        = if ("pass_rt" %in% names(.)) sum(pass_rt, na.rm = TRUE) else NA_integer_,
    n_keep_for_corr  = if ("keep_for_correction" %in% names(.)) sum(keep_for_correction, na.rm = TRUE) else NA_integer_,
    .groups = "drop"
  )

cat("  [Feature 级 QC 条件统计（按 batch）]\n")
print(feat_qc_by_batch)
cat("\n")

# 写出一个 QA 表方便以后看
qa_feat_qc_file <- file.path(qc_dir, "qa_feature_qc_by_batch.tsv")
readr::write_tsv(feat_qc_by_batch, qa_feat_qc_file)
cat("  [INFO] Feature QC 条件统计已写出: ", qa_feat_qc_file, "\n\n")

## 3.3 keep_for_correction 标记在长表中的一致性检查

cat("  [keep_for_correction 标记一致性检查]\n")

if (!"keep_for_correction" %in% names(dedup_flag)) {
  stop("[dedup_withflag] 中缺少 keep_for_correction 列。", call. = FALSE)
}

# 同一 (batch, lipid) 下，keep_for_correction 是否一致
keep_consistency <- dedup_flag %>%
  group_by(batch, lipidName, LipidIon, CalcMz) %>%
  summarise(
    n_rows = n(),
    n_keep_distinct = n_distinct(keep_for_correction),
    .groups = "drop"
  ) %>%
  filter(n_keep_distinct > 1)

n_inconsistent <- nrow(keep_consistency)
cat("    - (batch, lipid) 下 keep_for_correction 不一致的 feature 数:", n_inconsistent, "\n")
if (n_inconsistent > 0) {
  cat("    [WARN] 前几条不一致的 feature:\n")
  print(keep_consistency %>% head(10))
}
cat("\n")

# 对比 feature_qc 中的 keep_for_correction 与 dedup_flag 中的情况
feat_keep_summary <- feature_qc %>%
  group_by(batch, lipidName, LipidIon, CalcMz) %>%
  summarise(
    keep_feat = any(keep_for_correction %in% TRUE, na.rm = TRUE),
    .groups = "drop"
  )

keep_from_long <- dedup_flag %>%
  group_by(batch, lipidName, LipidIon, CalcMz) %>%
  summarise(
    keep_long = any(keep_for_correction %in% TRUE, na.rm = TRUE),
    .groups = "drop"
  )

join_keep <- feat_keep_summary %>%
  full_join(keep_from_long,
            by = c("batch", "lipidName", "LipidIon", "CalcMz"))

n_keep_feat <- sum(join_keep$keep_feat, na.rm = TRUE)
n_keep_long <- sum(join_keep$keep_long, na.rm = TRUE)
n_keep_disagree <- join_keep %>%
  filter((keep_feat %in% TRUE) != (keep_long %in% TRUE)) %>%
  nrow()

cat("    - feature_qc 中 keep_for_correction == TRUE 的 feature 数:", n_keep_feat, "\n")
cat("    - dedup_withflag 中 keep_for_correction == TRUE 的 feature 数:", n_keep_long, "\n")
cat("    - 两者 keep 标记不一致的 feature 数:", n_keep_disagree, "\n\n")

## 3.4 for_batch_correction 子集一致性检查

cat("  [for_batch_correction 子集检查]\n")

# 在 dedup_flag 中，所有 keep_for_correction == TRUE 的行
long_keep <- dedup_flag %>%
  filter(keep_for_correction %in% TRUE)

cat("    - dedup_withflag 中 keep_for_correction == TRUE 的行数:", nrow(long_keep), "\n")
cat("    - for_batch_correction 中行数:", nrow(for_corr), "\n")

# 检查 for_corr 是否是 long_keep 的子集
# 用所有列做 anti_join 看差异
if (nrow(for_corr) > 0) {
  anti1 <- for_corr %>%
    anti_join(long_keep, by = intersect(names(for_corr), names(long_keep)))
  anti2 <- long_keep %>%
    anti_join(for_corr, by = intersect(names(for_corr), names(long_keep)))

  cat("    - for_corr 中不能在 long_keep 中找到的行数:", nrow(anti1), "\n")
  cat("    - long_keep 中不能在 for_corr 中找到的行数:", nrow(anti2), "\n")

  if (nrow(anti1) > 0) {
    cat("    [WARN] for_corr 中有行不在 long_keep 中，前几条:\n")
    print(anti1 %>% head(5))
  }
  if (nrow(anti2) > 0) {
    cat("    [WARN] long_keep 中有行不在 for_corr 中（说明有漏），前几条:\n")
    print(anti2 %>% head(5))
  }
} else {
  cat("    [WARN] for_corr 表为空，无法检查子集关系。\n")
}

cat("\n")

cat("============================================================\n")
cat("[DONE] 05_qc_check_dedup_and_filter.R QA 流程完成。\n")
cat("  建议：结合终端输出和 qa_feature_qc_by_batch.tsv，人工评估阈值是否需要调整。\n")
cat("============================================================\n")