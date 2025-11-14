#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("用法: Rscript 04_dedup_and_qc_filter.R lipid_long_all_batches.tsv [out_dir]", call. = FALSE)
}

infile <- args[1]
out_dir <- ifelse(length(args) >= 2, args[2], "results/lipid")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_tables_dir <- file.path(out_dir, "tables")
out_qc_dir     <- file.path(out_dir, "qc")
dir.create(out_tables_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(out_qc_dir, showWarnings = FALSE, recursive = TRUE)

cat("==== 读取长表 ====\n")
df <- readr::read_tsv(infile, show_col_types = FALSE)
cat("行数:", nrow(df), "列数:", ncol(df), "\n\n")

##--------------------------------------------------
## 0. 基础检查 & 类型清洗
##--------------------------------------------------

required_cols <- c("batch", "sample_id", "intensity",
                   "lipidName", "LipidIon", "CalcMz", "is_qc")

missing_req <- setdiff(required_cols, names(df))
if (length(missing_req) > 0) {
  stop("缺少必要列: ", paste(missing_req, collapse = ", "), call. = FALSE)
}

# 转 numeric
numeric_cols <- intersect(c("intensity", "qc_sn", "qc_rt", "qc_height", "CalcMz"),
                          names(df))

df <- df %>%
  mutate(across(all_of(numeric_cols), as.numeric))

cat("[INFO] 数值列转为 numeric: ",
    paste(numeric_cols, collapse = ", "), "\n\n")

## lipidName == \"Group\" 理论上已经没有，这里再保险过滤一次
df <- df %>% filter(is.na(lipidName) | lipidName != "Group")

##--------------------------------------------------
## 1. Dedup: 检查 (batch, sample_id, lipidName, LipidIon, CalcMz) 键是否有重复
##--------------------------------------------------

key_cols <- c("batch", "sample_id", "lipidName", "LipidIon", "CalcMz")

cat("==== Dedup 检查 ====\n")

# 每个键的行数 & 数值多样性
dup_groups <- df %>%
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

n_dup_total    <- dup_groups %>% filter(is_dup) %>% nrow()
n_dup_equal    <- dup_groups %>% filter(all_equal) %>% nrow()
n_dup_conflict <- dup_groups %>% filter(is_dup & !all_equal) %>% nrow()

cat("  有重复 key 组数:", n_dup_total, "\n")
cat("    其中完全相同(可安全合并)组数:", n_dup_equal, "\n")
cat("    其中存在数值冲突(需要择优/人工核查)组数:", n_dup_conflict, "\n\n")

# 提取键集合
keys_dup_all_equal <- dup_groups %>%
  filter(all_equal) %>%
  select(all_of(key_cols))

keys_dup_conflict <- dup_groups %>%
  filter(is_dup & !all_equal) %>%
  select(all_of(key_cols))

keys_is_dup <- dup_groups %>%
  filter(is_dup) %>%
  select(all_of(key_cols))

# 1) 非重复行：直接保留
df_notdup <- df %>%
  anti_join(keys_is_dup, by = key_cols)

# 2) 重复但完全相同：每组取第一行
df_dup_equal <- df %>%
  inner_join(keys_dup_all_equal, by = key_cols) %>%
  group_by(across(all_of(key_cols))) %>%
  slice_head(n = 1) %>%
  ungroup()

# 3) 重复且存在数值冲突：全部导出到 conflict 文件 + 选一行“最佳峰”
df_dup_conflict_all <- df %>%
  inner_join(keys_dup_conflict, by = key_cols)

if (nrow(df_dup_conflict_all) > 0) {
  conflict_file <- file.path(out_qc_dir, "dedup_conflicts.tsv")
  readr::write_tsv(df_dup_conflict_all, conflict_file)
  cat("[WARN] 发现存在数值不一致的重复行，已输出到:\n  ",
      conflict_file, "\n")
}

# 选择策略：优先 qc_sn 高，其次 intensity 高
has_qc_sn <- "qc_sn" %in% names(df_dup_conflict_all)

df_dup_conflict_resolved <- df_dup_conflict_all %>%
  mutate(
    .score_sn  = if (has_qc_sn) qc_sn else NA_real_,
    .score_int = intensity
  ) %>%
  group_by(across(all_of(key_cols))) %>%
  arrange(
    desc(coalesce(.score_sn, -Inf)),
    desc(coalesce(.score_int, -Inf))
  ) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  select(-.score_sn, -.score_int)

# 合并得到 Dedup 后长表
df_dedup <- bind_rows(
  df_notdup,
  df_dup_equal,
  df_dup_conflict_resolved
) %>%
  arrange(batch, lipidName, LipidIon, CalcMz, sample_id)

cat("[INFO] Dedup 后行数:", nrow(df_dedup), "\n\n")

dedup_file <- file.path(out_tables_dir, "lipid_long_dedup.tsv")
readr::write_tsv(df_dedup, dedup_file)
cat("[INFO] Dedup 后长表已写出: ", dedup_file, "\n\n")

##--------------------------------------------------
## 2. 基于 QC 的 feature 级 QC 统计
##--------------------------------------------------

cat("==== 基于 QC 的 feature 级 QC 统计 ====\n")

if (!"is_qc" %in% names(df_dedup)) {
  stop("数据中没有 is_qc 列，无法进行 QC 统计。", call. = FALSE)
}

has_qc_sn     <- "qc_sn"     %in% names(df_dedup)
has_qc_rt     <- "qc_rt"     %in% names(df_dedup)
has_qc_height <- "qc_height" %in% names(df_dedup)

# 参数：阈值（之后你可以根据数据分布调）
min_n_qc_non_na   <- 2L
min_qc_mean       <- 1e5      # QC 平均强度太低的剔除
max_qc_cv         <- 0.3      # CV > 30% 视为不稳定
min_qc_sn_median  <- 10       # 中位 S/N 太低视为噪声
max_qc_rt_sd      <- 0.05     # RT 波动太大的剔除（单位视你的 RT 而定）

cat("  阈值设置:\n")
cat("    - QC 非 NA 个数 ≥", min_n_qc_non_na, "\n")
cat("    - QC 平均强度 ≥", format(min_qc_mean, scientific = TRUE), "\n")
cat("    - QC CV ≤", max_qc_cv, "\n")
if (has_qc_sn) cat("    - QC S/N 中位数 ≥", min_qc_sn_median, "\n")
if (has_qc_rt) cat("    - QC RT SD ≤", max_qc_rt_sd, "\n")
cat("\n")

# 2.1 只看 QC 行，按 (batch, lipid) 汇总
qc_summary <- df_dedup %>%
  filter(is_qc) %>%
  group_by(batch, lipidName, LipidIon, CalcMz) %>%
  summarise(
    n_qc_total   = n(),
    n_qc_non_na  = sum(!is.na(intensity)),
    qc_mean      = ifelse(n_qc_non_na > 0, mean(intensity, na.rm = TRUE), NA_real_),
    qc_sd        = ifelse(n_qc_non_na > 1,  sd(intensity, na.rm = TRUE), NA_real_),
    qc_cv        = ifelse(!is.na(qc_mean) & qc_mean > 0 & !is.na(qc_sd),
                          qc_sd / qc_mean, NA_real_),
    qc_missing_prop = ifelse(n_qc_total > 0,
                             1 - n_qc_non_na / n_qc_total,
                             NA_real_),
    qc_sn_median = if (has_qc_sn) stats::median(qc_sn, na.rm = TRUE) else NA_real_,
    qc_rt_mean   = if (has_qc_rt) mean(qc_rt, na.rm = TRUE) else NA_real_,
    qc_rt_sd     = if (has_qc_rt) stats::sd(qc_rt, na.rm = TRUE) else NA_real_,
    qc_height_mean = if (has_qc_height) mean(qc_height, na.rm = TRUE) else NA_real_,
    .groups = "drop"
  )

cat("  有 QC 记录的 (batch, lipid) 个数:", nrow(qc_summary), "\n")

# 2.2 生物样本检测率统计（含 QC 与否都算）
det_summary <- df_dedup %>%
  group_by(batch, lipidName, LipidIon, CalcMz) %>%
  summarise(
    n_samples_total   = n(),
    n_samples_non_na  = sum(!is.na(intensity)),
    detect_rate       = n_samples_non_na / n_samples_total,
    .groups = "drop"
  )

# 2.3 合并 QC summary + detect summary
feature_summary <- qc_summary %>%
  full_join(det_summary,
            by = c("batch", "lipidName", "LipidIon", "CalcMz"))

# 对没有 QC 的 feature，后面 keep_for_correction 一定会是 FALSE
feature_summary <- feature_summary %>%
  mutate(
    # 各单项通过与否
    pass_n_qc = !is.na(n_qc_non_na) & n_qc_non_na >= min_n_qc_non_na,
    pass_mean = !is.na(qc_mean)     & qc_mean     >= min_qc_mean,
    pass_cv   = !is.na(qc_cv)       & qc_cv       <= max_qc_cv,
    pass_sn   = if (has_qc_sn) (!is.na(qc_sn_median) & qc_sn_median >= min_qc_sn_median) else TRUE,
    pass_rt   = if (has_qc_rt) (is.na(qc_rt_sd) | qc_rt_sd <= max_qc_rt_sd) else TRUE,
    # 最终是否用于“批次内矫正”这一步
    keep_for_correction = pass_n_qc & pass_mean & pass_cv & pass_sn & pass_rt
  )

qc_summary_file <- file.path(out_qc_dir, "lipid_feature_qc_summary.tsv")
readr::write_tsv(feature_summary, qc_summary_file)
cat("[INFO] feature 级 QC summary 已写出: ", qc_summary_file, "\n\n")

cat("  其中 keep_for_correction == TRUE 的 feature 数量: ",
    sum(feature_summary$keep_for_correction, na.rm = TRUE), "\n\n")

##--------------------------------------------------
## 3. 把 keep_for_correction 标记回长表 & 输出子集
##--------------------------------------------------

cat("==== 回写标记并输出长表 ====\n")

df_annot <- df_dedup %>%
  left_join(
    feature_summary %>%
      select(batch, lipidName, LipidIon, CalcMz, keep_for_correction),
    by = c("batch", "lipidName", "LipidIon", "CalcMz")
  )

annot_file <- file.path(out_tables_dir, "lipid_long_dedup_with_qcflag.tsv")
readr::write_tsv(df_annot, annot_file)
cat("[INFO] 含 keep_for_correction 标记的长表已写出: ", annot_file, "\n")

df_for_corr <- df_annot %>% filter(keep_for_correction %in% TRUE)

corr_file <- file.path(out_tables_dir, "lipid_long_for_batch_correction.tsv")
readr::write_tsv(df_for_corr, corr_file)
cat("[INFO] 用于“批次内矫正”的子集长表已写出: ", corr_file, "\n\n")

cat("==== DONE: Dedup + QC filtering 完成 ====\n")