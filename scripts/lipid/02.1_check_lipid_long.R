#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("用法: Rscript 03_check_lipid_long.R <lipid_long_all_batches.tsv>", call. = FALSE)
}
in_file <- args[1]

cat("============================================================\n")
cat("[INFO] 检查文件:", in_file, "\n")
cat("============================================================\n\n")

# ---------- 1. 读取 ----------
cat("[STEP1] 读取数据...\n")
dat <- read_tsv(in_file, col_types = cols(), progress = FALSE)
cat("  - 行数:", nrow(dat), "\n")
cat("  - 列数:", ncol(dat), "\n\n")

# ---------- 2. 列信息 ----------
cat("[STEP2] 列名与类型概览:\n")
print(tibble(
  col = names(dat),
  class = vapply(dat, function(x) class(x)[1], character(1))
))
cat("\n")

# 尝试识别关键列名（尽量容错）
get_first <- function(cands) {
  cands[cands %in% names(dat)][1]
}

col_batch     <- get_first(c("batch", "Batch"))
col_sample_id <- get_first(c("sample_id", "Sample", "sample"))
col_lipid     <- get_first(c("lipidName", "Lipid", "feature"))
col_value     <- get_first(c("value", "intensity", "Intensity", "area", "Area"))
col_is_qc     <- get_first(c("is_qc", "qc_flag"))
col_qc_grade  <- get_first(c("qc_grade"))
col_qc_sn     <- get_first(c("qc_sn"))
col_qc_rt     <- get_first(c("qc_rt"))
col_qc_h      <- get_first(c("qc_height", "qc_h", "height"))

cat("[INFO] 识别出的关键信息列:\n")
cat("  batch:     ", col_batch,    "\n")
cat("  sample_id: ", col_sample_id,"\n")
cat("  lipidName: ", col_lipid,    "\n")
cat("  value:     ", col_value,    "\n")
cat("  is_qc:     ", ifelse(is.na(col_is_qc),  "NONE", col_is_qc), "\n")
cat("  qc_grade:  ", ifelse(is.na(col_qc_grade),"NONE", col_qc_grade), "\n")
cat("  qc_sn:     ", ifelse(is.na(col_qc_sn),   "NONE", col_qc_sn), "\n")
cat("  qc_rt:     ", ifelse(is.na(col_qc_rt),   "NONE", col_qc_rt), "\n")
cat("  qc_height: ", ifelse(is.na(col_qc_h),    "NONE", col_qc_h), "\n\n")

# ---------- 3. 基本完整性：缺失值 ----------
cat("[STEP3] 每列缺失值统计:\n")
na_stat <- colSums(is.na(dat))
na_tab <- tibble(
  col = names(na_stat),
  n_na = as.integer(na_stat),
  prop_na = round(na_stat / nrow(dat), 4)
) %>% arrange(desc(prop_na))
print(na_tab)
cat("\n")

# ---------- 4. batch / sample / lipid 结构 ----------
if (!is.na(col_batch)) {
  cat("[STEP4] batch 分布:\n")
  print(dat %>% count(.data[[col_batch]]) %>% arrange(.data[[col_batch]]))
  cat("\n")
}

if (!is.na(col_sample_id)) {
  cat("[STEP5] 每个 batch 内各 sample 数量:\n")
  if (!is.na(col_batch)) {
    print(
      dat %>%
        distinct(.data[[col_batch]], .data[[col_sample_id]]) %>%
        count(.data[[col_batch]]) %>%
        arrange(.data[[col_batch]])
    )
  } else {
    print(dat %>% distinct(.data[[col_sample_id]]) %>% count())
  }
  cat("\n")
}

if (!is.na(col_lipid)) {
  cat("[STEP6] 总 lipid 数量 (去重):\n")
  print(dat %>% distinct(.data[[col_lipid]]) %>% count())
  cat("\n")
}

# ---------- 5. 键唯一性检查 ----------
cat("[STEP7] (batch, sample_id, lipid) 是否唯一？\n")
if (!is.na(col_batch) && !is.na(col_sample_id) && !is.na(col_lipid)) {
  dup <- dat %>%
    count(.data[[col_batch]], .data[[col_sample_id]], .data[[col_lipid]]) %>%
    filter(n > 1)
  if (nrow(dup) == 0) {
    cat("  [OK] 没有重复键 (batch, sample_id, lipid)\n\n")
  } else {
    cat("  [WARN] 存在重复键, 前几条如下:\n")
    print(head(dup, 20))
    cat("\n")
  }
} else {
  cat("  [SKIP] 由于缺少 batch / sample_id / lipidName 无法检查键唯一性\n\n")
}

# ---------- 6. value 分布 ----------
cat("[STEP8] 强度列 value 分布检查:\n")
if (!is.na(col_value)) {
  v <- dat[[col_value]]
  cat("  - 非 NA 个数:", sum(!is.na(v)), "\n")
  cat("  - NA 个数:   ", sum(is.na(v)), "\n")
  cat("  - 最小值:    ", min(v, na.rm = TRUE), "\n")
  cat("  - 最大值:    ", max(v, na.rm = TRUE), "\n")
  cat("  - 零值个数:  ", sum(v == 0, na.rm = TRUE), "\n\n")
} else {
  cat("  [SKIP] 未识别到 value 列，无法检查强度分布\n\n")
}

# ---------- 7. is_qc 与 sample_id 一致性 ----------
cat("[STEP9] QC 标记检查:\n")
if (!is.na(col_is_qc) && !is.na(col_sample_id)) {
  cat("  - is_qc 分布:\n")
  print(table(dat[[col_is_qc]], useNA = "ifany"))
  cat("\n")

  # 用 sample_id 中是否含 QC 做一个参考
  has_qc_pattern <- str_detect(dat[[col_sample_id]], "(?i)QC")
  tab_cmp <- table(
    is_qc_col  = dat[[col_is_qc]],
    sample_has_QC_in_name = has_qc_pattern,
    useNA = "ifany"
  )
  cat("  - is_qc 列 vs sample_id 名称中是否含“QC”的交叉表:\n")
  print(tab_cmp)
  cat("\n")
} else if (!is.na(col_sample_id)) {
  cat("  [WARN] 没有 is_qc 列，只能用 sample_id 名称中含 'QC' 作为 QC 样本判定\n")
  print(table(str_detect(dat[[col_sample_id]], "(?i)QC"), useNA = "ifany"))
  cat("\n")
} else {
  cat("  [SKIP] 无法检查 QC 标记（缺少 sample_id）\n\n")
}

# ---------- 8. QC 附属信息完整性 ----------
cat("[STEP10] QC 附属列完整性（仅对 QC 样本统计）:\n")
if (!is.na(col_is_qc)) {
  dat_qc <- dat %>% filter(.data[[col_is_qc]] %in% TRUE)
} else if (!is.na(col_sample_id)) {
  dat_qc <- dat %>% filter(str_detect(.data[[col_sample_id]], "(?i)QC"))
} else {
  dat_qc <- NULL
}

if (!is.null(dat_qc) && nrow(dat_qc) > 0) {
  cat("  - QC 行数:", nrow(dat_qc), "\n")
  if (!is.na(col_qc_grade)) {
    cat("  - qc_grade 非 NA 占比:", 1 - mean(is.na(dat_qc[[col_qc_grade]])), "\n")
  }
  if (!is.na(col_qc_sn)) {
    cat("  - qc_sn 非 NA 占比:   ", 1 - mean(is.na(dat_qc[[col_qc_sn]])), "\n")
  }
  if (!is.na(col_qc_rt)) {
    cat("  - qc_rt 非 NA 占比:   ", 1 - mean(is.na(dat_qc[[col_qc_rt]])), "\n")
  }
  if (!is.na(col_qc_h)) {
    cat("  - qc_height 非 NA 占比:", 1 - mean(is.na(dat_qc[[col_qc_h]])), "\n")
  }
  cat("\n")
} else {
  cat("  [WARN] 没有检测到 QC 行（按 is_qc 或 sample_id 模式）。\n\n")
}

# ---------- 9. 每个 batch × sample 的缺失比例 ----------
cat("[STEP11] 每个 (batch, sample_id) 的检测率（非 NA 占比）:\n")
if (!is.na(col_batch) && !is.na(col_sample_id) && !is.na(col_value)) {
  det_rate <- dat %>%
    group_by(.data[[col_batch]], .data[[col_sample_id]]) %>%
    summarise(
      n = n(),
      n_non_na = sum(!is.na(.data[[col_value]])),
      detect_rate = n_non_na / n,
      .groups = "drop"
    ) %>%
    arrange(.data[[col_batch]], desc(detect_rate))
  print(det_rate)
  cat("\n")
} else {
  cat("  [SKIP] 无法计算检测率（缺少 batch / sample_id / value 其中之一）\n\n")
}

cat("============================================================\n")
cat("[DONE] 质量与完整性检查完成。\n")
cat("============================================================\n")