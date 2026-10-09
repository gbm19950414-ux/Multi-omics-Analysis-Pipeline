#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
})

msg <- function(...) {
  cat("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", sprintf(...), "\n", sep = "")
}

main <- function(args) {
  if (length(args) < 2) {
    stop("用法: Rscript 01_mergeResult_to_batch1_like.R mergeResult.txt output.csv", call. = FALSE)
  }

  in_file  <- args[1]
  out_file <- args[2]

  if (!file.exists(in_file)) {
    stop(sprintf("输入文件不存在: %s", in_file), call. = FALSE)
  }

  msg("[INFO] 读取 mergeResult: %s", in_file)

  ## 1) 解析文件头，得到 [c-1] ~ [c-9] 对应的 raw 文件名 --------------------
  header_lines <- readLines(in_file, n = 50)
  c_lines <- grep("^#\\[c-[0-9]+\\]:", header_lines, value = TRUE)

  if (length(c_lines) == 0) {
    stop("在文件头部未找到 '#[c-x]:xxx.raw' 行，检查是否为正确的 mergeResult 文件。", call. = FALSE)
  }

  mat <- stringr::str_match(c_lines, "^#\\[c-(\\d+)\\]:(.+)$")
  c_idx     <- as.integer(mat[, 2])
  raw_files <- trimws(mat[, 3])
  raw_base  <- sub("\\.raw$", "", basename(raw_files))   # 去掉路径和 .raw

  # 按 c-idx 排序
  ord      <- order(c_idx)
  c_idx    <- c_idx[ord]
  raw_base <- raw_base[ord]

  msg("[INFO] c-idx 与 raw 文件对应关系:")
  print(data.frame(c_idx = c_idx, raw = raw_base))

  ## 2) 读取主表，跳过注释行 ----------------------------------------------
  df <- readr::read_tsv(
    in_file,
    comment = "#",
    show_col_types = FALSE
  )
  msg("[INFO] 主表维度: %d 行 × %d 列", nrow(df), ncol(df))

  # 可选：只保留 Rej. == 0 的峰
  if ("Rej." %in% names(df)) {
    msg("[INFO] 发现 Rej. 列，仅保留 Rej. == 0 的行")
    df <- df %>% dplyr::filter(`Rej.` == 0)
    msg("[INFO] 过滤后维度: %d 行 × %d 列", nrow(df), ncol(df))
  }

  ## 3) 检查必要列 ---------------------------------------------------------
  needed_cols <- c("LipidIon", "CalcMz")
  miss <- setdiff(needed_cols, names(df))
  if (length(miss) > 0) {
    stop(sprintf("缺少必要列: %s", paste(miss, collapse = ", ")), call. = FALSE)
  }

  ## 4) 构建 lipidName / LipidIon / CalcMz / class ------------------------
  # lipidName: 去掉 "+..."，得到例如 AEA(12:0)
  lipidName <- stringr::str_replace(df$LipidIon, "\\+.*$", "")

  # class 一般来自 Class 列，如果没有就用 NA
  if ("Class" %in% names(df)) {
    class_vec <- df$Class
  } else {
    class_vec <- NA_character_
  }

  ## 5) 提取 Area[c-1..] → 样本列（列名用 raw_base，模仿文件2） ----------
  area_cols_all <- names(df)[startsWith(names(df), "Area[c-")]
  if (length(area_cols_all) == 0) {
    stop("未发现 Area[c-x] 列，检查 mergeResult 文件结构。", call. = FALSE)
  }

  area_idx  <- as.integer(stringr::str_match(area_cols_all, "Area\\[c-(\\d+)\\]")[, 2])
  area_cols <- area_cols_all[order(area_idx)]

  if (length(area_cols) != length(raw_base)) {
    msg("[WARN] Area[c-x] 列数量 (%d) 与 header 中 c-x 数量 (%d) 不一致",
        length(area_cols), length(raw_base))
  }

  area_mat <- df %>% dplyr::select(all_of(area_cols))
  colnames(area_mat) <- raw_base  # 用原始 raw 基名作为列名，和 20240617_0612_HH_cxy_data.csv 一致

  ## 6) 提取 Grade[c-x] → Grade.c.1. ... ---------------------------------
  grade_cols_all <- names(df)[startsWith(names(df), "Grade[c-")]
  grade_mat <- NULL
  if (length(grade_cols_all) > 0) {
    grade_idx  <- as.integer(stringr::str_match(grade_cols_all, "Grade\\[c-(\\d+)\\]")[, 2])
    grade_cols <- grade_cols_all[order(grade_idx)]
    grade_mat  <- df %>% dplyr::select(all_of(grade_cols))
    colnames(grade_mat) <- paste0("Grade.c.", seq_along(grade_cols), ".")
  }

  ## 7) 提取 S/N[c-x] → S.N.c.1. ... --------------------------------------
  sn_cols_all <- names(df)[startsWith(names(df), "S/N[c-")]
  sn_mat <- NULL
  if (length(sn_cols_all) > 0) {
    sn_idx  <- as.integer(stringr::str_match(sn_cols_all, "S/N\\[c-(\\d+)\\]")[, 2])
    sn_cols <- sn_cols_all[order(sn_idx)]
    sn_mat  <- df %>% dplyr::select(all_of(sn_cols))
    colnames(sn_mat) <- paste0("S.N.c.", seq_along(sn_cols), ".")
  }

  ## 8) 组装宽表（不含 Group 行） -----------------------------------------
  out <- tibble(
    lipidName = lipidName,
    LipidIon  = df$LipidIon,
    CalcMz    = df$CalcMz,
    class     = class_vec
  ) %>%
    bind_cols(area_mat)

  if (!is.null(grade_mat)) out <- bind_cols(out, grade_mat)
  if (!is.null(sn_mat))    out <- bind_cols(out, sn_mat)

  ## 8.1) 确保样本列类型为字符，避免与 Group 行类型冲突 -------------------
  sample_cols <- raw_base[raw_base %in% names(out)]
  if (length(sample_cols) > 0) {
    out <- out %>%
      mutate(across(all_of(sample_cols), as.character))
  }

  ## 9) 构建 Group 行（第一行），模仿文件2 -------------------------------
  # 原文件2中 Group 行的结构：
  # lipidName = "Group", LipidIon / CalcMz / class 为空；
  # 每个样本列填 S 或 Q，这里我们根据 raw_base 是否包含 "QC" 来判断。
  group_row <- out[1, , drop = FALSE]
  group_row[1, ] <- NA

  group_row$lipidName <- "Group"

  # 对样本列填 S / Q
  for (nm in raw_base) {
    if (!nm %in% names(group_row)) next
    if (grepl("QC", nm, ignore.case = TRUE)) {
      group_row[[nm]] <- "Q"
    } else {
      group_row[[nm]] <- "S"
    }
  }
  # Grade.c.* / S.N.c.* 在 Group 行保持 NA

  out_final <- bind_rows(group_row, out)

  msg("[INFO] 输出到: %s", out_file)
  readr::write_csv(out_final, out_file)

  msg("[DONE] mergeResult → 文件2风格宽表 完成。")
}

if (identical(environment(), globalenv())) {
  args <- commandArgs(trailingOnly = TRUE)
  main(args)
}