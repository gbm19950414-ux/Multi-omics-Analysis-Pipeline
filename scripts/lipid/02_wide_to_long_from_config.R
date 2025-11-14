#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(yaml)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(tibble)
})

# 小工具：带时间戳的 log -----------------------------------------------
msg <- function(...) {
  cat("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", sprintf(...), "\n", sep = "")
}

# 小工具：解析列名（支持大小写轻微不一致） ------------------------------
resolve_col <- function(df, name, label = name) {
  if (!is.null(name) && name %in% names(df)) {
    return(name)
  }
  # 大小写不敏感匹配
  idx <- which(tolower(names(df)) == tolower(name))
  if (length(idx) > 0) {
    real_name <- names(df)[idx[1]]
    msg("  [WARN] 列 '%s' 未找到，自动使用 '%s'", name, real_name)
    return(real_name)
  }
  stop(sprintf("FATAL: 无法在数据中找到列 '%s'（label=%s）", name, label), call. = FALSE)
}

# 识别样本列（非元信息，非 QC 前缀） ------------------------------------
guess_sample_cols <- function(df, meta_cols, qc_prefixes) {
  all_cols <- names(df)
  qc_prefixes <- qc_prefixes[!is.na(qc_prefixes) & qc_prefixes != ""]
  
  is_qc_col <- function(col) {
    any(startsWith(col, qc_prefixes))
  }
  
  sample_cols <- all_cols[!(all_cols %in% meta_cols) & !vapply(all_cols, is_qc_col, logical(1))]
  sample_cols
}

# 处理单个 batch：wide → long --------------------------------------------
process_one_batch <- function(batch_cfg, global_cfg) {
  project_dir <- global_cfg$project_dir %||% getwd()
  file_rel    <- batch_cfg$file
  file_path   <- file.path(project_dir, file_rel)
  
  msg("[INFO] reading [%s] %s", batch_cfg$name, file_path)
  
  if (!file.exists(file_path)) {
    stop(sprintf("FATAL: 文件不存在: %s", file_path), call. = FALSE)
  }
  
  # 读表（自动类型）
  df_raw <- readr::read_csv(file_path, show_col_types = FALSE)
  
  col_cfg   <- global_cfg$columns
  # Remove Group pseudo-row early
  df_raw <- df_raw[df_raw[[col_cfg$lipid]] != "Group", , drop = FALSE]
  
  # 解析固定列名
  lipid_col <- resolve_col(df_raw, col_cfg$lipid,  "lipid")
  ion_col   <- resolve_col(df_raw, col_cfg$ion,    "ion")
  mz_col    <- resolve_col(df_raw, col_cfg$mz,     "mz")
  # class 列可能是 Class / class，交给 resolve_col 处理
  class_col <- resolve_col(df_raw, col_cfg$class %||% "Class", "class")
  
  # Group 列（第二行 group 信息）
  has_group_col <- "Group" %in% names(df_raw)
  if (has_group_col) {
    group_col <- "Group"
  } else {
    group_col <- NULL
  }
  
  # meta 列：脂质信息 + Group（如果存在）
  meta_cols <- c(lipid_col, ion_col, class_col, mz_col)
  if (!is.null(group_col)) {
    meta_cols <- union(meta_cols, group_col)
  }
  
  # QC 前缀（目前只用于从样本列中排除）
  qc_cols_cfg <- batch_cfg$qc_cols %||% list()
  qc_prefixes <- unlist(qc_cols_cfg)
  qc_prefixes <- qc_prefixes[!is.na(qc_prefixes) & qc_prefixes != ""]
  qc_prefixes <- unique(c(qc_prefixes, "Grade.c.", "S.N.c.", "Rt.c.", "Height.c."))
  
  # 识别样本列
  sample_cols <- guess_sample_cols(df_raw, meta_cols, qc_prefixes)
  msg("  [INFO] 识别到样本列 %d 个：%s",
      length(sample_cols),
      paste(sample_cols, collapse = ", "))
  
  # # 提取 group 行（如果有）
  # if (!is.null(group_col)) {
  #   grp_idx <- which(df_raw[[lipid_col]] == "Group")
  #   if (length(grp_idx) > 0) {
  #     group_row <- df_raw[grp_idx[1], , drop = FALSE]
  #     df <- df_raw[-grp_idx, , drop = FALSE]
  #     
  #     # 为每个样本列提取 S/Q 等类型信息
  #     sample_type_map <- as.list(group_row[1, sample_cols])
  #     sample_type_map <- vapply(sample_type_map, function(x) {
  #       x <- as.character(x)
  #       ifelse(is.na(x) | x == "", NA_character_, x)
  #     }, character(1))
  #   } else {
  #     df <- df_raw
  #     sample_type_map <- NULL
  #   }
  # } else {
  #   df <- df_raw
  #   sample_type_map <- NULL
  # }
  df <- df_raw
  sample_type_map <- NULL

  # 为每一行添加内部行号，后续 QC 值按 (.row_id, 样本列序号) 匹配
  df <- df %>% mutate(.row_id = dplyr::row_number())
  
  # wide → long：一个 batch 的长表
  long_df <- df %>%
    select(all_of(c(lipid_col, ion_col, class_col, mz_col)), .row_id, all_of(sample_cols)) %>%
    pivot_longer(
      cols      = all_of(sample_cols),
      names_to  = "sample_id",
      values_to = "intensity"
    ) %>%
    mutate(
      batch      = batch_cfg$name,
      lipidName  = .data[[lipid_col]],
      LipidIon   = .data[[ion_col]],
      Class      = .data[[class_col]],
      CalcMz     = .data[[mz_col]]
    ) %>%
    select(batch, sample_id, intensity, lipidName, LipidIon, Class, CalcMz, .row_id, everything())

  # ===== 提取 QC 列（按列顺序匹配到 sample_id，返回的是每行数值而非列名） =====
  qc_cols_all <- names(df_raw)
  qc_grade_cols  <- qc_cols_all[startsWith(qc_cols_all, "Grade.c.")]
  qc_sn_cols     <- qc_cols_all[startsWith(qc_cols_all, "S.N.c.")]
  qc_rt_cols     <- qc_cols_all[startsWith(qc_cols_all, "Rt.c.")]
  qc_height_cols <- qc_cols_all[startsWith(qc_cols_all, "Height.c.")]

  # 构建一个小工具：把某一种 QC 矩阵转成长表 (.row_id, sample_id, value)
  build_qc_long <- function(df_wide, qc_cols, sample_cols, value_name) {
    qc_cols <- qc_cols[!is.na(qc_cols)]
    if (length(qc_cols) == 0) {
      return(NULL)
    }
    # 按列顺序建立 qc_col ↔ sample_id 的映射
    map_tbl <- tibble(
      qc_col    = qc_cols,
      sample_id = sample_cols[seq_along(qc_cols)]
    )
    df_wide %>%
      select(.row_id, all_of(qc_cols)) %>%
      pivot_longer(
        cols      = - .row_id,
        names_to  = "qc_col",
        values_to = value_name
      ) %>%
      left_join(map_tbl, by = "qc_col") %>%
      filter(!is.na(sample_id)) %>%
      select(.row_id, sample_id, !!value_name)
  }

  n_samp <- length(sample_cols)
  qc_grade_long  <- build_qc_long(df, qc_grade_cols,  sample_cols, "qc_grade")
  qc_sn_long     <- build_qc_long(df, qc_sn_cols,     sample_cols, "qc_sn")
  qc_rt_long     <- build_qc_long(df, qc_rt_cols,     sample_cols, "qc_rt")
  qc_height_long <- build_qc_long(df, qc_height_cols, sample_cols, "qc_height")

  # 将各类 QC 信息按 (.row_id, sample_id) 合并回长表
  if (!is.null(qc_grade_long)) {
    long_df <- long_df %>%
      left_join(qc_grade_long, by = c(".row_id", "sample_id"))
  } else {
    long_df <- long_df %>% mutate(qc_grade = NA_real_)
  }

  if (!is.null(qc_sn_long)) {
    long_df <- long_df %>%
      left_join(qc_sn_long, by = c(".row_id", "sample_id"))
  } else {
    long_df <- long_df %>% mutate(qc_sn = NA_real_)
  }

  if (!is.null(qc_rt_long)) {
    long_df <- long_df %>%
      left_join(qc_rt_long, by = c(".row_id", "sample_id"))
  } else {
    long_df <- long_df %>% mutate(qc_rt = NA_real_)
  }

  if (!is.null(qc_height_long)) {
    long_df <- long_df %>%
      left_join(qc_height_long, by = c(".row_id", "sample_id"))
  } else {
    long_df <- long_df %>% mutate(qc_height = NA_real_)
  }
  
  # 补充 sample_type（来自 Group 行）
  if (!is.null(sample_type_map)) {
    long_df <- long_df %>%
      mutate(sample_type = unname(sample_type_map[sample_id]))
  } else {
    long_df <- long_df %>%
      mutate(sample_type = NA_character_)
  }
  
  # 标记是否 QC（利用全局 patterns.qc_regex）
  qc_regex <- global_cfg$patterns$qc_regex %||% "QC"
  long_df <- long_df %>%
    mutate(is_qc = str_detect(sample_id, qc_regex))

  # 移除内部辅助列和无用列：.row_id、原始 class 列、sample_type
  long_df <- long_df %>%
    select(-any_of(c(".row_id", "class", "sample_type")))

  long_df
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

# 主函数 -----------------------------------------------------------------
main <- function(args) {
  if (length(args) < 1) {
    stop("用法: Rscript 02_wide_to_long_from_config.R path/to/00_lipid_config.yaml", call. = FALSE)
  }
  cfg_path <- args[1]
  if (!file.exists(cfg_path)) {
    stop(sprintf("FATAL: 找不到配置文件: %s", cfg_path), call. = FALSE)
  }
  
  msg("[INFO] 读取配置: %s", cfg_path)
  cfg <- yaml::read_yaml(cfg_path)
  
  project_dir <- cfg$project_dir %||% getwd()
  out_dir_rel <- cfg$out_dir %||% "results/lipid/tables"
  out_dir     <- file.path(project_dir, out_dir_rel)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # 需要处理的批次
  all_batch_cfgs <- cfg$input$batches
  include_names  <- cfg$input$include_batches %||% vapply(all_batch_cfgs, `[[`, character(1), "name")
  
  batch_cfgs <- keep(all_batch_cfgs, ~ .x$name %in% include_names)
  msg("[INFO] 将转长表的批次: %s", paste(include_names, collapse = ", "))
  
  # 逐批处理
  long_list <- lapply(batch_cfgs, process_one_batch, global_cfg = cfg)
  long_all  <- bind_rows(long_list)
  
  # 写出长表
  out_file <- file.path(out_dir, "lipid_long_all_batches.tsv")
  msg("[INFO] 写出长表: %s (n=%d 行, %d 列)", out_file, nrow(long_all), ncol(long_all))
  readr::write_tsv(long_all, out_file)
  
  msg("[DONE] 完成 wide→long 转换。")
}

# 运行 --------------------------------------------------------------------
if (identical(environment(), globalenv())) {
  args <- commandArgs(trailingOnly = TRUE)
  main(args)
}