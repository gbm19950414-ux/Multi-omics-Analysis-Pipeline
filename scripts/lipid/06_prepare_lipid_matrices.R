#!/usr/bin/env Rscript

## ============================================================
## 06_prepare_lipid_matrices.R
##
## 从矫正后的长表构建：
##   1) 样本信息表 sample_info_lipid.tsv
##   2) 全脂质矩阵  lipid_matrix_allfat_log2.tsv
##   3) CL 子集矩阵 lipid_matrix_CL_log2.tsv
##
## 使用 00_lipid_downstream_config.yaml 中的配置：
##   - paths: 输入/输出路径
##   - feature_filter: allfat / CL 的过滤规则
##   - universe: allfat / CL 的批次和 class 选择
##   - sample_parsing: group 的解析和 overrides
##
## 约定：
##   表达量列优先级：intensity_pqn_bb > intensity_pqn > intensity
## ============================================================

suppressPackageStartupMessages({
  library(yaml)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

log_msg <- function(...){
  cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), ..., "\n")
}

## ------------------------------------------------------------
## 0. 读配置
## ------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1) args[1] else "scripts/lipid/00_lipid_downstream_config.yaml"

if (!file.exists(config_path)) {
  stop("Config file not found: ", config_path,
       "\nUsage: Rscript 06_prepare_lipid_matrices.R path/to/00_lipid_downstream_config.yaml")
}

log_msg("[INFO] 使用配置文件: ", config_path)
cfg <- yaml::read_yaml(config_path)

paths <- cfg$paths

## ------------------------------------------------------------
## 1. 读取矫正后的长表，确定表达量列
## ------------------------------------------------------------
log_msg("[STEP1] 读取矫正后的长表: ", paths$long_table_pqn_bb)
long <- readr::read_tsv(paths$long_table_pqn_bb, show_col_types = FALSE)
log_msg("[INFO] long nrow = ", nrow(long), " ; ncol = ", ncol(long))

expr_col <- NULL
if ("intensity_pqn_bb" %in% colnames(long)) {
  expr_col <- "intensity_pqn_bb"
} else if ("intensity_pqn" %in% colnames(long)) {
  expr_col <- "intensity_pqn"
} else if ("intensity" %in% colnames(long)) {
  expr_col <- "intensity"
} else {
  stop("无法在长表中找到表达量列（intensity_pqn_bb / intensity_pqn / intensity）")
}
log_msg("[INFO] 使用表达量列: ", expr_col)

long <- long %>%
  filter(!is.na(.data[[expr_col]])) %>%
  mutate(global_sample_id = paste(batch, sample_id, sep = "_"))

## ------------------------------------------------------------
## 2. 构建 sample_info_lipid.tsv
##    - batch: 直接来自 long$batch（假定每个 sample 只属于一个 batch）
##    - group: 根据 sample_parsing$group_regex 从 sample_id 解析
##    - condition: 默认 NA，可通过 overrides 指定
## ------------------------------------------------------------
log_msg("[STEP2] 构建 sample_info 表 ...")

sample_info <- long %>%
  distinct(batch, sample_id, global_sample_id)

## 检查同名 sample_id 是否跨 batch 复用（仅提示，不作为错误）
multi_batch <- sample_info %>%
  count(sample_id) %>%
  filter(n > 1)
if (nrow(multi_batch) > 0) {
  warning("以下 sample_id 出现在多个 batch；已通过 global_sample_id 区分:\n",
          paste(multi_batch$sample_id, collapse = ", "))
}

## group 解析
group_regex <- cfg$sample_parsing$group_regex %||% list()
if (length(group_regex) > 0) {
  sample_info$group <- NA_character_
  for (g in names(group_regex)) {
    pattern <- group_regex[[g]]
    idx <- is.na(sample_info$group) & grepl(pattern, sample_info$sample_id)
    sample_info$group[idx] <- g
  }
} else {
  sample_info$group <- NA_character_
}

## condition 列（可被 overrides 改写）
sample_info$condition <- NA_character_

## overrides（可选）
overrides <- cfg$sample_parsing$overrides %||% list()
if (length(overrides) > 0) {
  log_msg("[INFO] 应用 sample overrides ...")
  for (sid in names(overrides)) {
    ov <- overrides[[sid]]

    # 优先按 global_sample_id 匹配；若不存在，再按 sample_id 匹配
    if (sid %in% sample_info$global_sample_id) {
      idx <- sample_info$global_sample_id == sid
    } else if (sid %in% sample_info$sample_id) {
      idx <- sample_info$sample_id == sid
    } else {
      warning("override 指定的 ID 在数据中不存在 (既不是 global_sample_id 也不是 sample_id): ", sid)
      next
    }

    if (!is.null(ov$batch))     sample_info$batch[idx]     <- ov$batch
    if (!is.null(ov$group))     sample_info$group[idx]     <- ov$group
    if (!is.null(ov$condition)) sample_info$condition[idx] <- ov$condition
  }
}

sample_info_out <- sample_info %>%
  select(global_sample_id, batch, sample_id, dplyr::everything())

dir.create(dirname(paths$sample_info_output), recursive = TRUE, showWarnings = FALSE)
readr::write_tsv(sample_info_out, paths$sample_info_output)
log_msg("[OK] 写出 sample_info: ", paths$sample_info_output)

## 一个方便的小函数：为后续统计准备 sample_info（避免重复 join）
sample_info_min <- sample_info %>%
  select(global_sample_id, sample_id, batch, group)

## ------------------------------------------------------------
## 3. 通用函数：根据 universe + filter 构建矩阵
## ------------------------------------------------------------
build_matrix_for_universe <- function(long,
                                      sample_info_min,
                                      expr_col,
                                      universe_cfg,
                                      filter_cfg,
                                      label_universe = "allfat",
                                      out_path_matrix) {

  log_msg("------------------------------------------------------------")
  log_msg("[INFO] 构建矩阵: ", label_universe)

  include_batches <- universe_cfg$include_batches %||% unique(long$batch)

  long_u <- long %>%
    filter(batch %in% include_batches)

  if (!is.null(universe_cfg$classes)) {
    long_u <- long_u %>% filter(Class %in% universe_cfg$classes)
  }
  if (!is.null(universe_cfg$exclude_classes) &&
      length(universe_cfg$exclude_classes) > 0) {
    long_u <- long_u %>% filter(!Class %in% universe_cfg$exclude_classes)
  }

  log_msg("[INFO] universe = ", label_universe,
          " ; nrow = ", nrow(long_u),
          " ; n_features(粗略) = ", n_distinct(long_u$lipidName))

  if (nrow(long_u) == 0) {
    warning("Universe ", label_universe, " 为空，跳过。")
    return(invisible(NULL))
  }

  ## 统一 feature_id
  long_u <- long_u %>%
    mutate(
      feature_id = paste0(
        lipidName, " | ", LipidIon, " | ",
        sprintf("%.4f", CalcMz)
      )
    )

  ## 为 feature 注释表构建唯一的 feature_id 层级注释
  feature_ann_raw <- long_u %>%
    distinct(
      feature_id, lipidName, LipidIon, Class, CalcMz
    )

  ## 检查是否存在同一个 feature_id 对应多条注释记录（通常是 CalcMz 精度差异导致）
  dup_ann <- feature_ann_raw %>%
    group_by(feature_id) %>%
    summarise(
      n_ann = n(),
      .groups = "drop"
    ) %>%
    filter(n_ann > 1)

  if (nrow(dup_ann) > 0) {
    log_msg(
      "[WARN] 在 universe ", label_universe,
      " 中发现 ", nrow(dup_ann),
      " 个 feature_id 拥有多条注释记录（可能是 CalcMz 精度差异导致）；将按 first() 合并注释。"
    )
  }

  feature_ann <- feature_ann_raw %>%
    group_by(feature_id) %>%
    summarise(
      lipidName = first(lipidName),
      LipidIon  = first(LipidIon),
      Class     = first(Class),
      CalcMz    = first(CalcMz),
      .groups   = "drop"
    )

  ## 3.1 检测率统计
  ## 全局检测率
  stats_global <- long_u %>%
    group_by(feature_id) %>%
    summarise(
      n_samples_total = n_distinct(global_sample_id),
      n_non_na_global = sum(!is.na(.data[[expr_col]])),
      detect_rate_global = n_non_na_global / n_samples_total,
      .groups = "drop"
    )

  ## 分 group 检测率（只对有 group 标记的样本）
  long_u_with_group <- long_u %>%
    left_join(sample_info_min %>% select(global_sample_id, group),
              by = "global_sample_id")

  stats_group <- long_u_with_group %>%
    filter(!is.na(group)) %>%
    group_by(feature_id, group) %>%
    summarise(
      n_samples_group = n_distinct(global_sample_id),
      n_non_na_group  = sum(!is.na(.data[[expr_col]])),
      detect_rate_group = n_non_na_group / n_samples_group,
      .groups = "drop"
    )

  stats_group_min <- stats_group %>%
    group_by(feature_id) %>%
    summarise(
      min_detect_rate_group = min(detect_rate_group),
      .groups = "drop"
    )

  ## 3.2 方差统计（log2 scale）
  stats_var <- long_u %>%
    group_by(feature_id) %>%
    summarise(
      var_log = var(log2(.data[[expr_col]] + 1), na.rm = TRUE),
      .groups = "drop"
    )

  ## 3.3 合并统计并应用过滤阈值
  thr_global <- filter_cfg$min_detect_rate_global %||% 0
  thr_group  <- filter_cfg$min_detect_rate_per_group %||% 0
  thr_var_pct <- filter_cfg$min_variance_percentile %||% 0

  stats_all <- stats_global %>%
    left_join(stats_group_min, by = "feature_id") %>%
    left_join(stats_var,      by = "feature_id")

  ## 方差百分位数阈值
  if (thr_var_pct > 0) {
    cutoff_var <- stats_all %>%
      filter(!is.na(var_log)) %>%
      pull(var_log) %>%
      quantile(probs = thr_var_pct, na.rm = TRUE)
  } else {
    cutoff_var <- -Inf
  }

  stats_all <- stats_all %>%
    mutate(
      pass_global = detect_rate_global >= thr_global,
      pass_group  = is.na(min_detect_rate_group) |
                    (min_detect_rate_group >= thr_group),
      pass_var    = is.na(var_log) | (var_log >= cutoff_var),
      keep_feature = pass_global & pass_group & pass_var
    )

  n_keep <- sum(stats_all$keep_feature, na.rm = TRUE)
  n_all  <- nrow(stats_all)
  log_msg("[INFO] feature 过滤结果 (", label_universe, "): keep = ",
          n_keep, " / ", n_all)

  ## 3.4 构建矩阵
  keep_ids <- stats_all %>%
    filter(keep_feature) %>%
    pull(feature_id)

  long_keep <- long_u %>%
    filter(feature_id %in% keep_ids)

  matrix_long <- long_keep %>%
    select(feature_id, global_sample_id, value = all_of(expr_col)) %>%
    group_by(feature_id, global_sample_id) %>%
    summarise(
      value = mean(value, na.rm = TRUE),
      .groups = "drop"
    )

  mat_wide <- matrix_long %>%
    tidyr::pivot_wider(
      id_cols = feature_id,
      names_from = global_sample_id,
      values_from = value
    )

  ## log2 transform
  numeric_cols <- setdiff(colnames(mat_wide), "feature_id")
  mat_wide[numeric_cols] <- lapply(mat_wide[numeric_cols],
                                   function(x) log2(x + 1))

  ## 添加注释列
  mat_out <- feature_ann %>%
    right_join(mat_wide, by = "feature_id") %>%
    relocate(feature_id, lipidName, LipidIon, Class, CalcMz)

  ## 去重：如果在 join 之后仍然存在重复的 feature_id（例如上游注释异常导致），进行合并
  if (any(duplicated(mat_out$feature_id))) {
    dup_tbl <- mat_out %>%
      count(feature_id) %>%
      filter(n > 1)

    log_msg(
      "[WARN] 在矩阵 ", label_universe,
      " 中发现 ", nrow(dup_tbl),
      " 个重复 feature_id；将对这些 feature 进行合并（first()）。"
    )

    mat_out <- mat_out %>%
      group_by(feature_id) %>%
      summarise(
        lipidName = first(lipidName),
        LipidIon  = first(LipidIon),
        Class     = first(Class),
        CalcMz    = first(CalcMz),
        across(-c(feature_id, lipidName, LipidIon, Class, CalcMz), first),
        .groups = "drop"
      ) %>%
      relocate(feature_id, lipidName, LipidIon, Class, CalcMz)
  }

  dir.create(dirname(out_path_matrix), recursive = TRUE, showWarnings = FALSE)
  readr::write_tsv(mat_out, out_path_matrix)
  log_msg("[OK] 写出矩阵 (", label_universe, "): ", out_path_matrix)

  ## 顺带输出一个 QA 统计文件
  qa_path <- file.path(
    dirname(out_path_matrix),
    paste0("feature_filter_stats_", label_universe, ".tsv")
  )
  qa_tbl <- stats_all %>%
    left_join(feature_ann, by = "feature_id") %>%
    relocate(feature_id, lipidName, LipidIon, Class, CalcMz)
  readr::write_tsv(qa_tbl, qa_path)
  log_msg("[OK] 写出过滤 QA: ", qa_path)

  invisible(list(
    stats = stats_all,
    matrix = mat_out
  ))
}

## ------------------------------------------------------------
## 4. 构建 allfat & CL 矩阵
## ------------------------------------------------------------
log_msg("[STEP3] 构建 allfat / CL 矩阵 ...")

# allfat
build_matrix_for_universe(
  long = long,
  sample_info_min = sample_info_min,
  expr_col = expr_col,
  universe_cfg = cfg$universe$allfat,
  filter_cfg = cfg$feature_filter$allfat,
  label_universe = "allfat",
  out_path_matrix = paths$matrix_allfat_filtered %||% paths$matrix_allfat
)

# CL
if (!is.null(cfg$universe$CL)) {
  build_matrix_for_universe(
    long = long,
    sample_info_min = sample_info_min,
    expr_col = expr_col,
    universe_cfg = cfg$universe$CL,
    filter_cfg = cfg$feature_filter$CL,
    label_universe = "CL",
    out_path_matrix = paths$matrix_CL_filtered %||% paths$matrix_CL
  )
} else {
  log_msg("[INFO] 配置中未定义 universe$CL，跳过 CL 矩阵。")
}

log_msg("============================================================")
log_msg("[DONE] 06_prepare_lipid_matrices.R 完成。")
log_msg("  下一步：差异分析 / PCA / heatmap 等 downstream。")
log_msg("============================================================")