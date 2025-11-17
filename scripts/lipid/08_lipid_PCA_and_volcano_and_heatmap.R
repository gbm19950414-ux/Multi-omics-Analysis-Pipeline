#!/usr/bin/env Rscript

## ============================================================
## 08_lipid_PCA_and_volcano_and_heatmap.R  （简化稳定版）
##
## 基于：
##   - 06_prepare_lipid_matrices.R 生成的矩阵 (allfat / CL)
##   - 07_lipid_DE_allfat_and_CL.R 生成的 DE 结果
##
## 生成：
##   1) PCA 图（allfat / CL）
##   2) 火山图（allfat / CL）
##   3) 热图（allfat：most-variable；CL：基于 DE）
##
## 统一 NA 处理策略（PCA + Heatmap）：
##   - 仅保留总体缺失比例 < max_missing_prop_overall (默认 0.3) 的 feature
##   - 对剩余 feature，用同一 group(HO/WT) 内中位数插补缺失；
##     若仍有缺失，则用该 feature 全局中位数插补。
##
## 关键简化：
##   - 热图不再使用 annotation_col，不画样本注释条，
##     从根本上绕开 pheatmap 在 annotation 上的各种坑。
##
## 依赖 00_lipid_downstream_config.yaml 中的 paths / de / na / plots 参数。
## ============================================================

suppressPackageStartupMessages({
  library(yaml)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(ggplot2)
})

## pheatmap 用于热图；如未安装，则跳过热图部分
has_pheatmap <- requireNamespace("pheatmap", quietly = TRUE)

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
       "\nUsage: Rscript 08_lipid_PCA_and_volcano_and_heatmap.R path/to/00_lipid_downstream_config.yaml")
}

log_msg("[INFO] 使用配置文件: ", config_path)
cfg   <- yaml::read_yaml(config_path)
paths <- cfg$paths

plots_dir <- paths$plots_dir %||% "results/lipid/plots"
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

## DE 参数（全局共用）
de_fdr_cutoff   <- cfg$de$fdr_cutoff   %||% 0.05
de_logFC_cutoff <- cfg$de$logFC_cutoff %||% 1

## NA / 缺失处理参数（全局共用）
max_missing_prop_overall <- cfg$na$max_missing_prop_overall %||% 0.3

## Heatmap 缺失率阈值（可以比整体更严格）
heatmap_max_missing_prop_allfat <- cfg$plots$heatmap_max_missing_prop_allfat %||% 0.1
heatmap_max_missing_prop_CL     <- cfg$plots$heatmap_max_missing_prop_CL     %||% max_missing_prop_overall

## Heatmap top N
heatmap_top_n_allfat <- cfg$plots$heatmap_top_n_allfat %||% 50
heatmap_top_n_CL     <- cfg$plots$heatmap_top_n_CL     %||% 50

log_msg("[INFO] Heatmap(allfat) 缺失率阈值: ", heatmap_max_missing_prop_allfat)
log_msg("[INFO] Heatmap(CL)     缺失率阈值: ", heatmap_max_missing_prop_CL)
log_msg("[INFO] DE 阈值: FDR < ", de_fdr_cutoff,
        " & |logFC| >= ", de_logFC_cutoff)
log_msg("[INFO] NA/缺失处理: 保留总体缺失 < ", max_missing_prop_overall,
        " 的 feature，并按 group 中位数插补 (不足时用全局中位数)。")

## ------------------------------------------------------------
## 1. 读 sample_info
## ------------------------------------------------------------
sample_info_path <- paths$sample_info_output %||% "results/lipid/tables/sample_info_lipid.tsv"
if (!file.exists(sample_info_path)) {
  stop("sample_info 文件不存在: ", sample_info_path)
}
sample_info <- readr::read_tsv(sample_info_path, show_col_types = FALSE)

if (!all(c("global_sample_id","batch","group") %in% colnames(sample_info))) {
  stop("sample_info 中必须包含 global_sample_id, batch, group 列。")
}

## ------------------------------------------------------------
## 2. 通用 helper：构建插补后的表达矩阵
## ------------------------------------------------------------

# 返回 numeric matrix，行 = feature_id，列 = global_sample_id
# 步骤：
#  1) 可选：按 feature_subset 子集 + 保持顺序
#  2) 计算总体缺失率 (NA/NaN/Inf)，保留 missing_prop <= max_missing_prop 的 feature
#  3) 按 group 中位数插补 NA / NaN / Inf，剩余用全局中位数补
build_imputed_expr <- function(mat,
                               sample_info,
                               universe_label    = "allfat",
                               max_missing_prop  = 0.3,
                               feature_subset    = NULL,
                               log_prefix        = "[NA]") {

  if (!"feature_id" %in% colnames(mat)) {
    stop(log_prefix, " universe = ", universe_label,
         " ; 矩阵中缺少 feature_id 列。")
  }

  ## 可选：只保留给定 feature_subset，并按给定顺序排序
  if (!is.null(feature_subset)) {
    mat <- mat %>%
      filter(feature_id %in% feature_subset) %>%
      mutate(feature_id = factor(feature_id, levels = feature_subset)) %>%
      arrange(feature_id)
    mat$feature_id <- as.character(mat$feature_id)
  }

  ann_cols <- c("feature_id","lipidName","LipidIon","Class","CalcMz")
  sample_cols <- setdiff(colnames(mat), ann_cols)

  if (length(sample_cols) == 0) {
    warning(log_prefix, " universe = ", universe_label,
            " ; 找不到样本列，返回空矩阵。")
    return(matrix(nrow = 0, ncol = 0))
  }

  expr <- as.matrix(mat[, sample_cols, drop = FALSE])
  mode(expr) <- "numeric"  # 确保 numeric
  rownames(expr) <- mat$feature_id

  n_before <- nrow(expr)

  ## 计算总体缺失比例（NA / NaN / Inf 都视为缺失）
  missing_mat  <- !is.finite(expr)
  missing_prop <- rowMeans(missing_mat)
  missing_prop[is.na(missing_prop)] <- 1

  keep_idx <- which(missing_prop <= max_missing_prop)
  n_keep   <- length(keep_idx)

  log_msg(log_prefix, " universe = ", universe_label,
          " ; feature 保留数(缺失比例<=", max_missing_prop, "): ",
          n_keep, " / ", n_before)

  if (n_keep == 0) {
    warning(log_prefix, " universe = ", universe_label,
            " ; 无 feature 满足缺失率阈值，返回空矩阵。")
    return(matrix(nrow = 0, ncol = 0))
  }

  expr_keep <- expr[keep_idx, , drop = FALSE]
  feat_ids_keep    <- rownames(expr_keep)
  sample_cols_keep <- colnames(expr_keep)

  ## 对应的 group 信息（按列顺序对齐）
  ann <- sample_info %>%
    filter(global_sample_id %in% sample_cols_keep) %>%
    distinct(global_sample_id, group) %>%
    arrange(match(global_sample_id, sample_cols_keep))

  group_vec <- ann$group
  uniq_groups <- sort(unique(group_vec[!is.na(group_vec)]))

  ## 按行做插补
  for (i in seq_len(nrow(expr_keep))) {
    x <- expr_keep[i, ]
    if (!any(!is.finite(x))) {
      next
    }

    x_imp <- x

    ## 先 group 中位数插补（只在有非 NA group 的样本上）
    for (g in uniq_groups) {
      idx_g  <- which(group_vec == g)
      vals_g <- x[idx_g]
      if (any(is.finite(vals_g))) {
        med_g <- median(vals_g[is.finite(vals_g)], na.rm = TRUE)
        fill_idx <- idx_g[!is.finite(vals_g)]
        if (length(fill_idx) > 0 && is.finite(med_g)) {
          x_imp[fill_idx] <- med_g
        }
      }
    }

    ## 仍有非 finite，则用全局 median 插补
    if (any(!is.finite(x_imp))) {
      med_all <- median(x_imp[is.finite(x_imp)], na.rm = TRUE)
      if (is.finite(med_all)) {
        x_imp[!is.finite(x_imp)] <- med_all
      }
    }

    expr_keep[i, ] <- x_imp
  }

  return(expr_keep)
}

## ------------------------------------------------------------
## 3. 通用函数：PCA / 火山图 / 热图
## ------------------------------------------------------------

## 3.1 PCA 函数
plot_pca_for_universe <- function(matrix_path,
                                  sample_info,
                                  universe_label = "allfat",
                                  out_png,
                                  max_missing_prop = 0.3) {

  log_msg("[PCA] universe = ", universe_label,
          " ; 使用矩阵: ", matrix_path)

  if (!file.exists(matrix_path)) {
    warning("[PCA] 矩阵文件不存在，跳过: ", matrix_path)
    return(invisible(NULL))
  }

  mat <- readr::read_tsv(matrix_path, show_col_types = FALSE)

  expr_imputed <- build_imputed_expr(
    mat               = mat,
    sample_info       = sample_info,
    universe_label    = universe_label,
    max_missing_prop  = max_missing_prop,
    feature_subset    = NULL,
    log_prefix        = "[PCA]"
  )

  if (nrow(expr_imputed) < 2) {
    warning("[PCA] 可用于 PCA 的 feature 数 < 2，跳过: ", universe_label)
    return(invisible(NULL))
  }

  expr_t <- t(expr_imputed)

  pca_res <- prcomp(expr_t, center = TRUE, scale. = TRUE)

  pca_df <- as.data.frame(pca_res$x)
  pca_df$global_sample_id <- rownames(pca_res$x)

  pca_df <- pca_df %>%
    left_join(sample_info, by = "global_sample_id")

  var_explained <- (pca_res$sdev^2) / sum(pca_res$sdev^2)
  pc1_var <- sprintf("%.1f%%", var_explained[1] * 100)
  pc2_var <- sprintf("%.1f%%", var_explained[2] * 100)

  p <- ggplot(pca_df, aes(x = PC1, y = PC2,
                          color = group, shape = batch)) +
    geom_point(size = 3, alpha = 0.9) +
    theme_bw() +
    labs(
      title = paste0("PCA (", universe_label, ")"),
      subtitle = "log2 表达矩阵（缺失率过滤 + group 中位数插补）",
      x = paste0("PC1 (", pc1_var, ")"),
      y = paste0("PC2 (", pc2_var, ")"),
      color = "Group",
      shape = "Batch"
    ) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.title = element_text(face = "bold")
    )

  dir.create(dirname(out_png), recursive = TRUE, showWarnings = FALSE)
  ggsave(out_png, p, width = 6, height = 5, dpi = 300)
  log_msg("[PCA] 写出 PCA 图: ", out_png)

  invisible(p)
}

## 3.2 火山图（不需要矩阵，不涉及 NA 插补）
plot_volcano_for_universe <- function(de_path,
                                      universe_label = "allfat",
                                      out_png,
                                      fdr_cutoff = 0.05,
                                      logFC_cutoff = 1) {

  log_msg("[Volcano] universe = ", universe_label,
          " ; 使用 DE 文件: ", de_path)

  if (!file.exists(de_path)) {
    warning("[Volcano] DE 文件不存在，跳过: ", de_path)
    return(invisible(NULL))
  }

  de <- readr::read_tsv(de_path, show_col_types = FALSE)

  required_cols <- c("feature_id","logFC","adj.P.Val")
  if (!all(required_cols %in% colnames(de))) {
    warning("[Volcano] DE 文件缺少必要列 (feature_id, logFC, adj.P.Val)，跳过: ", de_path)
    return(invisible(NULL))
  }

  de <- de %>%
    mutate(
      neg_log10_FDR = -log10(adj.P.Val),
      sig = case_when(
        adj.P.Val < fdr_cutoff & abs(logFC) >= logFC_cutoff ~ "sig",
        TRUE ~ "nonsig"
      ),
      direction = case_when(
        adj.P.Val < fdr_cutoff & logFC >  logFC_cutoff ~ "up",
        adj.P.Val < fdr_cutoff & logFC < -logFC_cutoff ~ "down",
        TRUE ~ "ns"
      )
    )

  n_sig <- sum(de$sig == "sig")
  log_msg("[Volcano] universe = ", universe_label,
          " ; sig feature 数(FDR<", fdr_cutoff,
          " & |logFC|>=", logFC_cutoff, "): ", n_sig)

  p <- ggplot(de, aes(x = logFC, y = neg_log10_FDR)) +
    geom_point(aes(color = direction), alpha = 0.7, size = 1.5) +
    scale_color_manual(
      values = c(
        "up"   = "red",
        "down" = "blue",
        "ns"   = "grey70"
      )
    ) +
    theme_bw() +
    geom_vline(xintercept = c(-logFC_cutoff, logFC_cutoff),
               linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(fdr_cutoff),
               linetype = "dashed", color = "black") +
    labs(
      title = paste0("Volcano plot (", universe_label, ")"),
      subtitle = paste0("FDR <", fdr_cutoff,
                        " & |logFC| >=", logFC_cutoff),
      x = "log2( HO / WT )",
      y = "-log10(FDR)",
      color = "Direction"
    ) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.title = element_text(face = "bold")
    )

  dir.create(dirname(out_png), recursive = TRUE, showWarnings = FALSE)
  ggsave(out_png, p, width = 6, height = 5, dpi = 300)
  log_msg("[Volcano] 写出火山图: ", out_png)

  invisible(p)
}

## 3.3 热图（allfat: most-variable；CL: DE-based）
plot_heatmap_for_universe <- function(matrix_path,
                                      de_path       = NULL,
                                      sample_info,
                                      universe_label = "allfat",
                                      out_png,
                                      top_n          = 50,
                                      fdr_cutoff     = 0.05,
                                      logFC_cutoff   = 1,
                                      max_missing_prop = 0.3,
                                      use_de         = TRUE) {

  if (!has_pheatmap) {
    warning("[Heatmap] 未安装 pheatmap 包，跳过热图绘制。")
    return(invisible(NULL))
  }

  if (!file.exists(matrix_path)) {
    warning("[Heatmap] 矩阵文件不存在，跳过: ", matrix_path)
    return(invisible(NULL))
  }

  if (use_de && !file.exists(de_path)) {
    warning("[Heatmap] DE 文件不存在，跳过: ", de_path)
    return(invisible(NULL))
  }

  mode_desc <- if (use_de) "DE-based" else "most-variable"
  log_msg("[Heatmap] universe = ", universe_label,
          " ; 模式 = ", mode_desc,
          " ; 矩阵: ", matrix_path,
          if (use_de) paste0(" ; DE: ", de_path) else "")

  mat <- readr::read_tsv(matrix_path, show_col_types = FALSE)

  if (!"feature_id" %in% colnames(mat)) {
    warning("[Heatmap] 矩阵中缺少 feature_id，跳过: ", universe_label)
    return(invisible(NULL))
  }

  feature_subset <- NULL

  if (use_de) {
    de <- readr::read_tsv(de_path, show_col_types = FALSE)
    if (!all(c("feature_id","logFC","adj.P.Val") %in% colnames(de))) {
      warning("[Heatmap] DE 文件中缺少 feature_id/logFC/adj.P.Val，跳过: ", universe_label)
      return(invisible(NULL))
    }

    de <- de %>%
      mutate(
        is_sig = adj.P.Val < fdr_cutoff & abs(logFC) >= logFC_cutoff
      ) %>%
      arrange(desc(is_sig), adj.P.Val)

    feature_subset <- de %>%
      slice_head(n = top_n) %>%
      pull(feature_id)

    log_msg("[Heatmap] universe = ", universe_label,
            " ; DE-based 候选 feature 数(top_n): ", length(feature_subset))
  }

  ## 构建插补后的表达矩阵：
  ##   - allfat: use_de = FALSE 时，先对全矩阵做缺失率过滤，然后选 topN variance
  ##   - CL: use_de = TRUE 时，先用 feature_subset 限定，再做缺失率过滤
  expr_imputed <- NULL

  if (use_de) {
    expr_imputed <- build_imputed_expr(
      mat               = mat,
      sample_info       = sample_info,
      universe_label    = universe_label,
      max_missing_prop  = max_missing_prop,
      feature_subset    = feature_subset,
      log_prefix        = "[Heatmap]"
    )

    if (nrow(expr_imputed) == 0) {
      warning("[Heatmap] universe = ", universe_label,
              " ; topN 中无 feature 满足缺失率阈值，跳过热图。")
      return(invisible(NULL))
    }

  } else {
    ## most-variable 模式
    expr_all <- build_imputed_expr(
      mat               = mat,
      sample_info       = sample_info,
      universe_label    = universe_label,
      max_missing_prop  = max_missing_prop,
      feature_subset    = NULL,
      log_prefix        = "[Heatmap]"
    )

    if (nrow(expr_all) == 0) {
      warning("[Heatmap] universe = ", universe_label,
              " ; 无 feature 满足缺失率阈值，跳过热图。")
      return(invisible(NULL))
    }

    vars <- apply(expr_all, 1, stats::var, na.rm = TRUE)
    ord  <- order(vars, decreasing = TRUE)
    n_sel <- min(top_n, length(ord))
    sel_idx <- ord[seq_len(n_sel)]
    expr_imputed <- expr_all[sel_idx, , drop = FALSE]

    log_msg("[Heatmap] universe = ", universe_label,
            " ; most-variable feature 行数: ", nrow(expr_imputed))
  }

  ## 强制 numeric matrix
  expr <- suppressWarnings(matrix(
    as.numeric(expr_imputed),
    nrow = nrow(expr_imputed),
    ncol = ncol(expr_imputed),
    dimnames = dimnames(expr_imputed)
  ))

  ## 显式行标准化（代替 pheatmap 的 scale="row"）
  row_means <- rowMeans(expr, na.rm = TRUE)
  row_sds   <- apply(expr, 1, stats::sd, na.rm = TRUE)
  expr_scaled <- expr
  for (i in seq_len(nrow(expr))) {
    if (is.finite(row_sds[i]) && row_sds[i] > 0) {
      expr_scaled[i, ] <- (expr[i, ] - row_means[i]) / row_sds[i]
    } else {
      expr_scaled[i, ] <- 0
    }
  }

  dir.create(dirname(out_png), recursive = TRUE, showWarnings = FALSE)

  ## 关键简化：不传 annotation_col，完全绕过注释颜色问题
  pheatmap::pheatmap(
    expr_scaled,
    scale = "none",
    show_rownames = FALSE,
    show_colnames = FALSE,
    clustering_method = "complete",
    main = paste0("Heatmap (", universe_label,
                  " ; ", mode_desc,
                  " ; ", nrow(expr_scaled), " features)"),
    filename = out_png,
    width = 8,
    height = 7
  )

  log_msg("[Heatmap] 写出热图: ", out_png)

  invisible(NULL)
}

## ------------------------------------------------------------
## 4. 针对 allfat / CL 分别跑 PCA / Volcano / Heatmap
## ------------------------------------------------------------

log_msg("============================================================")
log_msg("[INFO] 08_lipid_PCA_and_volcano_and_heatmap.R 开始执行")
log_msg("============================================================")

## allfat
matrix_allfat <- paths$matrix_allfat_filtered %||%
  paths$matrix_allfat %||%
  "results/lipid/tables/lipid_matrix_allfat_log2_filtered.tsv"

de_allfat <- paths$de_allfat %||%
  "results/lipid/de/de_allfat_HO_vs_WT.tsv"

pca_allfat_png      <- file.path(plots_dir, "pca_allfat_PC1_PC2.png")
volcano_allfat_png  <- file.path(plots_dir, "volcano_allfat_HO_vs_WT.png")
heatmap_allfat_png  <- file.path(plots_dir, "heatmap_allfat_topN_most_variable.png")

log_msg("------------------------------------------------------------")
log_msg("[INFO] 处理 universe: allfat")

plot_pca_for_universe(
  matrix_path        = matrix_allfat,
  sample_info        = sample_info,
  universe_label     = "allfat",
  out_png            = pca_allfat_png,
  max_missing_prop   = max_missing_prop_overall
)

plot_volcano_for_universe(
  de_path        = de_allfat,
  universe_label = "allfat",
  out_png        = volcano_allfat_png,
  fdr_cutoff     = de_fdr_cutoff,
  logFC_cutoff   = de_logFC_cutoff
)

plot_heatmap_for_universe(
  matrix_path        = matrix_allfat,
  de_path            = de_allfat,      # 这里虽然传了，但 use_de = FALSE，不会用
  sample_info        = sample_info,
  universe_label     = "allfat",
  out_png            = heatmap_allfat_png,
  top_n              = heatmap_top_n_allfat,
  fdr_cutoff         = de_fdr_cutoff,
  logFC_cutoff       = de_logFC_cutoff,
  max_missing_prop   = heatmap_max_missing_prop_allfat,
  use_de             = FALSE          # allfat: most-variable 模式
)

## CL
matrix_CL <- paths$matrix_CL_filtered %||%
  paths$matrix_CL %||%
  "results/lipid/tables/lipid_matrix_CL_log2_filtered.tsv"

de_CL <- paths$de_CL %||%
  "results/lipid/de/de_CL_HO_vs_WT.tsv"

pca_CL_png      <- file.path(plots_dir, "pca_CL_PC1_PC2.png")
volcano_CL_png  <- file.path(plots_dir, "volcano_CL_HO_vs_WT.png")
heatmap_CL_png  <- file.path(plots_dir, "heatmap_CL_topN_DE_based.png")

log_msg("------------------------------------------------------------")
log_msg("[INFO] 处理 universe: CL")

plot_pca_for_universe(
  matrix_path        = matrix_CL,
  sample_info        = sample_info,
  universe_label     = "CL",
  out_png            = pca_CL_png,
  max_missing_prop   = max_missing_prop_overall
)

plot_volcano_for_universe(
  de_path        = de_CL,
  universe_label = "CL",
  out_png        = volcano_CL_png,
  fdr_cutoff     = de_fdr_cutoff,
  logFC_cutoff   = de_logFC_cutoff
)

plot_heatmap_for_universe(
  matrix_path        = matrix_CL,
  de_path            = de_CL,
  sample_info        = sample_info,
  universe_label     = "CL",
  out_png            = heatmap_CL_png,
  top_n              = heatmap_top_n_CL,
  fdr_cutoff         = de_fdr_cutoff,
  logFC_cutoff       = de_logFC_cutoff,
  max_missing_prop   = heatmap_max_missing_prop_CL,
  use_de             = TRUE          # CL: DE-based 模式
)

log_msg("============================================================")
log_msg("[DONE] 08_lipid_PCA_and_volcano_and_heatmap.R 完成。")
log_msg("  请查看: ", plots_dir,
        " 中的 pca_* / volcano_* / heatmap_* 图。")
log_msg("============================================================")