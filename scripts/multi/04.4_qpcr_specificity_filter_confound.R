#!/usr/bin/env Rscript

## ============================================================
## 04.4_qpcr_specificity_filter_confound.R
##
## 目的：
##   对 results/multiomics/tables/qpcr_value_ranking/ 下所有
##   *_qpcr_value_ranking.tsv 进行 IFN/UPR confound 特异性筛选：
##     - 计算每个 batch 的 Confound_IFN_I_response 与 Confound_UPR_ER_stress score
##     - 对每个候选基因：计算其 signed LFC pattern 与 confound score pattern 的相关性
##     - 输出：*_with_confound.tsv, *_filtered.tsv, 以及总汇总表
##
## 默认输入：
##   ranking_dir: results/multiomics/tables/qpcr_value_ranking
##   LFC table  : results/rna/deseq2/interaction/interaction_summary_with_perBatchLFC.tsv
##   gene set Y : scripts/multi/geneset_for_qPCR.yaml   (可改成你的实际路径)
##   registry   : <yaml_base>_gene_registry.tsv (与 04.2 同逻辑)
##
## 运行：
##   Rscript 04.4_qpcr_specificity_filter_confound.R
##   Rscript 04.4_qpcr_specificity_filter_confound.R <ranking_dir> <geneset_yaml> <lfc_path> [corr_cutoff]
##
## ============================================================

suppressPackageStartupMessages({
  library(yaml)
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
  library(purrr)
  library(tools)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

safe_cor <- function(x, y, method = "pearson") {
  x <- as.numeric(x); y <- as.numeric(y)
  ok <- is.finite(x) & is.finite(y) & !is.na(x) & !is.na(y)
  x <- x[ok]; y <- y[ok]
  if (length(x) < 3) return(NA_real_)  # 3 batches 时才有点意义
  if (sd(x) == 0 || sd(y) == 0) return(NA_real_)
  suppressWarnings(cor(x, y, method = method))
}

sanitize_filename <- function(x) {
  x %>%
    as.character() %>%
    stringr::str_replace_all("[\\s/\\\\:+()\"'\t]+", "_") %>%
    stringr::str_replace_all("_+", "_") %>%
    stringr::str_replace_all("^_+|_+$", "")
}

## ----------------- 0) Args & defaults -----------------------

args <- commandArgs(trailingOnly = TRUE)

ranking_dir_default <- "results/multiomics/tables/qpcr_value_ranking"
geneset_yaml_default <- "scripts/multi/geneset_for_qPCR.yaml"
lfc_path_default <- "results/rna/deseq2/interaction/interaction_summary_with_perBatchLFC.tsv"
corr_cutoff_default <- 0.7

get_arg <- function(i) {
  if (length(args) >= i) {
    v <- args[[i]]
    if (is.null(v) || is.na(v) || v == "") return(NULL)
    return(v)
  }
  NULL
}

ranking_dir <- get_arg(1) %||% ranking_dir_default
geneset_yaml <- get_arg(2) %||% geneset_yaml_default
lfc_path <- get_arg(3) %||% lfc_path_default
corr_cutoff <- as.numeric(get_arg(4) %||% corr_cutoff_default)
if (is.na(corr_cutoff) || !is.finite(corr_cutoff)) corr_cutoff <- corr_cutoff_default

cat("============================================================\n")
cat("[INFO] 04.4 confound specificity filter\n")
cat("  ranking_dir : ", ranking_dir, "\n", sep = "")
cat("  geneset_yaml : ", geneset_yaml, "\n", sep = "")
cat("  lfc_path     : ", lfc_path, "\n", sep = "")
cat("  corr_cutoff  : ", corr_cutoff, "\n", sep = "")
cat("============================================================\n\n")

if (!dir.exists(ranking_dir)) stop("[ERROR] ranking_dir not found: ", ranking_dir)
if (!file.exists(geneset_yaml)) stop("[ERROR] geneset_yaml not found: ", geneset_yaml)
if (!file.exists(lfc_path)) stop("[ERROR] lfc_path not found: ", lfc_path)

## ----------------- 1) Load LFC wide table -------------------

lfc_wide <- readr::read_tsv(lfc_path, show_col_types = FALSE)

if (!"GeneID" %in% colnames(lfc_wide)) stop("[ERROR] LFC table missing GeneID")
lfc_cols <- grep("^LFC_", colnames(lfc_wide), value = TRUE)
if (length(lfc_cols) == 0) stop("[ERROR] No LFC_* columns found (e.g., LFC_batch1)")

cat("[STEP1] Found LFC cols: ", paste(lfc_cols, collapse = ", "), "\n\n", sep = "")

## ----------------- 2) Load gene sets & registry --------------

pw_cfg <- yaml::read_yaml(geneset_yaml)

# 兼容两种结构：
# A) 顶层是 pathways:
# B) 顶层就是各个 geneset key（你现在 geneset_for_qPCR.yaml 可能就是这种）
if (!is.null(pw_cfg$pathways)) {
  pathways_list <- pw_cfg$pathways
} else {
  pathways_list <- pw_cfg
}

# 取 confound gene sets：名字以 Confound_ 开头 或 axis == Confound_control
confound_keys <- names(pathways_list)
confound_keys <- confound_keys[grepl("^Confound_", confound_keys) |
                                 purrr::map_lgl(pathways_list, ~ (.x$axis %||% "") == "Confound_control")]

if (length(confound_keys) == 0) {
  stop("[ERROR] No confound gene sets found in YAML. Expect keys like Confound_IFN_I_response / Confound_UPR_ER_stress")
}

cat("[STEP2] Confound sets: ", paste(confound_keys, collapse = ", "), "\n", sep = "")

yaml_dir  <- dirname(geneset_yaml)
yaml_base <- tools::file_path_sans_ext(basename(geneset_yaml))
registry_path <- file.path(yaml_dir, paste0(yaml_base, "_gene_registry.tsv"))

if (!file.exists(registry_path)) {
  stop("[ERROR] registry not found: ", registry_path,
       "\nPlease generate it (same as 04.2 logic) before running 04.4.")
}

registry <- readr::read_tsv(registry_path, show_col_types = FALSE)
if (!all(c("gene_symbol", "ensembl_gene_id") %in% colnames(registry))) {
  stop("[ERROR] registry must include gene_symbol and ensembl_gene_id")
}

confound_gene_tbl <- purrr::map_dfr(confound_keys, function(k) {
  genes <- pathways_list[[k]]$genes %||% list()
  if (length(genes) == 0) return(tibble())
  purrr::map_dfr(genes, function(g) {
    tibble::tibble(
      confound_set = k,
      gene_symbol = g$name,
      gene_effect_sign = as.numeric(g$gene_effect_sign %||% 1),
      weight = as.numeric(g$weight %||% 1)
    )
  })
}) %>%
  left_join(registry %>% select(gene_symbol, ensembl_gene_id), by = "gene_symbol") %>%
  rename(GeneID = ensembl_gene_id) %>%
  mutate(
    gene_effect_sign = if_else(is.na(gene_effect_sign), 1, gene_effect_sign),
    weight = if_else(is.na(weight) | !is.finite(weight) | weight <= 0, 1, weight)
  )

if (any(is.na(confound_gene_tbl$GeneID))) {
  cat("[WARN] Some confound gene_symbol cannot map to GeneID; they will be ignored.\n")
}

## ----------------- 3) Compute confound scores per batch -------

# Join confound genes to LFC
confound_lfc <- confound_gene_tbl %>%
  filter(!is.na(GeneID)) %>%
  left_join(lfc_wide %>% select(GeneID, all_of(lfc_cols)), by = "GeneID")

# Compute signed LFC per batch for confound sets
for (cc in lfc_cols) {
  confound_lfc[[paste0("signed_", cc)]] <- as.numeric(confound_lfc[[cc]]) * confound_lfc$gene_effect_sign
}

signed_cols <- paste0("signed_", lfc_cols)

# Score per confound set per batch: weighted mean of signed LFC
confound_scores <- confound_lfc %>%
  pivot_longer(cols = all_of(signed_cols), names_to = "signed_batch", values_to = "signed_lfc") %>%
  mutate(batch = str_replace(signed_batch, "^signed_LFC_", "")) %>%
  filter(is.finite(signed_lfc) & !is.na(signed_lfc)) %>%
  group_by(confound_set, batch) %>%
  summarise(
    n_genes = n(),
    score = weighted.mean(signed_lfc, w = weight, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = batch, values_from = score)

cat("[STEP3] Confound scores computed.\n\n")

# Helper to extract score vector in the same order as LFC batches in ranking file
get_confound_vec <- function(confound_set, batch_names) {
  row <- confound_scores %>% filter(confound_set == !!confound_set)
  if (nrow(row) == 0) return(rep(NA_real_, length(batch_names)))
  v <- purrr::map_dbl(batch_names, function(bn) {
    vv <- row[[bn]]
    if (length(vv) == 0) return(NA_real_)
    as.numeric(vv)
  })
  v
}

## ----------------- 4) Scan ranking files & annotate -----------

files <- list.files(ranking_dir, pattern = "_qpcr_value_ranking\\.tsv$", full.names = TRUE)
if (length(files) == 0) stop("[ERROR] No *_qpcr_value_ranking.tsv found in: ", ranking_dir)

cat("[STEP4] Found ranking files: ", length(files), "\n", sep = "")

summary_rows <- list()

for (f in files) {
  cat("  -> ", basename(f), "\n", sep = "")
  df <- readr::read_tsv(f, show_col_types = FALSE)

  # Identify signed LFC columns (preferred), else fallback to LFC columns
  signed_cols_in <- grep("^signed_LFC_", colnames(df), value = TRUE)
  if (length(signed_cols_in) == 0) {
    stop("[ERROR] Ranking file has no signed_LFC_* columns: ", f,
         "\nPlease ensure it is produced by 04.2 (weight-upgraded version is ok).")
  }

  batch_names <- str_replace(signed_cols_in, "^signed_LFC_", "")
  # confound vectors aligned to batch order
  ifn_key <- if ("Confound_IFN_I_response" %in% confound_keys) "Confound_IFN_I_response" else confound_keys[1]
  upr_key <- if ("Confound_UPR_ER_stress" %in% confound_keys) "Confound_UPR_ER_stress" else confound_keys[min(2, length(confound_keys))]

  ifn_vec <- get_confound_vec(ifn_key, batch_names)
  upr_vec <- get_confound_vec(upr_key, batch_names)

  # Compute correlations per gene
  df2 <- df %>%
    rowwise() %>%
    mutate(
      corr_IFN = safe_cor(c_across(all_of(signed_cols_in)), ifn_vec, method = "pearson"),
      corr_UPR = safe_cor(c_across(all_of(signed_cols_in)), upr_vec, method = "pearson"),
      max_abs_confound_corr = max(abs(corr_IFN), abs(corr_UPR), na.rm = TRUE),
      confound_flag = is.finite(max_abs_confound_corr) & (max_abs_confound_corr >= corr_cutoff)
    ) %>%
    ungroup()

  # Recommend filtered topN: Tier1 but not confounded
  df2 <- df2 %>%
    mutate(
      recommended_topN_filtered = (Tier == "Tier1") & (!confound_flag)
    )

  # Write outputs
  base <- str_replace(basename(f), "_qpcr_value_ranking\\.tsv$", "")
  out_with <- file.path(ranking_dir, paste0(base, "_qpcr_value_ranking_with_confound.tsv"))
  out_filt <- file.path(ranking_dir, paste0(base, "_qpcr_value_ranking_filtered.tsv"))

  write_tsv(df2, out_with)

  # filtered table: keep all columns, but sort by original rank / or score_used_for_ranking if exists
  sort_col <- if ("score_used_for_ranking" %in% colnames(df2)) "score_used_for_ranking" else "score"
  df_filt <- df2 %>%
    filter(eligible) %>%
    arrange(desc(.data[[sort_col]])) %>%
    mutate(filtered_rank = row_number())

  write_tsv(df_filt, out_filt)

  # Summaries
  n_total <- nrow(df2)
  n_tier1 <- sum(df2$Tier == "Tier1", na.rm = TRUE)
  n_flag <- sum(df2$confound_flag, na.rm = TRUE)
  n_tier1_kept <- sum(df2$recommended_topN_filtered, na.rm = TRUE)

  top_keep <- df2 %>%
    filter(recommended_topN_filtered) %>%
    arrange(desc(.data[[sort_col]])) %>%
    slice_head(n = 10) %>%
    pull(gene_symbol) %>%
    paste(collapse = ",")

  summary_rows[[f]] <- tibble::tibble(
    file = basename(f),
    out_with_confound = basename(out_with),
    out_filtered = basename(out_filt),
    confound_IFN_key = ifn_key,
    confound_UPR_key = upr_key,
    n_total = n_total,
    n_tier1 = n_tier1,
    n_confound_flag = n_flag,
    n_tier1_filtered_kept = n_tier1_kept,
    top10_kept_gene_symbols = top_keep
  )
}

summary_tbl <- bind_rows(summary_rows)
summary_path <- file.path(ranking_dir, "qpcr_value_ranking_confound_filter_summary.tsv")
write_tsv(summary_tbl, summary_path)

cat("\n[OK] Wrote per-file outputs:\n")
cat("  *_qpcr_value_ranking_with_confound.tsv\n")
cat("  *_qpcr_value_ranking_filtered.tsv\n")
cat("[OK] Summary:\n  ", summary_path, "\n", sep = "")
cat("[DONE]\n")