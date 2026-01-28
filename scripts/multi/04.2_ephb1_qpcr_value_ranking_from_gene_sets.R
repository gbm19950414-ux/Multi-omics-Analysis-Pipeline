#!/usr/bin/env Rscript

## ============================================================
## 04.2_ephb1_qpcr_value_ranking_from_gene_sets.R
##
## 目的：
##   基于 per-gene per-batch log2FC（KO vs WT）表与下游通路 gene set YAML，
##   对每个 gene set 内的基因进行“qPCR 验证价值”排序，并输出：
##     - 每个 gene set 一个 *_qpcr_value_ranking.tsv
##     - Tier 标注（Tier1/Tier2/Tier3）
##     - 推荐 top N（Tier1，默认 N=10）
##
## 设计逻辑（可审稿的表述）：
##   对每个 gene set 内每个基因，在统一方向（gene_effect_sign）后，
##   以“效应强度 × 跨批一致性 × 稳定性（低异质性） × 统计证据”综合评分。
##
## 输入（默认，可通过 config 覆盖）：
##   1) per-gene per-batch LFC 表（宽表）：
##        results/rna/deseq2/interaction/interaction_summary_with_perBatchLFC.tsv
##      要求至少包含列：
##        GeneID, IC_stat, IC_pvalue, IC_padj, LFC_batch1, LFC_batch2, ...
##   2) gene set YAML：
##        scripts/multi/ephb1_downstream_signaling_sets.yaml
##      需要同目录下存在 registry：
##        ephb1_downstream_signaling_sets_gene_registry.tsv
##
## 输出：
##   outdir 下：
##     - <gene_set_name>_qpcr_value_ranking.tsv （每个通路一个）
##     - qpcr_value_ranking_topN_summary.tsv      （汇总，可选但默认生成）
##
## 运行：
##   Rscript 04.2_ephb1_qpcr_value_ranking_from_gene_sets.R [config.yaml]
##
## ============================================================

suppressPackageStartupMessages({
  library(yaml)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(tools)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

sanitize_filename <- function(x) {
  x %>%
    as.character() %>%
    stringr::str_replace_all("[\\s/\\\\:+()\"'\t]+", "_") %>%
    stringr::str_replace_all("_+", "_") %>%
    stringr::str_replace_all("^_+|_+$", "")
}

safe_p_evidence <- function(p) {
  # return -log10(p) with guards; NA -> 0
  p <- suppressWarnings(as.numeric(p))
  ifelse(is.na(p) | p <= 0, 0, -log10(p))
}

sign_consistency <- function(v) {
  # v: numeric vector of signed LFC across batches
  v <- v[is.finite(v) & !is.na(v)]
  if (length(v) == 0) return(NA_real_)
  m <- mean(v)
  if (!is.finite(m) || m == 0) return(0)
  s <- sign(m)
  mean(sign(v) == s)
}

## --------- 0. 解析命令行参数 & 配置 -------------------------

args <- commandArgs(trailingOnly = TRUE)

config_path <- if (length(args) >= 1) {
  args[1]
} else {
  "scripts/multiomics/04.2_qpcr_value_ranking_config.yaml"
}

if (file.exists(config_path)) {
  cat("[INFO] 使用配置文件: ", config_path, "\n", sep = "")
  cfg <- yaml::read_yaml(config_path)
} else {
  cat("[WARN] 未找到配置文件: ", config_path,
      "，将使用脚本内置默认路径。\n", sep = "")
  cfg <- list()
}

deseq2_outdir_default     <- "results/rna/deseq2/interaction"
lfc_table_default         <- "interaction_summary_with_perBatchLFC.tsv"
downstream_yaml_default   <- "scripts/multi/ephb1_downstream_signaling_sets.yaml"
outdir_default            <- "results/multiomics/tables/qpcr_value_ranking"

top_n_default             <- 10
tier2_n_default           <- 20
min_batches_present_default <- 2
min_effect_abs_default    <- 0.20
min_sign_consistency_default <- 0.67

deseq2_outdir   <- cfg$deseq2_outdir            %||% deseq2_outdir_default
lfc_table       <- cfg$lfc_table                %||% lfc_table_default
downstream_yaml <- cfg$downstream_gene_set_yaml %||% downstream_yaml_default
outdir          <- cfg$outdir                   %||% outdir_default

top_n           <- as.integer(cfg$top_n %||% top_n_default)
tier2_n         <- as.integer(cfg$tier2_n %||% tier2_n_default)
min_batches_present <- as.integer(cfg$min_batches_present %||% min_batches_present_default)
min_effect_abs  <- as.numeric(cfg$min_effect_abs %||% min_effect_abs_default)
min_sign_consistency <- as.numeric(cfg$min_sign_consistency %||% min_sign_consistency_default)

lfc_path <- file.path(deseq2_outdir, lfc_table)

cat("============================================================\n")
cat("[INFO] 04.2 qPCR 验证价值排序（gene set 内基因）\n")
cat("  LFC table path : ", lfc_path, "\n", sep = "")
cat("  gene set YAML  : ", downstream_yaml, "\n", sep = "")
cat("  outdir         : ", outdir, "\n", sep = "")
cat("  top_n / tier2_n: ", top_n, " / ", tier2_n, "\n", sep = "")
cat("  min_batches    : ", min_batches_present, "\n", sep = "")
cat("  min_abs_effect : ", min_effect_abs, "\n", sep = "")
cat("  min_consistency: ", min_sign_consistency, "\n", sep = "")
cat("============================================================\n\n")

## --------- 1. 读取 per-gene per-batch LFC -------------------

if (!file.exists(lfc_path)) {
  stop("[ERROR] 找不到 LFC 表: ", lfc_path)
}

cat("[STEP1] 读取 LFC 宽表 ...\n")
lfc_wide <- readr::read_tsv(lfc_path, show_col_types = FALSE)

if (!"GeneID" %in% colnames(lfc_wide)) {
  stop("[ERROR] LFC 表中缺少 GeneID 列。")
}

lfc_cols <- grep("^LFC_", colnames(lfc_wide), value = TRUE)
if (length(lfc_cols) == 0) {
  stop("[ERROR] 未找到任何 LFC_* 列（如 LFC_batch1）。")
}

cat("  [INFO] 检测到 LFC 列: ", paste(lfc_cols, collapse = ", "), "\n", sep = "")

need_cols <- c("IC_stat", "IC_pvalue", "IC_padj")
missing_need <- setdiff(need_cols, colnames(lfc_wide))
if (length(missing_need) > 0) {
  cat("  [WARN] LFC 表缺少列: ", paste(missing_need, collapse = ", "),
      "；将按存在列计算（缺失列置 NA）。\n", sep = "")
  for (cc in missing_need) lfc_wide[[cc]] <- NA
}

## --------- 2. 读取 gene set YAML 并映射到 GeneID -----------

if (!file.exists(downstream_yaml)) {
  stop("[ERROR] 找不到 gene set YAML: ", downstream_yaml)
}

cat("\n[STEP2] 读取 gene set YAML: ", downstream_yaml, " ...\n", sep = "")
pw_cfg <- yaml::read_yaml(downstream_yaml)
if (is.null(pw_cfg$pathways)) {
  stop("[ERROR] YAML 中未找到顶层字段 pathways。")
}

pathway_gene_sym <- purrr::imap_dfr(
  pw_cfg$pathways,
  function(pw, pw_name) {
    genes <- pw$genes
    if (is.null(genes)) {
      return(tibble::tibble(
        pathway = character(0), axis = character(0),
        gene_symbol = character(0), gene_effect_sign = numeric(0), weight = numeric(0)
      ))
    }
    purrr::map_dfr(genes, function(g) {
      tibble::tibble(
        pathway          = pw_name,
        axis             = pw$axis %||% NA_character_,
        gene_symbol      = g$name,
        gene_effect_sign = as.numeric(g$gene_effect_sign %||% 1),
        weight           = as.numeric(g$weight %||% 1)
      )
    })
  }
)

yaml_dir  <- dirname(downstream_yaml)
yaml_base <- tools::file_path_sans_ext(basename(downstream_yaml))
registry_path <- file.path(yaml_dir, paste0(yaml_base, "_gene_registry.tsv"))

if (!file.exists(registry_path)) {
  stop("[ERROR] 找不到基因注册表: ", registry_path,
       "\n请先运行 09a.1_update_gene_symbol_to_ensembl.R ", downstream_yaml, " 生成。")
}

registry <- readr::read_tsv(registry_path, show_col_types = FALSE)
if (!all(c("gene_symbol", "ensembl_gene_id") %in% colnames(registry))) {
  stop("[ERROR] registry 必须包含 gene_symbol 与 ensembl_gene_id 列。")
}

pathway_gene_tbl <- pathway_gene_sym %>%
  dplyr::left_join(registry %>% dplyr::select(gene_symbol, ensembl_gene_id), by = "gene_symbol") %>%
  dplyr::rename(GeneID = ensembl_gene_id) %>%
  dplyr::mutate(
    gene_effect_sign = dplyr::if_else(is.na(gene_effect_sign), 1, gene_effect_sign),
    weight           = dplyr::if_else(is.na(weight), 1, weight)
  )

cat("  [INFO] 通路数 = ", length(unique(pathway_gene_tbl$pathway)),
    " ; gene_symbol 数 = ", length(unique(pathway_gene_tbl$gene_symbol)),
    " ; 映射到 GeneID 数 = ", length(unique(pathway_gene_tbl$GeneID)), "\n", sep = "")

if (any(is.na(pathway_gene_tbl$GeneID))) {
  n_na <- sum(is.na(pathway_gene_tbl$GeneID))
  cat("  [WARN] 有 ", n_na, " 个 gene_symbol 未能映射到 GeneID（registry 中缺失）。\n", sep = "")
}

## --------- 3. 对每个通路：计算每基因的综合评分 ----------------

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

lfc_sub <- lfc_wide %>%
  dplyr::select(GeneID, dplyr::all_of(need_cols), dplyr::all_of(lfc_cols))

rank_one_pathway <- function(pw_name) {
  pg <- pathway_gene_tbl %>% dplyr::filter(pathway == pw_name)
  if (nrow(pg) == 0) return(NULL)

  # Map to LFC table (keep genes even if missing some batches; filter later)
  tbl <- pg %>%
    dplyr::left_join(lfc_sub, by = "GeneID")

  # signed LFC across batches
  for (cc in lfc_cols) {
    tbl[[paste0("signed_", cc)]] <- suppressWarnings(as.numeric(tbl[[cc]])) * tbl$gene_effect_sign
  }

  signed_cols <- paste0("signed_", lfc_cols)

  # per-gene metrics
  tbl2 <- tbl %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      n_batches = sum(is.finite(c_across(all_of(signed_cols))) & !is.na(c_across(all_of(signed_cols)))),
      mean_LFC_signed = ifelse(n_batches > 0, mean(c_across(all_of(signed_cols)), na.rm = TRUE), NA_real_),
      mean_abs_LFC = ifelse(n_batches > 0, mean(abs(c_across(all_of(signed_cols))), na.rm = TRUE), NA_real_),
      sd_LFC_signed = ifelse(n_batches > 1, sd(c_across(all_of(signed_cols)), na.rm = TRUE), NA_real_),
      se_LFC_signed = ifelse(n_batches > 1 & is.finite(sd_LFC_signed) & sd_LFC_signed > 0, sd_LFC_signed / sqrt(n_batches), NA_real_),
      z_meta = ifelse(is.finite(se_LFC_signed) & se_LFC_signed > 0, mean_LFC_signed / se_LFC_signed, NA_real_),
      sign_consistency = sign_consistency(as.numeric(c_across(all_of(signed_cols)))),
      min_abs_LFC = ifelse(n_batches > 0, min(abs(c_across(all_of(signed_cols))), na.rm = TRUE), NA_real_),
      p_use = dplyr::coalesce(as.numeric(IC_padj), as.numeric(IC_pvalue)),
      p_evidence = safe_p_evidence(p_use)
    ) %>%
    dplyr::ungroup()

  # composite score (transparent components)
  tbl3 <- tbl2 %>%
    dplyr::mutate(
      # stability penalty: smaller is better
      stability = 1 / (1 + dplyr::coalesce(sd_LFC_signed, 0)),
      # evidence factor: mild boost; 0 -> 1, 10 -> 2, 20 -> 3
      evidence_factor = 1 + (p_evidence / 10),
      score = dplyr::coalesce(mean_abs_LFC, 0) * dplyr::coalesce(sign_consistency, 0) * stability * evidence_factor
    )

  # ranking within pathway
  tbl4 <- tbl3 %>%
    dplyr::mutate(
      eligible = (n_batches >= min_batches_present) & is.finite(score),
      hard_fail_reason = dplyr::case_when(
        n_batches < min_batches_present ~ paste0("n_batches<", min_batches_present),
        TRUE ~ ""
      )
    ) %>%
    dplyr::arrange(dplyr::desc(eligible), dplyr::desc(score), dplyr::desc(mean_abs_LFC)) %>%
    dplyr::mutate(
      rank_in_set = dplyr::if_else(eligible, dplyr::row_number(), NA_integer_)
    )

  # Tiering rules (objective + minimal biological guardrails)
  tbl5 <- tbl4 %>%
    dplyr::mutate(
      pass_basic = eligible & (mean_abs_LFC >= min_effect_abs) & (sign_consistency >= min_sign_consistency),
      Tier = dplyr::case_when(
        pass_basic & !is.na(rank_in_set) & rank_in_set <= top_n ~ "Tier1",
        eligible & !is.na(rank_in_set) & rank_in_set <= tier2_n ~ "Tier2",
        eligible ~ "Tier3",
        TRUE ~ "NotEligible"
      ),
      recommended_topN = (Tier == "Tier1")
    )

  # output columns: keep original LFC columns + signed versions + metrics
  out <- tbl5 %>%
    dplyr::select(
      pathway, axis, gene_symbol, GeneID, gene_effect_sign, weight,
      dplyr::all_of(need_cols), dplyr::all_of(lfc_cols),
      dplyr::all_of(signed_cols),
      n_batches, mean_LFC_signed, mean_abs_LFC, sd_LFC_signed, se_LFC_signed, z_meta,
      sign_consistency, min_abs_LFC, p_use, p_evidence,
      stability, evidence_factor, score,
      eligible, hard_fail_reason,
      rank_in_set, Tier, recommended_topN
    )

  out
}

pathways <- sort(unique(pathway_gene_tbl$pathway))

cat("\n[STEP3] 逐通路计算 qPCR 价值排序 ...\n")
summary_rows <- list()

for (pw in pathways) {
  cat("  [", pw, "] ... ", sep = "")
  out_tbl <- rank_one_pathway(pw)
  if (is.null(out_tbl) || nrow(out_tbl) == 0) {
    cat("SKIP（无基因）\n")
    next
  }

  fname <- paste0(sanitize_filename(pw), "_qpcr_value_ranking.tsv")
  out_path <- file.path(outdir, fname)
  readr::write_tsv(out_tbl, out_path)

  # summary: top genes (Tier1 else top eligible)
  top_tbl <- out_tbl %>%
    dplyr::filter(eligible) %>%
    dplyr::arrange(dplyr::desc(score)) %>%
    dplyr::slice_head(n = top_n)

  top_genes <- paste(top_tbl$gene_symbol, collapse = ",")
  n_tier1 <- sum(out_tbl$Tier == "Tier1", na.rm = TRUE)
  n_elig  <- sum(out_tbl$eligible, na.rm = TRUE)

  summary_rows[[pw]] <- tibble::tibble(
    pathway = pw,
    axis = dplyr::first(out_tbl$axis),
    n_total = nrow(out_tbl),
    n_eligible = n_elig,
    n_tier1 = n_tier1,
    topN_gene_symbols = top_genes,
    ranking_file = out_path
  )

  cat("OK（Tier1=", n_tier1, ", eligible=", n_elig, "） -> ", fname, "\n", sep = "")
}

summary_tbl <- dplyr::bind_rows(summary_rows)
summary_out <- file.path(outdir, "qpcr_value_ranking_topN_summary.tsv")
readr::write_tsv(summary_tbl, summary_out)

cat("\n[OK] 写出汇总：", summary_out, "\n", sep = "")
cat("[DONE] 输出目录：", outdir, "\n", sep = "")
