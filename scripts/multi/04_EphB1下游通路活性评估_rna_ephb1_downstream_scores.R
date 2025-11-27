#!/usr/bin/env Rscript

## ============================================================
## 05_rna_ephb1_downstream_scores.R
##
## 目的：
##   基于基因层面 DESeq2 log2FC（KO vs WT），
##   使用 EphB1 下游信号通路基因集 YAML
##     (ephb1_downstream_signaling_sets.yaml)，
##   计算每条 EphB1 下游信号通路在每个 batch 的 HO vs WT 效应 Z 分数：
##     - mean_LFC（方向已按 gene_effect_sign 校正）
##     - sd_LFC / se_LFC
##     - z_pathway = mean_LFC / se_LFC
##
## 输入（默认，可通过 config 覆盖）：
##
##   1) per-gene per-batch LFC 表：
##        results/rna/deseq2/interaction/interaction_summary_with_perBatchLFC.tsv
##      要求至少包含列：
##        GeneID, 以及若干 "LFC_<batch>" 列（如 LFC_batch1）
##
##   2) EphB1 下游信号通路基因集 YAML：
##        scripts/multi/ephb1_downstream_signaling_sets.yaml
##
##   3) 可选配置文件（命令行第 1 个参数）：
##        scripts/multiomics/00_ephb1_downstream_scores_config.yaml
##
##      支持字段：
##        deseq2_outdir: "results/rna/deseq2/interaction"
##        downstream_gene_set_yaml: "scripts/multi/ephb1_downstream_signaling_sets.yaml"
##        lfc_table: "interaction_summary_with_perBatchLFC.tsv"
##        outdir: "results/multiomics/tables"
##
## 输出：
##
##   - results/multiomics/tables/ephb1_downstream_signaling_scores_per_batch.tsv
##
##     列：
##       pathway, axis, batch,
##       n_genes, mean_LFC, sd_LFC, se_LFC, z_pathway
##
##   z_pathway 的含义：
##     在该 EphB1 下游信号通路内所有基因的 HO vs WT 有效 log2FC
##     （按 gene_effect_sign 校正）的均值，用 pathway 内基因间变异
##     估计不确定性（SE），z_pathway = mean_LFC / se_LFC。
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

## --------- 0. 解析命令行参数 & 配置 -------------------------

args <- commandArgs(trailingOnly = TRUE)

config_path <- if (length(args) >= 1) {
  args[1]
} else {
  "scripts/multiomics/00_ephb1_downstream_scores_config.yaml"
}

if (file.exists(config_path)) {
  cat("[INFO] 使用配置文件: ", config_path, "\n", sep = "")
  cfg <- yaml::read_yaml(config_path)
} else {
  cat("[WARN] 未找到配置文件: ", config_path,
      "，将使用脚本内置默认路径。\n", sep = "")
  cfg <- list()
}

deseq2_outdir_default <- "results/rna/deseq2/interaction"
outdir_default        <- "results/multiomics/tables"
downstream_yaml_default <- "scripts/multi/ephb1_downstream_signaling_sets.yaml"
lfc_table_default     <- "interaction_summary_with_perBatchLFC.tsv"

deseq2_outdir   <- cfg$deseq2_outdir            %||% deseq2_outdir_default
outdir          <- cfg$outdir                   %||% outdir_default
downstream_yaml <- cfg$downstream_gene_set_yaml %||% downstream_yaml_default
lfc_table       <- cfg$lfc_table                %||% lfc_table_default

lfc_path <- file.path(deseq2_outdir, lfc_table)

cat("============================================================\n")
cat("[INFO] 05_rna_ephb1_downstream_scores.R\n")
cat("  LFC table path : ", lfc_path, "\n", sep = "")
cat("  pathway YAML   : ", downstream_yaml, "\n", sep = "")
cat("  outdir         : ", outdir, "\n", sep = "")
cat("============================================================\n\n")

## --------- 1. 读取 per-gene per-batch LFC -------------------

if (!file.exists(lfc_path)) {
  stop("[ERROR] 找不到 LFC 表: ", lfc_path)
}

cat("[STEP1] 读取 per-gene per-batch LFC 表 ...\n")

lfc_wide <- readr::read_tsv(lfc_path, show_col_types = FALSE)

cat("  [lfc_wide] nrow = ", nrow(lfc_wide),
    " ; ncol = ", ncol(lfc_wide), "\n", sep = "")

## 检查 GeneID 列
if (!"GeneID" %in% colnames(lfc_wide)) {
  stop("[ERROR] LFC 表中缺少 GeneID 列，请检查输入文件结构。")
}

## 找出所有 LFC_* 列
lfc_cols <- grep("^LFC_", colnames(lfc_wide), value = TRUE)

if (length(lfc_cols) == 0) {
  stop("[ERROR] 在 LFC 表中未找到任何以 'LFC_' 开头的列，请检查文件。")
}

cat("  [INFO] 检测到的 LFC 列: ",
    paste(lfc_cols, collapse = ", "), "\n", sep = "")

## 将 wide → long：每行 = GeneID × batch
lfc_long <- lfc_wide %>%
  dplyr::select(GeneID, dplyr::all_of(lfc_cols)) %>%
  tidyr::pivot_longer(
    cols      = dplyr::all_of(lfc_cols),
    names_to  = "batch",
    values_to = "log2FC"
  ) %>%
  dplyr::mutate(
    batch = stringr::str_replace(batch, "^LFC_", "")  # "LFC_batch1" -> "batch1"
  )

cat("  [lfc_long] nrow = ", nrow(lfc_long),
    " ; 批次数 = ", length(unique(lfc_long$batch)), "\n", sep = "")

## --------- 2. 读取 EphB1 下游通路基因集 YAML ---------------

if (!file.exists(downstream_yaml)) {
  stop("[ERROR] 找不到 EphB1 下游信号通路基因集 YAML: ", downstream_yaml)
}

cat("\n[STEP2] 读取 EphB1 下游通路 YAML: ", downstream_yaml, " ...\n", sep = "")
pw_cfg <- yaml::read_yaml(downstream_yaml)

if (is.null(pw_cfg$pathways)) {
  stop("[ERROR] 下游通路 YAML 中未找到 'pathways' 顶层字段。")
}

## 先从 YAML 中提取 symbol 层信息
pathway_gene_sym <- purrr::imap_dfr(
  pw_cfg$pathways,
  function(pw, pw_name) {
    genes <- pw$genes
    if (is.null(genes)) {
      return(tibble::tibble(
        pathway          = character(0),
        axis             = character(0),
        gene_symbol      = character(0),
        gene_effect_sign = numeric(0)
      ))
    }
    purrr::map_dfr(genes, function(g) {
      tibble::tibble(
        pathway          = pw_name,
        axis             = pw$axis %||% NA_character_,
        gene_symbol      = g$name,
        gene_effect_sign = as.numeric(g$gene_effect_sign %||% 1)
      )
    })
  }
)

cat("  [pathway_gene_sym] nrow = ", nrow(pathway_gene_sym),
    " ; 通路数 = ", length(unique(pathway_gene_sym$pathway)),
    " ; 涉及 gene_symbol 数 = ", length(unique(pathway_gene_sym$gene_symbol)), "\n", sep = "")

## 使用基因注册表 symbol -> Ensembl GeneID
yaml_dir  <- dirname(downstream_yaml)
yaml_base <- tools::file_path_sans_ext(basename(downstream_yaml))
registry_path <- file.path(yaml_dir, paste0(yaml_base, "_gene_registry.tsv"))

if (!file.exists(registry_path)) {
  stop("[ERROR] 找不到基因注册表文件: ", registry_path,
       "\n请先运行 09a.1_update_gene_symbol_to_ensembl.R ",
       downstream_yaml, " 生成该文件。")
}

registry <- readr::read_tsv(registry_path, show_col_types = FALSE)

if (!all(c("gene_symbol", "ensembl_gene_id") %in% colnames(registry))) {
  stop("[ERROR] 基因注册表中必须包含列: gene_symbol, ensembl_gene_id。实际列: ",
       paste(colnames(registry), collapse = ", "))
}

pathway_gene_tbl <- pathway_gene_sym %>%
  dplyr::left_join(registry %>% dplyr::select(gene_symbol, ensembl_gene_id),
                   by = "gene_symbol") %>%
  dplyr::rename(GeneID = ensembl_gene_id)

cat("  [pathway_gene_tbl] nrow = ", nrow(pathway_gene_tbl),
    " ; 通路数 = ", length(unique(pathway_gene_tbl$pathway)),
    " ; 涉及 GeneID 数 = ", length(unique(pathway_gene_tbl$GeneID)), "\n", sep = "")

## --------- 3. 将 pathway 基因集映射到 LFC -------------------

cat("\n[STEP3] 将 EphB1 下游 pathway 基因映射到 LFC 表 ...\n")

## inner_join：只保留在 LFC 表中有 log2FC 的基因
pw_lfc <- pathway_gene_tbl %>%
  dplyr::inner_join(lfc_long, by = "GeneID") %>%
  dplyr::mutate(
    gene_effect_sign = dplyr::if_else(
      is.na(gene_effect_sign), 1, gene_effect_sign
    ),
    effective_log2FC = log2FC * gene_effect_sign
  )

cat("  [pw_lfc] nrow = ", nrow(pw_lfc),
    " ; pathway×batch 组合数 = ",
    length(unique(paste(pw_lfc$pathway, pw_lfc$batch, sep = "_"))),
    "\n", sep = "")

## --------- 4. 在 pathway × batch 层面汇总 --------------------

cat("\n[STEP4] 计算每个 EphB1 下游 pathway × batch 的 mean_LFC / Z 分数 ...\n")

pathway_scores <- pw_lfc %>%
  dplyr::group_by(pathway, axis, batch) %>%
  dplyr::summarise(
    n_genes  = sum(!is.na(effective_log2FC)),
    mean_LFC = dplyr::if_else(
      n_genes > 0,
      mean(effective_log2FC, na.rm = TRUE),
      NA_real_
    ),
    sd_LFC = dplyr::if_else(
      n_genes > 1,
      sd(effective_log2FC, na.rm = TRUE),
      NA_real_
    ),
    se_LFC = dplyr::if_else(
      n_genes > 1 & !is.na(sd_LFC),
      sd_LFC / sqrt(n_genes),
      NA_real_
    ),
    z_pathway = dplyr::if_else(
      !is.na(se_LFC) & se_LFC > 0,
      mean_LFC / se_LFC,
      NA_real_
    ),
    .groups = "drop"
  )

cat("  [pathway_scores] nrow = ", nrow(pathway_scores),
    " ; pathway 数 = ", length(unique(pathway_scores$pathway)),
    " ; batch 数 = ", length(unique(pathway_scores$batch)), "\n", sep = "")

## --------- 5. 写出结果 --------------------------------------

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

out_path <- file.path(outdir, "ephb1_downstream_signaling_scores_per_batch.tsv")
readr::write_tsv(pathway_scores, out_path)

cat("\n[OK] 写出 EphB1 下游信号通路 Z 分数表: ", out_path, "\n", sep = "")
cat("============================================================\n")
cat("[DONE] 05_rna_ephb1_downstream_scores.R 完成。\n")
cat("  说明：\n")
cat("    - mean_LFC > 0 且 z_pathway > 0：该 EphB1 下游通路在 HO 中整体激活；\n")
cat("    - mean_LFC < 0 且 z_pathway < 0：该通路在 HO 中整体受抑；\n")
cat("    - 绝对值越大，效应越强，可与 upstream 机制轴 Z / Supply / Transport 进行方向匹配。\n")
cat("============================================================\n")