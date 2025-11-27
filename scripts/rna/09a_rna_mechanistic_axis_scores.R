#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(purrr)
  library(yaml)
})

# ============================================================
# 09a_rna_mechanistic_axis_scores.R
#
# Usage（可以传参，也可以不传）:
#   Rscript 09a_rna_mechanistic_axis_scores.R \
#       [deseq2_outdir] \
#       [rna_mechanistic_gene_sets.yaml] \
#       [outdir]
#
# 若不提供参数，则使用默认：
#   deseq2_outdir  = "results/rna/deseq2"
#   gene_set_yaml  = "scripts/rna/rna_mechanistic_gene_sets.yaml"
#   outdir         = "results/rna/mechanistic_axis"
#
# 依赖：
#   - 08_deseq2.R 已运行完毕
#   - deseq2_outdir/interaction/interaction_summary_with_perBatchLFC.tsv 存在
#   - gene_set_yaml: rna_mechanistic_gene_sets.yaml
#   - 若 YAML 中包含 mapping: 字段，则 ref/gene_symbol_to_ensembl.tsv 存在
# ============================================================

## ---------- 解析命令行参数 + 默认路径 ----------
args <- commandArgs(trailingOnly = TRUE)

deseq2_outdir_default <- "results/rna/deseq2"
gene_set_yaml_default <- "scripts/rna/rna_mechanistic_gene_sets.yaml"
outdir_default        <- "results/rna/mechanistic_axis"

deseq2_outdir <- if (length(args) >= 1 && nzchar(args[1])) args[1] else deseq2_outdir_default
gene_set_yaml <- if (length(args) >= 2 && nzchar(args[2])) args[2] else gene_set_yaml_default
outdir        <- if (length(args) >= 3 && nzchar(args[3])) args[3] else outdir_default

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

message("============================================================")
message("[INFO] 09a_rna_mechanistic_axis_scores.R")
message("  deseq2_outdir : ", deseq2_outdir)
message("  gene_set_yaml : ", gene_set_yaml)
message("  outdir        : ", outdir)
message("============================================================")

# ---------- STEP1: 读取 DE 源文件（per-batch LFC） ----------
ic_path <- file.path(deseq2_outdir, "interaction", "interaction_summary_with_perBatchLFC.tsv")
if (!file.exists(ic_path)) {
  stop(
    "[FATAL] interaction_summary_with_perBatchLFC.tsv not found: ", ic_path,
    "\n  请检查: 1) 08_deseq2.R 是否已运行; 2) deseq2_outdir 是否正确。"
  )
}

ic_tbl <- read_tsv(ic_path, show_col_types = FALSE)

if (!"GeneID" %in% colnames(ic_tbl)) {
  stop(
    "[FATAL] interaction_summary_with_perBatchLFC.tsv 缺少 GeneID 列。",
    "\n  当前脚本假定 DE 结果中存在 GeneID 列（通常为 Ensembl ID）。"
  )
}

# LFC_ 开头的列视为各 batch 的 KO vs WT log2FC
lfc_cols <- grep("^LFC_", colnames(ic_tbl), value = TRUE)
if (length(lfc_cols) == 0) {
  stop("[FATAL] 在 ", ic_path, " 中未找到以 'LFC_' 开头的列，请检查 08_deseq2.R 输出。")
}

batch_levels <- sub("^LFC_", "", lfc_cols)

message("[INFO] 读取 per-gene per-batch LFC：")
message("  n_gene = ", nrow(ic_tbl))
message("  LFC columns: ", paste(lfc_cols, collapse = ", "))
message("  batch levels: ", paste(batch_levels, collapse = ", "))

# ---------- STEP2: 读取机制基因集 YAML ----------
if (!file.exists(gene_set_yaml)) {
  stop("[FATAL] gene set YAML not found: ", gene_set_yaml)
}

gs_conf <- yaml::read_yaml(gene_set_yaml)

# 期望结构：axes: { Synthesis: { core_genes: [..], pc_pe_metabolism_genes: [..], ... }, ... }
if (!"axes" %in% names(gs_conf)) {
  stop("[FATAL] YAML 中未找到 'axes' 字段，请检查 rna_mechanistic_gene_sets.yaml。")
}

axes_list <- gs_conf$axes

# helper: 从一个“基因容器”（字符向量 / list / module$genes）里提取 gene_symbol 和 gene_effect_sign
extract_genes_from_container <- function(container) {
  if (is.null(container)) {
    return(tibble(gene_symbol = character(0), gene_effect_sign = numeric(0), weight = numeric(0)))
  }

  # 若是 module 结构，优先使用 $genes 字段
  if (is.list(container) && "genes" %in% names(container)) {
    entries <- container$genes
  } else {
    entries <- container
  }

  # entries 可以是字符向量，或者 list（其中每个元素是字符或带 name/gene_effect_sign 的列表）
  if (is.character(entries)) {
    return(
      tibble(
        gene_symbol      = entries,
        gene_effect_sign = rep(1, length(entries)),
        weight           = rep(1, length(entries))
      )
    )
  }

  if (is.list(entries)) {
    return(
      purrr::map_dfr(entries, function(e) {
        # 单个 gene entry 既可以直接是字符，也可以是 list(name, gene_symbol, gene_effect_sign)
        if (is.character(e)) {
          tibble(
            gene_symbol      = e,
            gene_effect_sign = 1,
            weight           = 1
          )
        } else if (is.list(e)) {
          sym <- e$name
          if (is.null(sym) && !is.null(e$gene_symbol)) {
            sym <- e$gene_symbol
          }
          if (is.null(sym)) {
            stop("[FATAL] gene entry 缺少 name / gene_symbol 字段，请检查 YAML。")
          }
          ges <- e$gene_effect_sign
          if (is.null(ges) || is.na(ges)) {
            ges <- 1
          }
          w <- e$weight
          if (is.null(w) || is.na(w)) {
            w <- 1
          }
          tibble(
            gene_symbol      = as.character(sym),
            gene_effect_sign = as.numeric(ges),
            weight           = as.numeric(w)
          )
        } else {
          stop("[FATAL] 无法解析 gene entry（既不是字符向量也不是带 name/gene_effect_sign 的列表），请检查 YAML。")
        }
      })
    )
  }

  stop("[FATAL] 无法解析 gene container（既不是字符向量也不是列表），请检查 YAML。")
}

# 将每个 axis 下的 gene 列表展平成一张表：axis × gene_symbol，并附带 gene_effect_sign
axis_gene_tbl <- purrr::map2_dfr(
  .x = names(axes_list),
  .y = axes_list,
  .f = function(axis_name, axis_def) {
    # axis_def: 可能包含 description / label，以及：
    #   - 新写法：modules: { module_name: { label / description / genes } }
    #   - 旧写法：直接挂 core_genes / extended_genes 等 gene 容器
    if ("modules" %in% names(axis_def)) {
      # 新写法：axis 下有 modules 层，每个 module 里面再有 genes
      gene_containers <- axis_def$modules
    } else {
      # 兼容旧写法：axis 直接挂 core_genes / extended_genes 等 gene 容器
      meta_fields <- c("description", "label")
      gene_containers <- axis_def[setdiff(names(axis_def), meta_fields)]
    }
    gene_containers <- gene_containers[!vapply(gene_containers, is.null, logical(1))]

    if (length(gene_containers) == 0) {
      warning("[WARN] axis ", axis_name, " 未找到任何基因列表（除 description/label 外为空）。")
      return(tibble(
        axis             = character(0),
        gene_symbol      = character(0),
        gene_effect_sign = numeric(0),
        weight           = numeric(0)
      ))
    }

    gene_tbl <- purrr::map_dfr(gene_containers, extract_genes_from_container)

    gene_tbl %>%
      distinct(gene_symbol, .keep_all = TRUE) %>%
      mutate(
        axis             = axis_name,
        gene_effect_sign = dplyr::if_else(is.na(gene_effect_sign), 1, gene_effect_sign)
      ) %>%
      select(axis, gene_symbol, gene_effect_sign, weight)
  }
)

message("[INFO] 机制基因集（按 axis × gene_symbol）初步统计：")
axis_gene_tbl %>%
  count(axis, name = "n_genes") %>%
  arrange(axis) %>%
  { print(.) }

# ---------- STEP2.1: 若 YAML 提供 mapping 配置，则使用 symbol -> GeneID 映射 ----------
# 期望 YAML 里有：
# mapping:
#   table: "ref/gene_symbol_to_ensembl.tsv"
#   symbol_column: "gene_symbol"
#   ensembl_column: "ensembl_gene_id"
#
# 且 DE 结果的 GeneID 列使用的是 ensembl_gene_id

if ("mapping" %in% names(gs_conf)) {
  map_conf <- gs_conf$mapping

  required_keys <- c("table", "symbol_column", "ensembl_column")
  missing_keys  <- setdiff(required_keys, names(map_conf))
  if (length(missing_keys) > 0) {
    stop(
      "[FATAL] YAML 中 mapping 配置不完整，缺少字段: ",
      paste(missing_keys, collapse = ", ")
    )
  }

  map_path <- "scripts/rna/rna_mechanistic_gene_sets_gene_registry.tsv"
  sym_col  <- map_conf$symbol_column
  ens_col  <- map_conf$ensembl_column

  if (!file.exists(map_path)) {
    stop("[FATAL] mapping 表不存在: ", map_path)
  }

  map_tbl <- read_tsv(map_path, show_col_types = FALSE)

  if (!all(c(sym_col, ens_col) %in% colnames(map_tbl))) {
    stop(
      "[FATAL] mapping 表中未找到指定列: ",
      sym_col, " 或 ", ens_col,
      "\n  实际列名为: ", paste(colnames(map_tbl), collapse = ", ")
    )
  }

  map_tbl2 <- map_tbl[, c(sym_col, ens_col)]
  colnames(map_tbl2) <- c("gene_symbol", "GeneID")

  axis_gene_tbl <- axis_gene_tbl %>%
    left_join(map_tbl2, by = "gene_symbol")

  n_total  <- nrow(axis_gene_tbl)
  n_mapped <- sum(!is.na(axis_gene_tbl$GeneID))

  message("[INFO] 机制基因映射到 GeneID（Ensembl）:")
  message("  total genes in gene sets = ", n_total)
  message("  mapped (non-NA GeneID)   = ", n_mapped)

  if (n_mapped == 0) {
    stop(
      "[FATAL] 机制基因集中没有任何基因在映射表中找到 GeneID，",
      "请检查 gene_symbol_to_ensembl.tsv 与 YAML 中的 gene_symbol 是否一致。"
    )
  }

  axis_gene_tbl <- axis_gene_tbl %>%
    filter(!is.na(GeneID))

  message("[INFO] 映射后各 axis 可用基因数：")
  axis_gene_tbl %>%
    count(axis, name = "n_genes_mapped") %>%
    arrange(axis) %>%
    { print(.) }

} else {
  message("[WARN] YAML 中未提供 mapping 字段，将假定 GeneID 与 gene_symbol 相同。")
  axis_gene_tbl <- axis_gene_tbl %>%
    rename(GeneID = gene_symbol)
}

# ---------- STEP3: 把 per-gene per-batch LFC 转成长表 ----------
lfc_long <- ic_tbl %>%
  select(GeneID, all_of(lfc_cols)) %>%
  pivot_longer(
    cols      = all_of(lfc_cols),
    names_to  = "batch",
    values_to = "log2FC"
  ) %>%
  mutate(batch = sub("^LFC_", "", batch))

message("[INFO] per-gene LFC long format: n = ", nrow(lfc_long))

# ---------- STEP4: 将基因映射到 axis ----------
axis_lfc_long <- axis_gene_tbl %>%
  inner_join(lfc_long, by = "GeneID") %>%
  mutate(
    gene_effect_sign = dplyr::if_else(is.na(gene_effect_sign), 1, gene_effect_sign),
    weight           = dplyr::if_else(is.na(weight), 1, weight),
    # NOTE: 这里使用 log2FC * gene_effect_sign * weight 作为该基因在该轴上的有效 log2FC，
    # weight 由 YAML 中的 weight 字段给出，未指定时默认为 1。
    effective_log2FC = log2FC * gene_effect_sign * weight
  )

if (nrow(axis_lfc_long) == 0) {
  stop(
    "[FATAL] axis_lfc_long 为空：机制基因集中的 GeneID 与 DE 结果中的 GeneID 完全不匹配。\n",
    "  请检查：1) DE 结果 GeneID 是否为 Ensembl ID；2) mapping 表；3) YAML 中 gene_symbol。"
  )
}

message("[INFO] axis × gene × batch LFC（用于轴级汇总）：")
axis_lfc_long %>%
  count(axis, batch, name = "n_genes_used") %>%
  arrange(axis, batch) %>%
  { print(.) }

# ---------- STEP5: 计算 axis per-batch 汇总值 ----------
axis_scores <- axis_lfc_long %>%
  group_by(axis, batch) %>%
  summarise(
    n_genes    = sum(!is.na(effective_log2FC)),
    mean_LFC   = ifelse(n_genes > 0, mean(effective_log2FC, na.rm = TRUE), NA_real_),
    median_LFC = ifelse(n_genes > 0, median(effective_log2FC, na.rm = TRUE), NA_real_),
    sd_LFC     = ifelse(n_genes > 1, sd(effective_log2FC, na.rm = TRUE), NA_real_),
    se_LFC     = ifelse(n_genes > 1, sd_LFC / sqrt(n_genes), NA_real_),
    z_score    = ifelse(!is.na(se_LFC) & se_LFC > 0, mean_LFC / se_LFC, NA_real_),
    .groups    = "drop"
  ) %>%
  arrange(axis, batch)

out_axis_scores <- file.path(outdir, "rna_mechanistic_axis_scores_per_batch.tsv")
write_tsv(axis_scores, out_axis_scores)
message("[OK] 写出轴级 per-batch 分数表: ", out_axis_scores)

# ---------- STEP6: 写出 per-gene axis LFC 长表 ----------
out_axis_gene <- file.path(outdir, "rna_mechanistic_axis_gene_LFC_long.tsv")
write_tsv(axis_lfc_long, out_axis_gene)
message("[OK] 写出 axis × gene × batch LFC 长表: ", out_axis_gene)

# ---------- SUMMARY ----------
message("============================================================")
message("[SUMMARY]")
message("  机制轴（axes）: ", paste(sort(unique(axis_scores$axis)), collapse = ", "))
message("  批次（batches）: ", paste(sort(unique(axis_scores$batch)), collapse = ", "))
message("  输出表：")
message("    - ", out_axis_scores)
message("    - ", out_axis_gene)
message("============================================================")