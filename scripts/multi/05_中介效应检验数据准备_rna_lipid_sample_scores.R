#!/usr/bin/env Rscript

## ============================================================
## 05_中介效应检验数据准备_rna_lipid_sample_scores.R
##
## 目的：
##   1) 从 RNA featureCounts 矩阵中，为每个样本计算：
##        - M_score_RNA ：指定 pathway 的样本级通路活性
##        - Y_score_RNA ：指定机制轴的样本级轴活性
##   2) （可选）从脂质表达矩阵中，为每个样本计算：
##        - Y_score_lipid ：指定脂质轴/feature 集的样本级表型
##   3) 利用 sample_matching 表，将 RNA 样本名与脂质样本名对齐，
##      输出一张适合中介分析的样本级汇总表。
##
## 输入（默认可被 config 覆盖）：
##   1) RNA counts:
##        data/processed/rna/all_batches_featureCounts.tsv
##      要求至少包含列：
##        GeneID, 以及若干 sample 列（与 sample_matching$rna_sample_id 对应）
##
##   2) 脂质矩阵（已 log2 处理）：
##        data/processed/lipid/lipid_matrix_allfat_log2_filtered.tsv
##      要求至少包含列：
##        feature_id, 以及若干 sample 列（与 sample_matching$lipid_sample_id 对应）
##
##   3) RNA–lipid 样本匹配表：
##        results/multiomics/sample_matching_batch1.tsv
##      需要列：
##        sample_id, group, batch, rna_sample_id, lipid_sample_id
##
##   4) RNA M（中介通路）的 gene set YAML：
##        默认：scripts/multi/ephb1_downstream_signaling_sets.yaml
##      顶层结构：
##        pathways:
##          <pathway_name>:
##            axis: ...
##            genes:
##              - name: GeneSymbol
##                gene_effect_sign: 1/-1
##      需要配套 gene registry：
##        同目录下 <yaml_basename>_gene_registry.tsv
##        包含列：gene_symbol, ensembl_gene_id
##
##   5) RNA Y（机制轴）的 gene set YAML：
##        默认：scripts/rna/rna_mechanistic_gene_sets.yaml
##      顶层结构：
##        axes:
##          <axis_name>:
##            modules:
##              <module_name>:
##                genes:
##                  - name: GeneSymbol
##                    gene_effect_sign: 1/-1
##      同样需要 gene registry：
##        <yaml_basename>_gene_registry.tsv
##
##   6) 脂质 Y 的 feature 集（可选）：
##        默认：scripts/lipid/lipid_supply_axis_feature_set.tsv
##      TSV 列：feature_id, axis, effect_sign
##
##   7) Config 文件（命令行第 1 个参数）：
##        scripts/multiomics/00_sample_scores_config.yaml
##
## 输出：
##   - results/multiomics/sample_scores/sample_scores_for_mediation.tsv
##     列：
##       sample_id, group, batch,
##       M_score_RNA, Y_score_RNA,
##       (可选) Y_score_lipid
##
## ============================================================

suppressPackageStartupMessages({
  library(yaml)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(tools)
  library(tibble)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

## ---------------- 0. 解析 config -----------------------------

args <- commandArgs(trailingOnly = TRUE)

## 1) 必须有一个 config 路径
config_path <- if (length(args) >= 1) {
  args[1]
} else {
  stop("[ERROR] 必须在命令行提供配置文件路径，例如：\n",
       "  Rscript scripts/multi/05_rna_lipid_sample_scores.R scripts/multi/sample_scores_config.yaml\n")
}

if (!file.exists(config_path)) {
  stop("[ERROR] 找不到配置文件: ", config_path)
}

cat("[INFO] 使用配置文件: ", config_path, "\n", sep = "")
cfg <- yaml::read_yaml(config_path)

## 2) 分析模式与所需字段
# analysis mode:
#   - "paired"  : original behavior, uses batch1 paired RNA+lipid
#   - "rna_only": compute RNA sample-level X/M/Z for all RNA samples (no lipid needed)
mode <- cfg$mode %||% "paired"

required_fields <- if (mode == "rna_only") {
  c(
    "rna_counts",
    "ephb1_M_yaml",
    "upstream_M_yaml",
    "outdir"
  )
} else {
  c(
    "rna_counts",
    "lipid_matrix",
    "sample_matching",
    "ephb1_M_yaml",
    "upstream_M_yaml",
    "lipid_feature_set_tsv",      # now used as lipid indicator def YAML
    "lipid_indicator_sign_yaml",  # new: sign matrix YAML
    "outdir"
  )
}

missing_fields <- required_fields[!required_fields %in% names(cfg)]
if (length(missing_fields) > 0) {
  stop("[ERROR] 配置文件中缺少以下字段：\n  ",
       paste(missing_fields, collapse = ", "),
       "\n请在配置文件中补全。")
}

## 3) 不再使用脚本内默认值，全部直接从 cfg 取
rna_counts      <- cfg$rna_counts

# paired mode only
lipid_matrix    <- if (mode == "rna_only") NA_character_ else cfg$lipid_matrix
sample_matching <- if (mode == "rna_only") NA_character_ else cfg$sample_matching

# RNA gene sets
ephb1_M_yaml      <- cfg$ephb1_M_yaml
upstream_M_yaml   <- cfg$upstream_M_yaml

# paired mode only (lipid indicator defs)
lipid_feature_set_tsv      <- if (mode == "rna_only") NA_character_ else cfg$lipid_feature_set_tsv
lipid_indicator_sign_yaml  <- if (mode == "rna_only") NA_character_ else cfg$lipid_indicator_sign_yaml

# optional (rna_only) sample meta: must include sample_id, batch, group (or genotype)
rna_sample_meta_tsv <- cfg$rna_sample_meta_tsv %||% ""

outdir           <- cfg$outdir

cat("============================================================\n")
cat("[INFO] 05_rna_lipid_sample_scores.R\n")
cat("  Mode            : ", mode, "\n", sep = "")
cat("  RNA counts       : ", rna_counts,      "\n", sep = "")
if (mode != "rna_only") {
  cat("  Lipid matrix     : ", lipid_matrix,    "\n", sep = "")
  cat("  Sample matching  : ", sample_matching, "\n", sep = "")
}
cat("  RNA M YAML       : ", ephb1_M_yaml,      "\n", sep = "")
cat("  RNA Y YAML       : ", upstream_M_yaml,      "\n", sep = "")
cat("  RNA M registry    : ", cfg$ephb1_M_registry %||% "(auto *_gene_registry.tsv)", "\n", sep = "")
cat("  RNA Y registry    : ", cfg$upstream_M_registry %||% "(auto *_gene_registry.tsv)", "\n", sep = "")
if (mode != "rna_only") {
  cat("  Lipid indicator def : ", lipid_feature_set_tsv,      "\n", sep = "")
  cat("  Lipid sign YAML     : ", lipid_indicator_sign_yaml,  "\n", sep = "")
}
cat("  Outdir           : ", outdir,          "\n", sep = "")
if (mode == "rna_only") {
  cat("  RNA sample meta  : ", ifelse(nchar(rna_sample_meta_tsv) == 0, "<infer>", rna_sample_meta_tsv), "\n", sep = "")
}
cat("============================================================\n\n")

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
## ---------------- helpers ------------------------------------

## 按行对矩阵做 z-score（每个 feature 在样本间标准化）
## 按行对矩阵做 z-score（每个 feature 在样本间标准化）
## 支持 batch 内标准化：同一基因在每个 batch 内单独做 z-score
row_zscore <- function(mat, batch_vec = NULL) {
  if (is.null(batch_vec)) {
    m <- rowMeans(mat, na.rm = TRUE)
    s <- apply(mat, 1, sd, na.rm = TRUE)
    s[s == 0 | is.na(s)] <- 1
    return((mat - m) / s)
  }

  if (length(batch_vec) != ncol(mat)) {
    stop("[ERROR] row_zscore: batch_vec length must equal ncol(mat).")
  }

  z <- mat
  for (b in unique(batch_vec)) {
    idx <- which(batch_vec == b)
    sub <- mat[, idx, drop = FALSE]
    m <- rowMeans(sub, na.rm = TRUE)
    s <- apply(sub, 1, sd, na.rm = TRUE)
    s[s == 0 | is.na(s)] <- 1
    z[, idx] <- (sub - m) / s
  }
  z
}

## 估计 RNA counts 的 size factors（DESeq2 median-of-ratios 思路的轻量实现）
## 输入：raw counts matrix [genes x samples]
## 输出：named numeric 向量，每个样本一个 size factor
estimate_size_factors_median_ratio <- function(count_mat) {
  # geometric mean per gene (ignore zeros)
  gm <- apply(count_mat, 1, function(x) {
    x <- x[x > 0]
    if (length(x) == 0) return(NA_real_)
    exp(mean(log(x)))
  })

  keep <- !is.na(gm) & gm > 0
  if (sum(keep) < 100) {
    warning("[WARN] 可用于 size factor 的基因过少（<100）。size factor 可能不稳定。")
  }

  ratios <- sweep(count_mat[keep, , drop = FALSE], 1, gm[keep], FUN = "/")
  sf <- apply(ratios, 2, function(r) stats::median(r[r > 0 & is.finite(r)], na.rm = TRUE))

  sf[is.na(sf) | sf <= 0] <- 1
  # normalize to geometric mean 1
  sf <- sf / exp(mean(log(sf)))
  sf
}
## 根据 gene 表 + RNA expr 矩阵计算 sample-level score
## expr_mat: matrix[genes x samples]，行名为 GeneID
## gene_tbl: 必须包含列 GeneID, effect_sign_col
compute_score_rna <- function(expr_mat,
                              gene_tbl,
                              effect_sign_col = "gene_effect_sign",
                              weight_col = "weight",
                              batch_vec = NULL) {
  gene_tbl <- gene_tbl %>%
    filter(!is.na(GeneID)) %>%
    distinct(GeneID, .keep_all = TRUE)

  common_genes <- intersect(rownames(expr_mat), gene_tbl$GeneID)
  if (length(common_genes) == 0) {
    stop("[ERROR] compute_score_rna: 当前 gene 集在表达矩阵中一个基因都没匹配上。")
  }

  expr_sub <- expr_mat[common_genes, , drop = FALSE]
  gene_tbl_sub <- gene_tbl %>%
    filter(GeneID %in% common_genes) %>%
    arrange(match(GeneID, rownames(expr_sub)))

  if (!all(rownames(expr_sub) == gene_tbl_sub$GeneID)) {
    stop("[ERROR] compute_score_rna: GeneID 对齐失败。")
  }

  sign_vec <- as.numeric(gene_tbl_sub[[effect_sign_col]])
  sign_vec[is.na(sign_vec)] <- 1

  w_vec <- if (weight_col %in% colnames(gene_tbl_sub)) {
    as.numeric(gene_tbl_sub[[weight_col]])
  } else {
    rep(1, nrow(gene_tbl_sub))
  }
  w_vec[is.na(w_vec) | w_vec <= 0] <- 1

  # 04 对齐口径：用 log2(norm+1) 的绝对表达做通路 raw score
  # 先乘 sign 与 weight，再对通路内基因做简单平均（而不是除以 sum(weight) 的加权平均）
  eff_expr <- sweep(expr_sub, 1, sign_vec * w_vec, FUN = "*")
  raw_score <- colMeans(eff_expr, na.rm = TRUE)

  # 04 对齐口径：对每条通路的 raw_score 做 batch 内 z-score（通路层面标准化）
  if (!is.null(batch_vec)) {
    if (length(batch_vec) != length(raw_score)) {
      stop("[ERROR] compute_score_rna: batch_vec length must equal ncol(expr_mat).")
    }
    z_score <- raw_score
    for (b in unique(batch_vec)) {
      idx <- which(batch_vec == b)
      m <- mean(raw_score[idx], na.rm = TRUE)
      s <- stats::sd(raw_score[idx], na.rm = TRUE)
      if (is.na(s) || s == 0) s <- 1
      z_score[idx] <- (raw_score[idx] - m) / s
    }
    return(z_score)
  }

  raw_score
}

## 根据 feature set + lipid 矩阵计算 sample-level score
## lipid_mat: matrix[features x samples]，行名为 feature_id
## feat_tbl: 必须包含列 feature_id, effect_sign
compute_score_lipid <- function(lipid_mat, feat_tbl) {
  feat_tbl <- feat_tbl %>%
    filter(!is.na(feature_id)) %>%
    distinct(feature_id, .keep_all = TRUE)

  common_feats <- intersect(rownames(lipid_mat), feat_tbl$feature_id)
  if (length(common_feats) == 0) {
    stop("[ERROR] compute_score_lipid: 当前 feature 集在脂质矩阵中一个 feature 都没匹配上。")
  }

  lipid_sub <- lipid_mat[common_feats, , drop = FALSE]
  feat_tbl_sub <- feat_tbl %>%
    filter(feature_id %in% common_feats) %>%
    arrange(match(feature_id, rownames(lipid_sub)))

  if (!all(rownames(lipid_sub) == feat_tbl_sub$feature_id)) {
    stop("[ERROR] compute_score_lipid: feature_id 对齐失败。")
  }

  sign_vec <- as.numeric(feat_tbl_sub$effect_sign)
  sign_vec[is.na(sign_vec)] <- 1

  z_mat <- row_zscore(lipid_sub)
  eff_mat <- z_mat * sign_vec

  score <- colMeans(eff_mat, na.rm = TRUE)
  score
}

## ---------------- 1. 读取 sample_matching / 构造 RNA-only meta -------------------

if (mode != "rna_only") {
  if (!file.exists(sample_matching)) {
    stop("[ERROR] 找不到 sample_matching 表: ", sample_matching)
  }

  cat("[STEP1] 读取样本匹配表: ", sample_matching, "\n", sep = "")
  sm <- readr::read_tsv(sample_matching, show_col_types = FALSE)

  # If only one column due to delimiter issues, retry with whitespace
  if (length(colnames(sm)) == 1 && grepl("\\s+", colnames(sm)[1])) {
    warning("[WARN] 检测到 sample_matching 只包含一个列名，疑似整行未被拆分，将使用 read_table() 以空白分隔重新读取。")
    sm <- readr::read_table(sample_matching, show_col_types = FALSE)
  }

  colnames(sm) <- gsub("\ufeff", "", colnames(sm))
  colnames(sm) <- trimws(colnames(sm))

  required_sm_cols <- c("sample_id", "group", "batch", "rna_sample_id", "lipid_sample_id")
  if (!all(required_sm_cols %in% colnames(sm))) {
    stop("[ERROR] sample_matching 必须包含列: ",
         paste(required_sm_cols, collapse = ", "),
         "\n实际列名: ", paste(colnames(sm), collapse = ", "))
  }

  cat("  [sample_matching] nrow = ", nrow(sm), "\n", sep = "")
  # paired-mode: only use batch1
  sm <- sm %>% dplyr::filter(.data$batch == "batch1")
  cat("  [sample_matching] after filter batch1, nrow = ", nrow(sm), "\n", sep = "")
  if (nrow(sm) == 0) {
    stop("[ERROR] sample_matching 过滤 batch1 后没有样本，请检查 batch 列是否为 batch1。")
  }

} else {
  cat("[STEP1] RNA-only 模式：不读取 sample_matching，直接从 counts/metadata 构造样本信息。\n")
  sm <- NULL
}
## ---------------- 2. 读取 RNA counts 并构建 expr_mat -------

if (!file.exists(rna_counts)) {
  stop("[ERROR] 找不到 RNA counts 表: ", rna_counts)
}

cat("\n[STEP2] 读取 RNA counts 矩阵: ", rna_counts, "\n", sep = "")
rna_counts_tbl <- readr::read_tsv(rna_counts, show_col_types = FALSE)

if (!"GeneID" %in% colnames(rna_counts_tbl)) {
  stop("[ERROR] RNA counts 表中缺少 GeneID 列。")
}

if (mode == "rna_only") {
  # all sample columns except GeneID
  rna_samples <- setdiff(colnames(rna_counts_tbl), "GeneID")
} else {
  rna_samples <- unique(sm$rna_sample_id)
}

missing_rna_samples <- setdiff(rna_samples, colnames(rna_counts_tbl))
if (length(missing_rna_samples) > 0) {
  stop("[ERROR] 下列 RNA sample 列在 counts 表中找不到: ",
       paste(missing_rna_samples, collapse = ", "))
}

rna_expr_log2 <- rna_counts_tbl %>%
  dplyr::select(GeneID, dplyr::all_of(rna_samples)) %>%
  column_to_rownames("GeneID") %>%
  as.matrix()
# batch 向量（与 rna_expr_log2 的列顺序一致）
if (mode != "rna_only") {
  rna_sample_batch <- sm %>%
    dplyr::select(rna_sample_id, batch) %>%
    distinct() %>%
    dplyr::filter(rna_sample_id %in% colnames(rna_expr_log2)) %>%
    dplyr::arrange(match(rna_sample_id, colnames(rna_expr_log2)))

  if (!all(rna_sample_batch$rna_sample_id == colnames(rna_expr_log2))) {
    stop("[ERROR] rna_sample_batch 与 rna_expr_log2 列对齐失败。")
  }
  rna_batch_vec <- rna_sample_batch$batch

} else {
  # build meta from optional file or infer from sample names
  if (nchar(rna_sample_meta_tsv) > 0 && file.exists(rna_sample_meta_tsv)) {
    meta <- readr::read_tsv(rna_sample_meta_tsv, show_col_types = FALSE)
    # accept either group or genotype
    if (!all(c("sample_id", "batch") %in% colnames(meta))) {
      stop("[ERROR] rna_sample_meta_tsv 必须至少包含列: sample_id, batch。")
    }
    if (!"group" %in% colnames(meta) && "genotype" %in% colnames(meta)) {
      meta <- meta %>% dplyr::rename(group = genotype)
    }
    if (!"group" %in% colnames(meta)) {
      meta$group <- NA_character_
    }
    meta <- meta %>% dplyr::filter(sample_id %in% colnames(rna_expr_log2))
  } else {
    meta <- tibble::tibble(sample_id = colnames(rna_expr_log2)) %>%
      dplyr::mutate(
        batch = stringr::str_extract(sample_id, "batch[0-9]+"),
        group = dplyr::case_when(
          stringr::str_detect(sample_id, "WT") ~ "WT",
          stringr::str_detect(sample_id, "HO") ~ "HO",
          stringr::str_detect(sample_id, "KO") ~ "HO",
          TRUE ~ NA_character_
        )
      )
  }

  if (any(is.na(meta$batch))) {
    warning("[WARN] RNA-only: 有样本未能推断 batch（batch\\d+），将赋值为 'batch1'（仅用于批内 z-score）。")
    meta$batch[is.na(meta$batch)] <- "batch1"
  }

  meta <- meta %>% dplyr::arrange(match(sample_id, colnames(rna_expr_log2)))
  rna_batch_vec <- meta$batch
  rna_only_meta <- meta
}

 # 使用 size factor 对 raw counts 做归一化，避免 KO/WT library size 差异导致 sample-level 分数翻号
sf <- estimate_size_factors_median_ratio(rna_expr_log2)
rna_expr_log2 <- sweep(rna_expr_log2, 2, sf, FUN = "/")
rna_expr_log2 <- log2(rna_expr_log2 + 1)

cat("  [rna_expr_log2] genes = ", nrow(rna_expr_log2),
    " ; samples = ", ncol(rna_expr_log2), "\n", sep = "")

## ---------------- 3. RNA M_score: pathway YAML --------------

cat("\n[STEP3] 计算 RNA M_score（通路级，中介变量）...\n")

if (!file.exists(ephb1_M_yaml)) {
  stop("[ERROR] 找不到 RNA M YAML: ", ephb1_M_yaml)
}

M_cfg <- yaml::read_yaml(ephb1_M_yaml)
if (is.null(M_cfg$pathways)) {
  stop("[ERROR] RNA M YAML 中未找到 'pathways' 顶层字段。")
}

M_pathways <- M_cfg$pathways
if (length(M_pathways) == 0) {
  stop("[ERROR] RNA M YAML 中 'pathways' 为空，至少需要一个通路。")
}

# 将所有 EphB1 下游信号通路的基因集展开为长表，保留 pathway 名称用于列名
M_gene_tbl_sym <- purrr::imap_dfr(
  M_pathways,
  function(pw, pw_name) {
    genes <- pw$genes
    if (is.null(genes)) {
      return(tibble::tibble(
        pathway          = character(0),
        gene_symbol      = character(0),
        gene_effect_sign = numeric(0)
      ))
    }
    purrr::map_dfr(genes, function(g) {
      tibble::tibble(
        pathway          = pw_name,
        gene_symbol      = g$name,
        gene_effect_sign = as.numeric(g$gene_effect_sign %||% 1),
        weight           = as.numeric(g$weight %||% 1)
      )
    })
  }
)

yaml_dir_M  <- dirname(ephb1_M_yaml)
yaml_base_M <- tools::file_path_sans_ext(basename(ephb1_M_yaml))
default_registry_M  <- file.path(yaml_dir_M, paste0(yaml_base_M, "_gene_registry.tsv"))
registry_M <- cfg$ephb1_M_registry %||% default_registry_M

if (!file.exists(registry_M)) {
  stop("[ERROR] 找不到 RNA M gene registry: ", registry_M,
       "\n请在配置文件中通过 ephb1_M_registry 指定，或在 YAML 同目录下生成 *_gene_registry.tsv。")
}

reg_M <- readr::read_tsv(registry_M, show_col_types = FALSE)
if (!all(c("gene_symbol", "ensembl_gene_id") %in% colnames(reg_M))) {
  stop("[ERROR] RNA M gene registry 必须包含 gene_symbol, ensembl_gene_id。")
}

M_gene_tbl <- M_gene_tbl_sym %>%
  left_join(reg_M %>% dplyr::select(gene_symbol, ensembl_gene_id),
            by = "gene_symbol") %>%
  rename(GeneID = ensembl_gene_id)

cat("  [RNA M gene set] n_pathway = ", length(unique(M_gene_tbl$pathway)),
    " ; total n_gene = ", nrow(M_gene_tbl),
    " ; matched Ensembl = ", sum(!is.na(M_gene_tbl$GeneID)), "\n", sep = "")

# 对每个 pathway 计算一个样本级 M_score 列，列名直接使用 YAML 中的 pathway 名称
M_scores_list <- M_gene_tbl %>%
  dplyr::filter(!is.na(GeneID)) %>%
  split(.$pathway) %>%
  purrr::map(~ compute_score_rna(rna_expr_log2, .x,
                  effect_sign_col = "gene_effect_sign",
                  weight_col = "weight",
                  batch_vec = rna_batch_vec))

M_score_df <- tibble::tibble(rna_sample_id = colnames(rna_expr_log2))
for (pw_name in names(M_scores_list)) {
  vec <- M_scores_list[[pw_name]]
  M_score_df[[pw_name]] <- vec[M_score_df$rna_sample_id]
}

## ---------------- 4. RNA Y_score: mechanistic axis ----------

cat("\n[STEP4] 计算 RNA Y_score（机制轴级）...\n")

if (!file.exists(upstream_M_yaml)) {
  stop("[ERROR] 找不到 RNA Y YAML: ", upstream_M_yaml)
}

Y_cfg <- yaml::read_yaml(upstream_M_yaml)

# 支持两种结构：
# 1) 旧的 mechanistic 轴：axes -> modules -> genes
# 2) 新的上游调控通路：pathways -> genes
if (!is.null(Y_cfg$axes)) {
  Y_gene_tbl_sym <- purrr::imap_dfr(
    Y_cfg$axes,
    function(axis_obj, axis_name) {
      modules <- axis_obj$modules %||% list()
      purrr::imap_dfr(
        modules,
        function(mod, mod_name) {
          genes <- mod$genes
          if (is.null(genes)) {
            return(tibble::tibble(
              axis             = character(0),
              module           = character(0),
              gene_symbol      = character(0),
              gene_effect_sign = numeric(0),
              weight           = numeric(0)
            ))
          }
          purrr::map_dfr(genes, function(g) {
            tibble::tibble(
              axis             = axis_name,
              module           = mod_name,
              gene_symbol      = g$name,
              gene_effect_sign = as.numeric(g$gene_effect_sign %||% 1),
              weight           = as.numeric(g$weight %||% 1)
            )
          })
        }
      )
    }
  )
} else if (!is.null(Y_cfg$pathways)) {
  Y_gene_tbl_sym <- purrr::imap_dfr(
    Y_cfg$pathways,
    function(pw, pw_name) {
      genes <- pw$genes
      if (is.null(genes)) {
        return(tibble::tibble(
          axis             = character(0),
          module           = character(0),
          gene_symbol      = character(0),
          gene_effect_sign = numeric(0),
          weight           = numeric(0)
        ))
      }
      purrr::map_dfr(genes, function(g) {
        tibble::tibble(
          axis             = pw_name,
          module           = NA_character_,
          gene_symbol      = g$name,
          gene_effect_sign = as.numeric(g$gene_effect_sign %||% 1),
          weight           = as.numeric(g$weight %||% 1)
        )
      })
    }
  )
} else {
  stop("[ERROR] RNA Y YAML 中既没有 'axes' 也没有 'pathways' 顶层字段，无法解析。")
}

yaml_dir_Y  <- dirname(upstream_M_yaml)
yaml_base_Y <- tools::file_path_sans_ext(basename(upstream_M_yaml))
default_registry_Y  <- file.path(yaml_dir_Y, paste0(yaml_base_Y, "_gene_registry.tsv"))
registry_Y <- cfg$upstream_M_registry %||% default_registry_Y

if (!file.exists(registry_Y)) {
  stop("[ERROR] 找不到 RNA Y gene registry: ", registry_Y,
       "\n请在配置文件中通过 upstream_M_registry 指定，或为对应 YAML 生成 *_gene_registry.tsv。")
}

reg_Y <- readr::read_tsv(registry_Y, show_col_types = FALSE)
if (!all(c("gene_symbol", "ensembl_gene_id") %in% colnames(reg_Y))) {
  stop("[ERROR] RNA Y gene registry 必须包含 gene_symbol, ensembl_gene_id。")
}

Y_gene_tbl <- Y_gene_tbl_sym %>%
  left_join(reg_Y %>% dplyr::select(gene_symbol, ensembl_gene_id),
            by = "gene_symbol") %>%
  rename(GeneID = ensembl_gene_id)

cat("  [RNA Y axis/pathway gene set] n_axis/pathway = ", length(unique(Y_gene_tbl$axis)),
    " ; total n_gene = ", nrow(Y_gene_tbl),
    " ; matched Ensembl = ", sum(!is.na(Y_gene_tbl$GeneID)), "\n", sep = "")

Y_scores_list <- Y_gene_tbl %>%
  dplyr::filter(!is.na(GeneID)) %>%
  split(.$axis) %>%
  purrr::map(~ compute_score_rna(rna_expr_log2, .x,
                  effect_sign_col = "gene_effect_sign",
                  weight_col = "weight",
                  batch_vec = rna_batch_vec))

Y_score_df <- tibble::tibble(rna_sample_id = colnames(rna_expr_log2))
for (axis_name in names(Y_scores_list)) {
  vec <- Y_scores_list[[axis_name]]
  Y_score_df[[axis_name]] <- vec[Y_score_df$rna_sample_id]
}

# Initialize lipid score flag
has_lipid_scores <- FALSE

## ---- STEP5: 计算 / 读取脂质 Y_score_lipid ---------------------------------

# RNA-only 模式下不应触碰任何 lipid 逻辑（sm 为 NULL），直接跳过
if (mode == "rna_only") {
  cat("\n[STEP5] RNA-only：跳过脂质 Y_score_lipid 读取/计算。\n")
  precomp_lipid_path <- NULL
} else {
  precomp_lipid_path <- cfg$precomputed_lipid_axis_scores
}

if (mode != "rna_only" && !is.null(precomp_lipid_path) && file.exists(precomp_lipid_path)) {

  message("[STEP5] 使用预计算好的脂质轴分数表: ", precomp_lipid_path)

  lipid_axis_df <- readr::read_tsv(precomp_lipid_path, show_col_types = FALSE)

  # 用户希望额外作为 Y 的脂质/指标列（原始列名）
  extra_y_raw <- c("sum_CL", "CL_fraction", "mito_lipid_mass",
                   "PC", "PE", "CL")

  extra_y_cols <- intersect(extra_y_raw, colnames(lipid_axis_df))
  if (length(extra_y_cols) > 0) {
    message("[INFO] 在预计算脂质表中发现以下将作为额外 Y 使用的列: ",
            paste(extra_y_cols, collapse = ", "))

    # 为了让后续中介脚本将其识别为 Y（outcome），统一加上 `_score` 后缀
    extra_rename <- stats::setNames(
      extra_y_cols,
      paste0(extra_y_cols, "_score")
    )

    lipid_axis_df <- lipid_axis_df %>%
      dplyr::rename(!!!extra_rename)
  } else {
    message("[INFO] 在预计算脂质表中未找到指定的额外 Y 列（sum_CL, CL_fraction, mito_lipid_mass, PC, PE, CL），",
            "将只使用已有的 *_score 轴分数列。")
  }

  # 自动识别所有 *_score 列，一次性 join 进来，作为 CL 相关多轴脂质表型
  score_cols <- grep("_score$", colnames(lipid_axis_df), value = TRUE)

  if (length(score_cols) == 0) {
    stop("[ERROR] 在 ", precomp_lipid_path, " 中找不到任何 *_score 列，无法作为脂质轴分数使用。")
  }

  message("  [INFO] 将以下脂质轴分数列 join 到 sample_matching: ",
          paste(score_cols, collapse = ", "))

  Y_score_lipid_df <- lipid_axis_df %>%
    dplyr::mutate(
      lipid_sample_id = if ("batch" %in% colnames(.)) {
        paste0(.data$batch, "_", .data$sample_id)
      } else {
        .data$sample_id
      }
    ) %>%
    dplyr::select(lipid_sample_id, dplyr::all_of(score_cols))

  sm <- sm %>%
    dplyr::left_join(Y_score_lipid_df, by = "lipid_sample_id")

  has_lipid_scores <- TRUE


} else if (mode != "rna_only") {

  # 没有预计算表，再看要不要走旧的“从 lipid_matrix 重算”逻辑
  lipid_matrix_path        <- cfg$lipid_matrix
  lipid_indicator_def_yaml <- cfg$lipid_feature_set_tsv

  if (is.null(lipid_matrix_path) || !file.exists(lipid_matrix_path) ||
      is.null(lipid_indicator_def_yaml) || !file.exists(lipid_indicator_def_yaml)) {

    message("[STEP5] 未提供预计算脂质轴分数，也未提供完整脂质矩阵 / 指标定义 YAML，跳过 Y_score_lipid 计算。")

  } else {

    message("[STEP5] 未提供预计算脂质轴分数，改为根据脂质矩阵重算 Y_score_lipid ...")

    ## 这里保留你原来的 Step5 计算代码（读 lipid_matrix、构造 indicator_mat、z-score、bottleneck_sign 等）
    ## ...
  }
}
#
# ---------------- RNA-only output and early exit ----------------
if (mode == "rna_only") {
  cat("\n[STEP6] RNA-only：汇总样本级 X/M/Z 并写出 sample_scores_rna_only.tsv ...\n")

  score_tbl_rna_only <- tibble::tibble(
    sample_id     = colnames(rna_expr_log2),
    rna_sample_id = colnames(rna_expr_log2)
  ) %>%
    dplyr::left_join(M_score_df, by = "rna_sample_id") %>%
    dplyr::left_join(Y_score_df, by = "rna_sample_id")

  # attach batch/group if available
  if (exists("rna_only_meta")) {
    score_tbl_rna_only <- score_tbl_rna_only %>%
      dplyr::left_join(rna_only_meta %>% dplyr::rename(rna_sample_id = sample_id), by = "rna_sample_id") %>%
      dplyr::select(sample_id, group, batch, rna_sample_id, dplyr::everything())
  }

  out_path_rna_only <- file.path(outdir, "sample_scores_rna_only.tsv")
  readr::write_tsv(score_tbl_rna_only, out_path_rna_only)

  cat("  [OK] 写出 RNA-only 样本级宽表: ", out_path_rna_only, "\n", sep = "")
  cat("============================================================\n")
  cat("[DONE] 05_rna_lipid_sample_scores.R (mode=rna_only) 完成。\n")
  cat("============================================================\n")
  quit(save = "no", status = 0)
}

## ---------------- 6. 汇总样本级得分并写出 -------------------

cat("\n[STEP6] 汇总样本级得分并写出...\n")

score_tbl <- sm %>%
  dplyr::mutate(
    group = group,
    batch = batch
  ) %>%
  dplyr::left_join(M_score_df, by = "rna_sample_id") %>%
  dplyr::left_join(Y_score_df, by = "rna_sample_id")

out_path <- file.path(outdir, "sample_scores_for_mediation.tsv")
readr::write_tsv(score_tbl, out_path)

cat("  [OK] 写出样本级得分表: ", out_path, "\n", sep = "")
cat("  列包括：sample_id, group, batch, rna_sample_id, lipid_sample_id，\n")
cat("          多列 RNA M 候选通路分数、RNA Y 候选轴/通路分数",
    if (has_lipid_scores) "，以及一组脂质轴得分列。\n" else "。\n", sep = "")

cat("============================================================\n")
cat("[DONE] 05_rna_lipid_sample_scores.R 完成。\n")
cat("  下一步：可将此表作为 08_mediation_analysis_minimal.R 的输入，\n")
cat("          选择适当的 M_score / Y_score 列进行中介分析。\n")
cat("============================================================\n")