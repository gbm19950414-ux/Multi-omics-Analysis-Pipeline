#!/usr/bin/env Rscript

## ============================================================
## 06_rna_path_rank_XMZ.R
##
## 目的（RNA-only Figure6 主计算）：
##   输入 sample_scores_rna_only.tsv（跨批次 RNA 样本级 X/M/Z 宽表），
##   计算并排序候选 EphB1 下游通路（M）：
##     A) X→M：各 batch 的 ΔM(HO−WT)、方向投票、汇总效应
##     B) M↔Z：各 batch 的耦合（默认 Z ~ M + X），以及 LOO 稳健性
##     C) 综合排序：一致性 + 耦合强度 + LOO
##
## 输入：
##   1) sample_scores_rna_only.tsv
##      必须列：sample_id, batch, group（WT/HO）
##      必须有 Z 列（你指定），以及多个 M 列（你指定或自动推断）
##
## 输出（outdir）：
##   - rank_X_to_M.tsv
##   - rank_M_to_Z.tsv
##   - rank_combined_M_for_Z.tsv
##   - loo_matrix_M_to_Z.tsv   (每个 M × 每个 leave-one-out 样本的系数)
##
## 用法：
##   Rscript scripts/multi/06_rna_path_rank_XMZ.R \
##     results/multiomics/tables/sample_scores/sample_scores_rna_only.tsv \
##     results/multiomics/tables/mediation \
##     "Membrane context" \
##     ""          # 可选：用逗号指定 M 列名；留空自动识别
##
##   # 传入 "__ALL__" 可对所有机制轴逐一执行
##
## ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(tibble)
})

args <- commandArgs(trailingOnly = TRUE)

# ---- defaults: allow running with no args ----
def_in_path <- "results/multiomics/tables/sample_scores/sample_scores_rna_only.tsv"
def_outdir  <- "results/multiomics/tables/mediation"
def_Z_col   <- "__ALL__"  # 默认对所有机制轴逐一执行

def_usage <- paste0(
  "[INFO] 未提供参数，使用默认值运行：\n",
  "  input : ", def_in_path, "\n",
  "  outdir: ", def_outdir,  "\n",
  "  Z_col : ", def_Z_col,   "\n",
  "  M_cols: <auto>\n",
  "如需自定义：Rscript 06_rna_path_rank_XMZ.R <sample_scores_rna_only.tsv> <outdir> <Z_col> [M_cols_csv]\n"
)

in_path <- def_in_path
outdir  <- def_outdir
Z_col   <- def_Z_col
M_cols_csv <- ""

if (length(args) == 0) {
  cat(def_usage)
} else {
  if (length(args) < 3) {
    stop("[ERROR] 用法：Rscript 06_rna_path_rank_XMZ.R <sample_scores_rna_only.tsv> <outdir> <Z_col> [M_cols_csv]\n")
  }
  in_path <- args[1]
  outdir  <- args[2]
  Z_col   <- args[3]
  M_cols_csv <- if (length(args) >= 4) args[4] else ""
}

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

cat("============================================================\n")
cat("[INFO] 06_rna_path_rank_XMZ.R\n")
cat("  input   : ", in_path, "\n", sep = "")
cat("  outdir  : ", outdir,  "\n", sep = "")
cat("  Z_col   : ", Z_col,   "\n", sep = "")
cat("  M_cols  : ", ifelse(nchar(M_cols_csv)==0, "<auto>", M_cols_csv), "\n", sep = "")
cat("============================================================\n\n")

if (!file.exists(in_path)) stop("[ERROR] 找不到输入表: ", in_path)
df <- readr::read_tsv(in_path, show_col_types = FALSE)

# ---- required columns ----
req <- c("sample_id", "batch", "group")
miss <- setdiff(req, colnames(df))
if (length(miss) > 0) stop("[ERROR] 输入表缺少列: ", paste(miss, collapse=", "))

# group normalize: WT / HO
df <- df %>%
  mutate(
    group = case_when(
      str_to_upper(group) %in% c("WT") ~ "WT",
      str_to_upper(group) %in% c("HO","KO") ~ "HO",
      TRUE ~ group
    )
  )

 # ---- helper: detect numeric-like columns robustly (avoid list subscript issues) ----
is_numeric_like <- function(x) {
  if (is.list(x)) return(FALSE)
  if (is.numeric(x)) return(TRUE)
  xx <- suppressWarnings(as.numeric(x))
  if (length(xx) == 0) return(FALSE)
  # numeric-like if not all NA after coercion
  !all(is.na(xx))
}

# ---- determine M columns ----
M_cols <- if (nchar(M_cols_csv) > 0) {
  str_split(M_cols_csv, ",", simplify = TRUE) %>% as.character() %>% trimws() %>% .[. != ""]
} else {
  # 自动识别：默认把 EphB1_ 前缀列视为 M（下游通路集合）
  # 注意：当 Z_col=="__ALL__" 时，不能把所有数值列都归为 M，否则 Z 会被排空。
  exclude <- c("sample_id","rna_sample_id","batch","group")
  cand <- setdiff(colnames(df), exclude)
  cand <- cand[vapply(df[cand], is_numeric_like, logical(1))]

  if (Z_col == "__ALL__") {
    cand <- cand[grepl("^EphB1_", cand)]
  } else {
    # 单一 Z 模式：排除指定的 Z 列，其余数值列默认当作 M
    cand <- setdiff(cand, Z_col)
  }
  cand
}
M_cols <- intersect(M_cols, colnames(df))
if (length(M_cols) == 0) stop("[ERROR] 未能识别任何 M 列。请在第4参数用逗号指定 M 列名。")

# ---- identify all Z columns if Z_col == "__ALL__" ----
all_Z_cols <- NULL
if (Z_col == "__ALL__") {
  # 候选 Z 列 = 所有数值列
  exclude_z <- c("sample_id","rna_sample_id","batch","group", M_cols)
  candidate_z <- setdiff(colnames(df), exclude_z)
  all_Z_cols <- candidate_z[vapply(df[candidate_z], is_numeric_like, logical(1))]
  if (length(all_Z_cols) == 0) stop("[ERROR] 未能识别任何 Z 列。请指定一个有效的 Z_col，或使用 __ALL__ 且数据中有数值列。")
}

cat("[INFO] n_samples = ", nrow(df), " ; n_batches = ", n_distinct(df$batch), " ; n_M = ", length(M_cols), "\n", sep="")

# ---- helper: delta (HO-WT) in a batch ----
delta_ho_wt <- function(x, g) {
  mu_ho <- mean(x[g=="HO"], na.rm = TRUE)
  mu_wt <- mean(x[g=="WT"], na.rm = TRUE)
  mu_ho - mu_wt
}

# ---- main function to run for one Z ----
run_for_one_Z <- function(Z_col_in) {
  if (!(Z_col_in %in% colnames(df))) {
    stop("[ERROR] 指定的 Z_col 不存在：", Z_col_in, "\n可用列：", paste(colnames(df), collapse=", "))
  }
  cat("\n[INFO] 处理中间变量 Z = ", Z_col_in, "\n", sep = "")

  # ---- A) X -> M : per-batch delta ----
  cat("\n[STEP1] X→M: 计算每个 batch 的 ΔM(HO−WT)...\n")
  rank_X_to_M_long <- df %>%
    select(sample_id, batch, group, all_of(M_cols)) %>%
    pivot_longer(all_of(M_cols), names_to = "M", values_to = "M_score") %>%
    group_by(M, batch) %>%
    summarise(
      n_WT = sum(group=="WT", na.rm=TRUE),
      n_HO = sum(group=="HO", na.rm=TRUE),
      delta_M = delta_ho_wt(M_score, group),
      .groups = "drop"
    )

  rank_X_to_M <- rank_X_to_M_long %>%
    group_by(M) %>%
    summarise(
      n_batches = n(),
      # 方向投票：delta_M > 0 视作 HO 更高
      n_pos = sum(delta_M > 0, na.rm=TRUE),
      n_neg = sum(delta_M < 0, na.rm=TRUE),
      frac_pos = n_pos / n_batches,
      # 汇总效应量：中位数更稳健
      delta_M_median = median(delta_M, na.rm=TRUE),
      delta_M_mean   = mean(delta_M, na.rm=TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(abs(delta_M_median)))

  # ---- B) M <-> Z coupling within batch: Z ~ M + X ----
  cat("\n[STEP2] M↔Z: 在每个 batch 内拟合 Z ~ M + X，并汇总...\n")

  fit_batch <- function(dat) {
    # dat: a single batch, expected columns: group, Z, M
    if (nrow(dat) < 4) {
      return(tibble(beta_M = NA_real_, p_M = NA_real_, beta_X = NA_real_, p_X = NA_real_, n = nrow(dat)))
    }

    d <- dat %>%
      dplyr::select(group, Z, M) %>%
      dplyr::mutate(group = factor(group, levels = c("WT", "HO")))

    d$Z <- as.numeric(d$Z)
    d$M <- as.numeric(d$M)

    mod <- tryCatch(stats::lm(Z ~ M + group, data = d), error = function(e) NULL)
    if (is.null(mod)) {
      return(tibble(beta_M = NA_real_, p_M = NA_real_, beta_X = NA_real_, p_X = NA_real_, n = nrow(d)))
    }

    sm <- summary(mod)$coefficients
    beta_M <- if ("M" %in% rownames(sm)) sm["M", "Estimate"] else NA_real_
    p_M    <- if ("M" %in% rownames(sm)) sm["M", "Pr(>|t|)"] else NA_real_
    beta_X <- if ("groupHO" %in% rownames(sm)) sm["groupHO", "Estimate"] else NA_real_
    p_X    <- if ("groupHO" %in% rownames(sm)) sm["groupHO", "Pr(>|t|)"] else NA_real_

    tibble(beta_M = beta_M, p_M = p_M, beta_X = beta_X, p_X = p_X, n = nrow(d))
  }

  rank_M_to_Z_long <- map_dfr(M_cols, function(mn) {
    df %>%
      select(sample_id, batch, group, Z = all_of(Z_col_in), M = all_of(mn)) %>%
      group_by(batch) %>%
      group_modify(~{
        fit_batch(.x)
      }) %>%
      ungroup() %>%
      mutate(M = mn) %>%
      select(M, batch, n, beta_M, p_M, beta_X, p_X)
  })

  rank_M_to_Z <- rank_M_to_Z_long %>%
    group_by(M) %>%
    summarise(
      n_batches = n(),
      n_beta_pos = sum(beta_M > 0, na.rm=TRUE),
      n_beta_neg = sum(beta_M < 0, na.rm=TRUE),
      frac_beta_pos = n_beta_pos / n_batches,
      beta_M_median = median(beta_M, na.rm=TRUE),
      beta_M_mean   = mean(beta_M, na.rm=TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(abs(beta_M_median)))

  # ---- C) LOO robustness on pooled data (not per-batch) ----
  # rationale: n per batch may be small; we use pooled LOO as a stability check,
  # and still report per-batch summaries above.
  cat("\n[STEP3] LOO: pooled 数据做 leave-one-out，评估 beta_M 稳健性...\n")

  pooled_fit_betaM <- function(dat) {
    # dat expected columns: group, Z, M
    d <- dat %>%
      dplyr::select(group, Z, M) %>%
      dplyr::mutate(group = factor(group, levels = c("WT", "HO")))

    d$Z <- as.numeric(d$Z)
    d$M <- as.numeric(d$M)

    mod <- tryCatch(stats::lm(Z ~ M + group, data = d), error = function(e) NULL)
    if (is.null(mod)) return(NA_real_)

    cf <- coef(mod)
    if (!"M" %in% names(cf)) NA_real_ else as.numeric(cf[["M"]])
  }

  loo_tbl <- map_dfr(M_cols, function(mn) {
    full_beta <- pooled_fit_betaM(df %>% dplyr::transmute(group, Z = .data[[Z_col_in]], M = .data[[mn]]))
    betas <- map_dbl(df$sample_id, function(sid) {
      d2 <- df %>% dplyr::filter(sample_id != sid) %>% dplyr::transmute(group, Z = .data[[Z_col_in]], M = .data[[mn]])
      pooled_fit_betaM(d2)
    })
    tibble(
      M = mn,
      sample_left_out = df$sample_id,
      beta_M_loo = betas,
      beta_M_full = full_beta
    )
  })

  # ---- file name suffix ----
  suffix <- ifelse(Z_col == "__ALL__", paste0("__", Z_col_in), "")

  readr::write_tsv(rank_X_to_M, file.path(outdir, paste0("rank_X_to_M", suffix, ".tsv")))
  cat("  [OK] 写出 rank_X_to_M", suffix, ".tsv\n", sep = "")

  readr::write_tsv(rank_M_to_Z, file.path(outdir, paste0("rank_M_to_Z", suffix, ".tsv")))
  cat("  [OK] 写出 rank_M_to_Z", suffix, ".tsv\n", sep = "")

  readr::write_tsv(loo_tbl, file.path(outdir, paste0("loo_matrix_M_to_Z", suffix, ".tsv")))
  cat("  [OK] 写出 loo_matrix_M_to_Z", suffix, ".tsv\n", sep = "")

  loo_summary <- loo_tbl %>%
    group_by(M) %>%
    summarise(
      full = unique(beta_M_full)[1],
      n_ok = sum(sign(beta_M_loo) == sign(full) & !is.na(beta_M_loo) & !is.na(full)),
      n_total = sum(!is.na(beta_M_loo)) ,
      frac_ok = ifelse(n_total==0, NA_real_, n_ok/n_total),
      beta_M_loo_min = min(beta_M_loo, na.rm=TRUE),
      beta_M_loo_max = max(beta_M_loo, na.rm=TRUE),
      .groups = "drop"
    )

  # ---- combined rank ----
  cat("\n[STEP4] 合成排序：一致性(X→M) + 耦合(M↔Z) + LOO...\n")

  rank_combined <- rank_X_to_M %>%
    rename(delta_M = delta_M_median, frac_XM_pos = frac_pos) %>%
    left_join(rank_M_to_Z %>% rename(beta_M = beta_M_median, frac_MZ_pos = frac_beta_pos), by="M") %>%
    left_join(loo_summary, by="M") %>%
    mutate(
      # 你可以按 Figure6 想要的“保护性”定义来改符号：
      # 这里给一个默认：希望 HO 导致 M 改变（delta_M）与 M→Z 耦合（beta_M）在方向上可解释。
      # 如果 Z 定义为“障碍越高越坏”，而你希望“通路激活越高越好”，通常希望：
      #   HO 若使通路下降(delta_M<0) 且 beta_M<0 (通路越高，障碍越低) => 解释为保护性链条。
      # 但不同模块可能相反；因此这里只做“强度+稳健性”排序，不强行定义保护性标签。
      score_strength = abs(delta_M) * abs(beta_M),
      score_robust   = ifelse(is.na(frac_ok), 0, frac_ok),
      score_final    = score_strength * (0.5 + 0.5*score_robust)
    ) %>%
    arrange(desc(score_final))

  readr::write_tsv(rank_combined, file.path(outdir, paste0("rank_combined_M_for_Z", suffix, ".tsv")))
  cat("  [OK] 写出 rank_combined_M_for_Z", suffix, ".tsv\n", sep = "")
}

# ---- main execution ----
if (Z_col == "__ALL__") {
  cat("[INFO] 识别到多机制轴模式，开始对所有候选 Z 逐一计算，共 ", length(all_Z_cols), " 个。\n", sep = "")
  for (i in seq_along(all_Z_cols)) {
    cat("\n============================================\n")
    cat("[INFO] 处理第 ", i, "/", length(all_Z_cols), " 个 Z: ", all_Z_cols[i], "\n", sep = "")
    run_for_one_Z(all_Z_cols[i])
  }
} else {
  run_for_one_Z(Z_col)
}

cat("\n============================================================\n")
cat("[DONE] 06_rna_path_rank_XMZ.R 完成。\n")
cat("============================================================\n")