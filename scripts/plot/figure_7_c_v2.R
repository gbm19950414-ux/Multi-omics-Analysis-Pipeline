#!/usr/bin/env Rscript

## ============================================================
## 06a_plot_rank_XMZ_panels_ABC.R
## Figure6 (RNA-only) plotting script: A/B/C panels
##
## Inputs:
##   1) sample_scores_rna_only.tsv (X/M/Z wide table)
##   2) rank_X_to_M.tsv            (from 06_rna_path_rank_XMZ.R)
##   3) rank_combined_M_for_Z.tsv  (from 06_rna_path_rank_XMZ.R)
##   4) loo_matrix_M_to_Z.tsv      (from 06_rna_path_rank_XMZ.R)
##
## Outputs (outdir):
##   - fig6A_deltaM_bar.png
##   - fig6B_scatter_topM.png
##   - fig6C_LOO_heatmap.png
##   - fig6_ABC_combined.png  (if patchwork available; otherwise skipped)
##
## Usage:
##   Rscript scripts/multi/06a_plot_rank_XMZ_panels_ABC.R \
##     results/multiomics/tables/sample_scores/sample_scores_rna_only.tsv \
##     results/multiomics/tables/mediation \
##     "Membrane_context_mito_mass_surface_area" \
##     12 \
##     3
##     # (top_n_A=12, top_k_scatter=3)
## ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)

# Defaults: allow running this script with no args.
# You can override any of them by passing the usual CLI arguments.
default_sample_scores_path <- "results/multiomics/tables/sample_scores/sample_scores_rna_only.tsv"
default_outdir <- "results/multiomics/tables/mediation"
## Prefer a Z that actually exists in your current sample_scores_rna_only.tsv
## (Change this default if you want a different Z axis.)
# Setting to "__ALL__" will plot all Z axes found in rank_combined_M_for_Z__*.tsv files in outdir.
default_Z_col <- "__ALL__"
default_top_n_A <- 12L
default_top_k_scatter <- 3L

if (length(args) >= 3) {
  sample_scores_path <- args[1]
  outdir <- args[2]
  Z_col <- args[3]
  top_n_A <- if (length(args) >= 4) as.integer(args[4]) else default_top_n_A
  top_k_scatter <- if (length(args) >= 5) as.integer(args[5]) else default_top_k_scatter
} else {
  sample_scores_path <- default_sample_scores_path
  outdir <- default_outdir
  Z_col <- default_Z_col
  top_n_A <- default_top_n_A
  top_k_scatter <- default_top_k_scatter
  cat("[INFO] No CLI args provided; using defaults:\n")
  cat("  sample_scores : ", sample_scores_path, "\n", sep = "")
  cat("  outdir        : ", outdir, "\n", sep = "")
  cat("  Z_col         : ", Z_col, "\n", sep = "")
  cat("  top_n_A       : ", top_n_A, "\n", sep = "")
  cat("  top_k_scatter : ", top_k_scatter, "\n", sep = "")
  cat("\n")
}

# If someone provides 1-2 args by mistake, show usage.
if (length(args) > 0 && length(args) < 3) {
  stop(
    "Usage: Rscript ", basename(commandArgs(trailingOnly = FALSE)[1]),
    " <sample_scores_rna_only.tsv> <mediation_outdir> <Z_col> [top_n_A] [top_k_scatter]\n",
    "(Or run with no args to use defaults.)\n"
  )
}

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ---- helper: safe filename suffix ----
sanitize_tag <- function(x) {
  x <- gsub("\\s+", "_", x)
  x <- gsub("[^A-Za-z0-9_\\-]", "_", x)
  x
}

# ---- helper: resolve per-Z input paths (prefer __Z suffix if exists) ----
resolve_paths_for_Z <- function(outdir, Zname) {
  tag <- sanitize_tag(Zname)
  paths <- list(
    rank_X = file.path(outdir, paste0("rank_X_to_M__", Zname, ".tsv")),
    rank_MZ = file.path(outdir, paste0("rank_M_to_Z__", Zname, ".tsv")),
    rank_comb = file.path(outdir, paste0("rank_combined_M_for_Z__", Zname, ".tsv")),
    loo = file.path(outdir, paste0("loo_matrix_M_to_Z__", Zname, ".tsv"))
  )

  # fallback to unsuffixed names (older mode)
  if (!file.exists(paths$rank_X)) paths$rank_X <- file.path(outdir, "rank_X_to_M.tsv")
  if (!file.exists(paths$rank_MZ)) paths$rank_MZ <- file.path(outdir, "rank_M_to_Z.tsv")
  if (!file.exists(paths$rank_comb)) paths$rank_comb <- file.path(outdir, "rank_combined_M_for_Z.tsv")
  if (!file.exists(paths$loo)) paths$loo <- file.path(outdir, "loo_matrix_M_to_Z.tsv")

  paths$tag <- tag
  paths
}

# ---- helper: discover all Z names from per-Z files ----
discover_Zs <- function(outdir) {
  f <- list.files(outdir, pattern = "^rank_combined_M_for_Z__.*\\.tsv$", full.names = FALSE)
  if (length(f) == 0) return(character(0))
  z <- sub("^rank_combined_M_for_Z__", "", f)
  z <- sub("\\.tsv$", "", z)
  z
}

# Check that sample_scores_path exists
if (!file.exists(sample_scores_path)) {
  stop("[ERROR] Missing input: ", sample_scores_path)
}

cat("============================================================\n")
cat("[INFO] figure_7_c_v2.R (ABC panels)\n")
cat("  sample_scores : ", sample_scores_path, "\n", sep="")
cat("  outdir        : ", outdir, "\n", sep="")
cat("  Z_col         : ", Z_col, "\n", sep="")
cat("  top_n_A       : ", top_n_A, "\n", sep="")
cat("  top_k_scatter : ", top_k_scatter, "\n", sep="")
cat("============================================================\n\n")

df <- readr::read_tsv(sample_scores_path, show_col_types = FALSE)

# normalize group labels
df <- df %>%
  mutate(group = case_when(
    str_to_upper(group) %in% c("WT") ~ "WT",
    str_to_upper(group) %in% c("HO","KO") ~ "HO",
    TRUE ~ group
  ))

# ------------------------- Plotting function ---------------------
plot_for_one_Z <- function(Z_col, paths) {
  # Read per-Z files
  rank_X <- readr::read_tsv(paths$rank_X, show_col_types = FALSE)
  rank_comb <- readr::read_tsv(paths$rank_comb, show_col_types = FALSE)
  loo <- readr::read_tsv(paths$loo, show_col_types = FALSE)

  if (!Z_col %in% colnames(df)) {
    stop("[ERROR] Z_col not found in sample_scores: ", Z_col, "\nAvailable: ", paste(colnames(df), collapse=", "))
  }
  req <- c("sample_id", "batch", "group")
  miss <- setdiff(req, colnames(df))
  if (length(miss) > 0) stop("[ERROR] sample_scores missing columns: ", paste(miss, collapse=", "))

  # ---------- Panel A: ΔM bar ----------
  # rank_X already contains delta_M_median etc
  rank_A <- rank_X %>%
    mutate(abs_delta = abs(delta_M_median)) %>%
    arrange(desc(abs_delta)) %>%
    slice_head(n = top_n_A) %>%
    mutate(M = factor(M, levels = rev(M)))

  pA <- ggplot(rank_A, aes(x = M, y = delta_M_median)) +
    geom_col() +
    coord_flip() +
    labs(
      title = "Figure6A  X→M (HO−WT) effect size",
      x = NULL,
      y = "ΔM median across batches (HO−WT)"
    ) +
    theme_bw(base_size = 11)

  ggsave(file.path(outdir, paste0("fig6A_deltaM_bar__", paths$tag, ".png")), pA, width = 7.5, height = 4.5, dpi = 300)
  cat("[OK] fig6A_deltaM_bar__", paths$tag, ".png\n", sep = "")

  # ---------- choose top mediators for B/C from combined rank ----------
  # rank_comb may include M names that are not present in the current sample_scores table
  # (e.g., if earlier ranking treated other numeric columns as mediators). We only keep
  # columns that exist in df, and prefer EphB1_ mediators.
  rank_comb2 <- rank_comb %>% arrange(desc(score_final))

  # prefer EphB1_ mediators if available
  rank_pref <- rank_comb2 %>% filter(grepl("^EphB1_", M))
  if (nrow(rank_pref) == 0) rank_pref <- rank_comb2

  candM <- rank_pref$M
  candM <- candM[candM %in% colnames(df)]

  if (length(candM) == 0) {
    stop("[ERROR] No mediators from rank_combined exist in sample_scores columns.\n",
         "Check that rank_combined_M_for_Z.tsv corresponds to the same sample_scores table and that mediator columns are present.")
  }

  if (length(candM) < top_k_scatter) {
    cat("[WARN] Only ", length(candM), " mediator columns found in sample_scores; using all of them for scatter/LOO.\n", sep="")
  }

  topM <- head(candM, top_k_scatter)

  # ---------- Panel B: scatter Z vs M (topK), colored by group, shaped by batch ----------
  df_longB <- df %>%
    select(sample_id, batch, group, Z = all_of(Z_col), all_of(topM)) %>%
    pivot_longer(all_of(topM), names_to = "M", values_to = "M_score") %>%
    mutate(M = factor(M, levels = topM))

  pB <- ggplot(df_longB, aes(x = M_score, y = Z, shape = batch)) +
    geom_point(aes(color = group), size = 2.6, alpha = 0.9) +
    scale_shape_manual(values = c(16, 15, 17, 18)[seq_along(unique(df_longB$batch))]) +
    facet_wrap(~ M, scales = "free_x") +
    labs(
      title = paste0("Figure6B  Z vs M (Z = ", Z_col, ")"),
      x = "Mediator score (M)",
      y = "Axis score (Z)"
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "right")

  ggsave(file.path(outdir, paste0("fig6B_scatter_topM__", paths$tag, ".png")), pB, width = 9.0, height = 4.8, dpi = 300)
  cat("[OK] fig6B_scatter_topM__", paths$tag, ".png\n", sep = "")

  # ---------- Panel C: LOO heatmap for topK mediators ----------
  loo_C <- loo %>%
    filter(M %in% topM) %>%
    mutate(
      sample_left_out = factor(sample_left_out, levels = unique(sample_left_out)),
      M = factor(M, levels = topM)
    )

  # add full beta per M (for annotation)
  full_beta_tbl <- loo_C %>%
    group_by(M) %>%
    summarise(beta_M_full = beta_M_full[1], .groups = "drop")

  pC <- ggplot(loo_C, aes(x = sample_left_out, y = M, fill = beta_M_loo)) +
    geom_tile(color = "white", linewidth = 0.2) +
    labs(
      title = "Figure6C  LOO robustness (pooled): beta_M_loo",
      x = "Left-out sample",
      y = "Mediator (top hits)"
    ) +
    theme_bw(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      panel.grid = element_blank()
    )

  ggsave(file.path(outdir, paste0("fig6C_LOO_heatmap__", paths$tag, ".png")), pC, width = 9.0, height = 3.8, dpi = 300)
  cat("[OK] fig6C_LOO_heatmap__", paths$tag, ".png\n", sep = "")

  # ---------- optional combined (patchwork) ----------
  if (requireNamespace("patchwork", quietly = TRUE)) {
    library(patchwork)
    pABC <- pA / pB / pC + patchwork::plot_layout(heights = c(1.0, 1.15, 0.9))
    ggsave(file.path(outdir, paste0("fig6_ABC_combined__", paths$tag, ".png")), pABC, width = 10.0, height = 12.0, dpi = 300)
    cat("[OK] fig6_ABC_combined__", paths$tag, ".png\n", sep = "")
  } else {
    cat("[SKIP] patchwork not installed; combined figure not generated.\n")
  }
}

# ------------------------- Main control flow ---------------------
if (Z_col == "__ALL__") {
  zs <- discover_Zs(outdir)
  if (length(zs) == 0) {
    stop("[ERROR] No per-Z files found in outdir; nothing to plot for __ALL__ mode.")
  }
  n_done <- 0L
  for (Zname in zs) {
    if (!Zname %in% colnames(df)) {
      cat("[WARN] Z axis '", Zname, "' not found in sample_scores table; skipping.\n", sep = "")
      next
    }
    paths <- resolve_paths_for_Z(outdir, Zname)
    needed <- c(paths$rank_X, paths$rank_comb, paths$loo)
    missing <- needed[!file.exists(needed)]
    if (length(missing) > 0) {
      cat("[WARN] Missing per-Z input file(s) for Z = ", Zname, ": ", paste(missing, collapse = ", "), "; skipping.\n", sep = "")
      next
    }
    plot_for_one_Z(Zname, paths)
    n_done <- n_done + 1L
  }
  cat("============================================================\n")
  cat("[DONE] Plotted ", n_done, " Z axes (per-Z mode).\n", sep = "")
  cat("============================================================\n")
} else {
  # If default/CLI Z_col doesn't exist, try common alternatives.
  if (!Z_col %in% colnames(df)) {
    z_candidates <- c(
      Z_col,
      "Membrane context",
      "Membrane_context",
      "Membrane context_score",
      "Membrane_context_score"
    )
    z_pick <- z_candidates[z_candidates %in% colnames(df)][1]
    if (!is.na(z_pick) && nchar(z_pick) > 0) {
      cat("[WARN] Z_col not found; using fallback Z_col = ", z_pick, "\n", sep = "")
      Z_col <- z_pick
    }
  }
  paths <- resolve_paths_for_Z(outdir, Z_col)
  needed <- c(paths$rank_X, paths$rank_comb, paths$loo)
  missing <- needed[!file.exists(needed)]
  if (length(missing) > 0) {
    stop("[ERROR] Missing input file(s) for Z = ", Z_col, ": ", paste(missing, collapse = ", "))
  }
  plot_for_one_Z(Z_col, paths)
  cat("============================================================\n")
  cat("[DONE] 06a plotting done.\n")
  cat("============================================================\n")
}