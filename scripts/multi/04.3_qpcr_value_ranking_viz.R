#!/usr/bin/env Rscript

# 04.3_qpcr_value_ranking_viz.R
# Visualize qPCR-value ranking table(s) produced by 04.2
#
# Usage:
#   Rscript 04.3_qpcr_value_ranking_viz.R <ranking.tsv|ranking_dir> [out_dir] [top_n]
#
# Examples:
#   Rscript 04.3_qpcr_value_ranking_viz.R results/multiomics/tables/qpcr_value_ranking
#   Rscript 04.3_qpcr_value_ranking_viz.R EphB1_PI3K_AKT_mTOR_qpcr_value_ranking.tsv
#   Rscript 04.3_qpcr_value_ranking_viz.R results/multiomics/tables/qpcr_value_ranking out/qpcr_viz 12
#
# Output:
#   One PDF per input TSV (multi-page):
#     <out_dir>/<tsv_stem>_vis.pdf

suppressPackageStartupMessages({
  pkgs <- c("ggplot2", "dplyr", "readr", "tidyr")
  missing <- pkgs[!vapply(pkgs, requireNamespace, FUN.VALUE = logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop(
      "Missing R packages: ", paste(missing, collapse = ", "),
      "\nPlease install them first, e.g.: install.packages(c(",
      paste(sprintf('"%s"', missing), collapse = ", "),
      "))",
      call. = FALSE
    )
  }
})

library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  cat("Usage: Rscript 04.3_qpcr_value_ranking_viz.R <ranking.tsv|ranking_dir> [out_dir] [top_n]\n")
  quit(status = 1)
}

in_path <- args[1]

out_dir <- if (length(args) >= 2) args[2] else NULL
top_n <- if (length(args) >= 3) as.integer(args[3]) else 10
if (is.na(top_n) || top_n < 1) top_n <- 10

if (dir.exists(in_path)) {
  tsv_files <- list.files(in_path, pattern = "_qpcr_value_ranking(.*)\\.tsv$", full.names = TRUE)
  # Skip summary tables
  # (e.g., qpcr_value_ranking_topN_summary.tsv)
  tsv_files <- tsv_files[!grepl("qpcr_value_ranking_topN_summary\\.tsv$", basename(tsv_files))]
  if (length(tsv_files) == 0) {
    stop("No ranking TSV files found in directory: ", in_path, call. = FALSE)
  }
} else if (file.exists(in_path)) {
  tsv_files <- c(in_path)
} else {
  stop("Input path does not exist: ", in_path, call. = FALSE)
}

save_one_pdf <- function(pdf_path, plots, w = 10, h = 7) {
  grDevices::pdf(pdf_path, width = w, height = h)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (p in plots) {
    print(p)
  }
}

for (in_tsv in tsv_files) {
  df <- readr::read_tsv(in_tsv, show_col_types = FALSE)

  required_cols <- c("gene_symbol", "score", "rank_in_set", "Tier", "recommended_topN", "eligible")
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop("Input TSV missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  signed_cols <- grep("^signed_LFC_batch\\d+$", colnames(df), value = TRUE)
  lfc_cols <- grep("^LFC_batch\\d+$", colnames(df), value = TRUE)

  if (length(signed_cols) == 0 && length(lfc_cols) > 0) {
    warning("No signed_LFC_batch* columns found. Heatmap will use LFC_batch* (unsigned).")
  }

  out_dir_use <- out_dir
  if (is.null(out_dir_use)) {
    out_dir_use <- dirname(normalizePath(in_tsv))
  }
  if (!dir.exists(out_dir_use)) dir.create(out_dir_use, recursive = TRUE, showWarnings = FALSE)

  stem <- tools::file_path_sans_ext(basename(in_tsv))
  pdf_path <- file.path(out_dir_use, paste0(stem, "_vis.pdf"))

  # Prepare df_plot with score_plot column
  df_plot <- df %>%
    mutate(
      Tier = factor(Tier, levels = c("Tier1", "Tier2", "Tier3", "NotEligible")),
      recommended_topN = as.logical(recommended_topN),
      score_plot = if ("score_used_for_ranking" %in% colnames(df)) score_used_for_ranking else score
    ) %>%
    arrange(rank_in_set)

  p1 <- ggplot(df_plot, aes(x = rank_in_set, y = score_plot, shape = recommended_topN, color = Tier)) +
    geom_point(size = 2.8, alpha = 0.9) +
    # Label only Tier1 and recommended genes to avoid clutter
    geom_text(
      data = df_plot %>% filter(Tier == "Tier1" | recommended_topN),
      aes(label = gene_symbol),
      vjust = -0.8,
      size = 3,
      show.legend = FALSE
    ) +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17)) +
    labs(
      title = "qPCR validation value ranking",
      subtitle = paste0("Input: ", basename(in_tsv)),
      x = "Rank in set (1 = highest)",
      y = "Composite score (higher = better)",
      shape = "Recommended topN",
      color = "Tier"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "right"
    )

  top_df <- df_plot %>%
    filter(eligible == TRUE) %>%
    arrange(desc(score_plot)) %>%
    slice_head(n = top_n) %>%
    mutate(gene_symbol = factor(gene_symbol, levels = rev(gene_symbol)))

  p2 <- ggplot(top_df, aes(x = gene_symbol, y = score_plot, fill = Tier)) +
    geom_col(width = 0.75) +
    coord_flip() +
    labs(
      title = paste0("Top ", top_n, " genes by qPCR validation value"),
      x = NULL,
      y = "Composite score"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "right"
    )

  heat_cols <- if (length(signed_cols) > 0) signed_cols else lfc_cols
  plots <- list(p1, p2)
  if (length(heat_cols) > 0) {
    heat_long <- top_df %>%
      select(gene_symbol, all_of(heat_cols)) %>%
      pivot_longer(cols = all_of(heat_cols), names_to = "batch", values_to = "LFC") %>%
      mutate(
        batch = gsub("^signed_", "", batch),
        batch = factor(batch, levels = unique(batch))
      )

    p3 <- ggplot(heat_long, aes(x = batch, y = gene_symbol, fill = LFC)) +
      geom_tile(color = "white", linewidth = 0.25) +
      labs(
        title = if (length(signed_cols) > 0) "TopN signed LFC consistency across batches" else "TopN LFC across batches",
        subtitle = "Helps see direction consistency and outlier batches",
        x = NULL,
        y = NULL
      ) +
      theme_minimal(base_size = 12) +
      theme(
        plot.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 0, hjust = 0.5)
      )

    plots[[length(plots) + 1]] <- p3
  }

  save_one_pdf(pdf_path, plots, w = 10, h = 7)
  cat("[OK] Wrote: ", pdf_path, "\n", sep = "")
}

cat("[DONE] Visualized ", length(tsv_files), " file(s).\n", sep = "")
