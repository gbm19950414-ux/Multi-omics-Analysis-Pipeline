#!/usr/bin/env Rscript

# 04.3_qpcr_value_ranking_viz.R
# Visualize qPCR-value ranking table produced by 04.2
#
# Usage:
#   Rscript 04.3_qpcr_value_ranking_viz.R <ranking.tsv> [out_prefix] [top_n]
#
# Examples:
#   Rscript 04.3_qpcr_value_ranking_viz.R EphB1_PI3K_AKT_mTOR_qpcr_value_ranking.tsv
#   Rscript 04.3_qpcr_value_ranking_viz.R ranking.tsv out/PI3K_rank 12
#
# Outputs (PNG + PDF):
#   <out_prefix>__score_vs_rank.(png|pdf)
#   <out_prefix>__topN_bar.(png|pdf)
#   <out_prefix>__topN_heatmap_signedLFC.(png|pdf)

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
  cat("Usage: Rscript 04.3_qpcr_value_ranking_viz.R <ranking.tsv> [out_prefix] [top_n]\n")
  quit(status = 1)
}

in_tsv <- args[1]

# Default out_prefix: same dir, file stem + "_vis"
if (length(args) >= 2) {
  out_prefix <- args[2]
} else {
  base_dir <- dirname(normalizePath(in_tsv))
  stem <- tools::file_path_sans_ext(basename(in_tsv))
  out_prefix <- file.path(base_dir, paste0(stem, "_vis"))
}

top_n <- if (length(args) >= 3) as.integer(args[3]) else 10
if (is.na(top_n) || top_n < 1) top_n <- 10

# --------- Load ----------
df <- readr::read_tsv(in_tsv, show_col_types = FALSE)

required_cols <- c("gene_symbol", "score", "rank_in_set", "Tier", "recommended_topN")
missing_cols <- setdiff(required_cols, colnames(df))
if (length(missing_cols) > 0) {
  stop("Input TSV missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
}

# Discover signed LFC batch columns (from 04.2)
signed_cols <- grep("^signed_LFC_batch\\d+$", colnames(df), value = TRUE)
lfc_cols <- grep("^LFC_batch\\d+$", colnames(df), value = TRUE)

if (length(signed_cols) == 0 && length(lfc_cols) > 0) {
  warning("No signed_LFC_batch* columns found. Heatmap will use LFC_batch* (unsigned).")
}

# Create output dir if needed
out_dir <- dirname(out_prefix)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Helper: save plot as PNG+PDF
save_both <- function(p, prefix, w = 10, h = 6, dpi = 300) {
  ggsave(paste0(prefix, ".png"), p, width = w, height = h, dpi = dpi, bg = "white")
  ggsave(paste0(prefix, ".pdf"), p, width = w, height = h, bg = "white")
}

# --------- Plot 1: score vs rank (overview) ----------
# Ensure stable ordering
df_plot <- df %>%
  mutate(
    Tier = factor(Tier, levels = c("Tier1", "Tier2", "Tier3", "NotEligible")),
    recommended_topN = as.logical(recommended_topN)
  ) %>%
  arrange(rank_in_set)

p1 <- ggplot(df_plot, aes(x = rank_in_set, y = score, shape = recommended_topN, color = Tier)) +
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

save_both(p1, paste0(out_prefix, "__score_vs_rank"), w = 10, h = 6)

# --------- Plot 2: topN bar ----------
top_df <- df_plot %>%
  filter(eligible == TRUE) %>%
  arrange(desc(score)) %>%
  slice_head(n = top_n) %>%
  mutate(gene_symbol = factor(gene_symbol, levels = rev(gene_symbol)))

p2 <- ggplot(top_df, aes(x = gene_symbol, y = score, fill = Tier)) +
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

save_both(p2, paste0(out_prefix, "__topN_bar"), w = 10, h = 7)

# --------- Plot 3: heatmap of (signed) LFC across batches for topN ----------
heat_cols <- if (length(signed_cols) > 0) signed_cols else lfc_cols
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

  save_both(p3, paste0(out_prefix, "__topN_heatmap_signedLFC"), w = 8, h = max(4, 0.35 * top_n + 2))
} else {
  warning("No LFC batch columns found; skipping heatmap plot.")
}

cat("[OK] Wrote outputs with prefix: ", out_prefix, "\n", sep = "")
