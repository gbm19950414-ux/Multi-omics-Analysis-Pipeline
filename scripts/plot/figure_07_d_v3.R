#!/usr/bin/env Rscript

## ============================================================
## figure_07_d_v3.R  (Figure6e panel)
##
## Input:
##   concordance_Z_lipid.tsv
##     columns required: batch, Z_effect, lipid_effect, direction_match
##
## Output (outdir):
##   - fig6e_concordance_slope.png
##   - fig6e_concordance_bar.png
##
## Usage:
##   Rscript scripts/multi/figure_07_d_v3.R <concordance_tsv> <outdir>
##
## Or run with no args (defaults):
##   concordance_tsv = results/multiomics/tables/mediation/concordance_Z_lipid.tsv
##   outdir          = results/multiomics/tables/mediation
## ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)

default_in <- "results/multiomics/tables/mediation/concordance_Z_lipid.tsv"
default_outdir <- "results/multiomics/tables/mediation"

if (length(args) == 0) {
  in_path <- default_in
  outdir <- default_outdir
  cat("[INFO] No CLI args provided; using defaults:\n")
  cat("  input : ", in_path, "\n", sep = "")
  cat("  outdir: ", outdir, "\n\n", sep = "")
} else if (length(args) == 1) {
  in_path <- args[1]
  outdir <- default_outdir
  cat("[INFO] One CLI arg provided; using:\n")
  cat("  input : ", in_path, "\n", sep = "")
  cat("  outdir: ", outdir, "\n\n", sep = "")
} else if (length(args) >= 2) {
  in_path <- args[1]
  outdir <- args[2]
} else {
  stop("Usage: Rscript figure_07_d_v3.R <concordance_tsv> <outdir>\n(Or run with no args for defaults.)\n")
}

if (!file.exists(in_path)) stop("[ERROR] Input not found: ", in_path)
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

df <- readr::read_tsv(in_path, show_col_types = FALSE)

req <- c("batch", "Z_effect", "lipid_effect", "direction_match")
miss <- setdiff(req, colnames(df))
if (length(miss) > 0) stop("[ERROR] Missing columns in concordance table: ", paste(miss, collapse=", "))

# Normalize batch ordering
df <- df %>%
  mutate(
    batch = as.character(batch),
    batch = factor(batch, levels = unique(batch)),
    direction_match = case_when(
      is.na(direction_match) ~ NA,
      direction_match %in% c(TRUE, "TRUE", "T") ~ TRUE,
      direction_match %in% c(FALSE, "FALSE", "F") ~ FALSE,
      TRUE ~ NA
    ),
    match_label = case_when(
      is.na(direction_match) ~ "NA",
      direction_match ~ "match",
      !direction_match ~ "mismatch"
    )
  )

# Long format for plotting
df_long <- df %>%
  select(batch, Z_effect, lipid_effect, direction_match, match_label) %>%
  pivot_longer(c(Z_effect, lipid_effect), names_to = "metric", values_to = "effect") %>%
  mutate(
    metric = recode(metric,
      Z_effect = "RNA Z (HO-WT)",
      lipid_effect = "Lipid (HO-WT)"
    ),
    metric = factor(metric, levels = c("RNA Z (HO-WT)", "Lipid (HO-WT)"))
  )

# ---------- Figure6e option 1: slopegraph (paired effects per batch) ----------
p_slope <- ggplot(df_long, aes(x = metric, y = effect, group = batch)) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  geom_line(aes(linetype = match_label), linewidth = 0.7, alpha = 0.9) +
  geom_point(aes(shape = batch), size = 2.6, alpha = 0.95) +
  facet_wrap(~ batch, nrow = 1, scales = "free_y") +
  labs(
    title = "Figure6e  Concordance of per-batch effects (HO-WT)",
    x = NULL,
    y = "Effect size (HO-WT)",
    linetype = "Direction"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

ggsave(file.path(outdir, "fig6e_concordance_slope.png"),
       p_slope, width = 11.5, height = 3.6, dpi = 300)

# ---------- Figure6e option 2: grouped bars (two bars per batch) ----------
# This is often easier to read if effect scales are comparable across batches.
p_bar <- ggplot(df_long, aes(x = batch, y = effect, fill = metric)) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  geom_col(position = position_dodge(width = 0.75), width = 0.68) +
  # annotate match/mismatch above batch
  geom_text(
    data = df %>% mutate(y_annot = pmax(abs(Z_effect), abs(lipid_effect), na.rm = TRUE) * 1.1),
    aes(x = batch, y = y_annot, label = match_label),
    inherit.aes = FALSE,
    size = 3.2,
    vjust = 0
  ) +
  labs(
    title = "Figure6e  Concordance summary by batch",
    x = "Batch",
    y = "Effect size (HO-WT)",
    fill = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

ggsave(file.path(outdir, "fig6e_concordance_bar.png"),
       p_bar, width = 7.2, height = 4.2, dpi = 300)

cat("[OK] wrote:\n")
cat("  ", file.path(outdir, "fig6e_concordance_slope.png"), "\n", sep = "")
cat("  ", file.path(outdir, "fig6e_concordance_bar.png"), "\n", sep = "")
cat("[DONE] figure_07_d_v3.R\n")