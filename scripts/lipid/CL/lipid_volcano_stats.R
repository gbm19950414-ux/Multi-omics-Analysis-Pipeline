#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)

# Inputs
in_csv  <- if (length(args) >= 1) args[[1]] else "data/raw/lipid/batch1_clean.csv"
out_csv <- if (length(args) >= 2) args[[2]] else "data/interim/lipid/batch1_volcano.csv"

# Options (can be overridden via args if you want later)
test_method <- if (length(args) >= 3) args[[3]] else "wilcox"  # "wilcox" or "t"
pseudo_count <- if (length(args) >= 4) as.numeric(args[[4]]) else 1.0  # for log2FC stability

# Ensure output dir exists
dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)

# ---- read ----
dat0 <- readr::read_csv(in_csv, show_col_types = FALSE)

# Drop common first row like "Group" (prevents sample cols being character)
if ("lipidName" %in% names(dat0)) {
  dat0 <- dat0 %>% filter(!is.na(lipidName), lipidName != "Group")
} else {
  stop("Column `lipidName` not found in input CSV.")
}

# sample columns: WT_1, HO_5, QC_1 ...
sample_cols <- names(dat0)[str_detect(names(dat0), "^(WT|HO|QC)_\\d+$")]
if (length(sample_cols) == 0) stop("No sample columns matched ^(WT|HO|QC)_\\d+$")

# coerce numeric (avoid character columns)
dat0 <- dat0 %>% mutate(across(all_of(sample_cols), ~ suppressWarnings(as.numeric(.x))))

# keep only WT/HO sample columns (exclude QC)
wt_cols <- sample_cols[str_detect(sample_cols, "^WT_")]
ho_cols <- sample_cols[str_detect(sample_cols, "^HO_")]
if (length(wt_cols) < 2 || length(ho_cols) < 2) {
  stop("Need >=2 WT and >=2 HO replicates. Found WT=", length(wt_cols), " HO=", length(ho_cols))
}

# optional columns
cls_col <- if ("class" %in% names(dat0)) "class" else NULL

# ---- volcano stats per lipid ----
# IMPORTANT: keep WT/HO columns until after row-wise p-value computation
vol <- dat0 %>%
  mutate(
    mean_WT = rowMeans(as.matrix(across(all_of(wt_cols))), na.rm = TRUE),
    mean_HO = rowMeans(as.matrix(across(all_of(ho_cols))), na.rm = TRUE),
    log2FC = log2((mean_HO + pseudo_count) / (mean_WT + pseudo_count))
  ) %>%
  rowwise() %>%
  mutate(
    p_value = {
      x <- c_across(all_of(wt_cols))
      y <- c_across(all_of(ho_cols))
      x <- x[!is.na(x)]; y <- y[!is.na(y)]
      if (length(x) < 2 || length(y) < 2) NA_real_
      else if (test_method == "t") t.test(x, y)$p.value else wilcox.test(x, y)$p.value
    }
  ) %>%
  ungroup() %>%
  transmute(
    lipidName = .data$lipidName,
    class = if (!is.null(cls_col)) .data[[cls_col]] else NA_character_,
    mean_WT,
    mean_HO,
    log2FC,
    p_value,
    padj = p.adjust(p_value, method = "BH"),
    neglog10p = -log10(p_value),
    neglog10padj = -log10(padj)
  )

# clean up infinities
vol <- vol %>%
  mutate(
    neglog10p = ifelse(is.finite(neglog10p), neglog10p, NA_real_),
    neglog10padj = ifelse(is.finite(neglog10padj), neglog10padj, NA_real_)
  )

# write
readr::write_csv(vol, out_csv)
message("[OK] wrote: ", normalizePath(out_csv))