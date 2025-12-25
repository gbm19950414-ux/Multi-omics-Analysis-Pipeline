#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)

in_csv  <- if (length(args) >= 1) args[[1]] else "data/raw/lipid/batch1_clean.csv"
out_dir <- if (length(args) >= 2) args[[2]] else "data/interim/lipid"
method  <- if (length(args) >= 3) args[[3]] else "wilcox"  # "wilcox" or "t"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------- read ----------
dat0 <- readr::read_csv(in_csv, show_col_types = FALSE)

# drop the "Group" row (common in your file)
if ("lipidName" %in% names(dat0)) {
  dat0 <- dat0 %>% filter(!is.na(lipidName), lipidName != "Group")
}

# sample columns: WT_1 HO_5 QC_1 etc
sample_cols <- names(dat0)[str_detect(names(dat0), "^(WT|HO|QC)_\\d+$")]
if (length(sample_cols) == 0) stop("No sample columns matched ^(WT|HO|QC)_\\d+$")

# coerce numeric
dat0 <- dat0 %>% mutate(across(all_of(sample_cols), ~ suppressWarnings(as.numeric(.x))))

# ---------- filter CL ----------
if (!("class" %in% names(dat0))) stop("Column `class` not found.")
cl <- dat0 %>% filter(class == "CL")

# ---------- long format (raw points) ----------
cl_long <- cl %>%
  select(lipidName, all_of(sample_cols)) %>%
  pivot_longer(cols = all_of(sample_cols), names_to = "sample", values_to = "value") %>%
  mutate(
    genetype = case_when(
      str_detect(sample, "^WT_") ~ "WT",
      str_detect(sample, "^HO_") ~ "HO",
      str_detect(sample, "^QC_") ~ "QC",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(genetype %in% c("WT", "HO")) %>%   # exclude QC
  filter(!is.na(value))

# ---------- stats per lipidName ----------
cl_stats <- cl_long %>%
  group_by(lipidName) %>%
  summarise(
    n_WT = sum(genetype == "WT"),
    n_HO = sum(genetype == "HO"),
    mean_WT = mean(value[genetype == "WT"], na.rm = TRUE),
    sd_WT   = sd(value[genetype == "WT"], na.rm = TRUE),
    mean_HO = mean(value[genetype == "HO"], na.rm = TRUE),
    sd_HO   = sd(value[genetype == "HO"], na.rm = TRUE),
    p_value = {
      x <- value[genetype == "WT"]
      y <- value[genetype == "HO"]
      if (length(na.omit(x)) < 2 || length(na.omit(y)) < 2) NA_real_
      else if (method == "t") t.test(x, y)$p.value else wilcox.test(x, y)$p.value
    },
    .groups = "drop"
  ) %>%
  arrange(p_value)

# ---------- write ----------
out_long  <- file.path(out_dir, "batch1_clean_cl_long.csv")
out_stats <- file.path(out_dir, "batch1_clean_cl_stats.csv")
readr::write_csv(cl_long,  out_long)
readr::write_csv(cl_stats, out_stats)

message("[OK] wrote: ", out_long)
message("[OK] wrote: ", out_stats)