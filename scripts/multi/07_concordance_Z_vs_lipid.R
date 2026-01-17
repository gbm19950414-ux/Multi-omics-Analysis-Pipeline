#!/usr/bin/env Rscript

## ============================================================
## 06e_concordance_Z_vs_lipid.R
##
## Figure6e: RNA Z axis vs lipid phenotype concordance (per-batch)
##
## Inputs:
##   1) RNA-only sample-level table: sample_scores_rna_only.tsv
##      required columns: batch, group (WT/HO), and Z_col
##   2) Lipid per-sample table: lipid_mechanistic_indicators_per_sample.tsv
##      required columns: batch, group (WT/HO), and lipid_col
##
## Output:
##   concordance_Z_lipid.tsv
##     columns: batch, Z_effect, lipid_effect, direction_match (T/F/NA)
##
## Effect definition per batch:
##   effect = mean(HO) - mean(WT)
##
## Usage:
##   Rscript scripts/multi/06e_concordance_Z_vs_lipid.R \
##     <rna_scores_tsv> <lipid_scores_tsv> <outdir> <Z_col> <lipid_col>
##
## Or run with no args to use defaults.
## ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
})

args <- commandArgs(trailingOnly = TRUE)

# defaults (adjust if your paths differ)
rna_scores_tsv_default   <- "results/multiomics/tables/sample_scores/sample_scores_rna_only.tsv"
lipid_scores_tsv_default <- "results/lipid/tables/lipid_mechanistic_indicators_per_sample.tsv"
outdir_default           <- "results/multiomics/tables/mediation"

# fallbacks if you don't pass Z/lipid col
Z_col_default_candidates <- c(
  "Membrane_context_score",
  "Membrane context_score",
  "Membrane_context_mito_mass_surface_area",
  "Membrane_context_OXPHOS_supercomplex_organization",
  # extra candidates for new axis naming
  "Membrane context",
  "Oxidation",
  "Remodeling",
  "Supply",
  "Synthesis",
  "Transport"
)

lipid_col_default_candidates <- c(
  "Membrane context_score",
  "Membrane_context_score",
  "sum_CL_score",
  "CL_score",
  "sum_CL",
  "CL"
)

if (length(args) >= 5) {
  rna_scores_tsv   <- args[1]
  lipid_scores_tsv <- args[2]
  outdir           <- args[3]
  Z_col            <- args[4]
  lipid_col        <- args[5]
} else if (length(args) == 0) {
  rna_scores_tsv   <- rna_scores_tsv_default
  lipid_scores_tsv <- lipid_scores_tsv_default
  outdir           <- outdir_default
  Z_col            <- NA_character_
  lipid_col        <- NA_character_
  cat("[INFO] No CLI args provided; using defaults.\n")
} else {
  stop("Usage: Rscript 06e_concordance_Z_vs_lipid.R <rna_scores_tsv> <lipid_scores_tsv> <outdir> <Z_col> <lipid_col>\n(Or run with no args for defaults.)\n")
}

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

cat("============================================================\n")
cat("[INFO] 06e_concordance_Z_vs_lipid.R\n")
cat("  RNA scores  : ", rna_scores_tsv, "\n", sep="")
cat("  Lipid table : ", lipid_scores_tsv, "\n", sep="")
cat("  outdir      : ", outdir, "\n", sep="")
cat("============================================================\n")

if (!file.exists(rna_scores_tsv)) stop("[ERROR] RNA scores file not found: ", rna_scores_tsv)
if (!file.exists(lipid_scores_tsv)) stop("[ERROR] Lipid table file not found: ", lipid_scores_tsv)

rna <- readr::read_tsv(rna_scores_tsv, show_col_types = FALSE)
lip <- readr::read_tsv(lipid_scores_tsv, show_col_types = FALSE)

# normalize group labels
norm_group <- function(x) {
  dplyr::case_when(
    str_to_upper(x) == "WT" ~ "WT",
    str_to_upper(x) %in% c("HO","KO") ~ "HO",
    TRUE ~ x
  )
}

# ensure batch/group exist (try minimal inference)
if (!"batch" %in% colnames(rna) && "sample_id" %in% colnames(rna)) rna$batch <- str_extract(rna$sample_id, "batch[0-9]+")
if (!"group" %in% colnames(rna) && "genotype" %in% colnames(rna)) rna <- rna %>% rename(group = genotype)
if (!all(c("batch","group") %in% colnames(rna))) stop("[ERROR] RNA table must have batch & group (or inferable).")

if (!"batch" %in% colnames(lip) && "sample_id" %in% colnames(lip)) lip$batch <- str_extract(lip$sample_id, "batch[0-9]+")
if (!"group" %in% colnames(lip) && "genotype" %in% colnames(lip)) lip <- lip %>% rename(group = genotype)
if (!all(c("batch","group") %in% colnames(lip))) stop("[ERROR] Lipid table must have batch & group (or inferable).")

rna <- rna %>% mutate(group = norm_group(group))
lip <- lip %>% mutate(group = norm_group(group))

## pick Z/lipid columns if not provided
if (is.na(Z_col) || nchar(Z_col) == 0) {
  Z_col <- Z_col_default_candidates[Z_col_default_candidates %in% colnames(rna)][1]
  if (is.na(Z_col) || nchar(Z_col) == 0) {
    # fallback 1: any "Membrane" axis-like column (exclude EphB1_)
    cand1 <- colnames(rna)[grepl("Membrane", colnames(rna)) & !grepl("^EphB1_", colnames(rna))]
    if (length(cand1) > 0) Z_col <- cand1[1]
  }
  if (is.na(Z_col) || nchar(Z_col) == 0) {
    # fallback 2: first non-id, non-EphB1_ column
    exclude <- c("sample_id","group","batch","rna_sample_id")
    cand2 <- setdiff(colnames(rna), exclude)
    cand2 <- cand2[!grepl("^EphB1_", cand2)]
    if (length(cand2) > 0) Z_col <- cand2[1]
  }
}
if (is.na(Z_col) || !Z_col %in% colnames(rna)) {
  stop("[ERROR] Z_col not found in RNA table. Provide as 4th arg.\nAvailable: ",
       paste(colnames(rna), collapse=", "))
}

if (is.na(lipid_col) || nchar(lipid_col) == 0) {
  lipid_col <- lipid_col_default_candidates[lipid_col_default_candidates %in% colnames(lip)][1]
}
if (is.na(lipid_col) || !lipid_col %in% colnames(lip)) {
  stop("[ERROR] lipid_col not found in lipid table. Provide as 5th arg.\nAvailable: ",
       paste(colnames(lip), collapse=", "))
}

effect_ho_wt <- function(x, g) {
  mean(as.numeric(x)[g=="HO"], na.rm = TRUE) - mean(as.numeric(x)[g=="WT"], na.rm = TRUE)
}

rna_eff <- rna %>%
  group_by(batch) %>%
  summarise(
    Z_effect = effect_ho_wt(.data[[Z_col]], group),
    n_WT = sum(group=="WT"), n_HO = sum(group=="HO"),
    .groups = "drop"
  )

lip_eff <- lip %>%
  group_by(batch) %>%
  summarise(
    lipid_effect = effect_ho_wt(.data[[lipid_col]], group),
    n_WT_lipid = sum(group=="WT"), n_HO_lipid = sum(group=="HO"),
    .groups = "drop"
  )

out <- rna_eff %>%
  inner_join(lip_eff, by = "batch") %>%
  mutate(
    direction_match = ifelse(
      is.na(Z_effect) | is.na(lipid_effect) | Z_effect == 0 | lipid_effect == 0,
      NA,
      sign(Z_effect) == sign(lipid_effect)
    )
  ) %>%
  select(batch, Z_effect, lipid_effect, direction_match,
         n_WT, n_HO, n_WT_lipid, n_HO_lipid)

out_path <- file.path(outdir, "concordance_Z_lipid.tsv")
readr::write_tsv(out, out_path)

cat("\n[OK] wrote: ", out_path, "\n", sep="")
cat("[INFO] Z_col used     : ", Z_col, "\n", sep="")
cat("[INFO] lipid_col used : ", lipid_col, "\n", sep="")
cat("============================================================\n")
cat("[DONE] 06e concordance done.\n")
cat("============================================================\n")