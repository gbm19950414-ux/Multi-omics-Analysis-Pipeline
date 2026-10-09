#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr); library(purrr)
})

# ------------------------------------------------------------------
# Usage:
#   07_merge_counts.R <processed_dir> <out_counts_tsv> [out_qc_tsv]
# Notes:
#   - Auto-discovers all batch count matrices under <processed_dir>
#     matching "featureCounts_matrix.tsv".
#   - Merges by GeneID with full join (no gene dropped).
#   - Cleans duplicate sample names by prefixing batch id.
#   - Optionally merges per-batch QC summaries (qc_basic.tsv) into one table
#     if a third argument is provided.
#   - Compatible with both STAR and HISAT2 pipelines (aligner-agnostic).
# ------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("用法: 07_merge_counts.R <processed_dir> <out_counts_tsv> [out_qc_tsv]")
}
proc_dir <- args[1]
out_counts <- args[2]
out_qc <- if (length(args) >= 3) args[3] else ""

message("[INFO] processed_dir = ", proc_dir)
message("[INFO] out_counts   = ", out_counts)
if (nzchar(out_qc)) message("[INFO] out_qc       = ", out_qc)

# Ensure output parent directories exist and guard against placeholder paths
oc_dir <- dirname(out_counts)
dir.create(oc_dir, recursive = TRUE, showWarnings = FALSE)
if (nzchar(out_qc)) {
  oq_dir <- dirname(out_qc)
  dir.create(oq_dir, recursive = TRUE, showWarnings = FALSE)
}

# Prevent accidental placeholder paths like "/.../" which will fail to write
if (stringr::str_detect(out_counts, "/\\.\\.\\./") || (nzchar(out_qc) && stringr::str_detect(out_qc, "/\\.\\.\\./"))) {
  stop("输出路径包含 '...' 占位符：请改成真实绝对路径。")
}

# ---------- helpers ----------
# infer batch name from a path like: .../data/processed/rna/<batch>/counts/featureCounts_matrix.tsv
infer_batch <- function(path) {
  m <- str_match(path, "/rna/([^/]+)/counts/")
  if (is.na(m[1,2])) return(NA_character_)
  m[1,2]
}

# Parse sample name from BAM-like path; prefer .../(hisat2|star)/SAMPLE/Aligned...
parse_sample_from_bam <- function(p) {
  pp <- str_replace_all(p, "\\\\", "/")
  if (str_detect(pp, "/(star|hisat2)/[^/]+/Aligned")) {
    return(str_replace(pp, ".*/(star|hisat2)/([^/]+)/Aligned.*$", "\\2"))
  }
  # fallback: basename without .bam
  str_remove(basename(pp), "\\.bam$")
}

# From a featureCounts.txt header, extract BAM paths and convert to sample names
# fc_txt: path to featureCounts.txt in the same counts directory
# n_cols: number of data columns (excluding GeneID) that matrix has
sample_names_from_fc_header <- function(fc_txt, n_cols) {
  if (!file.exists(fc_txt)) return(NULL)
  hdr <- read_lines(fc_txt, n_max = 2)
  if (length(hdr) < 2) return(NULL)
  # featureCounts header: first 6 columns are meta, following are BAM paths
  # Create enough splits: 6 meta + n_cols sample columns
  parts <- str_split_fixed(hdr[2], "\t", 6 + n_cols)
  if (ncol(parts) < (6 + n_cols)) return(NULL)
  bam_paths <- as.vector(parts[, -(1:6)])
  # Convert each BAM path to sample name
  vapply(bam_paths, parse_sample_from_bam, character(1))
}

# Make sample-like column names tidy and unique
tidy_colnames <- function(nms, path_vec = NULL) {
  # If columns already look like sample IDs (WT_1 etc.), keep them.
  # Else, try to extract sample from BAM-like paths if any.
  nn <- nms
  # Strip absolute paths, keep basename-ish tails
  nn <- str_replace_all(nn, "\\\\", "/")
  nn <- str_replace(nn, ".*/(star|hisat2)/([^/]+)/Aligned.*$", "\\2")
  # Remove common suffixes
  nn <- str_replace_all(nn, "\\.bam$", "")
  nn <- str_replace_all(nn, "\\.sortedByCoord\\.out$", "")
  # Ensure syntactically valid and unique
  nn <- make.unique(nn, sep = "__dup")
  nn
}

# ---------- discover matrices ----------
mat_files <- list.files(proc_dir, pattern = "featureCounts_matrix.tsv$", recursive = TRUE, full.names = TRUE)
if (length(mat_files) == 0) {
  stop("未找到任何 featureCounts_matrix.tsv 于: ", proc_dir)
}

message("[INFO] merging matrices:")
for (f in mat_files) message("  - ", f)

# read all, ensure GeneID is character
dfs <- lapply(mat_files, function(f) {
  df <- read_tsv(f, show_col_types = FALSE)
  if (!"GeneID" %in% names(df)) stop("文件缺少 GeneID 列: ", f)
  df$GeneID <- as.character(df$GeneID)

  # Determine batch from path
  b <- infer_batch(f)

  # Try to find featureCounts.txt in the same counts directory and derive sample names
  counts_dir <- dirname(f)
  fc_txt <- file.path(counts_dir, "featureCounts.txt")
  n_data_cols <- ncol(df) - 1
  samp <- sample_names_from_fc_header(fc_txt, n_data_cols)

  if (is.null(samp) || length(samp) != n_data_cols) {
    # Fallback: try to tidy existing column names
    samp <- tidy_colnames(names(df)[-1])
  }

  # Prefix with batch to avoid duplicates across batches when available
  if (!is.na(b) && nzchar(b)) {
    samp <- paste0(b, "_", samp)
  }

  names(df) <- c("GeneID", samp)

  # attach batch attribute for downstream (if needed)
  attr(df, "batch") <- b
  df
})

# progressive full join by GeneID
merged <- Reduce(function(x, y) full_join(x, y, by = "GeneID"), dfs)

# Column names already set to batch_sample; just ensure uniqueness
names(merged) <- make.unique(names(merged), sep = "__dup")

# write merged counts
write_tsv(merged, out_counts)
message("[OK ] merged counts -> ", out_counts)

# ---------- optional: merge QC ----------
if (nzchar(out_qc)) {
  qc_files <- list.files(proc_dir, pattern = "qc_basic.tsv$", recursive = TRUE, full.names = TRUE)
  if (length(qc_files) == 0) {
    warning("未找到任何 qc_basic.tsv，将跳过 QC 合并")
  } else {
    message("[INFO] merging QC summaries:")
    for (f in qc_files) message("  - ", f)
    qc_list <- lapply(qc_files, function(f) {
      b <- infer_batch(f)
      q <- suppressMessages(read_tsv(f, show_col_types = FALSE))
      # Normalize columns
      if (!all(c("sample","assigned","total","assigned_rate") %in% names(q))) {
        stop("QC 文件缺失列: ", f)
      }
      # Standardize sample name
      q <- q %>%
        mutate(
          batch = b,
          sample = str_replace_all(sample, "\\\\", "/"),
          sample = str_replace(sample, ".*/(star|hisat2)/([^/]+)/Aligned.*$", "\\2"),
          sample = str_replace(sample, "^.+/", ""),
          assigned = as.numeric(assigned),
          total = as.numeric(total),
          assigned_rate = as.numeric(assigned_rate)
        )
      q
    })
    qc_merged <- bind_rows(qc_list) %>%
      select(batch, sample, total, assigned, assigned_rate) %>%
      arrange(batch, sample)
    write_tsv(qc_merged, out_qc)
    message("[OK ] merged QC -> ", out_qc)

    # quick stats
    stats <- qc_merged %>%
      group_by(batch) %>%
      summarise(n = n(), mean_map = mean(assigned_rate, na.rm = TRUE)) %>%
      ungroup()
    message("[STATS] samples per batch & mean assigned_rate:")
    capture.output(print(stats, n = Inf)) %>% paste(collapse = "\n") %>% message()
  }
}

# final echo
message("[SUMMARY] matrices: ", length(mat_files),
        "; samples: ", ncol(merged) - 1,
        if (nzchar(out_qc)) paste0("; QC rows: ", if (exists("qc_merged")) nrow(qc_merged) else 0) else "")
