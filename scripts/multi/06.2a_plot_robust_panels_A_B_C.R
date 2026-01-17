#!/usr/bin/env Rscript

## ============================================================
## 06.2a_plot_robust_panels_A_B_C.R
##
## 输入：
##   1) robust_rank_M1_membrane_context.tsv
##   2) sample_scores_for_mediation.tsv
##
## 输出：
##   - robust_06.2a_panels_ABC.png / .pdf
##   - robust_06.2a_panel_A_deltaM.png
##   - robust_06.2a_panel_B_scatter.png
##   - robust_06.2a_panel_C_LOO_heatmap.png
##
## 三联图：
##   A: ΔM（HO - WT）条形图（按 rank_sum 排序）
##   B: scatter（Top K mediators：M vs Y，按组上色）
##   C: LOO heatmap（每次 leave-one-out 的 Spearman r）
##
## 用法：
##   Rscript scripts/multi/06.2a_plot_robust_panels_A_B_C.R \
##     results/multiomics/tables/mediation/robust_rank_M1_membrane_context.tsv \
##     results/multiomics/tables/sample_scores/sample_scores_for_mediation.tsv \
##     results/multiomics/tables/mediation \
##     3
##
## ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(ggplot2)
})

## ---------------- 参数 ----------------
default_rank_infile  <- "results/multiomics/tables/mediation/robust_rank_M1_membrane_context.tsv"
default_score_infile <- "results/multiomics/tables/sample_scores/sample_scores_for_mediation.tsv"
default_outdir       <- "results/multiomics/plots/mediation"
default_top_k_scatter <- 3

args <- commandArgs(trailingOnly = TRUE)
rank_infile   <- if (length(args) >= 1) args[1] else default_rank_infile
score_infile  <- if (length(args) >= 2) args[2] else default_score_infile
outdir        <- if (length(args) >= 3) args[3] else default_outdir
top_k_scatter <- if (length(args) >= 4) as.integer(args[4]) else default_top_k_scatter

cat("============================================================\n")
cat("[INFO] 06.2a_plot_robust_panels_A_B_C.R\n")
cat("  rank_table      : ", rank_infile, "\n", sep = "")
cat("  sample_scores   : ", score_infile, "\n", sep = "")
cat("  outdir          : ", outdir, "\n", sep = "")
cat("  top_k_scatter   : ", top_k_scatter, "\n", sep = "")
cat("============================================================\n\n")

if (!file.exists(rank_infile))  stop("[ERROR] 找不到 rank_table: ", rank_infile)
if (!file.exists(score_infile)) stop("[ERROR] 找不到 sample_scores: ", score_infile)
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

## ---------------- 读取数据 ----------------
rank_tbl <- readr::read_tsv(rank_infile, show_col_types = FALSE)
sc       <- readr::read_tsv(score_infile, show_col_types = FALSE)

needed_rank <- c("mediator_col","outcome_col","Y_axis","dM_HO_minus_WT","cor_spearman_M_Y",
                 "loo_both_ok_n","rank_sum")
miss_rank <- setdiff(needed_rank, colnames(rank_tbl))
if (length(miss_rank) > 0) stop("[ERROR] rank_table 缺少列：", paste(miss_rank, collapse = ", "))

needed_sc <- c("sample_id","group","batch")
miss_sc <- setdiff(needed_sc, colnames(sc))
if (length(miss_sc) > 0) stop("[ERROR] sample_scores 缺少列：", paste(miss_sc, collapse = ", "))

# X 编码（与 06.2a 一致）
is_focal <- function(g) str_detect(tolower(g), "ho|ko|knockout|ephb1-/-|ephb1_ko")
sc <- sc %>% mutate(X = ifelse(is_focal(group), 1, 0))

# 只画一个 outcome_col（通常就是 Membrane context_score）
Y_col <- unique(rank_tbl$outcome_col)
if (length(Y_col) != 1) {
  stop("[ERROR] rank_table 中 outcome_col 不唯一：", paste(Y_col, collapse = ", "))
}
Y_col <- Y_col[1]
if (!Y_col %in% colnames(sc)) stop("[ERROR] sample_scores 中找不到 Y_col：", Y_col)

## ---------------- Panel A：ΔM barplot ----------------
A_df <- rank_tbl %>%
  arrange(rank_sum) %>%
  mutate(mediator_col = factor(mediator_col, levels = mediator_col))

pA <- ggplot(A_df, aes(x = mediator_col, y = dM_HO_minus_WT)) +
  geom_col() +
  coord_flip() +
  labs(
    title = "A. ΔM (HO - WT) by mediator (sorted)",
    x = NULL,
    y = "ΔM = mean(M_HO) - mean(M_WT)"
  ) +
  theme_bw(base_size = 11)

## ---------------- Panel B：scatter（Top K mediators） ----------------
# dplyr::n() cannot be used here; compute a constant K first
K <- min(top_k_scatter, nrow(rank_tbl))

top_mediators <- rank_tbl %>%
  arrange(rank_sum) %>%
  slice_head(n = K) %>%
  pull(mediator_col)

B_long <- sc %>%
  select(sample_id, group, X, Y = all_of(Y_col), all_of(top_mediators)) %>%
  pivot_longer(cols = all_of(top_mediators), names_to = "mediator_col", values_to = "M") %>%
  mutate(
    mediator_col = factor(mediator_col, levels = top_mediators),
    group = factor(group)
  )

pB <- ggplot(B_long, aes(x = M, y = Y, color = group)) +
  geom_point(size = 2.4) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~ mediator_col, scales = "free_x") +
  labs(
    title = paste0("B. Scatter (Top ", length(top_mediators), " mediators): M vs Y"),
    x = "Mediator score (M)",
    y = paste0("Y = ", Y_col)
  ) +
  theme_bw(base_size = 11)

## ---------------- Panel C：LOO heatmap（每次 leave-one-out 的 Spearman r） ----------------
loo_spearman <- function(df, m_col, y_col) {
  n <- nrow(df)
  out <- vector("list", n)
  for (i in seq_len(n)) {
    sub <- df[-i, , drop = FALSE]
    r <- suppressWarnings(cor(sub[[m_col]], sub[[y_col]], method = "spearman", use = "pairwise.complete.obs"))
    out[[i]] <- tibble(left_out = df$sample_id[i], r_spearman = r)
  }
  bind_rows(out)
}

C_rows <- lapply(as.character(rank_tbl$mediator_col), function(mcol) {
  df <- sc %>% select(sample_id, X, Y = all_of(Y_col), M = all_of(mcol))
  loo <- loo_spearman(df, "M", "Y")
  loo %>% mutate(mediator_col = mcol)
})

C_df <- bind_rows(C_rows) %>%
  left_join(
    rank_tbl %>% select(mediator_col, rank_sum),
    by = "mediator_col"
  ) %>%
  mutate(
    mediator_col = factor(mediator_col, levels = rank_tbl %>% arrange(rank_sum) %>% pull(mediator_col)),
    left_out = factor(left_out, levels = unique(sc$sample_id))
  )

pC <- ggplot(C_df, aes(x = left_out, y = mediator_col, fill = r_spearman)) +
  geom_tile() +
  geom_text(aes(label = ifelse(is.na(r_spearman), "NA", sprintf("%.2f", r_spearman))), size = 3) +
  scale_fill_gradient2() +
  labs(
    title = "C. LOO robustness heatmap (Spearman r of M vs Y)",
    x = "Left-out sample",
    y = "Mediator (sorted by rank_sum)",
    fill = "Spearman r"
  ) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## ---------------- 合并三联图（优先 patchwork；否则分别保存） ----------------
out_A <- file.path(outdir, "robust_06.2a_panel_A_deltaM.png")
out_B <- file.path(outdir, "robust_06.2a_panel_B_scatter.png")
out_C <- file.path(outdir, "robust_06.2a_panel_C_LOO_heatmap.png")

ggsave(out_A, pA, width = 8, height = 4.8, dpi = 300)
ggsave(out_B, pB, width = 10, height = 5.2, dpi = 300)
ggsave(out_C, pC, width = 10, height = 5.5, dpi = 300)

cat("[OK] 写出单图：\n")
cat("  ", out_A, "\n", sep = "")
cat("  ", out_B, "\n", sep = "")
cat("  ", out_C, "\n", sep = "")

# 合并图（如果 patchwork 或 cowplot 可用）
out_ABC_png <- file.path(outdir, "robust_06.2a_panels_ABC.png")
out_ABC_pdf <- file.path(outdir, "robust_06.2a_panels_ABC.pdf")

made_combo <- FALSE
if (requireNamespace("patchwork", quietly = TRUE)) {
  cat("[INFO] 使用 patchwork 合并三联图。\n")
  library(patchwork)
  pABC <- pA / pB / pC + plot_layout(heights = c(1.1, 1.2, 1.25))
  ggsave(out_ABC_png, pABC, width = 11, height = 14, dpi = 300)
  ggsave(out_ABC_pdf, pABC, width = 11, height = 14)
  made_combo <- TRUE
} else if (requireNamespace("cowplot", quietly = TRUE)) {
  cat("[INFO] 使用 cowplot 合并三联图。\n")
  library(cowplot)
  pABC <- cowplot::plot_grid(pA, pB, pC, ncol = 1, rel_heights = c(1.1, 1.2, 1.25))
  ggsave(out_ABC_png, pABC, width = 11, height = 14, dpi = 300)
  ggsave(out_ABC_pdf, pABC, width = 11, height = 14)
  made_combo <- TRUE
} else {
  cat("[WARN] 未检测到 patchwork/cowplot：将只保留单图 PNG。\n")
}

if (made_combo) {
  cat("[OK] 三联图输出：\n")
  cat("  ", out_ABC_png, "\n", sep = "")
  cat("  ", out_ABC_pdf, "\n", sep = "")
}

cat("============================================================\n")
cat("[DONE] 06.2a panels A/B/C finished.\n")
cat("============================================================\n")