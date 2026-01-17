#!/usr/bin/env Rscript

## ============================================================
## 06.2a_mediation_rank_m1_membrane_context_robust.R
##
## 目的：
##   在样本量极小（n=5）前提下，以更可辩护的方式筛选/排序
##   EphB1 下游通路（M1）中“最有希望解释 membrane context 脂质轴障碍分数（Y）变化”的候选通路。
##
##   与 06.1 的区别：
##     - 不以 p 值/显著性为核心筛选（n 很小不可靠）
##     - 引入方向一致性 + LOO（leave-one-out）稳健性作为主判据
##     - 中介回归输出（a/b/c/c'/ab）仅作为辅助特征 join 进来
##
## 你的方向定义（本脚本默认遵循）：
##   X: KO/HO = 1, WT = 0
##   Y: membrane context_score（越高=障碍越严重）
##   M: pathway activation score（越高=激活越强/越好）
##
##   “保护性中介候选”在方向上的期望：
##     c_total  > 0   （KO 增加障碍）
##     a        < 0   （KO 降低通路激活）
##     b        < 0   （通路越激活，障碍越轻）
##     ab = a*b > 0   （间接效应为正）
##
## 输入：
##   1) mediation_summary.tsv（来自 06_...minimal.R）
##   2) sample_scores_for_mediation.tsv（来自 05_...sample_scores.R）
##
## 输出（默认 outdir = results/multiomics/tables/mediation）：
##   1) robust_rank_M1_membrane_context.tsv      （全量排序表）
##   2) robust_top_hits_M1_membrane_context.tsv  （通过方向+LOO 的 topN）
##
## 用法示例：
##   Rscript scripts/multi/06.2a_mediation_rank_m1_membrane_context_robust.R \
##     results/multiomics/tables/mediation/mediation_summary.tsv \
##     results/multiomics/tables/sample_scores/sample_scores_for_mediation.tsv \
##     results/multiomics/tables/mediation \
##     "Membrane context" \
##     5
##
## ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
})

## ---------------- 0. 默认路径 & 参数 ------------------------

default_med_infile   <- "results/multiomics/tables/mediation/mediation_summary.tsv"
default_score_infile <- "results/multiomics/tables/sample_scores/sample_scores_for_mediation.tsv"
default_outdir       <- "results/multiomics/tables/mediation"
default_y_axis       <- "Membrane context"
default_top_n        <- 5
default_min_loo_ok   <- 4  # 5 折 LOO 中至少 4 折方向一致

args <- commandArgs(trailingOnly = TRUE)

med_infile   <- if (length(args) >= 1) args[1] else default_med_infile
score_infile <- if (length(args) >= 2) args[2] else default_score_infile
outdir       <- if (length(args) >= 3) args[3] else default_outdir
y_axis_keep  <- if (length(args) >= 4) args[4] else default_y_axis
top_n        <- if (length(args) >= 5) as.integer(args[5]) else default_top_n

min_loo_ok <- default_min_loo_ok

cat("============================================================\n")
cat("[INFO] 06.2a_mediation_rank_m1_membrane_context_robust.R\n")
cat("  mediation_summary : ", med_infile, "\n", sep = "")
cat("  sample_scores     : ", score_infile, "\n", sep = "")
cat("  outdir            : ", outdir, "\n", sep = "")
cat("  Y_axis keep       : ", y_axis_keep, "\n", sep = "")
cat("  top_n             : ", top_n, "\n", sep = "")
cat("  min_loo_ok        : ", min_loo_ok, " / 5\n", sep = "")
cat("============================================================\n\n")

if (!file.exists(med_infile))   stop("[ERROR] 找不到 mediation_summary.tsv: ", med_infile)
if (!file.exists(score_infile)) stop("[ERROR] 找不到 sample_scores_for_mediation.tsv: ", score_infile)

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

## ---------------- 1. 读取数据 -------------------------------

cat("[STEP1] 读取 mediation_summary.tsv...\n")
med <- readr::read_tsv(med_infile, show_col_types = FALSE)

needed_med_cols <- c(
  "mediator_col", "outcome_col",
  "n", "group_ref", "group_focal",
  "a", "b", "c_total", "c_prime",
  "ab", "ab_boot_mean", "ab_boot_ci_low", "ab_boot_ci_high",
  "p_sobel"
)
missing_med <- setdiff(needed_med_cols, colnames(med))
if (length(missing_med) > 0) {
  stop("[ERROR] mediation_summary.tsv 缺少列：", paste(missing_med, collapse = ", "))
}

cat("[STEP2] 读取 sample_scores_for_mediation.tsv...\n")
sc <- readr::read_tsv(score_infile, show_col_types = FALSE)

needed_sc_cols <- c("sample_id", "group", "batch")
missing_sc <- setdiff(needed_sc_cols, colnames(sc))
if (length(missing_sc) > 0) {
  stop("[ERROR] sample_scores_for_mediation.tsv 缺少列：", paste(missing_sc, collapse = ", "))
}

## ---------------- 2. 定义 Y 列（只保留 membrane context）----

cat("[STEP3] 定位 Y 列（outcome_col）...\n")

sc_num_cols <- names(sc)[vapply(sc, is.numeric, logical(1))]

# Y 候选：以 _score 结尾的数值列
y_candidates <- sc_num_cols[str_detect(sc_num_cols, "_score$")]
if (length(y_candidates) == 0) stop("[ERROR] sample_scores_for_mediation.tsv 中未找到任何 *_score 数值列。")

# 在 mediation_summary 中，Y_axis = 去掉 _score 的 outcome_col
med2 <- med %>%
  mutate(
    Y_axis = str_remove(outcome_col, "_score$")
  )

# 只保留指定 Y_axis
med_y <- med2 %>% filter(Y_axis == y_axis_keep)

if (nrow(med_y) == 0) {
  stop("[ERROR] mediation_summary.tsv 中找不到 Y_axis == '", y_axis_keep, "' 的记录。\n",
       "  你可以检查 outcome_col 的具体列名，或将第 4 个参数改成实际的 Y_axis。")
}

# 在 sample_scores 中定位对应 outcome_col（必须存在）
outcome_col_keep <- unique(med_y$outcome_col)
if (length(outcome_col_keep) != 1) {
  stop("[ERROR] 期望指定 Y_axis 对应唯一 outcome_col，但发现：",
       paste(outcome_col_keep, collapse = ", "))
}
Y_col <- outcome_col_keep

if (!Y_col %in% colnames(sc)) {
  stop("[ERROR] sample_scores 表中找不到 outcome_col 列：", Y_col)
}

cat("  [OK] Y_col = ", Y_col, "\n", sep = "")
cat("  [INFO] 用于排序的 mediation 记录数（含不同 mediator）：", nrow(med_y), "\n", sep = "")

## ---------------- 3. 定义 M1 / M2（复用 06.1 的规则） ---------

axis_prefixes <- c(
  "Synthesis_",
  "Remodeling_",
  "Oxidation_",
  "Supply_",
  "Transport_",
  "Membrane_context_"
)
axis_pattern <- paste0("^(", paste(axis_prefixes, collapse = "|"), ")")

med_y <- med_y %>%
  mutate(
    M_is_axis = str_detect(mediator_col, axis_pattern),
    M_type = ifelse(M_is_axis, "M2_axis_upstream", "M1_ephb1_downstream_or_other")
  )

## ---------------- 4. 把 sample_scores 转成长表用于取 M 列 -----

cat("[STEP4] 准备 sample-level 的 X / Y / M 数据...\n")

# X 编码：尽可能沿用你此前的命名习惯（KO/HO 为 1，WT 为 0）
is_focal <- function(g) {
  # 你可以按需扩充关键字
  str_detect(tolower(g), "ho|ko|knockout|ephb1-/-|ephb1_ko")
}

sc2 <- sc %>%
  mutate(
    X = ifelse(is_focal(group), 1, 0)
  )

if (length(unique(sc2$X)) < 2) {
  warning("[WARN] X 只有一个水平（全是 0 或全是 1）。请检查 group 命名是否包含 KO/HO 关键词。")
}

# 提取所有 mediator 列（数值列，且不是 *_score）
m_candidates <- setdiff(sc_num_cols, y_candidates)
if (length(m_candidates) == 0) stop("[ERROR] sample_scores 表中未找到任何 mediator 数值列（非 *_score）。")

## ---------------- 5. 定义 LOO 计算函数 -----------------------

loo_metrics <- function(df, m_col, y_col) {
  # df 必须含 X, 以及 m_col / y_col 两列
  n <- nrow(df)
  if (n < 4) {
    return(list(
      loo_r_neg_n = NA_integer_,
      loo_dM_neg_n = NA_integer_,
      loo_both_ok_n = NA_integer_,
      loo_r_median = NA_real_,
      loo_r_minabs = NA_real_
    ))
  }

  r_list <- numeric(0)
  dm_list <- numeric(0)
  ok_r <- integer(0)
  ok_dm <- integer(0)
  ok_both <- integer(0)

  for (i in seq_len(n)) {
    sub <- df[-i, , drop = FALSE]
    r <- suppressWarnings(cor(sub[[m_col]], sub[[y_col]], method = "spearman", use = "pairwise.complete.obs"))
    dm <- mean(sub[[m_col]][sub$X == 1], na.rm = TRUE) - mean(sub[[m_col]][sub$X == 0], na.rm = TRUE)

    r_list <- c(r_list, r)
    dm_list <- c(dm_list, dm)

    ok_r <- c(ok_r, ifelse(!is.na(r)  && r < 0, 1L, 0L))
    ok_dm <- c(ok_dm, ifelse(!is.na(dm) && dm < 0, 1L, 0L))
    ok_both <- c(ok_both, ifelse(ok_r[length(ok_r)] == 1L && ok_dm[length(ok_dm)] == 1L, 1L, 0L))
  }

  list(
    loo_r_neg_n    = sum(ok_r),
    loo_dM_neg_n   = sum(ok_dm),
    loo_both_ok_n  = sum(ok_both),
    loo_r_median   = median(r_list, na.rm = TRUE),
    loo_r_minabs   = min(abs(r_list), na.rm = TRUE)
  )
}

## ---------------- 6. 对每个 M1 mediator 计算稳健性指标 ---------

cat("[STEP5] 计算每个 mediator 的 ΔM / 相关 / LOO 稳健性...\n")

# 只对 M1 做排序（你的目标是 EphB1 downstream）
med_m1 <- med_y %>% filter(M_type == "M1_ephb1_downstream_or_other")

if (nrow(med_m1) == 0) {
  stop("[ERROR] 在 Y_axis == '", y_axis_keep, "' 下没有任何 M1 mediator 记录。")
}

# 去重：每个 mediator_col 在同一 outcome_col 下应唯一
med_m1 <- med_m1 %>%
  distinct(mediator_col, outcome_col, .keep_all = TRUE)

# 确认这些 mediator 列在 sample_scores 里存在
missing_m_cols <- setdiff(med_m1$mediator_col, colnames(sc2))
if (length(missing_m_cols) > 0) {
  stop("[ERROR] 下列 mediator_col 在 sample_scores 表中不存在：\n  ",
       paste(missing_m_cols, collapse = ", "),
       "\n请检查 05 输出列名是否与 06 的 mediator_col 一致。")
}

# 计算每个 mediator 的样本级指标
rows <- lapply(med_m1$mediator_col, function(mcol) {
  df <- sc2 %>% select(sample_id, group, batch, X, Y = all_of(Y_col), M = all_of(mcol))

  dY <- mean(df$Y[df$X == 1], na.rm = TRUE) - mean(df$Y[df$X == 0], na.rm = TRUE)
  dM <- mean(df$M[df$X == 1], na.rm = TRUE) - mean(df$M[df$X == 0], na.rm = TRUE)

  r_spear <- suppressWarnings(cor(df$M, df$Y, method = "spearman", use = "pairwise.complete.obs"))
  r_pear  <- suppressWarnings(cor(df$M, df$Y, method = "pearson",  use = "pairwise.complete.obs"))

  loo <- loo_metrics(df, "M", "Y")

  tibble(
    mediator_col = mcol,
    outcome_col  = Y_col,
    Y_axis       = y_axis_keep,
    n_samples    = nrow(df),
    dY_HO_minus_WT = dY,
    dM_HO_minus_WT = dM,
    cor_spearman_M_Y = r_spear,
    cor_pearson_M_Y  = r_pear,
    loo_r_neg_n      = loo$loo_r_neg_n,
    loo_dM_neg_n     = loo$loo_dM_neg_n,
    loo_both_ok_n    = loo$loo_both_ok_n,
    loo_r_median     = loo$loo_r_median,
    loo_r_minabs     = loo$loo_r_minabs
  )
})

rob <- bind_rows(rows)

## ---------------- 7. join mediation（a/b/c/c'/ab）并打标签 ----

cat("[STEP6] join mediation a/b/c/c'/ab，并计算方向一致性标签...\n")

rob2 <- rob %>%
  left_join(
    med_m1 %>%
      select(
        mediator_col, outcome_col,
        group_ref, group_focal,
        a, b, c_total, c_prime,
        ab, ab_boot_mean, ab_boot_ci_low, ab_boot_ci_high,
        p_sobel
      ),
    by = c("mediator_col", "outcome_col")
  ) %>%
  mutate(
    # 你定义的“保护性候选”方向规则（注意：mutate 里要用向量化的 & / |，不能用 && / ||）
    dir_c_ok  = !is.na(c_total) & (c_total > 0),
    dir_a_ok  = !is.na(a)       & (a < 0),
    dir_b_ok  = !is.na(b)       & (b < 0),
    dir_ab_ok = !is.na(ab_boot_mean) & (ab_boot_mean > 0), # 用 bootstrap mean 作为稳定版

    loo_ok = !is.na(loo_both_ok_n) & (loo_both_ok_n >= min_loo_ok),

    # 解释削弱（仅作排序特征，不作因果比例宣称）
    shrink_ratio = ifelse(!is.na(c_total) & (c_total != 0) & !is.na(c_prime),
                          (c_total - c_prime) / c_total,
                          NA_real_),

    # 总体通过：方向 + LOO
    pass_protective = dir_c_ok & dir_a_ok & dir_b_ok & dir_ab_ok & loo_ok
  )

## ---------------- 8. 排序规则（rank-sum，稳健、易辩护） -------

cat("[STEP7] 生成排序分数（rank-sum）...\n")

# 排序偏好：
#   1) shrink_ratio 越大越好（解释力）
#   2) cor_spearman 越负越好（保护性关联）
#   3) dM 越负越好（KO 使通路变弱）
#   4) LOO 稳健性越强越好（loo_both_ok_n）
rob3 <- rob2 %>%
  mutate(
    rank_shrink = rank(-shrink_ratio, ties.method = "average", na.last = "keep"),
    rank_cor    = rank(cor_spearman_M_Y, ties.method = "average", na.last = "keep"), # 越小(更负)越好
    rank_dM     = rank(dM_HO_minus_WT,   ties.method = "average", na.last = "keep"), # 越小(更负)越好
    rank_loo    = rank(-loo_both_ok_n,   ties.method = "average", na.last = "keep"),
    rank_sum    = rank_shrink + rank_cor + rank_dM + rank_loo
  ) %>%
  arrange(rank_sum)

## ---------------- 9. 写出输出表 ------------------------------

out_all <- file.path(outdir, "robust_rank_M1_membrane_context.tsv")
out_top <- file.path(outdir, "robust_top_hits_M1_membrane_context.tsv")

cat("[STEP8] 写出结果...\n")
readr::write_tsv(rob3, out_all)

top_tbl <- rob3 %>%
  filter(pass_protective) %>%
  arrange(rank_sum) %>%
  slice_head(n = top_n)

readr::write_tsv(top_tbl, out_top)

cat("  [OK] 全量排序表: ", out_all, "\n", sep = "")
cat("  [OK] 通过方向+LOO 的 top hits: ", out_top, "\n", sep = "")

cat("\n[SUMMARY]\n")
cat("  Y_col: ", Y_col, "\n", sep = "")
cat("  n mediators (M1): ", nrow(rob3), "\n", sep = "")
cat("  pass_protective: ", sum(rob3$pass_protective, na.rm = TRUE), "\n", sep = "")
cat("  top_n written: ", nrow(top_tbl), "\n", sep = "")
cat("============================================================\n")
cat("[DONE] 06.2a robust ranking 完成。\n")
cat("============================================================\n")