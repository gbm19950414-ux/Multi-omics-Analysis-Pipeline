#!/usr/bin/env Rscript

## ============================================================
## 06_mediation_analysis_minimal.R
##
## 目的：
##   对 X -> M -> Y 进行一个极简但可扩展的线性中介分析。
##
##   - X: group（二分类，如 WT / HO 或 WT / KO）
##   - M: 中介变量（例如 RNA 通路分数、上游调控通路分数）
##   - Y: 结果变量（例如 CL 相关多轴脂质分数 *_score）
##
## 使用模式：
##
##   1) 单一 M/Y 模式（与旧版兼容）：
##        Rscript scripts/multi/06_mediation_analysis_minimal.R \\
##          sample_scores_for_mediation.tsv \\
##          M_col_name \\
##          Y_col_name \\
##          outdir
##
##      - 其中 M_col_name / Y_col_name 必须是表中的列名。
##
##   2) 自动网格模式（推荐探索阶段使用）：
##        Rscript scripts/multi/06_mediation_analysis_minimal.R \\
##          sample_scores_for_mediation.tsv \\
##          outdir
##
##      - 此时脚本会自动：
##          * 把所有数值型且不以 "_score" 结尾、且不属于
##            (sample_id, group, batch, rna_sample_id, lipid_sample_id)
##            的列当作候选 M（上游通路 / 轴等）。
##          * 把所有列名以 "_score" 结尾的数值型列当作候选 Y
##            （例如：Synthesis_score / Supply_score …）。
##          * 对每一个 (M, Y) 组合做一次中介分析。
##
## 输入（共同要求）：
##   - sample_scores_for_mediation.tsv 必须至少包含：
##       * sample_id
##       * group（二分类）
##       * 若干数值型 M 候选列
##       * 若干数值型 Y 候选列（对自动模式尤其重要）
##
## 输出（两种模式一致）：
##   - <outdir>/mediation_summary.tsv
##       每一行对应一组 (M, Y) 组合，包含：
##         * n
##         * group_ref, group_focal
##         * mediator_col, outcome_col
##         * a, se_a, p_a
##         * b, se_b, p_b
##         * c_total, se_c_total, p_c_total        （总效应）
##         * c_prime, se_c_prime, p_c_prime        （直接效应）
##         * ab, se_ab_sobel, z_sobel, p_sobel     （间接效应 + Sobel）
##         * ab_boot_mean, ab_boot_ci_low, ab_boot_ci_high
##
##   - <outdir>/mediation_bootstrap_ab.tsv
##       每一行对应一次 bootstrap 抽样：
##         * mediator_col
##         * outcome_col
##         * iter
##         * ab
##
## 统计方法说明：
##   - 线性模型：
##       * a:     M ~ X
##       * c:     Y ~ X
##       * b,c':  Y ~ X + M
##   - 间接效应：ab = a * b
##   - Sobel 检验：使用 se(a)、se(b) 组合估计 se(ab)
##   - Bootstrap：对样本编号做 B=5000 次有放回重采样，
##     每次重新估计 a、b 并记录 a*b，给出经验均值与 95% CI。
##
## 注意：
##   - n 很小（例如 5）时，所有 p 值与 CI 仅作“方向性提示”，
##     解读时必须依赖生物学合理性与多组学一致性。
## ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
})

default_sample_scores_tsv <- "results/multiomics/tables/sample_scores/sample_scores_for_mediation.tsv"
default_outdir           <- "results/multiomics/tables/mediation"

args <- commandArgs(trailingOnly = TRUE)

if (! (length(args) %in% c(0, 1, 2, 4)) ) {
  stop(
    "[ERROR] 参数数量不正确。\n",
    "用法一（单一 M/Y 模式，与旧版兼容）：\n",
    "  Rscript scripts/multi/06_mediation_analysis_minimal.R \\\n",
    "    sample_scores_for_mediation.tsv M_col_name Y_col_name outdir\n\n",
    "用法二（自动网格模式，显式指定输入表）：\n",
    "  Rscript scripts/multi/06_mediation_analysis_minimal.R \\\n",
    "    sample_scores_for_mediation.tsv outdir\n\n",
    "用法三（自动网格模式，默认使用固定输入表，输出也固定）：\n",
    "  Rscript scripts/multi/06_mediation_analysis_minimal.R\n",
    "  此时默认输入: ", default_sample_scores_tsv, "\n",
    "       默认输出: ", default_outdir, "\n",
    call. = FALSE
  )
}

if (length(args) == 4) {
  ## 单一 M/Y 模式
  mode               <- "single"
  sample_scores_tsv  <- args[1]
  mediator_cols      <- args[2]
  outcome_cols       <- args[3]
  outdir             <- args[4]
} else if (length(args) == 2) {
  ## 自动网格模式，显式指定输入表
  mode               <- "grid"
  sample_scores_tsv  <- args[1]
  outdir             <- args[2]
  mediator_cols      <- character(0)  # 占位，后面自动识别
  outcome_cols       <- character(0)
} else if (length(args) == 1) {
  ## 自动网格模式，默认输入表，手动指定输出目录
  mode               <- "grid"
  sample_scores_tsv  <- default_sample_scores_tsv
  outdir             <- args[1]
  mediator_cols      <- character(0)
  outcome_cols       <- character(0)
} else {
  ## length(args) == 0: 自动网格模式，输入和输出都使用固定默认路径
  mode               <- "grid"
  sample_scores_tsv  <- default_sample_scores_tsv
  outdir             <- default_outdir
  mediator_cols      <- character(0)
  outcome_cols       <- character(0)
}

cat("============================================================\n")
cat("[INFO] 06_mediation_analysis_minimal.R\n")
cat("  mode              : ", mode,              "\n", sep = "")
cat("  sample_scores_tsv : ", sample_scores_tsv, "\n", sep = "")
if (mode == "single") {
  cat("  mediator_col      : ", mediator_cols,     "\n", sep = "")
  cat("  outcome_col       : ", outcome_cols,      "\n", sep = "")
} else {
  cat("  mediator_col      : auto detect numeric non-_score columns\n")
  cat("  outcome_col       : auto detect numeric *_score columns\n")
}
cat("  outdir            : ", outdir,            "\n", sep = "")
cat("============================================================\n\n")

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

## ---------------- 1. 读取数据 & 基本检查 ---------------------

if (!file.exists(sample_scores_tsv)) {
  stop("[ERROR] 找不到样本级得分表: ", sample_scores_tsv)
}

cat("[STEP1] 读取样本级得分表...\n")
dat <- readr::read_tsv(sample_scores_tsv, show_col_types = FALSE)

required_cols <- c("sample_id", "group")
missing_cols  <- setdiff(required_cols, colnames(dat))
if (length(missing_cols) > 0) {
  stop("[ERROR] 数据表中缺少以下列： ",
       paste(missing_cols, collapse = ", "),
       "\n实际列名: ", paste(colnames(dat), collapse = ", "))
}

## 在 grid 模式下，自动识别 mediator_cols / outcome_cols --------------------

if (mode == "grid") {
  non_candidate <- c("sample_id", "group", "batch",
                     "rna_sample_id", "lipid_sample_id")
  num_cols <- names(dat)[vapply(dat, is.numeric, logical(1))]

  outcome_cols <- setdiff(num_cols[grepl("_score$", num_cols)], non_candidate)
  mediator_cols <- setdiff(num_cols[!grepl("_score$", num_cols)], non_candidate)

  if (length(mediator_cols) == 0) {
    stop("[ERROR] 自动模式下未识别到任何候选中介列 (M)。\n",
         "  规则：数值型 & 不以 '_score' 结尾 & 不在 ",
         paste(non_candidate, collapse = ", "), " 内。")
  }
  if (length(outcome_cols) == 0) {
    stop("[ERROR] 自动模式下未识别到任何候选结果列 (Y)。\n",
         "  规则：数值型 & 列名以 '_score' 结尾。")
  }

  cat("[INFO] 自动识别到候选 M 列 (", length(mediator_cols), ")：\n  ",
      paste(mediator_cols, collapse = ", "), "\n", sep = "")
  cat("[INFO] 自动识别到候选 Y 列 (", length(outcome_cols), ")：\n  ",
      paste(outcome_cols, collapse = ", "), "\n", sep = "")
}

## 若是 single 模式，检查指定列是否存在且为数值型 --------------------------

if (mode == "single") {
  if (!mediator_cols %in% colnames(dat)) {
    stop("[ERROR] mediator_col 在数据表中不存在: ", mediator_cols)
  }
  if (!outcome_cols %in% colnames(dat)) {
    stop("[ERROR] outcome_col 在数据表中不存在: ", outcome_cols)
  }

  if (!is.numeric(dat[[mediator_cols]])) {
    stop("[ERROR] mediator_col 不是数值型: ", mediator_cols)
  }
  if (!is.numeric(dat[[outcome_cols]])) {
    stop("[ERROR] outcome_col 不是数值型: ", outcome_cols)
  }
}

## ---------------- 2. 定义一个对单个 (M,Y) 组合做中介的函数 -----------------

run_one_mediation <- function(dat_full, mediator_col, outcome_col, B = 5000) {

  message("[PAIR] mediator_col = ", mediator_col,
          " ; outcome_col = ", outcome_col)

  ## 2.1 提取并清理数据
  if (!all(c(mediator_col, outcome_col) %in% colnames(dat_full))) {
    stop("[ERROR] 某个 (M,Y) 组合在数据表中找不到对应列: ",
         mediator_col, ", ", outcome_col)
  }

  dat_pair <- dat_full[, c("sample_id", "group", mediator_col, outcome_col)]
  colnames(dat_pair) <- c("sample_id", "group", "M", "Y")

  dat_clean <- dat_pair %>%
    dplyr::filter(!is.na(group), !is.na(M), !is.na(Y))

  n <- nrow(dat_clean)
  message("  [INFO] 有效样本数 n = ", n)

  if (n < 4) {
    warning("[WARN] 有效样本数 < 4（n = ", n, "），回归与 bootstrap 结果极不稳定，仅作探索性参考。\n")
  }

  ## 2.2 处理 group -> X (0/1)
  g_levels <- sort(unique(dat_clean$group))
  if (length(g_levels) != 2) {
    stop("[ERROR] group 列必须恰好有 2 个水平，用于构建二分类自变量 X。\n",
         "当前 group 水平: ", paste(g_levels, collapse = ", "))
  }
  ko_keywords <- c("HO","KO","knockout","EphB1_KO","EphB1KO")

  group_focal <- g_levels[which(g_levels %in% ko_keywords)]
  group_ref   <- g_levels[which(!(g_levels %in% ko_keywords))]

  if (length(group_focal) != 1 | length(group_ref) != 1) {
    stop("[ERROR] group 无法唯一识别 KO 与 WT，请检查命名。\n",
         "  当前水平: ", paste(g_levels, collapse = ", "))
  }

  message("  [INFO] group 映射：")
  message("    X = 0 -> ", group_ref,   " (reference)")
  message("    X = 1 -> ", group_focal, " (focal)")

  dat_clean <- dat_clean %>%
    mutate(X = ifelse(group == group_focal, 1, 0))

  ## 2.3 拟合三个回归模型：a, c, b/c'
  message("  [STEP] 拟合 a, b, c, c' 路径线性模型...")

  ## 模型 1: M ~ X  (a 路径)
  fit_a <- lm(M ~ X, data = dat_clean)
  sum_a <- summary(fit_a)
  coef_a <- coef(sum_a)

  a     <- unname(coef_a["X", "Estimate"])
  se_a  <- unname(coef_a["X", "Std. Error"])
  p_a   <- unname(coef_a["X", "Pr(>|t|)"])

  ## 模型 2: Y ~ X  (c 路径，总效应)
  fit_c <- lm(Y ~ X, data = dat_clean)
  sum_c <- summary(fit_c)
  coef_c <- coef(sum_c)

  c_tot   <- unname(coef_c["X", "Estimate"])
  se_c    <- unname(coef_c["X", "Std. Error"])
  p_c     <- unname(coef_c["X", "Pr(>|t|)"])

  ## 模型 3: Y ~ X + M  (b 路径 & c' 路径)
  fit_b <- lm(Y ~ X + M, data = dat_clean)
  sum_b <- summary(fit_b)
  coef_b <- coef(sum_b)

  b        <- unname(coef_b["M", "Estimate"])
  se_b     <- unname(coef_b["M", "Std. Error"])
  p_b      <- unname(coef_b["M", "Pr(>|t|)"])

  c_prime  <- unname(coef_b["X", "Estimate"])
  se_c_p   <- unname(coef_b["X", "Std. Error"])
  p_c_p    <- unname(coef_b["X", "Pr(>|t|)"])

  message("  [INFO] a = ", round(a, 4),
          " ; b = ", round(b, 4),
          " ; c = ", round(c_tot, 4),
          " ; c' = ", round(c_prime, 4))

  ## 2.4 计算中介效应 a*b + Sobel
  message("  [STEP] 计算中介效应 a*b 及 Sobel 检验...")

  ab <- a * b

  se_ab_sobel <- sqrt(b^2 * se_a^2 + a^2 * se_b^2)
  if (se_ab_sobel > 0) {
    z_sobel <- ab / se_ab_sobel
    p_sobel <- 2 * (1 - pnorm(abs(z_sobel)))
  } else {
    z_sobel <- NA_real_
    p_sobel <- NA_real_
    warning("[WARN] Sobel 标准误为 0，无法计算 z 值和 p 值。")
  }

  message("  [INFO] 间接效应 a*b = ", round(ab, 4),
          " ; Sobel z = ",
          ifelse(is.na(z_sobel), "NA", round(z_sobel, 3)),
          " ; p = ",
          ifelse(is.na(p_sobel), "NA", signif(p_sobel, 3)))

  ## 2.5 Bootstrap 中介效应
  message("  [STEP] Bootstrap 间接效应 (a*b)...")

  set.seed(12345)  # 使结果可复现
  boot_ab <- numeric(B)
  boot_a  <- numeric(B)
  boot_b  <- numeric(B)

  for (b_i in seq_len(B)) {
    idx <- sample.int(n, size = n, replace = TRUE)
    dat_b <- dat_clean[idx, ]

    ok <- TRUE
    a_bi <- NA_real_
    b_bi <- NA_real_

    ## a 路径
    fit_a_b <- try(lm(M ~ X, data = dat_b), silent = TRUE)
    if (inherits(fit_a_b, "try-error")) {
      ok <- FALSE
    } else {
      coef_a_b <- coef(summary(fit_a_b))
      if ("X" %in% rownames(coef_a_b)) {
        a_bi <- coef_a_b["X", "Estimate"]
      } else {
        ok <- FALSE
      }
    }

    ## b 路径
    if (ok) {
      fit_b_b <- try(lm(Y ~ X + M, data = dat_b), silent = TRUE)
      if (inherits(fit_b_b, "try-error")) {
        ok <- FALSE
      } else {
        coef_b_b <- coef(summary(fit_b_b))
        if ("M" %in% rownames(coef_b_b)) {
          b_bi <- coef_b_b["M", "Estimate"]
        } else {
          ok <- FALSE
        }
      }
    }

    if (!ok || is.na(a_bi) || is.na(b_bi)) {
      boot_a[b_i]  <- NA_real_
      boot_b[b_i]  <- NA_real_
      boot_ab[b_i] <- NA_real_
    } else {
      boot_a[b_i]  <- a_bi
      boot_b[b_i]  <- b_bi
      boot_ab[b_i] <- a_bi * b_bi
    }
  }

  boot_ab_valid <- boot_ab[!is.na(boot_ab)]
  boot_a_valid  <- boot_a[!is.na(boot_a)]
  boot_b_valid  <- boot_b[!is.na(boot_b)]

  if (length(boot_ab_valid) < length(boot_ab) * 0.5) {
    warning("[WARN] 有大量 bootstrap 样本无法拟合模型，bootstrap 结果可能极不稳定。\n")
  }

  ab_boot_mean <- mean(boot_ab_valid, na.rm = TRUE)
  ab_boot_ci   <- quantile(boot_ab_valid, probs = c(0.025, 0.975), na.rm = TRUE)

  a_boot_mean <- mean(boot_a_valid, na.rm = TRUE)
  a_boot_ci   <- quantile(boot_a_valid, probs = c(0.025, 0.975), na.rm = TRUE)

  b_boot_mean <- mean(boot_b_valid, na.rm = TRUE)
  b_boot_ci   <- quantile(boot_b_valid, probs = c(0.025, 0.975), na.rm = TRUE)

  message("  [INFO] Bootstrap 间接效应均值 = ", round(ab_boot_mean, 4),
          " ; 95% CI = [", round(ab_boot_ci[1], 4), ", ",
          round(ab_boot_ci[2], 4), "]")
  message("  [INFO] Bootstrap a 均值 = ", round(a_boot_mean, 4),
          " ; 95% CI = [", round(a_boot_ci[1], 4), ", ",
          round(a_boot_ci[2], 4), "]")
  message("  [INFO] Bootstrap b 均值 = ", round(b_boot_mean, 4),
          " ; 95% CI = [", round(b_boot_ci[1], 4), ", ",
          round(b_boot_ci[2], 4), "]")

  summary_row <- tibble::tibble(
    n               = n,
    group_ref       = group_ref,
    group_focal     = group_focal,
    mediator_col    = mediator_col,
    outcome_col     = outcome_col,

    a               = a,
    se_a            = se_a,
    p_a             = p_a,

    b               = b,
    se_b            = se_b,
    p_b             = p_b,

    c_total         = c_tot,
    se_c_total      = se_c,
    p_c_total       = p_c,

    c_prime         = c_prime,
    se_c_prime      = se_c_p,
    p_c_prime       = p_c_p,

    ab              = ab,
    se_ab_sobel     = se_ab_sobel,
    z_sobel         = z_sobel,
    p_sobel         = p_sobel,

    ab_boot_mean    = ab_boot_mean,
    ab_boot_ci_low  = unname(ab_boot_ci[1]),
    ab_boot_ci_high = unname(ab_boot_ci[2]),

    a_boot_mean     = a_boot_mean,
    a_boot_ci_low   = unname(a_boot_ci[1]),
    a_boot_ci_high  = unname(a_boot_ci[2]),

    b_boot_mean     = b_boot_mean,
    b_boot_ci_low   = unname(b_boot_ci[1]),
    b_boot_ci_high  = unname(b_boot_ci[2])
  )

  list(
    summary_row = summary_row,
    boot_a      = boot_a,
    boot_b      = boot_b,
    boot_ab     = boot_ab
  )
}

## ---------------- 3. 对所有 (M,Y) 组合运行中介分析 -------------------------

cat("\n[STEP2] 对所有 (M,Y) 组合运行中介分析...\n")

all_summary <- list()
all_boot    <- list()
pair_idx    <- 1L

B_global <- 5000

for (m_col in mediator_cols) {
  for (y_col in outcome_cols) {

    cat("\n------------------------------------------------------------\n")
    cat("[PAIR] (", pair_idx, ") mediator_col = ", m_col,
        " ; outcome_col = ", y_col, "\n", sep = "")
    cat("------------------------------------------------------------\n")

    res <- try(
      run_one_mediation(dat_full = dat,
                        mediator_col = m_col,
                        outcome_col  = y_col,
                        B = B_global),
      silent = TRUE
    )

    if (inherits(res, "try-error")) {
      warning("[WARN] (", m_col, " -> ", y_col,
              ") 中介分析失败，将以 NA 占位。\n")
      tmp_summary <- tibble::tibble(
        n               = NA_integer_,
        group_ref       = NA_character_,
        group_focal     = NA_character_,
        mediator_col    = m_col,
        outcome_col     = y_col,

        a               = NA_real_,
        se_a            = NA_real_,
        p_a             = NA_real_,

        b               = NA_real_,
        se_b            = NA_real_,
        p_b             = NA_real_,

        c_total         = NA_real_,
        se_c_total      = NA_real_,
        p_c_total       = NA_real_,

        c_prime         = NA_real_,
        se_c_prime      = NA_real_,
        p_c_prime       = NA_real_,

        ab              = NA_real_,
        se_ab_sobel     = NA_real_,
        z_sobel         = NA_real_,
        p_sobel         = NA_real_,

        ab_boot_mean    = NA_real_,
        ab_boot_ci_low  = NA_real_,
        ab_boot_ci_high = NA_real_,

        a_boot_mean     = NA_real_,
        a_boot_ci_low   = NA_real_,
        a_boot_ci_high  = NA_real_,

        b_boot_mean     = NA_real_,
        b_boot_ci_low   = NA_real_,
        b_boot_ci_high  = NA_real_
      )

      tmp_boot <- tibble::tibble(
        mediator_col = m_col,
        outcome_col  = y_col,
        iter         = integer(0),
        a            = numeric(0),
        b            = numeric(0),
        ab           = numeric(0)
      )
    } else {
      tmp_summary <- res$summary_row
      tmp_boot <- tibble::tibble(
        mediator_col = m_col,
        outcome_col  = y_col,
        iter         = seq_along(res$boot_ab),
        a            = res$boot_a,
        b            = res$boot_b,
        ab           = res$boot_ab
      )
    }

    all_summary[[pair_idx]] <- tmp_summary
    all_boot[[pair_idx]]    <- tmp_boot
    pair_idx <- pair_idx + 1L
  }
}

summary_tbl <- dplyr::bind_rows(all_summary)
boot_tbl    <- dplyr::bind_rows(all_boot)

## ---------------- 4. 写出结果 -------------------------------

cat("\n[STEP3] 写出结果表...\n")

out_summary <- file.path(outdir, "mediation_summary.tsv")
out_boot    <- file.path(outdir, "mediation_bootstrap_ab.tsv")

readr::write_tsv(summary_tbl, out_summary)
readr::write_tsv(boot_tbl,    out_boot)

cat("  [OK] 写出中介分析汇总表: ", out_summary, "\n", sep = "")
cat("  [OK] 写出 bootstrap 间接效应分布: ", out_boot, "\n", sep = "")
cat("============================================================\n")
cat("[DONE] 中介效应分析完成。\n")
cat("  提示：n 较小时（例如 5），p 值和 CI 仅供方向性参考，\n")
cat("        请结合生物学合理性和多组学证据综合解读。\n")
cat("============================================================\n")