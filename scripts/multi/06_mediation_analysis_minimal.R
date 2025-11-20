#!/usr/bin/env Rscript

## ============================================================
## 08_mediation_analysis_minimal.R
##
## 目的：
##   对 X -> M -> Y 进行一个极简的线性中介分析：
##     - X: group（二分类，如 WT / HO 或 WT / KO）
##     - M: 中介变量（例如 M_score_RNA）
##     - Y: 结果变量（例如 Y_score_lipid 或 Y_score_RNA）
##
## 输入（命令行参数）：
##   1) sample_scores_tsv : 样本级得分表
##        例如: results/multiomics/tables/sample_scores/sample_scores_for_mediation.tsv
##        必须包含：
##          - sample_id
##          - group    （二分类）
##          - 以及指定的 M 列、Y 列
##   2) mediator_col      : M 列名，例如 "M_score_RNA"
##   3) outcome_col       : Y 列名，例如 "Y_score_lipid" 或 "Y_score_RNA"
##   4) outdir            : 输出目录，例如 "results/multiomics/tables/mediation"
##
## 输出：
##   - <outdir>/mediation_summary.tsv
##       包含：
##         * n
##         * group_ref, group_focal
##         * a, se_a, p_a
##         * b, se_b, p_b
##         * c, se_c, p_c        （总效应）
##         * c_prime, se_c_prime, p_c_prime  （直接效应）
##         * ab, se_ab_sobel, z_sobel, p_sobel
##         * ab_boot_mean, ab_boot_ci_low, ab_boot_ci_high
##   - <outdir>/mediation_bootstrap_ab.tsv
##       每次 bootstrap 的间接效应 a*b
##
## 注意：
##   - 使用简单线性回归 + Sobel + bootstrap（默认 B=5000）
##   - n 较小时（例如 5），p 值和 CI 仅供参考，解释时要非常保守。
## ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
  stop(
    "[ERROR] 参数不足。\n",
    "用法示例：\n",
    "  Rscript scripts/multi/08_mediation_analysis_minimal.R \\\n",
    "    results/multiomics/tables/sample_scores/sample_scores_for_mediation.tsv \\\n",
    "    M_score_RNA \\\n",
    "    Y_score_lipid \\\n",
    "    results/multiomics/tables/mediation\n",
    call. = FALSE
  )
}

sample_scores_tsv <- args[1]
mediator_col      <- args[2]
outcome_col       <- args[3]
outdir            <- args[4]

cat("============================================================\n")
cat("[INFO] 08_mediation_analysis_minimal.R\n")
cat("  sample_scores_tsv : ", sample_scores_tsv, "\n", sep = "")
cat("  mediator_col      : ", mediator_col,      "\n", sep = "")
cat("  outcome_col       : ", outcome_col,       "\n", sep = "")
cat("  outdir            : ", outdir,            "\n", sep = "")
cat("============================================================\n\n")

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

## ---------------- 1. 读取数据 & 基本检查 ---------------------

if (!file.exists(sample_scores_tsv)) {
  stop("[ERROR] 找不到样本级得分表: ", sample_scores_tsv)
}

cat("[STEP1] 读取样本级得分表...\n")
dat <- readr::read_tsv(sample_scores_tsv, show_col_types = FALSE)

required_cols <- c("sample_id", "group", mediator_col, outcome_col)
missing_cols  <- setdiff(required_cols, colnames(dat))
if (length(missing_cols) > 0) {
  stop("[ERROR] 数据表中缺少以下列： ",
       paste(missing_cols, collapse = ", "),
       "\n实际列名: ", paste(colnames(dat), collapse = ", "))
}

## 去掉 NA
dat_clean <- dat %>%
  dplyr::select(sample_id, group,
                M = dplyr::all_of(mediator_col),
                Y = dplyr::all_of(outcome_col)) %>%
  dplyr::filter(!is.na(group), !is.na(M), !is.na(Y))

cat("  [INFO] 有效样本数 n = ", nrow(dat_clean), "\n", sep = "")

if (nrow(dat_clean) < 4) {
  warning("[WARN] 有效样本数 < 4，回归与 bootstrap 结果极不稳定，仅作探索性参考。\n")
}

## ---------------- 2. group -> X（二分类编码） -----------------

cat("\n[STEP2] 处理 group 列，映射为 X (0/1)...\n")

g_levels <- sort(unique(dat_clean$group))
if (length(g_levels) != 2) {
  stop("[ERROR] group 列必须恰好有 2 个水平，用于构建二分类自变量 X。\n",
       "当前 group 水平: ", paste(g_levels, collapse = ", "))
}
ko_keywords <- c("HO","KO","knockout","EphB1_KO","EphB1KO")

group_focal <- g_levels[ which(g_levels %in% ko_keywords) ]
group_ref   <- g_levels[ which(!(g_levels %in% ko_keywords)) ]

if (length(group_focal) != 1 | length(group_ref) != 1) {
    stop("[ERROR] group 无法唯一识别 KO 与 WT，请检查命名。")
}
cat("  [INFO] group 映射：\n")
cat("    X = 0 -> ", group_ref,   " (reference)\n", sep = "")
cat("    X = 1 -> ", group_focal, " (focal)\n", sep = "")

dat_clean <- dat_clean %>%
  mutate(
    X = ifelse(group == group_focal, 1, 0)
  )

## ---------------- 3. 拟合三个回归模型 ------------------------

cat("\n[STEP3] 拟合 a, b, c, c' 路径的线性模型...\n")

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

cat("  [INFO] a = ", round(a, 4), " ; b = ", round(b, 4),
    " ; c = ", round(c_tot, 4), " ; c' = ", round(c_prime, 4), "\n", sep = "")

## ---------------- 4. 中介效应 a*b + Sobel --------------------

cat("\n[STEP4] 计算中介效应 a*b 及 Sobel 检验...\n")

ab <- a * b

## Sobel 标准误
se_ab_sobel <- sqrt(b^2 * se_a^2 + a^2 * se_b^2)
if (se_ab_sobel > 0) {
  z_sobel <- ab / se_ab_sobel
  p_sobel <- 2 * (1 - pnorm(abs(z_sobel)))
} else {
  z_sobel <- NA_real_
  p_sobel <- NA_real_
  warning("[WARN] Sobel 标准误为 0，无法计算 z 值和 p 值。")
}

cat("  [INFO] 间接效应 a*b = ", round(ab, 4),
    " ; Sobel z = ", ifelse(is.na(z_sobel), "NA", round(z_sobel, 3)),
    " ; p = ", ifelse(is.na(p_sobel), "NA", signif(p_sobel, 3)), "\n", sep = "")

## ---------------- 5. Bootstrap 中介效应 ----------------------

cat("\n[STEP5] Bootstrap 间接效应 (a*b)...\n")

set.seed(12345)  # 使结果可复现
B <- 5000

boot_ab <- numeric(B)

n <- nrow(dat_clean)

for (b_i in seq_len(B)) {
  idx <- sample.int(n, size = n, replace = TRUE)
  dat_b <- dat_clean[idx, ]

  ## 尝试拟合，防止奇异情况报错
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
    boot_ab[b_i] <- NA_real_
  } else {
    boot_ab[b_i] <- a_bi * b_bi
  }
}

boot_ab_valid <- boot_ab[!is.na(boot_ab)]
if (length(boot_ab_valid) < B * 0.5) {
  warning("[WARN] 有大量 bootstrap 样本无法拟合模型，bootstrap 结果可能极不稳定。")
}

ab_boot_mean <- mean(boot_ab_valid, na.rm = TRUE)
ab_boot_ci   <- quantile(boot_ab_valid, probs = c(0.025, 0.975), na.rm = TRUE)

cat("  [INFO] Bootstrap 间接效应均值 = ", round(ab_boot_mean, 4),
    " ; 95% CI = [", round(ab_boot_ci[1], 4), ", ",
    round(ab_boot_ci[2], 4), "]\n", sep = "")

## 写出 bootstrap 抽样结果
boot_path <- file.path(outdir, "mediation_bootstrap_ab.tsv")
readr::write_tsv(
  tibble::tibble(iter = seq_along(boot_ab), ab = boot_ab),
  boot_path
)
cat("  [OK] 写出 bootstrap 间接效应分布: ", boot_path, "\n", sep = "")

## ---------------- 6. 汇总结果写出 ---------------------------

cat("\n[STEP6] 写出汇总表...\n")

summary_tbl <- tibble::tibble(
  n               = nrow(dat_clean),
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
  ab_boot_ci_high = unname(ab_boot_ci[2])
)

out_summary <- file.path(outdir, "mediation_summary.tsv")
readr::write_tsv(summary_tbl, out_summary)

cat("  [OK] 写出中介分析汇总表: ", out_summary, "\n", sep = "")
cat("============================================================\n")
cat("[DONE] 中介效应分析完成。\n")
cat("  提示：n 较小时（例如 5），p 值和 CI 仅供方向性参考，请结合生物学合理性解读。\n")
cat("============================================================\n")