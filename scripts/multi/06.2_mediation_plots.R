#!/usr/bin/env Rscript

## ============================================================
## 06.2_mediation_plots.R
##
## 目的：
##   可视化两个中介模型的结果：
##     1) RNA 内部：    X -> M_score_RNA -> Y_score_RNA
##     2) RNA -> 脂质： X -> M_score_RNA -> Y_score_lipid
##
## 输入（命令行参数）：
##   1) mediation_dir1 : 第一个中介模型结果目录
##        例如: results/multiomics/tables/mediation_rna_only
##   2) mediation_dir2 : 第二个中介模型结果目录
##        例如: results/multiomics/tables/mediation_rna_to_lipid
##   3) outdir         : 输出图像目录
##        例如: results/multiomics/plots/mediation
##
## 要求每个目录中包含：
##   - mediation_summary.tsv      （由 06_mediation_analysis_minimal.R 生成）
##   - mediation_bootstrap_ab.tsv （同上）
##
## 输出：
##   - <outdir>/mediation_effects_comparison.pdf
##       比较两个模型的 a / b / c / c' / 间接效应（ab）
##
##   - <outdir>/mediation_indirect_bootstrap.pdf
##       两个模型间接效应 a*b 的 bootstrap 分布对比
##
## ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(ggplot2)
})

## ---------------- 0. 解析命令行参数 --------------------------

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop(
    "[ERROR] 参数不足。\n",
    "用法示例：\n",
    "  Rscript scripts/multi/06.2_mediation_plots.R \\\n",
    "    results/multiomics/tables/mediation_rna_only \\\n",
    "    results/multiomics/tables/mediation_rna_to_lipid \\\n",
    "    results/multiomics/plots/mediation\n",
    call. = FALSE
  )
}

med_dir1 <- args[1]
med_dir2 <- args[2]
outdir   <- args[3]

cat("============================================================\n")
cat("[INFO] 06.2_mediation_plots.R\n")
cat("  mediation_dir1 : ", med_dir1, "\n", sep = "")
cat("  mediation_dir2 : ", med_dir2, "\n", sep = "")
cat("  outdir         : ", outdir,   "\n", sep = "")
cat("============================================================\n\n")

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

## 辅助函数：读取某个中介目录的 summary & bootstrap
read_mediation_dir <- function(med_dir, label = NULL) {
  summary_path <- file.path(med_dir, "mediation_summary.tsv")
  boot_path    <- file.path(med_dir, "mediation_bootstrap_ab.tsv")

  if (!file.exists(summary_path)) {
    stop("[ERROR] 在目录中找不到 mediation_summary.tsv: ", summary_path)
  }
  if (!file.exists(boot_path)) {
    stop("[ERROR] 在目录中找不到 mediation_bootstrap_ab.tsv: ", boot_path)
  }

  summary_tbl <- readr::read_tsv(summary_path, show_col_types = FALSE)
  boot_tbl    <- readr::read_tsv(boot_path, show_col_types = FALSE)

  if (!"ab" %in% colnames(boot_tbl)) {
    stop("[ERROR] mediation_bootstrap_ab.tsv 中必须包含列 'ab': ", boot_path)
  }

  if (!"outcome_col" %in% colnames(summary_tbl)) {
    stop("[ERROR] mediation_summary.tsv 中必须包含列 'outcome_col': ", summary_path)
  }

  if (is.null(label)) {
    label <- unique(summary_tbl$outcome_col)
    if (length(label) != 1) {
      label <- paste0(basename(med_dir))
    }
  }

  summary_tbl$Model <- label
  boot_tbl$Model    <- label

  list(summary = summary_tbl, boot = boot_tbl)
}

## ---------------- 1. 读取两个模型的结果 -----------------------

cat("[STEP1] 读取两个中介模型的结果...\n")

med1 <- read_mediation_dir(med_dir1, label = NULL)
med2 <- read_mediation_dir(med_dir2, label = NULL)

sum1 <- med1$summary
sum2 <- med2$summary
boot1 <- med1$boot
boot2 <- med2$boot

all_summary <- dplyr::bind_rows(sum1, sum2)
all_boot    <- dplyr::bind_rows(boot1, boot2)

## 根据 outcome_col 生成更易读的模型标签
## 例如：Y_score_RNA -> "RNA Y (M_score_RNA → Y_score_RNA)"
##      Y_score_lipid -> "Lipid Y (M_score_RNA → Y_score_lipid)"
model_labels <- all_summary %>%
  dplyr::distinct(Model, mediator_col, outcome_col) %>%
  dplyr::mutate(
    Model_pretty = dplyr::case_when(
      outcome_col == "Y_score_RNA"   ~ "RNA Y (M -> RNA axis)",
      outcome_col == "Y_score_lipid" ~ "Lipid Y (M -> lipid axis)",
      TRUE ~ paste0("M: ", mediator_col, " → Y: ", outcome_col)
    )
  )

all_summary <- all_summary %>%
  dplyr::left_join(model_labels %>% dplyr::select(Model, Model_pretty),
                   by = "Model")

all_boot <- all_boot %>%
  dplyr::left_join(model_labels %>% dplyr::select(Model, Model_pretty),
                   by = "Model")

## ---------------- 2. 组合效应尺寸表用于绘图 -------------------

cat("\n[STEP2] 组合 a / b / c / c' / 间接效应 的表...\n")

effects_long <- all_summary %>%
  dplyr::select(
    Model, Model_pretty,
    a, se_a,
    b, se_b,
    c_total, se_c_total,
    c_prime, se_c_prime,
    ab, ab_boot_ci_low, ab_boot_ci_high
  ) %>%
  tidyr::pivot_longer(
    cols = c(a, b, c_total, c_prime, ab),
    names_to = "effect",
    values_to = "estimate"
  ) %>%
  dplyr::mutate(
    se = dplyr::case_when(
      effect == "a"        ~ se_a,
      effect == "b"        ~ se_b,
      effect == "c_total"  ~ se_c_total,
      effect == "c_prime"  ~ se_c_prime,
      effect == "ab"       ~ NA_real_,   # ab 用 bootstrap CI，不用 se
      TRUE ~ NA_real_
    ),
    effect_label = dplyr::case_when(
      effect == "a"        ~ "a: X → M",
      effect == "b"        ~ "b: M → Y (控制 X)",
      effect == "c_total"  ~ "c: X → Y 总效应",
      effect == "c_prime"  ~ "c': X → Y 直接效应",
      effect == "ab"       ~ "间接效应 a×b",
      TRUE ~ effect
    )
  )

## 为误差条生成上下界：
effects_long <- effects_long %>%
  dplyr::mutate(
    ymin = dplyr::case_when(
      effect %in% c("a", "b", "c_total", "c_prime") ~ estimate - 1.96 * se,
      effect == "ab" ~ ab_boot_ci_low,
      TRUE ~ NA_real_
    ),
    ymax = dplyr::case_when(
      effect %in% c("a", "b", "c_total", "c_prime") ~ estimate + 1.96 * se,
      effect == "ab" ~ ab_boot_ci_high,
      TRUE ~ NA_real_
    )
  )

## 限制绘图顺序
effects_long$effect_label <- factor(
  effects_long$effect_label,
  levels = c("a: X → M",
             "b: M → Y (控制 X)",
             "c: X → Y 总效应",
             "c': X → Y 直接效应",
             "间接效应 a×b")
)

## ---------------- 3. 绘制效应对比图 --------------------------

cat("\n[STEP3] 绘制 a / b / c / c' / ab 对比图...\n")

p_effects <- ggplot(effects_long,
                    aes(x = effect_label,
                        y = estimate,
                        colour = Model_pretty)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(ymin = ymin, ymax = ymax),
                  position = position_dodge(width = 0.6)) +
  labs(
    x = "效应类型",
    y = "效应大小 (估计值 ± 95% CI)",
    colour = "模型",
    title  = "两个中介模型的效应对比",
    subtitle = "a / b / c / c' 为正态近似 95% CI，间接效应 a×b 为 bootstrap 95% CI"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 20, hjust = 1)
  )

out_effects <- file.path(outdir, "mediation_effects_comparison.pdf")
ggsave(out_effects, p_effects,
       width = 8, height = 5)
cat("  [OK] 写出效应对比图: ", out_effects, "\n", sep = "")

## ---------------- 4. 绘制间接效应 bootstrap 分布 --------------

cat("\n[STEP4] 绘制间接效应 a×b 的 bootstrap 分布...\n")

p_boot <- ggplot(all_boot, aes(x = ab, fill = Model_pretty, colour = Model_pretty)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_density(alpha = 0.3, adjust = 1) +
  labs(
    x = "间接效应 a×b (bootstrap 样本)",
    y = "密度",
    fill   = "模型",
    colour = "模型",
    title  = "间接效应 a×b 的 bootstrap 分布",
    subtitle = "虚线为 0；若分布整体偏离 0，则间接效应方向较稳定"
  ) +
  theme_bw(base_size = 12)

out_boot <- file.path(outdir, "mediation_indirect_bootstrap.pdf")
ggsave(out_boot, p_boot,
       width = 7, height = 5)
cat("  [OK] 写出间接效应 bootstrap 分布图: ", out_boot, "\n", sep = "")

cat("============================================================\n")
cat("[DONE] 中介效应可视化完成。\n")
cat("  - 效应对比: ", out_effects, "\n", sep = "")
cat("  - 间接效应分布: ", out_boot,   "\n", sep = "")
cat("============================================================\n")