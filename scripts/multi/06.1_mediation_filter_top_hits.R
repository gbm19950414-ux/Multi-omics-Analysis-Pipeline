#!/usr/bin/env Rscript

## ============================================================
## 07_mediation_filter_top_hits.R
##
## 目的：
##   从 06_mediation_analysis_minimal.R 生成的
##   mediation_summary.tsv 中，筛选出：
##     - M1（Ephb1 下游信号通路）里，
##       每个 CL 轴（Y_axis）“最有希望”的 (M, Y) 组合
##     - M2（CL 上游 RNA 轴通路）里，
##       每个 CL 轴（Y_axis）“最有希望”的 (M, Y) 组合
##
## 输入（默认）：
##   results/multiomics/tables/mediation/mediation_summary.tsv
##
## 输出（固定在 results/multiomics/tables/mediation）：
##   results/multiomics/tables/mediation/mediation_top_hits_M1.tsv
##   results/multiomics/tables/mediation/mediation_top_hits_M2.tsv
##
## “最有希望”的 operational 定义（当前版本）：
##   - 过滤掉 ab_boot_mean 为 NA 的组合
##   - 在同一 Y_axis 内：
##       1. 按 |ab_boot_mean| 从大到小排序
##       2. 再按 p_sobel 从小到大排序
##       3. 取每个 Y_axis 的前 1 条记录
##
##   （你想多保留几个候选，只要把 slice_head(n = 1)
##     改成 3 或 5 即可。）
##
## M1 / M2 的定义：
##   - mediator_col 以以下前缀开头的，视为 M2：
##       Synthesis_
##       Remodeling_
##       Oxidation_
##       Supply_
##       Transport_
##       Membrane_context_
##   - 其它 mediator 视为 M1（Ephb1 下游或其它）
##
## Y 轴名：
##   - 从 outcome_col 中去掉 "_score" 即为 Y_axis
##     例如 outcome_col = "Synthesis_score" -> Y_axis = "Synthesis"
##          outcome_col = "Membrane context_score" -> Y_axis = "Membrane context"
##
## ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
})

## ---------------- 0. 默认路径 & 参数 ------------------------

default_infile  <- "results/multiomics/tables/mediation/mediation_summary.tsv"
default_outdir  <- "results/multiomics/tables/mediation"

args <- commandArgs(trailingOnly = TRUE)

infile <- if (length(args) >= 1) args[1] else default_infile
outdir <- if (length(args) >= 2) args[2] else default_outdir

cat("============================================================\n")
cat("[INFO] 07_mediation_filter_top_hits.R\n")
cat("  infile  : ", infile, "\n", sep = "")
cat("  outdir  : ", outdir, "\n", sep = "")
cat("============================================================\n\n")

if (!file.exists(infile)) {
  stop("[ERROR] 找不到 mediation_summary.tsv: ", infile)
}

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

## ---------------- 1. 读入 mediation_summary.tsv ----------------

cat("[STEP1] 读取 mediation_summary.tsv...\n")
med <- readr::read_tsv(infile, show_col_types = FALSE)

## 简单检查一些关键列是否存在
needed_cols <- c(
  "mediator_col", "outcome_col",
  "ab", "ab_boot_mean", "ab_boot_ci_low", "ab_boot_ci_high",
  "z_sobel", "p_sobel",
  "c_total", "p_c_total"
)

missing <- setdiff(needed_cols, colnames(med))
if (length(missing) > 0) {
  stop("[ERROR] mediation_summary.tsv 中缺少以下列： ",
       paste(missing, collapse = ", "))
}

## ---------------- 2. 标注 M1 / M2 & 轴信息 --------------------

cat("[STEP2] 标注 mediator 类型 (M1/M2) 与轴信息...\n")

axis_prefixes <- c(
  "Synthesis_",
  "Remodeling_",
  "Oxidation_",
  "Supply_",
  "Transport_",
  "Membrane_context_"
)

axis_pattern <- paste0("^(", paste(axis_prefixes, collapse = "|"), ")")

med2 <- med %>%
  mutate(
    ## 是否是“轴级别上游通路”（M2）
    M_is_axis = str_detect(mediator_col, axis_pattern),

    M_type = ifelse(
      M_is_axis,
      "M2_axis_upstream",
      "M1_ephb1_downstream_or_other"
    ),

    ## 从 mediator 列名中提取 M 所属轴（只对 M2 有意义）
    M_axis = case_when(
      str_starts(mediator_col, "Synthesis_")         ~ "Synthesis",
      str_starts(mediator_col, "Remodeling_")        ~ "Remodeling",
      str_starts(mediator_col, "Oxidation_")         ~ "Oxidation",
      str_starts(mediator_col, "Supply_")            ~ "Supply",
      str_starts(mediator_col, "Transport_")         ~ "Transport",
      str_starts(mediator_col, "Membrane_context_")  ~ "Membrane context",
      TRUE                                           ~ NA_character_
    ),

    ## Y 轴名：去掉 "_score"
    Y_axis = str_remove(outcome_col, "_score$")
  )

## ---------------- 3. 定义一个通用的“每轴取最优一条”函数 ------------

pick_top_per_Y_axis <- function(df, label_M_type) {

  if (nrow(df) == 0) {
    warning("[WARN] 在 ", label_M_type, " 中没有任何记录，返回空表。\n")
    return(df)
  }

  df %>%
    filter(!is.na(ab_boot_mean)) %>%
    filter(ab_boot_mean < 0) %>%
    mutate(
      rank_metric = abs(ab_boot_mean)
    ) %>%
    group_by(Y_axis) %>%
    arrange(
      desc(rank_metric),
      p_sobel,
      .by_group = TRUE
    ) %>%
    slice_head(n = 1) %>%   # 如需每轴多条，把 1 改成 3 或 5 即可
    ungroup() %>%
    select(
      ## 轴信息
      Y_axis,
      outcome_col,
      M_type,
      M_axis,
      mediator_col,
      ## 基本样本信息
      n,
      group_ref,
      group_focal,
      ## 路径系数
      a, se_a, p_a,
      b, se_b, p_b,
      c_total, se_c_total, p_c_total,
      c_prime, se_c_prime, p_c_prime,
      ## 中介效应
      ab,
      ab_boot_mean, ab_boot_ci_low, ab_boot_ci_high,
      z_sobel, p_sobel
    ) %>%
    arrange(Y_axis)
}

## ---------------- 4. 分别对 M1 / M2 做筛选 -------------------

cat("[STEP3] 对 M1 (Ephb1 downstream) 做筛选...\n")
top_M1 <- med2 %>%
  filter(M_type == "M1_ephb1_downstream_or_other") %>%
  pick_top_per_Y_axis("M1_ephb1_downstream_or_other")

cat("[STEP4] 对 M2 (CL 上游 RNA 轴通路) 做筛选...\n")
top_M2 <- med2 %>%
  filter(M_type == "M2_axis_upstream") %>%
  pick_top_per_Y_axis("M2_axis_upstream")

## ---------------- 5. 写出精简结果表 -------------------------

out_M1 <- file.path(outdir, "mediation_top_hits_M1.tsv")
out_M2 <- file.path(outdir, "mediation_top_hits_M2.tsv")

cat("[STEP5] 写出精简结果表...\n")
readr::write_tsv(top_M1, out_M1)
readr::write_tsv(top_M2, out_M2)

cat("  [OK] M1 精简表: ", out_M1, "\n", sep = "")
cat("  [OK] M2 精简表: ", out_M2, "\n", sep = "")
cat("============================================================\n")
cat("[DONE] 筛选“每轴最有希望的中介组合”完成。\n")
cat("  若需每个轴保留多条候选，请修改函数 pick_top_per_Y_axis 中的 slice_head(n = 1)。\n")
cat("============================================================\n")