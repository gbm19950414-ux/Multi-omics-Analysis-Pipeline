#!/usr/bin/env Rscript

## ============================================================
## Figure 7D: Forest plot of mediation effects (Eph downstream → Membrane context)
##
## 输入：
##   results/multiomics/tables/mediation/mediation_figure7D_membrane_axis_eph_downstream.tsv
##
## 输出：
##   results/multiomics/figs/figure_7D_mediation_membrane_axis.pdf
##
## 展示内容：
##   - ab_boot_mean（间接效应点估计）
##   - ab_boot_ci_low / ab_boot_ci_high（95% bootstrap CI）
##   - 方向标签基于 a_boot_mean 和 b_boot_mean
##
## 图类型：
##   横向森林图（forest plot）
##
## ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(stringr)
})

## -------------------- Input / Output -------------------------

infile  <- "results/multiomics/tables/mediation/mediation_figure7D_membrane_axis_eph_downstream.tsv"
outfile <- "results/figs/figure_7_d_mediation_membrane_context.pdf"

dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)

cat("[INFO] 读取输入文件: ", infile, "\n")

df <- readr::read_tsv(infile, show_col_types = FALSE)

## -------------------- Direction label 使用 boot mean -----------

df <- df %>%
  mutate(
    direction = dplyr::case_when(
      a_boot_mean > 0 & b_boot_mean > 0 ~ "a+, b+",
      a_boot_mean > 0 & b_boot_mean < 0 ~ "a+, b−",
      a_boot_mean < 0 & b_boot_mean > 0 ~ "a−, b+",
      a_boot_mean < 0 & b_boot_mean < 0 ~ "a−, b−",
      TRUE ~ ""
    ),
    # 画图时按 ab_boot_mean 绝对值排序（脚本里已排，这里确保）
    mediator_col = factor(mediator_col, levels = df$mediator_col)
  )

## -------------------- 绘图 -----------------------------------

p <- ggplot(df,
            aes(x = ab_boot_mean,
                y = mediator_col)) +

  ## 95% bootstrap CI
  geom_errorbarh(aes(xmin = ab_boot_ci_low,
                     xmax = ab_boot_ci_high),
                 height = 0.3,
                 color = "grey50",
                 size = 0.8) +

  ## 点估计
  geom_point(aes(color = ab_boot_mean > 0),
             size = 3) +

  ## 零效应虚线
  geom_vline(xintercept = 0,
             linetype = "dashed",
             color = "grey40") +

  ## 方向标签（右侧）
  geom_text(aes(label = direction,
                x = ab_boot_ci_high + 0.03 * max(abs(df$ab_boot_mean))),
            hjust = 0,
            size = 3.3) +

  scale_color_manual(
    values = c("TRUE" = "#1f78b4",   # ab>0 blue
               "FALSE" = "#e31a1c"), # ab<0 red
    guide = "none"
  ) +

  labs(
    x = "Indirect effect (a × b), bootstrap mean ± 95% CI",
    y = "EphB1 downstream pathway",
    title = "Figure 7D — Mediation of EphB1 loss → Membrane context axis",
    subtitle = "Bootstrap-based indirect effects (ab) with direction labels determined by a_boot_mean and b_boot_mean"
  ) +

  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.text.y = element_text(size = 9),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(10, 50, 10, 10)   # 右边留出给方向标签
  )

## -------------------- 保存输出 --------------------------------

ggsave(outfile, p, width = 7, height = 6)

cat("[OK] 输出图像: ", outfile, "\n")