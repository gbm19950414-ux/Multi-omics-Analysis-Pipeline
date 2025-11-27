#!/usr/bin/env Rscript

## ============================================================
## figure_6_e.R
##
## Figure 6E – 跨组学“综合机制得分”的条形图
##
## 目标：
##   - 把每条机制轴在 RNA + Lipid 的综合 multi-omics Z
##     (Z_multiomics_overall) 画成一维条形图。
##   - 使用 p_multiomics_overall 添加显著性星号，
##     高亮“主角”机制轴。
##
## 依赖：
##   - scripts/multi/01_multiomics_axis_Z_xxx.R 的输出：
##       results/multiomics/tables/multiomics_axis_Z_overall.tsv
##
## 默认用法：
##   Rscript scripts/plot/figure_6_e.R \
##     scripts/plot/figure_6_e_config.yaml
##
## 若不提供配置文件，则使用脚本内置默认路径。
##
## ============================================================

suppressPackageStartupMessages({
  library(yaml)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(stringr)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

## -------- 0. 解析命令行参数与配置 ----------------------------

args <- commandArgs(trailingOnly = TRUE)

config_path <- if (length(args) >= 1) {
  args[1]
} else {
  "scripts/plot/figure_6_e_config.yaml"
}

if (file.exists(config_path)) {
  cat("[INFO] 使用配置文件: ", config_path, "\n", sep = "")
  cfg <- yaml::read_yaml(config_path)
} else {
  cat("[WARN] 未找到配置文件: ", config_path,
      " ，将使用脚本内置默认路径。\n", sep = "")
  cfg <- list()
}

tables_dir       <- cfg$tables_dir       %||% "results/multiomics/tables"
overall_file     <- cfg$overall_file     %||% "multiomics_axis_Z_overall.tsv"
plots_dir        <- cfg$plots_dir        %||% "results/figs"
output_file      <- cfg$output_file      %||% "figure_6_e_multiomics_axis_barplot.pdf"

## 主角轴高亮标准：默认按 p_multiomics_overall ≤ 0.05
main_axis_p_cutoff <- cfg$main_axis_p_cutoff %||% 0.05

## 轴顺序：如果配置文件中给了 axis_levels，则用该顺序
axis_levels_cfg <- cfg$axis_levels %||% NULL

overall_path <- file.path(tables_dir, overall_file)
out_plot_path <- file.path(plots_dir, output_file)

cat("============================================================\n")
cat("[INFO] figure_6_e.R\n")
cat("  input overall Z: ", overall_path,  "\n", sep = "")
cat("  output plot    : ", out_plot_path, "\n", sep = "")
cat("============================================================\n\n")

## -------- 1. 读取 multi-omics 轴汇总表 -----------------------

if (!file.exists(overall_path)) {
  stop("[ERROR] 找不到 multi-omics 轴汇总表: ", overall_path)
}

overall <- readr::read_tsv(overall_path, show_col_types = FALSE)

cat("[STEP1] 读入 multiomics_axis_Z_overall.tsv\n")
cat("  nrow = ", nrow(overall), " ; ncol = ", ncol(overall), "\n", sep = "")

## 必要列检查
required_cols <- c(
  "axis",
  "Z_multiomics_overall",
  "p_multiomics_overall"
)
miss <- setdiff(required_cols, colnames(overall))
if (length(miss) > 0) {
  stop("[ERROR] 汇总表缺少必要列: ", paste(miss, collapse = ", "))
}

## -------- 2. 整理数据：轴顺序 + 显著性标记 -------------------

cat("[STEP2] 整理轴顺序和显著性标记...\n")

## 轴顺序：优先使用配置文件中的 axis_levels，否则按文件出现顺序
if (!is.null(axis_levels_cfg)) {
  axis_levels <- axis_levels_cfg
} else {
  axis_levels <- overall$axis %>% unique()
}

dat <- overall %>%
  mutate(
    axis = factor(axis, levels = axis_levels),
    ## p 值星号
    p_multiomics_overall = as.numeric(p_multiomics_overall),
    p_label = case_when(
      is.na(p_multiomics_overall)            ~ "",
      p_multiomics_overall <= 0.001          ~ "***",
      p_multiomics_overall <= 0.01           ~ "**",
      p_multiomics_overall <= 0.05           ~ "*",
      TRUE                                   ~ ""
    ),
    ## 高亮“主角轴”：
    is_main = !is.na(p_multiomics_overall) &
              p_multiomics_overall <= main_axis_p_cutoff
  )

## -------- 3. 画图 ------------------------------------------------

cat("[STEP3] 绘制 Figure 6E 条形图...\n")

if (!dir.exists(plots_dir)) {
  dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
}

# y 轴范围稍微留出空间给星号
y_min <- min(dat$Z_multiomics_overall, na.rm = TRUE)
y_max <- max(dat$Z_multiomics_overall, na.rm = TRUE)
y_range <- y_max - y_min
y_star_offset <- 0.05 * y_range
y_star_pos <- ifelse(dat$Z_multiomics_overall >= 0,
                     dat$Z_multiomics_overall + y_star_offset,
                     dat$Z_multiomics_overall - y_star_offset)

p <- ggplot(dat, aes(x = axis, y = Z_multiomics_overall)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
  geom_col(
    aes(
      fill = is_main,
      color = is_main
    ),
    width = 0.6,
    linewidth = 0.5
  ) +
  ## 星号标记
  geom_text(
    aes(
      label = p_label,
      y = y_star_pos
    ),
    vjust = ifelse(dat$Z_multiomics_overall >= 0, 0, 1),
    size = 4,
    show.legend = FALSE
  ) +
  scale_fill_manual(
    values = c(
      `FALSE` = "grey75",   # 非主角轴：浅灰
      `TRUE`  = "grey40"    # 主角轴：深灰或可替换为带色
    )
  ) +
  scale_color_manual(
    values = c(
      `FALSE` = "grey40",
      `TRUE`  = "black"     # 主角轴：更明显边框
    )
  ) +
  labs(
    x = "Mechanistic axis",
    y = "Multi-omics mechanism score (Z)",
    title = "Figure 6E – Cross-omics mechanistic axis summary",
    subtitle = "Positive Z: mechanism enhanced in Ephb1−/−; Negative Z: mechanism suppressed"
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)
  )

ggsave(out_plot_path, p, width = 5.5, height = 3.8)

cat("[OK] 成功输出图像: ", out_plot_path, "\n", sep = "")
cat("============================================================\n")
cat("[DONE] Figure 6E 绘制完成。\n")
cat("============================================================\n")