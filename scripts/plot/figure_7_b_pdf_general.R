#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(gt)
})

# === 路径设置 ===
input_file  <- "results/figs/figure_7_b_ephb1_downstream_signaling_table.tsv"
output_file <- "results/figs/figure_7_b.pdf"

if (!file.exists(input_file)) {
  stop("找不到输入文件: ", input_file)
}

# === 读入 TSV ===
tab <- read_tsv(input_file, col_types = cols())

# 全部转成字符，避免某些列被科学计数法或其他格式化
tab_chr <- tab %>%
  mutate(across(everything(), as.character))

# === 用 gt 生成表格 ===
gt_tbl <- tab_chr %>%
  gt() %>%
  # 字体尽量小一点，行距也压缩一点，方便放下长文字
  tab_options(
    table.font.size      = 7,
    data_row.padding     = px(2),
    table.width          = pct(100),
    column_labels.font.weight = "bold"
  ) %>%
  cols_align(
    align = "left",
    columns = everything()
  )

# === 导出为 PDF ===
# 关键：vwidth / vheight 设得足够大，避免被裁掉
# expand 在周围多留一点边框
gtsave(
  gt_tbl,
  filename = output_file,
  expand   = 5,
  vwidth   = 2200,  # 水平方向视窗像素
  vheight  = 800    # 垂直方向视窗像素
)

cat("Input: ", normalizePath(input_file), "\n")
cat("Output:", normalizePath(output_file), "\n")