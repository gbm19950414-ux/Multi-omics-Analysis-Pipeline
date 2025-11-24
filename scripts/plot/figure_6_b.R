#!/usr/bin/env Rscript

## figure_6_b.R
## --------------
## 从 results/figs/figure_6_b.xlsx 读取数据，
## 使用 results/figs/table_config_lipid_axes.yaml 中的“外观配置”
## 生成 Word 表格（.docx）。
##
## 特点：
##   - 不规定列名，列名 & 顺序完全来自 xlsx；
##   - YAML 只控制外观：对齐方式、字体、字号、是否加粗表头等；
##   - 输出：results/figs/figure_6_b.docx（或 YAML 中 output.docx.file 指定的名字）。

suppressPackageStartupMessages({
  library(yaml)
  library(readxl)
  library(officer)
  library(flextable)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

# ---------- 根据 YAML 样式生成 flextable ----------

make_flextable <- function(df, config) {
  style_cfg <- config$style %||% list()

  # 对齐设置：style.align 是一个向量，如 ["l","c","l",...]
  align_cfg <- style_cfg$align
  n_col <- ncol(df)

  if (is.null(align_cfg)) {
    # 默认都左对齐
    aligns <- rep("left", n_col)
  } else {
    # l/c/r → left/center/right
    align_vec <- unlist(align_cfg)
    if (length(align_vec) >= n_col) {
      align_vec <- align_vec[seq_len(n_col)]
    } else {
      last_align <- align_vec[length(align_vec)]
      align_vec <- c(align_vec, rep(last_align, n_col - length(align_vec)))
    }
    aligns <- dplyr::recode(align_vec,
                            "l" = "left",
                            "c" = "center",
                            "r" = "right",
                            .default = "left")
  }

  font_family <- style_cfg$font_family %||% "Times New Roman"
  font_size   <- style_cfg$font_size   %||% 10

  header_cfg <- style_cfg$header %||% list()
  header_bold <- header_cfg$bold %||% TRUE
  header_bg   <- header_cfg$bg   %||% NA_character_

  # 构建 flextable
  ft <- flextable(df)

  # 对齐
  for (j in seq_len(n_col)) {
    ft <- align(ft, j = j, align = aligns[j], part = "all")
  }

  # 字体 & 字号
  ft <- font(ft, fontname = font_family, part = "all")
  ft <- fontsize(ft, size = font_size, part = "all")

  # 表头加粗
  if (isTRUE(header_bold)) {
    ft <- bold(ft, part = "header", bold = TRUE)
  }
  # 表头背景色（如果提供）
  if (!is.na(header_bg)) {
    ft <- bg(ft, part = "header", bg = header_bg)
  }

  # 边框：简单的 box + 内横线
  std_border <- fp_border(color = "black", width = 0.75)
  ft <- border_remove(ft)
  ft <- border_outer(ft, border = std_border, part = "all")
  ft <- border_inner_h(ft, border = std_border, part = "body")

  # 自动调整列宽
  ft <- autofit(ft)

  ft
}

# ---------- 主函数 ----------

main <- function() {
  config_path <- "results/figs/table_config_lipid_axes.yaml"
  data_path   <- "results/figs/figure_6_b.xlsx"
  output_dir  <- "results/figs"

  if (!file.exists(config_path)) {
    stop("[ERROR] 配置文件未找到: ", config_path)
  }
  if (!file.exists(data_path)) {
    stop("[ERROR] 数据文件未找到: ", data_path)
  }

  message("[INFO] 读取 YAML 配置: ", config_path)
  config <- yaml::read_yaml(config_path)

  # 读取 Excel
  sheet <- config$input$sheet %||% 1
  message("[INFO] 读取 Excel: ", data_path, " (sheet = ", sheet, ")")
  df <- readxl::read_excel(data_path, sheet = sheet)
  df <- as.data.frame(df, stringsAsFactors = FALSE)

  # 构建 flextable
  message("[INFO] 构建 flextable 表格...")
  ft <- make_flextable(df, config)

  # 输出路径
  docx_cfg   <- config$output$docx %||% list()
  out_file   <- docx_cfg$file %||% "figure_6_b.docx"
  output_path <- file.path(output_dir, out_file)

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # 新建 Word 文档并插入表格
  message("[INFO] 生成 Word 文件: ", output_path)
  doc <- read_docx()

  # 如果有 caption，就在表格前加一行文字，方便你在 Word 里手动改成“题注”
  caption <- config$caption %||% NULL
  if (!is.null(caption) && nzchar(caption)) {
    doc <- body_add_par(doc, caption, style = "Normal")
  }

  doc <- body_add_flextable(doc, ft)

  print(doc, target = output_path)

  message("[INFO] 完成 Word 表格生成。")
}

if (identical(environment(), globalenv())) {
  main()
}