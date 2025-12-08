#!/usr/bin/env Rscript

## figure_7_b_pdf_general.R
## 在 Nature 风格配置约束下，将 figure_7_b 表格导出为 PDF（自动换行 + 自动尺寸）

suppressPackageStartupMessages({
  library(yaml)
  library(readr)
  library(dplyr)
  library(gridExtra)
  library(grid)
  library(stringr)
})

## -------------------------------
## 1. 解析脚本路径 & 项目根目录
## -------------------------------

get_script_path <- function() {
  # 在 Rscript 下从 --file= 参数获取脚本路径
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 0) {
    # 如果在交互环境中 source()，则退回当前工作目录
    return(normalizePath(getwd()))
  }
  normalizePath(sub("^--file=", "", file_arg[1]))
}

script_path  <- get_script_path()
script_dir   <- dirname(script_path)
project_root <- normalizePath(file.path(script_dir, "..", ".."))

## -------------------------------
## 2. 路径设置
## -------------------------------

config_path <- "/Users/gongbaoming/Library/CloudStorage/OneDrive-个人/EphB1/02_protocols/figure_style_nature.yaml"

input_tsv  <- file.path(project_root, "results/figs/figure_7_b_ephb1_downstream_signaling_table.tsv")
output_pdf <- file.path(project_root, "results/figs/figure_7_b.pdf")

if (!file.exists(config_path)) {
  stop("配置文件未找到：", config_path)
}
if (!file.exists(input_tsv)) {
  stop("输入 TSV 文件未找到：", input_tsv)
}

dir.create(dirname(output_pdf), recursive = TRUE, showWarnings = FALSE)

## -------------------------------
## 3. 读取配置（Nature 风格）
## -------------------------------

cfg <- yaml::read_yaml(config_path)

font_family   <- cfg$typography$font_family_primary
size_body_pt  <- cfg$typography$sizes_pt$table_body_default
size_head_pt  <- cfg$typography$sizes_pt$table_header_default

mm_per_inch   <- cfg$units$mm_per_inch
max_width_mm  <- cfg$page_layout$max_width_mm

axis_line_pt  <- cfg$lines$axis_line_default_pt
pt_per_lwd    <- cfg$lines$r_lwd_scale$pt_per_lwd
border_lwd    <- axis_line_pt / pt_per_lwd

## -------------------------------
## 4. 读取数据，只保留 3 列 + 文本换行
## -------------------------------

df_raw <- readr::read_tsv(input_tsv, col_types = cols())

df_table <- df_raw %>%
  dplyr::select(label, n_genes)

# 长文本列换行（经验值，可按需要调整）
wrap_label_chars           <- 32



## -------------------------------
## 5. 构建 tableGrob 并应用样式
## -------------------------------

tt <- gridExtra::ttheme_minimal(
  core = list(
    fg_params = list(
      fontfamily = font_family,
      fontsize   = size_body_pt,
      lineheight = 1.1
    ),
    bg_params = list(
      lwd = border_lwd
    ),
    padding = unit(c(1.2, 1.2), "mm")
  ),
  colhead = list(
    fg_params = list(
      fontfamily = font_family,
      fontsize   = size_head_pt,
      fontface   = "bold",
      lineheight = 1.1
    ),
    bg_params = list(
      lwd = border_lwd
    ),
    padding = unit(c(1.2, 1.2), "mm")
  )
)

# 让表头和表体文字左对齐，避免长行向左溢出
tt$core$fg_params$hjust    <- 0
tt$core$fg_params$x        <- 0
tt$colhead$fg_params$hjust <- 0
tt$colhead$fg_params$x     <- 0

# 使用自动换行后的 df_table 构建表格
tbl <- gridExtra::tableGrob(
  df_table,
  rows  = NULL,  # 不显示左侧行名
  theme = tt
)

## 列宽使用相对宽度，保证整表在设备宽度内
# 2 列：label / n_genes（第一列更宽，留出更多空间放文字）
tbl$widths <- unit(c(0.88, 0.12), "null")

# 统一所有 cell 的边框线宽
for (i in seq_along(tbl$grobs)) {
  g <- tbl$grobs[[i]]
  if (inherits(g, "grob") && !is.null(g$gp)) {
    g$gp$lwd <- border_lwd
    tbl$grobs[[i]] <- g
  }
}

## -------------------------------
## 6. 输出 PDF（用户自定义尺寸）
## -------------------------------

# 用户自定义图像尺寸（单位：mm），可根据需要修改
width_mm  <- 95
height_mm <- 35

width_in  <- width_mm  / mm_per_inch
height_in <- height_mm / mm_per_inch

pdf(file = output_pdf,
    width  = width_in,
    height = height_in,
    family = font_family,
    useDingbats = FALSE)  # 避免 Illustrator 中字体替换问题

grid::grid.newpage()
# 使表格在页面中居中
tbl$vp <- grid::viewport(x = 0.5, y = 0.5, just = c("center", "center"))
grid::grid.draw(tbl)
dev.off()

message(
  sprintf(
    "figure_7_b 已输出: %s (%.2f x %.2f mm, %.2f x %.2f in)",
    output_pdf,
    width_mm, height_mm,
    width_in, height_in
  )
)