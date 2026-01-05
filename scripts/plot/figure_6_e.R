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
  library(grid)
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

# -------- Global figure style (Nature) --------
style_path <- "/Volumes/Samsung_SSD_990_PRO_2TB_Media/EphB1/02_protocols/figure_style_nature.yaml"
if (!file.exists(style_path)) {
  stop("[ERROR] 找不到 figure style 文件: ", style_path)
}
style <- yaml::read_yaml(style_path)

# 线宽（pt）
line_width_pt <- style$lines$line_width_pt %||% 0.5
axis_line_pt  <- style$lines$axis_line_default_pt %||% 0.25
pt_per_lwd     <- style$lines$r_lwd_scale$pt_per_lwd %||% 0.75

# ggplot2 的 linewidth 在不同元素/设备上会有“发丝线”观感差异；
# 这里采用 box_panel_from_yaml.R 的一致做法：用 pt_per_lwd 将 pt 映射为 lwd/linewidth 标度
axis_line_lwd <- axis_line_pt / pt_per_lwd
bar_line_lwd  <- line_width_pt / pt_per_lwd

# 单位换算
mm_per_inch <- style$units$mm_per_inch %||% 25.4
mm_to_pt <- function(mm) (72 / mm_per_inch) * mm

# ---- 布局：沿用 box_panel_from_yaml.R 的 axis_outer_frac + axis_gap_pt 逻辑 ----
style_layout <- style$layout %||% list()
axis_outer_frac <- style_layout$axis_outer_frac %||% list(
  left   = 0.01,  # 收紧左侧留白
  right  = 0.04,
  bottom = 0.01,  # 明显减少下方留白
  top    = 0.06
)
axis_gap_pt <- style_layout$axis_gap_pt %||% list(
  x_title_to_ticks = 1,
  x_ticks_to_axis  = 1,
  y_title_to_ticks = 3,
  y_ticks_to_axis  = 2
)

# 字体与字号（pt）
font_family <- style$typography$font_family_primary %||% "Helvetica"
axis_tick_pt  <- style$typography$sizes_pt$axis_tick_default  %||% 5.5
axis_label_pt <- style$typography$sizes_pt$axis_label_default %||% 6.5
legend_text_pt  <- style$typography$sizes_pt$legend_text_default  %||% 6
legend_title_pt <- style$typography$sizes_pt$legend_title_default %||% 6.5
panel_label_pt  <- style$typography$sizes_pt$panel_label_default %||% 8

# 布局边距（pt）
plot_margin_pt <- style$layout$plot_margin_pt %||% list(top = 0, right = 0, bottom = 0, left = 0)
axis_title_margin_pt <- style$layout$axis_title_margin_pt %||% list(x = 3, y = 3)

# 保证从 YAML 读取到的数值都是 length-1 的 numeric（避免 margin() 报错）
get_num1 <- function(x, default) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) return(default)
  as.numeric(x)[1]
}

pm_top    <- get_num1(plot_margin_pt$top,    0)
pm_right  <- get_num1(plot_margin_pt$right,  0)
pm_bottom <- get_num1(plot_margin_pt$bottom, 0)
pm_left   <- get_num1(plot_margin_pt$left,   0)

atm_x <- get_num1(axis_title_margin_pt$x, 3)
atm_y <- get_num1(axis_title_margin_pt$y, 3)

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

# 输出尺寸：宽度严格 84 mm；高度按轴数量自适应（与 box_panel_from_yaml.R 一致的布局计算需要用到）
width_mm <- 84
n_axis <- nlevels(dat$axis)
height_mm <- max(30, 6 * n_axis)  # 每个 axis 约 6 mm，下限 30 mm

# 将 axis_outer_frac 换算成 plot.margin 所需的 pt（基于当前 panel 的物理 size）
frac_or0 <- function(x) if (is.null(x)) 0 else x
margin_top_mm    <- frac_or0(axis_outer_frac$top)    * height_mm
margin_bottom_mm <- frac_or0(axis_outer_frac$bottom) * height_mm
margin_left_mm   <- frac_or0(axis_outer_frac$left)   * width_mm
margin_right_mm  <- frac_or0(axis_outer_frac$right)  * width_mm
plot_margin_cfg <- list(
  top    = mm_to_pt(margin_top_mm),
  right  = mm_to_pt(margin_right_mm),
  bottom = mm_to_pt(margin_bottom_mm),
  left   = mm_to_pt(margin_left_mm)
)

p <- ggplot(dat, aes(x = axis, y = Z_multiomics_overall)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = axis_line_lwd, color = "grey50") +
  geom_col(
    fill = "grey60",
    color = "grey40",
    width = 0.6,
    linewidth = bar_line_lwd
  ) +
  ## 星号标记
  geom_text(
    aes(
      label = p_label,
      y = y_star_pos
    ),
    vjust = ifelse(dat$Z_multiomics_overall >= 0, 0, 1),
    size = legend_text_pt / 2.5,
    color = "grey20",
    show.legend = FALSE
  ) +
  labs(
    x = NULL,  # 去除横坐标轴标题
    y = "Axis impairment (Z)",
    title = NULL,
    subtitle = NULL
  ) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = axis_tick_pt, base_family = font_family) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(
      family = font_family,
      size   = axis_label_pt,
      margin = margin(r = axis_gap_pt$y_title_to_ticks %||% 0)
    ),
    axis.text.x = element_text(
      angle = 30,
      hjust = 1,
      vjust = 1,
      family = font_family,
      size   = axis_tick_pt,
      margin = margin(t = axis_gap_pt$x_ticks_to_axis %||% 0)
    ),
    axis.text.y = element_text(
      family = font_family,
      size   = axis_tick_pt,
      margin = margin(r = axis_gap_pt$y_ticks_to_axis %||% 0)
    ),
    axis.line  = element_line(linewidth = axis_line_lwd, colour = "black"),
    axis.ticks = element_line(linewidth = axis_line_lwd, colour = "black"),
    plot.margin = margin(
      plot_margin_cfg$top    %||% 0,
      plot_margin_cfg$right  %||% 0,
      plot_margin_cfg$bottom %||% 0,
      plot_margin_cfg$left   %||% 0
    )
  )


ggsave(out_plot_path, p,
       width = width_mm, height = height_mm,
       units = "mm", device = grDevices::cairo_pdf)

cat("[OK] 成功输出图像: ", out_plot_path, "\n", sep = "")
cat("============================================================\n")
cat("[DONE] Figure 6E 绘制完成。\n")
cat("============================================================\n")