#!/usr/bin/env Rscript

# -------------
# 从原始 figure_6_b 表格中，抽取 axis / feature_type / name / bottleneck_sign / weight
# 然后按 feature_type + axis 聚合成一个“轴级摘要表”：
#   前几行为 feature_type = lipid、按 axis 排序
#   后几行为 feature_type = gene、按 axis 排序
# 每一行的内容 = 该 feature_type + axis 下所有指标
#   以 "name (bottleneck_sign, weight)" 形式，用逗号连接
#
# 生成表结构（用于排版进 Figure 6 的 PDF）：
#   Level | Axis | Indicators
#
# 默认输入：results/figs/figure_6_b.xlsx
# 输出：    results/figs/figure_6_b.pdf
#
# 导出逻辑（方案 A）：
#   1) gt::gtsave() → PNG（通过 vwidth 控制近似 83 mm 宽）
#   2) magick::image_read() 读 PNG
#   3) magick::image_trim() 裁掉所有多余白边
#   4) magick::image_write(..., format = "pdf") → 画布贴合内容的 PDF
# -------------

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(stringr)
  library(gt)
  library(magick)
})

# 1. 定义项目路径 & 输入输出路径 ----------------------------------------------

proj_dir <- getwd()

input_xlsx  <- file.path(proj_dir, "results", "figs", "figure_6_b.xlsx")
output_pdf  <- file.path(proj_dir, "results", "figs", "figure_6_b.pdf")

# 中间 PNG（放在临时目录，也可以改成 results/figs/ 下看得见）
tmp_png <- file.path(tempdir(), "figure_6_b_tmp.png")

message("Input:   ", input_xlsx)
message("Tmp PNG: ", tmp_png)
message("Output:  ", output_pdf)

# 2. 读入原始表，只保留需要的 5 列 -----------------------------------------

df_raw <- readxl::read_xlsx(input_xlsx)

# 确保列名存在，否则给出友好报错
required_cols <- c("axis", "feature_type", "name", "bottleneck_sign", "weight")
missing_cols  <- setdiff(required_cols, names(df_raw))
if (length(missing_cols) > 0) {
  stop(
    "下列列名在输入表中不存在，请检查：",
    paste(missing_cols, collapse = ", ")
  )
}

df <- df_raw %>%
  dplyr::select(all_of(required_cols))

# 3. 定义 axis 排序顺序（可按需要修改） -------------------------------------

preferred_axes <- c("Synthesis", "Oxidation", "Transport", "Supply", "Membrane context")

# 4. 定义一个函数：对某个 feature_type 进行 axis 维度的聚合 ------------------

aggregate_by_type <- function(dat, ft) {
  dat_ft <- dat %>% dplyr::filter(feature_type == ft)

  if (nrow(dat_ft) == 0L) {
    return(tibble())
  }

  dat_ft %>%
    dplyr::mutate(
      # 把每一行变成 "name (bottleneck_sign, weight)" 形式
      item = stringr::str_c(name, " (", bottleneck_sign, ", ", weight, ")")
    ) %>%
    dplyr::group_by(feature_type, axis) %>%
    dplyr::summarise(
      items = stringr::str_c(item, collapse = ", "),
      .groups = "drop"
    ) %>%
    # axis 排序：先按 preferred_axes 的顺序排，再把其他没列出的 axis 放后面
    dplyr::mutate(
      axis_factor = factor(axis, levels = preferred_axes)
    ) %>%
    dplyr::arrange(axis_factor, axis) %>%
    dplyr::select(-axis_factor)
}

# 5. 分别处理 lipid 和 gene，然后上下绑定 -------------------------------------

res_lipid <- aggregate_by_type(df, "lipid")
res_gene  <- aggregate_by_type(df, "gene")

# 确保顺序：先 lipid，后 gene
res_all <- dplyr::bind_rows(res_lipid, res_gene)

# 6. 为排版做一点列名与文本美化 ---------------------------------------------

plot_df <- res_all %>%
  dplyr::mutate(
    Level = dplyr::case_when(
      feature_type == "lipid" ~ "Lipid-level indicators",
      feature_type == "gene"  ~ "Gene-level indicators",
      TRUE                    ~ feature_type
    ),
    Axis = axis,
    Indicators = items
  ) %>%
  dplyr::select(Level, Axis, Indicators)

# 7. 用 gt 生成表格对象 ------------------------------------------------------

if (nrow(plot_df) == 0L) {
  stop("没有生成任何聚合结果，请检查原始表中的 feature_type / axis 是否正确。")
}

tab <- plot_df %>%
  gt::gt() %>%
  # 如果你不想在最终图里显示标题，可以把 tab_header 整块注释掉
  gt::tab_header(
    title = gt::md("**Figure 6B. Axis-level mechanistic indicators summary**")
  ) %>%
  # 列标题加粗
  gt::tab_style(
    style = list(gt::cell_text(weight = "bold")),
    locations = gt::cells_column_labels(gt::everything())
  ) %>%
  # 基本排版选项：字体、字号、行距
  gt::tab_options(
    table.font.names = "Arial",
    table.font.size  = 7,
    data_row.padding = gt::px(1)
  ) %>%
  # 缩小列宽，使总宽度适合 83 mm 画布（约 314 px）
  gt::cols_width(
    Level      ~ gt::px(70),
    Axis       ~ gt::px(70),
    Indicators ~ gt::px(170)
  ) %>%
  gt::cols_align(
    align = "center",
    columns = c(Level, Axis)
  ) %>%
  gt::cols_align(
    align = "left",
    columns = Indicators
  )

# 8. 第一步：gt → PNG（通过 webshot2 截图） ---------------------------------
# 说明：
# - gtsave(filename = *.png) 走的是 webshot2::webshot 路线；
# - vwidth/vheight 单位是像素；
# - 83 mm ≈ 314 px（按 96 dpi 估），这里 vwidth 取 320 px 稍微留一点余地；
# - vheight 取 600 px 足够容纳 12 行；expand = 0 减少额外白边（仍可能有少量）。

message("Step 1: save gt table to PNG ...")

gt::gtsave(
  data     = tab,
  filename = tmp_png,
  vwidth   = 320,   # 近似 83 mm 宽
  vheight  = 600,
  expand   = 0
)

# 9. 第二步：magick 读 PNG → 裁白边 → 写成 PDF -------------------------------

message("Step 2: trim PNG and convert to PDF ...")

img <- magick::image_read(tmp_png)

# 自动裁掉所有背景白边（以左上角像素为背景色）
img_trim <- magick::image_trim(img)

# 可选：设置一个合理的分辨率，方便在 AI 里看到物理尺寸
# 例如：300 dpi（不设置也可以，AI 里按像素缩放）
magick::image_write(
  image  = img_trim,
  path   = output_pdf,
  format = "pdf",
  density = "300x300"
)

message("Done: ", output_pdf)