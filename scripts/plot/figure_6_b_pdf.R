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
#   1) gt::gtsave() → PNG（通过 vwidth 控制近似 84 mm 宽）
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
  library(grid)
  library(gridExtra)
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

# --- 固定输出物理尺寸（严格 84 mm 宽）与字号（按 300 dpi 换算） -----------------
# 最终 PDF 使用 density = 300 dpi，因此：
#   px = inch * dpi；pt = inch * 72  ⇒  px = pt/72 * dpi
DPI_OUT <- 300
TARGET_WIDTH_MM <- 84
TARGET_WIDTH_PX <- round(TARGET_WIDTH_MM / 25.4 * DPI_OUT)  # 84 mm @ 300 dpi ≈ 992 px

BODY_PT  <- 3
LABEL_PT <- 4
BODY_PX  <- BODY_PT  / 72 * DPI_OUT   # 2 pt @ 300 dpi ≈ 8.33 px
LABEL_PX <- LABEL_PT / 72 * DPI_OUT   # 3 pt @ 300 dpi ≈ 12.5 px

# 线条粗细（按 300 dpi 换算）
THICK_PT <- 0.5
THIN_PT  <- 0.25
THICK_PX <- THICK_PT / 72 * DPI_OUT   # 0.5 pt @ 300 dpi ≈ 2.08 px
THIN_PX  <- THIN_PT  / 72 * DPI_OUT   # 0.25 pt @ 300 dpi ≈ 1.04 px

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

# 6.5 依据内容自适应列宽（Level / Axis），为 Indicators 留足空间 ---------------
# 经验：Arial 在 6 pt 下平均字符宽度约为字号的 ~0.55 倍（像素单位）
char_px <- BODY_PX * 0.55
pad_px  <- 28  # 预留左右内边距/分隔

level_max_chars <- max(nchar(plot_df$Level), na.rm = TRUE)
axis_max_chars  <- max(nchar(plot_df$Axis),  na.rm = TRUE)

# 先按内容估计宽度，并做上/下限约束
level_w <- level_max_chars * char_px + pad_px
axis_w  <- axis_max_chars  * char_px + pad_px

# 上下限：避免前两列过宽挤占 Indicators
level_w <- max(180, min(320, level_w))
axis_w  <- max(120, min(260, axis_w))

# 保障 Indicators 至少有 520 px
min_ind_w <- 520
ind_w <- TARGET_WIDTH_PX - level_w - axis_w
if (ind_w < min_ind_w) {
  # 优先压缩 Level / Axis 到各自下限，确保 Indicators
  deficit <- (min_ind_w - ind_w)
  shrink_level <- min(deficit * 0.6, level_w - 180)
  level_w <- level_w - shrink_level
  deficit <- deficit - shrink_level
  shrink_axis <- min(deficit, axis_w - 120)
  axis_w <- axis_w - shrink_axis
  ind_w <- TARGET_WIDTH_PX - level_w - axis_w
}

# 最终取整，避免小数像素带来的渲染抖动
level_w <- as.integer(round(level_w))
axis_w  <- as.integer(round(axis_w))
ind_w   <- as.integer(round(ind_w))

# 6.6 为 Indicators 自动换行（通过 \n 触发行高自适应） ---------------------------
# 估算每行可容纳字符数：Indicators 列宽 / 单字符像素宽度
ind_pad_px <- 24
wrap_chars <- max(20, floor((ind_w - ind_pad_px) / char_px))

plot_df2 <- plot_df %>%
  dplyr::mutate(
    Indicators = stringr::str_wrap(Indicators, width = wrap_chars)
  )

# 7. 生成矢量 PDF（可在 Illustrator 中保持为矢量文字/线条） -----------------------

if (nrow(plot_df2) == 0L) {
  stop("没有生成任何聚合结果，请检查原始表中的 feature_type / axis 是否正确。")
}

# 目标物理宽度（英寸）
W_IN <- TARGET_WIDTH_MM / 25.4

# 线宽（grid 的 lwd 单位是“线宽倍数”，这里直接用 pt 更直观：在 PDF 中仍为矢量）
# 通过 gpar(lwd=...) 传递，通常在 pdf 设备中表现为 pt 级线宽。

# 构建表格主题：
# - 正文 3 pt Helvetica
# - 表头 4 pt Helvetica Bold
# - 细线 0.25 pt（行分隔）
# - 粗线 0.5 pt（表头分隔/外框）

thm <- gridExtra::ttheme_minimal(
  core = list(
    fg_params = list(fontfamily = "Helvetica", fontsize = BODY_PT, hjust = 0, x = 0.02, y = 0.5),
    bg_params = list(fill = "white", col = NA)
  ),
  colhead = list(
    fg_params = list(fontfamily = "Helvetica", fontsize = LABEL_PT, fontface = "bold", hjust = 0.5, x = 0.5),
    bg_params = list(fill = "white", col = NA)
  )
)

# 生成 table grob（Indicators 已经包含换行符，会自动增高行高）
tg <- gridExtra::tableGrob(plot_df2, rows = NULL, theme = thm)

# 进一步压缩行高（对所有数据行生效，不影响换行逻辑）
# 系数 < 1 会让行更紧凑
compact_factor <- 0.85
# 第 1 行是表头，从第 2 行开始是正文
for (i in seq_along(tg$heights)) {
  if (i > 1) tg$heights[i] <- tg$heights[i] * compact_factor
}

# 对齐：Level/Axis 居中，Indicators 左对齐
# （tableGrob 的列从 1 开始；行 1 是表头）
# 核心单元格在 tg$grobs 中，使用 gtable 来微调对齐
suppressWarnings({
  tg <- gtable::gtable_add_grob(tg, tg$grobs, t = 1, l = 1)
})

# 设置列宽：使用“毫米”等效宽度，确保总体宽度 = 84 mm
# 注意：grid 的 unit 可以使用 "null"/"mm"，这里用 mm 直接锁物理尺寸
col_w_mm <- c(level_w, axis_w, ind_w) / DPI_OUT * 25.4
# 防御：避免极端情况下出现负数/0
col_w_mm[col_w_mm < 5] <- 5

tg$widths[1] <- grid::unit(col_w_mm[1], "mm")
tg$widths[2] <- grid::unit(col_w_mm[2], "mm")
tg$widths[3] <- grid::unit(col_w_mm[3], "mm")

# --- 自定义横线：仅行间细线 + 表头上下粗线（无竖线、无外框） -------------------
# 表格有 1 行表头 + N 行数据
n_data_rows <- nrow(plot_df2)

# 一条横线 grob（在所处单元格内的 y=0 或 y=1 位置画线，跨越所有列）
make_hline <- function(y_npc, lwd_pt) {
  grid::segmentsGrob(
    x0 = grid::unit(0, "npc"), x1 = grid::unit(1, "npc"),
    y0 = grid::unit(y_npc, "npc"), y1 = grid::unit(y_npc, "npc"),
    gp = grid::gpar(col = "black", lwd = lwd_pt, lineend = "butt")
  )
}

# 表头上下各一条 0.5 pt 横线（明确保证列标题行下方有粗线）
# 表头：上粗线（y=1）+ 下粗线（y=0）放在第 1 行
n_cols <- ncol(plot_df2)
header_top <- make_hline(1, THICK_PT)
# 放在表头行内稍微上移，避免刚好落在边界被渲染吞掉
header_bot <- make_hline(0.02, THICK_PT)

tg <- gtable::gtable_add_grob(tg, header_top, t = 1, l = 1, r = n_cols, z = Inf, clip = "off")
tg <- gtable::gtable_add_grob(tg, header_bot, t = 1, l = 1, r = n_cols, z = Inf, clip = "off")

# 进一步保险：在第一条数据行顶部再画一条 0.5 pt 粗线（与表头下边界重合）
if (n_data_rows > 0) {
  tg <- gtable::gtable_add_grob(
    tg,
    make_hline(0.98, THICK_PT),
    t = 2, l = 1, r = n_cols,
    z = Inf, clip = "off"
  )
}

# 数据行：每行底部一条细线（y=0）
if (n_data_rows > 0) {
  for (i in seq_len(n_data_rows)) {
    row_idx <- 1 + i  # tableGrob: 1=header, 2..=data
    tg <- gtable::gtable_add_grob(
      tg,
      make_hline(0, THIN_PT),
      t = row_idx, l = 1, r = n_cols,
      z = Inf, clip = "off"
    )
  }
}

# 根据 gtable 行高总和自动计算 PDF 高度（英寸），确保完整显示
# 使用 sum(tg$heights) 比 grobHeight() 更稳（避免在某些设备上低估高度）
H_IN <- grid::convertHeight(sum(tg$heights), "in", valueOnly = TRUE) + 0.20

message("Step: write vector PDF ...")

grDevices::pdf(
  file = output_pdf,
  width = W_IN,
  height = H_IN,
  family = "Helvetica",
  useDingbats = FALSE
)

grid::grid.newpage()
# 仅绘制表格（横线已作为 grob 叠加在 tg 内部）
grid::grid.draw(tg)

grDevices::dev.off()

message("Done (vector): ", output_pdf)