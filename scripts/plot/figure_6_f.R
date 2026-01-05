#!/usr/bin/env Rscript

## ------------------------------------------------------------
## Figure 6F – mechanistic axis effect vs CL phenotypes (per batch)
##
## 用途：
##   从多组学机制轴效应矩阵 + 每批次脂质指标 HO–WT 差值表，
##   根据 YAML 配置，画出：
##     X = 轴效应(每批次 HO vs WT 的 combined Z)
##     Y = CL 关键表型 diff_HO_WT (HO - WT)
##   每个点 = 一个 batch。
## ------------------------------------------------------------

suppressPackageStartupMessages({
  library(yaml)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(ggplot2)
  library(broom)
  library(patchwork)
  library(grid)
})

## ---------------- 1. 读取命令行参数 & 配置 ------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("用法: Rscript figure_6_f.R path/to/figure_6_f.yaml\n", call. = FALSE)
}
config_path <- args[1]
cat("[INFO] 使用配置文件: ", config_path, "\n", sep = "")

cfg <- yaml::read_yaml(config_path)
# ---------------- Global figure style (Nature) ----------------
style_path <- "/Volumes/Samsung_SSD_990_PRO_2TB_Media/EphB1/02_protocols/figure_style_nature.yaml"
if (!file.exists(style_path)) {
  stop("[ERROR] 找不到 figure style 文件: ", style_path)
}
style <- yaml::read_yaml(style_path)


get_num1 <- function(x, default) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) return(default)
  as.numeric(x)[1]
}

# NULL-coalescing operator (avoid requiring rlang)
`%||%` <- function(a, b) if (!is.null(a)) a else b

# Fonts
font_family <- style$typography$font_family_primary %||% "Helvetica"

# Sizes (pt)
axis_tick_pt    <- get_num1(style$typography$sizes_pt$axis_tick_default, 5.5)
axis_label_pt   <- get_num1(style$typography$sizes_pt$axis_label_default, 6.5)
legend_text_pt  <- get_num1(style$typography$sizes_pt$legend_text_default, 6)
legend_title_pt <- get_num1(style$typography$sizes_pt$legend_title_default, 6.5)

# Lines (pt)
axis_line_pt   <- get_num1(style$lines$axis_line_default_pt, 0.25)
major_grid_pt  <- get_num1(style$lines$major_grid_default_pt, 0.35)
minor_grid_pt  <- get_num1(style$lines$minor_grid_default_pt, 0.25)

# pt -> lwd for ggplot linewidth (align with box_panel_from_yaml.R)
pt_per_lwd <- 0.75 # 1 lwd = 0.75 pt
pt_to_lwd <- function(pt) pt / pt_per_lwd
axis_line_lwd   <- pt_to_lwd(axis_line_pt)
major_grid_lwd  <- pt_to_lwd(major_grid_pt)
minor_grid_lwd  <- pt_to_lwd(minor_grid_pt)

# Colors
col_dir_match <- style$colors$dir_match %||% list(match="#B2182B", mismatch="grey40", unknown="black")
col_strip_bg  <- style$colors$strip_background %||% "grey90"
col_grid_major<- style$colors$grid_major %||% "grey85"
col_grid_minor<- style$colors$grid_minor %||% "grey92"

# Layout
legend_pos <- style$layout$legend_position_default %||% "bottom"
plot_margin_pt <- style$layout$plot_margin_pt %||% list(top=0, right=0, bottom=0, left=0)
pm_top    <- get_num1(plot_margin_pt$top, 0)
pm_right  <- get_num1(plot_margin_pt$right, 0)
pm_bottom <- get_num1(plot_margin_pt$bottom, 0)
pm_left   <- get_num1(plot_margin_pt$left, 0)

# Marks
point_size_mm     <- get_num1(style$marks$point_size, 1.6)
smooth_line_pt    <- get_num1(style$marks$smooth_line_width_pt, 0.5)
smooth_line_lwd   <- pt_to_lwd(smooth_line_pt)
smooth_se_alpha   <- get_num1(style$marks$smooth_se_alpha, 0.15)
facet_strip_text_pt <- get_num1(style$marks$facet_strip_text_pt, legend_text_pt)
stat_label_text_pt  <- get_num1(style$marks$stat_label_text_pt, legend_text_pt)
## 期望的 YAML 结构示例：
## data_files:
##   axis_effects_matrix: "results/multiomics/mechanistic_axis_effects_matrix.tsv"
##   lipid_indicator_de:  "results/lipid/tables/lipid_mechanistic_indicator_DE_per_batch.tsv"
##   lipid_indicator_yaml: "scripts/lipid/lipid_mechanistic_indicators.yaml"
## axes_to_plot:
##   - Synthesis
##   - Oxidation
## axis_phenotype_map:
##   Synthesis: ["sumCL"]
##   Oxidation: ["oxCL_CL_ratio"]
## expected_direction:
##   Synthesis: "negative"
##   Oxidation: "positive"
## output:
##   pdf: "results/figs/figure_6_f.pdf"
##   png: "results/figs/figure_6_f.png"
##   stats_tsv: "results/stats/figure_6_f_all_stats.tsv"  # optional
## layout:
##   nrow: 2
##   ncol: 1
## stats:
##   method: "pearson"

get_cfg <- function(lst, ..., default = NULL) {
  # 小工具：安全读取嵌套字段
  keys <- c(...)
  cur <- lst
  for (k in keys) {
    if (is.null(cur[[k]])) return(default)
    cur <- cur[[k]]
  }
  cur
}

axis_mat_path <- get_cfg(cfg, "data_files", "axis_effects_matrix")
lipid_de_path <- get_cfg(cfg, "data_files", "lipid_indicator_de")
lipid_ind_yaml_path <- get_cfg(cfg, "data_files", "lipid_indicator_yaml")

if (is.null(axis_mat_path) || is.null(lipid_de_path)) {
  stop("[ERROR] 配置文件中的 data_files.axis_effects_matrix / lipid_indicator_de 未设置\n")
}

axes_to_plot <- get_cfg(cfg, "axes_to_plot", default = NULL)
axis_pheno_map <- get_cfg(cfg, "axis_phenotype_map", default = NULL)
expected_dir   <- get_cfg(cfg, "expected_direction", default = list())

out_pdf <- get_cfg(cfg, "output", "pdf", default = "results/figs/figure_6_f.pdf")
out_png <- get_cfg(cfg, "output", "png", default = "results/figs/figure_6_f.png")
out_stats_tsv <- get_cfg(cfg, "output", "stats_tsv", default = NULL)
out_full_pdf <- get_cfg(cfg, "output", "full_pdf", default = "results/figs/figure_6_f_full/figure_6_f_full.pdf")
out_full_png <- get_cfg(cfg, "output", "full_png", default = "results/figs/figure_6_f_full/figure_6_f_full.png")

layout_nrow <- get_cfg(cfg, "layout", "nrow", default = NULL)
layout_ncol <- get_cfg(cfg, "layout", "ncol", default = NULL)

stat_method <- get_cfg(cfg, "stats", "method", default = "pearson")
transform_y_mode <- get_cfg(cfg, "transform_y", "mode", default = "none")

width_cfg  <- get_cfg(cfg, "plot_size", "width",  default = 8)
height_cfg <- get_cfg(cfg, "plot_size", "height", default = 6)

# Main figure size: strict single-column width (84 mm). Height defaults to YAML inches converted to mm.
width_mm_main  <- get_cfg(cfg, "plot_size", "width_mm",  default = 84)
height_mm_main <- get_cfg(cfg, "plot_size", "height_mm", default = height_cfg * 25.4)

# 强制单栏宽度为 84 mm（忽略 config 中可能的误设），并转换为英寸用于 device 尺寸锁定
width_mm_main <- 84
width_in_main  <- width_mm_main / 25.4
height_in_main <- height_mm_main / 25.4
cat(sprintf("[INFO] Main plot size locked: %.2f mm (%.4f in) × %.2f mm (%.4f in)\n",
            width_mm_main, width_in_main, height_mm_main, height_in_main))

# After rendering, enforce PDF CropBox/MediaBox to exact dimensions (points)
mm_to_pt <- function(mm) mm * 72 / 25.4
page_w_pt <- mm_to_pt(width_mm_main)
page_h_pt <- mm_to_pt(height_mm_main)

fix_pdf_boxes <- function(pdf_path, w_pt, h_pt) {
  # Use Ghostscript if available to force an exact page size/crop box.
  gs_bin <- Sys.which("gs")
  if (!nzchar(gs_bin)) return(invisible(FALSE))

  tmp_out <- paste0(pdf_path, ".tmp.pdf")
  args <- c(
    "-o", tmp_out,
    "-sDEVICE=pdfwrite",
    "-dSAFER",
    "-dBATCH",
    "-dNOPAUSE",
    "-dUseCropBox",
    sprintf("-dDEVICEWIDTHPOINTS=%.2f", w_pt),
    sprintf("-dDEVICEHEIGHTPOINTS=%.2f", h_pt),
    "-dPDFFitPage",
    pdf_path
  )
  system2(gs_bin, args = args)
  if (file.exists(tmp_out) && file.info(tmp_out)$size > 0) {
    ok <- file.rename(tmp_out, pdf_path)
    return(invisible(isTRUE(ok)))
  }
  invisible(FALSE)
}

# full 图（每个 phenotype 单独图）的尺寸：同样锁定为单栏宽 84 mm；高度可单独配置。
full_width_mm  <- get_cfg(cfg, "plot_size", "full_width_mm",  default = 84)
full_height_mm <- get_cfg(cfg, "plot_size", "full_height_mm", default = full_width_mm * 1.25)  # 默认给一个相对合理的纵横比

# 强制单栏宽度 84 mm
full_width_mm <- 84
full_width_in  <- full_width_mm / 25.4
full_height_in <- full_height_mm / 25.4

## ---------------- 2. 读取数据表 -----------------------------

cat("[STEP1] 读取机制轴效应矩阵: ", axis_mat_path, "\n", sep = "")
axis_mat <- readr::read_tsv(axis_mat_path, col_types = cols())

## 预期结构（示意）：
## axis    RNA_batch1  RNA_batch2  Lipid_batch1  ...
## Oxidation   1.08      1.10        -0.51      ...
## ...
if (!"axis" %in% colnames(axis_mat)) {
  stop("[ERROR] mechanistic_axis_effects_matrix.tsv 中必须包含列 'axis'\n")
}

cat("[INFO] axis_mat: nrow = ", nrow(axis_mat), "; ncol = ", ncol(axis_mat), "\n", sep = "")

cat("[STEP2] 读取脂质指标 HO-WT 差值表: ", lipid_de_path, "\n", sep = "")
lipid_de <- readr::read_tsv(lipid_de_path, col_types = cols())

## 预期结构（示意）：
## batch  indicator          mean_HO  mean_WT  diff_HO_WT   p_value  FDR
## B1     sumCL              ...      ...      ...          ...      ...
## B1     oxCL_CL_ratio      ...      ...      ...          ...      ...

needed_cols <- c("batch", "indicator", "diff_HO_WT")
missing_cols <- setdiff(needed_cols, colnames(lipid_de))
if (length(missing_cols) > 0) {
  stop("[ERROR] lipid_mechanistic_indicator_DE_per_batch.tsv 缺少列: ",
       paste(missing_cols, collapse = ", "), "\n")
}

cat("[INFO] lipid_de: nrow = ", nrow(lipid_de), "; ncol = ", ncol(lipid_de), "\n", sep = "")

## ---------------- 3. 整理轴效应：axis × batch --------------

## 3.1 宽转长：把每个 omic × batch 列转成长格式
## 注意：这里假定除了 'axis' 之外的列名编码了 omic & batch 信息，
## 例如 "RNA_batch1", "Lipid_batch1" 等。
## 你可以根据实际列名改动下面的正则解析部分。

axis_long <- axis_mat %>%
  tidyr::pivot_longer(
    cols      = -axis,
    names_to  = "source",
    values_to = "effect_z"
  ) %>%
  mutate(
    ## 简单示例：根据列名前缀判断 omic
    omic = case_when(
      str_detect(source, regex("rna", ignore_case = TRUE))   ~ "RNA",
      str_detect(source, regex("lipid", ignore_case = TRUE)) ~ "Lipid",
      TRUE                                                   ~ "Other"
    ),
    ## 解析 batch 名：根据实际列名修改
    batch = str_extract(source, "batch[0-9A-Za-z_]+")
  )

cat("[INFO] axis_long: nrow = ", nrow(axis_long), "\n", sep = "")

## 如果列名格式比较规整，可以在这里检查一下
if (any(is.na(axis_long$batch))) {
  warning("[WARN] 有部分列无法解析出 batch，建议检查 mechanistic_axis_effects_matrix.tsv 的列名格式。\n")
}

## 3.2 对同一 axis × batch 的多组学效应做 combined Z
## 简单版：Stouffer Z = sum(z) / sqrt(n)
axis_combined <- axis_long %>%
  filter(!is.na(batch)) %>%
  group_by(axis, batch) %>%
  summarise(
    n_components = sum(!is.na(effect_z)),
    combined_z   = if (n_components > 0) sum(effect_z, na.rm = TRUE) / sqrt(n_components) else NA_real_,
    .groups      = "drop"
  )

cat("[INFO] axis_combined: nrow = ", nrow(axis_combined), "\n", sep = "")

## 用于构造 axis × phenotype 的全组合
axes_all <- sort(unique(axis_combined$axis))

## ---------------- 4. 读取脂质指标 YAML，构造 axis × phenotype 全组合 ---------------------
##
## 方案 a：每个 axis 与所有指标做回归分析：
##   - 从 YAML 中只提取“所有 mechanistic 指标的列表”，不再用 mechanistic_axes 来限制 axis–indicator 映射；
##   - 使用 expand_grid(axis = axes_all, phenotype = phenos_all) 构造全组合。
## -------------------------------------------------------------------

if (!is.null(lipid_ind_yaml_path)) {
  cat("[STEP3] 读取脂质指标 YAML: ", lipid_ind_yaml_path, "\n", sep = "")
  lipid_ind_yaml <- yaml::read_yaml(lipid_ind_yaml_path)
  # 结构示例：
  # mechanistic_axes:
  #   Synthesis:
  #     indicators:
  #       - name: sumCL
  #       - name: ...
  #   Oxidation:
  #     indicators:
  #       - name: oxCL_CL_ratio
  #       - name: ...
  mech_axes <- lipid_ind_yaml$mechanistic_axes
  if (is.null(mech_axes)) {
    stop("[ERROR] lipid_indicator_yaml 文件中缺少 mechanistic_axes 字段\n")
  }

  ## 从所有 axis 的 indicators 中提取去重的指标列表
  phenos_all <- purrr::imap_dfr(mech_axes, function(ax_info, ax_name) {
    inds <- ax_info$indicators
    if (is.null(inds)) {
      return(tibble(phenotype = character(0)))
    }
    tibble(
      phenotype = purrr::map_chr(inds, "name")
    )
  }) %>%
    distinct(phenotype) %>%
    dplyr::pull(phenotype)

  cat("[INFO] 从 lipid_indicator_yaml 中获得指标种类数 (去重): ", length(phenos_all), "\n", sep = "")

  ## 构造 axis_indicator_pairs_all：每个 axis 与所有 phenotype 的笛卡尔积
  axis_indicator_pairs_all <- tidyr::expand_grid(
    axis      = axes_all,
    phenotype = phenos_all
  )

  cat("[INFO] 按方案 a 使用 axis × phenotype 全组合对数: ", nrow(axis_indicator_pairs_all), "\n", sep = "")
} else {
  stop("[ERROR] 配置文件中缺少 data_files.lipid_indicator_yaml，无法获取全部 axis-indicator 对\n")
}

## ---------------- 5. 过滤 lipid_de，只保留 axis_indicator_pairs_all 中的phenotype ---------------------

phenos_all <- unique(axis_indicator_pairs_all$phenotype)
lipid_de_filt <- lipid_de %>%
  filter(indicator %in% phenos_all) %>%
  select(batch, indicator, diff_HO_WT)

## ---------------- 6. 构造绘图用 axis-phenotype 对 ---------------------

## 如果 axis_phenotype_map 提供，则用它来限制绘图对；否则绘图对 = 全部对
if (is.null(axis_pheno_map)) {
  axis_indicator_pairs_plot <- axis_indicator_pairs_all
} else {
  axis_indicator_pairs_plot <- purrr::imap_dfr(axis_pheno_map, function(phenos, ax) {
    tibble(
      axis = ax,
      phenotype = phenos
    )
  })
}

## ---------------- 7. 组装作图用数据 -------------------------

## 7.1 全部 axis_indicator_pairs_all 与 axis_combined 和 lipid_de_filt 连接，准备计算统计量
full_data <- axis_indicator_pairs_all %>%
  left_join(axis_combined, by = "axis") %>%
  left_join(lipid_de_filt, by = c("batch", "phenotype" = "indicator"))

## 7.2 过滤 NA
full_data <- full_data %>%
  filter(!is.na(combined_z), !is.na(diff_HO_WT))

cat("[INFO] 用于统计分析的点数: ", nrow(full_data), "\n", sep = "")

## 为 full 图和主图准备 y 变量名称
y_var_full <- "diff_HO_WT"
y_lab_full <- "CL phenotype difference (HO - WT)"

if (!is.null(transform_y_mode) && identical(transform_y_mode, "z_by_indicator")) {
  if (nrow(full_data) > 0) {
    full_data <- full_data %>%
      group_by(phenotype) %>%
      mutate(
        diff_HO_WT_z = ifelse(
          sd(diff_HO_WT, na.rm = TRUE) > 0,
          (diff_HO_WT - mean(diff_HO_WT, na.rm = TRUE)) / sd(diff_HO_WT, na.rm = TRUE),
          0
        )
      ) %>%
      ungroup()
    y_var_full <- "diff_HO_WT_z"
    y_lab_full <- "Standardized CL phenotype difference (z, HO - WT)"
  }
}

## ---------------- 8. 计算所有 axis-phenotype 对的相关系数/回归 ---------

calc_cor <- function(x, y, method = "pearson") {
  ## 小样本下这一步要 tryCatch，以免 cor.test 报错
  out <- tryCatch({
    ct <- suppressWarnings(cor.test(x, y, method = method))
    tibble(
      r      = unname(ct$estimate),
      p      = ct$p.value
    )
  }, error = function(e) {
    tibble(r = NA_real_, p = NA_real_)
  })
  out
}

stats_df_all <- full_data %>%
  group_by(axis, phenotype) %>%
  summarise(
    n_points = n(),
    ## 线性回归，用于得到 slope
    {
      fit <- tryCatch(lm(diff_HO_WT ~ combined_z), error = function(e) NULL)
      if (is.null(fit)) {
        slope <- NA_real_
      } else {
        slope <- coef(fit)[["combined_z"]]
      }
      tibble(slope = slope)
    },
    ## 相关
    calc_cor(combined_z, diff_HO_WT, method = stat_method),
    .groups = "drop"
  )

## ---------------- 9. 根据 expected_direction 加入方向匹配标签 -----------------

stats_df_all <- stats_df_all %>%
  mutate(
    expected_dir = expected_dir[axis] %||% NA_character_,
    dir_match = case_when(
      is.na(slope) | is.na(expected_dir) ~ "unknown",
      expected_dir == "positive" & slope > 0 ~ "match",
      expected_dir == "negative" & slope < 0 ~ "match",
      TRUE ~ "mismatch"
    ),
    label = sprintf("r = %.2f\np = %.2f", r, p)
  )

## ---------------- 9.1 为全轴全指标图准备数据（可选补充图） -----------------
##
## full_plot_data 使用 full_data（所有 axis-phenotype 对），
## 并合并 stats_df_all 以便显示 r/p 和 dir_match。
## 该数据集将用于可选的“全轴全指标”补充图。
## ------------------------------------------------------------

full_plot_data <- full_data %>%
  left_join(stats_df_all, by = c("axis", "phenotype")) %>%
  mutate(
    axis      = factor(axis, levels = unique(axis)),
    phenotype = factor(phenotype)
  )

## ---------------- 10. 准备绘图数据 -------------------------

## 10.1 只用 axis_indicator_pairs_plot 过滤 full_data 作为绘图数据
plot_data <- axis_indicator_pairs_plot %>%
  left_join(axis_combined, by = "axis") %>%
  left_join(lipid_de_filt, by = c("batch", "phenotype" = "indicator")) %>%
  filter(!is.na(combined_z), !is.na(diff_HO_WT))

## 10.2 过滤 stats_df_all 只保留绘图用的 axis-phenotype 对
stats_df_plot <- stats_df_all %>%
  semi_join(axis_indicator_pairs_plot, by = c("axis", "phenotype"))

## 10.3 将 stats_df_plot 合并到 plot_data，用于 panel 标签和颜色映射
plot_data <- plot_data %>%
  left_join(stats_df_plot, by = c("axis", "phenotype"))

## ---------------- 11. 可选：对 y 做标准化（按 phenotype）-------------------

# 默认使用原始 diff_HO_WT；如果 transform_y_mode == "z_by_indicator"，
# 则对每个 phenotype 在所有批次内做 z-score，用于绘图，使不同指标的 y 轴更可比。
y_var <- "diff_HO_WT"
y_lab <- "CL phenotype difference (HO - WT)"

if (!is.null(transform_y_mode) && identical(transform_y_mode, "z_by_indicator")) {
  if (nrow(plot_data) > 0) {
    plot_data <- plot_data %>%
      group_by(phenotype) %>%
      mutate(
        diff_HO_WT_z = ifelse(
          sd(diff_HO_WT, na.rm = TRUE) > 0,
          (diff_HO_WT - mean(diff_HO_WT, na.rm = TRUE)) / sd(diff_HO_WT, na.rm = TRUE),
          0
        )
      ) %>%
      ungroup()
    y_var <- "diff_HO_WT_z"
    y_lab <- "Standardized CL phenotype difference (z, HO - WT)"
  }
}

## ---------------- 11. 画图 -----------------------------------

## 如配置要求，先输出“全轴全指标”补充图：
## 方案 A（更新版）：每个 phenotype 单独一张图，分别保存到 results/figs/figure_6_f_full/ 目录下
if (!is.null(out_full_pdf) || !is.null(out_full_png)) {

  # 按 phenotype 拆分数据，生成每个 phenotype 的 ggplot 对象
  full_plots_list <- full_plot_data %>%
    split(.$phenotype) %>%
    purrr::imap(function(df, pheno_name) {
      # 对应 phenotype 的统计结果（r, p, slope 等）
      stats_sub <- stats_df_all %>%
        dplyr::filter(phenotype == pheno_name)

      ggplot(df, aes(x = combined_z, y = .data[[y_var_full]])) +
        geom_point(size = point_size_mm) +
        geom_smooth(
          method = "lm",
          se = TRUE,
          aes(color = dir_match),
          show.legend = TRUE,
          linewidth = smooth_line_lwd,
          alpha = smooth_se_alpha
        ) +
        facet_grid(axis ~ ., scales = "free") +
        geom_text(
          data = stats_sub,
          aes(
            x = -Inf, y = Inf,
            label = label
          ),
          hjust = -0.1, vjust = 1.1,
          inherit.aes = FALSE,
          size = stat_label_text_pt / ggplot2::.pt
        ) +
        scale_color_manual(
          values = c(
            match    = col_dir_match$match,
            mismatch = col_dir_match$mismatch,
            unknown  = col_dir_match$unknown
          )
        ) +
        labs(
          title = pheno_name,
          x = "Multi-omics mechanistic axis effect (HO vs WT, combined Z)",
          y = y_lab_full,
          color = "Direction match"
        ) +
        theme_bw(base_size = axis_tick_pt, base_family = font_family) +
        theme(
          strip.background = element_rect(fill = col_strip_bg, color = NA),
          strip.text = element_text(size = facet_strip_text_pt, family = font_family),

          panel.grid.major = element_line(linewidth = major_grid_lwd, color = col_grid_major),
          panel.grid.minor = element_line(linewidth = minor_grid_lwd, color = col_grid_minor),

          axis.line  = element_line(linewidth = axis_line_lwd, color = "black"),
          axis.ticks = element_line(linewidth = axis_line_lwd, color = "black"),

          axis.text.x  = element_text(size = axis_tick_pt, family = font_family),
          axis.text.y  = element_text(size = axis_tick_pt, family = font_family),
          axis.title.x = element_text(size = axis_label_pt, family = font_family),
          axis.title.y = element_text(size = axis_label_pt, family = font_family),

          legend.position = legend_pos,
          legend.title = element_text(size = legend_title_pt, family = font_family),
          legend.text  = element_text(size = legend_text_pt,  family = font_family),

          plot.margin = margin(pm_top, pm_right, pm_bottom, pm_left, unit = "pt")
        )
    })

  # 目标输出目录：优先使用 out_full_pdf/out_full_png 的目录，否则默认 results/figs/figure_6_f_full
  out_full_dir <- NULL
  if (!is.null(out_full_pdf)) {
    out_full_dir <- dirname(out_full_pdf)
  } else if (!is.null(out_full_png)) {
    out_full_dir <- dirname(out_full_png)
  } else {
    out_full_dir <- "results/figs/figure_6_f_full"
  }
  dir.create(out_full_dir, showWarnings = FALSE, recursive = TRUE)

  # 逐个 phenotype 保存单独文件
  pheno_names <- names(full_plots_list)
  purrr::iwalk(full_plots_list, function(p_single, pheno_name) {
    # 将 phenotype 名字中的特殊字符替换为下划线，防止文件名非法
    pheno_safe <- stringr::str_replace_all(pheno_name, "[^A-Za-z0-9_]+", "_")
    file_base <- file.path(out_full_dir, paste0("figure_6_f_full_", pheno_safe))

    if (!is.null(out_full_pdf)) {
      out_pdf_single <- paste0(file_base, ".pdf")
      grDevices::cairo_pdf(
        file   = out_pdf_single,
        width  = full_width_in,
        height = full_height_in,
        family = font_family
      )
      print(p_single)
      dev.off()

      # enforce exact page boxes for Illustrator placement
      fixed_full <- fix_pdf_boxes(out_pdf_single, mm_to_pt(full_width_mm), mm_to_pt(full_height_mm))
      if (!isTRUE(fixed_full)) {
        cat("[WARN] full PDF: could not enforce CropBox/MediaBox (gs not found): ", out_pdf_single, "\n", sep = "")
      }
    }
    if (!is.null(out_full_png)) {
      ggsave(paste0(file_base, ".png"), p_single,
             width  = full_width_in,
             height = full_height_in,
             units  = "in", dpi = 300)
    }
  })
}

## 设置 factor 顺序（可按需要调整）
plot_data <- plot_data %>%
  mutate(
    axis      = factor(axis, levels = axes_to_plot %||% unique(axis)),
    phenotype = factor(phenotype)
  )

p <- ggplot(plot_data, aes(x = combined_z, y = .data[[y_var]])) +
  geom_point(size = point_size_mm) +
  ## 用 dir_match 映射颜色（match/mismatch/unknown），可以在 scale_color_manual 里自定义
  geom_smooth(
    method = "lm",
    se = TRUE,
    aes(color = dir_match),
    show.legend = TRUE,
    linewidth = smooth_line_lwd,
    alpha = smooth_se_alpha
  ) +
  facet_grid(axis ~ phenotype, scales = "free") +
  ## 在每个 panel 右上角写 r, p
  geom_text(
    data = stats_df_plot,
    aes(
      x = -Inf, y = Inf,
      label = label
    ),
    hjust = -0.1, vjust = 1.1,
    inherit.aes = FALSE,
    size = stat_label_text_pt / ggplot2::.pt
  ) +
  scale_color_manual(
    values = c(
      match    = col_dir_match$match,
      mismatch = col_dir_match$mismatch,
      unknown  = col_dir_match$unknown
    )
  ) +
  labs(
    x = "Multi-omics mechanistic axis effect (HO vs WT, combined Z)",
    y = y_lab,
    color = "Direction match"
  ) +
  theme_bw(base_size = axis_tick_pt, base_family = font_family) +
  theme(
    strip.background = element_rect(fill = col_strip_bg, color = NA),
    strip.text = element_text(size = facet_strip_text_pt, family = font_family),

    panel.grid.major = element_line(linewidth = major_grid_lwd, color = col_grid_major),
    panel.grid.minor = element_line(linewidth = minor_grid_lwd, color = col_grid_minor),

    axis.line  = element_line(linewidth = axis_line_lwd, color = "black"),
    axis.ticks = element_line(linewidth = axis_line_lwd, color = "black"),

    axis.text.x  = element_text(size = axis_tick_pt, family = font_family),
    axis.text.y  = element_text(size = axis_tick_pt, family = font_family),
    axis.title.x = element_text(size = axis_label_pt, family = font_family),
    axis.title.y = element_text(size = axis_label_pt, family = font_family),

    legend.position = legend_pos,
    legend.title = element_text(size = legend_title_pt, family = font_family),
    legend.text  = element_text(size = legend_text_pt,  family = font_family),

    plot.margin = margin(pm_top, pm_right, pm_bottom, pm_left, unit = "pt")
  )

## ---------------- 12. 保存图像和统计结果 -------------------


## ---- 导出 PDF（cairo_pdf + print）----
dir.create(dirname(out_pdf), showWarnings = FALSE, recursive = TRUE)
grDevices::cairo_pdf(
  file   = out_pdf,
  width  = width_in_main,
  height = height_in_main,
  family = font_family
)
print(p)
dev.off()

# Enforce page boxes for Illustrator placement
fixed <- fix_pdf_boxes(out_pdf, page_w_pt, page_h_pt)
if (isTRUE(fixed)) {
  cat(sprintf("[INFO] PDF boxes fixed to %.2f×%.2f pt (%.2f×%.2f mm)\n",
              page_w_pt, page_h_pt, width_mm_main, height_mm_main))
} else {
  cat("[WARN] Could not enforce PDF CropBox/MediaBox (gs not found). PDF is still vector; Illustrator may place by bounding box.\n")
}

## ---- 导出 PNG（用于快速预览）----
dir.create(dirname(out_png), showWarnings = FALSE, recursive = TRUE)
ggsave(out_png, p,
       width = width_in_main, height = height_in_main,
       units = "in", dpi = 300)

if (!is.null(out_stats_tsv)) {
  dir.create(dirname(out_stats_tsv), showWarnings = FALSE, recursive = TRUE)
  readr::write_tsv(stats_df_all, out_stats_tsv)
  cat("[INFO] 全部 axis-phenotype 统计结果已保存至: ", out_stats_tsv, "\n", sep = "")
}

cat("[OK] Figure 6F 已输出到:\n  - ", out_pdf, "\n  - ", out_png, "\n", sep = "")