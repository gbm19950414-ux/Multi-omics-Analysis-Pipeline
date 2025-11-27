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
})

## ---------------- 1. 读取命令行参数 & 配置 ------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("用法: Rscript figure_6_f.R path/to/figure_6_f.yaml\n", call. = FALSE)
}
config_path <- args[1]
cat("[INFO] 使用配置文件: ", config_path, "\n", sep = "")

cfg <- yaml::read_yaml(config_path)

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
out_full_pdf <- get_cfg(cfg, "output", "full_pdf", default = NULL)
out_full_png <- get_cfg(cfg, "output", "full_png", default = NULL)

layout_nrow <- get_cfg(cfg, "layout", "nrow", default = NULL)
layout_ncol <- get_cfg(cfg, "layout", "ncol", default = NULL)

stat_method <- get_cfg(cfg, "stats", "method", default = "pearson")
transform_y_mode <- get_cfg(cfg, "transform_y", "mode", default = "none")

width_cfg  <- get_cfg(cfg, "plot_size", "width",  default = 8)
height_cfg <- get_cfg(cfg, "plot_size", "height", default = 6)

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

## 如有 axes_to_plot 限定，就在这里过滤
if (!is.null(axes_to_plot)) {
  axis_combined <- axis_combined %>%
    filter(axis %in% axes_to_plot)
}

## ---------------- 4. 读取脂质指标 YAML，构造所有 axis-phenotype 对 ---------------------

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
  # 构造 axis_indicator_pairs_all: data.frame with columns axis, phenotype
  axis_indicator_pairs_all <- purrr::imap_dfr(mech_axes, function(ax_info, ax_name) {
    inds <- ax_info$indicators
    if (is.null(inds)) {
      return(tibble(axis = character(0), phenotype = character(0)))
    }
    tibble(
      axis = ax_name,
      phenotype = purrr::map_chr(inds, "name")
    )
  })
  cat("[INFO] 从 lipid_indicator_yaml 中获得 axis-phenotype 对数: ", nrow(axis_indicator_pairs_all), "\n", sep = "")
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

## 如配置要求，先输出“全轴全指标”补充图（使用原始 diff_HO_WT 作为 y）
if (!is.null(out_full_pdf) || !is.null(out_full_png)) {
  p_full <- ggplot(full_plot_data, aes(x = combined_z, y = diff_HO_WT)) +
    geom_point() +
    geom_smooth(
      method = "lm",
      se = TRUE,
      aes(color = dir_match),
      show.legend = TRUE
    ) +
    facet_grid(axis ~ phenotype, scales = "free") +
    geom_text(
      data = stats_df_all,
      aes(
        x = -Inf, y = Inf,
        label = label
      ),
      hjust = -0.1, vjust = 1.1,
      inherit.aes = FALSE,
      size = 3
    ) +
    scale_color_manual(
      values = c(
        match    = "red",
        mismatch = "grey40",
        unknown  = "black"
      )
    ) +
    labs(
      x = "Multi-omics mechanistic axis effect (HO vs WT, combined Z)",
      y = "CL phenotype difference (HO - WT)",
      color = "Direction match"
    ) +
    theme_bw() +
    theme(
      strip.background = element_rect(fill = "grey90"),
      panel.grid = element_line(size = 0.2),
      legend.position = "bottom"
    )

  if (!is.null(out_full_pdf)) {
    dir.create(dirname(out_full_pdf), showWarnings = FALSE, recursive = TRUE)
    ggsave(out_full_pdf, p_full, width = width_cfg, height = height_cfg, units = "in")
  }
  if (!is.null(out_full_png)) {
    dir.create(dirname(out_full_png), showWarnings = FALSE, recursive = TRUE)
    ggsave(out_full_png, p_full, width = width_cfg, height = height_cfg, units = "in", dpi = 300)
  }
}

## 设置 factor 顺序（可按需要调整）
plot_data <- plot_data %>%
  mutate(
    axis      = factor(axis, levels = axes_to_plot %||% unique(axis)),
    phenotype = factor(phenotype)
  )

p <- ggplot(plot_data, aes(x = combined_z, y = .data[[y_var]])) +
  geom_point() +
  ## 用 dir_match 映射颜色（match/mismatch/unknown），可以在 scale_color_manual 里自定义
  geom_smooth(
    method = "lm",
    se = TRUE,
    aes(color = dir_match),
    show.legend = TRUE
  ) +
  facet_grid(axis ~ phenotype, scales = if (identical(transform_y_mode, "z_by_indicator")) "free_x" else "free") +
  ## 在每个 panel 右上角写 r, p
  geom_text(
    data = stats_df_plot,
    aes(
      x = -Inf, y = Inf,
      label = label
    ),
    hjust = -0.1, vjust = 1.1,
    inherit.aes = FALSE,
    size = 3
  ) +
  scale_color_manual(
    values = c(
      match    = "red",
      mismatch = "grey40",
      unknown  = "black"
    )
  ) +
  labs(
    x = "Multi-omics mechanistic axis effect (HO vs WT, combined Z)",
    y = y_lab,
    color = "Direction match"
  ) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "grey90"),
    panel.grid = element_line(size = 0.2),
    legend.position = "bottom"
  )

## ---------------- 12. 保存图像和统计结果 -------------------

ggsave(out_pdf, p, width = width_cfg, height = height_cfg, units = "in")
ggsave(out_png, p, width = width_cfg, height = height_cfg, units = "in", dpi = 300)

if (!is.null(out_stats_tsv)) {
  dir.create(dirname(out_stats_tsv), showWarnings = FALSE, recursive = TRUE)
  readr::write_tsv(stats_df_all, out_stats_tsv)
  cat("[INFO] 全部 axis-phenotype 统计结果已保存至: ", out_stats_tsv, "\n", sep = "")
}

cat("[OK] Figure 6F 已输出到:\n  - ", out_pdf, "\n  - ", out_png, "\n", sep = "")