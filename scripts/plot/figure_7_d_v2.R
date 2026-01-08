#!/usr/bin/env Rscript

## ============================================================
## Figure 7D: Forest plot of mediation effects
##   EphB1 loss (X) -> EphB1 downstream pathway score (M) -> Membrane context axis (Y)
##
## Input:
##   results/multiomics/tables/mediation/mediation_figure7D_membrane_axis_eph_downstream.tsv
##
## Output:
##   results/figs/figure_7_d_mediation_membrane_context.pdf
##
## Plot intent (logic):
##   - Put effect size + uncertainty first: ab_boot_mean ± 95% bootstrap CI
##   - Rank pathways by |ab| to surface top mediators (PI3K–AKT–mTOR typically ranks first)
##   - De-emphasize non-top pathways visually, while still showing the full candidate set
##   - Direction label summarizes signs of a and b (derived from bootstrap means)
##
## Figure style:
##   - Read Nature-like style from:
##     /Volumes/Samsung_SSD_990_PRO_2TB_Media/EphB1/02_protocols/figure_style_nature.yaml
##   - Save as single-column width (84 mm) vector PDF for Illustrator
## ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(stringr)
})

`%||%` <- function(x, y) if (!is.null(x)) x else y
pt_to_mm <- function(pt) pt * 0.3527778

## -------------------- Input / Output -------------------------

infile  <- "results/multiomics/tables/mediation/mediation_figure7D_membrane_axis_eph_downstream.tsv"
outfile <- "results/figs/figure_7_d_mediation_membrane_context.pdf"

dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)

cat("[INFO] 读取输入文件: ", infile, "\n")
df <- readr::read_tsv(infile, show_col_types = FALSE)

stopifnot(all(c(
  "mediator_col",
  "ab_boot_mean", "ab_boot_ci_low", "ab_boot_ci_high",
  "a_boot_mean",  "b_boot_mean"
) %in% colnames(df)))

## -------------------- Load style (Nature-like) ----------------

style_file <- "/Volumes/Samsung_SSD_990_PRO_2TB_Media/EphB1/02_protocols/figure_style_nature.yaml"
style <- NULL
if (file.exists(style_file) && requireNamespace("yaml", quietly = TRUE)) {
  style <- yaml::read_yaml(style_file)
} else {
  cat("[WARN] 未能读取样式文件或缺少 yaml 包，将使用内置默认样式。\n")
}

font_family   <- style$typography$font_family_primary %||% "Helvetica"
size_tick_pt  <- style$typography$sizes_pt$axis_tick_default %||% 5.5
size_label_pt <- style$typography$sizes_pt$axis_label_default %||% 6.5
size_legend_pt <- style$typography$sizes_pt$legend_text_default %||% 6

axis_line_pt  <- style$lines$axis_line_default_pt %||% 0.25
errorbar_pt   <- style$lines$errorbar_default_pt %||% 0.25
line_pt       <- style$lines$line_width_pt %||% 0.5

axis_title_margin_x_pt <- style$layout$axis_title_margin_pt$x %||% 3
axis_title_margin_y_pt <- style$layout$axis_title_margin_pt$y %||% 3

point_size_mm <- style$marks$point_size %||% 1.6

## -------------------- Load pathway rename map ----------------

rename_yaml <- "/Volumes/Samsung_SSD_990_PRO_2TB_Media/multiomics_mech/scripts/plot/figure_7_c.yaml"
pathway_rename <- NULL
if (file.exists(rename_yaml) && requireNamespace("yaml", quietly = TRUE)) {
  rename_cfg <- yaml::read_yaml(rename_yaml)
  # Ensure a named character vector for safe indexing
  pathway_rename <- unlist(rename_cfg$pathway_rename)
} else {
  cat("[WARN] 未能读取 figure_7_c.yaml 或缺少 yaml 包，将使用原始通路名称。\n")
}

## -------------------- Derive labels / ranking ----------------

# Direction labels are based on bootstrap means of a and b
# (a: X -> M, b: M -> Y | X)

df <- df %>%
  mutate(
    mediator_raw = mediator_col,
    # Vectorized rename: map known IDs to display names; keep original if not mapped
    mediator_col = {
      if (!is.null(pathway_rename)) {
        mapped <- unname(pathway_rename[mediator_col])
        dplyr::coalesce(mapped, mediator_col)
      } else {
        mediator_col
      }
    },
    direction = dplyr::case_when(
      a_boot_mean > 0 & b_boot_mean > 0 ~ "a+/b+",
      a_boot_mean > 0 & b_boot_mean < 0 ~ "a+/b−",
      a_boot_mean < 0 & b_boot_mean > 0 ~ "a−/b+",
      a_boot_mean < 0 & b_boot_mean < 0 ~ "a−/b−",
      TRUE ~ ""
    ),
    highlight_group = dplyr::case_when(
      str_detect(mediator_col, regex("PI3K", ignore_case = TRUE)) ~ "PI3K",
      str_detect(mediator_col, regex("cytoskeleton|FAK|Src|Rho", ignore_case = TRUE)) ~ "Cytoskeleton/FAK/Src/Rho",
      TRUE ~ "Other"
    )
  ) %>%
  arrange(desc(abs(ab_boot_mean))) %>%
  mutate(
    mediator_col = factor(mediator_col, levels = unique(mediator_col))
  )

## -------------------- Layout helpers -------------------------

x_min <- min(df$ab_boot_ci_low,  df$ab_boot_mean, na.rm = TRUE)
x_max <- max(df$ab_boot_ci_high, df$ab_boot_mean, na.rm = TRUE)
x_span <- max(1e-8, x_max - x_min)

# Place direction labels slightly to the right of CI
label_x <- df$ab_boot_ci_high + 0.04 * x_span

# Provide extra room on the right for labels
xlim_left  <- x_min - 0.02 * x_span
xlim_right <- max(label_x, na.rm = TRUE) + 0.02 * x_span

## -------------------- Plot -----------------------------------

p <- ggplot(df, aes(x = ab_boot_mean, y = mediator_col)) +

  # 95% bootstrap CI
  geom_errorbarh(
    aes(xmin = ab_boot_ci_low, xmax = ab_boot_ci_high),
    height   = 0.28,
    linewidth = pt_to_mm(errorbar_pt),
    color    = "grey50"
  ) +

  # Point estimate (shape encodes sign; color de-emphasizes non-top pathways)
  geom_point(
    aes(
      shape = ab_boot_mean > 0,
      color = highlight_group
    ),
    size   = point_size_mm,
    stroke = pt_to_mm(line_pt)
  ) +

  # Zero-effect reference
  geom_vline(
    xintercept = 0,
    linetype   = "dashed",
    linewidth  = pt_to_mm(axis_line_pt),
    color      = "grey40"
  ) +

  # Direction labels (right side)
  geom_text(
    aes(label = direction, x = label_x),
    hjust = 0,
    size  = size_tick_pt / ggplot2::.pt,
    family = font_family,
    color = "grey20"
  ) +

  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), guide = "none") +
  scale_color_manual(
    values = c(
      "PI3K" = "black",
      "Cytoskeleton/FAK/Src/Rho" = "black",
      "Other" = "grey65"
    ),
    guide = "none"
  ) +

  # No title/subtitle inside panel (Nature-style multi-panel figure)
  labs(
    x = "Indirect effect (a × b), bootstrap mean ± 95% CI",
    y = NULL
  ) +

  coord_cartesian(xlim = c(xlim_left, xlim_right), clip = "off") +

  theme_classic(base_size = size_label_pt, base_family = font_family) +
  theme(
    axis.text.y  = element_text(size = size_tick_pt),
    axis.text.x  = element_text(size = size_tick_pt),
    axis.title.x = element_text(
      size = size_label_pt,
      margin = margin(t = axis_title_margin_y_pt, unit = "pt")
    ),

    axis.line  = element_line(linewidth = pt_to_mm(axis_line_pt), color = "black"),
    axis.ticks = element_line(linewidth = pt_to_mm(axis_line_pt), color = "black"),

    # keep clean background; no grids
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),

    # leave minimal right margin for labels (mm)
    plot.margin = margin(t = 0, r = 6, b = 0, l = 0, unit = "mm")
  )

## -------------------- Save output ----------------------------

# Single-column width (84 mm). Height scales with number of pathways.
width_mm <- 84
n_pathway <- nrow(df)
height_mm <- max(45, min(120, 12 + n_pathway * 4))

cat(sprintf("[INFO] 保存 PDF: %s (%.1f mm × %.1f mm)\n", outfile, width_mm, height_mm))
ggsave(
  filename = outfile,
  plot     = p,
  width    = width_mm,
  height   = height_mm,
  units    = "mm",
  device   = cairo_pdf
)

cat("[OK] 输出图像: ", outfile, "\n")