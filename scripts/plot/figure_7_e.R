#!/usr/bin/env Rscript

## ============================================================
## Figure 7E: Mediation effects on lipid/indicator readouts (small multiples)
##   EphB1 loss (X) -> EphB1 downstream pathway score (M) -> lipid/indicator (Y)
##
## Input:
##   results/multiomics/tables/mediation/mediation_summary.tsv
##
## Output (recommended for Illustrator / Nature-style multi-panel figures):
##   - Per-outcome single-column panels (84 mm width):
##     results/figs/figure_7_e_single/figure_7_e_<outcome>.pdf
##   - Optional stacked overview (84 mm width):
##     results/figs/figure_7_e_mediation_lipid_readouts_stack.pdf
##
## Plot intent (logic):
##   - Effect size + uncertainty first: ab_boot_mean ± 95% bootstrap CI
##   - Within each outcome Y, rank pathways by |ab| to surface top mediators
##   - De-emphasize non-top pathways visually; highlight PI3K and Cytoskeleton/FAK/Src/Rho
##   - Avoid color semantics that imply “activation/inhibition”; direction shown as a+/b± label
##
## Figure style:
##   - Read Nature-like style from:
##     /Volumes/Samsung_SSD_990_PRO_2TB_Media/EphB1/02_protocols/figure_style_nature.yaml
##   - Pathway display names are unified with Figure 7C via:
##     /Volumes/Samsung_SSD_990_PRO_2TB_Media/multiomics_mech/scripts/plot/figure_7_c.yaml
##   - Save vector PDF; single-column width 84 mm
## ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(forcats)
  library(ggplot2)
  library(tidyr)
})

`%||%` <- function(x, y) if (!is.null(x)) x else y
pt_to_mm <- function(pt) pt * 0.3527778

## -------------------- Paths ----------------------------------

in_file <- "results/multiomics/tables/mediation/mediation_summary.tsv"
out_dir_single <- "results/figs/figure_7_e_single"
out_stack <- "results/figs/figure_7_e_mediation_lipid_readouts_stack.pdf"

dir.create(dirname(out_stack), recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(out_dir_single)) dir.create(out_dir_single, recursive = TRUE)

message("[Figure 7E] 读取中介汇总表: ", in_file)
med <- readr::read_tsv(in_file, show_col_types = FALSE)

## -------------------- Load style (Nature-like) ----------------

style_file <- "/Volumes/Samsung_SSD_990_PRO_2TB_Media/EphB1/02_protocols/figure_style_nature.yaml"
style <- NULL
if (file.exists(style_file) && requireNamespace("yaml", quietly = TRUE)) {
  style <- yaml::read_yaml(style_file)
} else {
  message("[WARN] 未能读取样式文件或缺少 yaml 包，将使用内置默认样式。")
}

font_family   <- style$typography$font_family_primary %||% "Helvetica"
size_tick_pt  <- style$typography$sizes_pt$axis_tick_default %||% 5.5
size_label_pt <- style$typography$sizes_pt$axis_label_default %||% 6.5

axis_line_pt  <- style$lines$axis_line_default_pt %||% 0.25
errorbar_pt   <- style$lines$errorbar_default_pt %||% 0.25
line_pt       <- style$lines$line_width_pt %||% 0.5

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
  message("[WARN] 未能读取 figure_7_c.yaml 或缺少 yaml 包，将使用原始通路名称。")
}

## -------------------- Select mediator / outcomes --------------

# 1) Limit to EphB1 downstream pathways
if ("M_type" %in% names(med)) {
  downstream_tag <- "M1_ephb1_downstream_or_other"
  message("[Figure 7E] 使用 M_type = ", downstream_tag, " 筛选 EphB1 下游通路。")
  med2 <- med %>% filter(M_type == downstream_tag)
  if (nrow(med2) == 0) stop("[ERROR] 过滤 M_type = ", downstream_tag, " 后没有行，请检查 mediation_summary.tsv。")
} else {
  if (!"mediator_col" %in% names(med)) stop("[ERROR] mediation_summary.tsv 中既没有 M_type 也没有 mediator_col。")
  message("[Figure 7E] mediation_summary.tsv 中缺少 M_type 列，退回使用 mediator_col 以 'EphB1_' 前缀筛选。")
  med2 <- med %>% filter(str_detect(mediator_col, "^EphB1_"))
  if (nrow(med2) == 0) stop("[ERROR] 使用 mediator_col 以 'EphB1_' 前缀筛选后没有行，请检查 mediation_summary.tsv。")
}

# 2) Select lipid/indicator outcomes (scored variables)
y_candidates <- c(
  "sum_CL_score",
  "CL_fraction_score",
  "mito_lipid_mass_score",
  "PC_score",
  "PE_score",
  "CL_score"
)

if (!"outcome_col" %in% names(med2)) stop("[ERROR] mediation_summary.tsv 缺少 outcome_col 列。")

y_available <- intersect(y_candidates, unique(med2$outcome_col))
if (length(y_available) == 0) {
  stop("[ERROR] 在 mediation_summary.tsv 中找不到任何指定的 outcome_col: ", paste(y_candidates, collapse = ", "))
}
message("[Figure 7E] 使用的脂质/指标 outcome: ", paste(y_available, collapse = ", "))

# 3) Must have bootstrap indirect effect
needed_cols <- c("mediator_col", "outcome_col", "ab_boot_mean", "ab_boot_ci_low", "ab_boot_ci_high")
missing_cols <- setdiff(needed_cols, names(med2))
if (length(missing_cols) > 0) stop("[ERROR] mediation_summary.tsv 缺少列: ", paste(missing_cols, collapse = ", "))

med3 <- med2 %>%
  filter(outcome_col %in% y_available) %>%
  filter(!is.na(ab_boot_mean))

if (nrow(med3) == 0) stop("[ERROR] 对指定 outcome 和 Eph 下游通路来说，没有 ab_boot_mean 非 NA 的行。")

## -------------------- Pretty labels ---------------------------

# Outcome pretty labels (panel headers)
pretty_y <- c(
  sum_CL_score          = "sum CL",
  CL_fraction_score     = "CL fraction (CL/total)",
  mito_lipid_mass_score = "Mito lipid mass (PC+PE+CL)",
  PC_score              = "Membrane PC",
  PE_score              = "Membrane PE",
  CL_score              = "Membrane CL"
)

# Desired outcome order for stacked overview (use pretty labels)
outcome_order_pretty <- unname(pretty_y[intersect(y_candidates, names(pretty_y))])

# Direction labels: prefer bootstrap a/b means if available
use_boot_ab <- all(c("a_boot_mean", "b_boot_mean") %in% names(med3))
if (!use_boot_ab) {
  message("[WARN] mediation_summary.tsv 中缺少 a_boot_mean / b_boot_mean，方向标签将退回使用 a / b 点估计。")
}

med3 <- med3 %>%
  mutate(
    a_use = if (use_boot_ab) a_boot_mean else a,
    b_use = if (use_boot_ab) b_boot_mean else b,
    direction = case_when(
      a_use > 0 & b_use > 0 ~ "a+/b+",
      a_use > 0 & b_use < 0 ~ "a+/b−",
      a_use < 0 & b_use > 0 ~ "a−/b+",
      a_use < 0 & b_use < 0 ~ "a−/b−",
      TRUE                  ~ ""
    ),
    outcome_pretty = dplyr::recode(outcome_col, !!!pretty_y, .default = outcome_col),

    # Pathway display name (unified with Figure 7C)
    mediator_raw = mediator_col,
    pathway = {
      if (!is.null(pathway_rename)) {
        mapped <- unname(pathway_rename[mediator_col])
        dplyr::coalesce(mapped, mediator_col)
      } else {
        mediator_col
      }
    },

    highlight_group = dplyr::case_when(
      str_detect(pathway, regex("PI3K", ignore_case = TRUE)) ~ "PI3K",
      str_detect(pathway, regex("cytoskeleton|FAK|Src|Rho", ignore_case = TRUE)) ~ "Cytoskeleton/FAK/Src/Rho",
      TRUE ~ "Other"
    )
  )

## -------------------- Reorder within each outcome -------------

# To allow per-facet ordering without extra packages, create a facet-specific y key
# and set levels in the desired order.
med3 <- med3 %>%
  group_by(outcome_pretty) %>%
  arrange(desc(abs(ab_boot_mean)), .by_group = TRUE) %>%
  mutate(
    pathway_key = paste0(pathway, "___", outcome_pretty),
    rank_within_outcome = row_number()
  ) %>%
  ungroup()

# Global factor levels preserve per-facet ordering because each facet uses a subset
med3 <- med3 %>%
  mutate(
    pathway_key = factor(pathway_key, levels = unique(pathway_key))
  )

# Direction label positions computed per outcome to keep consistent padding
med3 <- med3 %>%
  group_by(outcome_pretty) %>%
  mutate(
    x_min = min(ab_boot_ci_low,  ab_boot_mean, na.rm = TRUE),
    x_max = max(ab_boot_ci_high, ab_boot_mean, na.rm = TRUE),
    x_span = max(1e-8, x_max - x_min),
    label_x = ab_boot_ci_high + 0.04 * x_span,
    xlim_left  = x_min - 0.02 * x_span,
    xlim_right = max(label_x, na.rm = TRUE) + 0.02 * x_span
  ) %>%
  ungroup()

## -------------------- Plot function (single outcome) ----------

plot_one_outcome <- function(df_y) {
  # Use per-outcome x limits
  xlim_left  <- df_y$xlim_left[1]
  xlim_right <- df_y$xlim_right[1]

  ggplot(df_y, aes(x = ab_boot_mean, y = pathway_key)) +
    geom_errorbarh(
      aes(xmin = ab_boot_ci_low, xmax = ab_boot_ci_high),
      height   = 0.28,
      linewidth = pt_to_mm(errorbar_pt),
      color    = "grey50"
    ) +
    geom_point(
      aes(shape = ab_boot_mean > 0, color = highlight_group),
      size   = point_size_mm,
      stroke = pt_to_mm(line_pt)
    ) +
    geom_vline(
      xintercept = 0,
      linetype   = "dashed",
      linewidth  = pt_to_mm(axis_line_pt),
      color      = "grey40"
    ) +
    geom_text(
      aes(label = direction, x = label_x),
      hjust  = 0,
      size   = size_tick_pt / ggplot2::.pt,
      family = font_family,
      color  = "grey20"
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
    scale_y_discrete(labels = function(x) str_replace(x, "___.*$", "")) +
    coord_cartesian(xlim = c(xlim_left, xlim_right), clip = "off") +
    labs(x = "Indirect effect (a × b), bootstrap mean ± 95% CI", y = NULL) +
    theme_classic(base_size = size_label_pt, base_family = font_family) +
    theme(
      axis.text.y  = element_text(size = size_tick_pt),
      axis.text.x  = element_text(size = size_tick_pt),
      axis.title.x = element_text(
        size   = size_label_pt,
        margin = margin(t = axis_title_margin_y_pt, unit = "pt")
      ),
      axis.line  = element_line(linewidth = pt_to_mm(axis_line_pt), color = "black"),
      axis.ticks = element_line(linewidth = pt_to_mm(axis_line_pt), color = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(t = 0, r = 6, b = 0, l = 0, unit = "mm")
    )
}

## -------------------- Output per-outcome panels ---------------

width_mm <- 84

for (yy in unique(med3$outcome_pretty)) {
  df_y <- med3 %>% filter(outcome_pretty == yy)

  # Height scales with number of pathways in this outcome
  n_pathway <- nrow(df_y)
  height_mm <- max(45, min(120, 12 + n_pathway * 4))

  p_y <- plot_one_outcome(df_y)

  safe_yy <- yy %>%
    stringr::str_replace_all("[^A-Za-z0-9_\\-]+", "_") %>%
    stringr::str_replace_all("_+", "_")

  out_file_y <- file.path(out_dir_single, paste0("figure_7_e_", safe_yy, ".pdf"))
  message("[Figure 7E] 输出单独图像: ", out_file_y)
  ggsave(
    filename = out_file_y,
    plot     = p_y,
    width    = width_mm,
    height   = height_mm,
    units    = "mm",
    device   = cairo_pdf
  )
}

## -------------------- Optional: stacked overview --------------

# This is for quick inspection; for final layout, use per-outcome PDFs above.
# Stack facets vertically to keep 84 mm width.

med3_stack <- med3 %>%
  mutate(
    outcome_pretty = factor(
      outcome_pretty,
      levels = outcome_order_pretty[outcome_order_pretty %in% unique(outcome_pretty)]
    )
  )

# Compute global height (approx.)
facet_n <- length(unique(med3$outcome_pretty))
stack_height_mm <- max(120, min(280, 30 + facet_n * 40))

p_stack <- ggplot(med3_stack, aes(x = ab_boot_mean, y = pathway_key)) +
  geom_errorbarh(
    aes(xmin = ab_boot_ci_low, xmax = ab_boot_ci_high),
    height   = 0.25,
    linewidth = pt_to_mm(errorbar_pt),
    color    = "grey50"
  ) +
  geom_point(
    aes(shape = ab_boot_mean > 0, color = highlight_group),
    size   = point_size_mm,
    stroke = pt_to_mm(line_pt)
  ) +
  geom_vline(
    xintercept = 0,
    linetype   = "dashed",
    linewidth  = pt_to_mm(axis_line_pt),
    color      = "grey40"
  ) +
  geom_text(
    aes(label = direction, x = label_x),
    hjust  = 0,
    size   = size_tick_pt / ggplot2::.pt,
    family = font_family,
    color  = "grey20"
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
  scale_y_discrete(labels = function(x) str_replace(x, "___.*$", "")) +
  facet_wrap(~ outcome_pretty, ncol = 1, scales = "free_y") +
  labs(x = "Indirect effect (a × b), bootstrap mean ± 95% CI", y = NULL) +
  theme_classic(base_size = size_label_pt, base_family = font_family) +
  theme(
    axis.text.y  = element_text(size = size_tick_pt),
    axis.text.x  = element_text(size = size_tick_pt),
    axis.title.x = element_text(
      size   = size_label_pt,
      margin = margin(t = axis_title_margin_y_pt, unit = "pt")
    ),
    axis.line  = element_line(linewidth = pt_to_mm(axis_line_pt), color = "black"),
    axis.ticks = element_line(linewidth = pt_to_mm(axis_line_pt), color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text       = element_text(size = size_label_pt, face = "bold"),
    plot.margin = margin(t = 0, r = 6, b = 0, l = 0, unit = "mm")
  )

message("[Figure 7E] 输出叠加预览图: ", out_stack)
ggsave(
  filename = out_stack,
  plot     = p_stack,
  width    = width_mm,
  height   = stack_height_mm,
  units    = "mm",
  device   = cairo_pdf
)

message("[Figure 7E] 完成。")