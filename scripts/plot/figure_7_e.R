#!/usr/bin/env Rscript

# scripts/plot/figure_7_e.R
# 目的：从 mediation_summary.tsv 中提取 EphB1 下游通路 → 脂质指标 的中介结果，
#      生成 Figure 7E：脂质层面的 a×b 间接效应森林图
#
# 输入（固定路径，可按需改成命令行参数）：
#   results/multiomics/tables/mediation/mediation_summary.tsv
# 输出：
#   results/figs/figure_7_e.pdf

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(forcats)
  library(ggplot2)
  library(tidyr)
})

# -------- 1. 路径设置 --------
in_file  <- "results/multiomics/tables/mediation/mediation_summary.tsv"
out_pdf  <- "results/figs/figure_7_e.pdf"

message("[Figure 7E] 读取中介汇总表: ", in_file)
med <- readr::read_tsv(in_file, show_col_types = FALSE)

# -------- 2. 选择我们要的 mediator / outcome --------
# 2.1 限定 EphB1 下游通路
if ("M_type" %in% names(med)) {
  # 优先使用 06.1 脚本写入的 M_type 标签
  downstream_tag <- "M1_ephb1_downstream_or_other"
  message("[Figure 7E] 使用 M_type = ", downstream_tag, " 筛选 EphB1 下游通路。")

  med2 <- med %>%
    filter(M_type == downstream_tag)

  if (nrow(med2) == 0) {
    stop("[ERROR] 过滤 M_type = ", downstream_tag, " 后没有行，请检查 mediation_summary.tsv。")
  }
} else {
  # mediation_summary.tsv 没有 M_type 列时，退回用 mediator_col 前缀筛选
  if (!"mediator_col" %in% names(med)) {
    stop("[ERROR] mediation_summary.tsv 中既没有 M_type 也没有 mediator_col，无法识别 EphB1 下游通路。")
  }

  message("[Figure 7E] mediation_summary.tsv 中缺少 M_type 列，",
          "退回使用 mediator_col 以 'EphB1_' 前缀筛选下游通路。")

  med2 <- med %>%
    filter(str_detect(mediator_col, "^EphB1_"))

  if (nrow(med2) == 0) {
    stop("[ERROR] 使用 mediator_col 以 'EphB1_' 前缀筛选后没有行，请检查 mediation_summary.tsv。")
  }
}

# 2.2 关心的脂质 / 指标 Y（在 05 脚本中已经 rename 为 *_score）
y_candidates <- c(
  "sum_CL_score",
  "CL_fraction_score",
  "mito_lipid_mass_score",
  "PC_score",
  "PE_score",
  "CL_score"
)

y_available <- intersect(y_candidates, unique(med2$outcome_col))
if (length(y_available) == 0) {
  stop("[ERROR] 在 mediation_summary.tsv 中找不到任何指定的 outcome_col: ",
       paste(y_candidates, collapse = ", "))
}

message("[Figure 7E] 使用的脂质/指标 outcome : ",
        paste(y_available, collapse = ", "))

med3 <- med2 %>%
  filter(outcome_col %in% y_available) %>%
  # 必须有 bootstrap 间接效应
  filter(!is.na(ab_boot_mean))

if (nrow(med3) == 0) {
  stop("[ERROR] 对指定 outcome 和 Eph 下游通路来说，没有 ab_boot_mean 非 NA 的行。")
}

# -------- 3. 整理名字 & 方向标签 --------

# 3.1 把 mediator 名字弄漂亮一点
med3 <- med3 %>%
  mutate(
    pathway = mediator_col %>%
      str_replace("^EphB1_", "") %>%
      str_replace_all("_", " ")
  )

# 3.2 outcome 的漂亮标签
pretty_y <- c(
  sum_CL_score          = "sum CL (总 CL)",
  CL_fraction_score     = "CL fraction (CL/total)",
  mito_lipid_mass_score = "Mito lipid mass (PC+PE+CL)",
  PC_score              = "PC (membrane PC)",
  PE_score              = "PE (membrane PE)",
  CL_score              = "CL (membrane CL)"
)

med3 <- med3 %>%
  mutate(
    outcome_pretty = dplyr::recode(outcome_col, !!!pretty_y,
                                   .default = outcome_col)
  )

# 3.3 根据 bootstrap a/b 做方向标签
if (!all(c("a_boot_mean", "b_boot_mean") %in% names(med3))) {
  warning("[WARN] mediation_summary.tsv 中没有 a_boot_mean / b_boot_mean，",
          "方向标签将退回使用 a / b 的点估计。")
  med3 <- med3 %>%
    mutate(
      a_use = a,
      b_use = b
    )
} else {
  med3 <- med3 %>%
    mutate(
      a_use = a_boot_mean,
      b_use = b_boot_mean
    )
}

med3 <- med3 %>%
  mutate(
    label_ab = case_when(
      a_use > 0 & b_use > 0 ~ "a+, b+",
      a_use > 0 & b_use < 0 ~ "a+, b−",
      a_use < 0 & b_use > 0 ~ "a−, b+",
      a_use < 0 & b_use < 0 ~ "a−, b−",
      TRUE                  ~ "a?, b?"
    )
  )

# 3.4 计算 ab 的排序（每个 outcome 内按 |ab| 排，越大越靠上）
med3 <- med3 %>%
  group_by(outcome_pretty) %>%
  arrange(desc(abs(ab_boot_mean)), .by_group = TRUE) %>%
  mutate(
    pathway = fct_reorder(pathway, ab_boot_mean),
    # 给每个 outcome 内加一个 rank，后续想筛 Top N 时可以用
    rank_within_outcome = row_number()
  ) %>%
  ungroup()

# 如你想只看每个 outcome Top N（比如 6 条），可以把下面注释打开：
# top_n_per_outcome <- 6
# med3 <- med3 %>%
#   filter(rank_within_outcome <= top_n_per_outcome)

# -------- 3.x 额外输出：每个 outcome 单独成图 --------
out_dir_single <- "results/figs/figure_7_e_single"
if (!dir.exists(out_dir_single)) dir.create(out_dir_single, recursive = TRUE)

unique_outcomes <- unique(med3$outcome_pretty)

for (yy in unique_outcomes) {
  df_y <- med3 %>% dplyr::filter(outcome_pretty == yy)

  x_min_y <- min(df_y$ab_boot_ci_low,  na.rm = TRUE)
  x_max_y <- max(df_y$ab_boot_ci_high, na.rm = TRUE)
  pad_y   <- 0.05 * (x_max_y - x_min_y)

  p_y <- ggplot(df_y,
                aes(x = ab_boot_mean,
                    y = pathway,
                    xmin = ab_boot_ci_low,
                    xmax = ab_boot_ci_high,
                    color = label_ab)) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_errorbarh(height = 0) +
    geom_point(size = 2) +
    scale_x_continuous(limits = c(x_min_y - pad_y, x_max_y + pad_y)) +
    scale_color_discrete(drop = FALSE) +
    labs(
      x = "Indirect effect (a × b)",
      y = "EphB1 downstream pathway",
      color = "Direction\n(a, b)",
      title = paste0("Figure 7E – ", yy),
      subtitle = "EphB1 (KO vs WT) → downstream pathway (M) → lipid / indicator (Y)"
    ) +
    theme_bw(base_size = 9) +
    theme(legend.position = "bottom")

  safe_yy <- yy %>%
    stringr::str_replace_all("[^A-Za-z0-9_\\-]+", "_")

  out_file_y <- file.path(out_dir_single, paste0("figure_7_e_", safe_yy, ".pdf"))
  message("[Figure 7E] 输出单独图像: ", out_file_y)
  ggsave(out_file_y, p_y, width = 6, height = 4, device = cairo_pdf)
}

# -------- 4. 画图：多 outcome facet 的森林图 --------

# 确定 X 轴范围，让 0 在中间稍微好看
x_min <- min(med3$ab_boot_ci_low,  na.rm = TRUE)
x_max <- max(med3$ab_boot_ci_high, na.rm = TRUE)
pad   <- 0.05 * (x_max - x_min)

p <- ggplot(med3,
            aes(x = ab_boot_mean,
                y = pathway,
                xmin = ab_boot_ci_low,
                xmax = ab_boot_ci_high,
                color = label_ab)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_errorbarh(height = 0) +
  geom_point(size = 2) +
  scale_x_continuous(
    limits = c(x_min - pad, x_max + pad)
  ) +
  scale_color_discrete(drop = FALSE) +
  facet_wrap(~ outcome_pretty, ncol = 2, scales = "free_y") +
  labs(
    x = "Indirect effect (a × b), bootstrap mean ± 95% CI",
    y = "EphB1 downstream pathway",
    color = "Direction\n(a, b)",
    title = "Figure 7E – Mediation of lipid readouts by EphB1 downstream pathways",
    subtitle = "EphB1 (KO vs WT) → downstream pathway (M) → lipid / indicator (Y)"
  ) +
  theme_bw(base_size = 9) +
  theme(
    strip.background = element_rect(colour = NA),
    strip.text       = element_text(face = "bold"),
    legend.position  = "bottom"
  )

message("[Figure 7E] 导出图像: ", out_pdf)
ggsave(out_pdf, p, width = 8, height = 6, device = cairo_pdf)

message("[Figure 7E] 完成。")