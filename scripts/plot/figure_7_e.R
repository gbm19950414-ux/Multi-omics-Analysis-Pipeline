#!/usr/bin/env Rscript

# Figure 7E: 可视化 EphB1 → membrane context 轴的总效应 c 与直接效应 c'
# 输入: results/multiomics/tables/mediation/mediation_figure7D_membrane_axis_eph_downstream.tsv
# 输出: results/figs/figure_7_e.pdf

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(forcats)
  library(tidyr)
  library(ggplot2)
})

# ---------------- 1. 路径配置 ----------------
in_tsv  <- "results/multiomics/tables/mediation/mediation_figure7D_membrane_axis_eph_downstream.tsv"
out_pdf <- "results/figs/figure_7_e.pdf"

message("[Figure 7E] 读取中介结果: ", in_tsv)
dat <- readr::read_tsv(in_tsv, show_col_types = FALSE)

# 如果文件里已经是 membrane context 子集，这里 filter 只是保险
dat <- dat %>%
  filter(
    !is.na(c_total),
    !is.na(c_prime)
  )

# ---------------- 2. 整理通路名称 & 长表 ----------------
dat_long <- dat %>%
  mutate(
    # 清理通路名字，去掉统一前缀并把下划线改成空格
    pathway = mediator_col %>%
      str_replace("^EphB1_", "") %>%
      str_replace_all("_", " ")
  ) %>%
  # 按 c_total 排序，越负的越靠上
  mutate(pathway = fct_reorder(pathway, c_total)) %>%
  # 把 c / c' 展开为长表
  transmute(
    pathway,
    c_total,  se_c_total,
    c_prime,  se_c_prime
  ) %>%
  pivot_longer(
    cols = c(c_total, c_prime),
    names_to = "effect_type",
    values_to = "estimate"
  ) %>%
  mutate(
    se = dplyr::case_when(
      effect_type == "c_total" ~ se_c_total,
      effect_type == "c_prime" ~ se_c_prime,
      TRUE ~ NA_real_
    ),
    effect_type = dplyr::recode(
      effect_type,
      c_total = "Total effect (c)",
      c_prime = "Direct effect (c')"
    ),
    effect_type = factor(
      effect_type,
      levels = c("Total effect (c)", "Direct effect (c')")
    ),
    ci_low  = estimate - 1.96 * se,
    ci_high = estimate + 1.96 * se
  )

# ---------------- 3. 画图 ----------------
message("[Figure 7E] 绘图并导出: ", out_pdf)

pd <- position_dodge(width = 0.6)

p <- ggplot(dat_long, aes(x = estimate, y = pathway, color = effect_type)) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.3) +
  geom_errorbarh(
    aes(xmin = ci_low, xmax = ci_high),
    height = 0.25,
    position = pd
  ) +
  geom_point(
    position = pd,
    size = 2
  ) +
  labs(
    x = "Effect of EphB1 KO on membrane context axis\n(coefficient for X in Y~X 或 Y~X+M)",
    y = "EphB1 downstream pathway (mediator)",
    color = NULL,
    title = "Figure 7E – Total effect (c) vs direct effect (c')",
    subtitle = "Mediation of EphB1 loss → Membrane context axis via downstream pathways"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "top"
  )

ggsave(out_pdf, p, width = 7, height = 4.5, device = cairo_pdf)

message("[Figure 7E] 完成。")