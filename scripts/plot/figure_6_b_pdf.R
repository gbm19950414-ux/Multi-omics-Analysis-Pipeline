#!/usr/bin/env Rscript

# -------------
# Figure 6B: summary table -> PDF (for Illustrator)
# 读取 results/figs/figure_6_b_summary.xlsx
# 输出 results/figs/figure_6_b.pdf
# -------------

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(gt)
})

# 工作目录假定是项目根目录（你当前就是在 multiomics_mech 下跑 Rscript）
proj_dir <- getwd()

input_xlsx  <- file.path(proj_dir, "results", "figs", "figure_6_b_summary.xlsx")
output_pdf  <- file.path(proj_dir, "results", "figs", "figure_6_b.pdf")

message("Input:  ", input_xlsx)
message("Output: ", output_pdf)

# 1. 读入精简表格 -------------------------------------------------------------
df_raw <- readxl::read_xlsx(input_xlsx)

# ！！！关键一步：根据实际列名修改这里！！！
# 当前 summary 表的实际英文列名为：
#   axis               机制轴（例如：Synthesis, Oxidation, Transport ...）
#   feature_type       特征类型（此处为 lipid，可不在图中展示）
#   name               指标名（例如：CL_PG_ratio, oxCL_CL_ratio ...）
#   description        指标的中文解释（简要功能说明）
#   raw_increase_means Ephb1-/- 中原始“升高”含义（中文解释，可视为“预期变化”）
#   effect_on_CL       对 CL 的影响（此处先不放入 figure 表格）
#   bottleneck_sign    瓶颈方向（内部计算用，不放入 figure 表格）
#   weight             权重（内部计算用，不放入 figure 表格）
#   note               额外备注（内部笔记，用于描述为何选为“核心指标”等）
#
# 在 figure 6B 的展示中，我们只挑选 4 列：
#   机制轴、指标、预期变化、解释
# 其中：
#   机制轴   <- axis
#   指标     <- name
#   预期变化 <- raw_increase_means
#   解释     <- description

df <- df_raw %>%
  dplyr::rename(
    机制轴   = axis,               # 机制轴（Synthesis, Oxidation, Transport ...）
    指标     = name,               # 指标名（CL_PG_ratio, oxCL_CL_ratio ...）
    预期变化 = raw_increase_means, # Ephb1-/- 中预期变化（原始“升高/变化”解释）
    解释     = description         # 简要功能说明
  )

# 确保行顺序是你想要的（例如手动在 Excel 里排好，这里就不改）
# 如果要显式设置顺序，可以用 mutate + factor，这里先省略。

# 2. 用 gt 生成紧凑表格 -------------------------------------------------------

tab <- df %>%
  gt::gt() %>%
  # 标题可选；如果你准备在 Illustrator 里自己写标题，可以去掉这一段
  gt::tab_header(
    title = gt::md("**Figure 6B. 脂质组学机制指标摘要**")
  ) %>%
  # 列标题加粗
  gt::tab_style(
    style = list(gt::cell_text(weight = "bold")),
    locations = gt::cells_column_labels(gt::everything())
  ) %>%
  # 统一排版选项：字体、字号、边距等（精简兼容版）
  gt::tab_options(
    table.font.names = "Arial",  # 和其他 panel 保持一致
    table.font.size  = 7,        # 适合放进整版 figure
    data_row.padding = gt::px(1) # 稍微紧凑一点的行距
  ) %>%
  # 列宽：这里给一个初始，之后可视效果再调
  # 如果列太多或内容太长，可以适当调小某些列宽
  gt::cols_width(
    gt::everything() ~ gt::px(90)
  ) %>%
  # 对齐方式：指标左对齐，解释左对齐，其他居中（可按需要调整）
  gt::cols_align(
    align = "left",
    columns = c(指标, 解释)
  ) %>%
  gt::cols_align(
    align = "center",
    columns = c(机制轴, 预期变化)
  )

# 3. 保存为 PDF ---------------------------------------------------------------

# 注意：保存为 PDF 通常不需要 webshot2；如果你改成输出 PNG/JPEG，
# gt 会走 webshot2 路径，又会报 "webshot2 未安装"。
gt::gtsave(
  data = tab,
  filename = output_pdf
)

message("Done: ", output_pdf)