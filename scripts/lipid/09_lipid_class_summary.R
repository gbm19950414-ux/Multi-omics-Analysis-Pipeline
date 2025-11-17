#!/usr/bin/env Rscript

## ============================================================
## 10_lipid_class_summary.R
##
## 功能：
##  1) 把单个 lipid 聚合到 Class 层面，得到：
##     - 每个样本 × 每个 Class 的汇总强度 (sum / mean / log2mean)
##  2) 基于 Class × sample 的 log2(mean_intensity) 矩阵，做：
##     - Class 层面的 HO vs WT 差异分析（limma）
##  3) 出图：
##     - 按 Class 汇总的 WT vs HO 柱状图（均值 ± SEM）
##     - 按 Class 的样本层面箱线图（WT vs HO）
##
## 输入：
##   - results/lipid/tables/lipid_long_dedup_pqn_bb.tsv
##   - results/lipid/tables/sample_info_lipid.tsv
##
## 输出：
##   - results/lipid/tables/lipid_class_summary.tsv
##       每行：Class × global_sample_id（样本）
##       包含：sum_intensity, mean_intensity, log2_mean_intensity, group 等
##   - results/lipid/tables/lipid_class_matrix_log2mean.tsv
##       行 = sample(global_sample_id)，列 = Class，值 = log2(mean_intensity)
##   - results/lipid/tables/lipid_class_DE.tsv
##       每行 = Class，包含 limma 差异分析结果（logFC, P.Value, adj.P.Val, ...）
##
##   - results/lipid/plots/lipid_class_bar_WT_vs_HO.png
##   - results/lipid/plots/lipid_class_boxplot_by_group.png
##
## 依赖：
##   - 00_lipid_downstream_config.yaml 中的 paths、class 配置（若有）
##
## 注意：
##   - intensity 列默认使用 intensity_pqn_bb
##   - log2 之前加一个很小的 pseudo_count，避免 0 取 log2 出现 -Inf
## ============================================================

suppressPackageStartupMessages({
  library(yaml)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(ggplot2)
  library(limma)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

log_msg <- function(...){
  cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), ..., "\n")
}

## ------------------------------------------------------------
## 0. 读配置
## ------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1) args[1] else "scripts/lipid/00_lipid_downstream_config.yaml"

if (!file.exists(config_path)) {
  stop("Config file not found: ", config_path,
       "\nUsage: Rscript 10_lipid_class_summary.R path/to/00_lipid_downstream_config.yaml")
}

log_msg("[INFO] 使用配置文件: ", config_path)
cfg   <- yaml::read_yaml(config_path)
paths <- cfg$paths

tables_dir <- paths$tables_dir %||% "results/lipid/tables"
plots_dir  <- paths$plots_dir  %||% "results/lipid/plots"
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_dir,  recursive = TRUE, showWarnings = FALSE)

## intensity 列名 & log2 pseudo-count
intensity_col   <- cfg$class$intensity_col   %||% "intensity_pqn_bb"
pseudo_count    <- cfg$class$pseudo_count    %||% 1e-6

log_msg("[INFO] Class summary 使用 intensity 列: ", intensity_col)
log_msg("[INFO] log2(mean_intensity + pseudo_count), pseudo_count = ", pseudo_count)

## ------------------------------------------------------------
## 1. 读 long 表 & sample_info
## ------------------------------------------------------------

long_path <- paths$long_lipid_pqn_bb %||%
  file.path(tables_dir, "lipid_long_dedup_pqn_bb.tsv")

sample_info_path <- paths$sample_info_output %||%
  file.path(tables_dir, "sample_info_lipid.tsv")

if (!file.exists(long_path)) {
  stop("long 表不存在: ", long_path)
}
if (!file.exists(sample_info_path)) {
  stop("sample_info 文件不存在: ", sample_info_path)
}

log_msg("[INFO] 读取 long 表: ", long_path)
long <- readr::read_tsv(long_path, show_col_types = FALSE)

log_msg("[INFO] 读取 sample_info: ", sample_info_path)
sample_info <- readr::read_tsv(sample_info_path, show_col_types = FALSE)

## 基本检查
required_long_cols <- c("batch", "sample_id", "Class", intensity_col)
miss_long <- setdiff(required_long_cols, colnames(long))
if (length(miss_long) > 0) {
  stop("[ERROR] long 表缺少必要列: ", paste(miss_long, collapse = ", "))
}

required_si_cols <- c("global_sample_id", "batch", "sample_id", "group")
miss_si <- setdiff(required_si_cols, colnames(sample_info))
if (length(miss_si) > 0) {
  stop("[ERROR] sample_info 缺少必要列: ", paste(miss_si, collapse = ", "))
}

## 在 long 里构建 global_sample_id，并与 sample_info 合并 group 信息
long2 <- long %>%
  mutate(
    global_sample_id = paste0(batch, "_", sample_id)
  ) %>%
  left_join(
    sample_info %>%
      select(global_sample_id, group),
    by = "global_sample_id"
  )

log_msg("[INFO] long 表行数: ", nrow(long2))
log_msg("[INFO] 其中 group 非 NA 的行数: ",
        sum(!is.na(long2$group)), "（非 HO/WT 的如 QC 可能为 NA）")

## ------------------------------------------------------------
## 2. Class × sample 聚合
## ------------------------------------------------------------

log_msg("[STEP] 聚合到 Class × sample 层面 ...")

class_summary <- long2 %>%
  group_by(Class, global_sample_id, group) %>%  # group 可能为 NA（QC 等）
  summarise(
    sum_intensity   = sum(.data[[intensity_col]], na.rm = TRUE),
    mean_intensity  = mean(.data[[intensity_col]], na.rm = TRUE),
    .groups         = "drop"
  ) %>%
  mutate(
    log2_mean_intensity = log2(mean_intensity + pseudo_count)
  )

log_msg("[INFO] Class × sample 行数: ", nrow(class_summary))
log_msg("[INFO] Class 种类数: ", n_distinct(class_summary$Class))
log_msg("[INFO] 样本数(按 global_sample_id): ",
        n_distinct(class_summary$global_sample_id))

## 写出长表形式的 class summary
class_summary_path <- file.path(tables_dir, "lipid_class_summary.tsv")
readr::write_tsv(class_summary, class_summary_path)
log_msg("[OK] 写出 Class summary 表: ", class_summary_path)

## ------------------------------------------------------------
## 3. 构建 Class × sample 矩阵（log2 mean）
## ------------------------------------------------------------

log_msg("[STEP] 构建 Class × sample 矩阵 (log2_mean_intensity) ...")

class_matrix_df <- class_summary %>%
  select(global_sample_id, Class, log2_mean_intensity) %>%
  tidyr::pivot_wider(
    id_cols     = global_sample_id,
    names_from  = Class,
    values_from = log2_mean_intensity
  )

## 按 sample_info 排序（方便和其它分析对齐）
class_matrix_df <- class_matrix_df %>%
  arrange(match(global_sample_id, sample_info$global_sample_id))

class_matrix_path <- file.path(tables_dir, "lipid_class_matrix_log2mean.tsv")
readr::write_tsv(class_matrix_df, class_matrix_path)
log_msg("[OK] 写出 Class matrix 表: ", class_matrix_path)

## ------------------------------------------------------------
## 4. Class 层面的 HO vs WT 差异分析（limma）
## ------------------------------------------------------------

log_msg("[STEP] Class 层面 HO vs WT 差异分析 (limma) ...")

## 只保留有 group 信息 & group 为 HO / WT 的样本
si_for_de <- sample_info %>%
  filter(!is.na(group)) %>%
  filter(group %in% c("WT", "HO"))

if (nrow(si_for_de) < 2) {
  stop("[ERROR] 用于 DE 的 HO/WT 样本数不足。")
}

## 与 class_matrix_df 对齐
class_matrix_de_df <- class_matrix_df %>%
  filter(global_sample_id %in% si_for_de$global_sample_id)

## 确保行顺序与 si_for_de 一致
class_matrix_de_df <- class_matrix_de_df %>%
  arrange(match(global_sample_id, si_for_de$global_sample_id))

si_for_de <- si_for_de %>%
  arrange(match(global_sample_id, class_matrix_de_df$global_sample_id))

## 构建表达矩阵
expr_mat <- class_matrix_de_df %>%
  select(-global_sample_id) %>%
  as.matrix()

rownames(expr_mat) <- class_matrix_de_df$global_sample_id

## 去掉全 NA 的 class 列（如果某个 class 在所有样本都 NA）
keep_class_cols <- which(colSums(!is.na(expr_mat)) > 0)
expr_mat <- expr_mat[, keep_class_cols, drop = FALSE]

if (ncol(expr_mat) == 0) {
  stop("[ERROR] 所有 Class 在 HO/WT 样本中都是 NA，无法做 DE。")
}

log_msg("[INFO] 用于 DE 的样本数: ", nrow(expr_mat),
        " ; Class 数(列数): ", ncol(expr_mat))

## 构建设计矩阵（WT 为 baseline）
si_for_de$group <- factor(si_for_de$group, levels = c("WT","HO"))
design <- model.matrix(~ group, data = si_for_de)
colnames(design)  # Intercept, groupHO

## 用 limma 拟合
fit <- lmFit(t(expr_mat), design)  # 注意：行 = feature(Class)，列 = sample
fit <- eBayes(fit)

tt <- topTable(fit,
               coef      = "groupHO",   # HO 相对 WT
               number    = Inf,
               sort.by   = "none")

## 把行名视为 Class
tt$Class <- rownames(tt)

## 整理列顺序
class_de <- tt %>%
  select(Class, logFC, AveExpr, t, P.Value, adj.P.Val, B)

class_de_path <- file.path(tables_dir, "lipid_class_DE.tsv")
readr::write_tsv(class_de, class_de_path)
log_msg("[OK] 写出 Class 层面的 DE 结果: ", class_de_path)

## ------------------------------------------------------------
## 5. 出图：Class 柱状图（WT vs HO） & 箱线图
## ------------------------------------------------------------

## 仅 WT / HO 的样本用于作图
class_summary_de <- class_summary %>%
  filter(group %in% c("WT","HO"))

## 5.1 柱状图：按 Class 汇总 WT vs HO（log2_mean_intensity 的样本均值 ± SEM）

log_msg("[STEP] 绘制 Class 柱状图 (WT vs HO) ...")

class_bar_df <- class_summary_de %>%
  group_by(Class, group) %>%
  summarise(
    n_sample   = n(),
    mean_log2  = mean(log2_mean_intensity, na.rm = TRUE),
    sd_log2    = sd(log2_mean_intensity,   na.rm = TRUE),
    sem_log2   = sd_log2 / sqrt(n_sample),
    .groups    = "drop"
  )

## 为了图的可读性，按 WT 的 mean_log2 排序 Class
class_levels <- class_bar_df %>%
  filter(group == "WT") %>%
  arrange(desc(mean_log2)) %>%
  pull(Class)

class_bar_df$Class <- factor(class_bar_df$Class, levels = class_levels)

p_bar <- ggplot(class_bar_df,
                aes(x = Class, y = mean_log2, fill = group)) +
  geom_col(position = position_dodge(width = 0.8)) +
  geom_errorbar(
    aes(ymin = mean_log2 - sem_log2,
        ymax = mean_log2 + sem_log2),
    position = position_dodge(width = 0.8),
    width = 0.3
  ) +
  theme_bw() +
  labs(
    title = "Lipid class summary (WT vs HO)",
    subtitle = "log2(mean_intensity) ± SEM",
    x = "Lipid Class",
    y = "log2(mean_intensity)",
    fill = "Group"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title  = element_text(face = "bold"),
    axis.title  = element_text(face = "bold")
  )

class_bar_png <- file.path(plots_dir, "lipid_class_bar_WT_vs_HO.png")
ggsave(class_bar_png, p_bar, width = 8, height = 5, dpi = 300)
log_msg("[OK] 写出柱状图: ", class_bar_png)

## 5.2 箱线图：每个 Class 内 WT / HO 样本的分布

log_msg("[STEP] 绘制 Class 箱线图 (by group) ...")

## 使用同样的 Class 排序
class_summary_de$Class <- factor(class_summary_de$Class, levels = class_levels)

p_box <- ggplot(class_summary_de,
                aes(x = Class, y = log2_mean_intensity, fill = group)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.size = 0.8) +
  theme_bw() +
  labs(
    title = "Lipid class distribution by group",
    subtitle = "log2(mean_intensity) per sample",
    x = "Lipid Class",
    y = "log2(mean_intensity)",
    fill = "Group"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title  = element_text(face = "bold"),
    axis.title  = element_text(face = "bold")
  )

class_box_png <- file.path(plots_dir, "lipid_class_boxplot_by_group.png")
ggsave(class_box_png, p_box, width = 8, height = 5, dpi = 300)
log_msg("[OK] 写出箱线图: ", class_box_png)

## ------------------------------------------------------------
## 6. 结束
## ------------------------------------------------------------

log_msg("============================================================")
log_msg("[DONE] 10_lipid_class_summary.R 完成。")
log_msg("  - Class summary 表: ", class_summary_path)
log_msg("  - Class matrix 表:  ", class_matrix_path)
log_msg("  - Class DE 结果:    ", class_de_path)
log_msg("  - 柱状图:           ", class_bar_png)
log_msg("  - 箱线图:           ", class_box_png)
log_msg("============================================================")