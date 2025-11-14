#!/usr/bin/env Rscript

## ============================================================
## 05_between_batch_panel_scaling.R
##
## 目的：
##   在「批次内 PQN 矫正」之后，利用“公共 feature panel”
##   做「无 anchor、对称」的批次间 scaling。
##
## 核心思想：
##   - Panel A（全脂 panel）：batch1, batch2, batch3 之间的公共 lipid
##     · 使用 keep_for_correction == TRUE 且有 QC 的 feature
##     · 在 batch1/2/3 都有 QC 信号的 lipid 组成 Panel A
##     · 在 Panel A 上，对每个 lipid 构建“虚拟参考强度”：
##         ref_lipid_center[l] = median_b( qc_center_pqn[b, l] )
##       对每个 batch 计算 ratio[b,l] = ref / center，
##       batch_factor_A[b] = median_l( ratio[b,l] )
##
##   - Panel B（CL 子 panel）：batch2, batch3, batch4 之间的 “共同 CL”
##     · 只考虑 Class == "CL"
##     · 在 batch2/3/4 都有 QC 信号的 CL 组成 Panel B
##     · 同样构建 ref_lipid_center_CL 和 batch_factor_B[b]
##
##   - Bridge（桥接两个 panel 的标尺）：
##     · batch2 和 batch3 同时出现在 Panel A / Panel B
##     · 在这两个 batch 上，有 factor_A[b] 和 factor_B[b]
##     · 定义一个缩放常数：
##         k = median_b( factor_A[b] / factor_B[b] , b ∈ {2,3} )
##     · 最终：
##         · batch1–3 使用 factor_A[b]
##         · batch4 使用 factor_final[4] = k * factor_B[4]
##       => 全部 batch 最终都落在 “Panel A 标尺” 上（以 Panel A 为 global scale）
##
## 输入：
##   - results/lipid/tables/lipid_long_dedup_pqn.tsv
##     · 批次内 PQN 之后的长表
##     · 需要包含列：
##         batch, sample_id, lipidName, intensity_pqn,
##         is_qc, keep_for_correction, Class
##
## 输出：
##   - results/lipid/tables/lipid_long_dedup_pqn_bb.tsv
##       · 在 PQN 基础上加入 batch 间 scaling 的强度：
##         intensity_pqn_bb
##   - results/lipid/qc/lipid_between_batch_panel_scaling_factors.tsv
##       · 每个 batch 的最终 batch_factor_between_panel
##   - results/lipid/qc/qa_between_batch_panelA_batch_factors.tsv
##   - results/lipid/qc/qa_between_batch_panelB_batch_factors.tsv
##   - results/lipid/qc/qa_between_batch_panel_bridge.tsv
##       · 方便检查 Panel A / B、以及桥接的情况
## ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
  library(purrr)
})

tables_dir <- "results/lipid/tables"
qc_dir     <- "results/lipid/qc"

long_pqn_path <- file.path(tables_dir, "lipid_long_dedup_pqn.tsv")

cat("============================================================\n")
cat("[INFO] Between-batch PANEL scaling (after within-batch PQN)\n")
cat("  input : ", long_pqn_path, "\n", sep = "")
cat("  output: ", file.path(tables_dir, "lipid_long_dedup_pqn_bb.tsv"), "\n", sep = "")
cat("  Panel A: batch1–3 共同 lipid\n")
cat("  Panel B: batch2–4 共同 CL\n")
cat("============================================================\n\n")

## -------- 1. 读取 PQN 后长表，基础检查 -------------------------------

cat("[STEP1] 读取 lipid_long_dedup_pqn.tsv ...\n")

if (!file.exists(long_pqn_path)) {
  stop("[ERROR] 找不到文件: ", long_pqn_path)
}

long_pqn <- readr::read_tsv(long_pqn_path, guess_max = 1e6)

cat("  [long_pqn] nrow = ", nrow(long_pqn),
    " ; ncol = ", ncol(long_pqn), "\n", sep = "")

required_cols <- c(
  "batch", "sample_id", "lipidName",
  "intensity_pqn", "is_qc", "keep_for_correction"
)

missing_cols <- setdiff(required_cols, colnames(long_pqn))
if (length(missing_cols) > 0) {
  stop(
    "[ERROR] '", long_pqn_path, "' 缺少必要列: ",
    paste(missing_cols, collapse = ", ")
  )
}

## Class 不是强制，但如果要构造 CL panel，就尽量需要
if (!"Class" %in% colnames(long_pqn)) {
  warning("[WARN] 长表中没有 'Class' 列，Panel B (CL-only) 无法构造，将只使用 Panel A (batch1–3)。")
}

## 确保逻辑列为 logical
to_logical <- function(x) {
  if (is.logical(x)) return(x)
  if (is.numeric(x)) return(x != 0)
  if (is.character(x)) {
    return(toupper(x) %in% c("TRUE", "T", "1", "YES", "Y"))
  }
  as.logical(x)
}

long_pqn <- long_pqn %>%
  mutate(
    is_qc = to_logical(is_qc),
    keep_for_correction = to_logical(keep_for_correction)
  )

## -------- 2. 先构建 QC 层面的 (batch, lipid) 中心值 ------------------

cat("\n[STEP2] 计算 QC 层面的 pqn 中心值（per batch, per lipid）...\n")

qc_long <- long_pqn %>%
  filter(
    keep_for_correction,
    is_qc,
    !is.na(intensity_pqn),
    intensity_pqn > 0
  ) %>%
  group_by(batch, lipidName, LipidIon, CalcMz, Class) %>%
  summarise(
    qc_center_pqn = median(intensity_pqn, na.rm = TRUE),
    n_qc = sum(!is.na(intensity_pqn)),
    .groups = "drop"
  ) %>%
  filter(n_qc > 0, qc_center_pqn > 0)

cat("  [qc_long] 有 QC 记录的 (batch, lipid) 数量: ",
    nrow(qc_long), "\n", sep = "")

batches <- sort(unique(long_pqn$batch))

## ============================================================
## Panel A: batch1–3 共同 lipid
## ============================================================

cat("\n[STEP3] 构建 Panel A（batch1–3 公共 lipid）...\n")

panelA_batches <- c("batch1", "batch2", "batch3")

qc_A <- qc_long %>%
  filter(batch %in% panelA_batches)

if (nrow(qc_A) == 0) {
  stop("[ERROR] 在 batch1–3 中没有任何 QC 记录，无法构造 Panel A。")
}

panelA_lipids <- qc_A %>%
  group_by(lipidName, LipidIon, CalcMz) %>%
  summarise(
    n_batch = n_distinct(batch),
    .groups = "drop"
  ) %>%
  filter(n_batch == length(panelA_batches))

cat("  [Panel A] 在 batch1–3 都有 QC 的 lipid 数: ",
    nrow(panelA_lipids), "\n", sep = "")

if (nrow(panelA_lipids) < 10) {
  warning("[WARN] Panel A 中可用 lipid < 10，批间 scaling 可能不稳定，请人工检查。")
}

qc_panelA <- qc_A %>%
  inner_join(panelA_lipids,
             by = c("lipidName", "LipidIon", "CalcMz")) %>%
  select(-n_batch)

## 对 Panel A 上的每个 lipid 构建“虚拟参考强度”
ref_panelA <- qc_panelA %>%
  group_by(lipidName, LipidIon, CalcMz) %>%
  summarise(
    ref_center_pqn = median(qc_center_pqn, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(!is.na(ref_center_pqn), ref_center_pqn > 0)

## 计算 batch-level factor（Panel A）
ratios_A <- qc_panelA %>%
  inner_join(ref_panelA,
             by = c("lipidName", "LipidIon", "CalcMz")) %>%
  mutate(
    ratio = ref_center_pqn / qc_center_pqn
  ) %>%
  filter(is.finite(ratio), ratio > 0)

batch_factor_A <- ratios_A %>%
  group_by(batch) %>%
  summarise(
    n_features_used        = n(),
    factor_between_panelA  = median(ratio, na.rm = TRUE),
    factor_min             = min(ratio, na.rm = TRUE),
    factor_q1              = quantile(ratio, 0.25, na.rm = TRUE),
    factor_median          = median(ratio, na.rm = TRUE),
    factor_q3              = quantile(ratio, 0.75, na.rm = TRUE),
    factor_max             = max(ratio, na.rm = TRUE),
    .groups = "drop"
  )

cat("  [Panel A] batch-level factor 概要：\n")
print(batch_factor_A)

## ============================================================
## Panel B: batch2–4 共享的 CL panel（可选）
## ============================================================

cat("\n[STEP4] 构建 Panel B（batch2–4 公共 CL；如无 Class 则跳过）...\n")

panelB_batches <- c("batch2", "batch3", "batch4")
batch_factor_B <- NULL
panel_bridge   <- NULL

if (!"Class" %in% colnames(qc_long)) {
  cat("  [INFO] 未找到 Class 列，跳过 Panel B（CL-only） 构建。\n")
} else {
  qc_B <- qc_long %>%
    filter(batch %in% panelB_batches,
           Class == "CL")

  if (nrow(qc_B) == 0) {
    cat("  [INFO] 在 batch2–4 中没有 CL 的 QC 记录，跳过 Panel B。\n")
  } else {
    panelB_lipids <- qc_B %>%
      group_by(lipidName, LipidIon, CalcMz) %>%
      summarise(
        n_batch = n_distinct(batch),
        .groups = "drop"
      ) %>%
      filter(n_batch == length(panelB_batches))

    cat("  [Panel B] 在 batch2–4 都有 QC 的 CL lipid 数: ",
        nrow(panelB_lipids), "\n", sep = "")

    if (nrow(panelB_lipids) < 5) {
      warning("[WARN] Panel B 中可用 CL lipid < 5，批间 scaling 对 batch4 可能不稳定。")
    }

    qc_panelB <- qc_B %>%
      inner_join(panelB_lipids,
                 by = c("lipidName", "LipidIon", "CalcMz")) %>%
      select(-n_batch)

    ref_panelB <- qc_panelB %>%
      group_by(lipidName, LipidIon, CalcMz) %>%
      summarise(
        ref_center_pqn_CL = median(qc_center_pqn, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      filter(!is.na(ref_center_pqn_CL), ref_center_pqn_CL > 0)

    ratios_B <- qc_panelB %>%
      inner_join(ref_panelB,
                 by = c("lipidName", "LipidIon", "CalcMz")) %>%
      mutate(
        ratio_CL = ref_center_pqn_CL / qc_center_pqn
      ) %>%
      filter(is.finite(ratio_CL), ratio_CL > 0)

    batch_factor_B <- ratios_B %>%
      group_by(batch) %>%
      summarise(
        n_features_used        = n(),
        factor_between_panelB  = median(ratio_CL, na.rm = TRUE),
        factor_min             = min(ratio_CL, na.rm = TRUE),
        factor_q1              = quantile(ratio_CL, 0.25, na.rm = TRUE),
        factor_median          = median(ratio_CL, na.rm = TRUE),
        factor_q3              = quantile(ratio_CL, 0.75, na.rm = TRUE),
        factor_max             = max(ratio_CL, na.rm = TRUE),
        .groups = "drop"
      )

    cat("  [Panel B] batch-level factor 概要：\n")
    print(batch_factor_B)

    ## -------- Bridge: 用 batch2 & batch3 对齐 Panel B 到 Panel A 标尺 -----

    cat("\n[STEP5] Panel A / Panel B 桥接（以 Panel A 为 global 标尺）...\n")

    bridge_batches <- intersect(
      panelA_batches,
      panelB_batches
    )  ## 即 batch2, batch3

    bridge_A <- batch_factor_A %>%
      filter(batch %in% bridge_batches) %>%
      select(batch, factor_between_panelA)

    bridge_B <- batch_factor_B %>%
      filter(batch %in% bridge_batches) %>%
      select(batch, factor_between_panelB)

    panel_bridge <- inner_join(bridge_A, bridge_B, by = "batch") %>%
      mutate(
        ratio_A_over_B = factor_between_panelA / factor_between_panelB
      )

    if (nrow(panel_bridge) < 1) {
      warning("[WARN] Panel bridge 上 batch2/3 不完整，无法桥接 Panel B -> Panel A；batch4 将无法通过 Panel B 单独矫正。")
    } else {
      k_bridge <- median(panel_bridge$ratio_A_over_B,
                         na.rm = TRUE)

      cat("  [Bridge] batch2/3 上 A/B 比例：\n")
      print(panel_bridge)
      cat("  [Bridge] k_bridge = median( factor_A / factor_B ) = ",
          signif(k_bridge, 4), "\n", sep = "")

      panel_bridge <- panel_bridge %>%
        mutate(k_bridge = k_bridge)
    }
  }
}

## ============================================================
## 6. 合成最终的 batch-level factor
##    - batch1–3 使用 Panel A 的 factor
##    - batch4 使用 Panel B factor * k_bridge（如可用）
## ============================================================

cat("\n[STEP6] 合成最终 batch-level factor（panel-based）...\n")

## 先准备所有 batch 的骨架
batch_df <- tibble(batch = sort(unique(long_pqn$batch)))

## 从 Panel A 拿 batch1–3 的 factor
final_factors <- batch_df %>%
  left_join(
    batch_factor_A %>% select(batch, factor_between_panelA, n_features_used_panelA = n_features_used),
    by = "batch"
  )

## Panel B 的 factor（可能为空）
if (!is.null(batch_factor_B)) {
  final_factors <- final_factors %>%
    left_join(
      batch_factor_B %>%
        select(batch, factor_between_panelB, n_features_used_panelB = n_features_used),
      by = "batch"
    )
}

## 默认最终 factor = Panel A 的（对于 batch1–3）
final_factors <- final_factors %>%
  mutate(
    factor_final = factor_between_panelA,
    source_panel = ifelse(!is.na(factor_between_panelA), "PanelA", NA_character_)
  )

## 对 batch4，如 Panel B + 桥接 k 可用，则 factor_final = k * factor_B
if (!is.null(batch_factor_B) && !is.null(panel_bridge) &&
    "batch4" %in% final_factors$batch) {

  if ("k_bridge" %in% colnames(panel_bridge)) {
    k_bridge_val <- unique(panel_bridge$k_bridge)

    final_factors <- final_factors %>%
      mutate(
        factor_final = dplyr::if_else(
          batch == "batch4" & !is.na(factor_between_panelB),
          k_bridge_val * factor_between_panelB,
          factor_final
        ),
        source_panel = dplyr::if_else(
          batch == "batch4" & !is.na(factor_between_panelB),
          "PanelB_via_bridge",
          source_panel
        )
      )
  } else {
    warning("[WARN] Panel bridge 未成功构建，batch4 无法通过 Panel B 进行矫正；其 factor_final 将保持 NA。")
  }
}

## 如果还有 batch 没有 factor_final（例如 future 的其它 batch），默认设为 1
final_factors <- final_factors %>%
  mutate(
    factor_final = ifelse(is.na(factor_final), 1, factor_final),
    source_panel = ifelse(is.na(source_panel), "default_1", source_panel)
  )

cat("  [Final factors] 概要：\n")
print(final_factors)

## ============================================================
## 7. 应用 final factor 到长表，生成 intensity_pqn_bb
## ============================================================

cat("\n[STEP7] 将 batch-level factor 应用到 PQN 后长表，生成 intensity_pqn_bb ...\n")

long_pqn_bb <- long_pqn %>%
  left_join(
    final_factors %>%
      select(batch, batch_factor_between_panel = factor_final),
    by = "batch"
  ) %>%
  mutate(
    batch_factor_between_panel = ifelse(
      is.na(batch_factor_between_panel) |
        !is.finite(batch_factor_between_panel) |
        batch_factor_between_panel <= 0,
      1,
      batch_factor_between_panel
    ),
    intensity_pqn_bb = intensity_pqn * batch_factor_between_panel
  )

cat("  [long_pqn_bb] nrow = ", nrow(long_pqn_bb),
    " ; ncol = ", ncol(long_pqn_bb), "\n", sep = "")

## ============================================================
## 8. 写出 QC 表 & 最终长表
## ============================================================

cat("\n[STEP8] 写出 panel-based scaling factor & QA 表...\n")

dir.create(qc_dir, showWarnings = FALSE, recursive = TRUE)

## Panel A factor QA
readr::write_tsv(
  batch_factor_A,
  file.path(qc_dir, "qa_between_batch_panelA_batch_factors.tsv")
)

## Panel B factor QA（如有）
if (!is.null(batch_factor_B)) {
  readr::write_tsv(
    batch_factor_B,
    file.path(qc_dir, "qa_between_batch_panelB_batch_factors.tsv")
  )
}

## Bridge QA（如有）
if (!is.null(panel_bridge)) {
  readr::write_tsv(
    panel_bridge,
    file.path(qc_dir, "qa_between_batch_panel_bridge.tsv")
  )
}

## 最终 batch-level factor
readr::write_tsv(
  final_factors,
  file.path(qc_dir, "lipid_between_batch_panel_scaling_factors.tsv")
)

## 输出最终长表
dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)

out_path <- file.path(tables_dir, "lipid_long_dedup_pqn_bb.tsv")
readr::write_tsv(long_pqn_bb, out_path)

cat("  [OK] 写出 panel-based batch factor: ",
    file.path(qc_dir, "lipid_between_batch_panel_scaling_factors.tsv"), "\n", sep = "")
cat("  [OK] 写出 PQN + panel-based between-batch scaling 后的长表: ",
    out_path, "\n", sep = "")

cat("============================================================\n")
cat("[DONE] 05_between_batch_panel_scaling.R 完成。\n")
cat("  Panel A: batch1–3 全脂 panel；Panel B: batch2–4 CL panel（如有）\n")
cat("  建议：用 QC-only PCA / boxplot 对比 PQN vs PQN+panel scaling，评估批间矫正效果。\n")
cat("============================================================\n")