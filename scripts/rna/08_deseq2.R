#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(DESeq2); library(readr); library(dplyr); library(stringr)
  library(ggplot2); library(pheatmap)
})

# Usage:
# 08_deseq2.R <counts_tsv> <outdir> [contrast=condition,KO,WT]
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 2) stop("Usage: 08_deseq2.R <counts_tsv> <outdir> [contrast=condition,KO,WT]")
counts_tsv <- args[1]
outdir     <- args[2]
contrast   <- ifelse(length(args)>=3, args[3], "condition,KO,WT")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ---------- 读取计数 ----------
cts <- read_tsv(counts_tsv, show_col_types = FALSE)
stopifnot("GeneID" %in% names(cts))
rownames_cts <- cts$GeneID
mat <- as.matrix(cts[,-1])
mode(mat) <- "integer"
rownames(mat) <- rownames_cts

# ---------- 构建样本信息（自动从列名提取 batch 与 condition） ----------
# 约定列名形如: batch1_KO_1 / batch1_WT_2 / ...
samples <- colnames(mat)
batch <- str_extract(samples, "^batch\\d+")
cond  <- ifelse(str_detect(samples, "_KO_"), "KO",
         ifelse(str_detect(samples, "_WT_"), "WT", NA))
if (any(is.na(cond))) stop("无法从列名推断 KO/WT，请检查列名。")

coldata <- data.frame(
  sample = samples,
  batch  = factor(batch),
  condition = factor(cond, levels = c("WT","KO"))
)
rownames(coldata) <- coldata$sample

# ---------- 过滤低表达 ----------
keep <- rowSums(mat) >= 10
mat_f <- mat[keep, , drop=FALSE]
message("[FILTER] kept genes: ", nrow(mat_f), " / ", nrow(mat))

# ---------- DESeq2：批次校正设计 ----------
dds <- DESeqDataSetFromMatrix(countData = mat_f, colData = coldata, design = ~ batch + condition)
dds <- DESeq(dds)

# ---------- 结果输出 ----------
# contrast 解析，如 "condition,KO,WT"
cc <- strsplit(contrast, ",")[[1]]
res <- results(dds, contrast = cc)
res <- lfcShrink(dds, contrast = cc, res = res, type = "ashr")  # 需要 ashr；若无可改 type="apeglm"/"normal"
res_tbl <- as.data.frame(res) %>%
  tibble::rownames_to_column("GeneID") %>%
  arrange(padj, pvalue)

out_res <- file.path(outdir, "deseq2_results.tsv")
write_tsv(res_tbl, out_res)
message("[OK] results: ", out_res)

# ---------- 归一化与可视化 ----------
vsd <- vst(dds, blind = FALSE)

# PCA
pca <- plotPCA(vsd, intgroup = c("batch","condition"), returnData = TRUE)
percentVar <- round(100 * attr(pca, "percentVar"))
p <- ggplot(pca, aes(PC1, PC2, color=condition, shape=batch, label=name)) +
  geom_point(size=3) + theme_bw() +
  xlab(paste0("PC1: ", percentVar[1], "%")) +
  ylab(paste0("PC2: ", percentVar[2], "%"))
ggsave(file.path(outdir, "PCA_vst.png"), p, width=6, height=4, dpi=200)

# 样本相关性热图
mat_vst <- assay(vsd)
cor_mat <- cor(mat_vst, method="spearman")
png(file.path(outdir, "sample_correlation.png"), width=900, height=800, res=140)
pheatmap(cor_mat, clustering_distance_rows="correlation", clustering_distance_cols="correlation")
dev.off()

# MA 图
png(file.path(outdir, "MA_plot.png"), width=800, height=600, res=140)
plotMA(res, ylim = c(-5,5)); abline(h=0,lty=2,col="grey50"); dev.off()

# 顶基因表
top_sig <- res_tbl %>% filter(!is.na(padj)) %>% arrange(padj) %>% head(1000)
write_tsv(top_sig, file.path(outdir, "top_DE_genes.tsv"))

# ---------- 交互项（batch × condition）检测与分批效应估计 ----------
# 说明：
# 1) 使用 LRT 检测是否存在显著的 batch × condition 交互效应；
# 2) 在包含交互项的 WALD 模型中，估计每个批次内的 KO vs WT 的 log2FC；
# 3) 导出显著交互的基因，并为前 N 个基因绘制按批次分组的表达图。

dir_ic <- file.path(outdir, "interaction")
dir.create(dir_ic, showWarnings = FALSE, recursive = TRUE)

## (1) LRT：是否存在交互效应
dds_ic_lrt <- DESeqDataSetFromMatrix(countData = mat_f, colData = coldata,
                                     design = ~ batch + condition + batch:condition)
dds_ic_lrt <- DESeq(dds_ic_lrt, test = "LRT", reduced = ~ batch + condition)
res_ic <- results(dds_ic_lrt)
res_ic_tbl <- as.data.frame(res_ic) %>%
  tibble::rownames_to_column("GeneID") %>%
  arrange(padj, pvalue)
ic_path <- file.path(dir_ic, "interaction_LRT.tsv")
write_tsv(res_ic_tbl, ic_path)
message("[OK] interaction LRT: ", ic_path)

## (2) 在包含交互项的 WALD 模型中提取每批次的 KO vs WT 效应
dds_ic_wald <- DESeqDataSetFromMatrix(countData = mat_f, colData = coldata,
                                      design = ~ batch + condition + batch:condition)
dds_ic_wald <- DESeq(dds_ic_wald)  # Wald
coef_names <- resultsNames(dds_ic_wald)

# 找到主效应和各批次交互项列名
# 主效应（参考批次的 KO vs WT）
main_coef <- coef_names[grepl("^condition_.*_vs_.*$", coef_names)]
if (length(main_coef) != 1) {
  stop("未能唯一定位 condition 主效应系数，请检查列名: ", paste(coef_names, collapse = ", "))
}
# 交互项：以 ".condition" 或 "condition" + "KO" 的形式命名（不同版本命名略有差异）
int_coefs <- coef_names[grepl("condition.*:", coef_names) | grepl("\\.condition", coef_names)]

# 批次水平
batch_lvls <- levels(coldata$batch)
ref_batch  <- batch_lvls[1]

# 提取每个基因的系数矩阵
bet <- as.data.frame(coef(dds_ic_wald))
bet$GeneID <- rownames(bet)

# 计算各批次内的 KO vs WT LFC：参考批次 = 主效应；其他批次 = 主效应 + 对应交互项
calc_batch_lfc <- function(b) {
  if (b == ref_batch) {
    if (!main_coef %in% colnames(bet)) stop("找不到主效应列：", main_coef)
    return(bet[[main_coef]])
  }
  # 交互项列名可能类似 "batchb.conditionKO" 或 "batch_b_vs_ref.condition_KO_vs_WT"
  pat <- paste0(b, ".*condition")
  cand <- int_coefs[grepl(pat, int_coefs)]
  if (length(cand) == 0) {
    warning("未找到批次 ", b, " 的交互项系数，默认按主效应处理。")
    return(bet[[main_coef]])
  }
  if (!cand[1] %in% colnames(bet)) stop("找不到交互项列：", cand[1])
  bet[[main_coef]] + bet[[cand[1]]]
}

lfc_by_batch <- lapply(batch_lvls, calc_batch_lfc)
lfc_df <- as.data.frame(lfc_by_batch)
colnames(lfc_df) <- paste0("LFC_", batch_lvls)
lfc_df$GeneID <- bet$GeneID

# 合并交互 LRT 的显著性，导出汇总表
ic_merged <- res_ic_tbl %>%
  select(GeneID, stat, pvalue, padj) %>%
  rename(IC_stat = stat, IC_pvalue = pvalue, IC_padj = padj) %>%
  left_join(lfc_df, by = "GeneID")

ic_sum_path <- file.path(dir_ic, "interaction_summary_with_perBatchLFC.tsv")
write_tsv(ic_merged, ic_sum_path)
message("[OK] interaction summary: ", ic_sum_path)

## (3) 可视化：为交互显著的前 N 个基因画按批次分组的表达图
topN <- 12
sig_ic <- ic_merged %>% filter(!is.na(IC_padj), IC_padj < 0.05) %>% head(topN) %>% pull(GeneID)

if (length(sig_ic) > 0) {
  for (g in sig_ic) {
    dfc <- plotCounts(dds_ic_wald, gene = g, intgroup = c("batch","condition"), returnData = TRUE)
    p_ic <- ggplot(dfc, aes(x=batch, y=count, color=condition, group=interaction(batch,condition))) +
      geom_point(position = position_jitter(width = 0.15, height = 0), size = 2) +
      geom_line(aes(group=condition), alpha = 0.5) +
      scale_y_log10() + theme_bw() + ggtitle(paste0(g, " (batch×condition)"))
    ggsave(filename = file.path(dir_ic, paste0("ICplot_", g, ".png")), plot = p_ic, width = 6, height = 4, dpi = 200)
  }
  message("[OK] interaction plots written to: ", dir_ic)
} else {
  message("[INFO] no significant interaction genes at padj<0.05")
}

# 简要摘要
sum_ok <- res_tbl %>% filter(!is.na(padj), padj < 0.05) %>% nrow()
message("[SUMMARY] DE genes (padj<0.05): ", sum_ok)