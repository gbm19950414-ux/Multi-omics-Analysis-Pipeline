#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(stringr); library(ggplot2); 
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: 09_deg_filter.R <deseq2_results.tsv> <outdir> <counts_tsv> [padj=0.05] [lfc=1] [topN=100]\n",
       " - deseq2_results.tsv: 由 08_deseq2.R 产出的差异分析结果\n",
       " - outdir: 输出目录\n",
       " - counts_tsv: 合并后的原始计数矩阵（用于画热图）。若不画热图可填 'none'\n",
       " - padj: FDR 阈值（默认 0.05）\n",
       " - lfc: 绝对 log2FC 阈值（默认 1）\n",
       " - topN: 热图展示前 topN 个基因（按 padj 排序，默认 100）"
  )
}

res_tsv   <- args[1]
outdir    <- args[2]
counts_tsv <- args[3]
padj_thr  <- ifelse(length(args)>=4, as.numeric(args[4]), 0.05)
lfc_thr   <- ifelse(length(args)>=5, as.numeric(args[5]), 1)
topN      <- ifelse(length(args)>=6, as.integer(args[6]), 100)

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

message("[INFO] res_tsv=", res_tsv)
message("[INFO] outdir=", outdir)
message("[INFO] counts_tsv=", counts_tsv)
message("[INFO] padj_thr=", padj_thr, "  lfc_thr=", lfc_thr, "  topN=", topN)

# ---------------- 读取 DESeq2 结果 ----------------
res <- read_tsv(res_tsv, show_col_types = FALSE)
if (!all(c("GeneID","log2FoldChange","pvalue","padj") %in% names(res))) {
  stop("输入结果表缺少必要列：GeneID, log2FoldChange, pvalue, padj")
}

# 清理并排序
res_clean <- res %>%
  filter(!is.na(padj)) %>%
  arrange(padj, pvalue)

# ---------------- 筛选 DEGs ----------------
deg <- res_clean %>% filter(padj < padj_thr, abs(log2FoldChange) > lfc_thr)
deg_up   <- deg %>% filter(log2FoldChange >  lfc_thr)
deg_down <- deg %>% filter(log2FoldChange < -lfc_thr)

# 输出表
out_all   <- file.path(outdir, "DEG_all.tsv")
out_up    <- file.path(outdir, "DEG_up.tsv")
out_down  <- file.path(outdir, "DEG_down.tsv")
write_tsv(deg,      out_all)
write_tsv(deg_up,   out_up)
write_tsv(deg_down, out_down)

# 简单统计
stats <- tibble::tibble(
  n_all = nrow(res_clean),
  n_DE  = nrow(deg),
  n_up  = nrow(deg_up),
  n_down= nrow(deg_down),
  padj_thr = padj_thr,
  lfc_thr  = lfc_thr
)
write_tsv(stats, file.path(outdir, "DEG_stats.tsv"))
message("[STATS] n_all=", nrow(res_clean), "  n_DE=", nrow(deg),
        " (up=", nrow(deg_up), ", down=", nrow(deg_down), ")")

# ---------------- 火山图 ----------------
vol_df <- res_clean %>%
  mutate(sig = case_when(
    padj < padj_thr & log2FoldChange >  lfc_thr ~ "Up",
    padj < padj_thr & log2FoldChange < -lfc_thr ~ "Down",
    TRUE ~ "NS"
  ))

vol <- ggplot(vol_df, aes(x=log2FoldChange, y=-log10(pvalue), color=sig)) +
  geom_point(alpha=0.6, size=1) +
  scale_color_manual(values=c("Down"="#377eb8","NS"="grey70","Up"="#e41a1c")) +
  geom_vline(xintercept=c(-lfc_thr, lfc_thr), linetype="dashed", color="grey50") +
  geom_hline(yintercept = -log10(max(min(res_clean$pvalue, 1e-300), 1e-300)), alpha=0) + # 保持尺度稳定
  theme_bw() + labs(title="Volcano plot",
    x="log2FoldChange", y="-log10(p-value)")
ggsave(filename=file.path(outdir,"volcano.png"), plot=vol, width=6, height=5, dpi=180)

# ---------------- Top 基因热图（可选） ----------------
if (!identical(tolower(counts_tsv), "none")) {
  has_pheatmap <- requireNamespace("pheatmap", quietly = TRUE)
  has_deseq2   <- requireNamespace("DESeq2", quietly = TRUE)
  if (!has_pheatmap || !has_deseq2) {
    msg <- paste0("[WARN] 缺少依赖，跳过热图：",
                  ifelse(!has_pheatmap, "pheatmap ", ""),
                  ifelse(!has_deseq2, "DESeq2 ", ""))
    message(msg)
  } else {
    suppressPackageStartupMessages({ library(DESeq2) })

    # 读取计数并构建样本信息
    cts <- read_tsv(counts_tsv, show_col_types = FALSE)
    stopifnot("GeneID" %in% names(cts))
    mat <- as.matrix(cts[,-1]); mode(mat) <- "integer"; rownames(mat) <- cts$GeneID

    samples <- colnames(mat)
    batch <- str_extract(samples, "^batch\\d+")
    cond  <- ifelse(str_detect(samples, "_KO_"), "KO",
             ifelse(str_detect(samples, "_WT_"), "WT", NA))
    if (any(is.na(cond))) {
      message("[WARN] 无法从列名推断 KO/WT，热图仅按样本名聚类。")
    }
    coldata <- data.frame(
      sample = samples,
      batch  = factor(batch),
      condition = factor(cond, levels=c("WT","KO"))
    )
    rownames(coldata) <- coldata$sample

    # 用同样设计矩阵重建 VST
    dds <- DESeqDataSetFromMatrix(countData = mat, colData = coldata, design = ~ batch + condition)
    vsd <- vst(dds, blind = FALSE)
    vsd_mat <- SummarizedExperiment::assay(vsd)

    # 取 topN 基因（按 padj）
    top_genes <- head(res_clean$GeneID, topN)
    top_genes <- intersect(top_genes, rownames(vsd_mat))
    if (length(top_genes) >= 2) {
      mat_top <- vsd_mat[top_genes, , drop=FALSE]
      # 侧注释（若可用）
      ann <- NULL
      if (!any(is.na(cond))) {
        ann <- data.frame(batch = coldata$batch, condition = coldata$condition)
        rownames(ann) <- rownames(coldata)
      }
      pheatmap::pheatmap(
        mat_top, scale="row", show_rownames=FALSE, fontsize_col=8,
        annotation_col = ann,
        filename = file.path(outdir,"heatmap_top_vst.png"), width=8, height=6
      )
      message("[OK] heatmap_top_vst.png 已生成")
    } else {
      message("[WARN] 可用于热图的 top 基因数不足，跳过热图。")
    }
  }
}

message("[SUMMARY] 输出目录：", outdir)