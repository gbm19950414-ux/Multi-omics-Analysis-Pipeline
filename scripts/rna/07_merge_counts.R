#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr)
})

args <- commandArgs(trailingOnly=TRUE)
if(length(args) < 2){
  stop("用法: 07_merge_counts.R <processed_dir> <out_tsv>")
}

proc_dir <- args[1]
out_tsv  <- args[2]

# 收集各批次 featureCounts_matrix.tsv
files <- list.files(proc_dir, pattern="featureCounts_matrix.tsv$", recursive=TRUE, full.names=TRUE)
stopifnot(length(files) > 0)

cat("[INFO] merging matrices:\n")
print(files)

# 读入并按 GeneID 合并
df_list <- lapply(files, function(f){
  read_tsv(f, show_col_types=FALSE)
})
merged <- Reduce(function(x,y) full_join(x,y, by="GeneID"), df_list) %>% replace_na(list())

# 去掉重复列名（若不同批次样本名重复，前缀加批次路径名）
dups <- names(merged)[duplicated(names(merged))]
if(length(dups) > 0){
  for(d in dups){
    idx <- which(names(merged) == d)
    for(i in seq_along(idx)){
      names(merged)[idx[i]] <- paste0(d,"__dup",i)
    }
  }
}

write_tsv(merged, out_tsv)
cat("[OK] merged counts -> ", out_tsv, "\n", sep="")
