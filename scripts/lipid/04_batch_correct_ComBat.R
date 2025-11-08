#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(sva)       # 需安装 r-sva
  library(readr)
})

option_list <- list(
  make_option(c("--cfg"), type="character", help="config yaml"),
  make_option(c("--in"),  type="character", help="input matrix (tsv, features x samples)"),
  make_option(c("--meta"),type="character", help="sample metadata tsv (SampleID,Group,Batch)"),
  make_option(c("--out"), type="character", help="output corrected tsv")
)
opt <- parse_args(OptionParser(option_list=option_list))

# 读入
X <- as.data.frame(readr::read_tsv(opt$`in`))
rownames(X) <- X[[1]]; X[[1]] <- NULL

meta <- readr::read_tsv(opt$meta, show_col_types = FALSE)
stopifnot(all(colnames(X) %in% meta$SampleID))
meta <- meta[match(colnames(X), meta$SampleID), ]
batch <- as.factor(meta$Batch)
mod   <- model.matrix(~ Group, data=meta)

# ComBat（对行=特征做校正）
Xc <- ComBat(dat=as.matrix(X), batch=batch, mod=mod, par.prior=TRUE, prior.plots=FALSE)
Xc <- as.data.frame(Xc)
Xc <- cbind(Feature=rownames(Xc), Xc)
readr::write_tsv(Xc, opt$out)
cat(sprintf("[OK] ComBat corrected -> %s\n", opt$out))