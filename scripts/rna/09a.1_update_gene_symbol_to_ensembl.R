#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(AnnotationDbi)
  library(org.Mm.eg.db)
  library(yaml)
  library(purrr)
  library(tools)
})

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1) args[1] else "scripts/multi/rna_upstream_pathway_sets.yaml"

message("[INFO] 读取配置文件: ", config_path)
config <- yaml::read_yaml(config_path)

if (is.null(config$pathways) || length(config$pathways) == 0) {
  stop("[ERROR] 配置文件中未找到有效的 pathways 字段")
}

# 提取所有基因 symbol（从每个 pathway 的 genes$name 中）
symbols <- config$pathways %>%
  purrr::map("genes") %>%
  purrr::compact() %>%
  purrr::map(function(glist) {
    # glist 是该 pathway 下的一串基因条目
    purrr::map_chr(glist, function(g) g$name %||% NA_character_)
  }) %>%
  unlist(use.names = FALSE) %>%
  unique()

# 去除 NA 和空字符串
symbols <- symbols[!is.na(symbols) & symbols != ""]

message("[INFO] 找到的唯一基因 symbol 数: ", length(symbols))

mapping <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys = symbols,
  keytype = "SYMBOL",
  columns = c("SYMBOL", "ENSEMBL")
)

mapping <- as.data.frame(mapping, stringsAsFactors = FALSE)

message("[INFO] select() 返回行数: ", nrow(mapping))

colnames(mapping)[colnames(mapping) == "SYMBOL"] <- "gene_symbol"
colnames(mapping)[colnames(mapping) == "ENSEMBL"] <- "ensembl_gene_id"

out_dir <- dirname(config_path)
yaml_base <- tools::file_path_sans_ext(basename(config_path))
out_path <- file.path(out_dir, paste0(yaml_base, "_gene_registry.tsv"))

write_tsv(mapping, out_path)

message("[INFO] 映射中包含基因 symbol 数: ", length(unique(mapping$gene_symbol)))
message("[INFO] 已写出基因注册表文件: ", out_path)