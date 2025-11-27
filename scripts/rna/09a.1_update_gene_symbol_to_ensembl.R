#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# 使用方法（如何指定输入文件）
#
# 本脚本是一个通用的 Gene Registry 生成器：
#   - 接受任意 YAML 文件作为输入
#   - 在整个 YAML 结构中递归查找形如 `name: <SYMBOL>` 的字段
#   - 将所有找到的基因 symbol 映射到 Ensembl Gene ID
#   - 在输入 YAML 所在目录下输出一个 *_gene_registry.tsv 文件
#
# 运行示例：
#
#   Rscript scripts/rna/09a.1_update_gene_symbol_to_ensembl.R scripts/rna/rna_mechanistic_gene_sets.yaml
#
# 若未提供参数，则默认读取：
#   scripts/rna/rna_mechanistic_gene_sets.yaml
#
# 注意：
#   - 本脚本会递归扫描所有嵌套列表，只要字段名为 `name` 且值为非空字符串，
#     即视为候选基因 symbol。
# ------------------------------------------------------------------------------

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
config_path <- if (length(args) >= 1) args[1] else "scripts/rna/rna_mechanistic_gene_sets.yaml"

message("[INFO] 读取配置文件: ", config_path)
config <- yaml::read_yaml(config_path)

# 递归遍历任意 YAML 结构，收集所有 name: <SYMBOL> 字段
collect_names <- function(x) {
  symbols <- character(0)

  if (is.list(x)) {
    # 如果是具名列表，优先检查是否存在 name 字段
    if (!is.null(names(x))) {
      # 当前层级含有 name 字段且为非空字符
      if ("name" %in% names(x)) {
        val <- x[["name"]]
        if (is.character(val) && length(val) == 1 && !is.na(val) && nzchar(val)) {
          symbols <- c(symbols, val)
        }
      }
    }
    # 递归所有子元素
    for (elem in x) {
      symbols <- c(symbols, collect_names(elem))
    }
  }

  symbols
}

# 从整个 YAML 配置中抽取所有 name: <SYMBOL>
symbols <- collect_names(config) %>%
  unique()

# 去除 NA 和空字符串（保险起见）
symbols <- symbols[!is.na(symbols) & symbols != ""]

if (length(symbols) == 0) {
  stop("[ERROR] 在 YAML 中未找到任何有效的 `name: <SYMBOL>` 字段，请检查配置。")
}

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