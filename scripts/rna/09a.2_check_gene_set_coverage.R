#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(yaml)
  library(dplyr)
  library(purrr)
  library(readr)
  library(tidyr)
  library(tibble)
})

message("============================================================")
message("[INFO] 09a.2_check_gene_set_coverage.R")
message("  yaml : scripts/rna/rna_mechanistic_gene_sets.yaml")
message("  map  : ref/gene_symbol_to_ensembl.tsv")
message("============================================================")

yaml_path <- "scripts/rna/rna_mechanistic_gene_sets.yaml"
map_path  <- "ref/gene_symbol_to_ensembl.tsv"

#------------------------------------------------------------
# STEP1 读取 YAML，抽出所有机制基因（全局 unique）
#------------------------------------------------------------
cfg <- yaml::read_yaml(yaml_path)

axes <- cfg$axes

# 抽取成 tibble: axis, module, gene_symbol
gene_tbl <- imap_dfr(
  axes,
  function(axis_obj, axis_name) {
    modules <- axis_obj$modules %||% list()

    imap_dfr(
      modules,
      function(mod_obj, mod_name) {
        genes <- mod_obj$genes %||% character(0)

        tibble(
          axis        = axis_name,
          module      = mod_name,
          gene_symbol = as.character(genes)
        )
      }
    )
  }
)

# 去掉 NA / 空串
gene_tbl <- gene_tbl %>%
  filter(!is.na(gene_symbol), gene_symbol != "")

# 全局 unique 基因列表
all_genes_unique <- gene_tbl %>%
  distinct(gene_symbol) %>%
  arrange(gene_symbol)

message("[INFO] YAML 中机制基因（含重复）总数: ", nrow(gene_tbl))
message("[INFO] YAML 中机制基因（全局 unique）数: ", nrow(all_genes_unique))

# 统计每个 axis 内 unique gene 数
message("[INFO] 每个 axis 内 unique 基因数：")
gene_tbl %>%
  distinct(axis, gene_symbol) %>%
  count(axis, name = "n_genes_axis_unique") %>%
  print()

#------------------------------------------------------------
# STEP2 读取映射表并比对覆盖情况
#------------------------------------------------------------
mapping <- read_tsv(map_path, show_col_types = FALSE)

# 这里假定映射表有列: gene_symbol, ensembl_gene_id
if (!all(c("gene_symbol", "ensembl_gene_id") %in% colnames(mapping))) {
  stop(
    "[FATAL] 映射表中必须包含列: gene_symbol, ensembl_gene_id\n",
    "  实际列名: ", paste(colnames(mapping), collapse = ", ")
  )
}

# 合并：看看哪些 symbol 有匹配，哪些没有
coverage_tbl <- all_genes_unique %>%
  left_join(mapping, by = "gene_symbol")

n_total   <- nrow(coverage_tbl)
n_mapped  <- sum(!is.na(coverage_tbl$ensembl_gene_id))
n_unmapped <- sum(is.na(coverage_tbl$ensembl_gene_id))

message("------------------------------------------------------------")
message("[SUMMARY] YAML 机制基因与映射表覆盖情况")
message("  全局 unique 基因数   : ", n_total)
message("  有 Ensembl 映射的基因 : ", n_mapped)
message("  无 Ensembl 映射的基因 : ", n_unmapped)
message("------------------------------------------------------------")

#------------------------------------------------------------
# STEP3 输出缺失映射的基因列表，方便人工修正
#------------------------------------------------------------
unmapped_tbl <- coverage_tbl %>%
  filter(is.na(ensembl_gene_id)) %>%
  arrange(gene_symbol)

if (nrow(unmapped_tbl) > 0) {
  out_unmapped <- "ref/rna_mechanistic_genes_unmapped.tsv"
  write_tsv(unmapped_tbl, out_unmapped)
  message("[WARN] 有 ", nrow(unmapped_tbl), " 个机制基因在映射表中找不到 Ensembl ID。")
  message("       已写出到: ", out_unmapped)
  message("       请检查：大小写、人鼠物种、旧名/新名、拼写错误等。")
} else {
  message("[OK] 所有机制基因在映射表中都有对应的 Ensembl ID。")
}

#------------------------------------------------------------
# STEP4（可选）输出一个“带 axis/module 信息的完整表”，方便审阅
#------------------------------------------------------------
full_out <- "ref/rna_mechanistic_genes_with_mapping.tsv"

full_tbl <- gene_tbl %>%
  distinct(axis, module, gene_symbol) %>%
  left_join(mapping, by = "gene_symbol") %>%
  arrange(axis, module, gene_symbol)

write_tsv(full_tbl, full_out)
message("[OK] 写出 axis/module × gene × mapping 的完整表: ", full_out)

message("============================================================")
message("[DONE] 09a.2_check_gene_set_coverage.R 结束")
message("============================================================")