#!/usr/bin/env Rscript
# figure_7_b.R
# 将 EphB1 下游信号 YAML 转换为论文可用的 Word 图表

## ---------------- 0. 基本设置 ----------------

options(stringsAsFactors = FALSE, encoding = "UTF-8")

required_pkgs <- c("yaml", "dplyr", "purrr", "tibble",
                   "stringr", "readr", "officer", "flextable")

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("[ERROR] 需要 R 包 '", pkg,
         "'，请先安装，例如：install.packages('", pkg, "')\n",
         call. = FALSE)
  }
}

library(yaml)
library(dplyr)
library(purrr)
library(tibble)
library(stringr)
library(readr)
library(officer)
library(flextable)

## ---------------- 1. 路径设置 ----------------

yaml_path <- "scripts/multi/ephb1_downstream_signaling_sets.yaml"
out_dir   <- "results/figs"

if (!file.exists(yaml_path)) {
  stop("[ERROR] 找不到 YAML 文件：", yaml_path,
       "\n请确认工作目录在项目根目录（multiomics_mech）。",
       call. = FALSE)
}

if (!dir.exists(out_dir)) {
  message("[INFO] 创建输出目录: ", out_dir)
  dir.create(out_dir, recursive = TRUE)
}

tsv_path  <- file.path(out_dir,
                       "figure_7_b_ephb1_downstream_signaling_table.tsv")
docx_path <- file.path(out_dir,
                       "figure_7_b_ephb1_downstream_signaling_table.docx")

## ---------------- 2. 读取 YAML ----------------

message("[STEP] 读取 YAML: ", yaml_path)
cfg <- yaml::read_yaml(yaml_path)

if (is.null(cfg$pathways)) {
  stop("[ERROR] YAML 中未找到 'pathways' 字段。", call. = FALSE)
}

pathways_list <- cfg$pathways

## ---------------- 3. 展开为数据框 ----------------

message("[STEP] 将 pathways 列表展开为数据框...")

# helper: 安全提取字符
get_chr <- function(x) {
  if (is.null(x)) return(NA_character_)
  as.character(x)
}

pathways_df <- purrr::imap_dfr(
  pathways_list,
  function(pw, pw_name) {

    axis <- get_chr(pw$axis)

    target_axes <- if (!is.null(pw$target_mechanistic_axes)) {
      paste0(pw$target_mechanistic_axes, collapse = "; ")
    } else {
      NA_character_
    }

    label <- get_chr(pw$label)
    desc  <- get_chr(pw$description)

    genes <- pw$genes %||% list()

    gene_symbols <- purrr::map_chr(genes, ~ get_chr(.x$name))
    gene_signs   <- purrr::map_chr(genes, ~ {
      if (is.null(.x$gene_effect_sign)) return(NA_character_)
      as.character(.x$gene_effect_sign)
    })

    # 纯基因名列表
    genes_concat <- if (length(gene_symbols) > 0) {
      paste(gene_symbols, collapse = ", ")
    } else {
      NA_character_
    }

    # 带 sign 的基因列表：如 "Bax (+1)"
    genes_with_sign_concat <- if (length(gene_symbols) > 0) {
      paste0(gene_symbols,
             " (", gene_signs, ")") |>
        paste(collapse = ", ")
    } else {
      NA_character_
    }

    tibble::tibble(
      pathway_id        = pw_name,
      axis              = axis,
      target_axes       = target_axes,
      label             = label,
      description       = str_squish(desc),
      n_genes           = length(genes),
      genes             = genes_concat,
      genes_with_sign   = genes_with_sign_concat
    )
  }
)

# 按 axis + pathway_id 排一下，保证 Word 表格顺序好看
pathways_df <- pathways_df %>%
  arrange(axis, pathway_id)

## ---------------- 4. 导出 TSV ----------------

message("[STEP] 写出 TSV: ", tsv_path)
readr::write_tsv(pathways_df, tsv_path)

## ---------------- 5. 生成 Word 图表 ----------------

message("[STEP] 生成 Word 图表: ", docx_path)

# 为论文准备一个更紧凑的列子集（如需可再调整）
table_df <- pathways_df %>%
  transmute(
    `Pathway ID`          = pathway_id,
    `模块名称（label）`     = label,
    `所属大轴 (axis)`      = axis,
    `目标 CL 机制轴`        = target_axes,
    `基因数 (n)`           = n_genes,
    `基因列表`              = genes_with_sign
  )

ft <- flextable::flextable(table_df)

ft <- ft |>
  flextable::autofit() |>
  flextable::align(align = "left", part = "all") |>
  flextable::valign(valign = "top", part = "all") |>
  flextable::fontsize(part = "all", size = 9) |>
  flextable::set_header_labels(
    `Pathway ID`          = "Pathway ID",
    `模块名称（label）`     = "模块名称（label）",
    `所属大轴 (axis)`      = "所属大轴 (axis)",
    `目标 CL 机制轴`        = "目标 CL 机制轴",
    `基因数 (n)`           = "基因数 (n)",
    `基因列表`              = "基因列表（含 gene_effect_sign）"
  )

doc <- officer::read_docx()
doc <- officer::body_add_par(
  doc,
  "Figure 7B. EphB1 下游信号模块（用于中介分析与机制图）",
  style = "heading 1"
)
doc <- officer::body_add_par(doc, "", style = "Normal")
doc <- flextable::body_add_flextable(doc, ft)

print(doc, target = docx_path)

message("[DONE] 已生成：")
message("  - TSV : ", tsv_path)
message("  - Word: ", docx_path)