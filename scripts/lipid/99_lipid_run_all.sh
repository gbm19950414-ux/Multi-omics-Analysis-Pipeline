#!/usr/bin/env bash
set -euo pipefail

CFG="scripts/lipid/00_lipid_config.yaml"

echo "==== [01] 合并 ===="
python scripts/lipid/01_combine_lipid.py "$CFG"

echo "==== [02] 去空白/置信度/透视 ===="
python scripts/lipid/02_blank_confident_pivot.py "$CFG"

echo "==== [03] 缺失填补 + 归一化 + 变换 + 标准化 ===="
python scripts/lipid/03_impute_normalize_scale.py "$CFG"

echo "==== [04] 批次校正（ComBat） ===="
P=$(awk -F': *' '/^project_dir:/ {print $2}' "$CFG" | tr -d '"')
OUT_DIR=$(awk -F': *' '/^out_dir:/ {print $2}' "$CFG" | tr -d '"')
IN="$P/$OUT_DIR/03_imputed_normalized_scaled.tsv"
META=$(awk -F': *' '/^metadata_file:/ {print $2}' "$CFG" | tr -d '"')
OUT="$P/$OUT_DIR/04_batch_corrected.tsv"
Rscript scripts/lipid/04_batch_correct_ComBat.R --cfg "$CFG" --in "$IN" --meta "$P/$META" --out "$OUT"

echo "==== [05] CL 指标 ===="
python scripts/lipid/05_cl_metrics.py "$CFG"

echo "==== ALL DONE ===="