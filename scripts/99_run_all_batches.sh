#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/00_config.env"

# 参考索引（若已构建会跳过）
bash "$(dirname "$0")/03_star_index.sh"

for B in "${BATCHES[@]}"; do
  echo "======== $B ========"
  BATCH="$B" bash "$(dirname "$0")/01_fastqc.sh"
  BATCH="$B" bash "$(dirname "$0")/02_fastp.sh"
  BATCH="$B" bash "$(dirname "$0")/04_star_align.sh"
  BATCH="$B" bash "$(dirname "$0")/05_featureCounts.sh"
  BATCH="$B" bash "$(dirname "$0")/06_multiqc.sh"
done

# 合并三个批次的计数矩阵
OUT_MERGED="$PROCESSED/featureCounts_merged_all_batches.tsv"
Rscript "$(dirname "$0")/07_merge_counts.R" "$PROCESSED" "$OUT_MERGED"

echo "[DONE] 全部完成。合并矩阵：$OUT_MERGED"
