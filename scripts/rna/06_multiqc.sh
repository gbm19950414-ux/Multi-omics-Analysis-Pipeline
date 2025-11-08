#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/00_config.env"
BATCH="${BATCH:?请用 BATCH=batch1 指定批次}"
OUT="$QC_OUT/$BATCH/multiqc"
mkdir -p "$OUT"

TARGETS=(
  "$QC_OUT/$BATCH/fastqc_raw"
  "$QC_OUT/$BATCH/fastp"
  "$INTERIM/$BATCH/star"
  "$PROCESSED/$BATCH/counts"
)

echo "[INFO] MultiQC for batch: $BATCH"
multiqc -f -o "$OUT" "${TARGETS[@]}"
echo "[OK] MultiQC: $OUT/index.html"
