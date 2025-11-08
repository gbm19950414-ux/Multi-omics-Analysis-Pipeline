#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/00_config.env"
BATCH="${BATCH:?请用 BATCH=batch1 指定批次}"
IN="$RAW/$BATCH"
OUT="$QC_OUT/$BATCH/fastqc_raw"
mkdir -p "$OUT"

shopt -s nullglob
files=()
for pat in "${R1_PATTERNS[@]}"; do
  for r1 in "$IN"/$pat; do
    r2="${r1/_R1/_R2}"; r2="${r2/_1/_2}"
    files+=("$r1")
    [[ -f "$r2" ]] && files+=("$r2")
  done
done

[[ "${#files[@]}" -gt 0 ]] || { echo "[ERR] 未找到 $IN 下的 FASTQ"; exit 1; }
echo "[INFO] FastQC -> $OUT (${#files[@]} files)"
fastqc -t "$THREADS" -o "$OUT" "${files[@]}"
