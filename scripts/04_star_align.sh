#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/00_config.env"
BATCH="${BATCH:?请用 BATCH=batch1 指定批次}"
IN="$INTERIM/$BATCH/trimmed"
OUTROOT="$INTERIM/$BATCH/star"
mkdir -p "$OUTROOT"

[[ -d "$STAR_INDEX" ]] || { echo "[ERR] 未检测到 STAR 索引，请先运行 03_star_index.sh"; exit 1; }

shopt -s nullglob
for R1 in "$IN"/*_R1.trim.fastq.gz "$IN"/*.trim.fastq.gz; do
  [[ -e "$R1" ]] || continue
  sample="$(basename "$R1" | sed -E 's/_R1\.trim\.fastq\.gz$//; s/\.trim\.fastq\.gz$//')"
  R2="$IN/${sample}_R2.trim.fastq.gz"
  OUT="$OUTROOT/$sample"; mkdir -p "$OUT"

  if [[ -f "$R2" ]]; then
    echo "[INFO] STAR PE: $sample"
    STAR --runThreadN "$THREADS" \
         --genomeDir "$STAR_INDEX" \
         --readFilesIn "$R1" "$R2" \
         --readFilesCommand zcat \
         --outFileNamePrefix "$OUT/" \
         --outSAMtype BAM SortedByCoordinate \
         --quantMode GeneCounts
  else
    echo "[INFO] STAR SE: $sample"
    STAR --runThreadN "$THREADS" \
         --genomeDir "$STAR_INDEX" \
         --readFilesIn "$R1" \
         --readFilesCommand zcat \
         --outFileNamePrefix "$OUT/" \
         --outSAMtype BAM SortedByCoordinate \
         --quantMode GeneCounts
  fi
done
