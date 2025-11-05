#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/00_config.env"
BATCH="${BATCH:?请用 BATCH=batch1 指定批次}"
IN="$RAW/$BATCH"
TRIM="$INTERIM/$BATCH/trimmed"
QC="$QC_OUT/$BATCH/fastp"
mkdir -p "$TRIM" "$QC"

shopt -s nullglob
found=0
for pat in "${R1_PATTERNS[@]}"; do
  for R1 in "$IN"/$pat; do
    R2="${R1/_R1/_R2}"; R2="${R2/_1/_2}"
    # 样本名去掉 _R1/_1 和扩展名
    base="$(basename "$R1")"
    sample="${base%.*}"; sample="${sample%.*}"  # 去 .gz/.fastq
    sample="$(echo "$sample" | sed -E 's/_R1(_001)?$//; s/_1$//')"

    if [[ -f "$R2" ]]; then
      echo "[INFO] fastp PE: $sample"
      fastp --in1 "$R1" --in2 "$R2" \
            --out1 "$TRIM/${sample}_R1.trim.fastq.gz" \
            --out2 "$TRIM/${sample}_R2.trim.fastq.gz" \
            --detect_adapter_for_pe \
            --thread "$THREADS" \
            --html "$QC/${sample}.fastp.html" \
            --json "$QC/${sample}.fastp.json"
    else
      echo "[INFO] fastp SE: $sample"
      fastp --in1 "$R1" \
            --out1 "$TRIM/${sample}.trim.fastq.gz" \
            --thread "$THREADS" \
            --html "$QC/${sample}.fastp.html" \
            --json "$QC/${sample}.fastp.json"
    fi
    found=1
  done
done
[[ "$found" -eq 1 ]] || { echo "[ERR] 未找到 $IN 下的 FASTQ"; exit 1; }
