#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/00_config.env"
BATCH="${BATCH:?请用 BATCH=batch1 指定批次}"
BAMROOT="$INTERIM/$BATCH/star"
OUTDIR="$PROCESSED/$BATCH/counts"
mkdir -p "$OUTDIR"

[[ -s "$GTF" ]] || { echo "[ERR] GTF 未找到，请检查 00_config.env"; exit 1; }

mapfile -t BAMS < <(find "$BAMROOT" -type f -name "Aligned.sortedByCoord.out.bam" | sort)
[[ "${#BAMS[@]}" -gt 0 ]] || { echo "[ERR] 未找到 BAM：$BAMROOT"; exit 1; }

echo "[INFO] featureCounts on ${#BAMS[@]} BAMs"
featureCounts -T "$THREADS" \
  -a "$GTF" \
  -o "$OUTDIR/featureCounts.txt" \
  -p -B -C \
  -s "$STRANDNESS" \
  "${BAMS[@]}"

# 导出纯计数矩阵（gene_id + 每样本 counts）
awk 'BEGIN{FS=OFS="\t"} NR==1{for(i=1;i<=6;i++)hdr[i]=$i; next}
     NR==2{printf "GeneID"; for(i=7;i<=NF;i++){printf "\t%s",$i}; printf "\n"; next}
     NR>2{printf "%s",$1; for(i=7;i<=NF;i++) printf "\t%s",$i; printf "\n"}' \
     "$OUTDIR/featureCounts.txt" > "$OUTDIR/featureCounts_matrix.tsv"

echo "[OK] $OUTDIR/featureCounts_matrix.tsv"
