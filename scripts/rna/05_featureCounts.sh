#!/usr/bin/env bash
# featureCounts wrapper (STAR/HISAT2 auto-detect + compatibility + QC)
# Usage example:
#   BATCH=batch1 THREADS=4 STRANDNESS=0 bash scripts/rna/05_featureCounts.sh
# Notes:
#   - Prefers HISAT2 BAMs under $INTERIM/<batch>/hisat2, falls back to STAR under $INTERIM/<batch>/star
#   - Assumes paired-end; uses -p -B -C. Adjust if needed.

set -euo pipefail

# ---------- 0) Load config & basic params ----------
script_dir=$(cd -- "$(dirname -- "$0")" && pwd)
# shellcheck disable=SC1091
source "$script_dir/00_config.env"

: "${BATCH:?请用 BATCH=batch1 指定批次}"

THREADS="${THREADS:-4}"
STRANDNESS="${STRANDNESS:-0}"    # 0=unstranded, 1=forward, 2=reverse

# ---------- 1) Locate input BAMs (prefer HISAT2, else STAR) ----------
CAND_HISAT2="$INTERIM/$BATCH/hisat2"
CAND_STAR="$INTERIM/$BATCH/star"

if [[ -d "$CAND_HISAT2" ]]; then
  BAMROOT="$CAND_HISAT2"
  WHICH_ALIGNER="hisat2"
elif [[ -d "$CAND_STAR" ]]; then
  BAMROOT="$CAND_STAR"
  WHICH_ALIGNER="star"
else
  echo "[ERR] 未找到对齐目录：$CAND_HISAT2 或 $CAND_STAR" >&2
  exit 1
fi

OUTDIR="$PROCESSED/$BATCH/counts"
mkdir -p "$OUTDIR"

# GTF must exist (prefer uncompressed if you keep both)
if [[ -f "$GTF" ]]; then
  GTF_IN="$GTF"
elif [[ -f "$REF_DIR/unzipped/$(basename "${GTF%.gz}")" ]]; then
  GTF_IN="$REF_DIR/unzipped/$(basename "${GTF%.gz}")"
else
  echo "[ERR] GTF 未找到，请检查 00_config.env 的 GTF 路径" >&2
  exit 1
fi

# ---------- 2) Collect BAMs robustly (avoid mapfile dependency) ----------
# Prefer typical name from STAR/HISAT2 pipelines; fallback to all *.bam
collect_bams() {
  local root="$1"
  local -a list
  while IFS= read -r -d '' f; do list+=("$f"); done < <(find "$root" -type f -name "Aligned.sortedByCoord.out.bam" -print0 | sort -z)
  if (( ${#list[@]} == 0 )); then
    while IFS= read -r -d '' f; do list+=("$f"); done < <(find "$root" -type f -name "*.bam" -print0 | sort -z)
  fi
  printf '%s\n' "${list[@]}"
}

# Build array
BAMS=()
while IFS= read -r line; do
  [[ -n "$line" ]] && BAMS+=("$line")
done < <(collect_bams "$BAMROOT")

if (( ${#BAMS[@]} == 0 )); then
  echo "[ERR] 未找到 BAM 于 $BAMROOT" >&2
  exit 1
fi

# ---------- 3) Run featureCounts ----------
FC_TXT="$OUTDIR/featureCounts.txt"
FC_SUM="$OUTDIR/featureCounts.txt.summary"
FC_MAT="$OUTDIR/featureCounts_matrix.tsv"

echo "[INFO] Aligner=$WHICH_ALIGNER  Batch=$BATCH  BAM数=${#BAMS[@]}  Threads=$THREADS  Strandness=$STRANDNESS"

echo "[RUN ] featureCounts …"
featureCounts \
  -T "$THREADS" \
  -a "$GTF_IN" \
  -o "$FC_TXT" \
  -t exon -g gene_id \
  -p -B -C \
  -s "$STRANDNESS" \
  "${BAMS[@]}"

echo "[OK  ] featureCounts 完成：$FC_TXT"

# ---------- 4) Build compact counts matrix (GeneID + sample columns) ----------
# featureCounts.txt 的第1行为注释行，第2行为列名，从第7列起为样本
awk 'BEGIN{FS=OFS="\t"} \
     NR==1{next} \
     NR==2{printf "GeneID"; for(i=7;i<=NF;i++) printf "\t"$i; printf "\n"; next} \
     NR>2{printf $1; for(i=7;i<=NF;i++) printf "\t"$i; printf "\n"}' \
     "$FC_TXT" > "$FC_MAT"

echo "[OK  ] 计数矩阵：$FC_MAT"

# Optional: Clean column names to simple sample names if CLEAN_COLNAMES=1
if [[ "${CLEAN_COLNAMES:-0}" == "1" ]]; then
  # Extract header line
  header=$(head -n 1 "$FC_MAT")
  # Extract sample names from header (skip first "GeneID")
  IFS=$'\t' read -r -a cols <<< "$header"
  new_header="GeneID"
  for ((i=1; i<${#cols[@]}; i++)); do
    # Simplify sample name: remove path and extensions, keep basename without extensions
    sample="${cols[i]}"
    sample=$(basename "$sample")
    sample="${sample%%.*}"
    new_header+=$'\t'"$sample"
  done
  # Write cleaned header and rest of file
  { echo -e "$new_header"; tail -n +2 "$FC_MAT"; } > "$FC_MAT.tmp" && mv "$FC_MAT.tmp" "$FC_MAT"
  echo "[OK  ] 样本列名已清理：$FC_MAT"
fi

# ---------- 5) QC summary from *.summary ----------
# featureCounts 会生成 $FC_SUM，第一列是 Status，其余列是各样本对应数量。
if [[ -s "$FC_SUM" ]]; then
  # 提取 Assigned 与总数，计算 Assigned rate
  # 输出：sample  assigned  total  assigned_rate
  awk 'BEGIN{FS=OFS="\t"}
       NR==1{for(i=2;i<=NF;i++) s[i]=$i; next}
       {for(i=2;i<=NF;i++) A[i,$1]=$i}
       END{
         print "sample","assigned","total","assigned_rate"
         for(i=2;i in s; i++){
           assigned=A[i,"Assigned"]; if(assigned=="") assigned=0;
           total=0; for (k in A) if (k ~ "^"i"\034") total+=A[k];
           rate=(total>0)?assigned/total:0;
           printf "%s\t%d\t%d\t%.4f\n", s[i], assigned, total, rate
         }
       }' "$FC_SUM" > "$OUTDIR/qc_basic.tsv"
  echo "[OK  ] QC摘要：$OUTDIR/qc_basic.tsv"
else
  echo "[WARN] 未找到 $FC_SUM（不同版本featureCounts可能写到同名路径）。"
fi

# Optional: Automatic strandness detection if DETECT_STRANDNESS=1
if [[ "${DETECT_STRANDNESS:-0}" == "1" ]]; then
  echo "[INFO] 自动检测链特异性…"
  # Use first BAM for detection
  SAMPLE_BAM="${BAMS[0]}"
  # Run featureCounts with -s 2 (reverse) and -s 1 (forward) and compare assigned reads
  TMP_DETECT_DIR="$OUTDIR/strand_detect_tmp"
  mkdir -p "$TMP_DETECT_DIR"
  DETECT_RESULTS="$TMP_DETECT_DIR/detect_results.txt"
  for s in 0 1 2; do
    OUT_PREFIX="$TMP_DETECT_DIR/count_s${s}"
    featureCounts -T 2 -a "$GTF_IN" -o "${OUT_PREFIX}.txt" -t exon -g gene_id -p -B -C -s "$s" "$SAMPLE_BAM" > /dev/null 2>&1
    # Extract assigned reads count
    assigned=$(awk 'NR>2{sum+=$7}END{print sum}' "${OUT_PREFIX}.txt")
    echo -e "$s\t$assigned" >> "$DETECT_RESULTS"
  done
  # Determine strandness with max assigned
  best_strand=$(sort -k2,2nr "$DETECT_RESULTS" | head -n1 | cut -f1)
  case "$best_strand" in
    0) strand_desc="unstranded";;
    1) strand_desc="forward";;
    2) strand_desc="reverse";;
    *) strand_desc="unknown";;
  esac
  echo "[INFO] 检测结果：最佳链特异性参数 -s $best_strand ($strand_desc)"
  echo "[INFO] 您可据此调整 STRANDNESS 参数后重新运行本脚本。"
  rm -r "$TMP_DETECT_DIR"
fi

# ---------- 6) Console preview ----------
{
  echo "\n== 矩阵预览 =="
  head -n 5 "$FC_MAT" || true
  echo "\n== QC 摘要（前10行） =="
  head -n 10 "$OUTDIR/qc_basic.tsv" || true
} >&2

echo "\n[SUMMARY] Done. 输出目录：$OUTDIR"
