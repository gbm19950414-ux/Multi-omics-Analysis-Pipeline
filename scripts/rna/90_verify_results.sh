#!/usr/bin/env bash
# Minimal verifier for RNA-seq preprocessing & alignment
# Checks: references & indexes -> raw fastq presence -> FastQC -> fastp -> alignment BAMs
# Usage:
#   bash 90_verify_results.sh               # verify default batches (batch1-3)
#   BATCH=batch2 bash 90_verify_results.sh  # verify one batch
# Exit code: 0 (pass, possibly with WARN), 1 (fail)

set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
source "$HERE/00_config.env"

ok()   { printf "[OK ] %s\n"   "$*"; }
warn() { printf "[WARN] %s\n"  "$*" >&2; }
fail() { printf "[FAIL] %s\n"  "$*" >&2; FAIL=1; }
say()  { printf "%s\n" "$*"; }

FAIL=0
REPORT_DIR="${QC_OUT%/}/_verify_simple"
mkdir -p "$REPORT_DIR"

# ---------------------------------------------
# 0) Batches
# ---------------------------------------------
if [[ -n "${BATCH:-}" ]]; then
  batches=("$BATCH")
else
  batches=(batch1 batch2 batch3)
fi
say "== Batches: ${batches[*]} =="

# ---------------------------------------------
# 1) References & Indexes (STAR and/or HISAT2)
# ---------------------------------------------
say "\n== References & Indexes =="
for f in "$GENOME_FA" "$GTF"; do
  if [[ -s "$f" ]]; then ok "Found compressed ref: $f"; else warn "Missing ref: $f"; fi
done

# Uncompressed (optional)
FA_UNZ="${GENOME_FA%.gz}"; GTF_UNZ="${GTF%.gz}"
[[ -s "$FA_UNZ"  ]] && ok "Uncompressed fasta present: $FA_UNZ"  || warn "Uncompressed fasta not found (ok)"
[[ -s "$GTF_UNZ" ]] && ok "Uncompressed gtf present: $GTF_UNZ"    || warn "Uncompressed gtf not found (ok)"

# STAR index (optional)
if [[ -d "$STAR_INDEX" ]]; then
  need=0
  for k in Genome SA SAindex; do
    if [[ -s "$STAR_INDEX/$k" ]]; then ok "STAR index: $k"; else need=1; warn "STAR index missing: $STAR_INDEX/$k"; fi
  done
  if [[ -f "$STAR_INDEX/genomeParameters.txt" ]]; then ok "STAR genomeParameters.txt present"; fi
  (( need==1 )) && FAIL=1
else
  warn "STAR index dir not found (ok if using HISAT2): $STAR_INDEX"
fi

# HISAT2 index (optional): look under $REF_DIR/hisat2_index/*.ht2
H2_DIR="${REF_DIR%/}/hisat2_index"
if compgen -G "$H2_DIR/*.ht2" > /dev/null; then
  cnt=$(ls "$H2_DIR"/*.ht2 2>/dev/null | wc -l | tr -d ' ')
  ok "HISAT2 index files: $cnt (*.ht2) in $H2_DIR"
else
  warn "HISAT2 index not found in $H2_DIR (ok if using STAR)"
fi

# ---------------------------------------------
# 2) Raw FASTQ presence & pairing
# ---------------------------------------------
say "\n== Raw FASTQ presence & pairing =="
for B in "${batches[@]}"; do
  IN="${RAW%/}/$B"
  if [[ ! -d "$IN" ]]; then warn "Batch dir missing: $IN"; FAIL=1; continue; fi
  shopt -s nullglob
  r1=("$IN"/*_R1*.f*q.gz "$IN"/*_1*.f*q.gz)
  r2=("$IN"/*_R2*.f*q.gz "$IN"/*_2*.f*q.gz)
  shopt -u nullglob
  say "  $B: R1=${#r1[@]}  R2=${#r2[@]}"
  (( ${#r1[@]}==0 )) && { warn "No R1 files in $IN"; FAIL=1; }
  (( ${#r1[@]} != ${#r2[@]} )) && { warn "R1/R2 count mismatch in $IN"; FAIL=1; }

done

# ---------------------------------------------
# 3) FastQC artifacts (zip+html) & FAIL/WARN counts
# ---------------------------------------------
say "\n== FastQC reports =="
echo -e "batch\tsample\tfail\twarn" > "$REPORT_DIR/fastqc_summary.tsv"
for B in "${batches[@]}"; do
  OUT="${QC_OUT%/}/$B/fastqc_raw"
  [[ -d "$OUT" ]] || { warn "FastQC dir missing: $OUT"; FAIL=1; continue; }
  shopt -s nullglob
  for zipf in "$OUT"/*_fastqc.zip; do
    stem="$(basename "$zipf" "_fastqc.zip")"
    if ! unzip -p "$zipf" */summary.txt 2>/dev/null | tr -d '\000' > /tmp/_fastqc_sum.$$; then
      warn "summary.txt missing or unreadable in $zipf"
      continue
    fi
    failc=$(awk '$1=="FAIL"{c++} END{print c+0}' /tmp/_fastqc_sum.$$)
    warnc=$(awk '$1=="WARN"{c++} END{print c+0}' /tmp/_fastqc_sum.$$)
    rm -f /tmp/_fastqc_sum.$$
    echo -e "$B\t$stem\t$failc\t$warnc" >> "$REPORT_DIR/fastqc_summary.tsv"
    [[ -f "$OUT/${stem}_fastqc.html" ]] || warn "HTML missing for $stem"
  done
  shopt -u nullglob
 done
ok "Wrote $REPORT_DIR/fastqc_summary.tsv"

# ---------------------------------------------
# 4) fastp outputs (trimmed pairs + per-sample json/html)
# ---------------------------------------------
say "\n== fastp outputs & trimmed fq =="
echo -e "batch\tsample\tR1.trim\tR2.trim\tjson\thtml" > "$REPORT_DIR/fastp_files.tsv"
for B in "${batches[@]}"; do
  IN="${RAW%/}/$B"; TRIM="${INTERIM%/}/$B/trimmed"; FPOUT="${QC_OUT%/}/$B/fastp"
  [[ -d "$TRIM" ]] || { warn "Trim dir missing: $TRIM"; FAIL=1; continue; }
  shopt -s nullglob
  for r1 in "$IN"/*_R1*.f*q.gz "$IN"/*_1*.f*q.gz; do
    base=$(basename "$r1"); stem="${base%.fastq.gz}"; stem="${stem%.fq.gz}"
    sample=$(echo "$stem" | sed -E 's/(_R?[12])$//')
    o1="$TRIM/${sample}_R1.trim.fastq.gz"; o2="$TRIM/${sample}_R2.trim.fastq.gz"
    j="$FPOUT/${sample}.fastp.json"; h="$FPOUT/${sample}.fastp.html"
    status=("MISS" "MISS" "MISS" "MISS")
    [[ -s "$o1" ]] && { gzip -t "$o1" 2>/dev/null || true; status[0]="OK"; } || warn "missing $o1"
    [[ -s "$o2" ]] && { gzip -t "$o2" 2>/dev/null || true; status[1]="OK"; } || warn "missing $o2"
    [[ -f "$j"  ]] && status[2]="OK" || warn "missing $j"
    [[ -f "$h"  ]] && status[3]="OK" || warn "missing $h"
    echo -e "$B\t$sample\t${status[0]}\t${status[1]}\t${status[2]}\t${status[3]}" >> "$REPORT_DIR/fastp_files.tsv"
  done
  shopt -u nullglob
 done
ok "Wrote $REPORT_DIR/fastp_files.tsv"

# ---------------------------------------------
# 5) Alignment BAMs (STAR and/or HISAT2) + flagstat
# ---------------------------------------------
say "\n== Alignment BAMs & flagstat =="
echo -e "batch\taligner\tsample\tbam\tmapped%\tproperly_paired%" > "$REPORT_DIR/aln_flagstat.tsv"
for B in "${batches[@]}"; do
  for ALG in star hisat2; do
    ROOT="${INTERIM%/}/$B/$ALG"
    [[ -d "$ROOT" ]] || continue
    while IFS= read -r -d '' bam; do
      sample=$(basename "$(dirname "$bam")")
      idx="$bam.bai"; [[ -f "$idx" ]] || samtools index -@ 1 "$bam" >/dev/null 2>&1 || true
      fs=$(samtools flagstat -@ 1 "$bam" 2>/dev/null || true)
      # parse mapped and properly paired
      mapP=$(awk '/\(%: mapped\)/{print $5}' <<< "$fs" | tr -d '()%')
      propP=$(awk '/properly paired/{print $5}' <<< "$fs" | tr -d '()%')
      echo -e "$B\t$ALG\t$sample\t$(basename "$bam")\t${mapP:-NA}\t${propP:-NA}" >> "$REPORT_DIR/aln_flagstat.tsv"
      ok "$B/$ALG/$sample : BAM checked"
    done < <(find "$ROOT" -type f -name 'Aligned.sortedByCoord.out.bam' -print0)
  done
 done
ok "Wrote $REPORT_DIR/aln_flagstat.tsv"

# ---------------------------------------------
# Done
# ---------------------------------------------
if (( FAIL==0 )); then
  say "\n[SUMMARY] PASS (with possible WARN). Reports: $REPORT_DIR"
  exit 0
else
  say "\n[SUMMARY] FAIL detected. Check: $REPORT_DIR"
  exit 1
fi
