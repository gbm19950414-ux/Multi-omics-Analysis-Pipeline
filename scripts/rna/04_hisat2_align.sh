#!/usr/bin/env bash
# === HISAT2 比对（低内存、可续跑、临时解压到 SSD）===
# 输入：data/interim/rna/<batch>/trimmed/*_R1.trim.fastq.gz 及 *_R2.trim.fastq.gz
# 输出：data/interim/rna/<batch>/hisat2/<sample>/Aligned.sortedByCoord.out.bam (+ .bai)
# 注：保持与 STAR 产物命名尽量一致，便于沿用后续 05_featureCounts.sh

set -euo pipefail

# -------- 基线路径（与你当前保持一致） --------
ROOT="/Volumes/Samsung_SSD_990_PRO_2TB_Media/multiomics_mech"
INTERIM="$ROOT/data/interim/rna"
REF="$ROOT/ref"
H2_IDX_DIR="$REF/hisat2_index"
H2_BASE="$H2_IDX_DIR/grcm39_primary"

# 从 INDEX.info 读取 splice/exon 文件（若存在）
SS="$H2_IDX_DIR/splice_sites.txt"
EX="$H2_IDX_DIR/exons.txt"

# -------- 运行参数（可用环境变量覆盖） --------
THREADS="${THREADS:-2}"                 # hisat2 线程
SAMTOOLS_THREADS="${SAMTOOLS_THREADS:-1}"
SAMTOOLS_SORT_MEM="${SAMTOOLS_SORT_MEM:-400M}"
READS_N="${READS_N:-0}"                 # 仅调试：0=不限制；>0 截断
TMPDIR_SSD="$ROOT/tmp_hisat2_sys"
export TMPDIR="${TMPDIR_OVERRIDE:-$TMPDIR_SSD}"

# 选择批次与样本
BATCHES=(batch1 batch2 batch3)
if [[ -n "${BATCH:-}" ]]; then BATCHES=("$BATCH"); fi
SAMPLE_GLOB="${SAMPLE_GLOB:-*}"

# -------- 依赖检查 --------
command -v hisat2 >/dev/null 2>&1 || { echo "[FATAL] 未找到 hisat2"; exit 10; }
command -v samtools >/dev/null 2>&1 || { echo "[FATAL] 未找到 samtools"; exit 11; }

# 仅强制要求 HISAT2 索引存在；splice site 文件可选
[[ -s "$H2_BASE.1.ht2" ]] || { echo "[FATAL] 缺少 HISAT2 索引: $H2_BASE.1.ht2"; exit 12; }
SS_OPT=()
if [[ -s "$SS" ]]; then
  SS_OPT=(--known-splicesite-infile "$SS")
fi

mkdir -p "$TMPDIR"; chmod 1777 "$TMPDIR" || true

echo "== Space check =="
df -h "$ROOT" "$TMPDIR" | sed -n '1,5p'
echo "[PARAM] THREADS=$THREADS SAMTOOLS_THREADS=$SAMTOOLS_THREADS SAMTOOLS_SORT_MEM=$SAMTOOLS_SORT_MEM READS_N=$READS_N TMPDIR=$TMPDIR"

for B in "${BATCHES[@]}"; do
  IN="$INTERIM/$B/trimmed"
  OUTROOT="$INTERIM/$B/hisat2"
  mkdir -p "$OUTROOT"
  echo "== RUN $B =="

  shopt -s nullglob
  for R1 in "$IN"/${SAMPLE_GLOB}_R1.trim.fastq.gz; do
    fname="$(basename "$R1")"
    sample="${fname%_R1.trim.fastq.gz}"
    R2="$IN/${sample}_R2.trim.fastq.gz"
    OUT="$OUTROOT/$sample"
    mkdir -p "$OUT"

    BAM="$OUT/Aligned.sortedByCoord.out.bam"
    if [[ -s "$BAM" && -s "$BAM.bai" ]]; then
      echo "  [SKIP] $B/$sample 已完成"; continue
    fi

    # ---- 临时解压到 SSD ----
    TMPD="$(mktemp -d "$TMPDIR/hisat2_${B}_${sample}.XXXXXX")"
    echo "  [UNZIP] $B/$sample -> $TMPD"
    /usr/bin/gzip -dc "$R1" > "$TMPD/R1.fastq"
    PE=0
    if [[ -f "$R2" ]]; then
      /usr/bin/gzip -dc "$R2" > "$TMPD/R2.fastq"; PE=1
    fi

    # ---- 组织 hisat2 命令 ----
    # --known-splicesite-infile（可选）: 使用已知剪接位点（存在才加）
    # 低内存：不在内存中排序；用 samtools 外排排序
    echo "  [HISAT2] $B/$sample (PE=$PE)"
    if [[ "$READS_N" -gt 0 ]]; then
      head -n $((READS_N*4)) "$TMPD/R1.fastq" > "$TMPD/R1.sub.fastq"
      [[ "$PE" -eq 1 ]] && head -n $((READS_N*4)) "$TMPD/R2.fastq" > "$TMPD/R2.sub.fastq"
      IN1="$TMPD/R1.sub.fastq"; IN2="$TMPD/R2.sub.fastq"
    else
      IN1="$TMPD/R1.fastq";     IN2="$TMPD/R2.fastq"
    fi

    # 比对 → BAM → 排序
    if [[ "$PE" -eq 1 ]]; then
      hisat2 -p "$THREADS" --dta \
        "${SS_OPT[@]}" \
        -x "$H2_BASE" -1 "$IN1" -2 "$IN2" \
        --summary-file "$OUT/align.summary.txt" \
      | samtools view -@ "$SAMTOOLS_THREADS" -b - \
      | samtools sort -@ "$SAMTOOLS_THREADS" -m "$SAMTOOLS_SORT_MEM" -T "$TMPDIR/sort_${B}_${sample}" -o "$BAM" -
    else
      hisat2 -p "$THREADS" --dta \
        "${SS_OPT[@]}" \
        -x "$H2_BASE" -U "$IN1" \
        --summary-file "$OUT/align.summary.txt" \
      | samtools view -@ "$SAMTOOLS_THREADS" -b - \
      | samtools sort -@ "$SAMTOOLS_THREADS" -m "$SAMTOOLS_SORT_MEM" -T "$TMPDIR/sort_${B}_${sample}" -o "$BAM" -
    fi

    samtools index -@ "$SAMTOOLS_THREADS" "$BAM"
    [[ -s "$OUT/align.summary.txt" ]] || { echo "  [ERR] $B/$sample 对齐 summary 缺失（hisat2 可能失败）"; exit 21; }

    # 产物检查
    if [[ -s "$BAM" && -s "$BAM.bai" ]]; then
      echo "  [OK ] $B/$sample 完成"
      rm -rf "$TMPD"
    else
      echo "  [ERR] $B/$sample 产物缺失，保留临时目录：$TMPD"
    fi
  done
  shopt -u nullglob
done

echo "== DONE =="