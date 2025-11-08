#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/00_config.env"
STAR_INDEX="$(cd "$STAR_INDEX" && pwd)"
BATCH="${BATCH:?请用 BATCH=batch1 指定批次}"

IN="$INTERIM/$BATCH/trimmed"
OUTROOT="$INTERIM/$BATCH/star"
mkdir -p "$OUTROOT"

[[ -d "$STAR_INDEX" ]] || { echo "[ERR] 未检测到 STAR 索引，请先运行 03_star_index.sh"; exit 1; }
echo "[DEBUG] genomeDir=[$STAR_INDEX]"
if [[ -z "${STAR_INDEX:-}" ]]; then
  echo "[ERR] STAR_INDEX 为空"; exit 2
fi
if [[ ! -r "$STAR_INDEX/geneInfo.tab" ]]; then
  echo "[ERR] geneInfo.tab 不可读: $STAR_INDEX/geneInfo.tab"; exit 3
fi
# 使用本机 APFS 的 /tmp 作为 STAR 临时目录与 FIFO 目录，避免外置盘不支持 FIFO
APFS_TMP_BASE="/tmp/multiomics_mech/star_${BATCH}"
mkdir -p "$APFS_TMP_BASE"

shopt -s nullglob
for R1 in "$IN"/*_R1.trim.fastq.gz; do
  [[ -e "$R1" ]] || continue
  sample="$(basename "$R1" | sed -E 's/_R1\.trim\.fastq\.gz$//')"
  R2="$IN/${sample}_R2.trim.fastq.gz"
  OUT="$OUTROOT/$sample"; mkdir -p "$OUT"

  # 每个样本独立的 STAR 临时目录；STAR 要求该目录在运行前不存在
  STAR_TMP="$APFS_TMP_BASE/tmp.${sample}"
  [[ -d "$STAR_TMP" ]] && rm -rf "$STAR_TMP"

  # 使用命名管道(FIFO)进行流式解压，彻底绕开 --readFilesCommand 的 exec 问题
  FIFO_DIR="$APFS_TMP_BASE/fifo.${sample}"
  rm -rf "$FIFO_DIR"; mkdir -p "$FIFO_DIR"

  if [[ -f "$R2" ]]; then
    echo "[INFO] STAR PE (FIFO): $sample"
    FIFO_R1="$FIFO_DIR/${sample}_R1.fifo"
    FIFO_R2="$FIFO_DIR/${sample}_R2.fifo"
    mkfifo "$FIFO_R1" "$FIFO_R2"

    # 清理函数（无论成功失败都清理 FIFO 与目录）
    cleanup() { rm -f "$FIFO_R1" "$FIFO_R2"; rm -rf "$FIFO_DIR"; }

    # 后台解压写入 FIFO
    /usr/bin/gzip -dc "$R1" > "$FIFO_R1" & pid1=$!
    /usr/bin/gzip -dc "$R2" > "$FIFO_R2" & pid2=$!

    set +e
    set -x  # 在调用 STAR 之前打开
    STAR --runThreadN "$THREADS" \
      --runMode alignReads \
      --genomeDir "$STAR_INDEX" \
      ...其余参数...
    set +x  # 紧跟在 STAR 之后关闭
    STAR --runThreadN "$THREADS" \
         --runMode alignReads \
         --genomeDir "$STAR_INDEX" \
         --readFilesIn "$FIFO_R1" "$FIFO_R2" \
         --outTmpDir "$STAR_TMP" \
         --outFileNamePrefix "$OUT/" \
         --outSAMtype BAM SortedByCoordinate \
         --quantMode GeneCounts
    status=$?
    # 等待解压进程结束（容忍其先于 STAR 退出）
    wait $pid1 2>/dev/null || true
    wait $pid2 2>/dev/null || true
    cleanup
    set -e
    if [[ $status -ne 0 ]]; then
      echo "[ERR] STAR 对齐失败: $sample" >&2
      exit $status
    fi
  else
    echo "[INFO] STAR SE (FIFO): $sample"
    FIFO_R1="$FIFO_DIR/${sample}_R1.fifo"
    mkfifo "$FIFO_R1"

    cleanup() { rm -f "$FIFO_R1"; rm -rf "$FIFO_DIR"; }

    /usr/bin/gzip -dc "$R1" > "$FIFO_R1" & pid1=$!

    set +e
    set -x  # 在调用 STAR 之前打开
    STAR --runThreadN "$THREADS" \
      --runMode alignReads \
      --genomeDir "$STAR_INDEX" \
      ...其余参数...
    set +x  # 紧跟在 STAR 之后关闭
    STAR --runThreadN "$THREADS" \
         --runMode alignReads \
         --genomeDir "$STAR_INDEX" \
         --readFilesIn "$FIFO_R1" \
         --outTmpDir "$STAR_TMP" \
         --outFileNamePrefix "$OUT/" \
         --outSAMtype BAM SortedByCoordinate \
         --quantMode GeneCounts
    status=$?
    wait $pid1 2>/dev/null || true
    cleanup
    set -e
    if [[ $status -ne 0 ]]; then
      echo "[ERR] STAR 对齐失败: $sample" >&2
      exit $status
    fi
  fi

done
