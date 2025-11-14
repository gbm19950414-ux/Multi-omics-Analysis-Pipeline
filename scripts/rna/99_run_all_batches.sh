#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/00_config.env"

# 参考索引（若已构建会跳过）
bash "$(dirname "$0")/03_star_index.sh"

# 小工具：存在且非空
exists_nonempty() {
  local p="$1"
  [[ -s "$p" ]] && return 0 || return 1
}

# 小工具：目录存在且包含至少一个文件
nonempty_dir() {
  local d="$1"
  [[ -d "$d" ]] && find "$d" -mindepth 1 -maxdepth 1 -type f | grep -q .
}

# 自愈：清理 0 字节 BAM 与残留临时目录（避免外置盘/中断导致的假完成）
clean_star_residuals() {
  local BATCH_NAME="$1"
  local STAR_DIR="$INTERIM/$BATCH_NAME/star"
  # 删 0B BAM 的父目录（样本目录）
  if find "$STAR_DIR" -type f -name 'Aligned.sortedByCoord.out.bam' -size 0 -print 2>/dev/null | grep -q .; then
    echo "[FIX ] 检测到 0B BAM，清理对应样本目录..."
    while IFS= read -r f; do
      echo "       - rm -rf \"$(dirname "$f")\""
      rm -rf "$(dirname "$f")"
    done < <(find "$STAR_DIR" -type f -name 'Aligned.sortedByCoord.out.bam' -size 0 -print)
  fi
  # 删 /tmp 下残留的 STAR 临时目录
  local TMP_BASE="/tmp/multiomics_mech/star_${BATCH_NAME}"
  if [[ -d "$TMP_BASE" ]]; then
    echo "[FIX ] 清理 /tmp 残留: $TMP_BASE"
    rm -rf "$TMP_BASE"/tmp.* || true
  fi
}

# 预检：确保 readFilesCommand 包装脚本存在且可执行（04_star_align.sh 依赖它）
ensure_read_cmd() {
  local rf="/tmp/read_fastq.sh"
  if [[ ! -x "$rf" ]]; then
    cat >&2 <<EOF
[ERR] 未找到可执行的 /tmp/read_fastq.sh。请先创建：
  cat > /tmp/read_fastq.sh <<'EOT'
#!/bin/sh
exec /usr/bin/gzip -dc "$@"
EOT
  chmod +x /tmp/read_fastq.sh
然后重新运行本脚本。
EOF
    exit 1
  fi
}

ensure_read_cmd

# 额外健壮性：警告脚本中是否混入 CR（\r）
if LC_ALL=C grep -q $'\r' "$0"; then
  echo "[WARN] 检测到本脚本含有 CR(\r) 行尾，可能导致参数异常。建议执行： perl -i -pe 's/\r$//' \"$0\"" >&2
fi

for B in "${BATCHES[@]}"; do
  echo "======== $B ========"

  # 1) FastQC（原始）
  RAW_QC_DIR="$QC_OUT/$B/fastqc_raw"
  if nonempty_dir "$RAW_QC_DIR"; then
    echo "[SKIP] FastQC raw: 已存在 -> $RAW_QC_DIR"
  else
    echo "[RUN ] FastQC raw"
    BATCH="$B" bash "$(dirname "$0")/01_fastqc.sh"
  fi

  # 2) fastp（清洗）
  TRIM_DIR="$INTERIM/$B/trimmed"
  if compgen -G "$TRIM_DIR/*_R?.trim.fastq.gz" > /dev/null 2>&1; then
    echo "[SKIP] fastp: 已存在 -> $TRIM_DIR"
  else
    echo "[RUN ] fastp"
    BATCH="$B" bash "$(dirname "$0")/02_fastp.sh"
  fi

  # 3) STAR 比对（先做自愈清理，再判定是否需要运行）
  clean_star_residuals "$B"
  STAR_DIR="$INTERIM/$B/star"
  if find "$STAR_DIR" -type f -name "Aligned.sortedByCoord.out.bam" ! -size 0 -print 2>/dev/null | grep -q .; then
    echo "[SKIP] STAR: 已存在 BAM -> $STAR_DIR"
  else
    echo "[RUN ] STAR align"
    BATCH="$B" bash "$(dirname "$0")/04_star_align.sh"
  fi

  # 4) featureCounts 计数
  FC_DIR="$PROCESSED/$B/counts"
  FC_MAT="$FC_DIR/featureCounts_matrix.tsv"
  if exists_nonempty "$FC_MAT"; then
    echo "[SKIP] featureCounts: 已存在 -> $FC_MAT"
  else
    echo "[RUN ] featureCounts"
    BATCH="$B" bash "$(dirname "$0")/05_featureCounts.sh"
  fi

  # 5) MultiQC 汇总
  MQC_HTML="$QC_OUT/$B/multiqc/index.html"
  if exists_nonempty "$MQC_HTML"; then
    echo "[SKIP] MultiQC: 已存在 -> $MQC_HTML"
  else
    echo "[RUN ] MultiQC"
    BATCH="$B" bash "$(dirname "$0")/06_multiqc.sh"
  fi

done

# 6) 合并三个批次的计数矩阵（若不存在则生成）
OUT_MERGED="$PROCESSED/featureCounts_merged_all_batches.tsv"
if exists_nonempty "$OUT_MERGED"; then
  echo "[SKIP] 合并矩阵已存在 -> $OUT_MERGED"
else
  echo "[RUN ] 合并三个批次计数矩阵"
  Rscript "$(dirname "$0")/07_merge_counts.R" "$PROCESSED" "$OUT_MERGED"
fi

echo "[DONE] 流水线执行完成。合并矩阵：$OUT_MERGED"
