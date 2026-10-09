#!/usr/bin/env bash
set -euo pipefail

# === Config & args ===
source "$(dirname "$0")/00_config.env"
BATCH="${BATCH:?请用 BATCH=batch1 指定批次}"
OUT="$QC_OUT/$BATCH/multiqc"
mkdir -p "$OUT"

# 收集各步骤输出目录（含 HISAT2）
TARGETS=(
  "$QC_OUT/$BATCH/fastqc_raw"
  "$QC_OUT/$BATCH/fastp"
  "$INTERIM/$BATCH/star"
  "$INTERIM/$BATCH/hisat2"   # 新增：纳入 HISAT2 对齐摘要
  "$PROCESSED/$BATCH/counts"
)

echo "[INFO] MultiQC for batch: $BATCH"

# 过滤掉不存在或空目录，避免 MultiQC 噪声
valid_targets=()
for d in "${TARGETS[@]}"; do
  if [[ -d "$d" ]] && [[ $(find "$d" -type f -maxdepth 1 2>/dev/null | wc -l | tr -d ' ') -gt 0 ]] ; then
    echo "[OK  ] 包含目录: $d"
    valid_targets+=("$d")
  else
    echo "[WARN] 跳过空目录: $d"
  fi
done

if [[ ${#valid_targets[@]} -eq 0 ]]; then
  echo "[ERR ] 没有可用输入，退出" >&2
  exit 2
fi

# 给报告命名更清晰一些，并强制覆盖
multiqc -f -n "multiqc_${BATCH}" -o "$OUT" "${valid_targets[@]}"
echo "[OK ] MultiQC: $OUT/multiqc_${BATCH}.html"
