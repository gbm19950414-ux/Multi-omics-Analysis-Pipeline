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

#!/usr/bin/env bash
set -euo pipefail

# 加载配置与批次
source "$(dirname "$0")/00_config.env"
BATCH="${BATCH:?请用 BATCH=batch1 指定批次}"
IN="$RAW/$BATCH"
TRIM="$INTERIM/$BATCH/trimmed"
QC="$QC_OUT/$BATCH/fastp"
mkdir -p "$TRIM" "$QC"

# 说明：
# 1) 严格以 R1 文件为驱动，样本基名保留编号（如 KO_1、WT_3）
# 2) 同时支持 *_R1.{fq,fastq}.gz 与 *_1.{fq,fastq}.gz，并跳过 ._* AppleDouble 隐藏文件
# 3) 防重复：不同通配模式命中同一 R1 时仅处理一次

shopt -s nullglob

# 收集 R1 列表（四种常见命名）
mapfile -t R1_LIST < <(\
  {
    ls -1 "$IN"/*_R1.fastq.gz 2>/dev/null; \
    ls -1 "$IN"/*_R1.fq.gz     2>/dev/null; \
    ls -1 "$IN"/*_1.fastq.gz   2>/dev/null; \
    ls -1 "$IN"/*_1.fq.gz      2>/dev/null; \
  } | awk '!/\/\._/' | sort -u)

if [[ ${#R1_LIST[@]} -eq 0 ]]; then
  echo "[ERR] 未在 $IN 找到 R1 FASTQ（支持 *_R1.{fq,fastq}.gz 或 *_1.{fq,fastq}.gz）。" >&2
  exit 1
fi

# 记录已处理样本，避免重复
declare -A seen

for R1 in "${R1_LIST[@]}"; do
  base=$(basename "$R1")
  # 计算样本基名：去掉扩展名与 R1/1 尾缀，保留编号（如 KO_1）
  sample=${base}
  sample=${sample%.fastq.gz}
  sample=${sample%.fq.gz}
  sample=$(sed -E 's/_R1(_001)?$//; s/_1(_001)?$//' <<<"$sample")

  # 防重复
  if [[ -n ${seen["$sample"]+x} ]]; then
    continue
  fi
  seen["$sample"]=1

  # 推断 R2：按原始文件名一一对应替换
  R2="$R1"
  R2=${R2/%_R1.fastq.gz/_R2.fastq.gz}
  R2=${R2/%_R1.fq.gz/_R2.fq.gz}
  R2=${R2/%_1.fastq.gz/_2.fastq.gz}
  R2=${R2/%_1.fq.gz/_2.fq.gz}

  # 输出与报告路径（样本级，一对一）
  OUT1="$TRIM/${sample}_R1.trim.fastq.gz"
  OUT2="$TRIM/${sample}_R2.trim.fastq.gz"
  HTML="$QC/${sample}.fastp.html"
  JSON="$QC/${sample}.fastp.json"

  mkdir -p "$(dirname "$OUT1")" "$(dirname "$HTML")"

  if [[ -f "$R2" ]]; then
    echo "[INFO] fastp PE: $sample"
    fastp \
      --in1 "$R1" --in2 "$R2" \
      --out1 "$OUT1" --out2 "$OUT2" \
      --detect_adapter_for_pe \
      --thread "$THREADS" \
      --html "$HTML" --json "$JSON"
  else
    echo "[WARN] 未找到配对 R2，按单端处理: $base" >&2
    # 若确为单端数据，输出无 _R1/_R2 的统一命名；为保持一致也可保留 _R1
    OUT_SE="$TRIM/${sample}.trim.fastq.gz"
    fastp \
      --in1 "$R1" \
      --out1 "$OUT_SE" \
      --thread "$THREADS" \
      --html "$HTML" --json "$JSON"
  fi

done