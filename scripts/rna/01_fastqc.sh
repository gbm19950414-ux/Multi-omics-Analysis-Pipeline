#!/usr/bin/env bash
# FastQC runner with completeness + quality report
# Usage:
#   bash 01_fastqc.sh batch1                # 跑该批次下所有 R1/R2
#   bash 01_fastqc.sh batch1 KO_1_R2.fq.gz  # 仅跑指定文件（可多给，支持相对/绝对路径）
# 也可用环境变量：BATCH=batch1 bash 01_fastqc.sh

set -euo pipefail

# --- load config ---
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "${SCRIPT_DIR}/00_config.env"

# --- parse args / env ---
ARG_BATCH="${1:-}"
if [[ -n "${ARG_BATCH}" && "${ARG_BATCH}" != -* && "${ARG_BATCH}" != *.fq.gz && "${ARG_BATCH}" != *.fastq.gz ]]; then
  BATCH="${ARG_BATCH}"
  shift || true
else
  : "${BATCH:?请用 'bash 01_fastqc.sh batch1' 或 'BATCH=batch1 bash 01_fastqc.sh' 指定批次}"
fi

IN="${RAW%/}/${BATCH}"
OUT="${QC_OUT%/}/${BATCH}/fastqc_raw"
mkdir -p "${OUT}"

# --- collect input fastq files ---
declare -a INPUTS=()

# If user provided explicit file names after batch, use them.
if [[ "$#" -gt 0 ]]; then
  while [[ "$#" -gt 0 ]]; do
    f="$1"; shift
    if [[ "$f" != /* ]]; then
      # treat as relative to IN
      f="${IN%/}/$f"
    fi
    if [[ -f "$f" ]]; then
      INPUTS+=("$f")
    else
      echo "[WARN] 指定文件不存在，跳过: $f" >&2
    fi
  done
else
  # Auto-discover: include both R1 and R2 (支持多种命名)
  shopt -s nullglob
  declare -A SEEN=()
  # 如果配置里提供了 R1_PATTERNS，则基于它们搜集，再附带推断的配对 R2
  if [[ "${#R1_PATTERNS[@]:-0}" -gt 0 ]]; then
    for pat in "${R1_PATTERNS[@]}"; do
      for r1 in "${IN}"/$pat; do
        if [[ -f "$r1" ]]; then
          SEEN["$r1"]=1
          r2="${r1/_R1/_R2}"; r2="${r2/_1/_2}"
          [[ -f "$r2" ]] && SEEN["$r2"]=1
        fi
      done
    done
  fi
  # 兜底：把目录里所有 *_R[12]* 与 *_[12]* 的 fq/fq.gz 都纳入
  for pat in "*_R"[12]"*.f*q.gz" "*_"[12]"*.f*q.gz"; do
    for f in "${IN}"/$pat; do
      [[ -f "$f" ]] && SEEN["$f"]=1
    done
  done
  # 转成数组
  for k in "${!SEEN[@]}"; do
    INPUTS+=("$k")
  done
  # 稳定排序
  IFS=$'\n' INPUTS=($(printf "%s\n" "${INPUTS[@]}" | sort)) ; unset IFS
fi

if [[ "${#INPUTS[@]}" -eq 0 ]]; then
  echo "[ERR] 未在 ${IN} 找到 FASTQ 文件" >&2
  exit 1
fi

echo "[INFO] FastQC 输入(${#INPUTS[@]} 个)："
printf '  - %s\n' "${INPUTS[@]}"

# --- run fastqc ---
echo "[INFO] FastQC 输出目录: ${OUT}"
fastqc -t "${THREADS}" -o "${OUT}" "${INPUTS[@]}"

# -------------------------------
# Completeness & Quality reporting
# -------------------------------

# 1) 原始 R1 -> R2 配对完整性检查（仅对自动发现/目录内文件）
echo "== 配对完整性检查 =="
set +e
missing_cnt=0
shopt -s nullglob
for r1 in "${IN}"/*_R1*.f*q.gz "${IN}"/*_1*.f*q.gz; do
  [[ -e "$r1" ]] || continue
  r2="${r1/_R1/_R2}"; r2="${r2/_1/_2}"
  if [[ ! -f "$r2" ]]; then
    echo "  [MISSING RAW R2] $r2"
    ((missing_cnt++))
  fi
done
if [[ "$missing_cnt" -eq 0 ]]; then
  echo "  [OK] 所有 R1 均有对应 R2"
fi

# 2) FastQC 产物存在性检查
echo "== 产物存在性检查 =="
prod_missing=0
declare -a ALL_BASES=()
for f in "${INPUTS[@]}"; do
  base="$(basename "$f")"
  stem="${base%.fastq.gz}"; stem="${stem%.fq.gz}"
  html="${OUT}/${stem}_fastqc.html"
  zipf="${OUT}/${stem}_fastqc.zip"
  ALL_BASES+=("${stem}")
  [[ -f "$html" ]] || { echo "  [MISSING HTML] $html"; ((prod_missing++)); }
  [[ -f "$zipf" ]] || { echo "  [MISSING ZIP ] $zipf"; ((prod_missing++)); }
done
[[ "$prod_missing" -eq 0 ]] && echo "  [OK] 所有 FastQC HTML/ZIP 均生成"

# 3) FastQC 质量统计（FAIL/WARN 数）
echo "== FastQC 质量统计（FAIL & WARN 计数）=="
summary_tsv="${OUT}/_fastqc_summary.tsv"
printf "sample\tFAIL\tWARN\n" > "$summary_tsv"
for stem in "${ALL_BASES[@]}"; do
  zipf="${OUT}/${stem}_fastqc.zip"
  if [[ -f "$zipf" ]]; then
    sum="$(unzip -p "$zipf" */summary.txt 2>/dev/null || true)"
    fail=$(awk '$1=="FAIL"{c++} END{print c+0}' <<< "$sum")
    warn=$(awk '$1=="WARN"{c++} END{print c+0}' <<< "$sum")
    printf "  [%s] FAIL=%s  WARN=%s\n" "$stem" "$fail" "$warn"
    printf "%s\t%s\t%s\n" "$stem" "$fail" "$warn" >> "$summary_tsv"
  else
    echo "  [$stem] (缺少 zip，无法统计)"
  fi
done
echo "[OK] 汇总表: $summary_tsv"

# 4) 专查 KO_1 / WT_1 的 R2
echo "== 专查 KO_1 / WT_1 的 R2 =="
for tag in KO_1 WT_1; do
  stem="${tag}_R2"
  html="${OUT}/${stem}_fastqc.html"
  zipf="${OUT}/${stem}_fastqc.zip"
  if [[ -f "$html" && -f "$zipf" ]]; then
    sum="$(unzip -p "$zipf" */summary.txt 2>/dev/null || true)"
    fail=$(awk '$1=="FAIL"{c++} END{print c+0}' <<< "$sum")
    warn=$(awk '$1=="WARN"{c++} END{print c+0}' <<< "$sum")
    echo "  [$stem] FAIL=$fail  WARN=$warn"
  else
    echo "  [$stem] 报告缺失 (HTML/ZIP 未找到)"
  fi
done

echo "[DONE] 01_fastqc 完成：批次 ${BATCH}"
