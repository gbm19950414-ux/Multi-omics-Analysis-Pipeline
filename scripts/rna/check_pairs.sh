#!/usr/bin/env bash
set -euo pipefail

RAW="/Volumes/MyPassport/multiomics_mech/data/raw/rna"
for B in batch1 batch2 batch3; do
  echo "=== checking $B ==="
  # 收集 R1/R2：同时匹配 *_R1*.f*q.gz 和 *_1.f*q.gz（兼容 .fq.gz/.fastq.gz）
  mapfile -t R1S < <(find "$RAW/$B" -type f \( -name "*_R1*.f*q.gz" -o -name "*_1.f*q.gz" \) | sort)
  mapfile -t R2S < <(find "$RAW/$B" -type f \( -name "*_R2*.f*q.gz" -o -name "*_2.f*q.gz" \) | sort)

  # 归一化样本基名：去掉 _R1/_R2 或 _1/_2 以及后缀（兼容 fastq/fq）
  norm_r1() { sed -E 's/_R1[^/]*\.f(ast)?q\.gz$//; s/_1[^/]*\.f(ast)?q\.gz$//'; }
  norm_r2() { sed -E 's/_R2[^/]*\.f(ast)?q\.gz$//; s/_2[^/]*\.f(ast)?q\.gz$//'; }

  comm -23 \
    <(printf "%s\n" "${R1S[@]}" | norm_r1 | sort -u) \
    <(printf "%s\n" "${R2S[@]}" | norm_r2 | sort -u) \
    | sed 's#^#[WARN] R2 missing for: #'

  comm -13 \
    <(printf "%s\n" "${R1S[@]}" | norm_r1 | sort -u) \
    <(printf "%s\n" "${R2S[@]}" | norm_r2 | sort -u) \
    | sed 's#^#[WARN] R1 missing for: #'

  echo "[SUMMARY] $B -> R1: ${#R1S[@]}  R2: ${#R2S[@]}"
done
