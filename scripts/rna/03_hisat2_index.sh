#!/usr/bin/env bash
# === HISAT2 索引构建（小鼠 GRCm39）===
# 基于你现有的解压参考：/ref/unzipped/*.fa / *.gtf
# 输出到：/ref/hisat2_index

set -euo pipefail

REF_DIR="/Volumes/Samsung_SSD_990_PRO_2TB_Media/multiomics_mech/ref"
FA="$REF_DIR/unzipped/Mus_musculus.GRCm39.dna.primary_assembly.fa"
GTF="$REF_DIR/unzipped/Mus_musculus.GRCm39.109.gtf"
IDX_DIR="$REF_DIR/hisat2_index"
IDX_BASE="$IDX_DIR/grcm39_primary"

THREADS="${THREADS:-4}"

mkdir -p "$IDX_DIR"

# 基本检查
for f in "$FA" "$GTF"; do
  [[ -s "$f" ]] || { echo "[FATAL] 缺少参考文件: $f"; exit 2; }
done

echo "== 提取 splice sites / exons =="
SS="$IDX_DIR/splice_sites.txt"
EX="$IDX_DIR/exons.txt"
hisat2_extract_splice_sites.py "$GTF" > "$SS"
hisat2_extract_exons.py        "$GTF" > "$EX"

echo "== 构建索引 =="
hisat2-build -p "$THREADS" "$FA" "$IDX_BASE"

echo "== 索引完成，记录配置 =="
cat > "$IDX_DIR/INDEX.info" <<EOF
FA=$FA
GTF=$GTF
IDX_BASE=$IDX_BASE
SS=$SS
EX=$EX
DATE=$(date)
EOF

echo "[OK] HISAT2 索引在: $IDX_BASE.*"