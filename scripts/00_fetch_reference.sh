#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/00_config.env"

mkdir -p "$REF_DIR"
cd "$REF_DIR"

# 1) 下载小鼠 GRCm39 主装配与注释（Ensembl）
# 如有现成文件，可跳过
if [[ ! -s "$GENOME_FA" ]]; then
  echo "[DL] Genome FASTA..."
  wget -c https://ftp.ensembl.org/pub/release-109/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz -O "$(basename "$GENOME_FA")"
fi
if [[ ! -s "$GTF" ]]; then
  echo "[DL] GTF..."
  wget -c https://ftp.ensembl.org/pub/release-109/gtf/mus_musculus/Mus_musculus.GRCm39.109.gtf.gz -O "$(basename "$GTF")"
fi

# 2) 构建 STAR 索引
mkdir -p "$STAR_INDEX"
if [[ ! -s "$STAR_INDEX/SAindex" ]]; then
  echo "[IDX] Building STAR index..."
  # sjdbOverhang 建议 = 读长-1；若读长不确定，用100通用
  STAR --runThreadN "$THREADS" \
       --runMode genomeGenerate \
       --genomeDir "$STAR_INDEX" \
       --genomeFastaFiles "$GENOME_FA" \
       --sjdbGTFfile "$GTF" \
       --sjdbOverhang 100 \
       --readFilesCommand zcat
fi

echo "[OK] Reference ready at: $REF_DIR"
