#!/usr/bin/env bash
# === STAR 对齐（低内存/稳健版）：强制本地临时解压 + samtools 外排排序 ===
# 设计目标：
# 1) 规避 readFilesCommand 在某些环境下的“Failed spawning readFilesCommand”问题；
# 2) 控制内存峰值，避免 VSCode、系统被挤爆；
# 3) 幂等可续跑：已有产物则跳过；失败样本保留临时目录以便排查。

set -euo pipefail

# -------------------- 可调参数（环境变量覆盖） --------------------
THREADS="${THREADS:-2}"                 # STAR 线程（建议 ≤ 2）
READS_N="${READS_N:-0}"                # 仅调试时限制读取条数；0 表示不限制
STAR_SORT_RAM="${STAR_SORT_RAM:-0}"    # 若 >0，则启用 STAR 内建排序（不推荐在低内存机）
SAMTOOLS_SORT_MEM="${SAMTOOLS_SORT_MEM:-400M}" # samtools sort 单线程内存块
SAMTOOLS_THREADS="${SAMTOOLS_THREADS:-1}"      # samtools sort 线程（建议 1）
# 强制把临时文件放到 APFS SSD（避免系统 /var… 空间不足）
TMPDIR_SSD="/Volumes/Samsung_SSD_990_PRO_2TB_Media/multiomics_mech/tmp_star_sys"
export TMPDIR="${TMPDIR_OVERRIDE:-$TMPDIR_SSD}"

# -------------------- 路径基线 --------------------
STAR_INDEX="/Volumes/Samsung_SSD_990_PRO_2TB_Media/multiomics_mech/ref/star_index"
INTERIM="/Volumes/Samsung_SSD_990_PRO_2TB_Media/multiomics_mech/data/interim/rna"
WORK_TMP="/Volumes/Samsung_SSD_990_PRO_2TB_Media/multiomics_mech/tmp_star"   # 解压工作区

# 仅跑指定批次/样本：
#   例：BATCH=batch2  SAMPLE_GLOB="WT_*" bash scripts/rna/04_star_align.sh
BATCHES=(batch1 batch2 batch3)
if [[ -n "${BATCH:-}" ]]; then BATCHES=("$BATCH"); fi
SAMPLE_GLOB="${SAMPLE_GLOB:-*}"    # 匹配样本名的 glob，默认全部

# -------------------- 依赖检查 --------------------
if ! command -v STAR >/dev/null 2>&1; then
  echo "[FATAL] 未找到 STAR，请先在 rnaseq_env 中安装" >&2; exit 11; fi
if ! command -v samtools >/dev/null 2>&1; then
  echo "[FATAL] 未找到 samtools，请先安装（用于外排排序）" >&2; exit 12; fi

# -------------------- 目录与空间 --------------------
mkdir -p "$WORK_TMP" "$TMPDIR"
chmod 1777 "$WORK_TMP" || true
case "$TMPDIR" in
  /Volumes/Samsung_SSD_990_PRO_2TB_Media/multiomics_mech/tmp_star*) chmod 1777 "$TMPDIR" || true ;;
 esac

echo "== Space check =="
df -h "$WORK_TMP" /private/tmp /Volumes/Samsung_SSD_990_PRO_2TB_Media | sed -n '1,5p'
echo "[PARAM] THREADS=$THREADS READS_N=$READS_N TMPDIR=$TMPDIR SAMTOOLS_SORT_MEM=$SAMTOOLS_SORT_MEM SAMTOOLS_THREADS=$SAMTOOLS_THREADS"

# -------------------- 主循环：逐批次/样本 --------------------
for B in "${BATCHES[@]}"; do
  IN="${INTERIM%/}/$B/trimmed"
  OUTROOT="${INTERIM%/}/$B/star"
  mkdir -p "$OUTROOT"
  echo "== RUN $B =="

  shopt -s nullglob
  for R1 in "$IN"/${SAMPLE_GLOB}_R1.trim.fastq.gz; do
    fname="$(basename "$R1")"
    sample="${fname%_R1.trim.fastq.gz}"
    R2="$IN/${sample}_R2.trim.fastq.gz"
    OUT="$OUTROOT/$sample"
    mkdir -p "$OUT"

    # 幂等：BAM + 计数文件均存在则跳过
    if [[ -s "$OUT/Aligned.sortedByCoord.out.bam" && -s "$OUT/ReadsPerGene.out.tab" ]]; then
      echo "  [SKIP] $B/$sample 已完成"; continue; fi

    # ---------- 强制本地临时解压 ----------
    TMPD="$(mktemp -d "${WORK_TMP%/}/star_${B}_${sample}.XXXXXX")" || { echo "  [ERR] mktemp 失败"; exit 1; }
    echo "  [UNZIP] $B/$sample -> $TMPD"
    /usr/bin/gzip -dc "$R1" > "$TMPD/R1.fastq" || { echo "  [ERR] 解压 R1 失败"; rm -rf "$TMPD"; continue; }
    PE=0
    if [[ -f "$R2" ]]; then
      /usr/bin/gzip -dc "$R2" > "$TMPD/R2.fastq" || { echo "  [ERR] 解压 R2 失败"; rm -rf "$TMPD"; continue; }
      PE=1
    fi

    # ---------- 低内存策略：STAR 仅输出 Unsorted，再用 samtools 外排排序 ----------
    echo "  [STAR] $B/$sample (PE=$PE)"
    STAR \
      --runThreadN "$THREADS" \
      --genomeDir "$STAR_INDEX" \
      --readFilesIn "$TMPD/R1.fastq" ${PE:+"$TMPD/R2.fastq"} \
      --outFileNamePrefix "$OUT/" \
      --outSAMtype BAM Unsorted \
      --quantMode GeneCounts \
      ${READS_N:+--readMapNumber "$READS_N"} \
      ${STAR_SORT_RAM:+--limitBAMsortRAM "$STAR_SORT_RAM"}

    if [[ ! -s "$OUT/Aligned.out.bam" ]]; then
      echo "  [ERR] $B/$sample STAR 未生成 Aligned.out.bam，保留临时目录: $TMPD"; continue; fi

    : "${SAMTOOLS_THREADS:=1}"
    : "${SAMTOOLS_SORT_MEM:=400M}"
    echo "  [SORT] $B/$sample 由 samtools 外排排序（-m $SAMTOOLS_SORT_MEM -@ $SAMTOOLS_THREADS）"
    samtools sort -@ "$SAMTOOLS_THREADS" -m "$SAMTOOLS_SORT_MEM" \
      -T "$TMPDIR/sort_${B}_${sample}" \
      -o "$OUT/Aligned.sortedByCoord.out.bam" \
      "$OUT/Aligned.out.bam" && rm -f "$OUT/Aligned.out.bam"

    # 基本产物检查
    if [[ -s "$OUT/Aligned.sortedByCoord.out.bam" && -s "$OUT/ReadsPerGene.out.tab" ]]; then
      echo "  [OK ] $B/$sample 完成"
      rm -rf "$TMPD"
    else
      echo "  [ERR] $B/$sample 产物缺失（BAM/计数其一不存在），保留临时目录: $TMPD"
    fi
  done
  shopt -u nullglob
 done

echo "== DONE =="