#!/usr/bin/env bash
set -euo pipefail

# 1) 载入配置并回显关键路径
source scripts/00_config.env
echo "[CONFIG] PROJECT_DIR=$PROJECT_DIR"
echo "[CONFIG] RAW=$RAW"
echo "[CONFIG] INTERIM=$INTERIM"
echo "[CONFIG] QC_OUT=$QC_OUT"
echo "[CONFIG] REF_DIR=$REF_DIR"
echo "[CONFIG] STAR_INDEX=$STAR_INDEX"
echo

# 设置要检查的批次（按需改）
BATCHES=("batch1" "batch2" "batch3")

# 一些小工具函数
ok(){ echo "  [OK] $*"; }
warn(){ echo "  [WARN] $*"; }
fail(){ echo "  [FAIL] $*"; }

# 2) 参考文件与索引完整性（到 03 步）
echo "== REF / INDEX 检查 =="
# 压缩包是否存在
test -r "$GENOME_FA"   && ok "找到参考基因组.gz: $(basename "$GENOME_FA")"   || fail "缺少参考基因组.gz: $GENOME_FA"
test -r "$GTF"         && ok "找到注释.gtf.gz:    $(basename "$GTF")"         || fail "缺少注释.gtf.gz: $GTF"
# 是否存在解压后的纯文本（若没有，仅提示；03 索引正确构建则可不强制）
if [[ -r "$REF_DIR/Mus_musculus.GRCm39.dna.primary_assembly.fa" ]]; then ok "存在解压FA"; else warn "未检测到解压FA（可选）"; fi
if [[ -r "$REF_DIR/Mus_musculus.GRCm39.109.gtf" ]]; then ok "存在解压GTF"; else warn "未检测到解压GTF（可选）"; fi
# STAR 索引关键文件
for f in Genome SA SAindex chrLength.txt geneInfo.tab sjdbInfo.txt sjdbList.fromGTF.out.tab; do
  if [[ -r "$STAR_INDEX/$f" ]]; then ok "索引含 $f"; else warn "索引缺少 $f（若缺 geneInfo.tab/sjdb*，说明索引没带注释）"; fi
done
echo

# 3) 原始 FASTQ 自检（00→01 前置）
echo "== RAW FASTQ 配对与完整性 =="
for B in "${BATCHES[@]}"; do
  RAWB="$RAW/$B"
  if [[ ! -d "$RAWB" ]]; then warn "$B: 原始目录不存在，跳过"; continue; fi
  echo "-- $B --"
  # 清点 Apple 的 '._' 隐藏文件
  HIDDEN=$(find "$RAWB" -maxdepth 1 -type f -name '._*' | wc -l | tr -d ' ')
  if [[ "$HIDDEN" != "0" ]]; then warn "$B: 发现 $HIDDEN 个 '._*' 隐藏文件（建议删除或忽略）"; else ok "$B: 无 '._*' 隐藏文件"; fi

  # 配对检查（R1 vs R2）
  R1S=$(find "$RAWB" -maxdepth 1 -type f -E -regex ".*(_R1|_1)\.f(ast)?q\.gz" | sed -E 's_.*/__' | sort) || true
  R2S=$(find "$RAWB" -maxdepth 1 -type f -E -regex ".*(_R2|_2)\.f(ast)?q\.gz" | sed -E 's_.*/__' | sort) || true
  # 统一基名
  R1BASE=$(printf "%s\n" "$R1S" | sed -E 's/_R1\.f(ast)?q\.gz$//' | sed -E 's/_1\.f(ast)?q\.gz$//' | sort || true)
  R2BASE=$(printf "%s\n" "$R2S" | sed -E 's/_R2\.f(ast)?q\.gz$//' | sed -E 's/_2\.f(ast)?q\.gz$//' | sort || true)

  # 找孤儿
  ORPH_R1=$(comm -23 <(printf "%s\n" "$R1BASE") <(printf "%s\n" "$R2BASE") || true)
  ORPH_R2=$(comm -13 <(printf "%s\n" "$R1BASE") <(printf "%s\n" "$R2BASE") || true)
  if [[ -n "$ORPH_R1" || -n "$ORPH_R2" ]]; then
    fail "$B: R1/R2 不成对："
    [[ -n "$ORPH_R1" ]] && echo "   · R2 缺失: $(echo "$ORPH_R1" | tr '\n' ' ')"
    [[ -n "$ORPH_R2" ]] && echo "   · R1 缺失: $(echo "$ORPH_R2" | tr '\n' ' ')"
  else
    ok "$B: R1/R2 成对数 = $(printf "%s\n" "$R1BASE" | grep -c . || true)"
  fi

  # 快速 gzip 完整性抽检（每对抽 1 个）
  ONE_R1=$(printf "%s\n" "$R1S" | head -n 1 || true)
  if [[ -n "$ONE_R1" ]]; then
    if gzip -t "$RAWB/$ONE_R1" 2>/dev/null; then ok "$B: gzip 抽检通过 ($ONE_R1)"; else fail "$B: gzip 抽检失败 ($ONE_R1)"; fi
  fi
done
echo

# 4) FastQC 原始质控输出（01）
echo "== FastQC 原始质控输出 =="
for B in "${BATCHES[@]}"; do
  D="$QC_OUT/$B/fastqc_raw"
  if [[ -d "$D" ]]; then
    NHTML=$(ls -1 "$D"/*_fastqc.html 2>/dev/null | wc -l | tr -d ' ')
    NZIP=$(ls -1 "$D"/*_fastqc.zip  2>/dev/null | wc -l | tr -d ' ')
    echo "-- $B --"
    ok "fastqc_raw: html=$NHTML  zip=$NZIP"
  else
    warn "$B: 未找到 $D（可能 01 未跑或输出目录不同）"
  fi
done
echo

# 5) fastp 剪切产物（02）
echo "== fastp 剪切产物 =="
for B in "${BATCHES[@]}"; do
  T="$INTERIM/$B/trimmed"
  R="$QC_OUT/$B/fastp"
  if [[ ! -d "$T" ]]; then warn "$B: 未找到 $T"; continue; fi
  echo "-- $B --"
  N1=$(ls -1 "$T"/*_R1.trim.fastq.gz 2>/dev/null | wc -l | tr -d ' ')
  N2=$(ls -1 "$T"/*_R2.trim.fastq.gz 2>/dev/null | wc -l | tr -d ' ')
  ok "trimmed 配对计数: R1=$N1 R2=$N2（应相等）"
  # 抽检 gzip 完整性
  ONE_T=$(ls -1 "$T"/*_R1.trim.fastq.gz 2>/dev/null | head -n1 || true)
  if [[ -n "$ONE_T" ]]; then
    if gzip -t "$ONE_T" 2>/dev/null; then ok "trimmed 抽检通过 ($(basename "$ONE_T"))"; else fail "trimmed 抽检失败 ($(basename "$ONE_T"))"; fi
  fi
  # fastp 报告
  if [[ -d "$R" ]]; then
    NJ=$(ls -1 "$R"/*.json 2>/dev/null | wc -l | tr -d ' ')
    NH=$(ls -1 "$R"/*.html 2>/dev/null | wc -l | tr -d ' ')
    ok "fastp 报告: json=$NJ html=$NH"
  else
    warn "$B: 未找到 fastp 报告目录: $R"
  fi
done
echo

# 6) （可选）fastp JSON 快速汇总（需要 jq；没有就跳过）
if command -v jq >/dev/null 2>&1; then
  echo "== fastp JSON 快速汇总（通过/过滤/dup 率） =="
  for B in "${BATCHES[@]}"; do
    R="$QC_OUT/$B/fastp"
    [[ -d "$R" ]] || continue
    for J in "$R"/*.json; do
      [[ -r "$J" ]] || continue
      S=$(basename "$J" .fastp.json)
      TOTAL=$(jq -r '.summary.before_filtering.total_reads' "$J")
      PASS=$(jq -r '.summary.after_filtering.total_reads' "$J")
      DUP=$(jq -r '.duplication.rate // .duplication' "$J" 2>/dev/null || echo "NA")
      echo "  $B/$S: total=$TOTAL  pass=$PASS  dup=$DUP"
    done
  done
else
  warn "未安装 jq，跳过 fastp JSON 汇总（可用 micromamba 安装：micromamba install -y -n rnaseq_env jq）"
fi
echo

echo "== 自检完成：若以上均为 [OK]/[WARN] 且无关键 [FAIL]，说明 00→03 阶段产物基本正常，可着手修复 04。=="
