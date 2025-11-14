#!/usr/bin/env bash
# Purpose: Validate STAR reference index integrity & freshness (and optionally rebuild)
# Usage:
#   bash scripts/rna/03_star_index.sh           # 只做完整性/质量检测
#   REBUILD=1 bash scripts/rna/03_star_index.sh # 若检测失败或路径不一致则触发重建（调用 00_fetch_reference.sh）

set -euo pipefail

# --- load config ---
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
. "$SCRIPT_DIR/00_config.env"

REF_DIR=${REF_DIR%/}
GENOME_FA=${GENOME_FA%/}
GTF=${GTF%/}
STAR_INDEX=${STAR_INDEX%/}
UNZ_DIR="$REF_DIR/unzipped"
FA_UNZ="$UNZ_DIR/$(basename "${GENOME_FA%.gz}")"
GTF_UNZ="$UNZ_DIR/$(basename "${GTF%.gz}")"

PASS=1 # overall pass flag (1=pass until proven otherwise)

say() { printf '%s\n' "$*"; }
ok()  { printf ' [OK ] %s\n' "$*"; }
warn(){ printf ' [WARN] %s\n' "$*"; PASS=1; }
fail(){ printf ' [FAIL] %s\n' "$*"; PASS=0; }
size_of(){ du -h "$1" 2>/dev/null | awk '{print $1}'; }
exists_nonempty(){ [ -s "$1" ]; }
contains(){ grep -Fq -- "$2" "$1"; }

say "== 路径基线 =="
say "REF_DIR=$REF_DIR"
say "GENOME_FA=$GENOME_FA"
say "GTF=$GTF"
say "STAR_INDEX=$STAR_INDEX"

say "\n== 原始参考压缩件是否存在且非空 =="
exists_nonempty "$GENOME_FA" && ok "$GENOME_FA" || fail "$GENOME_FA 不存在或大小为0"
exists_nonempty "$GTF"       && ok "$GTF"       || fail "$GTF 不存在或大小为0"

say "\n== 解压件是否存在且非空 =="
exists_nonempty "$FA_UNZ" && ok "$FA_UNZ" || fail "$FA_UNZ 缺失，请执行 00_fetch_reference.sh"
exists_nonempty "$GTF_UNZ" && ok "$GTF_UNZ" || fail "$GTF_UNZ 缺失，请执行 00_fetch_reference.sh"

# 轻量质量检查
if exists_nonempty "$FA_UNZ"; then
  # 统计 fasta 条目数
  FA_N=$(grep -c '^>' "$FA_UNZ" || true)
  [ "${FA_N:-0}" -gt 0 ] && ok "FASTA 条目数: $FA_N" || fail "FASTA 条目数=0，文件可能损坏"
fi
if exists_nonempty "$GTF_UNZ"; then
  # 行数（注释行+记录行）
  GTF_N=$(wc -l < "$GTF_UNZ" 2>/dev/null || echo 0)
  [ "${GTF_N:-0}" -gt 100000 ] && ok "GTF 行数: $GTF_N" || warn "GTF 行数异常(=$GTF_N)，请复核是否为完整注释"
fi

say "\n== STAR 索引关键文件是否存在 =="
need=( SAindex Genome SA chrName chrNameLength genomeParameters.txt )
missing=0
for f in "${need[@]}"; do
  p="$STAR_INDEX/$f"
  if exists_nonempty "$p"; then
    ok "$p ($(size_of "$p"))"
  else
    fail "$p 缺失"
    missing=1
  fi
done

# 路径一致性检查（防止从旧磁盘迁移后仍指向 MyPassport）
say "\n== genomeParameters.txt 路径一致性检查 =="
GPAR="$STAR_INDEX/genomeParameters.txt"
if exists_nonempty "$GPAR"; then
  if grep -q "/Volumes/MyPassport/" "$GPAR"; then
    fail "genomeParameters.txt 仍指向 /Volumes/MyPassport/，建议重建索引"
  else
    ok "genomeParameters.txt 未发现旧路径残留"
  fi
  # 校验记录的 fasta/gtf 路径是否与现用一致（仅提示，不作为致命失败）
  rec_fa=$(awk '/genomeFastaFiles/{print $2}' "$GPAR" 2>/dev/null || true)
  rec_gtf=$(awk '/sjdbGTFfile/{print $2}' "$GPAR" 2>/dev/null || true)
  [ -n "$rec_fa" ] && [ "$rec_fa" != "${FA_UNZ}" ] && warn "记录的 genomeFastaFiles=[$rec_fa] 与当前=[$FA_UNZ] 不一致"
  [ -n "$rec_gtf" ] && [ "$rec_gtf" != "${GTF_UNZ}" ] && warn "记录的 sjdbGTFfile=[$rec_gtf] 与当前=[$GTF_UNZ] 不一致"
fi

# 新鲜度检查：索引是否比解压参考文件更新
say "\n== 索引新鲜度检查 =="
probe="$STAR_INDEX/SAindex"
if exists_nonempty "$probe" && exists_nonempty "$FA_UNZ" && exists_nonempty "$GTF_UNZ"; then
  if [ "$probe" -nt "$FA_UNZ" ] && [ "$probe" -nt "$GTF_UNZ" ]; then
    ok "索引文件更新时间新于参考序列/注释"
  else
    warn "索引可能过旧（比 FASTA/GTF 解压件旧）"
  fi
fi

# 汇总与可选重建
if [ "$PASS" -eq 1 ] && [ "$missing" -eq 0 ]; then
  say "\n[SUMMARY] STAR 索引完整性检查：PASS"
else
  say "\n[SUMMARY] STAR 索引完整性检查：FAIL/需处理"
  if [ "${REBUILD:-0}" = "1" ]; then
    say "触发重建：调用 00_fetch_reference.sh"
    bash "$SCRIPT_DIR/00_fetch_reference.sh"
  else
    say "提示：可使用 REBUILD=1 重新构建，或手动运行 $SCRIPT_DIR/00_fetch_reference.sh"
  fi
fi
