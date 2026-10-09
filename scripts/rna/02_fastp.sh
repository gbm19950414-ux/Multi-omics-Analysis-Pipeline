#!/usr/bin/env bash
set -euo pipefail

# 读取路径配置
HERE="$(cd "$(dirname "$0")" && pwd)"
source "$HERE/00_config.env"

THREADS="${THREADS:-8}"

# BATCH 可选；不设则全批次
if [[ -n "${BATCH:-}" ]]; then
  batches=("$BATCH")
else
  batches=(batch1 batch2 batch3)
fi

echo "== fastp 开始 =="
echo "RAW=$RAW"
echo "INTERIM=$INTERIM"
echo "QC_OUT=$QC_OUT"
echo "BATCHES=${batches[*]}"
echo

shopt -s nullglob

for B in "${batches[@]}"; do
  IN="${RAW%/}/$B"
  TRIM="${INTERIM%/}/$B/trimmed"
  OUT="${QC_OUT%/}/$B/fastp"
  mkdir -p "$TRIM" "$OUT"

  echo "== 处理 $B =="
  # 兼容 *_R1*.fastq.gz / *_1*.fq.gz 等
  r1=("$IN"/*_R1*.f*q.gz "$IN"/*_1*.f*q.gz)

  for R1 in "${r1[@]}"; do
    base="${R1##*/}"
    # 找配对的 R2
    R2="${base/_R1/_R2}"; R2="${R2/_1/_2}"
    if [[ ! -e "$IN/$R2" ]]; then
      echo "  [SKIP] 缺R2: $base"
      continue
    fi
    R2="$IN/$R2"

    # 样本名：只剥离结尾的读段标记
    stem="${base%.fastq.gz}"; stem="${stem%.fq.gz}"
    S="$(printf '%s\n' "$stem" | sed -E 's/(_R?[12])$//')"

    O1="$TRIM/${S}_R1.trim.fastq.gz"
    O2="$TRIM/${S}_R2.trim.fastq.gz"
    J="$OUT/${S}.fastp.json"
    H="$OUT/${S}.fastp.html"

    need=0
    [[ -s "$O1" ]] && gzip -t "$O1" &>/dev/null || need=1
    [[ -s "$O2" ]] && gzip -t "$O2" &>/dev/null || need=1
    [[ -s "$J"  ]] && grep -q "\"summary\"" "$J" 2>/dev/null || need=1
    [[ -s "$H"  ]] || need=1

    if (( need==0 )); then
      echo "  [OK ] $S"
      continue
    fi

    echo "  [RUN] $S"
    rm -f "$O1" "$O2" "$J" "$H"
    fastp --in1 "$R1" --in2 "$R2" \
          --out1 "$O1" --out2 "$O2" \
          --detect_adapter_for_pe --thread "$THREADS" \
          --html "$H" --json "$J"
  done
done

echo
echo "== 构建 fastp 汇总（JSON → TSV） =="
python - <<'PY'
import os,glob,json
qc=os.environ.get("QC_OUT","")
rows=[]
for j in sorted(glob.glob(os.path.join(qc,"batch*","fastp","*.fastp.json"))):
    batch=os.path.basename(os.path.dirname(os.path.dirname(j)))
    sample=os.path.basename(j).replace(".fastp.json","")
    try:
        d=json.load(open(j))
    except Exception:
        continue
    bef=d.get("summary",{}).get("before_filtering",{})
    aft=d.get("summary",{}).get("after_filtering",{})
    rb=bef.get("total_reads",0)
    ra=aft.get("total_reads",0)
    q30b=bef.get("q30_rate",None)
    q30a=aft.get("q30_rate",None)
    dup=d.get("duplication",{}).get("rate",None)
    trim=d.get("adapter_cutting",{}).get("total_trimmed_bases",None)
    passrate=(ra/rb) if rb else None
    rows.append([batch,sample,rb,ra,
                 None if q30b is None else round(q30b,3),
                 None if q30a is None else round(q30a,3),
                 None if passrate is None else round(passrate,3),
                 None if dup is None else round(dup,3),
                 "" if trim is None else trim])
out=os.path.join(qc,"fastp_summary.tsv")
with open(out,"w") as f:
    f.write("\t".join(["batch","sample","reads_before","reads_after","q30_before","q30_after","pass_rate","dup_rate","adapter_trimmed_bases"])+"\n")
    for r in rows:
        f.write("\t".join("" if v is None else str(v) for v in r)+"\n")
print("[OK] 写出:", out, "(样本数:", len(rows), ")")
PY

echo "== 完成 =="