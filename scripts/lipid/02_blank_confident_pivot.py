#!/usr/bin/env python3
import os, sys, re, yaml
import pandas as pd
import numpy as np

# --- compatibility shim for old bash lacking mapfile ---
def readarray_compat(name, lines):
    globals()[name] = [line.strip() for line in lines if line.strip()]

def main():
    if len(sys.argv) != 2:
        print("Usage: 02_blank_confident_pivot.py <config.yaml>", file=sys.stderr)
        sys.exit(1)

    cfg = yaml.safe_load(open(sys.argv[1]))
    P = cfg["project_dir"]
    out_dir = os.path.join(P, cfg["out_dir"])
    os.makedirs(out_dir, exist_ok=True)

    long_path = os.path.join(out_dir, "01_combined_long.csv")
    if not os.path.isfile(long_path):
        print(f"[FATAL] not found: {long_path}", file=sys.stderr)
        sys.exit(1)

    long = pd.read_csv(long_path)

    # 确保字符串列为 str（修复 “float64 + '|' + str” 报错）
    for col in ["SampleID","Group","Batch","LipidName","Class"]:
        if col in long.columns:
            long[col] = long[col].astype(str)

    # 置信度阈值
    score_min = cfg["input"].get("score_min", None)
    if score_min is not None and "Score" in long.columns:
        long = long[ long["Score"].isna() | (long["Score"] >= float(score_min)) ].copy()

    # 空白过滤：mask 在 merge 前固定为一列，避免 reindex 警告
    blank_pat = cfg["input"]["patterns"].get("blank_regex", "")
    if blank_pat:
        is_blank = long["SampleID"].str.match(re.compile(blank_pat))
        long["__is_blank__"] = is_blank

        blank_bg = (long[long["__is_blank__"]]
                    .groupby(["Class","LipidName"])["Intensity"]
                    .median().rename("BlankMed"))
        long = long.merge(blank_bg, on=["Class","LipidName"], how="left")
        long["BlankMed"] = long["BlankMed"].fillna(0)

        real_mean = (long[~long["__is_blank__"]]
                     .groupby(["Class","LipidName"])["Intensity"]
                     .mean().rename("RealMean"))
        long = long.merge(real_mean, on=["Class","LipidName"], how="left")
        long["RealMean"] = long["RealMean"].replace(0, np.nan)

        ratio_max = float(cfg["input"].get("blank_ratio_max", 0.2))
        keep_feat = (long["BlankMed"] / long["RealMean"]) < ratio_max
        long = long[ keep_feat | long["__is_blank__"] ].copy()
        # 真正输出前，去掉空白行
        long = long[ ~long["__is_blank__"] ].copy()
        long.drop(columns=["__is_blank__","BlankMed","RealMean"], errors="ignore", inplace=True)

    # 检出率过滤（按样本维度）
    detect_rate_min = float(cfg["input"].get("detect_rate_min", 0.0))
    if detect_rate_min > 0:
        # 先透视到宽表，按“>0”计检出
        wide_tmp = long.pivot_table(index=["Class","LipidName"],
                                    columns="SampleID", values="Intensity", aggfunc="mean")
        det = (wide_tmp.fillna(0) > 0).mean(axis=1)
        keep_idx = det[det >= detect_rate_min].index
        long = long.set_index(["Class","LipidName"]).loc[keep_idx].reset_index()

    # 透视为宽表：Feature = Class|LipidName
    long["Feature"] = long["Class"] + "|" + long["LipidName"]
    wide = long.pivot_table(index="Feature", columns="SampleID", values="Intensity", aggfunc="mean").sort_index()

    wide_path = os.path.join(out_dir, "02_wide_blank_confident.tsv")
    annot_path = os.path.join(out_dir, "02_feature_annotation.tsv")
    wide.to_csv(wide_path, sep="\t")
    (long[["Class","LipidName","Feature"]].drop_duplicates()
         .to_csv(annot_path, sep="\t", index=False))

    print(f"[OK] wide saved -> {wide_path}")
    print(f"[OK] feature annotation -> {annot_path}")

if __name__ == "__main__":
    main()