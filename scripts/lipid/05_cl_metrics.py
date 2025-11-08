#!/usr/bin/env python3
import os, sys, yaml, pandas as pd

def main():
    if len(sys.argv) != 2:
        print("Usage: 05_cl_metrics.py <config.yaml>", file=sys.stderr)
        sys.exit(1)
    cfg = yaml.safe_load(open(sys.argv[1]))
    P = cfg["project_dir"]; out_dir = os.path.join(P, cfg["out_dir"])

    # 取 ComBat 后矩阵；若未跑04，可改读 03_
    mat = os.path.join(out_dir, "04_batch_corrected.tsv")
    if not os.path.isfile(mat):
        mat = os.path.join(out_dir, "03_imputed_normalized_scaled.tsv")
    X = pd.read_csv(mat, sep="\t", index_col=0)

    annot = pd.read_csv(os.path.join(out_dir, "02_feature_annotation.tsv"), sep="\t")
    annot.index = annot["Feature"]

    # 计算各类脂总强度、CL占比
    cls = annot["Class"].unique()
    sums = {}
    for cl in cls:
        feats = annot[annot["Class"]==cl].index
        inter = [f for f in feats if f in X.index]
        if inter:
            sums[cl] = X.loc[inter].sum(axis=0)
    sums_df = pd.DataFrame(sums)
    sums_df.to_csv(os.path.join(out_dir, "05_class_sums.tsv"), sep="\t")

    if "CL" in sums_df.columns:
        total = sums_df.sum(axis=1)
        cl_pct = (sums_df["CL"] / total).fillna(0)
        cl_pct.to_csv(os.path.join(out_dir, "05_CL_fraction.tsv"), sep="\t", header=["CL_fraction"])

    print("[OK] CL metrics exported.")

if __name__ == "__main__":
    main()