#!/usr/bin/env python3
import os, sys, yaml
import numpy as np
import pandas as pd
from sklearn.impute import KNNImputer

def pqn_norm(df):
    # Probabilistic Quotient Normalization
    ref = df.median(axis=1)  # feature-wise median -> reference spectrum
    quot = df.div(ref, axis=0)
    factors = quot.median(axis=0)  # sample-wise factor
    return df.div(factors, axis=1)

def main():
    if len(sys.argv) != 2:
        print("Usage: 03_impute_normalize_scale.py <config.yaml>", file=sys.stderr)
        sys.exit(1)
    cfg = yaml.safe_load(open(sys.argv[1]))
    P = cfg["project_dir"]
    out_dir = os.path.join(P, cfg["out_dir"])
    os.makedirs(out_dir, exist_ok=True)

    wide_path = os.path.join(out_dir, "02_wide_blank_confident.tsv")
    if not os.path.isfile(wide_path):
        print(f"[FATAL] not found: {wide_path}", file=sys.stderr); sys.exit(1)
    X = pd.read_csv(wide_path, sep="\t", index_col=0)

    # 0 视为缺失
    X = X.replace(0, np.nan)

    # MNAR: half-min 按特征
    if cfg["input"]["impute"].get("mnar_method","halfmin") == "halfmin":
        mins = X.min(axis=1, skipna=True)
        halfmins = mins / 2.0
        X = X.apply(lambda s: s.fillna(halfmins[s.name]))

    # MAR: KNN
    if cfg["input"]["impute"].get("mar_method","knn") == "knn":
        imp = KNNImputer(n_neighbors=5)
        X[:] = imp.fit_transform(X)

    # Normalize
    norm = cfg["input"].get("normalize","TIC").upper()
    if norm == "TIC":
        totals = X.sum(axis=0)
        X = X.div(totals, axis=1) * totals.median()
    elif norm == "PQN":
        X = pqn_norm(X)
    # IS 可在此扩展

    # Transform
    if cfg["input"].get("transform","log2").lower() == "log2":
        X = np.log2(X + 1e-9)

    # Scale (Pareto)
    if cfg["input"].get("scale","pareto").lower() == "pareto":
        mu = X.mean(axis=1)
        sd = X.std(axis=1, ddof=1).replace(0, np.nan)
        X = X.sub(mu, axis=0).div(np.sqrt(sd), axis=0)

    out_path = os.path.join(out_dir, "03_imputed_normalized_scaled.tsv")
    X.to_csv(out_path, sep="\t")
    print(f"[OK] impute/normalize/scale -> {out_path}")

if __name__ == "__main__":
    main()