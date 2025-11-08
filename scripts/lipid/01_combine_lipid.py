#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
01_combine_lipid.py

目的：在不可靠定列名的情况下，自动兼容你的 4 个原始批次表，
将其统一为 **长表** 格式，并导出样本元数据。

输入：单一参数 <config.yaml>
  - project_dir: 项目根目录
  - out_dir: 结果输出相对目录（例如 results/lipid 或 data/lipid/processed 均可）
  - raw_dir: 原始数据相对目录（例如 data/raw/lipid）
  - metadata_file: 样本元数据输出相对路径（例如 metadata/lipid_samples.tsv）
  - （可选）input.files: 原始文件相对路径列表；若缺省则自动扫描 raw_dir/*.csv

容错：
  - 自动识别常见列名同义词（大小写/空格/下划线/中英混排均做最大容忍）
  - 兼容两种表结构：
      A) 已经是“长表”：每行一个(样本, 脂质)的强度
      B) “宽表”：样本列为多列，统一熔融为长表
  - 批次 Batch：优先用文件内列，否则用文件名推断（例如 20240617 → Batch=20240617）
  - 分组 Group：若无显式列，尝试从 SampleID 前缀（如 KO_1 → KO）推断，失败则标为 Unknown

输出：
  - <out_dir>/01_combined_long.csv  （列：SampleID, Group, Batch, Class, LipidName, Intensity, Score[可空], SourceFile）
  - metadata/lipid_samples.tsv      （唯一的样本/分组/批次）
"""

import os
import re
import sys
import csv
import glob
import unicodedata
from typing import Dict, List, Tuple

import pandas as pd

# ------------------------------ 小工具 ------------------------------

def norm(s: str) -> str:
    """标准化列名：去空白、统一大小写，移除非字母数字字符。"""
    s = unicodedata.normalize("NFKC", str(s))
    s = s.strip().lower()
    s = re.sub(r"[\s\-]+", "_", s)
    s = re.sub(r"[^0-9a-zA-Z_]+", "", s)
    return s

def try_read_csv(path: str) -> pd.DataFrame:
    """以多编码尝试读取 CSV。"""
    for enc in ("utf-8-sig", "utf-8", "gb18030", "cp936"):
        try:
            return pd.read_csv(path, encoding=enc)
        except Exception:
            continue
    # 最后再试：严格 csv 解析（防止意外分隔符）
    return pd.read_csv(path, engine="python")

# 候选列同义词
SYN = {
    "sample": ["sample", "sampleid", "sample_id", "sample_name", "name", "filename", "file", "rawfile"],
    "group":  ["group", "condition", "phenotype", "class_group", "组", "分组"],
    "batch":  ["batch", "batchid", "run", "plate", "批次"],
    "lipid":  ["lipid", "lipidname", "compound", "feature", "metabolite", "name", "id", "品名", "化合物"],
    "class":  ["class", "lipidclass", "category", "superclass", "类别", "分类"],
    "intensity": [
        "intensity", "peak_area", "area", "abundance", "height", "value", "quant", "normalized", "norm_intensity"
    ],
    "score":  ["score", "confidence", "id_score", "msms_score", "ms2_score", "谱图得分", "置信度"],
}

# 样本列判定：
#   给定全部列（标准化后），哪些列看起来像“样本强度列”（用于宽表→长表）
SUSPECT_SAMPLE_PAT = re.compile(r"^(ko|wt|ctl|ctrl|case|ho|wt\d*|ko\d*|[a-z]+_\d+)$", re.I)


def pick_first(df_cols_norm: List[str], choices: List[str]) -> str:
    for c in df_cols_norm:
        if c in choices:
            return c
    return ""


def map_known_columns(df: pd.DataFrame) -> Dict[str, str]:
    """从 DataFrame 中识别关键列名；返回 {role: 原始列名}。若未找到返回 ''."""
    name_map = {norm(c): c for c in df.columns}
    inv = {v: k for k, v in name_map.items()}

    out = {}
    for role, alts in SYN.items():
        out[role] = ""
        for alt in alts:
            if alt in name_map:
                out[role] = name_map[alt]
                break

    # 若 sample/group/batch 未识别，稍后再做推断
    return out


def detect_wide_sample_columns(df: pd.DataFrame, known: Dict[str, str]) -> List[str]:
    """检测“宽表”中的样本强度列集合。排除已识别的注释列。"""
    ann_cols = {c for c in known.values() if c}
    cols_norm = {c: norm(c) for c in df.columns}
    sample_like = []
    for c in df.columns:
        if c in ann_cols:
            continue
        cn = cols_norm[c]
        # 含有明显样本模式/或首尾像样本名/或列名中包含 _R\d 等序列
        if SUSPECT_SAMPLE_PAT.match(cn) or re.search(r"_\d$", cn) or re.search(r"^(ko|wt)_", cn):
            sample_like.append(c)
        # 全数字/大多数字符也可能是样本，但风险较大，先不启用
    return sample_like


def infer_group_from_sample(sample: str) -> str:
    m = re.match(r"([A-Za-z]+)", str(sample))
    return m.group(1).upper() if m else "Unknown"


def infer_batch_from_filename(path: str) -> str:
    base = os.path.basename(path)
    m = re.search(r"(20\d{6})", base)  # 20240617 这样的日期
    if m:
        return m.group(1)
    m = re.search(r"(\d{8})", base)
    if m:
        return m.group(1)
    stem = os.path.splitext(base)[0]
    return stem


def to_long(df: pd.DataFrame, path: str) -> Tuple[pd.DataFrame, Dict[str, str]]:
    """将任意表格转为统一长表。
    返回 (long_df, colmap)；colmap 给出识别到的关键列对应。
    """
    colmap = map_known_columns(df)

    # 判断是否已是长表：是否同时具备 [样本, 脂质, 强度]
    is_long = bool(colmap.get("intensity")) and (
        bool(colmap.get("sample")) and bool(colmap.get("lipid"))
    )

    if is_long:
        long = df.copy()
        sample_col = colmap.get("sample") or "SampleID"
        lipid_col  = colmap.get("lipid")  or "LipidName"
        class_col  = colmap.get("class")  or None
        inten_col  = colmap.get("intensity") or None
        score_col  = colmap.get("score") or None
        group_col  = colmap.get("group") or None
        batch_col  = colmap.get("batch") or None

        rename_map = {}
        if sample_col: rename_map[sample_col] = "SampleID"
        if lipid_col:  rename_map[lipid_col]  = "LipidName"
        if class_col:  rename_map[class_col]  = "Class"
        if inten_col:  rename_map[inten_col]  = "Intensity"
        if score_col:  rename_map[score_col]  = "Score"
        if group_col:  rename_map[group_col]  = "Group"
        if batch_col:  rename_map[batch_col]  = "Batch"

        long = long.rename(columns=rename_map)

        # 若缺 Group/Batch，尝试推断
        if "Group" not in long.columns:
            long["Group"] = long["SampleID"].map(infer_group_from_sample)
        if "Batch" not in long.columns:
            long["Batch"] = infer_batch_from_filename(path)

        # 保留关键列
        keep = ["SampleID","Group","Batch","LipidName","Intensity"]
        if "Class" in long.columns: keep.insert(3, "Class")
        if "Score" in long.columns: keep.append("Score")
        long = long[keep]
        long["SourceFile"] = os.path.basename(path)
        return long, colmap

    # 否则视为“宽表”，需要找出样本强度列集合
    sample_cols = detect_wide_sample_columns(df, colmap)
    if not sample_cols:
        raise ValueError(f"无法在 {os.path.basename(path)} 中识别样本强度列，请检查列名或提供 config 中的列映射。")

    # 注释列优先级：lipid/class/score
    lipid_col = colmap.get("lipid") or "Lipid"
    if lipid_col not in df.columns:
        # 尝试若干候选
        for k in ["Lipid name", "Compound", "Feature", "Name", "ID"]:
            if k in df.columns:
                lipid_col = k
                break
    class_col = colmap.get("class") or ("Class" if "Class" in df.columns else None)
    score_col = colmap.get("score") if colmap.get("score") in df.columns else None

    ann_cols = [c for c in [lipid_col, class_col, score_col] if c]

    long = df.melt(
        id_vars=ann_cols,
        value_vars=sample_cols,
        var_name="SampleID",
        value_name="Intensity",
    )
    long = long.rename(columns={lipid_col: "LipidName"})
    if class_col:
        long = long.rename(columns={class_col: "Class"})
    if score_col:
        long = long.rename(columns={score_col: "Score"})

    long["Group"] = long["SampleID"].map(infer_group_from_sample)
    long["Batch"] = infer_batch_from_filename(path)
    long["SourceFile"] = os.path.basename(path)

    # 仅保留核心列（可选 Score）
    keep = ["SampleID","Group","Batch","LipidName","Intensity","SourceFile"]
    if "Class" in long.columns: keep.insert(3, "Class")
    if "Score" in long.columns: keep.append("Score")
    long = long[keep]
    return long, colmap


# ------------------------------ 主流程 ------------------------------

def main():
    if len(sys.argv) != 2:
        print("Usage: 01_combine_lipid.py <config.yaml>", file=sys.stderr)
        sys.exit(1)

    import yaml
    cfg = yaml.safe_load(open(sys.argv[1], "r"))

    P = cfg.get("project_dir", ".")
    raw_rel = cfg.get("raw_dir", "data/raw/lipid")
    out_rel = cfg.get("out_dir", "results/lipid")
    meta_rel = cfg.get("metadata_file", "metadata/lipid_samples.tsv")

    RAW_DIR = os.path.join(P, raw_rel)
    OUT_DIR = os.path.join(P, out_rel)
    META_PATH = os.path.join(P, meta_rel)
    os.makedirs(OUT_DIR, exist_ok=True)
    os.makedirs(os.path.dirname(META_PATH), exist_ok=True)

    # 文件清单：优先用 config.input.files；否则扫描 RAW_DIR/*.csv
    files = []
    in_cfg_files = (cfg.get("input", {}) or {}).get("files", [])
    if in_cfg_files:
        files = [os.path.join(P, f) for f in in_cfg_files]
    else:
        files = sorted(glob.glob(os.path.join(RAW_DIR, "*.csv")))
    if not files:
        print(f"[FATAL] 未找到原始 CSV：{RAW_DIR}", file=sys.stderr)
        sys.exit(1)

    all_long = []
    colmaps = []
    for fp in files:
        if not os.path.isfile(fp):
            print(f"[WARN] 找不到文件，跳过: {fp}")
            continue
        try:
            df = try_read_csv(fp)
        except Exception as e:
            print(f"[WARN] 读取失败，跳过 {os.path.basename(fp)}: {e}")
            continue
        long, cmap = to_long(df, fp)
        all_long.append(long)
        colmaps.append((os.path.basename(fp), {k: v for k, v in cmap.items() if v}))

    if not all_long:
        print("[FATAL] 没有可合并的数据表，请检查文件内容/列名。", file=sys.stderr)
        sys.exit(1)

    long = pd.concat(all_long, ignore_index=True)

    # 规范类型
    for col in ["SampleID","Group","Batch","LipidName"]:
        if col in long.columns:
            long[col] = long[col].astype(str)
    if "Class" in long.columns:
        long["Class"] = long["Class"].astype(str)
    # 强度转数值
    long["Intensity"] = pd.to_numeric(long["Intensity"], errors="coerce")

    out_long = os.path.join(OUT_DIR, "01_combined_long.csv")
    long.to_csv(out_long, index=False)
    print(f"[OK] combined long -> {out_long}")

    # 输出样本元数据
    meta = long.loc[:, [c for c in ["SampleID","Group","Batch"] if c in long.columns]].dropna().drop_duplicates()
    meta.to_csv(META_PATH, sep='\t', index=False)
    print(f"[OK] sample metadata -> {META_PATH}")

    # 附带打印每个源文件识别到的列映射，便于回溯
    print("[INFO] 列映射（供排错）：")
    for fname, mp in colmaps:
        print("  -", fname, ":", mp)

if __name__ == "__main__":
    main()