import os
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

# ---------- 如果目录不存在，则自动创建 ----------
save_path = "results/figs/Figure_6_a.svg"
os.makedirs(os.path.dirname(save_path), exist_ok=True)

# ---------- 一些小工具函数 ----------
def add_box(ax, xy, width, height, text, fontsize=9, ha="center"):
    x, y = xy
    box = FancyBboxPatch(
        (x, y),
        width,
        height,
        boxstyle="round,pad=0.2,rounding_size=0.05",
        linewidth=1,
        edgecolor="black",
        facecolor="white"
    )
    ax.add_patch(box)
    ax.text(
        x + width/2,
        y + height/2,
        text,
        ha=ha,
        va="center",
        fontsize=fontsize
    )
    return box

def add_arrow(ax, start, end, text=None, fontsize=8, offset=0.03):
    arrow = FancyArrowPatch(start, end, arrowstyle="->", linewidth=1)
    ax.add_patch(arrow)
    if text:
        mx = (start[0] + end[0]) / 2
        my = (start[1] + end[1]) / 2
        ax.text(mx, my + offset, text, ha="center", va="bottom", fontsize=fontsize)
    return arrow

# ---------- 生成图形 ----------
fig, ax = plt.subplots(figsize=(8, 6))
ax.set_axis_off()
ax.set_xlim(0, 10)
ax.set_ylim(0, 7)

# 线粒体轮廓
mito = FancyBboxPatch(
    (1.0, 1.0),
    8.0,
    5.0,
    boxstyle="round,pad=0.4,rounding_size=0.8",
    linewidth=1.5,
    edgecolor="black",
    facecolor="none"
)
ax.add_patch(mito)
ax.text(5.0, 5.8, "Mitochondrion", ha="center", va="bottom", fontsize=12, fontweight="bold")

# --------- 各条机制轴 ---------
syn_box = add_box(ax, (1.3, 4.4), 3.3, 1.2,
    "Synthesis axis\nPA → CDP-DAG → PG → CL\n"
    "Key genes: CDS, PGS1, CRLS1, TAZ/LYCAT\n"
    "Lipid readouts: CL/PG, ΣCL", fontsize=8)

rem_box = add_box(ax, (5.4, 4.4), 3.3, 1.2,
    "Remodeling axis\nMLCL ↔ CL\n"
    "Key genes: TAZ/LYCAT, PLA2G6\n"
    "Lipid readouts: MLCL/CL", fontsize=8)

ox_box = add_box(ax, (1.3, 1.4), 3.3, 1.3,
    "Oxidation axis\nCL → oxCL (ROS)\n"
    "Key genes: PLA2 family, ALOX\n"
    "Lipid readouts: oxCL/CL", fontsize=8)

ts_box = add_box(ax, (5.4, 1.4), 3.3, 1.3,
    "Transport/Supply axis\nER→mito, PG/PA supply\n"
    "Key genes: PRELID, TRIAP1\n"
    "Lipid readouts: PG/PA, ΣPG", fontsize=8)

# 箭头示意
add_arrow(ax, (3.0, 4.4), (3.0, 3.7), "CL synthesis")
add_arrow(ax, (7.0, 4.4), (7.0, 3.7), "CL remodeling")

ax.text(3.0, 3.2, "CL", fontsize=9)
ax.text(3.8, 3.2, "→", fontsize=9)
ax.text(4.6, 3.2, "oxCL", fontsize=9)
add_arrow(ax, (3.0, 3.0), (4.6, 3.0), "oxidation")

er_box = add_box(ax, (0.2, 2.7), 1.0, 0.8, "ER\n(PA, PG)", fontsize=8)
add_arrow(ax, (1.2, 3.1), (1.9, 3.1), "PG/PA transport")

# 说明文字
ax.text(
    5.0,
    0.4,
    "Each axis is quantified by both transcriptional and lipidomic readouts,\n"
    "yielding a directional Z score (mechanistic effect).",
    ha="center",
    fontsize=9
)

# ---------- 保存为 SVG ----------
fig.savefig(save_path, dpi=300, bbox_inches="tight")
plt.close(fig)

print(f"图像已成功保存到: {save_path}")