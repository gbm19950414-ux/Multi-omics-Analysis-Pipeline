#!/usr/bin/env python3
"""
plot_figure_7a_model.py

生成 Figure 7A: "EphB1-dependent signaling pathways mediating cardiolipin bottleneck"
的概念模型图。

特点：
- 使用 matplotlib 绘制方框 + 箭头
- 保存为 SVG，所有文字保持为可编辑文本（svg.fonttype='none'）
- 通过 CONFIG 在顶部修改各个方框的文字内容
"""

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

# ===================== 1. 配置区域：你可以在这里改文字 =====================

CONFIG = {
    "font_family": "Arial",   # 如果你机器上没有 Arial，可以改成 "Helvetica" 或其他
    "output_path": "results/figs/figure7a_ephb1_signaling_model.svg",

    # 左侧：X（基因状态）
    "left_box": "EphB1 loss\n(HO vs WT)",

    # 中间：1–2 条关键信号通路
    "middle_boxes": [
        "PGC-1α /\nmitochondrial\nbiogenesis",
        "NF-κB /\ninflammatory ROS"
    ],

    # 右侧：对应机制轴（可以和 Figure 6 的命名保持一致）
    "right_boxes": [
        "Synthesis axis ↓\n(CL_PG_ratio ↓)",
        "Oxidation axis ↑\n(oxCL_CL_ratio ↑)"
    ],

    # 底部：CL phenotype / 终点表型
    "bottom_box": "CL availability ↓\n(tetra-18:2 CL ↓)\nIL-1β release ↓",

    # 一些视觉细节
    "box_facecolor": "white",
    "box_edgecolor": "black",
    "arrow_color": "black",
    "line_width": 1.2,
}


# ===================== 2. 一些小工具函数 =====================

def add_box(ax, xy, width, height, text, fontsize=11,
            facecolor="white", edgecolor="black", lw=1.2,
            boxstyle="round,pad=0.3,rounding_size=0.08"):
    """
    在给定坐标处画一个带圆角的方框，并在中间写字。

    Parameters
    ----------
    ax : matplotlib.axes.Axes
    xy : (x, y)
        方框中心坐标（数据坐标）
    width, height : float
        方框宽、高（数据坐标）
    text : str
        方框内文字（可换行）
    """
    x, y = xy
    # FancyBboxPatch 以左下角为锚点，这里做一个换算
    rect = FancyBboxPatch(
        (x - width / 2.0, y - height / 2.0),
        width, height,
        boxstyle=boxstyle,
        linewidth=lw,
        edgecolor=edgecolor,
        facecolor=facecolor
    )
    ax.add_patch(rect)

    ax.text(
        x, y, text,
        ha="center", va="center",
        fontsize=fontsize
    )

    return rect


from matplotlib.patches import Arrow

def add_arrow(ax, xy_start, xy_end, lw=1.2, color="black"):
    """
    使用最原始的 Arrow，不会产生 FancyArrowPatch 的控制线。
    """
    x0, y0 = xy_start
    x1, y1 = xy_end
    dx = x1 - x0
    dy = y1 - y0

    arr = Arrow(
        x0, y0, dx, dy,
        width=0.01,
        color=color
    )
    ax.add_patch(arr)
    return arr

# ===================== 3. 主函数：画 Figure 7A =====================

def draw_figure_7a(config: dict | None = None):
    if config is None:
        config = CONFIG

    # --- 全局字体与 SVG 文本设置 ---
    plt.rcParams["font.family"] = config["font_family"]
    # 关键：保持 SVG 中文本为 text 对象，而非路径
    plt.rcParams["svg.fonttype"] = "none"

    fig, ax = plt.subplots(figsize=(7.5, 4.5))  # 你可以视需要调整尺寸

    # 使用一个 0-1 的规范化坐标系，方便布局
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    # ------------- 布局坐标（可根据美观微调） -------------
    x_left = 0.12
    x_mid = 0.42
    x_right = 0.72
    x_bottom = 0.50

    # 左箱 & 中间两箱 & 右侧两箱 & 底部箱在 y 方向的布局
    y_center = 0.60
    y_top = 0.78
    y_bottom_mid = 0.42
    y_phenotype = 0.13

    box_w_left = 0.16
    box_h_left = 0.20

    box_w_mid = 0.18
    box_h_mid = 0.22

    box_w_right = 0.22
    box_h_right = 0.22

    box_w_bottom = 0.30
    box_h_bottom = 0.22

    fc = config["box_facecolor"]
    ec = config["box_edgecolor"]
    lw = config["line_width"]

    # ------------- 画左边 EphB1 方框 -------------
    left_box = add_box(
        ax,
        (x_left, y_center),
        box_w_left, box_h_left,
        config["left_box"],
        facecolor=fc,
        edgecolor=ec,
        lw=lw
    )

    # ------------- 画中间两条信号通路 -------------
    middle_boxes_text = config["middle_boxes"]
    # 如果你有 1 条通路，也可以只用第一个，脚本会忽略多余坐标
    mid_positions = [(x_mid, y_top), (x_mid, y_bottom_mid)]

    middle_boxes = []
    for i, text in enumerate(middle_boxes_text):
        if i >= len(mid_positions):
            break
        middle_boxes.append(
            add_box(
                ax,
                mid_positions[i],
                box_w_mid, box_h_mid,
                text,
                facecolor=fc,
                edgecolor=ec,
                lw=lw
            )
        )

    # ------------- 画右侧两条机制轴 -------------
    right_boxes_text = config["right_boxes"]
    right_positions = [(x_right, y_top), (x_right, y_bottom_mid)]

    right_boxes = []
    for i, text in enumerate(right_boxes_text):
        if i >= len(right_positions):
            break
        right_boxes.append(
            add_box(
                ax,
                right_positions[i],
                box_w_right, box_h_right,
                text,
                facecolor=fc,
                edgecolor=ec,
                lw=lw
            )
        )

    # ------------- 画底部 CL phenotype 方框 -------------
    bottom_box = add_box(
        ax,
        (x_bottom, y_phenotype),
        box_w_bottom, box_h_bottom,
        config["bottom_box"],
        facecolor=fc,
        edgecolor=ec,
        lw=lw
    )

    # ------------- 箭头：EphB1 → 中间通路 -------------
    arrow_color = config["arrow_color"]

    for mb in middle_boxes:
        # 左箱右侧中点
        x0 = x_left + box_w_left / 2.0
        y0 = y_center
        # 中间箱左侧中点
        bbox = mb.get_bbox().bounds  # (x, y, w, h)
        x1 = bbox[0]                 # 左边
        y1 = bbox[1] + bbox[3] / 2.0
        add_arrow(ax, (x0, y0), (x1, y1),
                  lw=lw, color=arrow_color)

    # ------------- 箭头：中间通路 → 机制轴 -------------
    for mb, rb in zip(middle_boxes, right_boxes):
        mb_box = mb.get_bbox().bounds
        rb_box = rb.get_bbox().bounds

        x0 = mb_box[0] + mb_box[2]      # 中间箱右侧
        y0 = mb_box[1] + mb_box[3] / 2.0
        x1 = rb_box[0]                  # 右侧箱左侧
        y1 = rb_box[1] + rb_box[3] / 2.0

        add_arrow(ax, (x0, y0), (x1, y1),
                  lw=lw, color=arrow_color)

    # ------------- 箭头：机制轴 → CL phenotype（从右侧两轴汇聚到下方） -------------
    # 简单做法：从两个机制轴的“中点”向下汇聚到 phenotype
    if right_boxes:
        x_mid_axes = sum(rb.get_bbox().bounds[0] + rb.get_bbox().bounds[2] / 2.0
                         for rb in right_boxes) / len(right_boxes)
        y_mid_axes = sum(rb.get_bbox().bounds[1]
                         for rb in right_boxes) / len(right_boxes)

        # 先从中点往下画一小段，再接到 bottom_box 顶部
        y_mid_down = (y_mid_axes + y_phenotype) / 2.0

        add_arrow(ax,
                  (x_mid_axes, y_mid_axes),
                  (x_mid_axes, y_mid_down),
                  lw=lw, color=arrow_color)

        bb = bottom_box.get_bbox().bounds
        x_bottom_center = bb[0] + bb[2] / 2.0
        y_bottom_top = bb[1] + bb[3]

        add_arrow(ax,
                  (x_mid_axes, y_mid_down),
                  (x_bottom_center, y_bottom_top),
                  lw=lw, color=arrow_color)

    # ------------- 标题（可选） -------------
    ax.text(
        0.5, 0.97,
        "EphB1-dependent signaling pathways mediating the cardiolipin bottleneck",
        ha="center", va="top",
        fontsize=12
    )

    # ------------- 保存图像 -------------
    output_path = config["output_path"]
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    print(f"Saved Figure 7A model to: {output_path}")


if __name__ == "__main__":
    draw_figure_7a()