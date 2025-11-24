from graphviz import Digraph

# ===== Create graph =====
g = Digraph(
    name="panel_1c_pipeline",
    format="svg"
)

g.attr(
    rankdir="LR",
    fontsize="10",
    fontname="Helvetica"
)

# ===== RNA-seq branch =====
with g.subgraph(name="cluster_rna") as c:
    c.attr(
        label="RNA-seq (3 batches)",
        color="lightgrey",
        style="rounded",
        fontname="Helvetica"
    )
    c.node("r1", "QC")
    c.node("r2", "Alignment")
    c.node("r3", "Count matrix")
    c.node("r4", "Normalization + Batch correction\n(DESeq2: ~ batch + group)")
    c.node("r5", "Gene-level log2FC")
    c.node("r6", "Mechanism gene-set Z")
    c.edges([("r1","r2"),("r2","r3"),("r3","r4"),("r4","r5"),("r5","r6")])

# ===== Lipidomics branch =====
with g.subgraph(name="cluster_lipid") as c:
    c.attr(
        label="Lipidomics (Batch1–3, untargeted)",
        color="lightgrey",
        style="rounded",
        fontname="Helvetica"
    )
    c.node("l1", "QC & Normalization")
    c.node("l2", "Batch correction (ComBat)")
    c.node("l3", "ALR / CLR transform")
    c.node("l4", "Derive:\nCL/PG, MLCL/CL, oxCL/CL,\nΣCL, ΣPG")
    c.node("l5", "Mechanism lipid Z")
    c.edges([("l1","l2"),("l2","l3"),("l3","l4"),("l4","l5")])

# ===== Targeted CL panel =====
with g.subgraph(name="cluster_cl") as c:
    c.attr(
        label="CL targeted panel (Batch4)",
        color="lightgrey",
        style="rounded",
        fontname="Helvetica"
    )
    c.node("c1", "Targeted CL species")
    c.node("c2", "Derived indices:\nTotal CL, Oxidation index")
    c.node("c3", "Validation phenotypes /\nMechanistic indicators")
    c.edges([("c1","c2"),("c2","c3")])

# ===== Integration node =====
g.node(
    "meta",
    "Meta-analysis / Multi-omics integration",
    shape="diamond",
    style="filled",
    color="lightblue",
    fontname="Helvetica"
)

g.edges([
    ("r6","meta"),
    ("l5","meta"),
    ("c3","meta")
])

# ===== Annotation =====
g.node(
    "note",
    "All analyses performed within\neach batch first, then combined\nusing meta-analytic approaches.",
    shape="note",
    fontname="Helvetica",
    fontsize="9"
)

g.edge("meta", "note", style="dashed", arrowhead="none")

# ===== Render output =====
output_path = g.render(filename="results/figs/figure_6_c")
print(f"✅ SVG generated at: {output_path}")