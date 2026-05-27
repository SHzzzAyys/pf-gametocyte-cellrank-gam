"""
pySCENIC regulon activity analysis.
The scenic_out_020223.h5ad already contains AUCell scores (X_aucell)
for 8 ApiAP2 regulons — no need to re-run SCENIC.
This script visualizes regulon activity per cluster and on UMAP.
Run with: .venv/Scripts/python.exe test4_scenic_regulon.py
"""
import os, warnings
warnings.filterwarnings("ignore")
import matplotlib; matplotlib.use("Agg")
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.gridspec import GridSpec

os.chdir("D:/projects/pf_gametocyte_reproduce")

# Known biology of these regulons
REGULON_LABELS = {
    "Regulon(PF3D7_0604100(-))": "PF3D7_0604100\n(AP2-G2-like, male↑)",
    "Regulon(PF3D7_0611200(+))": "PF3D7_0611200\n(early gametocyte)",
    "Regulon(PF3D7_0613800(+))": "PF3D7_0613800\n(intermediate)",
    "Regulon(PF3D7_1007700(-))": "PF3D7_1007700\n(C8-specific)",
    "Regulon(PF3D7_1305200(+))": "PF3D7_1305200\n(broad expression)",
    "Regulon(PF3D7_1342900(+))": "PF3D7_1342900\n(stage III-IV)",
    "Regulon(PF3D7_1350900(+))": "PF3D7_1350900\n(AP2-O4, male↑)",
    "Regulon(PF3D7_1429200(+))": "PF3D7_1429200\n(intermediate)",
}

CLUSTER_ORDER = ["C0", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9"]
CLUSTER_COLORS = {
    "C0": "#aec6cf", "C1": "#e74c3c", "C2": "#f39c12", "C3": "#27ae60",
    "C4": "#2980b9", "C5": "#8e44ad", "C6": "#16a085", "C7": "#d35400",
    "C8": "#7f8c8d", "C9": "#2c3e50",
}


def main():
    print("Loading pf_expr_ready.h5ad ...", flush=True)
    adata = sc.read_h5ad("outputs/pf_expr_ready.h5ad")

    regulon_names = list(adata.uns["aucell"]["regulon_names"])
    X_aucell = adata.obsm["X_aucell"]
    auc_df = pd.DataFrame(X_aucell, index=adata.obs_names, columns=regulon_names)
    auc_df["cluster"] = adata.obs["annotation"].values

    summary = auc_df.groupby("cluster")[regulon_names].mean()
    summary = summary.loc[CLUSTER_ORDER]

    # ── Fig 1: Heatmap — regulon activity per cluster ────────────────
    print("Plotting regulon activity heatmap ...", flush=True)
    short_labels = [REGULON_LABELS.get(r, r.split("(")[1].rstrip(")")) for r in regulon_names]

    fig, ax = plt.subplots(figsize=(11, 5))
    data = summary.values
    # Z-score each regulon across clusters for better contrast
    from scipy.stats import zscore
    data_z = zscore(data, axis=0)

    im = ax.imshow(data_z.T, aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)
    ax.set_xticks(range(len(CLUSTER_ORDER)))
    ax.set_xticklabels(CLUSTER_ORDER, fontsize=11)
    ax.set_yticks(range(len(regulon_names)))
    ax.set_yticklabels(short_labels, fontsize=8)
    ax.set_xlabel("Cluster", fontsize=12)
    ax.set_title("pySCENIC AUCell Regulon Activity\n(z-score across clusters)", fontsize=13)
    plt.colorbar(im, ax=ax, label="z-score", shrink=0.8)

    # Annotate female (C1) and male (C5)
    ax.axvline(x=0.5, color="gray", lw=0.5, ls="--")
    for xi, cl in enumerate(CLUSTER_ORDER):
        if cl == "C1":
            ax.add_patch(plt.Rectangle((xi-0.5, -0.5), 1, len(regulon_names),
                         fill=True, color="#e74c3c", alpha=0.08, zorder=0))
            ax.text(xi, len(regulon_names)-0.1, "♀", ha="center", fontsize=14, color="#e74c3c")
        elif cl == "C5":
            ax.add_patch(plt.Rectangle((xi-0.5, -0.5), 1, len(regulon_names),
                         fill=True, color="#8e44ad", alpha=0.08, zorder=0))
            ax.text(xi, len(regulon_names)-0.1, "♂", ha="center", fontsize=14, color="#8e44ad")

    plt.tight_layout()
    plt.savefig("figures/scenic_regulon_heatmap.png", dpi=150)
    plt.close()
    print("Saved figures/scenic_regulon_heatmap.png", flush=True)

    # ── Fig 2: Violin plots for top sex-specific regulons ────────────
    print("Plotting violin plots ...", flush=True)
    # Most sex-specific: PF3D7_0604100 (male C5 high) and PF3D7_1350900 (AP2-O4, male)
    top_regulons = [
        "Regulon(PF3D7_0604100(-))",
        "Regulon(PF3D7_1350900(+))",
        "Regulon(PF3D7_1305200(+))",
        "Regulon(PF3D7_0611200(+))",
    ]

    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    axes = axes.flatten()

    for ax, reg in zip(axes, top_regulons):
        scores_per_cluster = [auc_df.loc[auc_df["cluster"] == cl, reg].values
                              for cl in CLUSTER_ORDER]
        colors = [CLUSTER_COLORS[cl] for cl in CLUSTER_ORDER]
        parts = ax.violinplot(scores_per_cluster, positions=range(len(CLUSTER_ORDER)),
                              showmedians=True, showextrema=False)
        for pc, color in zip(parts["bodies"], colors):
            pc.set_facecolor(color); pc.set_alpha(0.7)
        parts["cmedians"].set_color("black"); parts["cmedians"].set_linewidth(1.5)

        ax.set_xticks(range(len(CLUSTER_ORDER)))
        ax.set_xticklabels(CLUSTER_ORDER, fontsize=9)
        ax.set_title(REGULON_LABELS.get(reg, reg), fontsize=9)
        ax.set_ylabel("AUCell score")

        # Highlight C1 / C5
        ax.axvspan(CLUSTER_ORDER.index("C1")-0.4, CLUSTER_ORDER.index("C1")+0.4,
                   alpha=0.1, color="#e74c3c")
        ax.axvspan(CLUSTER_ORDER.index("C5")-0.4, CLUSTER_ORDER.index("C5")+0.4,
                   alpha=0.1, color="#8e44ad")

    fig.suptitle("pySCENIC: Top Regulon Activity by Cluster\n(red=C1 female, purple=C5 male)",
                 fontsize=13)
    plt.tight_layout()
    plt.savefig("figures/scenic_regulon_violin.png", dpi=150)
    plt.close()
    print("Saved figures/scenic_regulon_violin.png", flush=True)

    # ── Fig 3: UMAP colored by regulon activity ──────────────────────
    print("Plotting UMAP ...", flush=True)
    umap = adata.obsm["X_umap"]

    fig, axes = plt.subplots(2, 4, figsize=(16, 8))
    axes = axes.flatten()

    for ax, reg in zip(axes, regulon_names):
        scores = auc_df[reg].values
        sc_plot = ax.scatter(umap[:, 0], umap[:, 1], c=scores, s=2,
                             cmap="YlOrRd", vmin=0, vmax=np.percentile(scores, 99))
        ax.set_title(REGULON_LABELS.get(reg, reg), fontsize=7)
        ax.axis("off")
        plt.colorbar(sc_plot, ax=ax, shrink=0.7)

    fig.suptitle("pySCENIC AUCell Regulon Activity on UMAP", fontsize=13)
    plt.tight_layout()
    plt.savefig("figures/scenic_regulon_umap.png", dpi=150)
    plt.close()
    print("Saved figures/scenic_regulon_umap.png", flush=True)

    # ── Summary table ────────────────────────────────────────────────
    print("\n=== Regulon Activity Summary ===")
    print(f"{'Regulon':<35} {'C1(female)':>12} {'C5(male)':>10} {'C5/C1 ratio':>12} {'Top cluster':>12}")
    print("-" * 85)
    for reg in regulon_names:
        c1_val = summary.loc["C1", reg]
        c5_val = summary.loc["C5", reg]
        ratio  = c5_val / (c1_val + 1e-6)
        top_cl = summary[reg].idxmax()
        short  = reg.replace("Regulon(", "").rstrip(")")
        print(f"{short:<35} {c1_val:>12.4f} {c5_val:>10.4f} {ratio:>12.2f} {top_cl:>12}")

    print("\n[DONE] Saved 3 figures to figures/", flush=True)


if __name__ == "__main__":
    main()
