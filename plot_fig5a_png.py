"""Re-run GAM plots from scratch to generate PNG versions of Fig. 5A."""
import os, warnings
warnings.filterwarnings("ignore")
import matplotlib; matplotlib.use("Agg")
import numpy as np, pandas as pd, scipy.sparse as sp, scanpy as sc
import cellrank as cr
from cellrank.tl.kernels import ConnectivityKernel
from cellrank.tl.estimators import GPCCA
from cellrank.ul.models import GAM
import matplotlib.pyplot as plt

APAP2_FEMALE = ["PF3D7_1222600", "PF3D7_0404100", "PF3D7_1408200"]
APAP2_MALE   = ["PF3D7_1350900", "PF3D7_1239200"]


def main():
    os.chdir("D:/projects/pf_gametocyte_reproduce")

    print("Loading cached AnnData...", flush=True)
    adata = sc.read_h5ad("outputs/pf_expr_ready.h5ad")
    print(adata)

    # Recompute DPT
    print("Computing DPT...", flush=True)
    early_idx = np.where(adata.obs["annotation"].isin(["C0"]))[0]
    pca_early = adata.obsm["X_pca"][early_idx, 0]
    root_local = np.argmin(np.abs(pca_early - np.median(pca_early)))
    adata.uns["iroot"] = int(early_idx[root_local])
    sc.tl.dpt(adata, n_dcs=10)
    adata.obs["latent_time"] = adata.obs["dpt_pseudotime"]
    print(f"DPT: {adata.obs['latent_time'].min():.3f} – {adata.obs['latent_time'].max():.3f}")

    # CellRank
    print("Running ConnectivityKernel + GPCCA...", flush=True)
    ck = ConnectivityKernel(adata)
    ck.compute_transition_matrix()
    g = GPCCA(ck)
    g.compute_schur(n_components=20)
    g.compute_macrostates(n_states=10, cluster_key="clusters")
    macros = g.macrostates.cat.categories.tolist()
    print("macrostates:", macros)

    c1 = next((m for m in macros if "C1" in m), None)
    c5 = next((m for m in macros if "C5" in m), None)
    terminals = [t for t in [c1, c5] if t]
    print("terminals:", terminals)
    g.set_terminal_states_from_macrostates(names=terminals)
    g.compute_absorption_probabilities(solver="direct")
    terminal_names = g.terminal_states.cat.categories.tolist()
    print("terminal states:", terminal_names)

    # GAM
    print("Fitting GAM models...", flush=True)
    model = GAM(adata, distribution="normal", link="identity")
    lin_f = next((t for t in terminal_names if "C1" in t), None)
    lin_m = next((t for t in terminal_names if "C5" in t), None)
    lineages = [t for t in [lin_f, lin_m] if t]
    all_genes = list(dict.fromkeys(APAP2_FEMALE + APAP2_MALE))

    # Combined Fig. 5A as PNG
    print("Saving fig5A_combined_apap2.png ...", flush=True)
    cr.pl.gene_trends(
        adata, model=model, genes=all_genes, lineages=lineages,
        time_key="latent_time", data_key="X",
        same_plot=True, hide_cells=True, n_test_points=200,
        save="fig5A_combined_apap2.png",
    )
    plt.close("all")

    # Female only
    if lin_f:
        print("Saving fig5A_female_apap2.png ...", flush=True)
        cr.pl.gene_trends(
            adata, model=model, genes=APAP2_FEMALE, lineages=[lin_f],
            time_key="latent_time", data_key="X",
            same_plot=False, hide_cells=False, n_test_points=200,
            save="fig5A_female_apap2.png",
        )
        plt.close("all")

    # Male only
    if lin_m:
        print("Saving fig5A_male_apap2.png ...", flush=True)
        cr.pl.gene_trends(
            adata, model=model, genes=APAP2_MALE, lineages=[lin_m],
            time_key="latent_time", data_key="X",
            same_plot=False, hide_cells=False, n_test_points=200,
            save="fig5A_male_apap2.png",
        )
        plt.close("all")

    print("[DONE] PNGs saved to figures/", flush=True)


if __name__ == "__main__":
    main()
