"""
WOT (Waddington Optimal Transport) fate analysis.
Uses DPT pseudotime bins as proxy for developmental time points.
Run with: D:/miniforge3/python.exe test3_wot_fate.py
Requires: pip install wot scanpy statsmodels (in miniforge3)
"""
import os, warnings, shutil
warnings.filterwarnings("ignore")
import matplotlib; matplotlib.use("Agg")
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt

os.chdir("D:/projects/pf_gametocyte_reproduce")
TMAP_DIR = "outputs/wot_tmaps"


def main():
    print("=== WOT Fate Analysis on P. falciparum Gametocytes ===\n", flush=True)
    print("Loading pf_expr_ready.h5ad ...", flush=True)
    adata = sc.read_h5ad("outputs/pf_expr_ready.h5ad")

    # ── Compute DPT pseudotime ───────────────────────────────────────
    print("Computing DPT pseudotime ...", flush=True)
    early_idx = np.where(adata.obs["annotation"].isin(["C0"]))[0]
    pca_early = adata.obsm["X_pca"][early_idx, 0]
    root_local = np.argmin(np.abs(pca_early - np.median(pca_early)))
    adata.uns["iroot"] = int(early_idx[root_local])
    sc.tl.dpt(adata, n_dcs=10)
    dpt = adata.obs["dpt_pseudotime"].values.copy()

    # ── Bin DPT into 5 discrete time points ─────────────────────────
    n_bins = 5
    bins = pd.qcut(dpt, n_bins, labels=list(range(n_bins)), duplicates="drop")
    adata.obs["wot_day"] = bins.astype(float)
    print("Cells per time bin:")
    print(adata.obs["wot_day"].value_counts().sort_index().to_string())

    # ── Build WOT input AnnData (needs 'day' in obs) ─────────────────
    X = adata.X
    if sp.issparse(X):
        X = X.toarray()

    wot_adata = ad.AnnData(
        X=X.astype(np.float32),
        obs=pd.DataFrame({"day": adata.obs["wot_day"].values,
                          "cluster": adata.obs["annotation"].values},
                         index=adata.obs_names),
        var=pd.DataFrame(index=adata.var_names),
    )

    # ── Run WOT transport maps ────────────────────────────────────────
    import wot
    # WOT 1.0.8 uses deprecated anndata.read — patch for anndata >= 0.10
    if not hasattr(ad, "read"):
        ad.read = ad.read_h5ad  # type: ignore

    print("\nFitting OT transport maps (epsilon=0.05, lambda1=1, lambda2=50) ...", flush=True)
    if os.path.exists(TMAP_DIR):
        shutil.rmtree(TMAP_DIR)
    os.makedirs(TMAP_DIR, exist_ok=True)

    ot_model = wot.ot.OTModel(
        wot_adata,
        day_field="day",
        epsilon=0.05,
        lambda1=1,
        lambda2=50,
        growth_iters=3,
    )
    ot_model.compute_all_transport_maps(tmap_out=os.path.join(TMAP_DIR, "tmaps"))
    print("Transport maps saved.", flush=True)

    # ── Load transport map model ─────────────────────────────────────
    tmap_model = wot.tmap.TransportMapModel.from_directory(os.path.join(TMAP_DIR, "tmaps"))

    # ── Define terminal populations: C1 (female) and C5 (male) ──────
    # Use cells at the latest time bin that belong to C1 or C5
    last_day = wot_adata.obs["day"].max()
    last_obs = wot_adata.obs[wot_adata.obs["day"] == last_day]

    c1_ids = last_obs.index[last_obs["cluster"] == "C1"].tolist()
    c5_ids = last_obs.index[last_obs["cluster"] == "C5"].tolist()

    print(f"\nTerminal cells at t={last_day}: C1={len(c1_ids)}, C5={len(c5_ids)}", flush=True)

    if len(c1_ids) == 0 or len(c5_ids) == 0:
        # Fall back: use ALL C1/C5 cells as terminal
        c1_ids = wot_adata.obs.index[wot_adata.obs["cluster"] == "C1"].tolist()
        c5_ids = wot_adata.obs.index[wot_adata.obs["cluster"] == "C5"].tolist()
        print(f"  (fallback) All C1={len(c1_ids)}, C5={len(c5_ids)}", flush=True)

    # population_from_ids takes each *arg as one population (list of cell ids)
    # returns list of Population objects
    pops_c1 = tmap_model.population_from_ids(c1_ids, at_time=last_day)
    pops_c5 = tmap_model.population_from_ids(c5_ids, at_time=last_day)
    pop_c1 = pops_c1[0]; pop_c1.name = "C1_female"
    pop_c5 = pops_c5[0]; pop_c5.name = "C5_male"

    # ── Compute fates: pass BOTH populations together → one AnnData ──
    print("Computing fate probabilities ...", flush=True)
    fates_ds = tmap_model.fates([pop_c1, pop_c5])
    # fates_ds.X shape: (all_cells_up_to_last_day, 2)
    # fates_ds.var_names: ["C1_female", "C5_male"]
    print(f"Fate AnnData shape: {fates_ds.shape}", flush=True)

    # Only keep cells from time 0 (earliest bin) for visualisation
    early_fates = fates_ds[fates_ds.obs["day"] == 0].to_df()
    early_fates.index = early_fates.index  # already cell barcodes

    fate_df = fates_ds.to_df()
    fate_df["cluster"] = wot_adata.obs.loc[fate_df.index, "cluster"]

    print("\nMean fate probabilities per cluster:")
    summary = fate_df.groupby("cluster")[["C1_female", "C5_male"]].mean()
    print(summary.round(4).to_string())

    # ── Plot ─────────────────────────────────────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # UMAP coloured by C1 fate
    umap = adata.obsm["X_umap"]
    fate_c1 = fate_df.reindex(adata.obs_names)["C1_female"].fillna(0).values
    fate_c5 = fate_df.reindex(adata.obs_names)["C5_male"].fillna(0).values

    sc1 = axes[0].scatter(umap[:, 0], umap[:, 1], c=fate_c1, s=3, cmap="Reds", vmin=0, vmax=fate_c1.max())
    axes[0].set_title("WOT: Fate Probability → C1 (Female)"); axes[0].axis("off")
    plt.colorbar(sc1, ax=axes[0])

    sc2 = axes[1].scatter(umap[:, 0], umap[:, 1], c=fate_c5, s=3, cmap="Blues", vmin=0, vmax=fate_c5.max())
    axes[1].set_title("WOT: Fate Probability → C5 (Male)"); axes[1].axis("off")
    plt.colorbar(sc2, ax=axes[1])

    plt.tight_layout()
    plt.savefig("figures/wot_fate_umap.png", dpi=150)
    print("\nSaved figures/wot_fate_umap.png", flush=True)

    # Bar chart
    fig2, ax = plt.subplots(figsize=(10, 4))
    summary.plot(kind="bar", ax=ax, color=["#e74c3c", "#2980b9"], edgecolor="black")
    ax.set_title("WOT: Mean Fate Probabilities per Cluster")
    ax.set_xlabel("Cluster"); ax.set_ylabel("Fate Probability")
    ax.set_xticklabels(summary.index, rotation=0)
    ax.legend(["C1_female", "C5_male"])
    plt.tight_layout()
    plt.savefig("figures/wot_fate_bar.png", dpi=150)
    print("Saved figures/wot_fate_bar.png", flush=True)

    print("\n[DONE]", flush=True)


if __name__ == "__main__":
    main()
