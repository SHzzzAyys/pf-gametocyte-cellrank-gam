"""
Full CellRank 1.x + scVelo + GAM pipeline — validation on pancreas dataset.

Mirrors the analysis architecture of Mohammed et al. (Nat Commun 2024)
on a different biological substrate (mouse pancreas endocrinogenesis).

Stages:
    1. Load pancreas (built-in, has spliced/unspliced)
    2. scVelo: moments -> velocity -> velocity_graph -> recover_latent_time
    3. CellRank 1.5.x: VelocityKernel -> GPCCA -> terminal states -> abs probs -> drivers
    4. GAM: cellrank.ul.models.GAM + cr.pl.gene_trends
"""
import os
import warnings
warnings.filterwarnings("ignore")

import matplotlib
matplotlib.use("Agg")  # non-interactive backend, avoid blocking on plt.show()

import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import cellrank as cr
from cellrank.tl.kernels import VelocityKernel
from cellrank.tl.estimators import GPCCA
from cellrank.ul.models import GAM
import matplotlib.pyplot as plt

OUTDIR = "outputs"
CACHE_VELO = os.path.join(OUTDIR, "pancreas_post_velocity.h5ad")


def log(msg):
    print(f"\n[STAGE] {msg}\n" + "=" * 60, flush=True)


def stage_velocity():
    """Stage 1+2: load data, run scVelo velocity + latent_time."""
    if os.path.exists(CACHE_VELO):
        log("1+2/4 Cached velocity result found, loading")
        return sc.read_h5ad(CACHE_VELO)

    log("1/4 Load pancreas dataset")
    adata = scv.datasets.pancreas()
    print(adata)
    print("clusters:", adata.obs["clusters"].cat.categories.tolist())

    log("2/4 scVelo: moments -> velocity -> latent_time")
    scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

    scv.tl.velocity(adata, mode="stochastic")
    scv.tl.velocity_graph(adata)

    # Latent time needs dynamical recovery; use n_jobs=1 on Windows to avoid
    # multiprocessing bootstrap errors
    scv.tl.recover_dynamics(adata, n_jobs=1)
    scv.tl.latent_time(adata)
    print(
        "latent_time range:",
        float(adata.obs["latent_time"].min()),
        float(adata.obs["latent_time"].max()),
    )

    adata.write_h5ad(CACHE_VELO)
    print(f"[checkpoint] saved {CACHE_VELO}", flush=True)

    try:
        scv.pl.velocity_embedding_stream(
            adata, basis="umap", color="clusters",
            save="01_velocity_stream.png", show=False, dpi=150,
        )
        plt.close("all")
    except Exception as e:
        print(f"[warn] stream plot failed (non-fatal): {e}", flush=True)

    return adata


def stage_cellrank(adata):
    """Stage 3: CellRank GPCCA, terminal states, abs probs, lineage drivers."""
    log("3/4 CellRank 1.5.x: GPCCA + lineage drivers")
    vk = VelocityKernel(adata)
    vk.compute_transition_matrix()

    g = GPCCA(vk)
    g.compute_schur(n_components=20)
    g.compute_macrostates(n_states=3, cluster_key="clusters")
    print("macrostates:", g.macrostates.cat.categories.tolist())

    try:
        g.plot_macrostates(discrete=True, basis="umap", show=False,
                           save="02_macrostates.png")
    except Exception as e:
        print(f"[warn] plot_macrostates failed (non-fatal): {e}", flush=True)

    # For pancreas, Alpha and Beta are the true endocrine terminal states;
    # Ductal is the progenitor. Pin terminals manually so absorption probs
    # measure fate towards these two lineages only.
    macrostate_names = g.macrostates.cat.categories.tolist()
    endpoints = [n for n in macrostate_names if n in {"Alpha", "Beta"}]
    if len(endpoints) >= 2:
        g.set_terminal_states_from_macrostates(names=endpoints)
    else:
        g.compute_terminal_states(method="stability")
    terminal_names = g.terminal_states.cat.categories.tolist()
    print("terminal states:", terminal_names, flush=True)
    try:
        g.plot_terminal_states(discrete=True, basis="umap", show=False,
                               save="03_terminal_states.png")
    except Exception as e:
        print(f"[warn] plot_terminal_states failed (non-fatal): {e}", flush=True)

    # direct solver avoids gmres numerical issues on a 3.7k-cell matrix
    g.compute_absorption_probabilities(solver="direct")
    try:
        g.plot_absorption_probabilities(basis="umap", show=False,
                                        save="04_abs_probs.png")
    except Exception as e:
        print(f"[warn] plot_absorption_probabilities failed (non-fatal): {e}",
              flush=True)

    drivers = g.compute_lineage_drivers(use_raw=False)
    drivers.to_csv(os.path.join(OUTDIR, "lineage_drivers.csv"))
    print("\ndriver columns:", drivers.columns.tolist())
    return adata, g, drivers, terminal_names


def stage_gam(adata, drivers, terminal_names):
    """Stage 4: GAM gene_trends along latent_time."""
    log("4/4 GAM gene_trends along latent_time")
    # default GAM uses gamma family + log link; log-normalized expression has
    # many zeros, which makes log(0) blow up. Use normal + identity link.
    model = GAM(adata, distribution="normal", link="identity")

    # cellrank 1.5.x returns flat columns like 'Alpha_corr', 'Alpha_pval', ...
    top_genes = []
    for lin in terminal_names:
        corr_col = f"{lin}_corr"
        if corr_col not in drivers.columns:
            print(f"[warn] no driver column for {lin}", flush=True)
            continue
        top = drivers.sort_values(by=corr_col, ascending=False).index[:3].tolist()
        top_genes.extend(top)
    top_genes = list(dict.fromkeys(top_genes))
    print("Plotting GAM trends for:", top_genes)

    cr.pl.gene_trends(
        adata,
        model=model,
        genes=top_genes,
        lineages=terminal_names,
        time_key="latent_time",
        data_key="X",
        same_plot=True,
        hide_cells=True,
        n_test_points=200,
        save="05_gene_trends_combined.pdf",
    )

    cr.pl.gene_trends(
        adata,
        model=model,
        genes=top_genes[:4],
        lineages=terminal_names,
        time_key="latent_time",
        data_key="X",
        same_plot=False,
        hide_cells=False,
        n_test_points=200,
        save="06_gene_trends_per_lineage.pdf",
    )


def main():
    adata = stage_velocity()
    adata, g, drivers, terminal_names = stage_cellrank(adata)
    stage_gam(adata, drivers, terminal_names)

    adata.write_h5ad(os.path.join(OUTDIR, "pancreas_cellrank.h5ad"))
    print("\n[DONE] outputs at:", os.path.abspath(OUTDIR), flush=True)


if __name__ == "__main__":
    main()
