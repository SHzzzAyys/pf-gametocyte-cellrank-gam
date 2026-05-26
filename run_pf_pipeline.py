"""
P. falciparum gametocyte — reproduce Fig. 5A (Mohammed et al., Nat Commun 2024).

Data sources:
  - scenic_out_020223.h5ad (Zenodo 7652581): UMAP / PCA / neighbors / cell annotations
  - filtered_matrix.tsv.gz (GEO GSE226145): raw UMI counts for expression modelling
  - cells_metadata.csv (GEO): cell-level metadata (annotation C0–C9, cell_class…)

Route taken (no velocity in deposited data):
  ConnectivityKernel → GPCCA → set C1/C5 terminal → absorption probs → GAM
"""
import os
import gzip
import warnings
warnings.filterwarnings("ignore")

import matplotlib
matplotlib.use("Agg")

import numpy as np
import pandas as pd
import scipy.sparse as sp
import scanpy as sc
import cellrank as cr
from cellrank.tl.kernels import ConnectivityKernel
from cellrank.tl.estimators import GPCCA
from cellrank.ul.models import GAM
import matplotlib.pyplot as plt

DATA_DIR = "data"
OUTDIR   = "outputs"
FIGDIR   = "figures"
os.makedirs(OUTDIR, exist_ok=True)
os.makedirs(FIGDIR, exist_ok=True)

H5AD    = os.path.join(DATA_DIR, "scenic_out_020223.h5ad")
MATRIX  = os.path.join(DATA_DIR, "filtered_matrix.tsv.gz")
META    = os.path.join(DATA_DIR, "cells_metadata.csv")
CACHE   = os.path.join(OUTDIR,  "pf_expr_ready.h5ad")

APAP2_FEMALE = ["PF3D7_1222600", "PF3D7_0404100", "PF3D7_1408200"]  # AP2-G, female ApiAP2
APAP2_MALE   = ["PF3D7_1350900", "PF3D7_1239200"]                   # AP2-O4, male ApiAP2


def log(msg):
    print(f"\n[STAGE] {msg}\n" + "=" * 60, flush=True)


# ── Step 1: build AnnData with log-normalised expression ────────────────────

def build_adata():
    if os.path.exists(CACHE):
        log("1. Loading cached AnnData")
        return sc.read_h5ad(CACHE)

    log("1. Building AnnData: raw counts + scenic metadata")

    # 1a. Load raw UMI counts (genes × cells) → cells × genes
    print("  Loading filtered_matrix.tsv.gz …", flush=True)
    with gzip.open(MATRIX, "rt") as f:
        mat_df = pd.read_csv(f, sep=",", index_col=0)
    mat_df = mat_df.T   # → cells × genes
    print(f"  Matrix: {mat_df.shape[0]} cells × {mat_df.shape[1]} genes")

    # 1b. Load scenic h5ad for structure (UMAP, PCA, neighbors, obs labels)
    print("  Loading scenic_out_020223.h5ad …", flush=True)
    scenic = sc.read_h5ad(H5AD)
    print(f"  scenic: {scenic.shape}")

    # Align cells (should be identical order, but be safe)
    common_cells = mat_df.index.intersection(scenic.obs_names)
    mat_df  = mat_df.loc[common_cells]
    scenic  = scenic[common_cells].copy()

    # 1c. Build AnnData: log-normalised counts as X
    adata = sc.AnnData(
        X    = sp.csr_matrix(mat_df.values.astype(np.float32)),
        obs  = scenic.obs.copy(),
        var  = scenic.var.copy(),
        obsm = dict(scenic.obsm),
        obsp = dict(scenic.obsp),
        uns  = dict(scenic.uns),
    )
    # Store raw counts, then normalise
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    print(f"  X after log-norm: min={adata.X.min():.3f}  max={adata.X.max():.3f}")

    # 1d. Ensure annotation is a string categorical (needed by GPCCA cluster_key)
    adata.obs["annotation"] = pd.Categorical(adata.obs["annotation"].astype(str))
    adata.obs["clusters"]   = pd.Categorical(adata.obs["annotation"])  # alias

    adata.write_h5ad(CACHE)
    print(f"  Saved → {CACHE}", flush=True)
    return adata


# ── Step 2: DPT pseudotime (replaces scVelo latent_time) ────────────────────

def compute_dpt(adata):
    log("2. Compute DPT pseudotime")
    # Root cell: an early gametocyte (C0) near the centre of the early cluster
    early_idx = np.where(adata.obs["annotation"].isin(["C0"]))[0]
    # Choose the cell whose PC1 is closest to the median of C0
    pca_early = adata.obsm["X_pca"][early_idx, 0]
    median_pc1 = np.median(pca_early)
    root_local = np.argmin(np.abs(pca_early - median_pc1))
    root_idx   = int(early_idx[root_local])
    adata.uns["iroot"] = root_idx
    print(f"  Root cell index: {root_idx}  ({adata.obs['annotation'].iloc[root_idx]})")

    sc.tl.dpt(adata, n_dcs=10)
    adata.obs["latent_time"] = adata.obs["dpt_pseudotime"]
    print(f"  DPT range: {adata.obs['latent_time'].min():.3f} – "
          f"{adata.obs['latent_time'].max():.3f}")
    return adata


# ── Step 3: CellRank — ConnectivityKernel + GPCCA ───────────────────────────

def run_cellrank(adata):
    log("3. CellRank: ConnectivityKernel → GPCCA → C1/C5 terminal states")

    ck = ConnectivityKernel(adata)
    ck.compute_transition_matrix()

    g = GPCCA(ck)
    g.compute_schur(n_components=20)

    # Try increasing n_states until C1 and C5 both appear as separate macrostates
    best_g = None
    for n in [10, 8, 6, 4, 3]:
        try:
            g.compute_macrostates(n_states=n, cluster_key="clusters")
            macros = g.macrostates.cat.categories.tolist()
            print(f"  n_states={n}: {macros}")
            c1_found = any("C1" in m for m in macros)
            c5_found = any("C5" in m for m in macros)
            if c1_found and c5_found:
                best_g = g
                break
        except Exception as e:
            print(f"  n_states={n} failed: {e}")

    if best_g is None:
        print("  [warn] Could not separate C1/C5; using last computed macrostates")
        best_g = g
        macros = g.macrostates.cat.categories.tolist()

    g = best_g
    macros = g.macrostates.cat.categories.tolist()

    # Plot macrostates
    try:
        g.plot_macrostates(discrete=True, basis="umap", show=False,
                           save="pf_macrostates.png")
        plt.close("all")
    except Exception as e:
        print(f"  [warn] macrostate plot: {e}")

    # Set terminal states: pick the macrostate names containing C1 / C5
    c1_name = next((m for m in macros if "C1" in m), None)
    c5_name = next((m for m in macros if "C5" in m), None)
    print(f"  C1 macrostate: {c1_name}")
    print(f"  C5 macrostate: {c5_name}")

    terminals = [t for t in [c1_name, c5_name] if t is not None]
    if len(terminals) >= 1:
        g.set_terminal_states_from_macrostates(names=terminals)
    else:
        print("  [warn] Falling back to stability-based terminal states")
        g.compute_terminal_states(method="stability")

    terminal_names = g.terminal_states.cat.categories.tolist()
    print(f"  Terminal states: {terminal_names}")

    try:
        g.plot_terminal_states(discrete=True, basis="umap", show=False,
                               save="pf_terminal_states.png")
        plt.close("all")
    except Exception as e:
        print(f"  [warn] terminal state plot: {e}")

    g.compute_absorption_probabilities(solver="direct")
    try:
        g.plot_absorption_probabilities(basis="umap", show=False,
                                        save="pf_abs_probs.png")
        plt.close("all")
    except Exception as e:
        print(f"  [warn] abs_prob plot: {e}")

    drivers = g.compute_lineage_drivers(use_raw=False)
    drivers.to_csv(os.path.join(OUTDIR, "pf_lineage_drivers.csv"))
    print("  Driver columns (first 8):", drivers.columns.tolist()[:8])
    return adata, g, drivers, terminal_names


# ── Step 4: GAM — ApiAP2 gene trends (Fig. 5A) ──────────────────────────────

def run_gam(adata, terminal_names):
    log("4. GAM gene_trends — ApiAP2 TF expression (Fig. 5A)")

    present_f = [g for g in APAP2_FEMALE if g in adata.var_names]
    present_m = [g for g in APAP2_MALE   if g in adata.var_names]
    print("  Female ApiAP2:", present_f)
    print("  Male ApiAP2:",   present_m)

    model = GAM(adata, distribution="normal", link="identity")

    lin_female = next((t for t in terminal_names if "C1" in t), None)
    lin_male   = next((t for t in terminal_names if "C5" in t), None)
    print(f"  Female lineage: {lin_female},  Male lineage: {lin_male}")

    lineages = [t for t in [lin_female, lin_male] if t is not None]
    all_genes = list(dict.fromkeys(present_f + present_m))

    if not all_genes:
        print("  [error] No ApiAP2 genes found in adata.var_names")
        return

    # ── Fig. 5A combined: both lineages, all ApiAP2 genes ─────────────────
    print("  Plotting combined trends (Fig. 5A style)…", flush=True)
    cr.pl.gene_trends(
        adata,
        model=model,
        genes=all_genes,
        lineages=lineages,
        time_key="latent_time",
        data_key="X",
        same_plot=True,
        hide_cells=True,
        n_test_points=200,
        save="fig5A_combined_apap2.pdf",
    )
    plt.close("all")

    # ── Female lineage only ────────────────────────────────────────────────
    if lin_female and present_f:
        print("  Plotting female trends…", flush=True)
        cr.pl.gene_trends(
            adata,
            model=model,
            genes=present_f,
            lineages=[lin_female],
            time_key="latent_time",
            data_key="X",
            same_plot=False,
            hide_cells=False,
            n_test_points=200,
            save="fig5A_female_apap2.pdf",
        )
        plt.close("all")

    # ── Male lineage only ──────────────────────────────────────────────────
    if lin_male and present_m:
        print("  Plotting male trends…", flush=True)
        cr.pl.gene_trends(
            adata,
            model=model,
            genes=present_m,
            lineages=[lin_male],
            time_key="latent_time",
            data_key="X",
            same_plot=False,
            hide_cells=False,
            n_test_points=200,
            save="fig5A_male_apap2.pdf",
        )
        plt.close("all")

    print("  [GAM] All figures saved to", os.path.abspath(FIGDIR))


# ── Main ────────────────────────────────────────────────────────────────────

def main():
    adata = build_adata()

    print(f"\nData overview:")
    print(f"  Cells: {adata.n_obs}  Genes: {adata.n_vars}")
    print(f"  Clusters: {sorted(adata.obs['annotation'].unique().tolist())}")
    print(f"  Cell classes: {sorted(adata.obs['cell_class'].unique().tolist())}")
    print(f"  UMAP: {adata.obsm['X_umap'].shape}")
    print(f"  Connectivity: {adata.obsp['connectivities'].shape}")

    adata = compute_dpt(adata)
    adata, g, drivers, terminal_names = run_cellrank(adata)
    run_gam(adata, terminal_names)

    print("\n[DONE] outputs →", os.path.abspath(OUTDIR))
    print("[DONE] figures →", os.path.abspath(FIGDIR))


if __name__ == "__main__":
    main()
