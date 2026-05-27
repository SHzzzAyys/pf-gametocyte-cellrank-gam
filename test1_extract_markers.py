"""
Step 1 of mLLMCelltype test:
Extract top marker genes per cluster from pf_expr_ready.h5ad
Save to markers.json for mLLMCelltype annotation.
Run with: .venv/Scripts/python.exe test1_extract_markers.py
"""
import os, json, warnings
warnings.filterwarnings("ignore")
import numpy as np
import scanpy as sc

os.chdir("D:/projects/pf_gametocyte_reproduce")


def main():
    print("Loading pf_expr_ready.h5ad ...", flush=True)
    adata = sc.read_h5ad("outputs/pf_expr_ready.h5ad")
    print(adata)
    print("Clusters:", sorted(adata.obs["clusters"].cat.categories.tolist()), flush=True)

    # Rank genes per cluster (Wilcoxon)
    print("Computing marker genes (Wilcoxon) ...", flush=True)
    sc.tl.rank_genes_groups(adata, groupby="clusters", method="wilcoxon", n_genes=20)

    clusters = adata.obs["clusters"].cat.categories.tolist()
    markers = {}
    for cl in clusters:
        genes = sc.get.rank_genes_groups_df(adata, group=cl)
        # top 10 by score, filter out low-scoring
        top = genes.head(10)["names"].tolist()
        markers[cl] = top
        print(f"  {cl}: {top[:5]} ...", flush=True)

    with open("outputs/cluster_markers.json", "w") as f:
        json.dump(markers, f, indent=2)
    print("\nSaved to outputs/cluster_markers.json", flush=True)

    # Also print a formatted summary for the annotation prompt
    print("\n=== Marker gene summary for LLM annotation ===")
    for cl, genes in markers.items():
        print(f"Cluster {cl}: {', '.join(genes)}")


if __name__ == "__main__":
    main()
