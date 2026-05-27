"""
Step 2 of mLLMCelltype test:
Use mLLMCelltype to annotate P. falciparum gametocyte clusters via LLM.
Run with: D:/miniforge3/python.exe test2_mllmcelltype.py
Requires: pip install mllmcelltype anthropic (in miniforge3)
"""
import os, json, warnings
warnings.filterwarnings("ignore")
os.chdir("D:/projects/pf_gametocyte_reproduce")

ANTHROPIC_API_KEY = os.environ.get("ANTHROPIC_API_KEY", "")


def main():
    # Load markers from Step 1
    with open("outputs/cluster_markers.json") as f:
        markers = json.load(f)

    print(f"Loaded markers for {len(markers)} clusters\n", flush=True)

    if not ANTHROPIC_API_KEY:
        print("ERROR: ANTHROPIC_API_KEY not set. Run:  $env:ANTHROPIC_API_KEY='sk-ant-...'", flush=True)
        return

    # Method A: use mllmcelltype package if installed
    try:
        from mllmcelltype import annotate_clusters
        print("=== Using mLLMCelltype package (v2.0.x) ===", flush=True)

        results = annotate_clusters(
            marker_genes=markers,
            species="Plasmodium falciparum",
            tissue="gametocyte",
            provider="anthropic",
            model="claude-haiku-4-5-20251001",
            api_key=ANTHROPIC_API_KEY,
        )
        print("\nAnnotation results:")
        for cluster, annotation in results.items():
            print(f"  {cluster}: {annotation}")

        print("\nAnnotation results:")
        for cluster, annotation in results.items():
            print(f"  {cluster}: {annotation}")

    except Exception as e:
        print(f"mllmcelltype failed ({e}), falling back to direct Claude API ...", flush=True)
        _annotate_direct(markers)


def _annotate_direct(markers: dict):
    """Direct Claude API fallback — same logic as mLLMCelltype but simplified."""
    import anthropic

    client = anthropic.Anthropic(api_key=ANTHROPIC_API_KEY)
    context = (
        "You are an expert in Plasmodium falciparum biology. "
        "The single-cell RNA-seq data contains gametocyte clusters C0–C9. "
        "Known biology: C0/C4/C9 = early uncommitted, C1 = mature female, "
        "C5 = mature male, C2/C3/C6/C7/C8 = intermediate. "
        "AP2-G (PF3D7_1222600) drives female, AP2-O4 (PF3D7_1350900) drives male."
    )

    results = {}
    for cluster, genes in markers.items():
        prompt = (
            f"{context}\n\n"
            f"Cluster {cluster} top marker genes: {', '.join(genes)}\n\n"
            "Based on these markers, predict: (1) cell type / developmental stage, "
            "(2) sex (female/male/uncommitted), (3) confidence (high/medium/low). "
            "Reply in one line: CLUSTER | CELL_TYPE | SEX | CONFIDENCE"
        )

        response = client.messages.create(
            model="claude-haiku-4-5-20251001",
            max_tokens=200,
            messages=[{"role": "user", "content": prompt}],
        )
        answer = response.content[0].text.strip()
        results[cluster] = answer
        print(f"  {answer}", flush=True)

    # Save results
    with open("outputs/mllmcelltype_annotations.json", "w") as f:
        json.dump(results, f, indent=2)
    print("\nSaved to outputs/mllmcelltype_annotations.json", flush=True)


if __name__ == "__main__":
    main()
