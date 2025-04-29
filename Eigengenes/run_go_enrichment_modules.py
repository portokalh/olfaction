
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GO enrichment pipeline for WGCNA modules.
- Loads gene_module_assignments.csv
- Runs GO:BP enrichment using gprofiler
- Outputs CSVs and a heatmap of enriched terms across modules
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from gprofiler import GProfiler

# === SETUP ===
base_dir = os.path.dirname(__file__)
results_dir = os.path.join(base_dir, "results_output")
module_file = os.path.join(results_dir, "gene_module_assignments.csv")
os.makedirs(results_dir, exist_ok=True)

gp = GProfiler(return_dataframe=True)

# === PARAMETERS ===
top_n_terms = 10
min_genes = 20
organism = "hsapiens"  # change to 'mmusculus' for mouse genes

# === LOAD MODULE ASSIGNMENTS ===
print("Loading gene-module assignments...")
df = pd.read_csv(module_file)
df.columns = ["Gene", "Module"]
modules = df["Module"].unique()
summary_rows = []
heatmap_data = {}

# === GO ENRICHMENT LOOP ===
for mod in modules:
    genes = df[df["Module"] == mod]["Gene"].tolist()
    if len(genes) < min_genes:
        continue

    results = gp.profile(organism=organism, query=genes, sources=["GO:BP"])
    if results.empty:
        continue

    # Save individual module enrichment
    results.to_csv(os.path.join(results_dir, f"GO_enrichment_module_{mod}.csv"), index=False)

    # Remove generic terms
    results = results[~results['name'].str.lower().isin(['biological process', 'cellular process'])]

    # Top term summary
    top_term = results.iloc[0]
    summary_rows.append({
        "Module": f"Module_{mod}",
        "Top_GO_Term": top_term['name'],
        "p_value": top_term['p_value'],
        "Term_ID": top_term['native']
    })

    # Build heatmap matrix
    for _, row in results.head(top_n_terms).iterrows():
        val = min(-np.log10(row['p_value']), 50) if row['p_value'] > 0 else np.nan
        heatmap_data.setdefault(row['name'], {})[f"Module_{mod}"] = val

# === SAVE SUMMARY ===
summary_df = pd.DataFrame(summary_rows)
summary_df.to_csv(os.path.join(results_dir, "GO_enrichment_summary.csv"), index=False)

# === PLOT HEATMAP ===
heatmap_df = pd.DataFrame(heatmap_data).T.fillna(0)
heatmap_df.index = heatmap_df.index.str.slice(0, 60)

if not heatmap_df.empty:
    plt.figure(figsize=(max(15, heatmap_df.shape[1] * 1.2), max(8, 0.7 * len(heatmap_df))))
    sns.heatmap(heatmap_df, cmap="viridis", linewidths=0.5, cbar_kws={"label": "-log10(p-value)"})
    plt.title("GO Term Enrichment Across Modules")
    plt.xlabel("Modules")
    plt.ylabel("GO Terms")
    plt.tight_layout()
    plt.savefig(os.path.join(results_dir, "GO_enrichment_heatmap.png"))
    plt.close()

print("✅ GO enrichment complete.")
