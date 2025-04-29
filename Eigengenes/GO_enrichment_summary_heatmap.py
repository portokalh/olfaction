
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Summarize GO enrichment results:
- Extract top GO terms per module
- Create a summary CSV
- Create a heatmap of -log10(p-values) across top GO terms
"""

import os
import glob
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# === PARAMETERS ===
input_dir = "/Users/alex/wgcna_from_expression"
summary_csv = os.path.join(input_dir, "GO_enrichment_summary.csv")
heatmap_png = os.path.join(input_dir, "GO_enrichment_heatmap.png")
top_n_terms = 10

# === LOAD ALL ENRICHMENT FILES ===
files = glob.glob(os.path.join(input_dir, "GO_enrichment_module_*.csv"))
summary_rows = []
heatmap_data = {}

for file in files:
    module_id = os.path.basename(file).split("_")[-1].split(".")[0]
    df = pd.read_csv(file)
    if df.empty:
        continue

    # Top GO term for summary
    top_term = df.iloc[0]
    summary_rows.append({
        "Module": module_id,
        "Top_GO_Term": top_term['name'],
        "p_value": top_term['p_value'],
        "Term_ID": top_term['native']
    })

    # Top N GO terms for heatmap
    top_terms = df.head(top_n_terms)
    for _, row in top_terms.iterrows():
        term = row['name']
        val = -np.log10(row['p_value']) if row['p_value'] > 0 else np.nan
        heatmap_data.setdefault(term, {})[f"Module_{module_id}"] = val

# === SAVE SUMMARY CSV ===
summary_df = pd.DataFrame(summary_rows)
summary_df.to_csv(summary_csv, index=False)

# === GENERATE HEATMAP ===
heatmap_df = pd.DataFrame(heatmap_data).T.fillna(0)

if not heatmap_df.empty:
    plt.figure(figsize=(min(20, heatmap_df.shape[1] * 0.8), 0.5 * len(heatmap_df)))
    sns.heatmap(heatmap_df, cmap="viridis", annot=False, linewidths=0.5)
    plt.title("GO Term Enrichment (-log10 p-values)")
    plt.xlabel("Modules")
    plt.ylabel("GO Terms")
    plt.tight_layout()
    plt.savefig(heatmap_png)
    plt.close()

print(f"✅ Summary saved to: {summary_csv}")
print(f"✅ Heatmap saved to: {heatmap_png}")
