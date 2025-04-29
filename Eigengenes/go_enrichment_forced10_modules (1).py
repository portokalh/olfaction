
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
WGCNA + GO Enrichment Pipeline (Final Version):
- Forces exactly 10 gene modules
- Filters low-variance genes
- Runs GO:BP enrichment
- Saves summary + improved heatmap
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from gprofiler import GProfiler

# === PARAMETERS ===
input_file = "/Users/alex/wgcna_from_expression/RNASEQ_normalized_mouseblood2.xlsx"
output_dir = "/Users/alex/wgcna_from_expression"
soft_power = 10
variance_threshold = 0.05
desired_module_count = 10
min_genes = 20
top_n_terms = 10
organism = "hsapiens"

# === STEP 1: Load and filter expression ===
print("Loading gene expression...")
df = pd.read_excel(input_file)
gene_names = df.iloc[:, 0].values
expression = df.iloc[:, 1:].T
expression.columns = gene_names

print("Filtering low-variance genes...")
variances = expression.var()
expression = expression.loc[:, variances > variance_threshold]
print(f"Remaining genes: {expression.shape[1]}")

# === STEP 2: Correlation, adjacency, clustering ===
print("Clustering genes into exactly 10 modules...")
correlation_matrix = expression.corr(method='pearson')
adjacency_matrix = correlation_matrix.abs() ** soft_power
dissimilarity = 1 - adjacency_matrix
dissimilarity_condensed = squareform(dissimilarity, checks=False)

linkage_matrix = linkage(dissimilarity_condensed, method='average')
module_labels = fcluster(linkage_matrix, t=desired_module_count, criterion='maxclust')
gene_to_module = pd.Series(module_labels, index=expression.columns)

# === STEP 3: Save gene-module assignments ===
module_df = pd.DataFrame({'Gene': expression.columns, 'Module': module_labels})
module_df.to_csv(os.path.join(output_dir, "gene_module_assignments_forced10.csv"), index=False)

# === STEP 4: GO enrichment per module ===
print("Running GO enrichment...")
gp = GProfiler(return_dataframe=True)
module_counts = module_df['Module'].value_counts()
valid_modules = module_counts[module_counts >= min_genes].index
filtered = module_df[module_df["Module"].isin(valid_modules)]

summary_rows = []
heatmap_data = {}

for mod in valid_modules:
    genes = filtered[filtered["Module"] == mod]["Gene"].tolist()
    results = gp.profile(organism=organism, query=genes, sources=["GO:BP"])
    if results.empty:
        continue

    # Save enrichment CSV
    csv_path = os.path.join(output_dir, f"GO_enrichment_module_{mod}.csv")
    results.to_csv(csv_path, index=False)

    # Filter out broad/ambiguous terms
    results = results[~results['name'].str.lower().isin(['biological process', 'cellular process'])]

    # Add top term to summary
    top_term = results.iloc[0]
    summary_rows.append({
        "Module": f"Module_{mod}",
        "Top_GO_Term": top_term['name'],
        "p_value": top_term['p_value'],
        "Term_ID": top_term['native']
    })

    # Build heatmap matrix (capped -log10(p))
    for _, row in results.head(top_n_terms).iterrows():
        val = min(-np.log10(row['p_value']), 50) if row['p_value'] > 0 else np.nan
        heatmap_data.setdefault(row['name'], {})[f"Module_{mod}"] = val

# === STEP 5: Save summary + improved heatmap ===
summary_df = pd.DataFrame(summary_rows)
summary_df.to_csv(os.path.join(output_dir, "GO_enrichment_summary.csv"), index=False)

heatmap_df = pd.DataFrame(heatmap_data).T.fillna(0)
heatmap_df.index = heatmap_df.index.str.slice(0, 60)  # truncate long GO terms

if not heatmap_df.empty:
    plt.figure(figsize=(max(15, heatmap_df.shape[1] * 1.2), max(8, 0.7 * len(heatmap_df))))
    sns.heatmap(heatmap_df, cmap="viridis", linewidths=0.5, cbar_kws={"label": "-log10(p-value)"})
    plt.title("GO Term Enrichment Across Modules")
    plt.xlabel("Modules")
    plt.ylabel("GO Terms")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "GO_enrichment_heatmap.png"))
    plt.close()

# === Module size histogram ===
module_counts[module_counts >= min_genes].sort_index().plot(kind="bar", figsize=(6, 4))
plt.ylabel("Number of Genes")
plt.title("Module Sizes (≥ 20 genes)")
plt.tight_layout()
plt.savefig(os.path.join(output_dir, "module_size_histogram_forced10.png"))
plt.close()

print("✅ GO enrichment pipeline complete.")
