
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Dynamic WGCNA-GO pipeline:
- Z-score normalization per gene
- Filters low-variance genes
- Uses dynamic cut height (distance-based)
- Runs GO:BP enrichment
- Saves heatmap and summary
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import squareform
from gprofiler import GProfiler

# === PARAMETERS ===
input_file = "/Users/alex/wgcna_from_expression/RNASEQ_normalized_mouseblood2.xlsx"


# ==== PARAMETERS ====
gene_expression_file = "/Users/alex/AlexBadea_MyCodes/Image2Omics/RNASEQ_normalized_mouseblood2.xlsx"
metadata_file = "/Users/alex/AlexBadea_MyCodes/Image2Omics/MouseMetaData.xlsx"
connectome_folder = "/Users/alex/AlexBadea_MyPapers/StevenWinter/old/MouseBrainNetworks/MouseBrainNetworks_PNAS/old/mouse_connectomes_080823/data/input/mouse_connectomes/"
input_file = gene_expression_file





output_dir = "/Users/alex/wgcna_from_expression"
soft_power = 6
variance_threshold = 0.1
module_cut_height = 0.8
min_genes = 20
top_n_terms = 10
organism = "hsapiens"

# === STEP 1: Load, normalize, and filter ===
print("Loading gene expression...")
df = pd.read_excel(input_file)
gene_names = df.iloc[:, 0].values
expression = df.iloc[:, 1:].T
expression.columns = gene_names

print("Z-score normalizing gene expression...")
expression = (expression - expression.mean()) / expression.std()

print("Filtering low-variance genes...")
variances = expression.var()
expression = expression.loc[:, variances > variance_threshold]
print(f"Remaining genes: {expression.shape[1]}")

# === STEP 2: Build adjacency and cluster ===
print("Building adjacency and performing clustering...")
correlation_matrix = expression.corr(method='pearson')
adjacency_matrix = correlation_matrix.abs() ** soft_power
dissimilarity = 1 - adjacency_matrix
dissimilarity_condensed = squareform(dissimilarity, checks=False)

linkage_matrix = linkage(dissimilarity_condensed, method='average')
module_labels = fcluster(linkage_matrix, t=module_cut_height, criterion='distance')
gene_to_module = pd.Series(module_labels, index=expression.columns)

# Save dendrogram
plt.figure(figsize=(15, 5))
dendrogram(linkage_matrix, no_labels=True)
plt.axhline(y=module_cut_height, color='red', linestyle='--')
plt.title("Gene Clustering Dendrogram")
plt.tight_layout()
plt.savefig(os.path.join(output_dir, "gene_dendrogram_dynamic.png"))
plt.close()

# === STEP 3: Save gene-module assignments ===
module_df = pd.DataFrame({'Gene': expression.columns, 'Module': module_labels})
module_df.to_csv(os.path.join(output_dir, "gene_module_assignments_dynamic.csv"), index=False)

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

    # Remove broad terms
    results = results[~results['name'].str.lower().isin(['biological process', 'cellular process'])]

    # Summary of top GO term
    top_term = results.iloc[0]
    summary_rows.append({
        "Module": f"Module_{mod}",
        "Top_GO_Term": top_term['name'],
        "p_value": top_term['p_value'],
        "Term_ID": top_term['native']
    })

    # Heatmap matrix: top N terms
    for _, row in results.head(top_n_terms).iterrows():
        val = min(-np.log10(row['p_value']), 50) if row['p_value'] > 0 else np.nan
        heatmap_data.setdefault(row['name'], {})[f"Module_{mod}"] = val

# === STEP 5: Save summary + heatmap ===
summary_df = pd.DataFrame(summary_rows)
summary_df.to_csv(os.path.join(output_dir, "GO_enrichment_summary_dynamic.csv"), index=False)

heatmap_df = pd.DataFrame(heatmap_data).T.fillna(0)
heatmap_df.index = heatmap_df.index.str.slice(0, 60)

if not heatmap_df.empty:
    plt.figure(figsize=(max(15, heatmap_df.shape[1] * 1.2), max(8, 0.7 * len(heatmap_df))))
    sns.heatmap(heatmap_df, cmap="viridis", linewidths=0.5, cbar_kws={"label": "-log10(p-value)"})
    plt.title("GO Term Enrichment Across Modules (Dynamic)")
    plt.xlabel("Modules")
    plt.ylabel("GO Terms")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "GO_enrichment_heatmap_dynamic.png"))
    plt.close()

# === Module size histogram ===
module_counts[module_counts >= min_genes].sort_index().plot(kind="bar", figsize=(6, 4))
plt.ylabel("Number of Genes")
plt.title("Module Sizes (≥ 20 genes)")
plt.tight_layout()
plt.savefig(os.path.join(output_dir, "module_size_histogram_dynamic.png"))
plt.close()

print("✅ Dynamic GO enrichment pipeline complete.")
