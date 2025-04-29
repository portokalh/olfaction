
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Refined WGCNA-style GO enrichment pipeline with stricter parameters:
- Soft power = 10
- Variance threshold = 0.05
- Module cut height = 1.2
- Minimum genes per module = 20
- GO enrichment for human (hsapiens)
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import squareform
from gprofiler import GProfiler

# === PARAMETERS ===
input_file = "/Users/alex/wgcna_from_expression/RNASEQ_normalized_mouseblood2.xlsx"
output_dir = "/Users/alex/wgcna_from_expression"


# === PARAMETERS ===
input_file = "/Users/alex/wgcna_from_expression/RNASEQ_normalized_mouseblood2.xlsx"
output_dir = "/Users/alex/wgcna_from_expression"

# ==== PARAMETERS ====
gene_expression_file = "/Users/alex/AlexBadea_MyCodes/Image2Omics/RNASEQ_normalized_mouseblood2.xlsx"
metadata_file = "/Users/alex/AlexBadea_MyCodes/Image2Omics/MouseMetaData.xlsx"
connectome_folder = "/Users/alex/AlexBadea_MyPapers/StevenWinter/old/MouseBrainNetworks/MouseBrainNetworks_PNAS/old/mouse_connectomes_080823/data/input/mouse_connectomes/"
input_file = gene_expression_file


soft_power = 10
variance_threshold = 0.05
module_cut_height = 1.2
min_genes = 20
organism = "hsapiens"

# === LOAD EXPRESSION DATA ===
print("Loading gene expression...")
df = pd.read_excel(input_file)
gene_names = df.iloc[:, 0].values
expression = df.iloc[:, 1:].T
expression.columns = gene_names

# === FILTER LOW-VARIANCE GENES ===
print("Filtering low-variance genes...")
variances = expression.var()
expression = expression.loc[:, variances > variance_threshold]
print(f"Remaining genes after filtering: {expression.shape[1]}")

# === BUILD ADJACENCY AND DISSIMILARITY ===
print("Computing adjacency matrix...")
correlation_matrix = expression.corr(method='pearson')
adjacency_matrix = correlation_matrix.abs() ** soft_power
dissimilarity = 1 - adjacency_matrix
dissimilarity_condensed = squareform(dissimilarity, checks=False)

# === CLUSTER GENES ===
print("Clustering genes with dynamic height cut...")
linkage_matrix = linkage(dissimilarity_condensed, method='average')
module_labels = fcluster(linkage_matrix, t=module_cut_height, criterion='distance')
gene_to_module = pd.Series(module_labels, index=expression.columns)

# === SAVE MODULE ASSIGNMENTS ===
module_df = pd.DataFrame({'Gene': expression.columns, 'Module': module_labels})
module_df.to_csv(os.path.join(output_dir, "gene_module_assignments_dynamic_strict.csv"), index=False)

# === RUN GO ENRICHMENT ===
print("Running GO enrichment for stable modules...")
gp = GProfiler(return_dataframe=True)
module_counts = module_df['Module'].value_counts()
valid_modules = module_counts[module_counts >= min_genes].index
filtered = module_df[module_df["Module"].isin(valid_modules)]

for mod in valid_modules:
    genes = filtered[filtered["Module"] == mod]["Gene"].tolist()
    results = gp.profile(organism=organism, query=genes, sources=["GO:BP"])
    results.to_csv(os.path.join(output_dir, f"GO_enrichment_module_{mod}.csv"), index=False)

    top = results.head(10)
    if not top.empty:
        plt.figure(figsize=(10, 5))
        plt.barh(top['name'], -np.log10(top['p_value']))
        plt.xlabel("-log10(p-value)")
        plt.title(f"Top GO terms (Module {mod})")
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f"GO_enrichment_module_{mod}.png"))
        plt.close()

# === PLOT MODULE SIZE HISTOGRAM ===
plt.figure(figsize=(6, 4))
module_counts[module_counts >= min_genes].sort_index().plot(kind="bar")
plt.ylabel("Number of Genes")
plt.title("Module Sizes (≥ 20 genes)")
plt.tight_layout()
plt.savefig(os.path.join(output_dir, "module_size_histogram_dynamic_strict.png"))
plt.close()

print("✅ GO enrichment complete.")
