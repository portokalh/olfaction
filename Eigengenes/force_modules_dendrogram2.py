
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Expanded WGCNA-style clustering:
- Forces number of modules
- Saves gene-to-module assignments
- Calculates module eigengenes (sample x module matrix)
- Saves gene loadings (gene contribution to eigengenes)
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import squareform

# === PARAMETERS ===
input_file = "module_eigengenes.csv"  # Input gene expression matrix
output_dir = "dendrogram_trials_expanded"
soft_power = 6
forced_module_count = 15

# === SETUP ===
os.makedirs(output_dir, exist_ok=True)

# === STEP 1: Load eigengenes and compute dissimilarity ===
print("Loading gene expression...")
gene_expression = pd.read_csv(input_file, index_col=0)

print("Computing adjacency and dissimilarity matrices...")
correlation_matrix = gene_expression.corr(method='pearson')
adjacency_matrix = correlation_matrix.abs() ** soft_power
dissimilarity_condensed = squareform(1 - adjacency_matrix, checks=False)

# === STEP 2: Hierarchical clustering ===
print("Performing hierarchical clustering...")
linkage_matrix = linkage(dissimilarity_condensed, method='average')

# === STEP 3: Assign modules ===
print(f"Forcing {forced_module_count} modules...")
module_labels = fcluster(linkage_matrix, t=forced_module_count, criterion='maxclust')
gene_to_module = pd.Series(module_labels, index=gene_expression.columns)

# Save gene-to-module assignments
gene_to_module_df = pd.DataFrame({'Gene': gene_expression.columns, 'Module': module_labels})
gene_to_module_df.to_csv(os.path.join(output_dir, f"gene_module_assignments_k{forced_module_count}.csv"), index=False)

# === STEP 4: Calculate eigengenes and gene loadings ===
print("Calculating module eigengenes and gene loadings...")
eigengene_list = []
loading_list = []

for module_id in np.unique(module_labels):
    module_genes = gene_to_module[gene_to_module == module_id].index.tolist()
    if len(module_genes) < 3:
        continue  # Skip very small modules
    
    module_expr = gene_expression[module_genes]
    
    pca = PCA(n_components=1)
    eigengene_values = pca.fit_transform(module_expr)
    
    # Save sample eigengene values
    eigengene_list.append(pd.DataFrame(
        {f"Module_{module_id}_Eigengene": eigengene_values.flatten()},
        index=gene_expression.index
    ))
    
    # Save gene loadings
    loadings = pd.DataFrame(
        {f"Module_{module_id}_Loading": pca.components_.flatten()},
        index=module_genes
    )
    loading_list.append(loadings)

# Combine all eigengenes and loadings
module_eigengenes = pd.concat(eigengene_list, axis=1)
module_loadings = pd.concat(loading_list, axis=0)

# Save outputs
module_eigengenes.index.name = 'Sample_ID'
module_eigengenes.to_csv(os.path.join(output_dir, f"module_eigengenes_k{forced_module_count}.csv"))
module_loadings.index.name = 'Gene'
module_loadings.to_csv(os.path.join(output_dir, f"module_gene_loadings_k{forced_module_count}.csv"))

print("✅ Done. Files saved in", output_dir)
