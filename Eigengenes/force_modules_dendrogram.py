
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generate gene dendrograms for various forced module counts using hierarchical clustering.
Includes fcluster(maxclust=...) to directly control the number of modules.
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import squareform

# === PARAMETERS ===
input_file = "module_eigengenes.csv"  # eigengene matrix: samples x Module_X_Eigengene
output_dir = "dendrogram_trials"
soft_power = 6
module_counts_to_try = [10, 12, 15, 18, 20]

# === SETUP ===
os.makedirs(output_dir, exist_ok=True)

# === STEP 1: Load eigengenes and compute dissimilarity ===
print("Loading eigengenes...")
gene_expression = pd.read_csv(input_file, index_col=0)

print("Computing adjacency matrix with soft thresholding...")
correlation_matrix = gene_expression.corr(method='pearson')
adjacency_matrix = correlation_matrix.abs() ** soft_power
dissimilarity_condensed = squareform(1 - adjacency_matrix, checks=False)

# === STEP 2: Hierarchical clustering ===
print("Performing hierarchical clustering...")
linkage_matrix = linkage(dissimilarity_condensed, method='average')

# === STEP 3: Plot dendrograms for various forced module counts ===
print("Generating dendrograms...")
fig, axes = plt.subplots(len(module_counts_to_try), 1, figsize=(12, 3 * len(module_counts_to_try)))

for ax, k in zip(axes, module_counts_to_try):
    labels = fcluster(linkage_matrix, t=k, criterion='maxclust')
    dendrogram(linkage_matrix, no_labels=True, ax=ax, color_threshold=None)
    ax.set_title(f"{k} Modules Forced (maxclust)")

plt.tight_layout()
plt.savefig(os.path.join(output_dir, "dendrograms_by_module_count.png"))
plt.close()
print("Saved dendrogram figure.")

# === STEP 4: Save module labels for k = 15 ===
k_final = 15
print(f"Saving module assignments for {k_final} modules...")
labels_k15 = fcluster(linkage_matrix, t=k_final, criterion='maxclust')
labels_df = pd.DataFrame({'Gene': gene_expression.columns, 'Module': labels_k15})
labels_df.to_csv(os.path.join(output_dir, f"gene_module_assignments_k{k_final}.csv"), index=False)

print("✅ Done.")
