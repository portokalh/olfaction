#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Proper weighted coexpression pipeline with progress bars and module dendrogram
Created on Fri Jan 31 10:23:39 2025
@author: alex
"""

# ==== IMPORTS ====
import os
import glob
import numpy as np
import pandas as pd
import igraph as ig
from sklearn.decomposition import PCA
from scipy.stats import spearmanr
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
from scipy.spatial.distance import squareform
from statsmodels.stats.multitest import multipletests
from tqdm import tqdm
import matplotlib.pyplot as plt

# ==== PARAMETERS ====
gene_expression_file = "/Users/alex/AlexBadea_MyCodes/Image2Omics/RNASEQ_normalized_mouseblood2.xlsx"
metadata_file = "/Users/alex/AlexBadea_MyCodes/Image2Omics/MouseMetaData.xlsx"
connectome_folder = "/Users/alex/AlexBadea_MyPapers/StevenWinter/old/MouseBrainNetworks/MouseBrainNetworks_PNAS/old/mouse_connectomes_080823/data/input/mouse_connectomes/"
# === PARAMETERS ===
base_dir = os.path.dirname(__file__)  # directory where your script is
results_dir = os.path.join(base_dir, "results_output")
os.makedirs(results_dir, exist_ok=True)

soft_power = 6
module_cut_height = .92
normalize_genes = True  # ⚠️ Set to False if already normalized (e.g., z-score or log2TPM)

# ==== SETUP ====
os.makedirs(results_dir, exist_ok=True)

# ==== Step 1: Load Gene Expression ====
print("Loading gene expression...")
gene_expression = pd.read_excel(gene_expression_file)
gene_names = gene_expression.iloc[:, 0].values
gene_expression = gene_expression.iloc[:, 2:].dropna().T  # Samples as rows


print(gene_expression.iloc[:, :5].describe())  # look at first 5 genes

# Remove genes expressed in <10% of samples
min_samples = int(0.1 * gene_expression.shape[0])
expression = gene_expression.loc[:, (gene_expression > 0).sum(axis=0) >= min_samples]

# Log transform safely
expression = np.log2(expression + 1)

# Z-score normalize
expression = (expression - expression.mean()) / expression.std()



print("Rechecking normalization...")
print("Mean after normalization:", expression.mean().mean())
print("Std after normalization:", expression.std().mean())



import matplotlib.pyplot as plt
gene_expression.values.flatten()
plt.hist(gene_expression.values.flatten(), bins=100)
plt.title("Distribution of expression values")
plt.xlabel("Expression")
plt.ylabel("Frequency")
plt.tight_layout()
plt.savefig(os.path.join(results_dir, "expression_value_distribution.png"))
plt.close()

# === Inspect normalization status ===
means = expression.mean()
stds = expression.std()

print("Mean across genes (should be ≈0 if z-scored):", means.mean())
print("Std across genes (should be ≈1 if z-scored):", stds.mean())

# Plot a histogram of values
import matplotlib.pyplot as plt
plt.hist(expression.values.flatten(), bins=100)
plt.title("Distribution of expression values")
plt.xlabel("Expression level")
plt.ylabel("Frequency")
plt.tight_layout()
plt.savefig(os.path.join(results_dir, "expression_histogram.png"))
plt.close()




if normalize_genes:
    print("Normalizing gene expression (z-score)...")
    gene_expression = (gene_expression - gene_expression.mean()) / gene_expression.std()

print(f"Final gene expression shape: {gene_expression.shape}")

# ==== Step 2: Load Metadata ====
metadata = pd.read_excel(metadata_file)
metadata = metadata[~metadata["Badea_ID"].duplicated(keep='first')]

# ==== Step 3: Correlation Matrix + Soft-Thresholding ====
print("Calculating gene–gene correlation and adjacency...")
correlation_matrix = gene_expression.corr(method='pearson')
adjacency_matrix = correlation_matrix.abs() ** soft_power

# ==== Step 4: Cluster Genes into Modules ====
print("Clustering genes into modules...")
dissimilarity_condensed = squareform(1 - adjacency_matrix, checks=False)
linkage_matrix = linkage(dissimilarity_condensed, method='average')

# Plot dendrogram
plt.figure(figsize=(12, 6))
dendro = dendrogram(linkage_matrix, no_labels=True, color_threshold=module_cut_height)
plt.axhline(y=module_cut_height, color='red', linestyle='--')
plt.title("Gene Clustering Dendrogram (Colored by Module)")
plt.savefig(os.path.join(results_dir, "dendrogram_modules.png"))
plt.close()

# Assign module labels
module_labels = fcluster(linkage_matrix, t=module_cut_height, criterion='distance')
gene_to_module = pd.Series(module_labels, index=gene_expression.columns)
gene_to_module.to_csv(os.path.join(results_dir, "gene_module_assignments.csv"))

print(f"Identified {len(np.unique(module_labels))} modules.")

# ==== Step 5: Calculate Module Eigengenes ====
print("Calculating module eigengenes...")
eigengene_list = []

for module_id in tqdm(np.unique(module_labels), desc="Modules"):
    module_genes = gene_to_module[gene_to_module == module_id].index.tolist()
    if len(module_genes) < 3:
        continue  # Skip tiny modules
    module_expr = gene_expression[module_genes]
    
    pca = PCA(n_components=1)
    eigengene = pca.fit_transform(module_expr)
    
    eigengene_list.append(pd.DataFrame(
        {f"Module_{module_id}_Eigengene": eigengene.flatten()},
        index=gene_expression.index
    ))

module_eigengenes = pd.concat(eigengene_list, axis=1)
module_eigengenes.index.name = 'Badea_ID'
module_eigengenes.to_csv(os.path.join(results_dir, "module_eigengenes.csv"))
print(f"Saved eigengenes for {len(eigengene_list)} modules.")

# ==== Step 6: Load Connectomes and Compute Features ====
print("Processing connectomes...")
connectome_files = glob.glob(f"{connectome_folder}/*plain.csv")
graph_metrics = {}

for file in tqdm(connectome_files, desc="Connectomes"):
    subject_id = os.path.basename(file)[:6]
    match = metadata.loc[metadata['DWI'].astype(str) == subject_id, 'Badea_ID']
    
    if match.empty:
        continue
    badea_id = match.values[0]
    
    conn_matrix = pd.read_csv(file, index_col=0).apply(pd.to_numeric, errors='coerce').fillna(0).values
    G = ig.Graph.Weighted_Adjacency(conn_matrix.tolist(), mode="undirected", attr="weight")
    
    avg_clustering = G.transitivity_avglocal_undirected()
    path_length = G.average_path_length()
    
    graph_metrics[badea_id] = {
        "DWI": subject_id,
        "Avg_Clustering": avg_clustering,
        "Path_Length": path_length
    }

graph_metrics_df = pd.DataFrame.from_dict(graph_metrics, orient="index")
graph_metrics_df.index = graph_metrics_df.index.astype(str)
graph_metrics_df.reset_index(inplace=True)
graph_metrics_df.rename(columns={'index': 'Badea_ID'}, inplace=True)
graph_metrics_df.to_csv(os.path.join(results_dir, "connectome_features.csv"), index=False)

# ==== Step 7: Merge Data ====
print("Merging eigengenes and connectome features...")
merged_df = pd.merge(module_eigengenes, graph_metrics_df, on="Badea_ID", how="inner")
merged_df.to_csv(os.path.join(results_dir, "merged_module_eigengene_connectome.csv"), index=False)

# ==== Step 8: Correlate Eigengenes with Connectome Metrics ====
print("Correlating eigengenes with connectome metrics...")
correlation_results = []

for eigengene in tqdm(module_eigengenes.columns, desc="Eigengene correlations"):
    for metric in ["Avg_Clustering", "Path_Length"]:
        temp_df = pd.concat([
            module_eigengenes[[eigengene]],
            graph_metrics_df.set_index('Badea_ID')[[metric]],
            metadata.set_index('Badea_ID')[['Sex', 'Age_Months']]
        ], axis=1, join='inner')
        
        temp_df = temp_df.dropna()
        if temp_df.empty:
            continue
        
        rho, pval = spearmanr(temp_df[eigengene], temp_df[metric])
        correlation_results.append({
            "Eigengene": eigengene,
            "Connectome_Metric": metric,
            "Spearman_rho": rho,
            "P-value": pval
        })

correlation_df = pd.DataFrame(correlation_results)

if not correlation_df.empty:
    correlation_df["FDR-corrected P-value"] = multipletests(correlation_df["P-value"], method='fdr_bh')[1]

correlation_df.to_csv(os.path.join(results_dir, "eigengene_connectome_correlations.csv"), index=False)

print("✅ All steps complete! Results saved to:", results_dir)
