
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Final WGCNA-style weighted coexpression pipeline
with proper normalization, dynamic clustering, and saving.
Created on 2025-01-31
Author: Alex
"""

# ==== IMPORTS ====
import os
import glob
import numpy as np
import pandas as pd
import igraph as ig
import seaborn as sns
from sklearn.decomposition import PCA
from scipy.stats import spearmanr
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
from scipy.spatial.distance import squareform
from statsmodels.stats.multitest import multipletests
from tqdm import tqdm
import matplotlib.pyplot as plt
from matplotlib.colors import to_hex
from scipy.cluster.hierarchy import fcluster

# ==== SETUP ====
base_dir = os.path.dirname(__file__)
results_dir = os.path.join(base_dir, "results_output")
os.makedirs(results_dir, exist_ok=True)

def save_to_results(filename):
    return os.path.join(results_dir, filename)

# ==== PARAMETERS ====
gene_expression_file = "/Users/alex/AlexBadea_MyCodes/Image2Omics/RNASEQ_normalized_mouseblood2.xlsx"
metadata_file = "/Users/alex/AlexBadea_MyCodes/Image2Omics/MouseMetaData.xlsx"
connectome_folder = "/Users/alex/AlexBadea_MyPapers/StevenWinter/old/MouseBrainNetworks/MouseBrainNetworks_PNAS/old/mouse_connectomes_080823/data/input/mouse_connectomes/"
soft_power = 4
module_cut_height = 0.5

# ==== Step 1: Load and Normalize Gene Expression ====
print("Loading and normalizing gene expression...")
gene_expression = pd.read_excel(gene_expression_file)
gene_names = gene_expression.iloc[:, 0].values
expression = gene_expression.iloc[:, 2:].T
expression.columns = gene_names

# Filter genes expressed in at least 10% of samples
min_samples = int(0.1 * expression.shape[0])
expression = expression.loc[:, (expression > 0).sum(axis=0) >= min_samples]

# Apply log2(x+1) normalization
expression = np.log2(expression + 1)

# Apply Z-score normalization
expression = (expression - expression.mean()) / expression.std()

# Confirm normalization
print("Mean across genes (should be ≈0):", expression.mean().mean())
print("Std across genes (should be ≈1):", expression.std().mean())

# Save distribution plot
plt.hist(expression.values.flatten(), bins=100)
plt.title("Distribution of Normalized Expression Values")
plt.xlabel("Expression Level")
plt.ylabel("Frequency")
plt.tight_layout()
plt.savefig(save_to_results("expression_histogram.png"))
plt.close()

# ==== Step 2: Load Metadata ====
metadata = pd.read_excel(metadata_file)
metadata = metadata[~metadata["Badea_ID"].duplicated(keep='first')]

# ==== Step 3: Build Correlation and Cluster ====
print("Calculating gene–gene correlation and adjacency...")
correlation_matrix = expression.corr(method='pearson')
adjacency_matrix = correlation_matrix.abs() ** soft_power

print("Clustering genes into modules...")
dissimilarity_condensed = squareform(1 - adjacency_matrix, checks=False)
linkage_matrix = linkage(dissimilarity_condensed, method='average')

# Assign module labels dynamically
#module_labels = fcluster(linkage_matrix, t=module_cut_height, criterion='distance')
module_labels = fcluster(linkage_matrix, t=10, criterion='maxclust')

gene_to_module = pd.Series(module_labels, index=expression.columns)
gene_to_module.to_csv(save_to_results("gene_module_assignments.csv"))

print(f"Identified {len(np.unique(module_labels))} modules.")

# Plot colored dendrogram
print("Plotting colored dendrogram...")
module_series = pd.Series(module_labels, index=expression.columns)
unique_modules = np.unique(module_labels)
color_palette = sns.color_palette("husl", len(unique_modules))
module_color_map = dict(zip(unique_modules, map(to_hex, color_palette)))
leaf_colors = module_series.map(module_color_map)

def leaf_color_func(leaf_id):
    gene_name = expression.columns[leaf_id]
    return leaf_colors[gene_name]

plt.figure(figsize=(15, 6))
dendrogram(
    linkage_matrix,
    no_labels=True,
    leaf_rotation=90,
    link_color_func=lambda k: 'black',
    leaf_font_size=4,
    color_threshold=0,
    above_threshold_color='black',
    get_leaves=True,
    #leaf_color_func=leaf_color_func
)
plt.title("Colored Dendrogram by Module")
plt.tight_layout()
plt.savefig(save_to_results("dendrogram_colored_by_module.png"))
plt.close()

# ==== Step 4: Calculate Module Eigengenes ====
print("Calculating module eigengenes...")
eigengene_list = []

for module_id in tqdm(np.unique(module_labels), desc="Modules"):
    module_genes = gene_to_module[gene_to_module == module_id].index.tolist()
    if len(module_genes) < 3:
        continue
    module_expr = expression[module_genes]
    pca = PCA(n_components=1)
    eigengene = pca.fit_transform(module_expr)
    eigengene_list.append(pd.DataFrame(
        {f"Module_{module_id}_Eigengene": eigengene.flatten()},
        index=expression.index
    ))

module_eigengenes = pd.concat(eigengene_list, axis=1)
module_eigengenes.index.name = 'Badea_ID'
module_eigengenes.to_csv(save_to_results("module_eigengenes.csv"))
print(f"Saved eigengenes for {len(eigengene_list)} modules.")

# ==== Step 5: Connectome Features ====
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
    graph_metrics[badea_id] = {
        "DWI": subject_id,
        "Avg_Clustering": G.transitivity_avglocal_undirected(),
        "Path_Length": G.average_path_length()
    }

graph_metrics_df = pd.DataFrame.from_dict(graph_metrics, orient="index")
graph_metrics_df.index = graph_metrics_df.index.astype(str)
graph_metrics_df.reset_index(inplace=True)
graph_metrics_df.rename(columns={'index': 'Badea_ID'}, inplace=True)
graph_metrics_df.to_csv(save_to_results("connectome_features.csv"), index=False)

# ==== Step 6: Merge and Correlate Eigengenes ====
print("Merging eigengenes and connectome features...")
merged_df = pd.merge(module_eigengenes, graph_metrics_df, on="Badea_ID", how="inner")
merged_df.to_csv(save_to_results("merged_module_eigengene_connectome.csv"), index=False)

print("Correlating eigengenes with connectome metrics...")
correlation_results = []

for eigengene in tqdm(module_eigengenes.columns, desc="Eigengene correlations"):
    for metric in ["Avg_Clustering", "Path_Length"]:
        temp_df = pd.concat([
            module_eigengenes[[eigengene]],
            graph_metrics_df.set_index('Badea_ID')[[metric]],
            metadata.set_index('Badea_ID')[['Sex', 'Age_Months']]
        ], axis=1, join='inner').dropna()
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
correlation_df.to_csv(save_to_results("eigengene_connectome_correlations.csv"), index=False)

print("✅ All steps complete! Results saved to:", results_dir)
