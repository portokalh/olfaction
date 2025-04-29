# ==== IMPORTS ====
import os
import glob
import numpy as np
import pandas as pd
import igraph as ig
from sklearn.decomposition import PCA
from scipy.stats import spearmanr
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from statsmodels.stats.multitest import multipletests


import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
from scipy.spatial.distance import squareform
from scipy.stats import spearmanr
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
from scipy.spatial.distance import squareform

# ---- Step 1: Load Raw Data ----
gene_expression_file = "/Users/alex/AlexBadea_MyCodes/Image2Omics/RNASEQ_normalized_mouseblood2.xlsx"
metadata_file = "/Users/alex/AlexBadea_MyCodes/Image2Omics/MouseMetaData.xlsx"
metadata = pd.read_excel(metadata_file)

connectome_folder = "/Users/alex/AlexBadea_MyPapers/StevenWinter/old/MouseBrainNetworks/MouseBrainNetworks_PNAS/old/mouse_connectomes_080823/data/input/mouse_connectomes/"

results_dir = "results_output"
os.makedirs(results_dir, exist_ok=True)


#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 31 10:23:39 2025

@author: alex
"""

import os
import glob
import numpy as np
import pandas as pd
import igraph as ig
from sklearn.decomposition import PCA
from scipy.stats import spearmanr
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from statsmodels.stats.multitest import multipletests

# ---- PARAMETERS ----
gene_expression_file = "/Users/alex/AlexBadea_MyCodes/Image2Omics/RNASEQ_normalized_mouseblood2.xlsx"
metadata_file = "/Users/alex/AlexBadea_MyCodes/Image2Omics/MouseMetaData.xlsx"
connectome_folder = "/Users/alex/AlexBadea_MyPapers/StevenWinter/old/MouseBrainNetworks/MouseBrainNetworks_PNAS/old/mouse_connectomes_080823/data/input/mouse_connectomes/"
results_dir = "results_output"
soft_power = 6  # Set soft thresholding power

os.makedirs(results_dir, exist_ok=True)

# ---- Step 1: Load and Normalize Gene Expression Data ----
print("Loading gene expression data...")
gene_expression = pd.read_excel(gene_expression_file)
gene_names = gene_expression.iloc[:, 0].values
gene_expression = gene_expression.iloc[:, 2:].dropna().T  # Samples as rows

# Standardize each gene
gene_expression = (gene_expression - gene_expression.mean()) / gene_expression.std()

print(f"Normalized gene expression data. Shape: {gene_expression.shape}")

# ---- Step 2: Load Metadata ----
print("Loading metadata...")
metadata = pd.read_excel(metadata_file)
metadata = metadata[~metadata["Badea_ID"].duplicated(keep='first')]

# ---- Step 3: Build Gene–Gene Correlation Matrix ----
print("Computing gene–gene correlation matrix...")
correlation_matrix = gene_expression.corr(method='pearson')

# ---- Step 4: Soft-Thresholding to Get Adjacency Matrix ----
adjacency_matrix = correlation_matrix.abs() ** soft_power

# ---- Step 5: Cluster Genes into Modules ----
print("Clustering genes into modules...")
dissimilarity = 1 - adjacency_matrix
#linkage_matrix = linkage(dissimilarity, method='average')
module_labels = fcluster(linkage_matrix, t=0.9, criterion='distance')  # Adjust 't' for tighter or looser modules


dissimilarity_condensed = squareform(1 - adjacency_matrix, checks=False)
linkage_matrix = linkage(dissimilarity_condensed, method='average')


gene_to_module = pd.Series(module_labels, index=gene_expression.columns)

# Save module assignment
gene_to_module.to_csv(os.path.join(results_dir, "gene_module_assignments.csv"))

# ---- Step 6: Calculate Module Eigengenes ----
print("Calculating module eigengenes...")
eigengene_list = []
module_ids = np.unique(module_labels)

for module_id in module_ids:
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

print(f"Eigengenes calculated for {len(module_ids)} modules. Shape: {module_eigengenes.shape}")

# ---- Step 7: Process Connectomes ----
print("Loading connectomes...")
connectome_files = glob.glob(f"{connectome_folder}/*plain.csv")
graph_metrics = {}

for file in connectome_files:
    subject_id = file.split("/")[-1][:6]  # Extract first 6 characters
    badea_id_match = metadata.loc[metadata['DWI'].astype(str) == subject_id, 'Badea_ID']
    
    if badea_id_match.empty:
        print(f"Warning: No matching Badea_ID found for {subject_id}")
        continue
    
    badea_id = badea_id_match.values[0]
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

# ---- Step 8: Merge Eigengenes with Connectome Features ----
print("Merging eigengenes and connectome features...")
merged_df = pd.merge(module_eigengenes, graph_metrics_df, on="Badea_ID", how="inner")
merged_df.to_csv(os.path.join(results_dir, "merged_module_eigengene_connectome.csv"), index=False)

# ---- Step 9: Correlate Eigengenes with Connectome Metrics ----
print("Running correlation analysis...")
correlation_results = []

for eigengene in module_eigengenes.columns:
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

print("✅ Analysis complete! All results saved.")

'''
File    Purpose
gene_module_assignments.csv Which genes belong to which module
module_eigengenes.csv   Subject-specific eigengene matrix
connectome_features.csv Connectome-derived metrics
merged_module_eigengene_connectome.csv  Combined file of eigengenes + connectome metrics
eigengene_connectome_correlations.csv   Correlations between module eigengenes and brain connectome metrics
'''