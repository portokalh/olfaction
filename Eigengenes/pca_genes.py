#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 31 10:23:39 2025

@author: alex
"""
import os
import numpy as np
import pandas as pd
import glob
import igraph as ig
from sklearn.decomposition import PCA
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests

# ---- Step 1: Load Raw Data ----
gene_expression_file = "/Users/alex/AlexBadea_MyCodes/Image2Omics/RNASEQ_normalized_mouseblood2.xlsx"
metadata_file = "/Users/alex/AlexBadea_MyCodes/Image2Omics/MouseMetaData.xlsx"
metadata = pd.read_excel(metadata_file)

connectome_folder = "/Users/alex/AlexBadea_MyPapers/StevenWinter/old/MouseBrainNetworks/MouseBrainNetworks_PNAS/old/mouse_connectomes_080823/data/input/mouse_connectomes/"

results_dir = "results_output"
os.makedirs(results_dir, exist_ok=True)









'''

trying to save loadings


'''

# ---- Step 3: Compute Subject-Specific Eigengenes using PCA ----
gene_expression = pd.read_excel(gene_expression_file)
print(f"Gene expression data loaded. Shape: {gene_expression.shape}")
#gene_expression.rename(columns={gene_expression.columns[0]: "Gene_Name"}, inplace=True)
gene_names = gene_expression.iloc[:, 0].values



gene_expression = gene_expression.iloc[:, 2:]
gene_expression = gene_expression.dropna().T  # Ensure subjects are rows



gene_expression = (gene_expression - gene_expression.mean()) / gene_expression.std()

# Apply PCA
pca = PCA(n_components=min(10, gene_expression.shape[1]))
eigengenes = pca.fit_transform(gene_expression)

# Store PCA results
eigengene_df = pd.DataFrame(eigengenes, index=gene_expression.index, columns=[f"Eigengene_{i+1}" for i in range(10)])
eigengene_df.reset_index(inplace=True)
eigengene_df.rename(columns={'index': 'Badea_ID'}, inplace=True)
eigengene_df.to_csv(os.path.join(results_dir, "eigengenes.csv"), index=False)

# Save PCA Loadings to Identify Gene Contributions
pca_loadings = pd.DataFrame(
    pca.components_.T,  # Genes as rows, eigengenes as columns
    index=gene_names,  # Restore gene names
    columns=[f"Eigengene_{i+1}" for i in range(pca.n_components_)]
)
pca_loadings.to_csv(os.path.join(results_dir, "pca_loadings.csv"))

# Extract top contributing genes for Eigengene_6
eigengene_6_top_genes = pca_loadings["Eigengene_6"].abs().sort_values(ascending=False).head(20)
print("Top contributing genes for Eigengene_6:")
print(eigengene_6_top_genes)

# Save results
eigengene_6_top_genes.to_csv(os.path.join(results_dir, "Eigengene6_Top20Genes.csv"))




'''
end trying to save loadings

'''

# ---- Step 3: Remove Duplicates and Add Badea_ID for Matching ----
metadata = metadata[~metadata["Badea_ID"].duplicated(keep='first')]



# ---- Step 4: Compute Subject-Specific Eigengenes using PCA ----
gene_expression = pd.read_excel(gene_expression_file)
print(f"Gene expression data loaded. Shape: {gene_expression.shape}")
gene_expression.rename(columns={gene_expression.columns[0]: "Badea_ID"}, inplace=True)
gene_expression = gene_expression.T  # Ensure subjects are rows

# Ensure all values are numeric
gene_expression = gene_expression.apply(pd.to_numeric, errors='coerce')

# Drop columns with all NaN values (likely non-numeric columns)
gene_expression = gene_expression.dropna(axis=1, how='all')

# Drop rows with NaN values
gene_expression = gene_expression.dropna()

# Recheck data types to ensure only numeric columns remain
print("Gene expression data types after conversion:\n", gene_expression.dtypes)

# Normalize data
gene_expression = (gene_expression - gene_expression.mean()) / gene_expression.std()


# Convert all values to numeric, forcing errors to NaN
gene_expression = gene_expression.apply(pd.to_numeric, errors='coerce')

# Drop any remaining non-numeric or NaN-only columns
gene_expression = gene_expression.dropna(axis=1, how='all')

# Ensure that no object (string) columns exist
print("Gene expression data types:\n", gene_expression.dtypes)

# Check final shape before PCA
print("Final shape of gene expression before PCA:", gene_expression.shape)

# Run PCA
if gene_expression.shape[1] < 10:
    raise ValueError(f"Not enough valid gene expression features ({gene_expression.shape[1]}) for PCA.")

pca = PCA(n_components=min(10, gene_expression.shape[1]))
eigengenes = pca.fit_transform(gene_expression)



eigengene_df = pd.DataFrame(eigengenes, index=gene_expression.index, columns=[f"Eigengene_{i+1}" for i in range(10)])

eigengene_df = eigengene_df.reset_index().rename(columns={'index': 'Badea_ID'})
if 'Badea_ID' in eigengene_df.columns and eigengene_df.columns.duplicated().any():
    eigengene_df = eigengene_df.loc[:, ~eigengene_df.columns.duplicated()]
print("Eigengene matrix shape:", eigengene_df.shape)  # Should be (num_subjects, 10)

eigengene_df.to_csv(os.path.join(results_dir, "eigengenes.csv"))

print(set(gene_expression.index) & set(metadata.index))  # Subjects present in both



# ---- Step 4: Load and Process Connectomes ----
connectome_files = glob.glob(f"{connectome_folder}/*plain.csv")
print(f"Total connectome files found: {len(connectome_files)}")

graph_metrics = {}
for idx, file in enumerate(connectome_files):
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
    graph_metrics[badea_id] = {"DWI": subject_id, "Avg_Clustering": avg_clustering, "Path_Length": path_length}

graph_metrics_df = pd.DataFrame.from_dict(graph_metrics, orient="index")
graph_metrics_df.index = graph_metrics_df.index.astype(str)
graph_metrics_df.reset_index(inplace=True)
graph_metrics_df.rename(columns={'index': 'Badea_ID'}, inplace=True)
graph_metrics_df.to_csv(os.path.join(results_dir, "connectome_features.csv"), index=False)

print(f"Processing complete! Results saved in {results_dir}/")

# ---- Step 5: Merge Eigengenes with Connectome Features ----
merged_df = pd.merge(eigengene_df, graph_metrics_df, on="Badea_ID", how="inner")
merged_df.to_csv(os.path.join(results_dir, "merged_eigengene_connectome.csv"), index=False)

# ---- Step 5: Merge Eigengenes with Connectome Features ----
merged_df = pd.merge(eigengene_df, graph_metrics_df, on="Badea_ID", how="inner")
merged_df.to_csv(os.path.join(results_dir, "merged_eigengene_connectome.csv"), index=False)


print(f"Merging complete! Results saved in {results_dir}/")

# ---- Step 6: Correlate Eigengenes with Avg_Clustering ----
correlation_results = []
for eigengene in eigengene_df.columns[1:]:  # Skip Badea_ID
    rho, pval = spearmanr(merged_df[eigengene], merged_df["Avg_Clustering"])
    correlation_results.append({"Eigengene": eigengene, "Metric": "Avg_Clustering", "Rho": rho, "P-value": pval})

correlation_df = pd.DataFrame(correlation_results)
correlation_df["FDR-corrected P-value"] = multipletests(correlation_df["P-value"], method='fdr_bh')[1]
correlation_df.to_csv(os.path.join(results_dir, "correlation_avg_clustering.csv"), index=False)

# ---- Step 7: Correlate Eigengenes with Path_Length ----
correlation_results = []
for eigengene in eigengene_df.columns[1:]:  # Skip Badea_ID
    rho, pval = spearmanr(merged_df[eigengene], merged_df["Path_Length"])
    correlation_results.append({"Eigengene": eigengene, "Metric": "Path_Length", "Rho": rho, "P-value": pval})

correlation_df = pd.DataFrame(correlation_results)
correlation_df["FDR-corrected P-value"] = multipletests(correlation_df["P-value"], method='fdr_bh')[1]
correlation_df.to_csv(os.path.join(results_dir, "correlation_path_length.csv"), index=False)

print(f"Analysis complete! Results saved in {results_dir}/")



# ---- Step 5: Correlate Eigengenes with Connectome Features ----
graph_metrics_df = graph_metrics_df.set_index("Badea_ID")
eigengene_df = eigengene_df.set_index("Badea_ID")

correlation_results = []
for eigengene in eigengene_df.columns:
    for metric in ["Avg_Clustering", "Path_Length"]:  # Only these two metrics
        merged_df = pd.concat([
            eigengene_df[[eigengene]],
            graph_metrics_df[[metric]],
            metadata[['Sex', 'Age_Months']]
        ], axis=1, join='inner')
        
        print(f"Processing correlation for {eigengene} and {metric}: {merged_df.shape}")
        
        merged_df = merged_df.dropna()
        if merged_df.empty:
            print(f"Skipping correlation for {eigengene} and {metric} due to insufficient data.")
            continue
        
        rho, pval = spearmanr(merged_df[eigengene], merged_df[metric])
        correlation_results.append({
            "Eigengene": eigengene, 
            "Connectome_Metric": metric, 
            "Rho": rho, 
            "P-value": float(pval)
        })

correlation_df = pd.DataFrame(correlation_results)
correlation_df.to_csv(os.path.join(results_dir, "eigengene_connectome_correlation.csv"), index=False)

print(f"Processing complete! Results saved in {results_dir}/")

# Extract PCA loadings (principal component contributions per gene)
pca_loadings = pd.DataFrame(
    pca.components_.T,  # Transpose so genes are rows, eigengenes are columns
    index=gene_expression.T.index,  # Use gene names as index
    columns=[f"Eigengene_{i+1}" for i in range(pca.n_components_)]
)

# Get top contributing genes for Eigengene_6
eigengene_6_genes = pca_loadings["Eigengene_6"].abs().sort_values(ascending=False).head(20)  # Top 20 genes

# Display the genes contributing most to Eigengene_6
print("Top contributing genes for Eigengene_6:")
print(eigengene_6_genes)

# Save to CSV for further analysis
eigengene_6_genes.to_csv(os.path.join(results_dir, "Eigengene6_TopGenes.csv"))




print(f"Processing complete! Results saved in {results_dir}/")


if "P-value" not in correlation_df.columns:
    print("Error: 'P-value' column is missing. Check correlation_results content:")
    print(correlation_df.head())
else:
    if not correlation_df.empty:
        correlation_df["FDR-corrected P-value"] = multipletests(correlation_df["P-value"], method='fdr_bh')[1]
    else:
        print("Warning: No valid correlations found. FDR correction skipped.")

correlation_df.to_csv(os.path.join(results_dir, "eigengene_connectome_correlation.csv"))
print(f"Analysis complete! Results saved in {results_dir}/")


