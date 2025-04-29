
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Post-WGCNA module analysis:
- GO enrichment
- Trait-eigengene heatmaps
- Summary plots for publication
"""

import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from gprofiler import GProfiler
from scipy.stats import spearmanr
from sklearn.preprocessing import LabelEncoder

# === SETUP ===
base_dir = os.path.dirname(__file__)
results_dir = os.path.join(base_dir, "results_output")
os.makedirs(results_dir, exist_ok=True)

def save_to_results(filename):
    return os.path.join(results_dir, filename)

# === PARAMETERS ===
min_genes = 3
top_n_terms = 10
organism = "hsapiens"  # switch to 'mmusculus' if needed

# === LOAD INPUT FILES ===
#modules = pd.read_csv(save_to_results("gene_module_assignments.csv"))
modules = pd.read_csv(save_to_results("gene_module_assignments.csv"), header=None, names=["Gene", "Module"])

eigengenes = pd.read_csv(save_to_results("module_eigengenes.csv"), index_col=0)
metadata = pd.read_excel("/Users/alex/AlexBadea_MyCodes/Image2Omics/MouseMetaData.xlsx")

metadata = pd.read_excel("/Users/alex/AlexBadea_MyCodes/Eigengenes/MouseMetaData092024AB.xlsx")



# === 1. GO ENRICHMENT ===
gp = GProfiler(return_dataframe=True)
summary_rows = []
heatmap_data = {}

print("Running GO enrichment...")

for mod_id in modules["Module"].unique():
    genes = modules[modules["Module"] == mod_id]["Gene"].tolist()
    if len(genes) < min_genes:
        continue
    result = gp.profile(organism=organism, query=genes, sources=["GO:BP"])
    if result.empty:
        continue

    result.to_csv(save_to_results(f"GO_enrichment_module_{mod_id}.csv"), index=False)
    result = result[~result["name"].str.lower().isin(["biological process", "cellular process"])]
    top = result.iloc[0]
    summary_rows.append({
        "Module": mod_id,
        "Top_GO_Term": top["name"],
        "p_value": top["p_value"]
    })

    for _, row in result.head(top_n_terms).iterrows():
        val = min(-np.log10(row["p_value"]), 50) if row["p_value"] > 0 else np.nan
        heatmap_data.setdefault(row["name"], {})[f"Module_{mod_id}"] = val

summary_df = pd.DataFrame(summary_rows)
summary_df.to_csv(save_to_results("GO_enrichment_summary.csv"), index=False)

# GO heatmap
heatmap_df = pd.DataFrame(heatmap_data).T.fillna(0)
heatmap_df.index = heatmap_df.index.str.slice(0, 60)
plt.figure(figsize=(max(15, heatmap_df.shape[1]*1.2), max(8, 0.7*len(heatmap_df))))
sns.heatmap(heatmap_df, cmap="viridis", linewidths=0.5, cbar_kws={"label": "-log10(p)"})
plt.title("GO Term Enrichment Heatmap")
plt.xlabel("Module")
plt.ylabel("GO Term")
plt.tight_layout()
plt.savefig(save_to_results("GO_enrichment_heatmap.png"))
plt.close()

# === 2. TRAIT-EIGENGENE HEATMAP ===
print("Generating trait-eigengene heatmap...")

# Standardize column names to avoid issues
metadata.columns = metadata.columns.str.strip().str.upper()

# Use 'BADEAID' as the sample ID
if "BADEAID" not in metadata.columns:
    raise ValueError("Column 'BADEAID' not found in metadata.")

metadata["BADEAID"] = metadata["BADEAID"].astype(str).str.strip()
eigengenes.index = eigengenes.index.astype(str).str.strip()

# Align metadata and eigengenes by shared IDs
shared_ids = set(metadata["BADEAID"]) & set(eigengenes.index)
metadata = metadata[metadata["BADEAID"].isin(shared_ids)].set_index("BADEAID")
eigengenes = eigengenes.loc[shared_ids]




traits = metadata[["Sex", "Age_Months", "APOE"]].copy()

# Encode categorical traits
for col in traits.select_dtypes(include="object").columns:
    traits[col] = LabelEncoder().fit_transform(traits[col].astype(str))

# Compute Spearman correlations
cor_mat = pd.DataFrame(index=traits.columns, columns=eigengenes.columns)
pval_mat = pd.DataFrame(index=traits.columns, columns=eigengenes.columns)

for trait in traits.columns:
    for eig in eigengenes.columns:
        rho, pval = spearmanr(traits[trait], eigengenes[eig])
        cor_mat.loc[trait, eig] = rho
        pval_mat.loc[trait, eig] = pval

cor_mat = cor_mat.astype(float)
plt.figure(figsize=(1.2*len(eigengenes.columns), 3.5))
sns.heatmap(cor_mat, cmap="vlag", center=0, annot=True, fmt=".2f", cbar_kws={"label": "Spearman ρ"})
plt.title("Trait–Eigengene Correlation")
plt.ylabel("Trait")
plt.xlabel("Module Eigengene")
plt.tight_layout()
plt.savefig(save_to_results("trait_eigengene_heatmap.png"))
plt.close()

print("✅ Post-analysis complete. Summary figures saved.")
