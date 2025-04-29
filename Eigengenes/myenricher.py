#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GO enrichment per module using g:Profiler
- Input: gene_module_assignments_k15.csv
- Output: enrichment tables + optional plots
"""
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np  # ← add this!
from gprofiler import GProfiler


import os
import pandas as pd
import matplotlib.pyplot as plt
from gprofiler import GProfiler

# === PARAMETERS ===
input_path = "/Users/alex/wgcna_from_expression/gene_module_assignments_k15.csv"
output_dir = os.path.dirname(input_path)  # save in same directory
min_genes = 5  # filter out tiny modules
organism = "hsapiens"

# === LOAD DATA ===
gene_modules = pd.read_csv(input_path)
module_sizes = gene_modules['Module'].value_counts()
valid_modules = module_sizes[module_sizes >= min_genes].index
filtered_modules = gene_modules[gene_modules['Module'].isin(valid_modules)]

# === INIT gProfiler ===
gp = GProfiler(return_dataframe=True)

# === RUN ENRICHMENT ===
for module_id in valid_modules:
    genes = filtered_modules[filtered_modules["Module"] == module_id]["Gene"].tolist()
    results = gp.profile(organism=organism, query=genes, sources=["GO:BP"])
    
    out_csv = os.path.join(output_dir, f"GO_enrichment_module_{module_id}.csv")
    results.to_csv(out_csv, index=False)
    
    # Optional: Save barplot of top GO terms
    top = results.head(10)
    if not top.empty:
        plt.figure(figsize=(10, 5))
        plt.barh(top['name'], -np.log10(top['p_value']))
        plt.xlabel("-log10(p-value)")
        plt.title(f"Top GO terms (Module {module_id})")
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f"GO_enrichment_module_{module_id}.png"))
        plt.close()

# === OPTIONAL: Save histogram of module sizes ===
plt.figure(figsize=(6, 4))
module_sizes[module_sizes >= min_genes].sort_index().plot(kind="bar")
plt.ylabel("Number of Genes")
plt.title("Module Sizes (≥ 5 genes)")
plt.tight_layout()
plt.savefig(os.path.join(output_dir, "module_size_histogram.png"))
plt.close()

print(f"✅ GO enrichment complete. Results saved to {output_dir}")
