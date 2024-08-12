import os
import json

import anndata
import pandas as pd
import anndata as ad
import sys
import numpy as np


## VIASH START
par = {
    'perturbation_data': 'resources/grn-benchmark/perturbation_data.h5ad',
    'grn_folder': 'resources/grn-benchmark/grn_models',
    'grns': 'figr.csv,scenicplus.csv',
    'output': 'resources/grn-benchmark/consensus-num-regulators.json'
}
## VIASH END


# Load perturbation data
adata_rna = anndata.read_h5ad(par['perturbation_data'])
try:
    gene_names = adata_rna.var.gene.to_numpy()
except:
    gene_names = adata_rna.var.index.to_numpy()

print(par['perturbation_data'])
print(par['grns'])

# Load inferred GRNs
grns = []
for filename in par['grns'].split(','):
    filepath = os.path.join(par['grn_folder'], filename)
    gene_dict = {gene_name: i for i, gene_name in enumerate(gene_names)}
    A = np.zeros((len(gene_names), len(gene_names)), dtype=float)
    df = pd.read_csv(filepath, sep=',', header='infer', index_col=0)
    for source, target, weight in zip(df['source'], df['target'], df['weight']):
        if (source not in gene_dict) or (target not in gene_dict):
            continue
        i = gene_dict[source]
        j = gene_dict[target]
        A[i, j] = float(weight)
    print(f'Sparsity of {filepath}: {np.mean(A == 0)}')
    grns.append(A)
grns = np.asarray(grns)

# Compute consensus number of putative regulators for each gene (and for each theta value)
thetas = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
n_tfs = {}
for theta in thetas:
    n_tfs[theta] = np.round(np.quantile(np.sum(grns != 0, axis=1), theta, axis=0)).astype(int)

# Store results
results = {}
for i, gene_name in enumerate(gene_names):
    results[gene_name] = {theta: int(n_tfs[theta][i]) for theta in thetas}
with open(par['output'], 'w') as f:
    json.dump(results, f)
