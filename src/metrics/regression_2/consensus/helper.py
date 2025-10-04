
import os
import json
import anndata
import pandas as pd
import anndata as ad
import sys
import numpy as np

from util import process_links
def main(par):
    print(par)
    # Load perturbation data
    adata_rna = anndata.read_h5ad(par['evaluation_data'])
    gene_names = adata_rna.var_names
    
    gene_dict = {gene_name: i for i, gene_name in enumerate(gene_names)}

    # Load inferred GRNs
    grns = []
    for filepath in par['predictions']:
        net = ad.read_h5ad(filepath)
        net = net.uns['prediction']

        A = np.zeros((len(gene_names), len(gene_names)), dtype=float)
        net = process_links(net, par)
        for source, target, weight in zip(net['source'], net['target'], net['weight']):
            if (source not in gene_dict) or (target not in gene_dict):
                continue
            i = gene_dict[source]
            j = gene_dict[target]
            A[i, j] = float(weight)
        print(f'Sparsity of {filepath}: {np.mean(A == 0)}')
        grns.append(A)
    grns = np.asarray(grns)

    assert len(grns) > 0, "No GRNs were loaded. Check the models directory and the models specified."

    # Compute consensus number of putative regulators for each gene (and for each theta value)
    thetas = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
    n_tfs = {}
    for theta in thetas:
        n_tfs[theta] = np.round(np.quantile(np.sum(grns != 0, axis=1), theta, axis=0)).astype(int)

    # Store results
    results = {}
    for i, gene_name in enumerate(gene_names):
        results[gene_name] = {theta: int(n_tfs[theta][i]) for theta in thetas}
    with open(par['regulators_consensus'], 'w') as f:
        json.dump(results, f)
        
    return results
