import os
import json

import anndata
import pandas as pd
import anndata as ad
import sys
import numpy as np


## VIASH START
par = {
    'evaluation_data': 'resources/grn_benchmark/evaluation_data//op.h5ad',
    'models_dir': 'resources/grn_models/op/',
    'models': ['pearson_corr', 'pearson_causal', 'portia', 'ppcor', 'genie3', 'grnboost2', 'scenic', 'scglue', 'celloracle'],
    'regulators_consensus': 'resources/grn_benchmark/prior/regulators_consensus_op.json'
}
## VIASH END
def main(par):
    print(par)
    # Load perturbation data
    adata_rna = anndata.read_h5ad(par['evaluation_data'])
    gene_names = adata_rna.var_names

    gene_dict = {gene_name: i for i, gene_name in enumerate(gene_names)}

    # Load inferred GRNs
    grns = []
    for model in par['models']:
        filepath = os.path.join(par['models_dir'], f'{model}.csv')
        if not os.path.exists(filepath):
            print(f"{filepath} didnt exist. Skipped.")
            continue 
        if model == 'collectri':
            print(f"skipping collectri")
            continue
        A = np.zeros((len(gene_names), len(gene_names)), dtype=float)
        df = pd.read_csv(filepath)

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
    with open(par['regulators_consensus'], 'w') as f:
        json.dump(results, f)
        
    return results


if __name__ == '__main__':
    main(par)
