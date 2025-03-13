import os
import json

import anndata
import pandas as pd
import anndata as ad
import sys
import numpy as np


def main(dataset):
    par = {
        'evaluation_data': f'resources/grn_benchmark/evaluation_data//{dataset}_bulk.h5ad',
        'models_dir': f'resources/grn_models/{dataset}/',
        'models': ['pearson_corr', 'portia', 'ppcor', 'grnboost2', 'scenic', 'scprint', 'figr', 'celloracle', 'scenic+', 'granie', 'scglue'],
        'regulators_consensus': f'resources/grn_benchmark/prior/regulators_consensus_{dataset}.json'
    }
    print(par)
    # Load perturbation data
    adata_rna = anndata.read_h5ad(par['evaluation_data'])
    gene_names = adata_rna.var_names

    gene_dict = {gene_name: i for i, gene_name in enumerate(gene_names)}

    # Load inferred GRNs
    grns = []
    for model in par['models']:
        filepath = os.path.join(par['models_dir'], f'{model}.h5ad')
        if not os.path.exists(filepath):
            print(f"{filepath} didnt exist. Skipped.")
            continue 
        if model == 'collectri':
            print(f"skipping collectri")
            continue
        A = np.zeros((len(gene_names), len(gene_names)), dtype=float)
        net = ad.read_h5ad(filepath)
        net = pd.DataFrame(net.uns['prediction'])
        net['weight'] = net['weight'].astype(float)

        for source, target, weight in zip(net['source'], net['target'], net['weight']):
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
    for dataset in ['replogle', 'norman', 'adamson', 'nakatake', 'op']:
        main(dataset)

