import os
import json

import anndata
import pandas as pd
import anndata as ad
import sys
import numpy as np
import scipy
import random


## VIASH START
par = {
    'evaluation_data': 'resources/grn-benchmark/evaluation_data.h5ad',
     'models_dir': 'resources/grn-benchmark/grn_models/d0_hvg',
     'models': [pearson_corr, pearson_causal, portia, ppcor, genie3, grnboost2, scenic, scglue, celloracle],
     'negative_samples': 100,
     'p_value_threshold': 0.01,
     'file': 'resources/prior/causalbench_evaluation.json'
}
## VIASH END
def main(par):
    print(par)
    # Load perturbation data
    adata_rna = anndata.read_h5ad(par['evaluation_data'])
    gene_names = adata_rna.var_names

    gene_dict = {gene_name: i for i, gene_name in enumerate(gene_names)}
    
    def get_observational(self, child: str) -> np.array:
        """
        Return all the samples for gene "child" in cells where there was no perturbations

        Args:
            child: Gene name of child to get samples for

        Returns:
            np.array matrix of corresponding samples
        """
        return tbd
    
    def get_interventional(self, child: str, parent: str) -> np.array:
        """
        Return all the samples for gene "child" in cells where "parent" was perturbed

        Args:
            child: Gene name of child to get samples for
            parent: Gene name of gene that must have been perturbed

        Returns:
            np.array matrix of corresponding samples
        """
        return tbd

    # Load inferred GRNs
    grns = {}
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

        network_as_dict = {}
        for source, target, weight in zip(df['source'], df['target'], df['weight']):
            if (source not in gene_dict) or (target not in gene_dict):
                continue
            i = gene_dict[source]
            j = gene_dict[target]
            A[i, j] = float(weight)
            if weight != 0:
                network_as_dict.setdefault(i, set()).add(j)
        print(f'Sparsity of {filepath}: {np.mean(A == 0)}')
        grns[model] = network_as_dict
    
    # Compute wasserstein distances
    def wasserstein_distances(network_as_dict):                   
        wasserstein_distances = []
        for parent in network_as_dict.keys():
            children = network_as_dict[parent]
            for child in children:
                observational_samples = get_observational(child)
                interventional_samples = get_interventional(child, parent)
                wasserstein_distance = scipy.stats.wasserstein_distance(
                    observational_samples, interventional_samples, 
                )
                wasserstein_distances.append(wasserstein_distance)
        return wasserstein_distances
    
    def false_omission_rate(network_as_dict):                  
        true_positive = 0
        tests = 0
        for parent in network_as_dict.keys():
            children = network_as_dict[parent]
            for child in children:
                observational_samples = get_observational(child)
                interventional_samples = get_interventional(child, parent)
                ranksum_result = scipy.stats.mannwhitneyu(
                    observational_samples, interventional_samples
                )
                tests +=1
                p_value = ranksum_result[1]
                if p_value < par['p_value_threshold']:
                 true_positive += 1
        return true_positive / tests
    
    def negative_sampling(network_as_dict, negative_samples = 100, seed = 0):
        edges = set()
        random.seed(seed)
        while len(edges) < negative_samples:
            edge = random.sample(range(len(gene_names)), 2)
            if edge[0] in network_as_dict and edge[1] in network_as_dict[edge[0]]:
                continue
            edges.add(edge)
        negative_dict = {}
        for a, b in edges:
            negative_dict.setdefault(a, set()).add(b)
            
        return negative_dict

    results = {}
    for model in grns.keys():
        dists = wasserstein_distances(grns[model])
        negative_dict = negative_sampling(grns[model], par['negative_samples'])
        neg_dists = wasserstein_distances(negative_dict)
        fo_rate = false_omission_rate(negative_dict)
        
        results[model] = {'mean_wasserstein': np.mean(dists),
                          'mean_negative_wasserstein': np.mean(neg_dists),
                          'fo_rate': fo_rate,
                          'wasserstein_dist': dists,
                          'negative_wasserstein_dist': neg_dists}

    # Store results
    with open(par['file'], 'w') as f:
        json.dump(results, f)


if __name__ == '__main__':
    main(par)
