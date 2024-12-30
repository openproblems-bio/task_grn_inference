import os
import json

import anndata
import pandas as pd
import anndata as ad
import sys
import numpy as np
import scipy
import random
from tqdm import tqdm
import scanpy as sc

## VIASH START
par = {
     'evaluation_data': 'resources/datasets_raw/norman_sc_counts.h5ad',
     'prediction': 'resources/grn_models/norman/grnboost2.csv',
     'tf_all': 'resources/prior/tf_all.csv',
     'max_n_links': 50_000
}
## VIASH END

def main(par):

    def get_observational(child: str) -> np.array:
        """
        Return all the samples for gene "child" in cells where there was no perturbations
        

        Args:
            child: Gene name of child to get samples for

        Returns:
            np.array matrix of corresponding samples
        """
        mask_gene = gene_names==child
        mask_sample =  adata_rna.obs['is_control']
        
        X = adata_rna[mask_sample, mask_gene].X.todense().A

        return  X 

    def get_interventional(child: str, parent: str) -> np.array:
        """
        Return all the samples for gene "child" in cells where "parent" was perturbed

        Args:
            child: Gene name of child to get samples for
            parent: Gene name of gene that must have been perturbed

        Returns:
            np.array matrix of corresponding samples
        """
        mask_gene = gene_names==child
        assert parent in adata_rna.obs['perturbation'].unique()


        mask_sample =  adata_rna.obs['perturbation']==parent
        # print(mask_sample.sum())
        
        X = adata_rna[mask_sample, mask_gene].X.todense().A
        return X
    # - read adata and normalize
    adata_rna = anndata.read_h5ad(par['evaluation_data'])
    sc.pp.normalize_total(adata_rna)
    sc.pp.log1p(adata_rna)
    gene_names = adata_rna.var_names

    # - read the net and retain only those tfs that are mutually available in both adata perturbation and net
    net = pd.read_csv(par['prediction'])
    tfs_pertubed = adata_rna[~adata_rna.obs['is_control']].obs['perturbation'].unique()
    tfs_common = np.intersect1d(tfs_pertubed, net['source'].unique())

    net = net[net['source'].isin(tfs_common)]
    net = net.nlargest(par['max_n_links'], 'weight')

    print('Remaining net size: ', net.shape, ' TF size: ', net['source'].nunique(), ' common TFs: ', tfs_common.shape)

    # - calculate the scores
    wasserstein_distances = []
    links = []
    for tf in tqdm(tfs_common):
        edges = net[net['source']==tf]
        for parent, child in zip(edges['source'],edges['target']): #multiple cuts
            
            get_observational(child)
            observational_samples = get_observational(child)
            interventional_samples = get_interventional(child, parent)

            wasserstein_distance = scipy.stats.wasserstein_distance(
                observational_samples.reshape(-1), interventional_samples.reshape(-1)
            )
            wasserstein_distances.append(wasserstein_distance)
            links.append(f'{parent}_{child}')
    mean_score = np.mean(wasserstein_distances)
    return mean_score, wasserstein_distances, links

if __name__ == '__main__':
    if False:
        mean_score, wasserstein_distances, links = main(par)
        #TODO: put this into right format
    else:
        output_dir = 'output'
        datasets = ['adamson', 'norman']
        models = ['pearson_corr', 'grnboost2','portia', 'ppcor','scenic']
        n_maxs = [500, 1000, 5000, 10000, 50000]
        

        for dataset in datasets:
            print(dataset)
            scores_all = []
            for model in models:
                par['evaluation_data'] = f'resources/datasets_raw/{dataset}_sc_counts.h5ad'
                par['prediction'] = f'resources/grn_models/{dataset}/{model}.csv'
                if not os.path.exists(par['prediction']):
                    print(f'Skip {dataset}-{model}')
                    continue
                for n_max in n_maxs:
                    par['max_n_links'] = n_max
                    
                    _, wasserstein_distances, links = main(par)
                    for score, link in zip(wasserstein_distances, links):
                        scores_all.append({'model':model, 'n_max':n_max, 'link':link, 'score':score})
            scores_all = pd.DataFrame(scores_all)
            scores_all.to_csv(f'{output_dir}/scores_{dataset}.csv')
                