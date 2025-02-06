import anndata as ad 
from tqdm import tqdm
import sys
import scipy
import pandas as pd
import numpy as np

from multiprocessing import Pool



# - create background dist. 
par = {
    'evaluation_data_sc': f'resources/grn_benchmark/evaluation_data//replogle_sc.h5ad',
    'background_distance': 'resources/grn_benchmark/prior/ws_distance_background_replogle.csv',
    'tf_all': 'resources/grn_benchmark/prior/tf_all.csv',
    'layer': 'X_norm',
    'max_workers': 100
}


def calculate_ws_distance(net, adata) -> pd.DataFrame:
    """
    Get a net with source and target columns, and adata with gene expression data for control and interventaion, 
    and returns the net with added column of wasserstein distance.
    """
    gene_names = adata.var_names
    ws_distances = []

    for parent, child in zip(net['source'],net['target']): 
        # - get observational X
        mask_gene = gene_names==child
        if 'is_control' not in adata.obs.columns:
            adata.obs['is_control'] = adata.obs['perturbation']=='non-targeting'
        
        mask_sample =  adata.obs['is_control']
        X_observ = adata[mask_sample, mask_gene].X.todense().A
        # - get interventional
        mask_gene = gene_names==child
        assert parent in adata.obs['perturbation'].unique()
        mask_sample =  adata.obs['perturbation']==parent        
        X_interv = adata[mask_sample, mask_gene].X.todense().A

        assert X_observ.shape[0] != 0
        assert X_interv.shape[0] != 0

        # - calculate the distance
        ws_distance = scipy.stats.wasserstein_distance(
            X_observ.reshape(-1), X_interv.reshape(-1)
        )
        ws_distances.append(ws_distance)
    net['ws_distance'] = ws_distances
    return net


# def main(par):
#     adata = ad.read_h5ad(par['evaluation_data_sc'])
#     adata.X = adata.layers[par['layer']]

#     tf_all = np.loadtxt(par['tf_all'], dtype='str')
#     available_tfs = np.intersect1d(adata.obs['perturbation'].unique(), tf_all)

#     # - for each tf, calculate the distances for all possible genes in the network
#     net_all_store = []
#     for tf in tqdm(available_tfs):
#         net_all = pd.DataFrame([{'source':tf, 'target': gene} for gene in adata.var_names]) 
#         net_all = calculate_ws_distance(net_all, adata)
#         net_all_store.append(net_all)
#     net_all_ws_distance = pd.concat(net_all_store)
#     net_all_ws_distance.to_csv(par['background_distance'])

def load_adata(par):
    adata = ad.read_h5ad(par['evaluation_data_sc'])
    adata.X = adata.layers[par['layer']]
    return adata

def process_tf(tf_adata):
    """Calculate ws distance for a given TF"""
    adata = load_adata(par)
    tf = tf_adata  # Unpack tuple
    net_all = pd.DataFrame([{'source': tf, 'target': gene} for gene in adata.var_names])
    net_all = calculate_ws_distance(net_all, adata)
    return net_all

def main(par):
    adata = ad.read_h5ad(par['evaluation_data_sc'], backed='r')

    tf_all = np.loadtxt(par['tf_all'], dtype=str)
    available_tfs = np.intersect1d(adata.obs['perturbation'].unique(), tf_all)

    # Set up multiprocessing pool
    num_workers = min(par['max_workers'], len(available_tfs))  # Use available cores
    with Pool(num_workers) as pool:
        # tqdm now tracks progress inside imap_unordered()
        results = list(tqdm(pool.imap(process_tf, [(tf) for tf in available_tfs]), 
                            total=len(available_tfs), desc="Processing TFs"))

    # Combine results
    net_all_ws_distance = pd.concat(results)
    net_all_ws_distance.to_csv(par['background_distance'])
if __name__ == '__main__':
    main(par)