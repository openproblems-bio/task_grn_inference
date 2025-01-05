import anndata as ad 
from tqdm import tqdm
import sys
import scipy
import pandas as pd
import numpy as np



# - create background dist. 
par = {
    'evaluation_data_sc': f'resources/datasets_raw/norman_sc_counts.h5ad',
    'background_distance': 'resources/prior/ws_distance_background_norman.csv',
    'tf_all': 'resources/prior/tf_all.csv',
    'layer': 'X_norm'
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


def main(par):
    adata = ad.read_h5ad(par['evaluation_data_sc'])
    adata.X = adata.layers[par['layer']]

    tf_all = np.loadtxt(par['tf_all'], dtype='str')
    available_tfs = np.intersect1d(adata.obs['perturbation'].unique(), tf_all)

    # - for each tf, select a random background net and calculate the distances

    random_grn_store = []
    for tf in tqdm(available_tfs):
        # random_genes = np.random.choice(adata.var_names, par['n_selection'], replace=False)
        random_grn = pd.DataFrame([{'source':tf, 'target': gene} for gene in adata.var_names]) 
        
        random_grn = calculate_ws_distance(random_grn, adata)
        random_grn_store.append(random_grn)
    random_grn_ws_distance = pd.concat(random_grn_store)
    random_grn_ws_distance.to_csv(par['background_distance'])
if __name__ == '__main__':
    main(par)