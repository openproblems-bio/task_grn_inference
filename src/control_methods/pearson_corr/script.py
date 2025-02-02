
import anndata as ad 
import pandas as pd
import os
import scanpy as sc 
import sys
import numpy as np

# - whatever is between viash start and end will be replaced by Viash with the parameters from the config. file

## VIASH START
par = {
    'rna': 'resources/evaluation_datasets/op_rna.h5ad',
    'tf_all': 'resources/prior/tf_all.csv',
    'cell_type_specific': False,
    'max_n_links': 50000,
    'prediction': 'output/pearson_net.h5ad',
    'apply_tf': True,
    'normalize': True}
## VIASH END

sys.path.append(meta["resources_dir"])
from util import corr_net


def infer_net(par: dict) -> pd.DataFrame:
    print(par)
    print('Read data')
    adata = ad.read_h5ad(par["rna"])
    try:
        X = adata.layers['X_norm'].todense().A
    except:
        X = adata.X

    # - remove genes with 0 standard deviation
    gene_std = np.std(X, axis=0)
    nonzero_std_genes = gene_std > 0
    X = X[:, nonzero_std_genes]
    # - get the net
    gene_names = adata[:, nonzero_std_genes].var_names.to_numpy()
    grn = corr_net(X, gene_names, par)    
    return grn

if __name__ == '__main__':
    net = infer_net(par)
    # - format of et
    '''
        the net is a pandas dataframe with the following columns:
            - source: the source gene of the interaction
            - target: the target gene of the interaction
            - weight: the weight of the interaction
    '''

    print('Output GRN')
    net['weight'] = net['weight'].astype(str)
    output = ad.AnnData(X=None, uns={"method_id": par['method_id'], "dataset_id": par['dataset_id'], "prediction": net[["source", "target", "weight"]]})
    output.write(par['prediction'])
