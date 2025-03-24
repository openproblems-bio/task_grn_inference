import os
import pandas as pd
import numpy as np
import sys
import anndata as ad
import scanpy as sc
from sklearn.impute import KNNImputer
import magic

## VIASH START
par = {
    "rna": "resources/grn_benchmark/inference_data/op_rna.h5ad",
    'layer': 'X_norm',
    'imputation': 'knn',
    "rna_imputed": "output/imputed.h5ad"
}
## VIASH END


def main(par):

    # - read inputs and cluster with differen resolutions
    rna = ad.read_h5ad(par['rna'])
    # rna.X = rna.layers[par['layer']]

    sc.pp.pca(rna, layer=par['layer'])
    sc.pp.neighbors(rna)

    X_original = rna.layers[par['layer']]
    try:
        X_original = X_original.todense().A 
    except:
        pass
    
    # - pseudobulkd and run per res
    gene_names = rna.var_names
    imputation = par['imputation']
    if imputation == 'original':
        X = rna.X.todense().A
    elif imputation == 'knn':
        X_original[X_original == 0] = np.nan
        knn_imputer = KNNImputer(n_neighbors=10) 
        X = knn_imputer.fit_transform(X_original)
    elif imputation == 'magic':
        magic_operator = magic.MAGIC()
        X = magic_operator.fit_transform(X_original)
    else:
        raise ValueError('define first')
    
    rna.layers[par['layer']] = X
    rna.write(par['rna_imputed'])
if __name__ == '__main__':
    main(par)