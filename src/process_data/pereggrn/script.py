import os 
import anndata as ad 
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys
import scanpy as sc
from sklearn.model_selection import train_test_split 

from scipy.sparse import csr_matrix

meta = {
    'resources_dir': 'src/utils/'
}

sys.path.append(meta["resources_dir"])



def process_dataset(file_name):
    # - get the data
    adata = ad.read_h5ad(f'resources/datasets_raw/{file_name}.h5ad')
    adata.var_names_make_unique()
    adata.var.index = adata.var.index.astype(str)
    # - clearn up 
    del adata.obsp 
    del adata.varm
    del adata.uns
    del adata.obsm
    if 'gene_name' in adata.var.columns:
        adata.var = adata.var[['gene_name']]
        adata.var = adata.var.set_index('gene_name')
    else:
        adata.var = adata.var[[]]
    adata.obs = adata.obs[['perturbation', 'is_control', 'perturbation_type']]
    
    # - data split 
    if file_name == 'replogle':
        ctr_samples = adata.obs.is_control
        samples = adata.obs.index[~ctr_samples] 
        _, test_samples = train_test_split(samples, test_size=.2, random_state=32)
        adata.obs['is_test'] = adata.obs.index.isin(test_samples)
    elif file_name == 'norman':
        ctr_samples = adata.obs.is_control
        samples = adata[adata.obs.index[~ctr_samples]].obs.perturbation.unique()
        _, test_samples = train_test_split(samples, test_size=.5, random_state=32)
        adata.obs['is_test'] = adata.obs.perturbation.isin(test_samples)
    elif file_name == 'nakatake':
        samples = adata.obs.perturbation.unique()
        _, test_samples = train_test_split(samples, test_size=.5, random_state=32)
        adata.obs['is_test'] = adata.obs.perturbation.isin(test_samples)
    elif file_name == 'adamson':
        ctr_samples = adata.obs.is_control
        samples = adata[adata.obs.index[~ctr_samples]].obs.perturbation.unique()
        _, test_samples = train_test_split(samples, test_size=.8, random_state=32)
        adata.obs['is_test'] = adata.obs.perturbation.isin(test_samples)

    adata_train = adata[~adata.obs['is_test']] # we use single cells for train (if not already bulked)
    
    # - pseudo bulk if needed
    if file_name in ['norman', 'adamson']: # these two are single cells. for norman, we have .counts but not for adamson -> different preprocessing
        if file_name == 'norman':
            adata.X = adata.layers['counts'] 
        adata_bulked = psedudobulk_fun(adata) 
    else:
        adata_bulked = adata
        # adata_bulked.layers['X_norm'] = adata_bulked.X.copy()


    adata_test_bulk = adata_bulked[adata_bulked.obs['is_test']] # we use bulked data for feature-based evaluation 

    # - duplicated gene names
    duplicates = adata_train.var_names[adata_train.var_names.duplicated()].unique()
    adata_train = adata_train[:, ~adata_train.var_names.isin(duplicates)]

    duplicates = adata_test_bulk.var_names[adata_test_bulk.var_names.duplicated()].unique()
    adata_test_bulk = adata_test_bulk[:, ~adata_test_bulk.var_names.isin(duplicates)]


    adata_test_bulk = adata_test_bulk.copy()  # Ensure it's a full AnnData object
    adata_train = adata_train.copy()  # Ensure it's a full AnnData object
    adata = adata.copy()  # Ensure it's a full AnnData object

    
    print(adata_train)
    if file_name == 'norman':
        adata_train.layers['counts'] = adata_train.X.copy()
        adata_train.X = adata_train.layers['X_norm']
    else:
        adata_train.layers['X_norm'] = adata_train.X.copy()
        adata_test_bulk.layers['X_norm'] = adata_test_bulk.X.copy()
        adata.layers['X_norm'] = adata.X.copy()

    if file_nam == 'norman': # only dataset with counts
        adata_train.X = adata_train.layers['counts']

    if file_name in ['norman', 'adamson']:
        adata.write(f'resources/grn_benchmark/evaluation_data/{file_name}_sc.h5ad')

    adata_bulked.write(f'resources/extended_data/{file_name}_bulk.h5ad')
    adata_train.write(f'resources/grn_benchmark/inference_data/{file_name}_rna.h5ad')
    adata_test_bulk.write(f'resources/grn_benchmark/evaluation_data/{file_name}_bulk.h5ad')
    

def main(par):
    for file_name in par['datasets']:
        process_dataset(file_name)

if __name__ == '__main__':
    
    par = {
        'datasets': ['norman', 'adamson', 'nakatake'] #'norman'
    }
    

    main(par)