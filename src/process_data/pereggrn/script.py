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

from util import sum_by

def psedudobulk_fun(adata: ad.AnnData) -> ad.AnnData:
    metadata = (
        adata.obs.groupby('perturbation')
        .agg(lambda x: x.mode()[0] if x.nunique() == 1 else x.iloc[0])  
    )

    pseudobulk_data = adata.to_df().groupby(adata.obs['perturbation']).sum()
    # Ensure the metadata index matches pseudobulk counts
    metadata = metadata.loc[pseudobulk_data.index]

    # Create a new AnnData object for pseudobulked data
    adata_bulked = sc.AnnData(
            X=pseudobulk_data.values,
            obs=metadata.reset_index(),
            var=adata.var.copy()
        )

    return adata_bulked

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
    
    if file_name in ['norman', 'adamson']: # these two are single cells. for norman, we have .counts but not for adamson -> different preprocessing
        adata_bulked = psedudobulk_fun(adata) # also normalize
    else:
        adata_bulked = adata
        # adata_bulked.layers['X_norm'] = adata_bulked.X.copy()

    adata_test = adata_bulked[adata_bulked.obs['is_test']] # we use bulked data for test 


    # - duplicated gene names
    duplicates = adata_train.var_names[adata_train.var_names.duplicated()].unique()
    adata_train = adata_train[:, ~adata_train.var_names.isin(duplicates)]

    duplicates = adata_test.var_names[adata_test.var_names.duplicated()].unique()
    adata_test = adata_test[:, ~adata_test.var_names.isin(duplicates)]


    adata_test = adata_test.copy()  # Ensure it's a full AnnData object
    adata_train = adata_train.copy()  # Ensure it's a full AnnData object
    adata = adata.copy()  # Ensure it's a full AnnData object

    adata_train.layers['X_norm'] = adata_train.X.copy()
    adata_test.layers['X_norm'] = adata_test.X.copy()
    adata.layers['X_norm'] = adata.X.copy()

    if file_name in ['norman', 'adamson']:
        adata.write(f'resources/grn_benchmark/evaluation_data/{file_name}_sc.h5ad')

    adata_bulked.write(f'resources/extended_data/{file_name}_bulk.h5ad')
    adata_train.write(f'resources/grn_benchmark/inference_datasets/{file_name}_rna.h5ad')
    adata_test.write(f'resources/grn_benchmark/evaluation_data/{file_name}_bulk.h5ad')
    

def main(par):
    for file_name in par['datasets']:
        process_dataset(file_name)

if __name__ == '__main__':
    
    par = {
        'datasets': ['norman', 'adamson', 'nakatake']
    }
    

    main(par)