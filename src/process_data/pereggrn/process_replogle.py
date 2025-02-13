
import os 
import anndata as ad 
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys
import scanpy as sc
from sklearn.model_selection import train_test_split 

from scipy.sparse import csr_matrix

par = {
    'input_data': f'resources/datasets_raw/replogle.h5ad',
    'adata_bulk': f'resources/extended_data/replogle_bulk.h5ad',
    'adata_test_bulk': f'resources/grn_benchmark/evaluation_data/replogle_bulk.h5ad',
    'adata_train_bulk': f'resources/grn_benchmark/inference_data/replogle_rna.h5ad',
    'test_perturbs': f'resources/grn_benchmark/prior/replogle_test_perturbs.csv',
    'tf_all': f'resources/grn_benchmark/prior/tf_all.csv',
}

meta = {
    'resources_dir': 'src/utils/'
}
sys.path.append(meta["resources_dir"])

from util import sum_by


# - get the data
adata = ad.read_h5ad(par['input_data'])
adata = adata[:, ~adata.var_names.duplicated()]
tf_all = np.loadtxt(par['tf_all'], dtype=str)
adata.var_names_make_unique()
duplicates = adata.var_names[adata.var_names.duplicated()].unique()
adata = adata[:, ~adata.var_names.isin(duplicates)]
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
adata_bulk = adata.copy()
del adata 
adata_bulk.var.index = adata_bulk.var.index.astype(str)
adata_bulk.obs = adata_bulk.obs[['perturbation', 'is_control', 'perturbation_type']]

# preprocess  
sc.pp.filter_cells(adata_bulk, min_genes=100)
sc.pp.filter_genes(adata_bulk, min_cells=10)


# - split to inference and evaluation datasets
if False:
    ctr_pertb = adata_bulk[adata_bulk.obs['is_control']].obs['perturbation'].unique()
    non_ctr_pertubs =adata_bulk[~adata_bulk.obs['is_control']].obs['perturbation'].unique()
    train_perturbs, test_perturbs = train_test_split(non_ctr_pertubs, test_size=.2, random_state=32)
    train_perturbs = np.concatenate([train_perturbs, ctr_pertb]) # add control perturbations to test set for ws_distance
    test_perturbs = np.concatenate([test_perturbs, ctr_pertb]) 
    adata_train_bulk = adata_bulk[adata_bulk.obs['perturbation'].isin(train_perturbs)] 
    adata_test_bulk = adata_bulk[adata_bulk.obs['perturbation'].isin(test_perturbs)] 

else:
    obs = adata_bulk.obs
    obs['is_tf'] = obs['perturbation'].isin(tf_all)
    unique_perts = obs['perturbation'].unique()
    
    tf_all = obs[obs['is_tf']]['perturbation'].unique()
    non_tfs = np.setdiff1d(unique_perts, tf_all)
    train_perturbs, test_perturbs = train_test_split(non_tfs, test_size=.2, random_state=32)

    n_tfs_in_test = int(len(tf_all)/2)

    np.random.seed(32)
    test_tfs = np.random.choice(tf_all, size=n_tfs_in_test, replace=False) 
    train_tfs = np.setdiff1d(tf_all, test_tfs)
    test_perturbs = np.concatenate([test_tfs, test_perturbs])

    train_perturbs = np.concatenate([train_tfs, train_perturbs])

    
    print(f"Test TFs: {n_tfs_in_test}, Train TFs: {len(tf_all) - n_tfs_in_test}")


    adata_test_bulk = adata_bulk[adata_bulk.obs['perturbation'].isin(test_perturbs)]
    adata_train_bulk = adata_bulk[adata_bulk.obs['perturbation'].isin(train_perturbs)]

    print(f"Train : {adata_train_bulk.shape}, {adata_train_bulk.obs['perturbation'].nunique()} \
            Test: {adata_test_bulk.shape}", {adata_test_bulk.obs['perturbation'].nunique()})




# - filter genes and cells 
sc.pp.filter_cells(adata_train_bulk, min_genes=100)
sc.pp.filter_genes(adata_train_bulk, min_cells=10)

sc.pp.filter_cells(adata_test_bulk, min_genes=100)
sc.pp.filter_genes(adata_test_bulk, min_cells=10)


# - normalize 
adata_bulk.layers['X_norm'] = adata_bulk.X.copy()
adata_test_bulk.layers['X_norm'] = adata_test_bulk.X.copy()
adata_train_bulk.layers['X_norm'] = adata_train_bulk.X.copy()

# - add metadata
adata_train_bulk.uns['dataset_id'] = 'replogle'
adata_test_bulk.uns['dataset_id'] = 'replogle'
adata_bulk.uns['dataset_id'] = 'replogle'

# - save 
adata_bulk.write(par['adata_bulk'])
adata_test_bulk.write(par['adata_test_bulk'])
adata_train_bulk.write(par['adata_train_bulk'])

np.savetxt(par['test_perturbs'], test_perturbs, fmt="%s",  delimiter=",")
