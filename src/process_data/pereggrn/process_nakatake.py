
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
    'input_data': f'resources/datasets_raw/nakatake.h5ad',
    'adata_bulk': f'resources/extended_data/nakatake_bulk.h5ad',
    'adata_test_bulk': f'resources/grn_benchmark/evaluation_data/nakatake_bulk.h5ad',
    'adata_train_bulk': f'resources/grn_benchmark/inference_data/nakatake_rna.h5ad'
}

meta = {
    'resources_dir': 'src/utils/'
}
sys.path.append(meta["resources_dir"])

from util import sum_by


# - get the data
adata = ad.read_h5ad(par['input_data'])
adata = adata[:, ~adata.var_names.duplicated()]
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
ctr_pertb = adata_bulk[adata_bulk.obs['is_control']].obs['perturbation'].unique()
non_ctr_pertubs =adata_bulk[~adata_bulk.obs['is_control']].obs['perturbation'].unique()
train_perturbs, test_perturbs = train_test_split(non_ctr_pertubs, test_size=.5, random_state=32)
train_perturbs = np.concatenate([train_perturbs, ctr_pertb]) # add control perturbations to test set for ws_distance
test_perturbs = np.concatenate([test_perturbs, ctr_pertb]) 

adata_train_bulk = adata_bulk[adata_bulk.obs['perturbation'].isin(train_perturbs)] 
adata_test_bulk = adata_bulk[adata_bulk.obs['perturbation'].isin(test_perturbs)] 


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
adata_train_bulk.uns['dataset_id'] = 'nakatake'
adata_test_bulk.uns['dataset_id'] = 'nakatake'
adata_bulk.uns['dataset_id'] = 'nakatake'

# - save 
adata_bulk.write(par['adata_bulk'])
adata_test_bulk.write(par['adata_test_bulk'])
adata_train_bulk.write(par['adata_train_bulk'])