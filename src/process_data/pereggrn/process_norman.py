
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
    'input_data': f'resources/datasets_raw/norman.h5ad',
    'adata_bulk': f'resources/extended_data/norman_bulk.h5ad',
    'adata_test_sc': f'resources/grn_benchmark/evaluation_data/norman_sc.h5ad',
    'adata_test_bulk': f'resources/grn_benchmark/evaluation_data/norman_bulk.h5ad',
    'adata_train_sc': f'resources/grn_benchmark/inference_data/norman_rna.h5ad'
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
adata.var.index = adata.var.index.astype(str)
adata.obs = adata.obs[['perturbation', 'is_control', 'perturbation_type']]

# preprocess  
sc.pp.filter_cells(adata, min_genes=100)
sc.pp.filter_genes(adata, min_cells=10)

# - unique to norman
adata.layers['X_norm'] = adata.X.copy()

# - split to inference and evaluation datasets: unique to norman 
ctr_pertb = adata[adata.obs['is_control']].obs['perturbation'].unique()
non_ctr_pertubs =adata[~adata.obs['is_control']].obs['perturbation'].unique()
train_perturbs, test_perturbs = train_test_split(non_ctr_pertubs, test_size=.5, random_state=32)
train_perturbs = np.concatenate([train_perturbs, ctr_pertb]) # add control perturbations to test set for ws_distance
test_perturbs = np.concatenate([test_perturbs, ctr_pertb]) 

adata_train_sc = adata[adata.obs['perturbation'].isin(train_perturbs)] 
adata_test_sc = adata[adata.obs['perturbation'].isin(test_perturbs)] 


# - filter genes and cells 
sc.pp.filter_cells(adata_train_sc, min_genes=100)
sc.pp.filter_genes(adata_train_sc, min_cells=10)

sc.pp.filter_cells(adata_test_sc, min_genes=100)
sc.pp.filter_genes(adata_test_sc, min_cells=10)

# - pseudo bulk: unique to norman 
adata_bulk = sum_by(adata, unique_mapping=True, col='perturbation') 
adata_test_bulk = sum_by(adata_test_sc, unique_mapping=True, col='perturbation') # summing over X_norm 


# - normalize evaluation data
# sc.pp.normalize_total(adata_test_bulk)
adata_test_bulk.layers['X_norm'] = adata_test_bulk.X.copy()

# - normalize adata_train_sc
adata_train_sc.layers['X_norm'] = adata_train_sc.X.copy()

# - add metadata
adata_train_sc.uns['dataset_id'] = 'norman'
adata_test_sc.uns['dataset_id'] = 'norman'
adata_test_bulk.uns['dataset_id'] = 'norman'
adata_bulk.uns['dataset_id'] = 'norman'
# - save 
adata_bulk.write(par['adata_bulk'])
adata_test_sc.write(par['adata_test_sc'])
adata_test_bulk.write(par['adata_test_bulk'])
adata_train_sc.write(par['adata_train_sc'])