import pandas as pd
import numpy as np
import anndata as ad
import scanpy as sc
import os

def add_metadata(adata):
    adata.uns['dataset_summary'] = '300BCG scRNA cohort with 38 donors, 2 time points and 2 perturbaions (control and LPS)'
    adata.uns['dataset_description'] = '300BCG scRNA cohort with 38 donors, 2 time points and 2 perturbaions (control and LPS)'
    adata.uns['data_reference'] = "Not published"
    adata.uns['data_url'] = ''
    adata.uns['dataset_id'] = '300BCG'
    adata.uns['dataset_name'] = '300BCG'
    adata.uns['dataset_organism'] = 'human'
    adata.uns['normalization_id'] = 'lognorm'
    return adata

adata = ad.read_h5ad('/vol/projects/CIIM/300BCG/300BCG_scRNA/bcg4-0712.h5ad')

def format_data(adata):
    adata.obs.rename({'stim': 'perturbation', 'ids': 'donor_id', 'clusters1': 'cell_type'}, axis=1, inplace=True)
    adata.obs = adata.obs[['age', 'perturbation', 'time', 'donor_id', 'batch', 'pool', 'cell_type']]
    adata.obs.loc[:, 'is_control'] = adata.obs['perturbation']=='RPMI'
    adata.obs['perturbation_type'] = 'chemical'
    return adata
adata = format_data(adata)

adata.obs['qc_group'] = adata.obs['perturbation'] + '_' + adata.obs['time'].astype(str) + '_' + adata.obs['donor_id'].astype(str) + '_' + adata.obs['cell_type'].astype(str)

adata.obs['split_group'] = adata.obs['perturbation'] + '_' + adata.obs['time'].astype(str)
adata.obs['split_group'].value_counts()

def split_data_func(adata: ad.AnnData):
    obs = adata.obs
    
    train_group = ['RPMI_T0', 'LPS_T0']
    test_group = np.setdiff1d(obs['split_group'].unique(), train_group)

    return train_group, test_group

from src.process_data.helper_data import wrapper_large_perturbation_data
wrapper_large_perturbation_data(adata, covariates=['perturbation', 'time', 'donor_id', 'cell_type', 'qc_group'], 
                                    add_metadata=add_metadata, 
                                    save_name='300BCG', 
                                    split_func=split_data_func, 
                                    qc_group='qc_group', 
                                    qc_perturbation_effect=False,
                                    group='split_group')