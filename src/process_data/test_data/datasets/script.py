# !pip install sctk anndata
# !aws s3 cp s3://openproblems-bio/public/neurips-2023-competition/sc_counts.h5ad  ./resources_raw/ --no-sign-request

import anndata as ad 
import pandas as pd
import numpy as np
from scipy import sparse
import scanpy as sc
import requests
import os

## VIASH START

par = {
    # - inference data
    'rna': 'resources/grn_benchmark/inference_data/op_rna.h5ad',
    'rna_test': 'resources_test/grn_benchmark/inference_data//op_rna.h5ad',
    'atac': 'resources/grn_benchmark/inference_data/op_atac.h5ad',
    'atac_test': 'resources_test/grn_benchmark/inference_data//op_atac.h5ad',
    # - evaluation data
    'evaluation_data': 'resources/grn_benchmark/evaluation_data//op_bulk.h5ad',
    'op_perturbation_test': 'resources_test/grn_benchmark/evaluation_data//op_bulk.h5ad',
    'evaluation_data_sc': 'resources/grn_benchmark/evaluation_data/norman_sc.h5ad',
    'evaluation_data_sc_test': 'resources_test/grn_benchmark/evaluation_data/norman_sc.h5ad'
}
## VIASH END

def shorten_evaluation_data_sc(par):
    evaluation_data_sc = ad.read_h5ad(par['evaluation_data_sc'])
    # - control 
    ctr = evaluation_data_sc[evaluation_data_sc.obs['is_control']].obs['perturbation'].unique()

    selected_perturbation = evaluation_data_sc.obs.groupby(['perturbation']).size().sort_values()[:2].index

    selected_perturbation = list(set(list(ctr)+list(selected_perturbation)))
    print(selected_perturbation)
    mask = evaluation_data_sc.obs['perturbation'].isin(selected_perturbation)
    evaluation_data_sc = evaluation_data_sc[mask]
    print(evaluation_data_sc)
    evaluation_data_sc.write(par['evaluation_data_sc_test'])
if __name__ == '__main__':
    shorten_inference_data(par)
    shorten_evaluation_data(par)
    shorten_evaluation_data_sc(par)

# # - test datasets for raw counts
# os.makedirs('resources_test/datasets_raw/', exist_ok=True)
# multiomics_counts = ad.read_h5ad(par['multiomics_counts'])
# # shorten multiomics_counts 
# mask = multiomics_counts.obs.donor_id=='CCL-240'
# multiomics_counts_s = multiomics_counts[mask]
# # only one chr and n_cells 
# if 'interval' not in multiomics_counts_s.var:
#     raise ValueError('location is not given in rna')
# chr_mask = multiomics_counts_s.var.interval.str.split(':').str[0] == 'chr1'
# multiomics_counts_s = multiomics_counts_s[multiomics_counts_s.obs.sample(2000).index, multiomics_counts_s.var.sample(2000).index]
# multiomics_counts_s.write(par['multiomics_counts_test'])
# shorten perturbation
# evaluation_data_sc = ad.read_h5ad(par['perturbation_counts'])
# perturbation_counts_s = evaluation_data_sc[evaluation_data_sc.obs.donor_id=='Donor 1']
# evaluation_data_sc = evaluation_data_sc[evaluation_data_sc.obs.sample(1000).index,evaluation_data_sc.var.sample(1000).index]
# evaluation_data_sc.write(par['evaluation_data_sc_test'])

# # shorten sc-counts 

# adata = ad.read_h5ad(par['perturbation_counts'])
# few_perturbations = adata.obs.groupby('perturbation').size().sort_values()[:5].index
# adata = adata[adata.obs['perturbation'].isin(few_perturbations)]
# adata.write(par['evaluation_data_sc_test'])
