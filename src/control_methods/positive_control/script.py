import os
import pandas as pd
import numpy as np
import anndata as ad
import scanpy as sc
from tqdm import tqdm
from scipy.stats import spearmanr
from sklearn.preprocessing import StandardScaler

## VIASH START
par = {
    'perturbation_data': 'resources/grn-benchmark/perturbation_data.h5ad',
    'tf_all': 'resources/prior/tf_all.csv',
    'causal': True,
    'cell_type_specific': True,
    'max_n_links': 50000,
    'prediction': 'resources/grn_models/positive_control.csv',
    "seed": 32,
    'normalize': False}
## VIASH END
import sys
sys.path.append(meta["resources_dir"])
from util import create_corr_net

print('Create causal corr net')
par['causal'] = True
par['multiomics_rna'] = par['perturbation_data']
par['only_hvgs'] = False

net = create_corr_net(par)

print('Output GRN')
net.to_csv(par['prediction'])
