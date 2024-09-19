import os
import pandas as pd
import numpy as np
import anndata as ad
import scanpy as sc
from tqdm import tqdm
from scipy.stats import spearmanr
from sklearn.preprocessing import StandardScaler
from utils.util import process_data, create_corr_net

## VIASH START
par = {
    'multiomics_rna': 'resources/grn-benchmark/multiomics_rna_0.h5ad',
    'tf_all': 'resources/prior/tf_all.csv',
    'causal': False,
    'cell_type_specific': True,
    'max_n_links': 50000,
    'prediction': 'resources/grn_models/donor_0_default/pearson.csv',
    "seed": 32,
    'only_hvgs': True,
    'normalize': False
}
## VIASH END
print(par)

print('Read data')
multiomics_rna = ad.read_h5ad(par["multiomics_rna"])
process_data(multiomics_rna)
    
gene_names = multiomics_rna.var_names.to_numpy()
tf_all = np.loadtxt(par['tf_all'], dtype=str)
groups = multiomics_rna.obs.cell_type
tf_all = np.intersect1d(tf_all, gene_names)
n_tfs = len(tf_all)


print('Create corr net')
net = create_corr_net(multiomics_rna.X, multiomics_rna.var_names, groups, par, tf_all=None)

print('Output GRN')
net.to_csv(par['prediction'])
