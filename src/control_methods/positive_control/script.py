import os
import pandas as pd
import numpy as np
import anndata as ad
import scanpy as sc
from tqdm import tqdm
from scipy.stats import spearmanr
from sklearn.preprocessing import StandardScaler
import scanpy as sc 
import sys

## VIASH START
par = {
    'rna': 'resources/grn-benchmark/rna.h5ad',
    'tf_all': 'resources/grn_benchmark/prior/tf_all.csv',
    'max_n_links': 50000,
    'prediction': 'resources/grn_models/positive_control.csv',
    "seed": 32,
    'temp_dir': 'output/positive_control',
}
meta = {
    "resources_dir": 'src/utils'
}
## VIASH END

sys.path.append(meta["resources_dir"])
from util import corr_net

adata = ad.read_h5ad(par["rna_all"])
tf_all = np.loadtxt(par["tf_all"], dtype=str)
dataset_id = adata.uns['dataset_id']
net = corr_net(adata, tf_all, par)

print('Output GRN')
net = net.astype(str)
output = ad.AnnData(
    X=None,
    uns={
        "method_id": meta['name'],
        "dataset_id": adata.uns['dataset_id'],
        "prediction": net[["source", "target", "weight"]]
    }
)
output.write_h5ad(par['prediction'])
