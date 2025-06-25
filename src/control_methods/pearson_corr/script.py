
import anndata as ad 
import pandas as pd
import os
import numba
numba.config.CACHE = False
import scanpy as sc 
import sys
import numpy as np

## VIASH START
par = {
    'rna': 'resources/grn_benchmark/inference_data//op_rna.h5ad',
    'tf_all': 'resources/grn_benchmark/prior/tf_all.csv',
    'cell_type_specific': False,
    'max_n_links': 50000,
    'prediction': 'output/pearson_net.h5ad',
    'layer': 'X_norm',
    'apply_tf_methods': True
}
meta = {
    "resources_dir": 'src/utils'
}
## VIASH END

sys.path.append(meta["resources_dir"])
from util import corr_net

adata = ad.read_h5ad(par["rna"])
tf_all = np.loadtxt(par["tf_all"], dtype=str)
dataset_id = adata.uns['dataset_id']
net = corr_net(adata, tf_all, par)

print('Output GRN')
print('Shape of the network:', net.shape)
print(net.sort_values('weight', ascending=False, key=abs).head(10))
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
