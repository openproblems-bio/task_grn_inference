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
    'rna_all': 'resources/grn-benchmark/rna.h5ad',
    'tf_all': 'resources/grn_benchmark/prior/tf_all.csv',
    'max_n_links': 50000,
    'prediction': 'resources/grn_models/positive_control.csv',
    "seed": 32,
    'temp_dir': 'output/positive_control',
    'apply_tf_methods': True,
}
## VIASH END
import argparse
argparser = argparse.ArgumentParser()
argparser.add_argument('--rna_all', type=str, help='Path to the input RNA data in h5ad format.')
argparser.add_argument('--prediction', type=str, help='Path to the output prediction in h5ad format.')
args = argparser.parse_args()
if args.rna_all is not None:
  par['rna_all'] = args.rna_all
if args.prediction is not None:
  par['prediction'] = args.prediction


meta = {
    "resources_dir": 'src/utils',
    "name": "positive_control"
}

sys.path.append(meta["resources_dir"])
from util import corr_net

adata = ad.read_h5ad(par["rna_all"])
par['layer'] = 'lognorm' if 'lognorm' in adata.layers.keys() else 'X_norm'
tf_all = np.loadtxt(par["tf_all"], dtype=str)
dataset_id = adata.uns['dataset_id']
net = corr_net(adata, tf_all, par)
print('Shape of the network:', net.shape)
print(net.sort_values('weight', ascending=False, key=abs).head(10))

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
