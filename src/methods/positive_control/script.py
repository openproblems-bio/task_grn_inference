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
try:
    sys.path.append(meta["resources_dir"])
except:
  meta = {
      "resources_dir": 'src/utils',
      "name": "positive_control"
  }
  sys.path.append(meta["resources_dir"])
from util import parse_args, process_links, corr_net
par = parse_args(par)

if __name__ == '__main__':
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
