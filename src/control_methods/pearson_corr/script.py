
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
    'apply_tf_methods': True
}
## VIASH END
import argparse
argparser = argparse.ArgumentParser()
argparser.add_argument('--rna', type=str, help='Path to the input RNA data in h5ad format.')
argparser.add_argument('--prediction', type=str, help='Path to the output prediction in h5ad format.')
args = argparser.parse_args()
if args.rna is not None:
  par['rna'] = args.rna
if args.prediction is not None:
  par['prediction'] = args.prediction


try:
    sys.path.append(meta["resources_dir"])
except:
    meta = {
        "resources_dir": 'src/utils',
        "name": "pearson_corr"
    }
    sys.path.append(meta["resources_dir"])
from util import corr_net

def main(par):
    adata = ad.read_h5ad(par["rna"])
    par['layer'] = 'lognorm' if 'lognorm' in adata.layers.keys() else 'X_norm'
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
    return output

if __name__ == '__main__':
    output = main(par)
    output.write_h5ad(par['prediction'])

