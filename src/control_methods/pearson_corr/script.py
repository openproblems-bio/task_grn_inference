
import anndata as ad 
import pandas as pd
import os
import scanpy as sc 
import sys
import numpy as np
import argparse


## VIASH START
par = {
    'rna': 'resources/grn_benchmark/inference_data//op_rna.h5ad',
    'tf_all': 'resources/grn_benchmark/prior/tf_all.csv',
    'cell_type_specific': False,
    'max_n_links': 50000,
    'prediction': 'output/pearson_net.h5ad',
    'apply_tf': True,
    'normalize': True}
## VIASH END

## LOCAL START

parser = argparse.ArgumentParser(description="Process multiomics RNA data.")
parser.add_argument('--rna', type=str, help='Path to the multiomics RNA file')
parser.add_argument('--prediction', type=str, help='Path to the prediction file')
parser.add_argument('--resources_dir', type=str, help='Path to the prediction file')
parser.add_argument('--tf_all', type=str, help='Path to the tf_all')
parser.add_argument('--num_workers', type=str, help='Number of cores')
parser.add_argument('--max_n_links', type=str, help='Number of top links to retain')
parser.add_argument('--dataset_id', type=str, help='Dataset id')
parser.add_argument('--normalize', action='store_true')
args = parser.parse_args()

par_local = vars(args)

for key, value in par_local.items():
    if value is not None:
        par[key] = value

## LOCAL END

try:
    sys.path.append(meta["resources_dir"])
except:
    meta = {
        "resources_dir": 'src/utils'
    }
    sys.path.append(meta["resources_dir"])
from util import corr_net


if __name__ == '__main__':
    dataset_id = ad.read_h5ad(par['rna'], backed='r').uns['dataset_id']
    net = corr_net(par)

    print('Output GRN')
    print('Shape of the network:', net.shape)
    print(net.sort_values('weight', ascending=False, key=abs).head(10))
    net = net.astype(str)
    output = ad.AnnData(X=None, uns={"method_id": 'pearson_corr', "dataset_id": dataset_id, "prediction": net[["source", "target", "weight"]]})
    output.write(par['prediction'])
