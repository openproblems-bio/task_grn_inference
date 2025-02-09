
import anndata as ad 
import pandas as pd
import os
import scanpy as sc 
import sys
import numpy as np
import argparse


## VIASH START
par = {
    'rna': 'resources/grn_benchmark/evaluation_data//op_rna.h5ad',
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
    net = corr_net(par)
    # - format of et
    '''
        the net is a pandas dataframe with the following columns:
            - source: the source gene of the interaction
            - target: the target gene of the interaction
            - weight: the weight of the interaction
    '''

    print('Output GRN')
    net['weight'] = net['weight'].astype(str)
    output = ad.AnnData(X=None, uns={"method_id": 'pearson_corr', "dataset_id": par['dataset_id'], "prediction": net[["source", "target", "weight"]]})
    output.write(par['prediction'])
