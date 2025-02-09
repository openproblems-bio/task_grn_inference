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
import argparse

## VIASH START
par = {
    'rna': 'resources/grn-benchmark/rna.h5ad',
    'tf_all': 'resources/grn_benchmark/prior/tf_all.csv',
    'max_n_links': 50000,
    'prediction': 'resources/grn_models/positive_control.csv',
    "seed": 32,
    'temp_dir': 'output/positive_control',
    }
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

    print('Output GRN')
    net['weight'] = net['weight'].astype(str)
    output = ad.AnnData(X=None, uns={"method_id": 'positive_control', "dataset_id": par['dataset_id'], "prediction": net[["source", "target", "weight"]]})
    output.write(par['prediction'])
