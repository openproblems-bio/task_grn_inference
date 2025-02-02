import os
import pandas as pd
import numpy as np
import anndata as ad
import scanpy as sc
from tqdm import tqdm
from scipy.stats import spearmanr
from sklearn.preprocessing import StandardScaler
import scanpy as sc 


## VIASH START
par = {
    'rna': 'resources/grn-benchmark/rna.h5ad',
    'tf_all': 'resources/prior/tf_all.csv',
    'max_n_links': 50000,
    'prediction': 'resources/grn_models/positive_control.csv',
    "seed": 32,
    'temp_dir': 'output/positive_control',
    'donor_specific': False,
    'cell_type_specific': False}
## VIASH END
import sys
import argparse
parser = argparse.ArgumentParser(description="Process multiomics RNA data.")
parser.add_argument('--rna', type=str, help='Path to the multiomics RNA file')
parser.add_argument('--prediction', type=str, help='Path to the prediction file')
parser.add_argument('--resources_dir', type=str, help='Path to the prediction file')
parser.add_argument('--tf_all', type=str, help='Path to the tf_all')
parser.add_argument('--num_workers', type=str, help='Number of cores')
parser.add_argument('--max_n_links', type=str, help='Number of top links to retain')
parser.add_argument('--normalize', action='store_true')
args = parser.parse_args()
if args.rna:
    par['rna'] = args.rna
if args.prediction:
    par['prediction'] = args.prediction
if args.tf_all:
    par['tf_all'] = args.tf_all
if args.num_workers:
    par['num_workers'] = args.num_workers

if args.max_n_links:
    par['max_n_links'] = int(args.max_n_links)
if args.resources_dir:
    meta['resources_dir'] = args.resources_dir  

if args.normalize:
    par['normalize'] = True
else:
    par['normalize'] = False
try:
    sys.path.append(meta["resources_dir"])
except:
    meta = {
        "resources_dir": 'src/utils'
    }
    sys.path.append(meta["resources_dir"])
from util import corr_net

def create_corr_net(par):
    print(par)
    print('Read data')
    adata = ad.read_h5ad(par["rna"])
    X = adata.layers['X_norm']
    # - corr
    gene_names = adata.var_names.to_numpy()
    grn = corr_net(X, gene_names, par)    
    return grn

net = create_corr_net(par)

print('Output GRN')
net['weight'] = net['weight'].astype(str)
output = ad.AnnData(X=None, uns={"method_id": par['method_id'], "dataset_id": par['dataset_id'], "prediction": net[["source", "target", "weight"]]})
output.write(par['prediction'])
