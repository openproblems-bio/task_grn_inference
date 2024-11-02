
import anndata as ad 
import pandas as pd
import os
import scanpy as sc 
## VIASH START
par = {
    'rna': 'resources/grn-benchmark/rna.h5ad',
    'tf_all': 'resources/prior/tf_all.csv',
    'cell_type_specific': False,
    'max_n_links': 50000,
    'prediction': 'output/pearson_net.csv',
    "seed": 32,
    'normalize': False,
    'donor_specific': False,
    'temp_dir': 'output/pearson_corr',
    'apply_tf': True,
    'normalize': True}
## VIASH END


import argparse
parser = argparse.ArgumentParser(description="Process multiomics RNA data.")
parser.add_argument('--rna', type=str, help='Path to the multiomics RNA file')
parser.add_argument('--prediction', type=str, help='Path to the prediction file')
parser.add_argument('--tf_all', type=str, help='Path to the tf_all')
parser.add_argument('--num_workers', type=str, help='Number of cores')
parser.add_argument('--normalize', action='store_true')
parser.add_argument('--max_n_links', type=str, help='Number of top links to retain')

args = parser.parse_args()

if args.max_n_links:
    par['max_n_links'] = int(args.max_n_links)
    
if args.rna:
    par['rna'] = args.rna

if args.normalize:
    par['normalize'] = True
else:
    par['normalize'] = False

if args.prediction:
    par['prediction'] = args.prediction
if args.tf_all:
    par['tf_all'] = args.tf_all
if args.num_workers:
    par['num_workers'] = args.num_workers

os.makedirs(par['temp_dir'], exist_ok=True)
import sys

try:
    sys.path.append(meta["resources_dir"])
except:
    meta = {
    'resources_dir': 'src/utils'
    }
    sys.path.append(meta["resources_dir"])
from util import corr_net


def create_corr_net(par):
    print(par)
    print('Read data')
    adata = ad.read_h5ad(par["rna"])
    if 'normalize' in par:
        if par['normalize']:
            print('normalizing')
            # - lognorm 
            sc.pp.normalize_total(adata)
            sc.pp.log1p(adata)
    if 'is_test' in adata.obs.columns:
        adata = adata[~adata.obs.is_test] # train on only non test data
    # - corr
    gene_names = adata.var_names.to_numpy()
    grn = corr_net(adata.X, gene_names, par)    
    return grn

net = create_corr_net(par)

print('Output GRN')
net.to_csv(par['prediction'])
