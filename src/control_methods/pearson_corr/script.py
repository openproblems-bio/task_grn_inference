
import anndata as ad 
import pandas as pd
import os
import scanpy as sc 
## VIASH START
par = {
    'multiomics_rna': 'resources/grn-benchmark/multiomics_rna.h5ad',
    'tf_all': 'resources/prior/tf_all.csv',
    'cell_type_specific': False,
    'max_n_links': 50000,
    'prediction': 'output/pearson_net.csv',
    "seed": 32,
    'normalize': False,
    'donor_specific': False,
    'temp_dir': 'output/pearson_corr',
    'causal': True,
    'normalize': True}
## VIASH END


import argparse
parser = argparse.ArgumentParser(description="Process multiomics RNA data.")
parser.add_argument('--multiomics_rna', type=str, help='Path to the multiomics RNA file')
parser.add_argument('--prediction', type=str, help='Path to the prediction file')
parser.add_argument('--tf_all', type=str, help='Path to the tf_all')
parser.add_argument('--num_workers', type=str, help='Number of cores')
parser.add_argument('--max_n_links', type=str, help='Number of top links to retain')
parser.add_argument('--causal', action='store_true', help='Enable causal mode')
parser.add_argument('--normalize', action='store_true')

args = parser.parse_args()

if args.multiomics_rna:
    par['multiomics_rna'] = args.multiomics_rna
if args.causal:
    par['causal'] = True
else:
    par['causal'] = False

if args.causal:
    par['normalize'] = True
else:
    par['normalize'] = False

if args.prediction:
    par['prediction'] = args.prediction
if args.tf_all:
    par['tf_all'] = args.tf_all
if args.num_workers:
    par['num_workers'] = args.num_workers
if args.max_n_links:
    par['max_n_links'] = int(args.max_n_links)

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
    adata = ad.read_h5ad(par["multiomics_rna"])
    if 'normalize' in par:
        if par['normalize']:
            # - lognorm 
            sc.pp.normalize_total(adata)
            sc.pp.log1p(adata)
    # - corr
    gene_names = adata.var_names.to_numpy()
    grn = corr_net(adata.X, gene_names, par)    
    return grn

if par['donor_specific']:
    adata = ad.read_h5ad(par['multiomics_rna'])
    # - new dir for donor specific adata
    par['multiomics_rna'] = f"{par['temp_dir']}/multiomics_rna.h5ad"
    donor_ids = adata.obs.donor_id.unique()
    for i, donor_id in enumerate(donor_ids): # run for each donor and concat
        adata_sub = adata[adata.obs.donor_id.eq(donor_id), :]
        adata_sub.write(par['multiomics_rna'])
        net_sub = create_corr_net(par)
        net_sub['donor_id'] = donor_id
        if i == 0:
            net = net_sub
        else:
            net = pd.concat([net, net_sub])
else:
    net = create_corr_net(par)

print('Output GRN')
net.to_csv(par['prediction'])
