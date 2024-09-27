import os

import anndata
import numpy as np
import pandas as pd
from arboreto.algo import genie3
from distributed import Client
import scipy.sparse as sp



## VIASH START
par = {
  'multiomics_rna': 'resources/grn-benchmark/multiomics_rna_d0_hvg.h5ad',
  "tf_all": 'resources/prior/tf_all.csv',
  'prediction': 'output/genie3_donor0_hvg.csv',
  'max_n_links': 50000,
  'normalize': False
}
## VIASH END

import argparse
parser = argparse.ArgumentParser(description="Process multiomics RNA data.")
parser.add_argument('--multiomics_rna', type=str, help='Path to the multiomics RNA file')
parser.add_argument('--prediction', type=str, help='Path to the prediction file')
parser.add_argument('--resources_dir', type=str, help='Path to the prediction file')
parser.add_argument('--tf_all', type=str, help='Path to the tf_all')
parser.add_argument('--num_workers', type=str, help='Number of cores')
args = parser.parse_args()

if args.multiomics_rna:
    par['multiomics_rna'] = args.multiomics_rna
if args.prediction:
    par['prediction'] = args.prediction
if args.tf_all:
    par['tf_all'] = args.tf_all
if args.num_workers:
    par['num_workers'] = args.num_workers
    
if args.resources_dir:
    meta['resources_dir'] = args.resources_dir  
    
# Load scRNA-seq data
print('Reading data')
adata_rna = anndata.read_h5ad(par['multiomics_rna'])

gene_names = adata_rna.var.gene_ids.index.to_numpy()
if sp.issparse(adata_rna.X):
    adata_rna.X = adata_rna.X.toarray()
X = adata_rna.X

# Load list of putative TFs
df = pd.read_csv(par["tf_all"], header=None, names=['gene_name'])
tfs = set(list(df['gene_name']))
tf_names = [gene_name for gene_name in gene_names if (gene_name in tfs)]


print('GRN inference')
client = Client(processes=False)
network = genie3(X, client_or_address=client, gene_names=gene_names, tf_names=tf_names)

# Keep only top links
network = network.head(par['max_n_links'])

# Rename columns and index
network.rename(columns={'TF': 'source', 'target': 'target', 'importance': 'weight'}, inplace=True)
network.reset_index(drop=True, inplace=True)

# Save inferred GRN
network.to_csv(par['prediction'], sep=',')

print('Finished.')
