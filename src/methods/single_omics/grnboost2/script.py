import os

import anndata
import numpy as np
import pandas as pd
from arboreto.algo import grnboost2
from distributed import Client, LocalCluster
from tqdm import tqdm
import subprocess 
import sys

## VIASH START
par = {
  'rna': 'resources/grn-benchmark/rna_d0_hvg.h5ad',
  "tf_all": 'resources/prior/tf_all.csv',
  'prediction': 'output/grnboost2_donor_0_hvg.csv',
  'max_n_links': 50000,
  'cell_type_specific': False,
  'normalize': False,
  'qc': False
}
## VIASH END
import argparse
parser = argparse.ArgumentParser(description="Process multiomics RNA data.")
parser.add_argument('--rna', type=str, help='Path to the multiomics RNA file')
parser.add_argument('--prediction', type=str, help='Path to the prediction file')
parser.add_argument('--resources_dir', type=str, help='Path to the prediction file')
parser.add_argument('--tf_all', type=str, help='Path to the tf_all')
parser.add_argument('--num_workers', type=str, help='Number of cores')
parser.add_argument('--qc', action='store_true', help='Whether to do QC or not')
parser.add_argument('--max_n_links', type=str, help='Number of links')

args = parser.parse_args()

if args.max_n_links:
    par['max_n_links'] = int(args.max_n_links)
if args.rna:
    par['rna'] = args.rna
if args.prediction:
    par['prediction'] = args.prediction
if args.tf_all:
    par['tf_all'] = args.tf_all
if args.num_workers:
    par['num_workers'] = args.num_workers
if args.qc:
    par['qc'] = args.qc
    
if args.resources_dir:
    meta['resources_dir'] = args.resources_dir   

print(par)

try:
    sys.path.append(meta["resources_dir"])
except:
    meta= {
    "resources_dir": 'src/utils/'
    }
    sys.path.append(meta["resources_dir"])
from util import process_links, basic_qc
# Load scRNA-seq data
print('Reading data')
adata_rna = anndata.read_h5ad(par['rna'])
if 'qc' in par:
    if par['qc']:
        print('Shape before QC: ', adata_rna.shape)
        adata_rna = basic_qc(adata_rna)
        print('Shape after QC: ', adata_rna.shape)

gene_names = adata_rna.var_names
X = adata_rna.X

# Load list of putative TFs
tfs = np.loadtxt(par["tf_all"], dtype=str)
tf_names = [gene_name for gene_name in gene_names if (gene_name in tfs)]

# GRN inference
client = Client(processes=False)

def infer_grn(X, par):
  print("Infer grn", flush=True)
  
  network = grnboost2(X, client_or_address=client, gene_names=gene_names, tf_names=tf_names)
  network.rename(columns={'TF': 'source', 'target': 'target', 'importance': 'weight'}, inplace=True)
  network.reset_index(drop=True, inplace=True)
  network = process_links(network, par)
  
  return network

grn = infer_grn(X, par)       

# Save inferred GRN
grn.to_csv(par['prediction'], sep=',')

print('Finished.')
