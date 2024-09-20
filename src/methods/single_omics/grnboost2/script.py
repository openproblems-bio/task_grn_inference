import os

import anndata
import numpy as np
import pandas as pd
from arboreto.algo import grnboost2
from distributed import Client, LocalCluster
from tqdm import tqdm
import subprocess 



## VIASH START
par = {
  'multiomics_rna': 'resources/grn-benchmark/multiomics_rna_d0_hvg.h5ad',
  "tf_all": 'resources/prior/tf_all.csv',
  'prediction': 'output/grnboost2_donor_0_hvg.csv',
  'max_n_links': 50000,
  'cell_type_specific': False,
  'normalize': False
}
## VIASH END

import sys
meta= {
  "resources_dir": 'src/utils/'
}
sys.path.append(meta["resources_dir"])
from util import process_links
# Load scRNA-seq data
print('Reading data')
adata_rna = anndata.read_h5ad(par['multiomics_rna'])

groups = adata_rna.obs.cell_type
gene_names = adata_rna.var.gene_ids.index.to_numpy()
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

# par['cell_type_specific'] = False
if par['cell_type_specific']:
    i = 0
    for group in tqdm(np.unique(groups), desc="Processing groups"):
        X_sub = X[groups == group, :]
        net = infer_grn(X_sub, par)
        net['cell_type'] = group
        if i==0:
            grn = net
        else:
            grn = pd.concat([grn, net], axis=0).reset_index(drop=True)
        i += 1
else:
    grn = infer_grn(X, par)       

# Save inferred GRN
grn.to_csv(par['prediction'], sep=',')

print('Finished.')
