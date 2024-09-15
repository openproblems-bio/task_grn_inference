import os

import anndata
import numpy as np
import pandas as pd
from arboreto.algo import grnboost2
from distributed import Client


## VIASH START
par = {
  'multiomics_rna': 'resources/resources_test/grn-benchmark/multiomics_rna.h5ad',
  "tf_all": 'resources/prior/tf_all.csv',
  'prediction': 'output/grnboost2/prediction.csv',
  'max_n_links': 50000
}
## VIASH END


# Load scRNA-seq data
adata_rna = anndata.read_h5ad(par['multiomics_rna'])
gene_names = adata_rna.var.gene_ids.index.to_numpy()
X = adata_rna.X.toarray()

if False:
  # Remove genes with >=90% of zeros
  mask = (np.mean(X == 0, axis=0) >= 0.9)
  X = X[:, ~mask]
  gene_names = gene_names[~mask]

  # Remove samples with >=90% of zeros
  mask = (np.mean(X == 0, axis=1) >= 0.9)
  adata_rna = X[~mask, :]

# Load list of putative TFs
df = pd.read_csv(par["tf_all"], header=None, names=['gene_name'])
tfs = set(list(df['gene_name']))
tf_names = [gene_name for gene_name in gene_names if (gene_name in tfs)]

# GRN inference
client = Client(processes=False)
network = grnboost2(X, client_or_address=client, gene_names=gene_names, tf_names=tf_names)

# Keep only top links
network = network.head(par['max_n_links'])

# Rename columns and index
network.rename(columns={'TF': 'source', 'target': 'target', 'importance': 'weight'}, inplace=True)
network.reset_index(drop=True, inplace=True)

# Save inferred GRN
network.to_csv(par['prediction'], sep=',')

print('Finished.')
