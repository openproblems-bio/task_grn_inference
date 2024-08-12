import os

import anndata
import numpy as np
import pandas as pd
from arboreto.algo import genie3
from distributed import Client


## VIASH START
par = {
  'multiomics_rna': 'resources/resources_test/grn-benchmark/multiomics_rna.h5ad',
  'prediction': 'output/genie3/prediction.csv',
  'max_n_links': 50000
}
## VIASH END


# Load scRNA-seq data
adata_rna = anndata.read_h5ad(par['multiomics_rna'])
gene_names = adata_rna.var.gene_ids.index.to_numpy()
X = adata_rna.X

# GRN inference
client = Client(processes=False)
network = genie3(X, client_or_address=client, gene_names=gene_names)

# Keep only top links
network = network.head(par['max_n_links'])

# Rename columns and index
network.rename(columns={'TF': 'source', 'target': 'target', 'importance': 'weight'}, inplace=True)
network.reset_index(drop=True, inplace=True)

# Save inferred GRN
network.to_csv(par['prediction'], sep=',')

print('Finished.')
