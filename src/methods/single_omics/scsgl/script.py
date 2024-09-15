import os
import sys
import contextlib
import subprocess
from itertools import permutations

import anndata
import numpy as np
import scipy.sparse
from scipy.spatial.distance import squareform
import pandas as pd


## VIASH START
par = {
  'multiomics_rna': 'resources/resources_test/grn-benchmark/multiomics_rna.h5ad',
  'prediction': 'output/scsgl/prediction.csv',
  'temp_dir': 'output/scsgl',
  'max_n_links': 50000
}
## VIASH END

# Load scRNA-seq data
adata_rna = anndata.read_h5ad(par['multiomics_rna'])
gene_names = adata_rna.var.gene_ids.index.to_numpy()
X = adata_rna.X.toarray() if scipy.sparse.issparse(adata_rna.X) else adata_rna.X

# Import pysrc locally
sys.path.append(os.environ['SCSGL_PATH'])
from pysrc.graphlearning import learn_signed_graph

# Remove genes with >=90% of zeros
# mask = (np.mean(X == 0, axis=0) >= 0.75)
# X = X[:, ~mask]
# gene_names = gene_names[~mask]

# # Remove samples with >=90% of zeros
# mask = (np.mean(X == 0, axis=1) >= 0.75)
# X = X[~mask, :]

# Run scSGL
print('Starting scSGL', flush=True)
df = learn_signed_graph(
    X.T,
    pos_density=0.45,
    neg_density=0.45,
    assoc='correlation',
    gene_names=gene_names
)

# Rename columns and index
df.rename(columns={'Gene1': 'source', 'Gene2': 'target', 'EdgeWeight': 'weight'}, inplace=True)
df.reset_index(inplace=True)
df.drop('index', axis=1, inplace=True, errors='ignore')

# Load list of putative TFs
df_tfs = pd.read_csv(par['tf_all'], header=None, names=['gene_name'])
tfs = set(list(df_tfs['gene_name']))

# Ensure first gene in pair is a putative TF
print(df)
mask = np.asarray([(gene_name in tfs) for gene_name in df['source']], dtype=bool)
df = df[mask]

# Sort values
df = df.sort_values(by='weight', key=abs, ascending=False)

# Keep only top links
df = df.head(par['max_n_links'])

# Reset index
df.reset_index(drop=True, inplace=True)

# Save inferred GRN
df.to_csv(par['prediction'], sep=',')

print('Finished.')
