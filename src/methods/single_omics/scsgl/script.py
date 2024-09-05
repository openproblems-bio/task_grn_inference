import os
import sys
import contextlib
import subprocess

import anndata
import numpy as np
import scipy.sparse
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

# Download scSGL
SCSGL_FOLDER = os.path.join(par['temp_dir'], 'scSGL')
os.makedirs(par['temp_dir'], exist_ok=True)
if not os.path.exists(os.path.join(par['temp_dir'], 'scSGL', 'pysrc')):
    subprocess.run(['git', 'clone', 'https://github.com/Single-Cell-Graph-Learning/scSGL', SCSGL_FOLDER])

# Import pysrc locally (from the cloned repository)
sys.path.append(SCSGL_FOLDER)
from pysrc.graphlearning import learn_signed_graph
from pysrc.evaluation import auc

# Run scSGL
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

# Sort values
df = df.sort_values(by='weight', key=abs, ascending=False)

# Keep only top links
df = df.head(par['max_n_links'])

# Reset index
df.reset_index(drop=True, inplace=True)

# Save inferred GRN
df.to_csv(par['prediction'], sep=',')

print('Finished.')
