
import anndata
import numpy as np
import pandas as pd
import subprocess
import scipy.sparse as sp
import sys
import os
import anndata as ad

def format_data(par):
  print('Read data')
  adata_rna = anndata.read_h5ad(par['rna'])  
  if False:
    if 'perturbation' in adata_rna.obs.columns: #test
      print('Subsample data... Testing:')
      pertubs = adata_rna.obs['perturbation'].unique()[:10]
      adata_rna = adata_rna[adata_rna.obs['perturbation'].isin(pertubs)]
      print(adata_rna.shape)
  gene_names = adata_rna.var_names
  
  layer = 'lognorm' if 'lognorm' in adata_rna.layers else 'X_norm'
  X = adata_rna.layers[layer]

  if sp.issparse(X):
    X = X.toarray()
  pd.DataFrame(X, columns=gene_names).to_csv(par['expression_data'], sep='\t', index=False)
  
def run_grn(par):
  print('Run grn')
  command = [
      "pyscenic", "grn",
      "--num_workers", str(par['num_workers']),
      "--seed", str(par['seed']),
      "-o", str(par['expr_mat_adjacencies']),
      "--method", "grnboost2", 
      str(par['expression_data']),
      par['tf_all']
  ]
  subprocess.run(command, check=True)