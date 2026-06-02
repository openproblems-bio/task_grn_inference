
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
  gene_names = adata_rna.var_names

  # Subsample cells if too many — pyscenic grn scales poorly beyond 25k cells
  MAX_CELLS = 25000
  if adata_rna.n_obs > MAX_CELLS:
      np.random.seed(42)
      idx = np.random.choice(adata_rna.n_obs, MAX_CELLS, replace=False)
      adata_rna = adata_rna[idx]
      print(f"Subsampled to {MAX_CELLS} cells for grnboost/scenic")

  from util import manage_layer
  layer = manage_layer(adata_rna, par)
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
  try:
    result = subprocess.run(command, check=True, capture_output=True, text=True)
    print(result.stdout)
  except subprocess.CalledProcessError as e:
    print(f"Error running pyscenic grn command")
    print(f"Command: {' '.join(command)}")
    print(f"Exit code: {e.returncode}")
    print(f"STDOUT:\n{e.stdout}")
    print(f"STDERR:\n{e.stderr}")
    raise RuntimeError(f"pyscenic grn failed with exit code {e.returncode}. See error details above.") from e