import os

import anndata
import numpy as np
import pandas as pd
import scipy.sparse
import portia as pt
from tqdm import tqdm


## VIASH START
par = {
  'multiomics_rna': 'resources/grn-benchmark/multiomics_rna_d0_hvg.h5ad',
  'tf_all': 'resources/prior/tf_all.csv',
  'prediction': 'output/portia_donor_0_hvgs.csv',
  'max_n_links': 50000,
  'cell_type_specific': False,
  'normalize': False,
  'only_hvgs': True
}
## VIASH END

import sys
meta= {
  "resources_dir": 'src/utils/'
}

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
    
sys.path.append(meta["resources_dir"])
from util import process_links
# Load scRNA-seq data
print('Reading data')
adata_rna = anndata.read_h5ad(par['multiomics_rna'])

gene_names = adata_rna.var_names
X = adata_rna.X.toarray() if scipy.sparse.issparse(adata_rna.X) else adata_rna.X

# Load list of putative TFs
tfs = np.loadtxt(par['tf_all'], dtype=str)
tf_names = [gene_name for gene_name in gene_names if (gene_name in tfs)]
tf_idx = np.asarray([i for i, gene_name in enumerate(gene_names) if gene_name in tf_names], dtype=int)


# GRN inference
def infer_grn(X, gene_names):
  print('Inferring grn')
  dataset = pt.GeneExpressionDataset()
  
  for exp_id, data in enumerate(X):
      dataset.add(pt.Experiment(exp_id, data))
  
  M_bar = pt.run(dataset, tf_idx=tf_idx, method='no-transform')
  ranked_scores = pt.rank_scores(M_bar, gene_names, limit=par['max_n_links'])
  sources, targets, weights = zip(*[(gene_a, gene_b, score) for gene_a, gene_b, score in ranked_scores])

  grn = pd.DataFrame({'source':sources, 'target':targets, 'weight':weights})
  print(grn.shape)
  grn = grn[grn.source.isin(tf_names)]

  grn = process_links(grn, par)
  return grn
# par['cell_type_specific'] = False
grn = infer_grn(X, gene_names)

grn.to_csv(par['prediction'])

print('Finished.')
