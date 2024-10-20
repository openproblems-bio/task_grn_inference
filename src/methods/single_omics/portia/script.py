import os

import anndata
import numpy as np
import pandas as pd
import scipy.sparse
import portia as pt
from tqdm import tqdm
import anndata as ad


## VIASH START
par = {
  'multiomics_rna': 'resources/grn-benchmark/multiomics_rna.h5ad',
  'tf_all': 'resources/prior/tf_all.csv',
  'prediction': 'output/portia.csv',
  'max_n_links': 50000,
  'donor_specific': False,
  'temp_dir': 'output/portia'
}
## VIASH END

import sys


import argparse
parser = argparse.ArgumentParser(description="Process multiomics RNA data.")
parser.add_argument('--multiomics_rna', type=str, help='Path to the multiomics RNA file')
parser.add_argument('--prediction', type=str, help='Path to the prediction file')
parser.add_argument('--tf_all', type=str, help='Path to the tf_all')
parser.add_argument('--num_workers', type=str, help='Number of cores')
parser.add_argument('--max_n_links', type=str, help='Number of top links to retain')
parser.add_argument('--causal', action='store_true', help='Enable causal mode')
parser.add_argument('--normalize', action='store_true')

args = parser.parse_args()

if args.multiomics_rna:
    par['multiomics_rna'] = args.multiomics_rna
if args.causal:
    par['causal'] = True
else:
    par['causal'] = False

if args.causal:
    par['normalize'] = True
else:
    par['normalize'] = False

if args.prediction:
    par['prediction'] = args.prediction
if args.tf_all:
    par['tf_all'] = args.tf_all
if args.num_workers:
    par['num_workers'] = args.num_workers
if args.max_n_links:
    par['max_n_links'] = int(args.max_n_links)

os.makedirs(par['temp_dir'], exist_ok=True)
import sys

try:
    sys.path.append(meta["resources_dir"])
except:
    meta = {
    'resources_dir': 'src/utils'
    }
    sys.path.append(meta["resources_dir"])
from util import process_links


def main(par):
  print('Reading data')
  adata_rna = anndata.read_h5ad(par['multiomics_rna'])

  gene_names = adata_rna.var_names
  X = adata_rna.X.toarray() if scipy.sparse.issparse(adata_rna.X) else adata_rna.X

  # Load list of putative TFs
  tfs = np.loadtxt(par['tf_all'], dtype=str)
  tf_names = [gene_name for gene_name in gene_names if (gene_name in tfs)]
  tf_idx = np.asarray([i for i, gene_name in enumerate(gene_names) if gene_name in tf_names], dtype=int)


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


if __name__ == '__main__':
  os.makedirs(par['temp_dir'], exist_ok=True)
  adata = ad.read_h5ad(par['multiomics_rna'])
  if ('donor_id' in adata.obs) & (par['donor_specific']):
      # - new dir for donor specific adata
      par['multiomics_rna'] = f"{par['temp_dir']}/multiomics_rna.h5ad"
      donor_ids = adata.obs.donor_id.unique()
      for i, donor_id in enumerate(donor_ids): # run for each donor and concat
          print('GRN inference for ', donor_id)
          adata_sub = adata[adata.obs.donor_id.eq(donor_id), :]
          adata_sub.write(par['multiomics_rna'])
          net_sub = main(par)
          net_sub['donor_id'] = donor_id
          if i == 0:
              net = net_sub
          else:
              net = pd.concat([net, net_sub])
  else:
      net = main(par)

  net.to_csv(par['prediction'])

  print('Finished.')
