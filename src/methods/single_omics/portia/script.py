import os
import anndata
import numpy as np
import pandas as pd
import scipy.sparse
import portia as pt
from tqdm import tqdm
import anndata as ad
import sys
import argparse


## VIASH START
par = {
  'rna': 'resources/grn-benchmark/rna.h5ad',
  'tf_all': 'resources/grn_benchmark/prior/tf_all.csv',
  'prediction': 'output/portia.csv',
  'max_n_links': 50000,
  'donor_specific': False,
  'temp_dir': 'output/portia'
}
## VIASH END

## LOCAL START
parser = argparse.ArgumentParser(description="Process multiomics RNA data.")
parser.add_argument('--rna', type=str, help='Path to the multiomics RNA file')
parser.add_argument('--prediction', type=str, help='Path to the prediction file')
parser.add_argument('--resources_dir', type=str, help='Path to the prediction file')
parser.add_argument('--tf_all', type=str, help='Path to the tf_all')
parser.add_argument('--num_workers', type=str, help='Number of cores')
parser.add_argument('--max_n_links', type=str, help='Number of top links to retain')
parser.add_argument('--dataset_id', type=str, help='Dataset id')
parser.add_argument('--normalize', action='store_true')
args = parser.parse_args()

par_local = vars(args)

for key, value in par_local.items():
    if value is not None:
        par[key] = value

## LOCAL END

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
  adata_rna = anndata.read_h5ad(par['rna'])

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
    adata = ad.read_h5ad(par['rna'])

    net = main(par)

    print('Output GRN')
    net['weight'] = net['weight'].astype(str)
    output = ad.AnnData(X=None, uns={"method_id": 'portia', "dataset_id": par['dataset_id'], "prediction": net[["source", "target", "weight"]]})
    output.write(par['prediction'])
