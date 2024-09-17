import os

import anndata
import numpy as np
import pandas as pd
import scipy.sparse
import portia as pt
from tqdm import tqdm


## VIASH START
par = {
  'multiomics_rna': 'resources/grn-benchmark/multiomics_rna_0.h5ad',
  'tf_all': 'resources/prior/tf_all.csv',
  'prediction': 'output/portia_celltype_0.csv',
  'max_n_links': 50000,
  'cell_type_specific': True
}
## VIASH END


# Load scRNA-seq data
print('Reading data')
adata_rna = anndata.read_h5ad(par['multiomics_rna'])
gene_names = adata_rna.var.gene_ids.index.to_numpy()
X = adata_rna.X.toarray() if scipy.sparse.issparse(adata_rna.X) else adata_rna.X

def process_links(net, par):
  net = net[net.source!=net.target]
  net_sorted = net.reindex(net['weight'].abs().sort_values(ascending=False).index)
  net = net_sorted.head(par['max_n_links']).reset_index(drop=True)
  return net

# Load list of putative TFs
tfs = np.loadtxt(par['tf_all'], dtype=str)
tf_names = [gene_name for gene_name in gene_names if (gene_name in tfs)]

# GRN inference
def infer_grn(X, gene_names):
  print('Inferring grn')
  dataset = pt.GeneExpressionDataset()
  for exp_id, data in enumerate(X):
      dataset.add(pt.Experiment(exp_id, data))
  M_bar = pt.run(dataset, method='no-transform')

  ranked_scores = pt.rank_scores(M_bar, gene_names)
  sources, targets, weights = zip(*[(gene_a, gene_b, score) for gene_a, gene_b, score in ranked_scores])

  grn = pd.DataFrame({'source':sources, 'target':targets, 'weight':weights})
  print(grn.shape)
  grn = grn[grn.source.isin(tf_names)]

  grn = process_links(grn, par)
  return grn

if par['cell_type_specific']:
  groups = adata_rna.obs.cell_type
  i = 0
  for group in tqdm(np.unique(groups), desc="Processing groups"):
    X_sub = X[groups == group, :]
    
    net = infer_grn(X_sub, gene_names)
    net['cell_type'] = group
    if i==0:
        grn = net
    else:
        grn = pd.concat([grn, net], axis=0).reset_index(drop=True)

    i += 1
        
else:
  grn = infer_grn(X, gene_names)

grn.to_csv(par['prediction'])

# # Save inferred GRN
# with open(par['prediction'], 'w') as f:
#     f.write(f',source,target,weight\n')
#     for i, (gene_a, gene_b, score) in enumerate(pt.rank_scores(M_bar, gene_names, limit=par['max_n_links'])):
#         f.write(f'{i},{gene_a},{gene_b},{score}\n')

print('Finished.')
