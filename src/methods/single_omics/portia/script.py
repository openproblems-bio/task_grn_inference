import os

import anndata
import numpy as np
import pandas as pd
import scipy.sparse
import portia as pt


## VIASH START
par = {
  'multiomics_rna': 'resources/resources_test/grn-benchmark/multiomics_rna.h5ad',
  'tfs': 'resources/prior/tf_all.csv',
  'prediction': 'output/portia/prediction.csv',
  'max_n_links': 50000
}
## VIASH END


# Load scRNA-seq data
adata_rna = anndata.read_h5ad(par['multiomics_rna'])
gene_names = adata_rna.var.gene_ids.index.to_numpy()
X = adata_rna.X.toarray() if scipy.sparse.issparse(adata_rna.X) else adata_rna.X

# Remove genes with >=90% of zeros
mask = (np.mean(X == 0, axis=0) >= 0.9)
X = X[:, ~mask]
gene_names = gene_names[~mask]

# Remove samples with >=90% of zeros
mask = (np.mean(X == 0, axis=1) >= 0.9)
adata_rna = X[~mask, :]

# Load list of putative TFs
df = pd.read_csv(par['tfs'], header=None, names=['gene_name'])
tfs = set(list(df['gene_name']))
tf_names = [gene_name for gene_name in gene_names if (gene_name in tfs)]

# GRN inference
dataset = pt.GeneExpressionDataset()
for exp_id, data in enumerate(X):
    dataset.add(pt.Experiment(exp_id, data))
M_bar = pt.run(dataset, method='no-transform')

print(M_bar.shape)

# Save inferred GRN
with open(par['prediction'], 'w') as f:
    f.write(f',source,target,weight\n')
    for i, (gene_a, gene_b, score) in enumerate(pt.rank_scores(M_bar, gene_names, limit=par['max_n_links'])):
        f.write(f'{i},{gene_a},{gene_b},{score}\n')

print('Finished.')
