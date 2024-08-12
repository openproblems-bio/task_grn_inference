import os

import anndata
import numpy as np
import scipy.sparse
import portia as pt


## VIASH START
par = {
  'multiomics_rna': 'resources/resources_test/grn-benchmark/multiomics_rna.h5ad',
  'prediction': 'output/portia/prediction.csv',
  'max_n_links': 50000
}
## VIASH END


# Load scRNA-seq data
adata_rna = anndata.read_h5ad(par['multiomics_rna'])
gene_names = adata_rna.var.gene_ids.index.to_numpy()
X = adata_rna.X.toarray() if scipy.sparse.issparse(adata_rna.X) else adata_rna.X

# GRN inference
dataset = pt.GeneExpressionDataset()
for exp_id, data in enumerate(X):
    dataset.add(pt.Experiment(exp_id, data))
M_bar = pt.run(dataset, method='no-transform')

# Save inferred GRN
with open(par['prediction'], 'w') as f:
    f.write(f',source,target,weight\n')
    for i, (gene_a, gene_b, score) in enumerate(pt.rank_scores(M_bar, gene_names, limit=par['max_n_links'])):
        f.write(f'{i},{gene_a},{gene_b},{score}\n')

print('Finished.')
