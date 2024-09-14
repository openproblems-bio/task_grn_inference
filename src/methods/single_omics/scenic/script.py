import os

import anndata
import numpy as np
import pandas as pd
from arboreto.algo import grnboost2
from distributed import Client


## VIASH START
par = {
  'multiomics_rna': 'resources_test/grn-benchmark/multiomics_rna.h5ad',
  "tf_all": 'resources/prior/tf_all.csv',
  'prediction': 'output/grnboost2/prediction.csv',
  'max_n_links': 50000
}
## VIASH END
os.makedirs(par['temp_dir'], exist_ok=True)

# Load scRNA-seq data
adata_rna = anndata.read_h5ad(par['multiomics_rna'])
gene_names = adata_rna.var.gene_ids.index.to_numpy()
X = adata_rna.X.toarray()

# Load list of putative TFs
# df = pd.read_csv(par["tf_all"], header=None, names=['gene_name'])
# tfs = set(list(df['gene_name']))
# tf_names = [gene_name for gene_name in gene_names if (gene_name in tfs)]

# format output
expression_data = f"{par['temp_dir']}/expression_data.tsv"
pd.DataFrame(X, columns=gene_names).to_csv(expression_data, sep='\t', index=False)

expr_mat_adjacencies = f"{par['temp_dir']}/expr_mat_adjacencies.tsv"
command = [
    "pyscenic", "grn",
    "--num_workers", par['max_workers'],
    "-o", expr_mat_adjacencies,
    expression_data,
    par['tf_all']
]

# Run grn
import subprocess
subprocess.run(command, check=True)


# Run prune
regulons = f"{par['temp_dir']}/regulons.csv"
annotations_fname = "/data/motifs-v9-nr.hgnc-m0.001-o0.0.tbl" 
ranking_1 = "/data/hg19-tss-centered-5kb-7species.mc9nr.genes_vs_motifs.rankings.feather "
ranking_2 = /data/hg19-tss-centered-10kb-7species.mc9nr.genes_vs_motifs.rankings.feather
command = [
    "pyscenic", "ctx",
    expr_mat_adjacencies, ranking_1, ranking_2,
    "--annotations_fname", annotations_fname, 
    "--expression_mtx_fname", expression_data,
    "--mode", "custom_multiprocessing",
    "--output", regulons, 
    "--num_workers", par['max_workers']
]
subprocess.run(command, check=True)

# Save inferred GRN
print(expr_mat_adjacencies)
network = pd.read_csv(expr_mat_adjacencies,  sep='\t')
network.to_csv(par['prediction'], sep=',')

print('Finished.')


