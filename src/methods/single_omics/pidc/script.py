import os
import subprocess
import pathlib

import anndata
import numpy as np
import pandas as pd


## VIASH START
par = {
    'multiomics_rna': 'resources/resources_test/grn-benchmark/multiomics_rna.h5ad',
    'prediction': 'output/pidc/prediction.csv',
    'temp_dir': 'output/pidc',
    'max_n_links': 50000
}
## VIASH END

# Load scRNA-seq data
adata_rna = anndata.read_h5ad(par['multiomics_rna'])
gene_names = adata_rna.var.gene_ids.index.to_numpy()
df_rna = pd.DataFrame(adata_rna.X.toarray(), index=adata_rna.obs.index, columns=adata_rna.var.index)

# Remove genes with >=90% of zeros
mask = (np.mean(df_rna.to_numpy() == 0, axis=0) >= 0.75)
df_rna = df_rna.loc[:, ~mask]

# Remove samples with >=90% of zeros
mask = (np.mean(df_rna.to_numpy() == 0, axis=1) >= 0.75)
df_rna = df_rna.loc[~mask, :]

# (genes x samples) format is needed
df_rna = df_rna.transpose()

# Save expression data
if not os.path.exists(os.path.join(par['temp_dir'], 'multiomics_rna.tsv')):
    os.makedirs(par['temp_dir'], exist_ok=True)
    df_rna.to_csv(os.path.join(par['temp_dir'], 'multiomics_rna.tsv'), sep='\t', index=True, header=True)

# Run PIDC
JL_SCRIPT_PATH = os.path.join(meta['resources_dir'], 'pidc.jl')
subprocess.run([
    'julia',
    JL_SCRIPT_PATH,
    os.path.join(par['temp_dir'], 'multiomics_rna.tsv'),
    os.path.join(par['temp_dir'], 'output.tsv'),
])

# Load list of putative TFs
df = pd.read_csv(par['tf_all'], header=None, names=['gene_name'])
tfs = set(list(df['gene_name']))

# Re-format output
df = pd.read_csv(os.path.join(par['temp_dir'], 'output.tsv'), header=None, names=['source', 'target', 'weight'], sep='\t')
mask = np.asarray([(gene_name in tfs) for gene_name in df['source']], dtype=bool)
df = df[mask]
df = df.head(par['max_n_links'])
df.to_csv(par['prediction'], header=True, sep=',')

print('Finished.')
