# !pip install sctk anndata
# !aws s3 cp s3://openproblems-bio/public/neurips-2023-competition/sc_counts.h5ad  ./resources_raw/ --no-sign-request

import anndata as ad 
import pandas as pd
import numpy as np
from scipy import sparse
import scanpy as sc
import requests

## VIASH START

par = {
    'multiomics_rna': 'resources/grn-benchmark/multiomics_rna.h5ad',
    'multiomics_rna_test': 'resources_test/grn-benchmark/multiomics_rna.h5ad',

    'multiomics_atac': 'resources/grn-benchmark/multiomics_atac.h5ad',
    'multiomics_atac_test': 'resources_test/grn-benchmark/multiomics_atac.h5ad',

    'perturbation_data': 'resources/grn-benchmark/perturbation_data.h5ad',
    'perturbation_data_test': 'resources_test/grn-benchmark/perturbation_data.h5ad',
}
## VIASH END

n_cells = 2000
n_peaks = 10000

adata_rna = ad.read_h5ad(par['multiomics_rna'])
adata_atac = ad.read_h5ad(par['multiomics_atac'])

# shorten rna 
mask = adata_rna.obs.donor_id=='donor_0'
adata_rna_s = adata_rna[mask]
# only one chr and n_cells 
if 'interval' not in adata_rna.var:
    raise ValueError('location is not given in rna')
chr_mask = adata_rna_s.var.interval.str.split(':').str[0] == 'chr1'
adata_rna_s = adata_rna_s[:n_cells, chr_mask]

# shorten atac
mask = adata_atac.obs.donor_id=='donor_0'
adata_atac_s = adata_atac[mask]
chr_mask = adata_atac_s.var.index.str.split(':').str[0] == 'chr1'
adata_atac_s = adata_atac_s[adata_rna_s.obs_names, chr_mask]

# adata_atac_s = adata_atac_s[:, adata_atac_s.var.sample(n=n_peaks, random_state=42).index]
total_counts = adata_atac_s.X.sum(axis=0).A1  # .A1 converts the sparse matrix to a dense array

# Create a DataFrame with peak indices and their total counts
peaks_df = pd.DataFrame({
    'peak': adata_atac_s.var.index,
    'total_counts': total_counts
})

# Sort peaks by total counts in descending order and select the top n_peaks
top_peaks = peaks_df.nlargest(n_peaks, 'total_counts')['peak']

# Subset the ATAC data to keep only the top peaks
adata_atac_s = adata_atac_s[:, top_peaks]

# Ensure the ATAC variable information is updated
adata_atac_s.var = adata_atac_s.var.loc[top_peaks]

print(adata_rna_s)
print(adata_atac_s)
adata_rna_s.write(par['multiomics_rna_test'])
adata_atac_s.write(par['multiomics_atac_test'])

# shorten perturbation
adata_bulk = ad.read_h5ad(par['perturbation_data'])
adata_bulk[:200, adata_bulk.var_names.isin(adata_rna_s.var_names)].write(par['perturbation_data_test'])