
import anndata as ad 
import pandas as pd
import numpy as np
import scanpy as sc
import os

par = {
    'multiomics_counts': 'resources/datasets_raw/multiome_counts.h5ad',
    'multiomics_counts_test': 'resources_test/datasets_raw/multiome_counts.h5ad',

    'perturbation_counts': 'resources/datasets_raw/perturbation_counts.h5ad',
    'perturbation_counts_test': 'resources_test/datasets_raw/perturbation_counts.h5ad',

}

os.makedirs('resources_test/datasets_raw/', exist_ok=True)
multiomics_counts = ad.read_h5ad(par['multiomics_counts'])
# shorten multiomics_counts 
mask = multiomics_counts.obs.donor_id=='CCL-240'
multiomics_counts_s = multiomics_counts[mask]
# only one chr and n_cells 
if 'interval' not in multiomics_counts_s.var:
    raise ValueError('location is not given in rna')
chr_mask = multiomics_counts_s.var.interval.str.split(':').str[0] == 'chr1'
multiomics_counts_s = multiomics_counts_s[multiomics_counts_s.obs.sample(2000).index, multiomics_counts_s.var.sample(2000).index]
multiomics_counts_s.write(par['multiomics_counts_test'])
# shorten perturbation
perturbation_counts = ad.read_h5ad(par['perturbation_counts'])
perturbation_counts_s = perturbation_counts[perturbation_counts.obs.donor_id=='Donor 1']
perturbation_counts_s = perturbation_counts_s[perturbation_counts_s.obs.sample(1000).index,perturbation_counts_s.var.sample(1000).index]
perturbation_counts_s.write(par['perturbation_counts_test'])