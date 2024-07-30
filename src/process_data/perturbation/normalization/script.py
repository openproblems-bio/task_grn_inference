# !pip install sctk anndata
# !aws s3 cp s3://openproblems-bio/public/neurips-2023-competition/sc_counts.h5ad  ./resources_raw/ --no-sign-request

import anndata as ad 
import pandas as pd
import numpy as np
import sctk
from scipy import sparse
import scanpy as sc

par = {
    'perturbation_data': 'resources/grn-benchmark/perturbation_data.h5ad',
}

def normalize_func(bulk_adata):
    # normalize data based on mean pseudobulk
    bulk_adata.X = bulk_adata.layers['counts'].copy()

    bulk_adata_c = bulk_adata.copy()
    sc.experimental.pp.normalize_pearson_residuals(bulk_adata_c)
    bulk_adata.layers['pearson'] = bulk_adata_c.X

    bulk_adata_c = bulk_adata.copy()
    sc.pp.normalize_total(bulk_adata_c)
    sc.pp.log1p(bulk_adata_c)
    sc.pp.scale(bulk_adata_c)
    bulk_adata.layers['lognorm'] = bulk_adata_c.X
    
    return bulk_adata
    

bulk_adata_n = normalize_func(bulk_adata_filtered)
bulk_adata_n.write(par['perturbation_data'])
