# !pip install sctk anndata
# !aws s3 cp s3://openproblems-bio/public/neurips-2023-competition/sc_counts.h5ad  ./resources_raw/ --no-sign-request

import anndata as ad 
import pandas as pd
import numpy as np
import sctk
from scipy import sparse
import scanpy as sc

par = {
    'sc_counts': 'resources/raw-data/sc_counts.h5ad',
    'perturbation_data': 'resources/raw-data/perturbation_data.h5ad',

}
