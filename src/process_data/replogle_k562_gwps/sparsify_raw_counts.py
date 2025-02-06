

import anndata as ad
from scipy import sparse
adata = ad.read_h5ad('resources/datasets_raw/replogle_K562_gwps_raw_singlecell.h5ad')
adata.X = sparse.csr_matrix(adata.X)
adata.write_h5ad('resources/datasets_raw/replogle_K562_gwps_raw_singlecell.h5ad')

