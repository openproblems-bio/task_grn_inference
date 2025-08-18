import os
import pandas as pd
import numpy as np
import sys
import anndata as ad
import scanpy as sc


## VIASH START
par = {
    "rna": "resources/grn_benchmark/inference_data/op_rna.h5ad",
    'layer': 'X_norm',
    'granularity': 1,
    "rna_pseudobulked": "output/imputed.h5ad"
}
## VIASH END

sys.path.append(meta["resources_dir"])
from util import sum_by

def main(par):

    # - read inputs and cluster with differen resolutions
    rna = ad.read_h5ad(par['rna'])
    granularity = par['granularity']

    if granularity == -1:
        pass
    else:
        sc.pp.pca(rna, layer=par['layer'])
        sc.pp.neighbors(rna)
        sc.tl.leiden(rna, resolution=granularity, key_added=f'leiden_{granularity}')
        rna_bulk = sum_by(rna, f'leiden_{granularity}', unique_mapping=True)
        rna_bulk.layers[par['layer']] = rna_bulk.X
        for key in rna.uns.keys():
            rna_bulk.uns[key] = rna.uns[key]
        rna = rna_bulk

    rna.write(par['rna_pseudobulked'])
if __name__ == '__main__':
    main(par)