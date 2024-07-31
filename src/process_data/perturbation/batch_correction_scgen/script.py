import anndata as ad 
import pandas as pd
import numpy as np
import scanpy as sc
import scgen

## VIASH START
par = {
    'perturbation_data_n': 'resources_test/grn-benchmark/perturbation_data.h5ad',
    "perturbation_data_bc": 'resources_test/grn-benchmark/perturbation_data.h5ad'
}
## VIASH END

bulk_adata = ad.read_h5ad(par['perturbation_data_n'])
print(bulk_adata)

for norm_name in ['lognorm', 'pearson']:
    train = bulk_adata.copy()
    train.X = train.layers[norm_name]
    sc.pp.neighbors(train)
    sc.tl.umap(train)

    scgen.SCGEN.setup_anndata(train, batch_key=batch_key, labels_key=label_key)
    model = scgen.SCGEN(train)
    model.train(
        max_epochs=100,
        batch_size=64,
        early_stopping=True,
        early_stopping_patience=25,
    )

    corrected_adata = model.batch_removal()

    bulk_adata.layers[f'scgen_{norm_name}'] = corrected_adata.X
bulk_adata.write_h5ad(par['perturbation_data_bc'])