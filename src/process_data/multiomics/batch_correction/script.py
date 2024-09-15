import anndata as ad 
import pandas as pd
import numpy as np
import scanpy as sc
import scgen
import torch 
print(torch.cuda.is_available())

## VIASH START
par = {
    'multiomics_rna': 'resources/grn-benchmark/multiomics_rna_qc.h5ad',
    "multiomics_rna_bc": 'resources/grn-benchmark/multiomics_rna_qc_bc.h5ad'
}
## VIASH END
batch_key = 'donor_id'
label_key = 'cell_type'
bulk_adata = ad.read_h5ad(par['multiomics_rna'])

# normalize 
sc.experimental.pp.normalize_pearson_residuals(bulk_adata)

#  scgen pipeline
train = bulk_adata.copy()
sc.pp.neighbors(train)
sc.tl.umap(train)

scgen.SCGEN.setup_anndata(train, batch_key=batch_key, labels_key=label_key)
model = scgen.SCGEN(train)
model.train(
    max_epochs=100,
    batch_size=64,
    early_stopping=True,
    early_stopping_patience=25
)

corrected_adata = model.batch_removal()

bulk_adata.X = corrected_adata.X
print(f"Saving the file")
bulk_adata.write_h5ad(par['multiomics_rna_bc'])