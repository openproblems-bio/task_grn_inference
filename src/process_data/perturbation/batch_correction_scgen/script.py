import anndata as ad 
import pandas as pd
import numpy as np
import scanpy as sc
import scgen

## VIASH START
par = {
    'perturbation_data_n': 'resources/grn-benchmark/perturbation_data.h5ad',
    "perturbation_data_bc": 'resources/grn-benchmark/perturbation_data.h5ad',
    'batch_key': 'plate_name',
    'label_key': 'cell_type'
}
## VIASH END
import argparse
parser = argparse.ArgumentParser(description="Batch correction")
parser.add_argument('--perturbation_data_n', type=str, help='Path to the anndata file')
parser.add_argument('--perturbation_data_bc', type=str, help='Path to the anndata file')
parser.add_argument('--batch_key', type=str, help='Batch name')
parser.add_argument('--label_key', type=str, help='label name')

args = parser.parse_args()

if args.perturbation_data_n:
    par['perturbation_data_n'] = args.perturbation_data_n
if args.perturbation_data_bc:
    par['perturbation_data_bc'] = args.perturbation_data_bc
if args.label_key:
    par['label_key'] = args.label_key
if args.batch_key:
    par['batch_key'] = args.batch_key

print(par)

bulk_adata = ad.read_h5ad(par['perturbation_data_n'])
print(bulk_adata)

for norm_name in ['lognorm', 'pearson']:
    train = bulk_adata.copy()
    train.X = train.layers[norm_name]
    sc.pp.neighbors(train)
    sc.tl.umap(train)

    scgen.SCGEN.setup_anndata(train, batch_key=par['batch_key'], labels_key=par['label_key'])
    model = scgen.SCGEN(train)
    model.train(
        max_epochs=100,
        batch_size=64,
        early_stopping=True,
        early_stopping_patience=25,
    )

    corrected_adata = model.batch_removal()

    bulk_adata.layers[f'scgen_{norm_name}'] = corrected_adata.X
    print(f"{norm_name} finished batch correction")
bulk_adata.write_h5ad(par['perturbation_data_bc'])