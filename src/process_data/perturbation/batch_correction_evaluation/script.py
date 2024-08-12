import anndata as ad 
import pandas as pd
import numpy as np
import scanpy as sc
import sys
import warnings
warnings.filterwarnings("ignore")
warnings.filterwarnings('ignore', category=FutureWarning)
warnings.filterwarnings('ignore', category=DeprecationWarning)

## VIASH START
par = {
    'perturbation_data': 'resources/grn-benchmark/perturbation_data.h5ad',
    'output': 'output/batch_correction_metrics.csv'
}
## VIASH END

meta = {
    'resources_dir': './'
}

sys.path.append(meta['resources_dir'])
from helper import run_scib, run_classifier

bulk_adata = ad.read_h5ad(par['perturbation_data'])
print(bulk_adata)

baseline_layer = 'n_counts'
layers = ['n_counts', 'pearson', 'lognorm', 'seurat_lognorm', 'seurat_pearson', 'scgen_lognorm', 'scgen_pearson']
batch_key = 'plate_name'
label_key = 'cell_type' 


def run_metrics(bulk_adata, layer='lognorm', batch_key='plate_name', label_key='cell_type'):
    rr = run_scib(bulk_adata, layer=layer, layer_baseline=baseline_layer, batch_key=batch_key, label_key=label_key)
    print("classifier")
    rr_classifier = run_classifier(bulk_adata, layer, batch_key)
    rr = pd.concat([rr_scib, rr_classifier], axis=1)
    rr.index = [layer]
    return rr


for i, layer in enumerate(layers):
    print('\n', layer)
    rr = run_metrics(bulk_adata, layer=layer, batch_key=batch_key, label_key=label_key)
    if i == 0:
        rr_all = rr 
    else:
        rr_all = pd.concat([rr_all, rr], axis=0)
    print(rr_all)
    rr_all.to_csv(par["output"])