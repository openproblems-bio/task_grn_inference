import os, binascii
import scipy
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import scanpy as sc
import pandas as pd
import numpy as np
import random
import anndata as ad
from matplotlib.patches import Patch
import warnings
from scipy import sparse
import lightgbm as lgb
from sklearn.model_selection import cross_validate
from sklearn.linear_model import RidgeClassifier
from sklearn.metrics import r2_score, make_scorer, accuracy_score

batch_key = 'plate_name'
label_key = 'cell_type'

work_dir = '../output'
input_dir = '../input'

bulk_adata = ad.read_h5ad(f"{input_dir}/bulk_adata_integrated.h5ad")

import scib
def run_scib(bulk_adata, layer='lognorm', layer_baseline='n_counts', batch_key='plate_name', label_key='cell_type'):
    bulk_adata.X = bulk_adata.layers[layer_baseline].copy()

    bulk_adata_c = bulk_adata.copy()
    bulk_adata_c.X = bulk_adata_c.layers[layer].copy()

    scib.pp.reduce_data(
        bulk_adata_c, n_top_genes=None, batch_key=batch_key, pca=True, neighbors=True
    )
    rr = scib.metrics.metrics(bulk_adata, bulk_adata_c, batch_key, label_key, organism='human', 
                            # biological conservation (label)
                            nmi_=True, 
                            ari_=False,
                            silhouette_=True,
                            isolated_labels_f1_=False, # there is no isolated cell type
                            isolated_labels_asw_=False, # there is no isolated cell type
                            # biological conservation (label free)
                            cell_cycle_=True,
                            hvg_score_=False,
                            trajectory_=False,
                            # batch correction
                            pcr_=False, 
                            graph_conn_=False,
                            kBET_=True,
                            ilisi_=False,
                            clisi_=False,
                            # Not sure what they are
                            isolated_labels_=False,  # backwards compatibility
                            n_isolated=None,
                            lisi_graph_=False,
                            )
    rr = rr.dropna().T
    return rr 
def run_classifier(adata, layer, batch_key):

    # model = lgb.LGBMClassifier()
    model = RidgeClassifier()
    X = adata.layers[layer].copy()
    y = adata.obs[batch_key]
    scoring = {
        'accuracy_score': make_scorer(accuracy_score)
    }
    score = 1 - cross_validate(model, X, y, cv=5, scoring=scoring, return_train_score=False)['test_accuracy_score'].mean()
    
    return pd.DataFrame({'Batch classifier':[score]})


def run_metrics(bulk_adata, layer='lognorm', batch_key='plate_name', label_key='cell_type'):
    rr_scib = run_scib(bulk_adata, layer=layer, layer_baseline='n_counts', batch_key=batch_key, label_key=label_key)
    rr_classifier = run_classifier(bulk_adata, layer, batch_key)
    rr = pd.concat([rr_scib, rr_classifier], axis=1)
    rr.index = [layer]
    return rr

for i, layer in enumerate(['n_counts','lognorm','pearson','seurat_lognorm', 'seurat_pearson', 'scgen_lognorm', 'scgen_pearson']):
    rr = run_metrics(bulk_adata, layer=layer, batch_key=batch_key, label_key=label_key)
    if i == 0:
        rr_all = rr 
    else:
        rr_all = pd.concat([rr_all, rr], axis=0)
    
rr_all.to_csv(f'{results_dir}/preprocess/batch_correction_metrics.csv')
