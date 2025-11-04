
import os
import pandas as pd
import numpy as np
import anndata as ad
import scanpy as sc
import sys
import os


from regression.helper import main as main_reg2
from ws_distance.helper import main as main_ws_distance
from sem.helper import main as main_sem
from tf_recovery.helper import main as main_tf_rec
from tf_binding.helper import main as main_tf_binding
from replica_consistency.helper import main as main_replica_consistency
from metrics_config import datasets_metrics


def sem_metric(par, dataset_id):
    if dataset_id in datasets_metrics:
        if 'sem' in datasets_metrics[dataset_id]:
            output = main_sem(par)
            return output
    return None

def tf_rec_metric(par, dataset_id):
    if dataset_id in datasets_metrics:
        if 'tf_recovery' in datasets_metrics[dataset_id]:
            output = main_tf_rec(par)
            return output
    return None

def tf_binding_metric(par, dataset_id):
    if dataset_id in datasets_metrics:
        if 'tf_binding' in datasets_metrics[dataset_id]:
            output = main_tf_binding(par)
            return output
    return None

def replica_consistency_metric(par, dataset_id):
    if dataset_id in datasets_metrics:
        if 'replica_consistency' in datasets_metrics[dataset_id]:
            try:
                output = main_replica_consistency(par)
            except:
                output = None
            return output
    return None

def reg2_metric(par, dataset_id):
    if dataset_id in datasets_metrics:
        if 'regression' in datasets_metrics[dataset_id]:
            output = main_reg2(par)
            return output
    return None

def ws_distance_metric(par, dataset_id):
    if dataset_id in datasets_metrics:
        if 'ws_distance' in datasets_metrics[dataset_id]:
            _, output = main_ws_distance(par)
            return output
    return None


def main(par):
    dataset_id = ad.read_h5ad(par['evaluation_data'], backed='r').uns['dataset_id']
    rr_store = []
    metrics = [reg2_metric, ws_distance_metric, sem_metric, tf_rec_metric, replica_consistency_metric]
    for metric in metrics:
        print(f"Computing metric: {metric.__name__}")
        rr = metric(par, dataset_id)
        if rr is not None:
            if 'key' in rr.columns:
                if rr['key']=="None":
                    print(f"Skipping metric {metric.__name__} due to None output")
                    continue
            rr_store.append(rr)

    rr_all = pd.concat(rr_store, axis=1)
    assert rr_all.shape[1] >0
    return rr_all
