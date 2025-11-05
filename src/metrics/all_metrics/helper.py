
import os
import pandas as pd
import numpy as np
import anndata as ad
import scanpy as sc
import sys
import os

try:
    from regression_helper import main as main_reg
except:
    from regression.helper import main as main_reg

try:
    from ws_helper import main as main_ws_distance
except:
    from ws_distance.helper import main as main_ws_distance

try:
    from sem_helper import main as main_sem
except:
    from sem.helper import main as main_sem


try:
    from tf_recovery_helper import main as main_tf_rec
except:
    from tf_recovery.helper import main as main_tf_rec


try:
    from tf_binding_helper import main as main_tf_binding
except:
    from tf_binding.helper import main as main_tf_binding


try:
    from replica_consistency_helper import main as main_replica_consistency
except:
    from replica_consistency.helper import main as main_replica_consistency

from config import DATASETS_METRICS


def sem_metric(par, dataset_id):
    if dataset_id in DATASETS_METRICS:
        if 'sem' in DATASETS_METRICS[dataset_id]:
            output = main_sem(par)
            return output
    return None

def tf_rec_metric(par, dataset_id):
    if dataset_id in DATASETS_METRICS:
        if 'tf_recovery' in DATASETS_METRICS[dataset_id]:
            output = main_tf_rec(par)
            return output
    return None

def tf_binding_metric(par, dataset_id):
    if dataset_id in DATASETS_METRICS:
        if 'tf_binding' in DATASETS_METRICS[dataset_id]:
            output = main_tf_binding(par)
            return output
    return None

def replica_consistency_metric(par, dataset_id):
    if dataset_id in DATASETS_METRICS:
        if 'replica_consistency' in DATASETS_METRICS[dataset_id]:
            try:
                output = main_replica_consistency(par)
            except:
                output = None
            return output
    return None

def reg2_metric(par, dataset_id):
    if dataset_id in DATASETS_METRICS:
        if 'regression' in DATASETS_METRICS[dataset_id]:
            output = main_reg(par)
            return output
    return None

def ws_distance_metric(par, dataset_id):
    if dataset_id in DATASETS_METRICS:
        if 'ws_distance' in DATASETS_METRICS[dataset_id]:
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
