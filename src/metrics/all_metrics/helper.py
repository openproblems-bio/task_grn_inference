
import os
import pandas as pd
import numpy as np
import anndata as ad
import scanpy as sc
import sys
import os


from helper_regression_2 import main as main_reg2
from helper_ws_distance import main as main_ws_distance
from helper_sem import main as main_sem
from helper_tf_recovery import main as main_tf_rec
from helper_tf_binding import main as main_tf_binding
from helper_replica_consistency import main as main_replica_consistency

metrics_dict = {
    'replogle': ['regression_2', 'ws_distance', 'sem', 'tf_recovery', 'tf_binding'],
    'adamon': ['regression_2', 'ws_distance', 'sem', 'tf_recovery', 'tf_binding'],
    'norman': ['regression_2', 'ws_distance', 'sem', 'tf_recovery', 'tf_binding'],
    'nakatake': ['regression_2', 'sem'],
    'op': ['regression_2',  'sem',  'tf_binding', 'replica_consistency'],
    '300BCG': ['regression_2', 'sem',  'tf_binding', 'replica_consistency'],
    'ibd': ['regression_2', 'sem',  'tf_binding', 'replica_consistency'],
    'parsebioscience': ['regression_2','sem', 'tf_binding', 'replica_consistency'],
    'xaira_hek293T': ['regression_2', 'ws_distance', 'sem', 'tf_recovery', 'tf_binding'],
    'xaira_hct116': ['regression_2', 'ws_distance', 'sem', 'tf_recovery', 'tf_binding'],
}

def sem_metric(par, dataset_id):
    if dataset_id in metrics_dict:
        if 'sem' in metrics_dict[dataset_id]:
            output = main_sem(par)
            return output
    return None

def tf_rec_metric(par, dataset_id):
    if dataset_id in metrics_dict:
        if 'tf_recovery' in metrics_dict[dataset_id]:
            output = main_tf_rec(par)
            return output
    return None

def tf_binding_metric(par, dataset_id):
    if dataset_id in metrics_dict:
        if 'tf_binding' in metrics_dict[dataset_id]:
            output = main_tf_binding(par)
            return output
    return None

def replica_consistency_metric(par, dataset_id):
    if dataset_id in metrics_dict:
        if 'replica_consistency' in metrics_dict[dataset_id]:
            output = main_replica_consistency(par)
            return output
    return None

def reg2_metric(par, dataset_id):
    if dataset_id in metrics_dict:
        if 'regression_2' in metrics_dict[dataset_id]:
            output = main_reg2(par)
            return output
    return None

def ws_distance_metric(par, dataset_id):
    if dataset_id in metrics_dict:
        if 'ws_distance' in metrics_dict[dataset_id]:
            _, output = main_ws_distance(par)
            return output
    return None


def main(par):
    dataset_id = ad.read_h5ad(par['evaluation_data'], backed='r').uns['dataset_id']
    rr_store = []
    metrics = [reg2_metric, ws_distance_metric, sem_metric, tf_rec_metric, replica_consistency_metric]
    for metric in metrics:
        rr = metric(par, dataset_id)
        if rr is not None:
            rr_store.append(rr.set_index('key')['value'].rename(metric.__name__.replace('_metric','')))

    rr_all = pd.concat(rr_store, axis=1)
    assert rr_all.shape[1] >0
    return rr_all
