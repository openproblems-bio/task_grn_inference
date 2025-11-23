
import os
import pandas as pd
import numpy as np
import anndata as ad
import scanpy as sc
import sys
import os

from regression.helper import main as regression
from ws_distance.helper import main as ws_distance
from sem.helper import main as sem
from anchor_regression.helper import main as ar
from tf_recovery.helper import main as tf_recovery
from tf_binding.helper import main as tf_binding
from rc_tf_act.helper import main as rc_tf_act
from vc.helper import main as vc
from gs_recovery.helper import main as gs_recovery

from config import DATASETS_METRICS

METRIC_FUNCTIONS = {
    'regression': regression,
    'ws_distance': ws_distance,
    'sem': sem,
    'ar': ar,
    'tf_recovery': tf_recovery,
    'tf_binding': tf_binding,
    'rc_tf_act': rc_tf_act,
    'vc': vc,
    'gs_recovery': gs_recovery,
}

def main(par):
    dataset_id = ad.read_h5ad(par['evaluation_data'], backed='r').uns['dataset_id']
    rr_store = []
    metrics = DATASETS_METRICS[dataset_id]
    # metrics = ['gs_recovery', 'tf_binding']

    for metric_name in metrics:
        print(f"Computing metric: {metric_name}")
        metric_func = METRIC_FUNCTIONS.get(metric_name)
        if metric_func is None:
            print(f"Warning: No function found for metric '{metric_name}'")
            continue
        
        rr = metric_func(par)

        if rr is None:
            raise ValueError(f"Metric function for '{metric_name}' returned None")

        if len(rr)>1:
            rr = rr[1]
        print(rr)
        rr_store.append(rr)

    rr_all = pd.concat(rr_store, axis=1)
    assert rr_all.shape[1] >0
    return rr_all
