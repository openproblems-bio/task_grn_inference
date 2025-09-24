import pandas as pd
import anndata as ad
import sys
import numpy as np
import argparse
import os



## VIASH START
par = {
    'layer': 'lognorm',
    'prediction': f'resources/results/replogle/replogle.negative_control.negative_control.prediction.h5ad',
    # 'evaluation_data_sc': f'resources/grn_benchmark/evaluation_data/replogle_sc.h5ad',
    'ws_consensus': f'resources/grn_benchmark/prior/ws_consensus_replogle.csv',
    'ws_distance_background': f'resources/grn_benchmark/prior/ws_distance_background_replogle.csv',
    'silent_missing_dependencies': False
}
## VIASH END

## LOCAL START
parser = argparse.ArgumentParser()
parser.add_argument('--run_local', action='store_true', help='Run locally')

args = parser.parse_args()
var_local = vars(args)

## LOCAL END

try:
    sys.path.append(meta["resources_dir"])
except:
    meta = {
        "resources_dir":'src/metrics/wasserstein/',
        "util_dir":'src/utils'
    }
    sys.path.append(meta["resources_dir"])
    sys.path.append(meta["util_dir"])

from helper import main 
from util import format_save_score

if __name__ == '__main__':
    # if par.get('ws_consensus') is None:
    #     if par['silent_missing_dependencies']:
    #         dataset_id = 'missing'
    #         metric_ids =['ws']
    #         metric_values = ['']
    #     else:
    #         raise FileNotFoundError(f"Dependencies missing {par['ws_consensus']}. Please check the paths of the dependencies")
    # else:
    _, output = main(par)
    method_id = ad.read_h5ad(par['prediction'], backed='r').uns['method_id']
    dataset_id = ad.read_h5ad(par['evaluation_data'], backed='r').uns['dataset_id']
    format_save_score(output, method_id, dataset_id, par['score'])