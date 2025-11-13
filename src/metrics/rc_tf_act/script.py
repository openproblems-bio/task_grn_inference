import os
import sys
import anndata as ad
import numpy as np
import pandas as pd


## VIASH START
par = {
    'prediction': 'resources/results/op/op.grnboost.grnboost.prediction.h5ad',
    'evaluation_data': 'resources/grn_benchmark/evaluation_data/op_bulk.h5ad',
    'layer': 'lognorm',
    'max_n_links': 50000,
    'num_workers': 20,
    'tf_all': 'resources/grn_benchmark/prior/tf_all.csv',    
    'score': 'output/score.h5ad',
    'min_targets': 5
}
## VIASH END

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--prediction', type=str, help='Path to the predicted GRN in h5ad format')
parser.add_argument('--evaluation_data', type=str, help='Path to the evaluation data in h5ad format')
parser.add_argument('--score', type=str, help='Output path for the score h5ad file')
parser.add_argument('--layer', type=str, default='lognorm', help='Layer in the h5ad file to use')
parser.add_argument('--num_workers', type=int, default=20, help='Number of workers to use')
parser.add_argument('--min_targets', type=int, default=5, help='Minimum number of targets per TF')

try:
    sys.path.append(meta["resources_dir"])
except:
    meta = {
        "resources_dir": 'src/metrics/rc_tf_act/',
        "util_dir": 'src/utils',
        'helper_dir': 'src/metrics/rc_tf_act/'
    }
    sys.path.append(meta["resources_dir"])
    sys.path.append(meta["util_dir"])
    sys.path.append(meta["helper_dir"])

from helper import main as main_rc_tf_act
from util import format_save_score


args = parser.parse_args()
for key, value in vars(args).items():
    if value is not None:
        par[key] = value

if __name__ == "__main__":
    output = main_rc_tf_act(par)
    
    dataset_id = ad.read_h5ad(par['evaluation_data'], backed='r').uns['dataset_id']
    method_id = ad.read_h5ad(par['prediction'], backed='r').uns['method_id']
    format_save_score(output, method_id, dataset_id, par['score'])
