import os
import sys
import anndata as ad
import numpy as np
import pandas as pd


## VIASH START
par = {
    'prediction': f'resources/results/op/op.pearson_corr.pearson_corr.prediction.h5ad',
    'evaluation_data': f'resources/grn_benchmark/evaluation_data/op_bulk.h5ad',
    'layer': 'lognorm',
    "max_n_links": 50000,
    'num_workers': 20,
    'tf_all': 'resources/grn_benchmark/prior/tf_all.csv',    
    'score': 'output/score.h5ad',
    'genes_n': 5000,
    'use_improved': True  # Flag to use improved implementation
}
## VIASH END

try:
    sys.path.append(meta["resources_dir"])
except:
    meta = {
    "resources_dir":'src/metrics/vc/',
    "util_dir": 'src/utils',
    'helper_dir': 'src/metrics/vc/'
    }
    sys.path.append(meta["resources_dir"])
    sys.path.append(meta["util_dir"])
    sys.path.append(meta["helper_dir"])

# Import appropriate helper based on flag
if par.get('use_improved', False):
    from improved_helper import main as main_vc
    print("Using improved VC implementation with numerical stability fixes")
else:
    from helper import main as main_vc
    print("Using original VC implementation")

from util import format_save_score, parse_args

args = parse_args(par)


if __name__ == "__main__":
    output = main_vc(par)

    dataset_id = ad.read_h5ad(par['evaluation_data'], backed='r').uns['dataset_id']
    method_id = ad.read_h5ad(par['prediction'], backed='r').uns['method_id']
    format_save_score(output, method_id, dataset_id, par['score'])