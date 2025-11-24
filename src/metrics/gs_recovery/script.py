import os
import sys
import anndata as ad
import numpy as np
import pandas as pd
import argparse


## VIASH START
par = {
    'prediction': f'resources/results/op/op.pearson_corr.pearson_corr.prediction.h5ad',
    'evaluation_data': f'resources/grn_benchmark/evaluation_data/op_bulk.h5ad',
    'layer': 'lognorm',
    "max_n_links": 50000,
    'num_workers': 20,
    'tf_all': 'resources/grn_benchmark/prior/tf_all.csv',    
    'score': 'output/score.h5ad',
    'fdr_threshold': 0.05,
    'min_pathway_size': 5,
    'max_pathway_size': 500,
    'min_targets': 10,
    'max_targets': 100,  # Top K edges by absolute weight
    'ulm_baseline_method': 'zero_centered'
}
## VIASH END

run_local=False
try:
    sys.path.append(meta["resources_dir"])
except:
    run_local=True
    meta = {
    "resources_dir":'src/metrics/experimental/annotation/',
    "util_dir": 'src/utils',
    }
    sys.path.append(meta["resources_dir"])
    sys.path.append(meta["util_dir"])
from helper import main as main_helper
from util import format_save_score, parse_args

if run_local:
    par = parse_args(par)

if __name__ == "__main__":
    print(par)
    output = main_helper(par)
    print(output)

    dataset_id = ad.read_h5ad(par['evaluation_data'], backed='r').uns['dataset_id']
    method_id = ad.read_h5ad(par['prediction'], backed='r').uns['method_id']
    format_save_score(output, method_id, dataset_id, par['score'])
    print(par['score'])
