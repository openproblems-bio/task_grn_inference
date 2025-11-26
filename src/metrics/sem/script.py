import os
import sys
import anndata as ad
import numpy as np
import pandas as pd


## VIASH START
par = {
    # 'prediction': f'resources/results/replogle/replogle.pearson_corr.pearson_corr.prediction.h5ad',
    'prediction': 'resources/results/experiment/global_grns/Gtex:Whole blood.csv.h5ad',
    'evaluation_data': f'resources/grn_benchmark/evaluation_data/replogle_bulk.h5ad',
    'layer': 'lognorm',
    "max_n_links": 50000,
    'num_workers': 20,
    'tf_all': 'resources/grn_benchmark/prior/tf_all.csv',    
    'score': 'output/score.h5ad',
    'genes_n': 5000
}
## VIASH END
run_local = False
try:
    sys.path.append(meta["resources_dir"])
except:
    run_local=True
    meta = {
    "resources_dir":'src/metrics/sem/',
    "util_dir": 'src/utils',
    'helper_dir': 'src/metrics/sem/'
    }
    sys.path.append(meta["resources_dir"])
    sys.path.append(meta["util_dir"])
    sys.path.append(meta["helper_dir"])
from helper import main as main_sem 
from util import format_save_score, parse_args

if run_local:
    par = parse_args(par)


if __name__ == "__main__":
    
    
    output = main_sem(par)
    
    dataset_id = ad.read_h5ad(par['evaluation_data'], backed='r').uns['dataset_id']
    method_id = ad.read_h5ad(par['prediction'], backed='r').uns['method_id']
    format_save_score(output, method_id, dataset_id, par['score'])