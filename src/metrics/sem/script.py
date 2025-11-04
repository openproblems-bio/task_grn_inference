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

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--prediction', type=str, help='Path to the predicted GRN in h5ad format')
parser.add_argument('--evaluation_data', type=str, help='Path to the evaluation data in h5ad format')
parser.add_argument('--score', type=str)
parser.add_argument('--layer', type=str, default='lognorm', help='Layer in the h5ad file to use')
parser.add_argument('--num_workers', type=int, default=20, help='Number of workers to use')

try:
    sys.path.append(meta["resources_dir"])
except:
    meta = {
    "resources_dir":'src/metrics/sem/',
    "util_dir": 'src/utils',
    'helper_dir': 'src/metrics/sem/'
    }
    sys.path.append(meta["resources_dir"])
    sys.path.append(meta["util_dir"])
    sys.path.append(meta["helper_dir"])
from helper import main as main_sem 
from util import format_save_score



args = parser.parse_args()
for key, value in vars(args).items():
    if value is not None:
        par[key] = value

if __name__ == "__main__":
    try:
        output = main_sem(par)
    except Exception as e:
        print(f"Error in SEM evaluation: {e}")
        output = pd.DataFrame({
            'key': ["None"],
            'value': ["None"]
        })

    dataset_id = ad.read_h5ad(par['evaluation_data'], backed='r').uns['dataset_id']
    method_id = ad.read_h5ad(par['prediction'], backed='r').uns['method_id']
    format_save_score(output, method_id, dataset_id, par['score'])