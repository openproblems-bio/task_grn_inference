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

DATASET_GROUPS = {
    "op": {
        "match": ["plate_name", "donor_id", "cell_type", 'well'],
        "loose_match": ["donor_id", "cell_type", "plate_name"],
        "cv": ["perturbation", "cell_type"],
    },
    "parsebioscience": {
        "match": ["donor_id", "cell_type", "well"],
        "loose_match": ["donor_id", "cell_type"],
        "cv": ["perturbation", "cell_type"],
    },
    "300BCG": {
        "match": ["donor_id",  "cell_type"],
        "loose_match": ["cell_type"],
        "cv": ["perturbation", "cell_type"],
    },
}
if __name__ == "__main__":
    dataset_id = ad.read_h5ad(par['evaluation_data'], backed='r').uns['dataset_id']
    method_id = ad.read_h5ad(par['prediction'], backed='r').uns['method_id']

    par['cv_groups'] = DATASET_GROUPS[dataset_id]['cv']
    par['match'] = DATASET_GROUPS[dataset_id]['match']
    par['loose_match'] = DATASET_GROUPS[dataset_id]['loose_match']
    par['method_id'] = method_id
    output = main_sem(par)

    method_id = ad.read_h5ad(par['prediction'], backed='r').uns['method_id']
    dataset_id = ad.read_h5ad(par['evaluation_data'], backed='r').uns['dataset_id']
    format_save_score(output, method_id, dataset_id, par['score'])