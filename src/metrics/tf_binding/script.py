import os
import sys
import anndata as ad
import numpy as np
import pandas as pd


## VIASH START
par = {
    'prediction': f'resources/results/op/op.pearson_corr.pearson_corr.prediction.h5ad',
    'evaluation_data': f'resources/grn_benchmark/evaluation_data/op_bulk.h5ad',
    'ground_truth': f'resources/grn_benchmark/ground_truth/pbmc.csv',
    'layer': 'lognorm',
    "max_n_links": 50000,
    'num_workers': 20,
    'tf_all': 'resources/grn_benchmark/prior/tf_all.csv',    
    'score': 'output/score.h5ad',
    'genes_n': 5000,
    'output_detailed_metrics': True
}
## VIASH END

try:
    sys.path.append(meta["resources_dir"])
except:
    meta = {
    "resources_dir":'src/metrics/experimental/tf_binding/',
    "util_dir": 'src/utils',
    }
    sys.path.append(meta["resources_dir"])
    sys.path.append(meta["util_dir"])
from helper import main as main 
from util import format_save_score, parse_args

args = parse_args(par)


if __name__ == "__main__":
    output = main(par)
    print(output)

    dataset_id = ad.read_h5ad(par['evaluation_data'], backed='r').uns['dataset_id']
    method_id = ad.read_h5ad(par['prediction'], backed='r').uns['method_id']
    
    score_adata = format_save_score(output, dataset_id, method_id, par['score'])
    print('Completed', flush=True)