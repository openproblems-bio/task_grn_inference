import os
import sys
import anndata as ad
import numpy as np
import pandas as pd


## VIASH START
par = {
    'prediction': f'resources/results/op/op.pearson_corr.pearson_corr.prediction.h5ad',
    'evaluation_data_de': f'resources/grn_benchmark/evaluation_data/op_de.csv',
    'layer': 'lognorm',
    "max_n_links": 50000,
    'num_workers': 20,
    'tf_all': 'resources/grn_benchmark/prior/tf_all.csv',    
    'score': 'output/score.h5ad',
}
## VIASH END

try:
    sys.path.append(meta["resources_dir"])
except:
    meta = {
    "util_dir": 'src/utils',
    'helper_dir': 'src/metrics/recovery/'
    }
    sys.path.append(meta["util_dir"])
    sys.path.append(meta["helper_dir"])
from helper import main as main 
from util import format_save_score, parse_args

args = parse_args(par)


if __name__ == "__main__":
    print(par)
    try:
        output = main(par)
        method_id = ad.read_h5ad(par['prediction'], backed='r').uns['method_id']
        dataset_id = ad.read_h5ad(par['evaluation_data'], backed='r').uns['dataset_id']
    except:
        print(f"Error in TF recovery evaluation")
        output = pd.DataFrame({
            'key': ["None"],
            'value': ["None"]
        })
        method_id = ad.read_h5ad(par['prediction'], backed='r').uns['method_id']
        dataset_id = "None"
        
    format_save_score(output, method_id, dataset_id, par['score'])