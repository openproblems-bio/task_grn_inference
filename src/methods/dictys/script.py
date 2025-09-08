
import pandas as pd
import numpy as np
import os
import argparse
import warnings
import sys
warnings.filterwarnings("ignore")

## VIASH START
par = {
    'rna': 'resources_test/grn_benchmark/inference_data/op_rna.h5ad',
    'atac': 'resources_test/grn_benchmark/inference_data/op_atac.h5ad',
    'temp_dir': 'output/temp/',
    'prediction': 'output/temp/prediction.h5ad',
}
## VIASH END
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--rna', type=str, help='Path to RNA data in h5ad format', required=False)
parser.add_argument('--atac', type=str, help='Path to ATAC data in h5ad format', required=False)
parser.add_argument('--prediction', type=str, help='Path to save the prediction h5ad file', required=False)
args = parser.parse_args()

for k, v in vars(args).items():
    if v is not None:
        par[k] = v


try: 
    sys.path.append(meta["resources_dir"])
except:
    meta = {
        'utils_dir': 'src/utils',
        'helper_dir': 'src/methods/dictys/'
    }
    sys.path.append(meta['utils_dir'])
    sys.path.append(meta['helper_dir'])

from helper import define_vars, format_inputs, export_net 
# from util import process_links

if __name__ == '__main__':
    define_vars(par)
    format_inputs(par)
    # export_net


