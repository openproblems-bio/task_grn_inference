import pandas as pd
import anndata as ad
import sys
import numpy as np
import sys 

## VIASH START
par = {
  "rna": "resources/grn-benchmark/rna.h5ad",
  "prediction": "resources/grn_models/default/negative_control.csv",
  "tf_all": "resources/grn_benchmark/prior/tf_all.csv",
  "max_n_links": 50000
}
## VIASH END

import argparse
parser = argparse.ArgumentParser(description="Process multiomics RNA data.")
parser.add_argument('--rna', type=str, help='Path to the multiomics RNA file')
parser.add_argument('--prediction', type=str, help='Path to the prediction file')
parser.add_argument('--resources_dir', type=str, help='Path to the prediction file')
parser.add_argument('--tf_all', type=str, help='Path to the tf_all')
parser.add_argument('--num_workers', type=str, help='Number of cores')
parser.add_argument('--max_n_links', type=str, help='Number of top links to retain')
parser.add_argument('--normalize', action='store_true')
args = parser.parse_args()
if args.rna:
    par['rna'] = args.rna
if args.prediction:
    par['prediction'] = args.prediction
if args.tf_all:
    par['tf_all'] = args.tf_all
if args.num_workers:
    par['num_workers'] = args.num_workers
if args.max_n_links:
    par['max_n_links'] = int(args.max_n_links)
if args.resources_dir:
    meta['resources_dir'] = args.resources_dir  
print(par)

try:
    sys.path.append(meta['resources_dir'])
except:
    meta = {
    'resources_dir':'src/control_methods/negative_control'
    }
    sys.path.append(meta['resources_dir'])
from main import main

if __name__ == '__main__':
    dataset_id = ad.read_h5ad(par['inference_data'], backed='r').uns['dataset_id']

    net = main(par)

    print('Output GRN')
    net['weight'] = net['weight'].astype(str)
    output = ad.AnnData(X=None, uns={"method_id": 'negative_control', "dataset_id": dataset_id, "prediction": net[["source", "target", "weight"]]})
    output.write(par['prediction'])

