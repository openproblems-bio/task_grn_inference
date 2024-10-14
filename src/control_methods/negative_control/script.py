import pandas as pd
import anndata as ad
import sys
import numpy as np
import sys 

## VIASH START
par = {
  "perturbation_data": "resources/grn-benchmark/perturbation_data.h5ad",
  "prediction": "resources/grn_models/default/negative_control.csv",
  "tf_all": "resources/prior/tf_all.csv",
  "max_n_links": 50000
}
## VIASH END
meta = {
    'resources_dir': 'src/utils'
    }
import argparse
parser = argparse.ArgumentParser(description="Process multiomics RNA data.")
parser.add_argument('--multiomics_rna', type=str, help='Path to the multiomics RNA file')
parser.add_argument('--multiomics_atac', type=str, help='Path to the multiomics RNA file')
parser.add_argument('--prediction', type=str, help='Path to the prediction file')
parser.add_argument('--resources_dir', type=str, help='Path to the prediction file')
parser.add_argument('--tf_all', type=str, help='Path to the tf_all')
parser.add_argument('--num_workers', type=str, help='Number of cores')
parser.add_argument('--max_n_links', type=str, help='Number of top links to retain')
args = parser.parse_args()
if args.multiomics_rna:
    par['multiomics_rna'] = args.multiomics_rna
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
meta = {
  'resources_dir':'src/control_methods/negative_control'
}
sys.path.append(meta['resources_dir'])
from main import main
prediction = main(par)
prediction.to_csv(par["prediction"])

