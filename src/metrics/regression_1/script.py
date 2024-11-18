import pandas as pd
import anndata as ad
import sys
import numpy as np


file_name = 'op' #nakatake

## VIASH START
par = {
  "evaluation_data": f"resources/evaluation_datasets/{file_name}_perturbation.h5ad",
  "tf_all": "resources/prior/tf_all.csv",
  # "prediction": "output/pearson_net.csv",
  "prediction": f"resources/grn_models/{file_name}/grnboost2.csv",
  "method_id": "scenic",
  "max_n_links": 50000,
  "apply_tf": True,
  'score': 'output/score.h5ad',
  'reg_type': 'ridge',
  'layer': 'pearson',
  'subsample': -1,
  'num_workers': 4,
  'skeleton': 'resources/prior/skeleton.csv',
  'apply_skeleton': False,
  'verbose': 4,
  'binarize': True
}
## VIASH END

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--evaluation_data', type=str, help='Path to the evaluation_data file')
parser.add_argument('--prediction', type=str, help='Path to the prediction file')
parser.add_argument('--tf_all', type=str, help='Path to the tf_all')
parser.add_argument('--num_workers', type=str, help='Number of cores')
parser.add_argument('--layer', type=str, help='Layer to use')
parser.add_argument('--causal', action='store_true', help='Enable causal mode')
parser.add_argument('--normalize', action='store_true')

args = parser.parse_args()

if args.evaluation_data:
    par['evaluation_data'] = args.evaluation_data
if args.layer:
  par['layer'] = args.layer
if args.causal:
    par['causal'] = True
else:
    par['causal'] = False

if args.prediction:
    par['prediction'] = args.prediction
if args.tf_all:
    par['tf_all'] = args.tf_all
if args.num_workers:
    par['num_workers'] = args.num_workers

try:
  sys.path.append(meta["resources_dir"])
except:
  meta = {
    "resources_dir":'src/metrics/regression_1/',
    "util_dir":'src/utils'
  }
  sys.path.append(meta["resources_dir"])
  sys.path.append(meta["util_dir"])

from main import main 
print(par)
output = main(par) 
print(output)

metric_ids = output.columns.to_numpy()
metric_values = output.values[0]

print(metric_ids.shape, metric_values.shape)
output = ad.AnnData(
    X=np.empty((0, 0)),
    uns={
        "dataset_id": par["layer"],
        "method_id": f"reg1-{par['method_id']}",
        "metric_ids": metric_ids,
        "metric_values": metric_values
    }
)

output.write_h5ad(par["score"], compression="gzip")