import pandas as pd
import anndata as ad
import sys
import numpy as np
import argparse

## VIASH START
par = {
    'layer': 'X_norm',
    "max_n_links": 50000,
    'reg_type': 'ridge',
    'static_only': True,
    'subsample': -1,
    'num_workers': 4,
    'apply_tf': True,
    'clip_scores': True,
    'method_id': 'grnboost',
    'apply_skeleton': False,
    'skeleton': 'resources/grn_benchmark/prior/skeleton.csv',
    'tf_all': 'resources/grn_benchmark/prior/tf_all.csv',
    'verbose': 2
    
}
## VIASH END

## LOCAL START
parser = argparse.ArgumentParser()
parser.add_argument('--run_local', action='store_true', help='Run locally')
parser.add_argument('--evaluation_data', type=str)
parser.add_argument('--regulators_consensus', type=str)
parser.add_argument('--prediction', type=str, help='Path to the prediction file')
parser.add_argument('--method_id', type=str, help='Method id')
parser.add_argument('--dataset_id', type=str, help='Dataset id')
parser.add_argument('--score', type=str, help='score file')
parser.add_argument('--reg_type', type=str)
parser.add_argument('--apply_skeleton', action='store_true')

args = parser.parse_args()
var_local = vars(args)

## LOCAL END

if args.run_local:
    for key in var_local:
        if var_local[key] is not None:
            par[key] = var_local[key]
    meta = {
      "resources_dir":'src/metrics/regression_2/',
      "util_dir":'src/utils'
    }
    sys.path.append(meta["resources_dir"])
    sys.path.append(meta["util_dir"])

else:
  sys.path.append(meta["resources_dir"])

from main import main


if __name__ == '__main__':
    method_id = ad.read_h5ad(par['prediction'], backed='r').uns['method_id']
    dataset_id = ad.read_h5ad(par['evaluation_data'], backed='r').uns['dataset_id']
    print(f"Method id: {method_id}, Dataset id: {dataset_id}")

    # - Main function
    output = main(par)

    print('Write output to file', flush=True)
    print(output)
    metric_ids = output.columns.to_numpy()
    metric_values = output.values[0]

    output = ad.AnnData(
        X=np.empty((0, 0)),
        uns={
            "dataset_id": dataset_id,
            "method_id": method_id,
            "metric_ids": metric_ids,
            "metric_values": metric_values
        }
    )
    output.write_h5ad(par['score'], compression='gzip')
    print('Completed', flush=True)
