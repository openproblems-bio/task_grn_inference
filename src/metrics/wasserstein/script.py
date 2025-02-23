import pandas as pd
import anndata as ad
import sys
import numpy as np
import argparse
import os


## VIASH START
par = {
    'layer': 'X_norm'
}
## VIASH END

## LOCAL START
parser = argparse.ArgumentParser()
parser.add_argument('--run_local', action='store_true', help='Run locally')
parser.add_argument('--evaluation_data_sc', type=str)
parser.add_argument('--ws_consensus', type=str)
parser.add_argument('--ws_distance_background', type=str)
parser.add_argument('--prediction', type=str, help='Path to the prediction file')
parser.add_argument('--method_id', type=str, help='Method id')
parser.add_argument('--dataset_id', type=str, help='Dataset id')
parser.add_argument('--score', type=str, help='score file')

args = parser.parse_args()
var_local = vars(args)

## LOCAL END

if args.run_local:
    for key in var_local:
        if var_local[key] is not None:
            par[key] = var_local[key]
    meta = {
      "resources_dir":'src/metrics/wasserstein/',
      "util_dir":'src/utils'
    }
    sys.path.append(meta["resources_dir"])
    sys.path.append(meta["util_dir"])

else:
  sys.path.append(meta["resources_dir"])


from main import main 

if __name__ == '__main__':
    method_id = ad.read_h5ad(par['prediction'], backed='r').uns['method_id']

    # - check dependencies
    if par.get('ws_consensus') is None:
        if par['silent_missing_dependencies']:
            dataset_id = 'missing'
            metric_ids =['ws']
            metric_values = ['']
        else:
            raise FileNotFoundError(f"Dependencies missing {par['ws_consensus']}. Please check the paths of the dependencies")
    else:
        method_id = ad.read_h5ad(par['prediction'], backed='r').uns['method_id']
        dataset_id = ad.read_h5ad(par['evaluation_data_sc'], backed='r').uns['dataset_id']
        print(f"Method id: {method_id}, Dataset id: {dataset_id}")
        # - main function
        _, mean_scores = main(par)
        print(mean_scores)
        metric_ids = mean_scores.columns.values
        metric_values = mean_scores.values[0]

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