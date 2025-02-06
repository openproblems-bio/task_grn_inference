import pandas as pd
import anndata as ad
import sys
import numpy as np
import argparse


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
    _, mean_scores = main(par)
    print(mean_scores)
    output = ad.AnnData(
        X=np.empty((0, 0)),
        uns={
            "dataset_id": str(par["dataset_id"]),
            "method_id": f"reg2-{par['method_id']}",
            "metric_ids": mean_scores.columns.values,
            "metric_values": mean_scores.values[0]
        }
    )
    output.write_h5ad(par['score'], compression='gzip')
    print('Completed', flush=True)