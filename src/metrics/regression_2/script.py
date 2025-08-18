import pandas as pd
import anndata as ad
import sys
import numpy as np
import argparse

## VIASH START
par = {
    'prediction': 'resources/results/replogle/replogle.negative_control.negative_control.prediction.h5ad',
    'evaluation_data': 'resources/grn_benchmark/evaluation_data/replogle_bulk.h5ad',
    'regulators_consensus': 'resources/grn_benchmark/prior/regulators_consensus_replogle.json',
    'layer': 'lognorm',
    "max_n_links": 50000,
    'reg_type': 'ridge',
    'static_only': True,
    'subsample': -1,
    'num_workers': 4,
    'apply_tf': True,
    'clip_scores': True,
    'apply_skeleton': False,
    'skeleton': 'resources/grn_benchmark/prior/skeleton.csv',
    'tf_all': 'resources/grn_benchmark/prior/tf_all.csv',
    'verbose': 2
    
}
## VIASH END

## LOCAL START
parser = argparse.ArgumentParser()
parser.add_argument('--run_local', action='store_true', help='Run locally')

args = parser.parse_args()
var_local = vars(args)

## LOCAL END
    
try:
  sys.path.append(meta["resources_dir"])
except:
    meta = {
      "resources_dir":'src/metrics/regression_2/',
      "util_dir":'src/utils'
    }
    sys.path.append(meta["resources_dir"])
    sys.path.append(meta["util_dir"])
from helper import main

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