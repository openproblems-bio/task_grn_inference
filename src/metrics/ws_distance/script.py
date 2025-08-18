import pandas as pd
import anndata as ad
import sys
import numpy as np
import argparse
import os



## VIASH START
par = {
    'layer': 'lognorm',
    'prediction': f'resources/results/replogle/replogle.negative_control.negative_control.prediction.h5ad',
    # 'evaluation_data_sc': f'resources/grn_benchmark/evaluation_data/replogle_sc.h5ad',
    'ws_consensus': f'resources/grn_benchmark/prior/ws_consensus_replogle.csv',
    'ws_distance_background': f'resources/grn_benchmark/prior/ws_distance_background_replogle.csv',
    'silent_missing_dependencies': False
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
        "resources_dir":'src/metrics/wasserstein/',
        "util_dir":'src/utils'
    }
    sys.path.append(meta["resources_dir"])
    sys.path.append(meta["util_dir"])

from helper import main 

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
        # dataset_id = ad.read_h5ad(par['evaluation_data_sc'], backed='r').uns['dataset_id']
        if 'replogle' in par['ws_consensus']:
            dataset_id = 'replogle'
        elif 'adamson' in par['ws_consensus']:
            dataset_id = 'adamson'
        elif 'norman' in par['ws_consensus']:
            dataset_id = 'norman'
        elif 'xaira_HCT116' in par['ws_consensus']:
            dataset_id = 'xaira_HCT116'
        elif 'xaira_HEK293T' in par['ws_consensus']:
            dataset_id = 'xaira_HEK293T'
        else:
            raise ValueError(f"Dataset name did not match expected datasets")
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