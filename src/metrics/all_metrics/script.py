import pandas as pd
import anndata as ad
import sys
import numpy as np
import os 
import argparse
import scanpy as sc

argparser = argparse.ArgumentParser()
argparser.add_argument('--run_local', action='store_true', help='Run locally')
argparser.add_argument('--prediction', help='Path to the GRN prediction file')
argparser.add_argument('--evaluation_dataset')
argparser.add_argument('--evaluation_dataset_sc')
argparser.add_argument('--regulators_consensus')
argparser.add_argument('--ws_consensus')
argparser.add_argument('--ws_distance_background')
argparser.add_argument('--method_id', help='Method ID')
argparser.add_argument('--dataset_id', help='Dataset ID')
argparser.add_argument('--score', help='Where to store scores')

par_local = vars(argparser.parse_args())

## VIASH START
par = {
    "tf_all": f"resources/grn_benchmark/prior/tf_all.csv",
    'skeleton': f'resources/grn_benchmark/prior/skeleton.csv', 
    'layer': 'X_norm',
    "apply_tf": True,
    'subsample': -1,
    'verbose': 4,
    'num_workers': 20,
    'binarize': False, 
    'max_n_links': 50000, 
    'apply_skeleton': False, 
    'reg_type':'ridge', 
}
## VIASH END
if par_local['run_local']:
    for key in par_local:
        par[key] = par_local[key]
    
    meta = {
        "resources_dir": 'src/metrics/',
    }
    sys.path.append(meta["resources_dir"])
    from regression_1.script import main as main_reg1
    from regression_2.main import main as main_reg2
    from wasserstein.script import main as main_ws
else:
    sys.path.append(meta["resources_dir"])
    from reg1_main import main as main_reg1
    from reg2_main import main as main_reg2
    from ws_main import main as main_ws

# try:
#     sys.path.append(meta["resources_dir"])
#     from reg1_main.main import main as main_reg1
#     from reg2_main.main import main as main_reg2
#     from ws_main.main import main as main_ws
    
# except:
#     meta = {
#     "resources_dir": 'src/metrics/',
#     "util": 'src/utils'
#     }
#     sys.path.append(meta["resources_dir"])
#     sys.path.append(meta["util"])
#     from regression_1.main import main as main_reg1
#     from regression_2.main import main as main_reg2
#     from wasserstein.script import main as main_ws



def main(par):
    """
        Calculate all scores for a given model and dataset.
    """
    assert par['dataset_id']
    dataset = par['dataset_id']

    # par['evaluation_data'] = f'resources/grn_benchmark/evaluation_datasets/{dataset}_perturbation.h5ad'
    # par['evaluation_data_sc'] = f'resources/datasets_raw/{dataset}_sc_counts.h5ad'
    # par['regulators_consensus'] = f'resources/grn_benchmark/prior/regulators_consensus_{dataset}.json'
    # par['ws_consensus'] = f'resources/grn_benchmark/prior/ws_consensus_{dataset}.csv'
    # par['ws_distance_background'] = f'resources/grn_benchmark/prior/ws_distance_background_{dataset}.csv'
    
    scores_all = []

    scores_reg1 = main_reg1(par)
    scores_all.append(scores_reg1)
    scores_reg2 = main_reg2(par)
    scores_all.append(scores_reg2)
    if dataset in ['norman', 'adamson']:
        print(par)
        _, scores_ws = main_ws(par)
        scores_all.append(scores_ws)

    scores_all = pd.concat(scores_all, axis=1)

    return scores_all
if __name__ == '__main__':
    scores_all = main(par)

    print(scores_all)

    output = ad.AnnData(
        X=np.empty((0, 0)),
        uns={
            "dataset_id": par["dataset_id"],
            "method_id": par['method_id'],
            "metric_ids": scores_all.columns.values,
            "metric_values": scores_all.values[0]
        }
    )
    output.write_h5ad(par['score'], compression='gzip')
    print('Completed', flush=True)
