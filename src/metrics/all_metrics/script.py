import pandas as pd
import anndata as ad
import sys
import numpy as np
import os 


## VIASH START
par = {
    'prediction': f'resources/grn_models/norman/grnboost2.csv',
    'method_id': 'grnboost2',

    "tf_all": f"resources/prior/tf_all.csv",
    'skeleton': f'resources/prior/skeleton.csv', 
    'dataset_id': 'norman',
    'layer': 'X_norm',
    "apply_tf": True,
    'subsample': -1,
    'verbose': 4,
    'num_workers': 20,
    'binarize': False, 
    'max_n_links': 50000, 
    'apply_skeleton': False, 
    'reg_type':'ridge', 
    'score': 'output/score.h5ad'
}
## VIASH END

meta = {
  "resources_dir": 'src/metrics/',
  "util": 'src/utils'
}
sys.path.append(meta["resources_dir"])
sys.path.append(meta["util"])
from regression_1.main import main as main_reg1
from regression_2.main import main as main_reg2
from wasserstein.script import main as main_ws



def main(par):
    """
        Calculate all scores for a given model and daatset.
    """
    assert par['dataset_id']
    dataset = par['dataset_id']

    par['evaluation_data'] = f'resources/evaluation_datasets/{dataset}_perturbation.h5ad'
    par['evaluation_data_sc'] = f'resources/datasets_raw/{dataset}_sc_counts.h5ad'
    par['regulators_consensus'] = f'resources/prior/regulators_consensus_{dataset}.json'
    par['ws_consensus'] = f'resources/prior/ws_consensus_{dataset}.csv'
    par['ws_distance_background'] = f'resources/prior/ws_distance_background_{dataset}.csv'
    
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

    output = ad.AnnData(
        X=np.empty((0, 0)),
        uns={
            "dataset_id": par["dataset_id"],
            "method_id": par['method_id'],
            "metric_ids": scores_all.columns.values,
            "metric_values": scores_all.values[0]
        }
    )
    print(output)
    output.write_h5ad(par['score'], compression='gzip')
    print('Completed', flush=True)