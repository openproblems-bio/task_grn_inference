import pandas as pd
import anndata as ad
import sys
import numpy as np
from sklearn.metrics import r2_score
from sklearn.pipeline import Pipeline
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler
import lightgbm
import random 
from concurrent.futures import ThreadPoolExecutor
from tqdm import tqdm
from sklearn.multioutput import MultiOutputRegressor
import os
import warnings

warnings.filterwarnings("ignore", category=FutureWarning)


## VIASH START
par = {
  "evaluation_data": f"resources/grn_benchmark/evaluation_datasets//op.h5ad",
  "tf_all": "resources/grn_benchmark/prior/tf_all.csv",
  "prediction": f"resources/grn_models/op/grnboost2.csv",
  "method_id": "scenic",
  "max_n_links": 50000,
  "apply_tf": True,
  'score': 'output/score.h5ad',
  'reg_type': 'ridge',
  'layer': 'pearson',
  'subsample': -1,
  'num_workers': 4,
  'skeleton': 'resources/grn_benchmark/prior/skeleton.csv',
  'apply_skeleton': False,
  'verbose': 4,
  'binarize': True
}
## VIASH END

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--run_local', action='store_true', help='Run locally')
parser.add_argument('--evaluation_data', type=str, help='Path to the evaluation_data file')
parser.add_argument('--prediction', type=str, help='Path to the prediction file')
parser.add_argument('--tf_all', type=str, help='Path to the tf_all')
parser.add_argument('--num_workers', type=str, help='Number of cores')
parser.add_argument('--layer', type=str, help='Layer to use')
parser.add_argument('--causal', action='store_true', help='Enable causal mode')
parser.add_argument('--normalize', action='store_true')

args = parser.parse_args()

var_local = vars(args)

if args.run_local:
    for key in var_local:
        par[key] = var_local[key]
    meta = {
      "resources_dir":'src/metrics/regression_1/',
      "util_dir":'src/utils'
    }
    sys.path.append(meta["resources_dir"])
    sys.path.append(meta["util_dir"])

else:
  sys.path.append(meta["resources_dir"])

from util import verbose_print, process_links, verbose_tqdm, process_net, binarize_weight
from helper import set_global_seed

def read_net(par):
    net = ad.read_h5ad(par['prediction'])
    net = pd.DataFrame(net.uns['prediction'])
    if par['binarize']:
      net['weight'] = net['weight'].apply(binarize_weight)  
    if par['apply_skeleton']: #apply skeleton
        print('Before filtering with skeleton:', net.shape)
        skeleton = pd.read_csv(par['skeleton'])
        skeleton['link'] = skeleton['source'].astype(str) + '_' + skeleton['target'].astype(str)
        skeleton = skeleton['link'].values.flatten()
        
        net['link'] = net['source'].astype(str) + '_' + net['target'].astype(str)
        net = net[net['link'].isin(skeleton)]
        print('After filtering with skeleton:', net.shape)

    if par['apply_tf']:
        net = net[net.source.isin(tf_all)]
    if 'cell_type' in net.columns:
        print('Taking mean of cell type specific grns')
        net.drop(columns=['cell_type'], inplace=True)
        net = net.groupby(['source', 'target']).mean().reset_index()
      
    if (len(net)>par['max_n_links']) and (par['max_n_links']!=-1):
        net = process_links(net, par) #only top links are considered
        verbose_print(par['verbose'], f"Number of links reduced to {par['max_n_links']}", 2)
  

    return net
    
def main(par):
    random_state = 42
    set_global_seed(random_state)

    # -- read input data
    net = read_net(par)
    evaluation_data = ad.read_h5ad(par['evaluation_data'])
    evaluation_data.X = evaluation_data.layers[par["layer"]]
    gene_names = evaluation_data.var.index.to_numpy()
    tf_all = np.loadtxt(par['tf_all'], dtype=str)
    
    # -- calclate scores 
    results = cross_validation(net, evaluation_data, par)

    return results

if __name__ == '__main__':
  print(par)
  results = main(par) 

  # - formatize results  
  results = dict(all=[results['r2score-aver-all']], grn=[results['r2score-aver-grn']])
  results = pd.DataFrame(results)
  print(results)

  metric_ids = results.columns.to_numpy()
  metric_values = results.values[0]

  print(metric_ids.shape, metric_values.shape)
  results = ad.AnnData(
      X=np.empty((0, 0)),
      uns={
          "dataset_id": par["dataset_id"],
          "method_id": f"reg1-{par['method_id']}",
          "metric_ids": metric_ids,
          "metric_values": metric_values
      }
  )

  results.write_h5ad(par["score"], compression="gzip")