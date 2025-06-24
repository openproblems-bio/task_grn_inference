import pandas as pd
import anndata as ad
import sys
import numpy as np
from sklearn.metrics import r2_score
from sklearn.pipeline import Pipeline
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler
import random 
from concurrent.futures import ThreadPoolExecutor
from tqdm import tqdm
from sklearn.multioutput import MultiOutputRegressor
import os
import warnings
warnings.simplefilter("ignore")

## VIASH START
par = {
  "evaluation_data": f"resources/grn_benchmark/evaluation_data/op_bulk.h5ad",
  "tf_all": "resources/grn_benchmark/prior/tf_all.csv",
  "prediction": "resources/grn_models/op/grnboost2.h5ad",
  "max_n_links": 50000,
  "apply_tf": True,
  'score': 'output/score.h5ad',
  'reg_type': 'ridge',
  'layer': 'X_norm',
  'num_workers': 1,
  'skeleton': 'resources/grn_benchmark/prior/skeleton.csv',
  'apply_skeleton': False,
  'verbose': 4,
  'binarize': False
}
## VIASH END

## LOCAL START
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--run_local', action='store_true', help='Run locally')
parser.add_argument('--evaluation_data', type=str, help='Path to the evaluation_data file')
parser.add_argument('--prediction', type=str, help='Path to the prediction file')
parser.add_argument('--dataset_id', type=str, help='Dataset id')
parser.add_argument('--score', type=str, help='score file')
parser.add_argument('--reg_type', type=str)
parser.add_argument('--apply_skeleton', action='store_true')

args = parser.parse_args()


var_local = vars(args)

if args.run_local:
    for key in var_local:
        if var_local[key] is not None:
            par[key] = var_local[key]

## LOCAL END
try:
    sys.path.append(meta["resources_dir"])
except:
    meta = {
      "resources_dir":'src/metrics/regression_1/',
      "util_dir":'src/utils'
    }
    sys.path.append(meta["resources_dir"])
    sys.path.append(meta["util_dir"])
  

from util import verbose_print, process_links, verbose_tqdm, binarize_weight
from helper import set_global_seed, process_net, cross_validation

def read_net(par):
    print(par['prediction'], flush=True)
    net = ad.read_h5ad(par['prediction'])
    net = pd.DataFrame(net.uns['prediction'])
    net['weight'] = net['weight'].astype(float)
    if par['binarize']:
      net['weight'] = net['weight'].apply(binarize_weight)  
    if par['apply_skeleton']: 
        print('Before filtering with skeleton:', net.shape)
        skeleton = pd.read_csv(par['skeleton'])
        skeleton['link'] = skeleton['source'].astype(str) + '_' + skeleton['target'].astype(str)
        skeleton = skeleton['link'].values.flatten()
        
        net['link'] = net['source'].astype(str) + '_' + net['target'].astype(str)
        net = net[net['link'].isin(skeleton)]
        print('After filtering with skeleton:', net.shape)

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
    tf_all = np.loadtxt(par['tf_all'], dtype=str)
    net = read_net(par)
    if par['apply_tf']:
        net = net[net.source.isin(tf_all)]
    
    assert net.shape[0]>0, 'No links found in the network'
    assert all(column in net.columns for column in ['source', 'target', 'weight']), 'Columns in the network should be source, target, weight'
    
    evaluation_data = ad.read_h5ad(par['evaluation_data'])
    evaluation_data.X = evaluation_data.layers[par["layer"]]
    gene_names = evaluation_data.var.index.to_numpy()
    
    # -- calclate scores 
    results = cross_validation(net, evaluation_data, par)

    return results

if __name__ == '__main__':
    method_id = ad.read_h5ad(par['prediction'], backed='r').uns['method_id']
    dataset_id = ad.read_h5ad(par['evaluation_data'], backed='r').uns['dataset_id']
    print(f"Method id: {method_id}, Dataset id: {dataset_id}")

    # - run main function
    results = main(par) 
    # - formatize results  
    results = dict(r1_all=[results['r2score-aver-all']], r1_grn=[results['r2score-aver-grn']])
    results = pd.DataFrame(results)
    print(results)

    
    results = ad.AnnData(
        X=np.empty((0, 0)),
        uns={
            "dataset_id": dataset_id,
            "method_id": method_id,
            "metric_ids": results.columns.to_numpy(),
            "metric_values": results.values[0]
        }
    )
    print(results.uns)

    results.write_h5ad(par["score"], compression="gzip")