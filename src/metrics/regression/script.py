import pandas as pd
import anndata as ad
import sys
import numpy as np
import argparse

## VIASH START
par = {
    'prediction': 'resources/results/300BCG/300BCG.pearson_corr.pearson_corr.prediction.h5ad',
    'evaluation_data': 'resources/grn_benchmark/evaluation_data/300BCG_bulk.h5ad',
    'regulators_consensus': 'resources/grn_benchmark/prior/regulators_consensus_300BCG.json',
    'layer': 'lognorm',
    "max_n_links": 50000,
    'reg_type': 'ridge',
    'static_only': True,
    'subsample': -1,
    'num_workers': 20,
    'apply_tf': True,
    'clip_scores': True,
    'apply_skeleton': False,
    'skeleton': 'resources/grn_benchmark/prior/skeleton.csv',
    'tf_all': 'resources/grn_benchmark/prior/tf_all.csv',
    'verbose': 2
    
}
## VIASH END
try:
  sys.path.append(meta["resources_dir"])
except:
    meta = {
      "resources_dir":'src/metrics/regression/',
      "util_dir":'src/utils'
    }
    sys.path.append(meta["resources_dir"])
    sys.path.append(meta["util_dir"])
    
from helper import main
from util import format_save_score, parse_args

par = parse_args(par)

if __name__ == '__main__':
  print(par)
  detailed_output, output = main(par)
  print(output)
  method_id = ad.read_h5ad(par['prediction'], backed='r').uns['method_id']
  dataset_id = ad.read_h5ad(par['evaluation_data'], backed='r').uns['dataset_id']
  format_save_score(output, method_id, dataset_id, par['score'])

