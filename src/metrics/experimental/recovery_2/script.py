import pandas as pd
import anndata as ad
import sys
import numpy as np

## VIASH START
par = {
    'prediction': 'resources/results/300BCG/300BCG.pearson_corr.pearson_corr.prediction.h5ad',
    'evaluation_data': 'resources/grn_benchmark/evaluation_data/300BCG_bulk.h5ad',
    'layer': 'lognorm',
    "max_n_links": 50000,
    'num_workers': 20
}
## VIASH END
try:
  sys.path.append(meta["resources_dir"])
except:
    meta = {
      "resources_dir":'src/metrics/recovery_2/',
      "util_dir":'src/utils'
    }
    sys.path.append(meta["resources_dir"])
    sys.path.append(meta["util_dir"])

from helper import main
from util import format_save_score, parse_args
# from dataset_config import DATASET_GROUPS

par = parse_args(par)


if __name__ == '__main__':
  method_id = ad.read_h5ad(par['prediction'], backed='r').uns['method_id']
  dataset_id = ad.read_h5ad(par['evaluation_data'], backed='r').uns['dataset_id']
  if dataset_id == 'op':
    par['group'] = ['perturbation', 'cell_type']
  elif dataset_id in '300BCG':
    par['group'] = ['cell_type', 'time', 'perturbation']
  else:
    raise ValueError(f"Dataset {dataset_id} not supported yet")
  output = main(par)
  
  format_save_score(output, method_id, dataset_id, par['score'])