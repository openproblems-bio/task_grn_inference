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

par = parse_args(par)



DATASET_GROUPS = {
    "op": ["perturbation", "cell_type"],
    "parsebioscience": ["perturbation", "cell_type"],
    "300BCG": ["perturbation", "cell_type"]
}
if __name__ == '__main__':
    par['group'] = DATASET_GROUPS[ad.read_h5ad(par['evaluation_data'], backed='r').uns['dataset_id']]
    output = main(par)
    format_save_score(output, par)