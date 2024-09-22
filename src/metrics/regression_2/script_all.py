import pandas as pd
import anndata as ad
import sys
import numpy as np
import os 

## VIASH START
par = {
  'reg_type': 'ridge',
  'read_dir': "resources/grn_models/d0_hvgs",
  'write_dir': "resources/results/d0_hvgs_ridge",
  'methods': [ 'collectri', 'negative_control', 'positive_control', 'pearson_corr', 'pearson_causal',  'portia', 'ppcor', 'genie3', 'grnboost2', 'scenic', 'scglue', 'celloracle'],


  "perturbation_data": "resources/grn-benchmark/perturbation_data.h5ad",
  "tf_all": "resources/prior/tf_all.csv",
  "min_tf": False,
  "max_n_links": 50000,
  "apply_tf": "true",
  'layer': 'scgen_pearson',
  'subsample': -2,
  'num_workers': 4,
  'verbose': 1,
  'binarize': True,
  'num_workers': 20,
  'consensus': 'resources/prior/consensus-num-regulators.json',
  'static_only': True,
  'clip_scores': True
}
# VIASH END
meta = {
  "resources_dir": 'src/metrics/regression_2/',
  "util": 'src/utils'
}
sys.path.append(meta["resources_dir"])
sys.path.append(meta["util"])
from main import main 


os.makedirs(par['write_dir'], exist_ok=True)


for i, method in enumerate(par['methods']):
  par['prediction'] = f"{par['read_dir']}/{method}.csv"
  prediction = main(par)
  prediction.index = [method]
  if i==0:
    df_all = prediction
  else:
    df_all = pd.concat([df_all, prediction])
  df_all.to_csv(f"{par['write_dir']}/reg2.csv")
  print(df_all)
  
