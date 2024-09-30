import pandas as pd
import anndata as ad
import sys
import numpy as np
import os 

## VIASH START
par = {
  'reg_type': 'ridge',
  'models_dir': "resources/grn_models/",
  'scores_dir': "resources/scores/",
  
  'methods': [ 'collectri', 'negative_control', 'positive_control', 'pearson_corr', 'portia', 'ppcor', 'grnboost2', 'scenic', 'granie', 'scglue', 'celloracle'],
  # 'layers': ['scgen_pearson', 'lognorm', 'pearson', 'seurat_lognorm', 'seurat_pearson', 'scgen_lognorm'],
  'layers': ['pearson', 'seurat_lognorm', 'seurat_pearson', 'scgen_lognorm'],

  "perturbation_data": "resources/grn-benchmark/perturbation_data.h5ad",
  "tf_all": "resources/prior/tf_all.csv",
  "max_n_links": 50000,
  "apply_tf": "true",
  'subsample': -2,
  'verbose': 1,
  'binarize': True,
  'num_workers': 20,
  'consensus': 'resources/prior/consensus-num-regulators.json',
  'static_only': True,
  'clip_scores': True
}
# VIASH END

import argparse
parser = argparse.ArgumentParser(description="Process multiomics RNA data.")
parser.add_argument('--multiomics_rna', type=str, help='Path to the multiomics RNA file')
parser.add_argument('--num_workers', type=str, help='Number of cores')
parser.add_argument('--models_dir', type=str, help='where to read the models')
parser.add_argument('--scores_dir', type=str, help='where to write the model')

args = parser.parse_args()

if args.multiomics_rna:
  par['multiomics_rna'] = args.multiomics_rna
if args.num_workers:
  par['num_workers'] = args.num_workers
if args.models_dir:
  par['models_dir'] = args.models_dir
if args.scores_dir:
  par['scores_dir'] = args.scores_dir
    
meta = {
  "resources_dir": 'src/metrics/',
  "util": 'src/utils'
}
sys.path.append(meta["resources_dir"])
sys.path.append(meta["util"])

os.makedirs(par['scores_dir'], exist_ok=True)


for layer in par['layers']:
  par['layer'] = layer
  for i, method in enumerate(par['methods']):
    par['prediction'] = f"{par['models_dir']}/{method}.csv"
    from regression_1.main import main 
    reg1 = main(par)
    from regression_2.main import main 
    reg2 = main(par)
    score = pd.concat([reg1, reg2], axis=1)
    score.index = [method]
    if i==0:
      df_all = score
    else:
      df_all = pd.concat([df_all, score])
    df_all.to_csv(f"{par['scores_dir']}/{layer}-{par['reg_type']}.csv")
    print(df_all)
  
