import pandas as pd
import anndata as ad
import sys
import numpy as np
import os 

## VIASH START
par = {
  'reg_type': 'RF',
  'models_dir': "resources/grn_models/",
  'scores_dir': "resources/scores/",
  
  'methods': [ 'collectri', 'negative_control', 'positive_control', 'pearson_corr', 'portia', 'ppcor', 'grnboost2', 'scenic', 'granie', 'scglue', 'celloracle', 'figr', 'scenicplus'],
  # 'layers': ['scgen_pearson', 'lognorm', 'pearson', 'seurat_lognorm', 'seurat_pearson', 'scgen_lognorm'],
  'layers': ['lognorm', 'pearson'],

  "perturbation_data": "resources/grn-benchmark/perturbation_data.h5ad",
  "tf_all": "resources/prior/tf_all.csv",
  'skeleton': 'resources/prior/skeleton.csv', 
  "max_n_links": 50000,
  "apply_tf": "true",
  'subsample': -1,
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
parser.add_argument('--apply_skeleton', type=str, help='where to apply the skeleton')

args = parser.parse_args()

if args.multiomics_rna:
  par['multiomics_rna'] = args.multiomics_rna
if args.num_workers:
  par['num_workers'] = args.num_workers
if args.models_dir:
  par['models_dir'] = args.models_dir
if args.scores_dir:
  par['scores_dir'] = args.scores_dir
if args.apply_skeleton:
  par['apply_skeleton'] = args.apply_skeleton
    
meta = {
  "resources_dir": 'src/metrics/',
  "util": 'src/utils'
}
sys.path.append(meta["resources_dir"])
sys.path.append(meta["util"])

os.makedirs(par['scores_dir'], exist_ok=True)

for binarize in [True]:
  par['binarize'] = binarize
  for max_n_links in [50000]:
    par['max_n_links'] = max_n_links
    for apply_skeleton in [True]:
      par['apply_skeleton'] = apply_skeleton
      for layer in par['layers']:
        par['layer'] = layer
        i = 0
        for method in par['methods']:
          print(method)
          par['prediction'] = f"{par['models_dir']}/{method}.csv"
          if not os.path.exists(par['prediction']):
            print(f"{par['prediction']} doesnt exist. Skipped.")
            continue
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
          df_all.to_csv(f"{par['scores_dir']}/{max_n_links}-skeleton_{apply_skeleton}-binarize_{binarize}_{layer}-{par['reg_type']}.csv")
          print(df_all)
          i+=1
  
