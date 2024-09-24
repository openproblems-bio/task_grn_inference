import pandas as pd
import anndata as ad
import sys
import numpy as np
import os 

## VIASH START
par = {
  'reg_type': 'ridge',
  'read_dir': "resources/grn_models/d0_hvgs",
  'write_dir': "resources/results/scores",
  'methods': [ 'collectri', 'negative_control', 'positive_control', 'pearson_corr', 'portia', 'ppcor', 'genie3', 'grnboost2', 'scenic', 'scglue', 'celloracle'],
  'layers': ['lognorm', 'pearson', 'seurat_lognorm', 'seurat_pearson', 'scgen_lognorm', 'scgen_pearson'],

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
meta = {
  "resources_dir": 'src/metrics/',
  "util": 'src/utils'
}
sys.path.append(meta["resources_dir"])
sys.path.append(meta["util"])

os.makedirs(par['write_dir'], exist_ok=True)

for layer in par['layers']:
  par['layer'] = layer
  for i, method in enumerate(par['methods']):
    par['prediction'] = f"{par['read_dir']}/{method}.csv"
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
    df_all.to_csv(f"{par['write_dir']}/{layer}-{par['reg_type']}.csv")
    print(df_all)
  
