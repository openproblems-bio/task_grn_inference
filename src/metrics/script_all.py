import pandas as pd
import anndata as ad
import sys
import numpy as np
import os 

## VIASH START
par = {
  'reg_type': 'ridge',
  'models_dir': "resources/grn_models/replogle2",
  'scores_dir': "output/",
  
  # 'models': [ 'collectri', 'negative_control', 'positive_control', 'pearson_corr', 'portia', 'ppcor', 'grnboost2', 'scenic', 'granie', 'scglue', 'celloracle', 'figr', 'scenicplus'],
  'models': [ 'negative_control', 'positive_control', 'pearson_corr', 'portia', 'grnboost2', 'ppcor', 'scenic'],
  # 'models': [ 'negative_control', 'positive_control', 'pearson_corr', 'portia'],

  "evaluation_data": "resources/evaluation_datasets/replogle2.h5ad",
  'consensus': 'resources/prior/replogle2_consensus-num-regulators.json',
  # "evaluation_data": "resources/evaluation_datasets/frangieh_IFNg_v2.h5ad",
  # 'consensus': 'resources/prior/frangieh_IFNg_v2_consensus-num-regulators.json',

  'layers': ['X'],
  
  "tf_all": "resources/prior/tf_all.csv",
  'skeleton': 'resources/prior/skeleton.csv', 
  "apply_tf": True,
  'subsample': -1,
  'verbose': 4,
  'num_workers': 20,
  'clip_scores': True
}
# VIASH END

meta = {
  "resources_dir": 'src/metrics/',
  "util": 'src/utils'
}
sys.path.append(meta["resources_dir"])
sys.path.append(meta["util"])

os.makedirs(par['scores_dir'], exist_ok=True)

# - run consensus 
from consensus.script import main 
main(par)

# - run metrics 
for binarize in [False]:
  par['binarize'] = binarize
  for max_n_links in [50000]:
    par['max_n_links'] = max_n_links
    for apply_skeleton in [False]:
      par['apply_skeleton'] = apply_skeleton
      for layer in par['layers']:
        par['layer'] = layer
        i = 0
        for model in par['models']:
          print(model)
          par['prediction'] = f"{par['models_dir']}/{model}.csv"
          if not os.path.exists(par['prediction']):
            print(f"{par['prediction']} doesnt exist. Skipped.")
            continue
          from regression_1.main import main 
          reg1 = main(par)
          from regression_2.main import main 
          reg2 = main(par)
          score = pd.concat([reg1, reg2], axis=1)
          score.index = [model]
          if i==0:
            df_all = score
          else:
            df_all = pd.concat([df_all, score])
          df_all.to_csv(f"{par['scores_dir']}/{max_n_links}-skeleton_{apply_skeleton}-binarize_{binarize}_{layer}-{par['reg_type']}.csv")
          print(df_all)
          i+=1
  
