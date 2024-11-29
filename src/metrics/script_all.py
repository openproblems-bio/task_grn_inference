import pandas as pd
import anndata as ad
import sys
import numpy as np
import os 


def define_par(dataset):

  par = {
      'reg_type': 'ridge',
      'models_dir': f"resources/grn_models/{dataset}",
      'scores_dir': f"resources/scores/{dataset}",
      
      'models': [ 'negative_control', 'positive_control', 'pearson_corr', 'portia', 'ppcor', 'grnboost2', 'scenic', 'granie', 'scglue', 'celloracle', 'figr', 'scenicplus'],

      'global_models': [
                        'collectri',
                        'Ananse:Lung',
                        'Ananse:Stomach',
                        'Ananse:Heart',
                        'Ananse:Bone marrow',
                        'Gtex:Whole blood',
                        'Gtex:Brain amygdala',
                        'Gtex:Breast mammary tissue',
                        'Gtex:Lung',
                        'Gtex:Stomach',
                        'Cellnet:Bcell',
                        'Cellnet:Tcell',
                        'Cellnet:Skin',
                        'Cellnet:Neuron',
                        'Cellnet:Heart'
                        ],
      'global_models_dir': 'resources/grn_models/global/',

      "evaluation_data": f"resources/evaluation_datasets/{dataset}_perturbation.h5ad",
      'consensus':  f'resources/prior/{dataset}_consensus-num-regulators.json',

      'layer': 'X_norm',
      
      "tf_all": "resources/prior/tf_all.csv",
      'skeleton': 'resources/prior/skeleton.csv', 
      "apply_tf": True,
      'subsample': -1,
      'verbose': 4,
      'num_workers': 20
  }

  return par


meta = {
  "resources_dir": 'src/metrics/',
  "util": 'src/utils'
}
sys.path.append(meta["resources_dir"])
sys.path.append(meta["util"])

from regression_1.main import main as main_reg1
from regression_1.main import binarize_weight
from regression_2.main import main as main_reg2
from util import process_links
# - run consensus 
from consensus.script import main as main_consensus

# - run general models
global_models = False

# - run metrics 
for dataset in ['op']: #'op', 'replogle2', 'nakatake', 'norman', 'adamson'
  print('------ ', dataset, '------')
  par = define_par(dataset)
  os.makedirs(par['scores_dir'], exist_ok=True)
  main_consensus(par)
  if global_models:
    par['models'] = par['global_models']
    par['models_dir'] = par['global_models_dir']
  for binarize in [False]:
    par['binarize'] = binarize
    for max_n_links in [50000]:
      par['max_n_links'] = max_n_links
      for apply_skeleton in [False]:
        par['apply_skeleton'] = apply_skeleton
        # - determines models to run 
        grn_files_dict = {}
        # - add models
        for model in par['models']:
          print(model)
          grn_file = f"{par['models_dir']}/{model}.csv"
          if not os.path.exists(grn_file):
            print(f"{grn_file} doesnt exist. Skipped.")
            continue
          grn_files_dict[model] = grn_file
          
        # - actual runs 
        i = 0
        for model, grn_file in grn_files_dict.items():
          par['prediction'] = grn_file
          reg1 = main_reg1(par)
          reg2 = main_reg2(par)
          score = pd.concat([reg1, reg2], axis=1)
          score.index = [model]
          if i==0:
            df_all = score
          else:
            df_all = pd.concat([df_all, score])
          df_all.to_csv(f"{par['scores_dir']}/{par['layer']}-{max_n_links}-skeleton_{apply_skeleton}-binarize_{binarize}-{par['reg_type']}-global-{global_models}.csv")
          print(df_all)
          i+=1
  
