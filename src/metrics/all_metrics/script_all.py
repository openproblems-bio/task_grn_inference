import pandas as pd
import anndata as ad
import sys
import numpy as np
import os 

meta = {
  "resources_dir": 'src/metrics/',
  "util": 'src/utils'
}
sys.path.append(meta["resources_dir"])
sys.path.append(meta["util"])

from all_metrics.helper import run_consensus, run_ws_distance_background, run_scores_all

par = {
      'layer': 'X_norm',
      "tf_all": "resources/grn_benchmark/prior/tf_all.csv",
      'skeleton': 'resources/grn_benchmark/prior/skeleton.csv', 
      "apply_tf": True,
      'subsample': -1,
      'verbose': 4,
      'num_workers': 20,
      'binarize': False, 
      'max_n_links': 50000, 
      'apply_skeleton': False, 
      'reg_type':'ridge'
  }


def run_evaluation(dataset, models, models_dir, save_file_name):
  print('------ ', dataset, '------')
   
  # - determines models to run 
  grn_files_dict = {}
  # - add models
  for model in models:
    print(model)
    grn_file = f"{models_dir}/{model}.csv"
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
    df_all.to_csv(save_file_name)
    print(df_all)
    i+=1

if __name__ == '__main__':
  run_scores_flag = True
  run_consensus_flag = False 
  run_ws_distance_background_flag = False 
  datasets = ['op', 'replogle', 'nakatake', 'norman', 'adamson']

  if run_consensus_flag: # run consensus
    run_consensus(datasets)
  
  if run_ws_distance_background_flag: # run background scores for ws distance
    run_ws_distance_background(datasets)
  
  if run_scores_flag:
    models = ['negative_control', 'positive_control', 'pearson_corr', 'portia', 'ppcor', 'grnboost2', 'scenic', 'granie', 'scglue', 'celloracle', 'figr', 'scenicplus']

    run_scores_all(datasets, models=models)


  aaa

  if False: # default run 
    for dataset in dataset: 
      models_dir = f"resources/grn_models/{dataset}" 
      scores_dir = f"resources/scores/{dataset}"
      run_consensus(dataset)
      save_file_name = f"{scores_dir}/default_scores.csv" 

      run_evaluation(dataset, models, models_dir, scores_dir, save_file_name)

  if True: # subsample
    # for dataset in ['op', 'replogle', 'nakatake', 'norman', 'adamson']: #'op', 'replogle', 'nakatake', 'norman', 'adamson'
    for dataset in ['op']:
      if dataset == 'op':
        models_subsampled = [f'{model}_{subsample}' for subsample in [1, 2] for model in models]
      else:
        models_subsampled = [f'{model}_{subsample}' for subsample in [0.2, 0.5] for model in models]
      models_dir = f"resources/grn_models/{dataset}"
      scores_dir = f"resources/scores/{dataset}"
      
      save_file_name = f"{scores_dir}/subsampled.csv" 

      run_evaluation(dataset, models_subsampled, models_dir, scores_dir, save_file_name)
  

  
  if False: # run global models
    models = ['pearson_corr']
    dataset = 'op'

    models_dir = "resources/grn_models/global/"  
    scores_dir = f"resources/scores/{dataset}"
    # run_consensus(dataset)
    save_file_name = f"{scores_dir}/X_norm-50000-skeleton_False-binarize_False-ridge-global-True.csv" 

    run_evaluation(dataset, models, models_dir, scores_dir, run_global_models, save_file_name)
  
  if False: # run skeleton 
    models = ['negative_control', 'positive_control', 'pearson_corr', 'portia', 'ppcor', 'grnboost2', 'scenic', 'granie', 'scglue', 'celloracle', 'figr', 'scenicplus']

    dataset = 'op'

    models_dir = f"resources/grn_models/{dataset}"
    scores_dir = f"resources/scores/{dataset}"
    save_file_name = f"{scores_dir}/X_norm-50000-skeleton_True-binarize_False-ridge-global-False.csv" 

    # run_consensus(dataset)
    run_evaluation(dataset, models, models_dir, scores_dir, save_file_name, apply_skeleton=True)
    
  if False: # run GB 
    models = ['negative_control', 'positive_control', 'pearson_corr', 'portia', 'ppcor', 'grnboost2', 'scenic', 'granie', 'scglue', 'celloracle', 'figr', 'scenicplus']

    dataset = 'op'

    models_dir = f"resources/grn_models/{dataset}"
    scores_dir = f"resources/scores/{dataset}"
    save_file_name = f"{scores_dir}/X_norm-50000-skeleton_True-binarize_False-GB-global-False.csv" 

    # run_consensus(dataset)
    run_evaluation(dataset, models, models_dir, scores_dir, save_file_name, apply_skeleton=True, reg_type='GB')
    
  
    



# def define_par(dataset):

#   par = {
#       "evaluation_data": f"resources/grn_benchmark/evaluation_data//{dataset}.h5ad",
#       'consensus':  f'resources/grn_benchmark/prior/{dataset}_consensus-num-regulators.json',

#       'layer': 'X_norm',
      
#       "tf_all": "resources/grn_benchmark/prior/tf_all.csv",
#       'skeleton': 'resources/grn_benchmark/prior/skeleton.csv', 
#       "apply_tf": True,
#       'subsample': -1,
#       'verbose': 4,
#       'num_workers': 20
#   }

#   return par
