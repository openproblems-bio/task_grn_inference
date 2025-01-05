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

from regression_1.main import main as main_reg1
from regression_1.main import binarize_weight
from regression_2.main import main as main_reg2
from util import process_links
# - run consensus 
from consensus.script import main as main_consensus

global_models = [
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



def define_par(dataset):

  par = {
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

def run_consensus(dataset):
  par = define_par(dataset)
  main_consensus(par)

def run_evaluation(dataset, models, models_dir, scores_dir, save_file_name, binarize=False, max_n_links=50000, apply_skeleton=False, run_global_models=False, reg_type='ridge'):
  print('------ ', dataset, '------')
  par = define_par(dataset)
  os.makedirs(scores_dir, exist_ok=True)
  
  par['binarize'] = binarize
  par['max_n_links'] = max_n_links
  par['apply_skeleton'] = apply_skeleton
  par['reg_type'] = reg_type
      
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
  # - define settings
  # models = ['negative_control', 'positive_control', 'pearson_corr', 'portia', 'ppcor', 'grnboost2', 'scenic', 'granie', 'scglue', 'celloracle', 'figr', 'scenicplus']
  models = ['pearson_corr', 'portia', 'grnboost2']

  if True: # subsample
    # for dataset in ['op', 'replogle2', 'nakatake', 'norman', 'adamson']: #'op', 'replogle2', 'nakatake', 'norman', 'adamson'
    for dataset in ['op']:
      if dataset == 'op':
        models_subsampled = [f'{model}_{subsample}' for subsample in [1, 2] for model in models]
      else:
        models_subsampled = [f'{model}_{subsample}' for subsample in [0.2, 0.5] for model in models]
      models_dir = f"resources/grn_models/{dataset}"
      scores_dir = f"resources/scores/{dataset}"
      
      save_file_name = f"{scores_dir}/subsampled.csv" 

      run_evaluation(dataset, models_subsampled, models_dir, scores_dir, save_file_name)
  

  if False: # default run 
    models = ['negative_control', 'positive_control', 'pearson_corr', 'portia', 'ppcor', 'grnboost2', 'scenic', 'granie', 'scglue', 'celloracle', 'figr', 'scenicplus']

    for dataset in ['op', 'replogle2', 'nakatake', 'norman', 'adamson']: 
      models_dir = f"resources/grn_models/{dataset}" 
      scores_dir = f"resources/scores/{dataset}"
      run_consensus(dataset)
      save_file_name = f"{scores_dir}/X_norm-50000-skeleton_False-binarize_False-ridge-global-False.csv" 

      run_evaluation(dataset, models, models_dir, scores_dir, run_global_models, save_file_name)
  
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
    
  
    

# - repo
# models_dir = f"../ciim/output/grns/"
  # scores_dir = f"../ciim/output/scores/"
  # models = [
  #       'pearson_corr', 
  #       'grnboost2', 'celloracle', 'scenicplus',
  #       'net_all_celltypes_young_all_batches',
  #       'net_all_celltypes_young_batch_1',
  #       'net_all_celltypes_old_all_batches',
  #       'net_all_celltypes_old_batch_1', 
  #       'net_B cells_all_ages_all_batches', 
  #       'net_T cells_all_ages_all_batches', 
  #       'net_Myeloid cells_all_ages_all_batches', 
  #       'net_NK cells_all_ages_all_batches'],