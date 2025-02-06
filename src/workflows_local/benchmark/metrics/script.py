import pandas as pd
import anndata as ad
import sys
import numpy as np
import os 
import subprocess

import argparse
argparser = argparse.ArgumentParser()

argparser.add_argument('--datasets', nargs='+', help='List of datasets to include', required=True)
argparser.add_argument('--methods', nargs='+', help='List of methods to include', required=True)
argparser.add_argument('--run_consensus_flag', action='store_true', help='Run consensus')
argparser.add_argument('--save_scores_file',  help='Save file name', required=True)
args = argparser.parse_args()
par = vars(args)

par['temp_dir'] = 'output/temp'
os.makedirs(par['temp_dir'], exist_ok=True)

meta = {
  "resources_dir": './',
  "util": 'src/utils'
}
sys.path.append(meta["resources_dir"])
sys.path.append(meta["util"])

dependencies={
  'all_metrics': 'src/metrics/all_metrics/script.py'
}

# - run consensus 
from src.metrics.regression_2.consensus.script import main as main_consensus


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

  par_local = {
      "evaluation_dataset": f"resources/grn_benchmark/evaluation_datasetsets/{dataset}_bulk.h5ad",
      "evaluation_dataset_sc": f"resources/grn_benchmark/evaluation_datasetsets/{dataset}_sc.h5ad",
      'regulators_consensus':  f'resources/grn_benchmark/prior/regulators_consensus_{dataset}.json',
      'ws_consensus': f'resources/grn_benchmark/prior/ws_consensus_{dataset}.json',
      'ws_distance_background': f'resources/grn_benchmark/prior/ws_distance_background_{dataset}.json',

      'layer': 'X_norm',
      "tf_all": "resources/grn_benchmark/prior/tf_all.csv",
      'skeleton': 'resources/grn_benchmark/prior/skeleton.csv', 
      "apply_tf": True,
      'subsample': -1,
      'verbose': 4,
      'num_workers': 20
  }

  return par_local

def run_consensus(dataset):
  par_local = define_par(dataset)
  main_consensus(par_local)

def run_metrics(par):
    par['score'] = f"output/score_{par['dataset_id']}_{par['method_id']}.h5ad"
    args = f"--run_local --prediction {par['prediction']} \
            --dataset_id {par['dataset_id']} \
            --method_id {par['method_id']} \
            --evaluation_dataset {par['evaluation_dataset']} \
            --evaluation_dataset_sc {par['evaluation_dataset_sc']} \
            --regulators_consensus {par['regulators_consensus']} \
            --ws_consensus {par['ws_consensus']} \
            --ws_distance_background {par['ws_distance_background']} \
            --score {par['score']} "
    

    command = f"python {dependencies['all_metrics']} {args}"

    print(command)

    result = subprocess.run(command, shell=True, capture_output=True, text=True)

    if result.returncode != 0:
        print("STDOUT:", result.stdout)
        print("STDERR:", result.stderr)
        raise RuntimeError(f"Error: command proprocess failed with exit code {result.returncode}")
    print('preprocess completed')

    score = ad.read_h5ad(par['score']).uns['score']
  
    return score


def run_evaluation(dataset, binarize=False, max_n_links=50000, apply_skeleton=False,  reg_type='ridge'):
  print('------ ', dataset, '------')
  par_local = define_par(dataset)
  par.update(par_local)
  
  par['binarize'] = binarize
  par['max_n_links'] = max_n_links
  par['apply_skeleton'] = apply_skeleton
  par['reg_type'] = reg_type
      
  # - determines models to run 
  grn_files_dict = {}
  # - add models
  for model in par['methods']:
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
    par['method_id'] = model
    par['dataset_id'] = dataset

  
    score = run_metrics(par)


    score.index = [model]
    if i==0:
      df_all = score
    else:
      df_all = pd.concat([df_all, score])
    print(df_all)
    
    i+=1
  return df_all 


if __name__ == '__main__':
  # - define settings
  # models = ['negative_control', 'positive_control', 'pearson_corr', 'portia', 'ppcor', 'grnboost2', 'scenic', 'granie', 'scglue', 'celloracle', 'figr', 'scenicplus']

  for i, dataset in enumerate(par['datasets']): 
    models_dir = f"resources/grn_models/{dataset}" 
    if par['run_consensus_flag']:
      run_consensus(dataset)
    
    df_dataset = run_evaluation(dataset)
    df_dataset['dataset'] = dataset
    if i==0:
      df_all = df_dataset
    else:
      df_all = pd.concat([df_all, df_dataset])
  df_all.to_csv(par['save_scores_file'])
  

  # if False: # subsample
  #   # for dataset in ['op', 'replogle2', 'nakatake', 'norman', 'adamson']: #'op', 'replogle2', 'nakatake', 'norman', 'adamson'
  #   save_scores_file = f"{scores_dir}/subsampled.csv" 
  #   for i, dataset in enumerate(['op']):
  #     if dataset == 'op':
  #       models_subsampled = [f'{model}_{subsample}' for subsample in [1, 2] for model in models]
  #     else:
  #       models_subsampled = [f'{model}_{subsample}' for subsample in [0.2, 0.5] for model in models]
  #     models_dir = f"resources/grn_models/{dataset}"
      
    
  #     df_dataset = run_evaluation(dataset, models, models_di)
  #     df_dataset['dataset'] = dataset
  #     if i==0:
  #       df_all = df_dataset
  #     else:
  #       df_all = pd.concat([df_all, df_dataset])
  #   df_all.to_csv(save_scores_file)


  # if False: # run global models
  #   models = ['pearson_corr']
  #   dataset = 'op'

  #   models_dir = "resources/grn_models/global/"  
  #   scores_dir = f"resources/scores/{dataset}"
  #   # run_consensus(dataset)
  #   save_scores_file = f"{scores_dir}/X_norm-50000-skeleton_False-binarize_False-ridge-global-True.csv" 

  #   run_evaluation(dataset, models, models_dir, scores_dir, save_scores_file)
  
  # if False: # run skeleton 
  #   models = ['negative_control', 'positive_control', 'pearson_corr', 'portia', 'ppcor', 'grnboost2', 'scenic', 'granie', 'scglue', 'celloracle', 'figr', 'scenicplus']

  #   dataset = 'op'

  #   models_dir = f"resources/grn_models/{dataset}"
  #   scores_dir = f"resources/scores/{dataset}"
  #   save_scores_file = f"{scores_dir}/X_norm-50000-skeleton_True-binarize_False-ridge-global-False.csv" 

  #   # run_consensus(dataset)
  #   run_evaluation(dataset, models, models_dir, scores_dir, save_scores_file, apply_skeleton=True)
    
  # if False: # run GB 
  #   models = ['negative_control', 'positive_control', 'pearson_corr', 'portia', 'ppcor', 'grnboost2', 'scenic', 'granie', 'scglue', 'celloracle', 'figr', 'scenicplus']

  #   dataset = 'op'

  #   models_dir = f"resources/grn_models/{dataset}"
  #   scores_dir = f"resources/scores/{dataset}"
  #   save_scores_file = f"{scores_dir}/X_norm-50000-skeleton_True-binarize_False-GB-global-False.csv" 

  #   # run_consensus(dataset)
  #   run_evaluation(dataset, models, models_dir, scores_dir, save_scores_file, apply_skeleton=True, reg_type='GB')
    
  
    

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