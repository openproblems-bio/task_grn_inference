import pandas as pd
import anndata as ad
import sys
import numpy as np
import os 
import random

par = {
  'reg_type': 'ridge',
  'read_dir': "resources/grn_models/op/",
  'write_dir': "resources/results/robustness_analysis",
  'degrees': [0, 10, 20, 50, 100],
  # 'degrees': [50],
  'noise_types': ["net", "sign"],
  # 'noise_types': ['weight'],
  'methods': ['negative_control', 'positive_control', 'pearson_corr', 'portia', 'ppcor', 'grnboost2', 'scenic', 'granie', 'scglue', 'celloracle', 'figr', 'scenicplus'],

  "evaluation_data": "resources/evaluation_datasets/op_perturbation.h5ad",
  "tf_all": "resources/prior/tf_all.csv",
  "max_n_links": 50000,
  "apply_tf": True,
  'binarize': False, 
  'subsample': -1,
  'verbose': 0,
  'num_workers': 20,
  'consensus': 'resources/prior/op_consensus-num-regulators.json',
  'static_only': True,
  'layer': 'X_norm',
  'apply_skeleton': False,
  'skeleton': 'resources/prior/skeleton.csv'
}

meta = {
  "resources_dir": 'src/',
  "util": 'src/utils'
}
sys.path.append(meta["resources_dir"])
sys.path.append(meta["util"])

os.makedirs(par['write_dir'], exist_ok=True)
os.makedirs(f"{par['write_dir']}/tmp/", exist_ok=True)

def run_reg(par):
  from metrics.regression_1.main import main 
  reg1 = main(par)
  from metrics.regression_2.main import main 
  reg2 = main(par)
  score = pd.concat([reg1, reg2], axis=1)
  return score 
  
#------ noise types and degrees ------#
if True:
  for noise_type in par['noise_types']: # run for each noise type (net, sign, weight)
    for degree in par['degrees']: # run for each degree
      for i, method in enumerate(par['methods']): # run for each method
        par['prediction'] = f"{par['read_dir']}/{method}.csv"
        par['prediction_n'] = f"{par['write_dir']}/tmp/{method}.csv"
        par['degree'] = degree
        par['noise_type'] = noise_type
        
        # permute
        from robustness_analysis.permute_grn.main import main 
        prediction_n = main(par)
        prediction_n.to_csv(par['prediction_n'])
        
        # run regs 
        par['prediction'] = par['prediction_n']
        score = run_reg(par)
        score.index = [method]
        if i==0:
          df_all = score
        else:
          df_all = pd.concat([df_all, score])
        df_all.to_csv(f"{par['write_dir']}/{noise_type}-{degree}-scores.csv")
        print(df_all)

#------ causal vs corr ------#
if False:
  from util import create_corr_net
  par = {
    'reg_type': 'ridge',
    'write_dir': "resources/results/robustness_analysis",
    ## base corr
    "perturbation_data": "resources/grn-benchmark/perturbation_data.h5ad",
    'cell_type_specific': False,
    'normalize': False,
    ## metric 
    'multiomics_rna': 'resources/grn-benchmark/multiomics_rna_d0_hvg.h5ad',
    "tf_all": "resources/prior/tf_all.csv",
    "max_n_links": 50000,
    "apply_tf": False, #this has to be false
    'subsample': -2,
    'verbose': 2,
    'binarize': True,
    'num_workers': 20,
    'consensus': 'resources/prior/consensus-num-regulators.json',
    'static_only': True,
    'clip_scores': True,
    'layer': 'scgen_pearson',
    'seed': 32
  }

  # run for corr 
  os.makedirs(f"{par['write_dir']}/corr/", exist_ok=True)
  par['causal'] = False
  for i in range(100):
    par['causal']
    par['prediction'] = f"{par['write_dir']}/corr/corr.csv"
    par['seed'] = i
    random.seed(par['seed'])
    print('seed :', par['seed'])
    
    net = create_corr_net(par)
    net.to_csv(par['prediction'])
    score = run_reg(par)
    if i == 0:
      scores_corr = score 
    else:
      scores_corr = pd.concat([score, scores_corr], axis=0)
    print(scores_corr)
    scores_corr.to_csv(f"{par['write_dir']}/corr/scores_corr.csv")
    
  # run for causal corr
  par['prediction'] = f"{par['write_dir']}/corr/corr_causal.csv"
  par['causal'] = True
  net = create_corr_net(par)

  net.to_csv(par['prediction'])
  score = run_reg(par)
  score.to_csv(f"{par['write_dir']}/corr/scores_causal.csv")