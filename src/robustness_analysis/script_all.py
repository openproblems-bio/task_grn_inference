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
  # 'degrees': [20, 50, 100],
  # 'degrees': [50],
  # 'analysis_types': ["net", "sign", 'weight'],
  'analysis_types': ['direction'],
  'methods': ['negative_control', 'positive_control', 'pearson_corr', 'portia', 'ppcor', 'grnboost2', 'scenic', 'granie', 'scglue', 'celloracle', 'figr', 'scenicplus'],
  # 'methods': ['pearson_corr'],

  "evaluation_data": "resources/grn_benchmark/evaluation_datasets//op_perturbation.h5ad",
  "tf_all": "resources/grn_benchmark/prior/tf_all.csv",
  "max_n_links": 50000,
  "apply_tf": False,
  'binarize': False, 
  'subsample': -1,
  'verbose': 0,
  'num_workers': 20,
  'consensus': 'resources/grn_benchmark/prior/op_consensus-num-regulators.json',
  'static_only': True,
  'layer': 'X_norm',
  'apply_skeleton': False,
  'skeleton': 'resources/grn_benchmark/prior/skeleton.csv'
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
  for noise_type in par['analysis_types']: # run for each noise type (net, sign, weight)
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
        print(noise_type, degree, df_all)
        df_all.to_csv(f"{par['write_dir']}/{noise_type}-{degree}-scores.csv")
