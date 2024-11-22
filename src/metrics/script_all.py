import pandas as pd
import anndata as ad
import sys
import numpy as np
import os 


def define_par(dataset):

  par = {
      'reg_type': 'ridge',
      'models_dir': f"resources/grn_models/{dataset}",
      'scores_dir': f"output/temp/{dataset}",
      
      'models': [ 'collectri', 'negative_control', 'positive_control', 'pearson_corr', 'portia', 'ppcor', 'grnboost2', 'scenic', 'granie', 'scglue', 'celloracle', 'figr', 'scenicplus'],

      # 'models': [ 'positive_control', 'pearson_corr'],
      'global_models': [
                        'ANANSE_tissue/networks/lung.parquet',
                        'ANANSE_tissue/networks/stomach.parquet', 
                        'ANANSE_tissue/networks/heart.parquet',
                        'ANANSE_tissue/networks/bone_marrow.parquet',
                        
                        'gtex_rna/networks/Whole_Blood.parquet',
                        'gtex_rna/networks/Brain_Amygdala.parquet', 
                        'gtex_rna/networks/Breast_Mammary_Tissue.parquet', 
                        'gtex_rna/networks/Lung.parquet',
                        'gtex_rna/networks/Stomach.parquet',


                        'cellnet_human_Hg1332/networks/bcell.parquet',
                        'cellnet_human_Hg1332/networks/tcell.parquet',
                        'cellnet_human_Hg1332/networks/skin.parquet',
                        'cellnet_human_Hg1332/networks/neuron.parquet',
                        'cellnet_human_Hg1332/networks/heart.parquet',
                        ],
      'global_models_dir': '../eric/network_collection/networks/',

      "evaluation_data": f"resources/evaluation_datasets/{dataset}_perturbation.h5ad",
      'consensus': f'resources/prior/{dataset}_consensus-num-regulators.json',

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
for dataset in ['adamson']: #'op', 'replogle2', 'nakatake', 'norman', 'adamson'
  print('------ ', dataset, '------')
  par = define_par(dataset)
  os.makedirs(par['scores_dir'], exist_ok=True)
  main_consensus(par)
  for binarize in [True]:
    par['binarize'] = binarize
    for max_n_links in [50000]:
      par['max_n_links'] = max_n_links
      for apply_skeleton in [False]:
        par['apply_skeleton'] = apply_skeleton
        # - determines models to run 
        grn_files_dict = {}
        # - add global models
        if global_models:
          for model in par['global_models']:
            temp_dir = f"{par['scores_dir']}/nets/"
            os.makedirs(temp_dir, exist_ok=True)
            net = pd.read_parquet(f"{par['global_models_dir']}/{model}")
            net.columns = ['source','target','weight']
            net = process_links(net, par)
            if par['binarize']:
                net['weight'] = net['weight'].apply(binarize_weight) 
            model = model.replace('/','_')
            grn_file = f'{temp_dir}/{model}.csv'
            net.to_csv(grn_file)
            grn_files_dict[model] = grn_file
        else:
          # - add actual models
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
  
