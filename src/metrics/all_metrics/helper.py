import os
import pandas as pd
from regression_2.consensus.script import main as main_consensus_reg2
from wasserstein.consensus.script import main as main_consensus_ws
from wasserstein.background_distance.script import main as main_ws_background_distance
from all_metrics.script import main as main_scores
from all_metrics.script import par as main_par



def run_scores_all(datasets, models):
    scores_dir = 'resources/scores/'
    save_file_name = f"{scores_dir}/default_scores.csv" 

    scores_store = []
    for dataset in datasets:
        for model in models:
            par = main_par.copy()
            # - adjust par
            par['dataset_id'] = dataset
            par['prediction'] = f'resources/grn_models/{dataset}/{model}.csv'
            if not os.path.exists(par['prediction']):
                print('Skipping ', par['prediction'])
                continue
            # - run
            scores_model = main_scores(par)
            scores_model['model'] = model
            scores_model['dataset'] = dataset

            scores_store.append(scores_model)
            scores_all = pd.concat(scores_store)
            scores_all.to_csv(save_file_name)

def run_consensus(datasets):
    models = ['positive_control', 'pearson_corr', 'portia', 'ppcor', 'grnboost2', 'scenic', 'granie', 'scglue', 'celloracle', 'figr', 'scenicplus']

    for dataset in datasets:
        par = {
            'models': models,
            'evaluation_data': f'resources/grn_benchmark/evaluation_data//{dataset}.h5ad',
            'evaluation_data_sc': f'resources/datasets_raw/{dataset}_sc_counts.h5ad',
            'models_dir': f'resources/grn_models/{dataset}/',
            'regulators_consensus': f'resources/grn_benchmark/prior/regulators_consensus_{dataset}.json',
            'ws_consensus': f'resources/grn_benchmark/prior/ws_consensus_{dataset}.csv',
            'tf_all': 'resources/grn_benchmark/prior/tf_all.csv',

        }
        # - reg2 consensus 
        print(f'--determining consensus for reg2--{dataset}')
        main_consensus_reg2(par)

        # - ws consensus
        print(f'--determining consensus for ws--{dataset}')
        if dataset in ['norman', 'adamson']:
            main_consensus_ws(par)
def run_ws_distance_background(datasets):
  for dataset in datasets:
    par = {
          'evaluation_data_sc': f'resources/datasets_raw/{dataset}_sc_counts.h5ad',
          'background_distance': f'resources/grn_benchmark/prior/ws_distance_background_{dataset}.csv',
          'tf_all': 'resources/grn_benchmark/prior/tf_all.csv',
          'layer': 'X_norm'
    }
    print(f'--run ws distance background --{dataset}')
    if dataset in ['norman', 'adamson']:
      main_ws_background_distance(par)
