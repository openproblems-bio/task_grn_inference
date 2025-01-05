import pandas as pd
import anndata as ad
import sys
import numpy as np
import os

from script import main
from consensus.script import main as main_consensus 
from background_scores.script import main as main_background_score

# acutal evaluation -> all 


par = {
    'datasets': ['adamson', 'norman'],
    'models': ['pearson_corr', 'grnboost2', 'portia', 'ppcor', 'scenic'],
    'evaluation_data_sc': 'resources/datasets_raw/adamson_sc_counts.h5ad',
    'mean_scores_all': 'resources/scores/ws_distance_mean.csv',
    'scores_all': 'resources/scores/ws_distance.csv',
    'consensus': 'resources/prior/consensus_ws_distance_adamson.csv',
    'tf_all': 'resources/prior/tf_all.csv',
    'background_distance':'resources/prior/ws_distance_background_adamson.csv',
    'layer': 'X_norm',
}

def main(par):
    mean_scores_store = []
    scores_store = []
    for dataset in par['datasets']:
        par['grns_dir'] = f'resources/grn_models/{dataset}'

        par['evaluation_data_sc'] = f'resources/datasets_raw/{dataset}_sc_counts.h5ad'
        par['consensus'] = f'resources/prior/consensus_ws_distance_{dataset}.csv'
        par['background_distance'] = f'resources/prior/ws_distance_background_{dataset}.csv'

        if True:
            main_consensus(par)
        if False:
            main_background_score(par)

        for model in par['models']:
            par['prediction'] = f'resources/grn_models/{dataset}/{model}.csv'
            if not os.path.exists(par['prediction']):
                continue
            scores_model = main(par)
            mean_score = scores_model.groupby('theta')['ws_distance_pc'].mean().to_frame().T.reset_index(drop=True)
            mean_score['dataset'] = dataset
            mean_score['model'] = model
            mean_scores_store.append(mean_score)

            # - also store raw scores
            scores_model['dataset'] = dataset
            scores_model['model'] = model
            scores_store.append(scores_model)
    mean_scores_all = pd.concat(mean_scores_store)
    scores_all = pd.concat(scores_store)
    
    return mean_scores_all, scores_all

if __name__ == '__main__':
    mean_scores_all, scores_all = main(par)
    mean_scores_all.to_csv(par['mean_scores_all'])
    scores_all.to_csv(par['scores_all'])