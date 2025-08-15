import pandas as pd
import anndata as ad
import sys
import numpy as np

# For reproducibility
seed = 42
np.random.seed(seed)

def main(par):
    prediction = ad.read_h5ad(par['prediction'])
    prediction = pd.DataFrame(prediction.uns['prediction'])
    prediction['weight'] = prediction['weight'].astype(float)
    assert prediction.shape[0]>0, 'No links found in the network'
    assert all(column in prediction.columns for column in ['source', 'target', 'weight']), 'Columns in the network should be source, target, weight'
    
    consensus = pd.read_csv(par['ws_consensus']) #  ['source', 'theta', 'value']
    background_distance = pd.read_csv(par['ws_distance_background']) # ['source', 'target', 'ws_distance']
    background_tfs = background_distance['source'].unique()
    consensus_tfs = consensus['source'].unique()
    consensus_tfs_common = np.intersect1d(consensus_tfs, background_tfs) # only keep tfs that are in the background distance
    print('Number of tfs in the consensus:', len(consensus_tfs), 
          ', Number of tfs in the background distance:', len(background_tfs),
          ', Number of common tfs:', len(consensus_tfs_common))
    prediction_tfs = prediction['source'].unique()
    # evaluation_data = ad.read_h5ad(par['evaluation_data_sc']) 
    # evaluation_data.X = evaluation_data.layers[par['layer']]

    # - for each theta, and each tf: 
    scores_model = []
    for theta in consensus['theta'].unique():
        consensus_theta = consensus[consensus['theta'] == theta]
        for tf in consensus_tfs_common:
            # - get the background distance for the given tf
            background_distance_tf = background_distance[background_distance['source']==tf]
            n_edges = consensus_theta[consensus_theta['source'] == tf]['value'].values[0]
            if tf in prediction_tfs:             
                # - subset the prediction to the given tf: choose the top edges based on n_edges
                prediction_tf = prediction[prediction['source']==tf]
                prediction_tf = prediction_tf.nlargest(n_edges, 'weight')
                # - get the ws distance
                ws_distance = background_distance_tf[background_distance_tf['target'].isin(prediction_tf['target'])].copy()
            else: # fill the missing TF with random scores
                ws_distance = background_distance_tf.sample(n_edges, replace=True, random_state=seed).copy()
            # - fill the missing links with random scores
            n_missing = n_edges - len(ws_distance)
            ws_distance_missing = background_distance_tf.sample(n_missing, replace=True, random_state=seed)
            ws_distance = pd.concat([ws_distance, ws_distance_missing])
            ws_distance['present_edges_n'] = len(ws_distance)
            # - normalize to the background dist -> percentile rank 
            background_distance_random = np.random.choice(background_distance_tf['ws_distance'].values, 1000)
            ws_distance_pc = ws_distance['ws_distance'].apply(lambda val: (val>background_distance_random).sum())/len(background_distance_random)
            ws_distance['ws_distance_pc'] = ws_distance_pc
            
            ws_distance['theta'] = theta
            scores_model.append(ws_distance)
    scores_model = pd.concat(scores_model).reset_index(drop=True)
    mean_scores = scores_model.groupby('theta')['ws_distance_pc'].mean().to_frame().T.reset_index(drop=True)
    return scores_model, mean_scores