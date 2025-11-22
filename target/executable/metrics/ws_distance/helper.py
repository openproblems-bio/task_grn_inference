import pandas as pd
import anndata as ad
import sys
import numpy as np
from tqdm import tqdm

from util import read_prediction
# For reproducibility
seed = 42
np.random.seed(seed)

def main(par):
    prediction = read_prediction(par)
    
    assert prediction.shape[0]>0, 'No links found in the network'
    
    consensus = pd.read_csv(par['ws_consensus']) #  ['source', 'theta', 'value']
    background_distance = pd.read_csv(par['ws_distance_background']) # ['source', 'target', 'ws_distance']
    background_tfs = background_distance['source'].unique()
    consensus_tfs = consensus['source'].unique()
    consensus_tfs_common = np.intersect1d(consensus_tfs, background_tfs) # only keep tfs that are in the background distance
    print('Number of tfs in the consensus:', len(consensus_tfs), 
          ', Number of tfs in the background distance:', len(background_tfs),
          ', Number of common tfs:', len(consensus_tfs_common))
    prediction_tfs = prediction['source'].unique()
    # - for each theta, and each tf: 
    scores_model = []
    if True: # raw scores
        print('Processing raw (no theta):', flush=True)
        scores_raw = []
        for tf in tqdm(prediction_tfs, desc='Processing tfs (raw)', leave=False):
            # Only process TFs that are in the background distance
            if tf not in background_tfs:
                continue
            
            # Get all edges from the prediction for this TF
            prediction_tf = prediction[prediction['source']==tf]
            
            # Get the background distance for the given tf
            background_distance_tf = background_distance[background_distance['source']==tf]
            
            # Get the ws distance for predicted edges that exist in background
            ws_distance = background_distance_tf[background_distance_tf['target'].isin(prediction_tf['target'])].copy()
            
            if len(ws_distance) > 0:
                ws_distance['present_edges_n'] = len(ws_distance)
                # Use the raw ws_distance values directly (no percentile ranking)
                ws_distance['ws_distance_pc'] = ws_distance['ws_distance']
                ws_distance['tf'] = tf
                ws_distance['theta'] = 'ws_raw'
                scores_raw.append(ws_distance)
        if len(scores_raw)>0:
            scores_raw = pd.concat(scores_raw).reset_index(drop=True)
            scores_model.append(scores_raw)
    for theta in consensus['theta'].unique():
        print('Processing theta:', theta, flush=True)
        consensus_theta = consensus[consensus['theta'] == theta]
        for tf in tqdm(consensus_tfs_common, desc='Processing tfs', leave=False):
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
            ws_distance['tf'] = tf
            
            ws_distance['theta'] = theta
            scores_model.append(ws_distance)
    
    
    scores_model = pd.concat(scores_model).reset_index(drop=True)
    mean_scores = scores_model.groupby('theta')['ws_distance_pc'].mean().to_frame().T.reset_index(drop=True)
    
    # Rename theta columns to precision/recall
    column_mapping = {
        0.25: 'ws_precision',
        0.75: 'ws_recall',
        'ws_raw': 'ws_raw'
    }
    
    # Rename columns if they exist
    new_columns = {}
    for old_col in mean_scores.columns:
        if old_col in column_mapping:
            new_columns[old_col] = column_mapping[old_col]
        else:
            new_columns[old_col] = old_col
    mean_scores = mean_scores.rename(columns=new_columns)
    
    # Calculate F1 score
    if 'ws_precision' in mean_scores.columns and 'ws_recall' in mean_scores.columns:
        precision = mean_scores['ws_precision'].values[0]
        recall = mean_scores['ws_recall'].values[0]
        f1 = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0
        mean_scores['ws_f1'] = f1
    
    return scores_model, mean_scores