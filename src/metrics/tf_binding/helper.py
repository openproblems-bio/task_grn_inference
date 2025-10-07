
import pandas as pd
import anndata as ad
import sys
import numpy as np
from sklearn.metrics import average_precision_score
from tqdm import tqdm
from util import read_prediction

# For reproducibility
seed = 42
np.random.seed(seed)

def main(par):
    prediction = read_prediction(par)
    tf_all = np.loadtxt(par['tf_all'], dtype=str, delimiter=',', skiprows=1)
    true_graph = pd.read_csv(par['ground_truth'])
    true_graph = true_graph[true_graph['source'].isin(tf_all)]
    assert prediction.shape[0] > 0, 'No links found in the network'
    assert true_graph.shape[0] > 0, 'No links found in the ground truth'
    
    scores_model = []
    for tf in tqdm(true_graph['source'].unique()):
        true_edges = true_graph[true_graph['source'] == tf]
        if tf in prediction['source'].unique():
            pred_edges = prediction[prediction['source'] == tf]
            true_labels = true_edges['target'].isin(pred_edges['target']).astype(int)
            pred_scores = pred_edges.set_index('target').reindex(true_edges['target'])['weight'].fillna(0)
            ap = average_precision_score(true_labels, pred_scores)
        else:
            ap = float('nan')
        
        scores_model.append({'source': tf, 'ap': ap})
    
    scores_df = pd.DataFrame(scores_model)
    print(scores_df)

    # Precision: mean over available TFs (ignoring NaNs)
    precision = scores_df['ap'].mean(skipna=True)

    # Recall: mean over all TFs, punishing NaNs as 0
    recall = scores_df['ap'].fillna(0).mean()

    # One-row summary DataFrame
    summary_df = pd.DataFrame([{'precision': precision, 'recall': recall}])

    return summary_df