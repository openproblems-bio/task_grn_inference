import pandas as pd
import anndata as ad
import sys
import numpy as np
from sklearn.metrics import average_precision_score

from util import process_links
# For reproducibility
seed = 42
np.random.seed(seed)

def main(par):
    prediction = ad.read_h5ad(par['prediction'])
    prediction = pd.DataFrame(prediction.uns['prediction'])
    prediction = process_links(prediction, par)
    
    ground_truth_files = {
        'replogle': 'tbd/replogle_grn_gt.csv',
        'adamson': 'tbd/adamson_grn_gt.csv',
        'norman': 'tbd/norman_grn_gt.csv',
        'xaira_HCT116': 'tbd/xaira_HCT116_grn_gt.csv',
        'xaira_HEK293T': 'tbd/xaira_HEK293T_grn_gt.csv'
    }
    dataset_id = par.get('dataset_id')
    if dataset_id not in ground_truth_files:
        raise ValueError(f"No ground truth file for dataset {dataset_id}")
    true_graph = pd.read_csv(ground_truth_files[dataset_id])
    
    assert prediction.shape[0] > 0, 'No links found in the network'
    
    scores_model = []
    for tf in true_graph['source'].unique():
        true_edges = true_graph[true_graph['source'] == tf]
        if tf in prediction['source'].unique():
            pred_edges = prediction[prediction['source'] == tf]
            true_labels = true_edges['target'].isin(pred_edges['target']).astype(int)
            pred_scores = pred_edges.set_index('target').reindex(true_edges['target'])['weight'].fillna(0)

            ap = average_precision_score(true_labels, pred_scores)
        else:
            ap = float('nan')
        scores_model.append({'source': tf, 'ap': ap})
    
    scores_model = pd.DataFrame(scores_model)
    mean_scores = scores_model['ap'].mean()
    return scores_model, mean_scores