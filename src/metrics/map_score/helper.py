import pandas as pd
import anndata as ad
import sys
import numpy as np
from sklearn.metrics import average_precision_score
from fetch_ground_truth import build_celltype_grn
from util import process_links
# For reproducibility
seed = 42
np.random.seed(seed)

def main(par):
    prediction = ad.read_h5ad(par['prediction'])
    prediction = pd.DataFrame(prediction.uns['prediction'])
    prediction = process_links(prediction, par)
    cell_type_lookup = {
        'replogle': 'PBMC',
        'adamson': 'PBMC',
        'norman': 'PBMC',
        'xaira_HCT116': 'PBMC',
        'xaira_HEK293T': 'PBMC'
    }
    ground_truth_file = build_celltype_grn(cell_type_lookup[par.get('dataset_id')])
    true_graph = pd.read_csv(ground_truth_file)
    
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