
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
    test_data = ad.read_h5ad(par['evaluation_data'], backed='r')
    evaluation_genes = test_data.var_names.tolist()
    n_targets_total = len(evaluation_genes)

    tf_all = np.loadtxt(par['tf_all'], dtype=str, delimiter=',', skiprows=1)
    true_graph = pd.read_csv(par['ground_truth'])
    true_graph = true_graph[(true_graph['source'].isin(tf_all)) & (true_graph['target'].isin(evaluation_genes))]
    assert prediction.shape[0] > 0, 'No links found in the network'
    assert true_graph.shape[0] > 0, 'No links found in the ground truth'
    
    scores_model = []
    for tf in tqdm(true_graph['source'].unique()):
        true_edges = true_graph[true_graph['source'] == tf]
        if tf in prediction['source'].unique():
            pred_edges = prediction[prediction['source'] == tf]
            true_labels = true_edges['target'].isin(pred_edges['target']).astype(int)
            pred_scores = pred_edges.set_index('target').reindex(true_edges['target'])['weight'].fillna(0)
            if true_labels.sum() == 0:  # no positives
                ap = 0.0
            else:
                ap = average_precision_score(true_labels, pred_scores)
        else:
            ap = float('nan')
        n_targets = len(true_edges)

        # ----- Analytic random baseline -----
        # Extend true edges to all evaluation genes
        true_labels_random = np.zeros(n_targets_total)
        idx = [evaluation_genes.index(t) for t in true_edges['target']]
        true_labels_random[idx] = 1
        ap_random = true_labels_random.sum() / len(true_labels_random)

        scores_model.append({'source': tf, 'ap': ap, 'n_targets': n_targets, 'ap_random': ap_random})
    
    scores_df = pd.DataFrame(scores_model)
    print('Number of TFs in GRN:', len(scores_df[scores_df['ap'].notna()]))

    # Compute weighted mean (ignoring NaNs)
    valid = scores_df.dropna(subset=['ap'])
    weighted_precision = np.average(valid['ap'], weights=valid['n_targets'])

    # Compute unweighted means (for reference)
    precision = scores_df['ap'].mean(skipna=True)
    precision_random = scores_df['ap_random'].mean(skipna=True)
    recall = scores_df['ap'].fillna(0).mean()

    summary_df = pd.DataFrame([{
        'precision': precision,
        'precision_random': precision_random,
        'recall': recall,
        'weighted_precision': weighted_precision
    }])

    return summary_df