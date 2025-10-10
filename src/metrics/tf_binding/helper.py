import math
import pandas as pd
import numpy as np
from sklearn.metrics import average_precision_score
from tqdm import tqdm
from util import read_prediction


# For reproducibility
seed = 42
np.random.seed(seed)

def main(par):

    # Load predictions
    prediction = read_prediction(par)
    assert prediction.shape[0] > 0, 'No links found in the network'
    if not {'source', 'target', 'weight'}.issubset(set(prediction.columns)):
        raise ValueError("prediction must have columns: 'source', 'target', 'weight'")

    # Load TF list and ground truth edges
    tf_all = np.loadtxt(par['tf_all'], dtype=str, delimiter=',', skiprows=1)
    true_graph = pd.read_csv(par['ground_truth'])
    if not {'source', 'target'}.issubset(set(true_graph.columns)):
        raise ValueError("ground_truth must have columns: 'source', 'target'")

    # Keep only TFs that are in the provided TF list
    true_graph = true_graph[true_graph['source'].isin(tf_all)]
    assert true_graph.shape[0] > 0, 'No links found in the ground truth after filtering to TF list'

    # Use the union of all targets that appear in either predictions or ground-truth
    pred_targets = set(prediction['target'].astype(str).unique())
    true_targets = set(true_graph['target'].astype(str).unique())
    evaluation_genes = sorted(pred_targets.union(true_targets))
    if len(evaluation_genes) == 0:
        raise ValueError("Empty evaluation target set (no targets in predictions or ground truth).")

    # Precompute lookup to avoid O(n^2)
    gene_to_idx = {g: i for i, g in enumerate(evaluation_genes)}
    n_targets_total = len(evaluation_genes)

    scores_model = []
    # Iterate over TFs that actually have at least one positive in ground-truth
    for tf in tqdm(true_graph['source'].unique()):
        # Positives for this TF
        true_edges = true_graph[true_graph['source'] == tf]
        y_true = np.zeros(n_targets_total, dtype=int)
        # Mark positives
        for t in true_edges['target'].astype(str):
            y_true[gene_to_idx[t]] = 1

        # Scores over all candidate targets
        y_score = np.zeros(n_targets_total, dtype=float)
        if tf in prediction['source'].unique():
            pred_edges = prediction[prediction['source'] == tf]

            # Fill scores where predicted
            for tgt, w in zip(pred_edges['target'].astype(str), pred_edges['weight'].astype(float)):
                idx = gene_to_idx.get(tgt, None)
                if idx is not None:
                    y_score[idx] = w

            # Ensure numerical safety
            y_score = np.nan_to_num(y_score, nan=0.0, posinf=0.0, neginf=0.0)

        # AP requires at least one positive
        ap = average_precision_score(y_true, y_score) if (y_true.sum()) > 0 else 0.0

        # Random baseline AUPRC is just the prevalence of positives
        prevalence = y_true.mean()

        scores_model.append({
            'source': tf,
            'ap': ap,
            'n_targets_pos': int(y_true.sum()),
            'prevalence': prevalence
        })

    scores_df = pd.DataFrame(scores_model)
    valid = scores_df[scores_df['n_targets_pos'] > 0]
    macro_ap = valid['ap'].mean() if len(valid) > 0 else np.nan
    macro_ap_w = np.average(valid['ap'], weights=valid['n_targets_pos']) if len(valid) > 0 else np.nan
    mean_prev = valid['prevalence'].mean() if len(valid) > 0 else np.nan
    lift = macro_ap / mean_prev if (mean_prev is not None and mean_prev and not np.isnan(mean_prev)) else np.nan

    summary_df = pd.DataFrame([{
        # 'macro_ap': macro_ap,
        #'macro_ap_weighted': macro_ap_w,
        # 'mean_prevalence': mean_prev,
        'lift_over_random': lift
    }])

    return summary_df
