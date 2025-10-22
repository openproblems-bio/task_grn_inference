
import math
import pandas as pd
import numpy as np
from sklearn.metrics import average_precision_score
from tqdm import tqdm
from util import read_prediction
from scipy.stats import wilcoxon

# Wilcoxon signed-rank test: AP vs prevalence
def wilcoxon_logp(x, y):
    try:
        stat, p = wilcoxon(x, y, zero_method='wilcox', alternative='greater')
        return -np.log10(p) if p > 0 else np.inf
    except Exception:
        return np.nan


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
        assert y_true.sum() > 0, f'TF {tf} has no positive targets in ground truth'

        # Scores over all candidate targets
        y_score = np.zeros(n_targets_total, dtype=float)
        if tf in prediction['source'].unique():
            pred_edges = prediction[prediction['source'] == tf]

            # Fill scores where predicted
            for tgt, w in zip(pred_edges['target'].astype(str), pred_edges['weight'].abs().astype(float)):
                idx = gene_to_idx.get(tgt, None)
                if idx is not None:
                    y_score[idx] = w

            # Ensure numerical safety
            y_score = np.nan_to_num(y_score, nan=0.0, posinf=0.0, neginf=0.0)
            # AP requires at least one positive
            ap = average_precision_score(y_true, y_score)
        else:
            ap = np.nan

        # Random baseline AUPRC is just the prevalence of positives
        prevalence = y_true.mean()

        # Calculate recall@k for this TF
        
        if tf in prediction['source'].unique() and not np.isnan(ap):
            gt_targets = set(true_edges['target'].astype(str))
            k = len(gt_targets)  # Top k predictions to consider
            # Get top k predicted targets for this TF
            genes_in_grn = set(pred_edges.sort_values(by='weight', key=abs, ascending=False).head(k)['target'].astype(str).unique() )
            
            # Calculate recall@k
            true_positives_in_top_k = len(gt_targets & genes_in_grn)
            recall_at_k = true_positives_in_top_k / len(gt_targets)
            
            # Calculate random recall baseline
            random_recall = len(genes_in_grn) / n_targets_total
            lift_recall = recall_at_k / (random_recall + 1e-10)
        else:
            recall_at_k = np.nan
            lift_recall = np.nan
        
        scores_model.append({
            'source': tf,
            'ap': ap,
            'n_targets_pos': int(y_true.sum()),
            'prevalence': prevalence,
            'recall_at_k': recall_at_k,
            'lift_recall': lift_recall
        })

        from scipy.stats import wilcoxon

    scores_df = pd.DataFrame(scores_model)
    ap_values_pred = scores_df['ap'].dropna().values
    preval_all = scores_df['prevalence'].values
    preval_pred = scores_df[~scores_df['ap'].isna()]['prevalence'].values
    ap_pred_mean = ap_values_pred.mean() if len(ap_values_pred) > 0 else 0.0
    preval_pred_mean = preval_pred.mean() if len(preval_pred) > 0 else 0.0
    lift_pred = ap_pred_mean / (preval_pred_mean + 1e-10)

    # Calculate mean recall metrics
    recall_values_pred = scores_df['recall_at_k'].dropna().values
    lift_recall_values_pred = scores_df['lift_recall'].dropna().values
    lift_recall_values_all = scores_df['lift_recall'].fillna(1).values

    mean_recall_at_k = recall_values_pred.mean() if len(recall_values_pred) > 0 else 0.0
    mean_lift_recall = lift_recall_values_pred.mean() if len(lift_recall_values_pred) > 0 else 1
    mean_lift_recall_all = lift_recall_values_all.mean() if len(lift_recall_values_all) > 0 else 1

    # Calculate TF coverage
    n_tfs_in_gt = true_graph['source'].nunique()
    n_tfs_in_grn = set(prediction['source'].unique()) & set(true_graph['source'].unique())
    n_tfs_in_grn = len(n_tfs_in_grn)

    recall = mean_lift_recall_all * mean_lift_recall

    summary_df = pd.DataFrame([{
        'tf_binding_precision': lift_pred,
        'tf_binding_recall': recall
    }])

    return summary_df