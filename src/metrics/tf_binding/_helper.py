
import math
import pandas as pd
import numpy as np
from tqdm import tqdm
import sys
import os

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
    
    # Get TFs that are both in ground truth and TF_all list
    tfs_in_gt = set(true_graph['source'].unique())
    tfs_to_evaluate = tfs_in_gt & set(tf_all)
    
    # Iterate over TFs that should be evaluated
    for tf in tqdm(tfs_to_evaluate):
        # Positives for this TF
        true_edges = true_graph[true_graph['source'] == tf]
        gt_targets = set(true_edges['target'].astype(str))
        k = len(gt_targets)  # Number of true targets for this TF
        
        if tf in prediction['source'].unique():
            # TF is predicted - calculate actual metrics
            pred_edges = prediction[prediction['source'] == tf]
            
            # Get top k predicted targets for this TF
            top_k_pred_targets = set(pred_edges.sort_values(by='weight', key=abs, ascending=False).head(k)['target'].astype(str))
            
            # Calculate precision and recall for predicted edges
            true_positives = len(gt_targets & top_k_pred_targets)
            predicted_positives = len(top_k_pred_targets)
            actual_positives = len(gt_targets)
            
            # Precision and recall for predictions
            precision_pred = true_positives / predicted_positives if predicted_positives > 0 else 0.0
            recall_pred = true_positives / actual_positives if actual_positives > 0 else 0.0
            f1_pred = 2 * (precision_pred * recall_pred) / (precision_pred + recall_pred) if (precision_pred + recall_pred) > 0 else 0.0
            
            # Calculate baseline (random selection of k targets)
            # Expected precision for random selection
            precision_random = actual_positives / n_targets_total if n_targets_total > 0 else 0.0
            recall_random = min(k / actual_positives, 1.0) if actual_positives > 0 else 0.0  # Can't exceed 1.0
            f1_random = 2 * (precision_random * recall_random) / (precision_random + recall_random) if (precision_random + recall_random) > 0 else 0.0
            
            # Calculate lifts
            precision_lift = precision_pred / (precision_random + 1e-10)
            f1_lift = f1_pred / (f1_random + 1e-10)
            
        else:
            # TF is missing from predictions - assign random baseline scores (lift = 1.0)
            precision_lift = 1.0  # Random performance
            f1_lift = 1.0  # Random performance
        
        scores_model.append({
            'source': tf,
            'precision_lift': precision_lift,
            'f1_lift': f1_lift,
            'n_targets_pos': k,
            'is_predicted': tf in prediction['source'].unique()
        })

    scores_df = pd.DataFrame(scores_model)
    
    # Calculate the 4 scores as requested:
    
    # 1. Precision lift - available TFs only (TFs that are predicted)
    available_tfs = scores_df[scores_df['is_predicted'] == True]
    precision_lift_available = available_tfs['precision_lift'].mean() if len(available_tfs) > 0 else 1.0
    
    # 2. F1 lift - available TFs only (TFs that are predicted)
    f1_lift_available = available_tfs['f1_lift'].mean() if len(available_tfs) > 0 else 1.0
    
    # 3. Precision lift - all TFs (including missing ones with score 1.0)
    precision_lift_all = scores_df['precision_lift'].mean()
    
    # 4. F1 lift - all TFs (including missing ones with score 1.0)
    f1_lift_all = scores_df['f1_lift'].mean()
    
    # Calculate TF coverage statistics
    n_tfs_in_gt = len(tfs_in_gt)
    n_tfs_to_evaluate = len(tfs_to_evaluate)
    tfs_in_grn = set(prediction['source'].unique())
    tfs_evaluated = tfs_to_evaluate & tfs_in_grn
    n_tfs_in_grn = len(tfs_in_grn)
    n_tfs_evaluated = len(tfs_evaluated)
    
    # Print debug info
    # print(f"\n=== TF COVERAGE DEBUG ===")
    # print(f"TFs in ground truth: {n_tfs_in_gt}")
    # print(f"TFs to evaluate (in GT & TF_all): {n_tfs_to_evaluate}")
    # print(f"TFs in prediction: {n_tfs_in_grn}")
    # print(f"TFs evaluated (predicted & should evaluate): {n_tfs_evaluated}")
    # print(f"TF coverage: {n_tfs_evaluated/n_tfs_to_evaluate*100:.1f}%")
    
    # print(f"\n=== METRIC DEBUG ===")
    # print(f"Precision lift (available TFs): {precision_lift_available:.4f}")
    # print(f"F1 lift (available TFs): {f1_lift_available:.4f}")
    # print(f"Precision lift (all TFs): {precision_lift_all:.4f}")
    # print(f"F1 lift (all TFs): {f1_lift_all:.4f}")

    summary_df = pd.DataFrame([{
        'tf_binding_precision_available': precision_lift_available,
        'tf_binding_f1_available': f1_lift_available,
        'tf_binding_precision_all': precision_lift_all,
        'tf_binding_f1_all': f1_lift_all,
        
        'n_tfs_in_gt': n_tfs_in_gt,
        'n_tfs_to_evaluate': n_tfs_to_evaluate,
        'n_tfs_in_grn': n_tfs_in_grn,
        'n_tfs_evaluated': n_tfs_evaluated,
        'tf_coverage_pct': n_tfs_evaluated/n_tfs_to_evaluate*100 if n_tfs_to_evaluate > 0 else 0
    }])

    return summary_df