
import math
import pandas as pd
import numpy as np
import anndata as ad
from tqdm import tqdm
import sys
import os

from util import read_prediction


# For reproducibility
seed = 42
np.random.seed(seed)

def main_sub(par):

    evaluation_data = ad.read_h5ad(par['evaluation_data'], backed='r')
    genes = evaluation_data.var_names.tolist()
    n_targets = len(genes)
    prediction = read_prediction(par)
    assert prediction.shape[0] > 0, 'No links found in the network'
    tf_all = np.loadtxt(par['tf_all'], dtype=str, delimiter=',', skiprows=1)
    true_graph = pd.read_csv(par['ground_truth'])
    if not {'source', 'target'}.issubset(set(true_graph.columns)):
        raise ValueError("ground_truth must have columns: 'source', 'target'")

    true_graph = true_graph[true_graph['source'].isin(tf_all)]
    assert true_graph.shape[0] > 0, 'No links found in the ground truth after filtering to TF list'

    scores_model = []
    tfs_in_gt = set(true_graph['source'].unique())
    tfs_to_evaluate = tfs_in_gt & set(tf_all)  # TFs in both GT and TF_all list
    
    for tf in tqdm(tfs_to_evaluate):
        true_edges = true_graph[true_graph['source'] == tf]
        gt_targets = set(true_edges['target'].astype(str))
        gt_targets = set(gt_targets) & set(genes)
        k = len(gt_targets)  # Number of true targets for this TF
        if k == 0:
            continue
        
        if tf in prediction['source'].unique():
            pred_edges = prediction[prediction['source'] == tf]
            
            top_k_t = set(pred_edges.sort_values(by='weight', key=abs, ascending=False).head(k)['target'].astype(str))
            
            tp_n = len(gt_targets & top_k_t)  
            top_k_t_n = len(top_k_t)  # Actual number returned (might be less than k)          
            
            # Calculate precision and recall scores
            precision_score = tp_n / top_k_t_n if top_k_t_n > 0 else 0.0             
            precision_baseline = k / n_targets
            
        else:
            precision_score = 0.0
        
        scores_model.append({
            'source': tf,
            'precision_score': precision_score,
            'n_targets_pos': k,
            'is_predicted': tf in prediction['source'].unique()
        })

    scores_df = pd.DataFrame(scores_model)
    available_tfs = scores_df[scores_df['is_predicted'] == True]
    precision_available = available_tfs['precision_score'].mean() if len(available_tfs) > 0 else 0.0
    
    precision_all = scores_df['precision_score'].mean()
    
    
    # Calculate TF coverage statistics
    n_tfs_in_gt = len(tfs_in_gt)
    n_tfs_to_evaluate = len(tfs_to_evaluate)
    tfs_in_grn = set(prediction['source'].unique())
    tfs_evaluated = tfs_to_evaluate & tfs_in_grn
    n_tfs_in_grn = len(tfs_in_grn)
    n_tfs_evaluated = len(tfs_evaluated)
    
    summary_df = pd.DataFrame([{
        'tfb_grn': precision_available,
        'tfb_all': precision_all,

        'n_tfs_in_gt': n_tfs_in_gt,
        'n_tfs_to_evaluate': n_tfs_to_evaluate,
        'n_tfs_in_grn': n_tfs_in_grn,
        'n_tfs_evaluated': n_tfs_evaluated,
        'tf_coverage_pct': n_tfs_evaluated/n_tfs_to_evaluate*100 if n_tfs_to_evaluate > 0 else 0
    }])

    return summary_df

def main(par):
    gts = ['unibind', 'chipatlas', 'remap']

    rr_store = []
    for gt in gts:
        par['ground_truth'] = par[f'ground_truth_{gt}']
        df = main_sub(par)
        df['gt'] = gt
        rr_store.append(df)
    result_df = pd.concat(rr_store, ignore_index=True)
    
    # Weighted average across all ground truth sources
    weights = result_df['n_tfs_to_evaluate'].values
    total_weight = weights.sum()
    
    tfb_grn_weighted = (result_df['tfb_grn'] * weights).sum() / total_weight
    tfb_all_weighted = (result_df['tfb_all'] * weights).sum() / total_weight
    
    
    result_df = pd.DataFrame([{
        'tfb_grn': tfb_grn_weighted,
        'tfb_all': tfb_all_weighted
    }])
    
    return result_df
