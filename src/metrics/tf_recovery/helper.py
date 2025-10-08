import pandas as pd
import numpy as np
import decoupler as dc
import sys
import os
import re
import anndata as ad
from scipy.stats import ttest_rel


    
    
from util import read_prediction

def calculate_ulm_score(tf_mat, tf_grn):
    """Calculate ULM score for a TF and its network."""
    act, pval = dc.run_ulm(
        mat=tf_mat,
        net=tf_grn,
        weight='weight',
        min_n=3,
    )
    return act.values[0, 0], pval.values[0, 0]

def create_random_network(tf, n_targets, all_genes):
    """Create a random network for a TF."""
    random_targets = np.random.choice(all_genes, size=n_targets, replace=False)
    return pd.DataFrame({
        'source': [tf] * n_targets,
        'target': random_targets,
        'weight': [1.0] * n_targets
    })

def calculate_precision_recall(acts, pvals, total_tfs):
    """Calculate precision and recall scores."""
    padj = dc.p_adjust_fdr(pvals)
    tp = np.sum((acts < 0) & (padj < 0.05))
    prc = tp / acts.size if tp > 0 else 0.0
    rcl = tp / total_tfs if tp > 0 else 0.0
    return prc, rcl

def main(par):
    df_de = ad.read_h5ad(par['evaluation_data_de'])
    df_de = df_de.to_df()  # convert X to DataFrame
    net = read_prediction(par)

    # filter TFs with less than 3 targets
    tf_counts = net['source'].value_counts()
    tfs_to_keep = tf_counts[tf_counts >= 3].index
    net = net[net['source'].isin(tfs_to_keep)]
    n_tfs = net['source'].nunique()

    all_genes = list(df_de.columns)
    acts_in_net, pvals_in_net, acts_random_in_net, pvals_random_in_net = [], [], [], []
    
    # Process TFs that are in the network -> precision 
    for tf in df_de.index:
        if tf not in net['source'].unique():
            continue

        tf_mat = df_de[df_de.index == tf]
        tf_grn = net[net['source'] == tf]
        shared_genes = set(tf_grn['target']) & set(df_de.columns)
        if len(shared_genes) < 3:
            continue 
        
        # Calculate actual network scores
        act, pval = calculate_ulm_score(tf_mat, tf_grn)
        acts_in_net.append(act)
        pvals_in_net.append(pval)
        
        # Calculate random network scores
        tf_grn_random = create_random_network(tf, len(shared_genes), all_genes)
        act_random, pval_random = calculate_ulm_score(tf_mat, tf_grn_random)
        acts_random_in_net.append(act_random)
        pvals_random_in_net.append(pval_random)
    
    median_random_act = np.median(acts_random_in_net)
    
    # Process all TFs in df_de ->recall
    acts_all, pvals_all, acts_random_all, pvals_random_all = [], [], [], []
    
    for tf in df_de.index:
        tf_mat = df_de[df_de.index == tf]
        
        if tf in net['source'].unique():
            tf_grn = net[net['source'] == tf]
            shared_genes = set(tf_grn['target']) & set(df_de.columns)
            if len(shared_genes) < 3:
                continue
                
            act, pval = calculate_ulm_score(tf_mat, tf_grn)
            acts_all.append(act)
            pvals_all.append(pval)
            
            tf_grn_random = create_random_network(tf, len(shared_genes), all_genes)
            act_random, pval_random = calculate_ulm_score(tf_mat, tf_grn_random)
            acts_random_all.append(act_random)
            pvals_random_all.append(pval_random)
        else:
            # TF not in network - assign median random act
            acts_all.append(median_random_act)
            pvals_all.append(1.0)
            acts_random_all.append(median_random_act)
            pvals_random_all.append(1.0)

    # Analyze both scenarios
    datasets = {
        'in_net': (acts_in_net, pvals_in_net, acts_random_in_net, pvals_random_in_net),
        'all': (acts_all, pvals_all, acts_random_all, pvals_random_all)
    }
    
    results = {}
    for name, (acts, pvals, acts_random, pvals_random) in datasets.items():
        acts, pvals = np.array(acts), np.array(pvals)
        acts_random, pvals_random = np.array(acts_random), np.array(pvals_random)
        
        print(f"\n=== Analysis: {name} ===")
        
        # Calculate precision/recall (side scores)
        prc, rcl = calculate_precision_recall(acts, pvals, df_de.shape[0])
        prc_random, rcl_random = calculate_precision_recall(acts_random, pvals_random, df_de.shape[0])
        
        print(f"Side scores - Actual: precision={prc:.3f}, recall={rcl:.3f}")
        print(f"Side scores - Random: precision={prc_random:.3f}, recall={rcl_random:.3f}")
        
        # Paired t-test on absolute values
        t_stat, p_value = ttest_rel(np.abs(acts), np.abs(acts_random))
        results[name] = {'t_stat': t_stat, 'p_value': p_value, 'count': len(acts)}
    print({
        'tfs_in_net_count': [results['in_net']['count']],
        'tfs_all_count': [results['all']['count']],
        't_stat_in_net': [results['in_net']['t_stat']],
        'p_value_in_net': [results['in_net']['p_value']],
        't_stat_all': [results['all']['t_stat']],
        'p_value_all': [results['all']['p_value']],
        'median_random_act': [median_random_act]
    })
    df = pd.DataFrame({
        't_rec_precision': [results['in_net']['t_stat']],
        't_rec_recall': [results['all']['t_stat']]
    })
    print(f"\n=== Final Results ===")
    print(df)
    return df

