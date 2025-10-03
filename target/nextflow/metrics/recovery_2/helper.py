import scipy
import numpy as np
import pandas as pd
import anndata as ad
import warnings
import scipy
warnings.filterwarnings("ignore")
from util import format_save_score, read_prediction

from scipy.stats import spearmanr
np.random.seed(0)
from util import manage_layer
def actual_calculation(adata, target_genes, group):
    adata = adata[adata.obs['is_control'] == False] # only use non-control cells
    target_genes = [g for g in target_genes if g in adata.var_names]
    target_adata = adata[:, target_genes]
    target_expr = pd.DataFrame(
        target_adata.X.toarray() if scipy.sparse.issparse(target_adata.X) else target_adata.X,
        index=target_adata.obs.index,
        columns=target_adata.var_names
    )
    for g in group:
        target_expr[g] = target_adata.obs[g].values
    group_consistencies = {}
    for g_vals, df_group in target_expr.groupby(group):
        expr_data = df_group.drop(columns=group)
        if len(expr_data) < 2:
            continue  # skip if only one sample

        rows = expr_data.values
        corr_matrix = np.corrcoef(rows)
        correlation = corr_matrix[np.triu_indices_from(corr_matrix, k=1)] 
        correlation = correlation[~np.isnan(correlation)]
        correlation = np.mean(correlation)
        assert isinstance(correlation, float)
        group_consistencies[g_vals] = correlation
    return group_consistencies

def calculate_target_gene_consistency(adata, net, group, min_targets_t=10, max_targets_t=20):
    """
    Calculate consistency score for each TF based on target gene expression across groups.
    Returns a dictionary with TF as key and consistency score as value.
    """
    from scipy.stats import ttest_rel
    
    adata = adata[adata.obs['is_control'] == False]
    all_genes = adata.var_names.tolist()
    tf_size = net.groupby('source').size()
    tfs = tf_size[(tf_size > min_targets_t)].index
    net = net[net['source'].isin(tfs)]
    
    all_values = []
    all_random_values = []
    
    for tf in tfs:
        target_genes = net[net['source'] == tf]
        target_genes = target_genes.sort_values(by='weight', ascending=False, key=abs).head(max_targets_t)
        target_genes = target_genes['target'].tolist()
            
        actual_consistencies = actual_calculation(adata, target_genes, group)
        
        n_targets = len(target_genes)
        random_targets = np.random.choice(all_genes, size=n_targets, replace=False)
        random_consistencies = actual_calculation(adata, random_targets, group)

        # Get consistency values for groups that exist in both actual and random
        common_groups = actual_consistencies.keys()
        if len(common_groups) < 2:
            continue  # Need at least 2 groups for paired t-test
            
        actual_values = [actual_consistencies[g] for g in common_groups]
        random_values = [random_consistencies[g] for g in common_groups]
        assert len(actual_values) == len(random_values)
        all_values.append(np.mean(actual_values))
        all_random_values.append(np.mean(random_values))

    return all_values, all_random_values

def evaluate_setting(evaluation_adata, net, group, setting_name, n_tfs=None, max_targets_per_tf=None):
    """
    Evaluate a specific setting (recall, balanced, or precision)
    """
    from scipy.stats import wilcoxon
    
    # Filter TFs if specified
    if n_tfs is not None:
        tf_counts = net.groupby('source').size().sort_values(ascending=False)
        top_tfs = tf_counts.head(n_tfs).index
        filtered_net = net[net['source'].isin(top_tfs)]
    else:
        filtered_net = net.copy()
    
    # Set max_targets parameter
    if max_targets_per_tf is None:
        max_targets_per_tf = 50  # default value
    
    actual_values, random_values = calculate_target_gene_consistency(
        evaluation_adata, filtered_net, group, max_targets_t=max_targets_per_tf
    )
    
    # Convert to numpy arrays for arithmetic operations
    actual_values = np.array(actual_values)
    random_values = np.array(random_values)
    
    # Statistical test
    res = wilcoxon(actual_values - random_values, zero_method='wilcox', alternative='greater')
    p_value = res.pvalue
    if p_value == 0:
        p_value = 1e-300
    final_score = -np.log10(p_value)
    
    result = {
        'setting': setting_name,
        'value_min': np.min(actual_values),
        'value_random_min': np.min(random_values),
        'value_median': np.median(actual_values),
        'value_random_median': np.median(random_values),
        'value_max': np.max(actual_values),
        'value_random_max': np.max(random_values),
        'value_mean': np.mean(actual_values),
        'value_random_mean': np.mean(random_values),
        'value_std': np.std(actual_values),
        'value_random_std': np.std(random_values),
        'n_values': len(actual_values),
        'p_value': p_value,
        'final_score': final_score,
        'mean_difference': np.mean(actual_values - random_values)
    }
    
    print(f"{setting_name}: {result}")
    return result

def main(par):
    evaluation_adata = ad.read_h5ad(par['evaluation_data'])
    layer = manage_layer(evaluation_adata, par)
    evaluation_adata.X = evaluation_adata.layers[layer]
    net = read_prediction(par) # [source	target	weight]

    group = par['group']
    for g in group:
        assert g in evaluation_adata.obs.columns

    # Evaluate three settings as specified in TODO
    settings_results = []
    
    # 1. Recall setting: all TFs, all target genes (no filtering)
    recall_result = evaluate_setting(
        evaluation_adata, net, group, 
        setting_name="recall", 
        n_tfs=None,  # Use all TFs
        max_targets_per_tf=None  # Use all target genes (up to default limit)
    )
    settings_results.append(recall_result)
    
    # 2. Balanced setting: top 100 TFs, top 100 target genes
    balanced_result = evaluate_setting(
        evaluation_adata, net, group, 
        setting_name="balanced", 
        n_tfs=100, 
        max_targets_per_tf=100
    )
    settings_results.append(balanced_result)
    
    # 3. Precision setting: top 20 TFs, top 20 target genes
    precision_result = evaluate_setting(
        evaluation_adata, net, group, 
        setting_name="precision", 
        n_tfs=20, 
        max_targets_per_tf=20
    )
    settings_results.append(precision_result)
    
    # Create results DataFrame with all three scores
    results = pd.DataFrame({
        'recovery_recall': [recall_result['final_score']],
        'recovery_balanced': [balanced_result['final_score']],
        'recovery_precision': [precision_result['final_score']],
    })
    
    return results


