import scipy
import numpy as np
import pandas as pd
import anndata as ad
import warnings
warnings.filterwarnings("ignore")
from util import format_save_score, read_prediction

from scipy.stats import spearmanr

from util import manage_layer
USE_MAE = False
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
    group_consistencies = []
    for g_vals, df_group in target_expr.groupby(group):
        expr_data = df_group.drop(columns=group)
        if len(expr_data) < 2:
            continue  # skip if only one sample
        sims = []
        rows = expr_data.values
        for i in range(len(rows)):
            for j in range(i + 1, len(rows)):
                if USE_MAE:
                    # pairwise absolute mean error per gene
                    mae = np.abs(rows[i] - rows[j])
                    normed = mae / (np.mean([rows[i], rows[j]], axis=0) + 1e-8)
                    score = np.mean(normed)
                else:
                    # correlation, _ = spearmanr(rows[i], rows[j])
                    correlation = np.corrcoef(rows[i], rows[j])[0, 1]
                    assert isinstance(correlation, float)
                    score = correlation
                sims.append(score)
        group_consistencies.append(np.mean(sims))
    return group_consistencies

def calculate_target_gene_consistency(adata, net, group, min_targets_t=10, max_targets_t=50):
    """
    Calculate consistency score for each TF based on target gene expression across groups.
    Returns a dictionary with TF as key and consistency score as value.
    """
    adata = adata[adata.obs['is_control'] == False]
    all_genes = adata.var_names.tolist()
    tf_size = net.groupby('source').size()
    tfs = tf_size[(tf_size > min_targets_t)].index
    net = net[net['source'].isin(tfs)]
    
    tf_consistency_scores = []
    tf_consistency_scores_random = []
    
    for tf in tfs:
        target_genes = net[net['source'] == tf]
        target_genes = target_genes.sort_values(by='weight', ascending=False, key=abs).head(max_targets_t)
        target_genes = target_genes['target'].tolist()
            
        actual_consistencies = actual_calculation(adata, target_genes, group)
        tf_consistency_scores.append(np.mean(actual_consistencies))

        n_targets = len(target_genes)
        random_targets = np.random.choice(all_genes, size=n_targets, replace=False)
        random_consistencies = actual_calculation(adata, random_targets, group)
        tf_consistency_scores_random.append(np.mean(random_consistencies))
    
    return tf_consistency_scores, tf_consistency_scores_random

def main(par):
    evaluation_adata = ad.read_h5ad(par['evaluation_data'])
    layer = manage_layer(evaluation_adata, par)
    evaluation_adata.X = evaluation_adata.layers[layer]
    net = read_prediction(par) # [source	target	weight]

    if True:
        n_tfs = 20
        c = net.groupby('source').size().sort_values(ascending=False).head(n_tfs)
        tfs = c.index
        net = net[net['source'].isin(tfs)]

    group = par['group']
    for g in group:
        assert g in evaluation_adata.obs.columns

    
    scores, scores_baseline = calculate_target_gene_consistency(
        evaluation_adata, net, group
    )
    
    from scipy.stats import ttest_rel, spearmanr, wilcoxon
    if USE_MAE:
        res = wilcoxon(np.array(scores_baseline)- np.array(scores), zero_method='wilcox', alternative='greater')
    else:
        res = wilcoxon(np.array(scores_baseline) - np.array(scores), zero_method='wilcox', alternative='less')
    p_value = res.pvalue

    if USE_MAE:
        print({
            'consistency_error': np.mean(scores),
            'consistency_error_random': np.median(scores_baseline),
            'p_value': p_value,
        })
    else:
        print({
            'consistency_score': np.mean(scores),
            'consistency_score_random': np.mean(scores_baseline),
            'p_value': p_value,
        })

    results = pd.DataFrame({
        'score': [-np.log10(p_value)],
    })
    
    return results


