import scipy
import numpy as np
import pandas as pd
import anndata as ad
import warnings
warnings.filterwarnings("ignore")
from util import format_save_score, read_prediction


def tf_activity_local(net, adata, tf_all=None, shuffle=False):
    net = net.pivot(index='source', columns='target', values='weight').fillna(0)
    if shuffle:
        # Shuffle the weights across the matrix (null model)
        flat = net.values.flatten()
        np.random.shuffle(flat)
        net = pd.DataFrame(
            flat.reshape(net.shape),
            index=net.index,
            columns=net.columns
        )
    net = net[[g for g in adata.var_names if g in net.columns]]
    if tf_all is not None:
        tfs_present = np.intersect1d(net.index, tf_all)
    else:
        tfs_present = net.index
    net = net[net.index.isin(tfs_present)]
    adata = adata[:, adata.var_names.isin(net.columns)]
    X = adata.X
    X = X.toarray() if scipy.sparse.issparse(X) else X
    mat = X.T
    tf_acts = np.dot(net, mat)
    tf_acts = pd.DataFrame(tf_acts, index=net.index, columns=adata.obs.index)
    tf_acts = tf_acts.reset_index().melt(id_vars='source', var_name='sample', value_name='activity')
    cols = adata.obs.columns
    tf_acts = tf_acts.set_index('sample').merge(adata.obs[cols], left_index=True, right_index=True).reset_index(drop=False)
    if 'sample' not in tf_acts.columns:
        tf_acts['sample'] = tf_acts['index']
    return tf_acts
def calculate_tf_activity(adata, net, tf_all=None, n_targets_t=1, **kwargs):  
    import pandas as pd  
    if tf_all is not None:
        net = net[net['source'].isin(tf_all)]
    tf_size = net.groupby('source').size()
    tfs = tf_size[tf_size > n_targets_t].index
    net = net[net['source'].isin(tfs)]
    if False:
        tf_acts = tf_activity_local(net, adata, tf_all, **kwargs)
    else:
        import decoupler
        adata.obs['sample'] = adata.obs.index.astype(str)
        mat = pd.DataFrame(
            data=adata.X.todense(),  
            columns=adata.var_names,  
            index=adata.obs.index  
        )

        tf_acts, tf_pvals = decoupler.run_ulm(mat, net, source='source', target='target', weight='weight', use_raw=False)
        # - formatize
        tf_acts = tf_acts.reset_index().melt(id_vars='index', var_name='source', value_name='activity')
     
        # obs = adata.obs[cols]
        obs = adata.obs.copy()
        obs = obs.reset_index()
        
        tf_acts['index'] = tf_acts['index'].astype(str)
        obs['index'] = obs['index'].astype(str)
        tf_acts = tf_acts.merge(obs, on='index', how='left').drop('index', axis=1)
        assert tf_acts.shape[0]==tf_acts.shape[0]
    if 'index' in tf_acts.columns:
        tf_acts = tf_acts.drop('index', axis=1)
    
    assert tf_acts.shape[0] != 0, 'Empty'
    assert 'sample' in tf_acts.columns, 'sample not in tf_acts' 
    if False:
        tf_acts['sample'] = tf_acts['sample'].astype(int)
    X_df = tf_acts.pivot(index='sample', columns='source', values='activity')
    obs_df = tf_acts.drop_duplicates(subset='sample').set_index('sample')[[c for c in tf_acts.columns if c not in ['source', 'activity', 'sample']]]
    obs_df = obs_df.loc[X_df.index]
    var_df = pd.DataFrame(index=X_df.columns)
    var_df['source'] = var_df.index
    tf_acts_adata = ad.AnnData(X=X_df.values, obs=obs_df, var=var_df)
    return tf_acts_adata

from scipy.stats import spearmanr
def recovery_consistency_score(tf_activity, group):
    df_tfs = tf_activity.to_df().join(tf_activity.obs[group])
    results = []
    for g_vals, df_group in df_tfs.groupby(group):
        df_feats = df_group.drop(columns=group)
        if len(df_feats) < 2:
            continue  # skip if only one sample
        sims = []
        rows = df_feats.values
        for i in range(len(rows)):
            for j in range(i + 1, len(rows)):
                rho, _ = spearmanr(rows[i], rows[j])
                sims.append(rho)
        mean_sim = np.nanmean(sims) if sims else np.nan
        results.append(list(g_vals) + [mean_sim])
    return pd.DataFrame(results, columns=group + ["consistency_score"])
    

def main(par):
    n_tfs = 100
    evaluation_adata = ad.read_h5ad(par['evaluation_data'])
    net = read_prediction(par) # [source	target	weight]

    c = net.groupby('source').size().sort_values(ascending=False).head(n_tfs)
    tfs = c.index
    net = net[net['source'].isin(tfs)]
    net_shuffle = net.copy()
    np.random.seed(42)
    net_shuffle['target'] = np.random.permutation(net_shuffle['target'].values)
    net_shuffle = net_shuffle.drop_duplicates(subset=['source', 'target'])


    tf_act = calculate_tf_activity(evaluation_adata, net, n_targets_t=10)
    tf_act_random = calculate_tf_activity(evaluation_adata, net_shuffle, n_targets_t=10)

    group = par['group']
    for g in group:
        assert g in tf_act.obs.columns

    scores = recovery_consistency_score(tf_activity=tf_act, group=group)['consistency_score']    
    scores_baseline = recovery_consistency_score(tf_activity=tf_act_random, group=group)['consistency_score']


    from scipy.stats import ttest_rel, spearmanr, wilcoxon
    try:
        res = wilcoxon(scores - scores_baseline, zero_method='wilcox', alternative='greater')
    except Exception as e:
        print('Error in wilcoxon test:', e)
        res = None
    if res is not None:
        score = np.median(scores)
        score_baseline = np.median(scores_baseline)
        print({'res.pvalue': res.pvalue, 'score:': score , 'scores_baseline': score_baseline})

        eps = 1e-300  # very small number to avoid log(0)
        pval_clipped = max(res.pvalue, eps)

        effect = score - score_baseline
        score = effect * (-np.log10(pval_clipped))

        results = {
            'recovery_2': [float(score)]
        }
        results = pd.DataFrame(results)
    else:
        results = {
            'recovery_2': [float('nan')]
        }
        results = pd.DataFrame(results)
    return results


