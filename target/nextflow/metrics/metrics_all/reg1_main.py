import numpy as np
from sklearn.metrics import r2_score
from sklearn.pipeline import Pipeline
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler
import lightgbm
import random 
from concurrent.futures import ThreadPoolExecutor
import pandas as pd
import anndata as ad
from tqdm import tqdm
from util import verbose_print, process_links, verbose_tqdm
from sklearn.multioutput import MultiOutputRegressor
import os
import warnings

warnings.filterwarnings("ignore", category=FutureWarning)



class lightgbm_wrapper:
    def __init__(self, params, max_workers=None):
        self.params =  params
        self.max_workers = max_workers
        
    def fit(self, X_train, Y_train):
        self.n_sample = Y_train.shape[1]
        self.regr_samples = [None] * self.n_sample
        
        def fit_sample(i):
            regr = lightgbm.LGBMRegressor(**self.params)
            regr.fit(X_train, Y_train[:, i])
            self.regr_samples[i] = regr
        
        with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            executor.map(fit_sample, range(self.n_sample))
            
    def predict(self, X_test):
        def predict_sample(i):
            regr = self.regr_samples[i]
            return regr.predict(X_test)
        
        with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            y_pred_list = list(executor.map(predict_sample, range(self.n_sample)))
        
        return np.stack(y_pred_list, axis=1)
    def get_feature_importance(self):
        # Ensure models have been fitted
        if not all(self.regr_samples):
            raise ValueError("Models are not yet trained. Call `fit` first.")

        # Collect feature importances from each regressor
        importances = np.array([regr.feature_importances_ for regr in self.regr_samples])
        
        # Optionally, return the importance for each sample as well
        return importances


def cv_5(genes_n):
    '''5 fold standard'''
    np.random.seed(32)
    num_groups = 5
    group_size = genes_n // num_groups
    groups = np.repeat(np.arange(num_groups), group_size)
    if genes_n % num_groups != 0:
        groups = np.concatenate((groups, np.arange(genes_n % num_groups)))
    np.random.shuffle(groups)
    return groups

def cross_validation(net, prturb_adata, par:dict):
    np.random.seed(32)
    
    gene_names = prturb_adata.var_names
    
    net = process_net(net.copy(), gene_names)

    gene_names_grn = net.index.to_numpy()   
    
    n_tfs = net.shape[1]
    # construct feature and target space
    X_df = pd.DataFrame(np.zeros((len(gene_names), n_tfs)), index=gene_names)
    
    try:
        train_df = pd.DataFrame(prturb_adata.X, columns=gene_names).T
    except:
        train_df = pd.DataFrame(prturb_adata.X.todense().A, columns=gene_names).T
    Y_df = train_df.loc[gene_names,:]

    mask_gene_net = X_df.index.isin(net.index)
    
    # fill the actual regulatory links
    X_df.loc[mask_gene_net, :] = net.values
    X = X_df.values.copy()
    

    # run cv 
    groups = cv_5(len(gene_names))
    # initialize y_pred with the mean of gene expressed across all samples
    y_pred = Y_df.copy()
    means = Y_df.mean(axis=0)
    y_pred[:] = means
    y_pred = y_pred.values

    # initialize y_true
    Y = Y_df.values
    y_true = Y.copy()

    # print(r2_score(y_true, y_pred, multioutput='variance_weighted'))

    unique_groups = np.unique(groups)

    # determine regressor 
    reg_type = par['reg_type']
    if reg_type=='ridge':
        regr =  Ridge(**dict(random_state=32))
    elif reg_type=='GB':
        params = dict(random_state=32, 
                    n_estimators=100, min_samples_leaf=2, min_child_samples=1, 
                    feature_fraction=0.05, verbosity=-1
        )
        regr = lightgbm_wrapper(params, max_workers=par['num_workers'])
    elif reg_type=='RF':
        params = dict(boosting_type='rf',random_state=32, n_estimators=100,  
        feature_fraction=0.05, verbosity=-1)
        regr = lightgbm_wrapper(params, max_workers=par['num_workers'])
    else:
        raise ValueError(f"{reg_type} is not defined.")  
    reg_models = []

    n_baselines = 100
    y_pred_baselines = [y_pred.copy() for i in range(n_baselines)]
    for group in verbose_tqdm(unique_groups, "Processing groups", 2, par['verbose']):
        mask_va = groups == group
        mask_tr = ~mask_va
        X_tr = X[mask_tr & mask_gene_net, :]
        Y_tr = Y[mask_tr & mask_gene_net, :]

        X_va = X[mask_va & mask_gene_net, :]

        if X_tr.shape[0]<2:
            continue 
        regr.fit(X_tr, Y_tr)
        Y_pr = regr.predict(X_va)

        y_pred[mask_va & mask_gene_net, :] = Y_pr.copy()
        
        # Shuffle each column independently and keep the same shape
        for i in range(n_baselines):
            for col in range(Y_pr.shape[1]):
                np.random.shuffle(Y_pr[:, col])
            y_pred_baselines[i][mask_va & mask_gene_net, :] = Y_pr
        reg_models.append(regr)
    
    # - normalized r2 score
    r2score = r2_score(y_true, y_pred, multioutput='variance_weighted')
    r2score_baselines = []
    for y_pred_baseline in y_pred_baselines:
        r2score_baselines.append(r2_score(y_true, y_pred_baseline, multioutput='variance_weighted'))
    r2score_n = r2score - np.mean(r2score_baselines)

    # - normalized r2 score - S2
    r2score = r2_score(y_true[mask_gene_net, :], y_pred[mask_gene_net, :], multioutput='variance_weighted')
    r2score_baselines = []
    for y_pred_baseline in y_pred_baselines:
        r2score_baselines.append(r2_score(y_true[mask_gene_net, :], y_pred_baseline[mask_gene_net, :], multioutput='variance_weighted'))
    r2score_n_s2 = r2score - np.mean(r2score_baselines)
    
    # - normalized r2 scores per sample
    r2score_samples = []
    for i_sample in range(y_true.shape[1]):
        score_sample = r2_score(y_true[:, i_sample], y_pred[:, i_sample])
        r2score_baselines = []
        for y_pred_baseline in y_pred_baselines:
            r2score_baselines.append(r2_score(y_true[:, i_sample], y_pred_baseline[:, i_sample]))
        r2score_samples.append(score_sample - np.mean(r2score_baselines))
    

    results = {
        'r2score-aver': r2score_n,
        'r2score-aver-s2': r2score_n_s2,
        'r2scores': r2score_samples,
        'reg_models' : reg_models
    }
    
    
    return results

def regression_1(
            net: pd.DataFrame, 
            prturb_adata: ad.AnnData,
            par:dict,
            verbose: int = 0,
            max_workers: int = 4) -> None: 
    """
    net: a df with index as genes and columns as tfs
    prturb_adata: a adata 
    """
       
    gene_names = prturb_adata.var_names
    if 'donor_id' in prturb_adata.obs.columns:
        donor_ids = prturb_adata.obs.donor_id.unique()
    else:
        donor_ids = ['default']
    score_list = []
    score_s2_list = []
    for donor_id in donor_ids:
        verbose_print(par['verbose'], f'----cross validate for {donor_id}----', 2)
        # check if net is donor specific 
        if 'donor_id' in net.columns:
            if donor_id not in net.donor_id.unique():
                raise ValueError(f'{donor_id} is not present in grn.')
            net_sub = net[net.donor_id==donor_id]
        else:
            net_sub = net
        if (len(net_sub)>par['max_n_links']) and (par['max_n_links']!=-1):
            net_sub = process_links(net_sub, par) #only top links are considered
            verbose_print(par['verbose'], f"Number of links reduced to {par['max_n_links']}", 2)
        if par['binarize']:
            net['weight'] = net['weight'].apply(binarize_weight)        
        if 'donor_id' in prturb_adata.obs.columns:
            prturb_adata_sub = prturb_adata[prturb_adata.obs.donor_id==donor_id,:]
        else:
            prturb_adata_sub = prturb_adata
        results = cross_validation(net_sub, prturb_adata_sub, par)

        
        
        score_list.append(results['r2score-aver'])
        score_s2_list.append(results['r2score-aver-s2'])
    s1 = [np.mean(score_list)]
    s2 = [np.mean(score_s2_list)]
    
    output = dict(S1=s1, S2=s2)
    return output

def set_global_seed(seed):
    np.random.seed(seed)
    random.seed(seed)
    lightgbm.LGBMRegressor().set_params(random_state=seed)


def pivot_grn(net):
    ''' make net to have gene*tf format'''
    net = net.drop_duplicates(subset=['target', 'source'])
    df_tmp = net.pivot(index='target', columns='source', values='weight')
    return df_tmp.fillna(0)


def process_net(net, gene_names):
    # Remove self-regulations
    net = net.drop_duplicates()
    net = net[net['source'] != net['target']]
    # pivot
    net = pivot_grn(net)
    # subset 
    net = net[net.index.isin(gene_names)]
    return net

def binarize_weight(weight):
    if weight > 0:
        return 1
    elif weight < 0:
        return -1
    else:
        return 0
    
def main(par):
    random_state = 42
    set_global_seed(random_state)
    par['theta'] = 1 # no subsetting based on theta
    manipulate = None 
    ## read and process input data
    verbose_print(par['verbose'], 'Reading input files', 3)
    net = ad.read_h5ad(par['prediction'])
    net = pd.DataFrame(net.uns['prediction'])

    evaluation_data = ad.read_h5ad(par['evaluation_data'])

    if 'is_test' in evaluation_data.obs.columns:
        evaluation_data = evaluation_data[evaluation_data.obs['is_test']]
        print(evaluation_data.shape)

    tf_all = np.loadtxt(par['tf_all'], dtype=str)
    gene_names = evaluation_data.var.index.to_numpy()
    
    
    if par['apply_skeleton']: #apply skeleton
        print('Before filtering with skeleton:', net.shape)
        skeleton = pd.read_csv(par['skeleton'])
        skeleton['link'] = skeleton['source'].astype(str) + '_' + skeleton['target'].astype(str)
        skeleton = skeleton['link'].values.flatten()
        
        net['link'] = net['source'].astype(str) + '_' + net['target'].astype(str)
        net = net[net['link'].isin(skeleton)]
        print('After filtering with skeleton:', net.shape)

    
    # net['weight'] = net.weight.abs()
    # subset to keep only those links with source as tf
    
    if par['apply_tf']:
        net = net[net.source.isin(tf_all)]

    if 'cell_type' in net.columns:
        # print('Taking mean of cell type specific grns')
        # net.drop(columns=['cell_type'], inplace=True)
        # net = net.groupby(['source', 'target']).mean().reset_index()
        raise ValueError('define this')

    subsample = par['subsample']
    max_workers = par['num_workers']
    layer = par["layer"]
    if subsample == -1:
        pass
    elif subsample == -2: # one combination of cell_type, perturbation
        sampled_obs = evaluation_data.obs.groupby(['perturbation', 'cell_type'], observed=False).apply(lambda x: x.sample(1)).reset_index(drop=True)
        obs = evaluation_data.obs
        mask = []
        for _, row in obs.iterrows():
            mask.append((sampled_obs==row).all(axis=1).any())  
        evaluation_data = evaluation_data[mask,:]
    elif subsample == -3: #negative control
        mask = evaluation_data.obs.perturbation == 'Dimethyl Sulfoxide'
        evaluation_data = evaluation_data[mask,:]
    elif subsample == -4: #positive control
        mask = evaluation_data.obs.perturbation.isin(['Dabrafenib', 'Belinostat'])
        evaluation_data = evaluation_data[mask,:]
    else:
        evaluation_data = evaluation_data[np.random.choice(evaluation_data.n_obs, subsample, replace=False), :]
    verbose_print(par['verbose'], evaluation_data.shape,4)   

    if layer=='X':
        pass 
    else:
        evaluation_data.X = evaluation_data.layers[layer]

    verbose_print(par['verbose'], f'Compute metrics for layer: {layer}',  3)    


    results = regression_1(net, evaluation_data, par=par, max_workers=max_workers)


    df_results = pd.DataFrame(results)
    return df_results