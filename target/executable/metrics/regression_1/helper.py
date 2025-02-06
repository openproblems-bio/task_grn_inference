import pandas as pd
import anndata as ad
import sys
import numpy as np
from sklearn.metrics import r2_score
from sklearn.pipeline import Pipeline
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler
import lightgbm
import random 
from concurrent.futures import ThreadPoolExecutor
from tqdm import tqdm
from sklearn.multioutput import MultiOutputRegressor
import os
import warnings

def set_global_seed(seed):
    np.random.seed(seed)
    random.seed(seed)
    lightgbm.LGBMRegressor().set_params(random_state=seed)


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

    # initialize baselines
    n_baselines = 100
    unique_groups = np.unique(groups)
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
    r2score_n_all = r2score - np.mean(r2score_baselines)

    # - normalized r2 score - S2
    r2score = r2_score(y_true[mask_gene_net, :], y_pred[mask_gene_net, :], multioutput='variance_weighted')
    r2score_baselines = []
    for y_pred_baseline in y_pred_baselines:
        r2score_baselines.append(r2_score(y_true[mask_gene_net, :], y_pred_baseline[mask_gene_net, :], multioutput='variance_weighted'))
    r2score_n_grn = r2score - np.mean(r2score_baselines)
    
    # - normalized r2 scores per sample
    r2score_samples = []
    for i_sample in range(y_true.shape[1]):
        score_sample = r2_score(y_true[:, i_sample], y_pred[:, i_sample])
        r2score_baselines = []
        for y_pred_baseline in y_pred_baselines:
            r2score_baselines.append(r2_score(y_true[:, i_sample], y_pred_baseline[:, i_sample]))
        r2score_samples.append(score_sample - np.mean(r2score_baselines))
    
    results = {
        'r2score-aver-all': r2score_n_all,
        'r2score-aver-grn': r2score_n_grn,
        'r2scores': r2score_samples,
        'reg_models' : reg_models
    }
    
    return results

