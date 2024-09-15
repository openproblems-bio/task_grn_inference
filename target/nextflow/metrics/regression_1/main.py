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

import os


def select_top_links(net, par):
    print("Number of links reduced to ", par['max_n_links'])
    net_sorted = net.reindex(net['weight'].abs().sort_values(ascending=False).index)
    net = net_sorted.head(par['max_n_links']).reset_index(drop=True)
    return net

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


def cv_5(genes_n):
    '''5 fold standard'''
    num_groups = 5
    group_size = genes_n // num_groups
    groups = np.repeat(np.arange(num_groups), group_size)
    if genes_n % num_groups != 0:
        groups = np.concatenate((groups, np.arange(genes_n % num_groups)))
    np.random.shuffle(groups)
    return groups

def cross_validation(net, prturb_adata, par:dict):
    gene_names = prturb_adata.var_names
    net = process_net(net.copy(), gene_names)
    net_subset = net.copy()
    # Subset TFs 
    if par['tf_n'] == -1:
        degrees = net.abs().sum(axis=0)
        net = net.loc[:, degrees >= degrees.quantile((1 - par['theta']))]
    else:
        degrees = net.abs().sum(axis=0)
        net = net[degrees.nlargest(tf_n).index]

    gene_names_grn = net.index.to_numpy()   
    
    n_tfs = net.shape[1]
    # construct feature and target space
    if par['exclude_missing_genes']:
        included_genes = gene_names_grn
    else:
        included_genes = gene_names
    
    X_df = pd.DataFrame(np.zeros((len(included_genes), n_tfs)), index=included_genes)
    train_df = pd.DataFrame(prturb_adata.X, columns=gene_names).T
    Y_df = train_df.loc[included_genes,:]

    mask_shared_genes = X_df.index.isin(net.index)
    
    # fill the actuall regulatory links
    X_df.loc[mask_shared_genes, :] = net.values
    X = X_df.values.copy()

    # run cv 
    groups = cv_5(len(included_genes))
    # initialize y_pred with the mean of gene expressed across all samples
    means = Y_df.mean(axis=0)
    y_pred = Y_df.copy()
    y_pred[:] = means
    y_pred = y_pred.values

    # initialize y_true
    Y = Y_df.values
    y_true = Y.copy()

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
        regr = lightgbm_wrapper(params, max_workers=max_workers)
    elif reg_type=='RF':
        params = dict(boosting_type='rf',random_state=32, n_estimators=100,  
        feature_fraction=0.05, verbosity=-1)
        regr = lightgbm_wrapper(params, max_workers)
    else:
        print(f'{reg_type} is not defined')
        raise ValueError("define first")  

    for group in tqdm(unique_groups, desc="Processing groups"):
        mask_va = groups == group
        mask_tr = ~mask_va
        X_tr = X[mask_tr & mask_shared_genes, :]
        Y_tr = Y[mask_tr & mask_shared_genes, :]
        regr.fit(X_tr, Y_tr)
        y_pred[mask_va & mask_shared_genes, :] = regr.predict(X[mask_va & mask_shared_genes, :])
    return y_true, y_pred

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
    cell_types = prturb_adata.obs.cell_type.unique()
    score_list = []
    for cell_type in cell_types:
        print(f'----cross validate for {cell_type}----')
        # check if net is cell type specific 
        if 'cell_type' in net.columns:
            if cell_type not in net.cell_type.unique():
                raise ValueError(f'{cell_type} is not present in grn.')
            net_sub = net[net.cell_type==cell_type]
        else:
            net_sub = net
        if len(net_sub)>par['max_n_links']:
            net_sub = select_top_links(net_sub, par) #only top links are considered
                
        prturb_adata_sub = prturb_adata[prturb_adata.obs.cell_type==cell_type,:]
        y_true_sub, y_pred_sub = cross_validation(net_sub, prturb_adata_sub, par)

        score = r2_score(y_true_sub, y_pred_sub, multioutput='variance_weighted')
        score_list.append(score)
    mean_score_r2 = np.mean(score_list)
    output = dict(mean_score_r2=mean_score_r2)
    return output

def set_global_seed(seed):
    np.random.seed(seed)
    random.seed(seed)
    lightgbm.LGBMRegressor().set_params(random_state=seed)


def pivot_grn(net):
    ''' make net to have gene*tf format'''
    df_tmp = net.pivot(index='target', columns='source', values='weight')
    return df_tmp.fillna(0)

def degree_centrality(net, source='source', target='target', normalize=False):
    '''calculates centrality score of source'''
    counts = net.groupby(source)[target].nunique().values
    if normalize:
        total_targets = net[target].nunique()
        counts = counts / total_targets
    return counts


def process_net(net, gene_names):
    # Remove self-regulations
    net = net[net['source'] != net['target']]
    # pivot
    net = pivot_grn(net)
    # subset 
    net = net[net.index.isin(gene_names)]
    return net


def main(par):
    random_state = 42
    set_global_seed(random_state)
    par['theta'] = 1 # no subsetting based on theta
    manipulate = None 
    ## read and process input data

    print('Reading input files', flush=True)
    
    perturbation_data = ad.read_h5ad(par['perturbation_data'])
    tf_all = np.loadtxt(par['tf_all'], dtype=str)
    gene_names = perturbation_data.var.index.to_numpy()
    net = pd.read_csv(par['prediction'])
    # net['weight'] = net.weight.abs()
    # subset to keep only those links with source as tf
    if par['apply_tf']:
        net = net[net.source.isin(tf_all)]

    subsample = par['subsample']
    max_workers = par['max_workers']
    layer = par["layer"]
    if subsample == -1:
        pass
    elif subsample == -2: # one combination of cell_type, sm_name
        sampled_obs = perturbation_data.obs.groupby(['sm_name', 'cell_type'], observed=False).apply(lambda x: x.sample(1)).reset_index(drop=True)
        obs = perturbation_data.obs
        mask = []
        for _, row in obs.iterrows():
            mask.append((sampled_obs==row).all(axis=1).any())  
        perturbation_data = perturbation_data[mask,:]
    elif subsample == -3: #negative control
        mask = perturbation_data.obs.sm_name == 'Dimethyl Sulfoxide'
        perturbation_data = perturbation_data[mask,:]
    elif subsample == -4: #positive control
        mask = perturbation_data.obs.sm_name.isin(['Dabrafenib', 'Belinostat'])
        perturbation_data = perturbation_data[mask,:]
    else:
        perturbation_data = perturbation_data[np.random.choice(perturbation_data.n_obs, subsample, replace=False), :]
    
    print(perturbation_data.shape)
    
    perturbation_data.X = perturbation_data.layers[layer]

    

    print(f'Compute metrics for layer: {layer}', flush=True)
    tfs_cases = [-1]
    if par['min_tf']:
        tfs_cases += 30
    layer_results = {}  # Store results for this layer
    for exclude_missing_genes in [False, True]:  # two settings on target gene
        for tf_n in tfs_cases:  # two settings on tfs
            par['exclude_missing_genes'] = exclude_missing_genes
            par['tf_n'] = tf_n
            run_key = f'ex({exclude_missing_genes})_tf({tf_n})'
            print(run_key)
            
            output = regression_1(net, perturbation_data, par=par, max_workers=max_workers)
            
            layer_results[run_key] = [output['mean_score_r2']]

    # Convert results to DataFrame
    df_results = pd.DataFrame(layer_results)
    # if par['min_tf']:
    #     if 'ex(True)_tf(140)' not in df_results.columns:
    #         df_results['ex(True)_tf(140)'] = df_results['ex(True)_tf(-1)']
    #     if 'ex(False)_tf(140)' not in df_results.columns:
    #         df_results['ex(False)_tf(140)'] = df_results['ex(False)_tf(-1)']
        
    df_results['Mean'] = df_results.mean(axis=1)
    
    return df_results