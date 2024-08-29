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


def regression_1(
            net: pd.DataFrame, 
            train_df: pd.DataFrame,
            reg_type: str = 'GB',
            exclude_missing_genes: bool = False,
            verbose: int = 0,
            max_workers: int = 4) -> None: 
    """
    net: a df with index as genes and columns as tfs
    train_df: a df with index as genes and columns as samples
    """
    gene_names = train_df.index.to_numpy()
    gene_names_grn = net.index.to_numpy()
    # determine regressor 
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
    
    n_tfs = net.shape[1]
    # construct feature and target space
    if exclude_missing_genes:
        included_genes = gene_names_grn
    else:
        included_genes = gene_names
    
    X_df = pd.DataFrame(np.zeros((len(included_genes), n_tfs)), index=included_genes)
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

    for group in tqdm(unique_groups, desc="Processing groups"):
        mask_va = groups == group
        mask_tr = ~mask_va
        # Use logical AND to combine masks correctly
        X_tr = X[mask_tr & mask_shared_genes, :]
        Y_tr = Y[mask_tr & mask_shared_genes, :]

        regr.fit(X_tr, Y_tr)
        y_pred[mask_va & mask_shared_genes, :] = regr.predict(X[mask_va & mask_shared_genes, :])

    score_r2  = r2_score(y_true, y_pred, multioutput='variance_weighted') #uniform_average', 'variance_weighted

    mean_score_r2 = r2_score(y_true, y_pred, multioutput='variance_weighted')
    # gene_scores_r2 = r2_score(y_true.T, y_pred.T, multioutput='raw_values')

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


def process_net(net, gene_names, manipulate):
    # Remove self-regulations
    net = net[net['source'] != net['target']]
    # pivot
    net = pivot_grn(net)
    # subset 
    net = net[net.index.isin(gene_names)]
    # sign or shuffle
    if manipulate=='signed':
        net = net.map(lambda x: 1 if x>0 else (-1 if x<0 else 0))
    return net


def main(par):
    random_state = 42
    set_global_seed(random_state)
    theta = 1 # no subsetting based on theta
    manipulate = None 
    ## read and process input data

    print('Reading input files', flush=True)
    
    perturbation_data = ad.read_h5ad(par['perturbation_data'])
    tf_all = np.loadtxt(par['tf_all'], dtype=str)
    gene_names = perturbation_data.var.index.to_numpy()
    net = pd.read_csv(par['prediction'])
    # subset to keep only those links with source as tf
    net = net[net.source.isin(tf_all)]

    subsample = par['subsample']
    reg_type = par['reg_type']
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

    pert_df = pd.DataFrame(perturbation_data.layers[layer], columns=gene_names)
    pert_df = pert_df.T  # make it gene*sample

    # process net
    net_processed = process_net(net.copy(), gene_names, manipulate)

    print(f'Compute metrics for layer: {layer}', flush=True)
    layer_results = {}  # Store results for this layer
    for exclude_missing_genes in [False]:  # two settings on target gene
        for tf_n in [-1]:  # two settings on tfs
            run_key = f'ex({exclude_missing_genes})_tf({tf_n})'
            print(run_key)
            net_subset = net_processed.copy()

            # Subset TFs 
            if tf_n == -1:
                degrees = net_subset.abs().sum(axis=0)
                net_subset = net_subset.loc[:, degrees >= degrees.quantile((1 - theta))]
            else:
                if tf_n > net_subset.shape[1]:
                    print(f'Skip running because tf_n ({tf_n}) is bigger than net.shape[1] ({net_subset.shape[1]})')
                    continue
                degrees = net_subset.abs().sum(axis=0)
                net_subset = net_subset[degrees.nlargest(tf_n).index]

            output = regression_1(net_subset, pert_df, exclude_missing_genes=exclude_missing_genes, reg_type=reg_type, max_workers=max_workers)
            
            layer_results[run_key] = [output['mean_score_r2']]

    # Convert results to DataFrame
    df_results = pd.DataFrame(layer_results)
    # if 'ex(True)_tf(140)' not in df_results.columns:
    #     df_results['ex(True)_tf(140)'] = df_results['ex(True)_tf(-1)']
    # if 'ex(False)_tf(140)' not in df_results.columns:
    #     df_results['ex(False)_tf(140)'] = df_results['ex(False)_tf(-1)']
    
    df_results['Mean'] = df_results.mean(axis=1)
    
    return df_results