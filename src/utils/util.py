import pandas as pd 
import anndata as ad 
import numpy as np 
from tqdm import tqdm
import scanpy as sc 

def process_links(net, par):
    net = net[net.source!=net.target]
    net_sorted = net.reindex(net['weight'].abs().sort_values(ascending=False).index)
    net = net_sorted.head(par['max_n_links']).reset_index(drop=True)
    return net

def corr_net(X, gene_names, par, tf_all=None):
    X = StandardScaler().fit_transform(X)
    net = np.dot(X.T, X) / X.shape[0]
    net = pd.DataFrame(net, index=gene_names, columns=gene_names)  
    if tf_all is None:  
        net = net.sample(n_tfs, axis=1, random_state=par['seed'])
    else:
        net = net[tf_all]
    net = net.reset_index()
    index_name = net.columns[0]
    net = net.melt(id_vars=index_name, var_name='source', value_name='weight')
    
    net.rename(columns={index_name: 'target'}, inplace=True)
    net = process_links(net, par)
    
    return net

def create_corr_net(X, gene_names, groups, par, tf_all):
    if par['cell_type_specific']:
        i = 0
        for group in tqdm(np.unique(groups), desc="Processing groups"):
            X_sub = X[groups == group, :]
            net = corr_net(X_sub, gene_names, par)
            net['cell_type'] = group
            if i==0:
                grn = net
            else:
                grn = pd.concat([grn, net], axis=0).reset_index(drop=True)
            i += 1
    else:
        grn = corr_net(X, gene_names, par)    
    return grn
def process_data(adata, par):
    if par['normalize']:
        print('Noramlize data')
        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)
        sc.pp.scale(adata)
    adata.X = adata.X.todense()
    if par['only_hvgs']:
        print('Subsetting data to hvgs')
        adata = adata[:, adata.var.hvg_counts]
        print('New dimension of data: ', adata.shape)