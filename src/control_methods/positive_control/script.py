import os
import pandas as pd
import numpy as np
import anndata as ad
import scanpy as sc
from tqdm import tqdm
from scipy.stats import spearmanr
from sklearn.preprocessing import StandardScaler

## VIASH START
par = {
    'perturbation_data': 'resources/grn-benchmark/perturbation_data.h5ad',
    'tf_all': 'resources/prior/tf_all.csv',
    'causal': True,
    'cell_type_specific': True,
    'max_n_links': 50000,
    'prediction': 'resources/grn_models/positive_control.csv',
    "seed": 32
}
## VIASH END
print(par)

def process_links(net, par):
    net = net[net.source!=net.target]
    net_sorted = net.reindex(net['weight'].abs().sort_values(ascending=False).index)
    net = net_sorted.head(par['max_n_links']).reset_index(drop=True)
    return net

def corr_net(X, gene_names, par):
    X = StandardScaler().fit_transform(X)
    net = np.dot(X.T, X) / X.shape[0]
    net = pd.DataFrame(net, index=gene_names, columns=gene_names)
    
    net = net[tf_all]
    net = net.reset_index()
    index_name = net.columns[0]
    net = net.melt(id_vars=index_name, var_name='source', value_name='weight')
    
    net.rename(columns={index_name: 'target'}, inplace=True)
    net = process_links(net, par)
    
    return net

def create_corr_net(X, gene_names, groups, par):
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
print('Read data')
multiomics_rna = ad.read_h5ad(par["perturbation_data"])
 

if par['normalize']:
    print('Noramlize data')
    sc.pp.normalize_total(multiomics_rna)
    sc.pp.log1p(multiomics_rna)
    sc.pp.scale(multiomics_rna)
    
gene_names = multiomics_rna.var_names.to_numpy()
tf_all = np.loadtxt(par['tf_all'], dtype=str)
groups = multiomics_rna.obs.cell_type
tf_all = np.intersect1d(tf_all, gene_names)


print('Create corr net')
net = create_corr_net(multiomics_rna.X, multiomics_rna.var_names, groups, par)

print('Output GRN')
net.to_csv(par['prediction'])
