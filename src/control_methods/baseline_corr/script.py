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
    'multiomics_rna': 'resources/grn-benchmark/multiomics_rna.h5ad',
    'tf_all': 'resources/prior/tf_all.csv',
    'causal': False,
    'metacell': False,
    'cell_type_specific': False,
    'impute': False,
    'max_n_links': 50000,
    'corr_method': 'pearson',
    'prediction': 'resources/grn_models/alldonors_default/pearson.csv',
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
    if par['corr_method'] == "pearson":
        X = StandardScaler().fit_transform(X)
        net = np.dot(X.T, X) / X.shape[0]
    elif par['corr_method'] == "spearman":
        net = spearmanr(X).statistic
    net = pd.DataFrame(net, index=gene_names, columns=gene_names)

    if par['causal']:
        net = net[tf_all]
    else:
        net = net.sample(len(tf_all), axis=1, random_state=par['seed'])
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
        # grn.drop(columns=['cell_type'], inplace=True)
        # grn = grn.groupby(['source', 'target']).mean().reset_index()
        # grn = process_links(grn, par)        
    return grn
print('Read data')
multiomics_rna = ad.read_h5ad(par["multiomics_rna"])
# multiomics_rna = multiomics_rna[:,:5000] #TODO: togo
 
if par['metacell']:
    print('metacell')
    def create_meta_cells(df, n_cells=15):
        meta_x = []
        for i in range(0, df.shape[0], n_cells):
            meta_x.append(df.iloc[i:i+n_cells, :].sum(axis=0).values)
        df = pd.DataFrame(meta_x, columns=df.columns)
        return df
            
    adata_df = pd.DataFrame(multiomics_rna.X.todense(), columns=multiomics_rna.var_names)
    adata_df['cell_type'] = multiomics_rna.obs['cell_type'].values
    adata_df['donor_id'] = multiomics_rna.obs['donor_id'].values
    df = adata_df.groupby(['cell_type','donor_id']).apply(lambda df: create_meta_cells(df))
    X = df.values
    var = pd.DataFrame(index=df.columns)
    obs = df.reset_index()[['cell_type','donor_id']]
    multiomics_rna = ad.AnnData(X=X, obs=obs, var=var)

gene_names = multiomics_rna.var_names.to_numpy()
tf_all = np.loadtxt(par['tf_all'], dtype=str)
groups = multiomics_rna.obs.cell_type
tf_all = np.intersect1d(tf_all, gene_names)

print('Noramlize data')
sc.pp.normalize_total(multiomics_rna)
sc.pp.log1p(multiomics_rna)
sc.pp.scale(multiomics_rna)

if par['impute']:
    print("imputing")
    import magic
    import scprep
    
    magic_operator = magic.MAGIC()
    
    multiomics_rna = magic_operator.fit_transform(multiomics_rna)
    
    print('zero ration: ', (multiomics_rna.X==0).sum()/multiomics_rna.X.size)
print('Create corr net')
net = create_corr_net(multiomics_rna.X, multiomics_rna.var_names, groups, par)

print('Output GRN')
net.to_csv(par['prediction'])
