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
}

print(par)

## VIASH END

def process_links(net, par):
    net = net[net.source!=net.target]
    net_sorted = net.reindex(net['weight'].abs().sort_values(ascending=False).index)
    net = net_sorted.head(par['max_n_links']).reset_index(drop=True)
    return net

def create_corr_net(X: np.ndarray, groups: np.ndarray, method="pearson"):
    i = 0
    for group in tqdm(np.unique(groups), desc="Processing groups"):
        X_sub = X[groups == group, :]
        if method == "dotproduct":
            X_sub = StandardScaler().fit_transform(X_sub)
            net = np.dot(X_sub.T, X_sub) / X_sub.shape[0]
        elif method == "pearson":
            net = np.corrcoef(X_sub.T)
            # net = pd.DataFrame(X_sub).transpose().corr().values.to_numpy()
            net = np.nan_to_num(net, nan=0.0, posinf=0.0, neginf=0.0)
        elif method == "spearman":
            net = spearmanr(X_sub).statistic

    
        net = pd.DataFrame(net, index=gene_names, columns=gene_names)
        if par['causal']:
            net = net[tf_all]
        else:
            net = net.sample(len(tf_all), axis=1, random_state=par['seed'])
        # flatten to source-target-weight    
        net = net.reset_index().melt(id_vars='index', var_name='source', value_name='weight')
        net.rename(columns={'index': 'target'}, inplace=True)
        
        net = process_links(net, par)
        net['cell_type'] = group
        if i==0:
            grn = net
        else:
            grn = pd.concat([grn, net], axis=0).reset_index(drop=True)

        i += 1
            
    if par['cell_type_specific']==False:
        grn.drop(columns=['cell_type'], inplace=True)
        grn = grn.groupby(['source', 'target']).mean().reset_index()
        grn = process_links(grn, par)        
    return grn
print('Read data')
multiomics_rna = ad.read_h5ad(par["multiomics_rna"])
multiomics_rna = multiomics_rna[:,:5000] #TODO: togo
 
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
    # import magic
    # import scprep
    
    # magic_operator = magic.MAGIC()
    
    # multiomics_rna = magic_operator.fit_transform(multiomics_rna)
    from sklearn.impute import KNNImputer
    import numpy as np

    print("Imputing with KNN")

    # Convert to dense if the matrix is sparse
    if sc.sparse.issparse(multiomics_rna.X):
        multiomics_rna_dense = multiomics_rna.X.toarray()
    else:
        multiomics_rna_dense = multiomics_rna.X

    # Apply KNN imputation
    knn_imputer = KNNImputer(n_neighbors=5)  # You can adjust the number of neighbors
    multiomics_rna_imputed = knn_imputer.fit_transform(multiomics_rna_dense)

    # Update the AnnData object with the imputed values
    multiomics_rna.X = multiomics_rna_imputed
    print('zero ration: ', (multiomics_rna.X==0).sum()/multiomics_rna.size)
print('Create corr net')
net = create_corr_net(multiomics_rna.X, groups, par['corr_method'])



print('Output GRN')
net.to_csv(par['prediction'])
