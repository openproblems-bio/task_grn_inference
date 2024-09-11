import os
import pandas as pd
import numpy as np
import anndata as ad
import scanpy as sc
from tqdm import tqdm
from scipy.stats import spearmanr

## VIASH START
par = {
}
## VIASH END
def create_corr_net(X: np.ndarray, groups: np.ndarray, method="pearson"):
    grns = []
    for group in tqdm(np.unique(groups), desc="Processing groups"):
        X_sub = X[groups == group, :]
        if method is "pearson":
            grn = np.corrcoef(X_sub.T)
        elif method is "spearman":
            grn = spearmanr(X_sub).statistic
        grns.append(grn)
    return np.mean(grns, axis=0)
print('Read data')
multiomics_rna = ad.read_h5ad(par["multiomics_rna"])
# print('subsetting: remove this')
# multiomics_rna = multiomics_rna[:5000, :5000]
gene_names = multiomics_rna.var_names.to_numpy()
tf_all = np.loadtxt(par['tf_all'], dtype=str)
groups = multiomics_rna.obs.cell_type
tf_all = np.intersect1d(tf_all, gene_names)

print('Noramlize data')
sc.pp.normalize_total(multiomics_rna)
sc.pp.log1p(multiomics_rna)
sc.pp.scale(multiomics_rna)

print('Create corr net')
net = create_corr_net(multiomics_rna.X, groups)
net = pd.DataFrame(net, index=gene_names, columns=gene_names)

if par['causal']:
    net = net[tf_all]
else:
    net = net.sample(len(tf_all), axis=1, random_state=par['seed'])
    
net = net.reset_index().melt(id_vars='index', var_name='source', value_name='weight')
net.rename(columns={'index': 'target'}, inplace=True)


print('Output GRN')
net.to_csv(par['prediction'])
