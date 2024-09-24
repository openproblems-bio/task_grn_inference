import pandas as pd 
import anndata as ad 
import numpy as np 
from tqdm import tqdm
from sklearn.preprocessing import StandardScaler
import scipy.sparse as sp

colors_blind = [
    '#E69F00',  # Orange
    '#56B4E9',  # Sky Blue
    '#009E73',  # Bluish Green
    '#F0E442',  # Yellow
    '#0072B2',  # Blue
    '#D55E00',  # Vermillion
    '#CC79A7']  # Reddish Purple


def verbose_print(verbose_level, message, level):
    if level <= verbose_level:
        print(message)
def verbose_tqdm(iterable, desc, level, verbose_level):  
    if level <= verbose_level:
        return tqdm(iterable, desc=desc)
    else:
        return iterable  # Return the iterable without a progress bar

def basic_qc(adata, min_genes_per_cell = 200, max_genes_per_cell = 5000, max_mt_frac = 0.05, min_cells_per_gene = 10):
    mt = adata.var_names.str.startswith('MT-')

    # 1. stats
    total_counts = adata.X.sum(axis=1)
    n_genes_by_counts = (adata.X > 0).sum(axis=1)
    mt_frac = adata[:, mt].X.sum(axis=1) / total_counts

    # 2. Filter cells
    mask_cells=  (n_genes_by_counts > min_genes_per_cell)& \
                 (n_genes_by_counts < max_genes_per_cell)&\
                 (mt_frac < max_mt_frac)
    # 3. Filter genes
    n_cells = (adata.X!=0).sum(axis=0)
    mask_genes = n_cells>min_cells_per_gene

    return adata[mask_cells, mask_genes]
def process_links(net, par):
    net = net[net.source!=net.target]
    net_sorted = net.reindex(net['weight'].abs().sort_values(ascending=False).index)
    net = net_sorted.head(par['max_n_links']).reset_index(drop=True)
    return net

def corr_net(X, gene_names, par, tf_all, causal=False):
    X = StandardScaler().fit_transform(X)
    net = np.dot(X.T, X) / X.shape[0]
    net = pd.DataFrame(net, index=gene_names, columns=gene_names)  
    if causal:  
        print('Causal subsetting')
        net = net[tf_all]
    else:
        print('Random subsetting')
        net = net.sample(len(tf_all), axis=1, random_state=par['seed'])
    net = net.reset_index()
    index_name = net.columns[0]
    net = net.melt(id_vars=index_name, var_name='source', value_name='weight')
    
    net.rename(columns={index_name: 'target'}, inplace=True)
    net = process_links(net, par)
    
    return net

def process_data(adata, par):
    if par['normalize']:
        import scanpy as sc
        print('Noramlize data')
        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)
        sc.pp.scale(adata)
    if sp.issparse(adata.X):
        # Convert to dense if it's sparse
        adata.X = adata.X.toarray()  # You can also use .todense(), but .toarray() gives a NumPy array directly
    else:
        print("adata.X is already dense.")
        
def create_corr_net(par):
    print(par)

    print('Read data')
    multiomics_rna = ad.read_h5ad(par["multiomics_rna"])
    process_data(multiomics_rna, par)
    gene_names = multiomics_rna.var_names.to_numpy()
    tf_all = np.loadtxt(par['tf_all'], dtype=str)
    groups = multiomics_rna.obs.cell_type
    tf_all = np.intersect1d(tf_all, gene_names)
    
    X = multiomics_rna.X
    if par['cell_type_specific']:
        print('cell_type_specific')
        i = 0
        for group in tqdm(np.unique(groups), desc="Processing groups"):
            X_sub = X[groups == group, :]
            net = corr_net(X_sub, gene_names, par, tf_all, par['causal'])
            net['cell_type'] = group
            if i==0:
                grn = net
            else:
                grn = pd.concat([grn, net], axis=0).reset_index(drop=True)
            i += 1
    else:
        grn = corr_net(X, gene_names, par, tf_all, par['causal'])    
    return grn
