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

def basic_qc(adata, min_genes_per_cell = 200, max_genes_per_cell = 5000, min_cells_per_gene = 10):
    mt = adata.var_names.str.startswith('MT-')
    print('shape before ', adata.shape)
    # 1. stats
    total_counts = adata.X.sum(axis=1)
    n_genes_by_counts = (adata.X > 0).sum(axis=1)
    # mt_frac = adata[:, mt].X.sum(axis=1) / total_counts
    
    low_gene_filter = (n_genes_by_counts < min_genes_per_cell)
    high_gene_filter = (n_genes_by_counts > max_genes_per_cell)
    # mt_filter = mt_frac > max_mt_frac

    # 2. Filter cells
    # print(f'Number of cells removed: below min gene {low_gene_filter.sum()}, exceed max gene {high_gene_filter.sum()}')
    mask_cells=  (~low_gene_filter)& \
                 (~high_gene_filter)
                #  (~mt_filter)
    # 3. Filter genes
    n_cells = (adata.X!=0).sum(axis=0)
    mask_genes = n_cells>min_cells_per_gene
    adata_f = adata[mask_cells, mask_genes]
    print('shape after ', adata_f.shape)
    return adata_f

def process_links(net, par):
    net = net[net.source!=net.target]
    
    if par['max_n_links'] != -1:
        net_sorted = net.reindex(net['weight'].abs().sort_values(ascending=False).index)
        net = net_sorted.head(par['max_n_links']).reset_index(drop=True)
    return net
def efficient_melting(net, gene_names, par):
    '''to replace pandas melting'''
    upper_triangle_indices = np.triu_indices_from(net, k=1)

    # Extract the source and target gene names based on the indices
    sources = np.array(gene_names)[upper_triangle_indices[0]]
    targets = np.array(gene_names)[upper_triangle_indices[1]]

    # Extract the corresponding correlation values
    weights = net[upper_triangle_indices]
    # Create a structured array
    data = np.column_stack((targets, sources, weights))

    # Convert to DataFrame
    print('conver to df')
    net = pd.DataFrame(data, columns=['target', 'source', 'weight'])
    return net
    

def corr_net(X, gene_names, par):
    print('calculate correlation')
    if hasattr(X, 'todense'):
        net = np.corrcoef(X.todense().T)
    else:
        net = np.corrcoef(X.T)
    print('process corr results')
    
    tf_all = np.loadtxt(par['tf_all'], dtype=str)
    tf_all = np.intersect1d(tf_all, gene_names)

    # # Convert to a DataFrame with gene names as both row and column indices
    net = pd.DataFrame(net, index=gene_names, columns=gene_names)
    
    net = net.values

    net = efficient_melting(net, gene_names, par)
    if par['causal']:  
        print('TF subsetting')
        net = net[net.source.isin(tf_all)]
    net = process_links(net, par)
    print('Corr results are ready')
    return net


def read_gmt(file_path:str) -> dict[str, list[str]]:
    '''Reas gmt file and returns a dict of gene'''
    gene_sets = {}
    with open(file_path, 'r') as file:
        for line in file:
            parts = line.strip().split('\t')
            gene_set_name = parts[0]
            gene_set_description = parts[1]
            genes = parts[2:]
            gene_sets[gene_set_name] = {
                'description': gene_set_description,
                'genes': genes
            }
    return gene_sets
def quantile_transformation(values, one_sided=False, log1p_scale=True):
    from sklearn.preprocessing import QuantileTransformer
    if log1p_scale:
        log_data = np.log1p(values)  # log(x + 1) to avoid log(0)
    if one_sided:
        output_distribution = 'uniform'
    else:
        output_distribution = 'normal'
    quantile_transformer = QuantileTransformer(output_distribution=output_distribution)
    transformed_data = quantile_transformer.fit_transform(log_data.reshape(-1, 1)).reshape(len(log_data))
    return transformed_data
def zscore_transformation(values, one_sided=False, log1p_scale=True):
    if log1p_scale:
        log_data = np.log1p(values)  # log(x + 1) to avoid log(0)
    if one_sided:
        mean = 0
    else:
        mean = np.mean(values)
    std = np.std(values)
    transformed_data = (log_data-mean)/std
    return transformed_data
