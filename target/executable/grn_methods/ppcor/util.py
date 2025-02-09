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

def binarize_weight(weight):
    if weight > 0:
        return 1
    elif weight < 0:
        return -1
    else:
        return 0
def verbose_print(verbose_level, message, level):
    if level <= verbose_level:
        print(message)
def verbose_tqdm(iterable, desc, level, verbose_level):  
    if level <= verbose_level:
        return tqdm(iterable, desc=desc)
    else:
        return iterable  # Return the iterable without a progress bar
def efficient_melting(net, gene_names):
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
    # print('convert to df')
    net = pd.DataFrame(data, columns=['source', 'target', 'weight'])
    return net

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

def corr_net(par: dict) -> pd.DataFrame:
    # - read data
    adata = ad.read_h5ad(par["rna"])
    tf_all = np.loadtxt(par['tf_all'], dtype=str)
    if 'X_norm' in adata.layers.keys():
        X = adata.layers['X_norm']
    else:
        X = adata.X
    # - remove genes with 0 standard deviation
    gene_std = np.std(X, axis=0)
    nonzero_std_genes = gene_std > 0
    X = X[:, nonzero_std_genes]
    gene_names = adata[:, nonzero_std_genes].var_names.to_numpy()
    # - calculate correlation
    if hasattr(X, 'todense'):
        net = np.corrcoef(X.todense().T)
    else:
        net = np.corrcoef(X.T)    
    # - melt the matrix
    net = pd.DataFrame(net, index=gene_names, columns=gene_names).values
    net = efficient_melting(net, gene_names)
    # - subset to known TFs
    tf_all = np.intersect1d(tf_all, gene_names)
    net = net[net.source.isin(tf_all)]
    # - process links: size control
    net = process_links(net, par)
  
    return net


def plot_heatmap(scores, ax=None, name='', fmt='0.02f', cmap="viridis"):
    import matplotlib.pyplot as plt
    import numpy as np
    import seaborn

    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(4, 4), sharey=True)

    # Normalize each column individually
    scores_normalized = scores.apply(lambda x: (x - np.nanmin(x)) / (np.nanmax(x) - np.nanmin(x)), axis=0)
    scores_normalized = scores_normalized.round(2)
    # scores_normalized['Rank'] = scores['Rank'].max()-scores['Rank']
    # scores_normalized['Rank'] = scores_normalized['Rank']/scores_normalized['Rank'].max()

    # Define the color scale range for each column (0 to 1 after normalization)
    vmin = 0
    vmax = 1



    # Plot the heatmap with normalized scores
    seaborn.heatmap(scores_normalized, ax=ax, square=False, cbar=False, annot=True, fmt=fmt, vmin=vmin, vmax=vmax, cmap=cmap)
    # Overlay the original (unnormalized) scores as annotations
    # scores['Rank'] = scores['Rank'].astype(int)
    # print(scores['Rank'])
    # Overlay the original (unnormalized) scores as annotations
    for text, (i, j) in zip(ax.texts, np.ndindex(scores.shape)):
        value = scores.iloc[i, j]
        if isinstance(value, np.int64):  # Check if the value is an integer for 'Rank'
            text.set_text(f'{value:d}')
        else:
            text.set_text(f'{value:.2f}')

    # Customize the axes and title
    ax.tick_params(left=False, bottom=False)
    ax.xaxis.set_tick_params(width=0)
    ax.yaxis.set_tick_params(width=0)
    ax.set_title(name, pad=10)

    ax.xaxis.set_label_position('top')
    ax.xaxis.tick_top()
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='left')
    

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
def sum_by(adata: ad.AnnData, col: str, unique_mapping: bool=True) -> ad.AnnData:
    """
    Adapted from this forum post: 
    https://discourse.scverse.org/t/group-sum-rows-based-on-jobs-feature/371/4
    """
    
    assert pd.api.types.is_categorical_dtype(adata.obs[col])
    from scipy import sparse
    

    # sum `.X` entries for each unique value in `col`
    cat = adata.obs[col].values

    indicator = sparse.coo_matrix(
        (
            np.broadcast_to(True, adata.n_obs),
            (cat.codes, np.arange(adata.n_obs))
        ),
        shape=(len(cat.categories), adata.n_obs),
    )
  
    sum_adata = ad.AnnData(
        indicator @ adata.X,
        var=adata.var,
        obs=pd.DataFrame(index=cat.categories),
    )
    
    # copy over `.obs` values that have a one-to-one-mapping with `.obs[col]`
    obs_cols = adata.obs.columns
    obs_cols = list(set(adata.obs.columns) - set([col]))
    
    if unique_mapping:
        one_to_one_mapped_obs_cols = []
        nunique_in_col = adata.obs[col].nunique()
        for other_col in obs_cols:
            if len(adata.obs[[col, other_col]].drop_duplicates()) == nunique_in_col:
                one_to_one_mapped_obs_cols.append(other_col)
    else:
        one_to_one_mapped_obs_cols = obs_cols

    joining_df = adata.obs[[col] + one_to_one_mapped_obs_cols].drop_duplicates().set_index(col)
    assert (sum_adata.obs.index == sum_adata.obs.join(joining_df).index).all()
    sum_adata.obs = sum_adata.obs.join(joining_df)
    sum_adata.obs.index.name = col
    sum_adata.obs = sum_adata.obs.reset_index()
    sum_adata.obs.index = sum_adata.obs.index.astype('str')

    return sum_adata


def fetch_gene_info():

    from pybiomart import Server
    # Connect to Ensembl server
    server = Server(host='http://www.ensembl.org')

    # Select the dataset for human genes
    dataset = server.marts['ENSEMBL_MART_ENSEMBL'].datasets['hsapiens_gene_ensembl']

    # Query relevant attributes
    df = dataset.query(
        attributes=[
            'ensembl_gene_id',        # Ensembl Gene ID
            'external_gene_name',     # Gene Name
            'chromosome_name',        # Chromosome
            'start_position',         # Start site
            'end_position'            # End site
        ]
    )
    df.columns = ['gene_id', 'gene_name', 'chr', 'start', 'end']
    # df = df[df['chr'].isin(['X', 'Y']+list(map(str, range(1, 23))))]
    # - keep those genes with one mapping
    df = df.groupby('gene_name').filter(lambda x: len(x) == 1)
    df.set_index('gene_name', inplace=True)

    return df