import numpy as np
import scanpy as sc
import anndata as ad
import pandas as pd
import sys
from scipy.sparse import csr_matrix
import warnings
import gc
warnings.filterwarnings("ignore")


def sum_by(adata: ad.AnnData, col: str, unique_mapping: bool = True) -> ad.AnnData:
    """
    Adapted from this forum post:
    https://discourse.scverse.org/t/group-sum-rows-based-on-jobs-feature/371/4
    """

    assert pd.api.types.is_categorical_dtype(adata.obs[col])
    from scipy import sparse

    # sum `.X` entries for each unique value in `col`
    cat = adata.obs[col].values

    indicator = sparse.coo_matrix(
        (np.broadcast_to(True, adata.n_obs), (cat.codes, np.arange(adata.n_obs))),
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

    joining_df = (
        adata.obs[[col] + one_to_one_mapped_obs_cols].drop_duplicates().set_index(col)
    )

    assert (sum_adata.obs.index == sum_adata.obs.join(joining_df).index).all()
    sum_adata.obs = sum_adata.obs.join(joining_df)
    sum_adata.obs.index.name = col
    sum_adata.obs = sum_adata.obs.reset_index()
    sum_adata.obs.index = sum_adata.obs.index.astype("str")

    return sum_adata


def normalize_func(adata, log_norm=True, pearson_residual=True):
    import scanpy as sc

    # sc.pp.filter_cells(adata, min_genes=100)
    # sc.pp.filter_genes(adata, min_cells=10)

    if 'counts' not in adata.layers:
        adata.layers['counts'] = adata.X.copy()
    if pearson_residual:
        adata.layers['pearson_residual'] = sc.experimental.pp.normalize_pearson_residuals(adata)['X']
    if log_norm:
        X_norm = sc.pp.normalize_total(adata, layer='counts',inplace=False)['X']
        X_norm = sc.pp.log1p(X_norm, copy=False)
        adata.layers['lognorm'] = X_norm
    adata.X = adata.layers['counts']
    del adata.layers['counts']

    return adata

def pseudobulk_sum_func(adata, group):
    
    adata.obs[group] = adata.obs[group].astype('category')
    bulk_adata = sum_by(adata, group)
    bulk_adata.obs['cell_count'] = adata.obs.groupby(group).size().values
    print('ratio of missingness' , (bulk_adata.X==0).sum()/bulk_adata.X.size)

    bulk_adata.X = np.nan_to_num(bulk_adata.X, nan=0)
    return bulk_adata

def bulkify_main(adata, cell_count_t=10, covariates=['cell_type', 'donor_id', 'age']):
    adata.obs['sum_by'] = ''
    for covariate in covariates:
        adata.obs['sum_by'] += '_' + adata.obs[covariate].astype(str)
    adata.obs['sum_by'] = adata.obs['sum_by'].astype('category')
    adata_bulk = sum_by(adata, 'sum_by', unique_mapping=True)
    cell_count_df = adata.obs.groupby('sum_by').size().reset_index(name='cell_count')
    adata_bulk.obs = adata_bulk.obs.merge(cell_count_df, on='sum_by')
    adata_bulk = adata_bulk[adata_bulk.obs['cell_count']>=cell_count_t]
    adata_bulk.X = np.nan_to_num(adata_bulk.X, nan=0)
    return adata_bulk

def split_data_gene_perturbation(adata: ad.AnnData, train_share):
    tf_all = np.loadtxt('resources/grn_benchmark/prior/tf_all.csv', dtype=str)
    obs = adata.obs
    obs['is_tf'] = obs['perturbation'].isin(tf_all)
    
    unique_perts = obs['perturbation'].unique()
    tf_perturbs = obs[obs['is_tf']]['perturbation'].unique()
    non_tf_perturbs = np.setdiff1d(unique_perts, tf_perturbs)

    # calculate how many TFs and non-TFs to put in train
    n_train_tfs = int(train_share * len(tf_perturbs))
    n_train_non_tfs = int(train_share * len(non_tf_perturbs))

    control_perturbs = obs[obs['is_control']]['perturbation'].unique()

    # sample for train
    np.random.seed(32)
    train_tfs = np.random.choice(tf_perturbs, size=n_train_tfs, replace=False)
    train_non_tfs = np.random.choice(non_tf_perturbs, size=n_train_non_tfs, replace=False)

    # the rest go to test
    test_tfs = np.setdiff1d(tf_perturbs, train_tfs)
    test_non_tfs = np.setdiff1d(non_tf_perturbs, train_non_tfs)

    train_perturbs = np.concatenate([train_tfs, train_non_tfs])
    test_perturbs = np.concatenate([test_tfs, test_non_tfs])

    train_perturbs = np.concatenate([train_perturbs, control_perturbs])
    test_perturbs = np.concatenate([test_perturbs, control_perturbs])

    print(f"Train TFs: {len(train_tfs)}, Train non-TFs: {len(train_non_tfs)}")
    print(f"Test TFs: {len(test_tfs)}, Test non-TFs: {len(test_non_tfs)}")
    print(f"Train total: {len(train_perturbs)}, Test total: {len(test_perturbs)}")

    return train_perturbs.astype(str), test_perturbs.astype(str)
def wrapper_large_perturbation_data(adata, covariates, add_metadata, save_name, split_func, train_share=.5, cell_count_t=10):

    del adata.obsp 
    del adata.varm
    del adata.uns
    del adata.obsm
    del adata.var 
    del adata.uns

    adata.X = csr_matrix(adata.X) if not isinstance(adata.X, csr_matrix) else adata.X

    expected_cols = np.unique(['perturbation_type', 'perturbation', 'is_control'] + covariates)
    assert all(col in adata.obs.columns for col in expected_cols), "Expected columns not found in adata.obs"
    print('Processing large perturbation data...', flush=True)
    print('Pseudobulking data...', flush=True)
    print('Cell count threshold:', cell_count_t, flush=True)
    print('Covariates:', covariates, flush=True)

    print('QC on single cell number of perturbation ...', flush=True)
    single_cell_counts = adata.obs.groupby('perturbation').size().sort_values(ascending=False)
    passed_perturbations = single_cell_counts[single_cell_counts >= cell_count_t].index.tolist()
    adata = adata[adata.obs['perturbation'].isin(passed_perturbations)]

    print(f'Number of perturbations after QC: {len(passed_perturbations)}', flush=True)
    assert adata.shape[0] > 0, "No data found after QC"

    print('QC on number of genes and cells ...', flush=True)
    sc.pp.filter_genes(adata, min_cells=10)
    sc.pp.filter_cells(adata, min_genes=100)
    print(f'Shape of adata after QC: {adata.shape}', flush=True)

    print('Splitting data...', flush=True)
    train_perturbs, test_perturbs = split_func(adata, train_share)
    
    # - adata train sc full
    print('Processing adata train sc ...', flush=True)
    adata_train_sc = adata[adata.obs['perturbation'].isin(train_perturbs)].copy()
    assert adata_train_sc.obs['is_control'].sum()> 0, "No control data found in training set"
    assert adata_train_sc.shape[0] > 0, "No training data found after splitting"
    adata_train_sc = normalize_func(adata_train_sc, pearson_residual=False)
    adata_train_sc = add_metadata(adata_train_sc)
    adata_train_sc.write(f'resources/extended_data/{save_name}_train_sc.h5ad', compression='gzip')

    # - adata train sc subset
    print('Processing adata train sc subset (main inference data)...', flush=True)
    if save_name == 'parsebioscience':
        shortlisted_perts = train_perturbs[:10]
    else:
        tf_all = np.loadtxt('resources/grn_benchmark/prior/tf_all.csv', dtype=str)
        shortlisted_perts = np.intersect1d(tf_all, train_perturbs)
    
    adata_train_sc_subset = adata_train_sc[adata_train_sc.obs['perturbation'].isin(shortlisted_perts), :]
    adata_train_sc_subset.write(f'resources/grn_benchmark/inference_data/{save_name}_rna.h5ad', compression='gzip')
    

    del adata_train_sc
    gc.collect()

    print('Processing adata test sc ...', flush=True)
    adata_test_sc = adata[adata.obs['perturbation'].isin(test_perturbs)].copy()
    assert adata_test_sc.obs['is_control'].sum()> 0, "No control data found in test set"
    adata_test_sc = normalize_func(adata_test_sc, pearson_residual=False)
    assert adata_test_sc.shape[0] > 0, "No test data found after splitting"
    adata_test_sc = add_metadata(adata_test_sc)
    adata_test_sc.write( f'resources/grn_benchmark/evaluation_data/{save_name}_sc.h5ad', compression='gzip')
    del adata_test_sc
    gc.collect()

    print('Processing adata sc ...', flush=True)
    adata = normalize_func(adata, pearson_residual=False)
    adata = add_metadata(adata)
    adata.write(f'resources/processed_data/{save_name}_sc.h5ad', compression='gzip')
    

    print('Process adata bulk ...', flush=True)
    adata_bulk = bulkify_main(adata, cell_count_t=cell_count_t, covariates=covariates)
    adata_bulk = adata_bulk.copy()
    del adata
    gc.collect()
    
    assert adata_bulk.shape[0] > 0, "No bulk data found after pseudobulking"
    assert 'perturbation' in adata_bulk.obs.columns, "Expected 'perturbation' column not found in adata_bulk.obs"
    adata_bulk = normalize_func(adata_bulk, pearson_residual=False)
    adata_bulk = add_metadata(adata_bulk)
    adata_bulk.write(f'resources/extended_data/{save_name}_bulk.h5ad', compression='gzip')

    print('Process adata bulk test...', flush=True)
    adata_test_bulk = adata_bulk[adata_bulk.obs['perturbation'].isin(test_perturbs)].copy()
    assert adata_test_bulk.shape[0] > 0, "No test data found after splitting"
    adata_test_bulk = normalize_func(adata_test_bulk, pearson_residual=False)
    adata_test_bulk = add_metadata(adata_test_bulk)
    adata_test_bulk.write(f'resources/grn_benchmark/evaluation_data/{save_name}_bulk.h5ad', compression='gzip')
    
    print('Process adata bulk train...', flush=True)
    adata_train_bulk = adata_bulk[adata_bulk.obs['perturbation'].isin(train_perturbs)].copy()
    assert adata_train_bulk.shape[0] > 0, "No training bulk data found after splitting"
    adata_train_bulk = normalize_func(adata_train_bulk, pearson_residual=False)
    adata_train_bulk = add_metadata(adata_train_bulk)
    adata_train_bulk.write(f'resources/extended_data/{save_name}_train_bulk.h5ad', compression='gzip')
    
    