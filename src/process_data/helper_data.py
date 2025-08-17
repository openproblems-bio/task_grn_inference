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

import numpy as np
from scipy.stats import ttest_ind
from tqdm import tqdm

def compare_perturbed_to_control_sc(adata, gene):
    from scipy.sparse import issparse
    if gene not in adata.var_names:
        raise ValueError(f"Gene '{gene}' not found in the dataset.")
    
    mask_gene = adata.var_names == gene
    control = adata[adata.obs['is_control']]
    perturbed = adata[~adata.obs['is_control']]
    assert len(control) > 0, "No control samples found in the data."
    assert len(perturbed) > 0, "No perturbed samples found in the data."
    
    # Assuming control[:, mask_gene].X is sparse
    control_exprs = control[:, mask_gene].X
    if hasattr(control_exprs, "toarray"):
        control_exprs = control_exprs.toarray().ravel()
    else:
        control_exprs = np.ravel(control_exprs)

    perturbed_exprs = perturbed[:, mask_gene].X
    if hasattr(perturbed_exprs, "toarray"):
        perturbed_exprs = perturbed_exprs.toarray().ravel()
    else:
        perturbed_exprs = np.ravel(perturbed_exprs)

    eps = 1e-8
    if np.all(np.abs(control_exprs) < eps):
        print(f"All control expressions are zero for {gene}.")
    
    # Compute log2 fold change of means
    mean_control = np.mean(control_exprs)
    mean_perturbed = np.mean(perturbed_exprs)
    logf2c = np.log2((mean_perturbed + eps) / (mean_control + eps))
    
    # Statistical significance: t-test
    if (np.all(np.abs(control_exprs) < eps))& (np.all(np.abs(perturbed_exprs) < eps)):
        p_val = 1.0  # No variation, no significance
    else:
        t_stat, p_val = ttest_ind(perturbed_exprs, control_exprs, equal_var=False)
    if p_val < 0:
        print("P value is not valid: ", p_val)
        print(control_exprs)
        print(perturbed_exprs)
        raise ValueError("P value is not valid, check the data.")
    return {
        "logf2c": logf2c,
        "mean_control": float(mean_control),
        "mean_perturbed": float(mean_perturbed),
        "p_val": float(p_val)
    }
def qc_perturbation_effect(adata):
    from tqdm import tqdm
    p_val_t=0.05

    perturbation_type = adata.obs['perturbation_type'].unique()[0]
    if perturbation_type=='knockout':
        logf2c_t = -0.54
    elif perturbation_type=='knockdown':
        logf2c_t = -0.54
    elif perturbation_type=='overexpression':
        raise ValueError("Overexpression perturbation type is not supported for this analysis.")
    elif perturbation_type=='activation':
        raise ValueError("Activation perturbation type is not supported for this analysis.")
    elif perturbation_type=='cytokine':
        print('Cytokine perturbation type detected, no filtering for qc pertubration effect.')
        return adata

    # Collect results
    results = {}
    if 'lognorm' in adata.layers:
        adata.X = adata.layers['lognorm']
    elif 'X_norm' in adata.layers:
        adata.X = adata.layers['X_norm']
    else:
        raise ValueError("No suitable normalized layer found in adata.layers. Expected 'lognorm' or 'X_norm'.")

    for perturbation in tqdm(adata.obs['perturbation'].unique(), desc='Processing perturbations'):
        if perturbation not in adata.var_names:
            continue
        if perturbation != 'control':
            perturbed_gene = perturbation
            rr = compare_perturbed_to_control_sc(
                adata[(adata.obs['perturbation'] == perturbation) | adata.obs['is_control']], 
                perturbed_gene
            )
            results[perturbation] = rr

    # Create DataFrame
    df_results = pd.DataFrame.from_dict(results, orient="index")
    
    df_results.index.name = "perturbation"
    df_results = df_results.reset_index()

    # - Filter results based on significance and fold change
    from statsmodels.stats.multitest import multipletests
    reject, pvals_corrected, _, _ = multipletests(df_results["p_val"].values, method="fdr_bh")
    df_results["p_adj"] = pvals_corrected
    print(df_results)
    aaa
    df_filtered = df_results[(df_results['logf2c']<logf2c_t)&(df_results['p_adj']<p_val_t)]
    perturbations_qc_passed = df_filtered['perturbation'].tolist()

    print(f'Number of perturbations after QC: {len(perturbations_qc_passed)}', flush=True)

    if len(perturbations_qc_passed) == 0:
        raise ValueError("No perturbations passed the QC criteria. Please adjust the thresholds or check the data.")
    # - Update adata with passed perturbations
    adata = adata[(adata.obs['perturbation'].isin(perturbations_qc_passed)) | adata.obs['is_control']]

    return adata


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

    print('QC on perturbation effect ...', flush=True)
    adata = normalize_func(adata, pearson_residual=False)
    adata = qc_perturbation_effect(adata)

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
    # adata_train_sc = normalize_func(adata_train_sc, pearson_residual=False)
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
    # adata_test_sc = normalize_func(adata_test_sc, pearson_residual=False)
    assert adata_test_sc.shape[0] > 0, "No test data found after splitting"
    adata_test_sc = add_metadata(adata_test_sc)
    adata_test_sc.write( f'resources/processed_data/{save_name}_evaluation_sc.h5ad', compression='gzip')
    del adata_test_sc
    gc.collect()

    print('Processing adata sc ...', flush=True)
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
    # adata_test_bulk = normalize_func(adata_test_bulk, pearson_residual=False)
    adata_test_bulk = add_metadata(adata_test_bulk)
    adata_test_bulk.write(f'resources/grn_benchmark/evaluation_data/{save_name}_bulk.h5ad', compression='gzip')
    
    print('Process adata bulk train...', flush=True)
    adata_train_bulk = adata_bulk[adata_bulk.obs['perturbation'].isin(train_perturbs)].copy()
    assert adata_train_bulk.shape[0] > 0, "No training bulk data found after splitting"
    # adata_train_bulk = normalize_func(adata_train_bulk, pearson_residual=False)
    adata_train_bulk = add_metadata(adata_train_bulk)
    adata_train_bulk.write(f'resources/extended_data/{save_name}_train_bulk.h5ad', compression='gzip')
    
    print('Data processing completed successfully!', flush=True)
    