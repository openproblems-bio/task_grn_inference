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

def normalize_func(adata, log_norm=True, pearson_residual=False, target_sum=1e4):
    import scanpy as sc
    import scipy.sparse as sp
    if 'counts' not in adata.layers:
        adata.layers['counts'] = adata.X.copy()
    if pearson_residual:
        adata.layers['pearson_residual'] = sc.experimental.pp.normalize_pearson_residuals(adata)['X']
    if log_norm:
        X_norm = sc.pp.normalize_total(adata, layer='counts', inplace=False, target_sum=target_sum)['X']
        X_norm = sc.pp.log1p(X_norm, copy=False)
        adata.layers['lognorm'] = sp.csr_matrix(X_norm)
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

def split_data_gene_perturbation(adata: ad.AnnData, train_share=.5):
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

    train_group = np.concatenate([train_tfs, train_non_tfs])
    test_group = np.concatenate([test_tfs, test_non_tfs])

    train_group = np.concatenate([train_group, control_perturbs])
    test_group = np.concatenate([test_group, control_perturbs])

    print(f"Train TFs: {len(train_tfs)}, Train non-TFs: {len(train_non_tfs)}")
    print(f"Test TFs: {len(test_tfs)}, Test non-TFs: {len(test_non_tfs)}")
    print(f"Train total: {len(train_group)}, Test total: {len(test_group)}")

    return train_group.astype(str), test_group.astype(str)

import numpy as np
from scipy.stats import ttest_ind
from tqdm import tqdm

def compare_perturbed_to_control_sc(adata, gene, layer):
    from scipy.sparse import issparse
    adata.X = adata.layers[layer]
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

        # Guard 1: check for nans / infs
    if (np.any(~np.isfinite(control_exprs)) or np.any(~np.isfinite(perturbed_exprs))):
        raise ValueError(f"Non-finite values found in expressions for {gene}.")

    # Guard 2: too few samples
    if len(control_exprs) < 3 or len(perturbed_exprs) < 3:
        print(f"Not enough samples to perform statistical test for {gene}. ")
        return {
            "logf2c": 0,
            "mean_control": 0,
            "mean_perturbed": 0,
            "p_val": 1

        }
    # Guard 3: zero variance in both groups
    elif np.std(control_exprs) < eps and np.std(perturbed_exprs) < eps:
        return {
            "logf2c": 0,
            "mean_control": 0,
            "mean_perturbed": 0,
            "p_val": 1

        }
    else:
        t_stat, p_val = ttest_ind(perturbed_exprs, control_exprs, equal_var=False)
        if not np.isfinite(p_val):
            p_val = 1.0
    
    # Compute log2 fold change of means
    mean_control = np.mean(control_exprs)
    mean_perturbed = np.mean(perturbed_exprs)
    logf2c = np.log2((mean_perturbed + eps) / (mean_control + eps))
    
    # Statistical significance: t-test
    if (np.all(np.abs(control_exprs) < eps))& (np.all(np.abs(perturbed_exprs) < eps)):
        p_val = 1.0  # No variation, no significance
    else:
        t_stat, p_val = ttest_ind(perturbed_exprs, control_exprs, equal_var=False)
    if np.isnan(p_val):
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


def qc_perturbation(
    adata, 
    col='perturbation', 
    min_cells_per_pert=10, 
    min_cells_per_gene=10, 
    min_genes_per_cell=200
):
    """
    Filter AnnData object:
      1. Remove perturbations with too few cells.
      2. Remove genes with too few cells per condition.
      3. Remove cells with too few genes.
    """
    import scipy.sparse as sp
    # 1. Remove perturbations with too few cells
    pert_counts = adata.obs[col].value_counts()
    keep_perts = pert_counts[pert_counts >= min_cells_per_pert].index
    adata = adata[adata.obs[col].isin(keep_perts)]
    print(f"Kept {len(keep_perts)} perturbations (>= {min_cells_per_pert} cells each) out of {len(pert_counts)} total")

    # 2. Remove genes with too few cells in total (not per condition)
    if sp.issparse(adata.X):
        cell_counts = np.array((adata.X > 0).sum(axis=0)).flatten()
    else:
        cell_counts = (adata.X > 0).sum(axis=0)

    gene_mask = cell_counts >= min_cells_per_gene
    adata = adata[:, gene_mask]
    print(f"Kept {gene_mask.sum()} genes (>= {min_cells_per_gene} cells total) out of {len(gene_mask)} total")
    # 3. Remove cells with too few expressed genes
    if sp.issparse(adata.X):
        gene_counts = np.array((adata.X > 0).sum(axis=1)).flatten()
    else:
        gene_counts = (adata.X > 0).sum(axis=1)
    keep_cells = gene_counts >= min_genes_per_cell
    adata = adata[keep_cells]
    print(f"Kept {keep_cells.sum()} cells (>= {min_genes_per_cell} genes per cell) out of {len(keep_cells)} total")

    return adata


def qc_perturbation_effect_func(adata, n_jobs=-1):
    """
    Perform QC on perturbation effects with parallel processing.
    
    Parameters
    ----------
    adata : AnnData
        AnnData object with 'perturbation', 'perturbation_type', 'is_control' in .obs
    n_jobs : int
        Number of parallel jobs (-1 = use all cores)
        
    Returns
    -------
    adata : AnnData
        Subsetted AnnData with only QC-passed perturbations and control cells
    """
    from joblib import Parallel, delayed
    import pandas as pd
    from statsmodels.stats.multitest import multipletests
    from tqdm import tqdm
    p_val_t = 0.05

    perturbation_type = adata.obs['perturbation_type'].unique()[0]
    if perturbation_type in ['knockout', 'knockdown']:
        logf2c_t = -0.54
    elif perturbation_type in ['overexpression', 'activation']:
        raise ValueError(f"{perturbation_type} perturbation type is not supported for this analysis.")
    elif perturbation_type == 'cytokine':
        print('Cytokine perturbation type detected, no filtering for QC perturbation effect.')
        return adata

    # Determine which layer to use
    if 'lognorm' in adata.layers:
        layer = 'lognorm'
    elif 'X_norm' in adata.layers:
        layer = 'X_norm'
    else:
        raise ValueError("No suitable normalized layer found in adata.layers. Expected 'lognorm' or 'X_norm'.")

    # Precompute masks and sets
    perturbations = adata.obs['perturbation'].unique()
    var_names_set = set(adata.var_names)
    is_control_mask = adata.obs['is_control'].values

    # Function to process a single perturbation
    def process_perturbation(pert):
        if pert == 'control' or pert not in var_names_set:
            return (pert, None)
        pert_mask = adata.obs['perturbation'].values == pert
        subset_mask = pert_mask | is_control_mask
        rr = compare_perturbed_to_control_sc(
            adata[subset_mask].copy(),
            gene=pert,
            layer=layer
        )
        return (pert, rr)

    # Parallel computation
    results_list = Parallel(n_jobs=n_jobs)(
        delayed(process_perturbation)(p) for p in tqdm(perturbations, desc='Processing perturbations')
    )

    # Collect results into a dict
    results = {pert: rr for pert, rr in results_list if rr is not None}

    # Create DataFrame
    df_results = pd.DataFrame.from_dict(results, orient="index").reset_index()
    df_results.rename(columns={'index': 'perturbation'}, inplace=True)

    # Adjust p-values
    reject, pvals_corrected, _, _ = multipletests(df_results["p_val"].values, method="fdr_bh")
    df_results['p_adj'] = pvals_corrected

    # Filter based on log fold change and significance
    df_filtered = df_results[(df_results['logf2c'] < logf2c_t) & (df_results['p_adj'] < p_val_t)]
    perturbations_qc_passed = df_filtered['perturbation'].tolist()
    print(f'Number of perturbations after QC: {len(perturbations_qc_passed)}', flush=True)

    if len(perturbations_qc_passed) == 0:
        raise ValueError("No perturbations passed the QC criteria. Please adjust thresholds or check data.")

    # Subset adata to only include passed perturbations and controls
    adata = adata[(adata.obs['perturbation'].isin(perturbations_qc_passed)) | adata.obs['is_control']]

    return adata


def wrapper_large_perturbation_data(adata, covariates, add_metadata, save_name, split_func, 
                                    qc_group='perturbation', 
                                    qc_perturbation_effect=True,
                                    group='perturbation',
                                    save_data=True,
                                    cell_count_t=10):

    del adata.obsp 
    del adata.varm
    del adata.uns
    del adata.obsm
    del adata.var 
    del adata.uns
    if 'counts' in adata.layers:
        adata.X = adata.layers['counts']
    adata.X = csr_matrix(adata.X) if not isinstance(adata.X, csr_matrix) else adata.X
    adata = normalize_func(adata, pearson_residual=False)
    if "_index" in adata.var.columns:
        adata.var = adata.var.drop(columns=["_index"])

    del adata.raw

    expected_cols = np.unique(['perturbation_type', 'perturbation', 'is_control'] + covariates)
    assert all(col in adata.obs.columns for col in expected_cols), "Expected columns not found in adata.obs"
    print('Processing large perturbation data...', flush=True)
    print('Pseudobulking data...', flush=True)
    print('Covariates:', covariates, flush=True)

    if qc_perturbation_effect:
        print('QC on perturbation effect ...', flush=True)
        adata = qc_perturbation_effect_func(adata)

    print('QC on single cell number of perturbation ...', flush=True)
    adata = qc_perturbation(adata, col=qc_group, min_cells_per_pert=cell_count_t)

    print(f'Shape of adata after qc: {adata.shape}', flush=True)
    assert adata.shape[0] > 0, "No data found after QC"

    print('Splitting data...', flush=True)
    train_group, test_group = split_func(adata)
    
    # - adata train sc full
    print('Processing adata train sc ...', flush=True)
    adata_train_sc = adata[adata.obs[group].isin(train_group)].copy()
    assert adata_train_sc.obs['is_control'].sum()> 0, "No control data found in training set"
    assert adata_train_sc.shape[0] > 0, "No training data found after splitting"
    adata_train_sc = add_metadata(adata_train_sc)
    print('Shape of adata_train_sc: ', adata_train_sc.shape, flush=True)
    if save_data:
        adata_train_sc.write(f'resources/extended_data/{save_name}_train_sc.h5ad', compression='gzip')
    del adata_train_sc
    gc.collect()

    print('Processing adata test sc ...', flush=True)
    adata_test_sc = adata[adata.obs[group].isin(test_group)].copy()
    print('Shape of adata_test_sc: ', adata_test_sc.shape, flush=True)
    assert adata_test_sc.obs['is_control'].sum()> 0, "No control data found in test set"
    assert adata_test_sc.shape[0] > 0, "No test data found after splitting"
    adata_test_sc = add_metadata(adata_test_sc)
    if save_data:
        adata_test_sc.write( f'resources/processed_data/{save_name}_evaluation_sc.h5ad', compression='gzip')
    del adata_test_sc
    gc.collect()

    print('Processing adata sc ...', flush=True)
    adata = add_metadata(adata)
    if save_data:
        adata.write(f'resources/processed_data/{save_name}_sc.h5ad', compression='gzip')
    

    print('Processing adata bulk ...', flush=True)
    adata_bulk = bulkify_main(adata, covariates=covariates)
    adata_bulk = adata_bulk.copy()
    print('Shape of adata_bulk: ', adata_bulk.shape, flush=True)
    del adata
    gc.collect()
    
    assert adata_bulk.shape[0] > 0, "No bulk data found after pseudobulking"
    assert group in adata_bulk.obs.columns, f"Expected {group} column not found in adata_bulk.obs"
    adata_bulk = normalize_func(adata_bulk, pearson_residual=False)
    adata_bulk = add_metadata(adata_bulk)
    if save_data:
        adata_bulk.write(f'resources/extended_data/{save_name}_bulk.h5ad', compression='gzip')

    print('Process adata bulk test...', flush=True)
    adata_test_bulk = adata_bulk[adata_bulk.obs[group].isin(test_group)].copy()
    assert adata_test_bulk.shape[0] > 0, "No test data found after splitting"
    print('Shape of adata_test_bulk: ', adata_test_bulk.shape, flush=True)
    adata_test_bulk = add_metadata(adata_test_bulk)
    if save_data:
        adata_test_bulk.write(f'resources/grn_benchmark/evaluation_data/{save_name}_bulk.h5ad', compression='gzip')
    
    print('Process adata bulk train...', flush=True)
    adata_train_bulk = adata_bulk[adata_bulk.obs[group].isin(train_group)].copy()
    assert adata_train_bulk.shape[0] > 0, "No training bulk data found after splitting"
    print('Shape of adata_train_bulk: ', adata_train_bulk.shape, flush=True)
    adata_train_bulk = add_metadata(adata_train_bulk)
    if save_data:
        adata_train_bulk.write(f'resources/grn_benchmark/inference_data/{save_name}_rna.h5ad', compression='gzip')
    
    print('Data processing completed successfully!', flush=True)


def bulkify_func(adata, cell_count_t=10, covariates=['cell_type', 'donor_id', 'age']):
    from task_grn_inference.src.process_data.helper_data import sum_by
    adata.obs['sum_by'] = ''
    for covariate in covariates:
        adata.obs['sum_by'] += '_' + adata.obs[covariate].astype(str)
    # adata.obs['sum_by'] = '_' + adata.obs['cell_type'].astype(str) + '_' + adata.obs['donor_id'].astype(str) + '_' + adata.obs['age'].astype(str) 
    adata.obs['sum_by'] = adata.obs['sum_by'].astype('category')
    adata_bulk = sum_by(adata, 'sum_by', unique_mapping=True)
    cell_count_df = adata.obs.groupby('sum_by').size().reset_index(name='cell_count')
    adata_bulk.obs = adata_bulk.obs.merge(cell_count_df, on='sum_by')
    adata_bulk = adata_bulk[adata_bulk.obs['cell_count']>=cell_count_t]
    return adata_bulk

