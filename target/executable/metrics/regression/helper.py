from typing import Dict, List, Set, Tuple, Any, Union
import random
import json

import tqdm
import numpy as np
import lightgbm
import pandas as pd
import anndata as ad
from sklearn.preprocessing import LabelEncoder, RobustScaler, StandardScaler
from sklearn.model_selection import GroupKFold, LeaveOneGroupOut
from sklearn.linear_model import Ridge
from sklearn.metrics import r2_score, mean_squared_error

from util import verbose_print, verbose_tqdm, read_prediction, manage_layer


SEED = 0xCAFE
N_POINTS_TO_ESTIMATE_BACKGROUND = 20

def net_to_matrix(net, gene_names: np.ndarray) -> np.ndarray:
    gene_dict = {gene_name: i for i, gene_name in enumerate(gene_names)}
    # convert to matrix
    A = np.zeros((len(gene_names), len(gene_names)), dtype=float)
    for source, target, weight in zip(net['source'], net['target'], net['weight']):
        if (source not in gene_dict) or (target not in gene_dict):
            continue
        i = gene_dict[source]
        j = gene_dict[target]
        A[i, j] = float(weight)
    return A


def fill_zeros_in_grn(A: np.ndarray, eps: float = 1e-10) -> np.ndarray:
    A = np.copy(A)
    A[A > 0] = A[A > 0] + eps
    A[A < 0] = A[A < 0] - eps
    A[A == 0] = np.random.rand(*A[A == 0].shape) * 2 * eps - eps
    A += np.random.rand(*A.shape) * eps
    return A


def cross_validate_gene(
        reg_type: str,
        X: np.ndarray,
        groups: np.ndarray,
        grn: np.ndarray,
        j: int,
        n_features: int = 10,
        random_state: int = 0xCAFE,
        n_jobs: int = 10
) -> Dict[str, float]:
    
    results = {'r2': 0.0, 'avg-r2': 0.0}
    
    if n_features == 0:
        return results
    
    # Feature selection
    regulatory_importance = np.abs(grn[:, j])
    regulatory_importance[j] = -1
    selected_features = np.argsort(regulatory_importance)[-n_features:]
    selected_features = selected_features[regulatory_importance[selected_features] > 0]
    if len(selected_features) == 0:
        return results
    assert j not in selected_features
    X_ = X[:, selected_features]

    # Define labels
    y_ = X[:, j]

    y_pred, y_target, r2s, models = [], [], [], []
    for t, (train_index, test_index) in enumerate(LeaveOneGroupOut().split(X_, y_, groups)):

        if reg_type == 'ridge':
            model = Ridge(random_state=random_state)
        elif reg_type == 'GB':
            model = lightgbm.LGBMRegressor(verbosity=-1, n_estimators=100, n_jobs=n_jobs, random_state=random_state)
        elif reg_type == 'RF':
            model = lightgbm.LGBMRegressor(boosting_type='rf', feature_fraction=0.05, verbosity=-1, n_estimators=100, n_jobs=n_jobs, random_state=random_state)
        else:
            raise NotImplementedError(f'Unknown model "{reg_type}"')

        X_train = X_[train_index, :]
        X_test = X_[test_index, :]
        y_train = y_[train_index]
        y_test = y_[test_index]
        
        X_train = X_train.toarray() if hasattr(X_train, "toarray") else X_train
        X_test = X_test.toarray() if hasattr(X_test, "toarray") else X_test
        y_train = y_train.toarray().ravel() if hasattr(y_train, "toarray") else np.ravel(y_train)
        y_test = y_test.toarray().ravel() if hasattr(y_test, "toarray") else np.ravel(y_test)
        model.fit(X_train, y_train)

        y_pred.append(model.predict(X_test))
        y_target.append(y_test)
        r2s.append(r2_score(y_target[-1], y_pred[-1]))
        models.append(model)

    y_pred = np.concatenate(y_pred, axis=0)
    y_target = np.concatenate(y_target, axis=0)
    
    results['r2'] = float(np.clip(r2_score(y_target, y_pred), 0, 1))
    results['avg-r2'] = float(np.clip(np.mean(r2s), 0, 1))
    results['models'] = models
    
    return results


def cross_validate(
        reg_type: str,
        gene_names: np.ndarray,
        tf_names: Set[str],
        X: np.ndarray,
        groups: np.ndarray,
        grn: np.ndarray,
        n_features: np.ndarray,
        n_jobs: int
) -> List[Dict[str, float]]:
    n_genes = len(grn)

    grn = fill_zeros_in_grn(grn)

    # Remove interactions when first gene in pair is not in TF list
    mask = np.isin(gene_names, list(tf_names))
    grn[~mask, :] = 0
    
    from joblib import Parallel, delayed
    from tqdm.auto import tqdm
    from tqdm_joblib import tqdm_joblib
    with tqdm_joblib(tqdm(total=n_genes, desc=f"{reg_type} CV")):
        results = Parallel(n_jobs=n_jobs, backend="loky")(
            delayed(cross_validate_gene)(
                reg_type, X, groups, grn, j, int(n_features[j]), n_jobs
            )
            for j in range(n_genes) if n_features[j] > 0
        )
    included_gene_names = [gene_names[j] for j in range(n_genes) if n_features[j] > 0]

    return {
        'gene_names': list(included_gene_names),
        'results': list(results)
    }

def static_approach(
        grn: np.ndarray,
        n_features: np.ndarray,
        X: np.ndarray,
        groups: np.ndarray,
        gene_names: List[str],
        tf_names: Set[str],
        reg_type: str,
        n_jobs:int
) -> float:

    # Cross-validate each gene using the inferred GRN to define select input features
    res = cross_validate(reg_type, gene_names, tf_names, X, groups, grn, n_features, n_jobs=n_jobs)
    
    r2_scores = np.asarray([res['results'][j]['avg-r2'] for j in range(len(res['results']))])
    mean_r2_scores = np.mean(r2_scores)

    return mean_r2_scores


def cross_validate_gene_raw(
        reg_type: str,
        X: np.ndarray,
        groups: np.ndarray,
        grn: np.ndarray,
        j: int,
        gene_names: np.ndarray,
        tf_names: Set[str],
        random_state: int = 0xCAFE,
        n_jobs: int = 10
) -> Dict[str, float]:
    
    results = {'r2': 0.0, 'avg-r2': 0.0}
    
    # Get all TFs that regulate gene j from the GRN (without feature selection)
    regulators = grn[:, j]
    selected_features = np.where(regulators != 0)[0]
    
    # Filter to only keep TFs
    selected_features = selected_features[np.isin(gene_names[selected_features], list(tf_names))]
    
    # Remove self-regulation if present
    selected_features = selected_features[selected_features != j]
    
    if len(selected_features) == 0:
        return results
    
    X_ = X[:, selected_features]
    y_ = X[:, j]

    y_pred, y_target, r2s, models = [], [], [], []
    for t, (train_index, test_index) in enumerate(LeaveOneGroupOut().split(X_, y_, groups)):

        if reg_type == 'ridge':
            model = Ridge(random_state=random_state)
        elif reg_type == 'GB':
            model = lightgbm.LGBMRegressor(verbosity=-1, n_estimators=100, n_jobs=n_jobs, random_state=random_state)
        elif reg_type == 'RF':
            model = lightgbm.LGBMRegressor(boosting_type='rf', feature_fraction=0.05, verbosity=-1, n_estimators=100, n_jobs=n_jobs, random_state=random_state)
        else:
            raise NotImplementedError(f'Unknown model "{reg_type}"')

        X_train = X_[train_index, :]
        X_test = X_[test_index, :]
        y_train = y_[train_index]
        y_test = y_[test_index]
        
        X_train = X_train.toarray() if hasattr(X_train, "toarray") else X_train
        X_test = X_test.toarray() if hasattr(X_test, "toarray") else X_test
        y_train = y_train.toarray().ravel() if hasattr(y_train, "toarray") else np.ravel(y_train)
        y_test = y_test.toarray().ravel() if hasattr(y_test, "toarray") else np.ravel(y_test)
        model.fit(X_train, y_train)

        y_pred.append(model.predict(X_test))
        y_target.append(y_test)
        r2s.append(r2_score(y_target[-1], y_pred[-1]))
        models.append(model)

    y_pred = np.concatenate(y_pred, axis=0)
    y_target = np.concatenate(y_target, axis=0)
    
    results['r2'] = float(np.clip(r2_score(y_target, y_pred), 0, 1))
    results['avg-r2'] = float(np.clip(np.mean(r2s), 0, 1))
    results['models'] = models
    
    return results


def raw_approach(
        grn: np.ndarray,
        X: np.ndarray,
        groups: np.ndarray,
        gene_names: np.ndarray,
        tf_names: Set[str],
        reg_type: str,
        n_jobs: int
) -> float:
    """
    Evaluate GRN using raw network structure without theta-based standardization.
    For each gene, use all its TFs from the network as predictors.
    """
    n_genes = len(gene_names)
    
    # Don't fill zeros or modify the GRN - use it as is
    # Just filter by TF constraint
    mask = np.isin(gene_names, list(tf_names))
    grn_filtered = grn.copy()
    grn_filtered[~mask, :] = 0
    
    from joblib import Parallel, delayed
    from tqdm.auto import tqdm
    from tqdm_joblib import tqdm_joblib
    
    with tqdm_joblib(tqdm(total=n_genes, desc=f"{reg_type} CV (raw)")):
        results = Parallel(n_jobs=n_jobs, backend="loky")(
            delayed(cross_validate_gene_raw)(
                reg_type, X, groups, grn_filtered, j, gene_names, tf_names, SEED, n_jobs
            )
            for j in range(n_genes)
        )
    
    r2_scores = np.asarray([results[j]['avg-r2'] for j in range(n_genes)])
    mean_r2_scores = np.mean(r2_scores)
    
    return mean_r2_scores


def main(par: Dict[str, Any]) -> pd.DataFrame:
    random_state = SEED
    np.random.seed(random_state)
    random.seed(random_state)
    lightgbm.LGBMRegressor().set_params(random_state=random_state)
    
    # Load data
    prturb_adata = ad.read_h5ad(par['evaluation_data'])
    layer = manage_layer(prturb_adata, par)
    gene_names = prturb_adata.var.index.to_numpy()
    with open(par['regulators_consensus'], 'r') as f:
        data = json.load(f)
    print(len(data), len(gene_names))
    
    n_features_theta_min = np.asarray([data[gene_name]['0'] for gene_name in gene_names], dtype=int)
    n_features_theta_median = np.asarray([data[gene_name]['0.5'] for gene_name in gene_names], dtype=int)
    n_features_theta_max = np.asarray([data[gene_name]['1'] for gene_name in gene_names], dtype=int)
    n_genes = len(gene_names)
    
    net = read_prediction(par)
    
    n_cells = prturb_adata.shape[0]
    random_groups = np.random.choice(range(1, 5+1), size=n_cells, replace=True) # random sampling
    groups = LabelEncoder().fit_transform(random_groups)

    # Load and standardize perturbation data    
    X = prturb_adata.layers[layer]
    X = RobustScaler(with_centering=False).fit_transform(X)

    tf_names = np.loadtxt(par['tf_all'], dtype=str)
    if par['apply_tf']==False:
        tf_names = gene_names

    # Check if group-specific evaluation is requested
    if 'group_specific' in par and par['group_specific'] is not None:
        # Group-specific evaluation
        group_column = par['group_specific']
        unique_groups = prturb_adata.obs[group_column].unique()
        
        scores_min_by_group = []
        scores_median_by_group = []
        scores_max_by_group = []
        scores_raw_by_group = []
        
        for group in unique_groups:
            print(f'Evaluating group: {group}', flush=True)
            
            # Subset evaluation data by group
            group_mask = prturb_adata.obs[group_column] == group
            X_group = X[group_mask, :]
            groups_group = groups[group_mask]
            
            # Subset network if group name is in net (check if net has group-specific data)
            net_group = net.copy()
            if 'group' in net.columns and group in net['group'].values:
                net_group = net[net['group'] == group]
            
            net_matrix_group = net_to_matrix(net_group, gene_names)
            
            # Evaluate GRN for this group
            if (n_features_theta_min!=0).any()==False:
                score_static_min = np.nan
            else:
                score_static_min = static_approach(net_matrix_group, n_features_theta_min, X_group, groups_group, gene_names, tf_names, par['reg_type'], n_jobs=par['num_workers'])
            
            score_static_median = static_approach(net_matrix_group, n_features_theta_median, X_group, groups_group, gene_names, tf_names, par['reg_type'], n_jobs=par['num_workers'])
            score_static_max = static_approach(net_matrix_group, n_features_theta_max, X_group, groups_group, gene_names, tf_names, par['reg_type'], n_jobs=par['num_workers'])
            score_raw = raw_approach(net_matrix_group, X_group, groups_group, gene_names, tf_names, par['reg_type'], n_jobs=par['num_workers'])
            
            scores_min_by_group.append(score_static_min)
            scores_median_by_group.append(score_static_median)
            scores_max_by_group.append(score_static_max)
            scores_raw_by_group.append(score_raw)
        
        # Take the mean across groups
        score_static_min = np.nanmean(scores_min_by_group)
        score_static_median = np.nanmean(scores_median_by_group)
        score_static_max = np.nanmean(scores_max_by_group)
        score_raw = np.nanmean(scores_raw_by_group)
        
    else:
        net_matrix = net_to_matrix(net, gene_names)
        # Original evaluation (no group-specific)
        print(f'Raw approach (no theta):', flush=True)
        score_raw = raw_approach(net_matrix, X, groups, gene_names, tf_names, par['reg_type'], n_jobs=par['num_workers'])

        if (n_features_theta_min!=0).any()==False:
            score_static_min = np.nan
        else:
            print(f'Static approach (theta=0):', flush=True)
            score_static_min = static_approach(net_matrix, n_features_theta_min, X, groups, gene_names, tf_names, par['reg_type'], n_jobs=par['num_workers'])
        print(f'Static approach (theta=.5):', flush=True)
        score_static_median = static_approach(net_matrix, n_features_theta_median, X, groups, gene_names, tf_names, par['reg_type'], n_jobs=par['num_workers'])
        print(f'Static approach (theta=1):', flush=True)
        score_static_max = static_approach(net_matrix, n_features_theta_max, X, groups, gene_names, tf_names, par['reg_type'], n_jobs=par['num_workers'])
        
    results = {
        'r2-theta-0.0': [float(score_static_min)],
        'r2-theta-0.5': [float(score_static_median)],
        'r2-theta-1.0': [float(score_static_max)],
        'R_raw': [float(score_raw)],
    }

    df_results = pd.DataFrame(results)

    return df_results
