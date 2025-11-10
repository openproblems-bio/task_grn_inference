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


def create_regression_model(reg_type: str, random_state: int = 0xCAFE, n_jobs: int = 10):
    """Create a regression model based on the specified type."""
    if reg_type == 'ridge':
        return Ridge(random_state=random_state)
    elif reg_type == 'GB':
        return lightgbm.LGBMRegressor(verbosity=-1, n_estimators=100, n_jobs=n_jobs, random_state=random_state)
    elif reg_type == 'RF':
        return lightgbm.LGBMRegressor(boosting_type='rf', feature_fraction=0.05, verbosity=-1, n_estimators=100, n_jobs=n_jobs, random_state=random_state)
    else:
        raise NotImplementedError(f'Unknown model "{reg_type}"')


def to_dense_array(arr):
    """Convert sparse matrix to dense array if needed."""
    return arr.toarray() if hasattr(arr, "toarray") else arr


def perform_cross_validation(
        reg_type: str,
        X_: np.ndarray,
        y_: np.ndarray,
        groups: np.ndarray,
        random_state: int = 0xCAFE,
        n_jobs: int = 10
) -> Dict[str, Any]:
    """Perform leave-one-group-out cross-validation."""
    y_pred, y_target, r2s, models = [], [], [], []
    
    for train_index, test_index in LeaveOneGroupOut().split(X_, y_, groups):
        model = create_regression_model(reg_type, random_state, n_jobs)
        
        X_train = to_dense_array(X_[train_index, :])
        X_test = to_dense_array(X_[test_index, :])
        y_train = np.ravel(to_dense_array(y_[train_index]))
        y_test = np.ravel(to_dense_array(y_[test_index]))
        
        model.fit(X_train, y_train)
        
        y_pred.append(model.predict(X_test))
        y_target.append(y_test)
        r2s.append(r2_score(y_target[-1], y_pred[-1]))
        models.append(model)
    
    y_pred = np.concatenate(y_pred, axis=0)
    y_target = np.concatenate(y_target, axis=0)
    
    return {
        'r2': float(np.clip(r2_score(y_target, y_pred), 0, 1)),
        'avg-r2': float(np.clip(np.mean(r2s), 0, 1)),
        'models': models
    }


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
    y_ = X[:, j]
    
    return perform_cross_validation(reg_type, X_, y_, groups, random_state, n_jobs)


def cross_validate(
        reg_type: str,
        gene_names: np.ndarray,
        tf_names: Set[str],
        X: np.ndarray,
        groups: np.ndarray,
        grn: np.ndarray,
        n_features: np.ndarray,
        n_jobs: int,
        theta: str = None
) -> pd.DataFrame:
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
    included_n_features = [int(n_features[j]) for j in range(n_genes) if n_features[j] > 0]
    
    # Create detailed DataFrame
    detailed_results = []
    for gene_name, result, n_tf in zip(included_gene_names, results, included_n_features):
        detailed_results.append({
            'gene': gene_name,
            'r2': result['avg-r2'],
            'n_regulators': n_tf,
            'theta': theta
        })
    
    return pd.DataFrame(detailed_results)

def static_approach(
        grn: np.ndarray,
        n_features: np.ndarray,
        X: np.ndarray,
        groups: np.ndarray,
        gene_names: List[str],
        tf_names: Set[str],
        reg_type: str,
        n_jobs:int,
        theta: str = None
) -> pd.DataFrame:

    # Cross-validate each gene using the inferred GRN to define select input features
    detailed_df = cross_validate(reg_type, gene_names, tf_names, X, groups, grn, n_features, n_jobs=n_jobs, theta=theta)
    
    return detailed_df


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
    
    return perform_cross_validation(reg_type, X_, y_, groups, random_state, n_jobs)


def raw_approach(
        grn: np.ndarray,
        X: np.ndarray,
        groups: np.ndarray,
        gene_names: np.ndarray,
        tf_names: Set[str],
        reg_type: str,
        n_jobs: int
) -> pd.DataFrame:
    """
    Evaluate GRN using raw network structure without theta-based standardization.
    For each gene, use all its TFs from the network as predictors.
    Only evaluates genes with more than 2 regulators (excluding self-regulation).
    """
    n_genes = len(gene_names)
    
    # Don't fill zeros or modify the GRN - use it as is
    # Just filter by TF constraint
    mask = np.isin(gene_names, list(tf_names))
    grn_filtered = grn.copy()
    grn_filtered[~mask, :] = 0
    
    # Count regulators for each gene to determine which genes to evaluate
    n_regulators_per_gene = []
    genes_to_evaluate = []
    for j in range(n_genes):
        regulators = grn_filtered[:, j]
        # Count non-zero regulators excluding self-regulation
        regulator_indices = np.where(regulators != 0)[0]
        regulator_indices = regulator_indices[regulator_indices != j]  # Exclude self
        n_regulators = len(regulator_indices)
        
        if n_regulators > 2:
            n_regulators_per_gene.append(n_regulators)
            genes_to_evaluate.append(j)
    
    from joblib import Parallel, delayed
    from tqdm.auto import tqdm
    from tqdm_joblib import tqdm_joblib
    
    with tqdm_joblib(tqdm(total=len(genes_to_evaluate), desc=f"{reg_type} CV (raw)")):
        results = Parallel(n_jobs=n_jobs, backend="loky")(
            delayed(cross_validate_gene_raw)(
                reg_type, X, groups, grn_filtered, j, gene_names, tf_names, SEED, n_jobs
            )
            for j in genes_to_evaluate
        )
    
    # Create detailed DataFrame only for genes with regulators
    detailed_results = []
    for j, result, n_tf in zip(genes_to_evaluate, results, n_regulators_per_gene):
        detailed_results.append({
            'gene': gene_names[j],
            'r2': result['avg-r2'],
            'n_regulators': n_tf,
            'theta': 'r2_raw'
        })
    
    return pd.DataFrame(detailed_results)



def evaluate_static_approach(
        net_matrix: np.ndarray,
        n_features: np.ndarray,
        X: np.ndarray,
        groups: np.ndarray,
        gene_names: np.ndarray,
        tf_names: Set[str],
        reg_type: str,
        n_jobs: int,
        theta: str
) -> pd.DataFrame:
    """Evaluate GRN using static approach with given theta."""
    print(f'Static approach (theta={theta}):', flush=True)
    return static_approach(net_matrix, n_features, X, groups, gene_names, tf_names, reg_type, n_jobs=n_jobs, theta=theta)


def evaluate_network_for_group(
        net_matrix: np.ndarray,
        n_features_theta_min: np.ndarray,
        n_features_theta_median: np.ndarray,
        n_features_theta_max: np.ndarray,
        X: np.ndarray,
        groups: np.ndarray,
        gene_names: np.ndarray,
        tf_names: Set[str],
        reg_type: str,
        n_jobs: int
) -> List[pd.DataFrame]:
    """Evaluate a network using multiple theta values and raw approach."""
    detailed_scores = []
    
    # Evaluate with different theta values
    assert n_features_theta_min.sum() > 0, "No gene has regulators with theta=0.1"
    detailed_scores.append(evaluate_static_approach(
        net_matrix, n_features_theta_min, X, groups, gene_names, tf_names, reg_type, n_jobs, theta='r2-theta-0.1'
    ))
    detailed_scores.append(evaluate_static_approach(
        net_matrix, n_features_theta_median, X, groups, gene_names, tf_names, reg_type, n_jobs, theta='r2-theta-0.5'
    ))
    detailed_scores.append(evaluate_static_approach(
        net_matrix, n_features_theta_max, X, groups, gene_names, tf_names, reg_type, n_jobs, theta='r2-theta-1.0'
    ))
    
    # Evaluate with raw approach
    print(f'Raw approach (no theta):', flush=True)
    detailed_scores.append(raw_approach(net_matrix, X, groups, gene_names, tf_names, reg_type, n_jobs=n_jobs))
    
    return detailed_scores


def main(par: Dict[str, Any]) -> Tuple[pd.DataFrame, pd.DataFrame]:
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
    
    n_features_theta_min = np.asarray([data[gene_name]['0.1'] for gene_name in gene_names], dtype=int)
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

    detailed_scores = []
    
    # Check if group-specific evaluation is requested
    if 'group_specific' in par and par['group_specific'] is not None:
        # Group-specific evaluation
        group_column = par['group_specific']
        unique_groups = prturb_adata.obs[group_column].unique()
        
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
            
            # Evaluate network for this group
            group_scores = evaluate_network_for_group(
                net_matrix_group, n_features_theta_min, n_features_theta_median, n_features_theta_max,
                X_group, groups_group, gene_names, tf_names, par['reg_type'], par['num_workers']
            )
            detailed_scores.extend(group_scores)
        
    else:
        # Original evaluation (no group-specific)
        net_matrix = net_to_matrix(net, gene_names)
            
        # Evaluate network
        all_scores = evaluate_network_for_group(
            net_matrix, n_features_theta_min, n_features_theta_median, n_features_theta_max,
            X, groups, gene_names, tf_names, par['reg_type'], par['num_workers']
        )
        detailed_scores.extend(all_scores)
    
    # Combine all detailed scores
    detailed_df = pd.concat(detailed_scores, ignore_index=True)
    
    # Compute mean scores per theta
    mean_scores = detailed_df.groupby('theta')['r2'].mean().to_frame().T.reset_index(drop=True)

    return detailed_df, mean_scores
