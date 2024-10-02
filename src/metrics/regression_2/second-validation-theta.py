import sys
sys.path.append('/lustre1/project/stg_00019/research/Antoine/dependencies')

import os
from typing import Dict, List, Tuple, Any, Union
import random
import json
import argparse

import tqdm
import numpy as np
import lightgbm
import pandas as pd
import anndata as ad
from sklearn.preprocessing import LabelEncoder, RobustScaler, StandardScaler
from sklearn.model_selection import GroupKFold, LeaveOneGroupOut
from sklearn.linear_model import Ridge
from sklearn.metrics import r2_score, mean_squared_error


DATA_DIR = os.path.join('../../', 'output', 'preprocess')
GRN_DIR = os.path.join('..', 'output', 'benchmark', 'grn_models')
BASELINE_GRN_DIR = os.path.join('..', 'output', 'benchmark', 'baseline_models')
RESULTS_DIR = os.path.join('..', 'output', 'benchmark', 'second-validation')
RESOURCES_DIR = os.path.join('..', 'resources')
os.makedirs(RESULTS_DIR, exist_ok=True)



SEED = 0xCAFE


parser = argparse.ArgumentParser()
parser.add_argument('estimator', type=str)
parser.add_argument('norm', type=str)
args = parser.parse_args()

estimator_t = args.estimator
norm_t = args.norm
override_ = False


def load_grn(filepath: str, gene_names: np.ndarray) -> np.ndarray:
    gene_dict = {gene_name: i for i, gene_name in enumerate(gene_names)}
    A = np.zeros((len(gene_names), len(gene_names)), dtype=float)
    df = pd.read_csv(filepath, sep=',', header='infer')
    for source, target, weight in zip(df['source'], df['target'], df['weight']):
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
        estimator_t: str,
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
    scores = np.abs(grn[:, j])
    scores[j] = -1
    selected_features = np.argsort(scores)[-n_features:]
    assert j not in selected_features
    X_ = X[:, selected_features]

    # Define labels
    y_ = X[:, j]

    y_pred, y_target, r2s = [], [], []
    for t, (train_index, test_index) in enumerate(LeaveOneGroupOut().split(X_, y_, groups)):

        if estimator_t == 'ridge':
            model = Ridge(random_state=random_state)
        elif estimator_t == 'GB':
            model = lightgbm.LGBMRegressor(verbosity=-1, n_estimators=100, n_jobs=n_jobs, random_state=random_state)
        else:
            raise NotImplementedError(f'Unknown model "{estimator_t}"')

        X_train = X_[train_index, :]
        X_test = X_[test_index, :]
        y_train = y_[train_index]
        y_test = y_[test_index]
        
        model.fit(X_train, y_train)

        y_pred.append(model.predict(X_test))
        y_target.append(y_test)
        r2s.append(r2_score(y_target[-1], y_pred[-1]))

    y_pred = np.concatenate(y_pred, axis=0)
    y_target = np.concatenate(y_target, axis=0)
    
    results['r2'] = float(r2_score(y_target, y_pred))
    results['avg-r2'] = float(np.mean(r2s))

    
    return results


def cross_validate(
        estimator_t: str,
        gene_names: np.ndarray,
        X: np.ndarray,
        groups: np.ndarray,
        grn: np.ndarray,
        n_features: np.ndarray,
        n_jobs: int
) -> List[Dict[str, float]]:
    n_genes = len(grn)
    
    grn = fill_zeros_in_grn(grn)
    
    results = []
    #for j in tqdm.tqdm(range(n_genes), desc=f'{estimator_t} CV'):
    for j in range(n_genes):
        res = cross_validate_gene(estimator_t, X, groups, grn, j, n_features=int(n_features[j]),n_jobs=n_jobs)
        results.append(res)
    
    return {
        'gene_names': list(gene_names),
        'scores': list(results)
    }


def static_approach(grn: np.ndarray, n_features: np.ndarray, X: np.ndarray, groups: np.ndarray, gene_names: List[str], reg_type: str, n_jobs:int) -> float:

    # Cross-validate each gene using the inferred GRN to define select input features
    res = cross_validate(reg_type, gene_names, X, groups, grn, n_features, n_jobs=n_jobs)
    return res


# Set global seed for reproducibility purposes
random_state = SEED
np.random.seed(random_state)
random.seed(random_state)
lightgbm.LGBMRegressor().set_params(random_state=random_state)

print('Reading input files', flush=True)

# Load perturbation data
perturbation_data = ad.read_h5ad(os.path.join(DATA_DIR, 'bulk_adata_integrated.h5ad'))
gene_names = perturbation_data.var.index.to_numpy()
n_genes = len(gene_names)
groups = LabelEncoder().fit_transform(perturbation_data.obs.plate_name)


for theta in ['0', '0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9', '1']:
    for method in ['ananse', 'celloracle', 'collectRI', 'figr', 'granie', 'negative_control', 'positive_control', 'scenicplus', 'scglue']:
        print(method, theta)

        # Load inferred GRN
        filepath = os.path.join('..', 'output', 'benchmark', 'grn_models', f'{method}.csv')
        grn = load_grn(filepath, gene_names)

        # Load and standardize perturbation data
        X = perturbation_data.layers[norm_t]

        X = RobustScaler().fit_transform(X)

        # Load consensus numbers of putative regulators
        with open(os.path.join(RESOURCES_DIR, 'prior', 'consensus-num-regulators.json'), 'r') as f:
            data = json.load(f)

        n_features = np.asarray([data[gene_name][theta] for gene_name in gene_names], dtype=int)

        # Evaluate GRN
        res = static_approach(grn, n_features, X, groups, gene_names, estimator_t, n_jobs=4)

        folder = os.path.join('..', 'output', 'benchmark', 'second-validation', norm_t, estimator_t, theta)
        os.makedirs(folder, exist_ok=True)
        with open(os.path.join(folder, f'{method}.results.json'), 'w') as f:
            json.dump(res, f)
