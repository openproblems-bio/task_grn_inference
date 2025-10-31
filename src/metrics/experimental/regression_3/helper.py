from typing import Tuple
import os
import sys
import traceback
import random
import h5py
import numpy as np
import pandas as pd
import tqdm
import xgboost
from scipy.sparse.linalg import LinearOperator
from scipy.stats import pearsonr, spearmanr, wilcoxon, ConstantInputWarning
from scipy.sparse import csr_matrix
from sklearn.preprocessing import StandardScaler, LabelEncoder, OneHotEncoder
from sklearn.linear_model import Ridge
import anndata as ad
import warnings
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=ConstantInputWarning)

# For reproducibility purposes
seed = 0xCAFE
os.environ['PYTHONHASHSEED'] = str(seed)
random.seed(seed)
np.random.seed(seed)

from util import read_prediction, manage_layer
from dataset_config import DATASET_GROUPS
from baseline import create_grn_baseline


def encode_obs_cols(adata, cols):
    """Encode observation columns to integer codes."""
    encoded = []
    for col in cols:
        if col in adata.obs:
            codes = LabelEncoder().fit_transform(adata.obs[col].values)
            encoded.append(codes)
    return encoded


def combine_multi_index(*arrays) -> np.array:
    """Combine parallel label arrays into a single integer label per position."""
    A = np.stack(arrays)
    n_classes = tuple(int(A[i].max()) + 1 for i in range(A.shape[0]))
    return np.ravel_multi_index(A, dims=n_classes, order='C')


def compute_residual_correlations(
        X_train: np.ndarray,
        y_train: np.ndarray,
        X_test: np.ndarray,
        y_test: np.ndarray,
        Z_test: np.ndarray
) -> np.ndarray:
    model = xgboost.XGBRegressor(n_estimators=10)
    #model = xgboost.XGBRegressor(n_estimators=10)
    model = Ridge(alpha=1)
    model.fit(X_train, y_train)
    y_hat = model.predict(X_test)
    residuals = y_test - y_hat
    coefs = pearsonr(residuals[:, np.newaxis], Z_test, axis=0)[0]
    coefs = np.nan_to_num(coefs, nan=0)
    assert coefs.shape[0] == Z_test.shape[1]
    return np.abs(coefs)


def main(par):
    # Load evaluation data
    adata = ad.read_h5ad(par['evaluation_data'])
    dataset_id = adata.uns['dataset_id']
    method_id = ad.read_h5ad(par['prediction'], backed='r').uns['method_id']
    
    # Get dataset-specific anchor variables
    if dataset_id not in DATASET_GROUPS:
        raise ValueError(f"Dataset {dataset_id} not found in DATASET_GROUPS")
    
    anchor_cols = DATASET_GROUPS[dataset_id]['anchors']
    print(f"Using anchor variables: {anchor_cols}")

    # Manage layer
    layer = manage_layer(adata, par)
    X = adata.layers[layer]
    if isinstance(X, csr_matrix):
        X = X.toarray()
    X = X.astype(np.float32)

    gene_names = adata.var_names
    gene_dict = {gene_name: i for i, gene_name in enumerate(gene_names)}

    # Encode anchor variables
    anchor_variables = encode_obs_cols(adata, anchor_cols)
    anchor_encoded = combine_multi_index(*anchor_variables)
    
    if len(anchor_variables) == 0:
        raise ValueError(f"No anchor variables found in dataset for columns: {anchor_cols}")
    
    # One-hot encode anchor variables
    Z = OneHotEncoder(sparse_output=False, dtype=np.float32).fit_transform(anchor_encoded.reshape(-1, 1))
    print(f"Anchor matrix Z shape: {Z.shape}")

    # Load inferred GRN
    df = read_prediction(par)
    sources = df["source"].to_numpy()
    targets = df["target"].to_numpy()
    weights = df["weight"].to_numpy()

    A = np.zeros((len(gene_names), len(gene_names)), dtype=X.dtype)
    for source, target, weight in zip(sources, targets, weights):
        if (source in gene_dict) and (target in gene_dict):
            i = gene_dict[source]
            j = gene_dict[target]
            A[i, j] = float(weight)

    # Only consider the genes that are actually present in the inferred GRN,
    # and keep only the most-connected genes (for speed).
    gene_mask = np.logical_or(np.any(A, axis=1), np.any(A, axis=0))
    in_degrees = np.sum(A != 0, axis=0)
    out_degrees = np.sum(A != 0, axis=1)
    idx = np.argsort(np.maximum(out_degrees, in_degrees))[:-1000]
    gene_mask[idx] = False
    X = X[:, gene_mask]
    X = X.toarray() if isinstance(X, csr_matrix) else X
    A = A[gene_mask, :][:, gene_mask]
    gene_names = gene_names[gene_mask]

    # Remove self-regulations
    np.fill_diagonal(A, 0)
    print(f"Evaluating {X.shape[1]} genes with {np.sum(A != 0)} regulatory links")

    # Create baseline model
    A_baseline = create_grn_baseline(A)

    scores, baseline_scores = [], []
    for group in np.unique(anchor_encoded):

        # Train/test split
        mask = (anchor_encoded != group)
        X_train = X[mask, :]
        X_test = X[~mask, :]

        # Standardize features
        scaler = StandardScaler()
        X_train = scaler.fit_transform(X_train)
        X_test = scaler.transform(X_test)

        for j in tqdm.tqdm(range(X_train.shape[1])):

            # Evaluate inferred GRN
            selected = (A[:, j] != 0)
            unselected = ~np.copy(selected)
            unselected[j] = False
            if (not np.any(selected)) or (not np.any(unselected)):
                continue
            else:
                coefs = compute_residual_correlations(
                    X_train[:, selected],
                    X_train[:, j],
                    X_test[:, selected],
                    X_test[:, j],
                    X_test[:, ~selected]
                )
                scores.append(np.mean(coefs))

            # Evaluate baseline GRN
            selected = (A_baseline[:, j] != 0)
            unselected = ~np.copy(selected)
            unselected[j] = False
            coefs = compute_residual_correlations(
                X_train[:, selected],
                X_train[:, j],
                X_test[:, selected],
                X_test[:, j],
                X_test[:, ~selected]
            )
            baseline_scores.append(np.mean(coefs))
    scores = np.array(scores)
    baseline_scores = np.array(baseline_scores)
    reg3_lift = np.mean(scores) / (np.mean(baseline_scores) + 1e-6)
    p_value = wilcoxon(baseline_scores, scores, alternative="greater").pvalue
    p_value = max(p_value, 1e-300)

    # Calculate final score
    final_score = -np.log10(p_value)
    print(f"Anchor Regression Score: {final_score:.6f}")
    print(f"Method: {method_id}")

    # Return results as DataFrame
    results = {
        'reg3_precision': [reg3_lift],
        'reg3_balanced': [final_score]
    }

    df_results = pd.DataFrame(results)
    return df_results