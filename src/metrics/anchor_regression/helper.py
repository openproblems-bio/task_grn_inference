from typing import Tuple
import os
import sys
import traceback
import random
import h5py
import numpy as np
import pandas as pd
import tqdm
from scipy.sparse.linalg import LinearOperator
from scipy.stats import pearsonr, wilcoxon
from scipy.sparse import csr_matrix
from sklearn.model_selection import KFold
from sklearn.preprocessing import StandardScaler, LabelEncoder, OneHotEncoder
import anndata as ad
import warnings
warnings.filterwarnings("ignore", category=UserWarning)

# For reproducibility purposes
seed = 0xCAFE
os.environ['PYTHONHASHSEED'] = str(seed)
random.seed(seed)
np.random.seed(seed)

from util import read_prediction, manage_layer, create_grn_baseline
from config import DATASET_GROUPS


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


def get_whitening_transform(Z: np.ndarray, gamma: float = 1.0) -> LinearOperator:
    """Get whitening transformation.

    Args:
        Z: Anchor variables.
        gamma: Anchor strength.

    Returns:
        Sparse linear operator corresponding to a whitening transform.
    """
    n, k = Z.shape

    # Precompute Gram inverse with jitter for stability
    ZTZ = Z.T @ Z + 1e-10 * np.eye(k)
    ZTZ_inv = np.linalg.inv(ZTZ)
    sqrt_gamma = np.sqrt(gamma)

    def matvec(v: np.ndarray) -> np.ndarray:
        """Matrix-vector multiplication."""
        v = np.atleast_2d(v)
        if v.shape[0] != n:
            v = v.T
        Pv = Z @ (ZTZ_inv @ (Z.T @ v))
        out = v + (sqrt_gamma - 1.0) * Pv
        return out if out.shape[1] > 1 else out.ravel()

    return LinearOperator((n, n), matvec=matvec, rmatvec=matvec)


def anchor_regression(
        X: np.ndarray,
        Z: np.ndarray,
        Y: np.ndarray,
        l2_reg: float = 1e-6,
        gamma: float = 1.0
) -> np.ndarray:
    """Anchor regression for causal inference under confounding.

    Args:
        X: input features, of shape (n, d).
        Z: anchor variables. The model is required to be invariant to
            these environmental variables. Shape (n, z).
        Y: predicted variables, of shape (n, u).
        l2_reg: L2 regularization strength.
        gamma: Strength of anchor regularization. 
            1 = standard regression, higher = more anchor regularization.

    Returns:
        Inferred parameters, of shape (d, u).
    """

    # Whitening transformation
    W = get_whitening_transform(Z, gamma=gamma)
    X_t = W @ X
    Y_t = W @ Y

    # Ridge regression on the whitened data
    sigma_xx = X_t.T @ X_t / X_t.shape[0]
    sigma_xy = X_t.T @ Y_t / Y_t.shape[0]
    theta = np.linalg.solve(sigma_xx + l2_reg * np.eye(X_t.shape[1]), sigma_xy)

    return theta


def compute_stabilities(
        X: np.ndarray,
        y: np.ndarray,
        Z: np.ndarray,
        is_selected: np.ndarray,
        eps: float = 1e-50
) -> float:

    theta0 = np.abs(anchor_regression(X, Z, y, gamma=1))
    theta0 /= np.sum(theta0)
    theta = np.abs(anchor_regression(X, Z, y, gamma=1.2))
    theta /= np.sum(theta)
    s1 = theta[is_selected] * theta0[is_selected]
    s2 = theta[~is_selected] * theta0[~is_selected]

    s1 = np.mean(s1)
    s2 = np.mean(s2)

    stability = (s1 - s2) / (s1 + s2 + eps)

    return stability


def evaluate_gene_stability(
        X: np.ndarray,
        Z: np.ndarray,
        A: np.ndarray,
        j: int,
        eps: float = 1e-50
) -> float:
    """Evaluate stability of regulatory relationships for a single target gene.
    
    Args:
        X: Gene expression matrix (n_samples, n_genes).
        Z: Anchor variables matrix (n_samples, n_anchors).
        A: GRN adjacency matrix (n_genes, n_genes).
        j: Target gene index.
        eps: Small epsilon for numerical stability.
        
    Returns:
        Stability score for gene j.
    """
    is_selected = np.array(A[:, j] != 0)
    if (not np.any(is_selected)) or np.all(is_selected):
        return np.nan
    assert not is_selected[j]

    # Exclude target gene from features
    mask = np.ones(X.shape[1], dtype=bool)
    mask[j] = False

    return compute_stabilities(X[:, mask], X[:, j], Z, is_selected[mask])


def main(par):
    """Main anchor regression evaluation function."""

    # Load evaluation data
    adata = ad.read_h5ad(par['evaluation_data'])
    dataset_id = adata.uns['dataset_id']
    method_id = ad.read_h5ad(par['prediction'], backed='r').uns['method_id']

    # Manage layer
    layer = manage_layer(adata, par)
    X = adata.layers[layer]
    if isinstance(X, csr_matrix):
        X = X.toarray()
    X = X.astype(np.float32)

    gene_names = adata.var_names
    gene_dict = {gene_name: i for i, gene_name in enumerate(gene_names)}

    # Get dataset-specific anchor variables
    if dataset_id not in DATASET_GROUPS:
        raise ValueError(f"Dataset {dataset_id} not found in DATASET_GROUPS")
    anchor_cols = DATASET_GROUPS[dataset_id].get('anchors', ['donor_id', 'plate_name'])
    print(f"Using anchor variables: {anchor_cols}")

    # Encode anchor variables
    anchor_variables = encode_obs_cols(adata, anchor_cols)
    anchor_encoded = combine_multi_index(*anchor_variables)

    # Get CV groups
    if "cell_type" in adata.obs:
        cv_groups = LabelEncoder().fit_transform(adata.obs["cell_type"].values)
    else:
        np.random.randint(0, 5)
        cv_groups = np.random.shuffle()
    
    if len(anchor_variables) == 0:
        raise ValueError(f"No anchor variables found in dataset for columns: {anchor_cols}")
    
    # One-hot encode anchor variables
    Z = OneHotEncoder(drop="first", sparse_output=False, dtype=np.float32).fit_transform(anchor_encoded.reshape(-1, 1))
    print(f"Anchor matrix Z shape: {Z.shape}")

    # Add intercept
    Z = np.concatenate((Z, np.ones((len(Z), 1), dtype=np.float32)), axis=1)

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

    # Remove self-regulations from GRN
    np.fill_diagonal(A, 0)
    print(f"Evaluating {X.shape[1]} genes with {np.sum(A != 0)} regulatory links")

    # Whether or not to take into account the regulatory modes (enhancer/inhibitor)
    signed = np.any(A < 0)

    # Center and scale dataset
    scaler = StandardScaler()
    X = scaler.fit_transform(X)

    # Create baseline model
    A_baseline = create_grn_baseline(A)

    # Compute gene stabilities
    scores, scores_baseline = [], []
    for j in tqdm.tqdm(range(X.shape[1]), desc="Evaluating gene stability"):
        scores.append(evaluate_gene_stability(X, Z, A, j))
        scores_baseline.append(evaluate_gene_stability(X, Z, A_baseline, j))
    scores = np.array(scores)
    scores_baseline = np.array(scores_baseline)

    # Skip NaNs
    mask = ~np.logical_or(np.isnan(scores), np.isnan(scores_baseline))
    scores = scores[mask]
    scores_baseline = scores_baseline[mask]

    # Calculate final score
    p_value = wilcoxon(scores, scores_baseline).pvalue
    p_value = np.clip(p_value, 1e-300, 1)
    final_score = -np.log10(p_value)
    print(f"Method: {method_id}")
    print(f"Anchor Regression Score: {final_score:.6f}")

    # Return results as DataFrame
    results = {
        'anchor_regression': [final_score]
    }

    df_results = pd.DataFrame(results)
    return df_results