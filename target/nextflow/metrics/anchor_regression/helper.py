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
        eps: float = 1e-50,
        W1=None,
        W2=None
) -> float:
    """Compute stability score with optional pre-computed whitening transforms."""
    
    # Use pre-computed whitening transforms if available
    if W1 is None:
        W1 = get_whitening_transform(Z, gamma=1.0)
    if W2 is None:
        W2 = get_whitening_transform(Z, gamma=1.2)
    
    # Compute anchor regression coefficients for both gamma values
    X_t1 = W1 @ X
    y_t1 = W1 @ y
    sigma_xx1 = X_t1.T @ X_t1 / X_t1.shape[0]
    sigma_xy1 = X_t1.T @ y_t1 / y_t1.shape[0]
    theta0 = np.linalg.solve(sigma_xx1 + 1e-6 * np.eye(X_t1.shape[1]), sigma_xy1)
    theta0 = np.abs(theta0.ravel())
    theta0 /= np.sum(theta0) + eps
    
    X_t2 = W2 @ X
    y_t2 = W2 @ y
    sigma_xx2 = X_t2.T @ X_t2 / X_t2.shape[0]
    sigma_xy2 = X_t2.T @ y_t2 / y_t2.shape[0]
    theta = np.linalg.solve(sigma_xx2 + 1e-6 * np.eye(X_t2.shape[1]), sigma_xy2)
    theta = np.abs(theta.ravel())
    theta /= np.sum(theta) + eps
    
    s1 = np.mean(theta[is_selected] * theta0[is_selected])
    s2 = np.mean(theta[~is_selected] * theta0[~is_selected])

    stability = (s1 - s2) / (s1 + s2 + eps)

    return stability


def evaluate_gene_stability(
        X: np.ndarray,
        Z: np.ndarray,
        A: np.ndarray,
        j: int,
        eps: float = 1e-50,
        W1=None,
        W2=None
) -> float:
    """Evaluate stability of regulatory relationships for a single target gene.
    
    Args:
        X: Gene expression matrix (n_samples, n_genes).
        Z: Anchor variables matrix (n_samples, n_anchors).
        A: GRN adjacency matrix (n_genes, n_genes).
        j: Target gene index.
        eps: Small epsilon for numerical stability.
        W1: Pre-computed whitening transform for gamma=1.0 (optional).
        W2: Pre-computed whitening transform for gamma=1.2 (optional).
        
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

    return compute_stabilities(X[:, mask], X[:, j], Z, is_selected[mask], eps=eps, W1=W1, W2=W2)


def main(par):
    adata = ad.read_h5ad(par['evaluation_data'])
    dataset_id = adata.uns['dataset_id']
    method_id = ad.read_h5ad(par['prediction'], backed='r').uns['method_id']
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
    anchor_cols = DATASET_GROUPS[dataset_id]['anchors']
    print(f"Using anchor variables: {anchor_cols}")
    
    # Validate anchor columns exist in adata.obs
    missing_cols = [col for col in anchor_cols if col not in adata.obs.columns]
    if missing_cols:
        raise ValueError(f"Anchor columns {missing_cols} not found in adata.obs. Available columns: {adata.obs.columns.tolist()}")

    # Encode anchor variables
    anchor_variables = encode_obs_cols(adata, anchor_cols)
    
    if len(anchor_variables) == 0:
        raise ValueError(f"No anchor variables could be encoded from columns: {anchor_cols}")
    
    anchor_encoded = combine_multi_index(*anchor_variables)
    
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
    
    # Pre-compute whitening transforms (shared across all genes)
    print("Pre-computing whitening transforms...")
    W1 = get_whitening_transform(Z, gamma=1.0)
    W2 = get_whitening_transform(Z, gamma=1.2)

    # Create baseline model
    A_baseline = create_grn_baseline(A)

    # Compute gene stabilities in parallel
    def eval_gene(j):
        return (evaluate_gene_stability(X, Z, A, j, W1=W1, W2=W2),
                evaluate_gene_stability(X, Z, A_baseline, j, W1=W1, W2=W2))
    
    n_jobs = par['num_workers']
    from joblib import Parallel, delayed
    
    results = Parallel(n_jobs=n_jobs, backend='loky')(
        delayed(eval_gene)(j) for j in tqdm.tqdm(range(X.shape[1]), desc="Evaluating gene stability")
    )
    
    scores = np.array([r[0] for r in results])
    scores_baseline = np.array([r[1] for r in results])

    # Skip NaNs
    mask = ~np.logical_or(np.isnan(scores), np.isnan(scores_baseline))
    scores = scores[mask]
    scores_baseline = scores_baseline[mask]
    
    if len(scores) < 10:
        raise ValueError(f"Too few valid genes ({len(scores)}) for anchor regression evaluation")

    # Calculate final score
    p_value = wilcoxon(scores, scores_baseline).pvalue
    p_value = np.clip(p_value, 1e-300, 1)
    
    # Calculate raw score (mean difference)
    raw_score = float(np.mean(scores) - np.mean(scores_baseline))
    
    # Set to 0 if not significant (p >= 0.05)
    if p_value >= 0.05:
        final_score = 0.0
        print(f"Method: {method_id}")
        print(f"p-value: {p_value:.6f} (not significant, p >= 0.05)")
        print(f"Anchor Regression Score: {final_score:.6f} (set to 0)")
        print(f"Anchor Regression Raw Score: {raw_score:.6f}")
    else:
        final_score = -np.log10(p_value)
        print(f"Method: {method_id}")
        print(f"p-value: {p_value:.6f} (significant)")
        print(f"Anchor Regression Score: {final_score:.6f}")
        print(f"Anchor Regression Raw Score: {raw_score:.6f}")

    # Return results as DataFrame
    results = {
        'anchor_regression': [final_score],
        'anchor_regression_raw': [raw_score]
    }

    df_results = pd.DataFrame(results)
    return df_results