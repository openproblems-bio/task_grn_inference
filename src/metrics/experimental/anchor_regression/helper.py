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
from sklearn.preprocessing import StandardScaler, LabelEncoder, OneHotEncoder
import anndata as ad
import warnings
warnings.filterwarnings("ignore", category=UserWarning)

# For reproducibility purposes
seed = 0xCAFE
os.environ['PYTHONHASHSEED'] = str(seed)
random.seed(seed)
np.random.seed(seed)

from util import read_prediction, manage_layer
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
        l2_reg: float = 1e-2,
        anchor_strength: float = 1.0
) -> np.ndarray:
    """Anchor regression for causal inference under confounding.

    Args:
        X: input features, of shape (n, d).
        Z: anchor variables. The model is required to be invariant to
            these environmental variables. Shape (n, z).
        Y: predicted variables, of shape (n, u).
        l2_reg: L2 regularization strength.
        anchor_strength: Strength of anchor regularization. 
                        0 = standard regression, higher = more anchor regularization.

    Returns:
        Inferred parameters, of shape (d, u).
    """

    # Whitening transformation
    W = get_whitening_transform(Z, gamma=anchor_strength)
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
        A: np.ndarray,
        is_selected: np.ndarray,
        eps: float = 1e-50
) -> float:
    theta0_signed = anchor_regression(X, Z, y, anchor_strength=1)
    theta0_signed = theta0_signed[is_selected]
    theta0 = np.abs(theta0_signed)
    theta0 /= np.sum(theta0)

    theta_signed = anchor_regression(X, Z, y, anchor_strength=20)
    theta_signed = theta_signed[is_selected]
    theta = np.abs(theta_signed)
    theta /= np.sum(theta)

    stabilities = np.clip((theta0 - theta) / (theta0 + eps), 0, 1)
    stabilities[np.sign(theta0_signed) != np.sign(theta_signed)] = 0

    return stabilities



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
        Stability score for gene j
    """
    is_selected = np.array(A[:, j] != 0)
    if (not np.any(is_selected)) or np.all(is_selected):
        return 0.0
    assert not is_selected[j]

    # Exclude target gene from features
    mask = np.ones(X.shape[1], dtype=bool)
    mask[j] = False

    stabilities_selected = np.mean(compute_stabilities(X[:, mask], X[:, j], Z, A, is_selected[mask], eps=eps))
    #stabilities_non_selected = np.mean(compute_stabilities(X[:, mask], X[:, j], Z, A, ~is_selected[mask], eps=eps))

    #score = (stabilities_selected - stabilities_non_selected) / (stabilities_selected + stabilities_non_selected + eps)
    score = np.mean(stabilities_selected)

    return score


def main(par):
    """Main anchor regression evaluation function."""
    # Load evaluation data
    adata = ad.read_h5ad(par['evaluation_data'])
    dataset_id = adata.uns['dataset_id']
    method_id = ad.read_h5ad(par['prediction'], backed='r').uns['method_id']
    
    # Get dataset-specific anchor variables
    if dataset_id not in DATASET_GROUPS:
        raise ValueError(f"Dataset {dataset_id} not found in DATASET_GROUPS")
    
    anchor_cols = DATASET_GROUPS[dataset_id].get('anchors', ['donor_id', 'plate_name'])
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

    # Remove self-regulations from GRN
    np.fill_diagonal(A, 0)
    print(f"Evaluating {X.shape[1]} genes with {np.sum(A != 0)} regulatory links")

    # Center and scale dataset
    scaler = StandardScaler()
    X = scaler.fit_transform(X)

    # Create baseline GRN
    A_baseline = np.copy(A)
    for j in range(A_baseline.shape[1]):
        np.random.shuffle(A_baseline[:j, j])
        np.random.shuffle(A_baseline[j+1:, j])
    assert np.any(A != A_baseline)

    scores, scores_baseline = [], []
    for j in tqdm.tqdm(range(X.shape[1]), desc="Evaluating gene stability"):
        scores.append(evaluate_gene_stability(X, Z, A, j))
    scores = np.array(scores)

    # Calculate final score
    final_score = np.mean(scores)
    print(f"Method: {method_id}")
    print(f"Anchor Regression Score: {final_score:.6f}")

    # Return results as DataFrame
    results = {
        'anchor_regression': [final_score]
    }

    df_results = pd.DataFrame(results)
    return df_results
