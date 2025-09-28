import os
import sys
import traceback
import random
import h5py
import numpy as np
import pandas as pd
import tqdm
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
from dataset_config import DATASET_GROUPS


def encode_obs_cols(adata, cols):
    """Encode observation columns to integer codes."""
    encoded = []
    for col in cols:
        if col in adata.obs:
            codes = LabelEncoder().fit_transform(adata.obs[col].values)
            encoded.append(codes)
    return encoded


def anchor_regression(
        X: np.ndarray,
        Z: np.ndarray,
        Y: np.ndarray,
        l2_reg: float = 1e-3,
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
    # Projector onto the anchor subspace
    pi = Z @ np.linalg.pinv(Z.T @ Z) @ Z.T

    # Whitening transformation
    W = np.eye(len(Z)) + (np.sqrt(1.0 + anchor_strength) - 1.0) * pi
    X_t = W @ X
    Y_t = W @ Y

    # Ridge regression on the whitened data
    theta = np.linalg.solve(X_t.T @ X_t + l2_reg * np.eye(X_t.shape[1]), X_t.T @ Y_t)

    return theta


def evaluate_gene_stability(X, Z, A, j, eps=1e-10):
    """Evaluate stability of regulatory relationships for a single target gene.
    
    Args:
        X: Gene expression matrix (n_samples, n_genes)
        Z: Anchor variables matrix (n_samples, n_anchors)  
        A: GRN adjacency matrix (n_genes, n_genes)
        j: Target gene index
        eps: Small epsilon for numerical stability
        
    Returns:
        Stability score for gene j
    """
    is_selected = np.array(A[:, j] != 0)
    if (not np.any(is_selected)) or np.all(is_selected):
        return 0.0

    # Exclude target gene from predictors
    mask = np.ones(X.shape[1], dtype=bool)
    mask[j] = False

    # Compare standard vs anchor-regularized regression
    theta0 = anchor_regression(X[:, mask], Z, X[:, j], anchor_strength=0)
    theta0 = np.abs(theta0)

    theta = anchor_regression(X[:, mask], Z, X[:, j], anchor_strength=10.0)
    theta = np.abs(theta)

    is_selected = is_selected[mask]

    # Calculate relative coefficient changes
    selected_diff = float(
        np.mean(np.clip((theta0[is_selected] - theta[is_selected]) / (theta0[is_selected] + eps), 0, 1)))
    unselected_diff = float(
        np.mean(np.clip((theta0[~is_selected] - theta[~is_selected]) / (theta0[~is_selected] + eps), 0, 1)))
    
    # Score: Good GRNs have stable selected regulators, unstable unselected ones
    score = np.clip((unselected_diff - selected_diff) / (unselected_diff + selected_diff + eps), -1, 1)
    
    if np.isnan(score):
        return 0.0
    else:
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
    
    if len(anchor_variables) == 0:
        raise ValueError(f"No anchor variables found in dataset for columns: {anchor_cols}")
    
    # One-hot encode anchor variables
    anchor_variables = np.asarray(anchor_variables).T
    Z = OneHotEncoder(sparse_output=False, dtype=np.float32).fit_transform(anchor_variables)
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

    # Only consider genes present in the inferred GRN
    gene_mask = np.logical_or(np.any(A, axis=1), np.any(A, axis=0))
    X = X[:, gene_mask]
    A = A[gene_mask, :][:, gene_mask]
    gene_names = gene_names[gene_mask]
    
    print(f"Evaluating {X.shape[1]} genes with {np.sum(A != 0)} regulatory links")

    # Center and scale dataset
    scaler = StandardScaler()
    scaler.fit(X)
    X = scaler.transform(X)

    # Evaluate stability for each target gene
    scores = []
    for j in tqdm.tqdm(range(X.shape[1]), desc="Evaluating gene stability"):
        score = evaluate_gene_stability(X, Z, A, j)
        scores.append(score)

    # Calculate final score
    final_score = float(np.mean(scores))
    print(f"Method: {method_id}")
    print(f"Anchor Regression Score: {final_score:.6f}")

    # Return results as DataFrame
    results = {
        'anchor_regression': [final_score]
    }

    df_results = pd.DataFrame(results)
    return df_results