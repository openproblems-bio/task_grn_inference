import os
import traceback

from sklearn.model_selection import GroupShuffleSplit, KFold

os.environ['CUBLAS_WORKSPACE_CONFIG'] = ':4096:8' # For reproducibility purposes (on GPU)

import random
import h5py
import numpy as np
import pandas as pd
import torch
import tqdm
from scipy.sparse import csr_matrix
from scipy.stats import ttest_rel, spearmanr, wilcoxon
from sklearn.linear_model import ElasticNet
from sklearn.preprocessing import StandardScaler, LabelEncoder


# Hyper-parameters
DEVICE = 'cuda' if torch.cuda.is_available() else 'cpu'
NUMPY_DTYPE = np.float32

# For reproducibility purposes
seed = 0xCAFE
os.environ['PYTHONHASHSEED'] = str(seed)
random.seed(seed)
np.random.seed(seed)
torch.manual_seed(seed)
torch.cuda.manual_seed(seed)
torch.cuda.manual_seed_all(seed)
torch.use_deterministic_algorithms(True)


def combine_multi_index(*arrays) -> np.array:
    """Combine parallel label arrays into a single integer label per position."""
    A = np.stack(arrays)
    n_classes = tuple(int(A[i].max()) + 1 for i in range(A.shape[0]))
    return np.ravel_multi_index(A, dims=n_classes, order='C')


# Load perturbation data
with h5py.File("../../../resources_test/grn_benchmark/evaluation_data/op_bulk.h5ad", "r") as f:

    # Get sample info
    print(f["obs"].keys())
    cv_groups, match_groups, exact_match_groups = [], [], []
    for obs_name in ["perturbation", "perturbation_type", "cell_type", "donor_id", "plate_name", "row"]:
        if obs_name in f["obs"]:
            if "codes" in f["obs"][obs_name]:
                codes = f["obs"][obs_name]["codes"][:]
            else:
                codes = LabelEncoder().fit_transform(f["obs"][obs_name][:])
            if obs_name in {"perturbation", "cell_type"}:
                cv_groups.append(codes)
            if obs_name not in {"row", "perturbation"}:
                match_groups.append(codes)
            if obs_name != "perturbation":
                exact_match_groups.append(codes)

    # Groups used for cross-validation
    cv_groups = combine_multi_index(*cv_groups)

    # Groups used for matching with negative controls
    match_groups = combine_multi_index(*match_groups)

    # Groups used for exact matching with negative controls
    # (e.g., samples should be from the same plate row)
    exact_match_groups = combine_multi_index(*exact_match_groups)

    are_controls = f["obs"]["is_control"][:].astype(bool)
    X = f["layers"]["X_norm"]
    if isinstance(X, h5py.Dataset):  # Dense array
        X = X[:].astype(np.float32)
    else:  # Sparse array
        group = X
        data = group['data'][...].astype(np.float32)
        indices = group['indices'][...]
        indptr = group['indptr'][...]
        shape = group.attrs['shape']
        sparse_matrix = csr_matrix((data, indices, indptr), shape=shape)
        X = sparse_matrix.toarray()
        print(X.shape)

    # Get gene names
    if "index" in f["var"]:
        gene_names = f["var"]["index"][:].astype(str)
    elif "Probe_id" in f["var"]:
        gene_names = f["var"]["Probe_id"][:].astype(str)
    else:
        gene_names = f["var"]["gene_name"][:].astype(str)
    print(gene_names)

gene_dict = {gene_name: i for i, gene_name in enumerate(gene_names)}

# Load inferred GRN
df = pd.read_csv("../../../output/benchmark/grn_models/scenicplus.csv")
sources = df["source"].to_numpy()
targets = df["target"].to_numpy()
weights = df["weight"].to_numpy()
#with h5py.File("../../../resources_test/grn_models/op/collectri.h5ad", "r") as f:
#    sources = f["uns"]["prediction"]["source"][:].astype(str)
#    targets = f["uns"]["prediction"]["target"][:].astype(str)
#    weights = f["uns"]["prediction"]["weight"][:]

A = np.zeros((len(gene_names), len(gene_names)), dtype=X.dtype)
for source, target, weight in zip(sources, targets, weights):
    if (source in gene_dict) and (target in gene_dict):
        i = gene_dict[source]
        j = gene_dict[target]
        A[i, j] = float(weight)

# Only consider the genes that are actually present in the inferred GRN
gene_mask = np.logical_or(np.any(A, axis=1), np.any(A, axis=0))
X = X[:, gene_mask]
A = A[gene_mask, :][:, gene_mask]
gene_names = gene_names[gene_mask]

# Check whether the inferred GRN contains signed predictions
use_signs = np.any(A < 0)

# Center and scale dataset
scaler = StandardScaler()
scaler.fit(X[are_controls, :])  # Use controls only to infer statistics (to avoid data leakage)
X = scaler.transform(X)

# Get negative controls
X_controls = X[are_controls, :]


def compute_perturbations(groups: np.array) -> np.ndarray:
    # Compute perturbations as differences between samples and their matched controls.
    # By matched control, we mean a negative control from the same donor, located on
    # the same plate and from the same cell type.
    # By construction, perturbations will be 0 for control samples.
    delta_X = np.copy(X)
    control_map = {}
    for i, (is_control, group_id) in enumerate(zip(are_controls, groups)):
        if is_control:
            control_map[group_id] = i
    for i, group_id in enumerate(groups):
        j = control_map[group_id]
        delta_X[i, :] -= delta_X[j, :]
    return delta_X


try:  # Try exact control matching first (same plate row, if info is available)
    delta_X = compute_perturbations(exact_match_groups)
except:  # If samples are missing, use less stringent control matching instead
    print(traceback.format_exc())
    delta_X = compute_perturbations(match_groups)


# Remove negative controls from downstream analysis
delta_X = delta_X[~are_controls, :]
cv_groups = cv_groups[~are_controls]
match_groups = match_groups[~are_controls]
exact_match_groups = exact_match_groups[~are_controls]


def anchor_regression(
        X: np.ndarray,
        Z: np.ndarray,
        Y: np.ndarray,
        l2_reg: float = 1e-3,
        anchor_strength: float = 5.0
) -> np.ndarray:
    """Anchor regression.

    Args:
        X: input features, of shape (n, d).
        Z: anchor variables. The model is required to be invariant to
            these environmental variables. Shape (n, z).
        Y: predicted variables, of shape (n, u).

    Returns:
        Inferred parameters, of shape (d, u).
    """

    # Projector onto the anchor subspace
    pi = Z @ np.linalg.pinv(Z.T @ Z) @ Z.T

    # Whitening
    W = np.eye(len(Z)) + (np.sqrt(1.0 + anchor_strength) - 1.0) * pi
    X_t = W @ X
    Y_t = W @ Y

    # Ridge regression on the whitened data
    theta = np.linalg.solve(X_t.T @ X_t + l2_reg * np.eye(X_t.shape[1]), X_t.T @ Y_t)

    #Y_hat = X @ Theta
    return theta



