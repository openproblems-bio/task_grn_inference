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


MAX_N_ITER = 2000


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
    print(list(f["layers"]))
    #X = f["layers"]["X_norm"]
    X = f["layers"]["lognorm"]
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


def neumann_series(A: torch.Tensor, k: int = 2) -> torch.Tensor:
    """Approximate the inverse of I - A using Neumann series.

    Args:
        A: the matrix for which to invert I - A.
        k: the number of terms in the series. The higher, the more accurate.

    Returns:
        Approximated inverse of I - A.
    """
    B = torch.eye(A.shape[0], device=A.device, dtype=A.dtype)
    for k in range(k):
        B = B + torch.mm(B, A)
    return B


def solve_sem(A: torch.Tensor, delta: torch.Tensor, stable: bool = False) -> torch.Tensor:
    """Compute the perturbation using a linear structural equation model (SEM).

    The SEM is defined as: perturbation^T = ((I - A)^{-1})^T delta^T.

    Args:
         A: the matrix for which to invert I - A.
            A_{ij} is the directed regulatory link between genes i and j.
         delta: the exogenous shocks. delta_{kj} is the impact of the drugs
            administrated to sample k on gene j.
        stable: Whether to use Neumann series for more stable gradients.

    Returns:
        The estimated perturbations. A matrix of same shape as `delta`.
    """
    if stable:  # Neumann series for more stable gradients
        B = neumann_series(A, k=2)
        return torch.mm(B.mT, delta.mT).mT
    else:  # Explicit matrix inversion
        I = torch.eye(A.size(1), dtype=A.dtype)
        return torch.linalg.solve((I - A).mT, delta.mT).mT


def regression_based_grn_inference(X: np.ndarray, base: np.ndarray, signed: bool = False) -> np.ndarray:
    """Infer GRN weights given a GRN topology.

    Args:
        X: gene expression from control samples. A matrix of shape `(n_samples, n_genes)`.
        base: an inferred GRN matrix. A matrix of shape `(n_genes, n_genes)`.
            Because `base` may contain values not necessarily suitable for regression,
            such as -log(p-values) or arbitrary scores, we re-fit the GRN weights
            for all the non-zero elements in `base`.
        signed: Whether to force the output regulatory links to have the same sign as
            the input regulatory links.

    Returns:
        A new GRN matrix of shape `(n_genes, n_genes)`.
    """
    mask = (base != 0)
    n_genes = X.shape[1]
    A = np.zeros((n_genes, n_genes), dtype=X.dtype)
    for j in tqdm.tqdm(range(n_genes)):
        model = ElasticNet(alpha=0.001, fit_intercept=False)
        if not np.any(mask[:, j]):
            continue  # If gene has no regulator, skip it.
        model.fit(X[:, mask[:, j]], X[:, j])
        A[mask[:, j], j] = model.coef_

    # For each parameter that has the wrong regulatory sign, set it to zero.
    if signed:
        A[base * A < 0] = 0

    return A


def evaluate_grn(
        X_controls: np.ndarray,
        delta_X: np.ndarray,
        is_train: np.ndarray,
        is_reporter: np.ndarray,
        A: np.ndarray,
        signed: bool = True
) -> np.ndarray:

    n_iter = MAX_N_ITER
    learning_rate = 0.0005

    signs = np.sign(A)
    mask = np.asarray(A != 0)
    print(f"Number of GRN edges: {np.sum(mask)}")

    # Learn initial GRN weights from the first set
    print("Infer GRN weights from training set, given the provided GRN topology")
    A = regression_based_grn_inference(X_controls, A, signed=signed)

    # Learn shocks from perturbations, using reporter genes only
    print("Learn shocks from perturbations, using reporter genes only")
    reporter_idx = np.where(is_reporter)[0]
    regulator_idx = np.where(mask.any(axis=1))[0]
    F = np.linalg.inv(np.eye(A.shape[0]) - A)
    F_CR = F[np.ix_(regulator_idx, reporter_idx)]
    Y_R = delta_X[:, reporter_idx]
    lam = 0.01
    I = np.eye(len(regulator_idx), dtype=delta_X.dtype)
    M = F_CR @ F_CR.T + lam * I
    delta = np.zeros_like(delta_X)
    delta[:, regulator_idx] = Y_R @ F_CR.T @ np.linalg.inv(M)

    # Learn perturbations from shocks, using evaluation genes only
    print("Learn perturbations from shocks, using evaluation genes only")
    delta_X_train = torch.from_numpy(delta_X[is_train, :])
    delta_train = torch.from_numpy(delta[is_train, :])
    signs = torch.from_numpy(signs).to(torch.float)
    mask = torch.from_numpy(mask).to(torch.bool)
    A = torch.nn.Parameter(torch.from_numpy(A))
    A_eff = torch.abs(A) * signs if signed else A * mask
    optimizer = torch.optim.Adam([A], lr=learning_rate)
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
        optimizer, mode='min', patience=5,
        min_lr=1e-5, cooldown=3, factor=0.8
    )
    best_loss = np.inf
    best_A_eff = A_eff.detach()
    pbar = tqdm.tqdm(range(n_iter))
    X_non_reporter = delta_X_train[:, ~is_reporter]
    for _ in pbar:
        optimizer.zero_grad()
        A_eff = torch.abs(A) * signs if signed else A * mask
        delta_X_hat = solve_sem(A_eff, delta_train)
        loss = torch.mean(torch.sum(torch.square(X_non_reporter - delta_X_hat[:, ~is_reporter]), dim=1))
        loss = loss + 0.00001 * torch.sum(torch.abs(A))
        loss = loss + 0.00001 * torch.sum(torch.square(A))

        # Keep track of best solution
        if loss.item() < best_loss:
            best_loss = loss.item()
            best_A_eff = A_eff.detach()

        loss.backward()
        optimizer.step()
        scheduler.step(loss.item())
        pbar.set_description(str(loss.item()))
    A = best_A_eff
    mask = mask.detach().cpu().numpy().astype(bool)

    # Predict perturbations in test set
    delta_test = torch.from_numpy(delta[~is_train, :])
    delta_X_test = torch.from_numpy(delta_X[~is_train, :])
    delta_X_hat = solve_sem(A, delta_test)

    # Evaluate predictions (using GRN) on the evaluation genes
    delta_X_test = delta_X_test.detach().cpu().numpy()
    delta_X_hat = delta_X_hat.detach().cpu().numpy()
    has_parent = (mask.any(axis=0))
    eval_mask = ((~is_reporter) & has_parent)
    coefficients = []
    for j in range(len(eval_mask)):
        if eval_mask[j]:
            if np.all(delta_X_test[:, j] == delta_X_test[0, j]):  # Constant
                coefficients.append(0.0)
            else:
                coefficients.append(spearmanr(delta_X_test[:, j], delta_X_hat[:, j]).correlation)
        else:
            coefficients.append(0.0)
    return np.nan_to_num(coefficients, nan=0)


# Create a split between training and test sets.
# Make sure that no compound ends up in both sets.
try:
    splitter = GroupShuffleSplit(test_size=0.5, n_splits=2, random_state=0xCAFE)
    train_idx, _ = next(splitter.split(delta_X, groups=cv_groups))
except ValueError:
    splitter = KFold(n_splits=2, random_state=0xCAFE, shuffle=True)
    train_idx, _ = next(splitter.split(delta_X))
is_train = np.zeros(len(delta_X), dtype=bool)
is_train[train_idx] = True

# Create a split between genes: reporter genes and evaluation genes.
# All TFs should be included in the reporter gene set.
n_genes = A.shape[1]
reg_mask = np.asarray(A != 0).any(axis=1)  # TF mask
is_reporter = np.copy(reg_mask)
print(f"Proportion of reporter genes: {np.mean(is_reporter)}")
print(f"Use regulatory modes/signs: {use_signs}")

# Create a symmetric (causally-wrong) baseline GRN
print(f"Creating baseline GRN")
A_baseline = np.copy(A).T
np.random.shuffle(A_baseline)
A_baseline = A_baseline.T

# Evaluate inferred GRN
print("\n======== Evaluate inferred GRN ========")
scores = evaluate_grn(X_controls, delta_X, is_train, is_reporter, A, signed=use_signs)
print(f"Final score: {np.mean(scores)}")

# Evaluate baseline GRN
print("\n======== Evaluate shuffled GRN ========")
n_repeats = 5
scores_baseline = np.zeros_like(scores)
for _ in range(n_repeats):  # Repeat for more robust estimation
    scores_baseline += evaluate_grn(X_controls, delta_X, is_train, is_reporter, A_baseline, signed=use_signs)
scores_baseline /= n_repeats

# Perform T-test between inferred GRN and baseline GRN
print(f"Average Spearman: {np.mean(scores)}")
print(f"Average Spearman (shuffled): {np.mean(scores_baseline)}")
res = wilcoxon(scores - scores_baseline, zero_method='wilcox', alternative='greater')
print(f"Wilcoxon signed-rank test: pvalue={res.pvalue}")

# Compute final score
steepness = 1.5
f = lambda p: (-np.log(p)) ** steepness
score = f(res.pvalue) / (f(res.pvalue) + f(1e-10))
print(f"Final score: {score}")


"""
Average Spearman: 0.20102521156254408
Average Spearman (shuffled): 0.10157913679830942
Wilcoxon signed-rank test: pvalue=7.65004668751997e-21
Final score: 0.7404715122547103
"""
