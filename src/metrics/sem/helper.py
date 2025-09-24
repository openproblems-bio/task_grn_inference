

import os
import traceback
import random
import h5py
import numpy as np
import pandas as pd
import torch
import anndata as ad
import tqdm
from scipy.sparse import csr_matrix
from scipy.stats import ttest_rel, spearmanr, wilcoxon
from sklearn.linear_model import ElasticNet
from sklearn.preprocessing import StandardScaler, LabelEncoder
import sys
import warnings
warnings.filterwarnings("ignore", category=UserWarning)
from sklearn.model_selection import GroupShuffleSplit, KFold
os.environ['CUBLAS_WORKSPACE_CONFIG'] = ':4096:8' # For reproducibility purposes (on GPU)
DEVICE = 'cuda' if torch.cuda.is_available() else 'cpu'
NUMPY_DTYPE = np.float32

# Hyper-parameters
MAX_N_ITER = 2000

# For reproducibility purposes
seed = 0xCAFE
os.environ['PYTHONHASHSEED'] = str(seed)
random.seed(seed)
np.random.seed(seed)
torch.manual_seed(seed)
torch.cuda.manual_seed(seed)
torch.cuda.manual_seed_all(seed)
torch.use_deterministic_algorithms(True)

from util import read_prediction, manage_layer



def encode_obs_cols(adata, cols):
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
def compute_perturbations(X, are_controls, groups: np.array) -> np.ndarray:
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
        if group_id not in control_map:
            raise ValueError(f"No control found for sample {i} (group {group_id}).")
            
    
        j = control_map[group_id]
        delta_X[i, :] -= delta_X[j, :]
    return delta_X

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
    # A = spearman_based_grn_inference(X_controls, A, signed=signed)

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
        # pbar.set_description(str(loss.item()))
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

DATASET_GROUPS = {
    "op": {
        "match": ["plate_name", "donor_id", "cell_type", 'well'],
        "loose_match": ["donor_id", "cell_type", "plate_name"],
        "cv": ["perturbation", "cell_type"],
    },
    "parsebioscience": {
        "match": ["donor_id", "cell_type", "well"],
        "loose_match": ["donor_id", "cell_type"],
        "cv": ["perturbation", "cell_type"],
    },
    "300BCG": {
        "match": ["donor_id",  "cell_type"],
        "loose_match": ["cell_type"],
        "cv": ["perturbation", "cell_type"],
    },
}

def main(par):
    dataset_id = ad.read_h5ad(par['evaluation_data'], backed='r').uns['dataset_id']
    method_id = ad.read_h5ad(par['prediction'], backed='r').uns['method_id']

    par['cv_groups'] = DATASET_GROUPS[dataset_id]['cv']
    par['match'] = DATASET_GROUPS[dataset_id]['match']
    par['loose_match'] = DATASET_GROUPS[dataset_id]['loose_match']

    adata = ad.read_h5ad(par['evaluation_data'])
    layer = manage_layer(adata, par)
    X = adata.layers[layer]
    gene_names = adata.var_names
    assert 'is_control' in adata.obs.columns, "Evaluation dataset must contain 'is_control'."
    are_controls = adata.obs['is_control'].values.astype(bool)
    cv_groups = encode_obs_cols(adata, par['cv_groups'])
    match_groups = encode_obs_cols(adata, par['match'])
    loose_match_groups = encode_obs_cols(adata, par['loose_match'])
    cv_groups = combine_multi_index(*cv_groups)
    match_groups = combine_multi_index(*match_groups)
    loose_match_groups = combine_multi_index(*loose_match_groups)
    gene_dict = {gene_name: i for i, gene_name in enumerate(gene_names)}

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

    # Only consider the genes that are actually present in the inferred GRN
    gene_mask = np.logical_or(np.any(A, axis=1), np.any(A, axis=0))
    if True:
        idx = np.argsort(np.maximum(np.sum(A!=0, axis=1), np.sum(A!=0, axis=0)))[:-5000]
        gene_mask[idx] = False
        
    X = X[:, gene_mask]
    X = X.toarray() if isinstance(X, csr_matrix) else X
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
    try:
        delta_X = compute_perturbations(X, are_controls, match_groups)
    except:
        print("Error in compute_perturbations with strict matching. Using loose matching instead.")
        traceback.print_exc()
        delta_X = compute_perturbations(X, are_controls, loose_match_groups)

    # Remove negative controls from downstream analysis
    delta_X = delta_X[~are_controls, :]
    cv_groups = cv_groups[~are_controls]
    match_groups = match_groups[~are_controls]
    loose_match_groups = loose_match_groups[~are_controls]

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
    # If the reporter gene set does not represent at least 50% of all genes,
    # then we randomly add target genes until the 50% threshold is reached.
    n_genes = A.shape[1]
    reg_mask = np.asarray(A != 0).any(axis=1)  # TF mask
    is_reporter = np.copy(reg_mask)
    print(f"Proportion of reporter genes: {np.mean(is_reporter)}")
    print(f"Use regulatory modes/signs: {use_signs}")

    # Create a symmetric (causally-wrong) baseline GRN
    print(f"Creating baseline GRN")
    mask = np.abs(A) > np.abs(A.T)
    #A_baseline = mask * A + (~mask) * A.T
    A_baseline = np.copy(A).T
    np.random.shuffle(A_baseline)
    A_baseline = A_baseline.T

    # Evaluate inferred GRN
    print("\n======== Evaluate inferred GRN ========")
    scores = evaluate_grn(X_controls, delta_X, is_train, is_reporter, A, signed=use_signs)
    # Evaluate baseline GRN
    print("\n======== Evaluate shuffled GRN ========")
    scores_baseline = evaluate_grn(X_controls, delta_X, is_train, is_reporter, A_baseline, signed=use_signs)

    # Keep only the genes for which both GRNs got a score
    #mask = ~np.logical_or(np.isnan(scores), np.isnan(scores_baseline))
    #scores = scores[mask]
    #scores_baseline = scores_baseline[mask]

    rr_all = {}
    # Perform rank test between actual scores and baseline
    rr_all['spearman'] = float(np.mean(scores))
    rr_all['spearman_shuffled'] = float(np.mean(scores_baseline))
    if np.std(scores - scores_baseline) == 0:
        df_results = pd.DataFrame({'sem': [0.0]})
    else:
        res = wilcoxon(scores - scores_baseline, zero_method='wilcox', alternative='greater')
        rr_all['Wilcoxon pvalue'] = float(res.pvalue)

        print(rr_all)
        
        eps = 1e-300  # very small number to avoid log(0)
        pval_clipped = max(res.pvalue, eps)
        # Compute final score
        if False:
            steepness = 1.5
            f = lambda p: (-np.log(p)) ** steepness
            score = f(pval_clipped) / (f(pval_clipped) + f(1e-10))
        else:
            effect = np.mean(scores) - np.mean(scores_baseline)
            score = effect * (-np.log10(pval_clipped))
        print(f"Final score: {score}")

        results = {
            'sem': [float(score)]
        }

        df_results = pd.DataFrame(results)
    return df_results