import os
import traceback
import random
import h5py
import numpy as np
import pandas as pd
import torch
import anndata as ad
import scanpy as sc
import tqdm
from scipy.sparse import csr_matrix
from scipy.stats import ttest_rel, spearmanr, pearsonr, wilcoxon
from sklearn.linear_model import ElasticNet
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.metrics import r2_score
import sys
import warnings
warnings.filterwarnings("ignore", category=UserWarning)
from sklearn.model_selection import GroupShuffleSplit, KFold


os.environ['CUBLAS_WORKSPACE_CONFIG'] = ':4096:8' # For reproducibility purposes (on GPU)
# Additional environment variables for determinism
os.environ['TF_DETERMINISTIC_OPS'] = '1'
os.environ['PYTHONHASHSEED'] = '0'

# Force NumPy to use single thread for determinism
os.environ['OMP_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1'
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['NUMEXPR_NUM_THREADS'] = '1'

DEVICE = 'cuda' if torch.cuda.is_available() else 'cpu'
NUMPY_DTYPE = np.float32

# Hyper-parameters
MAX_N_ITER = 500

# For reproducibility purposes - use consistent seed
seed = 42  # Changed to standard seed
os.environ['PYTHONHASHSEED'] = str(seed)
random.seed(seed)
np.random.seed(seed)
torch.manual_seed(seed)
torch.cuda.manual_seed(seed)
torch.cuda.manual_seed_all(seed)
torch.use_deterministic_algorithms(True, warn_only=True)  # warn_only to avoid crashes
torch.backends.cudnn.deterministic = True
torch.backends.cudnn.benchmark = False


from util import read_prediction, manage_layer, create_grn_baseline
from config import DATASET_GROUPS


# List of immediate early genes (IEGs) - scCustomize v3.2.0
IEG = [
    "CCL7"  , "CCL2"    , "TNFSF9"  , "SIK1B"   , "RASD1"  ,  "NUP98",
    "CEBPD" , "NFKBID"  , "NOCT"    , "FOS"     , "IER2"   ,  "HBEGF",
    "RHEB"  , "PLAU"    , "IFNB1"   , "PCDH8"   , "WEE1"   ,  "FBXO33",
    "PMAIP1", "DUSP1"   , "PLK2"    , "TSC22D1" , "MAP3K8" ,  "PIAS1",
    "KLF10" , "BDNF"    , "CCN2"    , "TRIB1"   , "SOD2"   ,  "IER3",
    "PLAT"  , "RCAN1"   , "ZFP36"   , "CCL5"    , "NFKBIA" ,  "NRN1",
    "KLF2"  , "SERPINE1", "MAFF"    , "NCOA7"   , "GDF15"  ,  "LDLR",
    "TNF"   , "CCRL2"   , "CCL18"   , "CCL3"    , "FOSL1"  ,  "DUSP2",
    "INHBA" , "JUND"    , "NR4A1"   , "EGR3"    , "IRF1"   ,  "IFIT3",
    "IFIT1B", "NR4A3"   , "PER2"    , "GADD45G" , "SOCS3"  ,  "TLR2",
    "PELI1" , "IL1A"    , "RGS2"    , "CSF2"    , "F3"     ,  "APOLD1",
    "ACKR4" , "BHLHE40" , "IL23A"   , "GBP2"    , "IL1B"   ,  "ARIH1",
    "GBP2"  , "JUN"     , "NPTX2"   , "FOSB"    , "EGR1"   ,  "MBNL2",
    "VCAM1" , "ZFP36L2" , "ARF4"    , "MMP13"   , "JUNB"   ,  "ZFP36L1",
    "CCN1"  , "NPAS4"   , "PPP1R15A", "EGR4"    , "NFKBIZ" ,  "ARC",
    "MCL1"  , "RGS1"    , "KLF6"    , "CLEC4E"  , "SGK1"   ,  "ARHGEF3",
    "EGR2"  , "IFIT2"   , "ID2"     , "DUSP6"   , "SIK1"   ,  "DUSP5",
    "NFIB"  , "SAA2"    , "SAA1"    , "THBS1"   , "FOSL2"  ,  "ICAM1",
    "CXCL10", "CSRNP1"  , "BCL3"    , "MARCKSL1", "HOMER1" ,  "TRAF1",
    "ATF3"  , "FLG"     , "SRF"     , "PIM1"    , "GADD45B",  "HES1",
    "GEM"   , "TNFAIP3" , "PER1"    , "CREM"    , "CD69"   ,  "IL6",
    "BTG2"  , "ACOD1"   , "CEBPB"   , "CXCL11"  , "IL12B"  ,  "NR4A2",
    "PTGS2" , "IKBKE"   , "TXNIP"   , "CD83"    , "IER5"   ,  "IL10",
    "MYC"   , "CXCL12"  , "SLC2A3"  , "CXCL1"
]


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

import numpy as np

def compute_perturbations(X, are_controls, groups: np.ndarray, loose_match_groups: np.ndarray) -> np.ndarray:
    """
    Compute perturbations as differences between samples and their matched controls.
    Each sample is compared to its matched control (from the same donor, plate, etc.).
    If no matching control is found (even for loose grouping), assign any available control.

    Parameters
    ----------
    X : np.ndarray
        Expression or feature matrix (samples × features)
    are_controls : np.ndarray (bool)
        Boolean array indicating which samples are controls
    groups : np.ndarray
        Strict matching groups (e.g., donor+plate+celltype)
    loose_match_groups : np.ndarray
        Looser matching groups (e.g., donor+plate)
    """

    delta_X = np.copy(X)
    control_map, loose_control_map = {}, {}

    # Build maps for strict and loose control matching
    for i, (is_control, group_id) in enumerate(zip(are_controls, groups)):
        if is_control:
            control_map[group_id] = i
    for i, (is_control, group_id) in enumerate(zip(are_controls, loose_match_groups)):
        if is_control:
            loose_control_map[group_id] = i

    # Identify any available control indices (for fallback)
    all_controls = np.where(are_controls)[0]
    if len(all_controls) == 0:
        raise RuntimeError("No control samples found at all — cannot compute perturbations.")

    # Compute perturbations
    for i, (group_id, loose_group_id) in enumerate(zip(groups, loose_match_groups)):
        try:
            if group_id in control_map:
                j = control_map[group_id]
            elif loose_group_id in loose_control_map:
                j = loose_control_map[loose_group_id]
            else:
                # Fall back to any available control (e.g., first one)
                j = all_controls[0]

            delta_X[i, :] -= X[j, :]

        except Exception as e:
            raise RuntimeError(
                f"Failed to compute perturbation for sample index {i}:\n"
                f"  group_id = {group_id}\n"
                f"  loose_group_id = {loose_group_id}\n"
                f"  Are controls found for this group? "
                f"Strict = {group_id in control_map}, Loose = {loose_group_id in loose_control_map}\n"
                f"  Fallback control used = {len(all_controls) > 0}\n"
                f"  Total controls available: {sum(are_controls)} / {len(are_controls)}\n"
                f"Hint: Loose grouping may be inconsistent; used fallback control."
            ) from e

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
    I = torch.eye(A.shape[0], device=A.device, dtype=A.dtype)
    term = I.clone()
    B = I.clone()
    for _ in range(k):
        term = term @ A
        B = B + term
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
        model = ElasticNet(alpha=0.001, fit_intercept=False, random_state=seed)  # Use global seed
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
    learning_rate = 0.001

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
    lam = 0.001
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
    
    # Set manual seed again before optimizer initialization for determinism
    torch.manual_seed(seed)
    np.random.seed(seed)
    random.seed(seed)
    
    optimizer = torch.optim.Adam([A], lr=learning_rate)
    # Use StepLR instead of ReduceLROnPlateau for determinism
    scheduler = torch.optim.lr_scheduler.StepLR(
        optimizer, step_size=20, gamma=0.8
    )
    best_loss = np.inf
    best_A_eff = A_eff.detach()
    best_iteration = 0
    pbar = tqdm.tqdm(range(n_iter))
    X_non_reporter = delta_X_train[:, ~is_reporter]
    for iteration in pbar:
        optimizer.zero_grad()
        A_eff = torch.abs(A) * signs if signed else A * mask
        delta_X_hat = solve_sem(A_eff, delta_train)
        loss = torch.mean(torch.sum(torch.square(X_non_reporter - delta_X_hat[:, ~is_reporter]), dim=1))
        loss = loss + 0.001 * torch.sum(torch.abs(A))
        loss = loss + 0.001 * torch.sum(torch.square(A))
        pbar.set_description(str(loss.item()))

        # Keep track of best solution
        if loss.item() < best_loss:
            best_loss = loss.item()
            best_A_eff = A_eff.detach()
            best_iteration = iteration

        loss.backward()
        optimizer.step()
        scheduler.step()  # Changed from scheduler.step(loss.item())
        # pbar.set_description(str(loss.item()))
    A = best_A_eff
    mask = mask.detach().cpu().numpy().astype(bool)
    print(f"Best iteration: {best_iteration + 1} / {n_iter}")

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
                coefficients.append(np.nan)
            else:
                coefficients.append(max(0, r2_score(delta_X_test[:, j], delta_X_hat[:, j])))
        else:
            coefficients.append(np.nan)
    coefficients = np.array(coefficients)
    #return np.nan_to_num(coefficients, nan=0)
    return coefficients
    


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
    assert len(cv_groups) > 0, "No cv_groups could be encoded."
    assert len(match_groups) > 0, "No match_groups could be encoded."
    assert len(loose_match_groups) > 0, "No loose_match_groups could be encoded."
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

    # Compute HVGs from full evaluation data (for HVG-based evaluation)
    print("\n======== Computing HVGs from full evaluation data ========")
    n_top_hvg = par['n_top_genes']
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_hvg, flavor='seurat', layer=layer)
    hvg_mask_full = adata.var['highly_variable'].values
    hvg_genes = gene_names[hvg_mask_full]
    print(f"Total HVGs identified: {hvg_mask_full.sum()}")

    # For GRN-based evaluation: keep only most-connected genes in the GRN
    print("\n======== Filtering genes for GRN-based evaluation ========")
    gene_mask_grn = np.logical_or(np.any(A, axis=1), np.any(A, axis=0))
    in_degrees = np.sum(A != 0, axis=0)
    out_degrees = np.sum(A != 0, axis=1)
    n_genes_grn = par['n_top_genes']
    idx = np.argsort(np.maximum(out_degrees, in_degrees))[:-n_genes_grn]
    gene_mask_grn[idx] = False
    
    X_grn = X[:, gene_mask_grn]
    X_grn = X_grn.toarray() if isinstance(X_grn, csr_matrix) else X_grn
    A_grn = A[gene_mask_grn, :][:, gene_mask_grn]
    gene_names_grn = gene_names[gene_mask_grn]
    print(f"Genes for GRN-based evaluation: {len(gene_names_grn)}")

    # Remove self-regulations
    np.fill_diagonal(A, 0)
    np.fill_diagonal(A_grn, 0)

    # Check whether the inferred GRN contains signed predictions
    if False:
        use_signs = np.any(A < 0)
    else:
        use_signs = False

    # Center and scale dataset for GRN-based evaluation
    scaler_grn = StandardScaler()
    X_grn_controls = X_grn[are_controls, :]
    scaler_grn.fit(X_grn_controls)
    X_grn_scaled = scaler_grn.transform(X_grn)
    X_grn_controls_scaled = X_grn_scaled[are_controls, :]
    delta_X_grn = compute_perturbations(X_grn_scaled, are_controls, match_groups, loose_match_groups)
    delta_X_grn = delta_X_grn[~are_controls, :]

    # Center and scale dataset for HVG-based evaluation (use all HVG genes, even if not in GRN)
    X_hvg = X[:, hvg_mask_full]
    X_hvg = X_hvg.toarray() if isinstance(X_hvg, csr_matrix) else X_hvg
    A_hvg = A[hvg_mask_full, :][:, hvg_mask_full]
    gene_names_hvg = gene_names[hvg_mask_full]
    scaler_hvg = StandardScaler()
    X_hvg_controls = X_hvg[are_controls, :]
    scaler_hvg.fit(X_hvg_controls)
    X_hvg_scaled = scaler_hvg.transform(X_hvg)
    X_hvg_controls_scaled = X_hvg_scaled[are_controls, :]
    delta_X_hvg = compute_perturbations(X_hvg_scaled, are_controls, match_groups, loose_match_groups)
    delta_X_hvg = delta_X_hvg[~are_controls, :]
    print(f"Genes for HVG-based evaluation: {len(gene_names_hvg)}")

    # Remove negative controls from downstream analysis
    cv_groups = cv_groups[~are_controls]
    match_groups = match_groups[~are_controls]
    loose_match_groups = loose_match_groups[~are_controls]

    # Create a split between training and test sets.
    # Make sure that no compound ends up in both sets.
    try:
        splitter = GroupShuffleSplit(test_size=0.5, n_splits=2, random_state=seed)  # Use consistent seed
        train_idx, _ = next(splitter.split(delta_X_grn, groups=cv_groups))
    except ValueError:
        print("Group k-fold failed. Using k-fold CV instead.")
        splitter = KFold(n_splits=2, random_state=seed, shuffle=True)  # Use consistent seed
        train_idx, _ = next(splitter.split(delta_X_grn))
    is_train = np.zeros(len(delta_X_grn), dtype=bool)
    is_train[train_idx] = True

    # ========== GRN-based evaluation ==========
    print("\n======== Evaluate inferred GRN (GRN-based: most connected genes) ========")
    n_genes_grn = A_grn.shape[1]
    reg_mask_grn = np.asarray(A_grn != 0).any(axis=1)
    ieg_mask_grn = np.asarray([gene_name in IEG for gene_name in gene_names_grn], dtype=bool)
    is_reporter_grn = np.logical_or(reg_mask_grn, ieg_mask_grn)
    print(f"Proportion of reporter genes (GRN): {np.mean(is_reporter_grn)}")


    scores_grn = evaluate_grn(X_grn_controls_scaled, delta_X_grn, is_train, is_reporter_grn, A_grn, signed=use_signs)
    valid_scores_grn = scores_grn[~np.isnan(scores_grn)]

    if len(valid_scores_grn) == 0:
        print("WARNING: No valid genes to evaluate for GRN-based!")
        sem_grn_score = 0.0
    else:
        sem_grn_score = float(np.mean(valid_scores_grn))
        print(f"SEM GRN score (mean R²): {sem_grn_score:.4f}")
        print(f"Valid genes evaluated: {len(valid_scores_grn)}/{len(scores_grn)}")
        print(f"SEM GRN score (min): {np.min(valid_scores_grn):.4f}")
        print(f"SEM GRN score (max): {np.max(valid_scores_grn):.4f}")

    # ========== HVG-based evaluation ==========
    print("\n======== Evaluate inferred GRN (HVG-based: highly variable genes) ========")
    n_genes_hvg = A_hvg.shape[1]
    reg_mask_hvg = np.asarray(A_hvg != 0).any(axis=1)
    ieg_mask_hvg = np.asarray([gene_name in IEG for gene_name in gene_names_hvg], dtype=bool)
    is_reporter_hvg = np.logical_or(reg_mask_hvg, ieg_mask_hvg)
    print(f"Proportion of reporter genes (HVG): {np.mean(is_reporter_hvg)}")

    scores_hvg = evaluate_grn(X_hvg_controls_scaled, delta_X_hvg, is_train, is_reporter_hvg, A_hvg, signed=use_signs)
    
    # For HVGs: genes with no GRN connections get score of 0 (penalize missing connections)
    has_parent_hvg = (np.asarray(A_hvg != 0).any(axis=0))
    eval_mask_hvg = ~is_reporter_hvg
    scores_hvg_penalized = scores_hvg.copy()
    for j in range(len(scores_hvg_penalized)):
        if eval_mask_hvg[j]:
            if not has_parent_hvg[j]:  # Gene has no connections in GRN
                scores_hvg_penalized[j] = 0.0  # Penalize by setting score to 0
            elif np.isnan(scores_hvg_penalized[j]):
                scores_hvg_penalized[j] = 0.0  # Also set NaN to 0
    
    valid_scores_hvg = scores_hvg_penalized[eval_mask_hvg]
    
    if len(valid_scores_hvg) == 0:
        print("WARNING: No valid HVG genes to evaluate!")
        sem_hvg_score = 0.0
    else:
        sem_hvg_score = float(np.mean(valid_scores_hvg))
        n_missing = np.sum(~has_parent_hvg[eval_mask_hvg])
        print(f"SEM HVG score (mean R²): {sem_hvg_score:.4f}")
        print(f"HVG genes evaluated: {len(valid_scores_hvg)}")
        print(f"HVG genes missing in GRN (penalized with 0): {n_missing}")
        print(f"SEM HVG score (min): {np.min(valid_scores_hvg):.4f}")
        print(f"SEM HVG score (max): {np.max(valid_scores_hvg):.4f}")
    
    results = {
        'sem_grn': [float(sem_grn_score)],
        'sem_hvg': [float(sem_hvg_score)],
        'sem': [float((sem_grn_score + sem_hvg_score) / 2)]
    }

    df_results = pd.DataFrame(results)
    return df_results