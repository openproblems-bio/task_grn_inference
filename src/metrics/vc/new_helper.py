import os
import traceback
from typing import Tuple, Dict
import sys
import tqdm
import random
import h5py
import numpy as np
import pandas as pd
import torch
import xgboost
from scipy.sparse import csr_matrix
from scipy.stats import wilcoxon
from sklearn.decomposition import PCA
from sklearn.linear_model import Ridge
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import StratifiedGroupKFold, StratifiedKFold
from sklearn.preprocessing import StandardScaler, LabelEncoder, TargetEncoder
from torch.utils.data import Dataset
import anndata as ad
import warnings
warnings.filterwarnings("ignore", category=UserWarning)

os.environ['CUBLAS_WORKSPACE_CONFIG'] = ':4096:8' # For reproducibility purposes (on GPU)

# Hyper-parameters
DEVICE = 'cuda' if torch.cuda.is_available() else 'cpu'
NUMPY_DTYPE = np.float32

# For reproducibility purposes
def set_seed():
    seed = 0xCAFE
    os.environ['PYTHONHASHSEED'] = str(seed)
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    torch.use_deterministic_algorithms(True)

from util import read_prediction, manage_layer
from dataset_config import DATASET_GROUPS
from scipy.spatial.distance import cityblock


def combine_multi_index(*arrays) -> np.array:
    """Combine parallel label arrays into a single integer label per position."""
    A = np.stack(arrays)
    n_classes = tuple(int(A[i].max()) + 1 for i in range(A.shape[0]))
    return np.ravel_multi_index(A, dims=n_classes, order='C')


def encode_obs_cols(adata, cols):
    encoded = []
    for col in cols:
        if col in adata.obs:
            codes = LabelEncoder().fit_transform(adata.obs[col].values)
            encoded.append(codes)
    return encoded


def create_control_matching(are_controls: np.ndarray, match_groups: np.ndarray) -> Tuple[Dict[int, int], np.ndarray]:
    control_indices = np.where(are_controls)[0]
    
    if len(control_indices) == 0:
        raise ValueError("No control samples found in dataset!")
    
    # First, try to create exact matching (original approach)
    control_map = {}
    for i, (is_control, group_id) in enumerate(zip(are_controls, match_groups)):
        if is_control:
            control_map[int(group_id)] = i

    return control_map, match_groups


class Model(object):

    def __init__(self, A: np.ndarray):
        self.n_genes: int = A.shape[1]
        self.A: np.ndarray = A
        self._te = [TargetEncoder() for _ in range(self.n_genes)]
        self._models = [xgboost.XGBRegressor(n_estimators=100) for _ in range(self.n_genes)]
        #self._models = [RandomForestRegressor() for _ in range(self.n_genes)]

    def fit(self, X_train: np.ndarray, Y_train: np.ndarray, p_train: np.ndarray) -> None:
        Z_train = []
        for j in range(self.n_genes):
            self._te[j].fit(p_train[:, np.newaxis], Y_train[:, j])
            Z_train.append(np.squeeze(self._te[j].transform(p_train[:, np.newaxis])))
        Z_train = np.asarray(Z_train).T
        for j in tqdm.tqdm(range(self.n_genes), desc="Fit"):
            mask = (self.A[:, j] != 0)
            mask[j] = True
            if np.any(mask):
                self._models[j].fit(np.concatenate((X_train[:, mask], Z_train[:, mask]), axis=1), Y_train[:, j])

    def predict(self, X_test: np.ndarray, p_test: np.ndarray) -> np.ndarray:
        Z_test = []
        for j in range(self.n_genes):
            Z_test.append(np.squeeze(self._te[j].transform(p_test[:, np.newaxis])))
        Z_test = np.asarray(Z_test).T
        Y_hat = []
        for j in tqdm.tqdm(range(self.n_genes), desc="Predict"):
            mask = (self.A[:, j] != 0)
            mask[j] = True
            if np.any(mask):
                Y_hat.append(self._models[j].predict(np.concatenate((X_test[:, mask], Z_test[:, mask]), axis=1)))
            else:
                Y_hat.append(np.zeros(len(X_test)))
        return np.asarray(Y_hat).T


def main(par):

    # Load evaluation data
    adata = ad.read_h5ad(par['evaluation_data'])
    dataset_id = adata.uns['dataset_id']
    method_id = ad.read_h5ad(par['prediction'], backed='r').uns['method_id']
    
    # Configure dataset-specific groups
    par['cv_groups'] = DATASET_GROUPS[dataset_id]['cv']
    par['match'] = DATASET_GROUPS[dataset_id]['match'] 
    par['loose_match'] = DATASET_GROUPS[dataset_id]['loose_match']

    layer = manage_layer(adata, par)
    X = adata.layers[layer]
    if isinstance(X, csr_matrix):
        X = X.toarray()
    X = X.astype(np.float32)
    
    gene_names = adata.var_names
    gene_dict = {gene_name: i for i, gene_name in enumerate(gene_names)}

    # Get sample info
    perturbations = None
    cv_groups, match_groups, loose_match_groups = [], [], []
    
    # Encode observation columns
    cv_groups = encode_obs_cols(adata, par['cv_groups'])
    match_groups = encode_obs_cols(adata, par['match'])
    loose_match_groups = encode_obs_cols(adata, par['loose_match'])

    # Get cell types
    N_FOLDS = 5
    try:
        cell_types = np.squeeze(encode_obs_cols(adata, ["cell_type"]))
    except Exception:
        print(traceback.format_exc())
        cell_types = np.random.randint(0, 5, size=len(X))

    # Set perturbations to first column (perturbation)
    perturbations = cv_groups[0]  # perturbation codes
    
    # Groups used for cross-validation
    cv_groups = combine_multi_index(*cv_groups)

    # Groups used for matching with negative controls
    match_groups = combine_multi_index(*match_groups)

    # Groups used for loose matching with negative controls
    loose_match_groups = combine_multi_index(*loose_match_groups)

    # Check for control samples - this is required for perturbation evaluation
    if "is_control" not in adata.obs.columns:
        raise ValueError("Dataset must contain 'is_control' column for perturbation evaluation")
    
    are_controls = adata.obs["is_control"].values.astype(bool)
    
    # Check if we have any control samples
    n_controls = np.sum(are_controls)
    n_perturbed = np.sum(~are_controls)
    
    if n_controls == 0:
        raise ValueError(f"No control samples found in dataset! Found {n_perturbed} perturbed samples but 0 controls. "
                        "Perturbation evaluation requires control samples for comparison.")
    
    print(f"Found {n_controls} control samples and {n_perturbed} perturbed samples")

    # Load inferred GRN
    net = read_prediction(par)
    sources = net["source"].to_numpy()
    targets = net["target"].to_numpy()
    weights = net["weight"].to_numpy()

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

    # Filter genes based on GRN instead of HVGs
    # Keep all genes that are present in the GRN (already filtered above)
    print(f"Using {len(gene_names)} genes present in the GRN")
    
    # Additional memory-aware gene filtering for very large GRNs
    MAX_GENES_FOR_MEMORY = 3000  # Reduced further to avoid memory issues
    if len(gene_names) > MAX_GENES_FOR_MEMORY:
        print(f"Too many genes ({len(gene_names)}) for memory. Selecting top {MAX_GENES_FOR_MEMORY} by GRN connectivity.")
        
        # Select genes with highest connectivity in the GRN
        gene_connectivity = np.sum(np.abs(A), axis=0) + np.sum(np.abs(A), axis=1)
        top_gene_indices = np.argsort(gene_connectivity)[-MAX_GENES_FOR_MEMORY:]
        
        X = X[:, top_gene_indices]
        A = A[top_gene_indices, :][:, top_gene_indices]
        gene_names = gene_names[top_gene_indices]
        
        print(f"Final: Using {len(gene_names)} most connected genes for evaluation")

    # Remove self-regulations
    np.fill_diagonal(A, 0)

    # Mapping between gene expression profiles and their matched negative controls
    control_map, _ = create_control_matching(are_controls, match_groups)
    loose_control_map, _ = create_control_matching(are_controls, loose_match_groups)

    # Build dataset
    X_, Y_ = [], []
    for i, (match_group, loose_match_group) in enumerate(zip(match_groups, loose_match_groups)):
        if match_group in control_map:
            i_prime = control_map[match_group]
        else:
            i_prime = loose_control_map[loose_match_group]
        X_.append(X[i_prime, :])
        Y_.append(X[i, :])
    X, Y = np.asarray(X_), np.asarray(Y_)
    Y = Y - X
    del X_, Y_

    eps = 1e-20
    r2, r2_baseline = [], []
    cv = StratifiedGroupKFold(n_splits=N_FOLDS, shuffle=True, random_state=0xCAFE)
    for i, (train_index, test_index) in enumerate(cv.split(X, perturbations, cell_types)):

        if (len(train_index) == 0) or (len(test_index) == 0):
            continue

        # Create baseline model
        A_baseline = np.copy(A)
        for j in range(A_baseline.shape[1]):
            np.random.shuffle(A_baseline[:j, j])
            np.random.shuffle(A_baseline[j+1:, j])
        assert np.any(A_baseline != A)

        # Center and scale dataset
        scaler = StandardScaler()
        X_train = scaler.fit_transform(X[train_index, :])
        X_test = scaler.transform(X[test_index, :])
        scaler = StandardScaler()
        Y_train = scaler.fit_transform(Y[train_index, :])
        Y_test = scaler.transform(Y[test_index, :])
        p_train = perturbations[train_index]
        p_test = perturbations[test_index]

        model = Model(A)
        model.fit(X_train, Y_train, p_train)
        Y_hat = model.predict(X_test, p_test)
        ss_res = np.sum(np.square(Y_test - Y_hat)) + eps
        ss_tot = np.sum(np.square(Y_test - np.mean(Y_test, axis=0)[np.newaxis, :])) + eps
        r2.append(1 - ss_res / ss_tot)

        model = Model(A_baseline)
        model.fit(X_train, Y_train, p_train)
        Y_hat = model.predict(X_test, p_test)
        ss_res = np.sum(np.square(Y_test - Y_hat)) + eps
        ss_tot = np.sum(np.square(Y_test - np.mean(Y_test, axis=0)[np.newaxis, :])) + eps
        r2_baseline.append(1 - ss_res / ss_tot)

        print("r2", np.mean(r2), np.mean(r2_baseline))

    r2 = np.asarray(r2).flatten()
    r2_baseline = np.asarray(r2_baseline).flatten()

    print("R2", np.mean(r2), np.mean(r2_baseline))

    if np.all(r2 == r2_baseline):
        final_score = 0
    else:
        print(np.mean(r2), np.mean(r2_baseline))
        p_value = wilcoxon(r2, r2_baseline, alternative="greater").pvalue
        final_score = -np.log10(p_value)

    print(f"Method: {method_id}")
    print(f"R2: {final_score}")

    results = {
        'vc': [float(final_score)]
    }

    df_results = pd.DataFrame(results)
    return df_results
