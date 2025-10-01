import os
from typing import Tuple, Dict
import sys
import tqdm
import random
import h5py
import numpy as np
import pandas as pd
import torch
from scipy.sparse import csr_matrix
from sklearn.model_selection import GroupKFold
from sklearn.preprocessing import StandardScaler, LabelEncoder
from torch.utils.data import Dataset
import anndata as ad
import warnings
warnings.filterwarnings("ignore", category=UserWarning)

os.environ['CUBLAS_WORKSPACE_CONFIG'] = ':4096:8' # For reproducibility purposes (on GPU)

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


def compute_pds_cell_level(X_true, X_pred, perturbations, gene_names, max_cells_per_pert=10):
    """
    Compute PDS at individual cell level (more challenging than mean-based PDS).
    
    For each individual predicted cell, find its rank when compared to 
    true mean profiles of all perturbations.
    """
    unique_perts = np.unique(perturbations)
    n_perts = len(unique_perts)
    
    print(f"Computing cell-level PDS for {n_perts} unique perturbations")
    
    if n_perts < 2:
        return 0.0
    
    # Compute mean true expression profiles per perturbation
    true_means = {}
    for pert in unique_perts:
        mask = (perturbations == pert)
        if np.sum(mask) == 0:
            continue
        true_means[pert] = np.mean(X_true[mask, :], axis=0)
    
    all_scores = []
    
    # For each perturbation, sample some cells and compute their PDS
    for pert in unique_perts:
        if pert not in true_means:
            continue
            
        # Get cells for this perturbation
        mask = (perturbations == pert)
        pert_indices = np.where(mask)[0]
        
        if len(pert_indices) == 0:
            continue
            
        # Sample max_cells_per_pert cells to avoid bias from perturbations with many cells
        # Use seeded random generator for reproducibility
        cell_rng = np.random.RandomState(seed + int(pert))  # Different seed per perturbation
        sampled_indices = cell_rng.choice(
            pert_indices, 
            size=min(max_cells_per_pert, len(pert_indices)), 
            replace=False
        )
        
        for cell_idx in sampled_indices:
            # Get predicted profile for this individual cell
            pred_vec = X_pred[cell_idx, :].copy()
            
            # Calculate distances to all true mean profiles
            dists = []
            for t in unique_perts:
                if t not in true_means:
                    continue
                true_vec = true_means[t].copy()
                
                # Remove target gene if it exists
                if str(pert) in gene_names:
                    gene_idx = np.where(gene_names == str(pert))[0]
                    if len(gene_idx) > 0:
                        true_vec = np.delete(true_vec, gene_idx)
                        pred_vec_temp = np.delete(pred_vec, gene_idx)
                    else:
                        pred_vec_temp = pred_vec
                else:
                    pred_vec_temp = pred_vec
                
                dist = cityblock(pred_vec_temp, true_vec)
                dists.append((t, dist))
            
            # Sort by distance and find rank of correct perturbation
            dists_sorted = sorted(dists, key=lambda x: x[1])
            true_rank = next((i for i, (t, _) in enumerate(dists_sorted) if t == pert), n_perts-1)
            
            # Cell-level PDS
            pds = 1 - (true_rank / (n_perts - 1)) if n_perts > 1 else 1.0
            all_scores.append(pds)
    
    mean_pds = np.mean(all_scores) if all_scores else 0.0
    print(f"Cell-level PDS scores: min={min(all_scores):.3f}, max={max(all_scores):.3f}, mean={mean_pds:.3f} (n_cells={len(all_scores)})")
    return mean_pds


def compute_pds(X_true, X_pred, perturbations, gene_names):
    """
    Compute both mean-level and cell-level PDS for comparison.
    """
    # Mean-level PDS (original approach)
    unique_perts = np.unique(perturbations)
    n_perts = len(unique_perts)
    
    print(f"Computing mean-level PDS for {n_perts} unique perturbations")
    
    if n_perts < 2:
        return 0.0, 0.0
    
    # Compute mean expression profiles per perturbation
    true_means = {}
    pred_means = {}
    
    for pert in unique_perts:
        mask = (perturbations == pert)
        if np.sum(mask) == 0:
            continue
        true_means[pert] = np.mean(X_true[mask, :], axis=0)
        pred_means[pert] = np.mean(X_pred[mask, :], axis=0)
    
    scores = {}
    for pert in unique_perts:
        if pert not in pred_means or pert not in true_means:
            continue
            
        pred_vec = pred_means[pert].copy()
        dists = []
        
        for t in unique_perts:
            if t not in true_means:
                continue
            true_vec = true_means[t].copy()
            
            if str(pert) in gene_names:
                gene_idx = np.where(gene_names == str(pert))[0]
                if len(gene_idx) > 0:
                    true_vec = np.delete(true_vec, gene_idx)
                    pred_vec_temp = np.delete(pred_vec, gene_idx)
                else:
                    pred_vec_temp = pred_vec
            else:
                pred_vec_temp = pred_vec
            
            dist = cityblock(pred_vec_temp, true_vec)
            dists.append((t, dist))
        
        dists_sorted = sorted(dists, key=lambda x: x[1])
        true_rank = next((i for i, (t, _) in enumerate(dists_sorted) if t == pert), n_perts-1)
        pds = 1 - (true_rank / (n_perts - 1)) if n_perts > 1 else 1.0
        scores[pert] = pds
    
    mean_level_pds = np.mean(list(scores.values())) if scores else 0.0
    print(f"Mean-level PDS: min={min(scores.values()):.3f}, max={max(scores.values()):.3f}, mean={mean_level_pds:.3f}")
    
    # Cell-level PDS (more challenging)
    cell_level_pds = compute_pds_cell_level(X_true, X_pred, perturbations, gene_names)
    
    return mean_level_pds, cell_level_pds


class GRNLayer(torch.nn.Module):

    def __init__(
            self,
            A_weights: torch.nn.Parameter,
            A_signs: torch.Tensor,
            signed: bool = True,
            inverse: bool = True,
            alpha: float = 1.0
    ):
        torch.nn.Module.__init__(self)
        self.n_genes: int = A_weights.size(1)
        self.A_weights: torch.nn.Parameter = A_weights
        self.register_buffer('A_signs', A_signs.to(A_weights.device))
        self.register_buffer('A_mask', (A_signs > 0).to(self.A_weights.dtype).to(A_weights.device))
        self.register_buffer('I', torch.eye(self.n_genes, dtype=A_weights.dtype, device=A_weights.device))
        self.signed: bool = signed
        self.inverse: bool = inverse
        self.alpha: float = alpha

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        if self.signed:
            A = torch.abs(self.A_weights) * self.A_signs
        else:
            A = self.A_weights * self.A_mask
        
        if self.inverse:
            # For inverse transformation, use iterative solve to avoid memory issues
            # Solve (I - alpha * A.t()) * y = x for y
            ia = self.I - self.alpha * A.t()
            
            # Add small regularization to diagonal for numerical stability
            ia = ia + 1e-6 * self.I
            
            # Use solve instead of inversion to save memory
            try:
                # Solve ia * y.t() = x.t() for y.t(), then transpose
                result = torch.linalg.solve(ia, x.t()).t()
                return result
            except torch.linalg.LinAlgError:
                # Fallback: simple linear transformation without inversion
                print("Warning: Matrix solve failed, using simplified GRN transformation")
                return torch.mm(x, A)
        else:
            # Forward transformation: apply GRN directly
            return torch.mm(x, A.t())


class Model(torch.nn.Module):
    def __init__(self, A: np.ndarray, n_perturbations: int, n_hidden: int = 64, signed: bool = True):
        torch.nn.Module.__init__(self)
        self.n_genes: int = A.shape[1]
        self.n_perturbations: int = n_perturbations
        self.n_hidden: int = n_hidden
        self.perturbation_embedding = torch.nn.Embedding(n_perturbations, n_hidden)

        A_signs = torch.from_numpy(np.sign(A).astype(NUMPY_DTYPE))
        A_weights = np.copy(A).astype(NUMPY_DTYPE)
        A_weights /= (np.sqrt(self.n_genes) * float(np.std(A_weights)))
        A_weights = torch.nn.Parameter(torch.from_numpy(A_weights))
        # Ensure A_signs is on the same device as A_weights
        A_signs = A_signs.to(A_weights.device)

        # First layer: GRN-informed transformation of control expression
        self.grn_input_layer = GRNLayer(A_weights, A_signs, inverse=False, signed=signed, alpha=0.1)
        
        # Middle layers: perturbation processing
        self.encoder = torch.nn.Sequential(
            torch.nn.PReLU(1),
            torch.nn.Linear(self.n_genes, self.n_hidden),
            torch.nn.PReLU(1),
            torch.nn.Linear(self.n_hidden, self.n_hidden),
            torch.nn.PReLU(1),
        )
        
        self.decoder = torch.nn.Sequential(
            torch.nn.Linear(self.n_hidden, self.n_hidden),
            torch.nn.PReLU(1),
            torch.nn.Linear(self.n_hidden, self.n_genes),
            torch.nn.PReLU(1),
        )
        
        # Last layer: GRN-informed transformation to final expression
        self.grn_output_layer = GRNLayer(A_weights, A_signs, inverse=True, signed=signed, alpha=0.1)

    def forward(self, x: torch.Tensor, pert: torch.LongTensor) -> torch.Tensor:
        # Apply GRN transformation to control expression
        x = self.grn_input_layer(x)
        y = self.encoder(x)
        z = self.perturbation_embedding(pert)
        y = y + z
        y = self.decoder(y)
        y = self.grn_output_layer(y)
        return y


class PerturbationDataset(Dataset):

    def __init__(
            self,
            X: np.ndarray,
            idx: np.ndarray,
            match_groups: np.ndarray,
            control_map: Dict[int, int],
            loose_match_groups: np.ndarray,
            loose_control_map: Dict[int, int],
            perturbations: np.ndarray
    ):
        self.X: np.ndarray = X.astype(NUMPY_DTYPE)
        self.idx: np.ndarray = idx
        self.match_groups: np.ndarray = match_groups
        self.control_map: Dict[int, int] = control_map
        self.loose_match_groups: np.ndarray = loose_match_groups
        self.loose_control_map: Dict[int, int] = loose_control_map
        self.perturbations: np.ndarray = perturbations.astype(int)

    def __len__(self) -> int:
        return len(self.idx)

    def __getitem__(self, i: int) -> Tuple[torch.Tensor, int, torch.Tensor]:
        i = self.idx[i]
        y = torch.from_numpy(self.X[i, :])
        group = int(self.match_groups[i])

        if group in self.control_map:
            j = int(self.control_map[group])
        else:
            group = int(self.loose_match_groups[i])
            j = int(self.loose_control_map[group])

        x = torch.from_numpy(self.X[j, :])
        p = int(self.perturbations[i])
        return x, p, y



def evaluate(A, train_data_loader, test_data_loader, n_perturbations: int) -> Tuple[float, float]:
    # Training

    signed = np.any(A < 0)
    model = Model(A, n_perturbations, n_hidden=16, signed=signed)
    model = model.to(DEVICE)
    optimizer = torch.optim.Adam(model.parameters(), lr=0.001, weight_decay=1e-5)
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
        optimizer, mode='min', patience=5,
        min_lr=1e-5, cooldown=3, factor=0.8
    )
    pbar = tqdm.tqdm(range(100))  # Reduced epochs for faster testing
    best_val_loss = float('inf')
    best_ss_res = None
    best_epoch = 0
    model.train()
    for epoch in pbar:
        total_loss = 0
        for x, pert, y in train_data_loader:
            x, pert, y = x.to(DEVICE), pert.to(DEVICE), y.to(DEVICE)
            optimizer.zero_grad()
            # Model now predicts full perturbed expression directly
            y_hat = model(x, pert)
            loss = torch.mean(torch.square(y - y_hat))
            loss.backward()
            optimizer.step()
            total_loss += loss.item() * len(x)
        pbar.set_description(str(total_loss))
        scheduler.step(total_loss)

        model.eval()
        ss_res = 0
        with torch.no_grad():
            for x, pert, y in test_data_loader:
                x, pert, y = x.to(DEVICE), pert.to(DEVICE), y.to(DEVICE)
                # Model predicts full perturbed expression
                y_hat = model(x, pert)
                residuals = torch.square(y - y_hat).cpu().data.numpy()
                ss_res += np.sum(residuals, axis=0)
        if np.sum(ss_res) < best_val_loss:
            best_val_loss = np.sum(ss_res)
            best_epoch = epoch
            best_ss_res = ss_res
        model.train()
    ss_res = best_ss_res

    # Final evaluation with PDS
    model.eval()

    ss_tot = 0

    with torch.no_grad():
        for x, pert, y in test_data_loader:
            x, pert, y = x.to(DEVICE), pert.to(DEVICE), y.to(DEVICE)
            y_hat = model(x, pert)

            residuals = torch.square(y - torch.mean(y, dim=0).unsqueeze(0)).cpu().data.numpy()
            ss_tot += np.sum(residuals, axis=0)

    print(f"Best epoch: {best_epoch} ({best_val_loss})")
    return best_ss_res, ss_tot



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

    # Create baseline model
    A_baseline = np.copy(A)
    for j in range(A.shape[1]):
        np.random.shuffle(A[:j, j])
        np.random.shuffle(A[j+1:, j])
    assert np.any(A_baseline != A)

    # Mapping between gene expression profiles and their matched negative controls

    control_map, _ = create_control_matching(are_controls, match_groups)
    loose_control_map, _ = create_control_matching(are_controls, loose_match_groups)

    ss_res = 0
    ss_tot = 0
    cv = GroupKFold(n_splits=5)
    
    results = []
    
    for i, (train_index, test_index) in enumerate(cv.split(X, X, cv_groups)):
        # Center and scale dataset
        scaler = StandardScaler()
        scaler.fit(X[train_index, :])
        X_standardized = scaler.transform(X)

        # Create data loaders
        n_perturbations = int(np.max(perturbations) + 1)
        train_dataset = PerturbationDataset(
            X_standardized,
            train_index,
            match_groups,
            control_map,
            loose_match_groups,
            loose_control_map,
            perturbations
        )
        train_data_loader = torch.utils.data.DataLoader(train_dataset, batch_size=512)
        test_dataset = PerturbationDataset(
            X_standardized,
            test_index,
            match_groups,
            control_map,
            loose_match_groups,
            loose_control_map,
            perturbations
        )
        test_data_loader = torch.utils.data.DataLoader(test_dataset, batch_size=512)

        # Evaluate inferred GRN
        res = evaluate(A, train_data_loader, test_data_loader, n_perturbations)
        ss_res = ss_res + res[0]
        ss_tot = ss_tot + res[1]

        # Evaluate baseline GRN (shuffled target genes)
        #ss_tot = ss_tot + evaluate(A_baseline, train_data_loader, test_data_loader, n_perturbations)

    r2 = 1 - ss_res / ss_tot

    final_score = np.mean(np.clip(r2, 0, 1))
    print(f"Method: {method_id}")
    print(f"R2: {final_score}")

    results = {
        'vc': [float(final_score)]

    }

    df_results = pd.DataFrame(results)
    return df_results