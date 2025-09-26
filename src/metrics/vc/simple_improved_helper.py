"""
Simplified improved VC metric that follows the original architecture but with stability fixes.
This version is more memory-efficient and suitable for large gene networks.
"""

import torch
import torch.nn as nn
import numpy as np
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.model_selection import GroupKFold
from torch.utils.data import Dataset, DataLoader
import logging
from typing import Dict, Tuple
import tqdm
import anndata as ad
import pandas as pd
from scipy.sparse import csr_matrix

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Constants
NUMPY_DTYPE = np.float32
DEVICE = torch.device('cpu')  # Force CPU to avoid memory issues

# Utility functions
def combine_multi_index(*arrays):
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
    control_map = {}
    for i, (is_control, group_id) in enumerate(zip(are_controls, match_groups)):
        if is_control:
            control_map[int(group_id)] = i
    for i in range(len(are_controls)):
        assert match_groups[i] in control_map
    return control_map, match_groups

class ImprovedGRNLayer(torch.nn.Module):
    """Improved GRN layer with better numerical stability."""
    
    def __init__(self, A_weights: torch.nn.Parameter, A_signs: torch.Tensor, 
                 signed: bool = True, inverse: bool = True, alpha: float = 0.1):
        super().__init__()
        self.n_genes: int = A_weights.size(1)
        self.A_weights: torch.nn.Parameter = A_weights
        self.register_buffer('A_signs', A_signs.to(A_weights.device))
        self.register_buffer('A_mask', (A_signs > 0).to(self.A_weights.dtype).to(A_weights.device))
        self.register_buffer('I', torch.eye(self.n_genes, dtype=A_weights.dtype, device=A_weights.device))
        self.signed: bool = signed
        self.inverse: bool = inverse
        self.alpha: float = alpha  # Reduced from 1.0 to 0.1
        self.ridge_reg: float = 1e-6  # Ridge regularization

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        # Scale weights more conservatively
        A = self.A_weights * self.A_mask
        if self.signed:
            A = A * self.A_signs
        
        # Apply tanh to limit weight range
        A = torch.tanh(A) * 0.5  # Limit to [-0.5, 0.5]
        
        if self.inverse:
            # Compute (I - alpha*A^T) with ridge regularization
            M = self.I - self.alpha * A.T + self.ridge_reg * self.I
            try:
                # Use solve instead of inverse for better numerical stability
                y = torch.linalg.solve(M, x.unsqueeze(-1)).squeeze(-1)
                
                # Check for numerical issues
                if torch.isnan(y).any() or torch.isinf(y).any():
                    logger.warning("Numerical instability detected in GRN layer, using fallback")
                    return x  # Fallback to input
                    
                return y
            except torch.linalg.LinAlgError:
                logger.warning("Matrix solve failed in GRN layer, using fallback")
                return x
        else:
            # Forward direction
            return x + self.alpha * torch.matmul(x, A.T)

class ImprovedModel(torch.nn.Module):
    """Improved model with the original architecture but better stability."""
    
    def __init__(self, A: np.ndarray, n_perturbations: int, n_hidden: int = 64, signed: bool = True):
        super().__init__()
        self.n_genes: int = A.shape[0]
        self.n_perturbations: int = n_perturbations
        self.n_hidden: int = min(n_hidden, self.n_genes // 2)  # Adaptive hidden size
        
        # Perturbation embedding
        self.perturbation_embedding = torch.nn.Embedding(n_perturbations, self.n_hidden)

        # Initialize A with better scaling
        A_signs = torch.from_numpy(np.sign(A).astype(NUMPY_DTYPE))
        A_weights = np.copy(A).astype(NUMPY_DTYPE)
        # More conservative scaling
        A_weights /= (2.0 * float(np.sqrt(self.n_genes)) * (float(np.std(A_weights)) + 1e-8))
        A_weights = torch.nn.Parameter(torch.from_numpy(A_weights))
        A_signs = A_signs.to(A_weights.device)

        # Encoder with improved GRN layer
        self.encoder = torch.nn.Sequential(
            ImprovedGRNLayer(A_weights, A_signs, inverse=False, signed=signed, alpha=0.1),
            torch.nn.PReLU(1),
            torch.nn.Linear(self.n_genes, self.n_hidden),
            torch.nn.Dropout(0.2),  # Add dropout
            torch.nn.PReLU(1),
            torch.nn.Linear(self.n_hidden, self.n_hidden),
            torch.nn.Dropout(0.2),
            torch.nn.PReLU(1),
        )
        
        # Decoder with improved GRN layer
        self.decoder = torch.nn.Sequential(
            torch.nn.Linear(self.n_hidden, self.n_hidden),
            torch.nn.Dropout(0.2),
            torch.nn.PReLU(1),
            torch.nn.Linear(self.n_hidden, self.n_genes),
            torch.nn.PReLU(1),
            ImprovedGRNLayer(A_weights, A_signs, inverse=True, signed=signed, alpha=0.1)
        )
        
        # Initialize weights
        self._initialize_weights()

    def _initialize_weights(self):
        """Initialize weights with smaller values."""
        for module in self.modules():
            if isinstance(module, torch.nn.Linear):
                torch.nn.init.xavier_uniform_(module.weight, gain=0.1)
                if module.bias is not None:
                    torch.nn.init.zeros_(module.bias)
            elif isinstance(module, torch.nn.Embedding):
                torch.nn.init.normal_(module.weight, mean=0.0, std=0.1)

    def forward(self, x: torch.Tensor, pert: torch.LongTensor) -> torch.Tensor:
        y = self.encoder(x)
        z = self.perturbation_embedding(pert)
        y = y + z
        y = self.decoder(y)
        return x + y

class PerturbationDataset(Dataset):
    """Compatible dataset class matching the original."""
    
    def __init__(self, X: np.ndarray, idx: np.ndarray, match_groups: np.ndarray,
                 control_map: Dict[int, int], perturbations: np.ndarray):
        self.X: np.ndarray = X.astype(NUMPY_DTYPE)
        self.idx: np.ndarray = idx
        self.match_groups: np.ndarray = match_groups
        self.control_map: Dict[int, int] = control_map
        self.perturbations: np.ndarray = perturbations.astype(int)

    def __len__(self) -> int:
        return len(self.idx)

    def __getitem__(self, i: int) -> Tuple[torch.Tensor, int, torch.Tensor]:
        i = self.idx[i]
        y = torch.from_numpy(self.X[i, :])
        group = int(self.match_groups[i])
        j = int(self.control_map[group])
        x = torch.from_numpy(self.X[j, :])
        p = int(self.perturbations[i])
        return x, p, y

def improved_evaluate(A, train_data_loader, test_data_loader, n_perturbations: int) -> float:
    """Improved evaluation with better training and early stopping."""
    # Training
    signed = np.any(A < 0)
    model = ImprovedModel(A, n_perturbations, n_hidden=64, signed=signed)
    model = model.to(DEVICE)
    
    # Improved optimizer settings
    optimizer = torch.optim.Adam(model.parameters(), lr=0.001, weight_decay=1e-5)
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
        optimizer, mode='min', patience=5, min_lr=1e-5, cooldown=3, factor=0.8
    )
    
    # Training loop with early stopping
    pbar = tqdm.tqdm(range(200))  # Reduced epochs
    best_val_loss = float('inf')
    patience_counter = 0
    patience = 10
    
    model.train()
    for epoch in pbar:
        total_loss = 0
        for x, pert, y in train_data_loader:
            x, pert, y = x.to(DEVICE), pert.to(DEVICE), y.to(DEVICE)
            optimizer.zero_grad()
            y_hat = model(x, pert)
            loss = torch.mean(torch.square(y - y_hat))
            loss.backward()
            
            # Gradient clipping
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
            
            optimizer.step()
            total_loss += loss.item() * len(x)
        
        pbar.set_description(f"Loss: {total_loss:.2f}")
        
        # Validation
        model.eval()
        val_loss = 0
        with torch.no_grad():
            for x, pert, y in test_data_loader:
                x, pert, y = x.to(DEVICE), pert.to(DEVICE), y.to(DEVICE)
                y_hat = model(x, pert)
                loss = torch.mean(torch.square(y - y_hat))
                val_loss += loss.item() * len(x)
        
        scheduler.step(val_loss)
        
        # Early stopping
        if val_loss < best_val_loss:
            best_val_loss = val_loss
            patience_counter = 0
        else:
            patience_counter += 1
            
        if patience_counter >= patience:
            logger.info(f"Early stopping at epoch {epoch}")
            break
            
        model.train()

    logger.info(f"Training completed. Best validation loss: {best_val_loss}")
    return best_val_loss

def create_improved_baseline_grn(A_true, method='weight_shuffled'):
    """Create improved baseline GRN."""
    if method == 'weight_shuffled':
        A_baseline = A_true.copy()
        non_zero_mask = A_baseline != 0
        non_zero_values = A_baseline[non_zero_mask]
        np.random.shuffle(non_zero_values)
        A_baseline[non_zero_mask] = non_zero_values
        return A_baseline
    else:  # random (original method)
        A_baseline = A_true.T.copy()
        np.random.shuffle(A_baseline)
        return A_baseline.T

def main(par):
    """Main function compatible with original interface."""
    import sys
    import os
    
    # Import necessary modules
    sys.path.append(os.path.join(os.path.dirname(__file__), '../../utils'))
    from dataset_config import DATASET_GROUPS
    from util import manage_layer, read_prediction
    
    # Load evaluation data
    adata = ad.read_h5ad(par['evaluation_data'])
    dataset_id = adata.uns['dataset_id']
    method_id = ad.read_h5ad(par['prediction'], backed='r').uns['method_id']
    
    logger.info(f"Processing dataset: {dataset_id}, method: {method_id}")
    
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
    cv_groups = encode_obs_cols(adata, par['cv_groups'])
    match_groups = encode_obs_cols(adata, par['match'])
    loose_match_groups = encode_obs_cols(adata, par['loose_match'])
    
    perturbations = cv_groups[0]  # perturbation codes
    cv_groups = combine_multi_index(*cv_groups)
    match_groups = combine_multi_index(*match_groups)
    loose_match_groups = combine_multi_index(*loose_match_groups)

    are_controls = adata.obs["is_control"].values.astype(bool)

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

    # Only consider genes present in the GRN
    gene_mask = np.logical_or(np.any(A, axis=1), np.any(A, axis=0))
    X = X[:, gene_mask]
    A = A[gene_mask, :][:, gene_mask]
    gene_names = gene_names[gene_mask]

    # Restrict to most variable genes
    std = np.std(X, axis=0)
    idx = np.argsort(std)[-min(2000, len(std)):]  # Reduced to 2000 for memory
    X = X[:, idx]
    A = A[idx, :][:, idx]
    gene_names = gene_names[idx]
    
    logger.info(f"Using {len(gene_names)} genes for evaluation")

    # Control matching
    try:
        control_map, match_groups = create_control_matching(are_controls, match_groups)
    except AssertionError:
        logger.info("Failed to match with controls exactly. Using less stringent matching.")
        control_map, match_groups = create_control_matching(are_controls, loose_match_groups)

    ss_res = 0
    ss_tot = 0
    cv = GroupKFold(n_splits=5)
    
    for i, (train_index, test_index) in enumerate(cv.split(X, X, cv_groups)):
        logger.info(f"Processing fold {i+1}/5")
        
        # Standardize data
        scaler = StandardScaler()
        scaler.fit(X[train_index, :])
        X_standardized = scaler.transform(X)

        # Create datasets
        n_perturbations = int(np.max(perturbations) + 1)
        train_dataset = PerturbationDataset(X_standardized, train_index, match_groups, control_map, perturbations)
        train_data_loader = torch.utils.data.DataLoader(train_dataset, batch_size=512)
        test_dataset = PerturbationDataset(X_standardized, test_index, match_groups, control_map, perturbations)
        test_data_loader = torch.utils.data.DataLoader(test_dataset, batch_size=512)

        # Evaluate inferred GRN
        ss_res += improved_evaluate(A, train_data_loader, test_data_loader, n_perturbations)

        # Evaluate baseline GRN
        A_baseline = create_improved_baseline_grn(A, method='weight_shuffled')
        ss_tot += improved_evaluate(A_baseline, train_data_loader, test_data_loader, n_perturbations)

    # Calculate R²
    r2 = 1 - ss_res / ss_tot if ss_tot > 1e-10 else 0.0
    
    logger.info(f"Method: {method_id}")
    logger.info(f"Improved R²: {r2}")
    logger.info(f"SS_res: {ss_res}, SS_tot: {ss_tot}")

    results = {'vc': [float(r2)]}
    return pd.DataFrame(results)