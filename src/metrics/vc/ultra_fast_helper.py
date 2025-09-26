"""
Ultra-fast VC implementation with smart optimizations:
1. Smart Gene Selection (only genes in GRN)
2. Pre-computed Matrix Operations 
3. Batch Matrix Operations
4. Sparse Matrix Representation
"""

import torch
import torch.nn as nn
import numpy as np
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.model_selection import GroupKFold
from torch.utils.data import Dataset, DataLoader
import logging
from typing import Dict, Tuple, Set
import tqdm
import anndata as ad
import pandas as pd
from scipy.sparse import csr_matrix
import scipy.sparse as sp

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Constants
NUMPY_DTYPE = np.float32
DEVICE = torch.device('cpu')  # Force CPU for now

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

def smart_gene_selection(df_grn, gene_names, max_genes=800):
    """
    Smart gene selection: only include genes actually present in the GRN.
    This dramatically reduces the problem size while keeping all relevant information.
    """
    # Get all genes mentioned in the GRN
    grn_genes = set(df_grn['source'].tolist() + df_grn['target'].tolist())
    
    # Filter to genes present in our data
    available_genes = set(gene_names)
    selected_genes = grn_genes.intersection(available_genes)
    
    logger.info(f"GRN contains {len(grn_genes)} unique genes")
    logger.info(f"Available genes: {len(available_genes)}")
    logger.info(f"Selected genes: {len(selected_genes)}")
    
    # If still too many, prioritize by degree centrality
    if len(selected_genes) > max_genes:
        # Count degree for each gene
        gene_degrees = {}
        for gene in selected_genes:
            in_degree = len(df_grn[df_grn['target'] == gene])
            out_degree = len(df_grn[df_grn['source'] == gene])
            gene_degrees[gene] = in_degree + out_degree
        
        # Select top genes by degree
        sorted_genes = sorted(gene_degrees.items(), key=lambda x: x[1], reverse=True)
        selected_genes = set([gene for gene, _ in sorted_genes[:max_genes]])
        logger.info(f"Reduced to {len(selected_genes)} highest-degree genes")
    
    return selected_genes

def create_sparse_grn_matrix(df_grn, gene_names, selected_genes):
    """
    Create sparse GRN matrix using only selected genes.
    """
    # Create mapping from gene names to indices (only for selected genes)
    selected_gene_list = sorted(selected_genes)
    gene_to_idx = {gene: i for i, gene in enumerate(selected_gene_list)}
    
    # Filter GRN to only include selected genes
    mask = (df_grn['source'].isin(selected_genes)) & (df_grn['target'].isin(selected_genes))
    filtered_grn = df_grn[mask].copy()
    
    if len(filtered_grn) == 0:
        logger.warning("No GRN edges found between selected genes!")
        return np.zeros((len(selected_gene_list), len(selected_gene_list))), selected_gene_list
    
    # Create sparse matrix
    n_genes = len(selected_gene_list)
    rows = [gene_to_idx[source] for source in filtered_grn['source']]
    cols = [gene_to_idx[target] for target in filtered_grn['target']]
    weights = filtered_grn['weight'].astype(NUMPY_DTYPE).values
    
    # Create sparse matrix and convert to dense for PyTorch
    sparse_A = sp.csr_matrix((weights, (rows, cols)), shape=(n_genes, n_genes))
    dense_A = sparse_A.toarray()
    
    logger.info(f"Created {n_genes}x{n_genes} GRN matrix with {len(weights)} edges ({100*len(weights)/(n_genes**2):.2f}% density)")
    
    return dense_A, selected_gene_list

class PrecomputedGRNLayer(nn.Module):
    """
    Pre-computed GRN layer that calculates matrix operations once and reuses them.
    Massive speedup by avoiding repeated matrix inversions.
    """
    
    def __init__(self, A_weights: torch.Tensor, alpha: float = 0.1, ridge_reg: float = 1e-6):
        super().__init__()
        self.n_genes = A_weights.shape[0]
        self.alpha = alpha
        self.ridge_reg = ridge_reg
        
        # Pre-compute the inverse matrix ONCE
        logger.info(f"Pre-computing GRN matrix operations for {self.n_genes} genes...")
        
        # Scale A_weights conservatively
        A_scaled = torch.tanh(A_weights) * 0.3  # Limit to [-0.3, 0.3] for stability
        
        # Compute (I - α*A^T) + ridge regularization
        I = torch.eye(self.n_genes, dtype=A_weights.dtype, device=A_weights.device)
        M = I - self.alpha * A_scaled.T + self.ridge_reg * I
        
        try:
            # Pre-compute the inverse matrix
            inv_M = torch.inverse(M)
            self.register_buffer('inv_matrix', inv_M)
            self.precomputed = True
            logger.info("✓ Matrix inversion successful - using pre-computed operations")
            
            # Check condition number for stability
            cond_num = torch.linalg.cond(M).item()
            if cond_num > 1e10:
                logger.warning(f"High condition number: {cond_num:.2e} - results may be unstable")
                
        except torch.linalg.LinAlgError:
            logger.warning("Matrix inversion failed - using approximation")
            # Fallback to identity + small perturbation
            self.register_buffer('inv_matrix', I + 0.1 * A_scaled.T)
            self.precomputed = False
    
    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """
        Ultra-fast forward pass using pre-computed matrix.
        No matrix inversion per sample - just matrix multiplication!
        """
        if x.dim() == 1:
            x = x.unsqueeze(0)
        
        # Single matrix multiplication instead of expensive inversion
        y = torch.matmul(x, self.inv_matrix.T)
        
        return y.squeeze(0) if y.shape[0] == 1 else y

class PerSamplePredictionModel(nn.Module):
    """
    Model that predicts all genes simultaneously for each sample (per-sample prediction).
    Uses the same GRN dynamics but processes entire gene expression vectors at once.
    """
    
    def __init__(self, A: np.ndarray, n_perturbations: int, hidden_dim: int = 64):
        super().__init__()
        self.n_genes = A.shape[0]
        self.n_perturbations = n_perturbations
        self.hidden_dim = min(hidden_dim, self.n_genes // 2)  # Adaptive sizing
        
        logger.info(f"Creating per-sample prediction model: {self.n_genes} genes, {self.n_perturbations} perturbations")
        
        # Perturbation embedding
        self.perturbation_embedding = nn.Embedding(n_perturbations, self.hidden_dim)
        
        # Pre-computed GRN layer - applies to entire gene expression vector
        A_tensor = torch.from_numpy(A.astype(NUMPY_DTYPE))
        self.grn_layer = PrecomputedGRNLayer(A_tensor, alpha=0.1, ridge_reg=1e-6)
        
        # Encoder: processes entire gene expression vector
        self.encoder = nn.Sequential(
            nn.Linear(self.n_genes, self.hidden_dim),
            nn.ReLU(),
            nn.Dropout(0.1),
            nn.Linear(self.hidden_dim, self.hidden_dim),
            nn.ReLU()
        )
        
        # Perturbation integration layer
        self.pert_integration = nn.Sequential(
            nn.Linear(self.hidden_dim * 2, self.hidden_dim),  # encoded genes + perturbation
            nn.ReLU(),
            nn.Dropout(0.1),
            nn.Linear(self.hidden_dim, self.hidden_dim),
            nn.ReLU()
        )
        
        # Decoder: outputs full gene expression vector
        self.decoder = nn.Sequential(
            nn.Linear(self.hidden_dim, self.hidden_dim),
            nn.ReLU(),
            nn.Dropout(0.1),
            nn.Linear(self.hidden_dim, self.n_genes)  # Output all genes simultaneously
        )
        
        # Initialize with small weights
        self._initialize_weights()
    
    def _initialize_weights(self):
        for module in self.modules():
            if isinstance(module, nn.Linear):
                nn.init.xavier_uniform_(module.weight, gain=0.1)
                if module.bias is not None:
                    nn.init.zeros_(module.bias)
            elif isinstance(module, nn.Embedding):
                nn.init.normal_(module.weight, mean=0.0, std=0.05)
    
    def forward(self, x: torch.Tensor, pert: torch.LongTensor) -> torch.Tensor:
        """
        Per-sample prediction: predicts all genes simultaneously for each sample.
        Single forward pass processes entire gene expression vector.
        """
        # Step 1: Apply GRN dynamics to entire gene expression vector
        x_grn = self.grn_layer(x)  # (batch_size, n_genes) -> (batch_size, n_genes)
        
        # Step 2: Encode gene expression to hidden representation  
        x_encoded = self.encoder(x_grn)  # (batch_size, n_genes) -> (batch_size, hidden_dim)
        
        # Step 3: Get perturbation embedding
        pert_emb = self.perturbation_embedding(pert)  # (batch_size,) -> (batch_size, hidden_dim)
        
        # Step 4: Combine gene and perturbation information
        combined = torch.cat([x_encoded, pert_emb], dim=1)  # (batch_size, hidden_dim*2)
        integrated = self.pert_integration(combined)  # (batch_size, hidden_dim*2) -> (batch_size, hidden_dim)
        
        # Step 5: Decode to expression changes for ALL genes simultaneously
        delta = self.decoder(integrated)  # (batch_size, hidden_dim) -> (batch_size, n_genes)
        
        # Step 6: Return baseline + predicted changes
        return x + delta  # (batch_size, n_genes)

class FastPerturbationDataset(Dataset):
    """Dataset optimized for fast loading with minimal memory footprint."""
    
    def __init__(self, X: np.ndarray, idx: np.ndarray, match_groups: np.ndarray,
                 control_map: Dict[int, int], perturbations: np.ndarray):
        self.X = X.astype(NUMPY_DTYPE)
        self.idx = idx
        self.match_groups = match_groups
        self.control_map = control_map
        self.perturbations = perturbations.astype(int)
        
        # Pre-compute control indices for speed
        self.control_indices = np.array([
            self.control_map[int(self.match_groups[i])] for i in self.idx
        ])

    def __len__(self) -> int:
        return len(self.idx)

    def __getitem__(self, i: int) -> Tuple[torch.Tensor, int, torch.Tensor]:
        idx = self.idx[i]
        y = torch.from_numpy(self.X[idx])
        x = torch.from_numpy(self.X[self.control_indices[i]])
        p = int(self.perturbations[idx])
        return x, p, y

def ultra_fast_evaluate(A, train_data_loader, test_data_loader, n_perturbations: int) -> float:
    """
    Ultra-fast evaluation with per-sample prediction (all genes simultaneously).
    """
    logger.info(f"Starting per-sample evaluation with {A.shape[0]} genes, {n_perturbations} perturbations")
    
    # Create per-sample prediction model
    model = PerSamplePredictionModel(A, n_perturbations, hidden_dim=64)
    model = model.to(DEVICE)
    
    # Optimized training settings
    optimizer = torch.optim.Adam(model.parameters(), lr=0.002, weight_decay=1e-4)
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer, patience=3, factor=0.7)
    
    best_val_loss = float('inf')
    patience_counter = 0
    max_epochs = 50  # Reduced epochs
    patience = 8
    
    # Training with early stopping
    model.train()
    for epoch in range(max_epochs):
        total_loss = 0
        n_batches = 0
        
        for x, pert, y in train_data_loader:
            x, pert, y = x.to(DEVICE), pert.to(DEVICE), y.to(DEVICE)
            
            optimizer.zero_grad()
            # Per-sample prediction: model predicts all genes for each sample
            y_hat = model(x, pert)  # (batch_size, n_genes)
            
            # Loss over all genes for all samples in batch
            loss = torch.mean((y - y_hat) ** 2)  # MSE across all genes and samples
            loss.backward()
            
            # Gradient clipping
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
            
            optimizer.step()
            total_loss += loss.item()
            n_batches += 1
        
        avg_train_loss = total_loss / n_batches
        
        # Validation
        model.eval()
        val_loss = 0
        n_val_batches = 0
        
        with torch.no_grad():
            for x, pert, y in test_data_loader:
                x, pert, y = x.to(DEVICE), pert.to(DEVICE), y.to(DEVICE)
                # Per-sample prediction for validation
                y_hat = model(x, pert)  # (batch_size, n_genes)
                loss = torch.mean((y - y_hat) ** 2)  # MSE across all genes and samples
                val_loss += loss.item()
                n_val_batches += 1
        
        avg_val_loss = val_loss / n_val_batches
        scheduler.step(avg_val_loss)
        
        # Early stopping
        if avg_val_loss < best_val_loss:
            best_val_loss = avg_val_loss
            patience_counter = 0
        else:
            patience_counter += 1
        
        if epoch % 10 == 0:
            logger.info(f"Epoch {epoch}: Train={avg_train_loss:.4f}, Val={avg_val_loss:.4f}")
        
        if patience_counter >= patience:
            logger.info(f"Early stopping at epoch {epoch}")
            break
        
        model.train()
    
    logger.info(f"Training completed. Best validation loss: {best_val_loss:.4f}")
    return best_val_loss

def create_optimized_baseline_grn(A_true, method='weight_shuffled'):
    """Create baseline with same optimizations."""
    if method == 'weight_shuffled':
        A_baseline = A_true.copy()
        non_zero_mask = A_baseline != 0
        if np.any(non_zero_mask):
            non_zero_values = A_baseline[non_zero_mask]
            np.random.shuffle(non_zero_values)
            A_baseline[non_zero_mask] = non_zero_values
        return A_baseline
    else:
        # Random permutation
        A_baseline = A_true.T.copy()
        np.random.shuffle(A_baseline)
        return A_baseline.T

def main(par):
    """
    Ultra-fast main function with all optimizations.
    """
    import sys
    import os
    import time
    
    start_time = time.time()
    
    # Import necessary modules
    sys.path.append(os.path.join(os.path.dirname(__file__), '../../utils'))
    from dataset_config import DATASET_GROUPS
    from util import manage_layer, read_prediction
    
    # Load data
    adata = ad.read_h5ad(par['evaluation_data'])
    dataset_id = adata.uns['dataset_id']
    method_id = ad.read_h5ad(par['prediction'], backed='r').uns['method_id']
    
    logger.info(f"Processing dataset: {dataset_id}, method: {method_id}")
    
    # Configure dataset-specific groups
    par['cv_groups'] = DATASET_GROUPS[dataset_id]['cv']
    par['match'] = DATASET_GROUPS[dataset_id]['match'] 
    par['loose_match'] = DATASET_GROUPS[dataset_id]['loose_match']

    # Load expression data
    layer = manage_layer(adata, par)
    X = adata.layers[layer]
    if isinstance(X, csr_matrix):
        X = X.toarray()
    X = X.astype(np.float32)
    
    # Load GRN prediction
    df_grn = read_prediction(par)
    logger.info(f"Loaded GRN with {len(df_grn)} edges")
    
    # OPTIMIZATION 1: Smart Gene Selection
    selected_genes = smart_gene_selection(df_grn, adata.var_names, max_genes=800)
    
    # Filter expression data to selected genes
    gene_mask = adata.var_names.isin(selected_genes)
    X_filtered = X[:, gene_mask]
    selected_gene_names = adata.var_names[gene_mask]
    
    logger.info(f"Reduced from {X.shape[1]} to {X_filtered.shape[1]} genes ({100*X_filtered.shape[1]/X.shape[1]:.1f}% of original)")
    
    # OPTIMIZATION 2: Create Sparse GRN Matrix
    A, final_gene_names = create_sparse_grn_matrix(df_grn, adata.var_names, selected_genes)
    
    # Re-filter expression data to match GRN genes exactly
    gene_mask_final = selected_gene_names.isin(final_gene_names)
    X_final = X_filtered[:, gene_mask_final]
    
    logger.info(f"Final gene set: {len(final_gene_names)} genes")
    logger.info(f"Final expression matrix: {X_final.shape}")
    
    # Get sample info
    cv_groups = encode_obs_cols(adata, par['cv_groups'])
    match_groups = encode_obs_cols(adata, par['match'])
    loose_match_groups = encode_obs_cols(adata, par['loose_match'])
    
    perturbations = cv_groups[0]
    cv_groups = combine_multi_index(*cv_groups)
    match_groups = combine_multi_index(*match_groups)
    loose_match_groups = combine_multi_index(*loose_match_groups)

    are_controls = adata.obs["is_control"].values.astype(bool)

    # Control matching
    try:
        control_map, match_groups = create_control_matching(are_controls, match_groups)
    except AssertionError:
        logger.info("Using loose matching for controls")
        control_map, match_groups = create_control_matching(are_controls, loose_match_groups)

    # Cross-validation with optimizations
    ss_res = 0
    ss_tot = 0
    cv = GroupKFold(n_splits=5)
    
    for i, (train_index, test_index) in enumerate(cv.split(X_final, X_final, cv_groups)):
        logger.info(f"Processing fold {i+1}/5")
        
        # Standardize (only once per fold)
        scaler = StandardScaler()
        scaler.fit(X_final[train_index, :])
        X_standardized = scaler.transform(X_final)

        # OPTIMIZATION 3: Fast Dataset Creation
        n_perturbations = int(np.max(perturbations) + 1)
        train_dataset = FastPerturbationDataset(X_standardized, train_index, match_groups, control_map, perturbations)
        test_dataset = FastPerturbationDataset(X_standardized, test_index, match_groups, control_map, perturbations)
        
        # OPTIMIZATION 4: Batch-optimized data loaders
        train_loader = DataLoader(train_dataset, batch_size=1024, shuffle=True, num_workers=0)
        test_loader = DataLoader(test_dataset, batch_size=1024, shuffle=False, num_workers=0)

        # Evaluate inferred GRN
        fold_start = time.time()
        ss_res += ultra_fast_evaluate(A, train_loader, test_loader, n_perturbations)
        
        # Evaluate baseline GRN
        A_baseline = create_optimized_baseline_grn(A)
        ss_tot += ultra_fast_evaluate(A_baseline, train_loader, test_loader, n_perturbations)
        
        fold_time = time.time() - fold_start
        logger.info(f"Fold {i+1} completed in {fold_time:.1f}s")

    # Calculate final R²
    r2 = 1 - ss_res / ss_tot if ss_tot > 1e-10 else 0.0
    
    total_time = time.time() - start_time
    logger.info(f"=== RESULTS ===")
    logger.info(f"Method: {method_id}")
    logger.info(f"Ultra-fast VC R²: {r2:.6f}")
    logger.info(f"SS_res: {ss_res:.4f}, SS_tot: {ss_tot:.4f}")
    logger.info(f"Total runtime: {total_time:.1f}s")
    logger.info(f"Genes used: {len(final_gene_names)} (vs {X.shape[1]} original)")

    results = {'vc': [float(r2)]}
    return pd.DataFrame(results)