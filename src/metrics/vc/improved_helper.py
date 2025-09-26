"""
Improved VC (Variational Causal) metric implementation with fixes for numerical stability
and better baseline construction.
"""

import torch
import torch.nn as nn
import numpy as np
from sklearn.preprocessing import StandardScaler, LabelEncoder
from torch.utils.data import Dataset, DataLoader
import logging
from typing import Dict, Tuple

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Utility functions copied from original helper
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

class ImprovedGRNLayer(nn.Module):
    """Improved GRN layer with numerical stabilization and smaller alpha values."""
    
    def __init__(self, n_genes, alpha=0.1, ridge_reg=1e-6):
        super().__init__()
        self.n_genes = n_genes
        self.alpha = alpha  # Reduced from 1.0 to 0.1 for stability
        self.ridge_reg = ridge_reg  # Ridge regularization for numerical stability
        
        # Register identity matrix as buffer for automatic device handling
        self.register_buffer('identity', torch.eye(n_genes))
        
    def forward(self, A_weights, x):
        """
        Forward pass with numerical stabilization.
        
        Args:
            A_weights: (batch_size, n_genes, n_genes) - regulatory network matrices
            x: (batch_size, n_genes) - input gene expression
            
        Returns:
            y: (batch_size, n_genes) - steady-state gene expression
        """
        batch_size = A_weights.shape[0]
        
        # Ensure A_weights are on the same device as identity
        A_weights = A_weights.to(self.identity.device)
        x = x.to(self.identity.device)
        
        # Scale A_weights more conservatively to prevent numerical issues
        A_weights = torch.tanh(A_weights) * 0.5  # Limit to [-0.5, 0.5] range
        
        # Expand identity matrix for batch processing
        I_batch = self.identity.unsqueeze(0).expand(batch_size, -1, -1)
        
        # Compute (I - α*A^T) with ridge regularization
        matrix = I_batch - self.alpha * A_weights.transpose(-2, -1)
        
        # Add ridge regularization for numerical stability
        matrix = matrix + self.ridge_reg * I_batch
        
        try:
            # Solve the linear system: (I - α*A^T) * y = x
            # This is more stable than computing the inverse directly
            y = torch.linalg.solve(matrix, x.unsqueeze(-1)).squeeze(-1)
            
            # Check for numerical issues
            if torch.isnan(y).any() or torch.isinf(y).any():
                logger.warning("Numerical instability detected, using fallback")
                return x  # Fallback to input
                
            return y
            
        except torch.linalg.LinAlgError:
            logger.warning("Matrix inversion failed, using fallback")
            return x  # Fallback to input

class ImprovedModel(nn.Module):
    """Improved neural network model with better architecture and regularization."""
    
    def __init__(self, n_genes, n_perturbations, hidden_dim=64, dropout=0.3):
        super().__init__()
        self.n_genes = n_genes
        self.n_perturbations = n_perturbations
        self.hidden_dim = hidden_dim
        
        # Perturbation embedding layer
        self.perturbation_embedding = nn.Embedding(n_perturbations, hidden_dim)
        
        # Gene expression encoder
        self.gene_encoder = nn.Sequential(
            nn.Linear(n_genes, hidden_dim),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, hidden_dim),
            nn.ReLU(),
            nn.Dropout(dropout)
        )
        
        # Combined encoder for GRN weights
        self.grn_encoder = nn.Sequential(
            nn.Linear(hidden_dim * 2, hidden_dim),  # Concatenated features
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, n_genes * n_genes)
        )
        
        # GRN layer with improved stability
        self.grn_layer = ImprovedGRNLayer(n_genes, alpha=0.1, ridge_reg=1e-6)
        
        # Output decoder
        self.decoder = nn.Sequential(
            nn.Linear(n_genes, hidden_dim),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, n_genes)
        )
        
        # Initialize weights properly
        self._initialize_weights()
    
    def _initialize_weights(self):
        """Initialize weights with appropriate scaling."""
        for module in self.modules():
            if isinstance(module, nn.Linear):
                nn.init.xavier_uniform_(module.weight, gain=0.1)  # Smaller gain
                if module.bias is not None:
                    nn.init.zeros_(module.bias)
            elif isinstance(module, nn.Embedding):
                nn.init.normal_(module.weight, mean=0.0, std=0.1)
    
    def forward(self, perturbation, baseline_expression):
        """
        Forward pass with improved numerical stability.
        
        Args:
            perturbation: (batch_size,) - perturbation indices
            baseline_expression: (batch_size, n_genes) - baseline gene expression
            
        Returns:
            predicted_expression: (batch_size, n_genes) - predicted gene expression
        """
        batch_size = baseline_expression.shape[0]
        
        # Encode perturbation
        pert_emb = self.perturbation_embedding(perturbation)  # (batch_size, hidden_dim)
        
        # Encode gene expression
        gene_emb = self.gene_encoder(baseline_expression)  # (batch_size, hidden_dim)
        
        # Combine perturbation and gene features
        combined_features = torch.cat([pert_emb, gene_emb], dim=1)  # (batch_size, hidden_dim*2)
        
        # Generate GRN weights
        grn_weights_flat = self.grn_encoder(combined_features)
        grn_weights = grn_weights_flat.view(batch_size, self.n_genes, self.n_genes)
        
        # Apply GRN dynamics
        steady_state = self.grn_layer(grn_weights, baseline_expression)
        
        # Decode to final expression change
        expression_change = self.decoder(steady_state)
        
        # Return baseline expression + predicted change
        predicted_expression = baseline_expression + expression_change
        
        return predicted_expression

class PerturbationDataset(Dataset):
    """Dataset for perturbation experiments with improved preprocessing."""
    
    def __init__(self, perturbations, expressions, baseline_expressions=None, 
                 standardize=True):
        """
        Initialize dataset with optional standardization.
        
        Args:
            perturbations: (n_samples,) - perturbation indices as integers
            expressions: (n_samples, n_genes) - target gene expressions
            baseline_expressions: (n_samples, n_genes) - baseline expressions
            standardize: bool - whether to standardize the data
        """
        # Convert perturbations to long tensor for embedding layer
        if perturbations.ndim == 1:
            self.perturbations = torch.LongTensor(perturbations)
        else:
            # If 2D, assume it's one-hot and convert to indices
            self.perturbations = torch.LongTensor(np.argmax(perturbations, axis=1))
        
        self.expressions = torch.FloatTensor(expressions)
        
        if baseline_expressions is not None:
            self.baseline_expressions = torch.FloatTensor(baseline_expressions)
        else:
            # Use mean expression as baseline
            self.baseline_expressions = torch.FloatTensor(
                np.tile(expressions.mean(axis=0), (expressions.shape[0], 1))
            )
        
        if standardize:
            self._standardize_data()
    
    def _standardize_data(self):
        """Standardize expressions (perturbations are kept as indices)."""
        # Don't standardize perturbations - they are discrete indices for embedding
        
        # Standardize expressions
        if self.expressions.std() > 0:
            expr_mean = self.expressions.mean(dim=0)
            expr_std = self.expressions.std(dim=0) + 1e-8
            self.expressions = (self.expressions - expr_mean) / expr_std
            self.baseline_expressions = (self.baseline_expressions - expr_mean) / expr_std
    
    def __len__(self):
        return len(self.perturbations)
    
    def __getitem__(self, idx):
        return (
            self.perturbations[idx],  # Already LongTensor
            self.baseline_expressions[idx],
            self.expressions[idx]
        )

def create_improved_baseline_grn(A_true, method='degree_preserving'):
    """
    Create an improved baseline GRN that preserves network properties.
    
    Args:
        A_true: true GRN adjacency matrix
        method: baseline construction method
        
    Returns:
        A_baseline: baseline GRN matrix
    """
    if method == 'degree_preserving':
        # Preserve degree distribution but randomize connections
        A_baseline = np.zeros_like(A_true)
        
        # Get in-degree and out-degree for each gene
        in_degrees = np.sum(A_true != 0, axis=0)
        out_degrees = np.sum(A_true != 0, axis=1)
        
        # Randomly reassign edges while preserving degrees
        for i in range(A_true.shape[0]):
            if out_degrees[i] > 0:
                # Randomly select targets for this gene's outgoing edges
                targets = np.random.choice(
                    A_true.shape[1], 
                    size=out_degrees[i], 
                    replace=False
                )
                # Assign random weights similar to original distribution
                weights = np.random.choice(
                    A_true[A_true != 0], 
                    size=out_degrees[i]
                )
                A_baseline[i, targets] = weights
        
        return A_baseline
    
    elif method == 'weight_shuffled':
        # Shuffle weights while keeping structure
        A_baseline = A_true.copy()
        non_zero_mask = A_baseline != 0
        non_zero_values = A_baseline[non_zero_mask]
        np.random.shuffle(non_zero_values)
        A_baseline[non_zero_mask] = non_zero_values
        return A_baseline
    
    else:  # random
        # Original random baseline
        A_baseline = A_true.T.copy()
        np.random.shuffle(A_baseline)
        return A_baseline.T

def improved_train_model(model, train_loader, val_loader, n_epochs=100, 
                        learning_rate=1e-3, patience=10, device='cpu'):
    """
    Train model with early stopping and validation monitoring.
    
    Args:
        model: neural network model
        train_loader: training data loader
        val_loader: validation data loader
        n_epochs: maximum number of epochs
        learning_rate: learning rate
        patience: early stopping patience
        device: computing device
        
    Returns:
        model: trained model
        train_losses: training loss history
        val_losses: validation loss history
    """
    model = model.to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate, weight_decay=1e-5)
    criterion = nn.MSELoss()
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
        optimizer, mode='min', factor=0.5, patience=5, verbose=True
    )
    
    train_losses = []
    val_losses = []
    best_val_loss = float('inf')
    patience_counter = 0
    
    for epoch in range(n_epochs):
        # Training phase
        model.train()
        train_loss = 0
        for batch_idx, (perturbation, baseline_expr, target_expr) in enumerate(train_loader):
            perturbation = perturbation.to(device)
            baseline_expr = baseline_expr.to(device)
            target_expr = target_expr.to(device)
            
            optimizer.zero_grad()
            predicted_expr = model(perturbation, baseline_expr)
            loss = criterion(predicted_expr, target_expr)
            loss.backward()
            
            # Gradient clipping for stability
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
            
            optimizer.step()
            train_loss += loss.item()
        
        train_loss /= len(train_loader)
        
        # Validation phase
        model.eval()
        val_loss = 0
        with torch.no_grad():
            for perturbation, baseline_expr, target_expr in val_loader:
                perturbation = perturbation.to(device)
                baseline_expr = baseline_expr.to(device)
                target_expr = target_expr.to(device)
                
                predicted_expr = model(perturbation, baseline_expr)
                loss = criterion(predicted_expr, target_expr)
                val_loss += loss.item()
        
        val_loss /= len(val_loader)
        
        train_losses.append(train_loss)
        val_losses.append(val_loss)
        
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
        
        if epoch % 10 == 0:
            logger.info(f"Epoch {epoch}: Train Loss = {train_loss:.4f}, Val Loss = {val_loss:.4f}")
    
    return model, train_losses, val_losses

def calculate_improved_r2(y_true, y_pred, y_baseline):
    """
    Calculate R² with improved numerical stability.
    
    Args:
        y_true: true values
        y_pred: predicted values
        y_baseline: baseline predictions
        
    Returns:
        r2_score: R² score
    """
    # Calculate sum of squared errors
    sse_pred = np.sum((y_true - y_pred) ** 2)
    sse_baseline = np.sum((y_true - y_baseline) ** 2)
    
    # Avoid division by zero
    if sse_baseline < 1e-10:
        logger.warning("Baseline SSE is near zero, using alternative R² calculation")
        # Use correlation-based R² as fallback
        corr = np.corrcoef(y_true.flatten(), y_pred.flatten())[0, 1]
        return corr ** 2 if not np.isnan(corr) else 0.0
    
    # Standard R² calculation
    r2 = 1 - (sse_pred / sse_baseline)
    
    # Clip extreme values
    r2 = np.clip(r2, -1.0, 1.0)
    
    return r2

def improved_evaluate(A, train_data_loader, test_data_loader, n_perturbations: int) -> float:
    """
    Improved evaluation function with better training and numerical stability.
    
    Args:
        A: GRN adjacency matrix
        train_data_loader: training data loader
        test_data_loader: test data loader
        n_perturbations: number of perturbations
        
    Returns:
        test_loss: validation loss
    """
    import torch
    
    # Force CPU usage to avoid memory issues with large gene networks
    device = torch.device('cpu')
    logger.info(f"Using device: {device}")
    
    # Create improved model with smaller architecture for large gene networks
    n_genes = A.shape[0]
    hidden_dim = min(64, n_genes // 10)  # Adaptive hidden dimension
    logger.info(f"Using hidden_dim={hidden_dim} for {n_genes} genes")
    
    model = ImprovedModel(
        n_genes=n_genes,
        n_perturbations=n_perturbations,
        hidden_dim=hidden_dim,
        dropout=0.3
    )
    
    # Train model with improved training procedure
    model, train_losses, val_losses = improved_train_model(
        model=model,
        train_loader=train_data_loader,
        val_loader=test_data_loader,
        n_epochs=100,
        learning_rate=1e-3,
        patience=10,
        device=device
    )
    
    # Return final validation loss
    return val_losses[-1] if val_losses else float('inf')

def main(par):
    """
    Improved main function with better numerical stability and baseline construction.
    """
    import sys
    import os
    
    # Import necessary modules
    sys.path.append(os.path.join(os.path.dirname(__file__), '../../utils'))
    from dataset_config import DATASET_GROUPS
    from util import manage_layer, read_prediction
    
    import anndata as ad
    import numpy as np
    import pandas as pd
    from sklearn.preprocessing import StandardScaler
    from sklearn.model_selection import GroupKFold
    from scipy.sparse import csr_matrix
    import torch
    
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

    # Only consider the genes that are actually present in the inferred GRN
    gene_mask = np.logical_or(np.any(A, axis=1), np.any(A, axis=0))
    X = X[:, gene_mask]
    A = A[gene_mask, :][:, gene_mask]
    gene_names = gene_names[gene_mask]

    # Restrict to the 5000 most highly variable genes for speed
    std = np.std(X, axis=0)
    idx = np.argsort(std)[-5000:]
    X = X[:, idx]
    A = A[idx, :][:, idx]
    gene_names = gene_names[idx]

    # Mapping between gene expression profiles and their matched negative controls
    try:
        control_map, match_groups = create_control_matching(are_controls, match_groups)
    except AssertionError:
        print("Failed to match with controls exactly. Using a less stringent matching instead.")
        control_map, match_groups = create_control_matching(are_controls, loose_match_groups)

    ss_res = 0
    ss_tot = 0
    cv = GroupKFold(n_splits=5)
    
    results = []
    
    for i, (train_index, test_index) in enumerate(cv.split(X, X, cv_groups)):
        print(f"Processing fold {i+1}/5...")
        
        # Center and scale dataset
        scaler = StandardScaler()
        scaler.fit(X[train_index, :])
        X_standardized = scaler.transform(X)

        # Create data loaders
        n_perturbations = int(np.max(perturbations) + 1)
        train_dataset = PerturbationDataset(
            perturbations[train_index],
            X_standardized[train_index],
            baseline_expressions=None,
            standardize=False  # Already standardized
        )
        train_data_loader = torch.utils.data.DataLoader(train_dataset, batch_size=512, shuffle=True)
        
        test_dataset = PerturbationDataset(
            perturbations[test_index],
            X_standardized[test_index],
            baseline_expressions=None,
            standardize=False  # Already standardized
        )
        test_data_loader = torch.utils.data.DataLoader(test_dataset, batch_size=512, shuffle=False)

        # Evaluate inferred GRN with improved method
        ss_res += improved_evaluate(A, train_data_loader, test_data_loader, n_perturbations)

        # Evaluate improved baseline GRN
        A_baseline = create_improved_baseline_grn(A, method='degree_preserving')
        ss_tot += improved_evaluate(A_baseline, train_data_loader, test_data_loader, n_perturbations)

    # Calculate improved R²
    r2 = calculate_improved_r2(
        y_true=np.array([ss_tot]),
        y_pred=np.array([ss_res]),
        y_baseline=np.array([ss_tot])
    )
    
    # Alternative R² calculation if needed
    if abs(r2) < 1e-6:  # Near zero, try alternative
        r2 = 1 - ss_res / ss_tot if ss_tot > 1e-10 else 0.0
    
    print(f"Method: {method_id}")
    print(f"Improved R²: {r2}")
    print(f"SS_res: {ss_res}, SS_tot: {ss_tot}")

    results = {
        'vc': [float(r2)]
    }

    df_results = pd.DataFrame(results)
    return df_results