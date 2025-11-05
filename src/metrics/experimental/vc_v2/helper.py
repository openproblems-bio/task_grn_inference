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
from config import DATASET_GROUPS
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
    
    # Create control mapping
    control_map = {}
    for i, (is_control, group_id) in enumerate(zip(are_controls, match_groups)):
        if is_control:
            control_map[int(group_id)] = i
    
    # If no controls were mapped (shouldn't happen but safety check), 
    # map group 0 to the first control
    if not control_map and len(control_indices) > 0:
        control_map[0] = control_indices[0]

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


class SparseGRNLayer(torch.nn.Module):
    """
    Efficient sparse linear layer using PyTorch's scatter_add for fast computation.
    Maps gene expression to TF activities based on GRN connections.
    """
    
    def __init__(self, grn_connections: torch.Tensor, n_genes: int, n_tfs: int):
        """
        Args:
            grn_connections: [n_connections, 2] tensor with [source_gene_idx, tf_idx] pairs
            n_genes: number of input genes
            n_tfs: number of output TFs
        """
        super().__init__()
        self.n_genes = n_genes
        self.n_tfs = n_tfs
        
        # Store connections - genes regulate TFs
        self.register_buffer('gene_indices', grn_connections[:, 0].long())  # Source gene indices
        self.register_buffer('tf_indices', grn_connections[:, 1].long())    # Target TF indices
        
        # Non-trainable weights for each connection (fixed GRN structure)
        self.register_buffer('weights', torch.randn(grn_connections.size(0)) * 0.1)
        
    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """
        Args:
            x: (batch_size, n_genes) input gene expression
        Returns:
            (batch_size, n_tfs) TF activities
        """
        batch_size = x.size(0)
        
        # Gather gene expressions for all connections: (batch_size, n_connections)
        gene_vals = x[:, self.gene_indices]  # Shape: (batch_size, n_connections)
        
        # Apply weights: (batch_size, n_connections)
        weighted_vals = gene_vals * self.weights  # Broadcast weights (non-trainable)
        
        # Initialize output tensor
        tf_output = torch.zeros(batch_size, self.n_tfs, dtype=x.dtype, device=x.device)
        
        # Use scatter_add to efficiently sum contributions to each TF
        # Expand tf_indices to match batch dimension
        tf_indices_expanded = self.tf_indices.unsqueeze(0).expand(batch_size, -1)  # (batch_size, n_connections)
        
        # Scatter add: for each TF, sum all weighted gene contributions
        tf_output.scatter_add_(1, tf_indices_expanded, weighted_vals)
        
        return tf_output


class Model(torch.nn.Module):
    def __init__(self, net: pd.DataFrame, gene_names: np.ndarray, n_perturbations: int, n_hidden: int = 64):
        """
        Args:
            net: DataFrame with columns ['source', 'target', 'weight'] 
            gene_names: array of gene names (expression data genes)
            n_perturbations: number of perturbation types
            n_hidden: hidden layer size
        """
        super().__init__()
        self.n_genes = len(gene_names)
        self.n_perturbations = n_perturbations
        self.n_hidden = n_hidden
        
        # Create gene name to index mapping
        gene_to_idx = {gene: idx for idx, gene in enumerate(gene_names)}
        
        # Extract all TFs (sources) from the network - don't filter by gene_names
        tf_genes = set(net['source'].unique())
        tf_list = sorted(list(tf_genes))
        self.n_tfs = len(tf_list)
        
        print(f"Identified {self.n_tfs} transcription factors")
        print(f"Using {self.n_genes} genes from expression data")
        
        # Create TF index mapping
        tf_to_idx = {tf: idx for idx, tf in enumerate(tf_list)}
        
        connections = []
        edge_weights = []
        for _, row in net.iterrows():
            source, target, weight = row['source'], row['target'], row['weight']
            # If target gene is in expression data, it can contribute to source TF activity
            if target in gene_to_idx and source in tf_to_idx:
                gene_idx = gene_to_idx[target]  # Target gene in expression data
                tf_idx = tf_to_idx[source]      # Source TF that regulates this gene
                connections.append([gene_idx, tf_idx])
                edge_weights.append(float(weight))
        
        if len(connections) == 0:
            raise ValueError("No valid GRN connections found!")
            
        connections = torch.tensor(connections, dtype=torch.long)
        edge_weights = torch.tensor(edge_weights, dtype=torch.float32)
        
        print(f"Created {len(connections)} gene->TF connections")
        print(f"Average edge weight: {edge_weights.mean().item():.4f}")
        
        # Create sparse GRN layer: genes -> TFs
        self.grn_layer = SparseGRNLayer(connections, self.n_genes, self.n_tfs)
        
        # Initialize weights with GRN edge weights (non-trainable)
        self.grn_layer.weights.data = edge_weights
        
        self.perturbation_embedding = torch.nn.Embedding(n_perturbations, self.n_hidden)
        
        # Processing layers
        self.encoder = torch.nn.Sequential(
            torch.nn.PReLU(1),
            torch.nn.Linear(self.n_tfs, self.n_hidden),
            torch.nn.PReLU(1),
            torch.nn.Linear(self.n_hidden, self.n_hidden),
            torch.nn.PReLU(1),
        )
        
        self.decoder = torch.nn.Sequential(
            torch.nn.Linear(self.n_hidden, self.n_hidden),
            torch.nn.PReLU(1),
            torch.nn.Linear(self.n_hidden, self.n_genes),  # Output back to gene space
        )

    def forward(self, x: torch.Tensor, pert: torch.LongTensor) -> torch.Tensor:
        # 1. Map gene expression to TF activities using sparse GRN layer
        tf_activities = self.grn_layer(x)
        
        # 2. Process TF activities through encoder
        encoded = self.encoder(tf_activities)
        
        # 3. Add perturbation embedding
        pert_embedding = self.perturbation_embedding(pert)
        encoded = encoded + pert_embedding
        
        # 4. Decode back to gene expression space
        gene_expression = self.decoder(encoded)
        
        return gene_expression


class BaselineModel(torch.nn.Module):
    """
    Baseline model without GRN layer - uses fully connected layers instead.
    This serves as a comparison to show the value of the GRN structure.
    """
    def __init__(self, n_genes: int, n_perturbations: int, n_hidden: int = 64):
        """
        Args:
            n_genes: number of genes (input/output dimension)
            n_perturbations: number of perturbation types
            n_hidden: hidden layer size
        """
        super().__init__()
        self.n_genes = n_genes
        self.n_perturbations = n_perturbations
        self.n_hidden = n_hidden
        
        self.perturbation_embedding = torch.nn.Embedding(n_perturbations, self.n_hidden)
        
        # Fully connected encoder (no GRN structure)
        self.encoder = torch.nn.Sequential(
            torch.nn.Linear(self.n_genes, self.n_hidden),
            torch.nn.PReLU(1),
            torch.nn.Linear(self.n_hidden, self.n_hidden),
            torch.nn.PReLU(1),
            torch.nn.Linear(self.n_hidden, self.n_hidden),
            torch.nn.PReLU(1),
        )
        
        self.decoder = torch.nn.Sequential(
            torch.nn.Linear(self.n_hidden, self.n_hidden),
            torch.nn.PReLU(1),
            torch.nn.Linear(self.n_hidden, self.n_genes),  # Output back to gene space
        )

    def forward(self, x: torch.Tensor, pert: torch.LongTensor) -> torch.Tensor:
        # 1. Process gene expression through fully connected encoder
        encoded = self.encoder(x)
        
        # 2. Add perturbation embedding
        pert_embedding = self.perturbation_embedding(pert)
        encoded = encoded + pert_embedding
        
        # 3. Decode back to gene expression space
        gene_expression = self.decoder(encoded)
        
        return gene_expression


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
        elif int(self.loose_match_groups[i]) in self.loose_control_map:
            group = int(self.loose_match_groups[i])
            j = int(self.loose_control_map[group])
        else:
            # Fallback: use any available control sample
            # This handles cases where no matching control exists (e.g., single control scenarios)
            available_controls = list(self.control_map.values()) + list(self.loose_control_map.values())
            if available_controls:
                j = available_controls[0]  # Use first available control
            else:
                raise ValueError("No control samples available for matching!")

        x = torch.from_numpy(self.X[j, :])
        p = int(self.perturbations[i])
        return x, p, y



def evaluate_model(model, train_data_loader, test_data_loader, model_name: str = "Model") -> Tuple[float, float]:
    """
    Generic evaluation function that works with any model.
    """
    model = model.to(DEVICE)
    optimizer = torch.optim.Adam(model.parameters(), lr=0.001, weight_decay=1e-5)
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
        optimizer, mode='min', patience=5,
        min_lr=1e-5, cooldown=3, factor=0.8
    )
    pbar = tqdm.tqdm(range(100), desc=f"Training {model_name}")  # Reduced epochs for faster testing
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
        pbar.set_description(f"{model_name} - Loss: {total_loss:.4f}")
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

    # Final evaluation
    model.eval()

    ss_tot = 0

    with torch.no_grad():
        for x, pert, y in test_data_loader:
            x, pert, y = x.to(DEVICE), pert.to(DEVICE), y.to(DEVICE)
            y_hat = model(x, pert)

            residuals = torch.square(y - torch.mean(y, dim=0).unsqueeze(0)).cpu().data.numpy()
            ss_tot += np.sum(residuals, axis=0)

    # print(f"{model_name} - Best epoch: {best_epoch} (Loss: {best_val_loss:.4f})")
    return best_ss_res, ss_tot


def evaluate(net: pd.DataFrame, gene_names: np.ndarray, train_data_loader, test_data_loader, n_perturbations: int) -> Tuple[Tuple[float, float], Tuple[float, float]]:
    """
    Evaluate both GRN-based model and baseline fully-connected model.
    Returns: ((grn_ss_res, grn_ss_tot), (baseline_ss_res, baseline_ss_tot))
    """
    print("Evaluating GRN-based model...")
    grn_model = Model(net, gene_names, n_perturbations, n_hidden=32)
    grn_results = evaluate_model(grn_model, train_data_loader, test_data_loader, "GRN Model")
    
    print("\nEvaluating baseline fully-connected model...")
    baseline_model = BaselineModel(len(gene_names), n_perturbations, n_hidden=32)
    baseline_results = evaluate_model(baseline_model, train_data_loader, test_data_loader, "Baseline Model")
    
    return grn_results, baseline_results



def main(par):
    """
    Main evaluation function that compares GRN-based model performance against a baseline.
    
    The evaluation now includes:
    1. GRN Model: Uses sparse GRN layer to map genes to TF activities
    2. Baseline Model: Fully connected network without GRN structure
    3. Comparison metrics to show the value of incorporating GRN structure
    
    Returns results with:
    - vc: GRN model performance
    - vc_baseline: Baseline model performance  
    - vc_improvement: Improvement of GRN over baseline
    """
    # Load evaluation data
    adata = ad.read_h5ad(par['evaluation_data'])
    assert 'is_control' in adata.obs.columns, "'is_control' column is required in the dataset for perturbation evaluation"
    assert adata.obs['is_control'].sum() > 0, "'is_control' column must contain at least one True value for control samples"
    dataset_id = adata.uns['dataset_id']
    method_id = ad.read_h5ad(par['prediction'], backed='r').uns['method_id']
    
    # Configure dataset-specific groups
    par['cv_groups'] = DATASET_GROUPS[dataset_id]['cv']
    par['match'] = DATASET_GROUPS[dataset_id]['match'] 
    par['loose_match'] = DATASET_GROUPS[dataset_id]['loose_match']

    # X and gene_names will be defined later after network loading

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
    common_genes = set(net['target']).intersection(set(adata.var_names))
    net = net[net['target'].isin(common_genes)]
    adata = adata[:, adata.var_names.isin(common_genes)].copy()
    print(f"Loaded GRN with {len(net)} edges")
    print(f"Unique sources (TFs): {len(net['source'].unique())}")
    print(f"Unique targets: {len(net['target'].unique())}")
    
    # Update X and gene_names to match current adata (after any filtering that might have occurred)
    layer = manage_layer(adata, par)
    X = adata.layers[layer]
    if isinstance(X, csr_matrix):
        X = X.toarray()
    X = X.astype(np.float32)
    
    gene_names = adata.var_names.to_numpy()
    gene_dict = {gene_name: i for i, gene_name in enumerate(gene_names)}
    print(f"Using {len(gene_names)} genes from expression data")

    # Mapping between gene expression profiles and their matched negative controls
    control_map, _ = create_control_matching(are_controls, match_groups)
    loose_control_map, _ = create_control_matching(are_controls, loose_match_groups)

    # Initialize accumulators for both models
    grn_ss_res = 0
    grn_ss_tot = 0
    baseline_ss_res = 0
    baseline_ss_tot = 0
    cv = GroupKFold(n_splits=5)
    
    results = []
    
    for i, (train_index, test_index) in enumerate(cv.split(X, X, cv_groups)):
        # print(f"\nCross-validation fold {i+1}/5")
        
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

        # Evaluate both GRN and baseline models
        grn_res, baseline_res = evaluate(net, gene_names, train_data_loader, test_data_loader, n_perturbations)
        
        # Accumulate results
        grn_ss_res = grn_ss_res + grn_res[0]
        grn_ss_tot = grn_ss_tot + grn_res[1]
        baseline_ss_res = baseline_ss_res + baseline_res[0]
        baseline_ss_tot = baseline_ss_tot + baseline_res[1]

    # Calculate R2 scores for both models
    grn_r2 = 1 - grn_ss_res / grn_ss_tot
    baseline_r2 = 1 - baseline_ss_res / baseline_ss_tot

    grn_final_score = np.mean(np.clip(grn_r2, 0, 1))
    baseline_final_score = np.mean(np.clip(baseline_r2, 0, 1))
    
    print(f"\n=== FINAL RESULTS ===")
    print(f"Method: {method_id}")
    print(f"GRN Model R2: {grn_final_score:.4f}")
    print(f"Baseline Model R2: {baseline_final_score:.4f}")
    print(f"GRN vs Baseline Improvement: {grn_final_score - baseline_final_score:.4f}")
    print(f"Relative Improvement: {((grn_final_score - baseline_final_score) / baseline_final_score * 100):.2f}%")

    results = {
        'vc': [float(grn_final_score)],
        'vc_baseline': [float(baseline_final_score)],
        'vc_improvement': [float(grn_final_score - baseline_final_score)]
    }

    df_results = pd.DataFrame(results)
    return df_results