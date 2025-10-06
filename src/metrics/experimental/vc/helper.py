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
from scipy.sparse import csr_matrix
from scipy.stats import wilcoxon
from sklearn.model_selection import StratifiedGroupKFold, StratifiedKFold
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


class GRNLayer(torch.nn.Module):

    def __init__(
            self,
            A_weights: torch.nn.Parameter,
            A_signs: torch.Tensor,
            signed: bool = True,
            inverse: bool = True,
            alpha: float = 0.2,
            stable: bool = True,
            bias: bool = True
    ):
        torch.nn.Module.__init__(self)
        self.n_genes: int = A_weights.size(1)
        self.A_weights: torch.nn.Parameter = A_weights
        dtype = A_weights.dtype
        device = A_weights.device
        if bias:
            self.b: torch.nn.Parameter = torch.nn.Parameter(torch.zeros((1, self.n_genes), dtype=dtype, device=device))
        else:
            self.b = None
        self.register_buffer('A_signs', A_signs.to(A_weights.device))
        self.register_buffer('A_mask', (A_signs != 0).to(self.A_weights.dtype).to(A_weights.device))
        self.register_buffer('I', torch.eye(self.n_genes, dtype=dtype, device=device))
        self.signed: bool = signed
        self.inverse: bool = inverse
        self.alpha: float = alpha
        self.stable: bool = stable

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        if self.signed:
            A = torch.abs(self.A_weights) * self.A_signs
            assert torch.any(A < 0)
        else:
            A = self.A_weights * self.A_mask
        
        if self.inverse:
            if self.stable:
                # Approximation using Neumann series
                B = GRNLayer.neumann_series(A.t(), self.alpha)
                y = torch.mm(x, B)
            else:
                # For inverse transformation, use iterative solve to avoid memory issues
                # Solve (I - alpha * A.t()) * y = x for y
                ia = self.I - self.alpha * A.t()
                
                try:
                    # Use solve instead of inversion to save memory
                    y = torch.linalg.solve(ia, x.t()).t()
                except torch.linalg.LinAlgError:
                    # Fallback: approximation using Neumann series
                    B = GRNLayer.neumann_series(A.t(), self.alpha)
                    y = torch.mm(x, B)
        else:
            # Forward transformation: apply GRN directly
            y = torch.mm(x, self.I - self.alpha * A.t())

        # Add bias term
        if self.b is not None:
            y = y + self.b

        return y

    @staticmethod
    def neumann_series(A: torch.Tensor, alpha: float, k: int = 2) -> torch.Tensor:
        """Approximate the inverse of I - A using Neumann series.

        Args:
            A: the matrix for which to invert I - A.
            k: the number of terms in the series. The higher, the more accurate.

        Returns:
            Approximated inverse of I - A.
        """
        I = torch.eye(A.shape[0], device=A.device, dtype=A.dtype)
        M = alpha * A
        term = I.clone()
        B = I.clone()
        for _ in range(k):
            term = term @ M
            B = B + term
        return B


class Model(torch.nn.Module):

    def __init__(self, A: np.ndarray, n_perturbations: int, n_hidden: int = 16, signed: bool = True):

        # n_hidden needs to be small enough to prevent the NN from arbitrarily shifting the learning task
        # from the GRN to the MLPs.

        torch.nn.Module.__init__(self)
        self.n_genes: int = A.shape[1]
        self.n_perturbations: int = n_perturbations
        self.n_hidden: int = n_hidden
        self.signed: bool = signed

        # Perturbation transformations defined in the latent space
        #self.perturbation_embedding = torch.nn.Embedding(n_perturbations, n_hidden)
        self.perturbation_embedding = torch.nn.Embedding(n_perturbations, n_hidden * n_hidden)

        # First layer: GRN-informed transformation of control expression
        A_signs = torch.from_numpy(np.sign(A).astype(NUMPY_DTYPE))
        A_weights = np.copy(A).astype(NUMPY_DTYPE)
        A_weights /= (np.sqrt(self.n_genes) * float(np.std(A_weights)))
        A_weights = torch.nn.Parameter(torch.from_numpy(A_weights))
        # Ensure A_signs is on the same device as A_weights
        A_signs = A_signs.to(A_weights.device)
        self.grn_input_layer = GRNLayer(A_weights, A_signs, inverse=False, signed=signed)
        
        # Middle layers: encode/decode between expression profiles and latent space
        self.encoder = torch.nn.Sequential(
            torch.nn.LayerNorm(self.n_genes),
            torch.nn.Dropout(p=0.4, inplace=True),
            torch.nn.PReLU(1),
            torch.nn.Linear(self.n_genes, self.n_hidden),
            torch.nn.PReLU(self.n_hidden),
            torch.nn.Linear(self.n_hidden, self.n_hidden),
            torch.nn.PReLU(self.n_hidden),
            torch.nn.Linear(self.n_hidden, self.n_hidden),
            torch.nn.PReLU(self.n_hidden),
        )
        
        self.decoder = torch.nn.Sequential(
            torch.nn.LayerNorm(self.n_hidden),
            torch.nn.Linear(self.n_hidden, self.n_hidden),
            torch.nn.PReLU(self.n_hidden),
            torch.nn.Linear(self.n_hidden, self.n_hidden),
            torch.nn.PReLU(self.n_hidden),
            torch.nn.Linear(self.n_hidden, self.n_genes),
            torch.nn.PReLU(1),
            torch.nn.Dropout(p=0.4, inplace=True),
        )
        
        # Last layer: GRN-informed transformation to final expression
        self.grn_output_layer = GRNLayer(A_weights, A_signs, inverse=True, signed=signed)

    def forward(self, x: torch.Tensor, pert: torch.LongTensor) -> torch.Tensor:

        # Encode each expression profile
        x = self.grn_input_layer(x)
        y = self.encoder(x)

        # Each perturbation is a linear transformation in the latent space.
        # Apply perturbation transform to the encoded profile.
        z = self.perturbation_embedding(pert).view(len(x), self.n_hidden, self.n_hidden)
        y = torch.einsum('ij,ijk->ik', y, z)
        #z = self.perturbation_embedding(pert)
        #y = y + z

        # Decode each expression profile
        y = self.decoder(y)
        y = self.grn_output_layer(y)
        return y

    def set_grn(self, A: np.ndarray) -> None:
        signed = np.any(A < 0)
        A_signs = torch.from_numpy(np.sign(A).astype(NUMPY_DTYPE))
        A_weights = np.copy(A).astype(NUMPY_DTYPE)
        A_weights /= (np.sqrt(self.n_genes) * float(np.std(A_weights)))
        A_weights = torch.nn.Parameter(torch.from_numpy(A_weights))
        # Ensure A_signs is on the same device as A_weights
        A_signs = A_signs.to(A_weights.device)
        self.grn_input_layer = GRNLayer(A_weights, A_signs, inverse=False, signed=self.signed)
        self.grn_output_layer = GRNLayer(A_weights, A_signs, inverse=True, signed=self.signed)


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

        # Find matched control
        group = int(self.match_groups[i])
        if group in self.control_map:
            j = int(self.control_map[group])
        else:
            group = int(self.loose_match_groups[i])
            j = int(self.loose_control_map[group])

        x = torch.from_numpy(self.X[j, :])
        p = int(self.perturbations[i])
        d_x = y - x
        return x, p, d_x


def coefficients_of_determination(y_target: np.ndarray, y_pred: np.ndarray, eps: float = 1e-20) -> np.ndarray:
    residuals = np.square(y_target - y_pred)
    ss_res = np.sum(residuals, axis=0) + eps
    mean = np.mean(y_target, axis=0)[np.newaxis, :]
    residuals = np.square(y_target - mean)
    ss_tot = np.sum(residuals, axis=0) + eps
    return 1 - ss_res / ss_tot


def evaluate(A, train_data_loader, test_data_loader, state_dict, n_perturbations: int, signed: bool = True) -> np.ndarray:
    set_seed()
    A = np.copy(A)
    model = Model(A, n_perturbations, signed=signed)
    model = model.to(DEVICE)
    model.load_state_dict(state_dict, strict=False)
    model.set_grn(A)
    optimizer = torch.optim.Adam(model.parameters(), lr=0.001, weight_decay=1e-6)
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
        optimizer, mode='min', patience=5,
        min_lr=1e-6, cooldown=3, factor=0.8
    )
    pbar = tqdm.tqdm(range(300))
    best_avg_r2, best_r2 = -np.inf, None
    model.train()
    for epoch in pbar:
        total_loss = 0
        y_target, y_pred = [], []
        for x, pert, d_x in train_data_loader:
            x, pert, d_x = x.to(DEVICE), pert.to(DEVICE), d_x.to(DEVICE)

            # Reset gradients
            optimizer.zero_grad()

            # Model now predicts full perturbed expression directly
            d_x_hat = model(x, pert)
            y_target.append(d_x.cpu().data.numpy())
            y_pred.append(d_x_hat.cpu().data.numpy())

            # Compute mean squared error
            loss = torch.mean(torch.square(d_x - d_x_hat))
            total_loss += loss.item() * len(x)

            # Compute gradients (clip them to prevent divergence)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1)

            # Update parameters
            optimizer.step()

        scheduler.step(total_loss)
        r2_train = coefficients_of_determination(np.concatenate(y_target, axis=0), np.concatenate(y_pred, axis=0))

        model.eval()
        y_target, y_pred = [], []
        with torch.no_grad():
            for x, pert, d_x in test_data_loader:
                x, pert, d_x = x.to(DEVICE), pert.to(DEVICE), d_x.to(DEVICE)
                d_x_hat = model(x, pert)
                y_target.append(d_x.cpu().data.numpy())
                y_pred.append(d_x_hat.cpu().data.numpy())
        r2_test = coefficients_of_determination(np.concatenate(y_target, axis=0), np.concatenate(y_pred, axis=0))
        avg_r2 = np.mean(r2_test)
        if avg_r2 > best_avg_r2:
            best_avg_r2 = avg_r2
            best_r2 = r2_test
        pbar.set_description(str(np.mean(r2_train)) + "  " + str(np.mean(r2_test)))
        model.train()
    model.eval()

    return best_r2


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
    
    # Validation strategy: evaluate on unseen (perturbation, cell type) pairs.
    cv_groups = combine_multi_index(cell_types, perturbations)

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

    # Add self-regulations
    np.fill_diagonal(A, 1)

    # Mapping between gene expression profiles and their matched negative controls
    control_map, _ = create_control_matching(are_controls, match_groups)
    loose_control_map, _ = create_control_matching(are_controls, loose_match_groups)

    r2 = []
    r2_baseline = []
    cv = StratifiedGroupKFold(n_splits=N_FOLDS, shuffle=True, random_state=0xCAFE)
    
    for i, (train_index, test_index) in enumerate(cv.split(X, perturbations, cv_groups)):

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
        scaler.fit(X[train_index, :])
        X_standardized = scaler.transform(X)

        # Create data loaders
        n_perturbations = int(np.max(perturbations) + 1)
        def create_data_loaders():
            train_dataset = PerturbationDataset(
                X_standardized,
                train_index,
                match_groups,
                control_map,
                loose_match_groups,
                loose_control_map,
                perturbations
            )
            train_data_loader = torch.utils.data.DataLoader(train_dataset, batch_size=512, shuffle=True)
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
            return train_data_loader, test_data_loader

        # For fair comparison, we first randomly initialize NN parameters, and use these same
        # parameters to build both models (only the GRN weights will differ).
        signed = np.any(A < 0)
        model_template = Model(A, n_perturbations, signed=signed).to(DEVICE)
        state_dict = model_template.state_dict()

        # Evaluate inferred GRN
        train_data_loader, test_data_loader = create_data_loaders()
        r2.append(evaluate(A, train_data_loader, test_data_loader, state_dict, n_perturbations, signed=signed))

        # Evaluate baseline GRN (shuffled target genes)
        train_data_loader, test_data_loader = create_data_loaders()
        r2_baseline.append(evaluate(A_baseline, train_data_loader, test_data_loader, state_dict, n_perturbations, signed=signed))

        break

    r2 = np.asarray(r2).flatten()
    r2_baseline = np.asarray(r2_baseline).flatten()
    print("Mean R2", np.mean(r2), np.mean(r2_baseline))
    if np.all(r2 == r2_baseline):
        final_score = 0
    else:
        p_value = wilcoxon(r2, r2_baseline, alternative="greater").pvalue
        final_score = -np.log10(p_value)

    print(f"Method: {method_id}")
    print(f"Final score: {final_score}")

    results = {
        'r2': [float(np.mean(r2))],
        'r2_baseline': [float(np.mean(r2_baseline))],
        'r2_diff': [float(np.mean(r2)) - float(np.mean(r2_baseline))],
        'vc': [float(final_score)]
    }

    df_results = pd.DataFrame(results)
    return df_results
