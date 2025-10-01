import os
from typing import Tuple, Dict

import tqdm
from sklearn.model_selection import GroupKFold
from torch.utils.data import Dataset

os.environ['CUBLAS_WORKSPACE_CONFIG'] = ':4096:8' # For reproducibility purposes (on GPU)

import random
import h5py
import numpy as np
import pandas as pd
import torch
from scipy.sparse import csr_matrix
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


def combine_multi_index(*arrays) -> np.array:
    """Combine parallel label arrays into a single integer label per position."""
    A = np.stack(arrays)
    n_classes = tuple(int(A[i].max()) + 1 for i in range(A.shape[0]))
    return np.ravel_multi_index(A, dims=n_classes, order='C')


for method_name in ["collectRI", "granie", "figr", "scglue", "celloracle", "scenicplus"]:


    # Load perturbation data
    with h5py.File("../../../resources_test/grn_benchmark/evaluation_data/op_bulk.h5ad", "r") as f:

        # Get sample info
        perturbations = None
        cv_groups, match_groups, exact_match_groups = [], [], []
        for obs_name in ["perturbation", "perturbation_type", "cell_type", "donor_id", "plate_name", "row"]:
            if obs_name in f["obs"]:
                if "codes" in f["obs"][obs_name]:
                    codes = f["obs"][obs_name]["codes"][:]
                else:
                    codes = LabelEncoder().fit_transform(f["obs"][obs_name][:])
                if perturbations is None:
                    perturbations = codes
                if obs_name in {"perturbation", "cell_type"}:
                    cv_groups.append(codes)
                if obs_name not in {"row", "perturbation"}:
                    match_groups.append(codes)
                exact_match_groups.append(codes)

        # Groups used for cross-validation
        cv_groups = combine_multi_index(*cv_groups)

        # Groups used for matching with negative controls
        match_groups = combine_multi_index(*match_groups)

        # Groups used for exact matching with negative controls
        # (e.g., samples should be from the same plate row)
        exact_match_groups = combine_multi_index(*exact_match_groups)

        are_controls = f["obs"]["is_control"][:].astype(bool)
        X = f["layers"]["X_norm"]
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

        # Get gene names
        if "index" in f["var"]:
            gene_names = f["var"]["index"][:].astype(str)
        elif "Probe_id" in f["var"]:
            gene_names = f["var"]["Probe_id"][:].astype(str)
        else:
            gene_names = f["var"]["gene_name"][:].astype(str)

    gene_dict = {gene_name: i for i, gene_name in enumerate(gene_names)}

    # Load inferred GRN
    df = pd.read_csv(f"../../../output/benchmark/grn_models/{method_name}.csv")
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

    # Restrict to the 1000 most highly variable genes for speed
    std = np.std(X, axis=0)
    idx = np.argsort(std)[-1000:]
    X = X[:, idx]
    A = A[idx, :][:, idx]
    gene_names = gene_names[idx]

    # Mapping between gene expression profiles and their matched negative controls
    def create_control_matching(are_controls: np.ndarray, match_groups: np.ndarray) -> Tuple[Dict[int, int], np.ndarray]:
        control_map = {}
        for i, (is_control, group_id) in enumerate(zip(are_controls, match_groups)):
            if is_control:
                control_map[int(group_id)] = i
        for i in range(len(are_controls)):
            assert match_groups[i] in control_map
        return control_map, match_groups
    try:
        control_map, match_groups = create_control_matching(are_controls, exact_match_groups)
    except AssertionError:
        print("Failed to match with controls exactly. Using a less stringent matching instead.")
        control_map, match_groups = create_control_matching(are_controls, match_groups)


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
            self.A_signs: torch.Tensor = A_signs
            self.A_mask: torch.Tensor = (self.A_signs > 0).to(self.A_weights.dtype)
            self.I = torch.eye(self.n_genes, dtype=A_weights.dtype, device=DEVICE)
            self.signed: bool = signed
            self.inverse: bool = inverse
            self.alpha: float = alpha

        def forward(self, x: torch.Tensor) -> torch.Tensor:
            if self.signed:
                A = torch.abs(self.A_weights) * self.A_signs
            else:
                A = self.A_weights * self.A_mask
            ia = self.I - self.alpha * A.t()
            if self.inverse:
                ia = torch.linalg.inv(ia)
            return torch.mm(ia, x.t()).t()


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

            self.encoder = torch.nn.Sequential(
                GRNLayer(A_weights, A_signs, inverse=False, signed=signed),
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
                GRNLayer(A_weights, A_signs, inverse=True, signed=signed)
            )

        def forward(self, x: torch.Tensor, pert: torch.LongTensor) -> torch.Tensor:
            y = self.encoder(x)
            z = self.perturbation_embedding(pert)
            y = y + z
            y = self.decoder(y)
            return x + y


    class PerturbationDataset(Dataset):

        def __init__(
                self,
                X: np.ndarray,
                idx: np.ndarray,
                match_groups: np.ndarray,
                control_map: Dict[int, int],
                perturbations: np.ndarray
        ):
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


    def evaluate(A, train_data_loader, test_data_loader, n_perturbations: int) -> float:

        # Training
        signed = np.any(A < 0)
        model = Model(A, n_perturbations, n_hidden=64, signed=signed)
        optimizer = torch.optim.Adam(model.parameters(), lr=0.001)
        scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
            optimizer, mode='min', patience=5,
            min_lr=1e-5, cooldown=3, factor=0.8
        )
        pbar = tqdm.tqdm(range(350))
        best_val_loss = float('inf')
        best_epoch = 0
        model.train()
        for epoch in pbar:

            # Train for one epoch
            total_loss = 0
            for x, pert, y in train_data_loader:
                optimizer.zero_grad()
                y_hat = model(x, pert)
                loss = torch.mean(torch.square(y - y_hat))
                loss.backward()
                optimizer.step()
                total_loss += loss.item() * len(x)
            pbar.set_description(str(total_loss))
            scheduler.step(total_loss)

            # Evaluation
            model.eval()
            ss_res = 0
            with torch.no_grad():
                for x, pert, y in test_data_loader:
                    y_hat = model(x, pert)
                    residuals = torch.square(y - y_hat).cpu().data.numpy()
                    ss_res += np.sum(residuals)
            if ss_res < best_val_loss:
                best_val_loss = ss_res
                best_epoch = epoch
            model.train()

        model.eval()

        print(f"Best epoch: {best_epoch} ({best_val_loss})")

        return best_val_loss


    ss_res = 0
    ss_tot = 0
    cv = GroupKFold(n_splits=5, random_state=seed, shuffle=True)
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
            perturbations
        )
        train_data_loader = torch.utils.data.DataLoader(train_dataset, batch_size=512)
        test_dataset = PerturbationDataset(
            X_standardized,
            test_index,
            match_groups,
            control_map,
            perturbations
        )
        test_data_loader = torch.utils.data.DataLoader(test_dataset, batch_size=512)

        # Evaluate inferred GRN
        ss_res += evaluate(A, train_data_loader, test_data_loader, n_perturbations)

        # Evaluate baseline GRN (shuffled target genes)
        A_baseline = np.copy(A).T
        np.random.shuffle(A_baseline)
        A_baseline = A_baseline.T
        ss_tot += evaluate(A_baseline, train_data_loader, test_data_loader, n_perturbations)

    r2 = 1 - ss_res / ss_tot
    print(method_name)
    print(f"R2: {r2}")
