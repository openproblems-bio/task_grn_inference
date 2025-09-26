import os
import random
from typing import Iterable, Tuple

import h5py
import numpy as np
import scipy.stats
import pandas as pd
import torch
import tqdm
from scipy.linalg import solve
from scipy.sparse import csr_matrix
from sklearn.preprocessing import StandardScaler, LabelEncoder, OneHotEncoder

# Hyper-parameters
DEVICE = 'cuda' if torch.cuda.is_available() else 'cpu'
NUMPY_DTYPE = np.float32

# For reproducibility purposes
seed = 0xCAFE
os.environ['PYTHONHASHSEED'] = str(seed)
random.seed(seed)
np.random.seed(seed)


for method_name in ["collectRI", "granie", "figr", "scglue", "celloracle", "scenicplus"]:

    # Load perturbation data
    with h5py.File("../../../resources_test/grn_benchmark/evaluation_data/op_bulk.h5ad", "r") as f:

        # Get anchor variables
        anchor_variables = []
        for obs_name in ["donor_id", "plate_name", "cell_type"]:
            if obs_name in f["obs"]:
                if "codes" in f["obs"][obs_name]:
                    codes = f["obs"][obs_name]["codes"][:]
                else:
                    codes = LabelEncoder().fit_transform(f["obs"][obs_name][:])
                if (obs_name == "cell_type") and (len(anchor_variables) == 0):
                    anchor_variables.append(codes)
                elif obs_name != "cell_type":
                    anchor_variables.append(codes)

        # Load expression data
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
    A = np.zeros((len(gene_names), len(gene_names)), dtype=X.dtype)
    for source, target, weight in zip(sources, targets, weights):
        if (source in gene_dict) and (target in gene_dict):
            i = gene_dict[source]
            j = gene_dict[target]
            A[i, j] = float(weight)

    # Only consider the most highly-variable genes
    gene_mask = np.zeros(len(gene_names), dtype=bool)
    idx = np.argsort(np.std(X, axis=0))[-1000:]
    gene_mask[idx] = True
    X = X[:, gene_mask]
    X = X.toarray() if isinstance(X, csr_matrix) else X
    A = A[gene_mask, :][:, gene_mask]
    gene_names = gene_names[gene_mask]

    # Center and scale dataset
    scaler = StandardScaler()
    scaler.fit(X)
    X = scaler.transform(X)

    # One-hot-encode anchor variables
    anchor_variables = np.asarray(anchor_variables).T
    Z = OneHotEncoder(sparse_output=False, dtype=NUMPY_DTYPE).fit_transform(anchor_variables)


    def anchor_regression(
            X: np.ndarray,
            Z: np.ndarray,
            Y: np.ndarray,
            l2_reg: float = 1e-3,
            anchor_strength: float = 1.0
    ) -> np.ndarray:
        """Anchor regression.

        Args:
            X: input features, of shape (n, d).
            Z: anchor variables. The model is required to be invariant to
                these environmental variables. Shape (n, z).
            Y: predicted variables, of shape (n, u).

        Returns:
            Inferred parameters, of shape (d, u).
        """

        # Projector onto the anchor subspace
        pi = Z @ np.linalg.pinv(Z.T @ Z) @ Z.T

        # Whitening
        W = np.eye(len(Z)) + (np.sqrt(1.0 + anchor_strength) - 1.0) * pi
        X_t = W @ X
        Y_t = W @ Y

        # Ridge regression on the whitened data
        theta = np.linalg.solve(X_t.T @ X_t + l2_reg * np.eye(X_t.shape[1]), X_t.T @ Y_t)

        return theta


    eps = 1e-10
    scores = []
    for j in tqdm.tqdm(range(X.shape[1])):

        is_selected = np.array(A[:, j] != 0)
        if (not np.any(is_selected)) or np.all(is_selected):
            scores.append(0.0)
            continue

        mask = np.ones(X.shape[1], dtype=bool)
        mask[j] = False

        theta0 = anchor_regression(X[:, mask], Z, X[:, j], anchor_strength=0)
        theta0 = np.abs(theta0)

        theta = anchor_regression(X[:, mask], Z, X[:, j], anchor_strength=10.0)
        theta = np.abs(theta)

        is_selected = is_selected[mask]

        selected_diff = float(np.mean(np.clip((theta0[is_selected] - theta[is_selected]) / (theta0[is_selected] + eps), 0, 1)))
        unselected_diff = float(np.mean(np.clip((theta0[~is_selected] - theta[~is_selected]) / (theta0[~is_selected] + eps), 0, 1)))
        score = np.clip((unselected_diff - selected_diff) / (unselected_diff + selected_diff + eps), -1, 1)
        if np.isnan(score):
            scores.append(0.0)
        else:
            scores.append(score)

    #print(np.mean(np.clip(differences, 0, 1)))
    print(method_name, np.mean(scores))


"""
granie 0.004012749918621728
collectRI 0.07245154180677192
figr 0.12552857501981485
scglue 0.14638436014327466
scenicplus 0.15850517365008884
celloracle 0.19987462678804796
"""
