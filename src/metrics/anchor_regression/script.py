import os
import traceback

from sklearn.model_selection import GroupShuffleSplit, KFold
import random
import h5py
import numpy as np
import pandas as pd
import torch
import tqdm
from scipy.sparse import csr_matrix
from scipy.stats import ttest_rel, spearmanr, wilcoxon
from sklearn.linear_model import ElasticNet
from sklearn.preprocessing import StandardScaler, LabelEncoder, OneHotEncoder
import sys


# Hyper-parameters
DEVICE = 'cuda' if torch.cuda.is_available() else 'cpu'
NUMPY_DTYPE = np.float32

# For reproducibility purposes
seed = 0xCAFE
os.environ['PYTHONHASHSEED'] = str(seed)
random.seed(seed)
np.random.seed(seed)

try:
    sys.path.append(meta["resources_dir"])
except:
    meta = {
    "resources_dir":'src/metrics/anchor_regression',
    "util_dir":'src/utils'
    }
    sys.path.append(meta["resources_dir"])
    sys.path.append(meta["util_dir"])
from util import read_prediction
save_dir = 'output/anchor_regression/'


def combine_multi_index(*arrays) -> np.array:
    """Combine parallel label arrays into a single integer label per position."""
    A = np.stack(arrays)
    n_classes = tuple(int(A[i].max()) + 1 for i in range(A.shape[0]))
    return np.ravel_multi_index(A, dims=n_classes, order='C')


for dataset in ['op']:
    rr_all = {}
    # for method_id in ['pearson_corr', 'grnboost', 'ppcor', 'scenicplus', 'scenic', 'celloracle', 'scprint']:
    for method_id in ['pearson_corr']:
        ## VIASH START
        par = {
            'prediction': f'resources_test/results/{dataset}/{dataset}.{method_id}.{method_id}.prediction.h5ad',
            'evaluation_data': f'resources_test/grn_benchmark/evaluation_data/{dataset}_bulk.h5ad',
            'layer': 'lognorm',
            "max_n_links": 50000,
            'num_workers': 20,
            'tf_all': 'resources/grn_benchmark/prior/tf_all.csv',    
        }
        ## VIASH END
        # Load perturbation data
        with h5py.File(par['evaluation_data'], "r") as f:
            assert 'is_control' in f["obs"].keys(), "The evaluation dataset should contain negative controls."

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
            X = f["layers"][par['layer']]
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
            else:
                gene_names = f["var"]["gene_name"][:].astype(str)
        gene_dict = {gene_name: i for i, gene_name in enumerate(gene_names)}

        # Load inferred GRN
        df = read_prediction(par['prediction'], par)
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

            selected_diff = float(
                np.mean(np.clip((theta0[is_selected] - theta[is_selected]) / (theta0[is_selected] + eps), 0, 1)))
            unselected_diff = float(
                np.mean(np.clip((theta0[~is_selected] - theta[~is_selected]) / (theta0[~is_selected] + eps), 0, 1)))
            score = np.clip((unselected_diff - selected_diff) / (unselected_diff + selected_diff + eps), -1, 1)
            if np.isnan(score):
                scores.append(0.0)
            else:
                scores.append(score)

        score = float(np.mean(scores))
        rr_all[method_id]['final'] = score


    print(rr_all)
    import json
    with open(f'{save_dir}/results_{dataset}.json', 'w') as f:
        json.dump(rr_all, f, indent=4)
