from typing import Tuple, List, Optional

import scipy
import numpy as np
from scipy.sparse import csr_matrix
from scipy.stats import spearmanr, pearsonr, wilcoxon
from sklearn.preprocessing import StandardScaler, LabelEncoder, OneHotEncoder
import pandas as pd
import anndata as ad
import warnings
warnings.filterwarnings("ignore")
from util import format_save_score, read_prediction

np.random.seed(0)
from util import manage_layer


def encode_obs_cols(adata, cols):
    """Encode observation columns to integer codes."""
    encoded = []
    for col in cols:
        if col in adata.obs:
            codes = LabelEncoder().fit_transform(adata.obs[col].values)
            encoded.append(codes)
    return encoded


def combine_multi_index(*arrays) -> np.array:
    """Combine parallel label arrays into a single integer label per position."""
    A = np.stack(arrays)
    n_classes = tuple(int(A[i].max()) + 1 for i in range(A.shape[0]))
    return np.ravel_multi_index(A, dims=n_classes, order='C')


def evaluate_grn(
        C: np.ndarray,
        A: np.ndarray,
        n_tfs: Optional[int] = None,
        max_targets_per_tf: Optional[int] = None
) -> float:
    """
    Evaluate a specific setting (recall, balanced, or precision)
    """

    # Compute the min and max correlation coefficients across groups
    C_min = np.min(C, axis=0)
    C_max = np.max(C, axis=0)
    
    # Filter TFs if specified (keep the TFs with the most target genes)
    if n_tfs is None:
        tf_idx = np.arange(A.shape[1])
    else:
        tf_idx = np.argsort(np.sum(A != 0, axis=1))[-n_tfs:]
    tf_mask = np.zeros(A.shape[1], dtype=bool)
    tf_mask[tf_idx] = True

    scores = []
    for j in range(A.shape[1]):

        if tf_mask[j]:

            # Filter TGs if specified (keep the TGs with highest regulatory links)
            if max_targets_per_tf is None:
                tg_idx = np.where(A[j, :] != 0)[0]
            else:
                tg_idx = np.argsort(np.abs(A[j, :]))[-max_targets_per_tf:]
                tg_idx = tg_idx[A[j, tg_idx] != 0]

            # Compute scores
            if len(tg_idx) == 0:
                scores.append(np.nan)
            else:
                # Intuition behind the scoring scheme:
                # our score is a weighted average of consistency scores,
                # where the weights are proportional to regulatory weight magnitudes,
                # and consistency scores are penalized by the max-min correlation spread.
                # The more consistent correlation coefficients across groups, the higher
                # the consistency scores.
                weights = np.abs(A[j, tg_idx])
                weights /= np.sum(weights)
                consistency_scores = 2 - (C_max[j, tg_idx] - C_min[j, tg_idx])
                tf_score = float(np.sum(weights * consistency_scores))
                scores.append(tf_score)

        else:
            scores.append(np.nan)

    return np.array(scores)


def evaluate_setting(C: np.ndarray, A: np.ndarray, setting_name: str, **kwargs) -> float:

    # Evaluate inferred GRN
    scores = evaluate_grn(C, A, **kwargs)

    # Create baseline GRN
    A_baseline = np.copy(A)
    for j in range(A.shape[1]):
        np.random.shuffle(A_baseline[:j, j])
        np.random.shuffle(A_baseline[j+1:, j])
    assert np.any(A_baseline != A)

    # Evaluate baseline GRN
    scores_baseline = evaluate_grn(C, A_baseline, **kwargs)

    # Skip NaNs
    mask = ~np.logical_or(np.isnan(scores), np.isnan(scores_baseline))
    scores = scores[mask]
    scores_baseline = scores_baseline[mask]

    # Compare results
    if len(scores) == 0:
        print(f"Mean corr. coef. for setting '{setting_name}': grn=nan, baseline=nan")
        return 0.0
    else:
        print(f"Mean corr. coef. for setting '{setting_name}': grn={np.mean(scores):.3f}, baseline={np.mean(scores_baseline):.3f}")
        try:
            p_value = wilcoxon(scores, scores_baseline, alternative="greater", zero_method="pratt").pvalue
            p_value = np.clip(p_value, 1e-300, 1)
            print(f"p-value={p_value}")
            #return -np.log10(p_value)
            return np.mean(np.clip(scores - scores_baseline, 0, None))
        except ValueError:
            return 0.0


def main(par):
    # Load evaluation data
    adata = ad.read_h5ad(par['evaluation_data'])
    dataset_id = adata.uns['dataset_id']
    method_id = ad.read_h5ad(par['prediction'], backed='r').uns['method_id']

    # Manage layer
    layer = manage_layer(adata, par)
    X = adata.layers[layer]
    if isinstance(X, csr_matrix):
        X = X.toarray()
    X = X.astype(np.float32)

    gene_names = adata.var_names
    gene_dict = {gene_name: i for i, gene_name in enumerate(gene_names)}

    # Get groups
    group_variables, group_names = [], []
    if ("donor_id" in adata.obs) or ("plate_name" in adata.obs):
        variables = encode_obs_cols(adata, ["donor_id", "plate_name"])
        group_variables.append(combine_multi_index(*variables))
        group_names.append("donor_id/plate_name")
    if (len(group_variables) == 0) and ("cell_type" in adata.obs):
        variables = encode_obs_cols(adata, ["cell_type"])
        group_variables.append(combine_multi_index(*variables))
        group_names.append("cell_type")
    if len(group_variables) == 0:
        variables = encode_obs_cols(adata, ["perturbation"])
        group_variables.append(combine_multi_index(*variables))
        group_names.append("perturbation")

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

    # Only consider the genes that are actually present in the inferred GRN,
    # and keep only the most-connected genes (for speed).
    gene_mask = np.logical_or(np.any(A, axis=1), np.any(A, axis=0))
    X = X[:, gene_mask]
    X = X.toarray() if isinstance(X, csr_matrix) else X
    A = A[gene_mask, :][:, gene_mask]
    gene_names = gene_names[gene_mask]

    # Remove self-regulations
    np.fill_diagonal(A, 0)
    print(f"Evaluating {X.shape[1]} genes with {np.sum(A != 0)} regulatory links")

    recall_results = []
    balanced_results = []
    precision_results = []
    for group_name, groups in zip(group_names, group_variables):
        print("\n" + "=" * 30)
        print(f"Grouping samples by {group_name}")

        # Compute correlation matrix for each group
        S, C = [], []
        for group in np.unique(groups):
            mask = (groups == group)
            if np.sum(mask) >= 3:
                assert not np.any(np.isnan(X[mask, :]))
                res = spearmanr(X[mask, :])
                S.append(res.statistic)
                C.append(res.pvalue)
        S = np.nan_to_num(S, nan=0)
        C = np.nan_to_num(C, nan=1)
        #C = -np.sign(S) * np.log10(C)
        #C = np.clip(C, -4, 4)
        C = S
        
        # 1. Recall setting: all TFs, all target genes (no filtering)
        recall_results.append(evaluate_setting(
            C, A, "recall",
            n_tfs=None,  # Use all TFs
            max_targets_per_tf=None  # Use all target genes (up to default limit)
        ))
        
        # 2. Balanced setting: top 100 TFs, top 100 target genes
        balanced_results.append(evaluate_setting(
            C, A, "balanced",
            n_tfs=100, 
            max_targets_per_tf=100
        ))
        
        # 3. Precision setting: top 40 TFs, top 40 target genes
        precision_results.append(evaluate_setting(
            C, A, "precision",
            n_tfs=40, 
            max_targets_per_tf=40
        ))

    recovery_recall = np.mean(recall_results)
    recovery_balanced = np.mean(balanced_results)
    recovery_precision = np.mean(precision_results)
    
    # Create results DataFrame with all three scores
    results = pd.DataFrame({
        'recovery_recall': [recovery_recall],
        'recovery_balanced': [recovery_balanced],
        'recovery_precision': [recovery_precision]
    })
    print(results)
    
    return results
