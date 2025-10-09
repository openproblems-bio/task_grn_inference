from typing import Tuple, List, Optional

import tqdm
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
        max_targets_per_tf: Optional[int] = None,
        tf_tg: bool = True,
        tg_tg: bool = True,
        signed: bool = True
) -> float:
    """
    Evaluate a specific setting (recall, balanced, or precision)
    """

    # Compute the min and max correlation coefficients across groups.
    # If evaluation is signed, then the min-max spread occuring in the negative (-1, 0) range
    # counts 3 times more if the inferred regulatory edge is positive, and vice versa.
    C_min = np.min(C, axis=0)
    C_max = np.max(C, axis=0)

    #C_min = np.quantile(C, 0.2, axis=0)
    #C_max = np.quantile(C, 0.8, axis=0)
    if signed:
        is_positive = (A > 0)
        is_negative = (A < 0)
        #pos_spread = np.clip(C_max, 0, 1) - np.clip(C_min, 0, 1)
        #neg_spread = np.clip(C_max, -1, 0) - np.clip(C_min, -1, 0)
        #spread = is_positive * (pos_spread + 3 * neg_spread)
        #spread += is_negative * (3 * pos_spread + neg_spread)
        #consistency_scores = 1 - 0.25 * spread
        consistency_scores = is_positive * C_min + is_negative * (-C_max)
        #C_mean = np.mean(C, axis=0)
        #consistency_scores = is_positive * C_mean + is_negative * (-C_mean)
        #consistency_scores = 1 - 0.5 * (C_max - C_min)
    else:
        consistency_scores = 1 - 0.5 * (C_max - C_min)
    #assert np.all(consistency_scores >= 0)
    #assert np.all(consistency_scores <= 1)
    
    # Filter TFs if specified (keep the TFs with the most target genes)
    if n_tfs is None:
        tf_idx = np.arange(A.shape[1])
    else:
        tf_idx = np.argsort(np.sum(A != 0, axis=1))[-n_tfs:]
    tf_mask = np.zeros(A.shape[1], dtype=bool)
    tf_mask[tf_idx] = True

    scores = []
    for j in tqdm.tqdm(range(A.shape[1])):

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
                overall_score = 0.0
                n_terms = 0
                if tf_tg:
                    overall_score += float(np.sum(weights * consistency_scores[j, tg_idx]))
                    n_terms += 1
                if tg_tg:
                    S = consistency_scores[tg_idx, :][:, tg_idx]
                    np.fill_diagonal(S, 0)
                    overall_score += float(np.sum(weights * np.sum(S, axis=0) / (len(S) - 1)))
                    n_terms += 1
                if n_terms > 0:
                    overall_score /= n_terms
                scores.append(overall_score)

        else:
            scores.append(np.nan)

    return np.array(scores)


def evaluate_setting(C: np.ndarray, A: np.ndarray, setting_name: str, **kwargs) -> float:

    # Whether or not to take into account the regulatory modes (enhancer/inhibitor)
    signed = kwargs.get("signed", True)

    # Evaluate inferred GRN
    scores = evaluate_grn(C, A, **kwargs)

    # Create baseline model: for each TG, shuffle the TFs
    A_baseline = np.copy(A)
    if signed:
        A_baseline *= (2 * (np.random.randint(0, 2) - 0.5))
    tf_mask = np.any(A_baseline != 0, axis=1)
    for j in range(A.shape[1]):
        mask = np.copy(tf_mask)
        mask[j] = False
        if np.any(mask):
            values = np.copy(A_baseline[mask, j])
            np.random.shuffle(values)
            A_baseline[mask, j] = values
    # assert np.any(A_baseline != A)

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
            return -np.log10(p_value)
            #return np.mean(np.clip(scores - scores_baseline, 0, None))
        except ValueError:
            return 0.0


def spearman_corrcoef(X: np.ndarray) -> np.ndarray:
    eps = np.random.uniform(0, 1e-30, size=X.shape)  # To break ties
    #X = np.argsort(X + eps, axis=0)
    C = np.corrcoef(X.T)
    np.fill_diagonal(C, 0)
    return C


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
    if len(group_variables) == 0:
        group_variables.append(np.random.randint(0, 6, size=len(X)))
        group_names.append("random_split")

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

    # Whether or not to take into account the regulatory modes (enhancer/inhibitor)
    signed = np.any(A < 0)

    # Remove self-regulations
    np.fill_diagonal(A, 0)
    print(f"Evaluating {X.shape[1]} genes with {np.sum(A != 0)} regulatory links")

    tftg_results = []
    tgtg_results = []
    tftg_tgtg_results = []
    for group_name, groups in zip(group_names, group_variables):
        print("\n" + "=" * 30)
        print(f"Grouping samples by {group_name}")

        # Compute correlation matrix for each group
        C = []
        for group in np.unique(groups):
            print(f"Computing correlations for group {group}")
            mask = (groups == group)
            if np.sum(mask) >= 3:
                assert not np.any(np.isnan(X[mask, :]))
                #C.append(np.corrcoef(X[mask, :].T))
                C.append(spearman_corrcoef(X[mask, :]))
        C = np.asarray(C)

        tftg_results.append(evaluate_setting(
            C, A, "tftg",
            tf_tg=True,
            tg_tg=False,
            n_tfs=100,
            max_targets_per_tf=100,
            signed=signed
        ))

        tgtg_results.append(evaluate_setting(
            C, A, "tgtg",
            tf_tg=False,
            tg_tg=True,
            n_tfs=40, 
            max_targets_per_tf=300,
            signed=signed
        ))
        
        tftg_tgtg_results.append(evaluate_setting(
            C, A, "tftg+tgtg",
            tf_tg=True,
            tg_tg=True,
            n_tfs=100, 
            max_targets_per_tf=300,
            signed=signed
        ))

    tftg_score = np.mean(tftg_results)
    tgtg_score = np.mean(tgtg_results)
    tftg_tgtg_score = np.mean(tftg_tgtg_results)
    final_score = (tftg_score + tgtg_score + tftg_tgtg_score) / 3.0
    
    # Create results DataFrame with all three scores
    results = pd.DataFrame({
        'tftg': [tftg_score],
        'tgtg': [tgtg_score],
        'tftg_tgtg': [tftg_tgtg_score],
        'final_score': [final_score]
    })
    print(results)
    
    return results
