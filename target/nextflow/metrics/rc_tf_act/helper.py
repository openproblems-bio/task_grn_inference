from typing import Dict, List, Tuple
import os
import random
import numpy as np
import pandas as pd
import anndata as ad
import decoupler as dc
from scipy.stats import wilcoxon
from scipy.sparse import csr_matrix
import warnings
warnings.filterwarnings("ignore")

from util import read_prediction, manage_layer, create_grn_baseline
from config import DATASET_GROUPS

# Set environment variables for deterministic NumPy operations
os.environ['PYTHONHASHSEED'] = '0'
os.environ['OMP_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1'
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['NUMEXPR_NUM_THREADS'] = '1'

# Set seeds for reproducibility
SEED = 42
random.seed(SEED)
np.random.seed(SEED)


def calculate_tf_degree_centrality(net: pd.DataFrame) -> pd.DataFrame:
    """
    Calculate degree centrality for TFs (number of target genes).
    
    Parameters:
    -----------
    net : pd.DataFrame
        GRN network with columns: source, target, weight
        
    Returns:
    --------
    tf_centrality : pd.DataFrame
        DataFrame with columns: tf, degree, rank
    """
    # Count number of targets per TF (out-degree)
    tf_degrees = net.groupby('source', sort=True).size().reset_index(name='degree')  # sort=True for determinism
    tf_degrees.columns = ['tf', 'degree']
    
    # Rank TFs by degree (higher degree = lower rank number)
    # Use stable sort for deterministic ordering when degrees are equal
    tf_degrees = tf_degrees.sort_values('degree', ascending=False, kind='stable').reset_index(drop=True)
    tf_degrees['rank'] = tf_degrees.index + 1
    
    return tf_degrees


def select_top_tfs_and_edges(net: pd.DataFrame, n_tfs: int, n_edges_per_tf: int) -> pd.DataFrame:
    """
    Select top N TFs by degree centrality and top M edges for each selected TF.
    
    Parameters:
    -----------
    net : pd.DataFrame
        GRN network with columns: source, target, weight
    n_tfs : int
        Number of top TFs to select
    n_edges_per_tf : int
        Number of top edges to select per TF (by absolute weight)
        
    Returns:
    --------
    net_subset : pd.DataFrame
        Subset of network with selected TFs and edges
    """
    # Calculate TF centrality
    tf_centrality = calculate_tf_degree_centrality(net)
    
    # Select top N TFs
    top_tfs = tf_centrality.head(n_tfs)['tf'].tolist()
    
    # Filter network to only include top TFs
    net_filtered = net[net['source'].isin(top_tfs)].copy()
    
    # For each TF, select top edges by absolute weight
    net_subset = []
    for tf in top_tfs:
        tf_edges = net_filtered[net_filtered['source'] == tf].copy()
        # Sort by absolute weight (stable sort for determinism)
        tf_edges['abs_weight'] = tf_edges['weight'].abs()
        tf_edges = tf_edges.sort_values('abs_weight', ascending=False, kind='stable')
        # Take top edges
        # tf_edges = tf_edges.head(n_edges_per_tf)
        tf_edges = tf_edges.drop(columns=['abs_weight'])
        net_subset.append(tf_edges)
    
    if len(net_subset) > 0:
        net_subset = pd.concat(net_subset, ignore_index=True)
    else:
        net_subset = pd.DataFrame(columns=['source', 'target', 'weight'])
    
    return net_subset


def encode_obs_cols(adata, cols):
    """Encode observation columns to get group identifiers."""
    from sklearn.preprocessing import LabelEncoder
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


def calculate_tf_activities(X: np.ndarray, gene_names: np.ndarray, net: pd.DataFrame, min_targets: int = 3) -> Tuple[np.ndarray, np.ndarray]:
    """
    Calculate TF activities using decoupler's ULM method.
    
    Parameters:
    -----------
    X : np.ndarray
        Expression matrix (samples × genes)
    gene_names : np.ndarray
        Gene names
    net : pd.DataFrame
        GRN network with columns: source, target, weight
    min_targets : int
        Minimum number of targets required for TF
        
    Returns:
    --------
    activities : np.ndarray
        TF activities (samples × TFs)
    tf_names : np.ndarray
        Names of TFs in the same order as activity columns
    """
    # Create AnnData object
    adata = ad.AnnData(X=X)
    adata.var_names = gene_names
    adata.obs_names = [f"sample_{i}" for i in range(X.shape[0])]
    
    # Filter network to only include genes present in data
    net_filtered = net[net['target'].isin(gene_names)].copy()
    
    # Run ULM with correct API (dc.mt.ulm modifies adata in place)
    # Note: ULM is deterministic (simple linear regression), but we set seed for safety
    np.random.seed(SEED)
    dc.mt.ulm(
        adata,
        net_filtered,
        tmin=min_targets
    )
    
    # Extract TF activities (stored in obsm with key 'score_ulm')
    if 'score_ulm' in adata.obsm:
        activities = adata.obsm['score_ulm'].values
    elif 'ulm_estimate' in adata.obsm:
        activities = adata.obsm['ulm_estimate'].values
    else:
        raise ValueError(f"ULM did not produce activity estimates. Available keys: {list(adata.obsm.keys())}")
    
    # Get TF names from the obsm DataFrame
    if 'score_ulm' in adata.obsm:
        tf_names = adata.obsm['score_ulm'].columns.to_numpy()
    elif 'ulm_estimate' in adata.obsm:
        tf_names = adata.obsm['ulm_estimate'].columns.to_numpy()
    else:
        tf_names = np.array([])
    
    return activities, tf_names


def calculate_activity_dispersion(activities: np.ndarray, groups: np.ndarray, tf_names: np.ndarray) -> Tuple[np.ndarray, np.ndarray, pd.DataFrame]:
    """
    Calculate dispersion of TF activities within each group using Mean Absolute Deviation (MAD).
    
    MAD is more robust than CV for small sample sizes (2-3 donors).
    For each TF and group: MAD = mean(|activity_i - median(activities)|)
    
    Parameters:
    -----------
    activities : np.ndarray
        TF activities (samples × TFs)
    groups : np.ndarray
        Group identifiers for each sample (perturbation × celltype)
    tf_names : np.ndarray
        Names of TFs corresponding to activity columns
        
    Returns:
    --------
    dispersions : np.ndarray
        Mean Absolute Deviation per TF (averaged across groups)
    consistency_scores : np.ndarray
        Consistency score per TF (1 - normalized_dispersion)
    tf_scores_df : pd.DataFrame
        DataFrame with columns: tf, mad, consistency
    """
    n_tfs = activities.shape[1]
    tf_dispersions = []
    
    # Sort groups for deterministic iteration order
    unique_groups = np.sort(np.unique(groups))
    
    for tf_idx in range(n_tfs):
        tf_activities = activities[:, tf_idx]
        group_dispersions = []
        
        for group_id in unique_groups:
            group_mask = (groups == group_id)
            group_activities = tf_activities[group_mask]
            
            # Skip groups with too few samples
            if len(group_activities) < 2:
                continue
            
            # Calculate Mean Absolute Deviation from median
            median_activity = np.median(group_activities)
            mad = np.mean(np.abs(group_activities - median_activity))
            group_dispersions.append(mad)
        
        # Average dispersion across groups for this TF
        if len(group_dispersions) > 0:
            tf_dispersions.append(np.mean(group_dispersions))
        else:
            tf_dispersions.append(np.nan)
    
    tf_dispersions = np.array(tf_dispersions)
    
    # Convert dispersion to consistency score (lower dispersion = higher consistency)
    # Normalize dispersions to [0, 1] range using 95th percentile
    valid_dispersions = tf_dispersions[~np.isnan(tf_dispersions)]
    if len(valid_dispersions) > 0:
        # Use 95th percentile as upper bound to avoid outliers
        max_dispersion = np.percentile(valid_dispersions, 95)
        
        # Normalize: 0 dispersion = 1.0 consistency, max = 0.0 consistency
        normalized_dispersions = np.clip(tf_dispersions / max_dispersion, 0, 1)
        consistency_scores = 1.0 - normalized_dispersions
    else:
        consistency_scores = np.full_like(tf_dispersions, np.nan)
    
    # Create DataFrame with TF-level scores
    tf_scores_df = pd.DataFrame({
        'tf': tf_names,
        'mad': tf_dispersions,
        'consistency': consistency_scores
    })
    
    return tf_dispersions, consistency_scores, tf_scores_df


def calculate_standardized_scores(
    tf_scores_df: pd.DataFrame,
    n_tfs: int
) -> float:
    """
    Calculate mean consistency score for top N TFs.
    
    If GRN has fewer TFs than n_tfs, pad with zeros and take mean.
    
    Parameters:
    -----------
    tf_scores_df : pd.DataFrame
        DataFrame with columns: tf, mad, consistency
    n_tfs : int
        Number of TFs to consider (top N by consistency)
        
    Returns:
    --------
    score : float
        Mean consistency score (padded with zeros if needed)
    """
    # Sort TFs by consistency score (descending)
    # Use stable sort and secondary sort key (tf name) for determinism
    tf_scores_sorted = tf_scores_df.sort_values(
        ['consistency', 'tf'], 
        ascending=[False, True],  # Higher consistency first, then alphabetically by TF
        na_position='last',
        kind='stable'
    )
    
    # Get top N TFs (or all if fewer than N)
    top_tfs = tf_scores_sorted.head(n_tfs)
    
    # Get consistency scores
    consistency_values = top_tfs['consistency'].values
    
    # Pad with zeros if needed
    if len(consistency_values) < n_tfs:
        n_missing = n_tfs - len(consistency_values)
        consistency_values = np.concatenate([consistency_values, np.zeros(n_missing)])
    
    # Remove NaN values and replace with 0
    consistency_values = np.nan_to_num(consistency_values, nan=0.0)
    
    # Return mean
    return np.mean(consistency_values)


def evaluate_grn(
    X: np.ndarray,
    gene_names: np.ndarray,
    net: pd.DataFrame,
    groups: np.ndarray,
    min_targets: int = 5
) -> Tuple[np.ndarray, np.ndarray, pd.DataFrame]:
    """
    Evaluate GRN based on TF activity consistency across replicas.
    
    Parameters:
    -----------
    X : np.ndarray
        Expression matrix (samples × genes)
    gene_names : np.ndarray
        Gene names
    net : pd.DataFrame
        GRN network
    groups : np.ndarray
        Group identifiers (perturbation × celltype combinations)
    min_targets : int
        Minimum targets per TF
        
    Returns:
    --------
    dispersions : np.ndarray
        TF activity dispersions
    consistency_scores : np.ndarray
        TF activity consistency scores
    tf_scores_df : pd.DataFrame
        DataFrame with TF-level scores
    """
    # Calculate TF activities
    print(f"Calculating TF activities for {len(net['source'].unique())} TFs...")
    activities, tf_names = calculate_tf_activities(X, gene_names, net, min_targets=min_targets)
    
    # Calculate dispersion
    print(f"Calculating activity dispersion across {len(np.unique(groups))} groups...")
    dispersions, consistency_scores, tf_scores_df = calculate_activity_dispersion(activities, groups, tf_names)
    
    # Print statistics
    valid_consistency = consistency_scores[~np.isnan(consistency_scores)]
    print(f"Consistency - Mean: {np.mean(valid_consistency):.3f} ± {np.std(valid_consistency):.3f}, Median: {np.median(valid_consistency):.3f}")
    
    return dispersions, consistency_scores, tf_scores_df


def main(par):
    # Load data
    adata = ad.read_h5ad(par['evaluation_data'])
    dataset_id = adata.uns['dataset_id']
    method_id = ad.read_h5ad(par['prediction'], backed='r').uns['method_id']
    
    # Manage layer
    layer = manage_layer(adata, par)
    X = adata.layers[layer]
    if isinstance(X, csr_matrix):
        X = X.toarray()
    X = X.astype(np.float32)
    
    gene_names = adata.var_names.to_numpy()
    
    # Get grouping variables (perturbation × celltype × donor)
    grouping_vars = DATASET_GROUPS[dataset_id]['rc_tf_ac']
    print(f"Grouping by: {grouping_vars}")
    
    # Create groups for perturbation × celltype combinations
    # We want to measure consistency across donors within each perturbation × celltype
    group_vars = encode_obs_cols(adata, grouping_vars)
    groups = combine_multi_index(*group_vars)
    
    # Load predicted GRN
    print("Loading predicted GRN...")
    net_full = read_prediction(par)
    
    # Configuration for three scores: precision, balanced, recall
    # Each uses different TF set size and edges per TF
    score_configs = [
        # {'name': 'precision', 'n_tfs': 10, 'n_edges_per_tf': 20},
        {'name': 'balanced', 'n_tfs': 100, 'n_edges_per_tf': 100},
        {'name': 'recall', 'n_tfs': 300, 'n_edges_per_tf': 1000}
    ]
    
    results = {}
    
    for config in score_configs:
        score_name = config['name']
        n_tfs = config['n_tfs']
        n_edges_per_tf = config['n_edges_per_tf']
        
        print("\n" + "="*60)
        print(f"Evaluating {score_name.upper()}: Top {n_tfs} TFs (with top {n_edges_per_tf} edges each)")
        print("="*60)
        
        # Select top TFs and edges
        net = select_top_tfs_and_edges(net_full, n_tfs, n_edges_per_tf)
        print(f"Selected {len(net)} edges across {len(net['source'].unique())} TFs")
        
        # Filter to minimum targets
        min_targets = par.get('min_targets', 5)
        tf_counts = net['source'].value_counts(sort=True)  # sort=True for determinism
        tfs_to_keep = tf_counts[tf_counts >= min_targets].index
        net = net[net['source'].isin(tfs_to_keep)]
        print(f"After filtering (min_targets={min_targets}): {len(net)} edges across {len(net['source'].unique())} TFs")
        
        if len(net) == 0:
            print(f"Warning: No TFs with at least {min_targets} targets")
            results[f'rc_tf_act_{score_name}'] = 0.0
            continue
        
        # Evaluate GRN
        dispersions, consistency_scores, tf_scores_df = evaluate_grn(X, gene_names, net, groups, min_targets)
        
        # Calculate score (mean consistency with padding)
        score = calculate_standardized_scores(tf_scores_df, n_tfs)
        
        print(f"\n{score_name.upper()} Score: {score:.4f}")
        results[f'rc_tf_act_{score_name}'] = score
    
    results_df = pd.DataFrame([results])
    
    print("\n" + "="*60)
    print("Final Summary:")
    print("="*60)
    print(results_df.to_string(index=False))
    
    return results_df
