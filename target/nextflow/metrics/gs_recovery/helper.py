

import math
import pandas as pd
import numpy as np
import anndata as ad
from tqdm import tqdm
from scipy.stats import fisher_exact, ttest_1samp
from typing import Dict, List, Tuple, Set
import sys
import os
import requests
from multiprocessing import Pool, cpu_count
from functools import partial

from util import read_prediction, read_gmt

import decoupler as dc


# For reproducibility
SEED = 42
np.random.seed(SEED)


def get_enrichr_library(library_name: str) -> Dict[str, List[str]]:
    """
    Download gene sets from Enrichr.
    
    Parameters
    ----------
    library_name : str
        Name of the Enrichr library (e.g., 'MSigDB_Hallmark_2020', 'KEGG_2021_Human')
    
    Returns
    -------
    Dict[str, List[str]]
        Dictionary mapping pathway names to gene lists
    """
    url = f'https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName={library_name}'
    response = requests.get(url)
    
    if response.status_code != 200:
        raise ValueError(f"Failed to download {library_name}: {response.status_code}")
    
    gene_sets = {}
    for line in response.text.strip().split('\n'):
        if not line:
            continue
        parts = line.split('\t')
        if len(parts) < 3:
            continue
        term = parts[0]
        genes = [g for g in parts[2:] if g]  # Filter empty strings
        if genes:  # Only add if there are genes
            gene_sets[term] = genes
    
    return gene_sets


def convert_enrichr_to_pathways(gene_sets: Dict[str, List[str]], 
                                min_size: int = 5, 
                                max_size: int = 500) -> Dict[str, Set[str]]:
    """
    Convert Enrichr gene sets to pathway dictionary with size filtering.
    
    Parameters
    ----------
    gene_sets : Dict[str, List[str]]
        Gene sets from Enrichr
    min_size : int
        Minimum pathway size
    max_size : int
        Maximum pathway size
    
    Returns
    -------
    Dict[str, Set[str]]
        Filtered pathway dictionary
    """
    pathways = {}
    for name, genes in gene_sets.items():
        gene_set = set(genes)
        if min_size <= len(gene_set) <= max_size:
            pathways[name] = gene_set
    return pathways


def get_canonical_pathways(pathway_file: str) -> pd.DataFrame:
    """
    Load canonical pathways from GMT file and return as DataFrame.
    
    Parameters
    ----------
    pathway_file : str
        Path to GMT file (e.g., MSigDB Hallmark pathways)
    
    Returns
    -------
    pd.DataFrame
        DataFrame with 'gene' as index and 'pathway' column containing pathway names
    """
    genesets_all = read_gmt(pathway_file)
    genesets_all = {key: gs['genes'] for key, gs in genesets_all.items()}

    # Create gene-to-pathway mappings
    gene_to_pathway_list = [
        (gene, pathway)
        for pathway, genes in genesets_all.items()
        for gene in genes
    ]

    df_pathway = pd.DataFrame(gene_to_pathway_list, columns=["gene", "pathway"])
    df_pathway = df_pathway.set_index("gene")
    
    # Clean pathway names
    df_pathway['pathway'] = df_pathway['pathway'].str.replace('HALLMARK_', '')
    df_pathway['pathway'] = df_pathway['pathway'].str.replace('_', ' ')
    df_pathway['pathway'] = df_pathway['pathway'].str.title()

    return df_pathway


def load_pathways_as_sets(pathway_file: str, 
                          min_size: int = 5, 
                          max_size: int = 500) -> Dict[str, Set[str]]:
    """
    Load pathways as dictionary of gene sets with size filtering.
    
    Parameters
    ----------
    pathway_file : str
        Path to GMT file
    min_size : int
        Minimum pathway size
    max_size : int
        Maximum pathway size
    
    Returns
    -------
    Dict[str, Set[str]]
        Dictionary mapping pathway names to sets of genes
    """
    genesets = read_gmt(pathway_file)
    pathways = {}
    
    for pathway_name, pathway_data in genesets.items():
        genes = set(pathway_data['genes'])
        if min_size <= len(genes) <= max_size:
            # Clean pathway name
            clean_name = pathway_name.replace('HALLMARK_', '').replace('_', ' ').title()
            pathways[clean_name] = genes
    
    return pathways


def fishers_exact_test(gene_set: Set[str], 
                       pathway_genes: Set[str], 
                       all_genes: Set[str]) -> Tuple[float, float]:
    """
    Perform Fisher's exact test for pathway enrichment.
    
    Contingency table:
                    In Pathway    Not in Pathway
    In Gene Set         a               b
    Not in Gene Set     c               d
    
    Parameters
    ----------
    gene_set : Set[str]
        Set of genes to test (e.g., TF targets)
    pathway_genes : Set[str]
        Set of genes in the pathway
    all_genes : Set[str]
        Background set of all genes
    
    Returns
    -------
    Tuple[float, float]
        (odds_ratio, p_value)
    """
    # Intersection
    a = len(gene_set & pathway_genes)  # In both
    b = len(gene_set - pathway_genes)  # In gene set but not pathway
    c = len(pathway_genes - gene_set)  # In pathway but not gene set
    d = len(all_genes - gene_set - pathway_genes)  # In neither
    
    # Fisher's exact test (one-sided, greater enrichment)
    odds_ratio, p_value = fisher_exact([[a, b], [c, d]], alternative='greater')
    
    return odds_ratio, p_value


def multiple_testing_correction(p_values: List[float], method: str = 'fdr_bh') -> np.ndarray:
    """
    Apply multiple testing correction.
    
    Parameters
    ----------
    p_values : List[float]
        List of p-values
    method : str
        Correction method ('fdr_bh' for Benjamini-Hochberg)
    
    Returns
    -------
    np.ndarray
        Corrected p-values
    """
    p_values = np.array(p_values)
    n = len(p_values)
    
    if method == 'fdr_bh':
        # Benjamini-Hochberg FDR correction
        ranked_p_values = p_values.argsort()
        fdr_values = np.zeros(n)
        
        for i, idx in enumerate(ranked_p_values):
            fdr_values[idx] = p_values[idx] * n / (i + 1)
        
        # Ensure monotonicity
        for i in range(n - 2, -1, -1):
            if fdr_values[ranked_p_values[i]] > fdr_values[ranked_p_values[i + 1]]:
                fdr_values[ranked_p_values[i]] = fdr_values[ranked_p_values[i + 1]]
        
        return np.minimum(fdr_values, 1.0)
    else:
        raise ValueError(f"Unknown correction method: {method}")


def calculate_pathway_enrichment(targets: Set[str],
                                 pathways: Dict[str, Set[str]],
                                 all_genes: Set[str],
                                 fdr_threshold: float = 0.05) -> pd.DataFrame:
    """
    Calculate pathway enrichment for a set of target genes.
    
    Parameters
    ----------
    targets : Set[str]
        Set of target genes
    pathways : Dict[str, Set[str]]
        Dictionary of pathway name to gene sets
    all_genes : Set[str]
        Background gene set
    fdr_threshold : float
        FDR threshold for significance
    
    Returns
    -------
    pd.DataFrame
        Enrichment results with columns: pathway, overlap, p_value, fdr, odds_ratio
    """
    results = []
    
    for pathway_name, pathway_genes in pathways.items():
        # Only test pathways with at least one overlap
        overlap = targets & pathway_genes
        if len(overlap) == 0:
            continue
        
        odds_ratio, p_value = fishers_exact_test(targets, pathway_genes, all_genes)
        
        results.append({
            'pathway': pathway_name,
            'overlap': len(overlap),
            'overlap_genes': ','.join(sorted(overlap)),
            'pathway_size': len(pathway_genes),
            'p_value': p_value,
            'odds_ratio': odds_ratio
        })
    
    if len(results) == 0:
        return pd.DataFrame(columns=['pathway', 'overlap', 'p_value', 'fdr', 'odds_ratio'])
    
    df = pd.DataFrame(results)
    
    # Multiple testing correction
    df['fdr'] = multiple_testing_correction(df['p_value'].values)
    
    # Sort by p-value
    df = df.sort_values('p_value')
    
    return df


def calculate_tf_pathway_score(tf: str,
                               prediction: pd.DataFrame,
                               pathways: Dict[str, Set[str]],
                               all_genes: Set[str],
                               fdr_threshold: float = 0.05,
                               min_targets: int = 10,
                               max_targets: int = 100) -> Dict:
    """
    Calculate pathway fidelity score for a single TF.
    
    Uses edge weight thresholding: if TF has >100 targets, select top 100 by absolute weight.
    This handles network density variation naturally.
    
    Parameters
    ----------
    tf : str
        TF name
    prediction : pd.DataFrame
        Predicted GRN with columns: source, target, weight
    pathways : Dict[str, Set[str]]
        Pathways dictionary
    all_genes : Set[str]
        Background genes
    fdr_threshold : float
        FDR threshold for significance
    min_targets : int
        Minimum number of targets to evaluate
    max_targets : int
        Maximum number of targets to use (selects top K by edge weight)
    
    Returns
    -------
    Dict
        Score dictionary with keys: tf, score_raw, score_norm, n_pathways, specificity, etc.
    """
    # Check if TF is in prediction
    if tf not in prediction['source'].unique():
        return {
            'tf': tf,
            'is_predicted': False,
            'n_targets': 0,
            'n_pathways_enriched': 0,
            'score_raw': 0.0,
            'score_norm': 0.0,
            'specificity': 0.0,
            'best_pathway': None,
            'best_pvalue': 1.0
        }
    
    # Get all edges for this TF
    tf_edges = prediction[prediction['source'] == tf].copy()
    
    # If TF has more than max_targets, select top edges by absolute weight
    if len(tf_edges) > max_targets:
        tf_edges['abs_weight'] = tf_edges['weight'].abs()
        tf_edges = tf_edges.nlargest(max_targets, 'abs_weight')
    
    # Extract targets
    targets = set(tf_edges['target'].astype(str))
    
    # Filter to genes in background
    targets = targets & all_genes
    n_targets = len(targets)
    
    # Skip if too few targets
    if n_targets < min_targets:
        return {
            'tf': tf,
            'is_predicted': True,
            'n_targets': n_targets,
            'n_pathways_enriched': 0,
            'score_raw': 0.0,
            'score_norm': 0.0,
            'specificity': 0.0,
            'best_pathway': None,
            'best_pvalue': 1.0,
            'skip_reason': 'too_few_targets'
        }
    
    # Calculate enrichment
    enrichment_df = calculate_pathway_enrichment(targets, pathways, all_genes, fdr_threshold)
    
    # Get significant pathways
    sig_pathways = enrichment_df[enrichment_df['fdr'] < fdr_threshold]
    n_sig = len(sig_pathways)
    
    if n_sig == 0:
        return {
            'tf': tf,
            'is_predicted': True,
            'n_targets': n_targets,
            'n_pathways_enriched': 0,
            'score_raw': 0.0,
            'score_norm': 0.0,
            'specificity': 0.0,
            'best_pathway': None,
            'best_pvalue': enrichment_df['p_value'].min() if len(enrichment_df) > 0 else 1.0
        }
    
    # Calculate scores
    best_pvalue = sig_pathways['p_value'].iloc[0]
    best_pathway = sig_pathways['pathway'].iloc[0]
    
    # Raw score: -log10(p-value) of best pathway
    score_raw = -np.log10(best_pvalue) if best_pvalue > 0 else 50.0
    
    # Normalized score: account for random expectation
    # Expected p-value by chance for k targets and N pathways
    n_pathways_tested = len(enrichment_df)
    expected_pvalue = min(1.0, n_pathways_tested * fdr_threshold)
    expected_score = -np.log10(expected_pvalue) if expected_pvalue > 0 else 0.0
    
    # Normalize: (observed - expected) / observed
    score_norm = (score_raw - expected_score) / score_raw if score_raw > 0 else 0.0
    score_norm = max(0.0, score_norm)  # Ensure non-negative
    
    # Specificity: penalize TFs with too many enriched pathways
    # specificity = 1 / (1 + log(n_sig))
    specificity = 1.0 / (1.0 + np.log1p(n_sig))
    
    return {
        'tf': tf,
        'is_predicted': True,
        'n_targets': n_targets,
        'n_pathways_enriched': n_sig,
        'score_raw': score_raw,
        'score_norm': score_norm,
        'specificity': specificity,
        'best_pathway': best_pathway,
        'best_pvalue': best_pvalue
    }


def identify_active_pathways_ulm(
    evaluation_data: ad.AnnData,
    pathways: Dict[str, Set[str]],
    prediction: pd.DataFrame,
    min_targets: int = 10,
    activity_threshold: float = 0.0,
    pvalue_threshold: float = 0.01,
    baseline_method: str = 'zero_centered'
) -> Dict[str, Dict]:
    """
    Identify active pathways in the dataset using ULM (Univariate Linear Model).
    
    Parameters
    ----------
    evaluation_data : ad.AnnData
        Expression data (bulk samples × genes)
    pathways : Dict[str, Set[str]]
        Dictionary of pathway name to gene sets
    prediction : pd.DataFrame
        GRN prediction (used to format network for ULM)
    min_targets : int
        Minimum number of genes per pathway for ULM
    activity_threshold : float
        Minimum mean activity score to consider pathway active
    pvalue_threshold : float
        P-value threshold for significance
    baseline_method : str
        Method for baseline calculation:
        - 'zero_centered': Test if activity > 0 (default)
        - 'permutation': Compare to permuted gene labels (slower)
        - 'random_genesets': Compare to random gene sets of same size
    
    Returns
    -------
    Dict[str, Dict]
        Dictionary of active pathways with statistics:
        {
            'Apoptosis': {
                'mean_activity': 0.45,
                'std_activity': 0.12,
                'pvalue': 0.001,
                'is_active': True,
                'n_samples_positive': 8,
                'total_samples': 10
            },
            ...
        }
    """
    if dc is None:
        raise ImportError("decoupler package required for ULM. Install with: pip install decoupler")
    
    print("  [ULM] Identifying active pathways in dataset...")
    
    # Convert pathways to network format for decoupler
    # Format: source (pathway), target (gene), weight (1.0)
    net_list = []
    for pathway_name, genes in pathways.items():
        for gene in genes:
            net_list.append({
                'source': pathway_name,
                'target': gene,
                'weight': 1.0
            })
    
    net_df = pd.DataFrame(net_list)
    print(f"  [ULM] Network for ULM: {len(net_df)} pathway-gene relationships")
    
    # Filter network to genes present in data
    genes_in_data = set(evaluation_data.var_names)
    net_filtered = net_df[net_df['target'].isin(genes_in_data)].copy()
    
    # Filter pathways with too few genes
    pathway_sizes = net_filtered.groupby('source').size()
    valid_pathways = pathway_sizes[pathway_sizes >= min_targets].index
    net_filtered = net_filtered[net_filtered['source'].isin(valid_pathways)]
    
    print(f"  [ULM] Filtered network: {len(net_filtered)} relationships, {len(valid_pathways)} pathways")
    
    # Read expression data
    adata = evaluation_data.to_memory() if hasattr(evaluation_data, 'to_memory') else evaluation_data
    
    # Run ULM to get pathway activities using the exact method you provided
    print(f"  [ULM] Running ULM on {adata.n_obs} samples...")
    dc.mt.ulm(
        adata,
        net_filtered
    )
    
    # Extract TF activities (stored in obsm with key 'score_ulm' or 'ulm_estimate')
    if 'score_ulm' in adata.obsm:
        activities = adata.obsm['score_ulm']
    elif 'ulm_estimate' in adata.obsm:
        activities = adata.obsm['ulm_estimate']
    else:
        raise ValueError(f"ULM did not produce activity estimates. Available keys: {list(adata.obsm.keys())}")
    
    # Get pathway names
    if 'score_ulm' in adata.uns:
        pathway_names = adata.uns['score_ulm']['names']
    elif 'ulm_estimate' in adata.uns:
        pathway_names = adata.uns['ulm_estimate']['names']
    else:
        # Fallback: use unique sources from network
        pathway_names = net_filtered['source'].unique()
    
    # Convert to DataFrame
    activities_df = pd.DataFrame(
        activities,
        index=adata.obs_names,
        columns=pathway_names
    )
    
    print(f"  [ULM] Activity matrix: {activities_df.shape} (samples × pathways)")
    
    # Determine which pathways are active
    active_pathways = {}
    
    for pathway in activities_df.columns:
        activity_scores = activities_df[pathway].values
        
        # Calculate statistics
        mean_activity = activity_scores.mean()
        std_activity = activity_scores.std()
        n_positive = (activity_scores > activity_threshold).sum()
        total_samples = len(activity_scores)
        
        # Determine if pathway is active based on baseline method
        if baseline_method == 'zero_centered':
            # One-sample t-test: is mean significantly > 0?
            t_stat, p_value = ttest_1samp(activity_scores, 0.0, alternative='greater')
            is_active = (p_value < pvalue_threshold) and (mean_activity > activity_threshold)
            
        elif baseline_method == 'permutation':
            # Permutation test (expensive, not implemented yet)
            # TODO: Implement permutation baseline
            raise NotImplementedError("Permutation baseline not yet implemented")
            
        elif baseline_method == 'random_genesets':
            # Compare to random gene sets (expensive, not implemented yet)
            # TODO: Implement random geneset baseline
            raise NotImplementedError("Random geneset baseline not yet implemented")
        else:
            raise ValueError(f"Unknown baseline_method: {baseline_method}")
        
        active_pathways[pathway] = {
            'mean_activity': mean_activity,
            'std_activity': std_activity,
            'pvalue': p_value,
            'is_active': is_active,
            'n_samples_positive': n_positive,
            'total_samples': total_samples
        }
    
    n_active = sum(1 for p in active_pathways.values() if p['is_active'])
    print(f"  [ULM] Active pathways: {n_active} / {len(active_pathways)}")
    
    return active_pathways


def _test_pathway_enrichment_worker(pathway_item: Tuple[str, Set[str]], 
                                     tfs_data: Dict[str, Set[str]], 
                                     all_genes: Set[str],
                                     fdr_threshold: float,
                                     n_pathways: int) -> Tuple[str, Dict]:
    """
    Worker function for parallel pathway enrichment testing.
    
    Parameters
    ----------
    pathway_item : Tuple[str, Set[str]]
        (pathway_name, pathway_genes)
    tfs_data : Dict[str, Set[str]]
        Dictionary mapping TF names to their target gene sets
    all_genes : Set[str]
        Background genes
    fdr_threshold : float
        FDR threshold
    n_pathways : int
        Total number of pathways for Bonferroni correction
    
    Returns
    -------
    Tuple[str, Dict]
        (pathway_name, enrichment_result)
    """
    pathway_name, pathway_genes = pathway_item
    
    min_pvalue = 1.0
    enriched_tfs = []
    
    for tf, tf_targets in tfs_data.items():
        if len(tf_targets) < 10:  # Skip TFs with too few targets
            continue
        
        # Fisher's exact test
        odds_ratio, p_value = fishers_exact_test(tf_targets, pathway_genes, all_genes)
        
        if p_value < min_pvalue:
            min_pvalue = p_value
        
        if p_value < fdr_threshold / n_pathways:  # Bonferroni correction
            enriched_tfs.append((tf, p_value))
    
    result = {
        'is_enriched': len(enriched_tfs) > 0,
        'min_pvalue': min_pvalue,
        'enriched_tfs': enriched_tfs,
        'n_tfs': len(enriched_tfs)
    }
    
    return pathway_name, result


def calculate_gene_set_recovery(
    prediction: pd.DataFrame,
    pathways: Dict[str, Set[str]],
    active_pathways: Dict[str, Dict],
    all_genes: Set[str],
    fdr_threshold: float = 0.01,
    max_workers: int = None
) -> Dict:
    """
    Calculate gene set recovery metrics (precision, recall, F1).
    
    For each pathway:
    - TP: Pathway is enriched in GRN AND active in dataset
    - FP: Pathway is enriched in GRN but NOT active in dataset
    - FN: Pathway is active in dataset but NOT enriched in GRN
    
    Parameters
    ----------
    prediction : pd.DataFrame
        Predicted GRN
    pathways : Dict[str, Set[str]]
        Pathway gene sets
    active_pathways : Dict[str, Dict]
        Active pathways from ULM analysis
    all_genes : Set[str]
        Background genes
    fdr_threshold : float
        FDR threshold for enrichment
    max_workers : int, optional
        Number of parallel workers (default: use all CPUs)
    
    Returns
    -------
    Dict
        {
            'TP': int,
            'FP': int,
            'FN': int,
            'precision': float,
            'recall': float,
            'f1': float,
            'n_active_pathways': int,
            'n_enriched_pathways': int,
            'pathway_details': List[Dict]
        }
    """
    if max_workers is None:
        max_workers = cpu_count()
    
    print(f"\n  [Gene Set Recovery] Testing pathway enrichment in GRN (using {max_workers} workers)...")
    
    # Get all TFs in GRN and prepare their target sets
    print("  [Gene Set Recovery] Preparing TF target sets...")
    tfs_data = {}
    for tf in prediction['source'].unique():
        tf_targets = set(prediction[prediction['source'] == tf]['target'])
        tf_targets = tf_targets & all_genes
        if len(tf_targets) >= 10:  # Only include TFs with enough targets
            tfs_data[tf] = tf_targets
    
    print(f"  [Gene Set Recovery] Testing {len(pathways)} pathways with {len(tfs_data)} TFs...")
    
    # Filter pathways to only those in active_pathways
    pathways_to_test = [(name, genes) for name, genes in pathways.items() 
                        if name in active_pathways]
    
    # Parallel pathway enrichment testing
    pathway_enrichment = {}
    
    if len(pathways_to_test) > 10 and max_workers > 1:
        # Use parallel processing for large pathway sets
        worker_fn = partial(_test_pathway_enrichment_worker,
                           tfs_data=tfs_data,
                           all_genes=all_genes,
                           fdr_threshold=fdr_threshold,
                           n_pathways=len(pathways))
        
        with Pool(processes=max_workers) as pool:
            results = list(tqdm(
                pool.imap(worker_fn, pathways_to_test),
                total=len(pathways_to_test),
                desc="  Testing pathways"
            ))
        
        pathway_enrichment = dict(results)
    else:
        # Serial processing for small sets
        for pathway_name, pathway_genes in tqdm(pathways_to_test, desc="  Testing pathways"):
            _, result = _test_pathway_enrichment_worker(
                (pathway_name, pathway_genes),
                tfs_data=tfs_data,
                all_genes=all_genes,
                fdr_threshold=fdr_threshold,
                n_pathways=len(pathways)
            )
            pathway_enrichment[pathway_name] = result
    
    # Calculate TP, FP, FN
    TP = 0
    FP = 0
    FN = 0
    TN = 0
    
    pathway_details = []
    
    for pathway_name in active_pathways.keys():
        is_active = active_pathways[pathway_name]['is_active']
        is_enriched = pathway_enrichment.get(pathway_name, {}).get('is_enriched', False)
        
        if is_enriched and is_active:
            TP += 1
            category = 'TP'
        elif is_enriched and not is_active:
            FP += 1
            category = 'FP'
        elif not is_enriched and is_active:
            FN += 1
            category = 'FN'
        else:
            TN += 1
            category = 'TN'
        
        pathway_details.append({
            'pathway': pathway_name,
            'is_active': is_active,
            'is_enriched': is_enriched,
            'category': category,
            'mean_activity': active_pathways[pathway_name]['mean_activity'],
            'activity_pvalue': active_pathways[pathway_name]['pvalue'],
            'enrichment_pvalue': pathway_enrichment.get(pathway_name, {}).get('min_pvalue', 1.0),
            'n_tfs_enriched': pathway_enrichment.get(pathway_name, {}).get('n_tfs', 0)
        })
    
    # Calculate metrics
    precision = TP / (TP + FP) if (TP + FP) > 0 else 0.0
    recall = TP / (TP + FN) if (TP + FN) > 0 else 0.0
    f1 = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0.0
    
    
    return {
        'TP': TP,
        'FP': FP,
        'FN': FN,
        'TN': TN,
        'precision': precision,
        'recall': recall,
        'f1': f1,
        'n_active_pathways': sum(1 for p in active_pathways.values() if p['is_active']),
        'n_enriched_pathways': sum(1 for p in pathway_enrichment.values() if p['is_enriched']),
        'pathway_details': pathway_details
    }


def main(par: dict) -> pd.DataFrame:
    """
    Main function to calculate annotation metrics.
    
    If pathway_files is provided (list of paths), calculate metrics for each.
    Otherwise, use single pathway_file.
    """
    
    # Load data
    print("\n[1/5] Loading data...")
    evaluation_data = ad.read_h5ad(par['evaluation_data'], backed='r')
    all_genes = set(evaluation_data.var_names.tolist())
    prediction = read_prediction(par)
    tf_all = np.loadtxt(par['tf_all'], dtype=str, delimiter=',', skiprows=1)
    pathway_files_dict = par['pathway_files']
    min_pw_size = par.get('min_pathway_size', 5)
    max_pw_size = par.get('max_pathway_size', 500)
    
    all_results = []
    
    for geneset_name, pathway_file in pathway_files_dict.items():
        pathways = load_pathways_as_sets(pathway_file, min_pw_size, max_pw_size)
        # Parameters
        fdr_threshold = par.get('fdr_threshold', 0.05)
        min_targets = par.get('min_targets', 10)
        max_targets = par.get('max_targets', 100)
        print(f"\n  [Gene Set Recovery] Running for {geneset_name}...")
        
        # Identify active pathways using ULM
        active_pathways = identify_active_pathways_ulm(
            evaluation_data=evaluation_data,
            pathways=pathways,
            prediction=prediction,
            min_targets=min_targets,
            activity_threshold=par.get('ulm_activity_threshold', 0.0),
            pvalue_threshold=par.get('ulm_pvalue_threshold', 0.01),
            baseline_method=par.get('ulm_baseline_method', 'zero_centered')
        )
        
        # Calculate gene set recovery
        max_workers = par.get('max_workers', None)
        gs_results = calculate_gene_set_recovery(
            prediction=prediction,
            pathways=pathways,
            active_pathways=active_pathways,
            all_genes=all_genes,
            fdr_threshold=fdr_threshold,
            max_workers=max_workers
        )
        
        # Store results with gene set prefix
        result_dict = {
            'geneset_name': geneset_name,
            'precision': gs_results['precision'],
            'recall': gs_results['recall'],
            'f1': gs_results['f1'],
            'n_active_pathways': gs_results['n_active_pathways'],
            'n_enriched_pathways': gs_results['n_enriched_pathways']
        }
        
        all_results.append(result_dict)
       
    
    final_dict = {}
    for result in all_results:
        geneset_name = result['geneset_name']
        final_dict[f'{geneset_name}_gs_precision'] = result['precision']
        final_dict[f'{geneset_name}_gs_recall'] = result['recall']
        final_dict[f'{geneset_name}_gs_f1'] = result['f1']
        final_dict[f'{geneset_name}_gs_n_active'] = result['n_active_pathways']
    
    # Calculate mean across all gene sets
    if all_results:
        final_dict['gs_precision'] = np.mean([r['precision'] for r in all_results])
        final_dict['gs_recall'] = np.mean([r['recall'] for r in all_results])
        final_dict['gs_f1'] = np.mean([r['f1'] for r in all_results])
        final_dict['gs_n_active'] = np.mean([r['n_active_pathways'] for r in all_results])
    
    summary_df = pd.DataFrame([final_dict])
    print(summary_df)
    return summary_df