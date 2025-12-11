"""
Analyze robustness of vc metric using degree-preserved permutations at different intensities.

This script:
1. Loads a GRN prediction
2. Creates degree-preserved baseline using create_grn_baseline
3. Mixes original and baseline networks at different intensities (0%, 20%, 50%, 100%)
4. Evaluates vc metric for each mixing level
5. Analyzes how the metric changes with increasing permutation
"""

import os
import sys
import pandas as pd
import numpy as np
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns

# Add paths
sys.path.append('src/metrics/vc')
sys.path.append('src/utils')

from helper import main as main_vc
from util import process_links, parse_args, create_grn_baseline


def convert_edgelist_to_matrix(prediction: pd.DataFrame, genes: list) -> np.ndarray:
    """
    Convert edge list to adjacency matrix.
    
    Parameters:
    -----------
    prediction : pd.DataFrame
        GRN with columns: source, target, weight
    genes : list
        List of all genes
        
    Returns:
    --------
    matrix : np.ndarray
        Adjacency matrix (genes x genes)
    """
    n_genes = len(genes)
    gene_to_idx = {gene: i for i, gene in enumerate(genes)}
    
    matrix = np.zeros((n_genes, n_genes))
    
    for _, row in prediction.iterrows():
        source_idx = gene_to_idx.get(row['source'])
        target_idx = gene_to_idx.get(row['target'])
        
        if source_idx is not None and target_idx is not None:
            matrix[source_idx, target_idx] = row['weight']
    
    return matrix


def convert_matrix_to_edgelist(matrix: np.ndarray, genes: list) -> pd.DataFrame:
    """
    Convert adjacency matrix to edge list.
    
    Parameters:
    -----------
    matrix : np.ndarray
        Adjacency matrix (genes x genes)
    genes : list
        List of all genes
        
    Returns:
    --------
    prediction : pd.DataFrame
        GRN with columns: source, target, weight
    """
    rows, cols = np.where(matrix != 0)
    
    edges = []
    for i, j in zip(rows, cols):
        edges.append({
            'source': genes[i],
            'target': genes[j],
            'weight': matrix[i, j]
        })
    
    return pd.DataFrame(edges)


def create_degree_preserved_permutation(prediction: pd.DataFrame, intensity: float) -> pd.DataFrame:
    """
    Create degree-preserved permutation at specified intensity.
    
    Parameters:
    -----------
    prediction : pd.DataFrame
        Original GRN with columns: source, target, weight
    intensity : float
        Mixing intensity (0.0 = original, 1.0 = fully permuted)
        
    Returns:
    --------
    prediction_mixed : pd.DataFrame
        Mixed GRN
    """
    if intensity == 0.0:
        return prediction.copy()
    
    # Get all genes
    all_genes = sorted(set(prediction['source'].unique()) | set(prediction['target'].unique()))
    
    # Convert to matrix
    A_original = convert_edgelist_to_matrix(prediction, all_genes)
    
    # Create degree-preserved baseline
    A_baseline = create_grn_baseline(A_original)
    
    if intensity == 1.0:
        # Return fully permuted network
        return convert_matrix_to_edgelist(A_baseline, all_genes)
    
    # Mix original and baseline
    # For partial mixing, we need to decide which edges come from which network
    # Strategy: randomly select edges from each network based on intensity
    
    # Get nonzero positions from both networks
    original_mask = (A_original != 0)
    baseline_mask = (A_baseline != 0)
    
    # Create mixed network
    A_mixed = np.zeros_like(A_original)
    
    # For each position that has an edge in either network
    all_edges_mask = original_mask | baseline_mask
    rows, cols = np.where(all_edges_mask)
    
    np.random.seed(42)  # For reproducibility
    
    for i, j in zip(rows, cols):
        # Decide whether to use baseline or original edge
        use_baseline = np.random.random() < intensity
        
        if use_baseline and baseline_mask[i, j]:
            A_mixed[i, j] = A_baseline[i, j]
        elif not use_baseline and original_mask[i, j]:
            A_mixed[i, j] = A_original[i, j]
        # If the edge doesn't exist in the chosen network, we skip it
    
    # Convert back to edge list
    prediction_mixed = convert_matrix_to_edgelist(A_mixed, all_genes)
    
    return prediction_mixed


def evaluate_degree_preserved_grn(par, intensity, seed=42):
    """
    Evaluate vc metric on degree-preserved permuted GRN.
    
    Parameters:
    -----------
    par : dict
        Parameters for evaluation
    intensity : float
        Permutation intensity (0.0 to 1.0)
    seed : int
        Random seed for reproducibility
        
    Returns:
    --------
    results : dict
        Metric scores
    """
    np.random.seed(seed)
    
    # Load original prediction
    net = ad.read_h5ad(par['prediction'])
    prediction = pd.DataFrame(net.uns['prediction'])
    prediction = process_links(prediction, par={'max_n_links': 50_000})
    
    # Create degree-preserved permutation
    prediction_permuted = create_degree_preserved_permutation(prediction, intensity)
    
    # Create temporary file with permuted GRN
    temp_file = f"{par['output_dir']}/tmp_degree_preserved_{int(intensity*100)}.h5ad"
    os.makedirs(os.path.dirname(temp_file), exist_ok=True)
    
    net_permuted = ad.AnnData(
        X=None,
        uns={
            "method_id": net.uns['method_id'],
            "dataset_id": net.uns['dataset_id'],
            "prediction": prediction_permuted[["source", "target", "weight"]]
        }
    )
    net_permuted.write_h5ad(temp_file)
    
    # Evaluate metric
    par_eval = par.copy()
    par_eval['prediction'] = temp_file
    
    try:
        scores = main_vc(par_eval)
        results = {
            'permutation_type': 'degree_preserved',
            'intensity': intensity * 100,  # Convert to percentage
            'n_edges_original': len(prediction),
            'n_edges_permuted': len(prediction_permuted)
        }
        # Add all metric columns dynamically
        for col in scores.columns:
            results[col] = scores[col].values[0]
    except Exception as e:
        print(f"Error evaluating degree_preserved at intensity {intensity}: {e}")
        import traceback
        traceback.print_exc()
        results = {
            'permutation_type': 'degree_preserved',
            'intensity': intensity * 100,
            'vc': np.nan,
            'n_edges_original': len(prediction),
            'n_edges_permuted': 0
        }
    finally:
        # Clean up temp file
        if os.path.exists(temp_file):
            os.remove(temp_file)
    
    return results


def analyze_degree_preserved_robustness(par):
    """
    Main function to analyze metric robustness with degree-preserved permutations.
    
    Parameters:
    -----------
    par : dict
        Configuration parameters with keys:
        - prediction: path to GRN prediction
        - evaluation_data: path to evaluation data
        - output_dir: directory to save results
        - intensities: list of mixing intensities (0.0 to 1.0)
    """
    # Default parameters
    if 'intensities' not in par:
        par['intensities'] = [0.0, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0]
    
    # Create output directory
    os.makedirs(par['output_dir'], exist_ok=True)
    
    print("="*60)
    print("VC Robustness Analysis - Degree-Preserved Permutations")
    print("="*60)
    print(f"GRN: {par['prediction']}")
    print(f"Evaluation data: {par['evaluation_data']}")
    print(f"Mixing intensities: {[i*100 for i in par['intensities']]}%")
    print("="*60)
    
    # Evaluate for each intensity
    all_results = []
    for intensity in par['intensities']:
        print(f"\nEvaluating intensity: {intensity*100:.0f}%")
        results = evaluate_degree_preserved_grn(par, intensity)
        all_results.append(results)
        # Print scores
        if 'vc' in results:
            print(f"  VC R²: {results['vc']:.4f}")
            print(f"  Edges: {results['n_edges_original']} -> {results['n_edges_permuted']}")
    
    # Convert to DataFrame
    df_results = pd.DataFrame(all_results)
    
    # Save results
    csv_file = f"{par['output_dir']}/degree_preserved_robustness.csv"
    df_results.to_csv(csv_file, index=False)
    print(f"\nResults saved to: {csv_file}")
    
    # Create visualization
    create_robustness_plot(df_results, par['output_dir'])
    
    # Print summary
    print("\n" + "="*60)
    print("Summary:")
    print("="*60)
    print(df_results.to_string(index=False))
    
    # Calculate robustness metrics
    if len(df_results) >= 2 and 'vc' in df_results.columns:
        initial = df_results[df_results['intensity'] == 0.0]['vc'].values[0]
        final = df_results[df_results['intensity'] == 100.0]['vc'].values[0]
        change = final - initial
        pct_change = (change / initial * 100) if initial != 0 else np.nan
        
        print(f"\nRobustness Summary:")
        print(f"="*60)
        print(f"Initial R² (0% permutation): {initial:.4f}")
        print(f"Final R² (100% permutation): {final:.4f}")
        print(f"Change: {change:+.4f} ({pct_change:+.1f}%)")
        
        # Calculate area under curve (normalized)
        auc = np.trapz(df_results['vc'].values, df_results['intensity'].values / 100)
        print(f"Area Under Curve: {auc:.4f}")
    
    return df_results


def create_robustness_plot(df_results, output_dir):
    """
    Create visualization of robustness analysis.
    
    Parameters:
    -----------
    df_results : pd.DataFrame
        Results from robustness analysis
    output_dir : str
        Directory to save plot
    """
    fig, ax = plt.subplots(1, 1, figsize=(10, 6))
    
    if 'vc' in df_results.columns:
        ax.plot(df_results['intensity'], df_results['vc'], 
               marker='o', linewidth=2.5, markersize=10, 
               color='#2E86AB', label='Degree-Preserved Permutation')
        
        # Add shaded region to show robustness
        ax.fill_between(df_results['intensity'], 
                        df_results['vc'], 
                        alpha=0.2, color='#2E86AB')
    
    ax.set_xlabel('Permutation Intensity (%)', fontsize=12)
    ax.set_ylabel('R² Score', fontsize=12)
    ax.set_title('VC Metric Robustness to Degree-Preserved Network Permutations', 
                fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.legend(loc='best', fontsize=10)
    ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
    
    # Add annotation for key points
    if len(df_results) >= 2 and 'vc' in df_results.columns:
        initial_val = df_results[df_results['intensity'] == 0.0]['vc'].values[0]
        final_val = df_results[df_results['intensity'] == 100.0]['vc'].values[0]
        
        ax.annotate(f'Original: {initial_val:.3f}', 
                   xy=(0, initial_val), xytext=(10, 10),
                   textcoords='offset points', fontsize=9,
                   bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.3))
        
        ax.annotate(f'Fully Permuted: {final_val:.3f}', 
                   xy=(100, final_val), xytext=(-10, 10),
                   textcoords='offset points', fontsize=9,
                   bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.3),
                   ha='right')
    
    plt.tight_layout()
    
    # Save plot
    plot_file = f"{output_dir}/degree_preserved_robustness_plot.png"
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    print(f"Plot saved to: {plot_file}")
    plt.close()


if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(description='Analyze vc metric robustness with degree-preserved permutations')
    parser.add_argument('--prediction', type=str, required=True,
                       help='Path to GRN prediction h5ad file')
    parser.add_argument('--evaluation_data', type=str, required=True,
                       help='Path to evaluation data h5ad file')
    parser.add_argument('--output_dir', type=str, required=True,
                       help='Output directory for results')
    parser.add_argument('--n_top_genes', type=int, default=3000)
    parser.add_argument('--intensities', type=float, nargs='+',
                       default=[0.0, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0],
                       help='Permutation intensities (0.0 to 1.0)')
    
    args = parser.parse_args()
    
    # Create parameter dictionary
    par = {
        'prediction': args.prediction,
        'evaluation_data': args.evaluation_data,
        'output_dir': args.output_dir,
        'intensities': args.intensities,
        'score': f"{args.output_dir}/vc_score.h5ad",
        'n_top_genes': args.n_top_genes
    }
    
    # Run analysis
    analyze_degree_preserved_robustness(par)
