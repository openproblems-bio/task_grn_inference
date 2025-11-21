"""
Analyze robustness of tf_binding metric by permuting GRN directions at different intensities.

This script:
1. Loads a GRN prediction
2. Permutes edge directions at different intensities (0%, 20%, 50%, 100%)
3. Evaluates tf_binding metric for each permutation
4. Analyzes how the metric changes with increasing permutation
"""

import os
import sys
import pandas as pd
import numpy as np
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns

# Add paths
sys.path.append('src/metrics/tf_binding')
sys.path.append('src/utils')

from helper import main as main_tf_binding
from util import process_links


def permute_grn_direction(prediction: pd.DataFrame, degree: float) -> pd.DataFrame:
    """
    Permute edge directions (swap source and target) with given intensity.
    
    Parameters:
    -----------
    prediction : pd.DataFrame
        GRN with columns: source, target, weight
    degree : float
        Fraction of edges to permute (0.0 to 1.0)
        
    Returns:
    --------
    prediction_permuted : pd.DataFrame
        GRN with permuted directions
    """
    prediction = prediction.reset_index(drop=True).copy()
    
    # Number of edges to permute
    n_rows_to_permute = int(len(prediction) * degree)
    
    if n_rows_to_permute == 0:
        return prediction
    
    # Random indices to permute
    indices_to_permute = np.random.choice(prediction.index, size=n_rows_to_permute, replace=False)
    
    # Swap source and target for selected edges
    prediction.loc[indices_to_permute, ['source', 'target']] = \
        prediction.loc[indices_to_permute, ['target', 'source']].values
    
    return prediction


def evaluate_permuted_grn(par, degree, seed=42):
    """
    Evaluate tf_binding metric on permuted GRN.
    
    Parameters:
    -----------
    par : dict
        Parameters for evaluation
    degree : float
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
    
    # Permute directions
    prediction_permuted = permute_grn_direction(prediction, degree)
    
    # Create temporary file with permuted GRN
    temp_file = f"{par['output_dir']}/tmp_permuted_{int(degree*100)}.h5ad"
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
    print(par_eval)
    par_eval['prediction'] = temp_file
    
    
    scores = main_tf_binding(par_eval)
    results = {
        'degree': degree * 100,  # Convert to percentage
        'tfb_grn': scores['tfb_grn'].values[0],
        'tfb_all': scores['tfb_all'].values[0],
        'n_edges_permuted': int(len(prediction) * degree)
    }
    
    
    return results


def analyze_robustness(par):
    """
    Main function to analyze metric robustness.
    
    Parameters:
    -----------
    par : dict
        Configuration parameters with keys:
        - prediction: path to GRN prediction
        - evaluation_data: path to evaluation data
        - output_dir: directory to save results
        - degrees: list of permutation intensities (0.0 to 1.0)
    """
    # Default degrees if not specified
    if 'degrees' not in par:
        par['degrees'] = [0.0, 0.2, 0.5, 1.0]
    
    # Create output directory
    os.makedirs(par['output_dir'], exist_ok=True)
    
    print("="*60)
    print("tf_binding Robustness Analysis")
    print("="*60)
    print(f"GRN: {par['prediction']}")
    print(f"Evaluation data: {par['evaluation_data']}")
    print(f"Permutation degrees: {[d*100 for d in par['degrees']]}%")
    print("="*60)
    
    # Evaluate for each degree
    all_results = []
    for degree in par['degrees']:
        print(f"\nEvaluating permutation degree: {degree*100:.0f}%")
        results = evaluate_permuted_grn(par, degree)
        all_results.append(results)
        print(f"  Precision: {results['tfb_grn']:.4f}")
        print(f"  Balanced: {results['tfb_all']:.2f}")
    
    # Convert to DataFrame
    df_results = pd.DataFrame(all_results)
    
    # Save results
    csv_file = f"{par['output_dir']}/robustness_analysis.csv"
    df_results.to_csv(csv_file, index=False)
    print(f"\nResults saved to: {csv_file}")
    
    # Create visualization
    
    # Print summary
    print("\n" + "="*60)
    print("Summary:")
    print("="*60)
    print(df_results.to_string(index=False))
    
    # Calculate robustness metrics
    if len(df_results) >= 2:
        precision_drop = df_results['tfb_grn'].iloc[0] - df_results['tfb_grn'].iloc[-1]
        balanced_drop = df_results['tfb_all'].iloc[0] - df_results['tfb_all'].iloc[-1]
        
        print(f"\nRobustness metrics:")
        print(f"  Precision drop (0% → 100%): {precision_drop:.4f}")
        print(f"  Balanced drop (0% → 100%): {balanced_drop:.2f}")
        print(f"  Relative precision change: {precision_drop/df_results['tfb_grn'].iloc[0]*100:.1f}%")
    
    return df_results



if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(description='Analyze tf_binding metric robustness')
    parser.add_argument('--prediction', type=str, required=True,
                       help='Path to GRN prediction h5ad file')
    parser.add_argument('--evaluation_data', type=str, required=True,
                       help='Path to evaluation data h5ad file')
    parser.add_argument('--output_dir', type=str, required=True,
                       help='Output directory for results')
    parser.add_argument('--degrees', type=float, nargs='+',
                       default=[0.0, 0.2, 0.5, 1.0],
                       help='Permutation degrees (0.0 to 1.0)')
    parser.add_argument('--layer', type=str, default='lognorm',
                       help='Expression layer to use')
    parser.add_argument('--min_targets', type=int, default=5,
                       help='Minimum targets per TF')
    parser.add_argument('--cell_type', type=str, default='PBMC',
                       help='Cell type to analyze')
    
    args = parser.parse_args()
    
    # Create parameter dictionary
    par = {
        'prediction': args.prediction,
        'evaluation_data': args.evaluation_data,
        'output_dir': args.output_dir,
        'degrees': args.degrees,
        'layer': args.layer,
        'min_targets': args.min_targets,
        'ground_truth_remap': f'resources/grn_benchmark/ground_truth/{args.cell_type}_remap.csv',
        'ground_truth_chipatlas': f'resources/grn_benchmark/ground_truth/{args.cell_type}_chipatlas.csv',
        'ground_truth_unibind': f'resources/grn_benchmark/ground_truth/{args.cell_type}_unibind.csv',
        'tf_all': 'resources/grn_benchmark/prior/tf_all.csv'
        
    }
    
    # Run analysis
    analyze_robustness(par)
