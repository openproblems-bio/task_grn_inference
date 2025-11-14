"""
Analyze robustness of vc metric by permuting GRN at different intensities.

This script:
1. Loads a GRN prediction
2. Permutes edge directions/signs/weights at different intensities (0%, 20%, 50%, 100%)
3. Evaluates vc metric for each permutation
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
sys.path.append('src/metrics/vc')
sys.path.append('src/utils')

from helper import main as main_vc
from util import process_links


def permute_grn(prediction: pd.DataFrame, degree: float, noise_type: str) -> pd.DataFrame:
    """
    Permute GRN with different noise types.
    
    Parameters:
    -----------
    prediction : pd.DataFrame
        GRN with columns: source, target, weight
    degree : float
        Fraction of perturbation (0.0 to 1.0)
    noise_type : str
        Type of perturbation: 'direction', 'sign', 'weight', 'net'
        
    Returns:
    --------
    prediction_permuted : pd.DataFrame
        Perturbed GRN
    """
    prediction = prediction.reset_index(drop=True).copy()
    
    if noise_type == 'direction':
        # Swap source and target for selected edges
        n_rows_to_permute = int(len(prediction) * degree)
        if n_rows_to_permute > 0:
            indices_to_permute = np.random.choice(prediction.index, size=n_rows_to_permute, replace=False)
            prediction.loc[indices_to_permute, ['source', 'target']] = \
                prediction.loc[indices_to_permute, ['target', 'source']].values
    
    elif noise_type == 'sign':
        # Flip sign of weights for selected edges
        num_rows = len(prediction)
        num_to_modify = int(num_rows * degree)
        if num_to_modify > 0:
            random_indices = np.random.choice(prediction.index, size=num_to_modify, replace=False)
            prediction.loc[random_indices, 'weight'] *= -1
    
    elif noise_type == 'weight':
        # Add noise to weights
        if 'weight' in prediction.columns and degree > 0:
            std_dev = prediction['weight'].std()
            noise = np.random.normal(loc=0, scale=degree * std_dev, size=len(prediction))
            prediction['weight'] += noise
    
    elif noise_type == 'net':
        # Shuffle network structure
        if degree > 0:
            # Group by target and source
            prediction = prediction.groupby(['target', 'source'], as_index=False)['weight'].mean()
            
            # Create pivot table
            pivot_df = prediction.pivot(index='target', columns='source', values='weight')
            pivot_df.fillna(0, inplace=True)
            
            # Flatten and shuffle
            matrix_flattened = pivot_df.values.flatten()
            n_elements = len(matrix_flattened)
            n_shuffle = int(n_elements * degree)
            
            if n_shuffle > 0:
                shuffle_indices = np.random.choice(n_elements, n_shuffle, replace=False)
                shuffle_values = matrix_flattened[shuffle_indices].copy()
                np.random.shuffle(shuffle_values)
                matrix_flattened[shuffle_indices] = shuffle_values
            
            # Reshape back
            pivot_df_shuffled = pd.DataFrame(
                matrix_flattened.reshape(pivot_df.shape),
                index=pivot_df.index,
                columns=pivot_df.columns
            )
            
            # Convert back to edge list
            flat_df = pivot_df_shuffled.reset_index()
            prediction = flat_df.melt(id_vars='target', var_name='source', value_name='weight')
            prediction = prediction[prediction['weight'] != 0].reset_index(drop=True)
    
    else:
        raise ValueError(f"Unknown noise_type: {noise_type}")
    
    return prediction


def evaluate_permuted_grn(par, degree, noise_type, seed=42):
    """
    Evaluate vc metric on permuted GRN.
    
    Parameters:
    -----------
    par : dict
        Parameters for evaluation
    degree : float
        Permutation intensity (0.0 to 1.0)
    noise_type : str
        Type of perturbation
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
    
    # Permute network
    prediction_permuted = permute_grn(prediction, degree, noise_type)
    
    # Create temporary file with permuted GRN
    temp_file = f"{par['output_dir']}/tmp_permuted_{noise_type}_{int(degree*100)}.h5ad"
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
            'noise_type': noise_type,
            'degree': degree * 100,  # Convert to percentage
            'n_edges': len(prediction),
            'n_edges_permuted': int(len(prediction) * degree)
        }
        # Add all metric columns dynamically
        for col in scores.columns:
            results[col] = scores[col].values[0]
    except Exception as e:
        print(f"Error evaluating {noise_type} at degree {degree}: {e}")
        import traceback
        traceback.print_exc()
        results = {
            'noise_type': noise_type,
            'degree': degree * 100,
            'vc': np.nan,
            'n_edges': len(prediction),
            'n_edges_permuted': int(len(prediction) * degree)
        }
    finally:
        # Clean up temp file
        if os.path.exists(temp_file):
            os.remove(temp_file)
    
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
        - noise_types: list of noise types ('direction', 'sign', 'weight', 'net')
    """
    # Default parameters
    if 'degrees' not in par:
        par['degrees'] = [0.0, 0.2, 0.5, 1.0]
    if 'noise_types' not in par:
        par['noise_types'] = ['direction', 'sign', 'weight', 'net']
    
    # Create output directory
    os.makedirs(par['output_dir'], exist_ok=True)
    
    print("="*60)
    print("VC Robustness Analysis")
    print("="*60)
    print(f"GRN: {par['prediction']}")
    print(f"Evaluation data: {par['evaluation_data']}")
    print(f"Permutation types: {par['noise_types']}")
    print(f"Permutation degrees: {[d*100 for d in par['degrees']]}%")
    print("="*60)
    
    # Evaluate for each noise type and degree
    all_results = []
    for noise_type in par['noise_types']:
        print(f"\n{'='*60}")
        print(f"Noise Type: {noise_type}")
        print(f"{'='*60}")
        
        for degree in par['degrees']:
            print(f"\nEvaluating degree: {degree*100:.0f}%")
            results = evaluate_permuted_grn(par, degree, noise_type)
            all_results.append(results)
            # Print scores
            if 'vc' in results:
                print(f"  VC R²: {results['vc']:.4f}")
    
    # Convert to DataFrame
    df_results = pd.DataFrame(all_results)
    
    # Save results
    csv_file = f"{par['output_dir']}/robustness_analysis.csv"
    df_results.to_csv(csv_file, index=False)
    print(f"\nResults saved to: {csv_file}")
    
    # Create visualization
    create_robustness_plot(df_results, par['output_dir'])
    
    # Print summary
    print("\n" + "="*60)
    print("Summary:")
    print("="*60)
    print(df_results.to_string(index=False))
    
    # Calculate robustness metrics for each noise type
    print(f"\nRobustness Summary by Noise Type:")
    print(f"="*60)
    for noise_type in par['noise_types']:
        df_noise = df_results[df_results['noise_type'] == noise_type]
        if len(df_noise) >= 2 and 'vc' in df_noise.columns:
            initial = df_noise[df_noise['degree'] == 0.0]['vc'].values[0]
            final = df_noise[df_noise['degree'] == 100.0]['vc'].values[0]
            change = final - initial
            pct_change = (change / initial * 100) if initial != 0 else np.nan
            print(f"\n{noise_type}:")
            print(f"  Initial R² (0%): {initial:.4f}")
            print(f"  Final R² (100%): {final:.4f}")
            print(f"  Change: {change:+.4f} ({pct_change:+.1f}%)")
    
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
    noise_types = df_results['noise_type'].unique()
    colors = {'direction': '#2E86AB', 'sign': '#A23B72', 'weight': '#F18F01', 'net': '#6A994E'}
    
    # Create single plot for VC metric
    fig, ax = plt.subplots(1, 1, figsize=(10, 6))
    
    if 'vc' in df_results.columns:
        for noise_type in noise_types:
            df_noise = df_results[df_results['noise_type'] == noise_type]
            ax.plot(df_noise['degree'], df_noise['vc'], 
                   marker='o', linewidth=2, markersize=8, 
                   color=colors.get(noise_type, 'gray'), label=noise_type)
    
    ax.set_xlabel('Permutation Degree (%)', fontsize=12)
    ax.set_ylabel('R² Score', fontsize=12)
    ax.set_title('VC Metric Robustness to GRN Perturbations', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.legend(loc='best', fontsize=10)
    ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
    
    plt.tight_layout()
    
    # Save plot
    plot_file = f"{output_dir}/robustness_plot.png"
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    print(f"Plot saved to: {plot_file}")
    plt.close()


if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(description='Analyze vc metric robustness')
    parser.add_argument('--prediction', type=str, required=True,
                       help='Path to GRN prediction h5ad file')
    parser.add_argument('--evaluation_data', type=str, required=True,
                       help='Path to evaluation data h5ad file')
    parser.add_argument('--output_dir', type=str, required=True,
                       help='Output directory for results')
    parser.add_argument('--degrees', type=float, nargs='+',
                       default=[0.0, 0.2, 0.5, 1.0],
                       help='Permutation degrees (0.0 to 1.0)')
    parser.add_argument('--noise_types', type=str, nargs='+',
                       default=['direction', 'sign', 'weight', 'net'],
                       help='Permutation types: direction, sign, weight, net')
    
    args = parser.parse_args()
    
    # Create parameter dictionary
    par = {
        'prediction': args.prediction,
        'evaluation_data': args.evaluation_data,
        'output_dir': args.output_dir,
        'degrees': args.degrees,
        'noise_types': args.noise_types,
        'score': f"{args.output_dir}/vc_score.h5ad"
    }
    
    # Run analysis
    analyze_robustness(par)
