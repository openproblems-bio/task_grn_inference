"""
Robustness Analysis for Pathway Annotation Metric

This script analyzes the robustness of the pathway annotation metric by:
1. Testing sensitivity to parameter choices (FDR, target ranges)
2. Comparing against randomized baseline GRNs
3. Stratifying results by TF properties
4. Analyzing pathway database choice effects
"""

import os
import sys
import pandas as pd
import numpy as np
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Dict, List

# Add paths
sys.path.append('src/metrics/experimental/annotation')
sys.path.append('src/utils')

from helper import main as main_pathway_annotation
from util import process_links, create_grn_baseline


def analyze_parameter_sensitivity(par: dict) -> pd.DataFrame:
    """
    Test sensitivity to FDR threshold and target count ranges.
    
    Returns
    -------
    pd.DataFrame
        Results for different parameter combinations
    """
    print("\n" + "="*80)
    print("PARAMETER SENSITIVITY ANALYSIS")
    print("="*80)
    
    # Parameter ranges to test
    fdr_thresholds = [0.01, 0.05, 0.1, 0.2]
    target_ranges = [(5, 200), (10, 500), (20, 1000)]
    
    results = []
    
    for fdr in fdr_thresholds:
        for min_t, max_t in target_ranges:
            print(f"\nTesting: FDR={fdr}, Targets={min_t}-{max_t}")
            
            par_copy = par.copy()
            par_copy['fdr_threshold'] = fdr
            par_copy['min_targets'] = min_t
            par_copy['max_targets'] = max_t
            
            try:
                result = main_pathway_annotation(par_copy)
                result['fdr_threshold'] = fdr
                result['min_targets'] = min_t
                result['max_targets'] = max_t
                results.append(result)
            except Exception as e:
                print(f"  ERROR: {e}")
                continue
    
    df = pd.concat(results, ignore_index=True) if results else pd.DataFrame()
    return df


def compare_to_baseline(par: dict, n_randomizations: int = 5) -> pd.DataFrame:
    """
    Compare real GRN scores to randomized baseline networks.
    
    Parameters
    ----------
    par : dict
        Parameters
    n_randomizations : int
        Number of randomized networks to generate
    
    Returns
    -------
    pd.DataFrame
        Comparison results showing real vs random performance
    """
    print("\n" + "="*80)
    print(f"BASELINE COMPARISON ({n_randomizations} randomizations)")
    print("="*80)
    
    # Get real score
    print("\n[1] Evaluating real GRN...")
    real_result = main_pathway_annotation(par)
    real_result['network_type'] = 'real'
    
    # Load prediction as matrix for randomization
    print("\n[2] Loading GRN for randomization...")
    net = ad.read_h5ad(par['prediction'])
    prediction = pd.DataFrame(net.uns['prediction'])
    prediction = process_links(prediction, par)
    
    evaluation_data = ad.read_h5ad(par['evaluation_data'], backed='r')
    genes = evaluation_data.var_names.tolist()
    
    # Create adjacency matrix
    gene_to_idx = {g: i for i, g in enumerate(genes)}
    n_genes = len(genes)
    A = np.zeros((n_genes, n_genes))
    
    for _, row in prediction.iterrows():
        if row['source'] in gene_to_idx and row['target'] in gene_to_idx:
            i = gene_to_idx[row['source']]
            j = gene_to_idx[row['target']]
            A[i, j] = row['weight']
    
    # Generate randomized networks
    print("\n[3] Generating randomized networks...")
    random_results = []
    
    for i in range(n_randomizations):
        print(f"  Randomization {i+1}/{n_randomizations}...")
        
        # Create baseline (degree-preserving randomization)
        A_random = create_grn_baseline(A)
        
        # Convert back to edge list
        sources, targets = np.where(A_random != 0)
        random_edges = pd.DataFrame({
            'source': [genes[s] for s in sources],
            'target': [genes[t] for t in targets],
            'weight': A_random[sources, targets]
        })
        
        # Save temporarily and score
        temp_pred_file = f'/tmp/random_pred_pathway_{i}.h5ad'
        net_random = ad.AnnData(
            X=None,
            uns={
                "method_id": net.uns['method_id'],
                "dataset_id": net.uns['dataset_id'],
                "prediction": random_edges
            }
        )
        net_random.write_h5ad(temp_pred_file)
        
        par_copy = par.copy()
        par_copy['prediction'] = temp_pred_file
        
        try:
            random_result = main_pathway_annotation(par_copy)
            random_result['network_type'] = f'random_{i}'
            random_results.append(random_result)
        except Exception as e:
            print(f"    ERROR: {e}")
            continue
        finally:
            if os.path.exists(temp_pred_file):
                os.remove(temp_pred_file)
    
    # Combine results
    all_results = pd.concat([real_result] + random_results, ignore_index=True)
    
    # Print summary
    print("\n" + "="*80)
    print("BASELINE COMPARISON RESULTS:")
    print("="*80)
    real_score = all_results[all_results['network_type'] == 'real']['pf_grn'].values[0]
    random_scores = all_results[all_results['network_type'].str.startswith('random')]['pf_grn']
    
    print(f"Real GRN score:        {real_score:.4f}")
    print(f"Random mean:           {random_scores.mean():.4f}")
    print(f"Random std:            {random_scores.std():.4f}")
    print(f"Improvement vs random: {real_score - random_scores.mean():.4f}")
    print("="*80)
    
    return all_results


def plot_results(sensitivity_df: pd.DataFrame, 
                baseline_df: pd.DataFrame,
                output_dir: str):
    """
    Create visualization plots for robustness analysis.
    """
    os.makedirs(output_dir, exist_ok=True)
    
    # 1. Parameter sensitivity heatmap
    if len(sensitivity_df) > 0:
        print(f"\nCreating sensitivity plots...")
        plt.figure(figsize=(15, 5))
        
        for i, metric in enumerate(['pf_grn', 'pf_all', 'coverage_pct']):
            plt.subplot(1, 3, i+1)
            pivot = sensitivity_df.pivot_table(
                values=metric,
                index='fdr_threshold',
                columns='max_targets',
                aggfunc='mean'
            )
            sns.heatmap(pivot, annot=True, fmt='.3f', cmap='viridis', cbar_kws={'label': metric})
            plt.title(f'{metric}')
            plt.xlabel('Max Targets')
            plt.ylabel('FDR Threshold')
        
        plt.tight_layout()
        plt.savefig(f'{output_dir}/parameter_sensitivity.png', dpi=300, bbox_inches='tight')
        plt.close()
        print(f"  Saved: {output_dir}/parameter_sensitivity.png")
    
    # 2. Baseline comparison
    if len(baseline_df) > 0:
        print(f"Creating baseline comparison plots...")
        metrics = ['pf_grn', 'pf_all', 'coverage_pct']
        
        fig, axes = plt.subplots(1, 3, figsize=(15, 5))
        
        for i, metric in enumerate(metrics):
            ax = axes[i]
            
            # Separate real vs random
            real_val = baseline_df[baseline_df['network_type'] == 'real'][metric].values
            random_vals = baseline_df[baseline_df['network_type'].str.startswith('random')][metric].values
            
            if len(random_vals) > 0:
                # Box plot for random
                bp = ax.boxplot([random_vals], positions=[0], widths=0.5, 
                          patch_artist=True, boxprops=dict(facecolor='lightblue', alpha=0.7))
            
            # Real value as red dot
            if len(real_val) > 0:
                ax.scatter([0], real_val, color='red', s=150, zorder=5, 
                          label='Real GRN', marker='*', edgecolors='darkred', linewidths=1.5)
            
            ax.set_ylabel(metric.replace('_', ' ').title(), fontsize=11)
            ax.set_xticks([0])
            ax.set_xticklabels(['Randomized GRNs'])
            ax.legend(loc='best')
            ax.grid(axis='y', alpha=0.3)
            
        plt.suptitle('Real GRN vs Degree-Preserving Randomized Baselines', fontsize=13, fontweight='bold')
        plt.tight_layout()
        plt.savefig(f'{output_dir}/baseline_comparison.png', dpi=300, bbox_inches='tight')
        plt.close()
        print(f"  Saved: {output_dir}/baseline_comparison.png")


def analyze_robustness(par: dict):
    """
    Main function to run comprehensive robustness analysis.
    
    Parameters
    ----------
    par : dict
        Configuration parameters with keys:
        - prediction: path to GRN prediction
        - evaluation_data: path to evaluation data
        - output_dir: directory to save results
        - tf_all: path to TF list
        - pathway_file: path to pathway GMT file
    """
    # Create output directory
    output_dir = par.get('output_dir', 'output/pathway_annotation_robustness')
    os.makedirs(output_dir, exist_ok=True)
    
    print("\n" + "="*80)
    print("PATHWAY ANNOTATION METRIC - ROBUSTNESS ANALYSIS")
    print("="*80)
    print(f"GRN:             {par['prediction']}")
    print(f"Evaluation data: {par['evaluation_data']}")
    print(f"Output dir:      {output_dir}")
    print("="*80)
    
    # 1. Parameter sensitivity analysis
    print("\n[ANALYSIS 1/2] Parameter Sensitivity")
    sensitivity_results = analyze_parameter_sensitivity(par)
    if len(sensitivity_results) > 0:
        csv_file = f"{output_dir}/parameter_sensitivity.csv"
        sensitivity_results.to_csv(csv_file, index=False)
        print(f"\nSaved: {csv_file}")
    
    # 2. Baseline comparison
    print("\n[ANALYSIS 2/2] Baseline Comparison")
    baseline_results = compare_to_baseline(par, n_randomizations=5)
    if len(baseline_results) > 0:
        csv_file = f"{output_dir}/baseline_comparison.csv"
        baseline_results.to_csv(csv_file, index=False)
        print(f"\nSaved: {csv_file}")
    
    # 3. Create plots
    print("\n[VISUALIZATION] Creating plots...")
    plot_results(sensitivity_results, baseline_results, output_dir)
    
    print("\n" + "="*80)
    print("ROBUSTNESS ANALYSIS COMPLETE!")
    print(f"Results saved to: {output_dir}")
    print("="*80 + "\n")
    
    return {
        'sensitivity': sensitivity_results,
        'baseline': baseline_results
    }



if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(description='Analyze pathway annotation metric robustness')
    parser.add_argument('--prediction', type=str, required=True,
                       help='Path to GRN prediction h5ad file')
    parser.add_argument('--evaluation_data', type=str, required=True,
                       help='Path to evaluation data h5ad file')
    parser.add_argument('--output_dir', type=str, default='output/pathway_annotation_robustness',
                       help='Output directory for results')
    parser.add_argument('--tf_all', type=str, 
                       default='resources/grn_benchmark/prior/tf_all.csv',
                       help='Path to TF list CSV file')
    parser.add_argument('--pathway_file', type=str,
                       default='/home/jnourisa/projs/ongoing/ciim/input/prior/h.all.v2024.1.Hs.symbols.gmt',
                       help='Path to pathway GMT file')
    parser.add_argument('--max_n_links', type=int, default=50000,
                       help='Maximum number of GRN links to use')
    
    args = parser.parse_args()
    
    # Create parameter dictionary
    par = {
        'prediction': args.prediction,
        'evaluation_data': args.evaluation_data,
        'output_dir': args.output_dir,
        'tf_all': args.tf_all,
        'pathway_file': args.pathway_file,
        'max_n_links': args.max_n_links
    }
    
    # Run analysis
    results = analyze_robustness(par)

