"""
Compare different permutation strategies for robustness analysis.

This script runs multiple permutation types and creates comparative visualizations.
"""

import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Add paths
sys.path.append('src/metrics/vc')
sys.path.append('src/utils')

from analyze_robustness import analyze_robustness
from analyze_robustness_degree_preserved import analyze_degree_preserved_robustness


def run_all_analyses(par):
    """
    Run all robustness analyses and create comparative plots.
    
    Parameters:
    -----------
    par : dict
        Configuration parameters
    """
    print("="*70)
    print("COMPREHENSIVE VC ROBUSTNESS ANALYSIS")
    print("="*70)
    
    # 1. Run standard permutations
    print("\n" + "="*70)
    print("PHASE 1: Standard Permutations")
    print("="*70)
    
    par_standard = par.copy()
    par_standard['output_dir'] = f"{par['base_output_dir']}/standard"
    par_standard['degrees'] = par.get('intensities', [0.0, 0.2, 0.5, 1.0])
    par_standard['noise_types'] = ['net', 'sign', 'direction']
    
    df_standard = analyze_robustness(par_standard)
    
    # 2. Run degree-preserved permutations
    print("\n" + "="*70)
    print("PHASE 2: Degree-Preserved Permutations")
    print("="*70)
    
    par_degree_preserved = par.copy()
    par_degree_preserved['output_dir'] = f"{par['base_output_dir']}/degree_preserved"
    par_degree_preserved['intensities'] = par.get('intensities', [0.0, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0])
    
    df_degree_preserved = analyze_degree_preserved_robustness(par_degree_preserved)
    
    # 3. Create comparative visualizations
    print("\n" + "="*70)
    print("PHASE 3: Comparative Analysis")
    print("="*70)
    
    create_comparison_plots(df_standard, df_degree_preserved, par['base_output_dir'])
    
    # 4. Generate summary report
    create_summary_report(df_standard, df_degree_preserved, par['base_output_dir'])
    
    print("\n" + "="*70)
    print("ANALYSIS COMPLETE")
    print("="*70)
    print(f"All results saved to: {par['base_output_dir']}")


def create_comparison_plots(df_standard, df_degree_preserved, output_dir):
    """
    Create comparative visualization of all permutation types.
    
    Parameters:
    -----------
    df_standard : pd.DataFrame
        Results from standard permutations
    df_degree_preserved : pd.DataFrame
        Results from degree-preserved permutations
    output_dir : str
        Directory to save plots
    """
    # Prepare data for plotting
    plot_data = []
    
    # Add standard permutations
    for _, row in df_standard.iterrows():
        if 'vc' in row:
            plot_data.append({
                'intensity': row['degree'],
                'vc': row['vc'],
                'type': row['noise_type'],
                'category': 'Standard'
            })
    
    # Add degree-preserved permutations
    for _, row in df_degree_preserved.iterrows():
        if 'vc' in row:
            plot_data.append({
                'intensity': row['intensity'],
                'vc': row['vc'],
                'type': 'degree_preserved',
                'category': 'Degree-Preserved'
            })
    
    df_plot = pd.DataFrame(plot_data)
    
    # Create comprehensive comparison plot
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    
    # Plot 1: All permutation types together
    ax1 = axes[0]
    colors = {
        'net': '#6A994E',
        'sign': '#A23B72',
        'direction': '#2E86AB',
        'degree_preserved': '#F18F01'
    }
    
    for perm_type in df_plot['type'].unique():
        df_type = df_plot[df_plot['type'] == perm_type]
        ax1.plot(df_type['intensity'], df_type['vc'], 
                marker='o', linewidth=2.5, markersize=8,
                color=colors.get(perm_type, 'gray'),
                label=perm_type.replace('_', ' ').title())
    
    ax1.set_xlabel('Permutation Intensity (%)', fontsize=12)
    ax1.set_ylabel('VC R² Score', fontsize=12)
    ax1.set_title('Comparison: All Permutation Types', fontsize=14, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    ax1.legend(loc='best', fontsize=10)
    ax1.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
    
    # Plot 2: Degree-preserved vs. Random net permutation
    ax2 = axes[1]
    
    df_net = df_plot[df_plot['type'] == 'net']
    df_deg_pres = df_plot[df_plot['type'] == 'degree_preserved']
    
    if not df_net.empty:
        ax2.plot(df_net['intensity'], df_net['vc'],
                marker='s', linewidth=2.5, markersize=8,
                color='#6A994E', label='Random Network Shuffle',
                linestyle='--')
    
    if not df_deg_pres.empty:
        ax2.plot(df_deg_pres['intensity'], df_deg_pres['vc'],
                marker='o', linewidth=2.5, markersize=8,
                color='#F18F01', label='Degree-Preserved Shuffle')
    
    ax2.set_xlabel('Permutation Intensity (%)', fontsize=12)
    ax2.set_ylabel('VC R² Score', fontsize=12)
    ax2.set_title('Focus: Topology vs. Degree Preservation', fontsize=14, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    ax2.legend(loc='best', fontsize=10)
    ax2.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
    
    plt.tight_layout()
    
    # Save plot
    plot_file = f"{output_dir}/comparison_all_permutations.png"
    os.makedirs(output_dir, exist_ok=True)
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    print(f"Comparison plot saved to: {plot_file}")
    plt.close()
    
    # Create detailed comparison plot with error bars/variability if needed
    create_detailed_comparison(df_plot, output_dir)


def create_detailed_comparison(df_plot, output_dir):
    """
    Create detailed comparison with subplots for each category.
    """
    fig = plt.figure(figsize=(14, 10))
    gs = fig.add_gridspec(2, 2, hspace=0.3, wspace=0.3)
    
    colors = {
        'net': '#6A994E',
        'sign': '#A23B72', 
        'direction': '#2E86AB',
        'degree_preserved': '#F18F01'
    }
    
    # Plot each permutation type separately
    perm_types = df_plot['type'].unique()
    
    for idx, perm_type in enumerate(perm_types[:4]):  # Max 4 subplots
        row = idx // 2
        col = idx % 2
        ax = fig.add_subplot(gs[row, col])
        
        df_type = df_plot[df_plot['type'] == perm_type]
        ax.plot(df_type['intensity'], df_type['vc'],
               marker='o', linewidth=2, markersize=8,
               color=colors.get(perm_type, 'gray'))
        
        # Add shaded region
        ax.fill_between(df_type['intensity'], df_type['vc'],
                       alpha=0.2, color=colors.get(perm_type, 'gray'))
        
        ax.set_xlabel('Permutation Intensity (%)', fontsize=10)
        ax.set_ylabel('VC R² Score', fontsize=10)
        ax.set_title(f'{perm_type.replace("_", " ").title()}', fontsize=12, fontweight='bold')
        ax.grid(True, alpha=0.3)
        ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
        
        # Add stats
        if len(df_type) >= 2:
            initial = df_type[df_type['intensity'] == df_type['intensity'].min()]['vc'].values[0]
            final = df_type[df_type['intensity'] == df_type['intensity'].max()]['vc'].values[0]
            change = final - initial
            
            stats_text = f'Δ R² = {change:+.3f}'
            ax.text(0.05, 0.95, stats_text, transform=ax.transAxes,
                   fontsize=9, verticalalignment='top',
                   bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.suptitle('Detailed Robustness Analysis by Permutation Type', 
                fontsize=16, fontweight='bold', y=0.995)
    
    plot_file = f"{output_dir}/detailed_comparison.png"
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    print(f"Detailed comparison plot saved to: {plot_file}")
    plt.close()


def create_summary_report(df_standard, df_degree_preserved, output_dir):
    """
    Create text summary report of all analyses.
    """
    report_file = f"{output_dir}/summary_report.txt"
    
    with open(report_file, 'w') as f:
        f.write("="*70 + "\n")
        f.write("VC METRIC ROBUSTNESS ANALYSIS - COMPREHENSIVE SUMMARY\n")
        f.write("="*70 + "\n\n")
        
        # Standard permutations summary
        f.write("STANDARD PERMUTATIONS\n")
        f.write("-"*70 + "\n")
        
        for noise_type in df_standard['noise_type'].unique():
            df_noise = df_standard[df_standard['noise_type'] == noise_type]
            if len(df_noise) >= 2 and 'vc' in df_noise.columns:
                initial = df_noise[df_noise['degree'] == df_noise['degree'].min()]['vc'].values[0]
                final = df_noise[df_noise['degree'] == df_noise['degree'].max()]['vc'].values[0]
                change = final - initial
                pct_change = (change / initial * 100) if initial != 0 else np.nan
                
                f.write(f"\n{noise_type.upper()}:\n")
                f.write(f"  Initial R²: {initial:.4f}\n")
                f.write(f"  Final R²: {final:.4f}\n")
                f.write(f"  Absolute change: {change:+.4f}\n")
                f.write(f"  Percent change: {pct_change:+.1f}%\n")
        
        # Degree-preserved summary
        f.write("\n" + "="*70 + "\n")
        f.write("DEGREE-PRESERVED PERMUTATIONS\n")
        f.write("-"*70 + "\n")
        
        if len(df_degree_preserved) >= 2 and 'vc' in df_degree_preserved.columns:
            initial = df_degree_preserved[df_degree_preserved['intensity'] == 0.0]['vc'].values[0]
            final = df_degree_preserved[df_degree_preserved['intensity'] == 100.0]['vc'].values[0]
            change = final - initial
            pct_change = (change / initial * 100) if initial != 0 else np.nan
            
            f.write(f"\nDEGREE-PRESERVED:\n")
            f.write(f"  Initial R² (0%): {initial:.4f}\n")
            f.write(f"  Final R² (100%): {final:.4f}\n")
            f.write(f"  Absolute change: {change:+.4f}\n")
            f.write(f"  Percent change: {pct_change:+.1f}%\n")
            
            # Calculate AUC
            auc = np.trapz(df_degree_preserved['vc'].values, 
                          df_degree_preserved['intensity'].values / 100)
            f.write(f"  Area under curve: {auc:.4f}\n")
        
        # Key insights
        f.write("\n" + "="*70 + "\n")
        f.write("KEY INSIGHTS\n")
        f.write("-"*70 + "\n")
        
        # Compare degree-preserved vs random net
        df_net = df_standard[df_standard['noise_type'] == 'net']
        if not df_net.empty and not df_degree_preserved.empty:
            net_change = (df_net[df_net['degree'] == 100.0]['vc'].values[0] - 
                         df_net[df_net['degree'] == 0.0]['vc'].values[0])
            deg_pres_change = (df_degree_preserved[df_degree_preserved['intensity'] == 100.0]['vc'].values[0] - 
                              df_degree_preserved[df_degree_preserved['intensity'] == 0.0]['vc'].values[0])
            
            f.write(f"\nRandom network shuffle impact: {net_change:+.4f}\n")
            f.write(f"Degree-preserved shuffle impact: {deg_pres_change:+.4f}\n")
            f.write(f"Difference (topology effect): {abs(net_change - deg_pres_change):.4f}\n")
            
            if abs(deg_pres_change) < abs(net_change):
                f.write("\n→ Degree distribution is IMPORTANT for VC metric\n")
                f.write("  (Metric is more robust when degree is preserved)\n")
            else:
                f.write("\n→ Network topology beyond degree is IMPORTANT for VC metric\n")
                f.write("  (Metric changes even when degree is preserved)\n")
        
        f.write("\n" + "="*70 + "\n")
    
    print(f"Summary report saved to: {report_file}")


if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(description='Comprehensive robustness analysis comparison')
    parser.add_argument('--prediction', type=str, required=True,
                       help='Path to GRN prediction h5ad file')
    parser.add_argument('--evaluation_data', type=str, required=True,
                       help='Path to evaluation data h5ad file')
    parser.add_argument('--output_dir', type=str, required=True,
                       help='Base output directory for all results')
    parser.add_argument('--n_top_genes', type=int, default=3000)
    parser.add_argument('--intensities', type=float, nargs='+',
                       default=[0.0, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0],
                       help='Permutation intensities to test')
    
    args = parser.parse_args()
    
    # Create parameter dictionary
    par = {
        'prediction': args.prediction,
        'evaluation_data': args.evaluation_data,
        'base_output_dir': args.output_dir,
        'intensities': args.intensities,
        'n_top_genes': args.n_top_genes
    }
    
    # Run comprehensive analysis
    run_all_analyses(par)
