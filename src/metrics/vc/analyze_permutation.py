"""
Analyze how vc metric scores change when GRN structure is permuted.
Tests robustness by shuffling TF-gene pairs in the network.
"""

import numpy as np
import pandas as pd
import anndata as ad
from helper import main as main_vc
import sys

def permute_grn_edges(net_df, permutation_degree=1.0, seed=42):
    """
    Permute the GRN by shuffling TF-gene pairs.
    
    Args:
        net_df: DataFrame with 'source', 'target', 'weight' columns
        permutation_degree: Fraction of edges to permute (0.0 to 1.0)
        seed: Random seed
        
    Returns:
        Permuted network DataFrame
    """
    np.random.seed(seed)
    net_permuted = net_df.copy()
    
    n_edges = len(net_df)
    n_to_permute = int(n_edges * permutation_degree)
    
    if n_to_permute == 0:
        return net_permuted
    
    # Select random edges to permute
    indices = np.random.choice(n_edges, n_to_permute, replace=False)
    
    # Shuffle the targets for selected edges
    targets = net_permuted.loc[indices, 'target'].values
    np.random.shuffle(targets)
    net_permuted.loc[indices, 'target'] = targets
    
    return net_permuted


def run_permutation_analysis(par, permutation_degrees=[0.0, 0.2, 0.5, 0.8, 1.0]):
    """
    Run vc metric with different levels of GRN permutation.
    
    Args:
        par: Parameters dict with paths to data and prediction
        permutation_degrees: List of permutation fractions to test
        
    Returns:
        DataFrame with results for each permutation level
    """
    results = []
    
    # Load original prediction - read as DataFrame from var
    pred_adata = ad.read_h5ad(par['prediction'])
    original_net = pred_adata.var.copy()
    
    if 'source' not in original_net.columns or 'target' not in original_net.columns:
        # Try to get from varm
        if 'skeleton' in pred_adata.varm:
            original_net = pd.DataFrame(pred_adata.varm['skeleton'], columns=['source', 'target', 'weight'])
        else:
            raise ValueError("Cannot find GRN data in prediction file")
    
    for degree in permutation_degrees:
        print(f"\n{'='*60}")
        print(f"Testing permutation degree: {degree:.0%}")
        print(f"{'='*60}")
        
        # Permute network
        if degree > 0:
            permuted_net = permute_grn_edges(original_net, degree, seed=42)
            
            # Create permuted anndata with same structure
            permuted_adata = ad.AnnData(uns=pred_adata.uns.copy())
            permuted_adata.var = permuted_net
            
            permuted_path = par['prediction'].replace('.h5ad', f'_perm{int(degree*100)}.h5ad')
            permuted_adata.write_h5ad(permuted_path)
            
            # Update parameter
            par_temp = par.copy()
            par_temp['prediction'] = permuted_path
        else:
            par_temp = par.copy()
        
        # Run evaluation
        try:
            df_result = main_vc(par_temp)
            r2_score = df_result['vc'].values[0]
            
            results.append({
                'permutation_degree': degree,
                'r2_score': r2_score
            })
            
            print(f"R² score: {r2_score:.4f}")
            
        except Exception as e:
            print(f"Error at permutation degree {degree}: {e}")
            import traceback
            traceback.print_exc()
            results.append({
                'permutation_degree': degree,
                'r2_score': np.nan
            })
    
    return pd.DataFrame(results)


if __name__ == '__main__':
    # Configuration
    dataset = 'op'  # Changed from 'replogle' to 'op'
    method = 'grnboost'
    
    par = {
        'prediction': f'resources/results/{dataset}/{dataset}.{method}.{method}.prediction.h5ad',
        'evaluation_data': f'resources/grn_benchmark/evaluation_data/{dataset}_bulk.h5ad',
        'score': f'output/vc/vc_{dataset}_{method}_permutation.h5ad'
    }
    
    print(f"Dataset: {dataset}")
    print(f"Method: {method}")
    
    # Run permutation analysis
    results_df = run_permutation_analysis(
        par,
        permutation_degrees=[0.0, 0.2, 0.5, 0.8, 1.0]
    )
    
    # Save results
    output_path = f'output/vc/permutation_analysis_{dataset}_{method}.csv'
    results_df.to_csv(output_path, index=False)
    
    print(f"\n{'='*60}")
    print("Permutation Analysis Results")
    print(f"{'='*60}")
    print(results_df.to_string(index=False))
    print(f"\nResults saved to: {output_path}")
    
    # Calculate score degradation
    if len(results_df) > 0 and not results_df['r2_score'].isna().all():
        original_score = results_df.loc[results_df['permutation_degree'] == 0.0, 'r2_score'].values[0]
        final_score = results_df.loc[results_df['permutation_degree'] == 1.0, 'r2_score'].values[0]
        
        if not np.isnan(original_score) and not np.isnan(final_score):
            degradation = (original_score - final_score) / original_score * 100
            print(f"\nScore degradation with 100% permutation: {degradation:.1f}%")
            print(f"Original R²: {original_score:.4f}")
            print(f"Permuted R²: {final_score:.4f}")
