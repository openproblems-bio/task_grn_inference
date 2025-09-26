#!/usr/bin/env python3
"""
Quick comparison between original and improved VC on key methods
"""

import sys
import os
import time
import pandas as pd

# Add paths
sys.path.append('src/metrics/vc')
sys.path.append('src/utils')

def test_method(method_name, use_simple_improved=False):
    """Test a single method and return the score and runtime."""
    
    test_par = {
        'prediction': f'resources/results/op/op.{method_name}.{method_name}.prediction.h5ad',
        'evaluation_data': 'resources/grn_benchmark/evaluation_data/op_bulk.h5ad',
        'layer': 'lognorm',
        'max_n_links': 50000,
        'num_workers': 20,
        'tf_all': 'resources/grn_benchmark/prior/tf_all.csv',    
        'score': f'output/test_{method_name}.h5ad',
        'genes_n': 1000  # Reduced for faster testing
    }
    
    if not os.path.exists(test_par['prediction']):
        return None, None, "File not found"
    
    try:
        start_time = time.time()
        
        if use_simple_improved:
            # Temporarily modify genes_n in the dataset loading
            from simple_improved_helper import main as test_main
        else:
            from helper import main as test_main
        
        result = test_main(test_par)
        runtime = time.time() - start_time
        
        if 'vc' in result.columns:
            score = float(result['vc'].iloc[0])
        else:
            score = None
            
        return score, runtime, "Success"
        
    except Exception as e:
        runtime = time.time() - start_time
        return None, runtime, str(e)[:100]

def main():
    print("=== Quick VC Implementation Comparison ===")
    print("Testing key methods with reduced gene set for speed...")
    print()
    
    # Test methods
    test_methods = [
        'positive_control',  # Should have high performance
        'negative_control',  # Should have low performance  
        'pearson_corr',      # Baseline method
    ]
    
    results = []
    
    for method in test_methods:
        print(f"Testing {method}...")
        
        # Test original
        print(f"  Running original VC...")
        orig_score, orig_time, orig_status = test_method(method, use_simple_improved=False)
        
        # Test improved (skip for now due to speed issues)
        # print(f"  Running improved VC...")
        # improved_score, improved_time, improved_status = test_method(method, use_simple_improved=True)
        improved_score, improved_time, improved_status = None, None, "Skipped (too slow)"
        
        results.append({
            'method': method,
            'original_score': orig_score,
            'original_time': orig_time,
            'original_status': orig_status,
            'improved_score': improved_score,
            'improved_time': improved_time,
            'improved_status': improved_status
        })
        
        print(f"    Original: score={orig_score}, time={orig_time:.1f}s, status={orig_status}")
        print(f"    Improved: score={improved_score}, time={improved_time}, status={improved_status}")
        print()
    
    # Summary
    print("=== Summary ===")
    df = pd.DataFrame(results)
    print(df[['method', 'original_score', 'original_time', 'original_status']])
    
    # Analysis
    print("\n=== Analysis ===")
    original_scores = [r['original_score'] for r in results if r['original_score'] is not None]
    if original_scores:
        print(f"Original VC scores: {original_scores}")
        print(f"Score range: {min(original_scores):.3f} to {max(original_scores):.3f}")
        
        # Check if positive_control > negative_control (as expected)
        pos_score = next((r['original_score'] for r in results if r['method'] == 'positive_control'), None)
        neg_score = next((r['original_score'] for r in results if r['method'] == 'negative_control'), None)
        
        if pos_score is not None and neg_score is not None:
            if pos_score > neg_score:
                print("✓ Positive control outperforms negative control (good sign)")
            else:
                print("✗ Positive control does NOT outperform negative control (bad sign)")
        
        # Check for extreme values
        if any(abs(score) > 100 for score in original_scores):
            print("⚠️  Extreme scores detected - indicates numerical instability")
        else:
            print("✓ Scores are in reasonable range")
    else:
        print("No valid scores obtained")

if __name__ == "__main__":
    main()