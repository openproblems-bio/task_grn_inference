#!/usr/bin/env python3
"""
Test per-sample prediction VC implementation
"""

import sys
import os
import time

# Add paths
sys.path.append('src/metrics/vc')
sys.path.append('src/utils')

def test_per_sample_vc(method_name):
    """Test the per-sample prediction VC implementation."""
    
    test_par = {
        'prediction': f'resources/results/op/op.{method_name}.{method_name}.prediction.h5ad',
        'evaluation_data': 'resources/grn_benchmark/evaluation_data/op_bulk.h5ad',
        'layer': 'lognorm',
        'max_n_links': 50000,
        'num_workers': 20,
        'tf_all': 'resources/grn_benchmark/prior/tf_all.csv',    
        'score': f'output/test_per_sample_{method_name}.h5ad',
        'genes_n': 5000
    }
    
    if not os.path.exists(test_par['prediction']):
        return None, None, "File not found"
    
    try:
        print(f"Testing {method_name} with per-sample prediction VC...")
        start_time = time.time()
        
        from ultra_fast_helper import main as per_sample_main
        result = per_sample_main(test_par)
        
        runtime = time.time() - start_time
        
        if 'vc' in result.columns:
            score = float(result['vc'].iloc[0])
        else:
            score = None
            
        return score, runtime, "Success"
        
    except Exception as e:
        runtime = time.time() - start_time
        return None, runtime, str(e)[:200]

def main():
    print("=== Per-Sample Prediction VC Test ===")
    print("Architecture: Predict all genes simultaneously for each sample")
    print("Key difference: Single forward pass per sample (not per gene)")
    print()
    
    # Test key methods
    test_methods = [
        'positive_control',  # Should score high
        'negative_control',  # Should score low  
        'pearson_corr',      # Baseline
    ]
    
    results = []
    total_start = time.time()
    
    for method in test_methods:
        score, runtime, status = test_per_sample_vc(method)
        results.append({
            'method': method,
            'score': score,
            'runtime': runtime,
            'status': status
        })
        
        print(f"{method:20s} | Score: {score:8.4f} | Time: {runtime:6.1f}s | Status: {status}")
    
    total_time = time.time() - total_start
    print(f"\nTotal time: {total_time:.1f}s")
    
    # Analysis
    print("\n=== Analysis ===")
    valid_results = [r for r in results if r['score'] is not None]
    
    if len(valid_results) >= 2:
        scores = [r['score'] for r in valid_results]
        times = [r['runtime'] for r in valid_results]
        
        print(f"Score range: {min(scores):.4f} to {max(scores):.4f}")
        print(f"Average time per method: {sum(times)/len(times):.1f}s")
        
        # Check ranking
        pos_score = next((r['score'] for r in results if r['method'] == 'positive_control'), None)
        neg_score = next((r['score'] for r in results if r['method'] == 'negative_control'), None)
        
        if pos_score is not None and neg_score is not None:
            if pos_score > neg_score:
                print("✓ Positive control > Negative control (GOOD)")
                improvement = "BETTER ranking than before!"
            else:
                print("✗ Positive control ≤ Negative control (BAD)")  
                improvement = "Same ranking issue as before"
                
            print(f"  Positive control: {pos_score:.4f}")
            print(f"  Negative control: {neg_score:.4f}")
            print(f"  Difference: {pos_score - neg_score:.4f}")
            print(f"  → {improvement}")
        
        # Score magnitude analysis
        max_abs_score = max(abs(s) for s in scores)
        if max_abs_score < 0.1:
            print("⚠️  Scores very close to zero - may indicate metric issues")
        elif max_abs_score > 1.0:
            print("⚠️  Large scores - may indicate numerical instability")
        else:
            print("✓ Scores in reasonable range")
    else:
        print("❌ No valid results obtained")

if __name__ == "__main__":
    main()