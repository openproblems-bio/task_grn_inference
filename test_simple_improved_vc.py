#!/usr/bin/env python3
"""
Test the simplified improved VC implementation
"""

import sys
import os

# Add paths
sys.path.append('src/metrics/vc')
sys.path.append('src/utils')

# Test parameters
test_par = {
    'prediction': 'resources/results/op/op.positive_control.positive_control.prediction.h5ad',
    'evaluation_data': 'resources/grn_benchmark/evaluation_data/op_bulk.h5ad',
    'layer': 'lognorm',
    'max_n_links': 50000,
    'num_workers': 20,
    'tf_all': 'resources/grn_benchmark/prior/tf_all.csv',    
    'score': 'output/test_simple_improved_vc.h5ad',
    'genes_n': 5000
}

if __name__ == "__main__":
    print("=== Testing Simple Improved VC Implementation ===")
    print(f"Prediction file: {test_par['prediction']}")
    print(f"Evaluation file: {test_par['evaluation_data']}")
    print()
    
    try:
        # Import and run simple improved VC
        from simple_improved_helper import main as simple_improved_main
        
        result = simple_improved_main(test_par)
        
        print("=== Results ===")
        print(result)
        
        # Save results
        from util import format_save_score
        import anndata as ad
        
        dataset_id = ad.read_h5ad(test_par['evaluation_data'], backed='r').uns['dataset_id']
        method_id = ad.read_h5ad(test_par['prediction'], backed='r').uns['method_id']
        
        format_save_score(result, method_id, dataset_id, test_par['score'])
        print(f"Results saved to: {test_par['score']}")
        
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()