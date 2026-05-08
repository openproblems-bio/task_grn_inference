"""
Combine trace files, score_uns.yaml, and dataset_uns.yaml from individual dataset result folders.
"""

import os
import shutil
import yaml
import pandas as pd
from pathlib import Path
from collections import OrderedDict
import sys
import argparse

# Add paths for imports
RESULTS_DIR = Path('resources/results')
BENCHMARK_DIR = Path('resources/results/benchmark')

from src.utils.config import DATASETS, METRICS

def combine_results(local_run=False):
    """Combine results from individual dataset folders into all_new folder."""
    
    save_dir = BENCHMARK_DIR / 'all_new'
    
    # Create output directory
    os.makedirs(save_dir, exist_ok=True)
    
    print(f"Processing {len(DATASETS)} datasets: {DATASETS}")
    
    # 1. Copy trace_merged.txt file from 'op' dataset only
    print("\n1. Copying trace_merged.txt file (using only 'op' dataset)...")
    
    trace_dataset = 'op'
    trace_path = BENCHMARK_DIR / trace_dataset / 'trace_merged.txt'
    
    if trace_path.exists():
        print(f"  Reading trace from {trace_dataset}...")
        trace_df = pd.read_csv(trace_path, sep='\t')
        trace_df = trace_df.drop_duplicates(subset=['name'])
        # Save trace
        output_trace = save_dir / 'trace.csv'
        trace_df.to_csv(output_trace, sep='\t', index=False)
        print(f"  Saved trace to {output_trace}")
        print(f"  Total unique entries: {len(trace_df)}")
    else:
        print(f"  Warning: {trace_path} not found!")
    
    # 2. Combine score files
    print("\n2. Combining score files...")
    
    if local_run:
        # Copy all_scores.csv directly (local run mode)
        print("  Using local run mode - copying all_scores.csv...")
        scores_path = BENCHMARK_DIR / 'all_scores.csv'
        assert scores_path.exists(), f"Error: {scores_path} not found!"
        shutil.copy(scores_path, 'docs/source/data/all_scores.csv')
        
    else:
        # Read from individual score_uns.yaml files (AWS mode)
        print("  Using AWS mode - reading from individual score_uns.yaml files...")
        merged_scores = []
        
        for dataset in DATASETS:
            score_path = BENCHMARK_DIR / dataset / 'score_uns.yaml'
            
            if not score_path.exists():
                print(f"  Warning: {score_path} not found, skipping...")
                continue
            
            print(f"  Reading scores from {dataset}...")
            with open(score_path, 'r') as f:
                data = yaml.safe_load(f)
                
                # Filter out None and missing entries
                if data:
                    if isinstance(data, dict):
                        if 'missing' not in str(data):
                            merged_scores.append(data)
                    elif isinstance(data, list):
                        valid_data = [d for d in data if d is not None and 'missing' not in str(d)]
                        merged_scores.extend(valid_data)
        
        if merged_scores:
            output_scores = save_dir / 'score_uns.yaml'
            with open(output_scores, 'w') as f:
                yaml.dump(merged_scores, f)
            print(f"  Saved combined scores to {output_scores}")
            print(f"  Total score entries: {len(merged_scores)}")
        else:
            print("  No score files found!")
    
    # 3. Combine dataset_uns.yaml files
    print("\n3. Combining dataset_uns.yaml files...")
    merged_datasets = []
    
    for dataset in DATASETS:
        dataset_path = BENCHMARK_DIR / dataset / 'dataset_uns.yaml'
        
        if not dataset_path.exists():
            print(f"  Warning: {dataset_path} not found, skipping...")
            continue
        
        print(f"  Reading dataset info from {dataset}...")
        with open(dataset_path, 'r') as f:
            data = yaml.safe_load(f)
            
            if data:
                if isinstance(data, dict):
                    merged_datasets.append(data)
                elif isinstance(data, list):
                    merged_datasets.extend(data)
    
    if merged_datasets:
        output_datasets = save_dir / 'dataset_uns.yaml'
        with open(output_datasets, 'w') as f:
            yaml.dump(merged_datasets, f)
        print(f"  Saved combined datasets to {output_datasets}")
        print(f"  Total dataset entries: {len(merged_datasets)}")
    else:
        print("  No dataset files found!")
    
    # 4. Copy method_configs.yaml and metric_configs.yaml from first dataset
    print("\n4. Copying config files...")
    config_files = ['method_configs.yaml', 'metric_configs.yaml']
    
    for config_file in config_files:
        # Find first dataset that has this file
        for dataset in DATASETS:
            src_path = BENCHMARK_DIR / dataset / config_file
            if src_path.exists():
                dst_path = save_dir / config_file
                # Just copy the file directly without parsing (to avoid Groovy YAML issues)
                shutil.copy(src_path, dst_path)
                print(f"  Copied {config_file} from {dataset}")
                break
        else:
            print(f"  Warning: {config_file} not found in any dataset folder")
    
    print("\n✓ All files combined successfully!")
    print(f"Output directory: {save_dir}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Combine results from individual dataset folders'
    )
    parser.add_argument(
        '--aws_run',
        dest='local_run',
        action='store_false',
        help='Use AWS run mode - read scores from individual score_uns.yaml files instead of all_scores.csv'
    )
    parser.set_defaults(local_run=True)
    
    args = parser.parse_args()
    combine_results(local_run=args.local_run)
