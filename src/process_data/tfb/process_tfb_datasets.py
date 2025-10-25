#!/usr/bin/env python3
"""
Process TF binding (TFB) datasets for multiple cell types to create Ground Truth GRNs.

This script uses the existing functions from build_grn.py to process
TF binding files in the resources/datasets_raw/tfb/ directory.
Supports multiple cell types including PBMC, K562, HEK293, and HCT116.
"""

import os
import pandas as pd
import argparse
from build_grn import (
    read_tfb_file, 
    build_grn, 
    print_grn_summary
)

def get_available_cell_types(base_dir=None):
    """Get available cell types from the TFB datasets."""
    if base_dir is None:
        # Try to find the correct base directory
        current_dir = os.getcwd()
        if 'task_grn_inference' in current_dir:
            # Navigate to the task_grn_inference root
            while not os.path.basename(current_dir) == 'task_grn_inference' and current_dir != '/':
                current_dir = os.path.dirname(current_dir)
            tfb_base_dir = os.path.join(current_dir, "resources", "datasets_raw", "tfb")
        else:
            tfb_base_dir = "resources/datasets_raw/tfb/"
    else:
        tfb_base_dir = base_dir
    
    cell_types = set()
    
    # Check each data source for available cell types
    for source in ["unibind", "chipatlas", "remap2022"]:
        source_dir = os.path.join(tfb_base_dir, source)
        if os.path.exists(source_dir):
            for file in os.listdir(source_dir):
                if file.endswith('.bed') and '_' in file:
                    # Extract cell type from filename (e.g., "unibind_pbmc.bed" -> "pbmc")
                    parts = file.replace('.bed', '').split('_')
                    if len(parts) > 1:
                        cell_type = parts[-1].lower()
                        cell_types.add(cell_type)
    
    return sorted(list(cell_types))

def process_tfb_datasets_for_cell_type(cell_type, base_dir=None):
    """Process all TFB datasets for a specific cell type and create GT GRNs."""
    
    # Configuration - handle relative paths
    if base_dir is None:
        current_dir = os.getcwd()
        if 'task_grn_inference' in current_dir:
            # Navigate to the task_grn_inference root
            while not os.path.basename(current_dir) == 'task_grn_inference' and current_dir != '/':
                current_dir = os.path.dirname(current_dir)
            tfb_base_dir = os.path.join(current_dir, "resources", "datasets_raw", "tfb")
            tss_file = os.path.join(current_dir, "resources", "supp_data", "tss_h38.bed")
            output_dir = os.path.join(current_dir, "resources", "grn_benchmark", "ground_truth")
        else:
            tfb_base_dir = "resources/datasets_raw/tfb/"
            tss_file = "resources/supp_data/tss_h38.bed"
            output_dir = "resources/grn_benchmark/ground_truth/"
    else:
        tfb_base_dir = os.path.join(base_dir, "resources", "datasets_raw", "tfb")
        tss_file = os.path.join(base_dir, "resources", "supp_data", "tss_h38.bed")
        output_dir = os.path.join(base_dir, "resources", "grn_benchmark", "ground_truth")
    
    max_workers = 10  # Reduce workers to avoid memory issues
    
    cell_type_lower = cell_type.lower()
    cell_type_upper = cell_type.upper()
    
    # TFB datasets to process for the specified cell type
    datasets = [
        {"name": f"remap2022_{cell_type_lower}", "file": f"remap2022/remap2022_{cell_type_lower}.bed", "source": "remap2022"},
        {"name": f"chipatlas_{cell_type_lower}", "file": f"chipatlas/chipatlas_{cell_type_lower}.bed", "source": "chipatlas"},
        {"name": f"unibind_{cell_type_lower}", "file": f"unibind/unibind_{cell_type_lower}.bed", "source": "unibind"},
    ]
    
    print(f"=== Processing TFB datasets for {cell_type_upper} ===")
    
    # Check if TSS file exists
    if not os.path.exists(tss_file):
        print(f"Error: TSS file not found: {tss_file}")
        return {}
    
    all_results = {}
    
    for dataset in datasets:
        tfb_file = os.path.join(tfb_base_dir, dataset["file"])
        
        if not os.path.exists(tfb_file):
            print(f"Warning: TFB file not found: {tfb_file}")
            continue
            
        print(f"\n--- Processing {dataset['name']} ---")
        print(f"File: {tfb_file}")
        
        # Read TFB file
        print("Reading TF binding data...")
        peaks_df = read_tfb_file(tfb_file, dataset["source"])
        print(f"Loaded {len(peaks_df)} TF binding sites")
        print(f"Unique TFs: {peaks_df['tf'].nunique()}")
        
        # Build GRN
        print("Building GRN...")
        grn = build_grn(peaks_df, tss_file, max_workers=max_workers)
        
        # Save results
        output_file = f"{output_dir}/{cell_type_upper}_{dataset['name']}.csv"
        os.makedirs(output_dir, exist_ok=True)
        grn.to_csv(output_file, index=False)
        
        # Store results
        all_results[dataset['name']] = {
            'grn': grn,
            'output_path': output_file,
            'num_edges': len(grn),
            'num_tfs': grn['source'].nunique(),
            'num_targets': grn['target'].nunique()
        }
        
        print(f"GRN saved to: {output_file}")
        print(f"Edges: {len(grn)}, TFs: {grn['source'].nunique()}, Targets: {grn['target'].nunique()}")
    
    # Print summary
    if all_results:
        print_grn_summary(all_results, f"TFB {cell_type_upper} Datasets")
        
        # Create a combined summary file
        summary_file = f"{output_dir}/tfb_{cell_type_lower}_summary.txt"
        with open(summary_file, 'w') as f:
            f.write(f"TFB {cell_type_upper} Datasets GRN Summary\n")
            f.write("=" * 50 + "\n")
            for name, result in all_results.items():
                f.write(f"\n{name}:\n")
                f.write(f"  Edges: {result['num_edges']}\n")
                f.write(f"  TFs: {result['num_tfs']}\n")
                f.write(f"  Target genes: {result['num_targets']}\n")
                f.write(f"  Output: {result['output_path']}\n")
        
        print(f"\nSummary saved to: {summary_file}")
    else:
        print(f"\nNo datasets were processed successfully for {cell_type_upper}.")
    
    print(f"\n=== Processing completed for {cell_type_upper} ===")
    return all_results

def process_all_tfb_datasets(cell_types=None, base_dir=None):
    """Process all TFB datasets for specified cell types and create GT GRNs."""
    
    if cell_types is None:
        # Get all available cell types
        available_cell_types = get_available_cell_types(base_dir)
        print(f"Available cell types: {', '.join(available_cell_types)}")
        cell_types = available_cell_types
    
    all_cell_type_results = {}
    
    for cell_type in cell_types:
        print(f"\n{'=' * 60}")
        print(f"Processing cell type: {cell_type.upper()}")
        print(f"{'=' * 60}")
        
        results = process_tfb_datasets_for_cell_type(cell_type, base_dir)
        if results:
            all_cell_type_results[cell_type] = results
    
    # Print overall summary
    if all_cell_type_results:
        print(f"\n{'=' * 60}")
        print("OVERALL SUMMARY")
        print(f"{'=' * 60}")
        for cell_type, results in all_cell_type_results.items():
            print(f"\n{cell_type.upper()}:")
            for dataset_name, result in results.items():
                print(f"  {dataset_name}: {result['num_edges']} edges, {result['num_tfs']} TFs, {result['num_targets']} targets")
    
    return all_cell_type_results


def main():
    """Main function to handle command line arguments."""
    parser = argparse.ArgumentParser(description='Process TF binding datasets for multiple cell types')
    parser.add_argument('--cell-types', nargs='*', 
                        help='Cell types to process (e.g., pbmc k562 hek293). If not specified, all available cell types will be processed.')
    parser.add_argument('--list-available', action='store_true',
                        help='List available cell types and exit')
    parser.add_argument('--base-dir', type=str,
                        help='Base directory containing the task_grn_inference project')
    
    args = parser.parse_args()
    
    if args.list_available:
        available_cell_types = get_available_cell_types(args.base_dir)
        print("Available cell types:")
        for cell_type in available_cell_types:
            print(f"  - {cell_type}")
        return
    
    cell_types = args.cell_types if args.cell_types else None
    process_all_tfb_datasets(cell_types, args.base_dir)

if __name__ == '__main__':
    main()