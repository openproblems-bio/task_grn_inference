#!/usr/bin/env python3
"""
Process Nextflow trace files to extract and merge successful tasks.

This script reads all trace_*.txt files from dataset directories, filters for successful
tasks (COMPLETED status), keeps only the most recent run per method (including
all related tasks from that run), and creates a merged output file for each dataset.

Usage:
    python process_trace_files.py [--dataset DATASET] [--input_base INPUT_BASE]

Args:
    --dataset: Specific dataset to process (e.g., 'op', 'replogle'). If not provided, processes all datasets.
    --input_base: Base directory containing dataset result folders
                  (default: /home/jnourisa/projs/ongoing/task_grn_inference/resources/results/)
"""

import os
import glob
import sys
from datetime import datetime
from pathlib import Path
from typing import List, Dict, Tuple
from collections import defaultdict
import re
import argparse

# Import datasets from config
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from utils.config import DATASETS


def parse_trace_file(filepath: str) -> Tuple[List[Dict], datetime]:
    """
    Parse a trace file and extract all tasks with COMPLETED status.
    
    Args:
        filepath: Path to the trace file
        
    Returns:
        Tuple of (list of completed task dictionaries, file date)
    """
    completed_tasks = []
    
    # Extract date from filename (trace_YYYY-MM-DD.txt or trace.txt)
    filename = os.path.basename(filepath)
    date_match = re.search(r'trace_(\d{4}-\d{2}-\d{2})\.txt', filename)
    if date_match:
        file_date = datetime.strptime(date_match.group(1), '%Y-%m-%d')
    else:
        # If no date in filename, use file modification time
        file_date = datetime.fromtimestamp(os.path.getmtime(filepath))
    
    try:
        with open(filepath, 'r') as f:
            # Read header
            header_line = f.readline().strip()
            if not header_line:
                return completed_tasks, file_date
                
            headers = header_line.split('\t')
            
            # Process data lines
            for line_num, line in enumerate(f, start=2):
                line = line.strip()
                if not line:
                    continue
                    
                fields = line.split('\t')
                
                # Ensure we have enough fields
                if len(fields) < len(headers):
                    continue
                
                # Create dictionary from fields
                task_dict = {headers[i]: fields[i] if i < len(fields) else '' 
                            for i in range(len(headers))}
                
                # Filter for COMPLETED tasks only
                status = task_dict.get('status', '')
                if status == 'COMPLETED':
                    task_dict['source_file'] = filename
                    task_dict['file_date'] = file_date
                    completed_tasks.append(task_dict)
                    
    except Exception as e:
        print(f"Error processing {filepath}: {e}", file=sys.stderr)
        
    return completed_tasks, file_date


def extract_method_name(task_name: str) -> str:
    """
    Extract the method/component name from a task name.
    
    Args:
        task_name: Full task name from the trace file
        
    Returns:
        Extracted method name (e.g., "scenic_op.regression") or empty string
    """
    # Pattern to extract the method from task names like:
    # "run_grn_evaluation:run_wf:runEachWf:regression:processWf:regression_process (scenic_op.regression)"
    # Extract the part in parentheses
    match = re.search(r'\(([^)]+)\)', task_name)
    if match:
        return match.group(1)
    return ''


def group_tasks_by_method_and_run(all_tasks: List[Dict]) -> Dict[str, Dict[str, List[Dict]]]:
    """
    Group tasks by method identifier and then by run (file_date + source_file).
    
    Args:
        all_tasks: List of all task dictionaries
        
    Returns:
        Nested dictionary: {method_id: {run_id: [tasks]}}
    """
    # Group by method, then by run
    method_runs = defaultdict(lambda: defaultdict(list))
    
    for task in all_tasks:
        task_name = task.get('name', '')
        method = extract_method_name(task_name)
        
        if not method:
            # Skip tasks without a method identifier
            continue
        
        # Create a run identifier from file date and source file
        run_id = f"{task['file_date'].strftime('%Y-%m-%d')}::{task['source_file']}"
        
        method_runs[method][run_id].append(task)
    
    return method_runs


def get_most_recent_run_per_method(method_runs: Dict[str, Dict[str, List[Dict]]]) -> List[Dict]:
    """
    Keep only the most recent run for each method (all tasks from that run).
    
    Args:
        method_runs: Nested dictionary {method_id: {run_id: [tasks]}}
        
    Returns:
        List of all tasks from the most recent run of each method
    """
    most_recent_tasks = []
    
    for method, runs in method_runs.items():
        # Sort runs by file date (most recent first)
        sorted_runs = sorted(
            runs.items(),
            key=lambda x: x[1][0]['file_date'],  # Use file_date from first task
            reverse=True
        )
        
        # Get all tasks from the most recent run
        most_recent_run_id, most_recent_tasks_list = sorted_runs[0]
        most_recent_tasks.extend(most_recent_tasks_list)
    
    # Sort the final list by method name and task name for consistency
    most_recent_tasks.sort(key=lambda x: (
        extract_method_name(x.get('name', '')),
        x.get('name', '')
    ))
    
    return most_recent_tasks


def write_merged_file(tasks: List[Dict], output_file: str, headers: List[str]):
    """
    Write the merged successful tasks to an output file.
    
    Args:
        tasks: List of task dictionaries to write
        output_file: Path to the output file
        headers: List of header field names
    """
    # Add source_file to headers if not present
    if 'source_file' not in headers:
        headers = headers + ['source_file']
    
    # Create output directory if it doesn't exist
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    with open(output_file, 'w') as f:
        # Write header
        f.write('\t'.join(headers) + '\n')
        
        # Write tasks
        for task in tasks:
            row = []
            for header in headers:
                value = task.get(header, '')
                # Skip file_date as it's not in original format
                if header == 'file_date':
                    continue
                row.append(str(value))
            f.write('\t'.join(row) + '\n')


def process_dataset(input_dir: str, output_file: str, dataset_name: str) -> Dict:
    """
    Process trace files for a single dataset.
    
    Args:
        input_dir: Directory containing trace files
        output_file: Path for the output merged file
        dataset_name: Name of the dataset being processed
        
    Returns:
        Dictionary with processing statistics
    """
    print(f"\n{'='*60}")
    print(f"Processing dataset: {dataset_name}")
    print(f"{'='*60}")
    print(f"Input directory: {input_dir}")
    print(f"Output file: {output_file}")
    
    # Find all trace files
    trace_pattern = os.path.join(input_dir, 'trace*.txt')
    trace_files = glob.glob(trace_pattern)
    
    if not trace_files:
        print(f"⚠ No trace files found matching pattern: {trace_pattern}")
        return {
            'dataset': dataset_name,
            'status': 'no_files',
            'trace_files': 0,
            'completed_tasks': 0,
            'output_tasks': 0
        }
    
    print(f"Found {len(trace_files)} trace file(s)")
    
    # Process all trace files
    all_completed_tasks = []
    headers = None
    
    for trace_file in sorted(trace_files):
        print(f"  Processing: {os.path.basename(trace_file)}")
        completed_tasks, file_date = parse_trace_file(trace_file)
        
        if completed_tasks:
            print(f"    → Found {len(completed_tasks)} completed task(s)")
            all_completed_tasks.extend(completed_tasks)
            
            # Get headers from first file with completed tasks
            if headers is None and completed_tasks:
                headers = list(completed_tasks[0].keys())
                # Remove file_date from headers as it's internal
                if 'file_date' in headers:
                    headers.remove('file_date')
    
    if not all_completed_tasks:
        print("⚠ No completed tasks found in any trace file.")
        return {
            'dataset': dataset_name,
            'status': 'no_completed',
            'trace_files': len(trace_files),
            'completed_tasks': 0,
            'output_tasks': 0
        }
    
    print(f"\nTotal completed tasks found: {len(all_completed_tasks)}")
    
    # Group by method and run
    method_runs = group_tasks_by_method_and_run(all_completed_tasks)
    print(f"Unique methods found: {len(method_runs)}")
    
    # Keep only most recent run per method
    most_recent_tasks = get_most_recent_run_per_method(method_runs)
    
    print(f"Tasks from most recent runs: {len(most_recent_tasks)}")
    
    # Write merged output
    write_merged_file(most_recent_tasks, output_file, headers)
    
    print(f"✓ Merged file written: {output_file}")
    
    # Get statistics
    source_files = defaultdict(int)
    for task in most_recent_tasks:
        source_files[task.get('source_file', 'unknown')] += 1
    
    return {
        'dataset': dataset_name,
        'status': 'success',
        'trace_files': len(trace_files),
        'completed_tasks': len(all_completed_tasks),
        'unique_methods': len(method_runs),
        'output_tasks': len(most_recent_tasks),
        'source_files': dict(source_files)
    }


def main():
    """Main function to process trace files."""
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description='Process Nextflow trace files and merge successful tasks'
    )
    parser.add_argument(
        '--dataset',
        type=str,
        help='Specific dataset to process (e.g., op, replogle). If not provided, processes all datasets.'
    )
    parser.add_argument(
        '--input_base',
        type=str,
        default='/home/jnourisa/projs/ongoing/task_grn_inference/resources/results/',
        help='Base directory containing dataset result folders'
    )
    
    args = parser.parse_args()
    
    # Determine which datasets to process
    if args.dataset:
        if args.dataset not in DATASETS:
            print(f"Error: Dataset '{args.dataset}' not found in config.")
            print(f"Available datasets: {', '.join(DATASETS)}")
            sys.exit(1)
        datasets_to_process = [args.dataset]
    else:
        datasets_to_process = DATASETS
    
    print(f"Base directory: {args.input_base}")
    print(f"Datasets to process: {', '.join(datasets_to_process)}\n")
    
    # Process each dataset
    results = []
    for dataset in datasets_to_process:
        input_dir = os.path.join(args.input_base, dataset)
        output_file = os.path.join(input_dir, 'trace_merged.txt')
        
        # Check if input directory exists
        if not os.path.exists(input_dir):
            print(f"\n{'='*60}")
            print(f"Dataset: {dataset}")
            print(f"{'='*60}")
            print(f"⚠ Directory not found: {input_dir}")
            print("Skipping...")
            results.append({
                'dataset': dataset,
                'status': 'dir_not_found',
                'trace_files': 0,
                'completed_tasks': 0,
                'output_tasks': 0
            })
            continue
        
        # Process the dataset
        result = process_dataset(input_dir, output_file, dataset)
        results.append(result)
    
    # Print summary
    print(f"\n{'='*60}")
    print("SUMMARY")
    print(f"{'='*60}\n")
    
    successful = [r for r in results if r['status'] == 'success']
    skipped = [r for r in results if r['status'] != 'success']
    
    if successful:
        print(f"✓ Successfully processed {len(successful)} dataset(s):\n")
        for result in successful:
            print(f"  {result['dataset']:20} → {result['output_tasks']:4} tasks "
                  f"(from {result['trace_files']} trace files)")
    
    if skipped:
        print(f"\n⚠ Skipped {len(skipped)} dataset(s):\n")
        for result in skipped:
            status_msg = {
                'dir_not_found': 'directory not found',
                'no_files': 'no trace files',
                'no_completed': 'no completed tasks'
            }.get(result['status'], result['status'])
            print(f"  {result['dataset']:20} → {status_msg}")
    
    print()


if __name__ == '__main__':
    main()
