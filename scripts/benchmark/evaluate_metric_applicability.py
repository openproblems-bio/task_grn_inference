"""
Evaluate which metrics to keep for each dataset based on:
1. std/mean ratio (coefficient of variation) >= threshold (e.g., 0.2)
2. max score obtained is meaningful based on metric-specific thresholds
"""

import os
import yaml
import pandas as pd
import numpy as np
from pathlib import Path
from collections import defaultdict
import argparse

from src.utils.config import DATASETS, METRICS, METRIC_THRESHOLDS, METRICS_DATASETS, surrogate_names

RESULTS_DIR = Path("resources/results")
BENCHMARK_DIR = Path("resources/results/benchmark")
DOCS_IMAGES_DIR = Path("docs/source/images")


parser = argparse.ArgumentParser(
    description='Evaluate metric quality for each dataset'
)
parser.add_argument(
    '--cv_threshold',
    type=float,
    default=0.2,
    help='Threshold for coefficient of variation (max - min) / mean'
)
parser.add_argument(
    '--output',
    type=str,
    default=None,
    help='Output CSV file path (optional)'
)
parser.add_argument(
    '--datasets',
    nargs='+',
    default=None,
    help='List of datasets to evaluate (default: all)'
)
parser.add_argument(
    '--metrics',
    nargs='+',
    default=None,
    help='List of metrics to evaluate (default: all METRICS)'
)

parser.add_argument(
    "--no_local_run",
    action="store_false",
    dest="local_run",
    help="Disable local mode (use score_uns.yaml instead of all_scores.csv)"
)
parser.set_defaults(local_run=True)


args = parser.parse_args()
local_run = args.local_run
# print(f"Arguments: {args}")

def load_scores_from_yaml(dataset):
    """Load scores from score_uns.yaml for a given dataset."""
    score_path = RESULTS_DIR / dataset / "score_uns.yaml"
    
    if not score_path.exists():
        print(f"Warning: {score_path} not found")
        return None
    
    with open(score_path, 'r') as f:
        data = yaml.safe_load(f)
    
    if not data:
        return None
    
    # Convert to DataFrame format: rows = methods, columns = metrics
    scores_dict = defaultdict(dict)
    
    for entry in data:
        if not isinstance(entry, dict):
            continue
        
        method_id = entry.get('method_id')
        metric_ids = entry.get('metric_ids', [])
        metric_values = entry.get('metric_values', [])
        
        if not method_id or not metric_ids:
            continue
        
        # Store metric values for this method
        for metric_id, metric_value in zip(metric_ids, metric_values):
            try:
                scores_dict[method_id][metric_id] = float(metric_value)
            except (ValueError, TypeError):
                continue
    
    if not scores_dict:
        return None
    
    # Convert to DataFrame
    df = pd.DataFrame.from_dict(scores_dict, orient='index')
    return df


def evaluate_metric_for_dataset(scores_df, metric, cv_threshold=0.2, nc_score=None, pc_score=None):
    """
    Evaluate if a metric should be kept for a dataset.

    Two criteria must both pass:
      (i)  CV = (max - min) / mean >= cv_threshold (variability across methods)
      (ii) max score >= max(nc_score, global_threshold) where global_threshold
           is the metric-specific floor from METRIC_THRESHOLDS. The global floor
           ensures metrics are not kept when all methods perform at a trivially
           low absolute level, even if one nominally outperforms the NC.

    Args:
        scores_df: DataFrame with methods as rows and metrics as columns
        metric: The metric name to evaluate
        cv_threshold: Minimum CV threshold (default 0.2)
        nc_score: Negative control score; used as per-dataset threshold when
                  stricter than the global floor.
        pc_score: Not used in this criterion (kept for API compatibility).

    Returns:
        dict with evaluation results
    """
    if metric not in scores_df.columns:
        return {
            'metric': metric,
            'present': False,
            'reason': 'Metric not computed for this dataset'
        }

    values = scores_df[metric].dropna()

    if len(values) < 2:
        return {
            'metric': metric,
            'present': True,
            'keep': False,
            'reason': 'Insufficient data points (<2 methods)',
            'n_methods': len(values),
            'mean': np.nan,
            'std': np.nan,
            'cv': np.nan,
            'max': np.nan
        }

    mean_val = values.mean()
    std_val = values.std()
    max_val = values.max()
    min_val = values.min()

    # (i) Coefficient of variation
    cv = (max_val - min_val) / mean_val if mean_val != 0 else np.inf
    cv_passes = cv >= cv_threshold

    # (ii) Max score must exceed both NC and the global floor
    global_threshold = METRIC_THRESHOLDS.get(metric, 0.1)
    if nc_score is not None and not np.isnan(nc_score):
        threshold = max(nc_score, global_threshold)
    else:
        threshold = global_threshold
    max_passes = max_val >= threshold

    # (iii) Positive control must outperform negative control
    pc_nc_valid = (
        pc_score is not None and not np.isnan(pc_score) and
        nc_score is not None and not np.isnan(nc_score)
    )
    pc_passes = (pc_score > nc_score) if pc_nc_valid else True  # skip if controls unavailable

    keep = cv_passes and max_passes and pc_passes

    reasons = []
    if not cv_passes:
        reasons.append(f'Low variability (CV={cv:.3f} < {cv_threshold})')
    if not max_passes:
        reasons.append(f'Low max score (max={max_val:.3f} < threshold={threshold:.3f})')
    if not pc_passes:
        reasons.append(f'PC not > NC (PC={pc_score:.4f}, NC={nc_score:.4f})')

    reason = 'Passes all criteria' if keep else '; '.join(reasons)

    return {
        'metric': metric,
        'present': True,
        'keep': keep,
        'reason': reason,
        'n_methods': len(values),
        'mean': mean_val,
        'std': std_val,
        'cv': cv,
        'max': max_val,
        'threshold': threshold,
        'global_threshold': global_threshold,
        'pc_score': pc_score if pc_nc_valid else np.nan,
        'nc_score': nc_score if (nc_score is not None and not np.isnan(nc_score)) else np.nan,
        'cv_threshold': cv_threshold
    }


def evaluate_all_datasets(datasets=None, metrics=None, cv_threshold=0.2, output_file=None):
    """
    Evaluate all metrics for all datasets.
    
    Args:
        datasets: List of datasets to evaluate (default: all DATASETS)
        metrics: List of metrics to evaluate (default: all METRICS)
        cv_threshold: Minimum CV = (max - min) / mean threshold (default 0.2)
        output_file: Path to save results CSV (optional)
    
    Returns:
        DataFrame with evaluation results
    """
    if datasets is None:
        datasets = DATASETS
    if metrics is None:
        metrics = METRICS
    
    all_results = []
    if local_run:
        scores_all = pd.read_csv(BENCHMARK_DIR / "all_scores.csv")
        scores_all = scores_all[METRICS + ['method', 'dataset']]
        scores_all.rename(columns={'method': 'model'}, inplace=True)
    
    # Pre-extract negative and positive control scores per dataset per metric
    nc_scores = {}  # nc_scores[dataset][metric] = nc_score
    pc_scores = {}  # pc_scores[dataset][metric] = pc_score
    if local_run:
        nc_df = scores_all[scores_all['model'] == 'negative_control'].copy()
        for _, row in nc_df.iterrows():
            ds = row['dataset']
            nc_scores.setdefault(ds, {})
            for metric in metrics:
                if metric in row and not pd.isna(row[metric]):
                    nc_scores[ds][metric] = row[metric]
        pc_df = scores_all[scores_all['model'] == 'positive_control'].copy()
        for _, row in pc_df.iterrows():
            ds = row['dataset']
            pc_scores.setdefault(ds, {})
            for metric in metrics:
                if metric in row and not pd.isna(row[metric]):
                    pc_scores[ds][metric] = row[metric]
    
    # print(f"Evaluating {len(metrics)} metrics across {len(datasets)} datasets")
    # print(f"CV threshold: {cv_threshold}")
    # print(f"\nMetric thresholds:")
    # for metric in metrics:
    #     threshold = METRIC_THRESHOLDS.get(metric, 0.1)
    #     print(f"  {metric}: {threshold}")

    
    for dataset in datasets:
        # print(f"\n{'='*60}")
        # print(f"Processing dataset: {dataset}")
        # print(f"{'='*60}")
        
        # Load scores
        if local_run:
            scores_df = scores_all[scores_all['dataset'] == dataset].drop(columns=['dataset']).set_index('model')
            scores_df = scores_df.loc[:, ~scores_df.isnull().all()]  # Drop columns with all NaNs
        else:
            scores_df = load_scores_from_yaml(dataset)

        
        if scores_df is None:
            # print(f"  No scores found for {dataset}")
            for metric in metrics:
                all_results.append({
                    'dataset': dataset,
                    'metric': metric,
                    'present': False,
                    'keep': False,
                    'reason': 'No scores file found'
                })
            continue
        

        # Evaluate each metric
        for metric in metrics:
            nc_score = nc_scores.get(dataset, {}).get(metric, None)
            pc_score = pc_scores.get(dataset, {}).get(metric, None)
            result = evaluate_metric_for_dataset(scores_df, metric, cv_threshold, nc_score=nc_score, pc_score=pc_score)
            result['dataset'] = dataset
            all_results.append(result)
            
            # Print summary
            if result['present']:
                status = "✓ KEEP" if result.get('keep', False) else "✗ DROP"
                # print(f"  {status} {metric:20s} - {result['reason']}")
            else:
                # print(f"  - SKIP {metric:20s} - {result['reason']}")
                pass
    
    # Convert to DataFrame
    results_df = pd.DataFrame(all_results)
    
    # Reorder columns for better readability
    cols_order = ['dataset', 'metric', 'keep', 'present', 'reason',
                  'n_methods', 'mean', 'std', 'cv', 'max', 'threshold', 'global_threshold', 'pc_score', 'nc_score', 'cv_threshold']
    cols_order = [c for c in cols_order if c in results_df.columns]
    results_df = results_df[cols_order]
    
    # Save if output file specified
    if output_file:
        output_path = Path(output_file)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        results_df.to_csv(output_path, index=False)
        print(f"\n{'='*60}")
        print(f"Results saved to: {output_path}")
        print(f"{'='*60}")
        
        # Also save a summary of metrics that passed per dataset for easy consumption
        # This will be used by create_overview_figure.py
        kept_metrics = results_df[results_df['keep'] == True].groupby('dataset')['metric'].apply(list).to_dict()
        summary_output = RESULTS_DIR / 'exp_analysis' / 'metrics_kept_per_dataset.yaml'
        import yaml
        with open(summary_output, 'w') as f:
            yaml.dump(kept_metrics, f)
        print(f"Metrics summary saved to: {summary_output}")
    
    # Print summary statistics
    # print(f"\n{'='*60}")
    # print("SUMMARY")
    # print(f"{'='*60}")
    
    # Overall summary
    total_combinations = len(datasets) * len(metrics)
    present_count = results_df['present'].sum()
    keep_count = results_df['keep'].sum()
    
    # print(f"\nTotal dataset-metric combinations: {total_combinations}")
    # print(f"Metrics computed: {present_count} ({100*present_count/total_combinations:.1f}%)")
    # print(f"Metrics to keep: {keep_count} ({100*keep_count/total_combinations:.1f}%)")
    
    # Per dataset summary
    # print(f"\nPer-dataset summary:")
    dataset_present = results_df[results_df['present']].copy()
    # Ensure 'keep' is treated as integer (True=1, False=0)
    dataset_present['keep'] = dataset_present['keep'].astype(int)
    
    dataset_summary = dataset_present.groupby('dataset')['keep'].agg(['sum', 'count'])
    dataset_summary.columns = ['Keep', 'Total']
    # dataset_summary['Drop'] = dataset_summary['Total'] - dataset_summary['Keep']
    dataset_summary['Keep %'] = 100 * dataset_summary['Keep'] / dataset_summary['Total']
    dataset_summary = dataset_summary[['Total', 'Keep', 'Keep %']]
    dataset_summary.sort_values(by='Keep %', ascending=False, inplace=True)
    print(dataset_summary.to_string())
    
    # Per metric summary
    # print(f"\nPer-metric summary:")
    metric_present = results_df[results_df['present']].copy()
    # Ensure 'keep' is treated as integer (True=1, False=0)
    metric_present['keep'] = metric_present['keep'].astype(int)
    metric_summary = metric_present.groupby('metric')['keep'].agg(['sum', 'count'])
    metric_summary.columns = ['Keep', 'Total']
    
    # Consolidate Keep/Total into one column
    metric_summary['Keep/Total'] = metric_summary.apply(
        lambda row: f"{int(row['Keep'])}/{int(row['Total'])}", axis=1
    )
    metric_summary['Keep %'] = 100 * metric_summary['Keep'] / metric_summary['Total']
    
    # Add min/max metric values and CV across datasets
    metric_stats = metric_present.groupby('metric').agg({
        'max': ['min', 'max'],
        'cv': ['min', 'max']
    })

    # Consolidate min/max into single columns with rounding
    metric_summary['Value (min/max)'] = metric_stats.apply(
        lambda row: f"{row[('max', 'min')]:.1f}/{row[('max', 'max')]:.1f}", axis=1
    )
    metric_summary['Variability (min/max)'] = metric_stats.apply(
        lambda row: f"{row[('cv', 'min')]:.1f}/{row[('cv', 'max')]:.1f}", axis=1
    )

    # Add threshold column: global floor per metric
    global_thresholds = metric_present.groupby('metric')['global_threshold'].first()
    nc_stats = metric_present.groupby('metric')['nc_score'].agg(['min', 'max'])
    metric_summary['Threshold'] = global_thresholds.map(lambda g: f"{g:.2f}")
    def format_nc(row):
        nc_min, nc_max = row['min'], row['max']
        if pd.isna(nc_min):
            return ""
        elif round(nc_min, 2) == round(nc_max, 2):
            return f"{nc_min:.2f}"
        else:
            return f"{nc_min:.2f}/{nc_max:.2f}"
    metric_summary['NC'] = nc_stats.apply(format_nc, axis=1)
    
    # Get valid datasets (where metric was actually computed AND kept)
    metric_kept = metric_present[metric_present['keep'] == 1].copy()
    valid_datasets_dict = metric_kept.groupby('metric')['dataset'].apply(
        lambda x: set(sorted(x.unique()))
    ).to_dict()
    
    # Add additional metadata from METRICS_DATASETS
    from src.utils.config import METRICS_DATASETS
  
    
    # Create metadata for each metric group
    metric_group_metadata = {
        'r_precision': {
            "Summary": "Evaluates a GRN by the ability of TFs to predict target gene expression. It only evaluates top regulators.",
            "Applicable Datasets": METRICS_DATASETS.get('regression', []),
            "Dataset Type Required": surrogate_names.get('bulk', 'bulk'),
            "Required Inputs": ", ".join([surrogate_names.get(i, i) for i in ['prediction', 'evaluation_data', 'tf_all', 'regulators_consensus']])
        },
        'r_recall': {
            "Summary": "Evaluates a GRN by the ability of TFs to predict target gene expression. It evaluates broader set of regulators.",
            "Applicable Datasets": METRICS_DATASETS.get('regression', []),
            "Dataset Type Required": surrogate_names.get('bulk', 'bulk'),
            "Required Inputs": ", ".join([surrogate_names.get(i, i) for i in ['prediction', 'evaluation_data', 'tf_all', 'regulators_consensus']])
        },
        'ws_precision': {
            "Summary": "Evaluates a GRN by the shift in gene expression of target gene in response to TF perturbation. It evaluates top TF-edge interactions.",
            "Applicable Datasets": METRICS_DATASETS.get('ws_distance', []),
            "Dataset Type Required": surrogate_names.get('sc', 'sc'),
            "Required Inputs": ", ".join([surrogate_names.get(i, i) for i in ['prediction', 'ws_consensus', 'ws_distance_background']])
        },
        'ws_recall': {
            "Summary": "Evaluates a GRN by the shift in gene expression of target gene in response to TF perturbation. It evaluates broader set of TF-edge interactions.",
            "Applicable Datasets": METRICS_DATASETS.get('ws_distance', []),
            "Dataset Type Required": surrogate_names.get('sc', 'sc'),
            "Required Inputs": ", ".join([surrogate_names.get(i, i) for i in ['prediction', 'ws_consensus', 'ws_distance_background']])
        },
        'sem': {
            "Summary": "Structural Equation Modeling for GRN evaluation.",
            "Applicable Datasets": METRICS_DATASETS.get('sem', []),
            "Dataset Type Required": surrogate_names.get('bulk', 'bulk'),
            "Required Inputs": ", ".join([surrogate_names.get(i, i) for i in ['prediction', 'evaluation_data', 'tf_all']])
        },
        'vc': {
            "Summary": "Evaluate GRNs by their ability in predicting gene expression through neural networks-based virtual cell model.",
            "Applicable Datasets": METRICS_DATASETS.get('vc', []),
            "Dataset Type Required": surrogate_names.get('bulk', 'bulk'),
            "Required Inputs": ", ".join([surrogate_names.get(i, i) for i in ['prediction', 'evaluation_data']])
        },
        't_rec_precision': {
            "Summary": "Measures ability to recover TFs using differential expression and GRN. It only evaluates available TFs in the GRN.",
            "Applicable Datasets": METRICS_DATASETS.get('tf_recovery', []),
            "Dataset Type Required": surrogate_names.get('de', 'de'),
            "Required Inputs": ", ".join([surrogate_names.get(i, i) for i in ['prediction', 'evaluation_data_de', 'tf_all']])
        },
        't_rec_recall': {
            "Summary": "Measures ability to recover TFs using differential expression and GRN. It considers all TFs in the evaluation data.",
            "Applicable Datasets": METRICS_DATASETS.get('tf_recovery', []),
            "Dataset Type Required": surrogate_names.get('de', 'de'),
            "Required Inputs": ", ".join([surrogate_names.get(i, i) for i in ['prediction', 'evaluation_data_de', 'tf_all']])
        },
        'replicate_consistency': {
            "Summary": "Measures consistency of GRN predictions across biological replicates/groups.",
            "Applicable Datasets": METRICS_DATASETS.get('replicate_consistency', []),
            "Dataset Type Required": surrogate_names.get('sc', 'sc'),
            "Required Inputs": ", ".join([surrogate_names.get(i, i) for i in ['prediction', 'evaluation_data']])
        },
        'tfb_f1': {
            "Summary": "Evaluates GRNs by their ability to match Chip-seq TF binding data",
            "Applicable Datasets": METRICS_DATASETS.get('tf_binding', []),
            "Dataset Type Required": surrogate_names.get('bulk', 'bulk'),
            "Required Inputs": ", ".join([surrogate_names.get(i, i) for i in ['prediction', 'evaluation_data', 'Ground Truth']])
        },
        'gs_f1': {
            "Summary": "Evaluates GRNs by their ability to recover known gene sets",
            "Applicable Datasets": METRICS_DATASETS.get('gs_recovery', []),
            "Dataset Type Required": surrogate_names.get('bulk', 'bulk'),
            "Required Inputs": ", ".join([surrogate_names.get(i, i) for i in ['prediction', 'evaluation_data', 'Gene Sets']])
        }
    }
    
    # Create metadata DataFrame with proper index matching
    metadata_rows = []
    for metric_id in metric_summary.index:
        # Check if metric_id exists directly in metadata
        if metric_id in metric_group_metadata:
            metadata = metric_group_metadata[metric_id].copy()
            
            # Show ALL applicable datasets, but add star only to those that passed threshold
            applicable_ds = set(metadata['Applicable Datasets'])
            valid_ds = valid_datasets_dict.get(metric_id, set())
            
            # Build list with all applicable datasets, marking valid ones with star
            all_datasets_list = []
            for ds in sorted(applicable_ds):
                ds_name = surrogate_names.get(ds, ds)
                if ds in valid_ds:
                    all_datasets_list.append(f"{ds_name}*")
                else:
                    all_datasets_list.append(ds_name)
            
            metadata['Datasets'] = ", ".join(all_datasets_list) if all_datasets_list else ""
            del metadata['Applicable Datasets']
            
            metadata_rows.append(metadata)
        else:
            # Fill with empty values if no metadata found
            valid_ds = valid_datasets_dict.get(metric_id, set())
            datasets_str = ", ".join([f"{surrogate_names.get(ds, ds)}*" for ds in sorted(valid_ds)])
            
            metadata_rows.append({
                "Summary": "",
                "Datasets": datasets_str,
                "Dataset Type Required": "",
                "Required Inputs": ""
            })
    
    metadata_df = pd.DataFrame(metadata_rows, index=metric_summary.index)
    
    # Add metric name as a column for display with surrogate names
    metric_summary_display = metric_summary.copy()
    metric_summary_display.insert(0, 'Metric', [
        surrogate_names.get(m, m)
        for m in metric_summary_display.index
    ])
    
    # Select and order final columns: Summary, Inputs, Datasets, then the rest
    final_columns = ['Metric', 'Summary', 'Required Inputs', 'Datasets', 'Threshold', 'Keep/Total', 'Keep %', 'Value (min/max)', 'Variability (min/max)']
    
    # Build the final dataframe
    metric_summary_final = pd.DataFrame()
    metric_summary_final['Metric'] = metric_summary_display['Metric']
    metric_summary_final['Summary'] = metadata_df['Summary'].values
    metric_summary_final['Inputs'] = metadata_df['Required Inputs'].values
    metric_summary_final['Datasets'] = metadata_df['Datasets'].values
    metric_summary_final['Keep/Total'] = metric_summary_display['Keep/Total'].values
    metric_summary_final['Keep %'] = metric_summary_display['Keep %'].values
    metric_summary_final['Value\n(min/max)'] = metric_summary_display['Value (min/max)'].values
    metric_summary_final['Variability\n(min/max)'] = metric_summary_display['Variability (min/max)'].values
    metric_summary_final['Threshold'] = metric_summary_display['Threshold'].values
    metric_summary_final['NC'] = metric_summary_display['NC'].values
    
    # Sort by Keep %
    metric_summary_final = metric_summary_final.sort_values(by='Keep %', ascending=False)
    
    # Round Keep % to 1 decimal
    metric_summary_final['Keep %'] = metric_summary_final['Keep %'].round(1)

    metric_summary_final = metric_summary_final.drop(columns=['Keep %'])
    
    # Word wrap long text fields for display
    def wrap_text(text, width=40):
        """Wrap text to specified width, preserving full content"""
        if not isinstance(text, str) or len(text) <= width:
            return text
        # Split into lines of max width
        lines = []
        while len(text) > width:
            # Try to break at comma or space
            break_at = text.rfind(',', 0, width)
            if break_at == -1:
                break_at = text.rfind(' ', 0, width)
            if break_at == -1:
                break_at = width
            else:
                break_at += 1  # Include comma/space
            lines.append(text[:break_at].strip())
            text = text[break_at:].strip()
        lines.append(text)
        return '\n'.join(lines)
    
    max_col_width = 40
    for col in ['Summary', 'Inputs', 'Datasets']:
        if col in metric_summary_final.columns:
            metric_summary_final[col] = metric_summary_final[col].apply(
                lambda x: wrap_text(x, max_col_width)
            )
    
    # Print with tabulate for better formatting
    try:
        from tabulate import tabulate
        print("\n" + tabulate(metric_summary_final, headers='keys', tablefmt='grid', showindex=False))
    except ImportError:
        # Fallback to pandas styling if tabulate not available
        print("\n" + metric_summary_final.to_string(index=False))
    
    # Save table as figure
    if output_file:
        import matplotlib.pyplot as plt
        import matplotlib
        matplotlib.use('Agg')  # Use non-interactive backend
        
        # Manual figure and table sizing - adjust these values as needed
        fig_width = 16
        fig_height = 8
    
        col_widths = [0.12, 0.22, 0.17, 0.20, 0.06, 0.08, 0.08, 0.07, 0.07]  # Must sum to ~1.0
        row_height = 0.08
        font_size = 7
        
        # Create figure
        fig, ax = plt.subplots(figsize=(fig_width, fig_height))
        ax.axis('tight')
        ax.axis('off')
        
        # Create table
        table = ax.table(
            cellText=metric_summary_final.values,
            colLabels=metric_summary_final.columns,
            cellLoc='left',
            loc='center',
            colWidths=col_widths
        )
        
        # Style the table
        table.auto_set_font_size(False)
        table.set_fontsize(font_size)
        
        # Set row heights — header row taller to fit two-line labels
        for i in range(len(metric_summary_final) + 1):  # +1 for header
            for j in range(len(metric_summary_final.columns)):
                table[(i, j)].set_height(row_height * 1.6 if i == 0 else row_height)
        
        # Style header
        for i in range(len(metric_summary_final.columns)):
            cell = table[(0, i)]
            cell.set_facecolor('#4CAF50')
            cell.set_text_props(weight='bold', color='white')
        
        # Alternate row colors
        for i in range(1, len(metric_summary_final) + 1):
            for j in range(len(metric_summary_final.columns)):
                cell = table[(i, j)]
                if i % 2 == 0:
                    cell.set_facecolor('#f0f0f0')
        
        # Save figure
        output_path_fig = str(Path(output_file).with_suffix('.png'))
        plt.savefig(output_path_fig, dpi=300, bbox_inches='tight', facecolor='white')
        print(f"Table figure saved to: {output_path_fig}")
        
        output_path_fig = DOCS_IMAGES_DIR / 'metric_quality_evaluation.png'
        print(f"Also saving copy to: {output_path_fig}")
        plt.savefig(output_path_fig, dpi=300, bbox_inches='tight', facecolor='white')
    
    # Recommended metrics per dataset
    if False:
        print(f"\n{'='*60}")

        for dataset in datasets:
            kept_metrics = results_df[(results_df['dataset'] == dataset) & 
                                    (results_df['keep'] == True)]['metric'].tolist()
            if kept_metrics:
                print(f"\n{dataset}:")
                print(f"  {kept_metrics}")
            else:
                print(f"\n{dataset}: No metrics recommended")
    
    return results_df


def plot_metric_applicability():
    """Plot a heatmap of metric applicability per dataset from config (DATASETS_METRICS)."""
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches

    # Map config categories → individual METRICS
    category_to_metrics = {
        'regression':           ['r_precision', 'r_recall'],
        'ws_distance':          ['ws_precision', 'ws_recall'],
        'tf_recovery':          ['t_rec_precision', 't_rec_recall'],
        'tf_binding':           ['tfb_f1'],
        'sem':                  ['sem'],
        'gs_recovery':          ['gs_f1'],
        'vc':                   ['vc'],
        'replicate_consistency':['replicate_consistency'],
    }

    # Build binary applicability matrix (metrics × datasets)
    dataset_order = list(DATASETS)
    metric_order = METRICS  # as defined in config

    from src.utils.config import DATASETS_METRICS
    applicable = pd.DataFrame(False, index=metric_order, columns=dataset_order)
    for dataset, categories in DATASETS_METRICS.items():
        if dataset not in dataset_order:
            continue
        for cat in categories:
            for m in category_to_metrics.get(cat, []):
                if m in applicable.index:
                    applicable.loc[m, dataset] = True

    # Sort metrics by descending applicability count
    applicable = applicable.loc[applicable.sum(axis=1).sort_values(ascending=True).index]

    # Display names
    applicable.index   = [surrogate_names.get(m, m) for m in applicable.index]
    applicable.columns = [surrogate_names.get(d, d) for d in applicable.columns]

    # Load empirically-passed metrics from metrics_kept_per_dataset.yaml
    results_dir = str(RESULTS_DIR)
    metrics_kept_file = Path(results_dir) / 'exp_analysis' / 'metrics_kept_per_dataset.yaml'
    empirical = pd.DataFrame(False, index=metric_order, columns=dataset_order)
    if metrics_kept_file.exists():
        with open(metrics_kept_file) as f:
            metrics_kept = yaml.safe_load(f)
        for dataset, mlist in metrics_kept.items():
            if dataset not in dataset_order:
                continue
            for m in mlist:
                if m in empirical.index:
                    empirical.loc[m, dataset] = True

    # Apply same metric sort order as applicable
    empirical = empirical.loc[applicable.index.map(
        {v: k for k, v in {m: surrogate_names.get(m, m) for m in metric_order}.items()}
    )]
    empirical.index   = applicable.index
    empirical.columns = applicable.columns

    n_metrics, n_datasets = applicable.shape
    cell = 0.15
    left_margin = 1.1
    fig_w = n_datasets * cell + left_margin + 0.6
    fig_h = n_metrics * cell + 0.4
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    ax.set_xlim(-0.5, n_datasets - 0.5)
    ax.set_ylim(-0.5, n_metrics - 0.5)
    ax.set_facecolor('#f8f8f8')

    # Grid lines
    for x in range(n_datasets):
        ax.axvline(x - 0.5, color='white', linewidth=1.2)
    for y in range(n_metrics):
        ax.axhline(y - 0.5, color='white', linewidth=1.2)

    # Fill cells: green checkmark (config applicable) + red star (empirically passed)
    for yi, metric in enumerate(applicable.index):
        for xi, dataset in enumerate(applicable.columns):
            val = applicable.loc[metric, dataset]
            passed = empirical.loc[metric, dataset]
            color = '#d4edda' if val else '#efefef'
            ax.add_patch(mpatches.FancyBboxPatch(
                (xi - 0.44, yi - 0.41), 0.88, 0.82,
                boxstyle='round,pad=0.04', facecolor=color, edgecolor='white', linewidth=0.5
            ))
            if val and not passed:
                ax.text(xi, yi, '✓', ha='center', va='center',
                        fontsize=5, color='#2e7d32', fontweight='bold')
            elif val and passed:
                ax.text(xi - 0.18, yi, '✓', ha='center', va='center',
                        fontsize=5, color='#2e7d32', fontweight='bold')
                ax.text(xi + 0.18, yi, '★', ha='center', va='center',
                        fontsize=4, color='#c0392b', fontweight='bold')

    ax.set_xticks(range(n_datasets))
    ax.set_xticklabels(applicable.columns, rotation=35, ha='left', fontsize=5)
    ax.set_yticks(range(n_metrics))
    ax.set_yticklabels(applicable.index, fontsize=5, ha='right')
    ax.xaxis.set_ticks_position('top')
    ax.xaxis.set_label_position('top')
    ax.tick_params(axis='both', length=0, pad=1)
    ax.set_xlabel('Dataset', fontsize=6, labelpad=4)
    ax.set_ylabel('Metric', fontsize=6, labelpad=4)

    # Fix left margin so long y labels don't get clipped
    fig.subplots_adjust(left=left_margin / fig_w)
    plt.tight_layout()

    out =  f'{results_dir}/exp_analysis/metric_applicability.png'
    plt.savefig(out, dpi=400, bbox_inches='tight', facecolor='white')
    print(f"Applicability plot saved to: {out}")

    docs_out = DOCS_IMAGES_DIR / 'metric_applicability.png'
    plt.savefig(docs_out, dpi=400, bbox_inches='tight', facecolor='white')
    print(f"Also saved to: {docs_out}")
    plt.close()


def main():
    # Set default output path if not specified
    if args.output is None:
        results_dir = str(RESULTS_DIR)
        args.output = f"{results_dir}/exp_analysis/metric_quality_evaluation.csv"
    
    # Run evaluation
    results = evaluate_all_datasets(
        datasets=args.datasets,
        metrics=args.metrics,
        cv_threshold=args.cv_threshold,
        output_file=args.output
    )

    plot_metric_applicability()

    return results


if __name__ == '__main__':
    main()
