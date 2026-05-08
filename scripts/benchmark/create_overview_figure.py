"""
Create overview figure from combined results.
Processes trace.csv and score files, then calls the R script to generate the summary figure.
"""

import subprocess
import pandas as pd
import numpy as np
import yaml
from pathlib import Path
import argparse

from src.utils.util import compute_overall_scores
from src.utils.config import DATASETS, METHODS, METRICS, DATASET_INFO, surrogate_names


def _rank_on_datasets(scores_all, dataset_list, final_metrics):
    """Median percentile rank across all (dataset, metric) pairs, normalized to [0,1]."""
    records = []
    sub = scores_all[scores_all['dataset'].isin(dataset_list)]
    for dataset, grp in sub.groupby('dataset'):
        grp = grp.set_index('model')
        for metric in final_metrics:
            if metric not in grp.columns:
                continue
            col = grp[metric].dropna()
            if col.empty:
                continue
            for method, r in col.rank(ascending=True, pct=True).items():
                records.append({'model': method, 'rank': r})
    if not records:
        return pd.Series(dtype=float)
    med = pd.DataFrame(records).groupby('model')['rank'].median().sort_values(ascending=False)
    r_min, r_max = med.min(), med.max()
    if r_max > r_min:
        med = (med - r_min) / (r_max - r_min)
    return med


def print_rankings(scores_all, final_metrics, modality_groups):
    multi_ds = modality_groups.get('Multiomics', [])
    trans_ds = modality_groups.get('Transcriptomics', [])

    multi_rank = _rank_on_datasets(scores_all, multi_ds, final_metrics)
    trans_rank = _rank_on_datasets(scores_all, trans_ds, final_metrics)

    all_methods = set(multi_rank.index) | set(trans_rank.index)
    overall = {}
    for m in all_methods:
        parts = [v for v in [multi_rank.get(m), trans_rank.get(m)] if pd.notna(v)]
        overall[m] = np.mean(parts) if parts else np.nan
    overall = pd.Series(overall).dropna().sort_values(ascending=False)
    r_min, r_max = overall.min(), overall.max()
    if r_max > r_min:
        overall = (overall - r_min) / (r_max - r_min)

    def _print(title, ranked):
        print(f"\n  {title}")
        print(f"  {'-' * 44}")
        for i, (m, v) in enumerate(ranked.items(), 1):
            name = surrogate_names.get(m, m)
            print(f"    {i:2d}. {name:20s}  {v:.3f}")

    print("\n" + "=" * 80)
    print("Rankings")
    print("=" * 80)
    _print(f"Overall (modality-balanced, {len(multi_ds)} multiomics + {len(trans_ds)} transcriptomics datasets)", overall)
    _print(f"Multiomics datasets {multi_ds}", multi_rank)
    _print(f"Transcriptomics datasets ({len(trans_ds)} datasets)", trans_rank)

    # Per-metric rankings
    all_ds = list(scores_all['dataset'].unique())
    print(f"\n  {'Per-metric rankings (across all applicable datasets)'}")
    print(f"  {'-' * 44}")
    for metric in final_metrics:
        records = []
        for dataset, grp in scores_all.groupby('dataset'):
            grp = grp.set_index('model')
            if metric not in grp.columns:
                continue
            col = grp[metric].dropna()
            if col.empty:
                continue
            for method, r in col.rank(ascending=True, pct=True).items():
                records.append({'model': method, 'rank': r})
        if not records:
            continue
        med = pd.DataFrame(records).groupby('model')['rank'].median().sort_values(ascending=False)
        r_min, r_max = med.min(), med.max()
        if r_max > r_min:
            med = (med - r_min) / (r_max - r_min)
        metric_label = surrogate_names.get(metric, metric)
        n_datasets = scores_all.groupby('dataset').apply(
            lambda g: metric in g.columns and g[metric].notna().any()
        ).sum()
        print(f"\n    {metric_label} ({n_datasets} datasets)")
        for i, (m, v) in enumerate(med.items(), 1):
            name = surrogate_names.get(m, m)
            print(f"      {i:2d}. {name:20s}  {v:.3f}")

RESULTS_DIR = Path("resources/results")
BENCHMARK_DIR = Path("resources/results/benchmark")
DOCS_IMAGES_DIR = Path("docs/source/images")
GENERNBI_SUPP_DIR = Path("../genernbi_supp")


def process_scores_from_yaml(score_file):
    with open(score_file, 'r') as f:
        scores_data = yaml.safe_load(f)
    rows = []
    for entry in scores_data:
        if entry is None or 'missing' in str(entry):
            continue
        dataset_id = entry.get('dataset_id', '')
        method_id = entry.get('method_id', '')
        for metric, value in zip(entry.get('metric_ids', []), entry.get('metric_values', [])):
            if value != "None":
                try:
                    rows.append({'dataset': dataset_id, 'model': method_id,
                                 'metric': metric, 'value': float(value)})
                except (ValueError, TypeError):
                    pass
    df = pd.DataFrame(rows)
    return df.pivot_table(index=['dataset', 'model'], columns='metric', values='value').reset_index()


def process_scores_from_csv(score_file):
    scores_all = pd.read_csv(score_file)
    if 'method' in scores_all.columns:
        scores_all = scores_all.rename(columns={'method': 'model'})
    return scores_all


def parse_duration(duration_str):
    if pd.isna(duration_str):
        return 0
    total_seconds = 0
    for part in str(duration_str).strip().split():
        if 'h' in part:
            total_seconds += float(part.replace('h', '')) * 3600
        elif 'm' in part:
            total_seconds += float(part.replace('m', '')) * 60
        elif 's' in part:
            total_seconds += float(part.replace('s', ''))
    return total_seconds / 3600


def convert_to_gb(value):
    if pd.isna(value):
        return 0
    value = str(value).strip()
    unit_to_bytes = {"KB": 1024, "MB": 1024**2, "GB": 1024**3, "B": 1}
    try:
        parts = value.split()
        if len(parts) == 2:
            num, unit = parts
            if unit.upper() in unit_to_bytes:
                return float(num) * unit_to_bytes[unit.upper()] / (1024**3)
        return float(value)
    except Exception:
        return 0


def process_trace_to_csv(trace_file):
    trace_df = pd.read_csv(trace_file, sep='\t')
    print(f"   Total trace entries: {len(trace_df)}")
    resource_data = []
    for method in METHODS:
        pattern = f'run_grn_inference.*:{method}:processWf:{method}_process \\(op\\.{method}\\)'
        matches = trace_df[trace_df['name'].str.contains(pattern, regex=True, na=False)]
        if len(matches) == 0:
            alt = f'{method}_process \\(op\\.{method}\\)$'
            matches = trace_df[trace_df['name'].str.contains(alt, regex=True, na=False)]
        if len(matches) == 0:
            print(f"   Warning: No inference trace found for {method}")
            continue
        if len(matches) > 1:
            print(f"   Warning: Multiple traces found for {method}, using first")
        row = matches.iloc[0]
        duration = parse_duration(row['duration'])
        peak_rss = convert_to_gb(row['peak_rss'])
        cpu_val = row['%cpu']
        cpu_usage = float(str(cpu_val).replace('%', '')) if pd.notna(cpu_val) else 0
        resource_data.append({'method': method, 'Duration (hour)': duration,
                               'Peak memory (GB)': peak_rss, 'CPU Usage (%)': cpu_usage})
        print(f"   {method}: duration={duration:.2f}h, memory={peak_rss:.1f}GB, cpu={cpu_usage:.1f}%")
    df_resources = pd.DataFrame(resource_data).set_index('method')
    if df_resources.index.duplicated().any():
        df_resources = df_resources.groupby(df_resources.index).max()
    print(f"\n   Processed resources for {len(df_resources)} methods")
    return df_resources


def main(local_run=True, methods=None, datasets=None):
    combined_dir = BENCHMARK_DIR / 'all_new'
    trace_file = combined_dir / 'trace.csv'
    score_file = BENCHMARK_DIR / 'all_scores.csv' if local_run else combined_dir / 'score_uns.yaml'

    print("=" * 80)
    print("Creating Overview Figure")
    print(f"Mode: {'LOCAL RUN' if local_run else 'AWS'} | Score file: {score_file}")
    print("=" * 80)

    print("\n1. Processing trace data (op dataset only)...")
    df_res = process_trace_to_csv(trace_file)

    print("\n2. Processing scores...")
    scores_all = process_scores_from_csv(score_file) if local_run else process_scores_from_yaml(score_file)
    if methods is not None:
        scores_all = scores_all[scores_all['model'].isin(methods)]
    if datasets is not None:
        scores_all = scores_all[scores_all['dataset'].isin(datasets)]
    print(f"   {len(scores_all)} method-dataset combinations")

    print("\n2.5. Loading metrics that passed applicability criteria...")
    metrics_kept_file = RESULTS_DIR / 'exp_analysis' / 'metrics_kept_per_dataset.yaml'
    if not metrics_kept_file.exists():
        raise FileNotFoundError(f"Metrics kept file not found: {metrics_kept_file}")
    with open(metrics_kept_file) as f:
        metrics_kept_per_dataset = yaml.safe_load(f)
    all_kept_metrics = {m for mlist in metrics_kept_per_dataset.values() for m in mlist}
    print(f"   {len(all_kept_metrics)} metrics passed")

    print("\n3. Creating summary dataframe...")
    final_metrics = [m for m in METRICS if m in all_kept_metrics]
    modality_groups = {}
    for dataset in DATASETS:
        if dataset not in scores_all['dataset'].unique():
            continue
        modality = DATASET_INFO.get(dataset, {}).get('Modality', 'Transcriptomics')
        modality_groups.setdefault(modality, []).append(dataset)

    df_scores = compute_overall_scores(scores_all, final_metrics, DATASETS, modality_groups=modality_groups)
    methods_in_scores = df_scores.index.tolist()

    print_rankings(scores_all, final_metrics, modality_groups)

    df_summary = pd.concat([df_scores, df_res], axis=1)
    df_summary = df_summary[df_summary.index.isin(methods_in_scores)].fillna(0)
    df_summary.index.name = 'method_name'
    df_summary = df_summary.reset_index().sort_values(by='overall_score', ascending=False)

    df_summary['method_name'] = df_summary['method_name'].map(lambda x: surrogate_names.get(x, x))
    df_summary['User-friendly'] = df_summary['method_name'].map({
        'Scenic+': 1, 'GRNBoost2': 8, 'Positive Ctrl': 10, 'Pearson Corr.': 10,
        'Spearman Corr.': 10, 'CellOracle': 6, 'Portia': 9, 'scGLUE': 6,
        'Scenic': 7, 'FigR': 6, 'PPCOR': 7, 'Negative Ctrl': 10, 'GRaNIE': 6,
        'scPRINT': 5, 'Geneformer': 5, 'scGPT': 3,
    })
    df_summary['Complexity'] = df_summary['User-friendly'].max() - df_summary['User-friendly']
    df_summary.columns = [surrogate_names.get(col, col) for col in df_summary.columns]
    df_summary = df_summary.fillna(0)

    for col, label in [('Peak memory (GB)', 'memory'), ('Complexity', 'complexity'), ('Duration (hour)', 'duration')]:
        df_summary[f'{label}_log'] = np.log(df_summary[col] + 1)
        df_summary[f'{label}_log'] = np.max(df_summary[f'{label}_log']) - df_summary[f'{label}_log']
    df_summary['duration_str'] = df_summary['Duration (hour)'].round(1).astype(str)
    df_summary['memory_str'] = df_summary['Peak memory (GB)'].round(1).astype(str)

    summary_file = RESULTS_DIR / "benchmark" / "summary.tsv"
    summary_file.parent.mkdir(parents=True, exist_ok=True)
    df_summary.to_csv(summary_file, sep='\t', index=False)
    print(f"\n4. Saved summary to: {summary_file}")

    print("\n5. Calling R script to create figure...")
    r_script = GENERNBI_SUPP_DIR / "src" / "summary_figure.R"
    summary_figure = RESULTS_DIR / "benchmark" / "summary_figure"
    result = subprocess.run(f"Rscript {r_script} {summary_file} {summary_figure}",
                            shell=True, capture_output=True, text=True)
    if result.returncode == 0:
        print(f"\n✓ Figure saved to: {summary_figure}.pdf / .png")
        subprocess.run(f"cp {summary_figure}.png {DOCS_IMAGES_DIR / 'summary_figure.png'}", shell=True)
    else:
        print(f"\n✗ R script error:\n{result.stderr}\n{result.stdout}")
        return 1

    print("\n" + "=" * 80 + "\nDone!\n" + "=" * 80)
    return 0


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Create overview figure from combined results')
    parser.add_argument('--aws_run', dest='local_run', action='store_false')
    parser.set_defaults(local_run=True)
    parser.add_argument('--methods', nargs='+', default=None)
    parser.add_argument('--datasets', nargs='+', default=None)
    args = parser.parse_args()
    exit(main(local_run=args.local_run,
              methods=args.methods or METHODS,
              datasets=args.datasets or DATASETS))
