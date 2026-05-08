"""
Aggregate *_score.h5ad files from local GRN evaluation output into all_scores.csv.
Run this after evaluation is complete (--run_grn_evaluation).
"""

import shutil
import argparse
import numpy as np
import anndata as ad
import pandas as pd
from pathlib import Path

from src.utils.config import DATASETS_METRICS

RESULTS_DIR = Path("resources/results")
BENCHMARK_DIR = Path("resources/results/benchmark")
DOCS_DATA_DIR = Path("docs/source/data")


def aggregate_score_files(results_file: Path, known_datasets: list) -> pd.DataFrame:
    """
    Collect *_score.h5ad files, parse dataset/method from filename,
    and write a wide-format CSV to results_file.

    Filename convention: {dataset}__{method}_score.h5ad
    Score files live in: resources/results/{dataset}/scores/
    """
    score_files = sorted(BENCHMARK_DIR.glob("*/scores/*__*_score.h5ad"))
    print(f"Looking for score files under: {BENCHMARK_DIR}/*/scores/")
    print(f"Will save results to: {results_file}")

    all_results = []
    for score_file in score_files:
        stem = score_file.stem
        if stem.endswith("_score"):
            stem = stem[:-6]

        dataset, method = None, None
        for ds in known_datasets:
            if stem.startswith(ds + "__"):
                dataset = ds
                method = stem[len(ds) + 2:]
                break

        if not dataset or not method:
            print(f"  Warning: could not parse filename {score_file.name} — skipping")
            continue

        try:
            adata = ad.read_h5ad(score_file)
        except Exception as e:
            print(f"  Error reading {score_file.name}: {e}")
            continue

        metric_ids = adata.uns.get("metric_ids")
        metric_values = adata.uns.get("metric_values")
        if metric_ids is None or metric_values is None:
            print(f"  Warning: no metric data in {score_file.name}")
            continue

        for mid, mval in zip(metric_ids, metric_values):
            try:
                score = float(mval)
            except (ValueError, TypeError):
                score = np.nan
            all_results.append({"dataset": dataset, "method": method, "metric": mid, "score": score})

        print(f"  Processed: {dataset} / {method}  ({len(metric_ids)} metrics)")

    if not all_results:
        print("No results found.")
        return pd.DataFrame()

    df = pd.DataFrame(all_results)
    df_wide = (
        df.pivot_table(index=["dataset", "method"], columns="metric", values="score")
        .reset_index()
    )
    df_wide.columns.name = ""

    results_file.parent.mkdir(parents=True, exist_ok=True)
    df_wide.to_csv(results_file, index=False)

    DOCS_DATA_DIR.mkdir(parents=True, exist_ok=True)
    shutil.copy(results_file, DOCS_DATA_DIR / "all_scores.csv")
    print(f"Copied scores to docs: {DOCS_DATA_DIR / 'all_scores.csv'}")

    print(f"\nResults saved to: {results_file}")
    print(f"  Evaluations : {len(df_wide)}")
    print(f"  Datasets    : {df_wide['dataset'].nunique()}")
    print(f"  Methods     : {df_wide['method'].nunique()}")
    print(f"  Metrics     : {len([c for c in df_wide.columns if c not in ['dataset', 'method']])}")
    return df_wide


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Aggregate *_score.h5ad files from local evaluation into all_scores.csv"
    )
    parser.add_argument(
        "--results_file",
        type=Path,
        default=BENCHMARK_DIR / "all_scores.csv",
        help="Output CSV path",
    )
    args = parser.parse_args()

    aggregate_score_files(args.results_file, list(DATASETS_METRICS.keys()))
