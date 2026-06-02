"""
Plot per-dataset raw score heatmaps from GRN benchmark results.
"""

import os
import argparse
import warnings
import yaml
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path

warnings.filterwarnings("ignore")
matplotlib.rcParams["font.family"] = "Arial"

from src.utils.util import plot_heatmap
from src.utils.config import DATASETS_METRICS, METHODS, METRICS, surrogate_names

RESULTS_DIR = Path("resources/results")
BENCHMARK_DIR = Path("resources/results/benchmark")
DOCS_IMAGES_DIR = Path("docs/source/images")

CATEGORY_TO_METRICS = {
    'regression':            ['r_precision', 'r_recall'],
    'ws_distance':           ['ws_precision', 'ws_recall'],
    'tf_recovery':           ['t_rec_precision', 't_rec_recall'],
    'tf_binding':            ['tfb_f1'],
    'sem':                   ['sem'],
    'gs_recovery':           ['gs_f1'],
    'vc':                    ['vc'],
    'replicate_consistency': ['replicate_consistency'],
}

def metrics_for_dataset(dataset: str) -> list:
    applicable = []
    for cat in DATASETS_METRICS.get(dataset, []):
        applicable.extend(CATEGORY_TO_METRICS.get(cat, []))
    return applicable


# ── Score loading ──────────────────────────────────────────────────────────────

def load_scores_from_csv(scores_path: str, metrics: list, methods: list) -> pd.DataFrame:
    df = pd.read_csv(scores_path)
    keep_cols = [c for c in metrics if c in df.columns] + ["method", "dataset"]
    df = df[keep_cols]
    df.rename(columns={"method": "model"}, inplace=True)
    return df[df["model"].isin(methods)]


def load_scores_from_yaml(score_file: str, metrics: list, methods: list) -> pd.DataFrame:
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
    scores_all = df.pivot_table(
        index=['dataset', 'model'], columns='metric', values='value'
    ).reset_index()
    scores_all.columns.name = ""
    return scores_all[scores_all["model"].isin(methods)]


# ── Plotting ───────────────────────────────────────────────────────────────────

def plot_raw_scores(scores_all, datasets, output_dirs, metrics):
    for out_dir in output_dirs:
        os.makedirs(out_dir, exist_ok=True)

    for dataset in datasets:
        dataset_metrics = [m for m in metrics_for_dataset(dataset) if m in metrics]

        scores = scores_all[scores_all["dataset"] == dataset].copy()
        scores = scores.loc[:, ~scores.isna().all()]
        if scores.empty:
            print(f"  No data for dataset {dataset} — skipping")
            continue

        scores = scores.set_index("model").drop(columns="dataset", errors="ignore")
        scores = scores[[c for c in dataset_metrics if c in scores.columns]]
        if scores.empty:
            print(f"  No metric columns for dataset {dataset} — skipping")
            continue

        scores.columns = scores.columns.map(lambda n: surrogate_names.get(n, n))
        scores.index = scores.index.map(lambda n: surrogate_names.get(n, n))
        scores = scores.astype(float)

        ranks = scores.rank(axis=0, ascending=False, method="min")
        scores["_rank"] = ranks.mean(axis=1, skipna=True)
        scores = scores.sort_values("_rank").drop(columns=["_rank"])

        n_rows, n_cols = scores.shape
        fig, ax = plt.subplots(1, 1, figsize=(max(n_cols * 0.6, 4), max(n_rows * 0.4, 3)))
        plot_heatmap(scores, ax=ax, cmap="viridis")
        ax.set_yticklabels(ax.get_yticklabels(), rotation=0, ha="right")
        ax.set_ylabel("")
        plt.suptitle(surrogate_names.get(dataset, dataset), y=1.01, fontsize=14, weight="bold")
        plt.tight_layout()

        for out_dir in output_dirs:
            fpath = os.path.join(out_dir, f"raw_scores_{dataset}.png")
            fig.savefig(fpath, dpi=150, transparent=True, bbox_inches="tight")
            print(f"  Saved: {fpath}", flush=True)

        plt.close(fig)


# ── Entry point ────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot per-dataset raw score heatmaps")
    parser.add_argument(
        "--aws_run", dest="local_run", action="store_false",
        help="AWS mode — read scores from all_new/score_uns.yaml instead of all_scores.csv",
    )
    parser.set_defaults(local_run=True)
    parser.add_argument("--scores_file", default=str(BENCHMARK_DIR / "all_scores.csv"))
    parser.add_argument("--output_dir_local", default=str(RESULTS_DIR / "figs" / "raw_scores"))
    parser.add_argument("--output_dir_docs", default=str(DOCS_IMAGES_DIR))
    args = parser.parse_args()

    print("=== Plotting raw scores ===")
    if args.local_run:
        print(f"Mode: LOCAL — reading from {args.scores_file}")
        scores_all = load_scores_from_csv(args.scores_file, METRICS, METHODS)
    else:
        yaml_file = str(RESULTS_DIR / "all_new" / "score_uns.yaml")
        print(f"Mode: AWS — reading from {yaml_file}")
        scores_all = load_scores_from_yaml(yaml_file, METRICS, METHODS)

    output_dirs = [args.output_dir_local, args.output_dir_docs]
    plot_raw_scores(scores_all, list(DATASETS_METRICS.keys()), output_dirs, METRICS)
    print("Done.")
