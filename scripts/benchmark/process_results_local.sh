#!/usr/bin/env bash
# Post-metric processing pipeline: aggregate scores, combine results, and plot.
# Run this after GRN evaluation (metric calculation) is complete.
#
# Usage:
#   bash scripts/benchmark/process_results_local.sh

# conda activate genernbi
set -euo pipefail

echo "=========================================="
echo "GRN Post-Metric Processing — LOCAL"
echo "=========================================="
RESULTS_DIR="resources/results/"
BENCHMARK_DIR="$RESULTS_DIR/benchmark"

echo ""
echo "── Step 1: Aggregating h5ad scores → all_scores.csv ──"
python "scripts/benchmark/aggregate_local_scores.py" \
    --results_file "${BENCHMARK_DIR}/all_scores.csv"

echo ""
echo "── Step 2: Combining results → all_new/ ──"
python scripts/benchmark/combine_results.py

echo ""
echo "── Step 3: Evaluating metric applicability per dataset ──"
python "scripts/benchmark/evaluate_metric_applicability.py" \
    --output "${RESULTS_DIR}/exp_analysis/metric_quality_evaluation.csv"

echo ""
echo "── Step 4: Plotting dataset heatmaps ──"
python "scripts/benchmark/plot_raw_scores.py" \
    --output_dir_docs "docs/source/images"

echo ""
echo "── Step 5: Generating summary figure ──"
python "scripts/benchmark/create_overview_figure.py"

echo ""
echo "Done."
