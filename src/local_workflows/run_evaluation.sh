#!/bin/bash
# Local GRN Evaluation Script
# This script evaluates all GRN predictions for all datasets using local SLURM jobs
# 
# Usage:
#   bash src/local_workflows/run_evaluation.sh --output_dir <dir> --run_metrics true --process_results false
#   After all jobs complete:
#   bash src/local_workflows/run_evaluation.sh --output_dir <dir> --run_metrics false --process_results true

# Default parameters
RUN_METRICS=false
PROCESS_RESULTS=true
OUTPUT_DIR="output/evaluation"
LAYER="lognorm"
NUM_WORKERS=20


echo "=========================================="
echo "GRN Evaluation Configuration"
echo "=========================================="
echo "Output directory: $OUTPUT_DIR"
echo "Run metrics: $RUN_METRICS"
echo "Process results: $PROCESS_RESULTS"
echo "Number of workers: $NUM_WORKERS"
echo "=========================================="

# Create output directory
mkdir -p "$OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR/scores"
mkdir -p "$OUTPUT_DIR/logs"

# Generate and source dataset configuration
echo "Generating dataset configuration..."
python src/utils/config.py --output src/utils/dataset_config.env
source src/utils/dataset_config.env

# Get list of datasets from config
DATASETS=($(python -c "from src.utils.config import DATASET_GROUPS; print(' '.join(DATASET_GROUPS.keys()))"))

echo "Datasets to evaluate: ${DATASETS[@]}"

# Method names to check
GRN_METHODS=(
    "positive_control"
    "pearson_corr"
    "negative_control"
    "spearman_corr"
    "scglue"
    "scenicplus"
    "celloracle"
    "granie"
    "figr"
    "grnboost"
    "portia"
    "scenic"
    "scprint"
    "geneformer"
    "scgpt"
    "ppcor"
)

# Function to submit a metric evaluation job
submit_metric_job() {
    local dataset=$1
    local method=$2
    local prediction_file=$3
    
    local job_name="${dataset}_${method}"
    local score_file="${OUTPUT_DIR}/scores/${dataset}_${method}_score.h5ad"
    
    # Skip if score file already exists
    if [[ -f "$score_file" ]]; then
        echo "  Skipping ${job_name} - score file already exists"
        return
    fi
    
    echo "  Submitting job: ${job_name}"
    
    sbatch \
        --job-name="${job_name}" \
        --output="${OUTPUT_DIR}/logs/${job_name}_%j.out" \
        --error="${OUTPUT_DIR}/logs/${job_name}_%j.err" \
        src/metrics/all_metrics/run_local.sh \
        --dataset "${dataset}" \
        --prediction "${prediction_file}" \
        --score "${score_file}" \
        --num_workers "${NUM_WORKERS}" || echo "  [ERROR] Failed to submit job for ${job_name}"
}

# Function to process all results into a single CSV
process_all_results() {
    echo ""
    echo "=========================================="
    echo "Processing Results"
    echo "=========================================="
    
    local results_file="${OUTPUT_DIR}/all_scores.csv"
    
    # Create Python script to aggregate results
    python << 'EOF'
import sys
import os
import pandas as pd
import anndata as ad
from pathlib import Path

output_dir = sys.argv[1]
scores_dir = Path(output_dir) / "scores"
results_file = Path(output_dir) / "all_scores.csv"

print(f"Looking for score files in: {scores_dir}")

all_results = []

# Iterate through all score files
for score_file in scores_dir.glob("*_score.h5ad"):
    try:
        # Parse filename: {dataset}_{method}_score.h5ad
        filename = score_file.stem  # Remove .h5ad
        parts = filename.replace("_score", "").rsplit("_", 1)
        
        if len(parts) != 2:
            print(f"Warning: Could not parse filename {score_file.name}")
            continue
        
        dataset, method = parts
        
        # Load score file
        adata = ad.read_h5ad(score_file)
        
        # Extract metric scores from uns
        if 'metric_ids' in adata.uns and 'metric_values' in adata.uns:
            metric_ids = adata.uns['metric_ids']
            metric_values = adata.uns['metric_values']
            
            # Create row for each metric
            for metric_id, metric_value in zip(metric_ids, metric_values):
                all_results.append({
                    'dataset': dataset,
                    'method': method,
                    'metric': metric_id,
                    'score': metric_value
                })
            
            print(f"Processed: {dataset} - {method} ({len(metric_ids)} metrics)")
        else:
            print(f"Warning: No metric data in {score_file.name}")
    
    except Exception as e:
        print(f"Error processing {score_file.name}: {e}")

# Create DataFrame
if all_results:
    df = pd.DataFrame(all_results)
    
    # Pivot to wide format: rows=dataset+method, columns=metrics
    df_wide = df.pivot_table(
        index=['dataset', 'method'],
        columns='metric',
        values='score'
    ).reset_index()
    
    # Save results
    df_wide.to_csv(results_file, index=False)
    print(f"\n{'='*50}")
    print(f"Results saved to: {results_file}")
    print(f"Total evaluations: {len(df_wide)}")
    print(f"Datasets: {df_wide['dataset'].nunique()}")
    print(f"Methods: {df_wide['method'].nunique()}")
    print(f"Metrics: {len([c for c in df_wide.columns if c not in ['dataset', 'method']])}")
    print(f"{'='*50}")
else:
    print("No results found to process!")

EOF
    
    python -c "$(cat)" "$OUTPUT_DIR"
}

# Main execution
if [[ "$RUN_METRICS" == "true" ]]; then
    echo ""
    echo "=========================================="
    echo "Submitting Metric Evaluation Jobs"
    echo "=========================================="
    
    job_count=0
    
    for dataset in "${DATASETS[@]}"; do
        echo ""
        echo "Dataset: $dataset"
        echo "----------------------------------------"
        
        models_folder="resources/results/${dataset}/"
        # echo "Looking in: $models_folder"
        
        # Check each method for this dataset
        for method in "${GRN_METHODS[@]}"; do
            prediction_file="${models_folder}/${dataset}.${method}.${method}.prediction.h5ad"
            
            if [[ -f "$prediction_file" ]]; then
                submit_metric_job "$dataset" "$method" "$prediction_file"
                # echo "  Submitting job: ${dataset}_${method}"
                ((job_count++))
            else
                echo "  [NOT FOUND] ${prediction_file}"
            fi
        done
    done
    
    echo ""
    echo "=========================================="
    echo "Summary"
    echo "=========================================="
    echo "Total jobs submitted: $job_count"
    echo "Output directory: $OUTPUT_DIR"
    echo ""
    echo "Monitor jobs with: squeue -u \$USER"
    echo "Check logs in: $OUTPUT_DIR/logs/"
    echo ""
    echo "Once all jobs complete, run:"
    echo "  bash $0 --output_dir=$OUTPUT_DIR --run_metrics=false --process_results=true"
    echo "=========================================="
fi

if [[ "$PROCESS_RESULTS" == "true" ]]; then
    process_all_results
fi

echo ""
echo "Done!"
