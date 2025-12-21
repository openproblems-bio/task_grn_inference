#!/bin/bash
# Local GRN Evaluation Script
# This script evaluates all GRN predictions for all datasets using local SLURM jobs
# 
# set -e
# Default parameters
RUN_CONSENSUS=false
RUN_METRICS=false
PROCESS_RESULTS=false
TEMP_DIR="output/evaluation"
RESULTS_FILE="resources/results/all_scores.csv"


LAYER="lognorm"
NUM_WORKERS=20

# get the arguments
while [[ $# -gt 0 ]]; do
    key="$1"

    case $key in
        --run_consensus)
        RUN_CONSENSUS=true
        shift
        ;;
        --no_run_consensus)
        RUN_CONSENSUS=false
        shift
        ;;
        --run_metrics)
        RUN_METRICS=true
        shift
        ;;
        --no_run_metrics)
        RUN_METRICS=false
        shift
        ;;
        --process_results)
        PROCESS_RESULTS=true
        shift
        ;;
        --no_process_results)
        PROCESS_RESULTS=false
        shift
        ;;
        --temp_dir)
        TEMP_DIR="$2"
        shift
        shift
        ;;
        --results_file)
        RESULTS_FILE="$2"
        shift
        shift
        ;;
        --num_workers)
        NUM_WORKERS="$2"
        shift
        shift
        ;;
        *)
        echo "Unknown option: $key"
        exit 1
        ;;
    esac
done


echo "=========================================="
echo "GRN Evaluation Configuration"
echo "=========================================="
echo "Run consensus: $RUN_CONSENSUS"
echo "Run metrics: $RUN_METRICS"
echo "Output directory: $TEMP_DIR"
echo "Results file: $RESULTS_FILE"
echo "Process results: $PROCESS_RESULTS"
echo "Number of workers: $NUM_WORKERS"
echo "=========================================="

# Create output directory
mkdir -p "$TEMP_DIR"

# Generate and source dataset configuration
echo "Generating dataset configuration..."
python src/utils/config.py --output src/utils/dataset_config.env
source src/utils/dataset_config.env

# Get list of datasets from config
DATASETS=($(python -c "from src.utils.config import DATASET_GROUPS; print(' '.join(DATASET_GROUPS.keys()))"))
# DATASETS=('norman')

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
    local score_file="${TEMP_DIR}/${dataset}_${method}_score.h5ad"
    
    # Skip if score file already exists
    if [[ -f "$score_file" ]]; then
        echo "  Skipping ${job_name} - score file already exists"
        return
    fi
    
    echo "  Submitting job: ${job_name}"
    
    
    sbatch \
        --job-name="${job_name}" \
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
    
    local results_file="${RESULTS_FILE}"
    
    # Get list of datasets for parsing
    local datasets_list=$(python -c "from src.utils.config import DATASET_GROUPS; print(','.join(DATASET_GROUPS.keys()))")
    
    # Create Python script to aggregate results
    python - "$TEMP_DIR" "$RESULTS_FILE" "$datasets_list" << 'EOF'
import sys
import os
import pandas as pd
import anndata as ad
from pathlib import Path

temp_dir = sys.argv[1]
results_file = sys.argv[2]
datasets_str = sys.argv[3]
scores_dir = Path(temp_dir)

# Parse datasets list
known_datasets = datasets_str.split(',')

print(f"Looking for score files in: {scores_dir}")
print(f"Will save results to: {results_file}")
print(f"Known datasets: {known_datasets[:5]}...")

all_results = []

# Iterate through all score files
for score_file in scores_dir.glob("*_score.h5ad"):
    try:
        # Parse filename: {dataset}_{method}_score.h5ad
        filename = score_file.stem  # Remove .h5ad
        
        # Remove _score suffix
        if filename.endswith('_score'):
            filename = filename[:-6]  # Remove last 6 characters '_score'
        
        # Try to match against known datasets - find which dataset this file belongs to
        dataset = None
        method = None
        
        for ds in known_datasets:
            if filename.startswith(ds + '_'):
                dataset = ds
                method = filename[len(ds) + 1:]  # Everything after dataset_
                break
        
        if not dataset or not method:
            print(f"Warning: Could not parse filename {score_file.name}")
            continue
        
        # Load score file
        adata = ad.read_h5ad(score_file)
        
        # Extract metric scores from uns
        if 'metric_ids' in adata.uns and 'metric_values' in adata.uns:
            metric_ids = adata.uns['metric_ids']
            metric_values = adata.uns['metric_values']
            
            # Create row for each metric
            for metric_id, metric_value in zip(metric_ids, metric_values):
                # Convert metric_value to float if it's a string
                try:
                    score_value = float(metric_value) if isinstance(metric_value, str) else metric_value
                except (ValueError, TypeError):
                    score_value = metric_value
                
                all_results.append({
                    'dataset': dataset,
                    'method': method,
                    'metric': metric_id,
                    'score': score_value
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
    
    # Save results (create directory if needed)
    results_path = Path(results_file)
    results_path.parent.mkdir(parents=True, exist_ok=True)
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
}

# Function to run consensus for a dataset
run_consensus() {
    local dataset=$1
    
    echo ""
    echo "=========================================="
    echo "Running Consensus for Dataset: $dataset"
    echo "=========================================="
    
    # Build list of available predictions
    local models_dir="resources/results/${dataset}"
    local predictions=()
    
    for method in "${GRN_METHODS[@]}"; do
        local file="${models_dir}/${dataset}.${method}.${method}.prediction.h5ad"
        if [[ -f "$file" ]]; then
            predictions+=("$file")
            echo "  Found: ${method}"
        fi
    done
    
    if [[ ${#predictions[@]} -eq 0 ]]; then
        echo "  [WARNING] No prediction files found for ${dataset}"
        return
    fi
    
    echo "  Total predictions: ${#predictions[@]}"
    
    # Run Regression consensus
    echo ""
    echo "Running Regression consensus..."
    python src/metrics/regression/consensus/script.py \
        --dataset "$dataset" \
        --regulators_consensus "resources/grn_benchmark/prior/regulators_consensus_${dataset}.json" \
        --evaluation_data "resources/grn_benchmark/evaluation_data/${dataset}_bulk.h5ad" \
        --predictions "${predictions[@]}"
    
    # Run WS Distance consensus (only for applicable datasets)
    local applicable_datasets=("norman" "adamson" "replogle" "xaira_HEK293T" "xaira_HCT116")
    local skip=true
    
    for d in "${applicable_datasets[@]}"; do
        if [[ "$dataset" == "$d" ]]; then
            skip=false
            break
        fi
    done
    
    if $skip; then
        echo ""
        echo "Skipping WS Distance consensus for ${dataset} (not applicable)"
    else
        echo ""
        echo "Running WS Distance consensus..."
        
        # Extract model names for ws consensus
        local models=()
        for method in "${GRN_METHODS[@]}"; do
            local file="${models_dir}/${dataset}.${method}.${method}.prediction.h5ad"
            if [[ -f "$file" ]]; then
                models+=("$method")
            fi
        done
        
        python src/metrics/ws_distance/consensus/script.py \
            --dataset "$dataset" \
            --models_dir "$models_dir" \
            --ws_consensus "resources/grn_benchmark/prior/ws_consensus_${dataset}.csv" \
            --tf_all "resources/grn_benchmark/prior/tf_all.csv" \
            --evaluation_data_sc "resources/processed_data/${dataset}_evaluation_sc.h5ad" \
            --models "${models[@]}"
    fi
    
    echo ""
    echo "Consensus completed for ${dataset}"
}

# Main execution
if [[ "$RUN_CONSENSUS" == "true" ]]; then
    echo ""
    echo "=========================================="
    echo "Running Consensus for All Datasets"
    echo "=========================================="
    
    for dataset in "${DATASETS[@]}"; do
        run_consensus "$dataset"
    done
fi

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
    
fi

if [[ "$PROCESS_RESULTS" == "true" ]]; then
    process_all_results
fi

echo ""
echo "Done!"
