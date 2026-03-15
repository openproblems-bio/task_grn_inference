#!/bin/bash
# Local GRN Evaluation Script
# This script evaluates all GRN predictions for all datasets using local SLURM jobs
# 
# set -e
# Default parameters
RUN_CONSENSUS=false
RUN_METRICS=false
TEMP_DIR="output/evaluation"


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
        --run_consensus)
        RUN_CONSENSUS=true
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
        --temp_dir)        TEMP_DIR="$2"
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
echo "Number of workers: $NUM_WORKERS"
echo "=========================================="

# Create output directory
mkdir -p "$TEMP_DIR"

# Generate and source dataset configuration
echo "Generating dataset configuration..."
python src/utils/config.py 
source src/utils/config.env

# Get list of datasets from config
DATASETS=(${DATASETS//,/ })
METHODS=(${METHODS//,/ })

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

# Function to run consensus for a dataset
run_consensus() {
    local dataset=$1
    bash scripts/prior/run_consensus.sh --dataset "$dataset"
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
        for method in "${METHODS[@]}"; do
            prediction_file="${models_folder}/${dataset}.${method}.${method}.prediction.h5ad"
            
            if [[ -f "$prediction_file" ]]; then
                submit_metric_job "$dataset" "$method" "$prediction_file"
                ((job_count++))
            else
                echo "  [NOT FOUND] ${prediction_file}"
            fi
        done
    done
    
fi

echo ""
echo "Done!"
