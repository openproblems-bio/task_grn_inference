#!/bin/bash
# Consensus Calculation Script
# This script runs consensus calculations for both Regression and WS distance metrics
# Usage: bash scripts/run_consensus.sh <dataset> [run_mode]
#   dataset: name of the dataset (e.g., replogle, op, norman)
#   run_mode: 'local' (default) or 'aws'

set -e

DATASET=$1

if [ -z "$DATASET" ]; then
  echo "Usage: bash scripts/run_consensus.sh <dataset> [run_mode]"
  echo "  dataset: name of the dataset (required)"
  echo "  run_mode: 'local' (default) or 'aws'"
  exit 1
fi

echo "=========================================="
echo "Running Consensus Calculation"
echo "Dataset: $DATASET"
echo "Run mode: $RUN_MODE"
echo "=========================================="

# Set paths based on run mode
resources_dir="./resources"
models_dir="${resources_dir}/results/$DATASET"

# Get available methods from config
echo "Checking available methods..."
available_methods=$(python -c "
from src.utils.config import METHODS
import os
methods = []
for method in METHODS:
    file = f'resources/results/$DATASET/$DATASET.{method}.{method}.prediction.h5ad'
    if os.path.exists(file):
        methods.append(method)
print(' '.join(methods))
")

if [ -z "$available_methods" ]; then
  echo "No prediction files found for dataset: $DATASET"
  exit 1
fi

echo "Available methods: $available_methods"

# Convert space-separated list to array
methods_array=($available_methods)

# Build predictions list
predictions=()
for method in "${methods_array[@]}"; do
    file="resources/results/${DATASET}/${DATASET}.${method}.${method}.prediction.h5ad"
    if [ -e "$file" ]; then
        predictions+=("$file")
    fi
done

if [ ${#predictions[@]} -eq 0 ]; then
  echo "No prediction files found for consensus calculation"
  exit 1
fi

echo "Found ${#predictions[@]} prediction files for consensus calculation"
printf '%s\n' "${predictions[@]}"

# Run Regression consensus
echo ""
echo "Running Regression consensus..."
python src/metrics/regression/consensus/script.py \
    --dataset "$DATASET" \
    --regulators_consensus "resources/grn_benchmark/prior/regulators_consensus_${DATASET}.json" \
    --evaluation_data "resources/grn_benchmark/evaluation_data/${DATASET}_bulk.h5ad" \
    --predictions "${predictions[@]}"

echo "Regression consensus completed successfully"

# Run WS distance consensus (only for applicable datasets)
applicable_datasets=("norman" "adamson" "replogle" "xaira_HEK293T" "xaira_HCT116")
skip_ws=true
for d in "${applicable_datasets[@]}"; do
    if [[ "$DATASET" == "$d" ]]; then
        skip_ws=false
        break
    fi
done

if [ "$skip_ws" = true ]; then
    echo ""
    echo "Skipping WS distance consensus (not applicable for dataset: $DATASET)"
else
    echo ""
    echo "Running WS distance consensus..."
    python src/metrics/ws_distance/consensus/script.py \
        --dataset "$DATASET" \
        --models_dir "resources/results/$DATASET" \
        --ws_consensus "resources/grn_benchmark/prior/ws_consensus_${DATASET}.csv" \
        --tf_all "resources/grn_benchmark/prior/tf_all.csv" \
        --evaluation_data_sc "resources/processed_data/${DATASET}_evaluation_sc.h5ad" \
        --models "${methods_array[@]}"
    
    echo "WS distance consensus completed successfully"
fi

# Sync results to AWS if needed
if [ "$RUN_MODE" = "aws" ]; then
    echo ""
    echo "Syncing consensus results to AWS..."
    aws s3 sync resources/grn_benchmark/prior s3://openproblems-data/resources/grn/grn_benchmark/prior
    echo "Sync completed"
fi

echo ""
echo "=========================================="
echo "Consensus calculation completed for $DATASET"
echo "=========================================="
