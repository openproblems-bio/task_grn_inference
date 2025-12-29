#!/bin/bash
#SBATCH --job-name=rc_tf_act_robustness
#SBATCH --output=logs/robustness_%j.out
#SBATCH --error=logs/robustness_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=10:00:00
#SBATCH --mem=128GB
#SBATCH --partition=cpu

set -euo pipefail

# Configuration
DATASET="op"
METHOD="grnboost"
OUTPUT_DIR="output/rc_tf_act/robustness"

# Input files
PREDICTION="resources/results/${DATASET}/${DATASET}.${METHOD}.${METHOD}.prediction.h5ad"
EVALUATION_DATA="resources/grn_benchmark/evaluation_data/${DATASET}_bulk.h5ad"

# Create output directory
mkdir -p "${OUTPUT_DIR}"

echo "Starting RC_TF_ACT Robustness Analysis"
echo "======================================"
echo "Dataset: ${DATASET}"
echo "Method: ${METHOD}"
echo "Prediction: ${PREDICTION}"
echo "Evaluation data: ${EVALUATION_DATA}"
echo "Output directory: ${OUTPUT_DIR}"
echo ""

# Run robustness analysis with different permutation types and degrees
python src/metrics/rc_tf_act/analyze_robustness.py \
    --prediction "${PREDICTION}" \
    --evaluation_data "${EVALUATION_DATA}" \
    --output_dir "${OUTPUT_DIR}" \
    --degrees 0.0 0.2 0.5 1.0 \
    --noise_types direction sign weight net \
    --layer lognorm \
    --min_targets 5

echo ""
echo "Analysis complete!"
echo "Results saved to: ${OUTPUT_DIR}"
