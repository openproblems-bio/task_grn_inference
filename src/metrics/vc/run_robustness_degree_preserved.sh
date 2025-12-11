#!/bin/bash
#SBATCH --job-name=vc_deg_pres_robust
#SBATCH --output=logs/degree_preserved_robustness_%j.out
#SBATCH --error=logs/degree_preserved_robustness_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=12:00:00
#SBATCH --mem=128GB
#SBATCH --partition=cpu

set -euo pipefail

# Configuration
DATASET="op"
METHOD="pearson_corr"
OUTPUT_DIR="output/vc/robustness_degree_preserved"

# Input files
PREDICTION="resources/results/${DATASET}/${DATASET}.${METHOD}.${METHOD}.prediction.h5ad"
EVALUATION_DATA="resources/grn_benchmark/evaluation_data/${DATASET}_bulk.h5ad"

# Create output directory
mkdir -p "${OUTPUT_DIR}"
mkdir -p logs

echo "Starting VC Robustness Analysis - Degree-Preserved Permutations"
echo "======================================"
echo "Dataset: ${DATASET}"
echo "Method: ${METHOD}"
echo "Prediction: ${PREDICTION}"
echo "Evaluation data: ${EVALUATION_DATA}"
echo "Output directory: ${OUTPUT_DIR}"
echo ""

# Activate environment
source ~/.bash_profile
conda activate py10

# Run degree-preserved robustness analysis
python src/metrics/vc/analyze_robustness_degree_preserved.py \
    --prediction "${PREDICTION}" \
    --evaluation_data "${EVALUATION_DATA}" \
    --output_dir "${OUTPUT_DIR}" \
    --intensities 0.0 0.1 0.2 0.3 0.5 0.7 1.0 \
    --n_top_genes 3000

echo ""
echo "Analysis complete!"
echo "Results saved to: ${OUTPUT_DIR}"
