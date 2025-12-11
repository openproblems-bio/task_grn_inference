#!/bin/bash
#SBATCH --job-name=vc_robust_compare
#SBATCH --output=logs/robustness_comparison_%j.out
#SBATCH --error=logs/robustness_comparison_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=24:00:00
#SBATCH --mem=128GB
#SBATCH --partition=cpu

set -euo pipefail

# Configuration
DATASET="op"
METHOD="pearson_corr"
OUTPUT_DIR="output/vc/robustness_comprehensive"

# Input files
PREDICTION="resources/results/${DATASET}/${DATASET}.${METHOD}.${METHOD}.prediction.h5ad"
EVALUATION_DATA="resources/grn_benchmark/evaluation_data/${DATASET}_bulk.h5ad"

# Create output directory
mkdir -p "${OUTPUT_DIR}"
mkdir -p logs

echo "Starting Comprehensive VC Robustness Analysis"
echo "======================================"
echo "Dataset: ${DATASET}"
echo "Method: ${METHOD}"
echo "Prediction: ${PREDICTION}"
echo "Evaluation data: ${EVALUATION_DATA}"
echo "Output directory: ${OUTPUT_DIR}"
echo ""
echo "This will run:"
echo "  1. Standard permutations (net, sign, direction)"
echo "  2. Degree-preserved permutations"
echo "  3. Comparative analysis and visualizations"
echo ""

# Activate environment
source ~/.bash_profile
conda activate py10

# Run comprehensive robustness comparison
python src/metrics/vc/analyze_robustness_comparison.py \
    --prediction "${PREDICTION}" \
    --evaluation_data "${EVALUATION_DATA}" \
    --output_dir "${OUTPUT_DIR}" \
    --intensities 0.0 0.1 0.2 0.3 0.5 0.7 1.0 \
    --n_top_genes 3000

echo ""
echo "Comprehensive analysis complete!"
echo "Results saved to: ${OUTPUT_DIR}"
echo ""
echo "Output structure:"
echo "  ${OUTPUT_DIR}/standard/                    - Standard permutation results"
echo "  ${OUTPUT_DIR}/degree_preserved/            - Degree-preserved results"
echo "  ${OUTPUT_DIR}/comparison_all_permutations.png - Main comparison plot"
echo "  ${OUTPUT_DIR}/detailed_comparison.png      - Detailed subplots"
echo "  ${OUTPUT_DIR}/summary_report.txt           - Text summary"
