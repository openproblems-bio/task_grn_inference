#!/bin/bash

# Default values
grn_name="test"
subsample="200"  # Default value for sample
reg_type="ridge"
score="out/score.csv"
layer="pearson"
max_workers="4"

# Parse arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --grn_model) grn_model="$2"; shift ;;
        --grn_name) grn_name="$2"; shift ;;
        --subsample) subsample="$2"; shift ;;
        --reg_type) reg_type="$2"; shift ;;
        --score) score="$2"; shift ;;
        --layer) layer="$2"; shift ;;
        --max_workers) max_workers="$2"; shift ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done
# Ensure required arguments are provided
if [ -z "$grn_model"  ]; then
    echo "Usage: $0 --grn_model <grn_file> [--subsample <sample_value>]"
    exit 1
fi

# Print parsed arguments (for debugging purposes)
echo "GRN name: $grn_name   Regression model: $reg_type    Samples: $subsample     layer: $layer     max_workers: $max_workers"


# # Run regression analysis 1
echo "Regression 1"
score_reg1="${score/.csv/_reg1.csv}"
viash run src/metrics/regression_1/config.vsh.yaml -- --perturbation_data resources/grn-benchmark/perturbation_data.h5ad \
    --reg_type $reg_type --prediction $grn_model --score $score_reg1 --layer $layer --subsample $subsample --max_workers $max_workers

# # Run regression analysis 2
echo "Regression 2"
score_reg2="${score/.csv/_reg2.csv}"
time viash run src/metrics/regression_1/config.vsh.yaml -- --perturbation_data resources/grn-benchmark/perturbation_data.h5ad \
    --reg_type $reg_type --prediction $grn_model --score $score_reg2 --layer $layer --subsample $subsample --max_workers $max_workers


paste -d ',' "$score_reg1" <(cut -d ',' -f2- "$score_reg2") > $score
