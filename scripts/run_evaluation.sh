#!/bin/bash

# Default values
grn=""
sample="200"  # Default value for sample
reg_type="ridge"
score="out/score.csv"

# Parse arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --grn) grn="$2"; shift ;;
        --sample) sample="$2"; shift ;;
        --reg_type) reg_type="$2"; shift ;;
        --score) score="$2"; shift ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

# Ensure required arguments are provided
if [ -z "$grn" ]; then
    echo "Usage: $0 --grn <grn_file> [--sample <sample_value>]"
    exit 1
fi

# Print parsed arguments (for debugging purposes)
echo "GRN file: $grn"
echo "Sample value: $sample"
echo "Sample value: $reg_type"

# Example operation using the arguments
echo "Running GRN benchmark with $grn and sample size $sample"
echo "Regression 1"
viash build src/metrics/regression_1/config.vsh.yaml -p docker -o bin 
bin/regression_1 --perturbation_data resources/grn-benchmark/perturbation_data.h5ad --reg_type $reg_type --prediction $grn --score $score
