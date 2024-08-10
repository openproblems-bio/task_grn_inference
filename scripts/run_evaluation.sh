#!/bin/bash

# Default values
grn=""
sample="200"  # Default value for sample
reg_type="ridge"
score="output/score.csv"

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
echo "Regression model: $reg_type"

# Clean bin/ folder
rm -r bin
mkdir bin

# Run regression analysis 1
echo "Running GRN benchmark with $grn and sample size $sample"
echo "Regression 1"
mkdir -p bin/regression_1
viash build src/metrics/regression_1/config.vsh.yaml -p docker -o bin/regression_1
bin/regression_1/regression_1 --perturbation_data resources/grn-benchmark/perturbation_data.h5ad --reg_type $reg_type --prediction $grn --score $score

# Run regression analysis 2
echo "Regression 2"
if [ ! -f resources/grn-benchmark/consensus-num-regulators.json ]; then
    viash build src/metrics/regression_2/consensus/config.vsh.yaml --platform docker -o bin/regression_2/consensus
    bin/regression_2/consensus/consensus_for_regression_2 --perturbation_data resources/grn-benchmark/perturbation_data.h5ad --output resources/grn-benchmark/consensus-num-regulators.json --grn_folder resources/grn-benchmark/grn_models/ --grns ananse.csv,celloracle.csv,figr.csv,granie.csv,scenicplus.csv,scglue.csv
fi
mkdir -p bin/regression_2
viash build src/metrics/regression_2/config.vsh.yaml -p docker -o bin/regression_2
bin/regression_2/regression_2 --perturbation_data resources/grn-benchmark/perturbation_data.h5ad --consensus resources/grn-benchmark/consensus-num-regulators.json --layer scgen_pearson --reg_type $reg_type --prediction $grn --score $score
