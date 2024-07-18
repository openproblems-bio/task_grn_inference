#!/bin/bash

models_folder="resources/grn-benchmark/grn_models"
grn_models=(
    "celloracle"
    "scenicplus"
    "figr"
    "granie"
    "scglue"
)

# Loop through each GRN model and run the benchmark script
for grn_model in "${grn_models[@]}"; do
    bash scripts/run_evaluation_grn.sh --grn "$models_folder/$grn_model.csv" --sample 200 --reg_type ridge --score "out/$grn_model.csv"
done
