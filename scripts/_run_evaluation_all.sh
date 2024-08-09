#!/bin/bash

models_folder="resources/grn-benchmark/grn_models"
grn_names=(
    "celloracle"
    "scenicplus"
    "figr"
    "granie"
    "scglue"
)
layers=("pearson" "lognorm" "scgen_pearson" "scgen_lognorm" "seurat_pearson" "seurat_lognorm")


reg_type="ridge" 
subsample="2"
max_workers="4"
folder="out"


# layer="pearson"
layers=("pearson" "lognorm" "scgen_pearson" "scgen_lognorm" "seurat_pearson" "seurat_lognorm")

for layer in "${layers[@]}"; do
    Loop through each GRN model and run the benchmark script
    layer_folder="${folder}/${layer}"
    mkdir $layer_folder
    for grn_name in "${grn_names[@]}"; do
        grn_model="${models_folder}/${grn_name}.csv"
        score="${layer_folder}/${grn_name}.csv"
        bash scripts/run_evaluation.sh --grn_model "$grn_model" \
            --grn_name "$grn_name" --subsample "$subsample" --reg_type "$reg_type" \
            --score "$score" --max_workers "$max_workers" --layer "$layer"
    done

    for control_model in "${controls[@]}"; do
        prediction="predictions/{$layer}_{$control_model}.csv"
        viash run src/control_methods/$control_model/config.vsh.yaml -- --perturbation_data resources/grn-benchmark/perturbation_data.h5ad \
            --layer $layer \
            --prediction $prediction

        score="${layer_folder}/${control_model}.csv"
        bash scripts/run_evaluation.sh --grn_model "$prediction" \
            --grn_name "$grn_name" --subsample "$subsample" --reg_type "$reg_type" \
            --score "$score" --max_workers "$max_workers" --layer "$layer"
    done
done

# aggregate the scores across layers (axis=0)
for grn_name in "${grn_names[@]}"; do
    # Define a file to store the aggregated scores
    aggregated_score="${folder}/${grn_name}.csv"

    # Remove the file if it exists (to avoid appending to old data)
    rm -f "$aggregated_score"

    for layer in "${layers[@]}"; do
        # Define the path to the score file for the current layer
        layer_folder="${folder}/${layer}"
        score="${layer_folder}/${grn_name}.csv"

        # Append the score to the aggregated file
        # Skip the header for all but the first file
        if [ -s "$aggregated_score" ]; then
            tail -n +2 "$score" >> "$aggregated_score"
        else
            cat "$score" > "$aggregated_score"
        fi
    done
done