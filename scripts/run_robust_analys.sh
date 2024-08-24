#!/bin/bash

# RUN_ID="run_$(date +%Y-%m-%d_%H-%M-%S)"

RUN_ID="robust_analy"
# resources_dir="resources"
# publish_dir="output/${RUN_ID}"

resources_dir="s3://openproblems-data/resources/grn"
publish_dir="s3://openproblems-data/resources/grn/results/${RUN_ID}"

grn_models_folder="${resources_dir}/grn_models"


reg_type=ridge
subsample=-2
max_workers=10
layer=pearson

param_file="./params/${RUN_ID}.yaml"


grn_names=(
    "collectri"
    "celloracle"
    "scenicplus"
    "figr"
    "granie"
    "scglue"
)

degrees=(10 20 50 100)
types=(links weight)

# Start writing to the YAML file
cat > $param_file << HERE
param_list:
HERE

append_entry() {
  cat >> $param_file << HERE
  - id: ${1}_${2}_${3}
    perturbation_data: ${resources_dir}/grn-benchmark/perturbation_data.h5ad
    layer: ${layer}
    reg_type: $reg_type
    method_id: $1
    subsample: $subsample
    max_workers: $max_workers
    consensus: ${resources_dir}/prior/consensus-num-regulators.json
    prediction: ${grn_models_folder}/$1.csv
    degree: ${3}
    type: ${2}
HERE
}
# Loop through grn_names and layers
for type in "${types[@]}"; do
    for degree in "${degrees[@]}"; do
        for grn_name in "${grn_names[@]}"; do
            append_entry "$grn_name" "$type" "$degree" 
        done
    done
done

# Append the remaining output_state and publish_dir to the YAML file
cat >> $param_file << HERE
output_state: "state.yaml"
publish_dir: "$publish_dir"
HERE

nextflow run . \
  -main-script  target/nextflow/workflows/run_robustness_analysis/main.nf \
  -profile docker \
  -with-trace \
  -c src/common/nextflow_helpers/labels_ci.config \
  -params-file ${param_file}

