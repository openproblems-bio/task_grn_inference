#!/bin/bash

# RUN_ID="run_$(date +%Y-%m-%d_%H-%M-%S)"

subsamples=(-2 -3 -4)

RUN_ID="robust_analy_reg2_$1"
# resources_dir="resources"
resources_dir="s3://openproblems-data/resources/grn"

publish_dir="${resources_dir}/${RUN_ID}"

grn_models_folder="${resources_dir}/grn_models"


reg_type=ridge
max_workers=10
layer=scgen_pearson

param_file="./params/${RUN_ID}.yaml"

grn_names=(
    "collectri"
    "celloracle"
    "scenicplus"
    "figr"
    "granie"
    "scglue"
)


# Start writing to the YAML file
cat > $param_file << HERE
param_list:
HERE

append_entry() {
  cat >> $param_file << HERE
  - id: ${1}_${2}
    perturbation_data: ${resources_dir}/grn-benchmark/perturbation_data.h5ad
    layer: ${layer}
    reg_type: $reg_type
    method_id: ${2}-${1}
    subsample: $2
    max_workers: $max_workers
    consensus: ${resources_dir}/prior/consensus-num-regulators.json
    prediction: ${grn_models_folder}/$1.csv
    degree: 0

HERE
}
# Loop through grn_names and layers
for subsample in "${subsamples[@]}"; do
    for grn_name in "${grn_names[@]}"; do
        append_entry "$grn_name" "$subsample" 
    done
done



# Append the remaining output_state and publish_dir to the YAML file
cat >> $param_file << HERE
output_state: "state.yaml"
publish_dir: "$publish_dir"
HERE

# nextflow run . \
#   -main-script  target/nextflow/workflows/run_robustness_analysis/main.nf \
#   -profile docker \
#   -with-trace \
#   -c src/common/nextflow_helpers/labels_ci.config \
#   -params-file ${param_file}

# ./tw-windows-x86_64.exe launch `
#     https://github.com/openproblems-bio/task_grn_inference.git `
#     --revision build/main `
#     --pull-latest `
#     --main-script target/nextflow/workflows/run_grn_evaluation/main.nf `
#     --workspace 53907369739130 `
#     --compute-env 6TeIFgV5OY4pJCk8I0bfOh `
#     --params-file ./params/scgen_pearson_gb_pcs.yaml `
#     --config src/common/nextflow_helpers/labels_tw.config


