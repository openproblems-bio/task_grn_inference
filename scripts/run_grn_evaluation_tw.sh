#!/bin/bash

# RUN_ID="run_$(date +%Y-%m-%d_%H-%M-%S)"
submit=false
read_results=false

RUN_ID="subsample_200_ridge"
resources_dir="s3://openproblems-data/resources/grn"
# resources_dir="resources/"
publish_dir="s3://openproblems-data/resources/grn/results/${RUN_ID}"
reg_type=ridge
subsample=200
max_workers=20

param_file="./params/${RUN_ID}.yaml"


grn_names=(
    "celloracle"
    "scenicplus"
    "figr"
    "granie"
    "scglue"
)
layers=("pearson" "lognorm" "scgen_pearson" "scgen_lognorm" "seurat_pearson" "seurat_lognorm")

# Start writing to the YAML file
cat > $param_file << HERE
param_list:
HERE

append_entry() {
  cat >> $param_file << HERE
  - id: ${layer}_${1}
    perturbation_data: ${resources_dir}/grn-benchmark/perturbation_data.h5ad
    layer: ${layer}
    reg_type: $reg_type
    method_id: $1
    subsample: $subsample
    max_workers: $max_workers
    consensus: ${resources_dir}/prior/consensus-num-regulators.json
    ${2:+tf_all: ${resources_dir}/prior/tf_all.csv}
    ${3:+prediction: ${resources_dir}/grn_models/$1.csv}
HERE
}
# Loop through grn_names and layers
for grn_name in "${grn_names[@]}"; do
  for layer in "${layers[@]}"; do
    append_entry "$grn_name" "" "true"
  done
done

# Append negative control
grn_name="negative_control"
append_entry "$grn_name"

# Append positive controls
grn_name="positive_control"
for layer in "${layers[@]}"; do
  append_entry "$grn_name" "true"
done

# Append the remaining output_state and publish_dir to the YAML file
cat >> $param_file << HERE
output_state: "state.yaml"
publish_dir: "$publish_dir"
HERE

# nextflow run . \
#   -main-script  target/nextflow/workflows/run_grn_evaluation/main.nf \
#   -profile docker \
#   -with-trace \
#   -c src/common/nextflow_helpers/labels_ci.config \
#   -params-file ${param_file}

./tw-windows-x86_64.exe launch `
    https://github.com/openproblems-bio/task_grn_benchmark.git `
    --revision build/main `
    --pull-latest `
    --main-script target/nextflow/workflows/run_grn_evaluation/main.nf `
    --workspace 53907369739130 `
    --compute-env 6TeIFgV5OY4pJCk8I0bfOh `
    --params-file ./params/subsample_200_ridge.yaml `
    --config src/common/nextflow_helpers/labels_tw.config


