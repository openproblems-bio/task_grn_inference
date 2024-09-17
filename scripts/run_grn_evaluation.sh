#!/bin/bash

# RUN_ID="run_$(date +%Y-%m-%d_%H-%M-%S)"
# reg_type=${1} #GB, ridge
# viash ns build --parallel
reg_type=ridge

RUN_ID="celltype_donor_0_subset_${reg_type}"
resources_dir="s3://openproblems-data/resources/grn"
# resources_dir="./resources"
publish_dir="${resources_dir}/results/${RUN_ID}"
grn_models_folder="${resources_dir}/grn_models/donor_0_celltype"

subsample=-2
max_workers=10
layer=scgen_pearson
metric_ids="[regression_1, regression_2]"

param_file="./params/${RUN_ID}.yaml"

grn_names=(
    # "scglue"
    # "scenicplus"
    # "celloracle"
    # "granie"
    # "figr"
    # "collectri"
    # "genie3"
    "grnboost2"
    # "ppcor"
    "portia"
    "positive_control"
    "pearson_causal"
    )

# Start writing to the YAML file
cat > $param_file << HERE
param_list:
HERE

append_entry() {
  cat >> $param_file << HERE
  - id: ${reg_type}_${1}
    metric_ids: ${metric_ids}
    perturbation_data: ${resources_dir}/grn-benchmark/perturbation_data.h5ad
    multiomics_rna: ${resources_dir}/grn-benchmark/multiomics_rna.h5ad
    reg_type: $reg_type
    method_id: $1
    subsample: $subsample
    max_workers: $max_workers
    tf_all: ${resources_dir}/prior/tf_all.csv
    layer: ${layer}
    consensus: ${resources_dir}/prior/consensus-num-regulators.json
    prediction: ${grn_models_folder}/$1.csv
HERE
}

# Loop through grn_names and layers
for grn_name in "${grn_names[@]}"; do
  append_entry "$grn_name"  
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
# subl resources/results/grn_evaluation_all_ridge/scores.yaml

# ./tw-windows-x86_64.exe launch `
#     https://github.com/openproblems-bio/task_grn_inference.git `
#     --revision build/main `
#     --pull-latest `
#     --main-script target/nextflow/workflows/run_grn_evaluation/main.nf `
#     --workspace 53907369739130 `
#     --compute-env 6TeIFgV5OY4pJCk8I0bfOh `
#     --params-file ./params/grn_evaluation_so_ridge.yaml `
#     --config src/common/nextflow_helpers/labels_tw.config


./tw launch https://github.com/openproblems-bio/task_grn_inference \
  --revision build/main \
  --pull-latest \
  --main-script target/nextflow/workflows/run_grn_evaluation/main.nf \
  --workspace 53907369739130 \
  --compute-env 6TeIFgV5OY4pJCk8I0bfOh \
  --params-file ${param_file} \
  --config src/common/nextflow_helpers/labels_tw.config
