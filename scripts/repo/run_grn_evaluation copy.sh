#!/bin/bash

# RUN_ID="run_$(date +%Y-%m-%d_%H-%M-%S)"
# reg_type=${1} #GB, ridge
viash ns build --parallel
reg_type=ridge

RUN_ID="grn_evaluation_all_${reg_type}"
resources_dir="s3://openproblems-data/resources/grn"
# resources_dir="./resources"
publish_dir="${resources_dir}/results/${RUN_ID}"
grn_models_folder="${resources_dir}/grn_models"

subsample=-2
max_workers=10
layer=scgen_pearson
metric_ids="[regression_1, regression]"

param_file="./params/${RUN_ID}.yaml"

grn_names=(
    "scglue"
    "scenicplus"
    "celloracle"
    "granie"
    "figr"
    "collectri"
    "genie3"
    "grnboost2"
    "ppcor"
    "portia"
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

append_entry_control() {
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
    causal: ${2}
    corr_method: ${3}
    prediction: ${resources_dir}/grn_models/collectri.csv
    cell_type_specific:  ${4}
    metacell:  ${5}
    impute: ${6}
HERE

}

Loop through grn_names and layers
for grn_name in "${grn_names[@]}"; do
  append_entry "$grn_name" 
done

## controls
append_entry_control "negative_control" "" "" "false" "false" "false"
append_entry_control "positive_control" "" "" "false" "false" "false"
append_entry_control "baseline_pearson" "false" "pearson" "false" "false" "false"
append_entry_control "baseline_dotproduct" "false" "dotproduct" "false" "false" "false"
append_entry_control "baseline_dotproduct_causal" "true" "dotproduct" "false" "false" "false"
append_entry_control "baseline_dotproduct_causal_cell_type" "true" "dotproduct" "true" "false" "false"
append_entry_control "baseline_dotproduct_causal_metacell" "true" "dotproduct" "false" "true" "false"
append_entry_control "baseline_dotproduct_causal_impute" "true" "dotproduct" "false" "false" "true"
append_entry_control "baseline_corr_causal_spearman" "true" "spearman"


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


