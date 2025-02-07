#!/bin/bash

RUN_ID="scores_workflow"
echo $RUN_ID
# resources_dir="s3://openproblems-data/resources/grn/"
resources_dir="./resources/"
publish_dir="${resources_dir}/results/${RUN_ID}"

dataset="norman"
grn_models_folder="${resources_dir}/grn_models/" #TODO: change this

reg_type="ridge"
num_workers=10
metric_ids="[regression_1, regression_2, ws_distance]"

param_file="./params/${RUN_ID}.yaml"

# grn_names=(
#     "scglue"
#     "scenicplus"
#     "celloracle"
#     "granie"
#     "figr"
#     "collectri"
#     "genie3"
#     "grnboost2"
#     "ppcor"
#     "portia"
#     "positive_control"
#     "pearson_causal"
#     "pearson"
#     )
grn_names=(
    # "scglue"
    # "celloracle"

    "grnboost2"
    # "genie3"
    # "ppcor"
    # "scenic"
    # "portia"
    
    # "negative_control"
    # "positive_control"
    # "pearson_causal"
    "pearson_corr"

    # "collectri"
    )
echo $grn_names
# Start writing to the YAML file
cat > $param_file << HERE
param_list:
HERE

append_entry() {
  cat >> $param_file << HERE
  - id: ${reg_type}_${1}
    metric_ids: ${metric_ids}
    evaluation_data: ${resources_dir}/grn_benchmark/evaluation_data/${dataset}.h5ad
    evaluation_data_sc: ${resources_dir}/grn_benchmark/evaluation_data/${dataset}_sc.h5ad
    reg_type: $reg_type
    method_id: $1
    dataset_id: $dataset
    num_workers: $num_workers
    tf_all: ${resources_dir}/grn_benchmark/prior/tf_all.csv
    regulators_consensus: ${resources_dir}/grn_benchmark/prior/regulators_consensus_${dataset}.json
    ws_consensus: ${resources_dir}/grn_benchmark/prior/ws_consensus_${dataset}.csv
    ws_distance_background: ${resources_dir}/grn_benchmark/prior/ws_distance_background_${dataset}.csv
    prediction: ${grn_models_folder}/${dataset}/$1.h5ad
    layer: "X_norm"
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
#   -c common/nextflow_helpers/labels_ci.config \
#   -params-file ${param_file}

nextflow run openproblems-bio/task_grn_inference -r build/main \
  -latest \
  -main-script  target/nextflow/workflows/run_grn_evaluation/main.nf \
  -profile singularity \
  -with-trace \
  -c common/nextflow_helpers/labels_ci.config \
  -params-file ${param_file}


# ./tw-windows-x86_64.exe launch `
#     https://github.com/openproblems-bio/task_grn_inference.git `
#     --revision build/main `
#     --pull-latest `
#     --main-script target/nextflow/workflows/run_grn_evaluation/main.nf `
#     --workspace 53907369739130 `
#     --compute-env 6TeIFgV5OY4pJCk8I0bfOh `
#     --params-file ./params/grn_evaluation_so_ridge.yaml `
#     --config src/common/nextflow_helpers/labels_tw.config


# ./tw launch https://github.com/openproblems-bio/task_grn_inference \
#   --revision build/main \
#   --pull-latest \
#   --main-script target/nextflow/workflows/run_grn_evaluation/main.nf \
#   --workspace 53907369739130 \
#   --compute-env 6TeIFgV5OY4pJCk8I0bfOh \
#   --params-file ${param_file} \
#   --config src/common/nextflow_helpers/labels_tw.config
