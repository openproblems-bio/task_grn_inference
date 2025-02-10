#!/bin/bash
#TODO: add atac to the inputs for op
dataset="norman"
RUN_ID=${dataset}_test
resources_dir="./resources/"
# resources_dir="s3://openproblems-data/resources/grn"
publish_dir="${resources_dir}/results/${RUN_ID}"

reg_type=ridge
num_workers=10
layer='X_norm'
metric_ids="[regression_1, regression_2, ws_distance]"

method_ids="[pearson_corr]"

param_file="./params/${RUN_ID}.yaml"

# Start writing to the YAML file
cat > $param_file << HERE
param_list:
  - id: ${reg_type}
    metric_ids: $metric_ids
    method_ids: $method_ids
    evaluation_data: ${resources_dir}/grn_benchmark/evaluation_data/${dataset}_bulk.h5ad
    evaluation_data_sc: ${resources_dir}/grn_benchmark/evaluation_data/${dataset}_sc.h5ad
    rna: ${resources_dir}/grn_benchmark/inference_data/${dataset}_rna.h5ad
    reg_type: $reg_type
    num_workers: $num_workers
    layer: $layer
    regulators_consensus: ${resources_dir}/grn_benchmark/prior/regulators_consensus_${dataset}.json
    ws_consensus: ${resources_dir}/grn_benchmark/prior/ws_consensus_${dataset}.csv
    ws_distance_background: ${resources_dir}/grn_benchmark/prior/ws_distance_background_${dataset}.csv
    prediction: ${grn_models_folder}/${dataset}/$1.h5ad
    tf_all: ${resources_dir}/grn_benchmark/prior/tf_all.csv

output_state: "state.yaml"
publish_dir: "$publish_dir"
HERE

nextflow run . \
  -main-script  target/nextflow/workflows/run_benchmark/main.nf \
  -profile docker \
  -with-trace \
  -c common/nextflow_helpers/labels_ci.config \
  -params-file ${param_file}

# ./tw-windows-x86_64.exe launch `
#     https://github.com/openproblems-bio/task_grn_inference.git `
#     --revision build/main `
#     --pull-latest `
#     --main-script target/nextflow/workflows/run_benchmark/main.nf `
#     --workspace 53907369739130 `
#     --compute-env 6TeIFgV5OY4pJCk8I0bfOh `
#     --params-file ./params/benchmark_donor_0_default.yaml `
#     --config src/common/nextflow_helpers/labels_tw.config

# tw launch https://github.com/openproblems-bio/task_grn_inference \
#   --revision build/main \
#   --pull-latest \
#   --main-script target/nextflow/workflows/run_benchmark/main.nf \
#   --workspace 53907369739130 \
#   --compute-env 7gRyww9YNGb0c6BUBtLhDP \
#   --params-file ${param_file} \
#   --config common/nextflow_helpers/labels_tw.config
