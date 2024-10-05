#!/bin/bash

RUN_ID="process_multiomics"
# resources_dir="s3://openproblems-data/resources/grn/"
resources_dir="resources"
publish_dir="${resources_dir}/results/${RUN_ID}"

cat > ./params/${RUN_ID}.yaml << HERE
param_list:
  - id: process_multiomics
    multiome_counts: $resources_dir/datasets_raw/multiome_counts.h5ad

output_state: "state.yaml"
publish_dir: "$publish_dir"
HERE


# ./tw-windows-x86_64.exe launch  https://github.com/openproblems-bio/task_grn_inference.git `
#     --revision build/main --pull-latest `
#     --main-script target/nextflow/workflows/process_multiomics/main.nf `
#     --workspace 53907369739130 --compute-env 6TeIFgV5OY4pJCk8I0bfOh `
#     --params-file ./params/process_multiomics.yaml `
#     --config src/common/nextflow_helpers/labels_tw.config


nextflow run .   \
  -main-script  target/nextflow/workflows/process_multiomics/main.nf  \
  -profile docker     -with-trace     -c src/common/nextflow_helpers/labels_ci.config  \
  -params-file params/${RUN_ID}.yaml
