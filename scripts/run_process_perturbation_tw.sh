#!/bin/bash

RUN_ID="run_$(date +%Y-%m-%d_%H-%M-%S)"
resources_dir="s3://openproblems-data/resources/grn/"
publish_dir="s3://openproblems-data/resources/grn/results/${RUN_ID}"

cat > /tmp/params.yaml << HERE
param_list:
  - id: test_process_perturatbion
    perturbation_counts: "$resources_dir/datasets_raw/perturbation_counts.h5ad",

output_state: "state.yaml"
publish_dir: "$publish_dir"
HERE

./tw-windows-x86_64.exe launch openproblems-bio/task_grn_benchmark \
  --revision build/main \
  --pull-latest \
  --main-script target/nextflow/workflows/process_perturbation/main.nf \
  --workspace 53907369739130 \
  --compute-env 6TeIFgV5OY4pJCk8I0bfOh \
  --params-file /tmp/params.yaml \
  --config src/common/nextflow_helpers/labels_tw.config



  ./tw-windows-x86_64.exe launch s3://openproblems-bio/task_grn_benchmark --revision build/main --pull-latest --main-script target/nextflow/workflows/process_perturbation/main.nf --workspace 53907369739130 --compute-env 6TeIFgV5OY4pJCk8I0bfOh --params-file /tmp/params.yaml --config src/common/nextflow_helpers/labels_tw.config
