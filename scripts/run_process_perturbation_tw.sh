#!/bin/bash

RUN_ID="process_perturbation"
resources_dir="s3://openproblems-data/resources/grn/"
publish_dir="s3://openproblems-data/resources/grn/results/${RUN_ID}"

cat > ./params/${RUN_ID}.yaml << HERE
param_list:
  - id: test_process_perturatbion
    perturbation_counts: $resources_dir/datasets_raw/perturbation_counts.h5ad

output_state: "state.yaml"
publish_dir: "$publish_dir"
HERE


  ./tw-windows-x86_64.exe launch  https://github.com/openproblems-bio/task_grn_inference.git `
     --revision build/main --pull-latest `
     --main-script target/nextflow/workflows/process_perturbation/main.nf `
     --workspace 53907369739130 --compute-env 6TeIFgV5OY4pJCk8I0bfOh `
     --params-file ./params/process_perturbation.yaml `
     --config src/common/nextflow_helpers/labels_tw.config


nextflow run .   \
  -main-script  target/nextflow/workflows/grn_inference_granie/main.nf  \
  -profile docker     -with-trace     -c src/common/nextflow_helpers/labels_ci.config  \
  -params-file params/granie.yaml
