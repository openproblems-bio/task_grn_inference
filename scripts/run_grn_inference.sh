#!/bin/bash

# RUN_ID="run_$(date +%Y-%m-%d_%H-%M-%S)"


RUN_ID="celloracle_test"
resources_dir="s3://openproblems-data/resources/grn/"
publish_dir="s3://openproblems-data/resources/grn/results/${RUN_ID}"
reg_type=ridge
subsample=200
max_workers=20


# Start writing to the YAML file
cat > ./params/params_${RUN_ID}.yaml << HERE
param_list:
  - id: ${RUN_ID}
    multiomics_rna: resources/grn-benchmark/perturbation_data.h5ad
    layer: ${layer}
    prediction: resources/grn_models/${grn_name}.csv 
    reg_type: $reg_type
    method_id: $grn_name
    subsample: $subsample
    max_workers: $max_workers
output_state: "state.yaml"
publish_dir: "$publish_dir"
HERE


nextflow run . \
  -main-script  target/nextflow/workflows/run_grn_inference/main.nf \
  -profile docker \
  -with-trace \
  -c src/common/nextflow_helpers/labels_ci.config \
  -params-file params/params_${RUN_ID}.yaml




# ./tw-windows-x86_64.exe launch https://github.com/openproblems-bio/task_grn_benchmark.git --revision build/main --pull-latest --main-script target/nextflow/workflows/run_grn_evaluation/main.nf --workspace 53907369739130 --compute-env 6TeIFgV5OY4pJCk8I0bfOh --params-file ./params/params.yaml --config src/common/nextflow_helpers/labels_tw.config