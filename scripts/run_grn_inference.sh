#!/bin/bash

# RUN_ID="run_$(date +%Y-%m-%d_%H-%M-%S)"


RUN_ID="celloracle_test"
# resources_dir="s3://openproblems-data/resources_test/grn/"
resources_dir="./resources_test"
publish_dir="s3://openproblems-data/resources/grn/results/${RUN_ID}"
num_workers=20


# Start writing to the YAML file
cat > ./params/${RUN_ID}.yaml << HERE
param_list:
  - id: ${RUN_ID}
    multiomics_rna: ${resources_dir}/grn-benchmark/multiomics_rna.h5ad
    multiomics_atac: ${resources_dir}/grn-benchmark/multiomics_atac.h5ad
    num_workers: $num_workers
    temp_dir: ./tmp/celloracle
output_state: "state.yaml"
publish_dir: "$publish_dir"
HERE


nextflow run . \
  -main-script  target/nextflow/workflows/run_grn_inference/main.nf \
  -profile docker \
  -with-trace \
  -c src/common/nextflow_helpers/labels_ci.config \
  -params-file params/${RUN_ID}.yaml




# ./tw-windows-x86_64.exe launch https://github.com/openproblems-bio/task_grn_benchmark.git --revision build/main --pull-latest --main-script target/nextflow/workflows/run_grn_evaluation/main.nf --workspace 53907369739130 --compute-env 6TeIFgV5OY4pJCk8I0bfOh --params-file ./params/params.yaml --config src/common/nextflow_helpers/labels_tw.config