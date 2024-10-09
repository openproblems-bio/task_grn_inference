#!/bin/bash

RUN_ID="scenicplus"
resources_dir="s3://openproblems-data/resources_test/grn"
publish_dir="s3://openproblems-data/resources_test/grn/results/${RUN_ID}"

# resources_dir="./resources_test"
# publish_dir="./output/${RUN_ID}"

num_workers=10

param_file="./params/${RUN_ID}.yaml"
# Start writing to the YAML file
cat > $param_file << HERE
param_list:
  - id: ${RUN_ID}
    multiomics_rna: ${resources_dir}/grn-benchmark/multiomics_rna.h5ad
    multiomics_atac: ${resources_dir}/grn-benchmark/multiomics_atac.h5ad
    num_workers: $num_workers
    temp_dir: tmp/grn
output_state: "state.yaml"
publish_dir: "$publish_dir"
HERE


# nextflow run . \
#   -main-script  target/nextflow/workflows/grn_inference_scenicplus/main.nf \
#   -profile docker \
#   -with-trace \
#   -c src/common/nextflow_helpers/labels_ci.config \
#   -params-file params/${RUN_ID}.yaml


# ./tw-windows-x86_64.exe launch `
#       https://github.com/openproblems-bio/task_grn_inference.git `
#       --revision build/main `
#       --pull-latest `
#       --main-script target/nextflow/workflows/grn_inference_scenicplus/main.nf `
#       --workspace 53907369739130 `
#       --compute-env 6TeIFgV5OY4pJCk8I0bfOh `
#       --params-file ./params/scenicplus.yaml `
#       --config src/common/nextflow_helpers/labels_tw.config
./tw launch https://github.com/openproblems-bio/task_grn_inference \
  --revision build/main \
  --pull-latest \
  --main-script target/nextflow/workflows/grn_inference_scenicplus/main.nf \
  --workspace 53907369739130 \
  --compute-env 6TeIFgV5OY4pJCk8I0bfOh \
  --params-file ${param_file} \
  --config src/common/nextflow_helpers/labels_tw.config

