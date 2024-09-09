#!/bin/bash

# RUN_ID="run_$(date +%Y-%m-%d_%H-%M-%S)"
RUN_ID="single_omics_scgpt"
# resources_dir="./resources_test/"
resources_dir="s3://openproblems-data/resources_test/grn"
publish_dir="${resources_dir}/results/${RUN_ID}"


reg_type=ridge
subsample=-2
max_workers=10
layer='pearson'

param_file="./params/${RUN_ID}.yaml"

# Start writing to the YAML file
cat > $param_file << HERE
param_list:
  - id: ${reg_type}
    perturbation_data: ${resources_dir}/grn-benchmark/perturbation_data.h5ad
    multiomics_rna: ${resources_dir}/grn-benchmark/multiomics_rna.h5ad
    reg_type: $reg_type
    subsample: $subsample
    max_workers: $max_workers
    layer: $layer
    consensus: ${resources_dir}/prior/consensus-num-regulators.json
    tf_all: ${resources_dir}/prior/tf_all.csv
    model_file: ${resources_dir}/supplementary/finetuned_scGPT_adamson/best_model.pt
    model_config_file: ${resources_dir}/supplementary/finetuned_scGPT_adamson/args.json
    vocab_file: ${resources_dir}/supplementary/finetuned_scGPT_adamson/vocab.json

output_state: "state.yaml"
publish_dir: "$publish_dir"
HERE

# nextflow run . \
#   -main-script  target/nextflow/workflows/run_benchmark_single_omics/main.nf \
#   -profile docker \
#   -with-trace \
#   -c src/common/nextflow_helpers/labels_ci.config \
#   -params-file ${param_file}

# ./tw-windows-x86_64.exe launch `
#     https://github.com/openproblems-bio/task_grn_inference.git `
#     --revision build/main `
#     --pull-latest `
#     --main-script target/nextflow/workflows/run_benchmark_single_omics/main.nf `
#     --workspace 53907369739130 `
#     --compute-env 6TeIFgV5OY4pJCk8I0bfOh `
#     --params-file ./params/single_omics_scgpt.yaml `
#     --config src/common/nextflow_helpers/labels_tw.config


