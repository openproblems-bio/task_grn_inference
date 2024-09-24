#!/bin/bash

# RUN_ID="run_$(date +%Y-%m-%d_%H-%M-%S)"
RUN_ID="scglue"
resources_dir="s3://openproblems-data/resources/grn"
publish_dir="s3://openproblems-data/resources/grn/results/${RUN_ID}"


param_file="./params/${RUN_ID}.yaml"

cat > $param_file << HERE
param_list:
  - id: ${RUN_ID}
    multiomics_rna: ${resources_dir}/grn-benchmark/multiomics_rna.h5ad
    multiomics_atac: ${resources_dir}/grn-benchmark/multiomics_atac.h5ad
    annotation_file: ${resources_dir}/supplementary/gencode.v45.annotation.gtf.gz
    motif_file: ${resources_dir}/supplementary/JASPAR2022-hg38.bed.gz
    num_workers: $num_workers
    temp_dir: ./tmp/grn
output_state: "state.yaml"
publish_dir: "$publish_dir"
HERE



# ./tw-windows-x86_64.exe launch `
#        https://github.com/openproblems-bio/task_grn_inference.git `
#        --revision build/main `
#        --pull-latest `
#        --main-script target/nextflow/workflows/grn_inference_scglue/main.nf `
#        --workspace 53907369739130 `
#        --compute-env 6TeIFgV5OY4pJCk8I0bfOh `
#        --params-file ./params/scglue.yaml `
#        --config src/common/nextflow_helpers/labels_tw.config