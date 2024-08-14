#!/bin/bash

# RUN_ID="run_$(date +%Y-%m-%d_%H-%M-%S)"
submit=false
read_results=true

RUN_ID="scglue"
resources_dir="s3://openproblems-data/resources/grn/"
publish_dir="s3://openproblems-data/resources/grn/results/${RUN_ID}"

num_workers=20

param_file="./params/${RUN_ID}.yaml"

cat > $param_file << HERE
param_list:
  - id: ${RUN_ID}
    multiomics_rna: ${resources_dir}/grn-benchmark/multiomics_rna.h5ad
    multiomics_atac: ${resources_dir}/grn-benchmark/multiomics_atac.h5ad
    annotation_file: ${resources_dir}/supplements/gencode.v45.annotation.gtf.gz
    motif_file: ${resources_dir}/supplements/JASPAR2022-hg38.bed.gz
    num_workers: $num_workers
    temp_dir: ./tmp/grn
output_state: "state.yaml"
publish_dir: "$publish_dir"
HERE

