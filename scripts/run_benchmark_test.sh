#!/bin/bash

export NXF_VER=23.04.2

resources_dir="resources"
publish_dir="output/test_run_benchmark"

cat > /tmp/params.yaml << HERE
param_list:
  - id: test_run_1
    multiomics_rna: "$resources_dir/grn-benchmark/multiomics_rna.h5ad"
    multiomics_atac: "$resources_dir/grn-benchmark/multiomics_atac.h5ad"
    perturbation_data: "$resources_dir/grn-benchmark/perturbation_data.h5ad"
    layer: lognorm

output_state: "state.yaml"
publish_dir: "$publish_dir"
HERE

nextflow run . \
  -main-script target/nextflow/workflows/run_benchmark/main.nf \
  -profile docker \
  -resume \
  -params-file /tmp/params.yaml