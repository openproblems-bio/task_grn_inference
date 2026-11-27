#!/bin/bash
# Run after: viash ns build --parallel --setup cachedbuild
# Requires resources_test/task_grn_inference (see scripts/create_resources/test_resources.sh)
set -e
cd "$(git rev-parse --show-toplevel)"

OUT=output/process_datasets_test
mkdir -p $OUT

IN=resources_test/task_grn_inference/inputs
nextflow run . \
  -main-script target/nextflow/workflows/process_datasets/main.nf \
  -profile docker \
  -c common/nextflow_helpers/labels_ci.config \
  --param_list $IN/param_list.yaml \
  --publish_dir $OUT
