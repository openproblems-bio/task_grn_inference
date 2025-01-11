#!/bin/bash

set -e

echo ">> Downloading resources"

viash run src/common/sync_test_resources/config.vsh.yaml -- \
  --input "s3://openproblems-data/resources/grn/inference_datasets/" \
  --output "resources/inference_datasets/" \
  --delete

viash run src/common/sync_test_resources/config.vsh.yaml -- \
  --input "s3://openproblems-data/resources/grn/evaluation_datasets/" \
  --output "resources/evaluation_datasets/" \
  --delete

viash run src/common/sync_test_resources/config.vsh.yaml -- \
  --input "s3://openproblems-data/resources/grn/prior" \
  --output "resources/prior" \
  --delete

viash run src/common/sync_test_resources/config.vsh.yaml -- \
  --input "s3://openproblems-data/resources/grn/grn_models/" \
  --output "resources/grn_models/" \
  --delete


