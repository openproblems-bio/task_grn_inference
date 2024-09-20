#!/bin/bash

set -e

echo ">> Downloading resources"

viash run src/common/sync_test_resources/config.vsh.yaml -- \
  --input "s3://openproblems-data/resources/grn/grn-benchmark" \
  --output "resources/grn-benchmark" \
  --delete

viash run src/common/sync_test_resources/config.vsh.yaml -- \
  --input "s3://openproblems-data/resources/grn/prior" \
  --output "resources/prior" \
  --delete

viash run src/common/sync_test_resources/config.vsh.yaml -- \
  --input "s3://openproblems-data/resources/grn/grn_models" \
  --output "resources/grn_models" \
  --delete
echo ">> Downloading resources test"
viash run src/common/sync_test_resources/config.vsh.yaml -- \
  --input "s3://openproblems-data/resources_test/grn/grn-benchmark" \
  --output "resources_test/grn-benchmark" \
  --delete

viash run src/common/sync_test_resources/config.vsh.yaml -- \
  --input "s3://openproblems-data/resources_test/grn/prior" \
  --output "resources_test/prior" \
  --delete

viash run src/common/sync_test_resources/config.vsh.yaml -- \
  --input "s3://openproblems-data/resources_test/grn/grn_models" \
  --output "resources_tests/grn_models" \
  --delete


