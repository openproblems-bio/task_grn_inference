#!/bin/bash

set -e

echo ">> Downloading resources"

viash run src/common/sync_test_resources/config.vsh.yaml -- \
  --input "s3://openproblems-data/resources/grn/" \
  --output "resources/" \
  --delete



