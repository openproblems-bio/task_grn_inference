#!/bin/bash

set -e

# common/scripts/sync_resources

aws s3 sync  s3://openproblems-data/resources_test/grn/grn_benchmark resources_test/grn_benchmark  --no-sign-request
