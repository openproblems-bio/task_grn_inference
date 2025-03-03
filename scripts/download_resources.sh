#!/bin/bash

set -e

common/scripts/sync_resources

# aws s3 sync  s3://openproblems-data/resources/grn/grn_benchmark resources/grn_benchmark --delete --no-sign-request


# aws s3 sync  s3://openproblems-data/resources/grn/grn_models resources/grn_models --delete --no-sign-request
