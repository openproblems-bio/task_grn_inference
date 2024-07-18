#!/bin/bash

set -e

echo ">> Uploading resources"

aws s3 sync  ./resources s3://openproblems-data/resources/grn/


