#!/bin/bash

set -e

# Test all components in a namespace (refer https://viash.io/reference/cli/ns_test.html)
viash ns test --parallel
