#!/bin/bash

set -e

# Build all components in a namespace (refer https://viash.io/reference/cli/ns_build.html)
viash ns build --parallel
