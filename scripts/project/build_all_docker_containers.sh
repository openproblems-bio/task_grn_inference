#!/bin/bash

set -e

# Build all components in a namespace (refer https://viash.io/reference/cli/ns_build.html)
# and set up the container via a cached build
viash ns build --parallel --setup cachedbuild
