#!/bin/bash

set -e

viash run src/common/create_task_readme/config.vsh.yaml -- \
  --task "grn_benchmark" \
  --task_dir "src" \
  --github_url "https://github.com/openproblems-bio/task_grn_inference/tree/main/" \
  --output "README.md"
