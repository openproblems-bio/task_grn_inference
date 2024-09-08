#!/bin/bash

set -e

[[ ! -d ../openproblems-v2 ]] && echo "You need to clone the openproblems repository next to this repository" && exit 1

../openproblems-v2/bin/create_task_readme \
  --task "grn_benchmark" \
  --task_dir "src" \
  --github_url "https://github.com/openproblems-bio/task_grn_inference/tree/main/" \
  --output "README.md"
