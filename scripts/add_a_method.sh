#!/bin/bash

echo "This script is not supposed to be run directly."
echo "Please run the script step-by-step."
# exit 1

# sync resources
# scripts/download_resources.sh

# create a new component
method_id="dummy"
method_lang="python" # change this to "r" if need be


bash common/scripts/create_component --type method --language ${method_lang} --name ${method_id}

# TODO: fill in required fields in src/methods/foo/config.vsh.yaml
# TODO: edit src/methods/foo/script.py/R

# test the component
viash test src/methods/$method_id/config.vsh.yaml

# rebuild the container (only if you change something to the docker platform)
# You can reduce the memory and cpu allotted to jobs in _viash.yaml by modifying .platforms[.type == "nextflow"].config.labels
viash run src/methods/$method_id/config.vsh.yaml -- \
  ---setup cachedbuild ---verbose

# run the inference using the method for op dataset using only RNA data. Add more aurguments if needed.
viash run src/methods/$method_id/config.vsh.yaml -- \
  --rna "resources/inference_datasets/op_rna.h5ad" \
  --prediction "output/prediction.h5ad"

# run evaluation metrics
bash scripts/calculate_score.sh output/prediction.h5ad op

# print the score
python -c 'import pandas as ad; print(ad.read_h5ad("output/score.h5ad"))'