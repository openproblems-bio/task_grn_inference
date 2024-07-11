#!/bin/bash

echo "This script is not supposed to be run directly."
echo "Please run the script step-by-step."
# exit 1

# sync resources
# scripts/download_resources.sh

# create a new component
method_id="scglue"
method_lang="python" # change this to "r" if need be

viash run src/common/create_component/config.vsh.yaml -- \
  --language "$method_lang" \
  --name "$method_id"

# TODO: fill in required fields in src/methods/foo/config.vsh.yaml
# TODO: edit src/methods/foo/script.py/R

# test the component
viash test src/methods/$method_id/config.vsh.yaml

# rebuild the container (only if you change something to the docker platform)
# You can reduce the memory and cpu allotted to jobs in _viash.yaml by modifying .platforms[.type == "nextflow"].config.labels
viash run src/methods/$method_id/config.vsh.yaml -- \
  ---setup cachedbuild ---verbose

# run the method (using h5ad as input)
viash run src/methods/$method_id/config.vsh.yaml -- \
  --multiomics_rna "resources/grn-benchmark/multiomics_rna.h5ad" \
  --multiomics_atac "resources/grn-benchmark/multiomics_atac.h5ad" \
  --output "output/prediction.csv"

# run evaluation metric
viash run src/metrics/regression_1/config.vsh.yaml -- \
  --perturbation_h5ad "resources/grn-benchmark/perturbation.h5ad" \
  --prediction "output/prediction.csv" \
  --output "output/score.h5ad"

# print score on kaggle test dataset
python -c 'import anndata; print(anndata.read_h5ad("output/score.h5ad").uns)'