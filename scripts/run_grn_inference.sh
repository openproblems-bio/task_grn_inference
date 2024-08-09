#!/bin/bash

models_folder="resources/control_models"
layers=("pearson" "lognorm" "scgen_pearson" "scgen_lognorm" "seurat_pearson" "seurat_lognorm")

# run for negative control
prediction="$models_folder/negative_control.csv"
echo $prediction
viash run src/control_methods/negative_control/config.vsh.yaml -- --perturbation_data resources/grn-benchmark/perturbation_data.h5ad \
    --prediction $prediction

# run for positive control
for layer in "${layers[@]}"; do
  prediction="$models_folder/${layer}_positive_control.csv"
  echo $prediction
  viash run src/control_methods/positive_control/config.vsh.yaml -- --perturbation_data resources/grn-benchmark/perturbation_data.h5ad \
      --layer $layer \
      --prediction $prediction
done






RUN_ID="run_figr_$(date +%Y-%m-%d_%H-%M-%S)"
resources_dir="s3://openproblems-data/resources/grn"
publish_dir="s3://openproblems-data/resources/grn/results/${RUN_ID}"

# layers=("pearson","lognorm")
# predictions=("scenicplus","celloracle")
cat > /params/params_2.yaml << HERE
param_list:
  - id: pearson_scenicplus
    layer: "lognorm"
    prediction: resources/grn_models/scenicplus.csv 
publish_dir: "./"
output_state: "state.yaml"
HERE

# tw launch openproblems-bio/task_perturbation_prediction \
#   --revision build/main \
#   --pull-latest \
#   --main-script target/nextflow/workflows/run_benchmark/main.nf \
#   --workspace 53907369739130 \
#   --compute-env 6TeIFgV5OY4pJCk8I0bfOh \
#   --params-file /tmp/params.yaml \
#   --config src/common/nextflow_helpers/labels_tw.config
