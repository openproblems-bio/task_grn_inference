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






# RUN_ID="run_figr_$(date +%Y-%m-%d_%H-%M-%S)"
# resources_dir="s3://openproblems-data/resources/grn"
# publish_dir="s3://openproblems-data/resources/grn/results/${RUN_ID}"

# cat > /tmp/params.yaml << HERE
# param_list:
#   - id: neurips-2023-data
#     de_train_h5ad: "$resources_dir/neurips-2023-data/de_train.h5ad"
#     de_test_h5ad: "$resources_dir/neurips-2023-data/de_test.h5ad"
#     id_map: "$resources_dir/neurips-2023-data/id_map.csv"
#     layer: clipped_sign_log10_pval

# output_state: "state.yaml"
# publish_dir: "$publish_dir"
# HERE

# tw launch openproblems-bio/task_perturbation_prediction \
#   --revision build/main \
#   --pull-latest \
#   --main-script target/nextflow/workflows/run_benchmark/main.nf \
#   --workspace 53907369739130 \
#   --compute-env 6TeIFgV5OY4pJCk8I0bfOh \
#   --params-file /tmp/params.yaml \
#   --config src/common/nextflow_helpers/labels_tw.config
