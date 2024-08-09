#!/bin/bash

# RUN_ID="run_$(date +%Y-%m-%d_%H-%M-%S)"
RUN_ID="subsample_200_ridge"
resources_dir="s3://openproblems-data/resources/grn/"
publish_dir="s3://openproblems-data/resources/grn/results/${RUN_ID}"
reg_type=ridge
subsample=200

grn_names=(
    "celloracle"
    "scenicplus"
    "figr"
    "granie"
    "scglue"
)
layers=("pearson" "lognorm" "scgen_pearson" "scgen_lognorm" "seurat_pearson" "seurat_lognorm")

# Start writing to the YAML file
cat > ./params/params_${RUN_ID}.yaml << HERE
param_list:
HERE

# Nested loops to iterate over grn_names and layers
for grn_name in "${grn_names[@]}"; do
  for layer in "${layers[@]}"; do
    cat >> ./params/params_${RUN_ID}.yaml << HERE
  - id: ${layer}_${grn_name}
    perturbation_data: resources/grn-benchmark/perturbation_data.h5ad
    layer: ${layer}
    prediction: resources/grn_models/${grn_name}.csv 
    reg_type: $reg_type
    method_id: $grn_name
    subsample: $subsample

HERE
  done
done

# append the positive controls
grn_name="positive_control"
for layer in "${layers[@]}"; do
  cat >> ./params/params_${RUN_ID}.yaml << HERE
  - id: ${layer}_${grn_name}
    perturbation_data: resources/grn-benchmark/perturbation_data.h5ad
    layer: ${layer}
    prediction: resources/control_models/${layer}_${grn_name}.csv 
    reg_type: $reg_type
    method_id: $grn_name
    subsample: $subsample

HERE
done

# append negative control
grn_name="negative_control"
cat >> ./params/params_${RUN_ID}.yaml << HERE
  - id: ${layer}_${grn_name}
    perturbation_data: resources/grn-benchmark/perturbation_data.h5ad
    layer: ${layer}
    prediction: resources/control_models/${grn_name}.csv 
    reg_type: $reg_type
    method_id: $grn_name
    subsample: $subsample

HERE
# Append the remaining output_state and publish_dir to the YAML file
cat >> ./params/params_${RUN_ID}.yaml << HERE
output_state: "state.yaml"
publish_dir: "$publish_dir"
HERE




# ./tw-windows-x86_64.exe launch openproblems-bio/task_grn_benchmark \
#   --revision build/main \
#   --pull-latest \
#   --main-script target/nextflow/workflows/process_perturbation/main.nf \
#   --workspace 53907369739130 \
#   --compute-env 6TeIFgV5OY4pJCk8I0bfOh \
#   --params-file ./params/params_${RUN_ID}.yaml \
#   --config src/common/nextflow_helpers/labels_tw.config



#   ./tw-windows-x86_64.exe launch openproblems-bio/task_grn_benchmark --revision build/main --pull-latest --main-script target/nextflow/workflows/process_perturbation/main.nf --workspace 53907369739130 --compute-env 6TeIFgV5OY4pJCk8I0bfOh --params-file /tmp/params.yaml --config src/common/nextflow_helpers/labels_tw.config
